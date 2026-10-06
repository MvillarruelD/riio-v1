"""
ligand_regulon.py -- the INPUT side (cognate ligand) + OUTPUT side (regulon) for the report, plus the
ligand-guided expansion: search the genome for ligand-associated genes whose promoters are ADDITIONAL
intergenic regions to assess for operators (the fix for non-autoregulatory TFs whose operator sits at a
distant regulated gene).

  infer_ligand   -- Ligify on the TF's operon neighbourhood -> ranked cognate ligand(s).
  ligand_genes   -- scan the genome annotation (GFF product names) for genes related to the ligand /
                    effector element (transport / efflux / reductase / the metal name) -> their upstream
                    intergenic windows = more candidate operator regions.
  regulon        -- operator PWM/counts -> genome-wide intergenic scan -> the operons the TF likely drives.

All best-effort and offline-capable (Ligify touches cached UniProt; regulon reuses motif_rescan).
"""
from __future__ import annotations

import re
from concurrent.futures import ProcessPoolExecutor

# effector element -> product-name keywords that flag a ligand-handling gene
_ELEMENT_KEYWORDS = {
    "ZN": ["zinc", "znt", "czc", "czr", "zur", "znu", "cad", "zia"],
    "CU": ["copper", "cop", "cue", "cus", "ctp", "multicopper"],
    "CO": ["cobalt", "nickel", "cnr", "nik", "rcn"],
    "NI": ["nickel", "nik", "rcn", "urease"],
    "FE": ["iron", "ferr", "fec", "feo", "sider"],
    "MN": ["manganese", "mnt", "mnth", "sit"],
    "AS": ["arsen", "ars", "aox", "arsenate", "arsenite"],
    "CD": ["cadmium", "cad", "czc"],
    "HG": ["mercur", "mer"],
}


def _element_key(effector) -> str | None:
    """Return a canonical element token only for an explicit ion/element name.

    Substring matching is unsafe here: ``CO`` occurs in "carbon monoxide" and ``FE`` in many organic
    ligand names, which previously routed non-metal effectors through cobalt/iron gene brushes.
    """
    raw = re.sub(r"[^A-Za-z]", "", str(effector or "")).upper()
    aliases = {
        "ZN": "ZN", "ZINC": "ZN", "CU": "CU", "COPPER": "CU", "CO": "CO", "COBALT": "CO",
        "NI": "NI", "NICKEL": "NI", "FE": "FE", "IRON": "FE", "MN": "MN", "MANGANESE": "MN",
        "AS": "AS", "ARSENIC": "AS", "ARSENITE": "AS", "ARSENATE": "AS",
        "CD": "CD", "CADMIUM": "CD", "HG": "HG", "MERCURY": "HG", "MERCURIC": "HG",
    }
    return aliases.get(raw)


def infer_ligand(ctx, tf_gene, *, max_kb: float = 10.0):
    """Top cognate-ligand candidates from the TF operon (Ligify). Best-effort; [] on any failure.
    These are REPORTED, never used to rank operators -- see the discovery note in pipeline.run_novel."""
    from predictor.effector import ligify
    try:
        effs = ligify.effectors_from_context(ctx, tf_gene, max_kb=max_kb)
        return [{"ligand": e.ligand, "rank": e.rank, "enzyme": e.enzyme, "rhea": e.rhea,
                 "equation": e.equation} for e in effs[:5]]
    except Exception:
        return []


def ligand_genes(ctx, *, effector=None, ligand_names=None, max_genes: int = 12):
    """Genome genes whose product matches the ligand/effector element -> candidate regions. The keyword
    set is built from the AF3 ion / element, any inferred ligand names, AND the substrate-resolved
    transporter/efflux/enzyme keywords for that ion (`effector.substrate_map.keywords_for_effector`) --
    so a Zn sensor flags zntA/cadA/czcD/znuA, not just the bare element word. Each flagged gene is also
    run through `substrate_map.classify_gene`, attaching its substrate class + the ion it handles when the
    map recognises it (a confirmed regulated-gene, vs a bare keyword brush)."""
    from predictor.effector import substrate_map as sm

    kws = []
    el = _element_key(effector)
    is_metal = el is not None
    if el:
        kws += _ELEMENT_KEYWORDS.get(el, [])
    kws += sm.keywords_for_effector(effector)              # substrate-resolved gene/element stems
    # The inferred ligand NAME brush is noise for metal sensors: Ligify often returns a spurious organic
    # call (e.g. CueR -> 'L-glutamine', CzrA -> 'D-fructose 6-phosphate') whose words then flag unrelated
    # genes. Use it ONLY for non-metal effectors (where the organic ligand name is the real signal).
    if not is_metal:
        for nm in (ligand_names or []):
            kws += [w.lower() for w in re.findall(r"[A-Za-z]{4,}", nm)]
    kws = sorted(set(kws))
    if not kws:
        return [], kws
    # PRECISION: split keywords into short gene-name stems (cop/cue/cad/znt...) and descriptive words
    # (copper/cadmium/multicopper/'copper chaperone'...). Descriptive words may match the product free text;
    # short stems match the gene NAME only as a token/prefix -- never as a substring of the product, which
    # made 'cop'/'mer'/'cad' brush dozens of unrelated genes (polymerase, dimer, cascade, chaperone...).
    gene_stems = sorted({k for k in kws if " " not in k and len(k) <= 4})
    desc_words = sorted({k for k in kws if k not in gene_stems})

    def _matches(prod: str, name: str) -> bool:
        if any(w in prod or w in name for w in desc_words):
            return True
        name_toks = re.findall(r"[a-z0-9]+", name)
        return any(t == s or t.startswith(s) for s in gene_stems for t in name_toks)

    out = []
    for g in ctx.genes:
        prod = (g.product or "").lower()
        name = (g.name or "").lower()
        if _matches(prod, name):
            # upstream intergenic window (strand-aware) = a candidate operator region
            up = (g.start - 160, g.start) if g.strand == "+" else (g.end, g.end + 160)
            rec = {"gene": g.name or g.protein_id, "product": g.product, "start": g.start,
                   "end": g.end, "strand": g.strand, "upstream_region": [max(0, up[0]), up[1]]}
            hit = sm.classify_gene(g.product or "", g.name or "")
            if hit is not None:
                rec["substrate_effector"] = hit.effector       # the ion this handler's class is specific for
                rec["substrate_class"] = hit.substrate_class
                rec["substrate_specificity"] = hit.specificity
            out.append(rec)
    # Rank globally before truncating. Genome order is not biological evidence and previously caused the
    # first 12 lexical matches to hide later, substrate-specific transporter/enzyme matches.
    out.sort(key=lambda r: (-int(r.get("substrate_specificity", -1)),
                            str(r.get("gene") or "").lower(), int(r.get("start") or 0)))
    return out[:max_genes], kws


def _local_genes(ctx, site, *, flank: int = 4000, max_genes: int = 9):
    """The genes flanking an operator site (name/start/end/strand), for a genome-neighborhood figure of a
    target operator. Window is +/-`flank` bp around the site centre, capped at `max_genes`."""
    if not site:
        return []
    c = (int(site[0]) + int(site[1])) // 2
    near = [g for g in ctx.genes if g.end > c - flank and g.start < c + flank]
    near.sort(key=lambda g: g.start)
    if len(near) > max_genes:                      # keep the max_genes closest to the site centre
        near = sorted(near, key=lambda g: abs((g.start + g.end) // 2 - c))[:max_genes]
        near.sort(key=lambda g: g.start)
    return [{"name": g.name or g.protein_id, "start": g.start, "end": g.end, "strand": g.strand,
             "product": getattr(g, "product", "")} for g in near]


def regulon(ctx, counts, *, regulator=None, pvalue_thresh: float = 1e-4, top: int | None = None,
            with_local_genes: bool = True, stats: dict | None = None):
    """Operator counts -> ranked genome-wide intergenic candidates and their putative operons.
    Best-effort; these are not called significant unless their reported q-values pass an explicit gate.
    Each operon carries `local_genes` (the flanking gene context around its site) so the report can draw a
    genome-neighborhood diagram for the top-ranking natural operators.

    `top=None` (the default) writes the ranked candidates retained by `reconstruct_regulon`.
    `stats`, when given, is filled with `n_operons_found` / `n_operons_reported` / `truncated`: the
    returned list is capped, so its length alone cannot say whether the cap was reached."""
    if counts is None:
        return []
    from predictor.signals import regulon as regmod
    try:
        rg = regmod.reconstruct_regulon(ctx, counts=counts, regulator=regulator,
                                        pvalue_thresh=pvalue_thresh)
        if stats is not None:
            kept_n = len(rg.operons if top is None else rg.operons[:top])
            stats.update({"n_operons_found": rg.n_operons_found, "n_operons_reported": kept_n,
                          "truncated": rg.n_operons_found > kept_n, "n_sites": rg.n_sites})
        out = []
        for o in (rg.operons if top is None else rg.operons[:top]):
            site = [o.site_start, o.site_end]
            row = {"first_gene": o.first_gene, "strand": o.strand, "genes": o.genes,
                   "site": site, "qvalue": o.qvalue, "tier": o.tier}
            if with_local_genes:
                row["local_genes"] = _local_genes(ctx, site)
            out.append(row)
        return out
    except Exception as e:
        return [{"error": f"{type(e).__name__}: {e}"}]


# --------------------------------------------------------------------------- per-hit regulation (WS4)
def _regulated_gene(ctx, op_center, *, max_dist: int = 400):
    """The gene whose promoter the operator sits in: nearest gene 5'-end downstream of the operator
    (strand-aware). Returns (gene, operator->gene bp) or (None, None)."""
    best, bestd = None, max_dist + 1
    for g in ctx.genes:
        d = (g.start - op_center) if g.strand == "+" else (op_center - g.end)
        if 0 <= d <= max_dist and d < bestd:
            best, bestd = g, d
    return (best, bestd) if best else (None, None)


def _operon_gene_names(ctx, first_gene) -> list:
    """Every gene in the operon headed by `first_gene`, reusing the same walk the operator route uses
    (`signals.regulon.operon_of`) so both routes agree on what an operon is."""
    if first_gene is None:
        return []
    from predictor.signals import regulon as regmod
    try:
        members = regmod.operon_of(sorted(ctx.genes, key=lambda g: g.start), first_gene)
        return [g.name or g.protein_id for g in members if (g.name or g.protein_id)]
    except Exception:
        return [first_gene.name or first_gene.protein_id]


def _classify_mode_job(args):
    """Process-pool worker for one independent operator-neighbourhood calculation."""
    region, op_start, op_end = args
    try:
        from ..signals import promoter as pm
        return pm.classify_mode(region, op_start, op_end)
    except Exception:
        return None


def _classify_modes(jobs, *, max_workers: int = 6):
    """Classify independent operator windows in parallel, with a safe serial fallback."""
    if len(jobs) < 2:
        return [_classify_mode_job(job) for job in jobs]
    try:
        with ProcessPoolExecutor(max_workers=min(max_workers, len(jobs))) as executor:
            # executor.map preserves input order, so output tables remain deterministic.
            return list(executor.map(_classify_mode_job, jobs))
    except Exception:
        # Some embedded Python hosts cannot spawn child processes. Correctness still wins there; the
        # smaller operator-centred windows below keep the serial fallback bounded.
        return [_classify_mode_job(job) for job in jobs]


def per_hit_regulation(ctx, hits, *, window: int = 120, max_hits: int = 40):
    """Run the promoter/TSS + mode-of-regulation classifier on EVERY rescan operator hit (not just the
    autoregulatory one) and map each to the gene it regulates -- the genome-wide per-target regulation
    dataset (the redesign's WS4). Each record = {operator, regulated_gene, operator_to_gene_bp, mode, tss,
    spacer_len, mode_confidence, flags}. Sorted by hit score; best-effort per hit."""
    glen = len(ctx.sequence)
    ranked = sorted(hits, key=lambda h: -(getattr(h, "score", 0) or 0))[:max_hits]
    prepared = []
    for h in ranked:
        center = int(getattr(h, "dyad_center", None) or ((h.start + h.end) // 2))
        lo, hi = max(0, center - window), min(glen, center + window)
        op0, op1 = max(0, h.start - lo), min(hi - lo, h.end - lo)
        prepared.append((h, center, lo, op0, op1, ctx.sequence[lo:hi]))
    modes = _classify_modes([(region, op0, op1)
                             for _h, _center, _lo, op0, op1, region in prepared])
    out = []
    for (h, center, _lo, op0, op1, _region), mc in zip(prepared, modes):
        g, dist = _regulated_gene(ctx, center)
        _loc = getattr(ctx, "locate", None)
        _rep_id, _rep_start = (_loc(h.start) if _loc
                               else (getattr(ctx, "accession", ""), h.start))
        tss = (mc.promoter.get("tss") if (mc and mc.promoter) else None)
        # signed operator -> TSS distances IN THE REGION FRAME (op0/op1 and tss are region-relative):
        # positive = 3' (downstream) of the TSS, negative = upstream. The FULL FOOTPRINT is reported (both
        # edges + which promoter elements the operator physically covers), not only the centre -- so the
        # occlusion is described mechanistically ("covers -10 + TSS"), and the operator-vs-TSS figure can
        # draw the operator as a span rather than a point.
        fp = (mc.footprint if mc else {}) or {}
        op_to_tss = fp.get("op_center_to_tss")
        if op_to_tss is None and tss is not None:
            op_to_tss = ((op0 + op1) // 2) - tss
        out.append({
            # p/q come straight off the FIMO-equivalent scan (motifs/pwm_scan: exact lattice null +
            # Benjamini-Hochberg). They used to be dropped here, which forced downstream analyses to
            # threshold the RAW log-odds score -- an uncalibrated cutoff. Carry them so an operator hit
            # can be defined at a stated false-discovery rate instead.
            # `start`/`end`/`dyad` are in the CONCATENATED genome frame -- the frame the scan and
            # every gene coordinate here live in, so they stay as they are. But a multi-replicon
            # assembly is one string, and reporting a plasmid hit under the chromosome's accession
            # is a mislabelling: `replicon*` resolves the frame back to the record it belongs to.
            "operator": {"start": h.start, "end": h.end, "strand": h.strand, "dyad": center,
                         "replicon": _rep_id, "replicon_start": _rep_start,
                         "replicon_end": _rep_start + (h.end - h.start),
                         "score": (float(h.score) if getattr(h, "score", None) is not None else None),
                         "pvalue": (float(h.pvalue) if getattr(h, "pvalue", None) is not None else None),
                         "qvalue": (float(h.qvalue) if getattr(h, "qvalue", None) is not None else None),
                         "seq": getattr(h, "seq", "") or ""},
            "regulated_gene": (g.name or g.protein_id) if g else None,
            # `regulated_gene` is the human gene name and is NOT unique within a genome, so it
            # cannot be joined against an annotation on its own. The locus tag can, and the
            # replicon above says which record it is on.
            "regulated_gene_locus": (getattr(g, "locus_tag", "") or None) if g else None,
            # The whole OPERON the hit drives, not only its first gene. An operator upstream of a
            # polycistronic unit regulates every gene in it, so reporting one gene per operator
            # systematically under-called every multi-gene target (nikABCDE, znuABC, feoAB ...).
            "regulated_operon": _operon_gene_names(ctx, g),
            "gene_product": (getattr(g, "product", None) if g else None),
            "operator_to_gene_bp": dist,
            # `mc.mode` is a geometric tag ('promoter-core-overlap', 'elongated-spacer-overlap',
            # 'no-promoter-overlap'), computed the same way for every family and asserting no mechanism.
            "mode": (mc.mode if mc else None),
            "tss": tss,
            "op_to_tss_bp": op_to_tss,
            "op_start_to_tss": fp.get("op_start_to_tss"),
            "op_end_to_tss": fp.get("op_end_to_tss"),
            "elements": fp.get("elements"),          # -35/spacer/-10 boxes relative to TSS (promoter zoom)
            "occludes": (mc.occludes if mc else []),
            "overlaps": (mc.overlaps if mc else {}),
            "spacer_len": (mc.spacer_len if mc else None),
            "mode_confidence": (mc.confidence if mc else None),
            "flags": (mc.flags if mc else []),
        })
    return out


def confirm_regulated_genes(per_hit_records, ligand_gene_candidates, *, min_score=None,
                            qvalue_thresh: float = 0.05):
    """Annotate inducer-inferred genes with operator support.

    A lexical candidate with a promoter hit is ``operator_supported``. It is ``confirmed`` only when the
    hit also passes the explicitly stated q-value threshold; ordinary ranked candidates are not promoted
    to independent confirmation.
    """
    hit_genes = {}
    for r in per_hit_records:
        if min_score is not None and (r["operator"].get("score") or 0) < min_score:
            continue
        # match on ANY gene of the driven operon, not just its first gene: an operator upstream of
        # nikABCDE confirms nikB as much as nikA.
        for g in (r.get("regulated_operon") or []) or [r.get("regulated_gene")]:
            if g:
                hit_genes.setdefault(g, r)
    out = []
    for c in ligand_gene_candidates:
        rec = hit_genes.get(c.get("gene"))
        cc = dict(c)
        q = ((rec or {}).get("operator") or {}).get("qvalue")
        cc["operator_supported"] = rec is not None
        cc["confirmed"] = bool(rec is not None and q is not None and q == q and q <= qvalue_thresh)
        cc["support_label"] = (
            f"operator hit passes q<={qvalue_thresh:g}" if cc["confirmed"]
            else ("ranked operator-supported candidate" if rec else "no operator support")
        )
        if rec:
            cc["confirming_operator"] = rec["operator"]
            cc["predicted_mode"] = rec["mode"]
        out.append(cc)
    return out


# --------------------------------------------------------------------------- self-test
def _demo():
    from types import SimpleNamespace as NS
    prom = "TTGACA" + "GCTAGCATCGATCGAT" + "TATAAT"
    op = "ATATGAACAAATATTCATAT"
    head = "A" * 200
    genome = head + prom + op + ("C" * 40) + ("ACGT" * 500)
    gA_start = len(head + prom + op + ("C" * 40))
    genes = [NS(name="geneA", protein_id="A", product="zinc transporter", start=gA_start,
                end=gA_start + 300, strand="+"),
             NS(name="geneB", protein_id="B", product="hypothetical", start=4000, end=4300, strand="-")]
    # A REAL GenomeContext, and a two-replicon one: `per_hit_regulation` must resolve each hit
    # back to the record it sits on. A bare namespace here would leave that path untested, which is
    # how every coordinate came to be reported under the first record's accession.
    from predictor.signals.motif_rescan import GenomeContext
    ctx = GenomeContext("REPL_A.1", genome, genes,
                        contig_offsets={"REPL_A.1": 0, "REPL_B.1": len(genome) - 400})
    assert ctx.locate(10) == ("REPL_A.1", 10)
    assert ctx.locate(len(genome) - 100) == ("REPL_B.1", 300), "hit past the offset is on replicon B"
    op_start = len(head + prom)
    center = op_start + len(op) // 2
    hits = [NS(start=op_start, end=op_start + len(op), dyad_center=center, strand="+", score=10.0, seq=op)]
    recs = per_hit_regulation(ctx, hits)
    r = recs[0]
    print(f"per_hit_regulation: op@{r['operator']['dyad']} -> gene={r['regulated_gene']} "
          f"({r['operator_to_gene_bp']}bp) mode={r['mode']} tss={r['tss']} conf={r['mode_confidence']}")
    assert len(recs) == 1 and r["regulated_gene"] == "geneA", "operator should map to the downstream geneA"
    assert r["operator"]["replicon"] == "REPL_A.1" and         r["operator"]["replicon_start"] == r["operator"]["start"], "replicon A is at offset 0"
    conf = confirm_regulated_genes(recs, [{"gene": "geneA", "product": "zinc transporter"},
                                          {"gene": "geneC", "product": "other"}])
    assert conf[0]["operator_supported"] is True and conf[0]["confirmed"] is False
    assert conf[1]["operator_supported"] is False and conf[1]["confirmed"] is False
    print(f"confirm_regulated_genes: geneA supported={conf[0]['operator_supported']} "
          f"(mode {conf[0].get('predicted_mode')}) geneC={conf[1]['operator_supported']}")
    print("OK: WS4 per-hit mode/TSS + best-hit regulated-gene confirmation.")


if __name__ == "__main__":
    _demo()
