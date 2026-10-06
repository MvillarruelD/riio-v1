"""
figure_data.py -- pure extractors that turn the dossier dict into the minimal rows/arrays each figure or
report section needs. No I/O, no matplotlib here -- keeps `figures.py` presentation-only and `render_report`
orchestration-only, and makes the extraction unit-testable.

Every function tolerates missing/partial dossiers (Stage A may lack Stage-B fields) and never raises.
Run `python -m predictor.report.figure_data` for a self-test on a synthetic dossier.
"""
from __future__ import annotations


# operator-source -> (human category, one-line rationale) for the FINAL operator table (§5)
_SOURCE_RATIONALE = {
    "conservation_natural": ("natural",
        "A genome site from the PWM rescan, among the highest site scores. Ranking uses the site "
        "score alone; a site at the regulator's own promoter (proximal) is the autoregulatory candidate."),
    "conservation_logo": ("designed (consensus)",
        "The per-column most frequent base of the reported motif. A designed starting point that "
        "need not occur in the genome; its score is the motif's mean information content."),
}
# `structural` and `deeppbs_iterative` were listed here until 2026-09-02. Neither could ever be produced:
# both were Stage-B outputs of a `tfop finalize` whose engines had been deleted, so every report carried
# two rows that were permanently `pending_af3`. Removed with Stage B -- see report/operators.py.

_TSS_WINDOW_DEFAULT = 300            # per_hit_regulation uses +/-300 bp; operator centre ~= this when interior


def inducer_display(inducers: dict | None) -> str:
    """Presentation label for new and legacy dossiers; never use this value for classification."""
    ind = inducers or {}
    stored = ind.get("top_display")
    if isinstance(stored, str) and stored.strip():
        return stored
    from predictor.effector import inducer_vocab
    return inducer_vocab.display_label(ind.get("top"), ind.get("candidates") or ())


def promoter_occlusion(dossier: dict) -> dict:
    """The autoregulatory operator's promoter-occlusion record.

    The pre-2026-09-02 key ``mode_of_regulation`` is no longer read: its one-release compatibility
    window has closed, and every bundle of the canonical runs carries only this key."""
    return dossier.get("promoter_occlusion") or {}


def _mean_ic(logo: dict) -> float:
    """Mean bits/column of a stored logo dict ({pwm, per_col_ic?, consensus, total_ic, width})."""
    if not logo:
        return 0.0
    tic = logo.get("total_ic")
    w = logo.get("width") or len(logo.get("consensus") or "") or (
        len(logo["pwm"][0]) if logo.get("pwm") else 0)
    return (float(tic) / w) if (tic and w) else 0.0


# --------------------------------------------------------------------------- §2 operator track
def operator_track_rows(dossier: dict) -> list[dict]:
    """One row per genome rescan hit for the operator-position track: dyad, score, strand, seq,
    is_autoregulatory (the hit nearest the TF's own-promoter anchor), regulated_gene."""
    hits = (dossier.get("rescan") or {}).get("hits") or []
    neigh = dossier.get("neighborhood") or {}
    center = neigh.get("center")
    # map dyad -> regulated gene from the per-hit dataset (if present)
    gene_by_dyad = {}
    for r in (dossier.get("per_hit_regulation") or []):
        op = r.get("operator") or {}
        if op.get("dyad") is not None:
            gene_by_dyad[op["dyad"]] = r.get("regulated_gene")
    # the autoregulatory hit = the rescan hit CLOSEST to the autoregulatory anchor, within a generous
    # window (dyad framing vs the anchor varies by a few tens of bp; the recovery threshold is ~80).
    auto_dyad, best = None, None
    nearest_bp = ((dossier.get("autoregulatory_recovery") or dossier.get("cognate_recovery")
                   or {}).get("nearest_bp"))
    tol = max(150, (nearest_bp + 1) if nearest_bp is not None else 0)
    if center is not None:
        for h in hits:
            dy = h.get("dyad")
            if dy is None:
                continue
            dd = abs(dy - center)
            if dd <= tol and (best is None or dd < best):
                best, auto_dyad = dd, dy
    rows = []
    for h in hits:
        dy = h.get("dyad")
        rows.append({"dyad": dy, "score": h.get("score"), "strand": h.get("strand"),
                     "seq": h.get("seq"), "is_autoregulatory": (dy is not None and dy == auto_dyad),
                     "regulated_gene": gene_by_dyad.get(dy)})
    return rows


# --------------------------------------------------------------------------- §3 operator-vs-TSS
def tss_offsets(dossier: dict) -> list[dict]:
    """Signed operator -> TSS distance (bp) for every per-hit record with a TSS. Carries the operator CENTRE
    offset (`offset`, for the histogram) PLUS the full FOOTPRINT span (`span`=[start,end] edges relative to
    the TSS) and which promoter elements the operator `occludes` -- so the figure can draw the operator as a
    span, not a point. Falls back to (window - tss) for older records that lack the exact fields."""
    out = []
    for r in (dossier.get("per_hit_regulation") or []):
        tss = r.get("tss")
        if tss is None:
            continue
        off = r.get("op_to_tss_bp")
        if off is None:
            dy = (r.get("operator") or {}).get("dyad")
            off = (min(dy, _TSS_WINDOW_DEFAULT) - tss) if dy is not None else (_TSS_WINDOW_DEFAULT - tss)
        s, e = r.get("op_start_to_tss"), r.get("op_end_to_tss")
        span = [int(s), int(e)] if (s is not None and e is not None) else None
        out.append({"offset": int(off), "span": span, "occludes": r.get("occludes") or [],
                    "mode": r.get("mode"), "gene": r.get("regulated_gene"),
                    "score": (r.get("operator") or {}).get("score")})
    return out


def cds_offsets(dossier: dict) -> list[dict]:
    """Signed operator -> assigned-CDS-start distance (bp) for every per-hit record that has one.

    The companion to `tss_offsets`, and the one a null model can be built against: `operator_to_gene_bp`
    is measured to the 5' end of the assigned gene as read from the genome annotation, so unlike the TSS
    it does not depend on the operator's own position. `report.ligand_regulon.per_hit_regulation` already
    computes it (nearest gene 5' end downstream of the operator, within 400 bp, strand-aware); it was
    simply never surfaced. `analysis/fig_operator_cds.py` rebuilt the manuscript's F3 around this anchor
    for exactly that reason."""
    out = []
    for r in (dossier.get("per_hit_regulation") or []):
        off = r.get("operator_to_gene_bp")
        if off is None:
            continue
        out.append({"offset": int(off), "gene": r.get("regulated_gene"),
                    "gene_locus": r.get("regulated_gene_locus"),
                    "score": (r.get("operator") or {}).get("score")})
    return out


# --------------------------------------------------------------------------- §2 genome neighborhoods
def _zoom_from(occludes, fp, mode) -> dict:
    """Build a promoter-zoom dict (operator span + -35/spacer/-10 boxes relative to TSS) from a footprint
    (either promoter_occlusion.footprint or a per_hit_regulation record)."""
    fp = fp or {}
    if fp.get("op_start_to_tss") is None or not fp.get("elements"):
        return {}
    return {"op": [fp["op_start_to_tss"], fp["op_end_to_tss"]], "elements": fp["elements"],
            "occludes": occludes or [], "mode": mode, "spacer_len": fp.get("spacer_len")}


def neighborhood_panels(dossier: dict, *, max_targets: int = 3) -> list[dict]:
    """Gene-neighborhood panels for the top-ranking natural operators: the TF's OWN locus (from
    `neighborhood`) first, then the best target operons (from `regulon`, lowest q-value) that carry a
    `local_genes` context. Each panel = {title, subtitle, genes:[...], operator:[s,e], zoom:{...}} where
    `zoom` is that operator's promoter schematic (empty {} if the promoter wasn't resolved)."""
    panels = []
    neigh = dossier.get("neighborhood") or {}
    # per-hit record per regulated gene -> the promoter zoom for that operator
    by_gene = {}
    for r in (dossier.get("per_hit_regulation") or []):
        g = r.get("regulated_gene")
        if g and g not in by_gene:
            by_gene[g] = r
    if neigh.get("genes"):
        mo = promoter_occlusion(dossier)                  # the autoregulatory operator's promoter
        panels.append({"title": f"{dossier.get('tf_id','TF')} — autoregulatory locus",
                       "subtitle": "the TF's own promoter", "genes": neigh["genes"],
                       "operator": neigh.get("operator"),
                       "zoom": _zoom_from(mo.get("occludes"), mo.get("footprint"), mo.get("mode"))})
    ops = [o for o in (dossier.get("regulon") or [])
           if isinstance(o, dict) and o.get("first_gene") and o.get("local_genes")]
    ops.sort(key=lambda o: (o.get("qvalue") if o.get("qvalue") is not None else 1.0))
    own_names = {g.get("name") for g in (neigh.get("genes") or [{}])}
    for o in ops:
        if o["first_gene"] in own_names:                 # skip the autoregulatory operon (already shown)
            continue
        q = o.get("qvalue")
        r = by_gene.get(o["first_gene"]) or {}
        panels.append({"title": f"target: {o['first_gene']} operon",
                       "subtitle": (f"q={q:.1e}" if isinstance(q, (int, float)) else "") +
                                   (f" · {'+'.join(o['genes'])}" if o.get("genes") else ""),
                       "genes": o["local_genes"], "operator": o.get("site"),
                       "zoom": _zoom_from(r.get("occludes"),
                                          {k: r.get(k) for k in ("op_start_to_tss", "op_end_to_tss",
                                                                 "elements", "spacer_len")}, r.get("mode"))})
        if len(panels) >= max_targets + 1:
            break
    return panels


# --------------------------------------------------------------------------- promoter zoom (footprint)
def promoter_zoom(dossier: dict) -> dict:
    """The autoregulatory operator's footprint AS A PROMOTER ZOOM: the operator span + the -35/spacer/-10
    element boxes, all relative to the TSS (x=0), for a schematic that shows which promoter element the
    operator physically occludes. {} when the promoter/footprint wasn't resolved."""
    mo = promoter_occlusion(dossier)
    return _zoom_from(mo.get("occludes"), mo.get("footprint"), mo.get("mode"))


# --------------------------------------------------------------------------- selectivity
def selectivity(dossier: dict) -> dict:
    """Operator-search selectivity: genome-wide hits per Mb at the fixed p-value cut. A high density flags a
    LOW-information / degenerate operator (matches much of the genome) -- a less selective TF, hence a
    leakier switch -- vs a sharp, selective operator with few genome matches. Interpretation, not pass/fail."""
    rescan = dossier.get("rescan") or {}
    n = rescan.get("n_hits")
    glen = rescan.get("region_len")
    if not n or not glen:
        return {}
    per_mb = round(n / (glen / 1e6), 1)
    if per_mb >= 120:
        band, note = "low", "high hit density → low-information/degenerate operator (less selective TF)"
    elif per_mb >= 40:
        band, note = "moderate", "moderate hit density"
    else:
        band, note = "high", "few genome matches → a sharp, selective operator"
    return {"n_hits": n, "genome_bp": glen, "hits_per_mb": per_mb, "selectivity": band, "note": note}


# --------------------------------------------------------------------------- §3 inducer evidence
def inducer_evidence_rows(dossier: dict) -> dict:
    """The inducer-inference evidence: top pick, agreement, coordination gate, and one row per source
    (ssn_cluster / coordination / metalnet / ligify_db / ligify) with its ligand + evidence
    dict. `metalnet` always shows ligand `—`: it votes on whether a metal site exists, not on which
    metal, so it never names an inducer (its metal-type guess is advisory, in the evidence dict)."""
    ind = dossier.get("inducers") or {}
    rows = []
    for c in (ind.get("calls") or []):
        rows.append({"source": c.get("source"), "ligand": c.get("ligand"),
                     "confidence": c.get("confidence"), "role": c.get("role"),
                     "evidence": c.get("evidence") or {}})
    return {"top": inducer_display(ind), "raw_top": ind.get("top"),
            "agreement": ind.get("agreement"),
            "coordination_gate": ind.get("coordination_gate"), "rows": rows,
            "notes": ind.get("notes") or []}


# --------------------------------------------------------------------------- §5 final operators
def final_operator_rows(dossier: dict) -> list[dict]:
    """The final operator list with source category + rationale + status (the headline §5 table)."""
    ops = dossier.get("operators") or {}
    rows = []
    for c in (ops.get("candidates") or []):
        src = c.get("source")
        cat, why = _SOURCE_RATIONALE.get(src, (src or "?", ""))
        rows.append({"source": src, "category": cat, "rationale": why,
                     "sequence": c.get("sequence"), "score": c.get("score"),
                     "score_type": c.get("score_type"), "status": c.get("status"),
                     "provenance": c.get("provenance") or {}})
    return rows


# --------------------------------------------------------------------------- verdict banner
def _chip(label, status, detail):
    return {"label": label, "status": status, "detail": detail}


def verdict(dossier: dict) -> list[dict]:
    """Per-stage green/amber/red sanity chips -- the 'did the pipeline work?' banner. Conservative:
    unknown -> 'na' (grey), never crashes on a partial dossier."""
    chips = []
    # genome + locus
    chips.append(_chip("genome", "green" if dossier.get("genome_accession") else "red",
                       dossier.get("genome_accession") or "not resolved"))
    chips.append(_chip("TF locus", "green" if dossier.get("tf_locus") else "red",
                       (lambda t: f"{t[0]}-{t[1]} ({t[2]})" if t else "not found")(dossier.get("tf_locus"))))
    # autoregulation: is there a strong operator at the TF's OWN promoter? (a finding, not a pass/fail --
    # a weak/non-autoregulator like CueR legitimately has "no", with its operator at target genes instead).
    cr = dossier.get("autoregulatory_recovery") or dossier.get("cognate_recovery") or {}
    nb = cr.get("nearest_bp")
    tgt = (dossier.get("target_operator") or {})
    if cr.get("recovered"):
        st, txt = "green", f"autoregulatory operator {nb} bp from the self-promoter anchor"
    elif tgt.get("n_sites_added"):
        st, txt = "amber", (f"no operator at the self-promoter; {tgt['n_sites_added']} operator(s) found "
                            f"at target genes (weak/non-autoregulator)")
    elif nb is not None and nb <= 200:
        st, txt = "amber", f"operator {nb} bp from the self-promoter anchor (borderline)"
    else:
        st, txt = "red", "no operator near the self-promoter and none at target genes"
    chips.append(_chip("autoregulation", st, txt))
    # logo informativeness (genomic operator logo)
    mic = _mean_ic(dossier.get("genomic_logo") or {})
    st = "green" if mic >= 1.0 else ("amber" if mic >= 0.5 else "red")
    # Truncate, not round: 0.496 printed as "0.50" beside a below-0.5 status reads as a contradiction.
    chips.append(_chip("motif informative", st,
                       f"reported motif {int(mic * 100) / 100:.2f} bits/col (≥0.5 check, ≥1.0 good)"))
    # fold QC
    s = dossier.get("structure") or {}
    if not s.get("folded"):
        chips.append(_chip("apo fold QC", "na", "no fold this run (folding disabled)"))
    elif s.get("qc_pass") is True:
        chips.append(_chip("apo fold QC", "green",
                           f"pLDDT {round(s.get('plddt_mean'),1) if s.get('plddt_mean') else '?'}, "
                           f"AFDB RMSD {s.get('qc_rmsd')} A"))
    elif s.get("qc_pass") is False:
        chips.append(_chip("apo fold QC", "red", f"chain-A vs AFDB RMSD {s.get('qc_rmsd')} A (high)"))
    else:
        chips.append(_chip("apo fold QC", "amber",
                           f"folded; QC inconclusive (pLDDT {s.get('plddt_mean')})"))
    # inducer agreement
    ind = dossier.get("inducers") or {}
    ag = ind.get("agreement")
    top = inducer_display(ind)
    if ind.get("top") and ag is not None:
        st = "green" if ag >= 0.5 else ("amber" if ag >= 0.3 else "red")
        chips.append(_chip("inducer", st, f"{top} (agreement {ag:.2f})"))
    else:
        chips.append(_chip("inducer", "na", "not inferred"))
    # operators + regulon counts
    nhits = len((dossier.get("rescan") or {}).get("hits") or [])
    chips.append(_chip("operators found", "green" if nhits else "red", f"{nhits} putative site(s)"))
    nreg = len([o for o in (dossier.get("regulon") or []) if isinstance(o, dict) and o.get("first_gene")])
    reg_stats = dossier.get("regulon_stats") or {}
    n_found = reg_stats.get("n_operons_found")
    reg_detail = (f"{nreg} shown of {n_found} mapped operons (report cap)"
                  if reg_stats.get("truncated") and n_found is not None else f"{nreg} operon(s) shown")
    chips.append(_chip("regulon", "green" if nreg else "amber", reg_detail))
    # operator selectivity: genome-wide hits per Mb (degenerate/low-info operator vs a sharp, selective one)
    sel = selectivity(dossier)
    if sel:
        st = {"high": "green", "moderate": "amber", "low": "red"}.get(sel["selectivity"], "na")
        chips.append(_chip("operator selectivity", st,
                           f"{sel['hits_per_mb']} hits/Mb — {sel['selectivity']} selectivity"))
    # TSS coverage + sanity (most hits should have a TSS, offsets not absurd)
    offs = tss_offsets(dossier)
    nph = len(dossier.get("per_hit_regulation") or [])
    if nph and offs:
        sane = sum(1 for o in offs if abs(o["offset"]) <= 250) / len(offs)
        st = "green" if (len(offs) >= 0.6 * nph and sane >= 0.5) else "amber"
        chips.append(_chip("TSS positioning", st, f"{len(offs)}/{nph} hits with a TSS"))
    else:
        chips.append(_chip("TSS positioning", "na", "no TSS calls"))
    return chips


# --------------------------------------------------------------------------- self-test
def _demo() -> None:
    doss = {
        "genome_accession": "NC_007795.1", "tf_locus": [2209006, 2209324, "+"],
        "autoregulatory_recovery": {"recovered": True, "nearest_bp": 48},
        "genomic_logo": {"consensus": "ATATGAACAGATAATCATAT", "total_ic": 24.0, "pwm": [[0]*20]*4},
        "structure": {"folded": True, "qc_pass": True, "plddt_mean": 86.8, "qc_rmsd": 2.48},
        "inducers": {"top": "Zn2+", "agreement": 0.33, "coordination_gate": True,
                     "calls": [{"source": "ssn_cluster", "ligand": "Zn2+", "confidence": 1.0,
                                "role": "metal", "evidence": {"cluster": "ArsR_c2"}}]},
        "rescan": {"region_len": 2_800_000, "n_hits": 2,
                   "hits": [{"dyad": 2208995, "score": 17.8, "strand": "+", "seq": "ATAT"},
                            {"dyad": 2300000, "score": 9.0, "strand": "-", "seq": "GGGG"}]},
        "neighborhood": {"center": 2208931,
                         "genes": [{"name": "czrA", "start": 2208700, "end": 2209000, "strand": "+"},
                                   {"name": "czrB", "start": 2209100, "end": 2210000, "strand": "+"}],
                         "operator": [2208980, 2209010]},
        "regulon": [{"first_gene": "czrA", "genes": ["czrA", "czrB"], "site": [2208980, 2209010],
                     "qvalue": 2e-8, "local_genes": [{"name": "czrA", "start": 2208700, "end": 2209000,
                                                      "strand": "+"}]},
                    {"first_gene": "SA_RS09665", "genes": ["SA_RS09665"], "site": [1_000_000, 1_000_030],
                     "qvalue": 9e-7, "local_genes": [{"name": "SA_RS09665", "start": 999500, "end": 1000400,
                                                      "strand": "-"}]}],
        "per_hit_regulation": [{"operator": {"dyad": 2208995, "score": 17.8}, "regulated_gene": "czrB",
                                "mode": "promoter-core-overlap", "tss": 319, "op_to_tss_bp": -19,
                                "op_start_to_tss": -29, "op_end_to_tss": -9, "occludes": ["-10", "TSS"],
                                "operator_to_gene_bp": -64}],
        "promoter_occlusion": {"mode": "promoter-core-overlap", "occludes": ["-10", "TSS"],
                               "footprint": {"op_start_to_tss": -16, "op_end_to_tss": 4, "spacer_len": 17,
                                             "elements": {"-35": [-35, -29], "spacer": [-29, -12],
                                                          "-10": [-12, -6]}}},
        "operators": {"candidates": [
            {"source": "conservation_natural", "sequence": "ATATGAAAGGATTGTCATAT", "score": 17.8,
             "score_type": "rescan_site_score", "status": "ready", "provenance": {"in_tf_neighborhood": True}},
            {"source": "conservation_logo", "sequence": "ATATGAAAGGATTGTCATAT", "score": 1.6,
             "score_type": "mean_ic_bits", "status": "ready", "provenance": {}}]},
    }
    tr = operator_track_rows(doss)
    assert tr[0]["is_autoregulatory"] is True and tr[1]["is_autoregulatory"] is False, tr
    assert tr[0]["regulated_gene"] == "czrB"
    off = tss_offsets(doss)
    assert off and off[0]["offset"] == -19 and off[0]["span"] == [-29, -9], off
    assert "TSS" in off[0]["occludes"], off
    # the CDS anchor is carried alongside, and is operator-INDEPENDENT (see cds_offsets)
    cds = cds_offsets(doss)
    assert cds and cds[0]["offset"] == -64 and cds[0]["gene"] == "czrB", cds
    # neighborhood panels: TF own locus first (with a promoter zoom), then a target operon
    panels = neighborhood_panels(doss)
    assert panels[0]["genes"] and "autoregulatory" in panels[0]["title"], panels
    assert panels[0]["zoom"]["elements"]["-10"] == [-12, -6], panels[0]["zoom"]     # per-panel zoom attached
    assert any("target" in p["title"] for p in panels), panels
    # selectivity: 2 hits / 2.8 Mb -> < 40 -> high selectivity
    sel = selectivity(doss)
    assert sel["hits_per_mb"] < 40 and sel["selectivity"] == "high", sel
    # promoter zoom: operator + element boxes relative to TSS
    pz = promoter_zoom(doss)
    assert pz["elements"]["-10"] == [-12, -6] and pz["op"] == [-16, 4], pz
    ev = inducer_evidence_rows(doss)
    assert ev["top"] == "Zn2+" and ev["rows"][0]["evidence"]["cluster"] == "ArsR_c2"
    fr = final_operator_rows(doss)
    assert fr[0]["category"] == "natural" and fr[1]["category"] == "designed (consensus)", fr
    # every row of the final table must be producible: no permanently-pending Stage-B placeholders
    assert all(r.get("status") != "pending_af3" for r in fr), fr
    v = {c["label"]: c["status"] for c in verdict(doss)}
    assert v["genome"] == "green" and v["autoregulation"] == "green" and v["apo fold QC"] == "green", v
    assert v["motif informative"] == "green", v
    # partial dossier must not crash and should flag red/na
    v2 = {c["label"]: c["status"] for c in verdict({})}
    assert v2["genome"] == "red" and v2["apo fold QC"] == "na", v2
    print("OK: figure_data extractors + verdict (track, tss_offsets footprint, neighborhood_panels, "
          "selectivity, inducer, final, verdict) verified.")


if __name__ == "__main__":
    _demo()
