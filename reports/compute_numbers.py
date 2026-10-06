"""compute_numbers.py -- every number the figures, reports and Methods quote, in one JSON.

Reads ONLY the canonical run named in `analysis/canonical_runs.toml` and the benchmark outputs
derived from it, and writes `numbers.json`. Figure and report code reads that file and never
recomputes or hardcodes a result, so re-pointing at a new run is: edit canonical_runs.toml, re-run
the benchmarks, re-run this, rebuild.

    python reports/compute_numbers.py
"""
from __future__ import annotations

import csv
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "analysis"))
from canonical_config import load_canonical  # noqa: E402

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
WORKSPACE = HERE.parent
ANALYSIS = WORKSPACE / "analysis"
ECOLI_KB = ANALYSIS / "benchmarking" / "ecoli_kb"
SALM = ANALYSIS / "benchmarking" / "salmonella"
sys.path.insert(0, str(ANALYSIS))
sys.path.insert(0, str(ANALYSIS / "benchmarking"))
sys.path.insert(0, str(ECOLI_KB))
from project_config import PREDICTOR_ROOT  # noqa: E402

sys.path.insert(0, str(PREDICTOR_ROOT))
import genome_scoring as GS  # noqa: E402
from study_families import STUDY_FAMILY_PFAM  # noqa: E402
from predictor.effector import inducer_vocab as V  # noqa: E402

#: the twelve studied families, in the fixed order every figure uses (zeros included)
FAMILIES = list(STUDY_FAMILY_PFAM)
#: display names, in the order figures use
GENOME_NAME = {"MtubH37Rv": "M. tuberculosis", "MaviHomi": "M. avium",
               "VchoRFB16": "V. cholerae", "VvulNBRC": "V. vulnificus"}
SURVEY_GENOMES = list(GENOME_NAME.values())
#: a metal call corroborated by >=2 independent evidence channels
CORROBORATED = 0.5
#: autoregulation = a top-k operator in the regulator's own flanking intergenic region
AUTOREG_K = 5


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


#: the inducer-evidence channels other than the SSN clade, in display order
OTHER_CHANNELS = ["coordination", "metalnet", "neighborhood", "ligify"]
#: half-width of the genomic window drawn around a regulator in the context figure
CONTEXT_FLANK_BP = 4000
#: the project's "operator on site" slop, reused to match a rescan hit to a regulon operon site
SITE_SLOP_BP = 25
#: Pfam accession -> the DNA-binding family it names, for TFs OUTSIDE the twelve studied families.
#: A label lookup (Pfam's own family names), not a result; the first accession in this order that a
#: protein carries names its family.
PFAM_FAMILY = {"PF00486": "OmpR/PhoB", "PF00158": "NtrC-type", "PF02954": "NtrC-type",
               "PF03551": "PadR", "PF03459": "ModE"}


def pfam_family(pfam: str) -> str:
    accs = {p for p in str(pfam).split(";") if p}
    return next((name for acc, name in PFAM_FAMILY.items() if acc in accs), "")


def short_locus(tag: str) -> str:
    """'F0316_RS10150' -> 'RS10150'; 'Rv1234' stays."""
    return str(tag).rsplit("_", 1)[-1]


def resolved_ion(top: str) -> str:
    """The named ion of a metal call, or '' when the call names none (see composition())."""
    top = str(top)
    if V.inducer_class_of(top) != "metal" or "unresolved" in top:
        return ""
    return V.normalise(top) or top


def own_region_hit(dossier: dict, k: int = AUTOREG_K) -> tuple[int, dict] | None:
    """(rank, hit) of the best-ranked top-k operator in the regulator's own flanking region."""
    for rank, h in enumerate(GS.ranked_hits(dossier)[:k], 1):
        if h.get("generator") == "rescan_tier1":
            return rank, h
    return None


def evidence_channels(dossier: dict) -> dict:
    """What each inducer-evidence channel said, reduced to what the SSN-clade figure needs."""
    calls = {c["source"]: c for c in (dossier.get("inducers") or {}).get("calls", [])}
    ssn = calls.get("ssn_cluster") or {}
    other = [s for s in OTHER_CHANNELS
             if (calls.get(s) or {}).get("role") == "metal"
             or (s == "ligify" and any(V.inducer_class_of(str((calls.get(x) or {}).get("ligand"))) == "metal"
                                       for x in ("ligify", "ligify_db")))]
    return {"ssn_role": ssn.get("role") or "", "ssn_ligand": ssn.get("ligand") or "",
            "metal_channels": other}


def load_predictions(root: Path) -> pd.DataFrame:
    df = pd.read_csv(root / "results" / "operator_predictions.csv")
    df["tag"] = df["run_name"].str.split("__").str[0]
    df["genome_label"] = df["tag"].map(GENOME_NAME).fillna(df["tag"])
    df["cls"] = df["inducer_top"].apply(lambda t: V.inducer_class_of(str(t)) or "undetermined")
    df["ion_named"] = df["inducer_top"].apply(resolved_ion)
    # The pipeline's own family call, top-k autoregulation, evidence channels and PWM, read from each
    # dossier and seed PWM.
    rows = []
    for name in df["run_name"]:
        bundle = root / "jobs" / name
        d = GS.load_dossier(bundle)
        own = own_region_hit(d)
        rows.append({"pred_family": d.get("family") or "",
                     "autoreg": own is not None,
                     "autoreg_rank": own[0] if own else None,
                     **evidence_channels(d),
                     **pwm_stats(bundle)})
    return pd.concat([df, pd.DataFrame(rows, index=df.index)], axis=1)


def pwm_stats(bundle: Path) -> dict:
    """Seed-PWM quality: width, mean information content, effective sequences, palindromicity.

    Palindromicity is the Pearson r between the PWM and its reverse complement (rows A,C,G,T), so 1
    is a perfect inverted repeat and 0 no dyad symmetry.
    """
    p = bundle / "motif" / "seed_pwm.json"
    if not p.exists():
        return {"pwm_width": None, "pwm_mean_ic": None, "pwm_n_eff": None, "pwm_palindromicity": None}
    s = json.loads(p.read_text(encoding="utf-8"))
    m = np.asarray(s["pwm"], dtype=float)
    width = m.shape[1]
    return {"pwm_width": int(width),
            "pwm_mean_ic": float(s["total_ic"]) / width,
            "pwm_n_eff": int(s["n_effective"]),
            "pwm_palindromicity": float(np.corrcoef(m.ravel(), m[::-1, ::-1].ravel())[0, 1])}


def ssn_contribution(metal: pd.DataFrame) -> dict:
    """How the SSN clade contributes to each metal call (the TF's final inducer call is metal).

    clade_and_other: the clade names a metal and >=1 other channel also calls metal;
    clade_only:      the clade names a metal and no other channel does;
    other_only:      the clade abstains (no clade, or a clade with no label);
    clade_disagrees: the clade names a non-metal role and other channels carry the call.
    """
    def status(r) -> str:
        if r.ssn_role == "metal":
            return "clade_and_other" if r.metal_channels else "clade_only"
        return "other_only" if not r.ssn_role else "clade_disagrees"

    st = metal.apply(status, axis=1) if len(metal) else pd.Series(dtype=str)
    ion = metal[metal.ion_named != ""]
    from_clade = ion.apply(lambda r: r.ssn_role == "metal" and bool(r.ssn_ligand)
                           and (V.normalise(str(r.ssn_ligand)) or r.ssn_ligand) == r.ion_named, axis=1) \
        if len(ion) else pd.Series(dtype=bool)
    combos = Counter(tuple((["SSN clade"] if r.ssn_role == "metal" else []) + list(r.metal_channels))
                     for r in metal.itertuples())
    return {"status": {k: int((st == k).sum()) for k in
                       ("clade_and_other", "clade_only", "other_only", "clade_disagrees")},
            "ion_named_by_clade": int(from_clade.sum()),
            "ion_named_otherwise": int(len(ion) - from_clade.sum()),
            "channel_combinations": [{"channels": list(k), "n": int(v)} for k, v in combos.most_common()],
            "channel_totals": {c: int(sum(c in r.metal_channels for r in metal.itertuples()))
                               for c in OTHER_CHANNELS} | {"SSN clade": int((metal.ssn_role == "metal").sum())}}


def composition(df: pd.DataFrame) -> dict:
    """The survey axes: family mix, inducer-class mix, the metal funnel, autoregulation."""
    metal = df[df.cls == "metal"]
    corr = metal[metal.inducer_agreement >= CORROBORATED]
    # A metal call is ION-RESOLVED when its top call names a species. The bundle's `ion` column is
    # the STRUCTURAL metal code (CU/FE/MN/ZN) and is empty for metalloid calls such as As(III), so
    # reading it would report those two ArsR_c3 calls as unresolved.
    ions = Counter(V.normalise(str(t)) or str(t) for t in metal.inducer_top
                   if V.inducer_class_of(str(t)) == "metal" and "unresolved" not in str(t))
    autoreg = df.autoreg.astype(bool)
    fam_counts = Counter(df.pred_family)
    return {
        "n": int(len(df)),
        "families": {f: int(fam_counts.get(f, 0)) for f in FAMILIES},
        "families_outside_study": {f: int(n) for f, n in fam_counts.items() if f and f not in FAMILIES},
        "class_counts": dict(Counter(df.cls)),
        "metal_funnel": {"metal_calls": int(len(metal)),
                         "corroborated": int(len(corr)),
                         "ion_resolved": int(sum(ions.values())),
                         "ions": dict(ions.most_common())},
        "autoregulation": {"n": int(autoreg.sum()), "pct": round(100 * autoreg.mean(), 1),
                           "top_k": AUTOREG_K},
        "operator_supports_metal_call": {
            "n": int((autoreg & (df.cls == "metal")).sum()),
            "of_metal_calls": int(len(metal)),
            "pct": round(100 * (autoreg & (df.cls == "metal")).sum() / len(metal), 1) if len(metal) else None},
        # D3: an operator can only correspond to an ion when the call names one, so operator support
        # is scored over the ion-resolved metal calls; stacked by ion in the survey figure.
        "operator_supports_ion_call": {
            "n": int((autoreg & (df.ion_named != "")).sum()),
            "of_ion_resolved": int((df.ion_named != "").sum()),
            "pct": round(100 * (autoreg & (df.ion_named != "")).sum() / (df.ion_named != "").sum(), 1)
            if (df.ion_named != "").any() else None,
            "ions": dict(Counter(df.ion_named[autoreg & (df.ion_named != "")]).most_common())},
        "ssn_contribution": ssn_contribution(metal),
        "pwm": {"median_mean_ic_bits": round(float(df.motif_mean_ic.median()), 3),
                "median_homologs": int(df.n_homologs.median()),
                "homolog_range": [int(df.n_homologs.min()), int(df.n_homologs.max())]},
        "pwm_quality": pwm_quality(df),
        "regulon_operons_capped_at": int(df.regulon_operons.max()),
    }


def pwm_quality(df: pd.DataFrame) -> dict:
    q = df.dropna(subset=["pwm_n_eff"])
    fam = q.groupby("pred_family")
    return {"n": int(len(q)),
            "n_effective_counts": {str(int(k)): int(v) for k, v in sorted(Counter(q.pwm_n_eff).items())},
            "median_palindromicity": round(float(q.pwm_palindromicity.median()), 3),
            "median_mean_ic_bits": round(float(q.pwm_mean_ic.median()), 3),
            "median_width": int(q.pwm_width.median()),
            "by_family": {f: {"n": int(len(g)),
                              "median_palindromicity": round(float(g.pwm_palindromicity.median()), 3),
                              "median_mean_ic_bits": round(float(g.pwm_mean_ic.median()), 3)}
                          for f, g in fam if f},
            "per_tf": [{"palindromicity": round(float(r.pwm_palindromicity), 3),
                        "mean_ic_bits": round(float(r.pwm_mean_ic), 3), "n_eff": int(r.pwm_n_eff),
                        "family": r.pred_family} for r in q.itertuples()]}


#: how many regulators the context figure shows
CONTEXT_N = 6


def context_selection(df: pd.DataFrame, root: Path) -> dict:
    """The regulators Figure 4 shows, and the genomic window drawn for each.

    Eligible: a corroborated metal call that names an ion AND has a top-k operator in the regulator's
    own flanking region -- the call and the operator point at the same locus. Ranked by channel
    agreement, then operator rank, then the number of metal channels; one regulator per SSN clade so
    the set is not six copies of one clade.
    """
    ok = df[(df.ion_named != "") & df.autoreg & (df.inducer_agreement >= CORROBORATED)].copy()
    ok["n_channels"] = ok.metal_channels.apply(len) + (ok.ssn_role == "metal")
    ok = ok.sort_values(["inducer_agreement", "autoreg_rank", "n_channels", "run_name"],
                        ascending=[False, True, False, True])
    picked, seen = [], set()
    for r in ok.itertuples():
        clade = r.ssn_cluster if isinstance(r.ssn_cluster, str) and r.ssn_cluster else r.pred_family
        if clade in seen:
            continue
        seen.add(clade)
        picked.append((r, clade))
        if len(picked) == CONTEXT_N:
            break
    out = []
    for r, clade in picked:
        bundle = root / "jobs" / r.run_name
        d = GS.load_dossier(bundle)
        rank, hit = own_region_hit(d)
        genes = pd.read_csv(bundle / "genome" / "genes.tsv", sep="\t", dtype={"protein_id": str})
        acc = r.run_name.split("__")[-1]
        tf = genes[genes.protein_id == acc]
        if tf.empty:
            continue
        tf = tf.iloc[0]
        # Rescan hits are in the concatenated SCAN frame (their `accession` names the first replicon
        # even on a second chromosome), so the window is built in scan coordinates and reported in
        # replicon coordinates through the TF's own offset.
        offset = int(tf.scan_start) - int(tf.replicon_start)
        lo = min(int(tf.scan_start), int(hit["start"])) - CONTEXT_FLANK_BP
        hi = max(int(tf.scan_end), int(hit["end"])) + CONTEXT_FLANK_BP
        win = genes[(genes.replicon == tf.replicon) & (genes.scan_end >= lo) & (genes.scan_start <= hi)]
        out.append({
            "run_name": r.run_name, "genome": r.genome_label, "family": r.pred_family, "clade": clade,
            "tf_gene": (r.gene if isinstance(r.gene, str) and r.gene else str(tf["name"])),
            "locus_tag": str(tf.locus_tag), "product": str(tf["product"]),
            "ion": r.ion_named, "agreement": float(r.inducer_agreement),
            "channels": (["SSN clade"] if r.ssn_role == "metal" else []) + list(r.metal_channels),
            "operator": {"rank": rank, "start": int(hit["start"]) - offset, "end": int(hit["end"]) - offset,
                         "strand": hit.get("strand"), "score": round(float(hit["score"]), 2),
                         "qvalue": float(hit["qvalue"]) if hit.get("qvalue") is not None else None,
                         "seq": hit.get("seq", "")},
            "replicon": str(tf.replicon), "window": [int(lo) - offset, int(hi) - offset],
            "tf_span": [int(tf.replicon_start), int(tf.replicon_end), str(tf.strand)],
            "genes": [{"name": str(g["name"]), "locus_tag": str(g.locus_tag), "product": str(g["product"]),
                       "start": int(g.replicon_start), "end": int(g.replicon_end), "strand": str(g.strand)}
                      for _, g in win.iterrows()],
        })
    return {"criteria": "corroborated metal call (inducer_agreement >= "
                        f"{CORROBORATED}) that names an ion, with a top-{AUTOREG_K} operator in the "
                        "regulator's own flanking intergenic region; ranked by agreement, operator rank, "
                        "then number of metal channels; one regulator per SSN clade",
            "n_eligible": int(len(ok)), "selected": out}


#: how the survey table names the gene an operator regulates
TARGET_RULE = (f"Regulated gene: the first gene of every predicted-regulon operon whose site overlaps the "
               f"primary operator within {SITE_SLOP_BP} bp ('self' = the regulator's own gene); "
               "otherwise the nearest protein-coding gene on the same replicon whose 5' end lies "
               "downstream of the operator (italic).")
CHANNEL_SHORT = {"SSN clade": "clade", "coordination": "coord.", "metalnet": "MetalNet2",
                 "neighborhood": "neighb.", "ligify": "Ligify"}


def _has_name(g) -> bool:
    name = str(g["name"])
    return bool(name) and name != "nan" and name != str(g.locus_tag)


def _gene_label(g) -> str:
    """A gene's name, or its short locus tag when it has none."""
    return str(g["name"]) if _has_name(g) else short_locus(g.locus_tag)


def regulated_gene(bundle: Path, hit: dict, tf_locus: str, genes: pd.DataFrame) -> dict:
    """The gene(s) the primary operator regulates (TARGET_RULE). Everything is in the SCAN frame."""
    hs, he = int(hit["start"]), int(hit["end"])
    reg = json.loads((bundle / "regulon_v2.json").read_text(encoding="utf-8"))
    by_name = {str(g["name"]): g for _, g in genes.iterrows()}
    names: list[str] = []
    for o in reg.get("operons") or []:
        s, e = o["site"]
        if s <= he + SITE_SLOP_BP and e >= hs - SITE_SLOP_BP and o["first_gene"] not in names:
            names.append(o["first_gene"])
    if names:
        out = []
        for n in names:
            g = by_name.get(n)
            self_ = g is not None and str(g.locus_tag) == tf_locus
            out.append({"label": "self" if self_ else (_gene_label(g) if g is not None else n),
                        "product": "" if g is None else str(g["product"]), "self": self_,
                        "named": g is not None and _has_name(g)})
        return {"rule": "operon", "genes": out}
    # the replicon holding the hit, from the concatenation offsets (the hit's accession is unreliable)
    offs = json.loads((bundle / "genome" / "contig_offsets.json").read_text(encoding="utf-8"))
    rep = max(((o, r) for r, o in offs.items() if o <= hs), default=(0, genes.replicon.iloc[0]))[1]
    mid = (hs + he) / 2
    pc = genes[(genes.replicon == rep) & genes.protein_id.fillna("").ne("")]
    plus = pc[(pc.strand == "+") & (pc.scan_start >= mid)]
    minus = pc[(pc.strand == "-") & (pc.scan_end <= mid)]
    cands = [(abs(g.scan_start - mid), g) for _, g in plus.iterrows()] + \
            [(abs(mid - g.scan_end), g) for _, g in minus.iterrows()]
    if not cands:
        return {"rule": "none", "genes": []}
    g = min(cands, key=lambda t: t[0])[1]
    return {"rule": "nearest", "genes": [{"label": "self" if str(g.locus_tag) == tf_locus else _gene_label(g),
                                          "product": str(g["product"]),
                                          "self": str(g.locus_tag) == tf_locus, "named": _has_name(g)}]}


def survey_table(df: pd.DataFrame, root: Path) -> dict:
    """One row per survey candidate: the call, the clade, the primary operator and what it regulates."""
    fam_rank = {f: i for i, f in enumerate(FAMILIES)}
    rows = []
    for r in df.itertuples():
        bundle = root / "jobs" / r.run_name
        d = GS.load_dossier(bundle)
        hits = GS.ranked_hits(d)
        genes = pd.read_csv(bundle / "genome" / "genes.tsv", sep="\t", dtype={"protein_id": str})
        acc = r.run_name.split("__")[-1]
        tf = genes[genes.protein_id == acc]
        tf_locus = str(tf.iloc[0].locus_tag) if not tf.empty else str(r.locus_tag)
        gene = r.gene if isinstance(r.gene, str) and r.gene else ""
        calls = {c["source"]: c for c in (d.get("inducers") or {}).get("calls", [])}
        # the clade's metal shortlist, only for an unresolved metal call (never inducers.candidates,
        # which mixes in Ligify's organic ligands)
        shortlist = []
        if r.cls == "metal" and not r.ion_named:
            shortlist = [V.normalise(str(c)) or str(c)
                         for c in (calls.get("ssn_cluster") or {}).get("candidates") or []
                         if V.inducer_class_of(str(c)) == "metal"]
        top = str(r.inducer_top)
        channels = (["SSN clade"] if r.ssn_role == "metal" else []) + list(r.metal_channels)
        op = hits[0] if hits else None
        rows.append({
            "genome": r.genome_label, "run_name": r.run_name,
            "tf": gene or short_locus(tf_locus), "locus_tag": tf_locus, "named": bool(gene),
            "family": r.pred_family,
            "clade": r.ssn_cluster if isinstance(r.ssn_cluster, str) else "",
            "cls": r.cls, "metal_call": r.cls == "metal", "ion": r.ion_named,
            "inducer": r.ion_named or top, "shortlist": shortlist,
            "metal_channels": [CHANNEL_SHORT.get(c, c) for c in channels],
            "agreement": float(r.inducer_agreement),
            "operator": {"seq": op.get("seq", ""), "own_region": op.get("generator") == "rescan_tier1",
                         "score": round(float(op["score"]), 1)} if op else None,
            "own_region_top_k": bool(r.autoreg),
            "own_region_rank": int(r.autoreg_rank) if r.autoreg else None,
            "target": regulated_gene(bundle, op, tf_locus, genes) if op else {"rule": "none", "genes": []},
            "_order": (SURVEY_GENOMES.index(r.genome_label) if r.genome_label in SURVEY_GENOMES else 9,
                       fam_rank.get(r.pred_family, 99), (gene or tf_locus).lower()),
        })
    rows.sort(key=lambda x: x.pop("_order"))
    return {"target_rule": TARGET_RULE, "top_k": AUTOREG_K,
            "n": len(rows), "n_metal_calls": sum(x["metal_call"] for x in rows), "rows": rows}


def main() -> int:
    cfg = load_canonical()
    runs = cfg["runs"]
    out: dict = {"provenance": cfg["provenance"] | {"runs": {k: {"tag": v["tag"], "n_candidates": v["n_candidates"]}
                                                             for k, v in runs.items()}},
                 "definitions": {
                     "primary_operator":
                         "The single operator the pipeline recommends for a TF: the highest-scoring "
                         "site found when that TF's operator PWM is rescanned across the intergenic "
                         "DNA of its own genome (p <= 1e-4), ranked by PWM match score. The remaining "
                         "ranked sites are the operator candidates.",
                     "precise_operator":
                         "Truth: a binding site with genomic coordinates whose strongest evidence is "
                         "classical binding (footprinting/EMSA-type) at RegulonDB confidence Strong or "
                         "Confirmed.",
                     "operator_on_site":
                         "A predicted operator whose span overlaps a truth site span within 25 bp.",
                     "proximal_vs_distal":
                         "Proximal = the hit lies in the regulator's own flanking intergenic region "
                         "(rescan locality tier 1), where the locality prior already ranks that window "
                         "first; distal = anywhere else, where the motif alone placed it.",
                     "corroborated_metal_call":
                         f"inducer_agreement >= {CORROBORATED}, i.e. >=2 independent evidence channels agree.",
                     "metal_call":
                         "A TF is a metal call when its FINAL inducer call (the fused consensus the "
                         "pipeline reports, inducer_top) is a metal ion, a metalloid oxyanion or an "
                         "unresolved metal. A TF for which only one evidence channel named a metal, but "
                         "whose final call is something else, is not a metal call.",
                     "ion_resolved":
                         "A metal call is ion-resolved when the final call names one species "
                         "(e.g. Zn2+, Cu+, As(III)); 'divalent metal (ion unresolved)' is a metal call "
                         "that is not ion-resolved.",
                     "operator_supports_call":
                         f"One of the regulator's top-{AUTOREG_K} operators lies in its own flanking "
                         "intergenic region (rescan locality tier 1), so the inducer call and the "
                         "operator point at the same locus. Scored over ion-resolved calls only: an "
                         "operator can correspond to an ion only when the call names one. The locality "
                         "prior favours this region, so this is support, not independent proof.",
                     "regulon_cap":
                         "Predicted regulons are truncated to 25 operons per TF for display; every "
                         "regulon count is a function of that cap.",
                 }}

    # ---------------- benchmarks ----------------
    eco_recall = json.loads((ECOLI_KB / "discovery_recall_summary.json").read_text(encoding="utf-8"))
    kb_sources = json.loads((ECOLI_KB / "KB_SOURCES.json").read_text(encoding="utf-8"))
    salm_kb = read_tsv(SALM / "salmonella_kb.tsv")
    out["families"] = FAMILIES
    # The operator-scoring definition is the scorer's own docstring, so the text a reader sees
    # cannot drift from the code that produced the numbers.
    out["definitions"]["operator_scoring"] = " ".join(GS.__doc__.split())
    out["benchmark"] = {
        "ecoli": {
            "sources": kb_sources["regulondb"] | {"uniprot": kb_sources["uniprot"]},
            "census": kb_sources["counts"],
            "census_source": "RegulonDB transcription factors plus UniProt K-12 proteins carrying a "
                             "studied family's signature Pfam (build_ecoli_kb.py)",
            "recall": eco_recall,
            "scores": json.loads((ECOLI_KB / "discovered_scores_summary.json").read_text(encoding="utf-8")),
            "per_tf": read_tsv(ECOLI_KB / "discovered_scores.tsv"),
            "criteria": "analysis/benchmarking/ecoli_kb/PANEL_CRITERIA.md",
        },
        "salmonella": {
            "source": salm_kb[0]["source"] if salm_kb else "",
            "scores": json.loads((SALM / "discovered_scores_summary.json").read_text(encoding="utf-8")),
            "per_tf": read_tsv(SALM / "discovered_scores.tsv"),
            "anchored_sensors": sum(k["anchored"] == "yes" for k in salm_kb),
        },
        "independence_caveat":
            "Every discovered E. coli metal sensor anchors a shipped SSN clade, and every Osman "
            "sensor is >=91% identical to one of those anchors, so neither genome benchmark tests the "
            "inducer call independently of the curated anchors. The independent measure is the "
            "knowledge-base leave-one-out (analysis/benchmarking/core/bench_loo.py).",
    }

    # ---------------- metal-sensor accounting ----------------
    # Every known metal sensor of both genomes, whether discovery reached it, and what was called;
    # plus every non-metal TF called metal. This is the "did we miss any" table.
    eco_tf = {r["b_number"]: r for r in eco_per} if (eco_per := out["benchmark"]["ecoli"]["per_tf"]) else {}
    census_tf = read_tsv(ECOLI_KB / "ecoli_tf.tsv")
    census_pfam = {t["b_number"]: t["pfam"] for t in census_tf}
    recall_rows = read_tsv(ECOLI_KB / "discovery_recall.tsv")
    acct = []
    for r in recall_rows:
        if r["is_metal_sensor"] != "yes":
            continue
        call = eco_tf.get(r["b_number"], {})
        acct.append({"genome": "E. coli", "tf": r["tf"], "family": r["study_family"] or "outside the 12",
                     "family_name": r["study_family"] or pfam_family(census_pfam.get(r["b_number"], ""))
                     or "unassigned",
                     "known_ions": r["metal_ions"], "discovered": r["discovered"],
                     "miss_reason": r["miss_reason"], "called": call.get("pred_inducer", ""),
                     "ion": call.get("ion", "")})
    for r in out["benchmark"]["salmonella"]["per_tf"]:
        if r.get("osman_sensor") == "yes":
            acct.append({"genome": "SL1344", "tf": r["tf"], "family": r["truth_family"],
                         "family_name": r["truth_family"], "known_ions": r["truth_ions"], "discovered": "yes", "miss_reason": "",
                         "called": r.get("pred_inducer", ""), "ion": r.get("ion", "")})
    false_metal = [{"genome": "E. coli", "tf": r["tf"], "known_class": r["truth_class"],
                    "called": r["pred_inducer"], "agreement": r["agreement"]}
                   for r in eco_per if r.get("truth_class") and r["truth_class"] != "metal"
                   and r.get("pred_class") == "metal"]
    eco_acct = [a for a in acct if a["genome"] == "E. coli"]
    in_fam = [a for a in eco_acct if a["family"] != "outside the 12"]
    out["benchmark"]["metal_sensor_accounting"] = {
        "sensors": acct, "false_metal": false_metal,
        # The "did we miss any" answer as counts, so no figure or text has to add up rows.
        "ecoli_recovery": {
            "in_study_families": len(in_fam),
            "found": sum(a["discovered"] == "yes" for a in in_fam),
            "missed": {a["tf"]: a["miss_reason"] for a in in_fam if a["discovered"] != "yes"},
            "outside_study_families": sorted(a["tf"] for a in eco_acct if a["family"] == "outside the 12"),
            "outside_found": sorted(a["tf"] for a in eco_acct
                                    if a["family"] == "outside the 12" and a["discovered"] == "yes")},
        "salmonella_recovery": {
            "sensors": sum(a["genome"] == "SL1344" for a in acct),
            "metal_class": sum(a["genome"] == "SL1344" and a["ion"] in ("hit", "shortlist", "class only")
                               for a in acct),
            "exact_ion": sum(a["genome"] == "SL1344" and a["ion"] == "hit" for a in acct)},
        "ecoli_census_totals": {"tfs": len(census_tf),
                                "in_study_families": sum(t["in_study_families"] == "yes" for t in census_tf),
                                "other_families": sum(t["in_study_families"] != "yes" for t in census_tf),
                                "metal_sensors": sum(t["is_metal_sensor"] == "yes" for t in census_tf)},
    }

    # ---------------- discovered TFs outside the census families ----------------
    # The census counts only the twelve studied families; Figure 2a shows every discovered TF, so the
    # ones the census cannot place are listed here under their real family where Pfam names one.
    up_sl = {u["Entry"]: u["Pfam"] for u in read_tsv(SALM / "raw" / "uniprot_sl1344.tsv")}
    outside: dict[str, dict[str, list[str]]] = {"ecoli": {}, "salmonella": {}}
    for r in recall_rows:
        if r["discovered"] == "yes" and r["in_study_families"] != "yes":
            fam = pfam_family(census_pfam.get(r["b_number"], "")) or "unassigned"
            outside["ecoli"].setdefault(fam, []).append(r["tf"])
    for t in eco_recall.get("discovered_not_in_census", []):
        outside["ecoli"].setdefault("not in the census", []).append(str(t))
    for r in out["benchmark"]["salmonella"]["per_tf"]:
        if not r["truth_family"]:
            fam = ((pfam_family(up_sl.get(r["uniprot"], "")) or "unassigned") if r["uniprot"]
                   else "no UniProt match")
            outside["salmonella"].setdefault(fam, []).append(r["tf"])
    for g, fams in outside.items():
        out["benchmark"][g]["discovered_outside_census"] = fams

    # ---------------- survey + per-genome ----------------
    survey = load_predictions(Path(runs["survey"]["root"]))
    out["survey"] = {"all": composition(survey),
                     "per_genome": {g: composition(survey[survey.genome_label == g])
                                    for g in survey.genome_label.unique()}}
    out["benchmark_genome_composition"] = {
        "ecoli": composition(load_predictions(Path(runs["ecoli"]["root"]))),
        "salmonella": composition(load_predictions(Path(runs["salmonella"]["root"]))),
    }
    out["context"] = context_selection(survey, Path(runs["survey"]["root"]))
    out["survey_table"] = survey_table(survey, Path(runs["survey"]["root"]))

    (HERE / "numbers.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    n_metal = out["survey"]["all"]["metal_funnel"]
    print(f"wrote numbers.json | survey n={out['survey']['all']['n']} "
          f"metal {n_metal['metal_calls']} -> corroborated {n_metal['corroborated']} -> ion {n_metal['ion_resolved']}")
    for g in ("ecoli", "salmonella"):
        s = out["benchmark"][g]["scores"]
        print(f"  {g}: in-family {s['discovery']['in_family_discovered']}/{s['discovery']['in_family_total']}, "
              f"family {s['family']['correct']}/{s['family']['n']}, class {s['class']['correct']}/{s['class']['n']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
