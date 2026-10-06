"""discovery_recall.py -- which E. coli TFs does genome-first discovery find, which does it miss, and why.

Reads the census (`ecoli_tf.tsv`, built by build_ecoli_kb.py) and the canonical E. coli run named in
`analysis/canonical_runs.toml`, plus that run's discovery search output (the BITACORA cache for the
genome), and answers three questions:

1. Recall over the census: of the TFs carrying a studied family's signature domain, which were
   discovered, and for each miss, WHICH step lost it, read from the search's own tables:
     profile_scope   no clade profile of the TF's OWN family reports it (hmmsearch reporting
                     cutoff E<=10); a weak hit from another family's profile does not count
     threshold       its own family's profile reports it, but only above the run's cut (1e-5)
     curation        it passed the E-value cut and BITACORA's downstream steps (blastp support /
                     coverage / trimming) dropped it
2. Enrichment: are metal sensors over-represented among the discovered TFs, against (a) all TFs in
   the studied families and (b) all TFs in the census? One-sided Fisher exact tests.
3. Precision: discovered proteins the census does not list as a TF.

Outputs: discovery_recall.tsv, discovery_recall_summary.json
"""
from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from canonical_config import load_canonical  # noqa: E402

from scipy.stats import fisher_exact

HERE = Path(__file__).resolve().parent
ANALYSIS = HERE.parents[1]
GENOME = "GCF_000005845.2"
EVALUE_CUT = 1e-5          # discovery/bitacora.py SEARCH_EVALUE

#: Which studied family each clade seed profile belongs to. A weak hit from ANOTHER family's
#: profile is background similarity between HTH domains, not a near-miss, so it counts as
#: profile scope; only a sub-threshold hit from the TF's own family's profile is a threshold miss.
PROFILE_FAMILY = {
    "BsCzrA": "ArsR/SmtB", "mtNmtR": "ArsR/SmtB", "ecZntR": "MerR", "PbrR": "MerR",
    "ecFur": "Fur", "ecZur": "Fur", "CopY_saMecI": "CopY", "CsoR_bsCsoR": "CsoR/FrmR",
    "DtxR_bsMntR": "DtxR/MntR", "GntR_LldR": "GntR", "LysR_ModE": "LysR-type (LTTR)",
    "MarR_AdcR": "MarR/SlyA", "NikR_NikR": "NikR", "Rrf2_ecIscR": "Rrf2", "TetR_SczA": "TetR/AcrR",
}


def bitacora_dir() -> Path:
    explicit = os.environ.get("TFOP_BITACORA_CACHE")
    if explicit:
        return Path(explicit) / GENOME
    sys.path.insert(0, str(ANALYSIS))
    from project_config import PREDICTOR_ROOT
    sys.path.insert(0, str(PREDICTOR_ROOT))
    from predictor import resources
    return Path(resources.cache_path("bitacora")) / GENOME


def best_hmm_hits(bdir: Path) -> dict[str, tuple[float, str]]:
    """protein_id -> (best full-sequence E-value, profile) across every clade profile's tblout."""
    best: dict[str, tuple[float, str]] = {}
    for tbl in bdir.glob("*/hmmer/*.tblout"):
        profile = tbl.parent.parent.name
        for ln in tbl.read_text(encoding="utf-8").splitlines():
            if ln.startswith("#") or not ln.strip():
                continue
            f = ln.split()
            pid, ev = f[0], float(f[4])
            if pid not in best or ev < best[pid][0]:
                best[pid] = (ev, profile)
    return best


def main() -> int:
    cfg = load_canonical()
    run = cfg["runs"]["ecoli"]
    with (ANALYSIS / run["manifest"]).open(encoding="utf-8") as fh:
        cands = {r["locus_tag"]: r for r in csv.DictReader(fh)}
    with (HERE / "ecoli_tf.tsv").open(encoding="utf-8") as fh:
        census = list(csv.DictReader(fh, delimiter="\t"))
    bdir = bitacora_dir()
    if not list(bdir.glob("*/hmmer/*.tblout")):
        raise SystemExit("Discovery diagnostics require BITACORA tblout cache; set TFOP_BITACORA_CACHE to its parent.")
    hits = best_hmm_hits(bdir)

    rows = []
    for t in census:
        b = t["b_number"]
        found = b in cands
        np_ids = [x for x in t["refseq_protein"].split(";") if x]
        hit = min((hits[i] for i in np_ids if i in hits), default=None)
        if found:
            reason = ""
        elif t["in_study_families"] != "yes":
            reason = "outside_studied_families"
        elif hit is None or PROFILE_FAMILY.get(hit[1]) not in t["study_family"].split(";"):
            reason = "profile_scope"
        elif hit[0] > EVALUE_CUT:
            reason = "threshold"
        else:
            reason = "curation"
        c = cands.get(b, {})
        rows.append({
            "tf": t["tf"], "b_number": b, "refseq_protein": t["refseq_protein"],
            "study_family": t["study_family"], "in_study_families": t["in_study_families"],
            "in_regulondb": t["in_regulondb"], "is_metal_sensor": t["is_metal_sensor"],
            "effector_class": t["effector_class"], "metal_ions": t["metal_ions"],
            "discovered": "yes" if found else "no",
            "discovered_family": c.get("family", ""), "discovery_profile": c.get("discovery_profile", ""),
            "run_name": c.get("run_name", ""),
            "best_profile_hit": hit[1] if hit else "", "best_profile_evalue": f"{hit[0]:.2g}" if hit else "",
            "miss_reason": reason,
        })
    census_b = {t["b_number"] for t in census}
    not_in_census = [{"b_number": b, "gene": c["gene"], "family": c["family"], "product": c["product"],
                      "run_name": c["run_name"]} for b, c in cands.items() if b not in census_b]

    fields = list(rows[0])
    with (HERE / "discovery_recall.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r["in_study_families"] != "yes", r["discovered"] != "yes",
                                                r["study_family"], r["tf"].lower())))

    def table(pop):
        """2x2 for Fisher: rows metal/non-metal, cols discovered/not, over population `pop`."""
        a = sum(r["is_metal_sensor"] == "yes" and r["discovered"] == "yes" for r in pop)
        b = sum(r["is_metal_sensor"] == "yes" and r["discovered"] != "yes" for r in pop)
        c = sum(r["is_metal_sensor"] != "yes" and r["discovered"] == "yes" for r in pop)
        d = sum(r["is_metal_sensor"] != "yes" and r["discovered"] != "yes" for r in pop)
        odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        return {"metal_found": a, "metal_missed": b, "nonmetal_found": c, "nonmetal_missed": d,
                "metal_fraction_among_found": round(a / (a + c), 3) if a + c else None,
                "metal_fraction_in_population": round((a + b) / len(pop), 3),
                "odds_ratio": None if odds == float("inf") else round(float(odds), 2),
                "p_one_sided": float(f"{p:.3g}")}

    fam = [r for r in rows if r["in_study_families"] == "yes"]
    by_family = {}
    for r in fam:
        f = by_family.setdefault(r["study_family"], {"n": 0, "found": 0, "metal": 0, "metal_found": 0})
        f["n"] += 1
        f["found"] += r["discovered"] == "yes"
        f["metal"] += r["is_metal_sensor"] == "yes"
        f["metal_found"] += r["is_metal_sensor"] == "yes" and r["discovered"] == "yes"
    summary = {
        "run": run["tag"], "genome": GENOME, "bitacora_cache": str(bdir),
        "n_candidates": len(cands), "n_census": len(rows), "n_in_study_families": len(fam),
        "recall_in_study_families": f"{sum(r['discovered'] == 'yes' for r in fam)}/{len(fam)}",
        "metal_sensors_in_study_families": [r["tf"] for r in fam if r["is_metal_sensor"] == "yes"],
        "metal_sensors_found": [r["tf"] for r in fam if r["is_metal_sensor"] == "yes" and r["discovered"] == "yes"],
        "metal_sensors_missed": {r["tf"]: r["miss_reason"] for r in fam
                                 if r["is_metal_sensor"] == "yes" and r["discovered"] != "yes"},
        "metal_sensors_outside_study_families": {r["tf"]: r["metal_ions"] for r in rows
                                                 if r["is_metal_sensor"] == "yes" and r["in_study_families"] != "yes"},
        "miss_reasons": {k: sum(r["miss_reason"] == k for r in fam) for k in ("profile_scope", "threshold", "curation")},
        "by_family": dict(sorted(by_family.items())),
        "enrichment_vs_study_families": table(fam),
        "enrichment_vs_all_census_tfs": table(rows),
        "discovered_not_in_census": not_in_census,
    }
    (HERE / "discovery_recall_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({k: v for k, v in summary.items() if k not in ("by_family",)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
