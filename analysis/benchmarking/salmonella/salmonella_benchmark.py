"""salmonella_benchmark.py -- the Osman et al. 2019 panel in the same shape as the E. coli benchmark.

Applies the Salmonella adaptation of PANEL_CRITERIA.md (../ecoli_kb/PANEL_CRITERIA.md) to the seven
SL1344 metal sensors and scores the canonical Salmonella run (analysis/canonical_runs.toml):

  R1-R4   family in scope, single protein, effector known, ion known -- from Osman 2019 Table 1
  R5'     >=1 target gene with `salmonella_experimental` provenance (operator truth = its promoter)

Operator recovery is read from `score_targets.py`'s per-target table for the same run
(`<run root>/results/salmonella_target_recovery.tsv`), which must be produced first:

    TFOP_RUN_TAG=salm_rerun20260917 TFOP_RUN_ROOT=<root> python score_targets.py
    python salmonella_benchmark.py

Anchor status is decided by SEQUENCE (MD5 of the sensor protein against every shipped anchor's
`seq_md5`), because the sensors are named by RefSeq WP_ accessions and the anchors by UniProt.

Outputs: salmonella_kb.tsv, salmonella_panel_scores.tsv, salmonella_panel_scores_summary.json
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from canonical_config import load_canonical  # noqa: E402

HERE = Path(__file__).resolve().parent
ANALYSIS = HERE.parents[1]
sys.path.insert(0, str(ANALYSIS))
from project_config import PREDICTOR_ROOT  # noqa: E402

sys.path.insert(0, str(PREDICTOR_ROOT))
from predictor.annotate import family_kb  # noqa: E402
from predictor.effector import inducer_vocab as V  # noqa: E402

# score_targets resolves run paths at import (run_paths needs a tag); point it at the canonical run.
_CFG = load_canonical()
os.environ.setdefault("TFOP_RUN_TAG", _CFG["runs"]["salmonella"]["tag"])
os.environ.setdefault("TFOP_RUN_ROOT", _CFG["runs"]["salmonella"]["root"])
sys.path.insert(0, str(HERE))
from score_targets import SENSOR_ACC  # noqa: E402

#: Osman 2019 Table 1 cognate metal, as canonical ions. RcnR senses Co(II) and Ni(II).
TRUTH_IONS = {"MntR": {"Mn2+"}, "Fur": {"Fe2+"}, "RcnR": {"Co2+", "Ni2+"}, "NikR": {"Ni2+"},
              "CueR": {"Cu+"}, "ZntR": {"Zn2+"}, "Zur": {"Zn2+"}}
TRUTH_FAMILY = {"MntR": "DtxR/MntR", "Fur": "Fur", "RcnR": "CsoR/FrmR", "NikR": "NikR",
                "CueR": "MerR", "ZntR": "MerR", "Zur": "Fur"}


def seq_md5(seq: str) -> str:
    """Same normalisation as analysis/family_kb/kb_io.seq_md5."""
    return hashlib.md5(re.sub(r"\s+", "", seq).upper().encode()).hexdigest()


def identity(a: str, b: str) -> float:
    """Identical positions over the shorter sequence, BLOSUM62 global alignment."""
    from Bio import Align
    from Bio.Align import substitution_matrices
    al = Align.PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score = -10, -0.5
    aln = al.align(a, b)[0]
    same = sum(a[x] == b[y] for i, j in zip(*aln.aligned) for x, y in zip(range(*i), range(*j)))
    return same / min(len(a), len(b))


def ecoli_sequences(cfg: dict) -> dict[str, str]:
    """gene -> protein sequence for the canonical E. coli run's candidates (all 9 E. coli metal
    sensors in it anchor a shipped clade -- see ecoli_kb/panel.tsv)."""
    run = cfg["runs"]["ecoli"]
    out = {}
    with (ANALYSIS / run["manifest"]).open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            fasta = ANALYSIS.parent / r["fasta"]
            out[r["gene"].lower()] = "".join(ln.strip() for ln in fasta.read_text().splitlines()
                                             if not ln.startswith(">"))
    return out


#: A sensor this close to an anchored E. coli ortholog inherits its clade label by near-lookup
#: even when its own sequence anchors nothing, so "non-anchored" is not "independent" here.
NEAR_IDENTICAL = 0.90


def anchor_md5s() -> dict[str, list[str]]:
    out = defaultdict(list)
    for cid, aset in family_kb._load().items():
        for a in aset.anchors:
            if a.seq_md5:
                out[a.seq_md5].append(cid)
    return out


def main() -> int:
    cfg = load_canonical()
    run = cfg["runs"]["salmonella"]
    root = Path(run["root"])
    with (ANALYSIS / run["manifest"]).open(encoding="utf-8") as fh:
        manifest = {r["run_name"].rsplit("__", 1)[-1]: r for r in csv.DictReader(fh)}
    rec_path = root / "results" / "salmonella_target_recovery.tsv"
    if not rec_path.is_file():
        raise SystemExit(f"run score_targets.py for {run['tag']} first: {rec_path} is missing")
    with rec_path.open(encoding="utf-8") as fh:
        recovery = list(csv.DictReader(fh, delimiter="\t"))
    anchors = anchor_md5s()
    eco = ecoli_sequences(cfg)

    kb, scores = [], []
    for acc, sensor in SENSOR_ACC.items():
        m = manifest.get(acc)
        seq = ""
        if m:
            fasta = ANALYSIS.parent / m["fasta"]
            seq = "".join(ln.strip() for ln in fasta.read_text().splitlines() if not ln.startswith(">"))
        anchored = anchors.get(seq_md5(seq), []) if seq else []
        eco_seq = eco.get(sensor.lower())
        ident = round(identity(seq, eco_seq), 3) if seq and eco_seq else None
        units = defaultdict(lambda: {"exp": False, "site": False, "reg": False})
        for t in recovery:
            if t["sensor"] != sensor:
                continue
            u = units[t["target_unit"]]
            u["exp"] |= t["provenance"] == "salmonella_experimental"
            u["site"] |= t["site_in_promoter"] == "True"
            u["reg"] |= t["in_regulon"] == "True"
        exp_units = [k for k, v in units.items() if v["exp"]]
        included = bool(m) and bool(exp_units)
        kb.append({"tf": sensor, "refseq_protein": acc, "locus_tag": (m or {}).get("locus_tag", ""),
                   "truth_family": TRUTH_FAMILY[sensor], "truth_class": "metal",
                   "truth_ions": ";".join(sorted(TRUTH_IONS[sensor])),
                   "n_target_units": len(units), "n_experimental_units": len(exp_units),
                   "included": "yes" if included else "no",
                   "reason": "all rules pass" if included else ("not discovered" if not m else "fails R5' (no experimental target)"),
                   "anchored": "yes" if anchored else "no", "anchor_clusters": ";".join(sorted(set(anchored))),
                   "identity_to_ecoli_ortholog": ident,
                   "near_identical_to_anchor": "yes" if (ident or 0) >= NEAR_IDENTICAL else "no",
                   "source": "Osman D et al. (2019) Nat Chem Biol 15:241-249, Table 1 + osman2019_targets.tsv"})
        if not m:
            scores.append({"tf": sensor, "discovered": "no", "in_panel": "no"})
            continue
        d = json.loads((root / "jobs" / m["run_name"] / "dossier.json").read_text(encoding="utf-8"))
        ind = d.get("inducers") or {}
        top, cands = ind.get("top"), [c for c in (ind.get("candidates") or []) if c]
        if V.normalise(top) in TRUTH_IONS[sensor]:
            ion = "hit"
        elif TRUTH_IONS[sensor] & {V.normalise(c) for c in cands}:
            ion = "shortlist"
        elif V.inducer_class_of(top) == "metal":
            ion = "class only" if top and "unresolved" in top else "wrong ion"
        else:
            ion = "miss"
        exp = {k: v for k, v in units.items() if v["exp"]}
        scores.append({
            "tf": sensor, "discovered": "yes", "in_panel": "yes" if included else "no",
            "anchored": "yes" if anchored else "no",
            "identity_to_ecoli_ortholog": ident,
            "independent": "no" if anchored or (ident or 0) >= NEAR_IDENTICAL else "yes",
            "pred_family": d.get("family"), "family_ok": str(d.get("family") == TRUTH_FAMILY[sensor]),
            "pred_inducer": top, "pred_candidates": ";".join(cands),
            "pred_class": V.inducer_class_of(top) or "undetermined",
            "class_ok": str(V.inducer_class_of(top) == "metal"), "ion": ion,
            "exp_units": len(exp), "exp_units_site_in_promoter": sum(v["site"] for v in exp.values()),
            "exp_units_in_regulon": sum(v["reg"] for v in exp.values()),
            "all_units": len(units), "all_units_site_in_promoter": sum(v["site"] for v in units.values()),
        })

    for name, rows in (("salmonella_kb.tsv", kb), ("salmonella_panel_scores.tsv", scores)):
        fields = list(max(rows, key=len))
        with (HERE / name).open("w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n", restval="")
            w.writeheader()
            w.writerows(rows)

    def summarise(rs):
        disc = [r for r in rs if r["discovered"] == "yes"]
        return {"n": len(rs), "discovered": f"{len(disc)}/{len(rs)}",
                "family": f"{sum(r['family_ok'] == 'True' for r in disc)}/{len(disc)}",
                "inducer_class": f"{sum(r['class_ok'] == 'True' for r in disc)}/{len(disc)}",
                "metal_ion": {k: sum(r["ion"] == k for r in disc) for k in ("hit", "shortlist", "class only", "wrong ion", "miss")},
                "experimental_units_site_in_promoter": f"{sum(r['exp_units_site_in_promoter'] for r in disc)}/{sum(r['exp_units'] for r in disc)}",
                "experimental_units_in_regulon": f"{sum(r['exp_units_in_regulon'] for r in disc)}/{sum(r['exp_units'] for r in disc)}",
                "all_units_site_in_promoter": f"{sum(r['all_units_site_in_promoter'] for r in disc)}/{sum(r['all_units'] for r in disc)}"}
    # The panel is what the rules admit. A sensor with no experimentally supported target in this
    # organism (NikR: nikA is absent from SL1344 and nixA is a motif transfer) has no operator truth,
    # so it is excluded by R5' and scored nowhere -- the same rule the E. coli panel applies.
    panel_scores = [r for r in scores if r.get("in_panel") == "yes"]
    summary = {"run": run["tag"], "predictor_commit": cfg["provenance"]["predictor_commit"],
               "promoter_window_bp": 400,
               "excluded_from_panel": {k["tf"]: k["reason"] for k in kb if k["included"] != "yes"},
               "all": summarise(panel_scores),
               "non_anchored": summarise([r for r in panel_scores if r.get("anchored") != "yes"]),
               "anchored": summarise([r for r in panel_scores if r.get("anchored") == "yes"]),
               "anchored_sensors": {k["tf"]: k["anchor_clusters"] for k in kb if k["anchored"] == "yes"},
               "identity_to_ecoli_ortholog": {k["tf"]: k["identity_to_ecoli_ortholog"] for k in kb},
               "independent_of_anchors": summarise([r for r in panel_scores if r.get("independent") == "yes"]),
               "note": (f"every sensor is >= {NEAR_IDENTICAL:.0%} identical to its E. coli ortholog, and all "
                        "nine E. coli metal sensors anchor a shipped clade, so no Salmonella sensor tests "
                        "the inducer call independently of the anchors")}
    (HERE / "salmonella_panel_scores_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
