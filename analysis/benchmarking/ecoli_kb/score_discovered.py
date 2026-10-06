"""score_discovered.py -- score EVERY discovered E. coli TF on every axis the census has truth for.

Replaces the all-axes panel (see the 2026-09-24 entry in PANEL_CRITERIA.md's change log): each axis is
scored over all discovered TFs that carry truth for THAT axis, so no TF is dropped from the family or
class axes because it lacks, say, a confirmed operator. Denominators therefore differ by axis and
are always printed.

  discovery    census TFs in the twelve studied families, discovered vs total, per family
  family       discovered TFs with a studied-family signature Pfam (truth = UniProt Pfam)
  class        discovered TFs with a known effector (RegulonDB or cited curation)
  metal ion    discovered metal sensors
  operator     discovered TFs with >=1 precise operator; top-1/5/10/all vs chance (genome_scoring)

Inputs : ecoli_tf.tsv, ecoli_sites.tsv (build_ecoli_kb.py); the canonical E. coli run.
Outputs: discovered_scores.tsv, discovered_scores_summary.json
"""
from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from canonical_config import load_canonical  # noqa: E402

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import bench_paths as BP  # noqa: E402
import genome_scoring as GS  # noqa: E402

sys.path.insert(0, str(BP.REPO))
from predictor.annotate import family_kb  # noqa: E402
from predictor.effector import inducer_vocab as V  # noqa: E402

from study_families import STUDY_FAMILY_PFAM  # noqa: E402

FAMILIES = list(STUDY_FAMILY_PFAM)


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def truth_class(effector_class: str) -> str | None:
    if effector_class in ("", "unknown"):
        return None
    return "metal" if "metal" in effector_class.split(";") else effector_class.split(";")[0]


def ion_outcome(top: str | None, cands: list[str], truth_ions: set[str]) -> str:
    if V.normalise(top) in truth_ions:
        return "hit"
    if truth_ions & {V.normalise(c) for c in cands}:
        return "shortlist"
    if V.inducer_class_of(top) == "metal":
        return "class only" if top and "unresolved" in top else "wrong ion"
    return "miss"


def main() -> int:
    cfg = load_canonical()
    run = cfg["runs"]["ecoli"]
    jobs = Path(run["root"]) / "jobs"
    with (BP.ANALYSIS / run["manifest"]).open(encoding="utf-8") as fh:
        manifest = {r["locus_tag"]: r for r in csv.DictReader(fh)}
    census = read_tsv(HERE / "ecoli_tf.tsv")
    by_b = {t["b_number"]: t for t in census}
    precise, ht = defaultdict(list), defaultdict(list)
    for s in read_tsv(HERE / "ecoli_sites.tsv"):
        if s["left"]:
            span = (int(s["left"]), int(s["right"]))
            if s["precise_operator"] == "yes":
                precise[s["b_number"]].append(span)
            elif s["evidence_tier"] == "ht_binding":
                ht[s["b_number"]].append(span)

    rng = np.random.default_rng(GS.SEED)
    gaps = None
    rows = []
    for b, m in sorted(manifest.items(), key=lambda kv: kv[1]["run_name"]):
        bundle = jobs / m["run_name"]
        d = GS.load_dossier(bundle)
        if gaps is None:
            gaps = GS.intergenic(GS.load_genes(bundle))
        t = by_b.get(b, {})
        ind = d.get("inducers") or {}
        top, cands = ind.get("top"), [c for c in (ind.get("candidates") or []) if c]
        pred_cls = V.inducer_class_of(top) or "undetermined"
        t_cls = truth_class(t.get("effector_class", ""))
        t_ions = {x for x in t.get("metal_ions", "").split(";") if x}
        fam_truth = t.get("study_family", "")
        row = {
            "tf": t.get("tf") or m.get("gene") or b, "b_number": b, "run_name": m["run_name"],
            "uniprot": t.get("uniprot", ""),
            "anchored": "yes" if t.get("uniprot") and family_kb.clusters_anchored_by(t["uniprot"]) else "no",
            "truth_family": fam_truth, "pred_family": d.get("family"),
            "family_ok": "" if not fam_truth else str(d.get("family") in fam_truth.split(";")),
            "truth_class": t_cls or "", "pred_class": pred_cls,
            "class_ok": "" if not t_cls else str(pred_cls == t_cls),
            "truth_ions": ";".join(sorted(t_ions)), "pred_inducer": top, "pred_candidates": ";".join(cands),
            "agreement": ind.get("agreement"),
            "ion": ion_outcome(top, cands, t_ions) if t_cls == "metal" else "",
            "n_precise_sites": len(precise[b]), "n_ht_sites": len(ht[b]),
            "autoregulated_top5": str(GS.autoregulated(d, 5)),
        }
        if precise[b]:
            op = GS.score_operator(bundle, d, precise[b], gaps, rng)
            row.update({f"op_{k}": v for k, v in op.items()})
        if ht[b]:
            op_ht = GS.score_operator(bundle, d, ht[b], gaps, rng)
            row.update({f"ht_{k}": v for k, v in op_ht.items()})
        rows.append(row)

    fields = []
    for r in rows:
        fields += [k for k in r if k not in fields]
    with (HERE / "discovered_scores.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(rows)

    def frac(sub, key):
        scored = [r for r in sub if r[key] in ("True", "False")]
        return {"correct": sum(r[key] == "True" for r in scored), "n": len(scored),
                "wrong": [r["tf"] for r in scored if r[key] == "False"]}

    def op(sub, prefix):
        with_truth = [r for r in sub if f"{prefix}_rec_1" in r]
        return GS.operator_block([{k[len(prefix) + 1:]: v for k, v in r.items() if k.startswith(prefix + "_")}
                                  for r in with_truth]) | {
            "tfs_recovered_top1": [r["tf"] for r in with_truth if r[f"{prefix}_rec_1"]],
            "tfs_recovered_top5": [r["tf"] for r in with_truth if r[f"{prefix}_rec_5"]],
            "primary_hits_by_locality": {loc: sum(1 for r in with_truth if r[f"{prefix}_rec_1"]
                                                  and r[f"{prefix}_primary_locality"] == loc)
                                         for loc in ("proximal", "distal")}}

    in_fam = [t for t in census if t["in_study_families"] == "yes"]
    disc_b = set(manifest)
    confusion = defaultdict(lambda: defaultdict(int))
    for r in rows:
        if r["truth_class"]:
            confusion[r["truth_class"]][r["pred_class"]] += 1
    metal_rows = [r for r in rows if r["truth_class"] == "metal"]
    summary = {
        "run": run["tag"], "predictor_commit": cfg["provenance"]["predictor_commit"],
        "census_source": "RegulonDB 14.5 transcription factors + UniProt K-12 proteins carrying a "
                         "studied family's signature Pfam (build_ecoli_kb.py)",
        "discovery": {
            "n_candidates": len(rows),
            "in_family_total": len(in_fam),
            "in_family_discovered": sum(t["b_number"] in disc_b for t in in_fam),
            "by_family": {f: {"total": sum(t["study_family"] == f for t in in_fam),
                              "discovered": sum(t["study_family"] == f and t["b_number"] in disc_b for t in in_fam),
                              "metal_total": sum(t["study_family"] == f and t["is_metal_sensor"] == "yes" for t in in_fam),
                              "metal_discovered": sum(t["study_family"] == f and t["is_metal_sensor"] == "yes"
                                                      and t["b_number"] in disc_b for t in in_fam)}
                          for f in FAMILIES},
            "metal_sensors_in_family": [t["tf"] for t in in_fam if t["is_metal_sensor"] == "yes"],
            "metal_sensors_discovered": [t["tf"] for t in in_fam if t["is_metal_sensor"] == "yes" and t["b_number"] in disc_b],
            "discovered_outside_family": [r["tf"] for r in rows if not r["truth_family"]],
        },
        "family": frac(rows, "family_ok"),
        "class": frac(rows, "class_ok") | {"confusion": {k: dict(v) for k, v in confusion.items()},
                                           "false_metal": [r["tf"] for r in rows if r["truth_class"] and r["truth_class"] != "metal"
                                                           and r["pred_class"] == "metal"],
                                           "false_metal_agreement": {r["tf"]: r["agreement"] for r in rows if r["truth_class"]
                                                                     and r["truth_class"] != "metal" and r["pred_class"] == "metal"}},
        "ion": {k: sum(r["ion"] == k for r in metal_rows) for k in ("hit", "shortlist", "class only", "wrong ion", "miss")}
               | {"n": len(metal_rows)},
        "operator": {"all": op(rows, "op"), "metal_sensors": op(metal_rows, "op")},
        "operator_ht": {"all": op(rows, "ht")},
        "anchored_metal_sensors": sum(r["anchored"] == "yes" for r in metal_rows),
    }
    (HERE / "discovered_scores_summary.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")
    s = summary
    print(f"discovered {s['discovery']['in_family_discovered']}/{s['discovery']['in_family_total']} in-family; "
          f"family {s['family']['correct']}/{s['family']['n']}; class {s['class']['correct']}/{s['class']['n']}; "
          f"ion {s['ion']}")
    for k in GS.TOP_K:
        a = s["operator"]["all"][str(k)]
        print(f"  operator top-{k}: {a['observed']}/{s['operator']['all']['n']} "
              f"(chance {a['expected_by_chance']}, P={a['p_vs_chance']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
