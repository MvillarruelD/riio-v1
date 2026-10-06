"""score_discovered.py -- score every discovered SL1344 candidate on the axes that have truth.

Truth sources, kept deliberately narrow:
  family        UniProtKB (release in raw/uniprot_sl1344_headers.txt), SL1344 taxon 216597: the
                same rule as E. coli: a studied family's signature Pfam
                (ecoli_kb/build_ecoli_kb.STUDY_FAMILY_PFAM), or membership of the paper's own family
                sequence set (build_ecoli_kb.paper_family_members -- this is what puts ModE in LysR).
                This also gives the family census (discovered vs total per family).
  class, ion    only the seven Osman et al. 2019 metal sensors.
  operator      only Osman sensors with an experimentally supported target: truth is the 400 bp
                upstream of each experimental target unit's first gene (prediction-independent).

Candidates are matched to UniProt by RefSeq accession, then by exact sequence.
Outputs: salmonella_family_census.tsv, discovered_scores.tsv, discovered_scores_summary.json
"""
from __future__ import annotations

import csv
import functools
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
from predictor.effector import inducer_vocab as V  # noqa: E402

sys.path.insert(0, str(HERE.parent / "ecoli_kb"))
from study_families import PFAM_TO_STUDY, STUDY_FAMILY_PFAM, paper_family_members  # noqa: E402

sys.path.insert(0, str(HERE))
from salmonella_benchmark import SENSOR_ACC, TRUTH_IONS, anchor_md5s, seq_md5  # noqa: E402

UNIPROT = HERE / "raw" / "uniprot_sl1344.tsv"
WINDOW = 400


def read_fasta(path: Path) -> str:
    return "".join(ln.strip() for ln in path.read_text().splitlines() if not ln.startswith(">"))


@functools.cache
def _paper_members() -> dict:
    return paper_family_members()


def study_family(u: dict) -> str:
    """Studied family of a UniProt row: signature Pfam, or the paper's own family sequence set."""
    fams = {PFAM_TO_STUDY[p] for p in u["Pfam"].split(";") if p in PFAM_TO_STUDY}
    return ";".join(sorted(fams | _paper_members().get(u["Entry"], set())))


def unit_windows(genes: list[dict], targets: list[dict]) -> dict[tuple[str, str], tuple[int, int]]:
    """(sensor, unit) -> the 400 bp promoter window of the unit's first gene, scan coordinates."""
    by_tag = {g["locus_tag"]: g for g in genes}
    units = defaultdict(list)
    for t in targets:
        if t["provenance"] == "salmonella_experimental" and t["locus_tag"] in by_tag:
            units[(t["sensor"], t["target_unit"])].append(by_tag[t["locus_tag"]])
    out = {}
    for key, gs in units.items():
        if gs[0]["strand"] == "+":
            first = min(gs, key=lambda g: g["s"])
            out[key] = (max(1, first["s"] - WINDOW), first["s"])
        else:
            first = max(gs, key=lambda g: g["e"])
            out[key] = (first["e"], first["e"] + WINDOW)
    return out


def main() -> int:
    cfg = load_canonical()
    run = cfg["runs"]["salmonella"]
    jobs = Path(run["root"]) / "jobs"
    with (BP.ANALYSIS / run["manifest"]).open(encoding="utf-8") as fh:
        manifest = list(csv.DictReader(fh))

    # --- UniProt family census -------------------------------------------------------------------
    with UNIPROT.open(encoding="utf-8") as fh:
        up = list(csv.DictReader(fh, delimiter="\t"))
    by_wp, by_seq = {}, {}
    for u in up:
        for wp in (x for x in u["RefSeq"].split(";") if x):
            by_wp[wp] = u
        by_seq[u["Sequence"]] = u
    census = [u for u in up if study_family(u)]
    with (HERE / "salmonella_family_census.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["uniprot", "locus", "gene", "protein", "pfam", "study_family", "refseq"])
        for u in census:
            w.writerow([u["Entry"], u["Gene Names (ordered locus)"], u["Gene Names (primary)"],
                        u["Protein names"], u["Pfam"], study_family(u), u["RefSeq"]])

    with (HERE / "osman2019_targets.tsv").open(encoding="utf-8") as fh:
        targets = [t for t in csv.DictReader(fh, delimiter="\t") if t["locus_tag"]]
    sensor_of = {acc: s for acc, s in SENSOR_ACC.items()}
    anchors = anchor_md5s()
    rng = np.random.default_rng(GS.SEED)

    rows, matched_uniprot, windows, gaps = [], set(), None, None
    for m in sorted(manifest, key=lambda r: r["run_name"]):
        acc = m["run_name"].rsplit("__", 1)[-1]
        bundle = jobs / m["run_name"]
        d = GS.load_dossier(bundle)
        if gaps is None:
            genes = GS.load_genes(bundle)
            gaps = GS.intergenic(genes)
            windows = unit_windows(genes, targets)
        seq = read_fasta(BP.WORKSPACE / m["fasta"])
        u = by_wp.get(acc) or by_seq.get(seq)
        how = "refseq" if acc in by_wp else ("sequence" if u else "none")
        if u is None and len(seq) >= 50:
            # BITACORA trims candidates to the aligned region, so an exact match can fail where the
            # trimmed protein is still a unique substring of one UniProt entry.
            hits = [x for x in up if seq in x["Sequence"]]
            if len(hits) == 1:
                u, how = hits[0], "contained sequence"
        if u:
            matched_uniprot.add(u["Entry"])
        fam_truth = study_family(u) if u else ""
        sensor = sensor_of.get(acc)
        ind = d.get("inducers") or {}
        top, cands = ind.get("top"), [c for c in (ind.get("candidates") or []) if c]
        pred_cls = V.inducer_class_of(top) or "undetermined"
        row = {"tf": sensor or m.get("gene") or acc, "refseq_protein": acc, "run_name": m["run_name"],
               "uniprot": u["Entry"] if u else "", "uniprot_match": how,
               "osman_sensor": "yes" if sensor else "no",
               "anchored": "yes" if anchors.get(seq_md5(seq)) else "no",
               "truth_family": fam_truth, "pred_family": d.get("family"),
               "family_ok": "" if not fam_truth else str(d.get("family") in fam_truth.split(";")),
               "truth_class": "metal" if sensor else "", "pred_class": pred_cls,
               "class_ok": str(pred_cls == "metal") if sensor else "",
               "truth_ions": ";".join(sorted(TRUTH_IONS[sensor])) if sensor else "",
               "pred_inducer": top, "pred_candidates": ";".join(cands), "agreement": ind.get("agreement"),
               "autoregulated_top5": str(GS.autoregulated(d, 5))}
        if sensor:
            if V.normalise(top) in TRUTH_IONS[sensor]:
                row["ion"] = "hit"
            elif TRUTH_IONS[sensor] & {V.normalise(c) for c in cands}:
                row["ion"] = "shortlist"
            else:
                row["ion"] = "class only" if top and "unresolved" in top else "wrong ion"
            spans = [w for (s, _), w in windows.items() if s == sensor]
            row["n_experimental_units"] = len(spans)
            if spans:
                row.update({f"op_{k}": v for k, v in GS.score_operator(bundle, d, spans, gaps, rng).items()})
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
                "wrong": [f"{r['tf']} (called {r['pred_family']}, Pfam {r['truth_family']})" if key == "family_ok"
                          else r["tf"] for r in scored if r[key] == "False"]}

    osman = [r for r in rows if r["osman_sensor"] == "yes"]
    with_op = [r for r in osman if "op_rec_1" in r]
    disc_up = matched_uniprot
    summary = {
        "run": run["tag"], "predictor_commit": cfg["provenance"]["predictor_commit"],
        "uniprot_release": next((ln.split(":", 1)[1].strip() for ln in
                                 (HERE / "raw" / "uniprot_sl1344_headers.txt").read_text().splitlines()
                                 if ln.lower().startswith("x-uniprot-release:")), "?"),
        "census_source": "UniProtKB SL1344 (taxon 216597) proteins carrying a studied family's signature Pfam",
        "uniprot_match": {k: sum(r["uniprot_match"] == k for r in rows)
                          for k in ("refseq", "sequence", "contained sequence", "none")},
        "discovery": {
            "n_candidates": len(rows),
            "in_family_total": len(census),
            "in_family_discovered": sum(u["Entry"] in disc_up for u in census),
            "by_family": {f: {"total": sum(study_family(u) == f for u in census),
                              "discovered": sum(study_family(u) == f and u["Entry"] in disc_up for u in census)}
                          for f in STUDY_FAMILY_PFAM},
            "osman_sensors_discovered": f"{len(osman)}/{len(SENSOR_ACC)}",
        },
        "family": frac(rows, "family_ok"),
        "class": frac(osman, "class_ok"),
        "ion": {k: sum(r.get("ion") == k for r in osman) for k in ("hit", "shortlist", "class only", "wrong ion")}
               | {"n": len(osman)},
        "operator": {"all": GS.operator_block([{k[3:]: v for k, v in r.items() if k.startswith("op_")} for r in with_op])
                     | {"tfs_recovered_top1": [r["tf"] for r in with_op if r["op_rec_1"]],
                        "tfs_recovered_top5": [r["tf"] for r in with_op if r["op_rec_5"]]},
                     "excluded": [r["tf"] for r in osman if "op_rec_1" not in r]},
        "truth_family_note": (f"UniProt holds {len(up)} SL1344 proteins against {len(genes)} annotated genes, so the "
                              "family census is a lower bound and candidates absent from UniProt carry no family "
                              "truth; class, ion and operator truth exist only for the Osman sensors"),
        "unmatched_candidates": [r["refseq_protein"] for r in rows if r["uniprot_match"] == "none"],
    }
    (HERE / "discovered_scores_summary.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")
    s = summary
    print(f"UniProt match {s['uniprot_match']}; in-family discovered {s['discovery']['in_family_discovered']}/"
          f"{s['discovery']['in_family_total']}; family {s['family']['correct']}/{s['family']['n']} wrong={s['family']['wrong']}")
    print(f"class {s['class']['correct']}/{s['class']['n']}; ion {s['ion']}")
    for k in GS.TOP_K:
        a = s["operator"]["all"][str(k)]
        print(f"  operator top-{k}: {a['observed']}/{s['operator']['all']['n']} (chance {a['expected_by_chance']}, P={a['p_vs_chance']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
