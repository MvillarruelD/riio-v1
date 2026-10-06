#!/usr/bin/env python
"""Score the operator/regulon predictions against Osman et al. 2019's curated Salmonella targets.

    TFOP_RUN_TAG=salmonella20260831 python analysis/benchmarking/salmonella/score_targets.py

This is the tier-2 external check. It runs on whatever bundles exist, so it can be used mid-batch.

**Two measures, and only one of them is a real claim.**

*Regulon membership* asks whether the curated target gene appears anywhere in the predicted regulon.
It is weak: a regulon of ~60 genes assembled into operons will contain a gene for reasons that have
nothing to do with a site being right, and a hit anywhere inside a multi-gene operon counts the whole
operon. Reported because it is the number a reader will otherwise compute themselves.

*Site proximity* asks whether a predicted operator site lies in the window where the target's
promoter is -- upstream of the gene, or upstream of the FIRST gene of its operon, since that is where
a repressor binds. This is the claim. Recovering `iroD` because a site sits beside it, when the Fur
box for `iroBCDE` belongs upstream of `iroB`, is not a recovery.

**Provenance is carried through.** Six of the 22 curated targets are *E. coli* knowledge transferred
to Salmonella (four by homology, two by searching with an *E. coli* NikR box), so they are not
independent of the organism this benchmark exists to get away from. They are scored separately.
`nikA` does not exist in SL1344 and is excluded; `nixA` is a probable rather than confident locus
assignment and is flagged.
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent   # analysis/benchmarking/<suite>
ANALYSIS = HERE.parents[1]               # analysis/
sys.path.insert(0, str(ANALYSIS))
import run_paths as RP  # noqa: E402

TARGETS = HERE / "osman2019_targets.tsv"

#: Osman Table 1: the cognate metal of each sensor, for the tier-1 column.
COGNATE = {"MntR": "Mn", "Fur": "Fe", "RcnR": "Co", "NikR": "Ni",
           "CueR": "Cu", "ZntR": "Zn", "Zur": "Zn"}
#: run_name fragment -> sensor
SENSOR_ACC = {"WP_000131699.1": "Fur", "WP_000416272.1": "Zur", "WP_001026760.1": "CueR",
              "WP_000285586.1": "ZntR", "WP_000533542.1": "MntR", "WP_001190057.1": "NikR",
              "WP_000019953.1": "RcnR"}
PROMOTER_WINDOW = 400      # bp upstream of a gene start that counts as its promoter


def load_targets() -> list[dict]:
    rows = [r for r in csv.DictReader(TARGETS.open(encoding="utf-8"), delimiter="\t")]
    return [r for r in rows if r["locus_tag"]]        # nikA has none: absent from the organism


def bundle_for(sensor: str) -> Path | None:
    for acc, s in SENSOR_ACC.items():
        if s == sensor:
            for p in RP.JOBS.glob(f"*__{acc}"):
                if (p / "dossier.json").is_file():
                    return p
    return None


def score_sensor(sensor: str, targets: list[dict]) -> dict | None:
    b = bundle_for(sensor)
    if b is None:
        return None
    d = json.loads((b / "dossier.json").read_text(encoding="utf-8"))

    regulon_genes = set()
    operon_first = {}                       # gene -> first gene of its predicted operon
    for op in d.get("regulon") or []:
        gs = op.get("genes") or []
        for g in gs:
            regulon_genes.add(str(g))
            operon_first[str(g)] = str(gs[0]) if gs else None

    sites = []
    for h in d.get("per_hit_regulation") or []:
        o = h.get("operator") or {}
        if o.get("start") is not None:
            sites.append((int(o["start"]), int(o.get("end") or o["start"])))

    # The promoter of a gene inside an operon is upstream of the operon's FIRST gene, not of the
    # gene itself. Scored per gene, a site sitting beside iroD counts as recovering iroD -- but the
    # Fur box for iroBCDE is upstream of iroB, and a hit in the middle of the operon is not a
    # recovery of anything. `target_unit` carries the operon name from Supplementary Table 1, so the
    # window is computed once per unit, from whichever member starts it.
    by_unit = collections.defaultdict(list)
    for t in targets:
        by_unit[t["target_unit"]].append(t)
    unit_window = {}
    for unit, members in by_unit.items():
        strand = members[0]["strand"]
        if strand == "+":
            first = min(members, key=lambda m: int(m["start"]))
            lo, hi = int(first["start"]) - PROMOTER_WINDOW, int(first["start"])
        else:
            first = max(members, key=lambda m: int(m["end"]))
            lo, hi = int(first["end"]), int(first["end"]) + PROMOTER_WINDOW
        unit_window[unit] = (lo, hi, first["gene"])

    out = []
    for t in targets:
        lt, gene = t["locus_tag"], t["gene"]
        in_regulon = lt in regulon_genes or gene in regulon_genes
        lo, hi, first_gene = unit_window[t["target_unit"]]
        near = any(lo <= s <= hi or lo <= e <= hi for s, e in sites)
        out.append({**t, "in_regulon": in_regulon, "site_in_promoter": near,
                    "operon_first_gene": first_gene})

    ind = (d.get("inducers") or {}).get("top") or ""
    want = COGNATE.get(sensor, "")
    return {"sensor": sensor, "bundle": b.name, "inducer_predicted": ind,
            "inducer_expected": want,
            "inducer_ok": bool(want and want.lower() in ind.lower()),
            "n_regulon_operons": len(d.get("regulon") or []),
            "n_regulon_genes": len(regulon_genes), "n_sites": len(sites), "targets": out}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--window", type=int, default=PROMOTER_WINDOW)
    a = ap.parse_args(argv)
    globals()["PROMOTER_WINDOW"] = a.window

    tg = load_targets()
    by_sensor = collections.defaultdict(list)
    for t in tg:
        by_sensor[t["sensor"]].append(t)

    results, missing = [], []
    for sensor in COGNATE:
        r = score_sensor(sensor, by_sensor.get(sensor, []))
        (results if r else missing).append(r if r else sensor)

    print(f"scored {len(results)} of {len(COGNATE)} sensors"
          + (f"; no bundle yet for {', '.join(missing)}" if missing else ""))
    print(f"promoter window: {PROMOTER_WINDOW} bp upstream\n")

    print(f"{'sensor':<7}{'inducer predicted':<34}{'cognate':<9}{'ok':<4}"
          f"{'operons':>8}{'sites':>7}")
    for r in results:
        print(f"{r['sensor']:<7}{r['inducer_predicted'][:32]:<34}{r['inducer_expected']:<9}"
              f"{'yes' if r['inducer_ok'] else 'NO':<4}{r['n_regulon_operons']:>8}{r['n_sites']:>7}")

    print(f"\n{'sensor':<7}{'target':<9}{'locus tag':<17}{'provenance':<28}"
          f"{'in regulon':<12}site in promoter")
    rows = []
    for r in results:
        for t in r["targets"]:
            print(f"{r['sensor']:<7}{t['gene']:<9}{t['locus_tag']:<17}"
                  f"{t['provenance'][:26]:<28}{('yes' if t['in_regulon'] else '-'):<12}"
                  f"{'YES' if t['site_in_promoter'] else '-'}")
            rows.append({"sensor": r["sensor"], **{k: t[k] for k in
                        ("gene", "locus_tag", "target_unit", "operon_first_gene",
                         "provenance", "in_regulon", "site_in_promoter")}})

    if rows:
        exp = [x for x in rows if x["provenance"] == "salmonella_experimental"]
        eco = [x for x in rows if x["provenance"] != "salmonella_experimental"]
        print()
        for label, s in (("Salmonella-experimental", exp),
                         ("E. coli-transferred", eco), ("ALL", rows)):
            if not s:
                continue
            # per TRANSCRIPTIONAL UNIT: sitABCD is one Fur/MntR target, not four
            units = {}
            for x in s:
                u = units.setdefault((x["sensor"], x["target_unit"]), [False, False])
                u[0] |= bool(x["in_regulon"]); u[1] |= bool(x["site_in_promoter"])
            reg = sum(1 for v in units.values() if v[0])
            sit = sum(1 for v in units.values() if v[1])
            print(f"  {label:<26} genes {len(s):>3}   units {len(units):>2}   "
                  f"in regulon {reg:>2}/{len(units):<3}  site in promoter {sit:>2}/{len(units)}")
        print("\n  'site in promoter' is the claim; 'in regulon' counts a whole operon from one hit")

        out = RP.RESULTS / "salmonella_target_recovery.tsv"
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader(); w.writerows(rows)
        print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
