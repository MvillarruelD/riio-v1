#!/usr/bin/env python
"""F6 -- what this run can claim, per family, and on what evidence.

    TFOP_RUN_TAG=<tag> python analysis/fig_evidence_ledger.py

The per-regulator reports say what each call is. Nothing said how much of the set is CORROBORATED,
where the pipeline is blind, or which families it never had the evidence to assess -- a reader had to
open 140 bundles to find out. This is that summary in one figure.

Panel A: every candidate, by family, split by what its inducer call rests on.
Panel B: for each family, how often each evidence axis actually SPOKE.

The two panels answer each other, and that pairing is the point. Six families return zero corroborated
metal calls, and panel B separates the two reasons. For GntR, TetR, MarR, CopY and LysR the
coordination gate implements no chemistry and the SSN holds almost no clade, so the only axis with
anything to say is MetalNet: those zeros are the absence of evidence. Rrf2's zero is the opposite --
its gate speaks for all five and fires NEGATIVE, and all five come out redox, which is the right
answer for a family that senses [Fe-S] cluster status rather than free iron.

A reader given panel A alone cannot tell those two cases apart, and they are different claims. Making
them distinguishable is the whole reason the figure has a second panel.

TIERS. A metal call is counted CORROBORATED when at least two protein-level sources independently
support it -- a positive coordination gate, a MetalNet site, a metal-sensing SSN clade. One such
source is SINGLE-SOURCE; none is CONTEXT-ONLY, meaning the call rests on the regulon, the genomic
neighbourhood or the family prior.

Three sources, three different kinds of evidence, which is what makes agreement between them worth
something: the gate reads coordinating residues in the fold, MetalNet predicts a site from a
co-evolution MSA, the clade places the sequence among characterised sensors. Context sources are
deliberately excluded from the tier -- a regulator sitting beside a metal transporter may regulate it
without ever binding the metal, so that is evidence about the genome and not about the protein.

`inducer_class_of` decides metal/redox/organic. It is not re-derived here: a second copy of a
classification drifts from the first and the two then disagree by a count nobody can explain.
"""
from __future__ import annotations

import collections
import csv
import json
import sys
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_paths as RP  # noqa: E402
from fig_style import INK, NEUTRAL_L, NEUTRAL_M, save  # noqa: E402

#: outcome -> (colour, label). Blues carry metal and darken with evidence; grey is "not a metal call".
TIERS = [
    ("metal_corroborated", "#0F4D92", "metal, corroborated\n(2+ protein-level sources)"),
    ("metal_single", "#5B9BD5", "metal, single source"),
    ("metal_context", "#BDD7EE", "metal, context only"),
    ("redox", "#9A4D8E", "redox"),
    ("organic", "#CFCECE", "organic / non-metal"),
    ("undetermined", "#F3F3F5", "no call"),
]
#: evidence axis -> column header. Kept to one short word: a matrix column is about 0.4 in wide and
#: "coordination chemistry" set across it overprints its neighbour. The subtitle carries the names.
# AXES and classify() live in `inducer_evidence` so that anything spanning MORE THAN ONE RUN
# can import the same rule. This module cannot be imported without a run tag (it imports
# run_paths), which is right for a figure and wrong for a cross-run table -- so the shared
# definition moved out and is imported back here. There is still exactly one of it.
from inducer_evidence import AXES, classify  # noqa: E402,F401


def load():
    rows = list(csv.DictReader(RP.MANIFEST.open(encoding="utf-8")))
    per = collections.defaultdict(collections.Counter)
    spoke = collections.defaultdict(collections.Counter)
    for r in rows:
        f = RP.JOBS / r["run_name"] / "dossier.json"
        if not f.is_file():
            continue
        out, ax = classify(json.loads(f.read_text(encoding="utf-8")))
        per[r["family"]][out] += 1
        per[r["family"]]["n"] += 1
        for k, v in ax.items():
            spoke[r["family"]][k] += bool(v)
    return per, spoke, len(rows)


def main() -> int:
    per, spoke, n_tot = load()
    if not per:
        print("no bundles found -- nothing to summarise")
        return 1
    fams = sorted(per, key=lambda f: -per[f]["n"])
    tot = collections.Counter()
    for f in fams:
        tot.update(per[f])

    h = 0.26 * len(fams) + 1.55
    fig, (axA, axB) = plt.subplots(
        1, 2, figsize=(7.1, h), gridspec_kw={"width_ratios": [2.45, 1.05], "wspace": 0.05})

    # ---------------------------------------------------------------- A: outcome per family
    for i, f in enumerate(fams):
        y = len(fams) - 1 - i
        x = 0
        for key, colour, _lab in TIERS:
            v = per[f][key]
            if not v:
                continue
            axA.barh(y, v, left=x, height=0.62, color=colour, edgecolor="white", linewidth=0.6)
            if v >= 4:                       # only label a segment wide enough to hold the number
                axA.text(x + v / 2, y, str(v), ha="center", va="center", fontsize=5.4,
                         color="white" if key in ("metal_corroborated", "redox") else INK)
            x += v
        axA.text(-0.8, y, f, ha="right", va="center", fontsize=6.6, color=INK)
        axA.text(x + 0.9, y, str(per[f]["n"]), ha="left", va="center", fontsize=6.0,
                 color=NEUTRAL_M)
    axA.set_xlim(0, max(per[f]["n"] for f in fams) * 1.10)
    axA.set_ylim(-0.7, len(fams) - 0.3)
    axA.set_yticks([])
    axA.set_xlabel("candidate regulators", fontsize=6.4, color=NEUTRAL_M)
    axA.tick_params(axis="x", labelsize=6.0, colors=NEUTRAL_M, length=2)
    for sp in ("left", "right", "top"):
        axA.spines[sp].set_visible(False)
    axA.spines["bottom"].set_color(NEUTRAL_L)

    # ---------------------------------------------------------------- B: did the axis speak?
    for i, f in enumerate(fams):
        y = len(fams) - 1 - i
        n = per[f]["n"]
        for j, (key, _lab) in enumerate(AXES):
            frac = spoke[f][key] / n if n else 0.0
            # shade by fraction, not by a yes/no: "2 of 32" and "32 of 32" are different claims
            axB.add_patch(Rectangle((j - 0.42, y - 0.31), 0.84, 0.62, linewidth=0.6,
                                    edgecolor="white",
                                    facecolor=plt.cm.Blues(0.10 + 0.72 * frac)))
            if frac:
                axB.text(j, y, f"{frac:.0%}" if frac < 0.995 else "all", ha="center", va="center",
                         fontsize=5.2, color="white" if frac > 0.55 else INK)
            else:
                axB.text(j, y, "–", ha="center", va="center", fontsize=6.0, color=NEUTRAL_M)
    axB.set_xlim(-0.6, len(AXES) - 0.4)
    axB.set_ylim(-0.7, len(fams) - 0.3)
    axB.set_yticks([])
    axB.set_xticks(range(len(AXES)))
    axB.set_xticklabels([lab for _k, lab in AXES], fontsize=5.4, color=NEUTRAL_M)
    axB.tick_params(length=0)
    for sp in ("left", "right", "top", "bottom"):
        axB.spines[sp].set_visible(False)

    axA.annotate("A · what each family's calls rest on", xy=(0, 1), xycoords="axes fraction",
                 xytext=(-52, 20), textcoords="offset points", fontsize=7.4, fontweight="bold",
                 color=INK, ha="left", va="baseline", annotation_clip=False)
    axA.annotate(f"{n_tot} candidates from 4 genomes · {tot['metal_corroborated']} metal calls "
                 f"carry two independent protein-level sources",
                 xy=(0, 1), xycoords="axes fraction", xytext=(-52, 9),
                 textcoords="offset points", fontsize=5.6, color=NEUTRAL_M, ha="left",
                 va="baseline", annotation_clip=False)
    axB.annotate("B · how often each axis spoke", xy=(0, 1), xycoords="axes fraction",
                 xytext=(-14, 20), textcoords="offset points", fontsize=7.4, fontweight="bold",
                 color=INK, ha="left", va="baseline", annotation_clip=False)
    axB.annotate("share of the family assessed at all, by the coordination gate, MetalNet,\n"
                 "an SSN clade assignment and a folded structure",
                 xy=(0, 1), xycoords="axes fraction", xytext=(-14, 3), textcoords="offset points",
                 fontsize=5.4, color=NEUTRAL_M, ha="left", va="baseline", annotation_clip=False)

    handles = [plt.Line2D([], [], marker="s", ls="", markersize=5.2, color=c, label=lab)
               for k, c, lab in TIERS if tot[k]]
    axA.legend(handles=handles, fontsize=5.4, ncol=3, loc="upper left",
               bbox_to_anchor=(-0.20, -0.135 - 0.30 / len(fams)), handletextpad=0.4,
               columnspacing=1.1, labelspacing=0.8)

    save(fig, str(RP.FIGURES / "F6_evidence_ledger"))
    plt.close(fig)

    out = RP.RESULTS / "evidence_ledger.tsv"
    keys = [k for k, _c, _l in TIERS]
    with out.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["family", "n"] + keys + [f"spoke_{k}" for k, _ in AXES])
        for f in fams:
            w.writerow([f, per[f]["n"]] + [per[f][k] for k in keys]
                       + [spoke[f][k] for k, _ in AXES])
        w.writerow(["ALL", n_tot] + [tot[k] for k in keys]
                   + [sum(spoke[f][k] for f in fams) for k, _ in AXES])

    n_metal = tot["metal_corroborated"] + tot["metal_single"] + tot["metal_context"]
    print(f"{n_tot} candidates · {len(fams)} families")
    print(f"   metal {n_metal}  = {tot['metal_corroborated']} corroborated "
          f"+ {tot['metal_single']} single-source + {tot['metal_context']} context-only")
    print(f"   redox {tot['redox']} · organic {tot['organic']} · no call {tot['undetermined']}")
    zero = [f for f in fams if not per[f]["metal_corroborated"]]
    print(f"   families with no corroborated metal call: {', '.join(zero)}")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
