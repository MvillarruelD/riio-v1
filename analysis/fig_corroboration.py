#!/usr/bin/env python
"""F2 -- metal predictions per family, and how many carry operator support.

Replaces the old `fig_central.py`, which read `regulon_validation.csv` -- a file produced by a
one-off 57-regulator study and by no pipeline stage, so on a fresh run it either failed or drew
last year's numbers. Everything here comes from THIS run's bundles.

    TFOP_RUN_TAG=<tag> python analysis/fig_corroboration.py

WHAT "CORROBORATED" MEANS, precisely, because the word invites over-reading: the regulator is called
metal-responsive, and at least one gene in its predicted regulon that handles that metal also carries
a predicted operator site. Both halves come from this pipeline, so this is INTERNAL consistency
between two of its axes -- not independent validation, and not experimental truth. A regulator with
no such gene in range is "no candidate gene": untestable by this criterion rather than refuted.

The old panel also had a "contradicted" class. It is not reported here because the bundles record
whether a candidate gene has operator support, not whether the operator evidence points at a
DIFFERENT metal. Inventing the category from data that cannot express it would be worse than
omitting it.

Design note -- why stacked bars rather than one pie per family. Eleven families cannot be compared
across eleven pies: angle is the weakest quantitative encoding, and the comparison the reader
actually needs ("which families corroborate, and how does that scale with n?") is a length
comparison along a shared axis. The layout also matches the survey-facing panels, which use
"filled = has the property" horizontal bars.
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

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_paths as RP  # noqa: E402
from fig_style import INK, NEUTRAL_L, NEUTRAL_M, save  # noqa: E402

sys.path.insert(0, str(RP.REPO))

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",     # editable text in SVG
    "pdf.fonttype": 42,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.8,
    "legend.frameon": False,
})

#: One colour per family, from the shared house palette so this panel sits beside the others.
FAMC = {
    "ArsR": "#0F4D92", "MerR": "#9A4D8E", "Fur": "#42949E", "CsoR": "#2E7D5B",
    "DtxR": "#7B52AB", "Rrf2": "#3775BA", "MarR": "#C9752B", "GntR": "#767676",
    "TetR": "#9E9E9E", "CopY": "#B0AEB4", "LysR": "#8A8F98",
}
CORROB, UNSUP, NOGENE, NONMETAL = "corroborated", "unsupported", "no candidate gene", "non-metal"


def classify() -> tuple[dict, collections.Counter]:
    from predictor.effector import inducer_vocab as vocab

    rows = list(csv.DictReader(RP.MANIFEST.open(encoding="utf-8")))
    by_fam: dict[str, collections.Counter] = collections.defaultdict(collections.Counter)
    total = collections.Counter()
    for r in rows:
        d = RP.JOBS / r["run_name"] / "dossier.json"
        if not d.is_file():
            continue
        j = json.loads(d.read_text(encoding="utf-8", errors="replace"))
        top = ((j.get("inducers") or {}).get("top") or "")
        crg = j.get("confirmed_regulated_genes") or []
        if not vocab.is_metal(top):
            k = NONMETAL
        elif any(g.get("operator_supported") for g in crg):
            k = CORROB
        elif crg:
            k = UNSUP
        else:
            k = NOGENE
        by_fam[r["family"]][k] += 1
        total[k] += 1
    return by_fam, total


def main() -> int:
    by_fam, total = classify()
    if not by_fam:
        print(f"no bundles under {RP.JOBS} for {RP.MANIFEST.name}")
        return 1

    n_all = sum(total.values())
    metal = total[CORROB] + total[UNSUP] + total[NOGENE]
    testable = total[CORROB] + total[UNSUP]

    def n_metal_of(c):
        return c[CORROB] + c[UNSUP] + c[NOGENE]

    # Order by METAL predictions, not by family size: the claim is about metal sensing, and sorting
    # by total puts the two largest non-metal families (GntR, TetR) at the top of a figure that is
    # not about them.
    fams = sorted(by_fam, key=lambda f: (-n_metal_of(by_fam[f]), -sum(by_fam[f].values())))
    # No ALL row: its bar would run to 52 while no family exceeds 16, so the aggregate dominates the
    # panel and compresses exactly the per-family comparison the figure is for. The totals are in the
    # title, where they cost no ink.
    rows = [(f, by_fam[f]) for f in fams]

    # The x axis carries METAL predictions only. Drawn against all 140 candidates the interesting
    # segments collapse into the left tenth of the panel -- GntR's 9 metal calls vanish beside its 23
    # non-metal ones -- so the denominator moves into the right-hand label where it can be read.
    xmax = max(n_metal_of(c) for _, c in rows)
    fig, ax = plt.subplots(figsize=(5.4, 0.235 * len(rows) + 1.05))

    for i, (fam, c) in enumerate(rows):
        y = len(rows) - 1 - i
        base = FAMC.get(fam, NEUTRAL_M) if fam != "ALL" else INK
        segs = [(c[CORROB], base, 1.00), (c[UNSUP], base, 0.32), (c[NOGENE], NEUTRAL_L, 1.00)]
        x = 0
        for val, colour, alpha in segs:
            if not val:
                continue
            ax.barh(y, val, left=x, height=0.60 if fam != "ALL" else 0.46,
                    color=colour, alpha=alpha, edgecolor="white", linewidth=0.6, zorder=3)
            if val >= 3:
                ax.text(x + val / 2, y, str(val), ha="center", va="center", fontsize=5.8,
                        color="white" if alpha == 1.0 and colour != NEUTRAL_L else INK, zorder=4)
            x += val
        nm, nt = n_metal_of(c), c[CORROB] + c[UNSUP]
        rate = f"{100 * c[CORROB] / nt:.0f}%" if nt else "–"
        ax.text(xmax * 1.03, y, f"{nm}/{sum(c.values())}", va="center", ha="left",
                fontsize=6.0, color=INK if fam == "ALL" else NEUTRAL_M)
        ax.text(xmax * 1.30, y, rate, va="center", ha="right",
                fontsize=6.0, color=INK if fam == "ALL" else NEUTRAL_M)

    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f for f, _ in rows][::-1], fontsize=6.4)
    for lbl, (fam, _) in zip(ax.get_yticklabels()[::-1], rows):
        lbl.set_color(FAMC.get(fam, NEUTRAL_M) if fam != "ALL" else INK)
        if fam == "ALL":
            lbl.set_fontweight("bold")
    ax.set_xlim(0, xmax * 1.32)
    ax.set_xticks(range(0, xmax + 1, 2))
    ax.set_xlabel("regulators called metal-responsive", fontsize=6.8, labelpad=2)
    ax.tick_params(axis="x", labelsize=6.0, pad=1.5)
    ax.tick_params(axis="y", length=0, pad=2)
    ax.set_ylim(-0.75, len(rows) - 0.25)
    ax.grid(axis="x", color="#F3F3F5", lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)

    ax.text(xmax * 1.03, len(rows) - 0.35, "metal/n", fontsize=5.7, color=NEUTRAL_M, ha="left")
    ax.text(xmax * 1.30, len(rows) - 0.35, "corrob.", fontsize=5.7, color=NEUTRAL_M, ha="right")

    ax.set_title(f"{total[CORROB]} of {metal} metal-responsive calls carry an operator at a gene for "
                 f"that metal  ·  {total[CORROB]}/{testable} of the testable ones",
                 loc="left", fontsize=7.4, color=INK, pad=7, fontweight="bold")

    handles = [plt.Rectangle((0, 0), 1, 1, color=INK, alpha=a) for a in (1.0, 0.32)]
    handles += [plt.Rectangle((0, 0), 1, 1, color=NEUTRAL_L)]
    ax.legend(handles, ["operator support at a gene for that metal", "no operator support",
                        "no candidate gene in range (untestable)"],
              fontsize=5.7, ncol=3, loc="upper left", bbox_to_anchor=(0, -0.13),
              handlelength=1.0, columnspacing=1.0, handletextpad=0.45)

    fig.text(0.012, -0.085, "Corroboration is internal consistency between this run's inducer and "
                            "operator axes, not experimental validation.",
             fontsize=5.5, color=NEUTRAL_M, style="italic")

    RP.FIGURES.mkdir(parents=True, exist_ok=True)
    save(fig, str(RP.FIGURES / "F2_metal_corroboration"))
    print(f"metal {metal}/{n_all} · corroborated {total[CORROB]} · testable {testable}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
