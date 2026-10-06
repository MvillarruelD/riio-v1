"""build_survey_figures.py -- Figure 3 (four-genome survey) and Figure S2 (own-promoter recovery).

The survey genomes have no operator truth set, so NOTHING here is an accuracy claim: these panels
describe what the pipeline returned. Accuracy lives only in the benchmark figure (Figure 2).
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from figure_style import (BENCH_GENOMES, FS_BODY, FS_SMALL, C, SHORT, SURVEY_ORDER, W2, apply_style, numbers,
                          panel_head, save)
from panels import class_legend, class_panel, family_heatmap, metal_key, metal_panels


def survey_rows(N: dict) -> list[tuple[str, dict]]:
    per = N["survey"]["per_genome"]
    return [(g, per[g]) for g in SURVEY_ORDER if g in per]


def fig_survey(N: dict) -> None:
    rows = survey_rows(N)
    short_rows = [(SHORT[g], c) for g, c in rows]
    fig = plt.figure(figsize=(W2, 4.7))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.15], hspace=0.66, wspace=0.45,
                          left=0.125, right=0.985, top=0.94, bottom=0.075)

    ax = fig.add_subplot(gs[0, :2])
    family_heatmap(ax, N["families"], rows)
    panel_head(ax, "a", "Candidates by family", pad_pt=12)

    ax = fig.add_subplot(gs[0, 2])
    class_panel(ax, short_rows)
    class_legend(ax, loc="upper left", bbox_to_anchor=(0.0, -0.30), ncol=2)
    panel_head(ax, "b", "Inducer class", pad_pt=12, letter_dx_pt=-28)

    axes, _ = metal_panels(fig, gs[1, :2], rows, wspace=0.42)
    metal_key(fig.add_subplot(gs[1, 2]), rows[0][1]["autoregulation"]["top_k"])
    panel_head(axes[0], "c", "Metal calls, the ions they name, and the operators that support them",
               pad_pt=14, letter_dx_pt=-28)

    save(fig, "fig3_survey")


def fig_own_promoter(N: dict) -> None:
    """Own promoter among the top-k operators, per genome: survey genomes, then benchmark genomes."""
    rows = [(SHORT[g], c) for g, c in survey_rows(N)]
    rows += [(SHORT[lab], N["benchmark_genome_composition"][k]) for k, lab in BENCH_GENOMES]
    fig, ax = plt.subplots(figsize=(W2 * 0.55, 2.3))
    x = np.arange(len(rows))
    vals = [c["autoregulation"]["n"] for _, c in rows]
    tots = [c["n"] for _, c in rows]
    pct = [100 * v / t if t else 0 for v, t in zip(vals, tots)]
    n_survey = len(rows) - len(BENCH_GENOMES)
    colors = [C["metal_mid"]] * n_survey + [C["neutral_mid"]] * len(BENCH_GENOMES)
    ax.bar(x, pct, width=0.62, color=colors, edgecolor="white")
    for xi, p, v, t in zip(x, pct, vals, tots):
        ax.text(xi, p + 2, f"{v}/{t}", ha="center", va="bottom", fontsize=FS_SMALL)
    ax.axvline(n_survey - 0.5, color=C["neutral_light"], lw=0.8, ls="--")
    ax.text(n_survey - 0.4, 97, "benchmark genomes", fontsize=FS_SMALL, color=C["neutral_dark"], va="top")
    ax.set_xticks(x)
    ax.set_xticklabels([lab for lab, _ in rows], fontsize=FS_BODY, rotation=20, ha="right")
    ax.set_ylim(0, 100)
    ax.set_ylabel("candidates (%)")
    k = rows[0][1]["autoregulation"]["top_k"]
    panel_head(ax, "", f"Own promoter among the top {k} operators", letter_dx_pt=-10)
    save(fig, "figS2_own_promoter")


def main() -> int:
    apply_style()
    N = numbers()
    fig_survey(N)
    fig_own_promoter(N)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
