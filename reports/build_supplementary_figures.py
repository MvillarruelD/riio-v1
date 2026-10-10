"""build_supplementary_figures.py -- Figure S3 (operator / PWM quality) and Figure S4 (the SSN clade
and the metal call).

Both read the composition blocks of numbers.json for all six genomes: the four survey genomes and
the two benchmark genomes. Neither makes an accuracy claim; accuracy is Figure 2.
"""
from __future__ import annotations

from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

from figure_style import (BENCH_GENOMES, FS_BODY, FS_SMALL, C, SHORT, SURVEY_ORDER, W2, apply_style, is_species, numbers,
                          panel_head, save)

#: the families whose studied clades are metal-sensing by construction (study design, not a result)
METAL_FAMILIES = {"ArsR/SmtB", "MerR", "Fur", "DtxR/MntR", "CsoR/FrmR", "CopY", "NikR"}


def genome_rows(N: dict) -> list[tuple[str, dict]]:
    per = N["survey"]["per_genome"]
    rows = [(g, per[g]) for g in SURVEY_ORDER if g in per]
    return rows + [(lab, N["benchmark_genome_composition"][k]) for k, lab in BENCH_GENOMES]


def _xt(ax, labels: list[str], **kw) -> None:
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels([SHORT.get(lab, lab) for lab in labels], **kw)
    for t, lab in zip(ax.get_xticklabels(), labels):
        t.set_style("italic" if is_species(lab) else "normal")


# ---------------------------------------------------------------- Figure S3
def fig_pwm(N: dict) -> None:
    rows = genome_rows(N)
    per_tf = [t for _, c in rows for t in c["pwm_quality"]["per_tf"]]
    fig = plt.figure(figsize=(W2, 2.9))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.0, 1.0], wspace=0.55,
                          left=0.14, right=0.985, top=0.9, bottom=0.26)

    # a: palindromicity by family, pooled over the six genomes
    ax = fig.add_subplot(gs[0, 0])
    fams: dict[str, list[float]] = {}
    for t in per_tf:
        if t["family"]:
            fams.setdefault(t["family"], []).append(t["palindromicity"])
    order = sorted(fams, key=lambda f: np.median(fams[f]))
    med = [float(np.median(fams[f])) for f in order]
    ax.barh(range(len(order)), med, height=0.7,
            color=[C["metal"] if f in METAL_FAMILIES else C["organic"] for f in order], edgecolor="white")
    for i, f in enumerate(order):
        ax.text(med[i] + 0.01, i, f"n={len(fams[f])}", va="center", fontsize=FS_SMALL, color=C["neutral_dark"])
    overall = float(np.median([t["palindromicity"] for t in per_tf]))
    ax.axvline(overall, ls="--", lw=0.7, color=C["neutral_dark"])
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([f.replace("LysR-type (LTTR)", "LysR") for f in order], fontsize=FS_BODY)
    ax.set_xlim(0, 1)
    ax.set_xlabel("median palindromicity (r)")
    ax.legend(handles=[Patch(fc=C["metal"], label="metal-sensing family"), Patch(fc=C["organic"], label="other")],
              loc="upper left", bbox_to_anchor=(-0.05, -0.22), ncol=2, fontsize=FS_SMALL, handlelength=1.0)
    panel_head(ax, "a", f"Dyad symmetry (all: median r = {overall:.2f})", letter_dx_pt=-48)

    # b: effective sequences behind each PWM, per genome
    ax = fig.add_subplot(gs[0, 1])
    keys = sorted({int(k) for _, c in rows for k in c["pwm_quality"]["n_effective_counts"]})
    shades = [C["incorrect"], C["partial"], C["correct"], C["metal"]][:len(keys)] + \
        [C["neutral_dark"]] * max(0, len(keys) - 4)
    bottoms = np.zeros(len(rows))
    for k, col in zip(keys, shades):
        v = np.array([100 * c["pwm_quality"]["n_effective_counts"].get(str(k), 0) / c["pwm_quality"]["n"]
                      for _, c in rows])
        ax.bar(range(len(rows)), v, bottom=bottoms, color=col, edgecolor="white", width=0.66,
               label=f"{k}" + (" sequence" if k == 1 else " sequences"))
        bottoms += v
    _xt(ax, [g for g, _ in rows], rotation=35, ha="right", fontsize=FS_BODY)
    ax.set_ylim(0, 100)
    ax.set_ylabel("PWMs (%)")
    ax.legend(loc="upper left", bbox_to_anchor=(-0.05, -0.27), ncol=2, fontsize=FS_SMALL, handlelength=0.9,
              columnspacing=0.8)
    panel_head(ax, "b", "Effective sequences per PWM", letter_dx_pt=-26)

    # c: does a PWM built from more sequences carry more dyad symmetry?
    ax = fig.add_subplot(gs[0, 2])
    rng = np.random.default_rng(0)  # jitter only; fixed so a rebuild is byte-identical
    for i, (k, col) in enumerate(zip(keys, shades)):
        v = [t["palindromicity"] for t in per_tf if t["n_eff"] == k]
        ax.scatter(i + rng.uniform(-0.18, 0.18, len(v)), v, s=4, color=col, alpha=0.7, lw=0)
        if v:
            ax.plot([i - 0.28, i + 0.28], [np.median(v)] * 2, color=C["ink"], lw=1.0)
            ax.text(i, 1.02, f"n={len(v)}", ha="center", va="bottom", fontsize=FS_SMALL)
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels([str(k) for k in keys])
    ax.set_xlabel("effective sequences in the PWM")
    ax.set_ylabel("palindromicity (r); bar = median")
    ax.set_ylim(min(-0.2, min(t["palindromicity"] for t in per_tf) - 0.05), 1.12)
    panel_head(ax, "c", "Symmetry vs support", letter_dx_pt=-26)

    save(fig, "figS3_pwm_quality")


# ---------------------------------------------------------------- Figure S4
STATUS = [("clade_and_other", "clade + other channels", C["metal"]),
          ("clade_only", "clade only", C["metal_soft"]),
          ("other_only", "clade silent; other channels", C["organic"]),
          ("clade_disagrees", "clade says non-metal", C["incorrect"])]
CHANNELS = ["SSN clade", "coordination", "metalnet", "neighborhood", "ligify"]
CHANNEL_LABEL = {"metalnet": "MetalNet2", "neighborhood": "neighbourhood", "ligify": "Ligify"}


def fig_ssn(N: dict) -> None:
    rows = genome_rows(N)
    fig = plt.figure(figsize=(W2, 4.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.1], width_ratios=[1.5, 1.0], hspace=0.78,
                          wspace=0.35, left=0.10, right=0.985, top=0.93, bottom=0.04)

    # a: how the clade contributes to each genome's metal calls
    ax = fig.add_subplot(gs[0, 0])
    bottoms = np.zeros(len(rows))
    for key, lab, col in STATUS:
        v = np.array([c["ssn_contribution"]["status"][key] for _, c in rows])
        ax.bar(range(len(rows)), v, bottom=bottoms, color=col, edgecolor="white", width=0.66, label=lab)
        bottoms += v
    for i, (_, c) in enumerate(rows):
        ax.text(i, bottoms[i] + 0.4, f"{int(bottoms[i])}", ha="center", va="bottom", fontsize=FS_SMALL)
    _xt(ax, [g for g, _ in rows], fontsize=FS_BODY)
    ax.set_ylabel("metal calls")
    ax.set_ylim(0, bottoms.max() * 1.18)
    ax.tick_params(axis="x", pad=1.5)
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, -0.13), ncol=4, fontsize=FS_SMALL, handlelength=1.0,
              columnspacing=0.9)
    panel_head(ax, "a", "What the SSN clade says for each metal call", letter_dx_pt=-26)

    # b: which channels called metal, pooled over the six genomes (an UpSet-style matrix)
    combos: Counter = Counter()
    for _, c in rows:
        for e in c["ssn_contribution"]["channel_combinations"]:
            combos[tuple(e["channels"])] += e["n"]
    items = combos.most_common()
    sub = gs[1, :].subgridspec(2, 1, height_ratios=[1.0, 0.9], hspace=0.05)
    axb = fig.add_subplot(sub[0, 0])
    axm = fig.add_subplot(sub[1, 0], sharex=axb)
    x = np.arange(len(items))
    axb.bar(x, [n for _, n in items], width=0.6,
            color=[C["metal"] if "SSN clade" in k else C["organic"] for k, _ in items])
    for xi, (_, n) in zip(x, items):
        axb.text(xi, n + 0.3, str(n), ha="center", va="bottom", fontsize=FS_SMALL)
    axb.set_ylabel("metal calls")
    axb.tick_params(axis="x", bottom=False, labelbottom=False)
    axb.set_ylim(0, max(n for _, n in items) * 1.25)
    for j, ch in enumerate(CHANNELS):
        on = [ch in k for k, _ in items]
        axm.scatter(x, [j] * len(x), s=14, c=[(C["metal"] if ch == "SSN clade" else C["ink"]) if o
                                                else C["neutral_light"] for o in on], zorder=3)
        for xi, (k, _) in zip(x, items):
            idx = [CHANNELS.index(c) for c in k if c in CHANNELS]
            if len(idx) > 1:
                axm.plot([xi, xi], [min(idx), max(idx)], color=C["ink"], lw=0.8, zorder=2)
    axm.set_yticks(range(len(CHANNELS)))
    axm.set_yticklabels([CHANNEL_LABEL.get(c, c) for c in CHANNELS], fontsize=FS_BODY)
    axm.set_ylim(len(CHANNELS) - 0.5, -0.5)
    axm.tick_params(axis="x", bottom=False, labelbottom=False)
    for sp in ("top", "right", "bottom"):
        axm.spines[sp].set_visible(False)
    panel_head(axb, "c", f"Which channels called metal, pooled over the {len(rows)} genomes (blue bar: the "
               "SSN clade is one of them)", letter_dx_pt=-26)

    # c: where the named ion comes from
    ax = fig.add_subplot(gs[0, 1])
    by = np.array([c["ssn_contribution"]["ion_named_by_clade"] for _, c in rows])
    other = np.array([c["ssn_contribution"]["ion_named_otherwise"] for _, c in rows])
    ax.bar(range(len(rows)), by, color=C["metal"], edgecolor="white", width=0.66, label="named by the clade")
    ax.bar(range(len(rows)), other, bottom=by, color=C["organic"], edgecolor="white", width=0.66,
           label="named by another channel")
    _xt(ax, [g for g, _ in rows], rotation=35, ha="right", fontsize=FS_BODY)
    ax.set_ylabel("ion-resolved calls")
    ax.set_ylim(0, (by + other).max() * 1.45)
    ax.legend(loc="upper right", fontsize=FS_SMALL, handlelength=1.0)
    panel_head(ax, "b", "Who names the ion", letter_dx_pt=-26)

    save(fig, "figS4_ssn_clade")


def main() -> int:
    apply_style()
    N = numbers()
    fig_pwm(N)
    fig_ssn(N)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
