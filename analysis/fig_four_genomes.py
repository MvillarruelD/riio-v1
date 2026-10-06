#!/usr/bin/env python
"""F22 -- the improved regulon reconstruction applied to the paper's four genomes.

Panels:
  a  candidate TFs per genome and family (the paper's Table S2 repertoire)
  b  predicted regulon size per TF, by family
  c  are the predictions TF-SPECIFIC? pairwise gene-set overlap within each genome
  d  operon structure recovered -- genes per operon, which is what the operon-assembly fix changed

Reads the current run's `regulon_v2_summary.tsv` (written by regenerate_regulons.py); the run is
selected by `TFOP_RUN_TAG` -- see `run_paths.py`.

    $TFOP_PY fig_four_genomes.py
"""
import csv
import itertools
import statistics as st
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
# fig_style selects the Agg backend, so it is imported BEFORE pyplot.
from fig_style import INK, NEUTRAL_L, NEUTRAL_M, PANEL, FAM, require, save   # noqa: E402
import matplotlib.pyplot as plt                                     # noqa: E402

# the report chain's output locations, so a figure is never drawn from an archived run
import run_paths as RP                                              # noqa: E402


def _manifest_path():
    """The manifest for THIS run. Defined once, in run_paths."""
    return RP.MANIFEST

SUMMARY = RP.REGULON_SUMMARY
MANIFEST = _manifest_path()
OUT = RP.FIGURES / "F22_four_genomes"

GENOME_LABEL = {
    "GCF_000195955.2": "M. tuberculosis H37Rv",
    "GCF_004345205.2": "M. avium hominissuis",
    "GCA_004345205.2": "M. avium hominissuis",
    "GCF_008369605.1": "V. cholerae RFB16",
    "GCF_002224265.1": "V. vulnificus NBRC 15645",
}
SHORT = {                                   # panel c has four categories on one axis -- keep them short
    "M. tuberculosis H37Rv": "M. tub.",
    "M. avium hominissuis": "M. avium",
    "V. cholerae RFB16": "V. chol.",
    "V. vulnificus NBRC 15645": "V. vuln.",
}
FAMC = dict(FAM)
FAMC.setdefault("(auto)", "#767676")
OTHER = "#B0B0B8"


def rows():
    return list(csv.DictReader(SUMMARY.open(encoding="utf-8"), delimiter="\t")) if SUMMARY.is_file() else []


def panel_a(ax):
    """The candidate repertoire actually processed, per genome and family."""
    rs = rows()
    by = defaultdict(Counter)
    for r in rs:
        by[GENOME_LABEL.get(r["genome"], r["genome"])][r["family"] or "(auto)"] += 1
    gens = sorted(by, key=lambda g: -sum(by[g].values()))
    fams = ["ArsR/SmtB", "MerR", "Fur", "(auto)"]
    y = np.arange(len(gens))
    left = np.zeros(len(gens))
    for f in fams:
        v = np.array([by[g].get(f, 0) for g in gens], float)
        ax.barh(y, v, left=left, height=0.66, label=f if f != "(auto)" else "other families",
                color=FAMC.get(f, OTHER))
        left += v
    ax.set_yticks(y); ax.set_yticklabels(gens, fontsize=6.2)
    ax.invert_yaxis()
    ax.set_xlabel("candidate TFs with a regenerated regulon")
    ax.legend(fontsize=6, loc="lower right")
    n_own = sum(1 for r in rs if r.get("pwm_source", "").startswith("per_tf"))
    ax.set_title(f"a   {len(rs)}/150 TFs — {n_own} on their own motif, {len(rs)-n_own} on a family prior",
                 loc="left", fontsize=7.5, color=INK)


def panel_b(ax):
    """Regulon size split by MOTIF SOURCE -- the two tiers must never be read as equivalent."""
    rs = rows()
    own = [int(r["n_genes"]) for r in rs if r.get("pwm_source", "").startswith("per_tf")]
    fam = [int(r["n_genes"]) for r in rs if r.get("pwm_source", "").startswith("family_prior")]
    data = [d for d in (own, fam) if d]
    labs = [l for l, d in zip([f"per-TF motif\n(n={len(own)})", f"family prior\n(n={len(fam)})"],
                              (own, fam)) if d]
    bp = ax.boxplot(data, vert=True, widths=0.5, patch_artist=True,
                    medianprops=dict(color=INK, lw=1.1),
                    flierprops=dict(marker="o", ms=2, mfc=NEUTRAL_M, mec="none", alpha=0.6))
    for patch, c in zip(bp["boxes"], ["#0F4D92", "#B0B0B8"]):
        patch.set_facecolor(c); patch.set_alpha(0.55); patch.set_edgecolor(INK); patch.set_linewidth(0.7)
    ax.set_xticklabels(labs, fontsize=6.2)
    ax.set_ylabel("predicted regulon genes per TF")
    ax.set_title("b   Two motif tiers — never pooled", loc="left", fontsize=7.5, color=INK)
    n_cap = sum(1 for r in rs if int(r["n_operons"]) >= 25)
    ax.text(0.98, 0.03, f"{n_cap}/{len(rs)} TFs hit the 25-operon cap\n"
                        f"→ ranked candidates, not FDR-significant",
            transform=ax.transAxes, ha="right", fontsize=6.0, color=NEUTRAL_M, style="italic")


def panel_c(ax):
    """Specificity: do different TFs in one genome predict DIFFERENT regulons?

    Computed on PER-TF-MOTIF rows only. Two TFs sharing a family-prior PWM in the same genome are scanning
    the identical motif against the identical sequence, so their regulons are identical *by construction* --
    including them would manufacture an overlap signal that says nothing about the tool.
    """
    rs = [r for r in rows() if r.get("pwm_source", "").startswith("per_tf")]
    by = defaultdict(list)
    for r in rs:
        if r["genes"]:
            by[GENOME_LABEL.get(r["genome"], r["genome"])].append(set(r["genes"].split(";")))
    gens = sorted(by, key=lambda g: -len(by[g]))
    vals, labs = [], []
    for g in gens:
        js = [len(a & b) / len(a | b) for a, b in itertools.combinations(by[g], 2) if a | b]
        if js:
            vals.append(js); labs.append(f"{SHORT.get(g, g)}\n(n={len(by[g])})")
    if not vals:
        ax.axis("off"); return
    parts = ax.violinplot(vals, showmedians=True, widths=0.7)
    # a handful of TF pairs share an identical gene set (Jaccard 1.0) and stretch the axis; clip the view
    # to the bulk of the distribution rather than dropping the points
    hi = max(np.percentile(v, 97) for v in vals)
    ax.set_ylim(0, max(0.25, min(1.0, hi * 1.4)))
    for pc in parts["bodies"]:
        pc.set_facecolor("#0F4D92"); pc.set_alpha(0.45); pc.set_edgecolor(INK); pc.set_linewidth(0.6)
    for k in ("cmedians", "cbars", "cmins", "cmaxes"):
        if k in parts:
            parts[k].set_color(INK); parts[k].set_linewidth(0.8)
    ax.set_xticks(range(1, len(labs) + 1)); ax.set_xticklabels(labs, fontsize=5.6)
    ax.set_ylabel("pairwise gene-set Jaccard")
    med = st.mean([st.mean(v) for v in vals])
    ax.set_title(f"c   TF-specific where a per-TF motif exists — mean overlap {med:.03f}",
                 loc="left", fontsize=7.5, color=INK)


def panel_d(ax):
    """Operon structure: genes per operon is exactly what the assembly fix changed."""
    rs = rows()
    per = [int(r["n_genes"]) / max(1, int(r["n_operons"])) for r in rs]
    if not per:
        ax.axis("off"); return
    ax.hist(per, bins=18, color="#42949E", edgecolor="white", linewidth=0.6, zorder=3)
    m = st.median(per)
    ax.axvline(m, color=INK, lw=1.0, ls="--", zorder=4)
    ax.text(m, ax.get_ylim()[1] * 0.92, f"  median {m:.2f}", fontsize=6.2, color=INK)
    ax.axvline(1.0, color=NEUTRAL_L, lw=1.0, zorder=2)
    ax.text(1.0, ax.get_ylim()[1] * 0.62, " 1.0 = no operon\n structure recovered",
            fontsize=6, color=NEUTRAL_M)
    ax.set_xlabel("genes per predicted operon")
    ax.set_ylabel("TFs")
    ax.set_title("d   Polycistronic targets are now expanded", loc="left", fontsize=7.5, color=INK)


def main():
    # Previously an early `return`, i.e. exit 0 -- the driver was told the figure stage had
    # succeeded while F22 was never drawn, and on a dirty output directory the report then picked
    # up the PREVIOUS run's F22 and looked complete.
    require(SUMMARY.is_file(), f"the regulon summary ({SUMMARY.name})",
            "the regulons stage (regenerate_regulons.py)")
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.2))
    for ax in axes.ravel():
        ax.set_facecolor("white")
        ax.grid(axis="both", color=PANEL, lw=0.7, zorder=0)
        ax.set_axisbelow(True)
    panel_a(axes[0, 0]); panel_b(axes[0, 1]); panel_c(axes[1, 0]); panel_d(axes[1, 1])
    fig.suptitle("Improved regulon reconstruction applied to the four genomes of the SSN paper",
                 fontsize=8.5, color=INK, x=0.005, ha="left", y=1.005, weight="bold")
    fig.tight_layout(h_pad=2.6, w_pad=2.2)
    for f in save(fig, OUT):
        print("wrote", f)


if __name__ == "__main__":
    main()
