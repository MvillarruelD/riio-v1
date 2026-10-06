#!/usr/bin/env python
"""Publication schematics for the genome-to-report prediction workflow.

The two figures use one visual vocabulary at two levels of detail:

* ``F0_pipeline_graphical_abstract`` is the four-step overview.
* ``F1_pipeline_schematic`` separates batch discovery from the three per-TF evidence routes.

Both diagrams describe the current executable workflow. They deliberately avoid performance claims,
run-specific counts and implied feedback from predicted regulons into motif or inducer inference.

    TFOP_RUN_TAG=<tag> python analysis/fig_pipeline.py
"""
from __future__ import annotations

import sys
from importlib import import_module
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

# fig_style selects Agg before pyplot and owns the publication export contract.
from fig_style import save                                 # noqa: E402
import matplotlib as mpl                                   # noqa: E402
import matplotlib.pyplot as plt                            # noqa: E402
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch  # noqa: E402

RP = import_module("run_paths")


# --------------------------------------------------------------------------- visual system
INK = "#17211E"
MUTED = "#66716C"
LINE = "#BFCAC5"
PAPER = "#FFFFFF"
DISC = "#2F6B56"
DISC_SOFT = "#E5EFEA"
PRED = "#315E7D"
PRED_SOFT = "#E7EEF3"
OPER = "#4A8C8A"
OPER_SOFT = "#E5F0EF"
OUTPUT = "#A9612E"
OUTPUT_SOFT = "#F4E9DF"
STRUCT = "#747C78"
PANEL = "#F5F7F6"

mpl.rcParams.update({
    "text.color": INK,
    "figure.facecolor": PAPER,
    "savefig.facecolor": PAPER,
})


def _export(fig, stem: str) -> list[str]:
    """Write editable SVG/PDF and 600-dpi PNG/TIFF through the shared contract."""
    return save(fig, RP.FIGURES / stem)


def _canvas(figsize):
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.set_aspect("auto")
    ax.axis("off")
    return fig, ax


def _round_box(ax, x, y, w, h, *, fc=PAPER, ec=LINE, radius=1.4, lw=0.8, z=2):
    patch = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0.01,rounding_size={radius}",
        fc=fc, ec=ec, lw=lw, zorder=z,
    )
    ax.add_patch(patch)
    return patch


def _arrow(ax, start, end, *, color=LINE, lw=1.2, dashed=False, z=1,
           connection="arc3,rad=0", scale=8):
    ax.add_patch(FancyArrowPatch(
        start, end, arrowstyle="-|>", mutation_scale=scale,
        lw=lw, color=color, linestyle=(0, (3, 2.4)) if dashed else "-",
        connectionstyle=connection, shrinkA=2, shrinkB=2, zorder=z,
    ))


def _node(ax, x, y, w, h, title, subtitle="", *, color=PRED, fill=PAPER,
          title_size=6.2, subtitle_size=4.9, align="center"):
    _round_box(ax, x, y, w, h, fc=fill, ec=color, radius=1.2, lw=0.9, z=3)
    tx = x + w / 2 if align == "center" else x + 1.4
    ha = "center" if align == "center" else "left"
    title_y = y + h * (0.62 if subtitle else 0.50)
    ax.text(tx, title_y, title, ha=ha, va="center", fontsize=title_size,
            fontweight="bold", color=INK, linespacing=1.08, zorder=4)
    if subtitle:
        ax.text(tx, y + h * 0.28, subtitle, ha=ha, va="center", fontsize=subtitle_size,
                color=MUTED, linespacing=1.16, zorder=4)


def _section_label(ax, x, y, letter, title, frequency, color):
    ax.text(x, y, letter, ha="left", va="center", fontsize=7.8,
            fontweight="bold", color=color)
    ax.text(x + 2.3, y, title, ha="left", va="center", fontsize=7.8,
            fontweight="bold", color=INK)
    ax.text(97, y, frequency, ha="right", va="center", fontsize=5.5,
            color=MUTED, style="italic")


def _pill(ax, x, y, w, label, *, color=PRED, size=4.7):
    _round_box(ax, x, y, w, 5.2, fc=PAPER, ec=color, radius=1.0, lw=0.75, z=4)
    ax.text(x + w / 2, y + 2.6, label, ha="center", va="center",
            fontsize=size, color=INK, zorder=5)


def _genome_icon(ax, cx, cy, color):
    ax.add_patch(Circle((cx, cy), 5.7, fc=PAPER, ec=color, lw=1.1, zorder=4))
    for yy, width in ((cy + 2.1, 6.5), (cy, 8.1), (cy - 2.1, 5.3)):
        ax.plot([cx - width / 2, cx + width / 2], [yy, yy], color=color,
                lw=1.45, solid_capstyle="round", zorder=5)


def _document_icon(ax, x, y, color):
    _round_box(ax, x, y, 8.0, 13.0, fc=PAPER, ec=color, radius=.8, lw=1.0, z=4)
    for yy, width in ((y + 9.7, 4.7), (y + 6.8, 5.5), (y + 3.9, 3.8)):
        ax.plot([x + 1.3, x + 1.3 + width], [yy, yy], color=color,
                lw=1.25, solid_capstyle="round", zorder=5)


# --------------------------------------------------------------------------- graphical abstract
def graphical_abstract() -> None:
    """Four-step overview for a graphical-abstract slot."""
    fig, ax = _canvas((7.2, 3.05))

    ax.text(4, 93, "GENOME-TO-REPORT WORKFLOW", fontsize=6.3, fontweight="bold",
            color=DISC, va="center")
    ax.text(4, 84, "From annotated genomes to traceable TF predictions",
            fontsize=13.0, fontweight="bold", color=INK, va="center")
    ax.text(4, 76.5,
            "Candidate discovery is followed by independent evidence routes for each regulator.",
            fontsize=6.1, color=MUTED, va="center")

    cards = [
        (3, 25, 20, 43, DISC_SOFT, DISC, "01", "Discover TFs",
         "genome + annotation"),
        (27.5, 25, 20, 43, PRED_SOFT, PRED, "02", "Resolve context",
         "family · clade · locus"),
        (52, 25, 20, 43, OPER_SOFT, OPER, "03", "Infer hypotheses",
         "operators · inducer · regulon"),
        (76.5, 25, 20.5, 43, OUTPUT_SOFT, OUTPUT, "04", "Package evidence",
         "ranked calls + run record"),
    ]
    for x, y, w, h, fc, color, number, title, subtitle in cards:
        _round_box(ax, x, y, w, h, fc=fc, ec="none", radius=2.1, lw=0, z=0)
        ax.text(x + 2.0, y + h - 4.7, number, fontsize=5.6,
                fontweight="bold", color=color, va="center")
        ax.text(x + w / 2, y + 11.0, title, fontsize=7.6, fontweight="bold",
                color=INK, ha="center", va="center")
        ax.text(x + w / 2, y + 5.4, subtitle, fontsize=5.6, color=MUTED,
                ha="center", va="center")

    _genome_icon(ax, 13, 52.3, DISC)

    for label, yy in (("family", 57.0), ("SSN clade", 50.5), ("TF locus", 44.0)):
        _pill(ax, 32.0, yy, 11.0, label, color=PRED, size=4.7)

    # Operator/regulon and inducer inference stay independent until reporting.
    for label, yy in (("operator → regulon", 56.2), ("inducer", 46.6)):
        ax.text(55.0, yy, label, fontsize=5.0, color=INK, va="center")
        ax.plot([64.2, 67.5], [yy, 52.0], color=OPER, lw=0.9, zorder=3)
    ax.add_patch(Circle((68.2, 52.0), 2.2, fc=PAPER, ec=OPER, lw=1.0, zorder=4))
    ax.text(68.2, 52.0, "+", ha="center", va="center", fontsize=7.0,
            fontweight="bold", color=OPER, zorder=5)

    _document_icon(ax, 82.8, 45.5, OUTPUT)

    for start, end, color in (
        ((23.2, 46.5), (27.1, 46.5), DISC),
        ((47.7, 46.5), (51.6, 46.5), PRED),
        ((72.2, 46.5), (76.1, 46.5), OPER),
    ):
        _arrow(ax, start, end, color=color, lw=1.7, scale=10)

    ax.text(50, 13.5,
            "Predictions remain linked to inputs, evidence, parameters and software.",
            fontsize=6.8, color=INK, ha="center", va="center", fontweight="bold")
    ax.text(50, 7.1, "Evidence sources may abstain; all reported calls remain unverified predictions.",
            fontsize=5.5, color=MUTED, ha="center", va="center")

    _export(fig, "F0_pipeline_graphical_abstract")


# --------------------------------------------------------------------------- detailed main-text figure
def detailed_flowchart() -> None:
    """Main-text schematic with batch and per-regulator execution separated."""
    fig, ax = _canvas((7.2, 5.05))

    ax.text(3, 97, "Genome-to-report prediction workflow", fontsize=11.5,
            fontweight="bold", color=INK, va="center")
    ax.text(3, 93.2,
            "Batch candidate discovery feeds independent, traceable predictions for each TF",
            fontsize=5.9, color=MUTED, va="center")

    # a. Candidate discovery: the analysis-layer batch front end.
    _round_box(ax, 2, 75.0, 96, 14.0, fc=DISC_SOFT, ec="none", radius=1.7, lw=0, z=0)
    _section_label(ax, 3.2, 87.0, "a", "Candidate discovery", "once per genome", DISC)
    discovery = [
        (4.5, "Genome + GFF", "source sequence"),
        (28.2, "Translated proteins", "retain accessions"),
        (51.9, "BLAST/HMM profiles", "find TF candidates"),
        (75.6, "Candidate manifest", "one row per TF"),
    ]
    for i, (x, title, subtitle) in enumerate(discovery):
        _node(ax, x, 77.0, 19.7, 7.0, title, subtitle, color=DISC,
              title_size=5.9, subtitle_size=4.5)
        if i < len(discovery) - 1:
            _arrow(ax, (x + 19.9, 80.5), (discovery[i + 1][0] - .2, 80.5),
                   color=DISC, lw=1.15, scale=8)

    # b. Per-regulator prediction: current predictor.pipeline.run_novel logic.
    _round_box(ax, 2, 22.0, 96, 49.0, fc=PANEL, ec="none", radius=1.7, lw=0, z=0)
    _section_label(ax, 3.2, 68.7, "b", "Per-regulator prediction", "one independent run per TF", PRED)

    _node(ax, 4.2, 40.0, 12.5, 16.0, "TF sequence\n+ source genome", "declared inputs",
          color=DISC, title_size=5.9, subtitle_size=4.5)
    _node(ax, 20.2, 40.0, 15.2, 16.0, "Resolve context", "family · SSN clade\nTF locus + promoters",
          color=PRED, title_size=5.9, subtitle_size=4.5)
    _arrow(ax, (16.9, 48.0), (19.9, 48.0), color=PRED, lw=1.2, scale=8)

    ax.text(39.0, 65.2, "OPTIONAL STRUCTURE ROUTE", fontsize=4.8,
            color=STRUCT, fontweight="bold", va="center")
    _node(ax, 39.0, 56.0, 42.0, 7.8, "Optional apo-dimer structure",
          "fold/QC; may inform operator design", color=STRUCT,
          title_size=5.4, subtitle_size=4.7)
    _arrow(ax, (35.7, 53.0), (38.7, 59.4), color=STRUCT, lw=.9,
           dashed=True, connection="arc3,rad=-.16", scale=7)

    ax.text(39.0, 52.4, "OPERATOR + REGULON ROUTE", fontsize=5.0,
            color=OPER, fontweight="bold", va="center")
    operator_nodes = [
        (39.0, "Homolog promoters", "SSN / MSA"),
        (54.0, "Seed motifs", "palindrome prior"),
        (69.0, "Genome rescans", "width trials\n→ selected sites"),
    ]
    for i, (x, title, subtitle) in enumerate(operator_nodes):
        _node(ax, x, 42.0, 12.5, 8.4, title, subtitle, color=OPER,
              title_size=5.15, subtitle_size=4.45)
        if i < len(operator_nodes) - 1:
            _arrow(ax, (x + 12.7, 46.2), (operator_nodes[i + 1][0] - .2, 46.2),
                   color=OPER, lw=1.0, scale=7)
    _arrow(ax, (35.7, 48.0), (38.7, 46.2), color=OPER, lw=1.0, scale=7)

    ax.text(39.0, 37.7, "INDUCER ROUTE", fontsize=5.0,
            color=PRED, fontweight="bold", va="center")
    _node(ax, 39.0, 27.5, 30.8, 8.5, "Independent inducer evidence",
          "SSN · coordination · MetalNet\nLigify · neighbourhood", color=PRED,
          title_size=5.3, subtitle_size=4.45)
    _node(ax, 72.3, 27.5, 8.7, 8.5, "Consensus", "sources may\nabstain",
          color=PRED, title_size=4.9, subtitle_size=4.4)
    _arrow(ax, (35.7, 43.0), (38.7, 32.0), color=PRED, lw=.9,
           connection="arc3,rad=.16", scale=7)
    _arrow(ax, (70.0, 31.7), (72.0, 31.7), color=PRED, lw=1.0, scale=7)

    # One deliberately dominant endpoint: the report directory is the audit unit.
    _round_box(ax, 84.0, 27.5, 12.5, 36.3, fc=OUTPUT_SOFT, ec=OUTPUT,
               radius=1.4, lw=1.0, z=3)
    ax.text(90.25, 59.8, "SELF-CONTAINED", fontsize=4.5, color=OUTPUT,
            ha="center", va="center", fontweight="bold")
    ax.text(90.25, 55.7, "TF bundle", fontsize=6.3, color=INK,
            ha="center", va="center", fontweight="bold")
    for label, yy in (
        ("motif + sites", 49.8),
        ("candidate regulon", 45.5),
        ("inducer call", 41.2),
        ("apo structure", 36.9),
        ("AF3 hand-off jobs", 34.1),
        ("inputs + run record", 31.3),
    ):
        ax.text(90.25, yy, label, fontsize=4.55, color=INK, ha="center", va="center")

    _arrow(ax, (81.7, 46.2), (83.7, 48.5), color=OPER, lw=1.0, scale=7)
    _arrow(ax, (81.3, 31.7), (83.7, 37.5), color=PRED, lw=1.0,
           connection="arc3,rad=-.12", scale=7)
    _arrow(ax, (81.3, 59.8), (83.7, 56.0), color=STRUCT, lw=.9,
           dashed=True, connection="arc3,rad=.10", scale=7)

    ax.text(4.2, 24.6,
            "No feedback: predicted regulons do not raise motif or inducer confidence.",
            fontsize=4.7, color=MUTED, va="center", style="italic")

    # c. Batch synthesis: only after all per-TF bundles are complete.
    _round_box(ax, 2, 3.5, 96, 14.0, fc=OUTPUT_SOFT, ec="none", radius=1.7, lw=0, z=0)
    _section_label(ax, 3.2, 15.5, "c", "Batch synthesis", "after all TF bundles complete", OUTPUT)
    outputs = [
        (4.5, "Per-TF bundles", "validated file sets"),
        (36.0, "Aggregate + regulon tables", "one row per regulator"),
        (67.5, "Figures + author report", "publish share package"),
    ]
    for i, (x, title, subtitle) in enumerate(outputs):
        _node(ax, x, 5.3, 27.8, 6.8, title, subtitle, color=OUTPUT,
              title_size=5.6, subtitle_size=4.6)
        if i < len(outputs) - 1:
            _arrow(ax, (x + 28.0, 8.7), (outputs[i + 1][0] - .2, 8.7),
                   color=OUTPUT, lw=1.05, scale=7)

    _export(fig, "F1_pipeline_schematic")


def main() -> int:
    graphical_abstract()
    detailed_flowchart()
    print("  F0_pipeline_graphical_abstract.svg/.pdf/.png/.tiff")
    print("  F1_pipeline_schematic.svg/.pdf/.png/.tiff")
    return 0


if __name__ == "__main__":
    sys.exit(main())
