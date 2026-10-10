#!/usr/bin/env python
"""Run-scoped F0/F1 exports using the shared RIIO schematic renderer.

TFOP_RUN_TAG=<tag> python analysis/fig_pipeline.py
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

# fig_style selects Agg and owns the run-scoped export contract.
from fig_style import save
from reports.build_pipeline_figure import draw_pipeline_figure
import run_paths as RP


def graphical_abstract() -> None:
    """Preserve the overview artifact contract with the four supported input modes."""
    save(draw_pipeline_figure(compact=True), RP.FIGURES / "F0_pipeline_graphical_abstract")


def detailed_flowchart() -> None:
    """Use the same full renderer as the manuscript Figure 1 builder."""
    save(draw_pipeline_figure(), RP.FIGURES / "F1_pipeline_schematic")


def main() -> int:
    graphical_abstract()
    detailed_flowchart()
    print("  F0_pipeline_graphical_abstract.svg/.pdf/.png/.tiff")
    print("  F1_pipeline_schematic.svg/.pdf/.png/.tiff")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
