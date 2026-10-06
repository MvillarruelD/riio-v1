"""Contracts binding the figure scripts, the driver and figure QA together.

The failure these prevent is drift between three lists that must agree: what `run_pipeline.py` is
configured to RUN, what `figure_qa.py` CHECKS, and what the scripts actually produce. A figure
configured to run but not checked is a figure nobody is checking, and it reaches the report
unverified.
"""
from __future__ import annotations

import ast
import os

import pytest

from conftest import ANALYSIS

# `run_paths` refuses to import without a tag, and both modules under test import it.
os.environ.setdefault("TFOP_RUN_TAG", "contracttest")

import figure_qa  # noqa: E402
import fig_style  # noqa: E402
import run_pipeline  # noqa: E402

#: Scripts that run in the figure stage without owning a primary stem: one writes the table F23 is
#: drawn from, the other builds the explicitly non-primary F24 composite.
TABLE_AND_COMPOSITE = {"fig5_complement.py", "fig5_compose.py"}


def test_every_configured_figure_script_exists():
    missing = [f for f in run_pipeline.FIGURE_SCRIPTS if not (ANALYSIS / f).is_file()]
    assert not missing, f"configured but absent: {missing}"


def test_the_checked_set_and_the_configured_set_agree():
    configured = set(run_pipeline.FIGURE_SCRIPTS)
    producers = set(figure_qa.PRIMARY.values())
    assert producers | TABLE_AND_COMPOSITE == configured, (
        f"configured but not QA-checked: {sorted(configured - producers - TABLE_AND_COMPOSITE)}; "
        f"QA-checked but not configured: {sorted(producers - configured)}"
    )


def test_qa_agrees_with_itself_about_the_configured_set():
    """`figure_qa.check_configured_set()` is the runtime form of the test above."""
    assert figure_qa.check_configured_set() == []


def test_f24_is_excluded_from_the_primary_set():
    """It embeds raster artwork, so holding it to the editable-vector contract is wrong."""
    assert "F24_fourgenome_ABCDE_composite" in figure_qa.NON_PRIMARY_BY_DESIGN
    assert "F24_fourgenome_ABCDE_composite" not in figure_qa.PRIMARY


def test_the_export_formats_are_the_documented_four():
    assert figure_qa.FORMATS == ("svg", "pdf", "png", "tiff")


def test_the_width_limit_agrees_between_style_and_qa():
    assert fig_style.MAX_WIDTH_MM == figure_qa.MAX_WIDTH_MM == 183.0


def test_editable_text_and_truetype_rules_are_set():
    import matplotlib.pyplot as plt
    assert plt.rcParams["svg.fonttype"] == "none"
    assert plt.rcParams["pdf.fonttype"] == 42


# --------------------------------------------------------------------------- per-script rules
FIGURE_SCRIPTS = sorted(set(run_pipeline.FIGURE_SCRIPTS))


@pytest.mark.parametrize("script", FIGURE_SCRIPTS)
def test_no_figure_script_keeps_a_private_export(script):
    """Exports go through `fig_style.save()`; a private copy is a second contract to maintain."""
    text = (ANALYSIS / script).read_text(encoding="utf-8")
    offenders = [line.strip() for line in text.splitlines()
                 if ".savefig(" in line and not line.strip().startswith("#")]
    assert not offenders, (f"{script} calls savefig directly:\n  " + "\n  ".join(offenders)
                           + "\n  Use fig_style.save() so there is one export contract.")


@pytest.mark.parametrize("script", FIGURE_SCRIPTS)
def test_the_backend_is_chosen_before_pyplot_is_imported(script):
    """`fig_style` (or an explicit `use("Agg")`) must come first; otherwise it is a backend switch."""
    tree = ast.parse((ANALYSIS / script).read_text(encoding="utf-8"), filename=script)
    # `ast.walk` is breadth-first, NOT line-ordered, so every candidate is collected and the
    # EARLIEST line taken. Recording the first node walk happened to reach compared arbitrary
    # lines and reported correctly ordered files as violations.
    pyplot_lines: list[int] = []
    backend_lines: list[int] = []
    for node in ast.walk(tree):
        line = getattr(node, "lineno", None)
        if line is None:
            continue
        if isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name == "matplotlib.pyplot":
                    pyplot_lines.append(line)
                if alias.name == "fig_style":
                    backend_lines.append(line)
        elif isinstance(node, ast.ImportFrom):
            if node.module == "matplotlib" and any(a.name == "pyplot" for a in node.names):
                pyplot_lines.append(line)
            if node.module == "fig_style":
                backend_lines.append(line)
        elif isinstance(node, ast.Call):
            func = node.func
            if isinstance(func, ast.Attribute) and func.attr == "use":
                backend_lines.append(line)

    pyplot_at = min(pyplot_lines) if pyplot_lines else None
    backend_at = min(backend_lines) if backend_lines else None
    if pyplot_at is None:
        pytest.skip(f"{script} does not import pyplot directly")
    assert backend_at is not None, f"{script} imports pyplot without selecting a backend first"
    assert backend_at < pyplot_at, (
        f"{script}: the backend is selected at line {backend_at}, after pyplot is imported at "
        f"line {pyplot_at}. Import fig_style (or call matplotlib.use) first."
    )


@pytest.mark.parametrize("script", sorted(set(figure_qa.PRIMARY.values())))
def test_data_driven_figures_validate_their_inputs(script):
    """Every primary figure that reads run data must fail closed on missing input.

    `fig_pipeline.py` is exempt: F0/F1 are schematics of the stage graph with no data inputs.
    """
    if script == "fig_pipeline.py":
        pytest.skip("F0/F1 are schematics with no data inputs")
    text = (ANALYSIS / script).read_text(encoding="utf-8")
    assert "require(" in text or "SystemExit" in text or "return 1" in text, (
        f"{script} has no fail-closed guard; use fig_style.require() so a missing input names "
        f"itself and the stage that produces it"
    )


@pytest.mark.parametrize("script", FIGURE_SCRIPTS)
def test_production_figures_do_not_read_the_benchmark_tree(script):
    """The production/benchmark boundary, enforced at the figure layer."""
    text = (ANALYSIS / script).read_text(encoding="utf-8")
    reads_bench = 'benchmarking"' in text or "benchmarking/" in text
    assert not reads_bench, f"{script} reads the benchmark tree from the production figure stage"


def test_require_raises_with_the_producing_stage_named():
    with pytest.raises(SystemExit) as exc:
        fig_style.require(False, "the regulon summary", "the regulons stage")
    message = str(exc.value)
    assert "the regulon summary" in message and "the regulons stage" in message


def test_require_passes_through_when_the_input_is_present():
    assert fig_style.require(True, "x", "y") is None


def test_production_report_has_no_benchmark_dependency():
    source = (ANALYSIS / "build_report_v5_docx.py").read_text(encoding="utf-8")
    tree = ast.parse(source)
    imports = {alias.name for node in ast.walk(tree) if isinstance(node, ast.Import)
               for alias in node.names}
    assert not imports.intersection({"report_cross_run", "cross_run_tables"})
    assert "F21_regulon_limits" not in source
    assert "benchmarking/" not in source
