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


def test_pipeline_builders_share_the_renderer_and_preserve_stems(monkeypatch):
    import fig_pipeline
    from reports import build_pipeline_figure
    assert fig_pipeline.draw_pipeline_figure is build_pipeline_figure.draw_pipeline_figure
    exports = []
    monkeypatch.setattr(fig_pipeline, "save", lambda fig, path: exports.append((fig, path.name)))
    fig_pipeline.graphical_abstract()
    fig_pipeline.detailed_flowchart()
    assert [stem for _, stem in exports] == ["F0_pipeline_graphical_abstract", "F1_pipeline_schematic"]
    import matplotlib.pyplot as plt
    for fig, _ in exports:
        plt.close(fig)


@pytest.mark.parametrize("compact", [False, True])
def test_pipeline_renderer_contains_modes_exports_and_unclipped_text(compact):
    from reports.build_pipeline_figure import draw_pipeline_figure
    import matplotlib.pyplot as plt
    fig = draw_pipeline_figure(compact=compact)
    try:
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        labels = {text.get_text() for ax in fig.axes for text in ax.texts}
        assert {"Single protein", "Protein list", "Genome accession", "Local genome",
                "Candidate discovery", "Per-TF inference",
                "Per-TF bundle", "Run report", "AlphaFold 3 hand-off"} <= labels
        assert "CLI · Python API · browser GUI" in labels
        # A compact schematic carries states; interface syntax stays in the legend.
        assert {"Homologue promoters", "Motif", "Genome rescan", "Ranked operators",
                "Resolve context", "Fused inducer call"} <= labels
        assert len(fig.axes[0].patches) > 40  # editable state drawings, not embedded images
        assert not fig.axes[0].images
        canvas = fig.bbox
        for ax in fig.axes:
            for text in ax.texts:
                bbox = text.get_window_extent(renderer)
                assert bbox.x0 >= canvas.x0 and bbox.x1 <= canvas.x1, text.get_text()
                assert bbox.y0 >= canvas.y0 and bbox.y1 <= canvas.y1, text.get_text()
                assert text.get_fontsize() >= 7, text.get_text()
        assert fig.get_figwidth() * 25.4 <= 183
        assert fig.get_figheight() * 25.4 <= 130
        # The TIFF saver adds matplotlib's default padding around tight content.
        # Reserve it in the canvas contract so the exported file also fits the page.
        assert (fig.get_figwidth() + 2 * plt.rcParams["savefig.pad_inches"]) * 25.4 <= 183
    finally:
        plt.close(fig)


def test_pipeline_schematic_and_legend_do_not_claim_unverified_neutrality():
    from reports.build_legends import legend_fig1
    from reports.build_pipeline_figure import draw_pipeline_figure
    import matplotlib.pyplot as plt
    fig = draw_pipeline_figure()
    try:
        text = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
        legend = legend_fig1({})  # schematic must not depend on run numbers
        for content in (text, legend):
            assert "unchanged from v1.0" not in content
            assert "never alter" not in content
            assert "never read each other's" not in content
            assert "changes no call" not in content
            assert "default" in content
        assert "Single-protein calls enter prediction directly" in legend
        assert "input and policy fingerprints" in legend
        assert "run and discovery reports" in legend
        assert "without boosting their ranking" in legend
        assert "drawings illustrate states" in legend
    finally:
        plt.close(fig)
