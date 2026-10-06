#!/usr/bin/env python
"""Fail closed when primary publication figures violate the shared export contract.

    TFOP_RUN_TAG=<tag> $TFOP_PY analysis/figure_qa.py

A failing report is a production stop condition: `run_pipeline.py` runs this after the figure stage
and refuses to build the report if it fails.

What it is actually defending against. A figure that is merely ABSENT is easy to notice. The
expensive failures are the ones that leave a plausible file on disk:

* a figure left over from an EARLIER RUN, which looks complete and is a different experiment;
* a figure whose source table changed after it was drawn, so the panel and the number in the text
  disagree;
* a raster that is 600 dpi in one dimension and 96 in the other;
* an SVG whose text was converted to paths, which cannot be edited at proof stage;
* a figure wider than the journal's column, discovered after submission.

So the checks below are not only "does the file exist": they compare the figure set against the
CONFIGURED production scripts, require every file to be newer than its inputs and to have been
written during this invocation, and measure the actual page geometry.
"""
from __future__ import annotations

import json
import os
import re
import sys
import time
from pathlib import Path

from PIL import Image

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_paths as RP  # noqa: E402

#: Every stem the production figure set must contain, and the script that produces it. This is the
#: same mapping `run_pipeline.FIGURE_SCRIPTS` drives, and `test_figure_contract.py` asserts the two
#: agree -- a figure configured to run but not checked here is a figure nobody is checking.
PRIMARY: dict[str, str] = {
    "F0_pipeline_graphical_abstract": "fig_pipeline.py",
    "F1_pipeline_schematic": "fig_pipeline.py",
    "F2_metal_corroboration": "fig_corroboration.py",
    "F3_operator_cds_distance": "fig_operator_cds.py",
    "F4a_operator_logo_per_family": "fig_logos.py",
    "F4b_cluster_motif_sharing": "fig_cluster_motifs.py",
    "F4c_operator_logo_per_tf": "fig_logos.py",
    "F6_evidence_ledger": "fig_evidence_ledger.py",
    "F22_four_genomes": "fig_four_genomes.py",
    "F23_panelE_inducer_per_regulator": "fig5_panelE.py",
}

#: `F24_fourgenome_ABCDE_composite` is deliberately NOT here. It embeds raster artwork supplied by
#: the authors, so it is a layout mock-up rather than a submission master; holding it to the
#: editable-vector contract would either fail forever or licence a raster into the primary set.
NON_PRIMARY_BY_DESIGN = {"F24_fourgenome_ABCDE_composite"}

FORMATS = ("svg", "pdf", "png", "tiff")
MIN_RASTER_PX = 1_000
MIN_DPI = 590
#: Nature-family double-column width. A figure wider than this is rejected unless the contract
#: documents it as a full-width supplement.
MAX_WIDTH_MM = 183.0
MM_PER_INCH = 25.4
#: Files older than this at the start of the run are leftovers, not products of it.
FRESH_WINDOW_S = float(os.environ.get("AN_QA_FRESH_WINDOW", 6 * 3600))

#: Files the FIGURE STAGE itself writes into the run's results directory. They are products of the
#: stage, not inputs to it, and counting them as inputs makes every figure drawn before them look
#: stale -- `fig5_complement.py` runs late in the stage and writes the per-regulator table, so F6,
#: F21 and F22 were all reported as older than "an input" that did not exist when they were drawn.
STAGE_PRODUCTS = frozenset({
    "figures_manifest.json", "figure_provenance.json",
    Path(RP.COMPLEMENT_TSV).name,
})

#: Imported, not redefined. This gate and fig_style.save()'s early warning have to be the SAME
#: measurement or the warning is not a warning -- when save() measured the figure object via
#: get_tightbbox and this measured the file, F4b passed there at "within limit" and failed here at
#: 184.8 mm on 2026-09-04.
from fig_style import svg_size_mm as _svg_size_mm  # noqa: E402


def _newest_input_mtime(since: float = 0.0) -> float:
    """When the run's inputs last changed.

    A figure older than its own source is a figure that does not show the current numbers. The
    comparison is deliberately coarse -- the run's tables and bundle manifests -- because that is
    the granularity at which a figure can actually go stale.

    `since` is when the figure STAGE began. Anything under these roots written at or after that
    moment was written BY the stage, so it is a product, not an input.

    Without it this check cannot pass. Figure scripts write side-car stats next to the tables they
    read (`cluster_motif_similarity.csv`, `F4b_cluster_motif_stats.json`,
    `F3_operator_cds_stats.json`, `evidence_ledger.tsv` ...), and `STAGE_PRODUCTS` names only three
    files. So the LAST figure to write a side-car retroactively made every figure generated before it
    "older than the run's inputs" -- 8 of 11 failed that way on 2026-09-04, and the three that passed
    did so only because they happened to run last. Enumerating side-cars in STAGE_PRODUCTS was
    rejected as the fix: it silently re-breaks the moment a new figure writes one.
    """
    newest = 0.0
    for root, patterns in ((RP.RESULTS, ("*.csv", "*.tsv", "*.json")),
                           (RP.REGULON, ("*.tsv", "*.csv")),
                           (RP.JOBS, ("*/bundle_manifest.json",))):
        if not root.is_dir():
            continue
        for pattern in patterns:
            for path in root.glob(pattern):
                if path.name in STAGE_PRODUCTS:
                    continue        # written by the figure stage itself; not an input to it
                try:
                    mtime = path.stat().st_mtime
                except OSError:
                    continue
                if since and mtime >= since:
                    continue        # written DURING the stage -> a product, not an input
                newest = max(newest, mtime)
    return newest


def check_stem(stem: str, *, newest_input: float, started: float,
               check_freshness: bool = True) -> list[str]:
    problems: list[str] = []
    for suffix in FORMATS:
        path = RP.FIGURES / f"{stem}.{suffix}"
        if not path.is_file() or path.stat().st_size == 0:
            problems.append(f"missing or empty {path.name}")
    if problems:
        return problems                      # geometry checks below need the files to exist

    for suffix in FORMATS:
        path = RP.FIGURES / f"{stem}.{suffix}"
        mtime = path.stat().st_mtime
        if newest_input and mtime < newest_input:
            problems.append(f"{path.name}: older than the run's inputs "
                            f"({time.strftime('%H:%M:%S', time.localtime(mtime))} < "
                            f"{time.strftime('%H:%M:%S', time.localtime(newest_input))}) -- "
                            "it does not show the current numbers")
        if check_freshness and started - mtime > FRESH_WINDOW_S:
            problems.append(f"{path.name}: not written by this run (last modified "
                            f"{(started - mtime) / 3600:.1f} h ago) -- a leftover from an earlier "
                            "run must never be presented as this run's figure")

    svg = RP.FIGURES / f"{stem}.svg"
    raw = svg.read_text(encoding="utf-8", errors="replace")
    if not re.search(r"<text(?:\s|>)", raw):
        problems.append(f"{svg.name}: text is not editable (converted to paths)")
    width_mm, height_mm = _svg_size_mm(raw)
    if width_mm is None or height_mm is None:
        problems.append(f"{svg.name}: no usable width/height on the <svg> element")
    else:
        if width_mm > MAX_WIDTH_MM + 0.5:
            problems.append(f"{svg.name}: {width_mm:.1f} mm wide, over the {MAX_WIDTH_MM:.0f} mm "
                            "publication limit")
        if height_mm <= 0:
            problems.append(f"{svg.name}: non-positive height")

    pdf = RP.FIGURES / f"{stem}.pdf"
    if pdf.read_bytes()[:5] != b"%PDF-":
        problems.append(f"{pdf.name}: not a PDF")

    for suffix in ("png", "tiff"):
        path = RP.FIGURES / f"{stem}.{suffix}"
        with Image.open(path) as im:
            width, height = im.size
            # BOTH dimensions: a panel exported at 600 dpi horizontally and 96 vertically passes a
            # `max(size)` test and is unusable.
            if min(width, height) < MIN_RASTER_PX:
                problems.append(f"{path.name}: raster is {width}x{height} px; both dimensions must "
                                f"be at least {MIN_RASTER_PX}")
            dpi = im.info.get("dpi")
            if not dpi or len(dpi) < 2:
                problems.append(f"{path.name}: no DPI metadata recorded")
            elif min(float(x) for x in dpi[:2]) < MIN_DPI:
                problems.append(f"{path.name}: recorded resolution is {dpi}, expected 600 dpi")
            if suffix == "tiff":
                compression = im.info.get("compression")
                if compression != "tiff_lzw":
                    problems.append(f"{path.name}: compression is {compression!r}, expected "
                                    "'tiff_lzw'")
    return problems


def check_configured_set() -> list[str]:
    """The checked set and the configured production set must be exactly the same figures."""
    problems: list[str] = []
    try:
        sys.path.insert(0, str(HERE))
        import run_pipeline
        configured = set(run_pipeline.FIGURE_SCRIPTS)
    except Exception as exc:                 # a driver that cannot be imported is its own problem
        return [f"cannot read the configured figure set from run_pipeline.py: {exc}"]

    producers = set(PRIMARY.values())
    # `fig5_complement.py` writes the table F23 is drawn from and `fig5_compose.py` builds the
    # non-primary composite; neither owns a primary stem, so neither is expected here.
    table_only = {"fig5_complement.py", "fig5_compose.py"}
    unchecked = configured - producers - table_only
    unknown = producers - configured
    if unchecked:
        problems.append(f"configured but not QA-checked: {', '.join(sorted(unchecked))}")
    if unknown:
        problems.append(f"QA-checked but not configured to run: {', '.join(sorted(unknown))}")
    return problems


def write_provenance(results: dict[str, list[str]], *, started: float) -> Path | None:
    """Record what was checked, from what, and with what result.

    A figure without provenance cannot be re-derived, and at proof stage the question is always
    "which run produced this panel?".
    """
    payload = {
        "run_tag": RP.TAG,
        "checked_at": time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(started)),
        "figure_dir": str(RP.FIGURES),
        "results_dir": str(RP.RESULTS),
        "formats": list(FORMATS),
        "max_width_mm": MAX_WIDTH_MM,
        "min_dpi": MIN_DPI,
        "non_primary_by_design": sorted(NON_PRIMARY_BY_DESIGN),
        "figures": {},
    }
    for stem, problems in results.items():
        entry: dict = {"script": PRIMARY[stem], "ok": not problems, "problems": problems,
                       "outputs": {}}
        for suffix in FORMATS:
            path = RP.FIGURES / f"{stem}.{suffix}"
            if not path.is_file():
                continue
            info: dict = {"bytes": path.stat().st_size,
                          "modified": time.strftime("%Y-%m-%dT%H:%M:%S",
                                                    time.localtime(path.stat().st_mtime))}
            if suffix == "svg":
                w, h = _svg_size_mm(path.read_text(encoding="utf-8", errors="replace"))
                info["width_mm"] = None if w is None else round(w, 2)
                info["height_mm"] = None if h is None else round(h, 2)
            elif suffix in ("png", "tiff"):
                with Image.open(path) as im:
                    info["pixels"] = list(im.size)
                    info["dpi"] = [float(x) for x in (im.info.get("dpi") or ())]
                    if suffix == "tiff":
                        info["compression"] = im.info.get("compression")
            entry["outputs"][suffix] = info
        payload["figures"][stem] = entry

    try:
        RP.RESULTS.mkdir(parents=True, exist_ok=True)
        out = RP.RESULTS / "figure_provenance.json"
        out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        return out
    except OSError as exc:
        print(f"  (could not write figure provenance: {exc})")
        return None


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--existing", action="store_true",
                    help="audit a run's EXISTING figures: check every export rule but not "
                         "freshness. For inspecting a finished or archived run. The production "
                         "figure stage never passes this -- there, a stale figure is the defect.")
    ap.add_argument("--since", type=float, default=0.0, metavar="EPOCH",
                    help="when the figure stage began. Files under the run's roots written at or "
                         "after this are the stage's own side-car outputs, not its inputs, and are "
                         "excluded from the freshness comparison. run_pipeline passes it "
                         "automatically; omit it when auditing by hand.")
    a = ap.parse_args(argv)

    started = time.time()
    newest_input = 0.0 if a.existing else _newest_input_mtime(a.since)
    print(f"figure directory: {RP.FIGURES}")
    print(f"run tag         : {RP.TAG}")
    if a.existing:
        print("mode            : --existing (export rules only; freshness not checked)")
    if newest_input:
        print(f"newest input    : {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(newest_input))}")

    results = {stem: check_stem(stem, newest_input=newest_input, started=started,
                                check_freshness=not a.existing)
               for stem in PRIMARY}
    for stem, problems in results.items():
        print(f"  {'FAIL' if problems else 'ok  '}  {stem}")

    all_problems = [f"{stem}: {p}" for stem, problems in results.items() for p in problems]
    all_problems.extend(check_configured_set())

    provenance = write_provenance(results, started=started)
    if provenance:
        print(f"  wrote {provenance}")

    if all_problems:
        print("\nFigure QA failed:")
        for problem in all_problems:
            print(f"  - {problem}")
        return 1
    print(f"\nFigure QA passed: {len(PRIMARY)} primary figures, "
          f"{len(PRIMARY) * len(FORMATS)} files, all <= {MAX_WIDTH_MM:.0f} mm wide")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
