"""run_paths.py -- where THIS run's figures, results and regulons live. One tag resolves all of it.

Every output path in the analysis layer hangs off a single run tag:

    TFOP_RUN_TAG=prod20260821 python analysis/fig_logos.py

    analysis/run_<tag>/jobs/        predictor bundles
    analysis/run_<tag>/logs/        per-candidate execution logs
    analysis/run_<tag>/figures/     figures
    analysis/run_<tag>/results/     tables and cached statistics
    analysis/run_<tag>/regulon/     regulon summary + report documents

**The tag is required. There is no default.** That is the whole point of this module, and it is
worth stating why, because a default looks harmless and is not.

Outputs from different runs share a naming scheme -- both the survey's frozen candidate list and a
genome-first run produce `{genome}__{family}__{accession}` -- and the two lists mostly overlap. So
pointing a script at the wrong run does not fail. It silently drops the candidates unique to the
right one, draws every figure, writes every table, and leaves a report that looks complete. A run
that half-succeeds is far more expensive than one that refuses to start, and the only reliable way
to refuse is to have nothing to fall back to.

Individual directories can still be overridden (`TFOP_RUN_FIGS`, `TFOP_RUN_STATS`,
`TFOP_RUN_REGULON`, `TFOP_JOBS`, `AN_MANIFEST`) for one-off inspection. Those are explicit acts; an
unset tag is not.
"""
from __future__ import annotations

import os
from pathlib import Path

from project_config import PREDICTOR_ROOT

ANALYSIS = Path(__file__).resolve().parent
REPO = PREDICTOR_ROOT


class RunTagNotSet(RuntimeError):
    """Raised at import when `TFOP_RUN_TAG` is unset. See the module docstring for why."""


def _require_tag() -> str:
    tag = os.environ.get("TFOP_RUN_TAG", "").strip()
    if not tag:
        available = sorted(p.name[4:] for p in ANALYSIS.glob("run_*") if p.is_dir())
        raise RunTagNotSet(
            "TFOP_RUN_TAG is not set, and this module has no default run.\n"
            "  Set it to the run you mean, e.g.  TFOP_RUN_TAG=prod20260821\n"
            + (f"  Runs present: {', '.join(available)}\n" if available else "")
            + "  A default would let a script read a DIFFERENT run's candidates and still produce\n"
              "  a complete-looking report, which is the failure this module exists to prevent."
        )
    from manifest import unsafe_run_name
    issue = unsafe_run_name(tag)
    if issue:
        raise ValueError(f"invalid TFOP_RUN_TAG: {issue}")
    return tag


#: Names the run. Every path below hangs off it, so one variable re-points the whole report chain.
TAG = _require_tag()
ROOT = Path(os.environ.get("TFOP_RUN_ROOT", ANALYSIS / f"run_{TAG}"))

FIGURES = Path(os.environ.get("TFOP_RUN_FIGS", ROOT / "figures"))
#: `results/`, NOT `stats/`: this is the layout the one-off figure scripts already expect from
#: `AN_OUT`, and matching it is what lets them be re-pointed by environment alone.
RESULTS = Path(os.environ.get("TFOP_RUN_STATS", ROOT / "results"))
STATS = RESULTS                      # older name, kept so nothing that imported it breaks
REGULON = Path(os.environ.get("TFOP_RUN_REGULON", ROOT / "regulon"))
VALIDATION = RESULTS / "regulon_validation.csv"

#: Predictor bundles and logs are run-scoped too. Keeping them beside the derived outputs prevents a
#: new tag from silently treating an older run's bundles as completed work. `TFOP_JOBS` remains an
#: explicit inspection/migration override; it is never the default production location.
JOBS = Path(os.environ.get("TFOP_JOBS", ROOT / "jobs"))
LOGS = Path(os.environ.get("TFOP_RUN_LOGS", ROOT / "logs"))

#: WHICH CANDIDATES THIS RUN IS ABOUT -- derived from the tag like every other path here.
#: This is the single definition; scripts import `RP.MANIFEST` rather than re-deriving it, because
#: eleven private copies of the same three lines is eleven places for the tag to be read wrongly.
MANIFEST = Path(os.environ.get("AN_MANIFEST") or ANALYSIS / f"run_manifest_{TAG}.csv")
if not MANIFEST.is_absolute():
    MANIFEST = ANALYSIS / MANIFEST

#: The per-regulator comparison against the survey. It used to live at a FIXED path shared by every
#: run, so two runs overwrote each other and -- the real danger -- if the producing script failed,
#: the report read the PREVIOUS run's numbers and still looked complete.
COMPLEMENT_TSV = RESULTS / "fig5_complement_per_regulator.tsv"

REGULON_SUMMARY = REGULON / "regulon_v2_summary.tsv"
REPORT_DOCX = REGULON / "metalloregulator_predictions_four_genomes.docx"


def ensure() -> None:
    """Create the run's output directories. Safe to call repeatedly."""
    for d in (JOBS, LOGS, FIGURES, RESULTS, REGULON):
        d.mkdir(parents=True, exist_ok=True)


def env() -> dict:
    """The environment that re-points EVERY analysis script at this run.

    Several one-off `fig_*.py` scripts read an environment variable before falling back to a path of
    their own. They do not need editing; they need the variable set. Five names carry such a
    fallback (`AN_RUN`, `AN_OUT`, `AN_FIG`, `AN_VALID`, `TFOP_JOBS`), and all five are set here.

    This is the difference between a figure stage that fails loudly and one that quietly draws an
    old statistic into a new report -- the second being the expensive kind, because the report still
    looks finished.
    """
    return {
        "TFOP_RUN_TAG": TAG,
        "TFOP_RUN_ROOT": str(ROOT),
        "TFOP_RUN_FIGS": str(FIGURES),
        "TFOP_RUN_STATS": str(RESULTS),
        "TFOP_RUN_REGULON": str(REGULON),
        "AN_RUN": str(ROOT),
        "AN_FIG": str(FIGURES),
        "AN_OUT": str(RESULTS),
        "AN_VALID": str(VALIDATION),
        "TFOP_JOBS": str(JOBS),
        # The installable predictor writes bundles below `<output>/jobs`; align that root with JOBS.
        "PREDICTOR_OUTPUT_DIR": str(JOBS.parent),
        "TFOP_RUN_LOGS": str(LOGS),
        "AN_MANIFEST": str(MANIFEST),
    }


def describe() -> str:
    return (f"run tag  : {TAG}\n"
            f"  jobs   : {JOBS}\n"
            f"  logs   : {LOGS}\n"
            f"  manifest: {MANIFEST}\n"
            f"  figures: {FIGURES}\n"
            f"  results: {RESULTS}\n"
            f"  regulon: {REGULON}")


if __name__ == "__main__":
    print(describe())
