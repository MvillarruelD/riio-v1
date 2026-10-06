"""bench_paths.py -- the one place a benchmark script learns where anything lives.

Every suite sits two levels below the analysis workspace::

    analysis/benchmarking/<suite>/<script>.py

so a script that derives its roots by counting `..` gets a different answer than it did when the
suites lived at `analysis/<suite>_bench/`. Counting is what broke when they moved, and counting is
what this module exists to stop: import the name you mean instead.

    SUITE = Path(__file__).resolve().parent      # analysis/benchmarking/<suite>
    BENCHMARKING = SUITE.parent                  # analysis/benchmarking
    ANALYSIS = BENCHMARKING.parent               # analysis
    WORKSPACE = ANALYSIS.parent                  # repository root

Two roots are NOT derivable from this file's location and are resolved instead:

`REPO`
    the tf-operator-predictor checkout or install, via `analysis/project_config.py`, so it honours
    `TFOP_REPO` and falls back to conventional sibling/home locations rather than one account name.

`PROMOTERS`
    an EXTERNAL curated operator/RegulonDB tree that is not distributed with this repository and is
    never read by production. It is the same root `analysis/family_kb/kb_paths.py` pins, resolved
    through the same `TFOP_PROMOTERS_DB` variable and the same workspace-relative default, so the
    two layers cannot drift onto different copies. Scripts that need it call `require_promoters()`
    and get a message naming the variable when it is absent, rather than a `FileNotFoundError` on a
    path from somebody else's machine.

Nothing here is imported by production code, and nothing here may be imported BY production code --
see this directory's README for the boundary rule.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

#: `analysis/benchmarking` -- this file's own directory, and the parent of every suite.
BENCHMARKING = Path(__file__).resolve().parent
ANALYSIS = BENCHMARKING.parent
WORKSPACE = ANALYSIS.parent

CORE = BENCHMARKING / "core"
DISCOVERY = BENCHMARKING / "discovery"
#: The hand-picked 17-TF development panel and its ablations, archived 2026-09-24; superseded by
#: the census-based `ecoli_kb/` benchmark. Kept addressable for the scripts that still read it.
REGULONDB = WORKSPACE / "project_archive" / "benchmarking_regulondb"
ECOLI_KB = BENCHMARKING / "ecoli_kb"
SALMONELLA = BENCHMARKING / "salmonella"
ENVIRONMENT = BENCHMARKING / "environment"

# `project_config` lives in `analysis/`, which is not on the path when a suite script is run
# directly (sys.path[0] is the suite directory). Add it once, here, so every suite inherits the
# same portable predictor-root discovery instead of hardcoding a checkout location.
if str(ANALYSIS) not in sys.path:
    sys.path.insert(0, str(ANALYSIS))

from project_config import PREDICTOR_ROOT  # noqa: E402

#: The predictor checkout/install. Honours `TFOP_REPO`; see `analysis/project_config.py`.
REPO = PREDICTOR_ROOT

#: External curated-operator tree. Same variable and same default as `family_kb/kb_paths.py`.
PROMOTERS_ENV = "TFOP_PROMOTERS_DB"
_promoters = os.environ.get(PROMOTERS_ENV, "").strip()
PROMOTERS = Path(_promoters).expanduser() if _promoters else WORKSPACE.parent / "5.8 Promoters"

#: Conventional layout inside that mirror, so a suite names a table rather than a path.
REGULONDB_SOURCES = "operators_db/data_sources/regulondb"
NETWORK_REGULATOR_GENE = f"{REGULONDB_SOURCES}/NetworkRegulatorGene.tsv"
ALL_OPERATORS = "operators_db/figures/fig3_mode/source_data/all_operators.tsv"


def require_promoters(what: str = "") -> Path:
    """The external RegulonDB mirror, or a message naming what to set and why it is not shipped.

    `what` names the table the caller wanted, so the error says which measurement is unavailable
    rather than only that a directory is missing.
    """
    if PROMOTERS is None:
        raise SystemExit(
            f"{PROMOTERS_ENV} is not set, so {what or 'this benchmark'} cannot run.\n"
            "  This suite reads an EXTERNAL RegulonDB mirror that is deliberately not distributed\n"
            "  with this repository and is never read by production code. Point the variable at\n"
            f"  your copy, e.g.  {PROMOTERS_ENV}=/path/to/'5.8 Promoters'"
        )
    if not PROMOTERS.is_dir():
        raise SystemExit(f"{PROMOTERS_ENV} does not exist: {PROMOTERS}")
    return PROMOTERS


def promoters_file(relative: str, what: str = "") -> Path:
    """One table inside the external mirror, verified to exist before it is handed back."""
    path = require_promoters(what) / relative
    if not path.exists():
        raise SystemExit(f"missing from the {PROMOTERS_ENV} mirror: {path}")
    return path


def describe() -> str:
    return (f"workspace   : {WORKSPACE}\n"
            f"  analysis  : {ANALYSIS}\n"
            f"  benchmarks: {BENCHMARKING}\n"
            f"  predictor : {REPO}\n"
            f"  promoters : {PROMOTERS}{'' if PROMOTERS.is_dir() else '   (absent)'}")


if __name__ == "__main__":
    print(describe())
