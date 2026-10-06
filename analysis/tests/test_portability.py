"""Portability contracts for the analysis layer.

A default that points at one developer's home directory does not fail on that machine, which is
exactly why it survives review and then misconfigures every other install. These tests make the
absence of such defaults a checked property rather than a habit.
"""
from __future__ import annotations

import ast
from pathlib import Path

import pytest

from conftest import ANALYSIS, WORKSPACE

#: Spellings of a developer-specific location. Checked case-insensitively for the Windows forms,
#: since `c:/users/...` reaches the same directory.
FORBIDDEN = ("c:/users/matia", "c:\\users\\matia", "/home/matia", "/users/matia")

#: Historical material, kept for reference and never executed. It is allowed to record the paths
#: that were true when it was written; that is what makes it a record.
EXCLUDED_DIRS = {"__pycache__", "project_archive", "snapshots", "out"}


#: This file holds every forbidden spelling as DATA -- it is the detector. Scanning it reports the
#: detector as the defect, which is the same trap `run_pipeline.stale_input_scripts()` documents.
SELF = Path(__file__).resolve()


def _sources() -> list[Path]:
    return [p for p in sorted(ANALYSIS.rglob("*.py"))
            if not EXCLUDED_DIRS.intersection(p.parts) and p.resolve() != SELF]


def test_there_are_sources_to_check():
    """Guard the guard: a glob that silently matches nothing passes every test below."""
    assert len(_sources()) >= 25


@pytest.mark.parametrize("path", _sources(), ids=lambda p: str(p.relative_to(ANALYSIS)))
def test_no_developer_specific_paths(path):
    offenders = []
    for lineno, line in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
        low = line.lower()
        for needle in FORBIDDEN:
            if needle in low:
                offenders.append(f"{lineno}: {line.strip()[:100]}")
                break
    assert not offenders, (
        f"{path.relative_to(WORKSPACE)} hardcodes a developer-specific path:\n  "
        + "\n  ".join(offenders)
        + "\n\nResolve it through analysis/project_config.py, family_kb/kb_paths.py, "
          "analysis/benchmarking/bench_paths.py, or an environment override."
    )


@pytest.mark.parametrize("path", _sources(), ids=lambda p: str(p.relative_to(ANALYSIS)))
def test_sources_parse(path):
    """Every script parses. Cheap, and it catches a bad mechanical edit across many files."""
    ast.parse(path.read_text(encoding="utf-8", errors="replace"), filename=str(path))


def test_benchmark_suites_derive_roots_by_name_not_by_counting():
    """A suite two levels below `analysis/` must not compute the workspace with `parents[2]`.

    This is the exact defect the `analysis/<suite>_bench/` -> `analysis/benchmarking/<suite>/` move
    introduced: every root computed by counting `..` silently resolved one directory too high, so
    scripts read `analysis/` as though it were the workspace and found nothing.
    """
    offenders = []
    for path in sorted((ANALYSIS / "benchmarking").rglob("*.py")):
        if "__pycache__" in path.parts or path.name == "bench_paths.py":
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        for lineno, line in enumerate(text.splitlines(), 1):
            stripped = line.strip()
            if stripped.startswith("#"):
                continue
            # A suite script is at analysis/benchmarking/<suite>/x.py, so the workspace is
            # parents[3] and `analysis/` is parents[2]. parents[2] AS the workspace is the bug.
            if "WORKSPACE" in stripped and "parents[2]" in stripped:
                offenders.append(f"{path.relative_to(ANALYSIS)}:{lineno}: {stripped[:90]}")
            if "PROJ" in stripped and ".parent.parent.parent" in stripped:
                offenders.append(f"{path.relative_to(ANALYSIS)}:{lineno}: {stripped[:90]}")
    assert not offenders, ("benchmark roots computed by counting:\n  " + "\n  ".join(offenders)
                           + "\n\nImport bench_paths instead.")


def test_no_stale_benchmark_directory_names():
    """The pre-move directory names must not survive anywhere in code or configuration."""
    stale = ("analysis/discovery_bench", "analysis/regulondb_bench", "analysis/salmonella_bench",
             'HERE / "regulondb_bench"', '"salmonella_bench"', '"discovery_bench"')
    offenders = []
    targets = _sources() + [WORKSPACE / ".gitignore"]   # _sources() already excludes SELF
    for path in targets:
        if not path.is_file():
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        for needle in stale:
            if needle in text:
                offenders.append(f"{path.relative_to(WORKSPACE)}: {needle}")
    assert not offenders, "stale benchmark paths:\n  " + "\n  ".join(offenders)
