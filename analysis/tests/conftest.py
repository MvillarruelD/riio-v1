"""Shared fixtures for the analysis-workspace tests.

`analysis/` is a directory of scripts, not an installed package, so the modules under test are
imported by putting that directory on `sys.path` exactly as the production drivers do.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

WORKSPACE = Path(__file__).resolve().parents[2]
ANALYSIS = WORKSPACE / "analysis"

if str(ANALYSIS) not in sys.path:
    sys.path.insert(0, str(ANALYSIS))


@pytest.fixture(scope="session")
def workspace() -> Path:
    return WORKSPACE


@pytest.fixture(scope="session")
def analysis_dir() -> Path:
    return ANALYSIS


@pytest.fixture
def protein_fasta(tmp_path):
    """Write a valid single-record protein FASTA and return its path."""
    def _make(name: str = "tf", seq: str = "MSTNPKPQRKTKRNTNRRPQDVKFPGGMPEPTIDESEQ") -> Path:
        p = tmp_path / f"{name}.fasta"
        p.write_text(f">{name}\n{seq}\n", encoding="utf-8")
        return p
    return _make
