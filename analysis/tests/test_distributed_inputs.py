"""Every distributed manifest must reference distributed FASTAs.

The defect this exists to prevent was real and silent: three manifests were committed while **none**
of their 490 candidate FASTAs were, so a fresh clone got manifests pointing at files that did not
exist. Nothing failed locally, because the FASTAs were sitting untracked in the working tree of the
one machine that had ever run the pipeline. `.gitignore` and `CANDIDATE_INPUTS.md` both stated the
files were versioned; only Git disagreed.

A manifest is only executable together with its inputs, so "distributed" has to mean both.
"""
from __future__ import annotations

import csv
import subprocess
from pathlib import Path

import pytest

from conftest import ANALYSIS, WORKSPACE


def _tracked(path: Path) -> bool:
    """True when Git has this file in the index (committed or staged)."""
    result = subprocess.run(["git", "ls-files", "--error-unmatch", str(path)],
                            cwd=str(WORKSPACE), capture_output=True, text=True)
    return result.returncode == 0


def _rows(manifest: Path) -> list[dict]:
    with manifest.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


@pytest.fixture(scope="module")
def git_available() -> bool:
    result = subprocess.run(["git", "rev-parse", "--is-inside-work-tree"],
                            cwd=str(WORKSPACE), capture_output=True, text=True)
    return result.returncode == 0


@pytest.fixture(scope="module")
def distributed_manifests() -> list[Path]:
    """Discover shipped manifest paths only after a released-data test is selected."""
    return sorted(
        list(ANALYSIS.glob("run_manifest*.csv"))
        + list((ANALYSIS / "benchmarking").rglob("run_manifest*.csv"))
    )


@pytest.mark.released_data
def test_there_are_manifests_to_check(distributed_manifests):
    """Guard the guard: an empty glob would make every test below vacuous."""
    assert distributed_manifests, "no run_manifest*.csv found"


@pytest.mark.released_data
def test_a_distributed_manifest_has_distributed_inputs(distributed_manifests, git_available):
    """If the manifest ships, every FASTA it names must ship too."""
    if not git_available:
        pytest.skip("not a git work tree")

    manifests = [manifest for manifest in distributed_manifests if _tracked(manifest)]
    if not manifests:
        pytest.skip("no discovered manifest is distributed")

    untracked = []
    for manifest in manifests:
        for row in _rows(manifest):
            raw = (row.get("fasta") or "").strip()
            if not raw:
                continue
            path = Path(raw)
            if not path.is_absolute():
                path = WORKSPACE / path
            if not _tracked(path):
                untracked.append(f"{manifest.name}: {raw}")
    assert not untracked, (
        f"distributed manifests reference {len(untracked)} untracked FASTA(s):\n  "
        + "\n  ".join(untracked[:10])
        + "\n\nA clone would receive a manifest pointing at files that do not exist."
    )


@pytest.mark.released_data
def test_manifest_paths_are_repository_relative(distributed_manifests):
    """An absolute path pins a manifest to one machine and is not distributable."""
    absolute = []
    for manifest in distributed_manifests:
        absolute.extend(
            f"{manifest.name}: {row['fasta']}" for row in _rows(manifest)
            if (row.get("fasta") or "").strip()
            and Path(row["fasta"].strip()).is_absolute()
        )
    assert not absolute, f"manifests carry absolute paths: {absolute[:5]}"


@pytest.mark.released_data
def test_every_referenced_fasta_exists_on_disk(distributed_manifests):
    missing = []
    for manifest in distributed_manifests:
        for row in _rows(manifest):
            raw = (row.get("fasta") or "").strip()
            if raw and not (WORKSPACE / raw).is_file():
                missing.append(f"{manifest.name}: {raw}")
    assert not missing, f"manifests reference {len(missing)} absent file(s): {missing[:5]}"


@pytest.mark.released_data
def test_no_production_manifest_points_into_the_benchmark_tree():
    """Production must never take an input from a benchmark suite.

    A row for a benchmark organism (`run_discovery.BENCHMARK_GENOMES`, e.g. Salmonella SL1344) IS a
    benchmark input, wherever its manifest lives; the rule is about survey rows.
    """
    from run_discovery import BENCHMARK_GENOMES

    offenders = []
    for manifest in ANALYSIS.glob("run_manifest*.csv"):
        for row in _rows(manifest):
            if row.get("organism_acc") in BENCHMARK_GENOMES:
                continue
            if "benchmarking/" in (row.get("fasta") or ""):
                offenders.append(f"{manifest.name}: {row['fasta']}")
    assert not offenders, "production manifests reading benchmark inputs:\n  " + "\n  ".join(offenders)
