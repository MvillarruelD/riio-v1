#!/usr/bin/env python
"""Publish a share directory from staging; never merge old and new deliverables.

Scientific interpretation belongs in the reports, not this packaging script.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import shutil
import sys
import tempfile
import uuid
from pathlib import Path

import manifest as MF
import run_paths as RP
from project_config import PREDICTOR_ROOT

if str(PREDICTOR_ROOT) not in sys.path:
    sys.path.insert(0, str(PREDICTOR_ROOT))
from predictor.report.bundle_contract import contained_file, inventory, validate_files

HERE = Path(__file__).resolve().parent
PROJ = HERE.parent
JOBS = RP.JOBS
OUT = PROJ / f"SHARE_FOR_AUTHORS_{RP.TAG}"


def _existing(*paths: Path) -> Path:
    for path in paths:
        if path.is_file():
            return path
    raise FileNotFoundError("missing run export: " + ", ".join(str(path) for path in paths))


def _manifest_path() -> Path:
    return RP.MANIFEST


def _rows(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError("candidate manifest is empty")
    names = set()
    for row in rows:
        name = row.get("run_name", "")
        issue = MF.unsafe_run_name(name)
        if issue or name.casefold() in names:
            raise ValueError(f"invalid or duplicate run_name {name!r}: {issue or 'duplicate'}")
        if not row.get("fasta"):
            raise ValueError(f"missing candidate input path: {name}")
        names.add(name.casefold())
    return rows


def _coverage_note(n: int) -> str:
    canonical = RP.ANALYSIS / f"run_manifest_{RP.TAG}.csv"
    if not canonical.is_file() or canonical.resolve() == _manifest_path().resolve():
        return ""
    canonical_names = {row["run_name"].casefold() for row in _rows(canonical)}
    selected = {row["run_name"].casefold() for row in _rows(_manifest_path())}
    missing, extra = canonical_names - selected, selected - canonical_names
    if not missing and not extra:
        return ""
    return (f"This package differs from the canonical manifest: {n} selected; "
            f"{len(missing)} canonical candidates absent; {len(extra)} additional candidates. "
            "It must not be described as the complete canonical run.\n")


def build_readme(*, no_bundles: bool = False) -> str:
    rows = _rows(_manifest_path())
    predictions = list(csv.DictReader(
        _existing(RP.COMPLEMENT_TSV, RP.RESULTS / "per_regulator_predictions.tsv")
        .open(encoding="utf-8"), delimiter="\t"))
    n_metal = sum(row.get("inducer_class") == "metal" for row in predictions)
    n_protein = sum(
        row.get("inducer_class") == "metal"
        and bool((row.get("ssn_ligand") or "").strip()
                 or (row.get("gate") or "").strip() == "True")
        for row in predictions
    )
    scope = ("Individual bundles were omitted by --no-bundles. This is a summary export."
             if no_bundles else "05_per_regulator/ contains the selected individual report directories.")
    return (
        f"# Metalloregulator prediction package\n\nPrepared {dt.date.today().isoformat()} for run {RP.TAG}. "
        "Private working material for the manuscript authors; do not redistribute without agreement.\n\n"
        f"The selected manifest contains {len(rows)} candidates. The published survey contains 150; "
        "the cohorts are close but not identical. "
        f"The comparable protein-level subset calls {n_protein}/{len(predictions)} "
        f"({100*n_protein/len(predictions):.0f}%) metal-responsive. The broader fused result is "
        f"{n_metal}/{len(predictions)} ({100*n_metal/len(predictions):.0f}%) and reuses MetalNet and "
        "genomic-context evidence, so it is not an independent replication of the survey result. "
        f"{_coverage_note(len(rows))}\n\n"
        "## Files\n\n"
        "- 01_reports/: main Word report and concise Markdown summary.\n"
        "- 02_figures/: available SVG, PDF, PNG and TIFF exports.\n"
        "- 03_data/: result tables, candidate manifest and locally copied candidate inputs.\n"
        "- 04_manuscript_integration/: available composition exports.\n"
        f"- {scope}\n"
        "- 06_benchmark_context/: historical benchmark panels, when present, kept separately.\n\n"
        "## Reuse\n\n"
        "Keep this directory intact. Individual REPORT.html files open offline. The manifest in "
        "03_data/ uses paths relative to 03_data/; the original manifest is retained alongside it. "
        "share_manifest.json records file sizes and SHA-256 digests, excluding its own digest.\n\n"
        "Predicted operators and regulons are ranked hypotheses for experimental prioritisation, "
        "not validated sites. These are recorded deliverables, not a complete offline replay environment. Read the "
        "individual bundle guides for input and provenance limitations. Packaging does not "
        "validate predictions or establish computational reproducibility.\n\n"
        "Benchmark results are maintained separately under analysis/benchmarking/. No benchmark "
        "accuracy figures are inserted into this guide by the packaging code. Historical reports "
        "may still contain benchmark discussion.\n"
    )


def build_summary() -> str:
    """A short, cautious index for readers who do not start with the Word report."""
    rows = list(csv.DictReader(
        _existing(RP.COMPLEMENT_TSV, RP.RESULTS / "per_regulator_predictions.tsv")
        .open(encoding="utf-8"), delimiter="\t"))
    n = len(rows)
    metal = sum(row.get("inducer_class") == "metal" for row in rows)
    protein = sum(
        row.get("inducer_class") == "metal"
        and bool((row.get("ssn_ligand") or "").strip()
                 or (row.get("gate") or "").strip() == "True")
        for row in rows
    )
    return (
        "# Report summary\n\n"
        f"Run **{RP.TAG}** contains **{n}** candidate regulators from four bacterial genomes and "
        f"{n}/{n} recorded per-regulator bundles. The published SSN survey contains 150 candidates; "
        "the two cohorts are not identical.\n\n"
        f"The pipeline calls **{protein}/{n} ({100*protein/n:.0f}%)** candidates metal-responsive "
        "when restricted to calls supported by the SSN clade or coordination gate. This is the "
        "closest comparison with the survey's **41/150 (27%)** MetalNet-site result, but the differing "
        "cohorts make it convergence rather than validation. The full fusion calls "
        f"**{metal}/{n} ({100*metal/n:.0f}%)** metal-responsive; it also uses MetalNet and genomic "
        "context and is not an independent comparison.\n\n"
        "Each regulator bundle records the fused inducer hypothesis, ranked operator candidates, a "
        "candidate regulon and source-level evidence. Operator and regulon results are prioritisation "
        "tools, not experimentally validated calls. Read the Word report for methods, evidence tiers "
        "and limitations before reusing a result.\n"
    )


def copy(src: Path, dst: Path, label: str) -> int:
    """Copy one required deliverable; a missing input fails publication."""
    if not src.is_file():
        raise FileNotFoundError(f"missing {label}: {src}")
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return 1


def assemble(destination: Path, *, no_bundles: bool = False) -> None:
    """Populate a new directory from the selected run."""
    rows = _rows(_manifest_path())
    copy(RP.REPORT_DOCX, destination / "01_reports" / RP.REPORT_DOCX.name, "report")
    (destination / "01_reports" / "report_summary.md").write_text(build_summary(), encoding="utf-8")
    for source, name in ((_existing(RP.COMPLEMENT_TSV, RP.RESULTS / "per_regulator_predictions.tsv"),
                          "per_regulator_predictions.tsv"),
                         (_existing(RP.REGULON_SUMMARY, RP.RESULTS / "regulon_summary.tsv"),
                          "regulon_summary.tsv")):
        copy(source, destination / "03_data" / name, name)
    ledger = RP.RESULTS / "evidence_ledger.tsv"
    if ledger.is_file():
        copy(ledger, destination / "03_data" / ledger.name, "evidence ledger")
    copy(_manifest_path(), destination / "03_data" / "candidate_manifest.original.csv", "manifest")
    copied_rows = []
    for row in rows:
        source = MF.resolve_fasta(row["fasta"], PROJ)
        relative = f"inputs/{row['run_name']}.fasta"
        copy(source, destination / "03_data" / relative, "candidate input")
        copied_rows.append({**row, "fasta": relative})
        if not no_bundles:
            bundle = contained_file(JOBS, row["run_name"])
            issues = validate_files(bundle)
            if issues:
                raise ValueError(f"cannot share {row['run_name']}: {'; '.join(issues)}")
            target = destination / "05_per_regulator" / row["run_name"]
            shutil.copytree(bundle, target)
            issues = validate_files(target)
            if issues:
                raise ValueError(f"copied bundle failed verification: {'; '.join(issues)}")
    with (destination / "03_data" / "candidate_manifest.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(copied_rows)
    figures = sorted(p for p in RP.FIGURES.iterdir() if p.suffix.lower() in {".png", ".svg", ".pdf", ".tiff"})
    if not figures:
        raise ValueError("no figure exports are available")
    for figure in figures:
        folder = "06_benchmark_context" if figure.name.startswith("F21_") else "02_figures"
        copy(figure, destination / folder / figure.name, "figure")
        if figure.name.startswith(("F23_", "F24_")):
            copy(figure, destination / "04_manuscript_integration" / figure.name, "composition")
    (destination / "README.md").write_text(build_readme(no_bundles=no_bundles), encoding="utf-8")
    (destination / "share_manifest.json").write_text(json.dumps({
        "schema_version": 1, "run_tag": RP.TAG, "candidate_count": len(rows),
        "individual_bundles_included": not no_bundles, "artifacts": inventory(destination),
    }, indent=2), encoding="utf-8")


def publish(*, no_bundles: bool = False) -> Path:
    """Publish after assembly succeeds; preserve any prior release as a sibling backup."""
    final = OUT.resolve()
    if final.parent != PROJ.resolve() or not final.name.startswith("SHARE_FOR_AUTHORS_"):
        raise ValueError("share output must be a direct child of the project root")
    staging = Path(tempfile.mkdtemp(prefix=".share-", dir=final.parent))
    backup = None
    try:
        assemble(staging, no_bundles=no_bundles)
        if final.exists():
            backup = final.with_name(f".{final.name}.previous-{uuid.uuid4().hex[:10]}")
            final.replace(backup)
        try:
            staging.replace(final)
        except BaseException:
            if backup is not None and backup.exists() and not final.exists():
                backup.replace(final)
            raise
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    print(f"Shared deliverables: {final}")
    if backup is not None:
        print(f"Previous release preserved: {backup}")
    return final


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-bundles", action="store_true", help="create a summary without individual bundles")
    args = parser.parse_args(argv)
    try:
        publish(no_bundles=args.no_bundles)
    except (OSError, ValueError) as exc:
        print(f"Share publication failed: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
