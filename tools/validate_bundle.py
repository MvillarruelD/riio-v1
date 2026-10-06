#!/usr/bin/env python3
"""Validate the internal consistency of one completed prediction bundle.

Usage:
    python tools/validate_bundle.py results/jobs/<name>
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

from predictor.structure.af3_export import validate_jobs
from predictor.report.bundle_contract import validate_files


def _load_json(path: Path):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read valid JSON from {path}: {exc}") from exc


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            return list(reader.fieldnames or []), list(reader)
    except OSError as exc:
        raise ValueError(f"cannot read {path}: {exc}") from exc


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def validate_bundle(bundle: Path) -> dict[str, int | str]:
    """Validate a prediction bundle and return a compact summary; raise ValueError on inconsistency."""
    bundle = bundle.resolve()
    _require(bundle.is_dir(), f"bundle directory does not exist: {bundle}")
    _require(not (issues := validate_files(bundle)), "; ".join(issues))

    manifest = _load_json(bundle / "bundle_manifest.json")
    _require(manifest.get("status") == "complete", "bundle manifest is not marked complete")
    _require(manifest.get("report_generated") is True, "bundle manifest says the report was not generated")
    _require(not manifest.get("warnings"), f"bundle manifest has warnings: {manifest.get('warnings')}")
    for rel in manifest.get("required_artifacts") or []:
        _require((bundle / rel).is_file(), f"required artifact is missing: {rel}")

    dossier = _load_json(bundle / "dossier.json")
    operators = _load_json(bundle / "operators_ranked.json")
    per_hit = _load_json(bundle / "regulation" / "per_hit_regulation.json").get("per_hit") or []
    af3_jobs = _load_json(bundle / "af3_jobs.json")
    validate_jobs(af3_jobs)

    binding_fields, binding_rows = _read_tsv(bundle / "binding_sites.tsv")
    _, regulon_rows = _read_tsv(bundle / "regulon.tsv")
    hits = (dossier.get("rescan") or {}).get("hits") or []
    operons = [row for row in (dossier.get("regulon") or [])
               if isinstance(row, dict) and "first_gene" in row]
    _require(len(binding_rows) == len(hits),
             f"binding-site row mismatch: TSV={len(binding_rows)}, dossier={len(hits)}")
    _require(len(regulon_rows) == len(operons),
             f"regulon row mismatch: TSV={len(regulon_rows)}, dossier={len(operons)}")
    for field in ("pvalue", "qvalue", "qvalue_kind", "generator"):
        _require(field in binding_fields, f"binding_sites.tsv is missing {field!r}")
    if hits:
        ranked_hits = sorted(hits, key=lambda hit: -(hit.get("score") or 0))
        for field in ("pvalue", "qvalue"):
            value = ranked_hits[0].get(field)
            if value is not None:
                # The TSV deliberately uses four significant digits (``.3e``); compare at that
                # display precision rather than requiring the rounded text to equal the JSON float.
                _require(math.isclose(float(binding_rows[0][field]), float(value), rel_tol=5e-4),
                         f"top binding-site {field} differs between TSV and dossier")

    ranked = operators.get("ranked") or []
    primary = operators.get("primary")
    _require(ranked and primary, "ranked operators or primary operator is empty")
    _require(primary == ranked[0], "primary operator is not the first full ranked record")
    for index, row in enumerate(ranked, 1):
        _require("provenance" in row and "status" in row,
                 f"ranked operator {index} has lost provenance or status")

    if per_hit:
        missing_modes = sum(not row.get("mode") for row in per_hit)
        missing_tss = sum(row.get("tss") is None for row in per_hit)
        _require(not missing_modes, f"{missing_modes}/{len(per_hit)} per-hit mode calls are missing")
        _require(not missing_tss, f"{missing_tss}/{len(per_hit)} per-hit TSS calls are missing")
    headline = dossier.get("promoter_occlusion") or {}
    _require(headline.get("mode"), "headline promoter-occlusion mode is missing")
    _require(headline.get("tss") is not None or (headline.get("promoter") or {}).get("tss") is not None,
             "headline promoter-occlusion TSS is missing")

    report_path = bundle / "REPORT.html"
    report = report_path.read_text(encoding="utf-8")
    local_refs = set(re.findall(r'''(?:src|href)=["'](?!https?://|data:|#)([^"']+)["']''', report))
    missing_refs = sorted(ref for ref in local_refs if not (bundle / ref.split("?", 1)[0]).is_file())
    _require(not missing_refs, f"REPORT.html has missing local assets: {', '.join(missing_refs)}")

    return {
        "status": "valid",
        "binding_sites": len(binding_rows),
        "regulon_operons": len(regulon_rows),
        "ranked_operators": len(ranked),
        "promoter_calls": len(per_hit),
        "af3_jobs": len(af3_jobs),
        "report_assets": len(local_refs),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path, help="completed results/jobs/<name> directory")
    parser.add_argument("--files-only", action="store_true", help="check file integrity without interpreting scientific records")
    args = parser.parse_args()
    try:
        if args.files_only:
            issues = validate_files(args.bundle)
            _require(not issues, "; ".join(issues))
            summary = {"status": "valid", "scope": "files only"}
        else:
            summary = validate_bundle(args.bundle)
    except ValueError as exc:
        parser.exit(1, f"ERROR: {exc}\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
