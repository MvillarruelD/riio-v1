#!/usr/bin/env python3
"""Regenerate the deterministic manifest for packaged TFOP reference data."""
from __future__ import annotations

import hashlib
import json
import argparse
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "predictor" / "data"
OUT = DATA / "MANIFEST.json"
DATABASE_RELEASE = "0.3.0"
RELEASE_DATE = "2026-08-31"


def _records(path: Path) -> int | None:
    if path.suffix in {".faa", ".fasta"}:
        with path.open(encoding="utf-8", errors="replace") as handle:
            return sum(line.startswith(">") for line in handle)
    if path.suffix == ".json" and path.parent.name != "schemas" and path.name != OUT.name:
        try:
            value = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return None
        return len(value) if isinstance(value, (list, dict)) else None
    return None


def build() -> dict:
    files = []
    for path in sorted(DATA.rglob("*")):
        if not path.is_file() or path == OUT or "__pycache__" in path.parts:
            continue
        raw = path.read_bytes()
        row = {
            "path": path.relative_to(DATA).as_posix(),
            "bytes": len(raw),
            "sha256": hashlib.sha256(raw).hexdigest(),
        }
        records = _records(path)
        if records is not None:
            row["records"] = records
        files.append(row)
    return {
        "schema_version": 1,
        "database_release": DATABASE_RELEASE,
        "release_date": RELEASE_DATE,
        "file_count": len(files),
        "total_bytes": sum(row["bytes"] for row in files),
        "files": files,
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="fail if the committed manifest is stale")
    args = parser.parse_args(argv)
    rendered = json.dumps(build(), indent=2, sort_keys=True) + "\n"
    if args.check:
        current = OUT.read_text(encoding="utf-8") if OUT.is_file() else ""
        if current != rendered:
            print(f"stale data manifest: regenerate with {Path(__file__).name}")
            return 1
        print(f"data manifest current: {OUT}")
        return 0
    OUT.write_text(rendered, encoding="utf-8", newline="\n")
    print(OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
