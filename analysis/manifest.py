"""manifest.py -- the single definition of what a candidate manifest must satisfy.

Every stage that reads `run_manifest_<tag>.csv` validates it THROUGH THIS MODULE. Before it existed,
`run_pipeline.py` checked columns, uniqueness and FASTA shape, while `run_phase1.py` and
`aggregate_phase1.py` opened the same file with a bare `csv.DictReader` and no checks at all. That
made the driver's validation an advisory rather than a contract: invoking `run_phase1.py` directly
-- which is exactly what an operator does when resuming a long batch by hand -- skipped all of it,
and a bad row surfaced hours later as one failed candidate among many.

The checks are the ones that are cheap now and expensive later:

column set          a missing column is a typo in the generator, not a data condition
run_name            non-empty, unique, and usable as a directory name on every platform.
                    Uniqueness is case-INSENSITIVE: two names differing only in case are ONE
                    directory on Windows and macOS, so such rows are one output target, not two,
                    and the second candidate would overwrite the first's bundle.
fasta               repository-relative, resolvable and readable. An absolute path pins the
                    manifest to one machine; relative paths are what make it distributable.
protein             validated by the PRODUCTION validator,
                    `predictor.input_validation.read_protein` -- one record, non-empty, protein
                    alphabet, not nucleotide. Using the same function means this gate accepts
                    exactly what the pipeline accepts, rather than holding a second, drifting
                    opinion about what a valid query is.
organism_acc        present, since the batch passes it to `--organism`; an empty one silently
                    changes which genome a candidate is scanned against.

Validation returns findings rather than raising, so a caller can print all of them at once. A
manifest with any finding is not executable; `load_or_exit()` is the entry point every stage uses.
"""
from __future__ import annotations

import csv
import re
from collections import Counter
from pathlib import Path

#: Columns every stage depends on. Extra columns are carried through untouched.
REQUIRED_COLUMNS: tuple[str, ...] = ("run_name", "fasta", "organism_acc", "family", "family_flag")

#: A run_name becomes a directory under `analysis/run_<tag>/jobs/`, so it must survive every
#: filesystem we ship to. Windows is the binding constraint: no <>:"/\|?* control characters, no
#: trailing dot or space, and no reserved device name.
_UNSAFE_CHARS = re.compile('[<>:"/\\\\|?*\\x00-\\x1f]')
_RESERVED = frozenset(
    ["CON", "PRN", "AUX", "NUL"]
    + [f"COM{i}" for i in range(1, 10)]
    + [f"LPT{i}" for i in range(1, 10)]
)


def _read_protein():
    """The production protein validator, resolved the same way every other stage resolves it.

    Imported lazily and through `project_config` so this module stays importable on its own, works
    whether the predictor is pip-installed or a sibling checkout, and honours `TFOP_REPO`.
    """
    import sys
    try:
        from predictor.input_validation import read_protein
        return read_protein
    except ModuleNotFoundError:
        pass
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from project_config import PREDICTOR_ROOT
    if str(PREDICTOR_ROOT) not in sys.path:
        sys.path.insert(0, str(PREDICTOR_ROOT))
    from predictor.input_validation import read_protein
    return read_protein


def unsafe_run_name(name: str) -> str | None:
    """Why `name` cannot be a directory, or None if it can."""
    if not name:
        return "empty run_name"
    if name != name.strip():
        return "leading or trailing whitespace"
    if _UNSAFE_CHARS.search(name):
        return "contains a character illegal in a Windows path"
    if name.endswith(".") or name.endswith(" "):
        return "ends with a dot or space (illegal on Windows)"
    if name in {".", ".."} or name.split(".")[0].upper() in _RESERVED:
        return "reserved filesystem name"
    if len(name) > 120:
        return f"too long for a directory name ({len(name)} chars)"
    return None


def resolve_fasta(raw: str, project_root: Path) -> Path:
    """The FASTA a row names, resolved against the repository root when it is relative."""
    path = Path(raw)
    return path if path.is_absolute() else project_root / path


def validate_rows(rows: list[dict], project_root: Path, *,
                  check_sequences: bool = True) -> list[str]:
    """Every structural problem with `rows`. Empty means the manifest is executable."""
    issues: list[str] = []
    if not rows:
        return ["manifest has no candidates"]

    names = [(r.get("run_name") or "").strip() for r in rows]
    folded = Counter(n.casefold() for n in names if n)
    for name, count in sorted(folded.items()):
        if count > 1:
            clash = sorted({n for n in names if n.casefold() == name})
            issues.append(f"duplicate output target ({count} rows): {', '.join(clash)}")

    read_protein = None
    if check_sequences:
        try:
            read_protein = _read_protein()
        except (ModuleNotFoundError, RuntimeError) as exc:
            # Not being able to reach the validator is itself a finding: it means this manifest was
            # NOT sequence-checked. Reporting that is honest; silently downgrading to a structural
            # check would let a nucleotide FASTA through as "validated".
            return issues + [f"cannot reach the production protein validator, so no FASTA was "
                             f"sequence-checked: {exc}"]
    for row_no, row in enumerate(rows, start=2):        # row 1 is the header
        name = (row.get("run_name") or "").strip()
        why = unsafe_run_name(name)
        if why:
            issues.append(f"row {row_no}: unusable run_name {name!r}: {why}")

        if not (row.get("organism_acc") or "").strip():
            issues.append(f"row {row_no}: empty organism_acc ({name or 'unnamed'})")

        raw = (row.get("fasta") or "").strip()
        if not raw:
            issues.append(f"row {row_no}: empty fasta path")
            continue
        if Path(raw).is_absolute():
            issues.append(f"row {row_no}: fasta must be repository-relative, got {raw}")
        path = resolve_fasta(raw, project_root)
        if not path.is_file():
            issues.append(f"row {row_no}: FASTA missing: {raw}")
            continue
        if check_sequences:
            try:
                read_protein(str(path))
            except ValueError as exc:
                issues.append(f"row {row_no}: {raw}: {exc}")
            except OSError as exc:
                issues.append(f"row {row_no}: FASTA unreadable: {exc}")
    return issues


def read(path: Path) -> tuple[list[dict], list[str]]:
    """(rows, issues). Issues at this level are about the FILE, not about its contents."""
    if not path.is_file():
        return [], [f"manifest does not exist: {path}"]
    try:
        with path.open(encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh)
            missing = sorted(set(REQUIRED_COLUMNS) - set(reader.fieldnames or []))
            if missing:
                return [], [f"missing column(s): {', '.join(missing)}"]
            return list(reader), []
    except (OSError, csv.Error) as exc:
        return [], [f"cannot read manifest: {exc}"]


def validate(path: Path, project_root: Path, *, check_sequences: bool = True) -> list[str]:
    """Every problem with the manifest at `path`. Empty means it is safe to execute."""
    rows, issues = read(path)
    if issues:
        return issues
    return validate_rows(rows, project_root, check_sequences=check_sequences)


def report(path: Path, project_root: Path, *, check_sequences: bool = True,
           limit: int = 10) -> bool:
    """Print a one-line verdict plus the first `limit` findings. True when executable."""
    issues = validate(path, project_root, check_sequences=check_sequences)
    if not issues:
        rows, _ = read(path)
        print(f"manifest validation   : OK ({len(rows)} unique, readable, single-record FASTAs)")
        return True
    print(f"manifest validation   : FAILED ({len(issues)} issue(s)) in {path.name}")
    for issue in issues[:limit]:
        print(f"    {issue}")
    if len(issues) > limit:
        print(f"    ... and {len(issues) - limit} more")
    return False


def load_or_exit(path: Path, project_root: Path, *, check_sequences: bool = True) -> list[dict]:
    """The rows, or `SystemExit(2)` with every finding printed.

    This is the entry point for stages, and it is what makes running a stage script DIRECTLY as
    safe as running it through `run_pipeline.py`.
    """
    rows, issues = read(path)
    if not issues:
        issues = validate_rows(rows, project_root, check_sequences=check_sequences)
    if issues:
        print(f"refusing to run: {len(issues)} problem(s) with {path}")
        for issue in issues:
            print(f"  {issue}")
        raise SystemExit(2)
    return rows


def _cli(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="Validate a candidate manifest.")
    ap.add_argument("manifest", help="path to run_manifest_<tag>.csv")
    ap.add_argument("--no-sequences", action="store_true",
                    help="check structure only; do not open or validate the FASTAs")
    a = ap.parse_args(argv)
    root = Path(__file__).resolve().parent.parent
    path = Path(a.manifest)
    if not path.is_absolute() and not path.exists():
        path = root / "analysis" / path
    return 0 if report(path.resolve(), root, check_sequences=not a.no_sequences) else 2


if __name__ == "__main__":
    raise SystemExit(_cli())
