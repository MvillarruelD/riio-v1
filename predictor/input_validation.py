"""Validation for the public single-protein prediction contract.

The production pipeline accepts exactly one bacterial TF protein.  Keeping this
at the boundary prevents a multi-record FASTA, nucleotide sequence, or typo in
the family name from flowing through protein and structure backends as though it
were a valid query.
"""
from __future__ import annotations

import json
import re
from pathlib import Path


_PROTEIN_ALPHABET = frozenset("ACDEFGHIKLMNPQRSTVWYBXZJUO")
_DNA_ALPHABET = frozenset("ACGTUN")
_VALIDATED_FAMILIES = frozenset({"MerR", "ArsR/SmtB"})
_PARTIAL_FAMILIES = frozenset({"Fur"})


def _source_text(source: str) -> str:
    value = str(source)
    # Raw protein sequences routinely exceed Windows' path-length limit.  Only
    # ask the filesystem about strings that can plausibly be a path.
    if "\n" not in value and "\r" not in value and len(value) < 240:
        path = Path(value)
        if path.exists():
            if not path.is_file():
                raise ValueError(f"sequence input is not a file: {path}")
            return path.read_text(encoding="utf-8", errors="replace")
    return value


def _single_fasta_sequence(text: str) -> str:
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        raise ValueError("TF protein sequence is empty")
    if not any(line.startswith(">") for line in lines):
        return "".join(lines)

    records: list[str] = []
    current: list[str] | None = None
    for line in lines:
        if line.startswith(">"):
            if current is not None:
                records.append("".join(current))
            current = []
        elif current is None:
            raise ValueError("FASTA sequence data appears before the first header")
        else:
            current.append(line)
    if current is not None:
        records.append("".join(current))
    if len(records) != 1:
        raise ValueError(f"production prediction requires exactly one FASTA record; found {len(records)}")
    return records[0]


def read_protein(source: str) -> str:
    """Return one validated, upper-case protein sequence from raw text or FASTA path."""
    seq = re.sub(r"\s+", "", _single_fasta_sequence(_source_text(source))).upper()
    if seq.endswith("*"):
        seq = seq[:-1]
    if not seq:
        raise ValueError("TF protein sequence is empty")
    invalid = sorted(set(seq) - _PROTEIN_ALPHABET)
    if invalid:
        raise ValueError(f"TF protein contains invalid residue symbol(s): {', '.join(invalid)}")
    if len(seq) >= 20 and sum(base in _DNA_ALPHABET for base in seq) / len(seq) > 0.90:
        raise ValueError("input appears to be nucleotide sequence; provide the translated TF protein")
    if len(seq) < 20:
        raise ValueError(f"TF protein is implausibly short ({len(seq)} aa; minimum 20 aa)")
    return seq


def known_families() -> frozenset[str]:
    """Canonical family labels accepted by a manual production override."""
    families = set(_VALIDATED_FAMILIES | _PARTIAL_FAMILIES)
    try:
        from predictor.annotate import family_db
        path = Path(family_db._PFAM_MAP)
        if path.exists():
            families.update(json.loads(path.read_text(encoding="utf-8")).values())
    except Exception:
        pass
    try:
        from predictor.annotate import ssn_clusters
        families.update(ssn_clusters.FAMILIES)
    except Exception:
        pass
    return frozenset(families)


def validate_family(family: str) -> str:
    family = str(family).strip()
    if family not in known_families():
        raise ValueError(
            f"unknown TF family {family!r}; use a canonical family label or omit --family for classification"
        )
    return family


def family_support(family: str) -> str:
    """Reporting-only capability label; it never changes the production route."""
    if family in _VALIDATED_FAMILIES:
        return "family_validated"
    if family in _PARTIAL_FAMILIES:
        return "partially_validated"
    try:
        from predictor.annotate import ssn_clusters
        if family in ssn_clusters.FAMILIES:
            return "ssn_supported_unvalidated"
    except Exception:
        pass
    return "generic_unvalidated"
