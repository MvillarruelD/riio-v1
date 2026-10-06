"""
ligify_db.py -- look a TF up in the PUBLISHED Ligify sensor database before re-deriving its effector.

`effector/ligify.py` reimplements Ligify's genome-context route locally: find the TF's operon, read the
reactions of the neighbouring enzymes from UniProt, and score each (TF, ligand) pair with the formula in
`effector/rank.py`. That reimplementation stays exactly as it is -- it is what handles the ~90 % of
queries the published database has never seen.

What was missing was their DATA. `groov-bio/ligify-ui` ships a static database of 3,164 predicted
biosensors, each keyed by RefSeq accession and carrying the protein sequence, the candidate ligands, the
Rhea reaction and Ligify's own rank. Our discovery candidates ARE RefSeq proteins, so for a record in
that set we can return the published answer instead of recomputing an approximation of it -- and, where
we do recompute, we can say whether we agree with the published tool.

Measured coverage on the two production manifests: **13 of 140** and **14 of 155** candidates are in the
database by accession. That is small, and it is concentrated where our own inducer call is weakest --
GntR, TetR and MarR, the organic-ligand families for which the SSN clade layer says little and the
coordination gate is silent.

Two join keys, in order:

  1. **RefSeq accession** -- exact, and the normal case for a candidate discovered from a genome.
  2. **Protein sequence (MD5)** -- for a query pasted as a bare sequence with no accession, but only
     when that sequence occurs in exactly one published record. Ligify derives its answer from genomic
     context, and the database contains identical proteins in different contexts with different ligand
     calls. Those ambiguous sequence-only joins must abstain; an exact accession still resolves them.

There is deliberately NO fuzzy match. A near-neighbour is not the same protein, and silently attributing
a published biosensor's ligand to a homolog is exactly the error this module exists to avoid making.

    python -m predictor.effector.ligify_db --acc WP_003963520.1
"""
from __future__ import annotations

import gzip
import hashlib
import json
from dataclasses import dataclass, field

from predictor import resources

DB_PATH = resources.REFS_DIR / "ligify_db.json.gz"

_CACHE: dict | None = None


@dataclass
class LigifyRecord:
    """One published Ligify prediction, as we read it."""

    refseq: str = ""
    uniprot_id: str = ""
    annotation: str = ""
    organism: str = ""
    equation: str = ""
    rhea_id: str = ""
    #: Ligify's own operon-context score, 0-100, on the same scale `effector.rank.calculate_rank`
    #: produces -- so ours and theirs are directly comparable rather than merely similar-looking.
    rank: int | None = None
    rank_metrics: dict = field(default_factory=dict)
    ligands: tuple = ()                 # ordered candidate ligand names
    ligand_detail: tuple = ()           # the full {name, smiles, iupac} records
    matched_by: str = ""                # "accession" | "sequence"

    @property
    def ligand(self) -> str | None:
        return self.ligands[0] if self.ligands else None


def _load() -> dict:
    """Load exact accession and unambiguous sequence indexes.

    Absent is a legitimate state: the database is optional package data and a source that cannot reach
    its data must abstain in the open rather than fail the run."""
    global _CACHE
    if _CACHE is not None:
        return _CACHE
    empty = {"by_acc": {}, "by_md5": {}, "ambiguous_md5": {}, "meta": {}}
    if not DB_PATH.is_file():
        _CACHE = empty
        return _CACHE
    try:
        with gzip.open(DB_PATH, "rt", encoding="utf-8") as fh:
            blob = json.load(fh)
    except Exception:
        _CACHE = empty
        return _CACHE
    by_acc, sequence_groups = {}, {}
    for r in blob.get("records") or []:
        acc = (r.get("refseq") or "").strip()
        if acc:
            by_acc.setdefault(acc, r)
            by_acc.setdefault(acc.split(".", 1)[0], r)      # version-insensitive
        seq = "".join(str(r.get("protein_seq") or "").split()).upper()
        if seq:
            digest = hashlib.md5(seq.encode()).hexdigest()
            sequence_groups.setdefault(digest, []).append(r)

    # Ligify predictions are context-derived. The same protein sequence can occur beside different
    # enzymes and therefore have different published ligand calls; choosing the first record would be
    # deterministic but scientifically arbitrary. Sequence-only lookup is safe only for singleton
    # groups. Accession lookup above remains available for every record.
    by_md5 = {digest: records[0] for digest, records in sequence_groups.items()
              if len(records) == 1}
    ambiguous_md5 = {digest: tuple(records) for digest, records in sequence_groups.items()
                     if len(records) > 1}
    _CACHE = {"by_acc": by_acc, "by_md5": by_md5, "ambiguous_md5": ambiguous_md5,
              "meta": {k: v for k, v in blob.items() if k != "records"}}
    return _CACHE


def available() -> bool:
    return bool(_load()["by_acc"])


def meta() -> dict:
    """Provenance of the vendored copy: source, licence, upstream sha256, retrieval date."""
    return dict(_load()["meta"])


def _to_record(raw: dict, how: str) -> LigifyRecord:
    ligs = raw.get("candidate_ligands") or []
    names = tuple(x.get("name") for x in ligs if isinstance(x, dict) and x.get("name"))
    rk = raw.get("rank") or {}
    return LigifyRecord(
        refseq=raw.get("refseq") or "", uniprot_id=raw.get("uniprot_id") or "",
        annotation=raw.get("annotation") or "",
        organism=(raw.get("organism") if isinstance(raw.get("organism"), str)
                  else " ".join(raw.get("organism") or [])),
        equation=raw.get("equation") or "", rhea_id=str(raw.get("rhea_id") or ""),
        rank=rk.get("rank") if isinstance(rk, dict) else None,
        rank_metrics=(rk.get("metrics") or {}) if isinstance(rk, dict) else {},
        ligands=names, ligand_detail=tuple(ligs), matched_by=how)


def lookup(seq: str | None = None, acc: str | None = None) -> LigifyRecord | None:
    """The published Ligify prediction for this TF, or None.

    Accession first (exact identity), then an unambiguous sequence hash. Nothing fuzzy and no arbitrary
    choice among context-dependent records: see the module docstring.
    """
    db = _load()
    if acc:
        a = str(acc).strip()
        raw = db["by_acc"].get(a) or db["by_acc"].get(a.split(".", 1)[0])
        if raw:
            return _to_record(raw, "accession")
    if seq:
        s = "".join(str(seq).split()).upper()
        raw = db["by_md5"].get(hashlib.md5(s.encode()).hexdigest())
        if raw:
            return _to_record(raw, "sequence")
    return None


# --------------------------------------------------------------------------- self-test
def _demo() -> None:
    db = _load()
    if not db["by_acc"]:
        print("SKIP: no vendored Ligify database (run tools/vendor_ligify_db.py)")
        return
    m = meta()
    print(f"Ligify DB: {len(db['by_md5'])} unambiguous sequences, "
          f"{len(db['ambiguous_md5'])} ambiguous, retrieved {m.get('_retrieved')}, "
          f"upstream {m.get('_upstream_records')} records")
    assert m.get("_license"), "the redistributed database must carry its MIT notice"
    assert m.get("_source", "").startswith("https://"), m.get("_source")

    # round-trip a real record by BOTH keys, and confirm they resolve to the same protein
    raw = next(iter(db["by_md5"].values()))
    by_seq = lookup(seq=raw["protein_seq"])
    by_acc = lookup(acc=raw["refseq"])
    assert by_seq and by_acc and by_seq.refseq == by_acc.refseq == raw["refseq"]
    assert by_acc.matched_by == "accession" and by_seq.matched_by == "sequence"
    print(f"  {by_acc.refseq}: ligand={by_acc.ligand!r} rank={by_acc.rank} "
          f"({by_acc.annotation[:44]})")

    # a version-less accession still resolves; an unknown one abstains rather than guessing
    assert lookup(acc=raw["refseq"].split(".", 1)[0]) is not None
    assert lookup(acc="WP_000000000.1") is None
    assert lookup(seq="MKKKNOTAREALPROTEINSEQUENCEATALL") is None
    assert lookup() is None
    print("OK: accession + sequence lookup, version-insensitive, no fuzzy matching, abstains cleanly.")


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--acc")
    ap.add_argument("--seq")
    a = ap.parse_args()
    if a.acc or a.seq:
        r = lookup(seq=a.seq, acc=a.acc)
        print(json.dumps(r.__dict__, indent=2, default=str) if r else "not in the Ligify database")
    else:
        _demo()
