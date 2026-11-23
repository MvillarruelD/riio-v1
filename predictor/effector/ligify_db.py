"""
ligify_db.py -- look a TF up in the PUBLISHED Ligify sensor database before re-deriving its effector.

`effector/operon_chemistry.py` derives candidates from a divergent partner's operon and neighbouring
enzyme reactions in UniProt. It uses the formula ported in `effector/rank.py`, with a different
neighbourhood model from the published tool. `effector/ligify.py` is now only a compatibility shim.
The published and derived routes can share underlying evidence; agreement is not automatically
independent corroboration.

What was missing was their DATA. `groov-bio/ligify-ui` ships a static database of 3,164 predicted
biosensors, each keyed by RefSeq accession and carrying the protein sequence, the candidate ligands, the
Rhea reaction and Ligify's own rank. Our discovery candidates ARE RefSeq proteins, so for a record in
that set we can retain the published prediction -- and, where we also derive candidates, report
agreement or disagreement. These are computational predictions, not experimentally validated truth.

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
import re
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
    #: Ligify's uncalibrated operon-context score (can be negative). The local route uses the same
    #: formula on a different gene set: equal scores do not establish equal accuracy or compound odds.
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


# --------------------------------------------------------------------------- near-identical homologues
#: EXP071: published records of NEAR-IDENTICAL proteins. Measured on this database over records that
#: carry ligands (EXP071 type-concordance table): homologous records share an exact published ligand
#: 89 % of the time at >= 99 % identity, 85 % at 80-90 %, 60 % at 50-60 %, against 5.5 % for random
#: pairs. Informative, but a context-derived prediction, not this protein's validated inducer.
#: (EXP070's "70 %" counted ligand-less records; corrected.) EXP071 kept it non-voting; EXP072 lets the
#: >= 99 % tier vote for organic ligands -- see below.
HOMOLOG_MIN_IDENTITY = 0.99
HOMOLOG_MIN_COVERAGE = 0.90
HOMOLOG_MAX = 3
#: EXP072 (owner decision 2026-11-09): a published homologue is informative for ORGANIC effectors, where
#: the SSN has no comparable handle; ions come from the SSN analysis. So
#:   * a record at >= HOMOLOG_MIN_IDENTITY may VOTE, for its organic ligands only;
#:   * a record between HOMOLOG_CONTEXT_MIN_IDENTITY and HOMOLOG_MIN_IDENTITY is context, never counted;
#:   * a metal / monatomic-ion ligand is dropped from both tiers (`inorganic_reason`), and the drop logged.
#: The 0.99 vote threshold was chosen deliberately high: lower values reach real cases on the E. coli
#: panel, but those hits had been inspected before the threshold was set, so that panel cannot choose it.
HOMOLOG_CONTEXT_MIN_IDENTITY = 0.50
#: Share of homologous record pairs at >= 99 % identity that carry an exact published ligand in common
#: (39 of 44; EXP071 `type_concordance.json`). Scales the vote of a homologue relative to an exact record.
HOMOLOG_VOTE_CONCORDANCE = 0.886

#: One bracketed atom with an optional charge: `[Na+]`, `[Fe+2]`, `[Hg]`, `[Cl-]`. Not `[SH-]`.
_MONATOMIC_SMILES = re.compile(r"\[\d*[A-Z][a-z]?(?:[+-]+|[+-]\d+)?\]")
#: Element symbol of every bracketed SMILES atom (`[Fe+3]`, `[Mo]`, `[Co-3]`); metals are always bracketed.
_BRACKET_ATOM = re.compile(r"\[\d*([A-Z][a-z]?)")


def inorganic_reason(name, smiles=None) -> str:
    """Why a published ligand is a metal or monatomic ion ('' when it is not one).

    Names in the database are prose (`Sodium ion`, `Ferrous cation`) that the identity vocabulary does
    not resolve, so the SMILES is the reliable test and the name class is the fallback."""
    from predictor.effector import chemical_identity as _ci
    if smiles and _MONATOMIC_SMILES.fullmatch(str(smiles).strip()):
        return "monatomic_smiles"
    # A metal anywhere in the structure: Fe-S clusters, hemes, cobalamins, molybdate, arsenite, and
    # metal salts such as sodium butyrate. Found by the EXP072 review -- single atoms alone let 187
    # database records keep a metal-containing "organic" ligand. Carbon-free non-metals (H2O2, nitrate)
    # are deliberately NOT dropped: the rule is about metals, which are the SSN's business.
    if smiles and any(a in _ci._METAL_ELEMENTS for a in _BRACKET_ATOM.findall(str(smiles))):
        return "metal_in_smiles"
    text = _ci.transport_name(name)
    if _ci.chemistry(name)[0] == "metal":
        return "metal_class"
    if text and (_ci._monatomic_ion(text) or _ci._bare_metal_symbol(text)):
        return "monatomic_ion_name"
    return ""

_RECORDS: list | None = None


def _all_records() -> list:
    global _RECORDS
    if _RECORDS is None:
        try:
            with gzip.open(DB_PATH, "rt", encoding="utf-8") as fh:
                _RECORDS = list(json.load(fh).get("records") or [])
        except Exception:
            _RECORDS = []
    return _RECORDS


def _records_fasta():
    """One FASTA of every record with a sequence, header `r<index>`; built once per vendored copy."""
    recs = _all_records()
    tag = (_load()["meta"].get("_upstream_sha256") or "nometa")[:12]
    path = resources.cache_path("ligify_db", f"ligify_db_records_{tag}.faa")
    if not (path.exists() and path.stat().st_size > 0):
        with resources.atomic_path(path) as temporary, temporary.open("w", encoding="ascii", errors="ignore") as fh:
            for i, r in enumerate(recs):
                s = "".join(str(r.get("protein_seq") or "").split()).upper()
                if s:
                    fh.write(f">r{i}\n{s}\n")
    return path, recs


def homologs(seq: str, *, exclude_refseq=(), min_identity: float = HOMOLOG_MIN_IDENTITY,
             min_coverage: float = HOMOLOG_MIN_COVERAGE, limit: int = HOMOLOG_MAX) -> list:
    """Published records of near-identical proteins: [(LigifyRecord, identity, qcov, tcov)], best first.

    Context, not identity. Records whose RefSeq (version-insensitive) is in `exclude_refseq` -- the
    query's own exact record -- are skipped. Empty when the database or MMseqs2 is unavailable."""
    if not seq or not _all_records():
        return []
    from predictor.annotate.homologs import homolog_sweep
    fasta, recs = _records_fasta()
    excl = {str(x).split(".", 1)[0] for x in exclude_refseq if x}
    sweep = homolog_sweep("".join(str(seq).split()).upper(), str(fasta), cluster=False,
                          evalue=1e-10, max_seqs=50)
    out, seen = [], set()
    for h in sweep.hits:
        if h.pident < min_identity or h.qcov < min_coverage or h.tcov < min_coverage:
            continue
        raw = recs[int(h.target.lstrip("r"))]
        acc = (raw.get("refseq") or "").split(".", 1)[0]
        if acc in excl or acc in seen:
            continue
        seen.add(acc)
        out.append((_to_record(raw, "homolog"), round(h.pident, 4), round(h.qcov, 3), round(h.tcov, 3)))
    out.sort(key=lambda t: (-t[1], -t[2]))
    return out[:limit]


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
