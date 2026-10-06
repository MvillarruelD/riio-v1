"""
esmfold2_inputs.py — Build StructurePredictionInput objects for the biohub SDK.

Mirrors the three job-type builders in af3_schema.py but produces input
objects for ESMFold2 instead of AF3 server JSON dicts. Reuses dna_pad,
ligand_map, and sequence_checks so the data path stays consistent across the
two backends.

Chain id convention (matches the AF3 path so script 11 contact analysis can
remain TF-vs-DNA aware without code changes):

  apo_dimer  : protein A, protein B
  tf_dna     : protein A, protein B, DNA C (sense), DNA D (reverse complement)
  holo_dimer : protein A, protein B, ligand L (chain id chosen by SDK)
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from typing import Any, List, Optional

from .ligand_map import LigandSpec
from .sequence_checks import validate_tf_sequence

log = logging.getLogger(__name__)

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")
_AA_PATTERN = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")
_DNA_PATTERN = re.compile(r"^[ACGT]+$", re.IGNORECASE)


def reverse_complement(seq: str) -> str:
    """Return the reverse complement (uppercase)."""
    return seq.translate(_COMPLEMENT)[::-1].upper()


@dataclass
class BuiltInput:
    """Wraps a StructurePredictionInput with the metadata we need downstream."""
    structure_input: Any                       # esm input_builder.StructurePredictionInput
    chain_ids: List[str]                       # ordered chain letters
    has_dna: bool = False
    has_ligand: bool = False
    notes: List[str] = None                    # informational warnings (e.g. "ligand omitted")

    def __post_init__(self) -> None:
        if self.notes is None:
            self.notes = []


# ── SDK helpers (deferred imports so module is importable without SDK) ──────

def _sdk():
    """Return the input_builder module from the biohub SDK."""
    from esm.utils.structure import input_builder  # type: ignore
    return input_builder


# ── Builders ────────────────────────────────────────────────────────────────

def build_apo_homomer_input(tf_id: str, sequence: str, copies: int = 2) -> BuiltInput:
    """Apo homo-oligomer with ``copies`` identical protein chains.

    ESMFold2's ProteinInput accepts a list of chain IDs to instantiate multiple
    copies of the same sequence. This matters for tetrameric families because
    their inter-dimer interfaces can contain functional sites.
    """
    validate_tf_sequence(tf_id, sequence)
    if not _AA_PATTERN.match(sequence):
        raise ValueError(f"{tf_id}: sequence contains non-standard amino acids")
    if not isinstance(copies, int) or not 1 <= copies <= 26:
        raise ValueError(f"{tf_id}: copies must be an integer in [1, 26], got {copies!r}")

    ib = _sdk()
    chain_ids = [chr(ord("A") + i) for i in range(copies)]
    protein = ib.ProteinInput(id=chain_ids, sequence=sequence)
    structure_input = ib.StructurePredictionInput(sequences=[protein])
    return BuiltInput(
        structure_input=structure_input,
        chain_ids=chain_ids,
        has_dna=False,
        has_ligand=False,
    )


def build_apo_dimer_input(tf_id: str, sequence: str) -> BuiltInput:
    """Backward-compatible two-chain wrapper."""
    return build_apo_homomer_input(tf_id, sequence, copies=2)


def build_tf_dna_input(tf_id: str, sequence: str, dna_top: str) -> BuiltInput:
    """Two protein chains + double-stranded DNA passed as explicit strands.

    Per locked decision: pass sense (chain C) and reverse complement (chain D)
    as two separate ``DNAInput`` entities. Do not rely on SDK auto-pairing.
    """
    validate_tf_sequence(tf_id, sequence)
    if not _AA_PATTERN.match(sequence):
        raise ValueError(f"{tf_id}: sequence contains non-standard amino acids")

    dna_top = dna_top.upper()
    if not _DNA_PATTERN.match(dna_top):
        raise ValueError(f"{tf_id}: DNA top strand is not pure ACGT: {dna_top!r}")
    dna_bottom = reverse_complement(dna_top)

    ib = _sdk()
    protein = ib.ProteinInput(id=["A", "B"], sequence=sequence)
    dna_sense = ib.DNAInput(id="C", sequence=dna_top)
    dna_anti = ib.DNAInput(id="D", sequence=dna_bottom)
    structure_input = ib.StructurePredictionInput(
        sequences=[protein, dna_sense, dna_anti]
    )
    return BuiltInput(
        structure_input=structure_input,
        chain_ids=["A", "B", "C", "D"],
        has_dna=True,
        has_ligand=False,
    )


def build_holo_dimer_input(
    tf_id: str,
    sequence: str,
    ligand_spec: LigandSpec,
) -> BuiltInput:
    """Two protein chains + cognate ligand.

    Maps the LigandSpec produced by ligand_map.resolve_ligand() to an
    ESMFold2 ``LigandInput``:

      - kind=ion   → LigandInput(ccd=code)      # bare metal ion via CCD code
      - kind=ccd   → LigandInput(ccd=code)
      - kind=smiles → LigandInput(smiles=…)
      - kind=omit  → no ligand entity (job degenerates to apo)

    Ion handling is provisional: if the SDK rejects the call the caller
    catches the exception, marks the ledger row, and continues. There is no
    silent AF3 fallback.
    """
    validate_tf_sequence(tf_id, sequence)
    if not _AA_PATTERN.match(sequence):
        raise ValueError(f"{tf_id}: sequence contains non-standard amino acids")

    ib = _sdk()
    sequences: List[Any] = [ib.ProteinInput(id=["A", "B"], sequence=sequence)]
    notes: List[str] = []
    has_ligand = False

    kind = ligand_spec.kind
    code = ligand_spec.code
    smiles = ligand_spec.smiles
    count = max(1, int(ligand_spec.count or 1))

    if kind == "omit":
        notes.append(
            f"holo job for {tf_id} has no ligand (kind=omit, "
            f"reason={ligand_spec.warning or 'covalent/redox effector'})"
        )
    elif kind in ("ion", "ccd") and code:
        # Both routes use a CCD code on the SDK side. Strip a leading "CCD_"
        # if the config map already added one (AF3 path keeps it; ESMFold2
        # SDK doesn't).
        ccd_code = code[4:] if str(code).upper().startswith("CCD_") else str(code)
        ligand_entity = _build_ligand_entity(ib, ccd=ccd_code, count=count)
        sequences.append(ligand_entity)
        has_ligand = True
        if ligand_spec.substitution_reason:
            notes.append(
                f"ligand substituted ({ligand_spec.substitution_reason})"
            )
    elif kind == "smiles" and smiles and "VERIFY_NEEDED" not in str(smiles):
        ligand_entity = _build_ligand_entity(ib, smiles=str(smiles), count=count)
        sequences.append(ligand_entity)
        has_ligand = True
    else:
        notes.append(
            f"holo job for {tf_id}: ligand kind={kind!r} not modelable; "
            f"degenerating to apo"
        )

    structure_input = ib.StructurePredictionInput(sequences=sequences)
    chain_ids = ["A", "B"] + (["L"] if has_ligand else [])
    return BuiltInput(
        structure_input=structure_input,
        chain_ids=chain_ids,
        has_dna=False,
        has_ligand=has_ligand,
        notes=notes,
    )


# ── LigandInput construction helper ─────────────────────────────────────────

def _build_ligand_entity(
    ib,
    *,
    ccd: Optional[str] = None,
    smiles: Optional[str] = None,
    count: int = 1,
):
    """Instantiate a LigandInput across plausible SDK signatures.

    The cookbook excerpt shows ``LigandInput`` exists but the exact constructor
    signature isn't pinned. Try ccd= / smiles= keyword forms in order; if all
    fail, raise so the driver records the row as an error.
    """
    LigandInput = getattr(ib, "LigandInput", None)
    if LigandInput is None:
        raise RuntimeError(
            "ESMFold2 SDK does not expose LigandInput in input_builder; "
            "cannot build holo job."
        )

    # Build positional-id list e.g. ["L"] * count so chain count matches
    # stoichiometry, matching how multimer proteins use id=[...].
    chain_ids = ["L"] if count == 1 else [f"L{i+1}" for i in range(count)]

    last_exc: Optional[Exception] = None
    candidates = []
    if ccd is not None:
        candidates += [
            dict(id=chain_ids, ccd=ccd),
            dict(id=chain_ids[0], ccd=ccd, count=count),
            dict(id=chain_ids[0], ccd_code=ccd, count=count),
        ]
    if smiles is not None:
        candidates += [
            dict(id=chain_ids, smiles=smiles),
            dict(id=chain_ids[0], smiles=smiles, count=count),
        ]
    for kw in candidates:
        try:
            return LigandInput(**kw)
        except TypeError as exc:
            last_exc = exc
            continue
        except Exception as exc:
            last_exc = exc
            break
    raise RuntimeError(
        f"Could not construct LigandInput (ccd={ccd!r}, smiles={smiles!r}); "
        f"last error: {last_exc}"
    )
