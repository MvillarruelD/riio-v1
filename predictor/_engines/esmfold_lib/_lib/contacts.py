"""
contacts.py — Geometric TF–DNA contact detection (BioPython back-end).

Uses BioPython NeighborSearch (KD-tree) to find protein *sidechain* polar atoms
within max_dist (default 3.5 Å) of DNA atoms.

Design matches the reference PyMOL script (Analyze_contacts_afArsR):
  - Only SIDECHAIN polar atoms (N, O, S) — backbone N/CA/C/O excluded.
  - DNA atoms classified as: "base" | "sugar" | "phosphate"
    (equivalent to the reference script's base/sugar/phosphate triage).
  - Relative position indexing: protein i+ (relative to recognition helix
    start), DNA pos (relative to operator center).

Contact dataclass
-----------------
  Matches the reference CSV column layout plus pLDDT fields added by the
  BioPython standalone pipeline.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

log = logging.getLogger(__name__)

# ── amino acid conversion ─────────────────────────────────────────────────────
_THREE_TO_ONE: Dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def three_to_one(resn: str) -> str:
    return _THREE_TO_ONE.get(str(resn).upper()[:3], "?")


# ── atom sets ─────────────────────────────────────────────────────────────────

# Backbone atoms — EXCLUDED from contact detection (only sidechain polars counted)
_BACKBONE_NAMES: Set[str] = {"N", "CA", "C", "O", "OXT"}

# Sidechain polar atoms (N, O, S family — by element, not by name)
# We filter by element in code; this set is kept for reference/documentation.
_SIDECHAIN_POLAR_ELEMENTS: Set[str] = {"N", "O", "S"}

# DNA base atoms (ring N/O/C that can form H-bonds — no prime in name)
_DNA_BASE_ATOMS: Set[str] = {
    "N1", "N3", "N7", "N9",          # purine/pyrimidine ring N
    "O2", "O4", "O6",                # carbonyl oxygens
    "N2", "N4", "N6",                # exocyclic amines
    "C2", "C4", "C5", "C6", "C8",   # ring carbons (for CH-π etc.)
}

# DNA sugar atoms
_DNA_SUGAR_ATOMS: Set[str] = {
    "C1'", "C2'", "C3'", "C4'", "C5'",
    "O4'", "O5'", "O3'", "O2'",
}

# DNA phosphate atoms
_DNA_PHOSPHATE_ATOMS: Set[str] = {"P", "OP1", "OP2", "O1P", "O2P"}

# Standard DNA residue names
_DNA_RESIDUES: Set[str] = {
    "DA", "DT", "DG", "DC", "A", "T", "G", "C",
    "DA3", "DA5", "DT3", "DT5", "DG3", "DG5", "DC3", "DC5",
}


# ── contact dataclass ─────────────────────────────────────────────────────────

@dataclass
class Contact:
    """
    A single protein (sidechain) – DNA contact.

    Column layout mirrors the reference CSV format:
      ProteinChain, ResidueAbs, Residue_i, ResidueAA_3, ResidueAA_1,
      ProtAtom, ProtAtomElem,
      DNAChain, DNAResName, DNAAbs, DNA_pos, DNAAtom, DNAAtomClass,
      Distance_A, pLDDT_protein, pLDDT_dna
    """
    # Protein
    chain_id: str               # protein chain id (e.g. "A", "B")
    residue_number: int         # absolute PDB residue number
    residue_i: Optional[int]    # relative to recognition helix start (i+ numbering); None if unknown
    tf_residue: str             # 3-letter residue name
    aa_one: str                 # 1-letter residue code
    atom: str                   # atom name
    atom_elem: str              # element (N/O/S)
    # DNA
    dna_chain: str
    dna_base: str               # residue name (DA, DT, DG, DC)
    dna_position: int           # absolute PDB sequence number
    dna_pos: Optional[int]      # relative to operator center; None if unknown
    dna_atom: str
    dna_atom_class: str         # "base" | "sugar" | "phosphate" | "unknown"
    # Geometry / quality
    interaction_type: str       # "direct_base" | "backbone" (sugar/phosphate)
    distance_A: float
    plddt_protein: Optional[float] = None
    plddt_dna: Optional[float] = None


# ── DNA atom classification ───────────────────────────────────────────────────

def classify_dna_atom(atom_name: str, resn: str = "", element: str = "") -> str:
    """
    Classify a DNA atom into 'base', 'sugar', or 'phosphate'.

    Matches the logic from the reference PyMOL script (classify_dna_atom).
    """
    an = atom_name.strip()

    if an in _DNA_BASE_ATOMS:
        return "base"
    if an in _DNA_SUGAR_ATOMS or "'" in an:
        return "sugar"
    if an in _DNA_PHOSPHATE_ATOMS or element.upper() == "P":
        return "phosphate"

    # Fallback: if it's a known DNA residue, use letter/prime heuristic
    r = (resn or "").upper()
    if r in _DNA_RESIDUES:
        if any(ch.isalpha() and ch.isupper() for ch in an) and "'" not in an:
            return "base"
        return "sugar"

    return "unknown"


def classify_dna_contact(dna_atom_name: str) -> str:
    """Return 'direct_base' if atom is a base atom, else 'backbone'."""
    cls = classify_dna_atom(dna_atom_name)
    return "direct_base" if cls == "base" else "backbone"


# ── recognition helix motif detection ────────────────────────────────────────

def find_motif_SxxLxxL(
    residues_in_range: List[Dict],
) -> Optional[int]:
    """
    Search for the ArsR/SmtB α3-recognition-helix motif SxxLxxL in a list
    of residue dicts (each with keys 'resi', 'resn', 'chain').

    Returns the absolute residue number of the Ser (start of motif), or None.
    """
    if not residues_in_range:
        return None

    seq1 = "".join(three_to_one(r["resn"]) for r in residues_in_range)
    pat = re.compile(r"S..L..L")
    m = pat.search(seq1)
    if m:
        start_idx = m.start()
        return residues_in_range[start_idx]["resi"]
    return None


# ── pLDDT helper ─────────────────────────────────────────────────────────────

def _plddt_for_residue(
    plddt_arr: Optional[np.ndarray],
    residue_idx: int,
) -> Optional[float]:
    if plddt_arr is None or residue_idx < 0 or residue_idx >= len(plddt_arr):
        return None
    return float(plddt_arr[residue_idx])


# ── main contact finder ───────────────────────────────────────────────────────

def find_contacts(
    structure,                           # Bio.PDB.Structure
    plddt_arr: Optional[np.ndarray] = None,
    max_dist: float = 3.5,
    plddt_min: float = 70.0,
    operator_center_pos: Optional[int] = None,    # 1-based DNA position for DNA_pos=0
    helix_start_resi: Optional[int] = None,       # 1-based protein resi for residue_i=0
) -> Tuple[List[Contact], List[Dict]]:
    """
    Detect protein–DNA contacts in an AF3 structure.

    Only SIDECHAIN polar atoms (element N, O, or S; not backbone N/CA/C/O)
    are considered — matching the reference PyMOL script's approach.

    Parameters
    ----------
    structure : Bio.PDB.Structure
    plddt_arr : per-residue pLDDT (same ordering as structure residues)
    max_dist : Å cutoff for contact detection
    plddt_min : minimum pLDDT for protein residues to be included
    operator_center_pos : 1-based DNA position used as DNA_pos=0 reference
    helix_start_resi : 1-based protein residue used as residue_i=0 reference

    Returns
    -------
    contacts, excluded_residues
    """
    try:
        from Bio.PDB import NeighborSearch
    except ImportError:
        raise ImportError("BioPython is required: pip install biopython")

    # ── build residue index for pLDDT lookup ──────────────────────────────────
    res_list: List = []
    res_index_map: Dict[Tuple[str, int, str], int] = {}

    protein_atoms: List[Tuple] = []   # (atom, res, chain_id)
    dna_atoms: List[Tuple] = []       # (atom, res, chain_id)

    for chain in structure.get_chains():
        for res in chain.get_residues():
            resname = res.resname.strip()
            res_key = (chain.id, res.id[1], resname)
            idx = len(res_list)
            res_list.append(res)
            res_index_map[res_key] = idx

            is_dna = resname in _DNA_RESIDUES

            for atom in res.get_atoms():
                atom_name = atom.name.strip()
                elem = getattr(atom, "element", None) or ""
                elem = elem.strip().upper() or atom_name[0].upper()

                if is_dna:
                    dna_atoms.append((atom, res, chain.id))
                else:
                    # Sidechain polar filter: element N/O/S AND NOT backbone
                    if elem in _SIDECHAIN_POLAR_ELEMENTS and atom_name not in _BACKBONE_NAMES:
                        protein_atoms.append((atom, res, chain.id))

    if not dna_atoms:
        log.info("No DNA residues found in structure; skipping contact detection")
        return [], []
    if not protein_atoms:
        log.info("No protein sidechain polar atoms found; skipping contact detection")
        return [], []

    log.debug(
        "Contact search: %d sidechain polar atoms vs %d DNA atoms (cutoff=%.1f Å)",
        len(protein_atoms), len(dna_atoms), max_dist,
    )

    # ── KD-tree on DNA atoms ───────────────────────────────────────────────────
    dna_atom_objects = [a for a, _, _ in dna_atoms]
    ns = NeighborSearch(dna_atom_objects)

    contacts: List[Contact] = []
    excluded: List[Dict] = []
    seen_low_plddt: Set[Tuple] = set()

    for p_atom, p_res, p_chain in protein_atoms:
        p_resname = p_res.resname.strip()
        p_res_key = (p_chain, p_res.id[1], p_resname)
        p_res_idx = res_index_map.get(p_res_key, -1)
        p_plddt = _plddt_for_residue(plddt_arr, p_res_idx)

        # pLDDT filter
        if p_plddt is not None and p_plddt < plddt_min:
            if p_res_key not in seen_low_plddt:
                excluded.append({
                    "tf_residue": p_resname,
                    "residue_number": p_res.id[1],
                    "chain_id": p_chain,
                    "plddt": p_plddt,
                    "reason": f"pLDDT {p_plddt:.1f} < threshold {plddt_min}",
                })
                seen_low_plddt.add(p_res_key)
            continue

        # Relative i+ numbering
        res_i: Optional[int] = None
        if helix_start_resi is not None:
            res_i = p_res.id[1] - helix_start_resi

        p_atom_name = p_atom.name.strip()
        p_elem = (getattr(p_atom, "element", None) or p_atom_name[0]).strip().upper()

        # Find nearby DNA atoms
        nearby = ns.search(p_atom.coord, max_dist, "A")
        for d_atom in nearby:
            d_res = d_atom.get_parent()
            d_chain_obj = d_res.get_parent()
            d_chain = d_chain_obj.id
            d_resname = d_res.resname.strip()
            d_res_key = (d_chain, d_res.id[1], d_resname)
            d_res_idx = res_index_map.get(d_res_key, -1)
            d_plddt = _plddt_for_residue(plddt_arr, d_res_idx)

            d_atom_name = d_atom.name.strip()
            d_elem = (getattr(d_atom, "element", None) or d_atom_name[0]).strip().upper()

            cls = classify_dna_atom(d_atom_name, d_resname, d_elem)
            itype = "direct_base" if cls == "base" else "backbone"

            # DNA relative position
            dna_pos: Optional[int] = None
            if operator_center_pos is not None:
                dna_pos = d_res.id[1] - operator_center_pos

            dist = round(float(p_atom - d_atom), 3)

            contacts.append(Contact(
                chain_id=p_chain,
                residue_number=p_res.id[1],
                residue_i=res_i,
                tf_residue=p_resname,
                aa_one=three_to_one(p_resname),
                atom=p_atom_name,
                atom_elem=p_elem,
                dna_chain=d_chain,
                dna_base=d_resname,
                dna_position=d_res.id[1],
                dna_pos=dna_pos,
                dna_atom=d_atom_name,
                dna_atom_class=cls,
                interaction_type=itype,
                distance_A=dist,
                plddt_protein=p_plddt,
                plddt_dna=d_plddt,
            ))

    log.info(
        "Found %d contacts (%d direct_base, %d backbone); "
        "excluded %d low-pLDDT protein residues (pLDDT_min=%.0f)",
        len(contacts),
        sum(1 for c in contacts if c.interaction_type == "direct_base"),
        sum(1 for c in contacts if c.interaction_type == "backbone"),
        len(excluded),
        plddt_min,
    )
    return contacts, excluded


# ── binding summary ───────────────────────────────────────────────────────────

def summarize_binding(
    contacts: List[Contact],
    operator_seq: str,
    pad_offset: int = 0,
) -> Dict:
    """
    Summarize binding from a contact list relative to an operator sequence.

    pad_offset : 0-based index of operator start in the 60-bp padded construct.
    """
    if not contacts or not operator_seq:
        return {}

    direct = [c for c in contacts if c.interaction_type == "direct_base"]
    dna_positions = sorted({c.dna_position for c in direct})

    op_len = len(operator_seq)
    op_start = pad_offset + 1   # 1-based
    op_end = pad_offset + op_len

    contacted_in_op = [p for p in dna_positions if op_start <= p <= op_end]
    overlap_pct = 100.0 * len(contacted_in_op) / op_len if op_len > 0 else 0.0

    # Chain A vs chain B (half-site symmetry)
    chain_a = {(c.residue_number, c.dna_position) for c in direct if c.chain_id == "A"}
    chain_b = {(c.residue_number, c.dna_position) for c in direct if c.chain_id == "B"}

    return {
        "n_direct_contacts": len(direct),
        "n_backbone_contacts": len(contacts) - len(direct),
        "dna_positions_contacted": ",".join(str(p) for p in dna_positions),
        "positions_in_operator": ",".join(str(p) for p in contacted_in_op),
        "operator_coverage_pct": round(overlap_pct, 1),
        "chain_a_contacts": len(chain_a),
        "chain_b_contacts": len(chain_b),
        "symmetric_contacts": len(chain_a & chain_b),
    }
