"""
operator_score.py — Composite "is this an operator?" score for one ESMFold2
TF-DNA fold, for the operator-scanner pipeline (scripts 30-33).

The premise (see plan): ESMFold2 docks a wHTH/MerR dimer onto almost any B-DNA,
so a *single* fold's contact count is not, by itself, evidence of an operator.
What should distinguish the true operator window from a decoy window is the
*quality* of the engagement:

  • more direct base-edge (readout) contacts, concentrated (localised) rather
    than smeared along the construct;
  • read out through the MAJOR groove (winged-HTH / MerR recognition helix);
  • coming from the recognition helix specifically;
  • symmetric across both protomers of the homodimer (a real dyad operator);
  • a confident protein-DNA interface (high ipTM, low interface-PAE);
  • appropriate DNA distortion (MerR kinks the operator).

This module turns one result folder into a transparent feature vector plus a
composite score. Weights are interpretable defaults and are meant to be
*validated / refit* on the Phase-A benchmark (script 33), never trusted blind.
Every score ships with its per-feature breakdown and (when a null is supplied)
a z-score against shuffled-DNA decoys.

Reuses: _lib.af3_parse (folder parsing, pLDDT, PAE), _lib.contacts
(geometric contact detection, SxxLxxL motif), _lib.dna_geometry (groove +
bend/kink).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

from . import af3_parse
from . import dna_geometry
from .contacts import find_contacts, find_motif_SxxLxxL

log = logging.getLogger(__name__)

_DNA_RESIDUES = {
    "DA", "DT", "DG", "DC", "A", "T", "G", "C",
    "DA3", "DA5", "DT3", "DT5", "DG3", "DG5", "DC3", "DC5",
}

# Recognition-helix span (residues counted from helix start) treated as the
# DNA-reading helix for "recognition-helix engagement".
RECOG_HELIX_SPAN = 12


# ── Feature container ────────────────────────────────────────────────────────

@dataclass
class FoldScore:
    tf_id: str
    folder: str
    construct_len: int

    # raw contacts
    n_base: int = 0
    n_backbone: int = 0
    n_total: int = 0

    # dimer symmetry
    chain_a_base: int = 0
    chain_b_base: int = 0
    dimer_symmetry: float = 0.0          # 0..1 (min/max of the two protomers' base contacts)

    # recognition-helix engagement
    recog_helix_base: int = 0
    recog_helix_base_frac: float = 0.0

    # groove readout (base contacts only)
    groove_major: int = 0
    groove_minor: int = 0
    groove_other: int = 0
    groove_major_frac: float = 0.0       # major / (major+minor)

    # DNA geometry (localised to contacted region)
    bend_deg: Optional[float] = None
    max_kink_deg: Optional[float] = None

    # interface confidence
    iptm: Optional[float] = None
    ptm: Optional[float] = None
    interface_pae: Optional[float] = None   # mean protein↔DNA PAE (lower=better)
    mean_plddt_contacts: Optional[float] = None

    # localisation
    contact_centroid: Optional[float] = None   # top-strand 0-based position
    locality: Optional[float] = None           # SD (bp) of contacted positions; small=focused

    # per-position evidence on the construct top strand (for footprint aggregation)
    per_pos_base: List[float] = field(default_factory=list)
    per_pos_backbone: List[float] = field(default_factory=list)

    composite: Optional[float] = None
    notes: List[str] = field(default_factory=list)

    def feature_dict(self) -> Dict[str, float]:
        """The numeric features used by the composite / benchmark (NaN-safe)."""
        def f(x):
            return float(x) if x is not None else float("nan")
        return {
            "n_base": f(self.n_base),
            "n_backbone": f(self.n_backbone),
            "dimer_symmetry": f(self.dimer_symmetry),
            "recog_helix_base_frac": f(self.recog_helix_base_frac),
            "groove_major_frac": f(self.groove_major_frac),
            "bend_deg": f(self.bend_deg),
            "max_kink_deg": f(self.max_kink_deg),
            "iptm": f(self.iptm),
            "interface_pae": f(self.interface_pae),
            "locality": f(self.locality),
        }


# Interpretable default weights (refit on benchmark). Sign convention: higher
# composite = more operator-like. interface_pae and locality are penalties
# (negative weights). Features are scaled inside ``composite_score``.
DEFAULT_WEIGHTS: Dict[str, float] = {
    "n_base": 1.0,                 # per direct base contact
    "groove_major_frac": 6.0,      # major-groove readout fraction (0..1)
    "recog_helix_base_frac": 5.0,  # recognition-helix fraction (0..1)
    "dimer_symmetry": 4.0,         # half-site symmetry (0..1)
    "iptm": 6.0,                   # interface confidence (0..1)
    "interface_pae": -0.20,        # Å of PAE (penalty)
    "bend_deg": 0.03,             # mild reward for localised distortion
    "locality": -0.15,            # SD in bp (penalty for smeared contacts)
}


def composite_score(feat: FoldScore, weights: Dict[str, float] = None) -> float:
    """Transparent weighted sum over feature_dict(); NaN features contribute 0."""
    w = weights or DEFAULT_WEIGHTS
    fd = feat.feature_dict()
    total = 0.0
    for k, wk in w.items():
        v = fd.get(k, float("nan"))
        if v == v:  # not NaN
            total += wk * v
    return float(total)


# ── helix detection (mirrors 11_dna_contacts.detect_recognition_helix) ───────

# Family-specific recognition-helix search defaults (residue numbering).
#   ArsR/SmtB: α3 reading helix carries the SxxLxxL motif (~res 30-55).
#   MerR:      α2 reading helix is near the N-terminus (~res 10-35); documented
#              base-contacting residues e.g. CueR K15/F19, ZntR similar.
_HELIX_DEFAULTS = {
    "arsr": (28, 60),
    "arsr_smtb": (28, 60),
    "merr": (8, 38),
    "merra": (8, 38),
}


def detect_helix_start(
    structure,
    protein_chains: Sequence[str],
    *,
    family: str = "auto",
    res_start: int = 40,
    res_end: int = 80,
) -> Optional[int]:
    """Return the recognition-helix start residue for relative i+ numbering.

    ArsR/SmtB: SxxLxxL α3 motif. MerR: the α2 recognition helix is near the
    N-terminus, so we anchor at the family-default start. 'auto' tries the
    SxxLxxL motif first and falls back to res_start.
    """
    # Apply family default search window unless the caller overrode the default.
    if family in _HELIX_DEFAULTS and (res_start, res_end) == (40, 80):
        res_start, res_end = _HELIX_DEFAULTS[family]
    if family in ("merr", "merra"):
        return res_start
    residues_in_range, seen = [], set()
    for chain in structure.get_chains():
        if chain.id not in protein_chains:
            continue
        for res in chain.get_residues():
            resi = res.id[1]
            if res_start <= resi <= res_end:
                key = (chain.id, resi)
                if key not in seen:
                    residues_in_range.append(
                        {"chain": chain.id, "resi": resi, "resn": res.resname.strip()})
                    seen.add(key)
    residues_in_range.sort(key=lambda r: (r["chain"], r["resi"]))
    motif = find_motif_SxxLxxL(residues_in_range)
    if motif is not None:
        return motif
    return res_start if family != "auto" else res_start


# ── interface PAE ────────────────────────────────────────────────────────────

def _interface_pae(model: Dict, protein_chains, dna_chains) -> Optional[float]:
    """Mean PAE over the protein↔DNA token block (both directions)."""
    pae = model.get("pae")
    if pae is None or not hasattr(pae, "shape") or pae.ndim != 2:
        return None
    folder = Path(model["source_path"])
    full_files = af3_parse._find_files(folder, af3_parse._FULL_GLOBS)
    tci: List[str] = []
    for fp in full_files:
        try:
            tci = af3_parse.extract_token_chain_ids(af3_parse.parse_full_data_json(fp))
        except Exception:
            tci = []
        if tci:
            break
    if not tci or len(tci) != pae.shape[0]:
        return None
    tci = np.asarray(tci)
    prot_mask = np.isin(tci, list(protein_chains))
    dna_mask = np.isin(tci, list(dna_chains))
    if not prot_mask.any() or not dna_mask.any():
        return None
    block1 = pae[np.ix_(prot_mask, dna_mask)]
    block2 = pae[np.ix_(dna_mask, prot_mask)]
    return float(np.mean(np.concatenate([block1.ravel(), block2.ravel()])))


# ── main entry point ─────────────────────────────────────────────────────────

def score_fold(
    folder: Path,
    tf_id: str,
    *,
    plddt_min: float = 50.0,
    max_dist: float = 3.5,
    family: str = "auto",
    res_start: int = 40,
    res_end: int = 80,
    weights: Dict[str, float] = None,
) -> Optional[FoldScore]:
    """Score one ESMFold2 tf_dna result folder. Returns None if unparseable."""
    from Bio.PDB import MMCIFParser

    models = af3_parse.parse_result_folder(folder, tf_id, "tf_dna")
    if not models:
        log.warning("score_fold: no models in %s", folder)
        return None
    model = models[0]
    cif_path = model["cif_path"]

    try:
        struct = MMCIFParser(QUIET=True).get_structure(tf_id, str(cif_path))
    except Exception as exc:  # noqa: BLE001
        log.warning("score_fold: CIF parse failed %s: %s", cif_path, exc)
        return None

    # chain inventory
    protein_chains, dna_chains = [], []
    for chain in struct.get_chains():
        is_dna = any(r.resname.strip() in _DNA_RESIDUES for r in chain.get_residues())
        (dna_chains if is_dna else protein_chains).append(chain.id)
    top_chain = "C" if "C" in dna_chains else (dna_chains[0] if dna_chains else "C")
    bot_chain = "D" if "D" in dna_chains else (dna_chains[1] if len(dna_chains) > 1 else None)

    plddt_arr = model.get("plddt_per_residue")
    plddt_chains = model.get("plddt_chain_ids") or []
    L = plddt_chains.count(top_chain) or 80

    helix_start = detect_helix_start(
        struct, protein_chains, family=family, res_start=res_start, res_end=res_end)

    contacts, _excl = find_contacts(
        struct, plddt_arr, max_dist=max_dist, plddt_min=plddt_min,
        operator_center_pos=None, helix_start_resi=helix_start)

    feat = FoldScore(tf_id=tf_id, folder=str(folder), construct_len=L)
    feat.per_pos_base = [0.0] * L
    feat.per_pos_backbone = [0.0] * L

    def top_pos(c) -> Optional[int]:
        """0-based top-strand position of a DNA contact."""
        n = int(c.dna_position)
        if c.dna_chain == top_chain:
            x = n - 1
        else:
            x = L - n          # bottom resi n ↔ top index L-n
        return x if 0 <= x < L else None

    base_positions: List[int] = []
    plddt_vals: List[float] = []
    for c in contacts:
        x = top_pos(c)
        if c.interaction_type == "direct_base":
            feat.n_base += 1
            if x is not None:
                feat.per_pos_base[x] += 1.0
                base_positions.append(x)
            if c.chain_id == "A":
                feat.chain_a_base += 1
            elif c.chain_id == "B":
                feat.chain_b_base += 1
            if c.residue_i is not None and 0 <= c.residue_i < RECOG_HELIX_SPAN:
                feat.recog_helix_base += 1
            g = dna_geometry.classify_groove(c.dna_atom, c.dna_base)
            if g == "major":
                feat.groove_major += 1
            elif g == "minor":
                feat.groove_minor += 1
            else:
                feat.groove_other += 1
        else:
            feat.n_backbone += 1
            if x is not None:
                feat.per_pos_backbone[x] += 1.0
        if c.plddt_protein is not None:
            plddt_vals.append(c.plddt_protein)

    feat.n_total = feat.n_base + feat.n_backbone

    # derived fractions
    hi, lo = max(feat.chain_a_base, feat.chain_b_base), min(feat.chain_a_base, feat.chain_b_base)
    feat.dimer_symmetry = round(lo / hi, 3) if hi else 0.0
    feat.recog_helix_base_frac = round(feat.recog_helix_base / feat.n_base, 3) if feat.n_base else 0.0
    gm = feat.groove_major + feat.groove_minor
    feat.groove_major_frac = round(feat.groove_major / gm, 3) if gm else 0.0
    feat.mean_plddt_contacts = round(float(np.mean(plddt_vals)), 2) if plddt_vals else None

    # localisation: centroid + SD over base-contacted positions (fallback: all)
    pos_for_loc = base_positions or [
        i for i, v in enumerate(feat.per_pos_backbone) if v > 0]
    if pos_for_loc:
        feat.contact_centroid = round(float(np.mean(pos_for_loc)), 2)
        feat.locality = round(float(np.std(pos_for_loc)), 2)

    # interface confidence
    feat.iptm = model.get("iptm")
    feat.ptm = model.get("ptm")
    feat.interface_pae = _interface_pae(model, protein_chains, dna_chains)

    # geometry, localised to the contacted region (±3 bp) to avoid floppy ends
    if pos_for_loc:
        c0, c1 = min(pos_for_loc), max(pos_for_loc)
        resi_window = (max(1, c0 - 3 + 1), min(L, c1 + 3 + 1))  # +1: 0-based→1-based resi
    else:
        resi_window = None
    try:
        geo = dna_geometry.axis_geometry(
            cif_path, top_chain=top_chain, bottom_chain=bot_chain or "D",
            resi_window=resi_window)
        feat.bend_deg = geo.bend_deg
        feat.max_kink_deg = geo.max_kink_deg
    except Exception as exc:  # noqa: BLE001
        feat.notes.append(f"geometry failed: {exc}")

    feat.composite = round(composite_score(feat, weights), 3)
    return feat


# ── null normalisation ───────────────────────────────────────────────────────

def zscore_vs_null(observed: float, null_values: Sequence[float]) -> Optional[float]:
    """z = (observed - mean(null)) / sd(null); None if null too small/degenerate."""
    vals = [v for v in null_values if v == v]
    if len(vals) < 2:
        return None
    mu, sd = float(np.mean(vals)), float(np.std(vals))
    if sd == 0:
        return None
    return round((observed - mu) / sd, 3)


def dinuc_shuffle(seq: str, rng: np.random.Generator) -> str:
    """Dinucleotide-preserving shuffle (Altschul-Erikson) — destroys motif,
    preserves base + nearest-neighbour composition for a fair decoy null."""
    seq = seq.upper()
    if len(seq) < 4:
        return seq
    # build edge lists per nucleotide
    from collections import defaultdict
    edges = defaultdict(list)
    for a, b in zip(seq[:-1], seq[1:]):
        edges[a].append(b)
    # shuffle each adjacency list
    for k in edges:
        rng.shuffle(edges[k])
    # Eulerian-walk reconstruction with retries for connectivity
    for _ in range(20):
        e = {k: list(v) for k, v in edges.items()}
        for k in e:
            rng.shuffle(e[k])
        out = [seq[0]]
        cur = seq[0]
        ok = True
        for _ in range(len(seq) - 1):
            if not e.get(cur):
                ok = False
                break
            nxt = e[cur].pop()
            out.append(nxt)
            cur = nxt
        if ok and len(out) == len(seq):
            return "".join(out)
    return seq  # give up → return original (caller still gets a valid string)
