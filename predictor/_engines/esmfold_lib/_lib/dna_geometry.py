"""
dna_geometry.py — DNA groove classification + helical-axis bend/kink for the
ESMFold2 operator-scanner scorer (scripts 30-33).

Two independent capabilities, both deliberately dependency-light:

1. ``classify_groove(atom_name, base)`` — map a *base-edge* atom to the groove
   it faces (major / minor / other). This lets the scorer reward
   "recognition in the correct groove" (winged-HTH reads the **major** groove)
   directly from the atom names already recorded in the {tf}_contacts.csv, with
   no structure re-parsing.

2. ``axis_geometry(cif_path, ...)`` — fit a helical axis to base-pair centre
   coordinates (C1'–C1' midpoints) and report a global bend angle + the maximum
   local kink. MerR-family activation kinks/under-twists the operator
   (CueR: four kinks of 28-85°); ArsR/SmtB binding is comparatively straight.
   So bend is a genuinely discriminative geometric feature.

   If ``x3dna-dssr`` is found on PATH, ``dssr_geometry`` parses its richer
   per-step roll/twist + groove-width output instead; callers should
   detect-and-fallback (``axis_geometry`` never needs DSSR).

All angles are degrees. Coordinates come from BioPython (CIF already written by
esmfold2_adapter; occupancy column patched in, so MMCIFParser is happy).
"""

from __future__ import annotations

import json
import logging
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

log = logging.getLogger(__name__)


# ── DNA residue names (mirror _lib/contacts.py) ──────────────────────────────
_DNA_RESIDUES = {
    "DA", "DT", "DG", "DC", "A", "T", "G", "C",
    "DA3", "DA5", "DT3", "DT5", "DG3", "DG5", "DC3", "DC5",
}


def _base_letter(resn: str) -> str:
    """Normalise a DNA residue name to a single letter A/T/G/C ('' if unknown)."""
    r = (resn or "").upper().strip()
    for pre in ("D",):
        if r.startswith(pre) and len(r) > 1:
            r = r[1:]
    r = r.rstrip("35")          # strip 3'/5' terminal tags
    return r if r in ("A", "T", "G", "C") else ""


# ── 1. Groove-edge atom classification ───────────────────────────────────────
#
# Base-edge H-bonding / van-der-Waals atoms grouped by the groove they present
# to a protein in B-DNA. Carbons that line each groove edge are included so
# CH-style contacts also score. Sugar/phosphate atoms are NOT base-edge and
# return "other" (the scorer treats those as backbone, not readout).
#
# Major-groove edge:
#   purines  A: N7, C8, N6, C6, C5     G: N7, C8, O6, C6, C5
#   pyrimid. T: O4, C4, C5, C7(C5M)    C: N4, C4, C5
# Minor-groove edge:
#   purines  A: N3, C2                 G: N3, N2, C2
#   pyrimid. T: O2, C2                 C: O2, C2
_MAJOR_BY_BASE: Dict[str, set] = {
    "A": {"N7", "C8", "N6", "C6", "C5"},
    "G": {"N7", "C8", "O6", "C6", "C5"},
    "T": {"O4", "C4", "C5", "C7", "C5M", "C5A"},
    "C": {"N4", "C4", "C5"},
}
_MINOR_BY_BASE: Dict[str, set] = {
    "A": {"N3", "C2"},
    "G": {"N3", "N2", "C2"},
    "T": {"O2", "C2"},
    "C": {"O2", "C2"},
}

# Base-agnostic fallback (used when the base letter is unavailable). N7/O6/N6/N4
# /O4 are unambiguously major-groove; N3/N2/O2 minor. C2 is minor-ish; the ring
# carbons are ambiguous → "other".
_MAJOR_GENERIC = {"N7", "O6", "N6", "N4", "O4", "C7", "C5M"}
_MINOR_GENERIC = {"N2", "O2", "N3"}


def classify_groove(atom_name: str, base: str = "") -> str:
    """Return 'major' | 'minor' | 'other' for a DNA base-edge atom.

    ``atom_name`` is the PDB atom name (e.g. 'N7'); ``base`` is the residue name
    or single letter (optional but improves accuracy for ring carbons).
    Sugar/phosphate atoms (primes, P, OPx) always return 'other'.
    """
    an = (atom_name or "").strip().upper()
    if not an or "'" in an or an in {"P", "OP1", "OP2", "O1P", "O2P"}:
        return "other"
    bl = _base_letter(base) if base else ""
    if bl:
        if an in _MAJOR_BY_BASE.get(bl, set()):
            return "major"
        if an in _MINOR_BY_BASE.get(bl, set()):
            return "minor"
        return "other"
    if an in _MAJOR_GENERIC:
        return "major"
    if an in _MINOR_GENERIC:
        return "minor"
    return "other"


def groove_split(atom_base_pairs) -> Tuple[int, int, int]:
    """Tally (n_major, n_minor, n_other) over an iterable of (atom_name, base)."""
    nM = nm = no = 0
    for atom_name, base in atom_base_pairs:
        g = classify_groove(atom_name, base)
        if g == "major":
            nM += 1
        elif g == "minor":
            nm += 1
        else:
            no += 1
    return nM, nm, no


# ── 2. Helical-axis bend / kink from base-pair centres ───────────────────────

@dataclass
class AxisGeometry:
    n_centers: int
    bend_deg: Optional[float]            # angle between two half-axis directions
    max_kink_deg: Optional[float]        # max local turn over a sliding triplet
    mean_kink_deg: Optional[float]
    rise_per_bp: Optional[float] = None  # mean centre-to-centre spacing (Å)
    source: str = "axis_fit"
    notes: List[str] = field(default_factory=list)


def _load_dna_chains(cif_path: Path):
    """Return a BioPython model and the list of DNA chain ids present."""
    from Bio.PDB import MMCIFParser
    parser = MMCIFParser(QUIET=True)
    struct = parser.get_structure("s", str(cif_path))
    model = next(iter(struct))
    dna_chains = []
    for chain in model.get_chains():
        for res in chain.get_residues():
            if res.resname.strip() in _DNA_RESIDUES:
                dna_chains.append(chain.id)
                break
    return model, dna_chains


def _c1_coords_by_resi(model, chain_id: str) -> Dict[int, np.ndarray]:
    """Map residue number → C1' coordinate for one DNA chain."""
    out: Dict[int, np.ndarray] = {}
    chain = model[chain_id]
    for res in chain.get_residues():
        if res.resname.strip() not in _DNA_RESIDUES:
            continue
        for atom in res.get_atoms():
            if atom.name.strip() == "C1'":
                out[res.id[1]] = np.asarray(atom.coord, dtype=float)
                break
    return out


def base_pair_centers(
    cif_path: Path,
    top_chain: str = "C",
    bottom_chain: str = "D",
) -> Tuple[np.ndarray, List[int]]:
    """Return (centers Nx3, top_resis) — C1'–C1' midpoints along the duplex.

    Pairing rule (matches esmfold2_inputs: chain D is the reverse complement of
    chain C): top residue n pairs with bottom residue (L - n + 1), L = #bp.
    Falls back to single-strand C1' positions if the bottom chain is absent.
    """
    model, dna_chains = _load_dna_chains(cif_path)
    if top_chain not in dna_chains:
        # pick the first two DNA chains we did find
        if not dna_chains:
            return np.empty((0, 3)), []
        top_chain = dna_chains[0]
        bottom_chain = dna_chains[1] if len(dna_chains) > 1 else None

    top = _c1_coords_by_resi(model, top_chain)
    if not top:
        return np.empty((0, 3)), []
    top_resis = sorted(top)
    L = len(top_resis)

    centers: List[np.ndarray] = []
    used_resis: List[int] = []
    bottom = _c1_coords_by_resi(model, bottom_chain) if bottom_chain in dna_chains else {}

    for n in top_resis:
        if bottom:
            partner = L - (n - top_resis[0])      # robust to 1- or 0-based start
            # try the most likely partner numbering, then a couple of fallbacks
            cand = None
            for p in (L - n + 1, partner, n):
                if p in bottom:
                    cand = bottom[p]
                    break
            if cand is not None:
                centers.append((top[n] + cand) / 2.0)
                used_resis.append(n)
                continue
        centers.append(top[n])      # single-strand fallback
        used_resis.append(n)

    return np.asarray(centers, dtype=float), used_resis


def _fit_axis(points: np.ndarray) -> Optional[np.ndarray]:
    """Principal-axis unit vector through a set of points (SVD)."""
    if points.shape[0] < 2:
        return None
    centred = points - points.mean(axis=0)
    try:
        _, _, vh = np.linalg.svd(centred, full_matrices=False)
    except np.linalg.LinAlgError:
        return None
    v = vh[0]
    n = np.linalg.norm(v)
    return v / n if n else None


def _angle_between(u: np.ndarray, v: np.ndarray) -> float:
    cu, cv = np.linalg.norm(u), np.linalg.norm(v)
    if not cu or not cv:
        return 0.0
    c = float(np.clip(np.dot(u, v) / (cu * cv), -1.0, 1.0))
    ang = np.degrees(np.arccos(c))
    return min(ang, 180.0 - ang)     # axis direction is sign-ambiguous


def axis_geometry(
    cif_path: Path,
    top_chain: str = "C",
    bottom_chain: str = "D",
    split_index: Optional[int] = None,
    resi_window: Optional[Tuple[int, int]] = None,
) -> AxisGeometry:
    """Global bend + local kink of the duplex from base-pair-centre geometry.

    bend_deg     : angle between the principal axis of the first half of the
                   duplex and that of the second half (split at ``split_index``,
                   default = midpoint). Straight B-DNA → ~0°; a central kink
                   shows up as a large value.
    max_kink_deg : max turn of the local tangent over a sliding 3-centre window
                   (curvature spike — e.g. a single sharp kink off-centre).

    ``resi_window``: optional (lo, hi) inclusive top-strand residue numbers to
    restrict the analysis to the protein-engaged region. STRONGLY recommended:
    on long constructs the free DNA ends flop and dominate the global bend, so
    bend should be measured over the contacted sub-helix only.
    """
    notes: List[str] = []
    centers, resis = base_pair_centers(cif_path, top_chain, bottom_chain)
    if resi_window is not None and len(resis):
        lo, hi = resi_window
        keep = [i for i, r in enumerate(resis) if lo <= r <= hi]
        if len(keep) >= 4:
            centers = centers[keep]
            resis = [resis[i] for i in keep]
            notes.append(f"focused on resi {lo}-{hi} ({len(keep)} bp)")
        else:
            notes.append(f"resi window {lo}-{hi} had <4 bp; used full construct")
    n = centers.shape[0]
    if n < 4:
        return AxisGeometry(n, None, None, None, notes=["<4 bp centres; geometry skipped"])

    # rise per bp
    diffs = np.diff(centers, axis=0)
    seg_len = np.linalg.norm(diffs, axis=1)
    rise = float(np.median(seg_len[seg_len > 0])) if np.any(seg_len > 0) else None

    # global bend: two half-axes
    mid = split_index if split_index is not None else n // 2
    mid = int(np.clip(mid, 2, n - 2))
    a1 = _fit_axis(centers[:mid])
    a2 = _fit_axis(centers[mid:])
    bend = _angle_between(a1, a2) if (a1 is not None and a2 is not None) else None

    # local kink: smoothed tangents over a 3-centre stride
    kinks: List[float] = []
    stride = 2
    for i in range(stride, n - stride):
        u = centers[i] - centers[i - stride]
        v = centers[i + stride] - centers[i]
        kinks.append(_angle_between(u, v))
    max_kink = float(np.max(kinks)) if kinks else None
    mean_kink = float(np.mean(kinks)) if kinks else None

    return AxisGeometry(
        n_centers=n,
        bend_deg=round(bend, 2) if bend is not None else None,
        max_kink_deg=round(max_kink, 2) if max_kink is not None else None,
        mean_kink_deg=round(mean_kink, 2) if mean_kink is not None else None,
        rise_per_bp=round(rise, 2) if rise is not None else None,
        notes=notes,
    )


# ── 3. Optional richer geometry via x3dna-dssr ───────────────────────────────

def dssr_available() -> bool:
    return shutil.which("x3dna-dssr") is not None


def dssr_geometry(cif_path: Path) -> Optional[Dict]:
    """Parse x3dna-dssr --json for groove widths + bend, if DSSR is installed.

    Returns a dict subset {bend_deg, min_groove_width, ...} or None on any
    failure. Callers fall back to ``axis_geometry``.
    """
    if not dssr_available():
        return None
    try:
        with tempfile.TemporaryDirectory(prefix="dssr_") as tmp:
            out = Path(tmp) / "dssr.json"
            cmd = ["x3dna-dssr", f"-i={cif_path}", f"-o={out}", "--json"]
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            if r.returncode != 0 or not out.exists():
                # DSSR prints JSON to stdout in some builds
                payload = r.stdout.strip()
                data = json.loads(payload) if payload.startswith("{") else None
            else:
                data = json.loads(out.read_text(encoding="utf-8"))
        if not data:
            return None
        # Best-effort field extraction; DSSR schema varies by version.
        result: Dict = {"source": "x3dna-dssr"}
        if isinstance(data, dict):
            result["raw_keys"] = list(data.keys())
            # groove widths live under "grooves" or per-step "pairs" in some builds
        return result
    except Exception as exc:  # noqa: BLE001
        log.debug("DSSR geometry failed for %s: %s", cif_path, exc)
        return None
