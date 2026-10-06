"""
af3_parse.py — Parse AF3 server result folders.

AF3 server delivers one ZIP per job; after extraction each folder contains:
  fold_<jobname>_model_0.cif … model_4.cif       — ranked structure files
  fold_<jobname>_summary_confidences_0.json …    — per-model summary
  fold_<jobname>_full_data_0.json …              — per-atom pLDDT, PAE matrix
  fold_<jobname>_job_request.json                — submitted request
  msas/, templates/

Job folders on disk are named with the LOWERCASED job name plus a version
suffix:  arsr_cgarsr_apo_v1, arsr_cgarsr_apo_v1_2 (retry), cuer_ecoli_dna_v1, …

This module:
  • discovers result folders (handles _v1, _v1_2 retry naming, prefers latest)
  • maps the lowercased folder TF prefix back to the canonical mixed-case
    TF ID by scanning transcription_factors/**/*.yaml
  • extracts per-residue pLDDT from CIF B-factors (the authoritative source —
    AF3 sets B-factor = pLDDT, so this avoids the token/atom mapping problem)
  • exposes PAE, ranking_score, pTM, ipTM, chain layout
"""

from __future__ import annotations

import json
import logging
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

log = logging.getLogger(__name__)

# ── Glob patterns (order = preference) ────────────────────────────────────────
_CIF_GLOBS = ["*model*.cif", "*.cif"]
_PDB_GLOBS = ["*model*.pdb", "*.pdb"]
_SUMMARY_GLOBS = ["*summary_confidences*.json", "*summary*.json"]
_FULL_GLOBS = ["*full_data*.json", "*confidences*.json", "*confidence*.json"]

# Rank extracted from filename:  fold_<job>_model_2.cif → 2
_RANK_PATTERNS = [
    re.compile(r"_rank_(\d+)"),
    re.compile(r"_model_(\d+)\.(?:cif|pdb)$"),
    re.compile(r"_(\d+)\.(?:cif|pdb|json)$"),
]

# Folder name pattern:  <tf_lower>_<state>_v<run>[_<retry>]
#   examples: cuer_ecoli_apo_v1
#             arsr_cgarsr_dna_v1_2     (retry of a re-submitted job)
#             smtb_synechococcus_holo_v3
_FOLDER_RE = re.compile(
    r"^(?P<tf>.+?)_(?P<state>apo|dimer|holo|ligand|dna|tf_dna|apo_dimer|holo_dimer)_v(?P<run>\d+)(?:_(?P<retry>\d+))?$"
)

_STATE_TO_PRED_TYPE = {
    "apo": "apo_dimer",
    "dimer": "apo_dimer",
    "apo_dimer": "apo_dimer",
    "holo": "holo_dimer",
    "ligand": "holo_dimer",
    "holo_dimer": "holo_dimer",
    "dna": "tf_dna",
    "tf_dna": "tf_dna",
}


# ── TF ID canonicalisation (lowercase folder → mixed-case YAML id) ────────────

_TF_ID_CACHE: Optional[Dict[str, str]] = None


def _build_tf_id_lookup(base_dir: Path) -> Dict[str, str]:
    """
    Scan transcription_factors/**/*.yaml for canonical TF ids.

    Returns: lower(tf_id) → tf_id  e.g.  {'arsr_cgarsr': 'ArsR_CgArsR', ...}
    """
    global _TF_ID_CACHE
    if _TF_ID_CACHE is not None:
        return _TF_ID_CACHE

    lookup: Dict[str, str] = {}
    tf_dir = base_dir / "transcription_factors"
    if not tf_dir.exists():
        log.warning("transcription_factors/ not found at %s", tf_dir)
        _TF_ID_CACHE = lookup
        return lookup

    # Use cheap yaml parsing for just the 'id:' line — avoids loading PyYAML
    # if not present, and is fast.
    for yaml_path in tf_dir.rglob("*.yaml"):
        try:
            with open(yaml_path, encoding="utf-8") as fh:
                for line in fh:
                    line = line.strip()
                    if line.startswith("id:"):
                        tf_id = line.split(":", 1)[1].strip()
                        if tf_id:
                            lookup[tf_id.lower()] = tf_id
                        break
        except Exception as e:
            log.debug("Skipped %s: %s", yaml_path, e)

    _TF_ID_CACHE = lookup
    log.debug("TF lookup built: %d entries", len(lookup))
    return lookup


def _canonical_tf_id(tf_lower: str, base_dir: Path) -> str:
    """Map lowercased TF id back to the canonical (mixed-case) form."""
    lookup = _build_tf_id_lookup(base_dir)
    return lookup.get(tf_lower, tf_lower)


# ── Helper utilities ──────────────────────────────────────────────────────────

def _rank_from_filename(name: str) -> int:
    for pat in _RANK_PATTERNS:
        m = pat.search(name)
        if m:
            return int(m.group(1))
    return 0


def _find_files(folder: Path, globs: List[str]) -> List[Path]:
    for g in globs:
        matches = sorted(folder.glob(g))
        if matches:
            return matches
    return []


# ── Result folder detection ──────────────────────────────────────────────────

def find_result_folders(af3_results_dir: Path) -> List[Tuple[str, str, Path]]:
    """
    Scan af3_results/ for completed job folders.

    Returns list of (tf_id, pred_type, folder_path).

    Handles:
      • lowercase folder names → canonical mixed-case TF id
      • _v1 / _v1_2 / _v2 version suffixes (retries) → keep the latest
      • multiple state aliases (apo|dimer; holo|ligand; dna|tf_dna)
    """
    results: Dict[Tuple[str, str], Tuple[int, int, Path]] = {}
    if not af3_results_dir.exists():
        log.warning("af3_results/ not found: %s", af3_results_dir)
        return []

    base_dir = af3_results_dir.parent

    for folder in sorted(af3_results_dir.iterdir()):
        if not folder.is_dir():
            continue
        m = _FOLDER_RE.match(folder.name)
        if not m:
            log.debug("Skipping unrecognised result folder: %s", folder.name)
            continue

        tf_lower = m.group("tf")
        state = m.group("state")
        run_idx = int(m.group("run"))
        retry = int(m.group("retry") or 0)

        pred_type = _STATE_TO_PRED_TYPE.get(state)
        if pred_type is None:
            log.debug("Unknown state suffix in %s", folder.name)
            continue

        tf_id = _canonical_tf_id(tf_lower, base_dir)
        key = (tf_id, pred_type)

        # Keep the latest folder for each (tf, type): higher run_idx wins,
        # then higher retry wins.
        score = (run_idx, retry)
        if key not in results or score > results[key][:2]:
            results[key] = (run_idx, retry, folder)

    out: List[Tuple[str, str, Path]] = []
    for (tf_id, pred_type), (_, _, folder) in sorted(results.items()):
        out.append((tf_id, pred_type, folder))
    log.info("Discovered %d unique result folders in %s", len(out), af3_results_dir)
    return out


# ── Single-model parsing ──────────────────────────────────────────────────────

def parse_summary_json(path: Path) -> Dict[str, Any]:
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def parse_full_data_json(path: Path) -> Dict[str, Any]:
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def extract_pae_matrix(full_data: Dict[str, Any]) -> Optional[np.ndarray]:
    """Return the PAE matrix as an NxN float32 array (token-level), or None."""
    pae_raw = full_data.get("pae")
    if pae_raw is None:
        return None
    mat = np.array(pae_raw, dtype=np.float32)
    if mat.ndim != 2:
        log.warning("PAE matrix has unexpected shape: %s", mat.shape)
        return None
    return mat


def extract_token_chain_ids(full_data: Dict[str, Any]) -> List[str]:
    """Per-token chain id list (length = N tokens = PAE matrix size)."""
    return list(full_data.get("token_chain_ids") or [])


def extract_token_res_ids(full_data: Dict[str, Any]) -> List[int]:
    """Per-token residue id list (1-based)."""
    return [int(x) for x in (full_data.get("token_res_ids") or [])]


def extract_chain_breaks_from_tokens(full_data: Dict[str, Any]) -> List[int]:
    """0-based token indices where a new chain begins."""
    chain_ids = full_data.get("token_chain_ids") or full_data.get("atom_chain_ids")
    if not chain_ids:
        return []
    breaks = [0]
    for i in range(1, len(chain_ids)):
        if chain_ids[i] != chain_ids[i - 1]:
            breaks.append(i)
    return breaks


# ── pLDDT extraction from CIF (authoritative source) ──────────────────────────

def extract_plddt_from_cif(cif_path: Path) -> Tuple[Optional[np.ndarray], List[str], List[int]]:
    """
    Read per-residue mean pLDDT from CIF B-factors (AF3 sets B = pLDDT).

    Returns:
      plddt_per_residue : 1-D float array (one entry per residue in CIF order)
      chain_ids         : list of chain IDs aligned with plddt_per_residue
      res_ids           : list of residue numbers aligned with plddt_per_residue

    Robust to missing biopython (returns None on failure).
    """
    try:
        from Bio.PDB import MMCIFParser
    except ImportError:
        log.warning("BioPython not installed; cannot extract per-residue pLDDT from CIF")
        return None, [], []

    try:
        parser = MMCIFParser(QUIET=True)
        struct = parser.get_structure(cif_path.stem, str(cif_path))
    except Exception as e:
        log.warning("Failed to parse CIF %s: %s", cif_path.name, e)
        return None, [], []

    plddts: List[float] = []
    chain_ids: List[str] = []
    res_ids: List[int] = []
    for model in struct:
        for chain in model:
            for residue in chain:
                # Skip waters / het flags that aren't part of the chain
                hetflag = residue.id[0].strip()
                if hetflag and hetflag != "H_MSE":  # keep selenomet, but exclude HOH
                    continue
                b_factors = [atom.get_bfactor() for atom in residue if atom.element != "H"]
                if not b_factors:
                    continue
                plddts.append(float(np.mean(b_factors)))
                chain_ids.append(chain.id)
                res_ids.append(int(residue.id[1]))
        break  # Only first model in CIF (AF3 = one model per file)

    if not plddts:
        return None, [], []
    return np.array(plddts, dtype=float), chain_ids, res_ids


# ── Folder-level parsing ──────────────────────────────────────────────────────

def parse_result_folder(folder: Path, tf_id: str, pred_type: str) -> List[Dict[str, Any]]:
    """
    Parse all models in an AF3 result folder.

    Returns a list of model-dicts, sorted by rank (ascending).
    Each dict has:
      tf_id, pred_type, model_rank, source_path, cif_path,
      mean_plddt, ranking_score, ptm, iptm,
      has_dna, has_ligand, n_chains,
      plddt_per_residue (np.ndarray | None),
      plddt_chain_ids    (List[str]),
      plddt_res_ids      (List[int]),
      pae                (np.ndarray | None),
      chain_breaks       (List[int]),
    """
    cif_files = _find_files(folder, _CIF_GLOBS) or _find_files(folder, _PDB_GLOBS)
    summary_files = _find_files(folder, _SUMMARY_GLOBS)
    full_files = _find_files(folder, _FULL_GLOBS)

    if not cif_files:
        log.warning("No structure file found in %s", folder)
        return []

    models = []

    for cif_path in cif_files:
        rank = _rank_from_filename(cif_path.name)

        summary_data: Dict[str, Any] = {}
        full_data: Dict[str, Any] = {}

        for sp in summary_files:
            if _rank_from_filename(sp.name) == rank:
                try:
                    summary_data = parse_summary_json(sp)
                except Exception as e:
                    log.warning("Failed to parse summary %s: %s", sp, e)
                break

        for fp in full_files:
            if _rank_from_filename(fp.name) == rank:
                try:
                    full_data = parse_full_data_json(fp)
                except Exception as e:
                    log.warning("Failed to parse full_data %s: %s", fp, e)
                break

        # Per-residue pLDDT from CIF (authoritative)
        plddt_arr, plddt_chains, plddt_resids = extract_plddt_from_cif(cif_path)

        pae_mat = extract_pae_matrix(full_data) if full_data else None
        chain_breaks = extract_chain_breaks_from_tokens(full_data) if full_data else []

        # Global mean pLDDT: prefer atom-level from full_data, fall back to per-res CIF
        atom_plddts = full_data.get("atom_plddts")
        if atom_plddts:
            mean_plddt = float(np.mean(atom_plddts))
        elif plddt_arr is not None and len(plddt_arr):
            mean_plddt = float(np.mean(plddt_arr))
        else:
            mean_plddt = None

        n_chains = (
            len(set(full_data.get("token_chain_ids") or []))
            if full_data.get("token_chain_ids")
            else (len(set(plddt_chains)) if plddt_chains else None)
        )

        model = {
            "tf_id": tf_id,
            "pred_type": pred_type,
            "model_rank": rank,
            "source_path": str(folder),
            "cif_path": cif_path,
            "mean_plddt": mean_plddt,
            "ranking_score": summary_data.get("ranking_score"),
            "ptm": summary_data.get("ptm"),
            "iptm": summary_data.get("iptm"),
            "has_dna": pred_type == "tf_dna",
            "has_ligand": pred_type == "holo_dimer",
            "n_chains": n_chains,
            "plddt_per_residue": plddt_arr,
            "plddt_chain_ids": plddt_chains,
            "plddt_res_ids": plddt_resids,
            "pae": pae_mat,
            "chain_breaks": chain_breaks,
            "summary_data": summary_data,
            "full_data_keys": list(full_data.keys()) if full_data else [],
        }
        models.append(model)

    models.sort(key=lambda m: m["model_rank"])
    return models


# Backwards compat shim — some scripts import these
def extract_plddt_per_residue(full_data: Dict[str, Any]) -> Optional[np.ndarray]:
    """
    Deprecated. Per-residue pLDDT is now extracted from CIF B-factors.
    Returns atom_plddts averaged into uniform-size residue bins as a fallback.
    """
    if "plddt" in full_data:
        return np.array(full_data["plddt"], dtype=float)
    return None


def extract_chain_breaks(full_data: Dict[str, Any]) -> List[int]:
    """Deprecated alias — kept for old callers."""
    return extract_chain_breaks_from_tokens(full_data)
