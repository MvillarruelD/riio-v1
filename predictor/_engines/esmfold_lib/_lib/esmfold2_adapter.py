"""
esmfold2_adapter.py — Write an ESMFold2Result into the AF3 result layout.

Goal: make ESMFold2 outputs invisible to scripts 7–12. Those scripts already
walk ``af3_results/<tf>_<state>_v<n>/`` and parse three filename patterns
(see _lib/af3_parse.py lines 36–40):

    fold_<job>_model_<rank>.cif
    fold_<job>_summary_confidences_<rank>.json
    fold_<job>_full_data_<rank>.json

We write exactly one model (rank 0) per fold call because ESMFold2 returns a
single structure (no AF3-style 5-model ensemble). The folder regex in
af3_parse handles single-rank folders without modification.

Conventions on the JSON contents:

  summary_confidences_0.json
    ranking_score : iptm if present, else ptm, else None
    ptm           : float | None
    iptm          : float | None
    has_clash     : false   (ESMFold2 doesn't expose this; default false)
    fraction_disordered : null

  full_data_0.json
    pae               : nested list, (n_tokens × n_tokens) floats
    token_chain_ids   : list[str], length n_tokens
    token_res_ids     : list[int], 1-based, length n_tokens
    atom_plddts       : list[float], per-residue pLDDT (AF3 emits per-atom;
                        af3_parse reads CIF B-factors as authoritative so the
                        per-residue list is sufficient)
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from .esmfold2_backend import ESMFold2Result

log = logging.getLogger(__name__)


# ── State name conventions (matches AF3 folder naming used by af3_parse) ────
_STATE_FOLDER_TAG = {
    "apo_dimer": "apo",
    "tf_dna": "dna",
    "holo_dimer": "holo",
}
_STATE_JOB_TAG = {
    "apo_dimer": "apo",
    "tf_dna": "dna",
    "holo_dimer": "holo",
}


def folder_name(tf_id: str, pred_type: str, run_idx: int) -> str:
    """Return the folder name (matches the AF3 server's lowercased convention)."""
    tag = _STATE_FOLDER_TAG.get(pred_type, pred_type)
    return f"{tf_id.lower()}_{tag}_v{run_idx}"


def job_name(tf_id: str, pred_type: str, run_idx: int) -> str:
    """Return the <job> stem used in fold_<job>_model_0.cif et al."""
    tag = _STATE_JOB_TAG.get(pred_type, pred_type)
    return f"{tf_id.lower()}_{tag}_v{run_idx}"


@dataclass
class WriteReport:
    folder: Path
    cif_path: Path
    summary_path: Path
    full_data_path: Path
    n_tokens: int
    plddt_mean: Optional[float]
    notes: List[str]


# ── CIF patching ─────────────────────────────────────────────────────────────

def _patch_cif_add_occupancy(cif_text: str) -> str:
    """Insert a dummy _atom_site.occupancy column (value 1.00) into a CIF.

    BioPython's MMCIFParser requires this field; ESMFold2's SDK-generated CIF
    omits it.  We insert the column header right after B_iso_or_equiv in the
    loop_ block, and prepend '1.00' to every ATOM / HETATM data line.

    This function is conservative: if the header insertion fails it returns
    the original text unchanged.
    """
    import re

    # 1. Insert the column name after _atom_site.B_iso_or_equiv
    header_pattern = r"(_atom_site\.B_iso_or_equiv\s*\n)"
    patched, n_subs = re.subn(
        header_pattern,
        r"\1_atom_site.occupancy\n",
        cif_text,
        count=1,
    )
    if n_subs == 0:
        return cif_text  # couldn't find the anchor; give up

    # 2. In each ATOM / HETATM data line, insert '1.00' at the position of
    #    B_iso_or_equiv.  We need to know the column *index* of B_iso_or_equiv
    #    in the loop block so we can insert after it.
    #
    #    Strategy: find the loop_ header block, count which 1-based column
    #    B_iso_or_equiv is, then rewrite each data line.
    loop_header = re.search(
        r"loop_\s*\n((?:_atom_site\.\S+\s*\n)+)",
        patched,
    )
    if loop_header is None:
        return cif_text

    col_names = [
        c.strip()
        for c in re.findall(r"_atom_site\.\S+", loop_header.group(0))
    ]
    try:
        b_idx = col_names.index("_atom_site.B_iso_or_equiv")  # 0-based
    except ValueError:
        return cif_text

    def _insert_occ(m: re.Match) -> str:
        parts = m.group(0).split()
        if len(parts) <= b_idx:
            return m.group(0)
        parts.insert(b_idx + 1, "1.00")
        return " ".join(parts)

    patched = re.sub(r"^(?:ATOM|HETATM)\s+.+$", _insert_occ, patched, flags=re.MULTILINE)
    return patched


# ── Public entry point ───────────────────────────────────────────────────────

def write_af3_layout(
    result: ESMFold2Result,
    *,
    tf_id: str,
    pred_type: str,
    run_idx: int,
    chain_ids: List[str],
    base_dir: Path,
    chain_lengths: Optional[Dict[str, int]] = None,
) -> WriteReport:
    """Write the CIF + two JSONs that scripts 7+ expect.

    Parameters
    ----------
    result
        Normalised SDK output from ``esmfold2_backend.fold``.
    tf_id, pred_type, run_idx
        Used to compose folder + filenames matching the AF3 layout.
    chain_ids
        Ordered chain letters for this job (e.g. ["A","B"] or ["A","B","C","D"]).
        Used to synthesise token_chain_ids / token_res_ids when the SDK does
        not expose them directly.
    base_dir
        Project root; the function writes into ``base_dir / "af3_results"``.
    chain_lengths
        Optional override: {chain_id: n_tokens}. When omitted we derive it
        from the result's pLDDT vector and chain_ids by splitting evenly,
        which is only valid for homodimers (apo_dimer); the driver should
        pass this for tf_dna and holo_dimer.
    """
    notes: List[str] = []

    out_dir = base_dir / "af3_results" / folder_name(tf_id, pred_type, run_idx)
    out_dir.mkdir(parents=True, exist_ok=True)

    job = job_name(tf_id, pred_type, run_idx)
    cif_path = out_dir / f"fold_{job}_model_0.cif"
    summary_path = out_dir / f"fold_{job}_summary_confidences_0.json"
    full_path = out_dir / f"fold_{job}_full_data_0.json"

    # ── 1. CIF ──
    if not result.cif_text:
        notes.append("SDK returned no CIF text; writing empty file (script 7 will skip)")
    cif_text = result.cif_text or ""
    # BioPython's MMCIFParser requires _atom_site.occupancy. ESMFold2's SDK
    # omits it; patch it in so downstream scripts can extract per-residue pLDDT
    # from B-factors without a parsing warning.
    if cif_text and "_atom_site.occupancy" not in cif_text:
        cif_text = _patch_cif_add_occupancy(cif_text)
        if "_atom_site.occupancy" in cif_text:
            notes.append("Patched CIF: added dummy _atom_site.occupancy 1.00")
        else:
            notes.append("CIF occupancy patch failed; BioPython may warn")
    cif_path.write_text(cif_text, encoding="utf-8")

    # Sanity: B-factors plausibly pLDDT?
    if result.plddt is not None and len(result.plddt):
        bf_range = (float(np.nanmin(result.plddt)), float(np.nanmax(result.plddt)))
        if not (0 <= bf_range[0] <= 100 and 0 <= bf_range[1] <= 100):
            notes.append(
                f"pLDDT out of expected 0–100 range: {bf_range}"
            )

    # ── 2. summary_confidences ──
    def _r4(x):
        return round(x, 4) if x is not None else None

    ptm = _r4(result.ptm)
    iptm = _r4(result.iptm)
    ranking = _r4(iptm if result.iptm is not None else result.ptm)
    summary_payload = {
        "ranking_score": ranking,
        "ptm": ptm,
        "iptm": iptm,
        "has_clash": False,
        "fraction_disordered": None,
        # bookkeeping — read by ledger / debugging, ignored by af3_parse
        "_esmfold2": {
            "model": result.model_name,
            "num_loops": result.num_loops,
            "num_sampling_steps": result.num_sampling_steps,
        },
    }
    summary_path.write_text(
        json.dumps(summary_payload, indent=2), encoding="utf-8"
    )

    # ── 3. full_data ──
    token_chain_ids, token_res_ids = _resolve_token_layout(
        result=result,
        chain_ids=chain_ids,
        chain_lengths=chain_lengths,
    )
    n_tokens = len(token_chain_ids)

    pae_list: List[List[float]] = []
    if result.pae is not None and result.pae.size:
        pae_list = result.pae.astype(float).tolist()
        if result.pae.shape[0] != n_tokens:
            notes.append(
                f"PAE size {result.pae.shape[0]} != n_tokens {n_tokens}; "
                "downstream chain-break logic may misalign"
            )

    atom_plddts: List[float] = []
    plddt_mean: Optional[float] = None
    if result.plddt is not None and len(result.plddt):
        raw_vals = [float(x) for x in result.plddt.tolist()]
        # ESMFold2 SDK may return pLDDT in 0–1 range while AF3 stores 0–100.
        # Rescale when all values are ≤ 1.05 (a small tolerance for fp noise).
        if raw_vals and max(raw_vals) <= 1.05:
            notes.append(
                f"pLDDT appears to be in 0–1 range (max={max(raw_vals):.4f}); "
                "rescaling ×100 to match AF3 0–100 convention"
            )
            raw_vals = [v * 100.0 for v in raw_vals]
        atom_plddts = raw_vals
        plddt_mean = float(np.mean(atom_plddts))

    full_payload = {
        "pae": pae_list,
        "token_chain_ids": token_chain_ids,
        "token_res_ids": token_res_ids,
        "atom_plddts": atom_plddts,
    }
    full_path.write_text(
        json.dumps(full_payload, indent=2), encoding="utf-8"
    )

    return WriteReport(
        folder=out_dir,
        cif_path=cif_path,
        summary_path=summary_path,
        full_data_path=full_path,
        n_tokens=n_tokens,
        plddt_mean=plddt_mean,
        notes=notes,
    )


# ── Token layout resolution ──────────────────────────────────────────────────

def _resolve_token_layout(
    *,
    result: ESMFold2Result,
    chain_ids: List[str],
    chain_lengths: Optional[Dict[str, int]],
) -> Tuple[List[str], List[int]]:
    """Return (token_chain_ids, token_res_ids).

    Preference order:
      1. Use the SDK-provided token_chain_ids / token_res_ids if both populated.
      2. Use the explicit chain_lengths mapping if provided.
      3. Fall back to splitting the pLDDT vector evenly across chain_ids
         (only correct for true homodimers — caller's responsibility).
    """
    sdk_chains = result.token_chain_ids
    sdk_res = result.token_res_ids
    if sdk_chains and sdk_res and len(sdk_chains) == len(sdk_res):
        return list(sdk_chains), list(sdk_res)

    total_tokens = 0
    if result.plddt is not None and len(result.plddt):
        total_tokens = int(len(result.plddt))
    elif result.pae is not None and result.pae.size:
        total_tokens = int(result.pae.shape[0])

    token_chain_ids: List[str] = []
    token_res_ids: List[int] = []

    if chain_lengths:
        for cid in chain_ids:
            n = int(chain_lengths.get(cid, 0))
            token_chain_ids.extend([cid] * n)
            token_res_ids.extend(range(1, n + 1))
        return token_chain_ids, token_res_ids

    # Last-ditch even split — only valid for homodimers
    if total_tokens and chain_ids:
        per = total_tokens // len(chain_ids)
        for cid in chain_ids:
            token_chain_ids.extend([cid] * per)
            token_res_ids.extend(range(1, per + 1))

    return token_chain_ids, token_res_ids
