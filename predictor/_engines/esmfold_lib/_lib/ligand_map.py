"""
ligand_map.py — Map TF.ligand.primary strings to AF3 server entity specs.

AF3 server JSON schema v1 entity formats (per README_AF3_JSON.md):
  Ion    : {"ion":    {"ion": "ZN", "count": 2}}
  Ligand : {"ligand": {"ligand": "CCD_B12", "count": 2}}   (CCD_ prefix required)
  SMILES : NOT supported in server JSON upload (only CCD codes allowed)

Matching is case-insensitive substring; first match wins (ordered config list).

Resolution precedence (per handoff design decisions):
  1. Direct AF3 ion  → {"ion": {"ion": CODE, "count": N}}
  2. Unsupported metal → closest supported ion (ZN/CO/CU) — substitution logged
  3. CCD cofactor    → {"ligand": {"ligand": "CCD_CODE", "count": N}}
  4. Omit            → no entity added; job tagged "no-effector"
  5. SMILES          → treated as omit (not supported in server JSON format)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

log = logging.getLogger(__name__)


@dataclass
class LigandSpec:
    """Resolved specification for a single ligand entity in an AF3 job."""
    kind: str                          # "ion" | "ccd" | "smiles" | "omit"
    code: Optional[str] = None         # ion code or CCD code (without CCD_ prefix)
    smiles: Optional[str] = None
    count: int = 2                     # copies per dimer
    substitution_reason: Optional[str] = None
    warning: Optional[str] = None


def _build_af3_entity(spec: LigandSpec) -> Optional[Dict[str, Any]]:
    """
    Convert a LigandSpec to a complete AF3 server sequence entity dict.

    The count field is embedded INSIDE the inner dict per the AF3 server JSON
    schema v1 (README reference):
      - Ion:     {"ion":    {"ion": "ZN",       "count": 2}}
      - Ligand:  {"ligand": {"ligand": "CCD_B12", "count": 2}}

    Returns None for 'omit' or unsupported kinds.
    """
    if spec.kind == "ion":
        return {"ion": {"ion": spec.code, "count": spec.count}}

    if spec.kind == "ccd":
        code = spec.code or ""
        # CCD_ prefix required; don't double-prefix if already present
        ccd_code = code if code.upper().startswith("CCD_") else f"CCD_{code}"
        return {"ligand": {"ligand": ccd_code, "count": spec.count}}

    if spec.kind == "smiles":
        # SMILES is NOT supported in the AF3 server JSON format (README v1).
        # The server only accepts the specific CCD codes listed as allowed ligands.
        log.warning(
            "SMILES ligands are not supported in the AF3 server JSON format; "
            "omitting ligand entity (kind='smiles'). Submit this job as no-effector."
        )
        return None

    return None  # omit


def resolve_ligand(
    primary: str,
    config_map: List[Dict[str, Any]],
    copies_per_dimer: int = 2,
) -> LigandSpec:
    """
    Resolve a ligand.primary string against the ordered config_map.

    Parameters
    ----------
    primary : str
        The raw TF.ligand.primary string (e.g. "Ni(II)", "H2S / persulfide").
    config_map : list of dicts
        The ligand_map list from af3_config.yaml (ordered; first match wins).
    copies_per_dimer : int
        Stoichiometry from TF.ligand.copies_per_dimer.

    Returns
    -------
    LigandSpec
    """
    primary_lower = primary.lower().strip()

    for entry in config_map:
        match_str = entry.get("match", "").lower()
        if not match_str:
            continue
        if match_str in primary_lower:
            kind = entry.get("kind", "omit")
            spec = LigandSpec(
                kind=kind,
                code=entry.get("code"),
                smiles=entry.get("smiles"),
                count=copies_per_dimer,
                substitution_reason=entry.get("substitution_reason"),
            )
            if kind == "smiles" and (not spec.smiles or "VERIFY_NEEDED" in str(spec.smiles)):
                spec.kind = "omit"
                spec.warning = (
                    f"SMILES for '{primary}' is unverified/missing; "
                    "holo-dimer will be generated without ligand (no-effector)"
                )
                log.warning(spec.warning)
            if entry.get("substitution_reason"):
                log.info(
                    "Ligand '%s' substituted: %s",
                    primary, entry["substitution_reason"]
                )
            return spec

    # Fallback: nothing matched
    log.warning(
        "No ligand mapping found for '%s'; treating as omit (no-effector)", primary
    )
    return LigandSpec(
        kind="omit",
        count=copies_per_dimer,
        warning=f"No mapping found for '{primary}'; holo-dimer will be no-effector",
    )


def build_ligand_entities(spec: LigandSpec) -> List[Dict[str, Any]]:
    """
    Return a list of AF3 sequence entities for this ligand.

    The count is embedded INSIDE the entity dict by _build_af3_entity —
    do NOT append it at this level.  Returns an empty list for omit.
    """
    if spec.kind == "omit":
        return []
    entity = _build_af3_entity(spec)
    if entity is None:
        return []
    return [entity]
