"""
ligify.py -- DEPRECATED compatibility shim. The derived genome-context route was renamed to
`effector/operon_chemistry.py` and now emits `source="operon_chem"` (EXP066), because EXP065 measured
that only Ligify's scoring FORMULA was theirs -- the neighbourhood was a strand-blind flat box, not
Ligify's promoter-gated operon walk. The published-database lookup keeps the name Ligify
(`effector/ligify_db.py`): that one is genuinely their answer.

This module re-exports the new API unchanged for one release so existing imports keep working. New code
should import `predictor.effector.operon_chemistry` directly. It will be removed once every on-disk
bundle and analysis script has moved off the `"ligify"` source string (readers still accept both names
during the transition -- see `effector/inducer.py` and `effector/ligify_db.py`).

NOTE: `effectors_from_context` no longer takes `max_kb` / `max_genes` (the flat window is gone); callers
that passed them will raise, which is deliberate -- those windows no longer exist.
"""
from __future__ import annotations

from predictor.effector.operon_chemistry import (                     # noqa: F401 - re-exported
    EffectorCandidate,
    MAX_PROMOTER_BP,
    OPERON_CHEM_MAX_GAP_BP,
    OPERON_CHEM_MAX_SPAN_BP,
    divergent_partner,
    operon_from_partner,
    effectors_from_context,
    predict_effector,
    _resolve_context,
)

__all__ = [
    "EffectorCandidate", "MAX_PROMOTER_BP", "OPERON_CHEM_MAX_GAP_BP", "OPERON_CHEM_MAX_SPAN_BP",
    "divergent_partner", "operon_from_partner", "effectors_from_context", "predict_effector",
]
