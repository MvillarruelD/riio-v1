"""THE inducer-evidence classification, in one place and importable without a run tag.

`classify()` used to live in `fig_evidence_ledger.py`. It is a pure function of one dossier -- it
touches no paths -- but that module imports `run_paths` at import time, which fail-closes unless
TFOP_RUN_TAG is set. That is correct for a figure script, which must belong to exactly one run, and
wrong for anything that legitimately spans runs: the cross-run report tables cover three tags at
once and cannot name a single one.

Moving the function here keeps the property that matters -- ONE definition of "corroborated" -- while
removing the accidental coupling to a single run. `fig_evidence_ledger` imports it from here, so the
figure and any table built on it cannot drift apart.

The rule, unchanged:

    metal  and >=2 of {coordination gate positive, MetalNet site, clade ligand}  -> corroborated
    metal  and  1 of those                                                       -> single-source
    metal  and  0                                                                -> context only
    redox / organic                                                              -> that class
    anything else                                                                -> undetermined

`inducer_class_of` is NOT re-derived here either; it comes from `predictor.effector.inducer_vocab`.
A second copy of the metal/redox/organic vocabulary is how a report ends up disagreeing with the
scorer it quotes.
"""
from __future__ import annotations

#: The four evidence axes reported alongside the outcome. "Spoke" means the axis produced an
#: opinion, not that it agreed -- a silent coordination gate (family chemistry not implemented) is
#: not the same as a negative one.
AXES = [("gate", "gate"), ("metalnet", "MetalNet"), ("clade", "clade"), ("fold", "structure")]


def display_inducer(inducers: dict | None) -> str:
    """Return presentation text for a new or legacy dossier's inducer consensus.

    New bundles persist ``top_display``. Older bundles can derive the same label from their raw
    ``top`` and ``candidates`` fields, which lets tables be refreshed without rerunning predictions.
    Classification and metrics must continue to use ``top`` directly.
    """
    from predictor.effector import inducer_vocab as V

    ind = inducers or {}
    stored = ind.get("top_display")
    if isinstance(stored, str) and stored.strip():
        return stored
    return V.display_label(ind.get("top"), ind.get("candidates") or ())


def classify(d) -> tuple:
    """(outcome, {axis: spoke?}) for one bundle's dossier dict."""
    from predictor.effector import inducer_vocab as V

    ind = d.get("inducers") or {}
    calls = {c["source"]: c for c in ind.get("calls", [])}
    gate_pos = ind.get("coordination_gate") is True
    # gate SILENT (None) is not the same as gate negative: the family has no chemistry implemented,
    # so the axis never spoke. Panel B counts speaking, not agreeing.
    gate_spoke = ind.get("coordination_gate") is not None
    mn = calls.get("metalnet") or {}
    site = bool((mn.get("evidence") or {}).get("has_site"))
    clade = bool((calls.get("ssn_cluster") or {}).get("ligand"))
    folded = bool((d.get("structure") or {}).get("folded"))

    cls = V.inducer_class_of(ind.get("top") or "")
    if cls == "metal":
        k = sum([gate_pos, site, clade])
        out = ("metal_corroborated" if k >= 2 else
               "metal_single" if k == 1 else "metal_context")
    elif cls in ("redox", "organic"):
        out = cls
    else:
        out = "undetermined"
    return out, {"gate": gate_spoke, "metalnet": bool(mn), "clade": clade, "fold": folded}
