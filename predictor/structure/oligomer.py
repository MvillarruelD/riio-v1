"""
oligomer.py -- the biologically relevant oligomeric state to fold a TF at.

Why this matters: bacterial one-component TFs are obligate **homodimers** (winged-HTH ArsR/SmtB, MerR,
TetR, Fur, HTH regulators), and several families are higher-order -- the **CsoR/RcnR** copper/nickel
sensors are **homotetramers**, **LysR-type (LTTR)** and LacI/GalR assemble as **tetramers**. Folding a
single chain throws away every effector/DNA contact that forms ACROSS the subunit interface (e.g. the
ArsR alpha5 inter-subunit metal site, the MerR dimer reading both operator half-sites). So the folding
backend must be told how many copies to assemble; AlphaFold DB (monomer-only) cannot supply this.

`oligomeric_state(family)` returns the copy number; unknown families default to a homodimer (2), the
overwhelmingly common case, and the call is flagged so it can be overridden from curated annotation.
"""
from __future__ import annotations

# family (as classified by annotate/family_db) -> homo-oligomer copy number. Keys matched case-folded
# and by substring, so "CsoR/RcnR", "CsoR-like", "LysR-type (LTTR)" all resolve.
_STATE = {
    "csor": 4, "rcnr": 4,                 # Cu(I)/Ni(II) sensors -- homotetramer (disc tetramer)
    "lysr": 4, "lttr": 4,                 # LysR-type transcriptional regulators -- tetramer
    "laci": 4, "galr": 4,                 # LacI/GalR -- tetramer (dimer of dimers)
    "arsr": 2, "smtb": 2,                 # winged-helix metal/RSS sensors -- homodimer
    "merr": 2, "fur": 2, "tetr": 2,       # classic homodimers
    "marr": 2, "padr": 2, "crp": 2, "fnr": 2, "gntr": 2, "arac": 2, "xyls": 2, "asnc": 2,
}
DEFAULT_STATE = 2                          # homodimer -- the default for an unrecognised family

_NAME = {1: "monomer", 2: "homodimer", 3: "homotrimer", 4: "homotetramer", 6: "homohexamer"}


def oligomeric_state(family: str | None) -> int:
    """Copy number to fold at for `family` (family_db label). Defaults to 2 (homodimer)."""
    if not family:
        return DEFAULT_STATE
    f = family.lower()
    for key, n in _STATE.items():
        if key in f:
            return n
    return DEFAULT_STATE


def is_default(family: str | None) -> bool:
    """True when the state is the fall-back homodimer (i.e. family not explicitly known)."""
    if not family:
        return True
    f = family.lower()
    return not any(key in f for key in _STATE)


def state_name(n: int) -> str:
    return _NAME.get(n, f"{n}-mer")


def _demo() -> None:
    assert oligomeric_state("CsoR/RcnR") == 4, "CsoR must be a tetramer"
    assert oligomeric_state("LysR-type (LTTR)") == 4, "LysR must be a tetramer"
    assert oligomeric_state("ArsR/SmtB") == 2 and oligomeric_state("MerR") == 2, "archetypes are dimers"
    assert oligomeric_state("TetR/AcrR") == 2, "TetR is a dimer"
    assert oligomeric_state("SomeNovelFamily") == 2 and is_default("SomeNovelFamily"), \
        "unknown family must default to homodimer and flag it"
    assert not is_default("CsoR/RcnR"), "known family should not be flagged default"
    assert oligomeric_state(None) == 2, "missing family -> dimer"
    print("oligomeric states: CsoR=4, LysR=4, ArsR=2, MerR=2, TetR=2, unknown->2(default-flagged)")
    print("OK: family-aware oligomeric state (dimer default, CsoR/LysR tetramer).")


if __name__ == "__main__":
    _demo()
