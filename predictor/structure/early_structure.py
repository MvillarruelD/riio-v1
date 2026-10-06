"""
early_structure.py -- the EARLY structure step of the redesigned pipeline (WS1).

The redesign folds the apo homodimer up front, for every TF, then cross-checks it against the free AFDB
monomer as a fold-integrity QC -- so the structure is available to the cluster/ligand steps and the
operator design, not bolted on at the end. This module is the quota-aware orchestrator:

  run_structure_step(seq, family, *, uniprot, name, cfg) -> StructureResult

Policy (graceful degradation, never hard-fails the pipeline):
  * fold the biologically selected apo homo-oligomer on biohub ESMFold2 ONLY when
    `cfg.allow_folds` AND a backend token is present;
    otherwise emit a *request descriptor* (seq, copies=oligomeric_state) and continue with AFDB-only QC.
  * AFDB-monomer QC is FREE (cached API): fetch the monomer for `uniprot`, superpose it onto the dimer's
    chain A (reuse `cluster_fold.superpose_rmsd`) -> RMSD + a pLDDT/agreement flag. When there is no dimer
    fold, the AFDB monomer alone still gives a monomer pLDDT proxy + serves the ligand step.
  * the resulting fold (cif) feeds `effector.inducer.infer_inducer(..., fold=...)` for the structural
    inducer signal -- wired by the caller (run_novel), not here.

Reuses `signals.structure_backend.fold_apo_homomer` + `.backend_available`, `structure.afdb`,
`structure/oligomer`, `structure/cluster_fold.superpose_rmsd`. Run `python -m predictor.structure.early_structure`
for a self-test (offline logic; the live apo fold + AFDB QC SKIP without a token / network).
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from pathlib import Path

_REPO = Path(__file__).resolve().parents[2]
from predictor.structure import afdb, oligomer                     # noqa: E402

# Fold-integrity QC is on the largest consistent RIGID CORE (see cluster_fold.core_superpose_rmsd), not a
# whole-chain RMSD -- the latter flags an expected inter-domain hinge (MerR coiled-coil) as a fold error.
# A fold passes when a substantial core superposes tightly onto the free monomer. Calibrated on the 37
# Phase-1 folds: real folds keep a 0.5-1.0 core at <2.5 A; genuine misfolds collapse to a <=0.16 core.
_QC_CORE_FRACTION = 0.50      # >= half the aligned residues in the rigid core
_QC_CORE_RMSD = 2.5           # core RMSD (A)
_LOW_PLDDT_FLAG = 70.0


@dataclass
class StructureResult:
    folded: bool                              # an apo dimer fold was produced this run (or cached)
    backend: str | None                      # "esmfold2" | None
    oligomeric_state: int = 2
    state_name: str = "homodimer"
    dimer_cif: str | None = None             # path to the apo dimer fold (chain A,B)
    plddt_mean: float | None = None          # dimer mean pLDDT (or AFDB monomer proxy)
    ptm: float | None = None
    afdb_monomer_cif: str | None = None      # free AFDB monomer used for QC / ligand
    afdb_entry: str | None = None
    qc_rmsd: float | None = None             # chain A vs AFDB monomer whole-chain Calpha RMSD (A, diagnostic)
    qc_n_matched: int | None = None
    qc_core_fraction: float | None = None    # fraction of aligned residues in the rigid core (robust QC)
    qc_core_rmsd: float | None = None        # RMSD over that core (A)
    qc_pass: bool | None = None              # None when QC could not run
    fold_seconds: float | None = None
    request: dict | None = None              # {seq, copies} when no fold was produced (descriptor)
    flags: list = field(default_factory=list)
    notes: list = field(default_factory=list)

    def as_dict(self) -> dict:
        return asdict(self)

    @property
    def cif(self) -> str | None:
        """The best available structure for the ligand step: the dimer fold, else the AFDB monomer."""
        return self.dimer_cif or self.afdb_monomer_cif


def run_structure_step(seq: str, family: str | None, *, uniprot: str | None = None,
                       name: str = "query", cfg=None, verbose: bool = True) -> StructureResult:
    """Fold the selected apo homo-oligomer + query AFDB-monomer QC. Never raises."""
    allow_folds = bool(getattr(cfg, "allow_folds", False))
    allow_online = bool(getattr(cfg, "allow_online", True))
    state = oligomer.oligomeric_state(family)
    res = StructureResult(folded=False, backend=None, oligomeric_state=state,
                          state_name=oligomer.state_name(state))

    # ---- 1. apo homo-oligomer fold (ESMFold2, quota-aware) ----------------------------------------
    from predictor.signals import structure_backend as sb
    # State-aware pre-check. The backend caps a COMPLEX at MAX_COMPLEX_RESIDUES, so the monomer
    # budget depends on how many copies this family folds as -- a 434 aa ArsR is fine alone and over
    # cap as a dimer. Checking here costs nothing; not checking it cost run A seven candidates at
    # ~30 min and one quota unit each, every one of them spent on a request the server was always
    # going to refuse with a 422. The failure is recorded exactly as a backend refusal would be, so
    # nothing downstream needs to distinguish them.
    _n_res = len(seq) * max(1, state)
    _over_cap = _n_res > sb.MAX_COMPLEX_RESIDUES
    if allow_folds and _over_cap:
        res.flags.append("apo_fold_failed")
        res.notes.append(
            f"apo {oligomer.state_name(state)} not folded: {len(seq)} aa x {state} copies = "
            f"{_n_res} residues, over the backend's {sb.MAX_COMPLEX_RESIDUES}-residue complex cap. "
            f"Not sent -- the request cannot succeed. Fold the DNA-binding domain alone to inspect "
            f"this one.")
    elif allow_folds and sb.backend_available():
        try:
            fr = sb.fold_apo_homomer(name, seq, copies=state)
            if fr and fr.get("cif"):
                res.folded = True
                res.backend = "esmfold2"
                res.dimer_cif = fr["cif"]
                res.plddt_mean = fr.get("plddt_mean")
                res.ptm = fr.get("ptm")
                res.fold_seconds = fr.get("fold_seconds")
                _cached = not fr.get("folded")
                _tstr = "cached" if _cached else f"{fr.get('fold_seconds')}s"
                res.notes.append(f"apo {oligomer.state_name(state)} folded (ESMFold2, {_tstr})")
        except Exception as e:
            res.flags.append("apo_fold_failed")
            res.notes.append(f"apo fold failed ({type(e).__name__}: {e}); continuing AFDB-only")
    if not res.folded:
        res.request = {"seq": seq, "copies": state}
        # Say which of the four it was. This used to be a two-way choice between "folds disabled"
        # and "no ESMFold2 backend token", reached by EVERY unfolded candidate -- including one
        # whose fold had just been attempted and failed a few lines above. Those runs emitted two
        # notes: the real error, then "no ESMFold2 backend token", which was false and contradicted
        # it. A reader looking at a failed fold was told to go check a token that was working.
        if not allow_folds:
            why = "folds disabled (cfg.allow_folds=False)"
        elif _over_cap:
            why = f"over the {sb.MAX_COMPLEX_RESIDUES}-residue complex cap"
        elif "apo_fold_failed" in res.flags:
            why = "the fold was attempted and failed -- see the preceding note for the error"
        else:
            why = "no ESMFold2 backend token"
        res.notes.append(f"apo {oligomer.state_name(state)} not folded ({why}); "
                         f"emitted a request descriptor")

    # ---- 2. AFDB-monomer integrity QC (FREE; cached API) ------------------------------------------
    # Same switch as the fold: AlphaFold-DB retrieval is part of the structural branch, not a
    # separate capability. Gating it on `allow_online` alone meant "structure off" still fetched a
    # model whenever the network was allowed.
    if uniprot and allow_online and allow_folds:
        try:
            meta = afdb.lookup(uniprot)
            mono = afdb.fetch_model(uniprot) if meta else None
            if mono and Path(mono).exists():
                res.afdb_monomer_cif = mono
                res.afdb_entry = (meta or {}).get("entry_id")
                if res.dimer_cif:
                    # Robust core metric, NOT whole-chain RMSD: chain A of a dimer and the free monomer can
                    # differ by an EXPECTED domain hinge (MerR's coiled-coil splays out to dimerise), which a
                    # rigid whole-chain overlay scores as a huge RMSD even for a perfectly folded chain. QC on
                    # the largest consistent rigid core instead -- a real misfold has no such core.
                    from predictor.structure.cluster_fold import core_superpose_rmsd
                    cm = core_superpose_rmsd(res.dimer_cif, mono)
                    if cm["global_rmsd"] != float("inf"):
                        res.qc_rmsd = cm["global_rmsd"]                 # keep global as the reported diagnostic
                        res.qc_n_matched = cm["n_aln"]
                        res.qc_core_fraction = cm["core_fraction"]
                        res.qc_core_rmsd = cm["core_rmsd"]
                        # PASS iff a substantial core superposes tightly (hinge artifact) rather than the
                        # whole chain (which fails for any hinged multi-domain fold). Thresholds calibrated
                        # on the 37 Phase-1 folds: good ArsR/MerR cores sit at fraction 0.5-1.0 / rmsd <2.5;
                        # genuine misfolds collapse to fraction <=0.16.
                        res.qc_pass = (cm["core_fraction"] >= _QC_CORE_FRACTION
                                       and cm["core_rmsd"] <= _QC_CORE_RMSD)
                        if not res.qc_pass:
                            res.flags.append("fold_integrity_low_core")
                        res.notes.append(
                            f"AFDB-monomer QC: chainA vs {res.afdb_entry} core "
                            f"{int(100*cm['core_fraction'])}% @ {cm['core_rmsd']} A "
                            f"(global {cm['global_rmsd']} A over {cm['n_aln']} CA) -> "
                            f"{'pass' if res.qc_pass else 'FLAG'}")
                else:
                    res.notes.append(f"AFDB monomer {res.afdb_entry} fetched (no dimer fold to QC against)")
        except Exception as e:
            res.notes.append(f"AFDB QC skipped ({type(e).__name__}: {e})")

    # ---- 3. confidence flags ----------------------------------------------------------------------
    if res.plddt_mean is not None and res.plddt_mean < _LOW_PLDDT_FLAG:
        res.flags.append("low_plddt")
    if verbose:
        st = ("FOLDED" if res.folded else "not folded")
        print(f"structure: {st} [{res.backend or '-'}] plddt={res.plddt_mean} "
              f"qc_rmsd={res.qc_rmsd} qc={res.qc_pass} flags={res.flags or '-'}")
    return res


def _demo() -> None:
    # offline: with no backend + no uniprot, the step degrades cleanly to a request descriptor, no raise
    class _Cfg:
        allow_folds = False
        allow_online = False
    r = run_structure_step("M" + "A" * 110, "ArsR/SmtB", uniprot=None, name="UnitTF", cfg=_Cfg(), verbose=False)
    assert r.folded is False and r.request == {"seq": "M" + "A" * 110, "copies": 2}, r
    assert r.oligomeric_state == 2 and r.cif is None
    assert r.qc_pass is None, "QC cannot run without a fold or AFDB monomer"
    d = r.as_dict()
    # the note names the state the way oligomer does ("homodimer"), so read it from there rather than
    # hard-coding a word that drifts when the state names change
    assert d["folded"] is False and f"apo {oligomer.state_name(2)} not folded" in " ".join(d["notes"]), d["notes"]
    print(f"OK (offline): structure step degrades to request descriptor; flags={r.flags} notes={len(r.notes)}")

    # live (SKIPs without network/token): AFDB QC at zero fold cost on a known ArsR (P0ACS5, E. coli ArsR)
    try:
        class _Cfg2:
            allow_folds = False
            allow_online = True
        r2 = run_structure_step("M" + "A" * 110, "ArsR/SmtB", uniprot="P0ACS5", name="ArsR_Ec",
                                cfg=_Cfg2(), verbose=False)
        if r2.afdb_monomer_cif:
            print(f"live: AFDB monomer {r2.afdb_entry} fetched (free); QC needs a dimer fold to compare.")
        else:
            print("SKIP (network): AFDB unreachable for P0ACS5.")
    except Exception as e:
        print(f"SKIP (network): AFDB QC not exercised ({type(e).__name__}: {e})")


if __name__ == "__main__":
    _demo()
