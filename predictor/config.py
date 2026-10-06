"""
config.py -- RunConfig: what this environment can reach, not which algorithm to run.

The pipeline has ONE path. It does not branch on configuration: every TF of every family goes through the
same discovery, motif, regulon and inducer steps. What can legitimately differ between machines is what
the run can *reach* -- is the network up, is a folding-backend token present, is ColabFold available -- and
that is all this object carries. There is no production/validation switch, no choice of homolog source and
no sampling strategy: a run is agnostic of pre-gathered local data by construction, and the homolog route
is the TF's SSN cluster with a BLAST fallback for families that have no SSN.

`allow_online` and `allow_folds` default to on; a caller with no network or no fold backend turns them off
(`predictor.api._resolve_folds` does this automatically from the token). `allow_colabfold` follows
`allow_online` unless set.

Every field here answers "what may this run SPEND, and what can it REACH" -- the expensive optional
engines (ESMFold2, ColabFold, MetalNet2) and whether the AF3 hand-off is written. That is why they are
legitimate while a run mode is not: turning MetalNet off does not choose a different algorithm, it
produces exactly what a machine without the engine installed produces, and the source abstains in the
open. Do not add a field here that selects a scoring rule, a family branch or a homolog route.

Run `python config.py` for a self-test.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict


@dataclass
class RunConfig:
    # reachability of the outside world (NCBI, AFDB, ColabFold, UniProt, AlphaFill)
    allow_online: bool = True
    # structure step: run the early ESMFold2 apo-dimer fold (needs a backend token)
    allow_folds: bool = True
    # deep MSA via ColabFold (None -> follows allow_online)
    allow_colabfold: bool | None = None
    # MetalNet2 metal-site gate: ~2 min/sequence (ESM-2 650M + an AutoGluon bag load per process).
    # Off -> the source abstains and the coordination gate alone speaks, exactly as on a machine where
    # the engine is not installed. A COST lever, not an algorithm switch: the path does not change.
    allow_metalnet: bool = True
    # Write the AlphaFold-3 job JSONs for Stage-B structural validation. Off when a run is not going
    # to be folded, so the bundle does not carry jobs nobody will submit.
    emit_af3_jobs: bool = True
    # bookkeeping
    write: bool = True
    verbose: bool = True

    def __post_init__(self):
        if self.allow_colabfold is None:
            self.allow_colabfold = self.allow_online

    def as_dict(self) -> dict:
        return asdict(self)


def _demo() -> None:
    c = RunConfig()
    assert c.allow_online is True and c.allow_folds is True and c.allow_colabfold is True
    # offline unit-test style: same logic, no network -- ColabFold follows the network by default
    off = RunConfig(allow_online=False)
    assert off.allow_online is False and off.allow_colabfold is False
    # an explicit ColabFold choice is not overridden by the network default
    both = RunConfig(allow_online=False, allow_colabfold=True)
    assert both.allow_colabfold is True
    assert RunConfig(allow_folds=False).allow_folds is False
    assert set(RunConfig().as_dict()) == {"allow_online", "allow_folds", "allow_colabfold",
                                          "allow_metalnet", "emit_af3_jobs", "write", "verbose"}
    # every lever is about COST or REACHABILITY -- what the run may spend and what it can reach.
    # None of them selects an algorithm, a family branch or a scoring rule, and adding one that did
    # would break the single-path rule this object exists to protect.
    assert RunConfig(allow_metalnet=False).allow_metalnet is False
    assert RunConfig(emit_af3_jobs=False).emit_af3_jobs is False
    print("OK: RunConfig carries environment reachability and cost only -- no algorithm switches.")


if __name__ == "__main__":
    _demo()
