"""
api.py -- the stable programmatic entry point for the TF -> operator pipeline.

This is the ONE contract that the CLI (`predictor.run_cli` / the `tfop` command) and the GUI both call, so
internal refactors of the pipeline never break callers. It wraps `predictor.pipeline.run_novel`, resolves
the folding policy, and returns a small `Result` object (convenient attributes + the full raw dict).

    from predictor.api import predict
    res = predict("mytf.fasta", name="MyTF")
    if res.ok:
        print(res.primary_operator, res.bundle, res.report_html)
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class Result:
    """A friendly view over a pipeline run. `raw` holds the complete run_novel result dict."""
    ok: bool
    tf_id: str
    error: str | None = None
    family: str | None = None
    family_support: str | None = None
    genome: str | None = None
    bundle: Path | None = None
    primary_operator: str | None = None
    primary_source: str | None = None
    motif_consensus: str | None = None
    inducer: str | None = None
    raw: dict = field(default_factory=dict)

    @property
    def report_html(self) -> Path | None:
        """Path to the human-readable REPORT.html (None if the run did not write a bundle)."""
        return (self.bundle / "REPORT.html") if self.bundle else None

    @classmethod
    def from_run(cls, res: dict) -> "Result":
        prim = res.get("operators_primary") or {}
        bundle = res.get("bundle")
        return cls(
            ok=bool(res.get("required_stages_complete", not res.get("error"))) and not res.get("error"),
            tf_id=res.get("tf_id", "query"),
            error=res.get("error"),
            family=res.get("family"),
            family_support=res.get("family_support"),
            genome=res.get("genome"),
            bundle=(Path(bundle) if bundle else None),
            primary_operator=prim.get("sequence"),
            primary_source=prim.get("source"),
            motif_consensus=res.get("motif_consensus"),
            inducer=res.get("inducer"),
            raw=res,
        )


def _resolve_folds(cfg, fold):
    """Folding policy. `fold=True` opts in; anything else leaves the structural branch OFF.

    The default used to be "fold if a backend token happens to exist", which made the contents of a
    bundle depend on the machine and on run history rather than on the input: of the last run's 140
    candidates, 95 folded and 45 did not, and successes were cached while failures were not. Whether
    a protein has a structure is not supposed to be a property of when you ran it.

    Nothing downstream reads the fold. It produces report figures and a `pending_af3` placeholder in
    the operator candidate list; the family call, the inducer call, the operators and the regulon are
    identical with it on or off. So OFF is the honest default and ON is an enrichment.
    """
    cfg.allow_folds = fold is True
    return cfg


def predict(seq: str, *, name: str = "query", family: str | None = None, organism: str | None = None,
            genome_acc: str | None = None, effector: str | None = None, fold: bool | None = None,
            cfg=None, verbose: bool = True, **kw) -> Result:
    """Run the novel-TF pipeline on `seq` (a protein sequence or a FASTA path).

    fold: True runs the optional structural branch (ESMFold2 apo oligomer + AlphaFold-DB QC).
        The default is OFF; it enriches the report and changes no prediction. A Biohub token is
        required. Pass a ready `cfg` (RunConfig) to declare what this environment can reach; otherwise
        a default one is built. Returns a `Result`.
    """
    if "_seed_construction" in kw:
        raise TypeError("seed construction overrides are not accepted by the public API")
    from predictor.config import RunConfig
    from predictor import pipeline
    if cfg is None:
        cfg = RunConfig()
    cfg = _resolve_folds(cfg, fold)
    cfg.verbose = verbose
    res = pipeline.run_novel(seq, name=name, family=family, organism=organism, genome_acc=genome_acc,
                             effector=effector, cfg=cfg, **kw)
    return Result.from_run(res)
