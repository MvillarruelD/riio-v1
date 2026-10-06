"""Pure, replayable construction and selection of motif-width genome-scan trials.

The production orchestrator remains :func:`predictor.pipeline.run_novel`: it generates candidates and
runs the genome scans.  This module only prepares the seed PWMs, turns already-computed scan hits into
serializable trial records, and chooses among those records without access to family labels, effectors,
SSN annotations, curated sites, or benchmark coordinates.
"""
from __future__ import annotations

import copy
from dataclasses import asdict, dataclass
import re

import numpy as np

from predictor.signals import operator_logo
from predictor.signals import coords


@dataclass(frozen=True)
class SeedConstruction:
    """Private construction inputs for production and benchmark source bundles.

    The tuning harness is the only caller expected to override these values.  They deliberately do not
    live in ``RunConfig`` and are not exposed by the public API or CLI.
    """

    #: PROMOTED 2026-08-19 (was 24). Raising top_k to 24 does buy seed-pool site coverage
    #: (39 % -> 63 % overall, 44 % -> 73 % on the metal families) but that coverage is not convertible
    #: into a better shipped top pick by ANY of the eight truth-blind selectors tested: at top_k=24 the
    #: modal length shifts, a different candidate builds the PWM, and IscR and MntR are LOST while only
    #: ArsR is gained. Held-out outcome is identical (3/10 non-anchored under both LOTF and LOFO).
    #: Full working: the 2026-08-19 seed decision, recorded in the analysis repo before commit 5815559.
    top_k: int = 8
    tol: int = 3
    pseudocount_per_base: float = 0.25
    #: How `motif_finder` scores a dyad. `"raw"` is production's `2m - h`, whose value for a perfect
    #: dyad is just the half-site length -- see `motif_finder.SCORE_MODES`. The alternatives exist to
    #: be MEASURED, not switched on.
    #:
    #: **Deliberately not read from the environment.** This selects a SCORING RULE, and an ambient
    #: variable that does so is exactly what `config.RunConfig` exists to prevent -- "every lever is
    #: about COST or REACHABILITY ... none of them selects an algorithm, a family branch or a scoring
    #: rule". An env-var default would also be bound once at import, so a run could be scored by a
    #: different rule than the one its operator was designed under, with the shell as the only
    #: record. A benchmark that wants another mode passes it here explicitly, where `to_dict()`
    #: carries it into the construction provenance and the bundle says which rule produced the seed.
    score_mode: str = "raw"
    source_run_id: str | None = None
    implementation_commit: str | None = None
    implementation_dirty_diff_sha256: str | None = None
    generator_source_sha256: str | None = None
    source_runner_sha256: str | None = None
    fixed_family_inputs_sha256: str | None = None
    mirror_genome_id: str | None = None
    expected_sequence_accession: str | None = None
    mirror_genome_sha256: str | None = None
    panel_fasta_sha256: str | None = None
    construction_grid_sha256: str | None = None
    generation_context_sha256: str | None = None

    def __post_init__(self):
        if isinstance(self.top_k, bool) or not isinstance(self.top_k, int) or self.top_k < 1:
            raise ValueError("top_k must be a positive integer")
        if isinstance(self.tol, bool) or not isinstance(self.tol, int) or self.tol < 0:
            raise ValueError("tol must be a non-negative integer")
        if self.pseudocount_per_base <= 0:
            raise ValueError("pseudocount_per_base must be positive")
        from predictor.signals.motif_finder import SCORE_MODES
        if self.score_mode not in SCORE_MODES:
            raise ValueError(f"score_mode must be one of {SCORE_MODES}, got {self.score_mode!r}")
        provenance = (
            self.source_run_id,
            self.implementation_commit,
            self.implementation_dirty_diff_sha256,
            self.generator_source_sha256,
            self.source_runner_sha256,
            self.fixed_family_inputs_sha256,
            self.mirror_genome_id,
            self.expected_sequence_accession,
            self.mirror_genome_sha256,
            self.panel_fasta_sha256,
            self.construction_grid_sha256,
            self.generation_context_sha256,
        )
        if any(value is not None for value in provenance) and not all(provenance):
            raise ValueError("source-run provenance must be supplied completely or omitted")
        if self.source_run_id is not None:
            if not re.fullmatch(r"[0-9a-f]{40,64}", self.implementation_commit or ""):
                raise ValueError("implementation_commit must be a lowercase git object id")
            for field, value in (
                ("implementation_dirty_diff_sha256", self.implementation_dirty_diff_sha256),
                ("generator_source_sha256", self.generator_source_sha256),
                ("source_runner_sha256", self.source_runner_sha256),
                ("fixed_family_inputs_sha256", self.fixed_family_inputs_sha256),
                ("mirror_genome_sha256", self.mirror_genome_sha256),
                ("panel_fasta_sha256", self.panel_fasta_sha256),
                ("construction_grid_sha256", self.construction_grid_sha256),
                ("generation_context_sha256", self.generation_context_sha256),
            ):
                if not re.fullmatch(r"[0-9a-f]{64}", value or ""):
                    raise ValueError(f"{field} must be a lowercase SHA-256")
            if not self.mirror_genome_id or not self.expected_sequence_accession:
                raise ValueError("source-run genome identities must be non-empty")

    def to_dict(self) -> dict:
        """Serialize with the exact field names consumed by ``seed_policy_sweep.py``."""
        return asdict(self)


DEFAULT_CONSTRUCTION = SeedConstruction()


@dataclass(frozen=True)
class SeedPolicy:
    """A named seed-width policy with construction and replay provenance."""

    #: PROMOTED 2026-08-19. `cluster0` -- take the trial built from the largest-support anchor width --
    #: was selected by the held-out harness in 139 of 139 applications (80 leave-one-TF-out, 59
    #: leave-one-family-out). The run-6 `mean_ic` selector it replaces was chosen in NONE of them, and is
    #: measured at 1/17 complete-panel primary-site recovery against this policy's 7/17. Mean IC is a
    #: SELF-FIT score -- a cluster's PWM creates the very hits its mean IC is computed from -- so sparse
    #: matrices manufacture tidy hit sets and win on a criterion that measures nothing external.
    #: The historical policies below are retained for the benchmark harness, not for production.
    policy_name: str = "cluster0_promoted_20260819"
    selection_method: str = "cluster0"
    support_basis: str = "input"
    min_support: int = 1
    min_hits: int = 3
    max_hits: int = 40
    keep_frac: float = 0.75
    trim_ic: float = 0.35
    pseudocount_per_base: float = 0.25
    audit_hit_limit: int = 200
    hit_ranking: str = "raw_pwm_score"
    fallback: str = "cluster_0"
    rescan_scope: str = "intergenic"
    rescan_pvalue_thresh: float = 1e-4

    def __post_init__(self):
        if not self.policy_name:
            raise ValueError("policy_name must be non-empty")
        if self.support_basis not in {"input", "effective", "independent"}:
            raise ValueError(f"unsupported support_basis: {self.support_basis!r}")
        if self.selection_method not in {"mean_ic", "cluster0"}:
            raise ValueError(f"unsupported selection_method: {self.selection_method!r}")
        if self.fallback not in {"all_trials_mean_ic", "cluster_0"}:
            raise ValueError(f"unsupported fallback: {self.fallback!r}")
        if self.hit_ranking != "raw_pwm_score":
            raise ValueError(f"unsupported hit_ranking: {self.hit_ranking!r}")
        if self.rescan_scope not in {"intergenic", "genome"}:
            raise ValueError(f"unsupported rescan_scope: {self.rescan_scope!r}")
        if self.min_support < 0 or self.min_hits < 1 or self.max_hits < self.min_hits:
            raise ValueError("support must be non-negative; hit minima positive; max_hits >= min_hits")
        if self.pseudocount_per_base <= 0 or self.audit_hit_limit < 0:
            raise ValueError("pseudocount_per_base must be positive and audit_hit_limit non-negative")
        if not 0 < self.keep_frac <= 1:
            raise ValueError("keep_frac must be in (0, 1]")
        if not 0 <= self.trim_ic <= 2:
            raise ValueError("trim_ic must be in [0, 2]")
        if not 0 < self.rescan_pvalue_thresh <= 1:
            raise ValueError("rescan_pvalue_thresh must be in (0, 1]")

    def parameters(self) -> dict:
        out = asdict(self)
        out.pop("policy_name")
        return out


#: The one policy production runs. Everything below it is history kept runnable for the benchmark.
DEFAULT_POLICY = SeedPolicy()
#: PRODUCTION USES EXACTLY ONE POLICY: `DEFAULT_POLICY` above. The four below are not a second code
#: path -- nothing in `predictor/` selects them, `choose_trial` merely takes a policy argument, and
#: their only callers are the regression tests that pin why the promoted selector was chosen over
#: each of them. They are kept for that coverage and for no other reason; the sweep that generated
#: them is retired.
V4_MEAN_IC_CONTROL = SeedPolicy(
    policy_name="run6_v4_mean_ic_control",
    selection_method="mean_ic",
    fallback="all_trials_mean_ic",
)
V5_RAW_INPUT_SUPPORT_3_CONTROL = SeedPolicy(
    policy_name="raw_input_support_3_v5_control",
    selection_method="mean_ic",
    support_basis="input",
    min_support=3,
    fallback="cluster_0",
)
EFFECTIVE_SUPPORT_3_PROPOSAL = SeedPolicy(
    policy_name="effective_support_3_proposal",
    selection_method="mean_ic",
    support_basis="effective",
    min_support=3,
    fallback="cluster_0",
)
#: Historical alias: this IS the promoted policy, under the name the 2026-08-14 replay recorded.
CLUSTER0_CONTROL = SeedPolicy(
    policy_name="cluster0_pre_run6_control",
    selection_method="cluster0",
    fallback="cluster_0",
)

_REPLAY_SELECTION_FIELDS = {
    "policy_name", "selection_method", "support_basis", "min_support", "fallback",
}


def _construction_parameters(policy: SeedPolicy) -> dict:
    return {
        key: value
        for key, value in asdict(policy).items()
        if key not in _REPLAY_SELECTION_FIELDS
    }


@dataclass(frozen=True)
class SeedCluster:
    """A prepared seed cluster and its non-serializable PWM runtime value."""

    cluster_id: int
    anchor_w: int
    members: tuple
    member_metadata: tuple[dict, ...]
    prepared: operator_logo.PreparedSequences
    n_independent: int
    independence_rule: dict
    build_policy: SeedPolicy
    counts: np.ndarray
    pwm: np.ndarray

    @property
    def n_input(self) -> int:
        return len(self.members)

    @property
    def n_effective(self) -> int:
        return self.prepared.n_effective


@dataclass(frozen=True)
class SeedTrial:
    """One cluster's complete replay record plus its in-memory downstream values."""

    cluster: SeedCluster
    hits: tuple
    audit: dict

    @property
    def cluster_id(self) -> int:
        return self.cluster.cluster_id

    @property
    def anchor_w(self) -> int:
        return self.cluster.anchor_w

    @property
    def n_effective(self) -> int:
        return self.cluster.n_effective

    @property
    def n_input(self) -> int:
        return self.cluster.n_input

    @property
    def n_independent(self) -> int:
        return self.cluster.n_independent

    @property
    def n_hits(self) -> int:
        return int(self.audit["n_hits"])

    @property
    def mean_ic(self) -> float:
        return float(self.audit["mean_ic"])

    def to_dict(self) -> dict:
        return copy.deepcopy(self.audit)


@dataclass(frozen=True)
class ReplayTrial:
    """A selection-only trial reconstructed from a JSON audit record."""

    audit: dict
    build_policy: SeedPolicy

    @property
    def cluster_id(self) -> int:
        return int(self.audit["cluster_id"])

    @property
    def anchor_w(self) -> int:
        return int(self.audit["anchor_w"])

    @property
    def n_effective(self) -> int:
        return int(self.audit["n_effective"])

    @property
    def n_input(self) -> int:
        return int(self.audit["n_input"])

    @property
    def n_independent(self) -> int:
        return int(self.audit["n_independent"])

    @property
    def n_hits(self) -> int:
        return int(self.audit["n_hits"])

    @property
    def mean_ic(self) -> float:
        return float(self.audit["mean_ic"])

    def to_dict(self) -> dict:
        return copy.deepcopy(self.audit)


@dataclass(frozen=True)
class SelectionDecision:
    """Chosen runtime trial and the serializable dossier payload explaining the choice."""

    chosen: SeedTrial | ReplayTrial
    policy: SeedPolicy
    eligible_cluster_ids: tuple[int, ...]
    fallback_reason: str | None

    def to_dict(self, trials) -> dict:
        if self.policy.selection_method == "cluster0":
            reason = "explicit cluster-0 control"
        else:
            reason = (
                self.fallback_reason
                or f"best mean IC among {len(self.eligible_cluster_ids)} eligible cluster(s)"
            )
        return {
            "policy_name": self.policy.policy_name,
            "parameters": self.policy.parameters(),
            "chosen_cluster_id": self.chosen.cluster_id,
            "chosen_w": self.chosen.anchor_w,
            "eligible_cluster_ids": list(self.eligible_cluster_ids),
            "fallback_used": self.fallback_reason is not None,
            "fallback_reason": self.fallback_reason,
            # Back-compatible human-readable alias retained for existing dossier consumers.
            "chosen_why": reason,
            "trials": [trial.to_dict() for trial in trials],
        }


def _candidate_seq(candidate) -> str:
    return str(getattr(candidate, "seq", "") or "").upper()


def _candidate_span(candidate) -> tuple[int, int]:
    start = int(getattr(candidate, "start", 0) or 0)
    end = int(getattr(candidate, "end", start + len(_candidate_seq(candidate))) or start)
    return start, end


def _candidate_center(candidate) -> float:
    center = getattr(candidate, "center", None)
    if center is not None:
        return float(center)
    start, end = _candidate_span(candidate)
    return (start + end) / 2.0


def _reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def _identity(a: str, b: str) -> float:
    if not a or not b:
        return 0.0
    denom = max(len(a), len(b))
    direct = sum(x == y for x, y in zip(a, b)) / denom
    reverse = sum(x == y for x, y in zip(a, _reverse_complement(b))) / denom
    return max(direct, reverse)


def _overlap_fraction(a, b) -> float:
    a0, a1 = _candidate_span(a)
    b0, b1 = _candidate_span(b)
    overlap = max(0, min(a1, b1) - max(a0, b0))
    return overlap / max(1, min(a1 - a0, b1 - b0))


def independent_member_indices(
    members,
    *,
    center_min_sep: float = 4.0,
    overlap_fraction: float = 0.8,
    identity_fraction: float = 0.9,
) -> tuple[int, ...]:
    """Return deterministic non-redundant member indices in the supplied (ranked) order.

    The first rule mirrors motif_finder's existing center deduplication.  The second explicitly merges
    near-identical, substantially overlapping calls that survived with centers just outside that radius.
    """
    kept: list[int] = []
    for index, candidate in enumerate(members):
        redundant = False
        for prior_index in kept:
            prior = members[prior_index]
            centers_duplicate = abs(_candidate_center(candidate) - _candidate_center(prior)) <= center_min_sep
            overlap_duplicate = (
                _overlap_fraction(candidate, prior) >= overlap_fraction
                and _identity(_candidate_seq(candidate), _candidate_seq(prior)) >= identity_fraction
            )
            if centers_duplicate or overlap_duplicate:
                redundant = True
                break
        if not redundant:
            kept.append(index)
    return tuple(kept)


def _summary(values) -> dict:
    vals = [float(value) for value in values]
    if not vals:
        return {"n": 0, "min": None, "max": None, "mean": None}
    return {
        "n": len(vals),
        "min": min(vals),
        "max": max(vals),
        "mean": float(sum(vals) / len(vals)),
    }


def _eligible_candidate_pool(candidates) -> list:
    """Reproduce motif_finder.seed_clusters' exact eligible pool and deterministic order."""
    candidates = [candidate for candidate in candidates if _candidate_seq(candidate)]
    dyads = [candidate for candidate in candidates if getattr(candidate, "kind", None) in ("IR", "DR")]
    base = dyads or candidates
    conserved = [candidate for candidate in base if getattr(candidate, "conservation", 0.0) > 0]
    return sorted(
        conserved or base,
        key=lambda candidate: -getattr(candidate, "final", getattr(candidate, "score", 0.0)),
    )


def _match_members(candidates, clustered_sequences) -> list[tuple]:
    """Map seed_clusters' sequence-only result back to its candidate objects, occurrence by occurrence."""
    available: dict[str, list] = {}
    for candidate in _eligible_candidate_pool(candidates):
        available.setdefault(_candidate_seq(candidate), []).append(candidate)
    matched = []
    for sequences, _anchor_w in clustered_sequences:
        members = []
        for sequence in sequences:
            key = str(sequence or "").upper()
            if not available.get(key):
                raise ValueError(f"cluster sequence has no matching candidate: {key!r}")
            members.append(available[key].pop(0))
        matched.append(tuple(members))
    return matched


def prepare_seed_clusters(candidates, clustered_sequences, *, policy: SeedPolicy = DEFAULT_POLICY):
    """Build deterministic cluster PWMs and auditable raw/effective/independent support metadata."""
    clustered_sequences = list(clustered_sequences)
    matched = _match_members(candidates, clustered_sequences)
    out = []
    for cluster_id, ((sequences, anchor_w), members) in enumerate(zip(clustered_sequences, matched)):
        prepared = operator_logo.prepare_sequences(sequences)
        counts = operator_logo.counts_from_prepared(prepared)
        if counts is None or prepared.modal_w is None:
            continue
        contributing = tuple(
            member for member in members if len(_candidate_seq(member)) == prepared.modal_w
        )
        independent_local = independent_member_indices(contributing)
        independent_ids = {id(contributing[index]) for index in independent_local}
        member_metadata = []
        for member_index, member in enumerate(members):
            start, end = _candidate_span(member)
            contributed = len(_candidate_seq(member)) == prepared.modal_w
            member_metadata.append({
                "member_index": member_index,
                "sequence": _candidate_seq(member),
                "width": len(_candidate_seq(member)),
                "start": start,
                "end": end,
                "center": _candidate_center(member),
                "type": str(getattr(member, "kind", "") or ""),
                "score": float(getattr(member, "score", 0.0) or 0.0),
                "final": float(getattr(member, "final", getattr(member, "score", 0.0)) or 0.0),
                "conservation": float(getattr(member, "conservation", 0.0) or 0.0),
                "contributed_to_pwm": contributed,
                "independent_contributor": contributed and id(member) in independent_ids,
            })
        prior = policy.pseudocount_per_base
        pwm = (counts + prior) / (counts + prior).sum(axis=0, keepdims=True)
        out.append(SeedCluster(
            cluster_id=cluster_id,
            anchor_w=int(anchor_w),
            members=members,
            member_metadata=tuple(member_metadata),
            prepared=prepared,
            n_independent=len(independent_local),
            independence_rule={
                "center_min_sep": 4.0,
                "overlap_fraction": 0.8,
                "identity_fraction": 0.9,
                "orientation": "forward_or_reverse_complement",
            },
            build_policy=policy,
            counts=counts,
            pwm=pwm,
        ))
    return out


def _hit_metadata(hit, rank: int) -> dict:
    generator = str(getattr(hit, "generator", "") or "")
    tier = getattr(hit, "tier", None)
    if tier is None and generator.startswith("rescan_tier"):
        suffix = generator.removeprefix("rescan_tier")
        tier = int(suffix) if suffix.isdigit() else None
    combined = getattr(hit, "combined_score", None)
    if combined is None:
        combined = getattr(hit, "combined", None)
    if combined is None:
        combined = getattr(hit, "sequence_score", None)
    return {
        "rank": rank,
        "accession": getattr(hit, "genome_accession", getattr(hit, "accession", None)),
        "start": int(getattr(hit, "start", 0)),
        "end": int(getattr(hit, "end", 0)),
        "strand": str(getattr(hit, "strand", "+")),
        "dyad": int(getattr(hit, "dyad_center", 0)),
        "sequence": str(getattr(hit, "seq", getattr(hit, "matched", "")) or ""),
        "score": _optional_float(getattr(hit, "score", None)),
        "pvalue": _optional_float(getattr(hit, "pvalue", None)),
        "qvalue": _optional_float(getattr(hit, "qvalue", None)),
        "tier": tier,
        "combined_score": _optional_float(combined),
    }


def _optional_float(value):
    return None if value is None else float(value)


def _normalize_hit(hit):
    """Return a canonical, motif-oriented hit; adapt native motif_rescan GenomeHit values exactly once."""
    if getattr(hit, "matched", "") and not getattr(hit, "seq", ""):
        return coords.from_genome_hit(hit)
    if hasattr(hit, "seq"):
        return hit
    raise ValueError("hit must be a canonical operator or motif_rescan GenomeHit")


def build_trial(cluster: SeedCluster, hits, *, policy: SeedPolicy | None = None) -> SeedTrial:
    """Turn one prepared cluster and its already-computed genome hits into a replayable trial."""
    policy = cluster.build_policy if policy is None else policy
    if policy != cluster.build_policy:
        raise ValueError("trial policy does not match the policy/prior that prepared this seed cluster")
    original_hits = tuple(_normalize_hit(hit) for hit in hits)
    ranked = tuple(sorted(original_hits, key=lambda hit: -(getattr(hit, "score", None) or 0.0)))
    ranked_sequences = [str(getattr(hit, "seq", "") or "") for hit in ranked if getattr(hit, "seq", "")]
    selected = operator_logo.select_by_motif(
        ranked_sequences,
        min_hits=policy.min_hits,
        max_hits=policy.max_hits,
        keep_frac=policy.keep_frac,
        trim_ic=policy.trim_ic,
    ) if ranked_sequences else operator_logo.MotifSelection()
    logo = selected.logo
    mean_ic = (
        logo.total_ic / logo.pwm.shape[1]
        if logo.pwm is not None and logo.pwm.shape[1]
        else 0.0
    )
    limit = max(0, int(policy.audit_hit_limit))
    recorded_hits = [_hit_metadata(hit, rank) for rank, hit in enumerate(ranked[:limit], start=1)]
    finals = [metadata["final"] for metadata in cluster.member_metadata]
    conservation = [metadata["conservation"] for metadata in cluster.member_metadata]
    audit = {
        "build_policy": {
            "policy_name": policy.policy_name,
            "parameters": policy.parameters(),
        },
        "cluster_id": cluster.cluster_id,
        "anchor_w": cluster.anchor_w,
        "member_widths": [metadata["width"] for metadata in cluster.member_metadata],
        "members": list(cluster.member_metadata),
        "n_input": cluster.n_input,
        "modal_w": cluster.prepared.modal_w,
        "n_effective": cluster.n_effective,
        "n_independent": cluster.n_independent,
        "independence_rule": copy.deepcopy(cluster.independence_rule),
        "pwm_prior": {
            "kind": "symmetric_dirichlet",
            "pseudocount_per_base": policy.pseudocount_per_base,
            "total_mass_per_column": 4.0 * policy.pseudocount_per_base,
        },
        "candidate_final_summary": _summary(finals),
        "candidate_conservation_summary": _summary(conservation),
        "selector": {
            "min_hits": policy.min_hits,
            "max_hits": policy.max_hits,
            "keep_frac": policy.keep_frac,
            "trim_ic": policy.trim_ic,
            "hit_ranking": policy.hit_ranking,
        },
        "n_hits": len(ranked),
        "n_keep": selected.n_keep,
        "peak_n": selected.peak_n,
        "ic_curve": list(selected.ic_curve),
        "mean_ic": round(float(mean_ic), 4),
        "consensus": logo.consensus,
        "ranked_hits": recorded_hits,
        "ranked_hits_limit": limit,
        "n_hits_recorded": len(recorded_hits),
        "hits_truncated": len(recorded_hits) < len(ranked),
    }
    # Preserve motif_rescan's combined-score ordering for every downstream consumer.  The separate
    # ``ranked_hits`` audit follows the raw-PWM ordering used by select_by_motif, exactly as run_novel did.
    return SeedTrial(cluster=cluster, hits=original_hits, audit=audit)


def trial_from_audit(audit) -> ReplayTrial:
    """Reconstruct a selection-only trial from one JSON-decoded audit record."""
    record = copy.deepcopy(dict(audit))
    required = {"build_policy", "cluster_id", "anchor_w", "n_input", "n_effective", "n_independent",
                "n_hits", "mean_ic", "ranked_hits", "ranked_hits_limit", "n_hits_recorded",
                "hits_truncated"}
    missing = sorted(required - set(record))
    if missing:
        raise ValueError(f"trial audit is missing required fields: {', '.join(missing)}")
    policy_payload = record["build_policy"]
    try:
        build_policy = SeedPolicy(
            policy_name=policy_payload["policy_name"],
            **policy_payload["parameters"],
        )
    except (KeyError, TypeError) as exc:
        raise ValueError("invalid build_policy payload in trial audit") from exc
    ranked_hits = record["ranked_hits"]
    if not isinstance(ranked_hits, list):
        raise ValueError("ranked_hits must be a list")
    if int(record["n_hits_recorded"]) != len(ranked_hits):
        raise ValueError("n_hits_recorded does not match ranked_hits")
    if len(ranked_hits) > int(record["ranked_hits_limit"]):
        raise ValueError("ranked_hits exceeds ranked_hits_limit")
    return ReplayTrial(audit=record, build_policy=build_policy)


def _trial_build_policy(trial) -> SeedPolicy:
    if isinstance(trial, ReplayTrial):
        return trial.build_policy
    return trial.cluster.build_policy


def _support(trial, basis: str) -> int:
    return int(getattr(trial, f"n_{basis}"))


def _selection_key(trial) -> tuple[float, int, int]:
    # seed_clusters emits canonical ascending IDs. Lower ID is therefore legacy stable-first on ties,
    # made explicit so JSON record order cannot alter the decision.
    return trial.mean_ic, trial.n_hits, -trial.cluster_id


def _validate_unique_cluster_ids(trials) -> None:
    cluster_ids = [trial.cluster_id for trial in trials]
    duplicates = sorted({cluster_id for cluster_id in cluster_ids if cluster_ids.count(cluster_id) > 1})
    if duplicates:
        duplicate_text = ", ".join(str(cluster_id) for cluster_id in duplicates)
        raise ValueError(f"trial cluster_id values must be unique; duplicate cluster_id(s): {duplicate_text}")


def _cluster_zero(trials, *, reason: str):
    cluster_zero = tuple(trial for trial in trials if trial.cluster_id == 0)
    if len(cluster_zero) != 1:
        raise ValueError(f"{reason} requires exactly one trial with cluster_id == 0")
    return cluster_zero[0]


def choose_trial(trials, *, policy: SeedPolicy = DEFAULT_POLICY) -> SelectionDecision:
    """Choose a trial using support and self-fit fields only."""
    trials = tuple(trials)
    if not trials:
        raise ValueError("at least one seed trial is required")
    _validate_unique_cluster_ids(trials)
    expected_construction = _construction_parameters(policy)
    for trial in trials:
        if _construction_parameters(_trial_build_policy(trial)) != expected_construction:
            raise ValueError("replay policy changes fields fixed during trial construction")
    if policy.selection_method == "cluster0":
        chosen = _cluster_zero(trials, reason="cluster0 selection")
        return SelectionDecision(
            chosen=chosen,
            policy=policy,
            eligible_cluster_ids=(chosen.cluster_id,),
            fallback_reason=None,
        )
    eligible = tuple(
        trial for trial in trials
        if _support(trial, policy.support_basis) >= policy.min_support and trial.n_hits >= policy.min_hits
    )
    if eligible:
        chosen = max(eligible, key=_selection_key)
        fallback_reason = None
    elif policy.fallback == "all_trials_mean_ic":
        # Exact committed run-6/v4 behaviour: prefer n_hits>=3, but if no cluster reaches it, compare
        # every trial by the same (mean_ic, n_hits) tuple rather than silently changing the selector.
        chosen = max(trials, key=_selection_key)
        fallback_reason = f"no cluster reached n_hits>={policy.min_hits}; compared all trials (v4 control)"
    else:
        chosen = _cluster_zero(trials, reason="cluster_0 fallback")
        fallback_reason = (
            f"no cluster reached n_{policy.support_basis}>={policy.min_support} and "
            f"n_hits>={policy.min_hits}; used cluster 0"
        )
    return SelectionDecision(
        chosen=chosen,
        policy=policy,
        eligible_cluster_ids=tuple(trial.cluster_id for trial in eligible),
        fallback_reason=fallback_reason,
    )


def choose_trial_records(records, *, policy: SeedPolicy | None = None) -> SelectionDecision:
    """Replay selection directly from JSON-decoded trial records, without genome/network access."""
    trials = tuple(trial_from_audit(record) for record in records)
    if not trials:
        raise ValueError("at least one trial audit record is required")
    _validate_unique_cluster_ids(trials)
    replay_policy = trials[0].build_policy if policy is None else policy
    return choose_trial(trials, policy=replay_policy)
