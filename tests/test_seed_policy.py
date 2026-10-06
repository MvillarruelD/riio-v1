"""Offline regression tests for replayable, family-agnostic seed-width trial selection."""
from __future__ import annotations

import inspect
import json
from types import SimpleNamespace

import numpy as np
import pytest

from predictor.signals import operator_logo
from predictor.signals import seed_selection as seed_policy


class Candidate:
    def __init__(self, sequence, *, start, score=1.0, kind="IR", conservation=0.0):
        self.seq = sequence
        self.start = start
        self.end = start + len(sequence)
        self.center = (self.start + self.end) / 2
        self.score = score
        self.final = score * (1.0 + conservation)
        self.kind = kind
        self.conservation = conservation


def _hit(sequence, score, *, start=100, tier=3):
    return SimpleNamespace(
        genome_accession="TEST_GENOME",
        start=start,
        end=start + len(sequence),
        strand="+",
        dyad_center=start + len(sequence) // 2,
        seq=sequence,
        score=score,
        pvalue=10 ** (-score),
        qvalue=min(1.0, 10 ** (-score + 1)),
        sequence_score=score - 0.2,
        combined_score=None,
        generator=f"rescan_tier{tier}",
    )


def _cluster(candidates, sequences, anchor_w, *, policy=seed_policy.DEFAULT_POLICY):
    return seed_policy.prepare_seed_clusters(
        candidates,
        [(sequences, anchor_w)],
        policy=policy,
    )[0]


def test_prepare_sequences_exposes_exact_contributors_without_changing_counts():
    sequences = ["AAAA", "CCCC", "GGGGG"]
    prepared = operator_logo.prepare_sequences(sequences, orient=False)
    counts = operator_logo.counts_from_seqs(sequences, orient=False)

    assert prepared.n_input == 3
    assert prepared.modal_w == 4
    assert prepared.n_effective == 2
    assert prepared.sequences == ("AAAA", "CCCC")
    expected = np.zeros((4, 4))
    expected[0, :] = 1
    expected[1, :] = 1
    np.testing.assert_array_equal(counts, expected)


def test_mixed_width_trial_records_raw_effective_and_independent_support():
    candidates = [
        Candidate("ACGT", start=0, score=3.0),
        Candidate("ACGA", start=2, score=2.0),
        Candidate("TTTTT", start=20, score=1.0),
    ]
    cluster = _cluster(candidates, [candidate.seq for candidate in candidates], 4)
    trial = seed_policy.build_trial(
        cluster,
        [_hit("ACGT", 6.0), _hit("ACGA", 5.0), _hit("ACGG", 4.0)],
    )

    assert trial.audit["n_input"] == 3
    assert trial.audit["modal_w"] == 4
    assert trial.audit["n_effective"] == 2
    # The two modal-width candidates overlap and have centers within the existing 4-bp dedupe radius.
    assert trial.audit["n_independent"] == 1
    assert [member["contributed_to_pwm"] for member in trial.audit["members"]] == [True, True, False]
    assert sum(member["independent_contributor"] for member in trial.audit["members"]) == 1
    np.testing.assert_array_equal(cluster.counts, operator_logo.counts_from_prepared(cluster.prepared))


def test_independence_is_deterministic_and_keeps_distinct_calls():
    candidates = [
        Candidate("ACGTAC", start=0, score=3.0),
        Candidate("ACGTAC", start=2, score=2.0),
        Candidate("TTGGCC", start=20, score=1.0),
    ]
    first = seed_policy.independent_member_indices(candidates)
    second = seed_policy.independent_member_indices(list(candidates))
    assert first == second == (0, 2)


def test_prepared_clusters_are_deterministic_disjoint_and_lossless():
    candidates = [
        Candidate("AAAA", start=0, score=4.0),
        Candidate("CCCC", start=10, score=3.0),
        Candidate("GGGGGGGG", start=20, score=2.0),
    ]
    seed_sets = [(["AAAA", "CCCC"], 4), (["GGGGGGGG"], 8)]
    first = seed_policy.prepare_seed_clusters(candidates, seed_sets)
    second = seed_policy.prepare_seed_clusters(candidates, seed_sets)

    assert [[member.seq for member in cluster.members] for cluster in first] == [
        ["AAAA", "CCCC"],
        ["GGGGGGGG"],
    ]
    assert [[member.seq for member in cluster.members] for cluster in first] == [
        [member.seq for member in cluster.members] for cluster in second
    ]
    pooled = [member.seq for cluster in first for member in cluster.members]
    assert pooled == [candidate.seq for candidate in candidates]
    assert len(pooled) == len(set(pooled))


def test_member_matching_uses_seed_clusters_eligible_dyad_pool():
    sequence = "ACGTACGT"
    filtered_pwm = Candidate(sequence, start=100, score=99.0, kind="pwm")
    eligible_ir = Candidate(sequence, start=10, score=3.0, kind="IR")
    cluster = _cluster([filtered_pwm, eligible_ir], [sequence], len(sequence))

    assert cluster.member_metadata[0]["type"] == "IR"
    assert cluster.member_metadata[0]["start"] == 10


def test_duplicate_occurrences_bind_once_each_in_eligible_order():
    first = Candidate("ACGTACGT", start=0, score=4.0, kind="IR")
    second = Candidate("ACGTACGT", start=30, score=3.0, kind="IR")
    cluster = _cluster([second, first], [first.seq, second.seq], len(first.seq))

    assert [member["start"] for member in cluster.member_metadata] == [0, 30]
    assert len({member["start"] for member in cluster.member_metadata}) == 2


def test_policy_uses_effective_support_not_raw_cluster_size_and_falls_back_to_cluster_zero():
    policy = seed_policy.EFFECTIVE_SUPPORT_3_PROPOSAL
    cluster0_candidates = [
        Candidate("AAAA", start=0, score=3.0),
        Candidate("AAAT", start=10, score=2.0),
        Candidate("GGGGG", start=20, score=1.0),
    ]
    cluster1_candidates = [
        Candidate("ACGTACGT", start=30, score=3.0),
        Candidate("ACGTACGA", start=50, score=2.0),
    ]
    candidates = cluster0_candidates + cluster1_candidates
    clusters = seed_policy.prepare_seed_clusters(
        candidates,
        [
            ([candidate.seq for candidate in cluster0_candidates], 4),
            ([candidate.seq for candidate in cluster1_candidates], 8),
        ],
        policy=policy,
    )
    trials = [
        seed_policy.build_trial(
            cluster,
            [_hit(cluster.prepared.sequences[0], 6.0, start=100 + 100 * cluster.cluster_id),
             _hit(cluster.prepared.sequences[0], 5.0, start=120 + 100 * cluster.cluster_id),
             _hit(cluster.prepared.sequences[0], 4.0, start=140 + 100 * cluster.cluster_id)],
            policy=policy,
        )
        for cluster in clusters
    ]
    decision = seed_policy.choose_trial(trials, policy=policy)
    recorded_v5 = seed_policy.choose_trial(trials, policy=seed_policy.V5_RAW_INPUT_SUPPORT_3_CONTROL)

    assert trials[0].audit["n_input"] == 3
    assert trials[0].audit["n_effective"] == 2
    assert decision.chosen.cluster_id == 0
    assert decision.eligible_cluster_ids == ()
    assert decision.fallback_reason is not None
    assert recorded_v5.eligible_cluster_ids == (0,)
    assert recorded_v5.fallback_reason is None
    report = decision.to_dict(trials)
    assert report["fallback_used"] is True
    assert report["parameters"]["support_basis"] == "effective"
    assert report["parameters"]["min_support"] == 3


def test_supported_trial_with_best_mean_ic_wins_without_fallback():
    candidates = [
        Candidate("AAAA", start=0, score=6.0),
        Candidate("AAAT", start=10, score=5.0),
        Candidate("AATT", start=20, score=4.0),
        Candidate("ACGTACGT", start=40, score=3.0),
        Candidate("ACGTACGA", start=60, score=2.0),
        Candidate("ACGTACGG", start=80, score=1.0),
    ]
    clusters = seed_policy.prepare_seed_clusters(
        candidates,
        [([candidate.seq for candidate in candidates[:3]], 4),
         ([candidate.seq for candidate in candidates[3:]], 8)],
    )
    low_ic = seed_policy.build_trial(
        clusters[0],
        [_hit("AAAA", 6.0), _hit("CCCC", 5.0), _hit("GGGG", 4.0)],
    )
    high_ic = seed_policy.build_trial(
        clusters[1],
        [_hit("ACGTACGT", 6.0), _hit("ACGTACGT", 5.0), _hit("ACGTACGT", 4.0)],
    )
    # Mean IC is no longer the production selector (SEED_DECISION_20260819.md), so this exercises the
    # retained control explicitly rather than through the default.
    decision = seed_policy.choose_trial([low_ic, high_ic], policy=seed_policy.V4_MEAN_IC_CONTROL)

    assert high_ic.mean_ic > low_ic.mean_ic
    assert decision.chosen.cluster_id == 1
    assert decision.eligible_cluster_ids == (0, 1)
    assert decision.fallback_reason is None


def test_explicit_cluster0_control_ignores_mean_ic_and_replays_deterministically():
    candidates = [Candidate("AAAA", start=0), Candidate("ACGTACGT", start=20)]
    clusters = seed_policy.prepare_seed_clusters(
        candidates,
        [([candidates[0].seq], 4), ([candidates[1].seq], 8)],
    )
    trials = [
        seed_policy.build_trial(
            clusters[0],
            [_hit("AAAA", 3.0), _hit("CCCC", 2.0), _hit("GGGG", 1.0)],
        ),
        seed_policy.build_trial(
            clusters[1],
            [_hit("ACGTACGT", 3.0), _hit("ACGTACGT", 2.0), _hit("ACGTACGT", 1.0)],
        ),
    ]
    assert trials[1].mean_ic > trials[0].mean_ic

    decision = seed_policy.choose_trial(trials, policy=seed_policy.CLUSTER0_CONTROL)
    replay = seed_policy.choose_trial_records(
        [trial.to_dict() for trial in reversed(trials)],
        policy=seed_policy.CLUSTER0_CONTROL,
    )

    assert decision.chosen.cluster_id == replay.chosen.cluster_id == 0
    assert decision.to_dict(trials)["chosen_why"] == "explicit cluster-0 control"
    assert decision.fallback_reason is None


def test_cluster0_control_rejects_missing_provenance_cluster_zero():
    candidates = [Candidate("AAAA", start=0), Candidate("ACGTACGT", start=20)]
    clusters = seed_policy.prepare_seed_clusters(
        candidates,
        [([candidates[0].seq], 4), ([candidates[1].seq], 8)],
    )
    cluster_one_only = seed_policy.build_trial(
        clusters[1],
        [_hit("ACGTACGT", 3.0), _hit("ACGTACGT", 2.0), _hit("ACGTACGT", 1.0)],
    )

    with pytest.raises(ValueError, match=r"cluster0 selection requires exactly one.*cluster_id == 0"):
        seed_policy.choose_trial([cluster_one_only], policy=seed_policy.CLUSTER0_CONTROL)


def test_cluster0_fallback_rejects_missing_provenance_cluster_zero():
    candidates = [Candidate("AAAA", start=0), Candidate("ACGTACGT", start=20)]
    clusters = seed_policy.prepare_seed_clusters(
        candidates,
        [([candidates[0].seq], 4), ([candidates[1].seq], 8)],
        policy=seed_policy.V5_RAW_INPUT_SUPPORT_3_CONTROL,
    )
    cluster_one_only = seed_policy.build_trial(
        clusters[1],
        [_hit("ACGTACGT", 2.0)],
        policy=seed_policy.V5_RAW_INPUT_SUPPORT_3_CONTROL,
    )

    with pytest.raises(ValueError, match=r"cluster_0 fallback requires exactly one.*cluster_id == 0"):
        seed_policy.choose_trial(
            [cluster_one_only],
            policy=seed_policy.V5_RAW_INPUT_SUPPORT_3_CONTROL,
        )


def test_duplicate_cluster_ids_are_rejected_at_choose_and_replay_boundaries():
    candidate = Candidate("AAAA", start=0)
    cluster = _cluster([candidate], [candidate.seq], 4)
    trial = seed_policy.build_trial(
        cluster,
        [_hit("AAAA", 3.0), _hit("AAAT", 2.0), _hit("AATT", 1.0)],
    )
    duplicate_message = r"cluster_id values must be unique; duplicate cluster_id\(s\): 0"

    with pytest.raises(ValueError, match=duplicate_message):
        seed_policy.choose_trial([trial, trial])
    with pytest.raises(ValueError, match=duplicate_message):
        seed_policy.choose_trial_records([trial.to_dict(), trial.to_dict()])


def test_default_policy_is_cluster0_and_the_legacy_controls_still_reproduce_v4_and_v5():
    candidates = [
        Candidate("AAAA", start=0, score=4.0),
        Candidate("ACGTACGT", start=20, score=3.0),
    ]
    clusters = seed_policy.prepare_seed_clusters(
        candidates,
        [([candidates[0].seq], 4), ([candidates[1].seq], 8)],
    )
    cluster0 = seed_policy.build_trial(
        clusters[0],
        [_hit("AAAA", 6.0), _hit("AAAT", 5.0), _hit("AATT", 4.0)],
    )
    cluster1 = seed_policy.build_trial(
        clusters[1],
        [_hit("ACGTACGT", 6.0), _hit("ACGTACGT", 5.0), _hit("ACGTACGT", 4.0)],
    )
    trials = [cluster0, cluster1]

    legacy_v4 = max(
        [trial for trial in trials if trial.n_hits >= 3] or trials,
        key=lambda trial: (trial.mean_ic, trial.n_hits),
    )
    default = seed_policy.choose_trial(trials)
    v4 = seed_policy.choose_trial(trials, policy=seed_policy.V4_MEAN_IC_CONTROL)
    support3 = seed_policy.choose_trial(trials, policy=seed_policy.V5_RAW_INPUT_SUPPORT_3_CONTROL)
    support3_replay = seed_policy.choose_trial_records(
        json.loads(json.dumps([trial.to_dict() for trial in trials])),
        policy=seed_policy.V5_RAW_INPUT_SUPPORT_3_CONTROL,
    )

    # Production is cluster0 as of 2026-08-19 (SEED_DECISION_20260819.md): selected 139/139 held out.
    assert seed_policy.DEFAULT_POLICY.policy_name == "cluster0_promoted_20260819"
    assert seed_policy.DEFAULT_POLICY.selection_method == "cluster0"
    assert seed_policy.DEFAULT_POLICY.fallback == "cluster_0"
    assert seed_policy.DEFAULT_POLICY.support_basis == "input"
    assert seed_policy.DEFAULT_POLICY.min_support == 1
    assert default.chosen.cluster_id == 0

    # This case is exactly the pathology that got mean IC demoted: cluster 1 is built from ONE sequence,
    # its three genome hits are near-exact copies of it, and it therefore wins on self-fit IC while
    # cluster 0 -- the wider-support anchor -- is the one production now keeps.
    assert seed_policy.V4_MEAN_IC_CONTROL.policy_name == "run6_v4_mean_ic_control"
    assert v4.chosen.cluster_id == legacy_v4.cluster_id == 1
    assert v4.chosen.cluster_id != default.chosen.cluster_id

    assert seed_policy.V5_RAW_INPUT_SUPPORT_3_CONTROL.support_basis == "input"
    assert support3.chosen.cluster_id == 0
    assert support3_replay.chosen.cluster_id == support3.chosen.cluster_id
    assert support3.fallback_reason is not None


def test_default_v4_control_compares_all_trials_when_none_has_three_hits():
    candidates = [Candidate("AAAA", start=0), Candidate("ACGTACGT", start=20)]
    clusters = seed_policy.prepare_seed_clusters(
        candidates,
        [([candidates[0].seq], 4), ([candidates[1].seq], 8)],
    )
    trials = [
        seed_policy.build_trial(clusters[0], [_hit("AAAA", 2.0), _hit("CCCC", 1.0)]),
        seed_policy.build_trial(clusters[1], [_hit("ACGTACGT", 2.0), _hit("ACGTACGT", 1.0)]),
    ]
    legacy_v4 = max(trials, key=lambda trial: (trial.mean_ic, trial.n_hits))
    decision = seed_policy.choose_trial(trials, policy=seed_policy.V4_MEAN_IC_CONTROL)

    assert decision.chosen.cluster_id == legacy_v4.cluster_id == 1
    assert decision.fallback_reason and "v4 control" in decision.fallback_reason


def test_replay_ties_choose_lower_cluster_id_regardless_of_record_order():
    candidates = [Candidate("AAAA", start=0), Candidate("CCCC", start=20)]
    clusters = seed_policy.prepare_seed_clusters(
        candidates,
        [([candidates[0].seq], 4), ([candidates[1].seq], 4)],
    )
    identical_hits = [_hit("ACGT", 4.0), _hit("ACGT", 3.0), _hit("ACGT", 2.0)]
    trials = [seed_policy.build_trial(cluster, identical_hits) for cluster in clusters]
    records = [trial.to_dict() for trial in reversed(trials)]

    assert trials[0].mean_ic == trials[1].mean_ic
    assert seed_policy.choose_trial_records(records).chosen.cluster_id == 0


def test_trial_audit_is_json_serializable_and_hit_bound_is_explicit():
    policy = seed_policy.SeedPolicy(audit_hit_limit=2)
    candidates = [
        Candidate("ACGT", start=0, score=3.0),
        Candidate("ACGA", start=10, score=2.0),
        Candidate("ACGG", start=20, score=1.0),
    ]
    cluster = _cluster(
        candidates,
        [candidate.seq for candidate in candidates],
        4,
        policy=policy,
    )
    trial = seed_policy.build_trial(
        cluster,
        [_hit("ACGT", 6.0, tier=1), _hit("ACGA", 5.0, tier=2), _hit("ACGG", 4.0, tier=3)],
        policy=policy,
    )
    decision = seed_policy.choose_trial([trial], policy=policy)
    payload = decision.to_dict([trial])

    json.dumps(payload)
    audit = payload["trials"][0]
    assert audit["n_hits"] == 3
    assert audit["n_hits_recorded"] == 2
    assert audit["ranked_hits_limit"] == 2
    assert audit["hits_truncated"] is True
    assert set(audit["ranked_hits"][0]) >= {
        "start", "end", "dyad", "pvalue", "qvalue", "tier", "combined_score",
    }
    assert audit["pwm_prior"]["pseudocount_per_base"] == 0.25
    assert audit["ic_curve"] == trial.audit["ic_curve"]

    replayed = seed_policy.choose_trial_records(json.loads(json.dumps(payload["trials"])))
    assert replayed.chosen.cluster_id == decision.chosen.cluster_id
    assert replayed.chosen.to_dict()["ranked_hits"] == audit["ranked_hits"]


def test_policy_provenance_rejects_trial_build_or_replay_mismatch():
    candidates = [Candidate("ACGT", start=0)]
    cluster = _cluster(candidates, [candidates[0].seq], 4)
    different_prior = seed_policy.SeedPolicy(pseudocount_per_base=0.5)

    with pytest.raises(ValueError, match="does not match"):
        seed_policy.build_trial(cluster, [_hit("ACGT", 4.0)], policy=different_prior)

    trial = seed_policy.build_trial(cluster, [_hit("ACGT", 4.0)])
    with pytest.raises(ValueError, match="fixed during trial construction"):
        seed_policy.choose_trial([trial], policy=different_prior)


def test_unsupported_report_only_policy_options_are_rejected():
    with pytest.raises(ValueError, match="unsupported selection_method"):
        seed_policy.SeedPolicy(selection_method="oracle")
    with pytest.raises(ValueError, match="unsupported fallback"):
        seed_policy.SeedPolicy(fallback="most_independent")
    with pytest.raises(ValueError, match="unsupported hit_ranking"):
        seed_policy.SeedPolicy(hit_ranking="combined_score")
    with pytest.raises(ValueError, match="unsupported support_basis"):
        seed_policy.SeedPolicy(support_basis="raw")
    with pytest.raises(ValueError, match="unsupported rescan_scope"):
        seed_policy.SeedPolicy(rescan_scope="promoters")
    with pytest.raises(ValueError, match="pseudocount_per_base"):
        seed_policy.SeedPolicy(pseudocount_per_base=0)
    with pytest.raises(ValueError, match="keep_frac"):
        seed_policy.SeedPolicy(keep_frac=1.1)
    with pytest.raises(ValueError, match="rescan_pvalue_thresh"):
        seed_policy.SeedPolicy(rescan_pvalue_thresh=0)


def test_native_genome_hits_are_oriented_and_build_nonempty_logo_metrics():
    candidate = Candidate("ACGT", start=0)
    cluster = _cluster([candidate], [candidate.seq], 4)
    plus_hit = SimpleNamespace(
        accession="TEST_GENOME",
        start=100,
        end=104,
        strand="+",
        dyad_center=102,
        matched="ACGT",
        score=7.0,
        pvalue=1e-7,
        qvalue=1e-5,
        tier=2,
        combined=6.778,
    )
    minus_hit = SimpleNamespace(
        accession="TEST_GENOME",
        start=200,
        end=204,
        strand="-",
        dyad_center=202,
        matched="AACG",
        score=6.0,
        pvalue=1e-6,
        qvalue=1e-4,
        tier=3,
        combined=5.523,
    )
    trial = seed_policy.build_trial(cluster, [plus_hit, minus_hit])

    assert trial.audit["ranked_hits"][0]["combined_score"] == 6.778
    assert [hit.seq for hit in trial.hits] == ["ACGT", "CGTT"]
    assert trial.audit["n_keep"] == 2
    assert trial.audit["mean_ic"] > 0
    assert trial.audit["consensus"]


def test_trial_preserves_rescan_order_while_auditing_selector_rank_order():
    candidates = [
        Candidate("ACGT", start=0, score=3.0),
        Candidate("ACGA", start=10, score=2.0),
        Candidate("ACGG", start=20, score=1.0),
    ]
    cluster = _cluster(candidates, [candidate.seq for candidate in candidates], 4)
    combined_first = _hit("ACGT", 2.0, start=100)
    raw_score_first = _hit("ACGA", 8.0, start=200)
    third = _hit("ACGG", 4.0, start=300)
    trial = seed_policy.build_trial(cluster, [combined_first, raw_score_first, third])

    assert trial.hits == (combined_first, raw_score_first, third)
    assert [hit["start"] for hit in trial.audit["ranked_hits"]] == [200, 300, 100]


def test_public_selector_api_has_no_family_effector_ssn_or_truth_inputs():
    banned = ("family", "effector", "ssn", "truth", "known_site")
    functions = (
        seed_policy.prepare_seed_clusters,
        seed_policy.build_trial,
        seed_policy.choose_trial,
        seed_policy.trial_from_audit,
        seed_policy.choose_trial_records,
        seed_policy.independent_member_indices,
    )
    for function in functions:
        parameters = inspect.signature(function).parameters
        assert not any(term in parameter.lower() for parameter in parameters for term in banned)
