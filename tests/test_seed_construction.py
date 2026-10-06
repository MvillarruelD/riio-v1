"""Offline contract tests for the private Phase-B construction bridge."""
from __future__ import annotations

import inspect
from types import SimpleNamespace

import pytest

from predictor import api, pipeline
from predictor.signals import seed_selection as seed


def test_default_construction_is_the_promoted_production_policy():
    """Pins the 2026-08-19 promotion (analysis/regulondb_bench/SEED_DECISION_20260819.md).

    top_k is 8, not run 6's 24: the extra candidates move the modal length, so a different candidate
    builds the PWM and IscR and MntR are lost for ArsR alone, at identical held-out recovery.
    """
    assert seed.DEFAULT_CONSTRUCTION == seed.SeedConstruction(
        top_k=8,
        tol=3,
        pseudocount_per_base=0.25,
    )
    assert seed.DEFAULT_POLICY.policy_name == "cluster0_promoted_20260819"
    assert seed.DEFAULT_POLICY.selection_method == "cluster0"


def test_private_construction_routes_top_k_tol_prior_and_no_cluster_cap(monkeypatch):
    calls = []
    candidates = [SimpleNamespace(seq="ACGT", kind="IR", conservation=1.0, final=1.0)]
    clustered = [(["ACGT"], 4)]

    def fake_predict(*args, **kwargs):
        calls.append(("predict", kwargs))
        return candidates

    def fake_seed_clusters(value, **kwargs):
        assert value is candidates
        calls.append(("seed_clusters", kwargs))
        return clustered

    def fake_prepare(value, sets, *, policy):
        assert value is candidates and sets is clustered
        calls.append(("prepare", policy.pseudocount_per_base))
        return ["prepared"]

    monkeypatch.setattr(pipeline.motif_finder, "predict", fake_predict)
    monkeypatch.setattr(pipeline.motif_finder, "seed_clusters", fake_seed_clusters)
    monkeypatch.setattr(pipeline.ss, "prepare_seed_clusters", fake_prepare)

    for top_k in (8, 24):
        construction = seed.SeedConstruction(
            top_k=top_k,
            tol=1,
            pseudocount_per_base=0.5,
        )
        _, _, clusters, policy = pipeline._construct_seed_clusters(
            "ACGT" * 20,
            family="MerR",
            homolog_regions=[],
            construction=construction,
        )
        assert clusters == ["prepared"]
        assert policy.pseudocount_per_base == 0.5

    predict_calls = [payload for kind, payload in calls if kind == "predict"]
    cluster_calls = [payload for kind, payload in calls if kind == "seed_clusters"]
    assert [call["top_k"] for call in predict_calls] == [8, 24]
    assert all(call["tol"] == 1 for call in cluster_calls)
    assert all(call["max_clusters"] is None for call in cluster_calls)


def test_source_provenance_serializes_with_harness_field_names():
    construction = seed.SeedConstruction(
        top_k=8,
        tol=3,
        pseudocount_per_base=0.25,
        source_run_id="offline-fixture",
        implementation_commit="a" * 40,
        implementation_dirty_diff_sha256="b" * 64,
        generator_source_sha256="c" * 64,
        source_runner_sha256="d" * 64,
        fixed_family_inputs_sha256="2" * 64,
        mirror_genome_id="GCF_TEST",
        expected_sequence_accession="NC_TEST",
        mirror_genome_sha256="e" * 64,
        panel_fasta_sha256="f" * 64,
        construction_grid_sha256="1" * 64,
        generation_context_sha256="3" * 64,
    )
    payload = construction.to_dict()
    assert payload == {
        "top_k": 8,
        "tol": 3,
        "pseudocount_per_base": 0.25,
        # the dyad scoring rule belongs in the provenance: it decides which candidate becomes the
        # seed, so a bundle that does not record it cannot say what produced its operator
        "score_mode": "raw",
        "source_run_id": "offline-fixture",
        "implementation_commit": "a" * 40,
        "implementation_dirty_diff_sha256": "b" * 64,
        "generator_source_sha256": "c" * 64,
        "source_runner_sha256": "d" * 64,
        "fixed_family_inputs_sha256": "2" * 64,
        "mirror_genome_id": "GCF_TEST",
        "expected_sequence_accession": "NC_TEST",
        "mirror_genome_sha256": "e" * 64,
        "panel_fasta_sha256": "f" * 64,
        "construction_grid_sha256": "1" * 64,
        "generation_context_sha256": "3" * 64,
    }

    decision = SimpleNamespace(to_dict=lambda trials: {"trials": list(trials)})
    report = pipeline._seed_selection_report(decision, ["trial"], construction)
    assert report["construction"] == payload


def test_partial_source_provenance_is_rejected():
    try:
        seed.SeedConstruction(source_run_id="incomplete")
    except ValueError as exc:
        assert "supplied completely" in str(exc)
    else:
        raise AssertionError("partial source provenance must fail closed")


def test_public_api_signature_has_no_construction_parameter():
    parameters = inspect.signature(api.predict).parameters
    assert "_seed_construction" not in parameters
    assert list(parameters) == [
        "seq", "name", "family", "organism", "genome_acc", "effector", "fold", "cfg", "verbose", "kw",
    ]


def test_public_api_rejects_private_construction_before_pipeline(monkeypatch):
    monkeypatch.setattr(api, "_resolve_folds", lambda cfg, fold: cfg)
    with pytest.raises(TypeError, match="not accepted by the public API"):
        api.predict("MPEPTIDE", _seed_construction=seed.DEFAULT_CONSTRUCTION)


def test_internal_benchmark_seam_routes_complete_construction(monkeypatch):
    construction = seed.SeedConstruction(
        source_run_id="fixture",
        implementation_commit="a" * 40,
        implementation_dirty_diff_sha256="b" * 64,
        generator_source_sha256="c" * 64,
        source_runner_sha256="d" * 64,
        fixed_family_inputs_sha256="2" * 64,
        mirror_genome_id="GCF_TEST",
        expected_sequence_accession="NC_TEST",
        mirror_genome_sha256="e" * 64,
        panel_fasta_sha256="f" * 64,
        construction_grid_sha256="1" * 64,
        generation_context_sha256="3" * 64,
    )
    captured = {}

    def fake_impl(seq, **kwargs):
        captured.update(seq=seq, **kwargs)
        return {"tf_id": kwargs["name"]}

    monkeypatch.setattr(pipeline, "_run_novel_impl", fake_impl)
    result = pipeline._run_seed_construction_benchmark(
        "MPEPTIDE",
        construction=construction,
        name="CueR_seed_sweep_top8",
        family="MerR",
        organism="GCF_TEST",
    )
    assert result == {"tf_id": "CueR_seed_sweep_top8"}
    assert captured["seed_construction"] is construction
