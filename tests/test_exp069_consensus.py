"""Synthetic finalization, persistence and offline-report regressions for EXP069 package 2."""
from dataclasses import asdict
import json
from types import SimpleNamespace as NS

import pytest

from predictor import pipeline
from predictor.effector import inducer as IND, ligify_db, operon_chemistry as OC
from predictor.report import figure_data, ligand_regulon as LR, outputs, render_report, run_page, run_report
from predictor.schema import InducerCall as Call, InducerConsensus as Consensus


@pytest.mark.parametrize("confidence", [0, -0.3, None])
def test_nonpositive_explicit_sets_never_become_operative(confidence):
    raw = Call("raw", "rejected", confidence, candidates=("rejected", "also rejected"),
               mixture_kind="unresolved")
    cons = Consensus.from_calls([raw, Call("usable", "citrate", 0.5)])
    assert cons.top == "citrate"
    assert cons.candidates == () and cons.mixture_kind == ""
    assert cons.agreement == 1 and cons.compatibility == 0
    assert cons.calls[0] == raw and cons.calls[0].candidates == ("rejected", "also rejected")
    none = Consensus.from_calls([raw])
    assert none.top is None and none.candidates == () and none.mixture_kind == ""
    assert none.agreement == none.compatibility == 0


def test_exact_support_is_distinct_from_unresolved_set_compatibility():
    unresolved = Call("set", "Ni2+", 0.9, "metal", candidates=("Ni2+", "Co2+"), mixture_kind="unresolved")
    exact = Call("exact", "Ni2+", 0.8, "metal")
    cons = IND._finalize([unresolved, exact], None, sensor_class="metal")
    assert cons.top == "Ni2+" and cons.candidates == ("Ni2+", "Co2+")
    assert cons.support_sources == ("exact",) and cons.compatible_sources == ("set",)
    assert cons.agreement == 0.5 and cons.compatibility == 0.5
    assert cons.mixture_kind == "unresolved" and cons.finalized
    # A positive set-only source retains alternatives without inventing a single compound.
    set_only = IND._finalize([Call("set", None, 0.5, "metal", candidates=("Ni2+", "Co2+"),
                                  mixture_kind="unresolved")], None, sensor_class="metal")
    assert set_only.top == "divalent metal (ion unresolved)"
    assert set_only.candidates == ("Ni2+", "Co2+") and set_only.agreement == 0


@pytest.mark.parametrize("gate", [None, True])
@pytest.mark.parametrize("source", ["operon_chem", "ligify"])
def test_guard_rebuilds_support_and_candidates_after_changing_headline(gate, source):
    raw = Call(source, "discarded local proposal", 0.8, evidence={"weight": 0.1 if gate else 1})
    accepted = Call("ligify_db", "citrate", 0.3)
    cons = IND._finalize([raw, accepted], gate, sensor_class="metal")
    assert cons.top == cons.top_display == "citrate"
    assert cons.candidates == () and cons.mixture_kind == ""
    assert cons.agreement == 0.5 and cons.support_sources == ("ligify_db",)
    assert cons.accepted_calls == [accepted]
    assert cons.calls == [raw, accepted] and cons.calls[0].confidence == 0.8


def test_positive_gate_cannot_revive_a_zero_confidence_named_alternative():
    calls = [Call("ssn_cluster", "raw nonmetal", 0.9, "non-metal", candidates=("raw nonmetal", "other"),
                  mixture_kind="unresolved"), Call("other", "Zn2+", 0, "metal")]
    cons = IND._finalize(calls, True, sensor_class="metal")
    assert cons.top == cons.top_display == "divalent metal (ion unresolved)"
    assert cons.candidates == () and cons.mixture_kind == ""
    assert cons.agreement == cons.compatibility == 0 and cons.accepted_calls == []
    assert cons.calls == calls


def test_redox_early_return_rebuilds_all_headline_fields():
    calls = [Call("ssn_cluster", "Fe2+", 0.9, "metal", candidates=("Fe2+", "Zn2+"),
                  mixture_kind="unresolved"), Call("coordination", None, 0, "redox")]
    cons = IND._finalize(calls, None, sensor_class="redox")
    assert cons.top == cons.top_display == "reactive oxygen/nitrogen species (redox)"
    assert cons.candidates == () and cons.mixture_kind == "" and cons.finalized
    assert cons.agreement == cons.compatibility == 0 and cons.accepted_calls == []
    assert cons.calls == calls


@pytest.mark.parametrize("top", [None, "divalent metal (ion unresolved)", "accepted alternative"])
def test_pipeline_uses_final_headline_instead_of_first_raw_source(monkeypatch, top):
    cons = Consensus(top=top, finalized=True, calls=[Call("ssn_cluster", "rejected raw", 0.9)])
    monkeypatch.setattr(IND, "infer_inducer", lambda *a, **k: cons)
    ctx = NS(genes=[])
    result, _, effector, error = pipeline._infer_inducer("SYNTHETIC", "test", ctx, 0, None,
                                                      allow_ncbi=False, verbose=False)
    assert result is cons and effector == top and error is None
    assert pipeline._infer_inducer("SYNTHETIC", "test", ctx, 0, "user override",
                                   allow_ncbi=False, verbose=False)[2] == "user override"


def _stub_inference(monkeypatch):
    monkeypatch.setattr(IND.hsel, "assign_cluster", lambda *a: NS(cluster_id=None))
    monkeypatch.setattr(IND, "curated_member_ion", lambda *a, **k: None)
    monkeypatch.setattr(ligify_db, "lookup", lambda **k: None)
    monkeypatch.setattr(IND.metal_site, "coordination_gate", lambda *a, **k: {"has_site": None})
    monkeypatch.setattr(IND.metal_site, "ligify_weight", lambda *a, **k: (1, "synthetic"))
    monkeypatch.setattr(IND.metal_site, "inducer_class", lambda *a: ("organic", "synthetic"))
    from predictor.annotate import neighborhood
    monkeypatch.setattr(neighborhood, "analyse", lambda *a: None)
    monkeypatch.setattr(neighborhood, "inducer_call", lambda *a: None)


def test_original_chemistry_snapshot_is_portable_and_reporting_never_fetches_again(monkeypatch, tmp_path):
    _stub_inference(monkeypatch)
    acquired = []

    def fetch(protein_id):
        acquired.append(protein_id)
        assert len(acquired) == 1, "a report must not fetch a second annotation"
        return {"comments": [{"reaction": "original acquired record"}], "nested": [1, None, True]}

    def derive(ctx, tf, *, entry_fetcher):
        entry_fetcher("SYNTHETIC_PROTEIN")
        return [OC.EffectorCandidate("citrate", 80, "", chebi=["raw reaction ID"],
                                     equation="original equation", rhea="RHEA:synthetic"),
                OC.EffectorCandidate("raw alternative", 70, "", equation="second equation")]

    monkeypatch.setattr(OC, "effectors_from_context", derive)
    cons = IND.infer_inducer("SYNTHETIC", "test", genome_ctx=NS(), tf_gene=NS(), allow_ncbi=True,
                            metalnet=False, entry_fetcher=fetch)
    recorded = json.loads(json.dumps(asdict(cons), allow_nan=False))
    snapshot = recorded["chemistry_snapshot"]
    assert snapshot["status"] == "recorded" and len(snapshot["candidates"]) == 2
    assert snapshot["candidates"][0]["chebi"] == ["raw reaction ID"]
    assert snapshot["annotations"][0]["entry"]["comments"][0]["reaction"] == "original acquired record"
    assert snapshot["annotations"][0]["availability"] == "injected_unverified"
    assert recorded["calls"] == json.loads(json.dumps([asdict(c) for c in cons.calls]))
    # Nothing downstream can use a changed response, even if an acquisition helper is called by mistake.
    monkeypatch.setattr(LR, "infer_ligand", lambda *a, **k: pytest.fail("second chemistry inference"))
    accepted = LR.accepted_ligands(recorded)
    assert [p["ligand"] for p in accepted] == ["citrate"]
    assert accepted[0]["equation"] == "original equation"
    dossier = {"tf_id": "synthetic", "family": "test", "inducers": recorded,
               "ligand": accepted, "operators": {"candidates": []}, "effector": "explicit override",
               "effector_origin": "user_override", "effector_override": "explicit override"}
    html = render_report.build_html({"dossier": dossier}, {})
    assert "Raw acquired chemistry proposals" in html and "raw alternative" in html
    primary = {"source": "test", "sequence": "", "score": 0, "provenance": {}, "status": "test"}
    operators = {"primary": primary, "ranked": [primary], "candidates": []}
    bundle = outputs.write_stage_a(dossier, operators=operators, root=tmp_path, render=False)
    saved = json.loads((bundle / "dossier.json").read_text())
    assert saved["inducers"]["chemistry_snapshot"] == snapshot
    assert saved["inducers"]["calls"] == recorded["calls"]
    inducer_export = json.loads((bundle / "ligand" / "inducer.json").read_text())
    assert inducer_export["effector"] == inducer_export["effector_override"] == "explicit override"
    assert inducer_export["effector_origin"] == "user_override"
    assert inducer_export["inducers"]["top"] == "citrate"
    assert acquired == ["SYNTHETIC_PROTEIN"]


def test_legacy_missing_snapshot_does_not_authorize_acquisition(monkeypatch):
    monkeypatch.setattr(OC, "effectors_from_context", lambda *a, **k: pytest.fail("legacy fetch"))
    legacy = {"top": "legacy ligand", "calls": []}
    assert LR.accepted_ligands(legacy)[0]["chemistry_status"] == "legacy_unavailable"
    html = render_report.build_html({"dossier": {"tf_id": "legacy", "inducers": legacy}}, {})
    assert "legacy / unavailable" in html and "legacy ligand" in html


def test_saved_final_empty_set_does_not_revive_rejected_raw_ssn_candidates():
    cons = IND._finalize([Call("ssn_cluster", "rejected label", 0.9, "non-metal",
                              candidates=("rejected label", "rejected alternative"),
                              mixture_kind="unresolved")], True, sensor_class="metal")
    dossier = {"tf_id": "synthetic", "inducers": json.loads(json.dumps(asdict(cons)))}
    html = render_report.build_html({"dossier": dossier}, {})
    assert "Raw source set: rejected label, rejected alternative" in html
    assert "advisory / rejected by final policy" in html
    assert "UNRESOLVED between rejected label and rejected alternative.</b>" not in html
    assert "Candidate set:</b> rejected label" not in html


@pytest.mark.parametrize("finalized", [None, True])
def test_explicit_final_null_is_authoritative_in_stored_run_and_page(tmp_path, finalized):
    directory = tmp_path / "jobs" / "synthetic"
    directory.mkdir(parents=True)
    dossier = {"tf_id": "synthetic", "effector": "stale effector", "family": "test",
               "inducers": {"top": None, "top_display": "stale raw display",
                            "calls": [{"source": "ssn_cluster", "ligand": "Zn2+", "confidence": 0,
                                       "role": "metal"}], "accepted_calls": []}}
    if finalized is not None:
        dossier["inducers"]["finalized"] = finalized
    (directory / "dossier.json").write_text(json.dumps(dossier))
    (tmp_path / "summary.tsv").write_text("run_name\tinducer_top\nsynthetic\tstale summary\n")
    rows, _ = run_report._run_rows(tmp_path)
    assert rows[0]["inducer_top"] is None and rows[0]["inducer_class"] == "unresolved"
    assert figure_data.inducer_top(dossier, "stale summary") is None
    assert figure_data.inducer_display(dossier["inducers"]) == ""
    assert "stale" not in run_page._inducer_cell(rows[0], dossier)
    assert "stale raw display" not in render_report.build_html({"dossier": dossier}, {})
    row = {"inducer_class": "metal", "inducer_top": "divalent metal (ion unresolved)"}
    assert "SSN" not in run_page._inducer_cell(row, {"inducers": {"top": row["inducer_top"],
                                                                 "calls": dossier["inducers"]["calls"]}})


def test_genuinely_absent_legacy_top_can_use_documented_fallbacks():
    assert figure_data.inducer_display({"top_display": "legacy display"}) == "legacy display"
    assert figure_data.inducer_top({"effector": "legacy effector"}, "legacy summary") == "legacy effector"
    assert figure_data.inducer_top({}, "legacy summary") == "legacy summary"


def test_legacy_set_membership_agreement_is_not_relabelled_as_exact_support():
    raw = Call("set", "Ni2+", 0.9, "metal", candidates=("Ni2+", "Co2+"), mixture_kind="unresolved")
    legacy = {"top": "Ni2+", "agreement": 1.0, "calls": [asdict(raw)],
              "candidates": ["Ni2+", "Co2+"], "mixture_kind": "unresolved"}
    legacy_html = render_report.build_html({"dossier": {"tf_id": "legacy", "inducers": legacy}}, {})
    assert "legacy source agreement 1.00" in legacy_html
    assert "exact source agreement 1.00" not in legacy_html
    assert "Exact support was not recorded separately" in legacy_html
    assert legacy["agreement"] == 1.0
    legacy_ev = figure_data.inducer_evidence_rows({"inducers": legacy})
    assert legacy_ev["agreement_label"] == "legacy source agreement"
    assert "legacy agreement" in next(c["detail"] for c in figure_data.verdict({"inducers": legacy})
                                      if c["label"] == "inducer")
    finalized = asdict(IND._finalize([raw], None, sensor_class="metal"))
    new_html = render_report.build_html({"dossier": {"tf_id": "new", "inducers": finalized}}, {})
    assert "exact source agreement 0.00" in new_html and "unresolved-set compatibility 0.00" in new_html
    assert figure_data.inducer_evidence_rows({"inducers": finalized})["agreement_label"] == "exact source agreement"


@pytest.mark.parametrize("top", [None, "divalent metal (ion unresolved)", "different headline"])
def test_singleton_accepted_alternative_survives_without_manufacturing_a_mixture(top):
    accepted = Call("set", None, 0.6, candidates=("citrate",), mixture_kind="unresolved")
    cons = Consensus(top=top, calls=[accepted]).refresh([accepted], finalized=True)
    assert cons.top == top and cons.candidates == ("citrate",)
    assert cons.mixture_kind == "" and cons.agreement == 0
    assert cons.support_sources == () and cons.top_display == (top or "unresolved: citrate")
    result = LR.accepted_ligands(json.loads(json.dumps(asdict(cons))))
    assert result[0]["ligand"] == "citrate" and result[0]["final_role"] == "retained_alternative"
    represented = Consensus.from_calls([Call("exact", "citrate", 0.6,
                                             candidates=("citrate",))])
    assert represented.candidates == (), "a named singleton headline still needs no duplicate candidate set"


def test_singleton_policy_alternative_survives_a_guard_json_and_report():
    raw_ssn = Call("ssn_cluster", "rejected organic", 0.9, "non-metal")
    raw_coord = Call("coordination", None, 0, "metal", candidates=("Zn2+",),
                     evidence={"implied_ions": {"Zn2+": ["synthetic motif"]},
                               "motif_ion": "Zn2+", "motifs": ["synthetic motif"]})
    cons = IND._finalize([raw_ssn, raw_coord], True, sensor_class="metal")
    recorded = json.loads(json.dumps(asdict(cons)))
    assert recorded["top"] == recorded["top_display"] == "divalent metal (ion unresolved)"
    assert recorded["candidates"] == ["Zn2+"] and recorded["mixture_kind"] == ""
    assert recorded["support_sources"] == [] and recorded["agreement"] == 0
    assert recorded["calls"] == json.loads(json.dumps([asdict(raw_ssn), asdict(raw_coord)]))
    accepted = LR.accepted_ligands(recorded)
    assert accepted[0]["ligand"] == "Zn2+" and accepted[0]["final_role"] == "retained_alternative"
    html = render_report.build_html({"dossier": {"tf_id": "synthetic", "inducers": recorded}}, {})
    assert "Retained alternative:</b> Zn2+" in html
    assert "not an exact headline assignment" in html
    assert "UNRESOLVED between" not in html and "Candidate set:</b> rejected organic" not in html


def test_motif_only_report_shows_raw_abstention_and_distinct_accepted_interpretation():
    raw = Call("coordination", None, 0, "metal", candidates=("Zn2+",),
               evidence={"implied_ions": {"Zn2+": ["synthetic motif"]},
                         "motif_ion": "Zn2+", "motifs": ["synthetic motif"]})
    cons = IND._finalize([raw], True, sensor_class="metal")
    recorded = json.loads(json.dumps(asdict(cons)))
    assert recorded["calls"][0]["ligand"] is None and recorded["calls"][0]["confidence"] == 0
    assert recorded["accepted_calls"][0]["ligand"] == "Zn2+"
    assert recorded["accepted_calls"][0]["confidence"] == 0.35
    html = render_report.build_html({"dossier": {"tf_id": "synthetic", "inducers": recorded}}, {})
    before, accepted_section = html.split("<h3>Accepted policy interpretation of recorded evidence</h3>", 1)
    assert "abstained" in before and "raw abstention" in before
    assert "accepted policy interpretation shown separately" in before
    assert "accepted call" in accepted_section and "Zn2+" in accepted_section and "0.35" in accepted_section
    assert "synthetic motif" in accepted_section
