"""Locked phase 3 policy tested on generic synthetic evidence only."""
from dataclasses import asdict
import json
from types import SimpleNamespace as NS

import pytest

from predictor.annotate.family_kb import Anchor, AnchorSet
from predictor import genome_run
from predictor.effector import chemical_identity as CI, inducer as IND, inducer_vocab as V
from predictor.effector import ligand as L, ligify, ligify_db as DB, operon_chemistry as OC
from predictor.report import figure_data as FD, ligand_regulon as LR, render_report, run_page, run_report
from predictor.schema import InducerCall as Call, InducerConsensus as Consensus
from predictor.signals.motif_rescan import Gene, GenomeContext


def entry(equation):
    return {"comments": [{"commentType": "CATALYTIC ACTIVITY", "reaction": {
        "name": equation, "reactionCrossReferences": [
            {"database": "Rhea", "id": "RHEA:SYNTHETIC"},
            {"database": "ChEBI", "id": "CHEBI:SYNTHETIC_1"},
            {"database": "ChEBI", "id": "CHEBI:SYNTHETIC_2"}]}}]}


def context():
    tf = Gene(100, 400, "-", "tf", protein_id="SYN_TF", product="regulator")
    enzyme = Gene(450, 1000, "+", "enzyme", protein_id="SYN_ENZ", product="isomerase")
    return GenomeContext("SYNTHETIC.1", "N" * 5000, [tf, enzyme]), tf


def stub(monkeypatch, assignment=None, database=None):
    monkeypatch.setattr(IND.hsel, "assign_cluster", lambda *a: assignment or NS(cluster_id=None))
    monkeypatch.setattr(IND, "curated_member_ion", lambda *a, **k: None)
    monkeypatch.setattr(DB, "lookup", lambda **k: database)
    monkeypatch.setattr(IND.metal_site, "coordination_gate", lambda *a, **k: {"has_site": None})


def test_reaction_ids_are_scoped_and_occurrences_are_lossless():
    raw = entry("2 citrate(in) + citrate(out) = citrate(out) + (R)-lactate 3-phosphate")
    rx = L.reactions_from_entry(raw)[0]
    assert rx.chebi == rx.reaction_chebi == ["CHEBI:SYNTHETIC_1", "CHEBI:SYNTHETIC_2"]
    assert rx.source_record == raw["comments"][0]["reaction"]
    participants = L.candidate_ligands(raw)
    assert len(participants) == 4
    assert all(p["chebi"] == [] and p["participant_identity_status"] == "unresolved" for p in participants)
    assert all(p["reaction_chebi"] == rx.reaction_chebi for p in participants)
    assert [p["side"] for p in participants] == ["left", "left", "right", "right"]
    assert [p["coefficient"] for p in participants] == ["2", "1", "1", "1"]
    assert [p["compartment"] for p in participants] == ["in", "out", "out", ""]
    assert participants[0]["raw_token"] == "2 citrate(in)"
    assert participants[-1]["name"] == "(R)-lactate 3-phosphate"


@pytest.mark.parametrize("route", [OC, ligify])
def test_transport_grouping_keeps_all_occurrences_and_both_routes(route):
    ctx, tf = context()
    candidates = route.effectors_from_context(ctx, tf, entry_fetcher=lambda _: entry("2 citrate(in) = citrate(out) + citrate(out)"))
    assert len(candidates) == 1 and candidates[0].chebi == []
    assert len(candidates[0].occurrences) == 3
    assert [o["side"] for o in candidates[0].occurrences] == ["left", "right", "right"]


@pytest.mark.parametrize("first,second", [("citrate", "salicylate"), ("novel A", "novel Z")])
def test_ties_are_invariant_to_rename_and_reorder(monkeypatch, first, second):
    stub(monkeypatch)
    ctx, tf = context()
    results = []
    for equation in (f"{first} = {second}", f"{second} = {first}"):
        candidates = OC.effectors_from_context(ctx, tf, entry_fetcher=lambda _: entry(equation))
        cons = IND.infer_inducer("SYNTHETIC", "GntR", ligify_fn=lambda *a: list(reversed(candidates)), metalnet=False)
        source = next(c for c in cons.calls if c.source == "operon_chem")
        assert source.ligand is None and source.mixture_kind == "unresolved"
        assert cons.top == "non-metal effector (undetermined)" and cons.sensor_class == "organic"
        assert cons.identity_status == "unresolved" and cons.agreement == 0
        assert set(cons.candidates) == {first, second}
        assert len(source.evidence["candidates"]) == len(cons.chemistry_snapshot["candidates"]) == 2
        results.append((cons.top, cons.sensor_class, set(cons.candidates)))
    assert results[0] == results[1]


def test_lower_candidates_stay_raw_but_only_operative_tie_is_accepted(monkeypatch):
    stub(monkeypatch)
    candidates = [OC.EffectorCandidate("citrate", 80, "", distance=1),
                  OC.EffectorCandidate("salicylate", 80, "", distance=1),
                  OC.EffectorCandidate("formaldehyde", 80, "", distance=2),
                  OC.EffectorCandidate("unmapped other", 70, "", distance=1)]
    cons = IND.infer_inducer("SYNTHETIC", "GntR", ligify_fn=lambda *a: list(reversed(candidates)), metalnet=False)
    assert set(cons.candidates) == {"citrate", "salicylate"}
    source = next(c for c in cons.calls if c.source == "operon_chem")
    assert len(source.evidence["candidates"]) == len(cons.chemistry_snapshot["candidates"]) == 4


def test_published_collection_is_unresolved_with_full_details(monkeypatch):
    record = DB._to_record({"candidate_ligands": [{"name": "citrate", "smiles": "synthetic-A", "iupac": "A"},
                                                  {"name": "salicylate", "smiles": "synthetic-B", "iupac": "B"}],
                            "rank": {"rank": 80}}, "sequence")
    stub(monkeypatch, database=record)
    cons = IND.infer_inducer("SYNTHETIC", "GntR", metalnet=False)
    call = next(c for c in cons.calls if c.source == "ligify_db")
    assert call.ligand is None and call.mixture_kind == "unresolved"
    assert call.evidence["ligand_detail"] == list(record.ligand_detail)
    assert set(cons.candidates) == {"citrate", "salicylate"}
    assert cons.sensor_class == "organic" and cons.identity_status == "unresolved"


@pytest.mark.parametrize("inducers,grades,expected", [
    (["Ni2+/Co2+"], ["literature"], "co_sensed"),
    (["Ni2+", "Co2+"], ["literature", "literature"], "unresolved"),
    (["Ni2+/Co2+"], ["omics"], "unresolved"),
    (["Ni2+", "Ni3+"], ["literature", "literature"], "unresolved"),
])
def test_cosensing_requires_literature_grade(monkeypatch, inducers, grades, expected):
    anchors = tuple(Anchor(f"SYN{i}", inducer=name, evidence_grade=grade, evidence_level="experimental")
                    for i, (name, grade) in enumerate(zip(inducers, grades, strict=True)))
    aset = AnchorSet("synthetic", n_anchors=len(anchors), anchors=anchors)
    assignment = NS(cluster_id="synthetic", support=0.9, top_identity=0.8,
                    info=NS(inducer=inducers[0], role="metal", representative="SYN", anchor_set=aset))
    stub(monkeypatch, assignment=assignment)
    cons = IND.infer_inducer("SYNTHETIC", "test", metalnet=False)
    call = next(c for c in cons.calls if c.source == "ssn_cluster")
    assert call.mixture_kind == expected
    assert call.ligand is None if expected == "unresolved" else call.ligand == inducers[0]
    assert cons.sensor_class == "metal"
    if "Ni3+" in inducers:
        assert set(cons.candidates) == {"Ni2+", "Ni3+"}


def test_blank_label_does_not_invent_legacy_default_charge_from_anchor(monkeypatch):
    anchor = Anchor("SYN", inducer="Ni3+", evidence_grade="literature")
    aset = AnchorSet("synthetic", n_anchors=1, anchors=(anchor,))
    assignment = NS(cluster_id="synthetic", support=0.9, top_identity=0.8,
                    info=NS(inducer=None, role="metal", representative="SYN", anchor_set=aset))
    stub(monkeypatch, assignment=assignment)
    cons = IND.infer_inducer("SYNTHETIC", "test", metalnet=False)
    assert cons.top == "Ni3+" and "Ni2+" not in cons.candidates


@pytest.mark.parametrize("name,cls,status", [
    ("Ni3+", "metal", "resolved"), ("Ni(III)", "metal", "resolved"),
    ("Cu(2+)", "metal", "resolved"), ("citrate", "organic", "resolved"),
    ("formaldehyde", "organic", "resolved"), ("H2O2", "redox", "resolved"),
    ("persulfide", "redox", "resolved"), ("unfamiliar molecule", None, "unmapped"),
    ("nitrogen", None, "unmapped"), ("Sodium ion", None, "unmapped"),
    ("steroid-CoA", None, "unmapped"), ("ATP", None, "excluded"),
    ("Na(+)", None, "excluded"), ("Na+", None, "excluded"), ("Mg(2+)", None, "excluded"),
    ("tRNA(Gly)", None, "excluded"), ("AH2", None, "generic"),
    ("AH2(out)", None, "generic"), ("AH2 derivative", None, "unmapped"),
])
def test_no_default_organic_or_broad_substring_classifier(name, cls, status):
    assert CI.chemistry(name)[:2] == (cls, status)


@pytest.mark.parametrize("name", ["unfamiliar molecule", "ATP", "Na(+)", "Mg(2+)", "tRNA(Gly)", "AH2"])
def test_ineligible_direct_calls_do_not_gain_exact_headlines(name):
    raw = Call("synthetic", name, 0.8)
    cons = Consensus.from_calls([raw])
    assert cons.top is None and cons.calls == cons.accepted_calls == [raw]
    assert cons.candidates == (name,) and cons.agreement == 0
    assert cons.sensor_class is None
    finalized = IND._finalize([raw], None)
    assert finalized.top is None and finalized.sensor_class is None


def test_charge_stereo_phosphate_and_one_terminal_transport_tag():
    assert CI.identity_key("Ni(III)") == CI.identity_key("Ni3+") != CI.identity_key("Ni2+")
    assert CI.identity_key("citrate ( IN ) ") == CI.identity_key("citrate(out)") == "citrate"
    assert CI.transport_name("citrate(in)(out)") == "citrate(in)"
    assert CI.identity_key("(R)-lactate") != CI.identity_key("(S)-lactate")
    assert CI.identity_key("D-ribulose") != CI.identity_key("D-ribulose 5-phosphate")
    assert CI.identity_key("ATP(analogue)") != CI.identity_key("ATP")
    assert V.normalise("Ni3+") == "Ni2+" and V.inducer_class_of("Sodium ion") == "organic"


def test_exact_agreement_and_display_preserve_distinct_ion_states():
    cons = Consensus.from_calls([Call("a", "Ni2+", 0.8, "metal"), Call("b", "Ni3+", 0.7, "metal")])
    assert cons.top == "Ni2+" and cons.agreement == 0.5
    assert cons.candidates == ("Ni2+", "Ni3+")
    ambiguous = IND._finalize([Call("set", None, 0.7, "metal", candidates=("Ni2+", "Ni(III)"),
                                   mixture_kind="unresolved")], None, sensor_class="metal")
    saved = json.loads(json.dumps(asdict(ambiguous)))
    assert "Ni2+" in FD.inducer_display(saved) and "Ni(III)" in FD.inducer_display(saved)
    saved.pop("top_display")
    assert "Ni(III)" in FD.inducer_display(saved)
    assert len(LR.accepted_ligands(saved)) == 2
    compatible = Consensus.from_calls([Call("a", "Ni3+", 0.8, "metal"),
                                      Call("set", None, 0.7, "metal", candidates=("Ni2+", "Ni(III)"),
                                           mixture_kind="unresolved")])
    assert compatible.agreement == compatible.compatibility == 0.5


@pytest.mark.parametrize("cls", ["organic", "metal", "redox"])
def test_positive_set_only_preserves_class_provenance_and_support(cls):
    raw = Call("synthetic", None, 0.5, candidates=("alternative A", "alternative B"),
               mixture_kind="unresolved", identity_status="unresolved", sensor_class=cls,
               class_status="supported", class_evidence={"basis": "explicit_synthetic_source"})
    cons = IND._finalize([raw], None)
    saved = json.loads(json.dumps(asdict(cons)))
    assert cons.sensor_class == cls and cons.class_status == "supported"
    assert cons.top is None and cons.identity_status == "unresolved"
    assert cons.candidates == ("alternative A", "alternative B") and cons.agreement == 0
    assert cons.class_evidence["records"][0]["basis"] == "explicit_synthetic_source"
    assert FD.inducer_display(saved) == "unresolved: alternative A or alternative B"
    assert CI.class_from_record(saved) == (cls, "explicit")
    assert CI.class_from_record({"top": "unmapped legacy"}) == ("organic", "legacy_label_fallback")
    assert CI.class_from_record({"top": "citrate", "sensor_class": None}) == (None, "explicit")
    assert len(LR.accepted_ligands(saved)) == 2
    assert "alternative A" in run_page._inducer_cell({"abstention": True}, {"inducers": saved})
    html = render_report.build_html({"dossier": {"tf_id": "synthetic", "inducers": saved}}, {})
    assert "alternative A" in html and f"sensor class {cls}" in html


def test_conflicting_class_sources_are_visible_and_unknown_is_not_a_class():
    cons = Consensus.from_calls([Call("a", None, 0.7, sensor_class="metal", candidates=("A",)),
                                Call("b", None, 0.6, sensor_class="organic", candidates=("B",))])
    assert cons.sensor_class is None and cons.class_status == "conflicting"
    assert {r["class"] for r in cons.class_evidence["records"]} == {"metal", "organic"}
    unknown = Consensus.from_calls([Call("a", "arbitrary string", 0.6)])
    assert unknown.top is None and unknown.sensor_class is None and unknown.class_status == "unknown"


@pytest.mark.parametrize("role", ["", "metal"])
def test_cross_class_candidate_collection_keeps_internal_conflict_visible(role):
    cons = Consensus.from_calls([Call("synthetic", None, 0.8, role,
                                      candidates=("Cu+", "persulfide"), mixture_kind="unresolved")])
    assert cons.top is None and cons.sensor_class is None and cons.class_status == "conflicting"
    assert {r["class"] for r in cons.class_evidence["records"]} == {"metal", "redox"}


def test_unmapped_name_retains_explicit_organic_source_class_without_exact_identity():
    raw = Call("ssn_cluster", "unmapped synthetic compound", 0.8, "non-metal")
    cons = IND._finalize([raw], None)
    assert cons.top is None and cons.identity_status == "unmapped"
    assert cons.sensor_class == "organic" and cons.class_status == "supported"
    assert cons.candidates == (raw.ligand,) and cons.calls == [raw]


def test_ah2_stays_raw_without_exact_headline_or_broad_exclusion(monkeypatch):
    stub(monkeypatch)
    ctx, tf = context()
    candidates = OC.effectors_from_context(ctx, tf, entry_fetcher=lambda _: entry("AH2 = AH2(out)"))
    assert len(candidates) == 1 and not candidates[0].headline_eligible
    assert len(candidates[0].occurrences) == 2
    cons = IND.infer_inducer("SYNTHETIC", "test", ligify_fn=lambda *a: candidates, metalnet=False)
    assert cons.top is None and cons.chemistry_snapshot["candidates"][0]["ligand"] == "AH2"
    assert L.candidate_ligands(entry("AH2 derivative = salicylate"))[0]["name"] == "AH2 derivative"


def test_structural_cofactor_and_neighbor_substrates_do_not_create_sensor_class(monkeypatch):
    stub(monkeypatch)
    candidates = [OC.EffectorCandidate("citrate", 80, "")]
    cons = IND.infer_inducer("SYNTHETIC", "test", ligify_fn=lambda *a: candidates, metalnet=False)
    assert cons.sensor_class is None and cons.class_status == "unknown"
    caveated = IND._finalize([Call("coordination", None, 0, evidence={"family": "GntR"}),
                              Call("metalnet", None, 0, "metal")], None)
    assert caveated.sensor_class == "organic" and caveated.top == "non-metal effector (undetermined)"
    assert all(r["class"] != "metal" for r in caveated.class_evidence["records"])


def test_old_and_new_schema_fields_remain_additive():
    old = Call(**{"source": "legacy", "ligand": "unmapped old name", "confidence": 0.7})
    assert old.identity_status == "" and old.sensor_class is None
    new = Call("a", "Ni3+", 0.8, "metal", identity_status="resolved", sensor_class="metal")
    loaded = Call(**json.loads(json.dumps(asdict(new))))
    assert loaded == new and CI.identity_key(loaded.ligand) != CI.identity_key("Ni2+")


def test_saved_run_rows_use_explicit_class_with_legacy_fallback(tmp_path, monkeypatch):
    for name, inducer in (("new", {"top": None, "candidates": ["A", "B"], "identity_status": "unresolved",
                                   "sensor_class": "organic", "finalized": True}),
                          ("legacy", {"top": "unmapped old name"})):
        path = tmp_path / name
        path.mkdir()
        (path / "dossier.json").write_text(json.dumps({"tf_id": name, "family": "test", "inducers": inducer}))
    rows = {r["run_name"]: r for r in run_report._run_rows(tmp_path)[0]}
    assert rows["new"]["inducer_class"] == "organic"
    assert rows["new"]["inducer_class_provenance"] == "explicit"
    assert rows["legacy"]["inducer_class_provenance"] == "legacy_label_fallback"
    path = tmp_path / "jobs" / "new"
    path.mkdir(parents=True)
    (path / "dossier.json").write_text((tmp_path / "new" / "dossier.json").read_text())
    monkeypatch.setattr(genome_run, "bundle_complete", lambda _: True)
    summary = genome_run.summarise([{"run_name": "new", "locus_tag": "SYN", "gene": "", "product": "", "family": "test"}],
                                  {"new": "ok"}, out_dir=tmp_path)
    header, data = [line.split("\t") for line in summary.read_text().splitlines()]
    assert data[header.index("inducer_class_provenance")] == "explicit"
    assert data[header.index("inducer_class")] == "organic"


@pytest.mark.parametrize("first,second", [("MoO4-", "MoO42-"), ("AsO4-", "AsO43-"),
                                           ("MoO4", "MoO42-"), ("CO", "Co")])
def test_explicit_formula_charge_is_preserved_before_legacy_aliases(first, second):
    assert CI.identity_key(first) != CI.identity_key(second)
    raw = [Call("a", first, 0.8, "metal"), Call("b", second, 0.7, "metal")]
    cons = Consensus.from_calls(raw)
    assert cons.agreement <= 0.5 and set(cons.candidates) == {first, second}
    saved = asdict(IND._finalize([Call("set", None, 0.8, "metal", candidates=(first, second),
                                      mixture_kind="unresolved")], None, sensor_class="metal"))
    assert first in FD.inducer_display(saved) and second in FD.inducer_display(saved)
    assert len(LR.accepted_ligands(saved)) == 2


@pytest.mark.parametrize("status", ["unresolved", "unmapped", "unknown", "generic", "excluded", "class_only"])
def test_explicit_nonidentified_status_is_never_an_exact_vote(status):
    raw = Call("uncertain", "citrate", 0.8, identity_status=status)
    only = Consensus.from_calls([raw])
    assert only.top is None and only.agreement == 0 and only.calls == [raw]
    with_exact = Consensus.from_calls([raw, Call("exact", "citrate", 0.7)])
    assert with_exact.top == "citrate" and with_exact.support_sources == ("exact",)
    assert with_exact.compatible_sources == ("uncertain",) and with_exact.agreement == 0.5
    assert raw.identity_status == status


@pytest.mark.parametrize("status,ligand", [("unresolved", "Ni2+"), ("unmapped", "unmapped metal proposal")])
@pytest.mark.parametrize("gate,sensor_class", [(True, None), (None, "metal")])
def test_guard_fallback_cannot_promote_an_unidentified_metal(status, ligand, gate, sensor_class):
    raw = Call("ssn_cluster", ligand, 0.8, "metal", candidates=(ligand, "Ni3+"),
               mixture_kind="unresolved" if status == "unresolved" else "", identity_status=status)
    local = Call("operon_chem", "citrate", 0.9, evidence={"weight": 0.1})
    cons = IND._finalize([raw, local], gate, sensor_class=sensor_class)
    assert cons.top == "divalent metal (ion unresolved)" and cons.agreement == 0
    assert cons.sensor_class == "metal" and set(cons.candidates) == {ligand, "Ni3+"}
    assert cons.calls == [raw, local] and raw in cons.accepted_calls


def test_other_named_guard_paths_use_the_same_exact_rule():
    uncertain_organic = Call("ssn_cluster", "citrate", 0.8, "non-metal", identity_status="unresolved",
                             candidates=("citrate", "salicylate"))
    local_metal = Call("operon_chem", "Zn2+", 0.9, "metal")
    organic = IND._finalize([uncertain_organic, local_metal], None, sensor_class="organic")
    assert organic.top == "non-metal effector (undetermined)" and organic.agreement == 0
    assert set(organic.candidates) == {"citrate", "salicylate"}
    uncertain_metal = Call("other", "Ni2+", 0.8, "metal", identity_status="unresolved",
                           candidates=("Ni2+", "Ni3+"))
    metal = IND._finalize([Call("ssn_cluster", "citrate", 0.9, "non-metal"), uncertain_metal], True)
    assert metal.top == "divalent metal (ion unresolved)" and metal.agreement == 0
    assert set(metal.candidates) == {"Ni2+", "Ni3+"}


def test_motif_promotion_respects_explicit_identity_uncertainty():
    raw = Call("coordination", None, 0, "metal", identity_status="unresolved",
               candidates=("Ni2+", "Ni3+"), evidence={"motif_ion": "Ni2+", "implied_ions": {"Ni2+": ["synthetic"]}})
    cons = IND._finalize([raw], True)
    assert cons.top == "divalent metal (ion unresolved)" and cons.agreement == 0
    assert cons.calls == [raw] and raw.confidence == 0


@pytest.mark.parametrize("class_status", ["unknown", "conflicting"])
def test_supplied_class_uncertainty_survives_consensus(class_status):
    raw = Call("a", None, 0.8, sensor_class="organic", class_status=class_status,
               class_evidence={"basis": "explicit_source_report", "conflicts": ["metal", "organic"]},
               candidates=("A", "B"))
    cons = Consensus.from_calls([raw])
    assert cons.sensor_class is None and cons.class_status == class_status
    assert cons.class_evidence["source_uncertainty"][0]["conflicts"] == ["metal", "organic"]
    assert cons.calls == [raw] and raw.class_status == class_status


@pytest.mark.parametrize("invalid", ["unknown", "conflicting", "unmapped", "synthetic invalid"])
def test_canonical_class_boundary_keeps_invalid_values_only_as_audit(invalid):
    raw = Call("a", None, 0.8, sensor_class=invalid)
    cons = Consensus.from_calls([raw])
    assert cons.sensor_class is None and cons.class_status == "unknown"
    assert cons.class_evidence["source_uncertainty"][0]["reported_class"] == invalid
    policy = IND._finalize([], None, sensor_class=invalid)
    assert policy.sensor_class is None and policy.class_status == "unknown"
    assert policy.class_evidence["policy_input"] == invalid
    assert raw.sensor_class == invalid


def test_existing_canonical_no_redox_signal_remains_exact_and_cross_class_visible():
    assert CI.chemistry("NO")[:2] == ("redox", "resolved")
    cons = Consensus.from_calls([Call("a", "NO", 0.8)])
    assert cons.top == "NO" and cons.sensor_class == "redox" and cons.agreement == 1
    mixed = Consensus.from_calls([Call("set", None, 0.8, candidates=("Cu+", "NO"), mixture_kind="unresolved")])
    assert mixed.sensor_class is None and mixed.class_status == "conflicting"


def test_synonymous_singleton_anchor_does_not_claim_cosensing(monkeypatch):
    anchor = Anchor("SYN", inducer="Ni(II)", evidence_grade="literature")
    aset = AnchorSet("synthetic", n_anchors=1, anchors=(anchor,))
    assignment = NS(cluster_id="synthetic", support=0.9, top_identity=0.8,
                    info=NS(inducer="Ni2+", role="metal", representative="SYN", anchor_set=aset))
    stub(monkeypatch, assignment=assignment)
    cons = IND.infer_inducer("SYNTHETIC", "test", metalnet=False)
    call = next(c for c in cons.calls if c.source == "ssn_cluster")
    assert call.candidates == () and call.mixture_kind == "" and call.identity_status == "resolved"
    assert call.evidence["raw_components"] == ["Ni(II)", "Ni2+"]
    assert cons.candidates == () and cons.mixture_kind == "" and cons.identity_status == "resolved"


@pytest.mark.parametrize("formula,ion", [("CO2+", "Co2+"), ("CO(2+)", "Co(2+)"),
                                       ("CO(II)", "Co(II)"), ("CO2+(in)", "Co2+(in)")])
def test_multielement_formula_cannot_be_classified_as_monatomic_cobalt(formula, ion):
    assert CI.identity_key(formula) != CI.identity_key(ion)
    assert CI.chemistry(formula)[:2] == (None, "unmapped")
    assert CI.chemistry(ion)[:2] == ("metal", "resolved")
    raw = Call("formula", formula, 0.9)
    exact = Call("ion", ion, 0.8)
    assert not CI.exact_eligible(raw) and CI.exact_eligible(exact)
    cons = Consensus.from_calls([raw, exact])
    assert cons.top == ion and cons.support_sources == ("ion",) and cons.agreement == 0.5
    assert set(cons.candidates) == {formula, ion} and cons.sensor_class == "metal"
