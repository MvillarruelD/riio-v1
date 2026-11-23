"""EXP070 acceptance tests, written from EXP070_PREREG_20261109.md before the code (synthetic evidence only)."""
from types import SimpleNamespace as NS

import pytest

from predictor.annotate.family_kb import Anchor, AnchorSet
from predictor.effector import chemical_identity as CI, inducer as IND, ligify_db as DB

S = CI.subsume_bare_elements


# --- arm A: bare element symbols ---------------------------------------------------------------
@pytest.mark.parametrize("names,expected", [
    (("Fe", "Fe2+"), ("Fe2+",)),
    (("Ni", "Ni2+"), ("Ni2+",)),
    (("Ni(2+)(in)", "Ni"), ("Ni(2+)(in)",)),
    (("Fe", "Fe3+"), ("Fe3+",)),                      # no oxidation state is invented
    (("Fe2+", "Fe3+"), ("Fe2+", "Fe3+")),             # two explicit states stay two
    (("Fe", "Fe2+", "Fe3+"), ("Fe2+", "Fe3+")),
    (("Cu+", "Au(I)"), ("Cu+", "Au(I)")),             # different elements
    (("As(III)", "antimonite"), ("As(III)", "antimonite")),
    (("Co", "CO"), ("Co", "CO")),                     # CO is carbon monoxide, not an ion of Co
    (("CO", "Co2+"), ("CO", "Co2+")),
    (("Fe",), ("Fe",)),
])
def test_bare_symbol_is_subsumed_only_by_an_explicit_ion_of_the_same_element(names, expected):
    assert S(names) == expected


def test_subsumption_against_a_wider_context():
    assert S(("Fe",), context=("Fe", "Fe2+")) == ()
    assert S(("Fe", "Zn2+"), context=("Fe2+",)) == ("Zn2+",)


def test_co_and_carbon_monoxide_stay_distinct():
    assert CI.identity_key("CO") != CI.identity_key("Co")
    assert CI._monatomic_ion("CO") is None
    assert CI.chemistry("CO")[0] != "metal"
    assert CI.chemistry("Co")[:2] == ("metal", "resolved")


def test_lone_bare_symbol_is_metal_without_an_asserted_oxidation_state():
    cls, status, ev = CI.chemistry("Fe")
    assert (cls, status) == ("metal", "resolved")
    assert ev == {"basis": "bare_element_symbol", "element": "Fe"}
    assert CI.identity_key("Fe") != CI.identity_key("Fe2+")
    assert CI.identity_key("Fe") != CI.identity_key("Fe3+")


@pytest.mark.parametrize("name", ["C", "N", "S", "P", "co", "FE"])
def test_non_metal_or_wrong_case_symbols_stay_unmapped(name):
    assert CI.chemistry(name)[1] == "unmapped"


@pytest.mark.parametrize("name", ["fe", "ni"])
def test_lower_case_spellings_keep_their_pre_exp070_legacy_alias_route(name):
    # Pinned as found at 75eef4ab: these reach the legacy `_CANONICAL` alias, not the symbol rule.
    assert CI.chemistry(name) == ("metal", "resolved", {"basis": "existing_explicit_vocabulary"})


@pytest.mark.parametrize("name", ["K", "Mg"])
def test_bare_currency_elements_are_excluded_like_their_ions(name):
    assert CI.chemistry(name)[1] == "excluded"


def test_bare_na_is_still_read_as_the_na_abstention_token():
    # Pre-existing (75eef4ab): `inducer_vocab.is_abstention("Na")` is True ("NA" = not available).
    assert CI.chemistry("Na")[1] == "unknown"


def _stub(monkeypatch, assignment, database=None):
    monkeypatch.setattr(IND.hsel, "assign_cluster", lambda *a: assignment)
    monkeypatch.setattr(IND, "curated_member_ion", lambda *a, **k: None)
    monkeypatch.setattr(DB, "lookup", lambda **k: database)
    monkeypatch.setattr(IND.metal_site, "coordination_gate", lambda *a, **k: {"has_site": None})


def _cluster(label, *anchor_inducers, role="metal"):
    anchors = tuple(Anchor(f"SYN{i}", inducer=x, evidence_grade="literature")
                    for i, x in enumerate(anchor_inducers))
    aset = AnchorSet("synthetic", n_anchors=len(anchors), anchors=anchors)
    return NS(cluster_id="synthetic", support=0.95, top_identity=0.9,
              info=NS(inducer=label, role=role, representative="SYN", anchor_set=aset))


@pytest.mark.parametrize("label,anchor,ion", [("Fe2+", "Fe", "Fe2+"), ("Ni2+", "Ni", "Ni2+")])
def test_ssn_source_names_the_ion_when_anchors_write_the_bare_symbol(monkeypatch, label, anchor, ion):
    _stub(monkeypatch, _cluster(label, anchor, anchor))
    cons = IND.infer_inducer("SYNTHETIC", "test", metalnet=False)
    call = next(c for c in cons.calls if c.source == "ssn_cluster")
    assert call.ligand == ion and call.candidates == () and call.mixture_kind == ""
    assert call.evidence["raw_components"] == [anchor, anchor, label]      # raw spelling retained
    assert cons.top == ion and cons.identity_status == "resolved"


@pytest.mark.parametrize("label,anchors", [("Cu+", ("Au(I)",)), ("As3+", ("As(III)", "antimonite")),
                                           ("Fe2+", ("Fe3+",))])
def test_ssn_source_keeps_genuinely_different_species_unresolved(monkeypatch, label, anchors):
    _stub(monkeypatch, _cluster(label, *anchors))
    cons = IND.infer_inducer("SYNTHETIC", "test", metalnet=False)
    call = next(c for c in cons.calls if c.source == "ssn_cluster")
    assert call.ligand is None and call.mixture_kind == "unresolved" and len(call.candidates) == 2


# --- arm A2: a generic participant does not widen the operative tie -----------------------------
from predictor.effector import operon_chemistry as OC  # noqa: E402


def _operon_call(monkeypatch, candidates):
    _stub(monkeypatch, NS(cluster_id=None))
    cons = IND.infer_inducer("SYNTHETIC", "GntR", ligify_fn=lambda *a: candidates, metalnet=False)
    return cons, next(c for c in cons.calls if c.source == "operon_chem")


def test_ah2_is_not_a_member_of_the_operative_tie(monkeypatch):
    ah2 = OC.EffectorCandidate("AH2", 80, "", distance=1, headline_eligible=False)
    cands = [ah2, OC.EffectorCandidate("glycolate", 80, "", distance=1),
             OC.EffectorCandidate("glyoxylate", 80, "", distance=1)]
    cons, call = _operon_call(monkeypatch, cands)
    assert set(call.candidates) == {"glycolate", "glyoxylate"} and call.mixture_kind == "unresolved"
    assert "AH2" in [c["ligand"] for c in call.evidence["candidates"]]          # raw evidence kept
    assert "AH2" not in cons.candidates


def test_ah2_alone_never_becomes_a_headline(monkeypatch):
    cons, call = _operon_call(monkeypatch, [OC.EffectorCandidate("AH2", 80, "", headline_eligible=False)])
    assert call.ligand is None and call.candidates == () and cons.top != "AH2"


def test_single_eligible_name_beside_ah2_is_named_even_when_ah2_sorts_first(monkeypatch):
    cands = [OC.EffectorCandidate("AH2", 80, "", distance=1, headline_eligible=False),
             OC.EffectorCandidate("salicylate", 80, "", distance=1)]
    cons, call = _operon_call(monkeypatch, cands)
    assert call.ligand == "salicylate" and call.candidates == ("salicylate",)


# --- arm C: the coordination gate splits a cross-class clade set ---------------------------------
from predictor.schema import InducerCall as Call  # noqa: E402


def _cross(cands=("Co2+", "Formaldehyde")):
    return Call("ssn_cluster", None, 0.5, "non-metal", evidence={"cluster": "syn"},
                candidates=cands, mixture_kind="unresolved", identity_status="unresolved",
                sensor_class="organic")


@pytest.mark.parametrize("gate,top,cls", [(True, "Co2+", "metal"), (False, "Formaldehyde", "organic")])
def test_gate_splits_a_cross_class_set(gate, top, cls):
    raw = _cross()
    cons = IND._finalize([raw], gate)
    assert cons.top == top and cons.sensor_class == cls and cons.class_status == "supported"
    assert cons.calls == [raw] and raw.candidates == ("Co2+", "Formaldehyde")      # raw retained
    assert any("EXP070" in n for n in cons.notes)


def test_silent_gate_leaves_a_cross_class_set_alone():
    cons = IND._finalize([_cross()], None)
    assert cons.top is None and cons.class_status == "conflicting"


@pytest.mark.parametrize("cands", [("Cu+", "Au(I)"), ("As(III)", "antimonite"),
                                   ("Fe2+", "redox[2Fe-2S]"), ("Co2+", "unmapped organic name")])
def test_gate_does_not_pick_inside_same_class_redox_or_unclassified_sets(cands):
    cons = IND._finalize([_cross(cands)], True)
    assert not (cons.top in cands)


def test_gate_split_is_idempotent_on_refusion():
    first = IND._finalize([_cross()], True)
    again = IND._finalize(list(first.calls), True)
    assert (first.top, first.candidates, first.sensor_class) == (again.top, again.candidates, again.sensor_class)


# --- arm D: the inducer stage gets the TF, not the overlapping upstream gene ----------------------
def test_infer_inducer_stage_uses_best_overlap_not_start_containment(monkeypatch):
    from predictor import pipeline
    from predictor.signals.motif_rescan import Gene, GenomeContext
    upstream = Gene(1000, 2061, "+", "yegV", protein_id="UP.1")
    tf = Gene(2034, 2781, "-", "yegW", protein_id="TF.1")          # overlaps upstream by 27 bp
    ctx = GenomeContext("SYN.1", "N" * 5000, [upstream, tf])
    seen = {}

    def fake(seq, family, *, tf_gene=None, protein_acc=None, **kw):
        seen.update(gene=tf_gene, acc=protein_acc)
        return NS(top="x", agreement=0.0, coordination_gate=None, calls=[])
    monkeypatch.setattr(IND, "infer_inducer", fake)
    pipeline._infer_inducer("SEQ", "GntR", ctx, 2037, None, tf_end=2781, allow_ncbi=False, verbose=False)
    assert seen["gene"] is tf and seen["acc"] == "TF.1"
    pipeline._infer_inducer("SEQ", "GntR", ctx, 2037, None, allow_ncbi=False, verbose=False)
    assert seen["gene"] is upstream                     # legacy rule kept only for callers without tf_end
