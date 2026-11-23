"""EXP066: the operon-chemistry route -- promoter gate, two-bound operon walk, macromolecule filter.

These pin the behaviour the rename introduced. The scoring formula is Ligify's and is covered by
`rank.py`'s own self-test; here we test the NEIGHBOURHOOD and the ligand filter, which are the parts
EXP066 changed.
"""
from __future__ import annotations

import pytest

from predictor.effector import operon_chemistry as OC
from predictor.effector.ligand import _is_macromolecule, candidate_ligands
from predictor.effector import inducer as IND
from predictor.schema import InducerCall
from predictor.signals.motif_rescan import Gene, GenomeContext


def _ctx(genes):
    return GenomeContext("SYN.1", "N" * 40000, genes), sorted(genes, key=lambda g: g.start)


def test_divergent_partner_plus_strand():
    # TF on '+' diverges from its LEFT neighbour, which must be on '-' and within the promoter gap.
    tf = Gene(2000, 2400, "+", "tf", protein_id="WP_tf", product="regulator")
    part = Gene(1500, 1900, "-", "p", protein_id="WP_p", product="dehydrogenase")   # gap 100
    _, s = _ctx([part, tf])
    assert OC.divergent_partner(s, tf).name == "p"

    # a co-directional left neighbour is NOT a divergent partner
    codir = Gene(1500, 1900, "+", "c", protein_id="WP_c", product="enzyme")
    _, s2 = _ctx([codir, tf])
    assert OC.divergent_partner(s2, tf) is None


def test_divergent_partner_minus_strand_and_promoter_gate():
    # TF on '-' diverges from its RIGHT neighbour, which must be on '+'.
    tf = Gene(2000, 2400, "-", "tf", protein_id="WP_tf", product="regulator")
    near = Gene(2500, 3000, "+", "p", protein_id="WP_p", product="enzyme")           # gap 100
    _, s = _ctx([tf, near])
    assert OC.divergent_partner(s, tf).name == "p"

    # beyond MAX_PROMOTER_BP -> the hard gate returns no partner
    far = Gene(4000, 4500, "+", "pf", protein_id="WP_pf", product="enzyme")          # gap 1600 > 1000
    _, s2 = _ctx([tf, far])
    assert OC.divergent_partner(s2, tf) is None


def test_operon_gap_stop_is_wider_than_regulon():
    # A 250 bp gap (between regulon's 150 and operon_chem's 300) is kept; a >300 gap stops the walk.
    partner = Gene(1000, 1500, "+", "g1", protein_id="WP_1", product="enzyme")
    g2 = Gene(1750, 2200, "+", "g2", protein_id="WP_2", product="enzyme")            # gap 250: kept
    g3 = Gene(2600, 3000, "+", "g3", protein_id="WP_3", product="enzyme")            # gap 400: stops
    _, s = _ctx([partner, g2, g3])
    names = [g.name for g in OC.operon_from_partner(s, partner)]
    assert names == ["g1", "g2"], names


def test_operon_span_cap_and_noncoding_termination():
    partner = Gene(1000, 1500, "+", "g1", protein_id="WP_1", product="enzyme")
    g2 = Gene(1600, 2000, "+", "g2", protein_id="WP_2", product="enzyme")            # gap 100: kept
    trna = Gene(2100, 2180, "+", "t", protein_id="", product="tRNA-Gly")             # ncRNA: terminates
    g3 = Gene(2200, 2600, "+", "g3", protein_id="WP_3", product="enzyme")
    _, s = _ctx([partner, g2, trna, g3])
    assert [g.name for g in OC.operon_from_partner(s, partner)] == ["g1", "g2"]

    # a same-strand gene past the 8 kb reach is excluded even with small gaps leading up to it
    chain = [Gene(1000, 1500, "+", "g1", protein_id="WP_1", product="e")]
    start = 1600
    for k in range(2, 40):                       # dense same-strand run marching past 8 kb
        chain.append(Gene(start, start + 400, "+", f"g{k}", protein_id=f"WP_{k}", product="e"))
        start += 500
    _, s2 = _ctx(chain)
    op = OC.operon_from_partner(s2, chain[0])
    assert all(abs(g.start - chain[0].start) <= OC.OPERON_CHEM_MAX_SPAN_BP for g in op)
    assert len(op) < len(chain), "the 8 kb reach must bound a runaway dense run"


def test_promoter_gate_blocks_the_divergent_route_but_not_the_tfs_own_run():
    # EXP071 changed this contract: with no divergent partner the DIVERGENT route is still closed, but
    # the regulator's own co-directional run (here the upstream gene 100 bp away) is now read.
    tf = Gene(2000, 2400, "+", "tf", protein_id="WP_tf", product="regulator")
    codir = Gene(1500, 1900, "+", "c", protein_id="WP_c", product="enzyme")
    ctx, _ = _ctx([codir, tf])
    fetched = []

    def _fetch(pid):
        fetched.append(pid)
        return None
    assert OC.effectors_from_context(ctx, tf, entry_fetcher=_fetch) == []
    assert fetched == ["WP_c"]


def test_no_divergent_partner_and_no_run_fetches_nothing():
    tf = Gene(2000, 2400, "+", "tf", protein_id="WP_tf", product="regulator")
    far = Gene(500, 900, "+", "far", protein_id="WP_far", product="enzyme")        # 1.1 kb gap
    opp = Gene(2700, 3100, "-", "opp", protein_id="WP_opp", product="enzyme")       # downstream, opposite
    ctx, _ = _ctx([far, tf, opp])
    fetched = []
    assert OC.effectors_from_context(ctx, tf, entry_fetcher=lambda p: fetched.append(p)) == []
    assert fetched == [], "neither route exists, so no enzyme may be fetched"


def test_macromolecule_filter_rejects_polymers_keeps_metabolites():
    for bad in ["tRNA(Gly)", "RNA(n)", "protein L-histidine", "N(6)-acetyl-L-lysyl-[protein]",
                "a [protein]-L-cysteine", "mRNA", "rRNA", "ncRNA", "polynucleotide"]:
        assert _is_macromolecule(bad), f"should be rejected: {bad}"
    for good in ["arsenite", "Cu(+)", "D-gluconate", "alpha,alpha-trehalose", "pyridoxal",
                 "(2-aminoethyl)phosphonate", "long-chain fatty acid", "arsenic triglutathione",
                 "molybdate", "beta-L-arabinopyranose", "L-ribulose", "adenosine", "salicylate"]:
        assert not _is_macromolecule(good), f"should be kept: {good}"


def test_macromolecule_filter_applies_in_candidate_ligands():
    entry = {"comments": [{"commentType": "CATALYTIC ACTIVITY", "reaction": {
        "name": "L-histidine + tRNA(His) = protein L-histidine + AMP",
        "reactionCrossReferences": [{"database": "Rhea", "id": "RHEA:00001"}]}}]}
    names = [l["name"] for l in candidate_ligands(entry)]
    assert "L-histidine" in names                      # the free amino acid survives
    assert "tRNA(His)" not in names and "protein L-histidine" not in names


@pytest.mark.parametrize("currency", ["H(+)(in)", "H(+)(out)", "H(+) ( IN )", "proton(out)",
                                      "ATP(in)", "water (Out)", "NAD(+)(in)", "H(+)"])
def test_compartment_tagged_currency_is_filtered(currency):
    entry = {"comments": [{"commentType": "CATALYTIC ACTIVITY", "reaction": {
        "name": f"{currency} + pyruvate = lactate"}}]}
    assert [c["name"] for c in candidate_ligands(entry)] == ["pyruvate", "lactate"]


def test_compartment_filter_preserves_real_metabolite_names_and_chemical_parentheses():
    entry = {"comments": [{"commentType": "CATALYTIC ACTIVITY", "reaction": {
        "name": "molybdate(out) + lactate(in) = (R)-lactate + ATP(analogue)"}}]}
    assert [c["name"] for c in candidate_ligands(entry)] == [
        "molybdate(out)", "lactate(in)", "(R)-lactate", "ATP(analogue)"]


def _rejected_operon_calls(source, gate):
    # The positive-site branch rejects a suppressed vote; the silent-site metal-family branch rejects
    # an unsuppressed organic vote. Both fixtures make chemistry the initial highest-confidence source.
    return [InducerCall(source, "salicylate", 0.8, "",
                        evidence={"weight": 0.1 if gate is True else 1.0})]


@pytest.mark.parametrize("source", ["operon_chem", "ligify"])
@pytest.mark.parametrize("gate", [True, None])
def test_rejected_operon_headline_falls_back_to_surviving_organic(source, gate):
    # EXP072: the fallback is reached only when the rejected operon call was the headline, i.e. when the
    # alternative is not in a higher-precedence tier. An unlisted source shares the local tier.
    calls = _rejected_operon_calls(source, gate)
    alternative = InducerCall("curated_note", "citrate", 0.3, "non-metal")
    calls.append(alternative)
    got = IND._finalize(calls, gate, sensor_class="metal")
    assert got.top == "citrate"
    assert got.calls == calls and got.calls[0].confidence == 0.8
    assert alternative.confidence == 0.3
    assert any("surviving organic" in n for n in got.notes)


@pytest.mark.parametrize("source", ["operon_chem", "ligify"])
@pytest.mark.parametrize("gate", [True, None])
def test_published_organic_outranks_a_louder_operon_call_by_precedence(source, gate):
    # Before EXP072 this headline was reached through the guard's fallback (the operon call was the
    # initial top on confidence and was then rejected). It is now the headline from the start, by source.
    calls = _rejected_operon_calls(source, gate)
    alternative = InducerCall("ligify_db", "citrate", 0.3, "")
    calls.append(alternative)
    got = IND._finalize(calls, gate, sensor_class="metal")
    assert got.top == "citrate"
    assert got.calls == calls and got.calls[0].confidence == 0.8 and alternative.confidence == 0.3
    assert any(n.startswith("source precedence") for n in got.notes)
    assert not any("surviving organic" in n for n in got.notes)


@pytest.mark.parametrize("source", ["operon_chem", "ligify"])
@pytest.mark.parametrize("gate", [True, None])
def test_independent_named_metal_retains_priority_over_organic_fallback(source, gate):
    calls = _rejected_operon_calls(source, gate) + [
        InducerCall("ligify_db", "glycolate", 0.3, ""),
        InducerCall("ssn_cluster", "Zn2+", 0.2, "metal")]
    assert IND._finalize(calls, gate, sensor_class="metal").top == "Zn2+"


@pytest.mark.parametrize("gate", [True, None])
@pytest.mark.parametrize("alternative", [
    InducerCall("ligify_db", "glycolate", 0.0, ""),
    InducerCall("ssn_cluster", "peroxide", 0.3, "redox"),
    InducerCall("operon_chem", "glycolate", 0.3, "non-metal"),
    InducerCall("ligify", "glycolate", 0.3, "non-metal"),
    InducerCall("untyped_source", "glycolate", 0.3, ""),
    InducerCall("ligify_db", None, 0.3, ""),
    InducerCall("ligify_db", "unknown", 0.3, ""),
    InducerCall("ligify_db", "non-metal effector (undetermined)", 0.3, ""),
    InducerCall("curated_note", "H2O2", 0.3, ""),
    InducerCall("curated_note", "Zn2+", 0.3, ""),
    InducerCall("ligify_db", "glycolate/lactate", 0.3, "", mixture_kind="unresolved"),
])
def test_guard_does_not_promote_ineligible_organic_fallback(gate, alternative):
    calls = _rejected_operon_calls("operon_chem", gate) + [alternative]
    assert IND._finalize(calls, gate, sensor_class="metal").top == "divalent metal (ion unresolved)"


@pytest.mark.parametrize("name", ["H2O2", "Zn2+"])
def test_fallback_itself_never_returns_a_published_metal_or_redox_name(name):
    # EXP072 moved the two `ligify_db` rows of the table above to an unlisted source: a NAMED published
    # call now outranks a louder operon call by source precedence, so it is the headline before any guard
    # runs (previously only when it was also the louder one). The fallback's own rule is unchanged.
    published = InducerCall("ligify_db", name, 0.3, "")
    assert IND._surviving_organic_call([published]) is None
    got = IND._finalize(_rejected_operon_calls("operon_chem", None) + [published], None, sensor_class="metal")
    assert got.top == name and any(n.startswith("source precedence") for n in got.notes)


@pytest.mark.parametrize("role", ["non-metal", "multidrug"])
def test_surviving_organic_fallback_uses_strongest_eligible_source_without_db_precedence(role):
    calls = _rejected_operon_calls("operon_chem", None) + [
        InducerCall("ligify_db", "glycolate", 0.3, ""),
        InducerCall("ssn_cluster", "salicylate", 0.5, role)]
    assert IND._finalize(calls, None, sensor_class="metal").top == "salicylate"


def test_surviving_organic_fallback_does_not_change_ssn_metal_boundary_policy():
    calls = [InducerCall("ssn_cluster", "formaldehyde", 1.0, "non-metal"),
             InducerCall("ligify_db", "glycolate", 0.3, "")]
    got = IND._finalize(calls, True, sensor_class="metal")
    assert got.top == "divalent metal (ion unresolved)"
    assert any("metal-boundary" in n for n in got.notes)
