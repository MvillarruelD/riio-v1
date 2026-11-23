"""EXP072 tests (EXP072_PREREG_20261109.md), synthetic evidence only.

Three contracts: the within-tie ordering of one enzyme's participants (implemented, measured, and
SHIPPED OFF -- `TIE_BREAK_SIGNALS == ()` is a recorded negative), the organic-only vote of a
near-identical published homologue, and headline precedence by SOURCE before confidence.
"""
from types import SimpleNamespace as NS

import pytest

from predictor import schema as S
from predictor.effector import chemical_identity as CI, inducer as IND, ligand as L, ligify_db as DB
from predictor.effector import operon_chemistry as OC
from predictor.schema import InducerCall, InducerConsensus
from predictor.signals.motif_rescan import Gene, GenomeContext


def entry(equation, direction=None):
    comment = {"commentType": "CATALYTIC ACTIVITY", "reaction": {
        "name": equation, "reactionCrossReferences": [{"database": "Rhea", "id": "RHEA:SYN"}]}}
    if direction:
        comment["physiologicalReactions"] = [{"directionType": d} for d in direction.split("+")]
    return {"comments": [comment]}


def ctx(genes):
    return GenomeContext("SYN.1", "N" * 40000, genes)


def occurrences(equation, direction=None, enzyme="enz"):
    return [{**o, "enzyme": enzyme} for o in L.candidate_ligands(entry(equation, direction))]


def named(cands, name):
    return next(c for c in cands if c.ligand.startswith(name))


@pytest.fixture
def signals_on(monkeypatch):
    """The three measured signals, applied -- the unshipped arm E5 of EXP072."""
    monkeypatch.setattr(OC, "TIE_BREAK_SIGNALS", OC.ALL_TIE_BREAK_SIGNALS)


# --- task 2: the signals --------------------------------------------------------------------------
def test_physiological_direction_is_read_only_when_it_is_unambiguous():
    assert L.reactions_from_entry(entry("a = b", "left-to-right"))[0].physiological_direction == "left-to-right"
    assert L.reactions_from_entry(entry("a = b", "right-to-left"))[0].physiological_direction == "right-to-left"
    assert L.reactions_from_entry(entry("a = b", "left-to-right+right-to-left"))[0].physiological_direction == ""
    assert L.reactions_from_entry(entry("a = b"))[0].physiological_direction == ""
    assert all(o["physiological_direction"] == "right-to-left"
               for o in L.candidate_ligands(entry("salicylate = catechol", "right-to-left")))


def test_import_follows_the_sourced_direction_not_the_written_order():
    written_import = [o for o in occurrences("lactate(out) = lactate(in)")]
    assert OC.tie_break_evidence(written_import)["import"] == 1
    written_export = [o for o in occurrences("lactate(in) = lactate(out)")]
    assert OC.tie_break_evidence(written_export)["import"] == 0
    # aroP in UniProt: written as an export, annotated right-to-left -- it is an importer
    sourced = OC.tie_break_evidence(occurrences("lactate(in) = lactate(out)", "right-to-left"))
    assert sourced["import"] == 1 and sourced["direction_basis"] == ["physiological"]
    assert OC.tie_break_evidence(written_export)["direction_basis"] == ["as_written"]


def test_substrate_is_the_effective_left_side_of_a_chemical_reaction_only():
    occ = occurrences("salicylate = catechol")
    by = {o["name"]: OC.tie_break_evidence([o]) for o in occ}
    assert by["salicylate"]["substrate"] == 1 and by["catechol"]["substrate"] == 0
    rev = {o["name"]: OC.tie_break_evidence([o]) for o in occurrences("salicylate = catechol", "right-to-left")}
    assert rev["salicylate"]["substrate"] == 0 and rev["catechol"]["substrate"] == 1
    # a transported species is not a substrate of a chemical reaction
    assert OC.tie_break_evidence(occurrences("lactate(out) = lactate(in)"))["substrate"] == 0


def test_recurrence_counts_distinct_enzymes():
    occ = occurrences("salicylate = catechol", enzyme="a") + occurrences("salicylate = gentisate", enzyme="b")
    sal = [o for o in occ if o["name"] == "salicylate"]
    assert OC.tie_break_evidence(sal)["recurrence"] == 2
    assert OC.tie_break_evidence(sal + sal)["recurrence"] == 2


def test_shipped_default_applies_no_signal():
    # EXP072: arm A1 failed its preregistered kill criterion, so the shipped prefix is empty.
    assert OC.TIE_BREAK_SIGNALS == ()
    assert OC.ALL_TIE_BREAK_SIGNALS == ("import", "substrate", "recurrence")


def test_key_is_lexicographic_in_the_declared_order(monkeypatch):
    ev = {"import": 0, "substrate": 1, "recurrence": 3}
    assert OC.tie_break_key(ev) == ()
    assert OC.tie_break_key(ev, ("substrate",)) == (-1,)
    monkeypatch.setattr(OC, "TIE_BREAK_SIGNALS", OC.ALL_TIE_BREAK_SIGNALS)
    assert OC.tie_break_key(ev) == (0, -1, -3)
    # import outranks any amount of the later signals
    assert OC.tie_break_key({"import": 1}, ("import", "substrate", "recurrence")) < (0, -1, -9)


# --- task 2: the operative tie --------------------------------------------------------------------
def test_operative_tie_orders_only_inside_one_rank_and_distance():
    a = NS(ligand="a", rank=100, distance=1, tiebreak=(0, -1), headline_eligible=True)
    b = NS(ligand="b", rank=100, distance=1, tiebreak=(0, 0), headline_eligible=True)
    far = NS(ligand="far", rank=90, distance=2, tiebreak=(-1, -1), headline_eligible=True)   # better key, worse rank
    top, rank_tie, tied = OC.operative_tie([a, b, far])
    assert top is a and rank_tie == [a, b] and tied == [a]
    assert OC.operative_tie([]) == (None, [], [])


def test_candidates_without_a_key_all_tie_and_a_generic_participant_never_wins():
    stub = [NS(ligand="x", rank=80, distance=1), NS(ligand="y", rank=80, distance=1)]      # injected stubs
    assert OC.operative_tie(stub)[2] == stub
    generic = NS(ligand="AH2", rank=100, distance=1, tiebreak=(-1,), headline_eligible=False)
    real = NS(ligand="glycolate", rank=100, distance=1, tiebreak=(0,), headline_eligible=True)
    top, rank_tie, tied = OC.operative_tie([generic, real])
    assert rank_tie == [generic, real] and tied == [real]


def _lactate_unit():
    # tf - dehydrogenase - permease, one co-transcribed run: lactate is imported AND consumed.
    tf = Gene(1000, 1400, "+", "tf", protein_id="TF", product="regulator")
    enz = Gene(1450, 2400, "+", "enz", protein_id="E", product="dehydrogenase")
    perm = Gene(2450, 3400, "+", "perm", protein_id="P", product="permease")
    entries = {"E": entry("lactate = pyruvate"), "P": entry("lactate(out) = lactate(in)")}
    return ctx([tf, enz, perm]), tf, entries.get


def test_imported_and_consumed_substrate_wins_its_enzyme_tie(signals_on):
    context, tf, fetch = _lactate_unit()
    cands = OC.effectors_from_context(context, tf, entry_fetcher=fetch)
    lactate, pyruvate = named(cands, "lactate"), named(cands, "pyruvate")
    assert (lactate.rank, lactate.distance) == (pyruvate.rank, pyruvate.distance)       # the Ligify tie
    assert lactate.tiebreak_evidence == {"import": 1, "substrate": 1, "recurrence": 2,
                                         "direction_basis": ["as_written"]}
    assert pyruvate.tiebreak_evidence["substrate"] == 0
    assert [c.ligand for c in OC.operative_tie(cands)[2]] == [lactate.ligand]
    assert len(cands) == 2                                                               # nothing deleted


def test_as_shipped_every_participant_of_the_best_enzyme_ties_and_evidence_is_still_recorded():
    context, tf, fetch = _lactate_unit()
    cands = OC.effectors_from_context(context, tf, entry_fetcher=fetch)
    assert all(c.tiebreak == () for c in cands)
    assert sorted(c.ligand for c in OC.operative_tie(cands)[2]) == sorted(c.ligand for c in cands)
    assert named(cands, "lactate").tiebreak_evidence["import"] == 1      # recorded, not applied


def test_as_shipped_the_call_keeps_the_whole_tie(monkeypatch):
    _quiet(monkeypatch)
    context, tf, fetch = _lactate_unit()
    cons = IND.infer_inducer("SYNTHETIC", "GntR", genome_ctx=context, tf_gene=tf, allow_ncbi=True,
                             metalnet=False, entry_fetcher=fetch)
    call = next(c for c in cons.calls if c.source == "operon_chem")
    assert len(call.candidates) == 2 and call.mixture_kind == "unresolved"
    assert call.evidence["tie_break"] == {"signals": [], "rank_tie": list(call.candidates), "dropped": []}


def _quiet(monkeypatch, homologs=()):
    monkeypatch.setattr(IND.hsel, "assign_cluster", lambda *a: NS(cluster_id=None, status="x"))
    monkeypatch.setattr(IND, "curated_member_ion", lambda *a, **k: None)
    monkeypatch.setattr(DB, "lookup", lambda **k: None)
    monkeypatch.setattr(DB, "homologs", lambda seq, **k: list(homologs))
    monkeypatch.setattr(IND.metal_site, "coordination_gate", lambda *a, **k: {"has_site": None})


def test_the_call_names_the_winner_and_keeps_the_losers_in_evidence(monkeypatch, signals_on):
    _quiet(monkeypatch)
    context, tf, fetch = _lactate_unit()
    cons = IND.infer_inducer("SYNTHETIC", "GntR", genome_ctx=context, tf_gene=tf, allow_ncbi=True,
                             metalnet=False, entry_fetcher=fetch)
    call = next(c for c in cons.calls if c.source == "operon_chem")
    assert len(call.candidates) == 1 and call.candidates[0].startswith("lactate")
    tb = call.evidence["tie_break"]
    assert tb["signals"] == ["import", "substrate", "recurrence"] and tb["dropped"] == ["pyruvate"]
    assert {c["ligand"] for c in call.evidence["candidates"]} >= {"pyruvate"}            # full list retained
    # an unmapped organic name is still not a headline: the identity policy is unchanged
    assert call.ligand is None and cons.top == "non-metal effector (undetermined)"


# --- task 3: published homologues, organic ligands only -------------------------------------------
def _record(*ligands, rank=90, refseq="WP_OTHER.1"):
    return DB.LigifyRecord(refseq=refseq, organism="Synthetica", ligands=tuple(n for n, _ in ligands),
                           ligand_detail=tuple({"name": n, "smiles": s} for n, s in ligands), rank=rank)


SALICYLATE = ("salicylate", "OC(=O)c1ccccc1O")


@pytest.mark.parametrize("name,smiles,reason", [
    ("Sodium ion", "[Na+]", "monatomic_smiles"), ("Ferrous cation", "[Fe+2]", "monatomic_smiles"),
    ("Mercury", "[Hg]", "monatomic_smiles"), ("Zn2+", None, "metal_class"), ("Fe", None, "metal_class"),
    ("Cl(-)", None, "monatomic_ion_name"),
    ("Iron-sulfur cluster", "[SH-].[SH-].[Fe].[Fe+3]", "metal_in_smiles"),
    ("Molybdate ion", "[O-][Mo](=O)(=O)[O-]", "metal_in_smiles"),
    ("Sodium butyrate", "CCCC(=O)[O-].[Na+]", "metal_in_smiles"),
    ("Hydrogen peroxide", "OO", ""), ("Nitrate", "[N+](=O)([O-])[O-]", ""),
    ("Hydrosulfide", "[SH-]", ""), ("salicylate", SALICYLATE[1], ""), ("CO", None, ""),
])
def test_inorganic_ligands_are_recognised_by_smiles_or_name(name, smiles, reason):
    assert DB.inorganic_reason(name, smiles) == reason


def test_near_identical_homologue_votes_for_its_organic_ligand(monkeypatch):
    _quiet(monkeypatch, [(_record(SALICYLATE), 0.995, 1.0, 1.0)])
    cons = IND.infer_inducer("SYNTHETIC", "GntR", metalnet=False)
    call = next(c for c in cons.calls if c.source == "ligify_db_homolog")
    weight = call.evidence["weight"]
    assert call.ligand == "salicylate" and call.candidates == ("salicylate",)
    assert call.confidence == round(weight * 0.90 * DB.HOMOLOG_VOTE_CONCORDANCE, 2) > 0
    assert call.evidence["votes"] is True and call.evidence["vote_from"] == "WP_OTHER.1"
    assert call.evidence["homologs"][0]["tier"] == "vote"
    assert cons.top == "salicylate"
    assert any(">= 99 %" in n for n in cons.notes)


def test_homologue_below_the_vote_tier_is_context_and_never_a_candidate(monkeypatch):
    _quiet(monkeypatch, [(_record(SALICYLATE), 0.60, 1.0, 1.0)])
    cons = IND.infer_inducer("SYNTHETIC", "GntR", metalnet=False)
    call = next(c for c in cons.calls if c.source == "ligify_db_homolog")
    assert call.confidence == 0.0 and call.ligand is None and call.candidates == ()
    assert call.evidence["votes"] is False and call.evidence["homologs"][0]["tier"] == "context"
    assert cons.top != "salicylate" and "salicylate" not in cons.candidates
    assert any("context only" in n for n in cons.notes)


def test_metal_ligands_of_a_homologue_are_dropped_and_logged(monkeypatch):
    mixed = _record(("Zinc(2+)", "[Zn+2]"), SALICYLATE)
    _quiet(monkeypatch, [(mixed, 1.0, 1.0, 1.0)])
    cons = IND.infer_inducer("SYNTHETIC", "GntR", metalnet=False)
    call = next(c for c in cons.calls if c.source == "ligify_db_homolog")
    assert call.candidates == ("salicylate",)
    assert call.evidence["dropped_inorganic"] == [
        {"refseq": "WP_OTHER.1", "name": "Zinc(2+)", "reason": "monatomic_smiles"}]
    assert all("Zinc" not in str(x) for x in cons.candidates) and cons.top == "salicylate"
    assert any("metal/ion ligands of published homologues are not used" in n for n in cons.notes)


def test_a_homologue_with_only_metal_ligands_does_not_vote(monkeypatch):
    _quiet(monkeypatch, [(_record(("Ferrous cation", "[Fe+2]"), ("Zn2+", None)), 1.0, 1.0, 1.0)])
    cons = IND.infer_inducer("SYNTHETIC", "GntR", metalnet=False)
    call = next(c for c in cons.calls if c.source == "ligify_db_homolog")
    assert call.confidence == 0.0 and call.candidates == () and call.evidence["votes"] is False
    assert len(call.evidence["dropped_inorganic"]) == 2
    assert not any(CI.chemistry(x)[0] == "metal" for x in cons.candidates)


def test_the_best_voting_homologue_is_the_one_that_votes(monkeypatch):
    best = _record(SALICYLATE, refseq="WP_BEST.1")
    other = _record(("formaldehyde", "C=O"), refseq="WP_NEXT.1")
    _quiet(monkeypatch, [(best, 1.0, 1.0, 1.0), (other, 0.992, 1.0, 1.0), (_record(SALICYLATE, refseq="WP_CTX.1"), 0.7, 1.0, 1.0)])
    cons = IND.infer_inducer("SYNTHETIC", "GntR", metalnet=False)
    call = next(c for c in cons.calls if c.source == "ligify_db_homolog")
    assert call.candidates == ("salicylate",) and call.evidence["vote_from"] == "WP_BEST.1"
    assert [h["tier"] for h in call.evidence["homologs"]] == ["vote", "vote", "context"]


def test_an_all_inorganic_best_hit_does_not_pass_the_vote_down(monkeypatch):
    metal_best = _record(("Ferrous cation", "[Fe+2]"), refseq="WP_BEST.1")
    _quiet(monkeypatch, [(metal_best, 1.0, 1.0, 1.0), (_record(SALICYLATE, refseq="WP_NEXT.1"), 0.995, 1.0, 1.0)])
    cons = IND.infer_inducer("SYNTHETIC", "GntR", metalnet=False)
    call = next(c for c in cons.calls if c.source == "ligify_db_homolog")
    assert call.confidence == 0.0 and call.candidates == () and call.evidence["vote_from"] is None
    assert cons.top != "salicylate"


def test_a_protein_with_its_own_record_does_not_also_vote_through_a_homologue(monkeypatch):
    _quiet(monkeypatch, [(_record(SALICYLATE, rank=100), 1.0, 1.0, 1.0)])
    own = DB.LigifyRecord(refseq="WP_SELF.1", ligands=("formaldehyde",),
                          ligand_detail=({"name": "formaldehyde", "smiles": "C=O"},), rank=50)
    monkeypatch.setattr(DB, "lookup", lambda **k: own)
    cons = IND.infer_inducer("SYNTHETIC", "GntR", metalnet=False)
    call = next(c for c in cons.calls if c.source == "ligify_db_homolog")
    assert call.confidence == 0.0 and call.candidates == ()
    assert call.evidence["vote_withheld_exact_record"] is True
    assert cons.top == "formaldehyde" and cons.support_sources == ("ligify_db",)


def test_homologue_lookup_reaches_the_context_floor_only_when_asked(monkeypatch):
    from predictor.annotate import homologs as H
    recs = [{"refseq": "WP_NEAR.1", "protein_seq": "MAAA", "candidate_ligands": []},
            {"refseq": "WP_MID.1", "protein_seq": "MAAC", "candidate_ligands": []},
            {"refseq": "WP_FAR.1", "protein_seq": "MCCC", "candidate_ligands": []}]
    monkeypatch.setattr(DB, "_all_records", lambda: recs)
    monkeypatch.setattr(DB, "_records_fasta", lambda: ("unused.faa", recs))
    hits = [H.Homolog("r0", 0.995, 0, 99, 1.0, 1.0), H.Homolog("r1", 0.62, 0, 80, 0.95, 0.95),
            H.Homolog("r2", 0.41, 0, 50, 1.0, 1.0)]
    monkeypatch.setattr(H, "homolog_sweep", lambda *a, **k: NS(hits=hits))
    assert [r.refseq for r, *_ in DB.homologs("MAAA")] == ["WP_NEAR.1"]                   # default: vote tier
    got = DB.homologs("MAAA", min_identity=DB.HOMOLOG_CONTEXT_MIN_IDENTITY)
    assert [r.refseq for r, *_ in got] == ["WP_NEAR.1", "WP_MID.1"]                       # 41 % is not reported
    assert (DB.HOMOLOG_MIN_IDENTITY, DB.HOMOLOG_CONTEXT_MIN_IDENTITY) == (0.99, 0.50)


# --- task 4: source precedence --------------------------------------------------------------------
def call(source, ligand, confidence, role=""):
    return InducerCall(source=source, ligand=ligand, confidence=confidence, role=role,
                       identity_status="resolved")


def test_tiers_are_clade_then_published_then_local():
    assert S.source_tier("ssn_cluster") < S.source_tier("ligify_db") < S.source_tier("operon_chem")
    assert S.source_tier("ligify_db") == S.source_tier("ligify_db_homolog")
    assert S.source_tier("operon_chem") == S.source_tier("ligify") == S.source_tier("some_new_source")


def test_a_louder_local_name_does_not_outrank_a_clade_call():
    # E. coli FrmR: the clade call sits at 0.50, the co-transcribed local tie at 0.77.
    calls = [call("ssn_cluster", "formaldehyde", 0.50, "non-metal"), call("operon_chem", "salicylate", 0.77)]
    cons = InducerConsensus.from_calls(calls)
    assert cons.top == "formaldehyde"
    assert cons.candidates == ("formaldehyde", "salicylate")          # the lower tier is kept, not dropped
    assert cons.support_sources == ("ssn_cluster",)


def test_published_outranks_local_and_clade_outranks_published():
    assert InducerConsensus.from_calls([call("operon_chem", "salicylate", 0.9),
                                        call("ligify_db", "formaldehyde", 0.3)]).top == "formaldehyde"
    assert InducerConsensus.from_calls([call("ligify_db_homolog", "salicylate", 0.9),
                                        call("ssn_cluster", "formaldehyde", 0.2, "non-metal")]).top == "formaldehyde"


def test_confidence_still_decides_inside_one_tier():
    assert InducerConsensus.from_calls([call("ligify_db", "salicylate", 0.4),
                                        call("ligify_db_homolog", "formaldehyde", 0.7)]).top == "formaldehyde"
    assert InducerConsensus.from_calls([call("operon_chem", "salicylate", 0.4),
                                        call("regulon", "formaldehyde", 0.7)]).top == "formaldehyde"


def test_a_non_voting_higher_tier_call_does_not_block_a_lower_one():
    silent = InducerCall(source="ssn_cluster", ligand=None, confidence=0.0)
    unresolved = InducerCall(source="ligify_db", ligand=None, confidence=0.8, candidates=("a", "b"),
                             mixture_kind="unresolved", identity_status="unresolved")
    cons = InducerConsensus.from_calls([silent, unresolved, call("operon_chem", "salicylate", 0.3)])
    assert cons.top == "salicylate"


def test_finalize_explains_a_precedence_decision_and_keeps_the_candidate():
    calls = [call("ssn_cluster", "formaldehyde", 0.50, "non-metal"),
             InducerCall(source="coordination", ligand=None, confidence=0.0, role="non-metal",
                         evidence={"gate": "CsoR-His/Cys", "family": "CsoR"}),
             InducerCall(source="operon_chem", ligand="salicylate", confidence=0.77, identity_status="resolved",
                         evidence={"weight": 1.0})]
    cons = IND._finalize(calls, False, gate_name="CsoR-His/Cys", sensor_class="organic")
    assert cons.top == "formaldehyde" and "salicylate" in cons.candidates
    assert any(n.startswith("source precedence") and "operon_chem" in n for n in cons.notes)
    # agreement between the tiers needs no note
    same = IND._finalize([calls[0], calls[1], call("operon_chem", "formaldehyde", 0.77)], False,
                         sensor_class="organic")
    assert same.top == "formaldehyde" and not any(n.startswith("source precedence") for n in same.notes)


def test_surviving_organic_fallback_follows_the_same_precedence():
    low_tier_loud = call("ligify_db_homolog", "salicylate", 0.9)
    clade_quiet = call("ssn_cluster", "formaldehyde", 0.3, "non-metal")
    assert IND._surviving_organic_call([low_tier_loud, clade_quiet]) is clade_quiet
    assert IND._surviving_organic_call([low_tier_loud, call("ligify_db", "formaldehyde", 0.95)]).ligand == "formaldehyde"
    assert IND._surviving_organic_call([call("operon_chem", "salicylate", 0.9)]) is None      # local is never the fallback
