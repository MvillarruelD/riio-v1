"""EXP071 tests (EXP071_PREREG_20261109.md), synthetic evidence only."""
from types import SimpleNamespace as NS

import pytest

from predictor.effector import chemical_identity as CI, inducer as IND, ligify_db as DB
from predictor.effector import operon_chemistry as OC
from predictor.signals.motif_rescan import Gene, GenomeContext


def entry(equation):
    return {"comments": [{"commentType": "CATALYTIC ACTIVITY", "reaction": {
        "name": equation, "reactionCrossReferences": [{"database": "Rhea", "id": "RHEA:SYN"}]}}]}


def ctx(genes):
    return GenomeContext("SYN.1", "N" * 40000, genes)


# --- change 1: the regulator's own transcription run ---------------------------------------------
def test_run_walks_both_ways_and_stops_at_strand_gap_and_noncoding():
    a = Gene(1000, 1500, "+", "a", protein_id="A")
    b = Gene(1550, 1900, "+", "b", protein_id="B")
    tf = Gene(2000, 2400, "+", "tf", protein_id="TF")
    c = Gene(2450, 2900, "+", "c", protein_id="C")
    rna = Gene(2950, 3050, "+", "rna")                                   # non-coding stops the walk
    d = Gene(3100, 3500, "+", "d", protein_id="D")
    far = Gene(200, 600, "+", "far", protein_id="FAR")                  # 400 bp gap to `a` stops it
    opp = Gene(1600, 1700, "-", "opp", protein_id="OPP")
    genes = sorted([far, a, b, tf, c, rna, d], key=lambda g: g.start)
    assert [g.name for g in OC.cotranscribed_run(genes, tf, ctx=ctx(genes))] == ["a", "b", "c"]
    genes2 = sorted([a, opp, tf], key=lambda g: g.start)
    assert OC.cotranscribed_run(genes2, tf) == []                        # opposite strand interrupts


@pytest.mark.parametrize("strand", ["+", "-"])
def test_tail_regulator_reads_its_upstream_operon_with_route_tag(strand):
    # frlABCD-R style: enzyme then TF on the same strand; no divergent partner.
    enz = Gene(1000, 1900, strand, "enz", protein_id="ENZ", product="kinase")
    tf = Gene(1950, 2400, strand, "tf", protein_id="TF", product="regulator")
    cands = OC.effectors_from_context(ctx([enz, tf]), tf,
                                      entry_fetcher=lambda p: entry("salicylate = catechol"))
    assert {c.ligand for c in cands} == {"salicylate", "catechol"}
    assert all(c.route == "cotranscribed" and c.distance == 1 and c.rank == 100 for c in cands)
    assert all(o["route"] == "cotranscribed" for c in cands for o in c.occurrences)


def test_both_routes_feed_one_list_and_each_protein_is_fetched_once():
    part = Gene(1000, 1900, "-", "part", protein_id="P", product="dehydrogenase")   # divergent
    tf = Gene(2000, 2400, "+", "tf", protein_id="TF", product="regulator")
    own = Gene(2450, 3000, "+", "own", protein_id="O", product="permease")          # co-transcribed
    fetched = []

    def fetch(pid):
        fetched.append(pid)
        return entry("glycolate = glyoxylate") if pid == "P" else entry("lactate(out) = lactate(in)")
    cands = OC.effectors_from_context(ctx([part, tf, own]), tf, entry_fetcher=fetch)
    routes = {c.ligand: c.route for c in cands}
    lactate = [n for n in routes if n.startswith("lactate")]
    assert len(lactate) == 1                                     # (in)/(out) are one identity
    assert routes["glycolate"] == "divergent" and routes[lactate[0]] == "cotranscribed"
    assert sorted(fetched) == ["O", "P"]


# --- change 2: published homologues below the vote tier are context, never a vote ----------------
# (EXP071 pinned this at 0.995; EXP072 lets the >= 99 % tier vote for organic ligands -- see
# tests/test_exp072.py -- so the non-voting contract is now pinned on the context tier.)
def test_homologue_call_never_votes_and_is_noted(monkeypatch):
    rec = DB.LigifyRecord(refseq="WP_OTHER.1", organism="Synthetica", ligands=("salicylate",),
                          ligand_detail=({"name": "salicylate", "smiles": "OC(=O)c1ccccc1O"},), rank=90)
    monkeypatch.setattr(IND.hsel, "assign_cluster", lambda *a: NS(cluster_id=None, status="x"))
    monkeypatch.setattr(IND, "curated_member_ion", lambda *a, **k: None)
    monkeypatch.setattr(DB, "lookup", lambda **k: None)
    monkeypatch.setattr(DB, "homologs", lambda seq, **k: [(rec, 0.97, 1.0, 1.0)])
    monkeypatch.setattr(IND.metal_site, "coordination_gate", lambda *a, **k: {"has_site": None})
    cons = IND.infer_inducer("SYNTHETIC", "GntR", metalnet=False)
    call = next(c for c in cons.calls if c.source == "ligify_db_homolog")
    assert call.confidence == 0.0 and call.ligand is None and call.candidates == ()
    assert call.evidence["homologs"][0]["refseq"] == "WP_OTHER.1"
    assert cons.top != "salicylate" and "salicylate" not in cons.candidates
    assert any("context only" in n for n in cons.notes)


def test_homologue_lookup_excludes_the_query_s_own_record(monkeypatch):
    from predictor.annotate import homologs as H
    recs = [{"refseq": "WP_SELF.1", "protein_seq": "MAAA", "candidate_ligands": []},
            {"refseq": "WP_OTHER.2", "protein_seq": "MAAA", "candidate_ligands": []},
            {"refseq": "WP_FAR.1", "protein_seq": "MCCC", "candidate_ligands": []}]
    monkeypatch.setattr(DB, "_all_records", lambda: recs)
    monkeypatch.setattr(DB, "_records_fasta", lambda: ("unused.faa", recs))
    hits = [H.Homolog("r0", 1.0, 0, 99, 1.0, 1.0), H.Homolog("r1", 0.995, 0, 98, 0.95, 0.95),
            H.Homolog("r2", 0.97, 0, 90, 1.0, 1.0)]
    monkeypatch.setattr(H, "homolog_sweep", lambda *a, **k: NS(hits=hits))
    got = DB.homologs("MAAA", exclude_refseq=["WP_SELF"])
    assert [(r.refseq, i) for r, i, _, _ in got] == [("WP_OTHER.2", 0.995)]       # self out, 97 % below cut


# --- change 3: a class placeholder lists its same-class shortlist explicitly ---------------------
def test_metal_shortlist_display_separates_other_proposals():
    shown = CI.display_label("divalent metal (ion unresolved)", ("Cu+", "Au(I)", "salicylate"))
    assert shown == "metal-ion shortlist: Cu+ or Au(I); other proposals: salicylate"
    # a lone same-class member is a declined alternative, not a shortlist
    assert not CI.display_label("divalent metal (ion unresolved)", ("Zn2+",)).startswith("metal-ion shortlist")
