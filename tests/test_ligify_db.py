"""The vendored Ligify database and the lookup source built on it.

`effector/rank.py` was already a faithful port of Ligify's scoring formula, so what we lacked was their
DATA. `predictor/data/refs/ligify_db.json.gz` is the trimmed, compressed `ligifyDB.json` from
`groov-bio/ligify-ui` (MIT). These tests pin the three properties that make it safe to consult: the
redistribution carries its licence, the lookup never matches fuzzily, and a published organic-ligand
answer is weighted by exactly the same metal/redox rule as our own derivation -- so it cannot walk a
metalloregulator's inducer call over to a neighbouring enzyme's metabolite.
"""
import pytest

from predictor.effector import ligify_db as LDB

pytestmark = pytest.mark.released_data


@pytest.fixture(autouse=True)
def require_vendored_ligify_db():
    """Delay record loading until a released-data test body was actually selected."""
    if not LDB.available():
        pytest.skip("no vendored Ligify database (tools/vendor_ligify_db.py)")


def test_redistribution_carries_its_provenance_and_licence():
    """MIT requires the notice to travel with a redistributed substantial portion."""
    m = LDB.meta()
    assert "MIT" in m["_license"]
    assert m["_source"].startswith("https://raw.githubusercontent.com/groov-bio/ligify-ui/")
    assert m["_home"].startswith("https://github.com/groov-bio/ligify-ui")
    # enough to trace the vendored subset back to one upstream state
    assert len(m["_upstream_sha256"]) == 64
    assert m["_upstream_records"] >= 3000
    assert m["_retrieved"] and m["_kept_fields"] and m["_dropped"]


def test_lookup_resolves_by_accession_and_by_sequence():
    db = LDB._load()
    raw = next(iter(db["by_md5"].values()))
    by_acc = LDB.lookup(acc=raw["refseq"])
    by_seq = LDB.lookup(seq=raw["protein_seq"])
    assert by_acc and by_seq
    assert by_acc.refseq == by_seq.refseq == raw["refseq"]
    assert by_acc.matched_by == "accession" and by_seq.matched_by == "sequence"
    # the accession is version-insensitive: WP_x.1 and WP_x are the same protein
    assert LDB.lookup(acc=raw["refseq"].split(".", 1)[0]) is not None


def test_lookup_never_matches_fuzzily():
    """Attributing a published biosensor's ligand to a homolog is the error this must not make."""
    db = LDB._load()
    raw = next(iter(db["by_md5"].values()))
    mutated = raw["protein_seq"][:-5] + "AAAAA"          # a near neighbour, not the same protein
    assert LDB.lookup(seq=mutated) is None
    assert LDB.lookup(acc="WP_000000000.1") is None
    assert LDB.lookup(seq="", acc="") is None
    assert LDB.lookup() is None


def test_sequence_lookup_abstains_when_genomic_context_is_ambiguous():
    """Identical proteins in different operons must not inherit an arbitrary first ligand call."""
    db = LDB._load()
    records = next(iter(db["ambiguous_md5"].values()))

    assert len(records) > 1
    assert LDB.lookup(seq=records[0]["protein_seq"]) is None
    assert all(LDB.lookup(acc=raw["refseq"]) is not None for raw in records)


def test_records_expose_their_rank_on_our_scale():
    """Their rank is the same 0-100 operon-context score `effector.rank` produces, so the two compare."""
    db = LDB._load()
    ranks = [r.get("rank", {}).get("rank") for r in list(db["by_md5"].values())[:200]]
    ranks = [x for x in ranks if x is not None]
    assert ranks and max(ranks) <= 100
    # negative ranks are real upstream (a big operon with competing regulators) and must survive the
    # round-trip, because the inducer source clamps them to a 0 confidence rather than hiding them
    rec = LDB.lookup(acc=next(iter(db["by_acc"])))
    assert rec is not None and isinstance(rec.rank_metrics, dict)


def test_published_answer_is_weighted_by_the_same_metal_rule_as_our_own():
    """A published organic ligand must not become a metalloregulator's inducer for free.

    `ligify_weight` returns ~0.1 when the coordination gate fired, and the `ligify_db` source multiplies
    its confidence by exactly that factor -- the same one `ligify` uses. Without this, consulting the
    published database would reintroduce the failure mode the weighting was built to stop.
    """
    from predictor.effector import inducer as IND
    db = LDB._load()
    raw = next(r for r in db["by_md5"].values()
               if (r.get("rank") or {}).get("rank", 0) > 60 and r.get("candidate_ligands"))

    # the same record, read as a metalloregulator (gate fires) and as an organic-ligand family
    metal = IND.infer_inducer(raw["protein_seq"], "MerR", protein_acc=raw["refseq"],
                              allow_ncbi=False, metalnet_fn=lambda _s: None)
    organic = IND.infer_inducer(raw["protein_seq"], "TetR/AcrR", protein_acc=raw["refseq"],
                                allow_ncbi=False, metalnet_fn=lambda _s: None)
    m_call = next((c for c in metal.calls if c.source == "ligify_db"), None)
    o_call = next((c for c in organic.calls if c.source == "ligify_db"), None)
    assert m_call is not None and o_call is not None, "the DB source must fire for both families"
    assert m_call.ligand == o_call.ligand, "same record, same published ligand"
    assert m_call.evidence["weight"] <= o_call.evidence["weight"], (
        "a metalloregulator family must not weight a published organic ligand ABOVE an organic family "
        f"({m_call.evidence['weight']} vs {o_call.evidence['weight']})")
    assert m_call.confidence <= o_call.confidence


def test_source_is_named_separately_from_our_own_derivation():
    """A reader must be able to tell a published answer from one we recomputed."""
    from predictor.effector import inducer as IND
    db = LDB._load()
    raw = next(iter(db["by_md5"].values()))
    cons = IND.infer_inducer(raw["protein_seq"], "TetR/AcrR", protein_acc=raw["refseq"],
                             allow_ncbi=False, metalnet_fn=lambda _s: None)
    call = next((c for c in cons.calls if c.source == "ligify_db"), None)
    assert call is not None
    assert call.source == "ligify_db" and call.source != "ligify"
    assert call.evidence["source_db"].startswith("groov-bio/ligify-ui")
    assert call.evidence["matched_by"] in ("accession", "sequence")
    assert call.evidence["refseq"] == raw["refseq"]
