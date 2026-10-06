"""Anchor sets and cluster-assignment abstention (Phase 5).

These pin two behaviours that were previously impossible and are easy to regress:

  * a cluster's label rests on a SET of anchors with support and conflict, not one representative;
  * `assign_cluster` can decline to answer, and says why.

Both are load-bearing for honest benchmarking, so the tests assert the *reasons* rather than just
the outcomes.
"""
from __future__ import annotations

import pytest

from predictor.annotate import family_kb, homolog_selection as hs, ssn_clusters as ssn

@pytest.fixture(autouse=True)
def require_packaged_kb_for_released_data(request):
    """Delay record loading until a released-data test body was actually selected."""
    if request.node.get_closest_marker("released_data") and not family_kb.available():
        pytest.skip("packaged anchor knowledge base is absent")


@pytest.mark.released_data
def test_anchor_sets_are_populated():
    """Populated and not truncated -- counted the way the KB now stores it.

    Records are one per distinct SEQUENCE (168 of them), but every accession that shares a sequence
    is still carried. Asserting on records alone would pass even if the accessions were thrown away,
    which is the failure that would break `is_anchor` and leak into leave-one-out; so the accession
    total is asserted too.
    """
    sets = family_kb._load()
    assert len(sets) >= 60, f"expected ~71 clusters, got {len(sets)}"
    records = sum(len(s.anchors) for s in sets.values())
    accessions = sum(len(a.accessions or (a.accession,))
                     for s in sets.values() for a in s.anchors)
    assert records >= 150, f"expected ~168 distinct-sequence records, got {records}"
    assert accessions > 2500, f"accessions dropped: {accessions} carried across {records} records"


@pytest.mark.released_data
def test_multiple_independent_anchors_exist():
    """The whole point of the redesign: labels backed by more than one protein."""
    multi = [s for s in family_kb._load().values() if s.n_independent_anchors > 1]
    assert len(multi) >= 20, f"only {len(multi)} clusters have >1 independent anchor"


@pytest.mark.released_data
def test_independent_anchors_count_sequences_not_rows():
    """64 identical TrEMBL copies of one protein are ONE piece of evidence."""
    for s in family_kb._load().values():
        assert s.n_independent_anchors <= s.n_anchors
        distinct = len({a.seq_md5 or a.accession for a in s.anchors})
        assert s.n_independent_anchors <= max(distinct, 1) or not s.anchors


@pytest.mark.released_data
def test_conflicts_are_surfaced_not_averaged():
    conflicted = [s for s in family_kb._load().values() if s.conflict]
    assert conflicted, "expected clusters whose anchors disagree on the inducer"
    for s in conflicted:
        assert s.n_distinct_inducers > 1
        assert s.agreement < 1.0, f"{s.cluster_id} reports full agreement despite a conflict"
        assert "mixed" in s.tier


@pytest.mark.released_data
def test_leave_one_anchor_out_drops_the_whole_sequence_group():
    """Naming ANY accession of a protein must remove that protein, not just one of its names.

    A protein reaches the KB under many accessions -- one sequence carries 213 of them. Holding out
    the representative while the others still anchored the cluster would leave the label standing on
    exactly the evidence being withheld. Since records are stored per sequence, the test that
    matters is that a NON-representative accession is just as effective a holdout as the
    representative one.
    """
    target, anchor = None, None
    for s in family_kb._load().values():
        for a in s.anchors:
            if len(a.accessions or ()) > 5:
                target, anchor = s, a
                break
        if target:
            break
    assert target is not None, "expected an anchor carrying several accessions"

    other = next(x for x in anchor.accessions if x.upper() != anchor.accession.upper())
    for held in (anchor.accession, other, other.lower()):
        red = target.excluding([held])
        assert red.n_anchors == target.n_anchors - 1,             f"holding out {held!r} removed {target.n_anchors - red.n_anchors} anchors, expected 1"
        assert not any(x.upper() == held.upper()
                       for a in red.anchors for x in (a.accessions or (a.accession,))),             f"{held!r} survives the holdout under another record"
    assert family_kb.is_anchor(other, target.cluster_id),         "a non-representative accession must still be recognised as an anchor"


@pytest.mark.released_data
def test_is_anchor_is_not_truncated():
    """A vendored anchor list must never be capped: a missed anchor reads as independent evidence."""
    for cid, s in family_kb._load().items():
        assert len(s.anchors) == s.n_anchors or s.n_anchors == 0, \
            f"{cid}: vendored {len(s.anchors)} anchors but n_anchors={s.n_anchors}"


@pytest.mark.released_data
def test_cluster_info_exposes_the_anchor_set():
    with_anchors = [c for c in ssn.load_clusters() if c.n_independent_anchors > 0]
    assert with_anchors, "no cluster exposes an anchor set"
    c = max(with_anchors, key=lambda x: x.n_independent_anchors)
    # checked against the vocabulary itself, so a tier added to the data (as `structural-only` was on
    # 2026-09-17) fails here only if it was not also declared in family_kb.TIERS
    assert c.label_confidence in family_kb.TIERS, c.label_confidence
    assert {x.label_confidence for x in ssn.load_clusters()} <= set(family_kb.TIERS)
    assert len(c.anchors) >= c.n_independent_anchors


# --------------------------------------------------------------------------- abstention
def test_identity_floor_is_measured_not_zero():
    assert hs.MIN_IDENTITY >= 0.5, "the floor is what stops out-of-clade proteins being labelled"
    assert hs.MIN_CLUSTER_SIZE == 0, "size floor measured inert; enabling it only costs recall"


def test_family_without_an_ssn_says_so():
    a = hs.assign_cluster("MKTAYIAKQRQISFVKSHFSRQ", "LacI/GalR")
    assert a.cluster_id is None and a.status == hs.STATUS_NO_SSN
    assert not a.assigned


@pytest.mark.released_data
def test_true_member_still_assigns():
    info = ssn.info_for_cluster("MerR_9")
    _acc, seq = next(ssn.iter_sequences(info))
    a = hs.assign_cluster(seq, ssn.MERR, top_n=30)
    assert a.cluster_id == "MerR_9" and a.status == hs.STATUS_ASSIGNED and a.assigned


@pytest.mark.released_data
def test_unassigned_reasons_are_distinguishable():
    """`cluster_id=None` used to conflate three different situations."""
    assert len({hs.STATUS_NO_SSN, hs.STATUS_NO_HIT, hs.STATUS_LOW_IDENTITY,
                hs.STATUS_LOW_SUPPORT, hs.STATUS_SMALL_CLUSTER}) == 5
    a = hs.assign_cluster("WWWWWWWWWWWWWWWWWWWWWWWWWWWWWW", ssn.MERR, top_n=10)
    assert a.cluster_id is None
    assert a.status.startswith("unassigned") or a.status == hs.STATUS_NO_SSN


def test_fur_structural_zn_not_called_as_effector():
    """Fur structural Zn(Cys4) must not override the sensed metal (site-typing correction).

    Regression guard for the two-site fix: a query whose nearest curated binder is a Fur
    structural-Zn anchor must still be called by the clade's sensed metal, and Fe-Fur must not
    carry a spurious Zn in its candidate shortlist.
    """
    import json
    from predictor import resources
    from predictor.annotate import family_kb

    binders = json.loads((resources.SSN_DATABASE / "curated_binders_Fur.json").read_text(encoding="utf-8"))
    by_cluster = {}
    for b in binders:
        by_cluster.setdefault(b["cluster"], b)

    # every Fur binder is site-typed
    assert all("site_type" in b for b in binders)
    # a Mn sensor clade whose only curated anchor binds Zn is typed structural
    baMur = [b for b in binders if b["cluster"] == "Fur_baMur"]
    assert baMur and all(b["site_type"] == "structural" for b in baMur)

    family_kb._load.cache_clear()
    # Fe-Fur consensus is no longer the structural Zn
    ec = family_kb.anchor_set("Fur_ecFur")
    assert ec is not None and not ec.real_conflict
    assert "Zn2+" not in ec.distinct_inducers
    family_kb._load.cache_clear()
