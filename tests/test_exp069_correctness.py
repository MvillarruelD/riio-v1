"""EXP069 phase 1: interval provenance, independent of operon distance heuristics."""
from __future__ import annotations

from types import SimpleNamespace

import pytest

from predictor.effector import operon_chemistry as OC
from predictor.effector.ligand import _synthetic_entry
from predictor.signals import regulon
from predictor.signals.motif_rescan import Gene, GenomeContext


def _layout(strand):
    member = Gene(100, 400, "-", "member", protein_id="WP_member", product="enzyme")
    partner = Gene(500, 800, "-", "partner", protein_id="WP_partner", product="enzyme")
    tf = Gene(900, 1200, "+", "tf", protein_id="WP_tf", product="regulator")
    genes = [member, partner, tf]
    if strand == "-":
        for gene in genes:
            gene.start, gene.end = 1500 - gene.end, 1500 - gene.start
            gene.strand = "+" if gene.strand == "-" else "-"
    return sorted(genes, key=lambda gene: gene.start), tf, partner, member


def _context(genes, boundary=None):
    offsets = {} if boundary is None else {"first": 0, "second": boundary}
    return GenomeContext("SYN.1", "N" * 4000, genes, contig_offsets=offsets)


def _fetching(calls):
    def fetch(protein_id):
        calls.append(protein_id)
        return _synthetic_entry()
    return fetch


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("boundary", [None, 50, 2500])
def test_same_replicon_controls_preserve_candidates_and_geometry(strand, boundary):
    genes, tf, partner, member = _layout(strand)
    baseline = OC.effectors_from_context(_context(genes), tf, entry_fetcher=lambda _: _synthetic_entry())
    ctx = _context(genes, boundary)
    fetched = []
    got = OC.effectors_from_context(ctx, tf, entry_fetcher=_fetching(fetched))
    assert got == baseline and got
    assert fetched == [partner.protein_id, member.protein_id]
    assert OC.divergent_partner(genes, tf, ctx=ctx) is partner
    assert OC.operon_from_partner(genes, partner, ctx=ctx) == [partner, member]
    assert (OC.MAX_PROMOTER_BP, OC.OPERON_CHEM_MAX_GAP_BP, OC.OPERON_CHEM_MAX_SPAN_BP,
            regulon.DEFAULT_OPERON_MAX_GAP) == (1000, 300, 8000, 150)


@pytest.mark.parametrize("strand,boundary", [("+", 850), ("-", 650)])
def test_cross_replicon_partner_rejected_before_annotation(strand, boundary):
    genes, tf, partner, _ = _layout(strand)
    ctx = _context(genes, boundary)
    fetched = []
    # Geometry alone would admit this partner: the promoter is only 100 bp wide.
    assert OC.divergent_partner(genes, tf) is partner
    assert OC.divergent_partner(genes, tf, ctx=ctx) is None
    assert OC.effectors_from_context(ctx, tf, entry_fetcher=_fetching(fetched)) == []
    assert fetched == []


@pytest.mark.parametrize("strand,boundary", [("+", 450), ("-", 1050)])
def test_cross_replicon_member_terminates_walk_before_annotation(strand, boundary):
    genes, tf, partner, member = _layout(strand)
    ctx = _context(genes, boundary)
    fetched = []
    assert OC.operon_from_partner(genes, partner) == [partner, member]
    assert OC.operon_from_partner(genes, partner, ctx=ctx) == [partner]
    assert OC.effectors_from_context(ctx, tf, entry_fetcher=_fetching(fetched))
    assert fetched == [partner.protein_id]


@pytest.mark.parametrize("strand,boundary", [("+", 1100), ("-", 450)])
def test_tf_interval_crossing_boundary_is_rejected_before_annotation(strand, boundary):
    genes, tf, _, _ = _layout(strand)
    fetched = []
    assert OC.effectors_from_context(_context(genes, boundary), tf,
                                     entry_fetcher=_fetching(fetched)) == []
    assert fetched == []


@pytest.mark.parametrize("strand,boundary", [("+", 650), ("-", 850)])
def test_partner_interval_crossing_boundary_is_rejected_before_annotation(strand, boundary):
    genes, tf, _, _ = _layout(strand)
    fetched = []
    assert OC.effectors_from_context(_context(genes, boundary), tf,
                                     entry_fetcher=_fetching(fetched)) == []
    assert fetched == []


@pytest.mark.parametrize("strand,boundary", [("+", 250), ("-", 1250)])
def test_member_interval_crossing_boundary_is_rejected_before_annotation(strand, boundary):
    genes, tf, partner, _ = _layout(strand)
    fetched = []
    assert OC.effectors_from_context(_context(genes, boundary), tf, entry_fetcher=_fetching(fetched))
    assert fetched == [partner.protein_id]


@pytest.mark.parametrize("strand", ["+", "-"])
def test_locate_checks_both_half_open_interval_endpoints(strand, monkeypatch):
    genes, tf, _, _ = _layout(strand)
    ctx = _context(genes, 2500)
    visited = []
    locate = ctx.locate

    def recording_locate(pos):
        visited.append(pos)
        return locate(pos)

    monkeypatch.setattr(ctx, "locate", recording_locate)
    assert OC.effectors_from_context(ctx, tf, entry_fetcher=lambda _: _synthetic_entry())
    assert {pos for gene in genes for pos in (gene.start, gene.end - 1)} <= set(visited)


def test_interval_ending_at_boundary_is_within_its_original_replicon():
    genes, tf, partner, _ = _layout("+")
    ctx = _context(genes, partner.end)
    # `end` is exclusive: the partner ends at this boundary but does not cross it.
    assert OC._gene_replicon(ctx, partner) == "first"
    assert OC._gene_replicon(ctx, tf) == "second"


def test_single_record_synthetic_context_without_locate_remains_usable():
    genes, tf, _, _ = _layout("+")
    ctx = SimpleNamespace(accession="SYN.1", sequence="N" * 4000, genes=genes)
    assert OC.effectors_from_context(ctx, tf, entry_fetcher=lambda _: _synthetic_entry())


@pytest.mark.parametrize("offsets", [{"first": 0, "second": 2500}, {"second": 50}])
def test_mapped_context_without_resolver_cannot_claim_same_replicon(offsets):
    genes, tf, _, _ = _layout("+")
    ctx = SimpleNamespace(accession="SYN.1", sequence="N" * 4000, genes=genes,
                          contig_offsets=offsets)
    fetched = []
    assert OC.effectors_from_context(ctx, tf, entry_fetcher=_fetching(fetched)) == []
    assert fetched == []


@pytest.mark.parametrize("strand", ["+", "-"])
def test_interval_before_first_mapped_offset_rejects_accession_fallback(strand):
    genes, tf, _, _ = _layout(strand)
    ctx = GenomeContext("r1", "N" * 4000, genes, contig_offsets={"r1": 2500})
    fetched = []
    assert ctx.locate(tf.start) == ("r1", tf.start), "exercise locate's accession fallback"
    assert OC.divergent_partner(genes, tf, ctx=ctx) is None
    assert OC.effectors_from_context(ctx, tf, entry_fetcher=_fetching(fetched)) == []
    assert fetched == []


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("endpoint", ["start", "end"])
def test_mapped_local_coordinates_must_agree_with_declared_offset(strand, endpoint, monkeypatch):
    genes, tf, _, _ = _layout(strand)
    ctx = _context(genes, 50)
    locate = ctx.locate
    target = tf.start if endpoint == "start" else tf.end - 1

    def inconsistent_locate(pos):
        replicon, local = locate(pos)
        return replicon, local + 1 if pos == target else local

    monkeypatch.setattr(ctx, "locate", inconsistent_locate)
    fetched = []
    assert OC.effectors_from_context(ctx, tf, entry_fetcher=_fetching(fetched)) == []
    assert fetched == []


@pytest.mark.parametrize("start,end", [(-1, 400), (900, 900), (900, 4001)])
def test_invalid_tf_interval_cannot_produce_an_annotation_request(start, end):
    genes, tf, _, _ = _layout("+")
    tf.start, tf.end = start, end
    fetched = []
    assert OC.effectors_from_context(_context(genes), tf, entry_fetcher=_fetching(fetched)) == []
    assert fetched == []
