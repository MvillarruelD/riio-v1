"""
coords.py -- lift every method's native coordinates into the canonical genome frame (Phase 0.1).

`schema.py` defines the frame `(genome_accession, strand, dyad_center, half_site_spans, offset_to_tss)`.
This module is the set of `to_genome_coords()` adapters the workplan mandates -- one per evidence type --
plus a Tomtom-style PWM aligner, because **motif PCC is undefined without an offset/orientation search**.

Adapters (native type -> CanonicalOperator):
  from_motif_finder   region-relative OperatorCandidate (+ region genome start)   [signals/motif_finder]
  from_genome_hit     GenomeHit (already genome coords)                            [signals/motif_rescan]
  from_regulon_operon RegulonOperon (operator site that drives an operon)          [signals/regulon]

Coordinate convention everywhere: 0-based, half-open, + strand. `region_gstart` is the genome 0-based
start of the extracted region whose + strand slice the region-relative offsets index (the `cli.predict`
convention: `region_gstart = region_start0 + region.start`).

Run `python coords.py` for a self-test (each adapter + a PWM-alignment round trip with shift & rc).
"""
from __future__ import annotations

import numpy as np

from predictor.schema import (CanonicalOperator, HalfSite,  # noqa: F401 - re-exported
                              OperatorPrediction, Provenance,
                              PREDICTED_UNVERIFIED)

_COMP = [3, 2, 1, 0]                       # A<->T, C<->G row permutation (rows are A,C,G,T)
_RC = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}


def _revcomp_str(s: str) -> str:
    return "".join(_RC.get(c, "N") for c in reversed((s or "").upper()))


# --------------------------------------------------------------------------- per-method adapters
def from_motif_finder(cand, *, genome_accession, region_gstart, strand="+", provenance=None):
    """signals.motif_finder.OperatorCandidate (region-relative, + strand) -> CanonicalOperator.

    `region_gstart` lifts region offsets to genome coords. IR/DR sites are dyad-symmetric so `strand`
    records the read strand only (default '+'); pass the candidate's PWM/importance through unchanged.
    """
    start = region_gstart + int(cand.start)
    end = region_gstart + int(cand.end)
    dyad = region_gstart + int(round(cand.center))
    half_sites = None
    if getattr(cand, "half_sites", None):
        (a0, a1), (b0, b1) = cand.half_sites
        half_sites = (HalfSite(region_gstart + a0, region_gstart + a1, strand),
                      HalfSite(region_gstart + b0, region_gstart + b1, strand))
    return CanonicalOperator(
        genome_accession=genome_accession, dyad_center=dyad, start=start, end=end, strand=strand,
        half_sites=half_sites, spacer_len=getattr(cand, "spacer", None), seq=cand.seq,
        pwm=getattr(cand, "pwm", None), per_base_importance=getattr(cand, "per_base_importance", None),
        sequence_score=getattr(cand, "final", None), score=getattr(cand, "score", None),
        confidence=(getattr(cand, "conservation", None)),
        generator=getattr(cand, "kind", ""),
        provenance=provenance or Provenance("motif_finder", PREDICTED_UNVERIFIED))


def from_genome_hit(hit, *, provenance=None):
    """signals.motif_rescan.GenomeHit (already genome coords) -> CanonicalOperator.

    `hit.matched` is ALWAYS the + strand substring (pwm_scan convention), even for a - strand hit. Store
    `seq` in the motif READING orientation (reverse-complement for - strand hits) so logos/PWMs built
    across hits of mixed strand stack in register (operator_logo.counts_from_seqs has no orientation
    handling). Genome coords (start/end/dyad) stay in the + strand frame.
    """
    matched = getattr(hit, "matched", "") or ""
    seq = _revcomp_str(matched) if getattr(hit, "strand", "+") == "-" else matched
    return CanonicalOperator(
        genome_accession=hit.accession, dyad_center=int(hit.dyad_center), start=int(hit.start),
        end=int(hit.end), strand=hit.strand, seq=seq,
        sequence_score=getattr(hit, "combined", None), score=getattr(hit, "score", None),
        pvalue=getattr(hit, "pvalue", None), qvalue=getattr(hit, "qvalue", None),
        qvalue_kind=getattr(hit, "qvalue_kind", None),
        generator=f"rescan_tier{getattr(hit, 'tier', '?')}",
        provenance=provenance or Provenance("motif_rescan", PREDICTED_UNVERIFIED))


def from_regulon_operon(op, accession, *, provenance=None):
    """signals.regulon.RegulonOperon (the operator site driving an operon) -> CanonicalOperator."""
    dyad = (int(op.site_start) + int(op.site_end)) // 2
    return CanonicalOperator(
        genome_accession=accession, dyad_center=dyad, start=int(op.site_start), end=int(op.site_end),
        strand=op.site_strand, score=getattr(op, "score", None),
        pvalue=getattr(op, "pvalue", None), qvalue=getattr(op, "qvalue", None),
        sequence_score=getattr(op, "score", None), generator=f"regulon:{op.first_gene}",
        provenance=provenance or Provenance("regulon", PREDICTED_UNVERIFIED))


# --------------------------------------------------------------------------- PWM alignment (Tomtom-style)
def _as_4xL(p):
    p = np.asarray(p, dtype=float)
    if p.shape[0] != 4 and p.shape[1] == 4:
        p = p.T
    return p


def rc_pwm(p):
    """Reverse-complement of a 4xL probability matrix (rows A,C,G,T)."""
    return _as_4xL(p)[_COMP][:, ::-1]


def align_pwms(a, b, *, min_overlap=4):
    """Tomtom-style alignment of two PWMs of possibly different widths. Slides `b` (and its
    reverse-complement) across `a`, scoring each register by Pearson correlation over the overlapping
    columns (flattened 4xk). Returns the best register.

    -> dict(pcc, offset, orientation, overlap): `offset` is b's start column relative to a (can be
    negative); `orientation` is '+' (b as-is) or '-' (b reverse-complemented).
    """
    a = _as_4xL(a)
    best = {"pcc": -2.0, "offset": 0, "orientation": "+", "overlap": 0}
    La = a.shape[1]
    for orient, bb in (("+", _as_4xL(b)), ("-", rc_pwm(b))):
        Lb = bb.shape[1]
        for off in range(-(Lb - min_overlap), La - min_overlap + 1):
            lo = max(0, off)
            hi = min(La, off + Lb)
            k = hi - lo
            if k < min_overlap:
                continue
            av = a[:, lo:hi].ravel()
            bv = bb[:, lo - off:hi - off].ravel()
            if av.std() < 1e-12 or bv.std() < 1e-12:
                pcc = 0.0
            else:
                pcc = float(np.corrcoef(av, bv)[0, 1])
            if pcc > best["pcc"]:
                best = {"pcc": pcc, "offset": off, "orientation": orient, "overlap": k}
    return best


def motif_pcc(a, b, *, min_overlap=4) -> float:
    """Best aligned per-column Pearson correlation between two PWMs (the motif-PCC metric)."""
    return align_pwms(a, b, min_overlap=min_overlap)["pcc"]


# --------------------------------------------------------------------------- self-test
def _demo() -> None:
    from .motif_finder import OperatorCandidate
    from .motif_rescan import GenomeHit
    from .regulon import RegulonOperon

    GSTART = 1000

    # 1) motif_finder region offsets -> genome
    mfc = OperatorCandidate(center=17.5, start=0, end=35, seq="A" * 35, score=9.0, kind="IR",
                            half=8, spacer=19, half_sites=((0, 8), (27, 35)),
                            conservation=0.83, final=16.5, pwm=np.full((4, 35), 0.25))
    op = from_motif_finder(mfc, genome_accession="NC_000964.3", region_gstart=GSTART)
    print(f"motif_finder -> {op.genome_accession}:{op.start}-{op.end} dyad={op.dyad_center} "
          f"spacer={op.spacer_len} conf={op.confidence}")
    assert op.start == GSTART and op.end == GSTART + 35 and op.dyad_center == GSTART + 18
    assert op.half_sites[1].start == GSTART + 27 and op.half_sites[1].width == 8
    assert op.provenance.method == "motif_finder" and not op.provenance.oracle

    # 2) GenomeHit already canonical
    gh = GenomeHit("NC_X.1", 500, 514, "-", dyad_center=507, score=8.4, pvalue=1e-9,
                   qvalue=1e-7, tier=1, prior=1.0, combined=8.42, matched="TTGACCTAGGTCAA")
    g = from_genome_hit(gh)
    assert g.dyad_center == 507 and g.strand == "-" and g.generator == "rescan_tier1"

    # a - strand hit must store its seq in READING orientation (revcomp of the + strand `matched`),
    # so logos built across mixed-strand hits stack in register (operator_logo has no orientation handling)
    gh_asym = GenomeHit("NC_X.1", 700, 706, "-", dyad_center=703, score=8.0, pvalue=1e-9,
                        qvalue=1e-7, tier=1, prior=1.0, combined=8.0, matched="AAATTC")  # rc -> GAATTT
    assert from_genome_hit(gh_asym).seq == "GAATTT", "minus-strand seq must be reverse-complemented"
    plus = from_genome_hit(GenomeHit("NC_X.1", 700, 706, "+", dyad_center=703, score=8.0, pvalue=1e-9,
                                     qvalue=1e-7, tier=1, prior=1.0, combined=8.0, matched="AAATTC"))
    assert plus.seq == "AAATTC", "plus-strand seq must be unchanged"

    # 3) RegulonOperon site
    ro = RegulonOperon("targA", "+", ["targA", "targB"], ["WP_A", "WP_B"], 1460, 1474, "+",
                       score=7.0, pvalue=1e-8, qvalue=5e-9, tier=1)
    r = from_regulon_operon(ro, "SYN.1")
    assert r.start == 1460 and r.dyad_center == 1467 and r.generator == "regulon:targA"

    # 4) PWM alignment: a peaked motif, recovered after a shift and after reverse-complement
    rng = np.random.default_rng(0)
    base = np.full((4, 8), 0.04)
    for j, ch in enumerate([0, 3, 3, 1, 2, 0, 2, 3]):     # a definite consensus
        base[ch, j] = 0.88
    a = np.concatenate([np.full((4, 3), 0.25), base, np.full((4, 3), 0.25)], axis=1)  # base at offset 3
    same = align_pwms(a, base)
    print(f"align same motif: pcc={same['pcc']:.3f} offset={same['offset']} orient={same['orientation']}")
    assert same["pcc"] > 0.95 and same["offset"] == 3 and same["orientation"] == "+"

    rc = align_pwms(a, rc_pwm(base))
    assert rc["pcc"] > 0.95 and rc["orientation"] == "-", "reverse-complement not detected"

    rand = align_pwms(a, rng.dirichlet(np.ones(4), size=8).T)
    assert rand["pcc"] < same["pcc"], "a random motif should not align as well as the true one"
    print(f"align random motif: pcc={rand['pcc']:.3f}  (< {same['pcc']:.3f}, as expected)")

    print("OK: every method lifts to the canonical frame; PWM alignment recovers shift + orientation.")


if __name__ == "__main__":
    _demo()
