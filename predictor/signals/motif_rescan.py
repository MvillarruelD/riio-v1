"""
motif_rescan.py -- PWM -> genome scan under a tiered locality prior (Phase 1.6A).

Decouples operator discovery from "the intergenic region immediately before the TF" -- Snowprint's
biggest failure mode (operator outside the extracted region, 48/80). Given ANY method's PWM (Snowprint
conservation, DeepPBS contact weights, our finder's consensus, or a family-library consensus), it scans
a genome and returns ranked hits, restricted to intergenic/promoter DNA and re-weighted by a locality
prior:

    tier 1  intergenic immediately flanking the regulator   prior 1.00
    tier 2  divergent promoters (head-to-head gene pair)     prior 0.60
    tier 3  any other intergenic region                      prior 0.35
    tier 4  coding DNA (only when scope="genome")            prior 0.10

The statistical engine is `motifs/pwm_scan.py` (our FIMO-equivalent: exact p-values). Combined rank =
-log10(p) + log10(prior), so a strong hit in a low-prior tier can still outrank a weak hit nearby, but
locality breaks ties -- matching the precision@k / FDR-at-prevalence metric the benchmark uses. The
structural verifier then folds only the top handful (closed loop, Appendix A).

Run `python motif_rescan.py` for a self-test (plants tier-1, distal, and coding sites).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field, replace

try:                                    # works as package and as `python motif_rescan.py`
    from ..motifs.pwm_scan import scan, background_from_seqs, pwm_from_counts
except ImportError:
    from predictor.motifs.pwm_scan import scan, background_from_seqs, pwm_from_counts


# --------------------------------------------------------------------------- genome model
@dataclass
class Gene:
    start: int          # 0-based half-open, + strand coordinates
    end: int
    strand: str         # '+' or '-'
    name: str = ""
    protein_id: str = ""    # CDS protein accession (from GenBank/GFF), used by effector inference
    product: str = ""       # free-text product description, used to flag enzymes vs regulators
    #: The GFF locus_tag. `name` prefers the human gene name, which is NOT unique -- SL1344 has 82
    #: duplicated names (rrf x8, tnpA x8), 21 of them on more than one replicon -- so a reported
    #: gene name alone cannot be resolved back to a locus. The locus_tag can.
    locus_tag: str = ""


@dataclass
class GenomeContext:
    accession: str
    sequence: str
    genes: list          # list[Gene]
    circular: bool = False
    origin_wrapped_features: tuple[str, ...] = ()
    #: replicon/contig id -> its offset in `sequence`. A multi-replicon assembly is concatenated
    #: (N-padded) into ONE coordinate space so the scan sees the whole genome, and GFF coordinates
    #: are lifted into it. Empty for a single-record genome. Carry it: without the map every
    #: coordinate reported downstream silently inherits `accession`, which is only the FIRST
    #: record -- a plasmid or a second chromosome then reports as the primary one.
    contig_offsets: dict = field(default_factory=dict)

    def locate(self, pos: int) -> tuple:
        """Concatenated coordinate -> (replicon id, coordinate within that replicon).

        Returns `(accession, pos)` unchanged for a single-record genome, so callers need no
        special case."""
        if not self.contig_offsets:
            return self.accession, pos
        best_id, best_off = self.accession, -1
        for cid, off in self.contig_offsets.items():
            if off <= pos and off > best_off:
                best_id, best_off = cid, off
        return (best_id, pos - best_off) if best_off >= 0 else (self.accession, pos)


@dataclass
class GenomeHit:
    accession: str
    start: int
    end: int
    strand: str
    dyad_center: int
    score: float
    pvalue: float
    qvalue: float
    tier: int
    prior: float
    combined: float
    matched: str = ""
    qvalue_kind: str = "conservative_bh_upper_bound"


DEFAULT_PRIORS = {1: 1.0, 2: 0.6, 3: 0.35, 4: 0.1}


# --------------------------------------------------------------------------- interval logic
def _merge(intervals):
    out = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def gene_spans(genes):
    return _merge([(g.start, g.end) for g in genes])


def intergenic_intervals(genes, length, *, circular: bool = False):
    """Gaps not covered by genes.

    For a circular replicon, the two terminal fragments are one origin-spanning
    interval.  The current canonical coordinate schema cannot represent a
    wrapped site, so that interval is omitted rather than misreported as two
    ordinary linear promoters.
    """
    gaps, prev = [], 0
    for s, e in gene_spans(genes):
        if s > prev:
            gaps.append((prev, s))
        prev = max(prev, e)
    if prev < length:
        gaps.append((prev, length))
    if circular:
        gaps = [(start, end) for start, end in gaps if start != 0 and end != length]
    return gaps


def _flanking(genes, gap):
    a, b = gap
    left = max((g for g in genes if g.end <= a), key=lambda g: g.end, default=None)
    right = min((g for g in genes if g.start >= b), key=lambda g: g.start, default=None)
    return left, right


def classify_tier(gap, genes, regulator):
    left, right = _flanking(genes, gap)
    if regulator is not None and (left is regulator or right is regulator):
        return 1
    if left is not None and right is not None and left.strand == "-" and right.strand == "+":
        return 2                        # divergent promoters point away from this gap
    return 3


def _resolve_regulator(genes, regulator):
    if regulator is None or isinstance(regulator, Gene):
        return regulator
    for g in genes:                     # allow passing a gene name or protein accession
        if g.name == regulator or g.protein_id == regulator:
            return g
    return None


# --------------------------------------------------------------------------- public API
def rescan_genome(ctx: GenomeContext, *, pwm=None, counts=None, regulator=None,
                  scope: str = "intergenic", pvalue_thresh: float = 1e-4,
                  both_strands: bool = True, prior_weights=None, boost_regions=None):
    """Scan `ctx` with a motif; return hits ranked by locality-weighted significance.

    scope="intergenic" (default) restricts to tiers 1-3; scope="genome" also scans coding DNA
    (tier 4) -- the loop's last-resort escalation when the upstream window fails.

    `boost_regions` (optional list of (start, end) genome intervals, e.g. promoters of inducer-consistent
    regulated genes from `report.ligand_regulon.ligand_genes`) lift a hit's locality prior toward tier-2
    when its dyad falls inside one -- a SOFT, additive re-rank for non-autoregulatory targets, never a
    gate (a hit outside all boost regions is scored exactly as before).
    """
    if pwm is None and counts is None:
        raise ValueError("provide pwm= or counts=")
    if pwm is None:
        pwm = pwm_from_counts(counts)
    W = pwm.shape[1] if pwm.shape[0] == 4 else pwm.shape[0]
    priors = prior_weights or DEFAULT_PRIORS
    reg = _resolve_regulator(ctx.genes, regulator)
    bg = background_from_seqs([ctx.sequence])
    boosts = [(int(lo), int(hi)) for lo, hi in (boost_regions or [])]

    gaps = intergenic_intervals(ctx.genes, len(ctx.sequence), circular=ctx.circular)
    intervals = [(s, e, classify_tier((s, e), ctx.genes, reg)) for s, e in gaps]
    if scope == "genome":
        intervals += [(s, e, 4) for s, e in gene_spans(ctx.genes)]

    hits = []
    n_tested = 0
    for s, e, tier in intervals:
        sub = ctx.sequence[s:e]
        if len(sub) < W:
            continue
        n_tested += (len(sub) - W + 1) * (2 if both_strands else 1)
        for h in scan(sub, prob=pwm, background=bg,
                      pvalue_thresh=pvalue_thresh, both_strands=both_strands):
            gs, ge = s + h.start, s + h.end
            dyad = (gs + ge) // 2
            prior = priors.get(tier, 0.1)
            if boosts and any(lo <= dyad < hi for lo, hi in boosts):
                prior = max(prior, priors.get(2, 0.6))    # regulated-gene promoter -> tier-2-like locality
            combined = -math.log10(max(h.pvalue, 1e-300)) + math.log10(prior)
            hits.append(GenomeHit(ctx.accession, gs, ge, h.strand, dyad,
                                  h.score, h.pvalue, h.qvalue, tier, prior, combined, h.matched))

    hits = _merge_adjacent(hits, W)
    hits = _genome_wide_qvalues(hits, n_tested)
    hits.sort(key=lambda x: x.combined, reverse=True)
    return hits


def _merge_adjacent(hits: list, width: int) -> list:
    """Collapse hits that describe the SAME site into one, keeping the best-scoring representative.

    `scan` is run per window and on both strands, so one real palindromic operator surfaces as several
    GenomeHits: the + and - strand reading of the same span, plus shifted windows overlapping it. Left
    unmerged they trebled the apparent site count, inflated false positives, and distorted the FDR
    denominator. The duplication is directly visible in the benchmark: true sites appear in adjacent-rank
    PAIRS (percentiles 0.00/0.01, 0.28/0.29, 0.62/0.63) -- two records of one site.

    Two hits are the same site when their dyads lie within half a motif width."""
    if not hits:
        return hits
    tol = max(1, width // 2)
    by_score = sorted(hits, key=lambda h: (-h.score, h.start))
    kept: list = []
    claimed: list = []                                    # dyads of sites already taken, sorted-insert
    for h in by_score:
        if any(abs(h.dyad_center - d) <= tol for d in claimed):
            continue
        claimed.append(h.dyad_center)
        kept.append(h)
    return kept


def _genome_wide_qvalues(hits: list, n_tested: int) -> list:
    """Conservative BH upper bounds over the whole genome scan.

    Two things were wrong before. (1) `scan()` applies BH inside each intergenic interval it is called on,
    so an interval yielding a single hit gets q == p; that made a downstream `qvalue_thresh` gate
    structurally incapable of doing anything, which is why an earlier benchmark found the q-value lever
    "inert". (2) The denominator must be the number of windows actually TESTED, not the number of hits
    that survived the p-value pre-filter -- BH over survivors alone assigns every survivor q == p and
    passes them all, since it silently drops the millions of tests that failed.

    The scanner retains only windows below ``pvalue_thresh``.  Without the discarded p-values an exact
    BH step-up adjustment cannot be reconstructed; using the total number tested and the retained ranks
    yields a conservative upper bound.  It is safe for fail-closed significance checks but must not be
    labelled an exact q-value.

    With the honest denominator the verdict is stark and worth stating plainly: scanning ~4.6 Mb on both
    strands is ~9.3e6 tests, so a hit needs p <= ~5e-9 to reach q <= 0.05, whereas the operator PWMs
    deliver p ~ 1e-4. **No site survives a genome-wide FDR.** That is a real property of the motifs, not a
    threshold to tune, and it is why site-level precision is ~2%. Callers should therefore select by RANK
    (the score ordering does carry signal -- true sites are ~4x enriched in the top decile) and report the
    q-value for what it is: an honest significance statement, not a usable filter."""
    n = len(hits)
    if n == 0:
        return hits
    N = max(int(n_tested), n)
    order = sorted(range(n), key=lambda i: hits[i].pvalue)
    q_prev, out = 1.0, [1.0] * n
    for rank, i in enumerate(reversed(order), start=1):    # walk worst -> best for the running minimum
        k = n - rank + 1                                   # BH rank among the retained hits
        q_prev = min(q_prev, hits[i].pvalue * N / k)
        out[i] = min(1.0, q_prev)
    return [replace(h, qvalue=out[i], qvalue_kind="conservative_bh_upper_bound") for i, h in enumerate(hits)]


class MotifRescan:
    """Thin wrapper matching the closed-loop call site (Appendix A.6)."""

    def __init__(self, ctx: GenomeContext):
        self.ctx = ctx

    def genome_hits(self, pwm, *, regulator=None, scope="intergenic", **kw):
        return rescan_genome(self.ctx, pwm=pwm, regulator=regulator, scope=scope, **kw)


# --------------------------------------------------------------------------- self-test
def _demo() -> None:
    import random
    import numpy as np

    motif = "TTGACCTAGGTCAA"            # distinctive 14-bp palindrome
    W = len(motif)
    counts = np.zeros((4, W))
    idx = {"A": 0, "C": 1, "G": 2, "T": 3}
    for j, ch in enumerate(motif):
        counts[idx[ch], j] += 1
    pwm = pwm_from_counts(counts)

    random.seed(3)
    L = 1400
    seq = list("".join(random.choice("ACGT") for _ in range(L)))

    def plant(at):
        seq[at:at + W] = list(motif)

    plant(320)      # tier-1 intergenic gap 250..400 (flanks regulator)
    plant(700)      # inside coding gene g2 600..800 -> excluded when scope=intergenic
    plant(1165)     # distal intergenic gap 1150..1200 -> tier 3
    seq = "".join(seq)

    genes = [
        Gene(100, 250, "-", "g0"),
        Gene(400, 550, "+", "reg"),     # the regulator: + strand; g0 is - -> divergent at 250..400
        Gene(600, 800, "+", "g2"),
        Gene(1000, 1150, "+", "g3"),
        Gene(1200, 1350, "+", "g4"),
    ]
    ctx = GenomeContext("synthetic", seq, genes)

    hits = rescan_genome(ctx, pwm=pwm, regulator="reg", scope="intergenic", pvalue_thresh=1e-4)
    print(f"intergenic scope: {len(hits)} hit(s)")
    for h in hits[:6]:
        print(f"  {h.start:>5}-{h.end:<5}({h.strand}) tier={h.tier} "
              f"p={h.pvalue:.1e} combined={h.combined:.2f}  {h.matched}")

    def covers(hits, pos):
        return [h for h in hits if h.start <= pos < h.end]

    assert covers(hits, 320), "tier-1 site missed"
    assert covers(hits, 1165), "distal tier-3 site missed"
    assert not covers(hits, 700), "coding-region site should be excluded under scope=intergenic"
    assert hits[0].start <= 320 < hits[0].end and hits[0].tier == 1, "tier-1 site should rank first"
    assert covers(hits, 320)[0].tier == 1 and covers(hits, 1165)[0].tier == 3, "tier mislabelled"

    genome_hits = rescan_genome(ctx, pwm=pwm, regulator="reg", scope="genome", pvalue_thresh=1e-4)
    assert any(h.tier == 4 and h.start <= 700 < h.end for h in genome_hits), \
        "coding site should appear as tier 4 under scope=genome"

    print("OK: tier-1 ranks first, coding site excluded (intergenic) / tier-4 (genome), "
          "distal site tier-3.")


if __name__ == "__main__":
    _demo()
