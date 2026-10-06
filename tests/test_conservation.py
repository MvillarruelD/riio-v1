"""Phylogenetic-conservation scoring: positional coherence + independence weighting.

`motif_finder.conservation` decides which candidate becomes the genome-rescan PWM (`seed_seqs` prefers
conserved candidates and ranks by `score * (1 + conservation)`), so it is the highest-leverage input to
operator recall. It used to score "the fraction of homolog promoters containing any significant hit
anywhere, on either strand", which credits two things that are not evidence of a conserved operator:
hits at unrelated positions, and repeated near-identical strains of one over-sequenced organism.

Every test here is built so the OLD score would return 1.0, and only a position-aware,
redundancy-aware score can tell the cases apart.
"""
import random

import pytest

from predictor.annotate import genome_resolver as gr
from predictor.signals import motif_finder as MF

MOTIF = "TTGACATGGCTAAAGACAATTACATAACATGTCAA"     # the module self-test's MerR-style dyad


def _filler(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def _region(motif, at, *, total=400, seed=0):
    """A promoter of `total` bp carrying `motif` starting at offset `at`."""
    left = _filler(at, seed)
    right = _filler(total - at - len(motif), seed + 10_000)
    return left + motif + right


class TestPositionalCoherence:
    def test_same_site_at_the_same_position_scores_high(self):
        regions = [_region(MOTIF, 150, seed=i) for i in range(8)]
        assert MF.conservation(MOTIF, regions) > 0.9

    def test_same_site_at_scattered_positions_scores_low(self):
        """Every region contains the motif, so the old score was 1.0; but the offsets disagree by far
        more than the drift a real orthologous site shows, so this is not one conserved site."""
        offsets = [20, 90, 160, 230, 300, 350, 60, 260]
        regions = [_region(MOTIF, off, seed=i) for i, off in enumerate(offsets)]
        scattered = MF.conservation(MOTIF, regions)
        aligned = MF.conservation(MOTIF, [_region(MOTIF, 150, seed=i) for i in range(8)])
        assert scattered < aligned
        assert scattered < 0.5, scattered

    def test_small_drift_is_tolerated(self):
        """Indels in the intergenic region shift a real site by a few bp; that must not be punished."""
        regions = [_region(MOTIF, 150 + d, seed=i)
                   for i, d in enumerate([0, 4, -6, 9, -3, 7, -8, 2])]
        assert MF.conservation(MOTIF, regions) > 0.9

    def test_position_is_measured_from_the_gene_proximal_end(self):
        """Windows clipped at a contig edge lose 5' sequence. Measuring from the 3' (gene) end keeps a
        conserved site coherent; measuring from the 5' end would scatter it by the clip amount."""
        full = [_region(MOTIF, 150, total=400, seed=i) for i in range(4)]
        clipped = [r[60:] for r in full]          # lose 60 bp off the 5' end
        assert MF.conservation(MOTIF, full + clipped) > 0.9


class TestIndependenceWeighting:
    def test_a_universal_site_stays_1_whether_or_not_genomes_are_redundant(self):
        """Redundancy must NOT be penalised when the site really is in every promoter: the score is a
        weighted FRACTION, so a duplicate group shrinks the numerator and the denominator together."""
        dup = _region(MOTIF, 150, seed=1)
        redundant = [dup] * 7 + [_region(MOTIF, 150, seed=99)]
        diverse = [_region(MOTIF, 150, seed=i) for i in range(8)]
        assert MF.conservation(MOTIF, redundant) == pytest.approx(1.0)
        assert MF.conservation(MOTIF, diverse) == pytest.approx(1.0)

    def test_redundancy_cannot_inflate_a_PARTIAL_signal(self):
        """Where it bites: the site is in 4 of 8 promoters either way, so the old fraction-of-regions
        score reported 0.5 for both. But in the redundant set those 4 are one re-sequenced strain --
        a single observation -- while in the diverse set they are four independent genomes."""
        dup = _region(MOTIF, 150, seed=1)
        redundant = [dup] * 4 + [_filler(400, 500 + i) for i in range(4)]
        diverse = [_region(MOTIF, 150, seed=i) for i in range(4)] + \
                  [_filler(400, 500 + i) for i in range(4)]
        assert MF.conservation(MOTIF, redundant) < MF.conservation(MOTIF, diverse)
        assert MF.conservation(MOTIF, redundant) == pytest.approx(0.2)
        assert MF.conservation(MOTIF, diverse) == pytest.approx(0.5)

    def test_weights_collapse_duplicate_groups(self):
        a, b = _region(MOTIF, 150, seed=1), _region(MOTIF, 150, seed=2)
        w = MF.independence_weights([a, a, a, b])
        assert w[:3] == [pytest.approx(1 / 3)] * 3
        assert w[3] == pytest.approx(1.0)
        assert sum(w) == pytest.approx(2.0)      # two independent observations

    def test_distinct_promoters_keep_full_weight(self):
        regions = [_region(MOTIF, 150, seed=i) for i in range(6)]
        assert MF.independence_weights(regions) == [1.0] * 6

    def test_a_site_in_one_over_sequenced_clade_does_not_look_universal(self):
        """The failure this guards: the motif is real in ONE species (sequenced many times) and absent
        elsewhere. Fraction-of-regions calls that highly conserved; weighted evidence does not."""
        clade = _region(MOTIF, 150, seed=5)
        others = [_filler(400, 200 + i) for i in range(4)]
        assert MF.conservation(MOTIF, [clade] * 6 + others) <= 0.5


class TestNoHits:
    def test_absent_motif_scores_near_zero(self):
        """A lone chance hit can still clear the per-region Bonferroni threshold (it did here, in 1 of 6
        random promoters -- the old score reported the same 1/6). What must not happen is a chance hit
        in ONE region looking like a conserved site: with nothing to agree with, it cannot."""
        assert MF.conservation(MOTIF, [_filler(400, 300 + i) for i in range(6)]) <= 1 / 6 + 1e-9

    def test_empty_inputs_are_safe(self):
        assert MF.conservation(MOTIF, []) == 0.0
        assert MF.conservation("ACGT", [_region(MOTIF, 150)]) == 0.0     # too short to score


class TestPromoterOrientationConvention:
    """Positional conservation is only meaningful because every extractor now returns the window in
    transcriptional orientation. Before that, a minus-strand promoter was mirrored AND complemented."""

    def test_minus_strand_window_is_reverse_complemented(self):
        seq = "AAAAGGGGCCCCTTTA"
        assert gr.orient_promoter(seq, "+") == seq
        assert gr.orient_promoter(seq, "-") == "TAAAGGGGCCCCTTTT"
        assert gr.orient_promoter(gr.orient_promoter(seq, "-"), "-") == seq

    def test_gene_start_sits_at_the_same_distance_from_the_3_prime_end_on_both_strands(self):
        """`conservation` compares hits by distance from the gene-proximal (3') end. After orientation
        that distance is exactly `downstream` for both strands -- which is what makes the coordinate
        comparable across homologs, and invariant to 5' clipping at a contig edge."""
        up, down = gr.PROMOTER_UPSTREAM, gr.PROMOTER_DOWNSTREAM
        tf_start, tf_end = 10_000, 10_300

        s1, e1 = gr.promoter_window_bounds(tf_start, tf_end, "+", upstream=up, downstream=down)
        assert e1 - tf_start == down                      # plus: 3' end is e1

        s1, e1 = gr.promoter_window_bounds(tf_start, tf_end, "-", upstream=up, downstream=down)
        assert tf_end - s1 == down                        # minus: after RC the 3' end is s1

    def test_all_three_extractors_share_the_convention(self):
        """A region set mixing conventions cannot be compared position-by-position, so the local,
        NCBI and SSN-cluster extractors must all route through `orient_promoter`."""
        import inspect
        from predictor.annotate import homolog_regions, msa_homologs
        for fn in (homolog_regions.regions_from_candidates,
                   gr.fetch_upstream_region,
                   msa_homologs._upstream):
            assert "orient_promoter" in inspect.getsource(fn), fn.__qualname__
