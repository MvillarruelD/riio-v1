"""The autoregulatory anchor uses the regulator's real intergenic gap on BOTH strands.

On the plus strand the old loop took the last gene starting before the TF's end -- the TF itself -- as the
upstream neighbour, so the gap came out negative and every plus-strand TF fell through to the fixed
160-bp fallback (anchor always 75 bp upstream). 61 of the 238 canonical anchors moved when this was
fixed; on 102 curated E. coli autoregulatory operators the |site - anchor| <= 80 bp test scored 72/102
before and after, so the fix is a correctness fix of the anchor, not a recovery gain.
"""
from __future__ import annotations

from predictor.annotate.context import locate_operator_window
from predictor.signals.motif_rescan import Gene, GenomeContext

SEQ = "A" * 3000


def test_plus_strand_window_is_the_gap_to_the_upstream_gene():
    ctx = GenomeContext("x", SEQ, [Gene(100, 900, "+", "up"), Gene(1000, 1400, "+", "tf"),
                                   Gene(1500, 2000, "-", "dn")])
    assert locate_operator_window(ctx, 1000, 1400, "+") == (950, 900, 1000)


def test_both_strands_are_mirror_images():
    plus = GenomeContext("x", SEQ, [Gene(100, 900, "+", "up"), Gene(1000, 1400, "+", "tf")])
    minus = GenomeContext("x", SEQ, [Gene(1000, 1400, "-", "tf"), Gene(1500, 2300, "-", "up")])
    c_p, lo_p, hi_p = locate_operator_window(plus, 1000, 1400, "+")
    c_m, lo_m, hi_m = locate_operator_window(minus, 1000, 1400, "-")
    assert hi_p - lo_p == hi_m - lo_m == 100
    assert 1000 - c_p == c_m - 1400 == 50


def test_a_gap_shorter_than_twenty_bp_still_uses_the_fixed_window():
    ctx = GenomeContext("x", SEQ, [Gene(100, 990, "+", "up"), Gene(1000, 1400, "+", "tf")])
    assert locate_operator_window(ctx, 1000, 1400, "+") == (925, 840, 1000)


def test_a_long_gap_caps_the_anchor_75_bp_from_the_start():
    ctx = GenomeContext("x", SEQ, [Gene(100, 400, "+", "up"), Gene(1000, 1400, "+", "tf")])
    assert locate_operator_window(ctx, 1000, 1400, "+")[0] == 925
