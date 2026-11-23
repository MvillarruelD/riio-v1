"""The regulator is the query's own gene, never another gene that shares its name.

The pipeline used to pass the regulator's gene NAME to the genome rescan, which resolved it back to the
first gene carrying that name. Names repeat (M. avium has several cmtR and ricR), so tier 1 -- the
intergenic gaps flanking the regulator -- was anchored on another copy for 2 of the 238 canonical
regulators. The pipeline now passes the Gene; a string is resolved only when it is unambiguous.
"""
from __future__ import annotations

import pytest

from predictor.signals import motif_rescan as MR


def _genes():
    return [MR.Gene(100, 400, "+", "cmtR", "WP_A", locus_tag="T_0001"),
            MR.Gene(900, 1200, "-", "cmtR", "WP_B", locus_tag="T_0002"),
            MR.Gene(2000, 2300, "+", "zntA", "WP_C", locus_tag="T_0003")]


def test_a_gene_object_is_used_as_given():
    genes = _genes()
    assert MR._resolve_regulator(genes, genes[1]) is genes[1]


def test_accession_and_locus_tag_resolve_to_their_own_copy():
    genes = _genes()
    assert MR._resolve_regulator(genes, "WP_B") is genes[1]
    assert MR._resolve_regulator(genes, "T_0002") is genes[1]


def test_a_unique_name_resolves_and_a_repeated_name_is_refused():
    genes = _genes()
    assert MR._resolve_regulator(genes, "zntA") is genes[2]
    with pytest.raises(ValueError, match="matches 2 genes"):
        MR._resolve_regulator(genes, "cmtR")
    assert MR._resolve_regulator(genes, "absent") is None


def test_tier_one_flanks_the_second_copy_when_that_copy_is_the_regulator():
    genes = _genes()
    flanks = MR._FlankIndex(genes)
    gap_next_to_second = (400, 900)          # between copy 1 and copy 2
    gap_after_second = (1200, 2000)          # only copy 2 flanks it
    assert MR.classify_tier(gap_after_second, genes, genes[1], flanks) == 1
    assert MR.classify_tier(gap_after_second, genes, genes[0], flanks) != 1
    assert MR.classify_tier(gap_next_to_second, genes, genes[1], flanks) == 1
