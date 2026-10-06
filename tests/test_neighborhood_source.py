"""The genomic-neighbourhood inducer source: what it may and may not do.

Measured behaviour behind these rules is in
`analysis/discovery_bench/NEIGHBORHOOD_DECISION_20260819.md`. The source is the weakest metal evidence
the pipeline carries, so its precedence and its silence matter more than its sensitivity.
"""
from __future__ import annotations

from predictor.annotate import neighborhood as NB
from predictor.effector import inducer as IND
from predictor.signals.motif_rescan import Gene, GenomeContext


def _ctx(*genes):
    return GenomeContext(accession="TEST", sequence="A" * 20000, genes=list(genes))


TF = Gene(1000, 1400, "+", "tfR", "WP_TF", "transcriptional regulator")
COPA = Gene(0, 900, "+", "copA", "WP_CU", "copper-translocating P-type ATPase")
FAR_COPA = Gene(9000, 9600, "+", "copA", "WP_CU", "copper-translocating P-type ATPase")
GYRA = Gene(1500, 2100, "+", "gyrA", "WP_GY", "DNA gyrase subunit A")


class TestClassification:
    def test_substring_traps_do_not_fire(self):
        """HANDOFF 5.4 -- these three strings contain 'co', 'ni' and 'metal' respectively."""
        assert NB.categorize("steroid-CoA dehydrogenase", "scdA") not in NB.METAL_CATEGORIES
        assert NB.categorize("nitrogen regulatory protein", "glnB") not in NB.METAL_CATEGORIES
        assert NB.categorize("non-metal effector binding protein", "yfoo") not in NB.METAL_CATEGORIES

    def test_short_symbols_need_a_word_boundary(self):
        """`cora`, `moda`, `dps` and `bfr` are short enough to appear inside unrelated words."""
        assert NB.categorize("decorated protein", "decorA") not in NB.METAL_CATEGORIES
        assert NB.categorize("modating enzyme", "remodase") not in NB.METAL_CATEGORIES

    def test_real_metal_machinery_is_still_caught(self):
        assert NB.categorize("copper-translocating P-type ATPase", "copA") == "metal_efflux_transport"
        assert NB.categorize("mercuric reductase", "merA") == "metal_resistance_enzyme"
        assert NB.categorize("TonB-dependent siderophore receptor", "fhuA") == "iron_siderophore"
        assert NB.categorize("bacterioferritin", "bfr") == "metal_storage_chelation"


class TestDistanceRule:
    def test_fires_inside_500bp(self):
        call = NB.analyse(_ctx(COPA, TF), TF)
        assert call.fired and call.metal_within_200
        assert call.nearest_metal.gene == "copA" and call.nearest_metal.gap == 100

    def test_does_not_fire_outside_500bp(self):
        call = NB.analyse(_ctx(TF, FAR_COPA), TF)
        assert call.nearest_metal is not None, "a distant metal gene is still recorded"
        assert not call.fired, "but 7600 bp away must not fire the 500 bp rule"
        assert NB.inducer_call(call) is None

    def test_no_metal_neighbour_abstains(self):
        call = NB.analyse(_ctx(TF, GYRA), TF)
        assert call.status == "ok" and not call.fired
        assert NB.inducer_call(call) is None, "abstain, do not guess"


class TestTheSourceNeverNamesAnIon:
    def test_ligand_is_always_none(self):
        call = NB.analyse(_ctx(COPA, TF), TF)
        ic = NB.inducer_call(call)
        assert ic.source == "neighborhood" and ic.role == "metal"
        assert ic.ligand is None, (
            "a copA neighbour says copper is handled nearby, NOT that this regulator senses copper"
        )

    def test_it_cannot_vote_for_a_headline_ligand(self):
        """`ligand=None` means the source is retained but never becomes `top` (schema.from_calls)."""
        call = NB.analyse(_ctx(COPA, TF), TF)
        cons = IND.InducerConsensus.from_calls([NB.inducer_call(call)])
        assert cons.top is None


class TestPrecedence:
    """It sits at the BOTTOM of the ladder: chemistry first, MetalNet next, co-location last."""

    @staticmethod
    def _nb_metal():
        return NB.inducer_call(NB.analyse(_ctx(COPA, TF), TF))

    def test_may_establish_metal_where_chemistry_is_silent(self):
        cons = IND._finalize([self._nb_metal()], None)
        assert cons.metal or "metal" in (cons.top or "").lower()

    def test_must_not_override_a_negative_family_gate(self):
        """The negative gate exists to catch the non-metal members of metal families (FrmR, TnrA)."""
        cons = IND._finalize([self._nb_metal()], False)
        top = (cons.top or "").lower()
        assert "non-metal" in top or not cons.metal, (
            "a gate-negative TF must stay non-metal however suggestive its neighbours are"
        )

    def test_a_positive_gate_keeps_the_verdict(self):
        cons = IND._finalize([self._nb_metal()], True)
        assert cons.metal or "metal" in (cons.top or "").lower()
