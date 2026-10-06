"""The report must distinguish an autoregulation call from a distal operator prediction.

Measured on the RegulonDB panel (HANDOFF 13.2): 12 of 13 primary-site recoveries were tier-1
proximal, and recovery given a distal top hit was 1/16. A single recovery number therefore conflates
"the motif found a site" with "the locality prior pointed at the regulator's own promoter". These
tests pin the split the report now reports.
"""
from __future__ import annotations

from predictor.report import operators as OPS


class TestLocalityOf:
    def test_tier1_is_proximal(self):
        assert OPS.locality_of("rescan_tier1") == "proximal"

    def test_tiers_beyond_one_are_distal(self):
        assert OPS.locality_of("rescan_tier2") == "distal"
        assert OPS.locality_of("rescan_tier3") == "distal"
        assert OPS.locality_of("rescan_tier4") == "distal"

    def test_untagged_is_neither(self):
        """Pre-refactor bundles carry no generator; counting them as distal would credit them with
        predictions they never made."""
        assert OPS.locality_of(None) == "unknown"
        assert OPS.locality_of("") == "unknown"
        assert OPS.locality_of("some_other_generator") == "unknown"


class TestLocalitySummary:
    @staticmethod
    def _dossier(*generators_and_scores):
        return {"rescan": {"hits": [{"generator": g, "score": s}
                                    for g, s in generators_and_scores]}}

    def test_counts_partition_the_hits(self):
        d = self._dossier(("rescan_tier1", 5.0), ("rescan_tier3", 9.0), (None, 1.0))
        s = OPS.locality_summary(d)
        assert s["proximal"] + s["distal"] + s["unknown"] == s["n"] == 3
        assert (s["proximal"], s["distal"], s["unknown"]) == (1, 1, 1)

    def test_primary_locality_follows_the_top_SCORING_hit_not_the_first(self):
        d = self._dossier(("rescan_tier1", 2.0), ("rescan_tier3", 9.0))
        assert OPS.locality_summary(d)["primary_locality"] == "distal"

    def test_empty_bundle_is_unknown_not_a_crash(self):
        s = OPS.locality_summary({})
        assert s["n"] == 0 and s["primary_locality"] == "unknown"


def test_in_tf_neighborhood_is_not_a_synonym_for_proximal():
    """HANDOFF 13.4: current bundles carry tier-3 hits that are still flagged in-neighborhood, so the
    looser flag must not be read as the locality claim."""
    d = {"rescan": {"hits": [{"generator": "rescan_tier3", "score": 9.0, "seq": "ACGTACGT",
                              "dyad": 100, "accession": "NC_1", "start": 90, "end": 110,
                              "strand": "+"}]},
         "neighborhood": {"window": [0, 1000]}}
    cands = OPS._natural_candidates(d)
    assert cands, "the hit should still be reported"
    prov = cands[0].provenance
    assert prov["in_tf_neighborhood"] is True, "it IS inside the looser window"
    assert prov["locality"] == "distal", "but its tier says distal, and the tier is what counts"
