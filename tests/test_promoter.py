"""Promoter positioning must be available in an ordinary package installation."""
from predictor.signals import promoter


def test_promoter_calculator_runtime_dependency_is_wired():
    sequence = (
        "GGGCGCGAACT" + "TTGACA" + "GCTAGCATCGATCGAT" + "TATAAT"
        + "GCATACTGGGCATGCATGCATGCGGGCCCAAATTTGGGACGT" * 2
    )

    promoters = promoter.find_promoters(sequence)
    assert promoters
    call = promoter.classify_mode(sequence, 18, 30)
    assert call.mode in {
        "elongated-spacer-overlap", "promoter-core-overlap", "no-promoter-overlap", "ambiguous",
    }
    assert call.promoter
