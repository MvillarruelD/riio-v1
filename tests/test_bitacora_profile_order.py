"""One canonical profile/trim annotation for identical raw engine hits on every OS."""
from predictor._engines.bitacora_runner import collect


def test_profile_order_and_first_trim_follow_historical_windows_order(tmp_path):
    # Create in the opposite order too; directory enumeration must not select the annotation.
    for profile, matched in (("PbrR", "MKLVTAAG"), ("ecZntR", "MKLV")):
        directory = tmp_path / profile
        directory.mkdir()
        (directory / f"{profile}_proteins_trimmed.fasta").write_text(f">protein\n{matched}\n")

    result = collect(tmp_path, {"protein": "MKLVTAAG"})
    assert len(result) == 1
    row = result[0]
    assert row["sequence"] == "MKLVTAAG"
    assert row["profiles"] == ["ecZntR", "PbrR"]
    assert (row["trim_start"], row["trim_end"]) == (0, 4)


def test_first_valid_trim_is_used_when_preceding_profile_has_no_matching_span(tmp_path):
    for profile, matched in (("ecZntR", "WWWW"), ("PbrR", "LVTA")):
        directory = tmp_path / profile
        directory.mkdir()
        (directory / f"{profile}_proteins_trimmed.fasta").write_text(f">protein\n{matched}\n")

    row = collect(tmp_path, {"protein": "MKLVTAAG"})[0]
    assert row["profiles"] == ["ecZntR", "PbrR"]
    assert row["sequence"] == "MKLVTAAG"
    assert (row["trim_start"], row["trim_end"]) == (2, 6)
