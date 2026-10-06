"""Unit contracts for `analysis/manifest.py`, the shared manifest gate.

The orchestration tests exercise these rules through the driver; these pin them at the module, so a
stage that calls the module directly is held to the same contract.
"""
from __future__ import annotations

import pytest

import manifest as MF

HEADER = "run_name,fasta,genome,organism_acc,family,family_flag\n"
GOOD_ROW = "cand_1,tf.fasta,Genome sp.,GCF_000005845.2,MerR,(auto)\n"


@pytest.fixture
def manifest_file(tmp_path, protein_fasta):
    protein_fasta("tf")

    def _write(body: str, header: str = HEADER):
        path = tmp_path / "run_manifest_t.csv"
        path.write_text(header + body, encoding="utf-8")
        return path
    return _write


def test_a_well_formed_manifest_has_no_findings(manifest_file, tmp_path):
    assert MF.validate(manifest_file(GOOD_ROW), tmp_path) == []


def test_a_missing_file_is_reported_not_raised(tmp_path):
    issues = MF.validate(tmp_path / "absent.csv", tmp_path)
    assert issues and "does not exist" in issues[0]


def test_missing_columns_are_named(manifest_file, tmp_path):
    path = manifest_file("cand_1,tf.fasta\n", header="run_name,fasta\n")
    issues = MF.validate(path, tmp_path)
    assert len(issues) == 1
    for column in ("organism_acc", "family", "family_flag"):
        assert column in issues[0]


def test_an_empty_manifest_is_not_executable(manifest_file, tmp_path):
    assert MF.validate(manifest_file(""), tmp_path) == ["manifest has no candidates"]


def test_an_empty_organism_accession_is_caught(manifest_file, tmp_path):
    path = manifest_file("cand_1,tf.fasta,Genome sp.,,MerR,(auto)\n")
    assert any("empty organism_acc" in i for i in MF.validate(path, tmp_path))


def test_an_absolute_fasta_path_is_rejected_as_unportable(manifest_file, tmp_path):
    abs_fasta = (tmp_path / "tf.fasta").as_posix()
    path = manifest_file(f"cand_1,{abs_fasta},Genome sp.,GCF_000005845.2,MerR,(auto)\n")
    assert any("repository-relative" in i for i in MF.validate(path, tmp_path))


@pytest.mark.parametrize("name, reason", [
    ("", "empty run_name"),
    ("has/slash", "illegal in a Windows path"),
    ("has:colon", "illegal in a Windows path"),
    ("trailing.", "dot or space"),
    ("CON", "reserved filesystem name"),
    ("x" * 200, "too long"),
])
def test_run_names_that_cannot_be_directories_are_rejected(name, reason):
    why = MF.unsafe_run_name(name)
    assert why is not None and reason in why


def test_ordinary_run_names_are_accepted():
    assert MF.unsafe_run_name("MtubH37Rv__ArsR__NP_214595.1") is None


def test_load_or_exit_raises_system_exit_2_with_every_finding(manifest_file, tmp_path, capsys):
    path = manifest_file("cand_1,missing.fasta,Genome sp.,,MerR,(auto)\n")
    with pytest.raises(SystemExit) as exc:
        MF.load_or_exit(path, tmp_path)
    assert exc.value.code == 2
    out = capsys.readouterr().out
    assert "empty organism_acc" in out and "FASTA missing" in out


def test_load_or_exit_returns_rows_when_valid(manifest_file, tmp_path):
    rows = MF.load_or_exit(manifest_file(GOOD_ROW), tmp_path)
    assert [r["run_name"] for r in rows] == ["cand_1"]


def test_structure_only_mode_skips_sequence_validation(manifest_file, tmp_path):
    """`--no-sequences` must still catch structural faults, and must not open the FASTAs."""
    bad = tmp_path / "nuc.fasta"
    bad.write_text(">n\n" + "ATGC" * 40 + "\n", encoding="utf-8")
    path = manifest_file(f"cand_1,{bad.name},Genome sp.,GCF_000005845.2,MerR,(auto)\n")
    assert MF.validate(path, tmp_path, check_sequences=False) == []
    assert any("nucleotide" in i for i in MF.validate(path, tmp_path, check_sequences=True))


@pytest.mark.released_data
def test_the_shipped_production_manifest_validates(workspace, analysis_dir):
    """The list the production run will actually execute."""
    path = analysis_dir / "run_manifest_survey_rerun20260924.csv"
    if not path.is_file():
        pytest.skip("final production manifest is not present in this checkout")
    assert MF.validate(path, workspace) == []
