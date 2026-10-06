"""The BITACORA adapter's seams, none of which the offline self-test could reach.

`find_candidates` was shipped passing `--genome` to a runner whose only required argument is
`--proteins`, so it failed 100 % of the time with the engine PRESENT while the self-test -- which
only ever exercises the engine-ABSENT path -- stayed green. Everything here tests the part that runs
when the engine exists, with the engine itself stubbed out.
"""
from __future__ import annotations

import subprocess
import types

import pytest

from predictor.discovery import bitacora


def _stub_genome_scan(monkeypatch, recs, kind):
    """Replace `genome_scan` where `_proteome_for` actually looks it up.

    Patching `sys.modules` is not enough: `from predictor.annotate import genome_scan` reads the
    ATTRIBUTE off the already-imported package, so a sys.modules stub is bypassed whenever anything
    earlier in the suite has imported it -- these tests passed alone and failed in the full run.
    """
    import predictor.annotate as _pkg

    monkeypatch.setattr(_pkg, "genome_scan",
                        types.SimpleNamespace(proteins_from_input=lambda *a, **k: (recs, kind)),
                        raising=False)


# --------------------------------------------------------------------- profile -> family
class TestFamilyOfProfiles:
    def test_single_profile(self):
        assert bitacora._family_of_profiles("BsCzrA") == "ArsR"

    def test_joined_profiles_keep_their_family(self):
        """The runner writes `A|B` when two seed profiles hit one protein.

        A bare dict lookup returns "" for those, dropping the family of exactly the proteins the
        search was most confident about -- two profiles agreeing is its strongest evidence.
        """
        assert bitacora._family_of_profiles("BsCzrA|mtNmtR") == "ArsR"
        assert bitacora._family_of_profiles("PbrR|ecZntR") == "MerR"

    def test_unknown_and_empty_are_blank_not_an_error(self):
        assert bitacora._family_of_profiles("") == ""
        assert bitacora._family_of_profiles("NotAProfile") == ""

    def test_first_known_family_wins_when_profiles_disagree(self):
        assert bitacora._family_of_profiles("NotAProfile|ecFur") == "Fur"


# --------------------------------------------------------------------- proteome derivation
class TestProteomeFor:
    def test_writes_full_length_records_with_unix_newlines(self, tmp_path, monkeypatch):
        recs = [("NP_000001.1", "MKTAYIAKQR"), ("NP_000002.1", "MSDQEAKPST")]
        _stub_genome_scan(monkeypatch, recs, "annotated_genome")

        path, kind = bitacora._proteome_for("g.fna", None, tmp_path)
        assert kind == "annotated_genome"
        raw = path.read_bytes()
        # BLAST+/HMMER read this inside WSL; CRLF makes the last residue of every line part of the
        # alphabet as far as they are concerned.
        assert b"\r\n" not in raw
        assert raw.decode() == ">NP_000001.1\nMKTAYIAKQR\n>NP_000002.1\nMSDQEAKPST\n"

    def test_empty_input_returns_none_rather_than_an_empty_fasta(self, tmp_path, monkeypatch):
        """Handing BITACORA an empty FASTA makes it report a failure that names nothing useful."""
        _stub_genome_scan(monkeypatch, [], "genome")
        path, kind = bitacora._proteome_for("g.fna", None, tmp_path)
        assert path is None and kind == ""

    def test_records_missing_a_field_do_not_produce_a_blank_record(self, tmp_path, monkeypatch):
        _stub_genome_scan(monkeypatch, [("NP_1", ""), ("", "MKT"), ("NP_2", "MKTA")], "proteome")
        path, _ = bitacora._proteome_for("g.fna", None, tmp_path)
        assert path.read_text() == ">NP_2\nMKTA\n"


# --------------------------------------------------------------------- the command it builds
class TestFindCandidatesCommand:
    """The regression that matters: the runner REQUIRES --proteins and ignores --genome."""

    @pytest.fixture
    def stubbed(self, tmp_path, monkeypatch):
        monkeypatch.setattr(bitacora, "engine_available", lambda: (True, ""))
        prof = tmp_path / "profiles"
        prof.mkdir()
        (prof / "BsCzrA_db.fasta").write_text(">x\nMKT\n")
        monkeypatch.setattr(bitacora, "_proteome_for",
                            lambda *a, **k: (tmp_path / "proteins.faa", "annotated_genome"))
        seen = {}

        def fake_run(cmd, **kw):
            seen["cmd"] = cmd
            (tmp_path / "out").mkdir(exist_ok=True)
            (tmp_path / "out" / "candidates.tsv").write_text(
                "protein_id\tprofile\tsequence\ttrim_start\ttrim_end\tsearch\n"
                "NP_1\tBsCzrA|mtNmtR\tMKTAYIAKQR\t2\t8\tbitacora_protein_mode\n")
            return subprocess.CompletedProcess(cmd, 0, "", "")

        monkeypatch.setattr(bitacora.subprocess, "run", fake_run)
        return tmp_path, prof, seen

    def test_passes_proteins_and_not_genome(self, stubbed):
        tmp_path, prof, seen = stubbed
        res = bitacora.find_candidates(tmp_path / "g.fna", profiles=prof, out_dir=tmp_path / "out")
        assert res.status == "ok", res.reason
        cmd = seen["cmd"]
        assert "--proteins" in cmd, "the runner's only required argument"
        assert "--genome" not in cmd, "protein mode ignores it; passing it alone was the bug"

    def test_invokes_a_real_interpreter_not_the_bare_name_python(self, stubbed):
        """`python` is not on PATH under that name here; it resolved to a Windows Store shim."""
        tmp_path, prof, seen = stubbed
        bitacora.find_candidates(tmp_path / "g.fna", profiles=prof, out_dir=tmp_path / "out")
        assert seen["cmd"][0] != "python"

    def test_carries_protein_source_onto_the_result(self, stubbed):
        tmp_path, prof, _ = stubbed
        res = bitacora.find_candidates(tmp_path / "g.fna", profiles=prof, out_dir=tmp_path / "out")
        # 'genome' would mean genes were CALLED, so the ids are ours and join to no external table.
        assert res.protein_source == "annotated_genome"

    def test_multi_profile_hit_keeps_its_family_end_to_end(self, stubbed):
        tmp_path, prof, _ = stubbed
        res = bitacora.find_candidates(tmp_path / "g.fna", profiles=prof, out_dir=tmp_path / "out")
        (c,) = res.candidates
        assert c.family == "ArsR"
        assert c.sequence == "MKTAYIAKQR"          # FULL length, never the trim
        assert c.trimmed_sequence == "TAYIAK"      # the trim is an annotation on it

    def test_clean_install_materializes_packaged_profiles(self, tmp_path, monkeypatch):
        """A user must not know that BITACORA profiles are derived cache files."""
        monkeypatch.setattr(bitacora, "engine_available", lambda: (True, ""))
        profiles = tmp_path / "new-profile-cache"
        built = []

        def fake_build(dest):
            built.append(dest)
            dest.mkdir(parents=True)
            (dest / "BsCzrA_db.fasta").write_text(">x\nMKT\n")
            return {"BsCzrA": {"status": "written"}}

        monkeypatch.setattr(bitacora, "build_profiles", fake_build)
        monkeypatch.setattr(bitacora, "_proteome_for",
                            lambda *a, **k: (tmp_path / "proteins.faa", "annotated_genome"))

        def fake_run(cmd, **kw):
            out = tmp_path / "out"
            out.mkdir(exist_ok=True)
            (out / "candidates.tsv").write_text(
                "protein_id\tprofile\tsequence\ttrim_start\ttrim_end\tsearch\n"
                "NP_1\tBsCzrA\tMKT\t1\t3\tbitacora_protein_mode\n")
            return subprocess.CompletedProcess(cmd, 0, "", "")

        monkeypatch.setattr(bitacora.subprocess, "run", fake_run)
        result = bitacora.find_candidates(
            tmp_path / "genome.fna", profiles=profiles, out_dir=tmp_path / "out")
        assert built == [profiles]
        assert result.status == "ok"
