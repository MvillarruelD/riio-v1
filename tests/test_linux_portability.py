"""Linux portability contracts: engine discovery and settings resolve the same way on every OS.

These run on any platform: the native (non-Windows) branches are exercised by switching the module's
`ON_WINDOWS` flag and stubbing the subprocess, so a Windows developer machine still guards the Linux
paths that CI runs for real.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from predictor import resources
from predictor.discovery import bitacora


def test_env_values_keep_spaces_and_drop_quotes_and_comments():
    txt = ('A=1\n'
           'BITACORA_DIR="/home/me/My Tools/bitacora"\n'
           "export MMSEQS='/opt/mm seqs/bin/mmseqs'\n"
           'BLAST_BIN=/opt/blast/bin   # the conda env\n'
           'EMPTY=\n')
    assert resources._env_value(txt, "A") == "1"
    assert resources._env_value(txt, "BITACORA_DIR") == "/home/me/My Tools/bitacora"
    assert resources._env_value(txt, "MMSEQS") == "/opt/mm seqs/bin/mmseqs"
    assert resources._env_value(txt, "BLAST_BIN") == "/opt/blast/bin"
    assert resources._env_value(txt, "EMPTY") == ""
    assert resources._env_value(txt, "MISSING") == ""


def test_every_engine_location_is_read_from_the_env_file():
    for key in ("MMSEQS", "BLAST_BIN", "DATASETS_EXE", "BITACORA_DIR", "BITACORA_BIN", "BITACORA_PY",
                "METALNET_PYTHON", "METALNET2_DIR"):
        assert key in resources._PATH_KEYS, key


def test_dotenv_override_precedence_and_all_engine_keys(monkeypatch, tmp_path):
    settings = tmp_path / "settings.env"
    keys = resources._PATH_KEYS + resources._TOKEN_KEYS
    settings.write_text("\n".join(f'{key}="value with spaces/{key}"' for key in keys), encoding="utf-8")
    for key in keys:
        monkeypatch.delenv(key, raising=False)
    monkeypatch.setenv("PREDICTOR_ENV_FILE", str(settings))
    monkeypatch.setenv("MMSEQS", "explicit/path")
    resources.load_tokens()
    assert __import__("os").environ["MMSEQS"] == "explicit/path"
    assert all(__import__("os").environ[key] == f"value with spaces/{key}" for key in keys if key != "MMSEQS")
    resources.load_tokens(override=True)
    assert __import__("os").environ["MMSEQS"] == "value with spaces/MMSEQS"


@pytest.mark.parametrize("value", ['"unterminated', '"closed" junk', "'unterminated"])
def test_dotenv_malformed_quotes_are_not_engine_paths(value):
    assert resources._env_value(f"BITACORA_DIR={value}\n", "BITACORA_DIR") == ""


def test_dotenv_keeps_windows_trailing_backslash_and_literal_hash():
    assert resources._env_value('BLAST_BIN="C:\\My Tools\\blast\\" # note\n', "BLAST_BIN") == "C:\\My Tools\\blast\\"
    assert resources._env_value('BLAST_BIN="/opt/hash#tools/bin"\n', "BLAST_BIN") == "/opt/hash#tools/bin"


def test_explicit_env_file_beats_checkout_and_user(monkeypatch, tmp_path):
    checkout = tmp_path / "checkout"
    checkout.mkdir()
    explicit, user = tmp_path / "explicit.env", tmp_path / "user.env"
    for path in (explicit, user, checkout / ".env"):
        path.touch()
    monkeypatch.setenv("PREDICTOR_ENV_FILE", str(explicit))
    monkeypatch.setattr(resources, "REPO", checkout)
    monkeypatch.setattr(resources, "USER_ENV_FILE", user)
    assert resources.env_file() == explicit
    monkeypatch.delenv("PREDICTOR_ENV_FILE")
    assert resources.env_file() == checkout / ".env"


def test_patch_check_is_read_only_and_repair_is_idempotent(monkeypatch, tmp_path):
    tracked, engine = tmp_path / "patches", tmp_path / "engine"
    tracked.mkdir()
    (engine / "Scripts").mkdir(parents=True)
    (tracked / "coverage.pl").write_text("new coverage rule\n", encoding="utf-8")
    target = engine / "Scripts" / "coverage.pl"
    target.write_text("old coverage rule\n", encoding="utf-8")
    monkeypatch.setattr(bitacora, "PATCH_DIR", tracked)
    monkeypatch.setattr(bitacora, "ENGINE_DIR", engine)
    assert bitacora.patch_status()[0] is False
    assert target.read_text() == "old coverage rule\n"
    assert bitacora.ensure_patched() == ["coverage.pl"]
    assert bitacora.patch_status()[0] is True
    assert bitacora.ensure_patched() == []


def test_unpublished_workspace_reads_require_opt_in(monkeypatch, tmp_path):
    from predictor.annotate import family_db, genome_db, genome_resolver, known_operators, tf_record
    monkeypatch.delenv("PREDICTOR_ALLOW_LEGACY_WORKSPACE", raising=False)
    monkeypatch.setattr(resources, "env_file", lambda: None)
    # Any glob of the neighbouring corpus would fail; generic tests never read it.
    monkeypatch.setattr(family_db.glob, "glob", lambda *a, **kw: pytest.fail("unpublished workspace read"))
    monkeypatch.setattr(known_operators.family_db, "_ALL_JSON", tmp_path / "absent.json")
    assert list(family_db._iter_tf_yamls()) == []
    assert genome_db.curated_genome_fastas() == []
    assert genome_resolver._env_file() is None
    assert known_operators.build_index(tmp_path / "index.json") == {"by_family": {}, "by_tf": {}}
    with pytest.raises(RuntimeError, match="PREDICTOR_ALLOW_LEGACY_WORKSPACE"):
        tf_record.build_reference(tmp_path / "missing.faa", tmp_path / "missing.json")
    monkeypatch.setenv("PREDICTOR_ALLOW_LEGACY_WORKSPACE", "true")
    assert resources.legacy_workspace_enabled()


def test_user_env_file_is_a_fallback_location(monkeypatch, tmp_path):
    user = tmp_path / ".env"
    user.write_text("NCBI_EMAIL=someone@example.org\n", encoding="utf-8")
    monkeypatch.delenv("PREDICTOR_ENV_FILE", raising=False)
    monkeypatch.setattr(resources, "REPO", tmp_path / "no_checkout")
    monkeypatch.setattr(resources, "USER_ENV_FILE", user)
    assert resources.env_file() == user


def test_bitacora_bin_goes_first_on_path_natively(monkeypatch, tmp_path):
    monkeypatch.setattr(bitacora, "ON_WINDOWS", False)
    monkeypatch.setenv("BITACORA_BIN", str(tmp_path))
    env = bitacora.native_env()
    assert env["PATH"].split(__import__("os").pathsep)[0] == str(tmp_path)
    assert bitacora._wsl_exe() is None, "Linux must never route through wsl.exe"


def test_missing_bitacora_bin_falls_back_to_path(monkeypatch, tmp_path):
    monkeypatch.setattr(bitacora, "ON_WINDOWS", False)
    monkeypatch.setenv("BITACORA_BIN", str(tmp_path / "absent"))
    assert bitacora.native_bin() is None
    assert bitacora.native_env()["PATH"] == __import__("os").environ.get("PATH", "")


def test_clade_build_uses_a_private_scratch_dir_and_reports_failure(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(bitacora, "ON_WINDOWS", False)
    (tmp_path / "P_db.fasta").write_text(">a\nMKV\n", encoding="utf-8")
    seen = {}

    def fake_run(cmd, **kw):
        seen["cmd"], seen["env"] = cmd, kw.get("env")
        return subprocess.CompletedProcess(cmd, 1, "", "hmmalign: command not found")

    monkeypatch.setattr(bitacora.subprocess, "run", fake_run)
    assert bitacora._build_clade_hmm("P", tmp_path / "P_db.hmm", tmp_path) == ""
    script = " ".join(seen["cmd"])
    assert "/tmp/_tfop_clade" not in script, "fixed scratch names collide between concurrent runs"
    assert "tfop_clade_" in script
    assert seen["env"] is not None and "PATH" in seen["env"]
    assert "clade HMM build failed" in capsys.readouterr().err


def test_fallback_profiles_are_tagged_in_candidate_provenance(tmp_path):
    (tmp_path / "TetR_SczA_db.hmm.source").write_text("pfam_family", encoding="utf-8")
    (tmp_path / "ecFur_db.hmm.source").write_text("clade", encoding="utf-8")
    cands = [bitacora.BitacoraCandidate("a", "MKV", profile="TetR_SczA"),
             bitacora.BitacoraCandidate("b", "MKV", profile="ecFur"),
             bitacora.BitacoraCandidate("c", "MKV", profile="ecFur|TetR_SczA")]
    assert bitacora.tag_fallback_profiles(cands, tmp_path) == 2
    assert cands[0].profile == "TetR_SczA" + bitacora.PFAM_FALLBACK_TAG
    assert cands[1].profile == "ecFur"
    assert cands[2].profile == "ecFur|TetR_SczA" + bitacora.PFAM_FALLBACK_TAG


def test_cached_fallback_models_are_retried(tmp_path):
    (tmp_path / "X_db.hmm").write_text("HMMER3\n//\n", encoding="utf-8")
    (tmp_path / "X_db.hmm.source").write_text("pfam_family", encoding="utf-8")
    assert bitacora.hmm_source(tmp_path, "X") == "pfam_family"
    assert bitacora.hmm_source(tmp_path, "never_built") == ""


def test_regulon_pool_spawns_rather_than_forks():
    src = (Path(bitacora.__file__).resolve().parents[1] / "report" / "ligand_regulon.py").read_text(
        encoding="utf-8")
    assert 'get_context("spawn")' in src


def test_patch_status_after_equal_size_equal_mtime_repair(monkeypatch, tmp_path):
    import filecmp
    import os
    source = tmp_path / "patches"
    engine = tmp_path / "engine"
    source.mkdir()
    (engine / "Scripts").mkdir(parents=True)
    original = source / "patch.pl"
    installed = engine / "Scripts" / original.name
    original.write_text("good1")
    installed.write_text("bad11")
    os.utime(installed, ns=(original.stat().st_atime_ns, original.stat().st_mtime_ns))
    filecmp.clear_cache()
    assert not filecmp.cmp(original, installed, shallow=False)
    monkeypatch.setattr(bitacora, "PATCH_DIR", source)
    monkeypatch.setattr(bitacora, "ENGINE_DIR", engine)
    assert bitacora.ensure_patched() == [original.name]
    assert bitacora.patch_status()[0]


def test_patch_status_reports_unreadable_installed_script(monkeypatch, tmp_path):
    source = tmp_path / "patches"
    engine = tmp_path / "engine"
    source.mkdir()
    (engine / "Scripts").mkdir(parents=True)
    original = source / "patch.pl"
    installed = engine / "Scripts" / original.name
    original.write_text("good1")
    installed.write_text("good1")
    monkeypatch.setattr(bitacora, "PATCH_DIR", source)
    monkeypatch.setattr(bitacora, "ENGINE_DIR", engine)
    read = Path.read_bytes

    def unreadable(path):
        if path == installed:
            raise PermissionError("script is locked")
        return read(path)

    monkeypatch.setattr(Path, "read_bytes", unreadable)
    assert not bitacora.patch_status()[0]
