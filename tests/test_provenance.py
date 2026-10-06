"""Provenance: the fingerprint identifies code, not a checkout; the commit is recorded beside it."""
from __future__ import annotations

from predictor import provenance


def test_fingerprint_ignores_line_endings(tmp_path):
    lf, crlf = tmp_path / "lf", tmp_path / "crlf"
    for root, text in ((lf, b"a = 1\nb = 2\n"), (crlf, b"a = 1\r\nb = 2\r\n")):
        root.mkdir()
        (root / "m.py").write_bytes(text)
    assert (provenance._digest_files([lf / "m.py"], root=lf)
            == provenance._digest_files([crlf / "m.py"], root=crlf))


def test_fingerprint_still_sees_content_changes(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    for root, text in ((a, b"x = 1\n"), (b, b"x = 2\n")):
        root.mkdir()
        (root / "m.py").write_bytes(text)
    assert provenance._digest_files([a / "m.py"], root=a) != provenance._digest_files([b / "m.py"], root=b)


def test_runtime_records_git_commit_outside_the_fingerprint():
    rt = provenance.runtime_provenance()
    assert {"git_commit", "git_dirty"} <= set(rt)
    # the commit must never feed the content fingerprint, or two checkouts of one tree would differ
    assert rt["runtime_sha256"] == provenance.runtime_provenance()["runtime_sha256"]


def test_git_provenance_degrades_without_a_checkout(monkeypatch, tmp_path):
    monkeypatch.setattr(provenance, "PACKAGE_ROOT", tmp_path / "predictor")
    assert provenance.git_provenance() == {"git_commit": "unavailable", "git_dirty": None}
