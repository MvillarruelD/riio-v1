"""The ColabFold MSA expansion reads its cache first, retries transient failures, and records a failure.

Measured on the 2026-11-10 and 2026-11-26 runs of the 238 canonical candidates: 6 ColabFold fetches
failed (read timeouts and refused connections), each leaving a sparse regulator on its 4-6 local homolog
promoters instead of about 65, and 3 of the 6 moved the primary operator. In one of them (E. coli GntR
YP_026245.1) the alignment and its 60-promoter record were already on disk. Nothing in the dossier or
the run record said an expansion had been attempted and failed.
"""
from __future__ import annotations

import io
import json
import tarfile
import urllib.error

import pytest

from predictor.annotate import msa_homologs as mh
from predictor.report.render_report import build_html

PROTEIN = "MKTLSELLAQHGIPVTRRQALVLAELARLGHASIAELSERLGLSRPTVSRLIDELEARG"
OTHER = "MSNIQEAAKRLGLSPRTVRRWEEEGLIKPERSESGYRLYSEEDIERLLEIKRLRELGL"
REGION = "ACGTTGACAATTAATCATCGGCTCGTATAATGTGTGGAATTGTGAGCGGATAACAATT"


def _a3m(query: str) -> str:
    return f">101\n{query}\n>UniRef90_P0A9E0\n{query}\n"


def _tar(members: dict) -> bytes:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tf:
        for name, text in members.items():
            data = text.encode()
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tf.addfile(info, io.BytesIO(data))
    return buf.getvalue()


class _Resp:
    def __init__(self, body: bytes):
        self._body = body

    def read(self):
        return self._body

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


class FakeServer:
    """urlopen stand-in for api.colabfold.com. `script[step]` is a list consumed one item per request:
    an exception instance is raised, anything else is the response."""

    def __init__(self, **script):
        self.script = {k: list(v) for k, v in script.items()}
        self.calls = {"submit": 0, "poll": 0, "download": 0}

    def __call__(self, req, timeout=None):
        url = req if isinstance(req, str) else req.full_url
        step = ("submit" if url.endswith("/ticket/msa") else
                "download" if "/result/download/" in url else "poll")
        self.calls[step] += 1
        queue = self.script[step]
        item = queue.pop(0) if len(queue) > 1 else queue[0]
        if isinstance(item, BaseException):
            raise item
        return _Resp(item if isinstance(item, bytes) else json.dumps(item).encode())


@pytest.fixture
def no_sleep(monkeypatch):
    waits = []
    monkeypatch.setattr(mh.time, "sleep", lambda s: waits.append(s))
    return waits


# ----------------------------------------------------------------------------------- cache first
def test_a_cached_alignment_is_used_before_the_network(tmp_path, monkeypatch):
    monkeypatch.setattr(mh, "CACHE", tmp_path)
    monkeypatch.setattr(mh, "_op_msa", lambda tf_id: None)
    (tmp_path / "TF1_colabfold_uniref.a3m").write_text(_a3m(PROTEIN))

    def network(*a, **k):
        raise AssertionError("the server must not be asked for an alignment that is on disk")

    monkeypatch.setattr(mh, "colabfold_a3m", network)
    monkeypatch.setattr(mh, "collect_msa_regions",
                        lambda a3m, **k: {"key": "k", "regions": [REGION], "n_orthologs": 1,
                                          "n_mapped": 1, "n_loci": 1, "sources": [{}]})
    regions, info = mh.expand_homolog_regions(["AAAA"], tf_id="TF1", seq=PROTEIN,
                                              allow_colabfold=True, verbose=False)
    assert regions == ["AAAA", REGION]
    assert info["source"] == "colabfold" and info["alignment_cached"] is True
    assert "failure" not in info


def test_an_alignment_of_another_protein_is_not_reused(tmp_path, monkeypatch):
    monkeypatch.setattr(mh, "CACHE", tmp_path)
    monkeypatch.setattr(mh, "_op_msa", lambda tf_id: None)
    (tmp_path / "TF1_colabfold_uniref.a3m").write_text(_a3m(OTHER))
    asked = []
    monkeypatch.setattr(mh, "colabfold_a3m", lambda *a, **k: asked.append(1) or None)
    mh.expand_homolog_regions(["AAAA"], tf_id="TF1", seq=PROTEIN, allow_colabfold=True, verbose=False)
    assert asked, "a file under this name whose query row is another protein is not this protein's MSA"


def test_a_rich_local_set_makes_no_request_and_records_no_failure(monkeypatch):
    monkeypatch.setattr(mh, "colabfold_a3m", lambda *a, **k: pytest.fail("no expansion was needed"))
    regions, info = mh.expand_homolog_regions(["A"] * 8, tf_id="TF1", seq=PROTEIN,
                                              allow_colabfold=True, verbose=False)
    assert len(regions) == 8 and "failure" not in info


# ----------------------------------------------------------------------------------- retries
def test_transient_failures_are_retried_at_every_step(tmp_path, monkeypatch, no_sleep):
    server = FakeServer(
        submit=[urllib.error.URLError("[WinError 10060] connection attempt failed"),
                {"status": "RATELIMIT"}, {"id": "T1", "status": "PENDING"}],
        poll=[TimeoutError("The read operation timed out"), {"status": "COMPLETE"}],
        download=[TimeoutError("The read operation timed out"),
                  _tar({"uniref.a3m": _a3m(PROTEIN), "bfd.mgnify30.a3m": _a3m(PROTEIN) * 3})])
    monkeypatch.setattr(mh.urllib.request, "urlopen", server)
    failure: dict = {}
    out = mh.colabfold_a3m(PROTEIN, tmp_path / "x.a3m", prefer="uniref", verbose=False, failure=failure)
    assert out == tmp_path / "x.a3m" and out.read_text() == _a3m(PROTEIN), "the UniRef member, written"
    assert server.calls == {"submit": 3, "poll": 2, "download": 2}
    assert failure == {}, "a fetch that succeeded on retry records no failure"


def test_a_persistent_failure_is_recorded_with_its_stage(tmp_path, monkeypatch, no_sleep):
    server = FakeServer(submit=[TimeoutError("The read operation timed out")])
    monkeypatch.setattr(mh.urllib.request, "urlopen", server)
    failure: dict = {}
    assert mh.colabfold_a3m(PROTEIN, tmp_path / "x.a3m", verbose=False, failure=failure) is None
    assert failure["stage"] == "submit" and failure["attempts"] == mh._CF_ATTEMPTS
    assert failure["error"].startswith("TimeoutError")
    assert no_sleep == [20, 40, 60, 80], "linear back-off between the five attempts"


def test_download_retries_are_bounded(tmp_path, monkeypatch, no_sleep):
    server = FakeServer(submit=[{"id": "T1", "status": "COMPLETE"}], poll=[{"status": "COMPLETE"}],
                        download=[TimeoutError("The read operation timed out")])
    monkeypatch.setattr(mh.urllib.request, "urlopen", server)
    failure: dict = {}
    assert mh.colabfold_a3m(PROTEIN, tmp_path / "x.a3m", verbose=False, failure=failure) is None
    assert server.calls["download"] == mh._CF_DOWNLOAD_ATTEMPTS == failure["attempts"]
    assert failure["stage"] == "download"


def test_a_client_error_is_not_retried(tmp_path, monkeypatch, no_sleep):
    server = FakeServer(submit=[urllib.error.HTTPError("u", 400, "bad request", {}, None)])
    monkeypatch.setattr(mh.urllib.request, "urlopen", server)
    failure: dict = {}
    assert mh.colabfold_a3m(PROTEIN, tmp_path / "x.a3m", verbose=False, failure=failure) is None
    assert server.calls["submit"] == 1 and failure["attempts"] == 1 and not no_sleep


def test_a_server_status_failure_is_recorded(tmp_path, monkeypatch, no_sleep):
    server = FakeServer(submit=[{"id": "T1", "status": "PENDING"}], poll=[{"status": "ERROR"}])
    monkeypatch.setattr(mh.urllib.request, "urlopen", server)
    failure: dict = {}
    assert mh.colabfold_a3m(PROTEIN, tmp_path / "x.a3m", verbose=False, failure=failure) is None
    assert failure["stage"] == "poll" and failure["error"] == "status=ERROR"


# ----------------------------------------------------------------------------------- recorded
def test_a_failed_expansion_says_so_in_its_record(tmp_path, monkeypatch, no_sleep):
    monkeypatch.setattr(mh, "CACHE", tmp_path)
    monkeypatch.setattr(mh, "_op_msa", lambda tf_id: None)
    monkeypatch.setattr(mh.urllib.request, "urlopen",
                        FakeServer(submit=[{"id": "T1", "status": "COMPLETE"}],
                                   poll=[{"status": "COMPLETE"}],
                                   download=[TimeoutError("The read operation timed out")]))
    regions, info = mh.expand_homolog_regions(["A"] * 6, tf_id="TF1", seq=PROTEIN,
                                              allow_colabfold=True, verbose=False)
    assert regions == ["A"] * 6 and info["source"] is None and info["msa_added"] == 0
    assert info["failure"] == {"source": "colabfold", "stage": "download",
                               "attempts": mh._CF_DOWNLOAD_ATTEMPTS,
                               "error": "TimeoutError: The read operation timed out"}


def _report(msa: dict, homolog_regions: dict | None = None) -> str:
    return build_html({
        "dossier": {"tf_id": "T", "family": "GntR",
                    "homolog_regions": {"pre_msa": [], "expanded": [], **(homolog_regions or {})}},
        "run_record": {"homolog_collection": {"attempts": [{"route": "blast_similarity", "selected": True}],
                                              "msa": msa}},
        "bundle_files": set(),
    }, {})


def test_the_report_states_a_failed_expansion():
    failure = {"source": "colabfold", "stage": "download", "attempts": 3,
               "error": "TimeoutError: The read operation timed out"}
    html = _report({"local": 6, "msa_added": 0, "source": None, "failure": failure},
                   {"expansion_failure": failure})
    assert "MSA expansion failed" in html and "download" in html and "TimeoutError" in html


def test_the_report_is_silent_when_no_expansion_was_needed():
    assert "MSA expansion failed" not in _report({"local": 12, "msa_added": 0, "source": None})
