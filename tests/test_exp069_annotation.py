"""Mock-only phase-4 mapping, immutable cache, and snapshot provenance checks."""
from copy import deepcopy
from dataclasses import asdict
from email.message import Message
import hashlib
import inspect
import io
import json
from types import SimpleNamespace as NS
import urllib.error
import urllib.parse

import pytest

from predictor.effector import annotation as A, inducer as IND, ligand as L, operon_chemistry as OC


def entry(primary="P00001", xref="NP_000001.2", *, taxon=123, sequence="MTEST"):
    return {"primaryAccession": primary, "secondaryAccessions": [],
            "organism": {"taxonId": taxon}, "sequence": {"value": sequence},
            "entryAudit": {"entryVersion": 7, "sequenceVersion": 3},
            "uniProtKBCrossReferences": [{"database": "RefSeq", "id": xref,
                                          "properties": [{"key": "ProteinId", "value": "preserved"}],
                                          "evidences": [{"evidenceCode": "synthetic"}]}],
            "comments": [{"commentType": "CATALYTIC ACTIVITY", "reaction": {
                "name": "citrate = isocitrate", "reactionCrossReferences": [
                    {"database": "Rhea", "id": "RHEA:synthetic"}]}}]}


class Response:
    def __init__(self, payload, *, headers=None, status=200, final_url=None):
        self.body = payload if isinstance(payload, bytes) else json.dumps(payload).encode()
        self.headers = Message()
        for key, value in (headers or {}).items():
            self.headers[key] = str(value)
        self.status, self.final_url = status, final_url

    def read(self):
        return self.body

    def geturl(self):
        return self.final_url

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


@pytest.fixture(autouse=True)
def forbid_network(monkeypatch, tmp_path):
    monkeypatch.setattr(A.urllib.request, "urlopen", lambda *a, **k: pytest.fail("unexpected network"))
    monkeypatch.setenv("PREDICTOR_EFFECTOR_CACHE", str(tmp_path / "cache"))
    monkeypatch.setattr(L, "_CACHE_DIR", tmp_path / "cache")


def transport(monkeypatch, *pages):
    calls = []
    pages = iter(pages)

    def fetch(url, *, timeout):
        calls.append((url, timeout))
        result = next(pages)
        if isinstance(result, BaseException):
            raise result
        result.final_url = result.final_url or url
        return result

    monkeypatch.setattr(A.urllib.request, "urlopen", fetch)
    return calls


def lookup(monkeypatch, results, **kwargs):
    transport(monkeypatch, Response({"results": results}, headers={"X-Total-Results": len(results)}))
    return A.annotation_record("NP_000001.2", use_cache=False, **kwargs)


@pytest.mark.parametrize("reverse", [False, True])
def test_exact_version_is_selected_by_identity_not_first_hit(monkeypatch, reverse):
    wrong, correct = entry("P00002", "NP_000001.1"), entry()
    results = [wrong, correct] if reverse else [correct, wrong]
    calls = transport(monkeypatch, Response({"results": results}, headers={
        "X-UniProt-Release": "2026_mock", "X-UniProt-Release-Date": "2026-01-01",
        "X-Total-Results": 2}))
    record = A.annotation_record("NP_000001.2", use_cache=False)
    assert record["state"] == "hit" and A.raw_entry(record) == correct
    assert record["selected"]["matched_cross_references"] == correct["uniProtKBCrossReferences"]
    assert record["selected"]["entry_version"] == 7 and record["selected"]["sequence_version"] == 3
    assert record["request"]["accession"] == "NP_000001.2" and record["request"]["version"] == "2"
    assert record["request"]["namespace"] == "RefSeq" and record["retrieved_utc"].endswith("+00:00")
    assert record["source"]["release"] == "2026_mock" and record["search"]["complete"]
    assert record["responses"][0]["sha256"] == hashlib.sha256(json.dumps({"results": results}).encode()).hexdigest()
    assert record["responses"][0]["final_url"] == calls[0][0]
    query = urllib.parse.parse_qs(urllib.parse.urlsplit(calls[0][0]).query)
    assert query == {"query": ["xref:RefSeq-NP_000001"], "format": ["json"], "size": ["500"]}
    assert calls[0][1] == 30 and "fields" not in query


@pytest.mark.parametrize("reverse", [False, True])
def test_multiple_acceptable_hits_remain_ambiguous(monkeypatch, reverse):
    entries = [entry(), entry("P00002")]
    record = lookup(monkeypatch, list(reversed(entries)) if reverse else entries)
    assert record["state"] == "ambiguous" and A.raw_entry(record) is None
    assert [c["primary_accession"] for c in record["candidates"]] == ["P00001", "P00002"]


@pytest.mark.parametrize("xref,state", [("NP_000001.1", "version_mismatch"),
                                       ("NP_000001", "version_unresolved"),
                                       ("NP_000002.2", "identity_unresolved")])
def test_mismatched_or_missing_protein_version_is_unresolved(monkeypatch, xref, state):
    record = lookup(monkeypatch, [entry(xref=xref)])
    assert record["state"] == state and A.raw_entry(record) is None


def test_refseq_ignores_nucleotide_or_other_database_matches(monkeypatch):
    data = entry(xref="NP_000002.2")
    data["uniProtKBCrossReferences"] += [{"database": "EMBL", "id": "NP_000001.2"}]
    assert lookup(monkeypatch, [data])["state"] == "identity_unresolved"


@pytest.mark.parametrize("has_hit", [False, True])
def test_unsupported_namespace_cannot_certify_nucleotide_or_arbitrary_xrefs(monkeypatch, has_hit):
    data = entry()
    data["uniProtKBCrossReferences"] = [{"database": "EMBL", "id": "AB000001.2"}]
    transport(monkeypatch, Response({"results": [data] if has_hit else []}))
    record = A.annotation_record("AB000001.2", use_cache=False)
    assert record["state"] == "unsupported_namespace" and A.raw_entry(record) is None
    if has_hit:
        assert record["candidates"][0]["state"] == "unsupported_namespace"
        assert record["candidates"][0]["matched_cross_references"] == data["uniProtKBCrossReferences"]


def test_cached_unsupported_namespace_hit_is_not_replayed_as_verified(tmp_path):
    identity = A._request_identity("AB000001.2", None, None)
    previous = A._new_record(identity)
    previous.update(state="hit", entry=entry(), entry_usable=True)
    previous = A._publish(tmp_path / A.CACHE_NAMESPACE, previous)
    assert A.annotation_record("AB000001.2", cache_dir=tmp_path,
                               pinned_record_id=previous["record_id"])["state"] == "pinned_record_unavailable"


def test_unsupported_namespace_legacy_entry_retains_explicit_unverified_compatibility(tmp_path):
    legacy = tmp_path / "AB000001.2.json"
    body = json.dumps(entry()).encode()
    legacy.write_bytes(body)
    record = A.annotation_record("AB000001.2", cache_dir=tmp_path)
    assert record["state"] == "legacy_unverified" and A.raw_entry(record) == entry()
    assert legacy.read_bytes() == body


@pytest.mark.parametrize("xref", ["NP_000001", "NP_000001.2"])
def test_unversioned_lookup_retains_unspecified_version_limitation(monkeypatch, xref):
    transport(monkeypatch, Response({"results": [entry(xref=xref)]}))
    record = A.annotation_record("NP_000001", use_cache=False)
    assert record["state"] == "hit" and record["request"]["version"] is None
    assert "unspecified" in record["rationale"]


@pytest.mark.parametrize("constraint,state", [({"expected_taxon_id": 999}, "organism_mismatch"),
                                             ({"expected_sequence_sha256": "0" * 64}, "sequence_mismatch")])
def test_expected_constraints_reject_conflicts(monkeypatch, constraint, state):
    record = lookup(monkeypatch, [entry()], **constraint)
    assert record["state"] == state and not record["entry_usable"]


@pytest.mark.parametrize("removed,constraint,state", [
    ("organism", {"expected_taxon_id": 123}, "organism_unresolved"),
    ("sequence", {"expected_sequence_sha256": "0" * 64}, "sequence_unresolved")])
def test_missing_identity_evidence_does_not_certify_constraints(monkeypatch, removed, constraint, state):
    data = entry()
    del data[removed]
    assert lookup(monkeypatch, [data], **constraint)["state"] == state


def test_constraints_can_disambiguate_only_on_verified_available_identity(monkeypatch):
    digest = hashlib.sha256(b"MTEST").hexdigest()
    record = lookup(monkeypatch, [entry(), entry("P00002", taxon=999)], expected_taxon_id=123,
                    expected_sequence_sha256=digest)
    assert record["state"] == "hit" and record["selected"]["sequence_sha256"] == digest


@pytest.mark.parametrize("secondary", [False, True])
def test_bare_uniprot_uses_full_json_and_verifies_primary_secondary(monkeypatch, secondary):
    data = entry("P00002" if secondary else "P00001")
    data["secondaryAccessions"] = ["P00001"] if secondary else []
    calls = transport(monkeypatch, Response(data))
    record = A.annotation_record("P00001", use_cache=False)
    assert record["state"] == "hit" and record["entry"] == data
    assert calls[0][0] == "https://rest.uniprot.org/uniprotkb/P00001.json"


def test_bare_uniprot_wrong_identity_is_not_a_hit(monkeypatch):
    transport(monkeypatch, Response(entry("P00002")))
    assert A.annotation_record("P00001", use_cache=False)["state"] == "identity_unresolved"


def test_valid_exhausted_empty_search_alone_is_no_hit(monkeypatch):
    record = lookup(monkeypatch, [])
    assert record["state"] == "no_hit" and record["search"] == {
        "complete": True, "consumed_results": 0, "total_results": 0}


@pytest.mark.parametrize("status", [400, 404, 410, 429, 500, 503])
def test_http_errors_remain_distinct_with_body_hash(monkeypatch, status):
    headers = Message()
    headers["Content-Type"] = "application/json"
    error = urllib.error.HTTPError("https://rest.uniprot.org/uniprotkb/search", status,
                                   "mock failure", headers, io.BytesIO(b'{"error":"synthetic"}'))
    transport(monkeypatch, error)
    record = A.annotation_record("NP_000001.2", use_cache=False)
    assert record["state"] == f"http_{status}" and not record["search"]["complete"]
    assert record["responses"][0]["http_status"] == status
    assert record["responses"][0]["sha256"] == hashlib.sha256(b'{"error":"synthetic"}').hexdigest()
    assert record["error_body"] == '{"error":"synthetic"}' and A.raw_entry(record) is None


@pytest.mark.parametrize("error,state", [(TimeoutError(), "timeout"),
                                        (urllib.error.URLError(TimeoutError()), "timeout"),
                                        (urllib.error.URLError("offline"), "transport_error")])
def test_transport_failure_never_means_no_hit(monkeypatch, error, state):
    transport(monkeypatch, error)
    assert A.annotation_record("NP_000001.2", use_cache=False)["state"] == state


def test_http_failure_without_headers_retains_status_and_response(monkeypatch):
    transport(monkeypatch, urllib.error.HTTPError("https://rest.uniprot.org/uniprotkb/search", 503,
                                                 "synthetic", None, io.BytesIO(b"synthetic error")))
    record = A.annotation_record("NP_000001.2", use_cache=False)
    assert record["state"] == "http_503" and record["responses"][0]["http_status"] == 503
    assert record["responses"][0]["content_type"] is None


@pytest.mark.parametrize("final_url", ["https://evil.invalid/uniprotkb/search", "https://[invalid"])
def test_invalid_final_response_url_remains_recorded_but_unusable(monkeypatch, final_url):
    transport(monkeypatch, Response({"results": [entry()]}, final_url=final_url))
    record = A.annotation_record("NP_000001.2", use_cache=False)
    assert record["state"] == "invalid_response" and A.raw_entry(record) is None
    assert record["responses"][0]["final_url"] == final_url and record["responses"][0]["sha256"]


@pytest.mark.parametrize("payload,state", [(b"broken", "malformed_json"),
                                          (b'{"results":[],"extra":NaN}', "malformed_json"),
                                          (None, "invalid_schema"),
                                          ({}, "invalid_schema"), ({"results": {}}, "invalid_schema"),
                                          ({"results": [{}]}, "invalid_schema"),
                                          ({"results": [{"primaryAccession": "P00001", "sequence": []}]}, "invalid_schema")])
def test_invalid_payloads_never_mean_no_hit(monkeypatch, payload, state):
    transport(monkeypatch, Response(payload))
    assert A.annotation_record("NP_000001.2", use_cache=False)["state"] == state


@pytest.mark.parametrize("accession", ["P00001", "NP_000001.2"])
@pytest.mark.parametrize("number", ["1e9999", "-1e9999"])
@pytest.mark.parametrize("location", ["audit", "nested"])
def test_nested_numeric_overflow_is_unavailable_for_envelope_and_raw_api(monkeypatch, accession, number, location):
    data = entry()
    if location == "audit":
        data["entryAudit"]["entryVersion"] = "OVERFLOW"
    else:
        data["extra"] = {"nested": [1, {"number": "OVERFLOW"}]}
    payload = {"results": [data]} if accession.startswith("NP_") else data
    body = json.dumps(payload).replace('"OVERFLOW"', number).encode()
    calls = transport(monkeypatch, Response(body, headers={"X-UniProt-Release": "synthetic"}), Response(body))
    record = A.annotation_record(accession, use_cache=False)
    assert record["state"] == "malformed_json" and A.raw_entry(record) is None
    assert not record["search"]["complete"] and record["retrieved_utc"]
    assert record["source"]["release"] == "synthetic"
    assert record["responses"][0]["sha256"] == hashlib.sha256(body).hexdigest()
    assert record["responses"][0]["http_status"] == 200
    assert record["responses"][0]["final_url"] == calls[0][0]
    assert L.uniprot_entry(accession, use_cache=False) is None
    assert len(calls) == 2 and not L._CACHE_DIR.exists()


def test_finite_nested_json_numbers_remain_usable(monkeypatch):
    data = entry()
    data["extra"] = {"nested": [1.5, {"number": 1e100}]}
    transport(monkeypatch, Response(data), Response(data))
    record = A.annotation_record("P00001", use_cache=False)
    assert record["state"] == "hit" and A.raw_entry(record) == data
    assert L.uniprot_entry("P00001", use_cache=False) == data


@pytest.mark.parametrize("accession", ["P00001", "NP_000001.2"])
@pytest.mark.parametrize("surrogate", ["\ud800", "\udfff"])
@pytest.mark.parametrize("nested", [False, True])
def test_unpaired_unicode_surrogates_are_unavailable_before_canonical_hashing(monkeypatch, accession, surrogate, nested):
    data = entry()
    data["extra"] = {"nested": [surrogate]} if nested else surrogate
    payload = {"results": [data]} if accession.startswith("NP_") else data
    body = json.dumps(payload).encode()
    calls = transport(monkeypatch, Response(body), Response(body))
    record = A.annotation_record(accession, use_cache=False)
    assert record["state"] == "malformed_json" and A.raw_entry(record) is None
    assert record["responses"][0]["sha256"] == hashlib.sha256(body).hexdigest()
    assert record["responses"][0]["http_status"] == 200
    assert L.uniprot_entry(accession, use_cache=False) is None
    assert len(calls) == 2 and not L._CACHE_DIR.exists()


def test_valid_unicode_and_paired_surrogate_escapes_remain_usable(monkeypatch):
    data = entry()
    data["extra"] = {"nested": ["α", "🧪"]}
    body = json.dumps(data).encode()
    transport(monkeypatch, Response(body), Response(body))
    record = A.annotation_record("P00001", use_cache=False)
    assert record["state"] == "hit" and A.raw_entry(record) == data
    assert L.uniprot_entry("P00001", use_cache=False) == data


def test_returned_pagination_is_consumed_before_unique_selection(monkeypatch):
    next_url = "https://rest.uniprot.org/uniprotkb/search?cursor=returned-opaque&format=json"
    calls = transport(monkeypatch,
                      Response({"results": [entry(xref="NP_000001.1")]}, headers={
                          "X-Total-Results": 2, "Link": f'<{next_url}>; rel="next"'}),
                      Response({"results": [entry("P00002")]}, headers={"X-Total-Results": 2}))
    record = A.annotation_record("NP_000001.2", use_cache=False)
    assert len(calls) == 2 and calls[1][0] == next_url
    assert record["state"] == "hit" and record["selected"]["primary_accession"] == "P00002"
    assert record["search"]["complete"] and len(record["responses"]) == 2


@pytest.mark.parametrize("link", ["https://evil.invalid/uniprotkb/search?cursor=x",
                                  "http://rest.uniprot.org/uniprotkb/search?cursor=x",
                                  "https://rest.uniprot.org/other?cursor=x"])
def test_foreign_or_invalid_next_links_are_not_requested(monkeypatch, link):
    calls = transport(monkeypatch, Response({"results": [entry()]}, headers={"Link": f'<{link}>; rel="next"'}))
    record = A.annotation_record("NP_000001.2", use_cache=False)
    assert record["state"] == "incomplete_search" and len(calls) == 1 and A.raw_entry(record) is None


def test_cycle_page_bound_truncation_and_changed_release_fail_closed(monkeypatch):
    initial = A._SEARCH + "?" + urllib.parse.urlencode({"query": "xref:RefSeq-NP_000001", "format": "json", "size": 500})
    calls = transport(monkeypatch, Response({"results": [entry()]}, headers={"Link": f'<{initial}>; rel="next"'}))
    assert A.annotation_record("NP_000001.2", use_cache=False)["state"] == "incomplete_search"
    assert len(calls) == 1
    monkeypatch.setattr(A, "_MAX_PAGES", 1)
    transport(monkeypatch, Response({"results": []}, headers={"Link": '<https://rest.uniprot.org/uniprotkb/search?cursor=2>; rel="next"'}))
    assert A.annotation_record("NP_000001.2", use_cache=False)["state"] == "incomplete_search"
    transport(monkeypatch, Response({"results": []}, headers={"X-Total-Results": 3}))
    assert A.annotation_record("NP_000001.2", use_cache=False)["state"] == "incomplete_search"
    monkeypatch.setattr(A, "_MAX_PAGES", 100)
    transport(monkeypatch,
              Response({"results": []}, headers={"X-UniProt-Release": "mock1", "Link": '<https://rest.uniprot.org/uniprotkb/search?cursor=2>; rel="next"'}),
              Response({"results": []}, headers={"X-UniProt-Release": "mock2"}))
    assert A.annotation_record("NP_000001.2", use_cache=False)["state"] == "incomplete_search"


@pytest.mark.parametrize("headers", [{"X-Total-Results": "invalid"}, {"Link": "broken"},
                                     {"Link": '<https://rest.uniprot.org/uniprotkb/search?a>; rel="next", <https://rest.uniprot.org/uniprotkb/search?b>; rel="next"'}])
def test_invalid_pagination_metadata_is_incomplete(monkeypatch, headers):
    transport(monkeypatch, Response({"results": []}, headers=headers))
    assert A.annotation_record("NP_000001.2", use_cache=False)["state"] == "incomplete_search"


@pytest.mark.parametrize("body,state", [(json.dumps(entry()).encode(), "legacy_unverified"),
                                       (b"{}", "legacy_empty"), (b"null", "legacy_empty"),
                                       (b"broken", "legacy_invalid"), (b'{"arbitrary":1}', "legacy_invalid")])
def test_historical_cache_bytes_never_change_or_silently_refetch(tmp_path, body, state):
    cache = tmp_path / "legacy"
    cache.mkdir()
    path = cache / "NP_000001.2.json"
    path.write_bytes(body)
    record = A.annotation_record("NP_000001.2", cache_dir=cache)
    assert record["state"] == state and record["retrieved_utc"] is None
    assert record["source"]["release"] is None and record["selected"] is None
    assert path.read_bytes() == body and not (cache / A.CACHE_NAMESPACE).exists()
    assert bool(A.raw_entry(record)) == (state == "legacy_unverified")


def test_constrained_legacy_remains_unverified_without_usable_entry(tmp_path):
    (tmp_path / "NP_000001.2.json").write_text(json.dumps(entry()))
    record = A.annotation_record("NP_000001.2", cache_dir=tmp_path, expected_taxon_id=123)
    assert record["state"] == "legacy_unverified" and record["entry"] == entry()
    assert A.raw_entry(record) is None


def test_cache_only_miss_never_writes(tmp_path):
    cache = tmp_path / "does_not_exist"
    assert A.annotation_record("NP_000001.2", cache_dir=cache, cache_only=True)["state"] == "cache_miss"
    assert not cache.exists()


def test_use_cache_false_bypasses_both_reads_and_writes(monkeypatch, tmp_path):
    legacy = tmp_path / "NP_000001.2.json"
    legacy.write_bytes(b"{}")
    calls = transport(monkeypatch, Response({"results": [entry()]}))
    record = A.annotation_record("NP_000001.2", cache_dir=tmp_path, use_cache=False)
    assert record["state"] == "hit" and len(calls) == 1 and "record_id" not in record
    assert legacy.read_bytes() == b"{}" and sorted(p.name for p in tmp_path.iterdir()) == [legacy.name]


def test_refresh_retains_attempts_and_pin_replays_deterministically(monkeypatch, tmp_path):
    legacy = tmp_path / "NP_000001.2.json"
    legacy.write_bytes(b"{}")
    transport(monkeypatch, Response({"results": [entry()]}), Response({"results": []}))
    first = A.annotation_record("NP_000001.2", cache_dir=tmp_path, refresh=True)
    second = A.annotation_record("NP_000001.2", cache_dir=tmp_path, refresh=True)
    assert first["state"] == "hit" and second["state"] == "no_hit"
    assert first["record_id"] != second["record_id"]
    assert A.annotation_record("NP_000001.2", cache_dir=tmp_path, cache_only=True) == second
    assert A.annotation_record("NP_000001.2", cache_dir=tmp_path, pinned_record_id=first["record_id"]) == first
    assert legacy.read_bytes() == b"{}" and len(list(tmp_path.rglob("records/*.json"))) == 2
    assert A.annotation_record("NP_000001.2", cache_dir=tmp_path, pinned_record_id="0" * 64)["state"] == "pinned_record_unavailable"


def test_constraint_policy_identity_and_pins_cannot_reuse_unconstrained_hit(monkeypatch, tmp_path):
    transport(monkeypatch, Response({"results": [entry()]}))
    first = A.annotation_record("NP_000001.2", cache_dir=tmp_path)
    assert A.annotation_record("NP_000001.2", cache_dir=tmp_path, expected_taxon_id=123,
                               cache_only=True)["state"] == "cache_miss"
    assert A.annotation_record("NP_000001.2", cache_dir=tmp_path, expected_taxon_id=123,
                               pinned_record_id=first["record_id"])["state"] == "pinned_record_unavailable"
    monkeypatch.setattr(A, "MAPPING_POLICY", "synthetic-new-policy")
    assert A.annotation_record("NP_000001.2", cache_dir=tmp_path, cache_only=True)["state"] == "cache_miss"


def test_corrupt_record_is_not_replayed_as_verified(monkeypatch, tmp_path):
    transport(monkeypatch, Response({"results": [entry()]}))
    record = A.annotation_record("NP_000001.2", cache_dir=tmp_path)
    path = next(tmp_path.rglob(f"{record['record_id']}.json"))
    path.write_text('{"schema_version":1,"entry_usable":true}')
    assert A.annotation_record("NP_000001.2", cache_dir=tmp_path, cache_only=True)["state"] == "cache_miss"
    assert A.annotation_record("NP_000001.2", cache_dir=tmp_path,
                               pinned_record_id=record["record_id"])["state"] == "pinned_record_unavailable"


def test_interrupted_publication_preserves_prior_latest_and_readable_record(monkeypatch, tmp_path):
    transport(monkeypatch, Response({"results": [entry()]}), Response({"results": []}))
    first = A.annotation_record("NP_000001.2", cache_dir=tmp_path)
    index = next(tmp_path.rglob("latest.json"))
    prior_bytes = index.read_bytes()
    original_replace = A.os.replace

    def fail_latest(source, destination):
        if destination.name == "latest.json":
            raise PermissionError("synthetic interruption")
        original_replace(source, destination)

    monkeypatch.setattr(A.os, "replace", fail_latest)
    monkeypatch.setattr(A.time, "sleep", lambda *a: None)
    failed = A.annotation_record("NP_000001.2", cache_dir=tmp_path, refresh=True)
    assert failed["state"] == "no_hit" and failed["cache_publication"]["state"] == "failed"
    assert index.read_bytes() == prior_bytes
    assert A.annotation_record("NP_000001.2", cache_dir=tmp_path, cache_only=True) == first
    assert not list(tmp_path.rglob("*.tmp"))
    assert len(list(tmp_path.rglob("records/*.json"))) == 2


def test_atomic_write_retries_bounded_permission_errors_and_cleans_temp(monkeypatch, tmp_path):
    original_replace = A.os.replace
    calls = []

    def replace(source, destination):
        calls.append((source, destination))
        assert source.parent == destination.parent and source.read_bytes()
        if len(calls) < 3:
            raise PermissionError("synthetic transient lock")
        original_replace(source, destination)

    monkeypatch.setattr(A.os, "replace", replace)
    monkeypatch.setattr(A.time, "sleep", lambda *a: None)
    destination = tmp_path / "record.json"
    A._atomic_json(destination, {"complete": True})
    assert len(calls) == 3 and json.loads(destination.read_text()) == {"complete": True}
    assert not list(tmp_path.glob("*.tmp"))


def test_wrapper_preserves_public_signature_and_returns_raw_entry_only(monkeypatch):
    assert str(inspect.signature(L.uniprot_entry)) == "(accession: 'str', *, use_cache: 'bool' = True, timeout: 'int' = 30) -> 'dict | None'"
    original = entry()
    calls = transport(monkeypatch, Response({"results": [original]}), Response({"results": []}))
    assert L.uniprot_entry("NP_000001.2", use_cache=False) == original
    assert L.uniprot_entry("NP_000001.2", use_cache=False) is None and len(calls) == 2
    assert not L._CACHE_DIR.exists()


def test_snapshot_records_default_provenance_once_and_keeps_availability_distinct(monkeypatch, tmp_path):
    from test_exp069_consensus import _stub_inference
    _stub_inference(monkeypatch)
    transport(monkeypatch, Response({"results": [entry()]}))
    recorded = []
    acquire = A.annotation_record

    def acquire_once(*args, **kwargs):
        recorded.append(args[0])
        assert len(recorded) == 1
        return acquire(*args, **kwargs)

    monkeypatch.setattr(A, "annotation_record", acquire_once)

    def derive(ctx, tf, *, entry_fetcher):
        raw = entry_fetcher("NP_000001.2")
        assert raw == entry() and "schema_version" not in raw
        return []

    monkeypatch.setattr(OC, "effectors_from_context", derive)
    cons = IND.infer_inducer("SYNTHETIC", "test", genome_ctx=NS(), tf_gene=NS(),
                            allow_ncbi=True, metalnet=False)
    snapshot = json.loads(json.dumps(asdict(cons)))["chemistry_snapshot"]
    annotation = snapshot["annotations"][0]
    assert annotation["availability"] == "hit" and annotation["annotation"]["record_id"]
    assert annotation["chemistry_status"] == "reactions_available" and snapshot["availability"] == "available"
    assert str(tmp_path) not in json.dumps(snapshot)
    from predictor.report import ligand_regulon
    ligand_regulon.accepted_ligands(asdict(cons))
    assert recorded == ["NP_000001.2"]


@pytest.mark.parametrize("data,expected", [(None, "unavailable"), ({}, "unavailable"),
                                           ({"primaryAccession": "P00001"}, "no_catalytic_annotation"),
                                           ({"primaryAccession": "P00001", "comments": [{"commentType": "CATALYTIC ACTIVITY"}]},
                                            "no_usable_reaction_participants")])
def test_snapshot_missing_annotation_does_not_claim_nonenzyme(monkeypatch, data, expected):
    from test_exp069_consensus import _stub_inference
    _stub_inference(monkeypatch)

    def derive(ctx, tf, *, entry_fetcher):
        entry_fetcher("SYNTHETIC")
        return []

    monkeypatch.setattr(OC, "effectors_from_context", derive)
    cons = IND.infer_inducer("SYNTHETIC", "test", genome_ctx=NS(), tf_gene=NS(), allow_ncbi=True,
                            metalnet=False, entry_fetcher=lambda _: deepcopy(data))
    assert cons.chemistry_snapshot["annotations"][0]["chemistry_status"] == expected
    assert cons.chemistry_snapshot["availability"] == ("unavailable" if not data else "available")


@pytest.mark.parametrize("invalid", [False, True])
def test_snapshot_callback_failures_keep_explicit_unavailability_without_local_paths(monkeypatch, invalid, tmp_path):
    from test_exp069_consensus import _stub_inference
    _stub_inference(monkeypatch)

    def fetch(_):
        if invalid:
            return "not a raw dictionary"
        raise OSError(str(tmp_path / "private-cache"))

    def derive(ctx, tf, *, entry_fetcher):
        entry_fetcher("SYNTHETIC")
        return []

    monkeypatch.setattr(OC, "effectors_from_context", derive)
    cons = IND.infer_inducer("SYNTHETIC", "test", genome_ctx=NS(), tf_gene=NS(), allow_ncbi=True,
                            metalnet=False, entry_fetcher=fetch)
    snapshot = cons.chemistry_snapshot
    assert snapshot["status"] == "unavailable" and snapshot["availability"] == "unavailable"
    assert snapshot["annotations"][0]["chemistry_status"] == "unavailable"
    assert str(tmp_path) not in json.dumps(snapshot)


def test_snapshot_partial_annotation_keeps_available_record_without_claiming_absent_chemistry(monkeypatch):
    from test_exp069_consensus import _stub_inference
    _stub_inference(monkeypatch)
    calls = []

    def fetch(pid):
        calls.append(pid)
        return entry() if pid == "FIRST" else None

    def derive(ctx, tf, *, entry_fetcher):
        entry_fetcher("FIRST")
        entry_fetcher("SECOND")
        return []

    monkeypatch.setattr(OC, "effectors_from_context", derive)
    cons = IND.infer_inducer("SYNTHETIC", "test", genome_ctx=NS(), tf_gene=NS(), allow_ncbi=True,
                            metalnet=False, entry_fetcher=fetch)
    assert calls == ["FIRST", "SECOND"]
    snapshot = cons.chemistry_snapshot
    assert snapshot["status"] == "recorded" and snapshot["availability"] == "partial"
    assert [r["chemistry_status"] for r in snapshot["annotations"]] == ["reactions_available", "unavailable"]
