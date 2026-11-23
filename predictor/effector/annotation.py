"""Versioned UniProt acquisition records; chemistry receives only ``raw_entry(record)``.

NEW API: ``annotation_record`` adds explicit cache-only/refresh and deterministic pin
replay. Its cache is independent of historical ``<accession>.json`` files. No cache
age expires automatically. Expected sequence hashes must describe the exact protein;
the version-stripped protein_seqs cache does not certify them. All identity choices
are recorded hypotheses about database mapping, not experimental reaction evidence.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import re
import socket
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from predictor import resources

SCHEMA_VERSION = 1
MAPPING_POLICY = "exact-protein-xref-v1"
CACHE_NAMESPACE = "annotation-v1"
_HOST = "rest.uniprot.org"
_SEARCH = f"https://{_HOST}/uniprotkb/search"
_MAX_PAGES = 100
_ACC = re.compile(r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})\Z")
_HEX = re.compile(r"[0-9a-f]{64}\Z")


def _bytes(value) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
                      allow_nan=False).encode("utf-8")


def _sha(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _valid_entry(entry) -> bool:
    if not isinstance(entry, dict) or not isinstance(entry.get("primaryAccession"), str) or not entry["primaryAccession"]:
        return False
    for name in ("secondaryAccessions", "comments", "uniProtKBCrossReferences"):
        if name in entry and not isinstance(entry[name], list):
            return False
    if any(not isinstance(acc, str) for acc in entry.get("secondaryAccessions", [])):
        return False
    if any(not isinstance(x, dict) or not isinstance(x.get("database"), str)
           or not isinstance(x.get("id"), str) for x in entry.get("uniProtKBCrossReferences", [])):
        return False
    if any(not isinstance(x, dict) for x in entry.get("comments", [])):
        return False
    for name in ("organism", "sequence", "entryAudit"):
        if name in entry and not isinstance(entry[name], dict):
            return False
    sequence = entry.get("sequence", {}).get("value")
    taxon = entry.get("organism", {}).get("taxonId")
    if sequence is not None and (not isinstance(sequence, str) or not sequence):
        return False
    if taxon is not None and (not isinstance(taxon, int) or isinstance(taxon, bool)):
        return False
    return True


def _request_identity(accession, expected_taxon_id, expected_sequence_sha256):
    accession = str(accession).strip()
    if not accession or not re.fullmatch(r"[A-Za-z0-9_.-]+", accession):
        raise ValueError("accession must be a nonempty accession token")
    base, dot, version = accession.rpartition(".")
    if not dot or not version.isdecimal():
        base, version = accession, None
    namespace = "UniProtKB" if _ACC.fullmatch(accession) else (
        "RefSeq" if re.match(r"(?:NP|WP|XP|YP|AP|ZP)_", base) else "cross_reference")
    if expected_taxon_id is not None and (isinstance(expected_taxon_id, bool)
                                         or not isinstance(expected_taxon_id, int) or expected_taxon_id <= 0):
        raise ValueError("expected_taxon_id must be a positive integer")
    if expected_sequence_sha256 is not None:
        expected_sequence_sha256 = str(expected_sequence_sha256).lower()
        if not _HEX.fullmatch(expected_sequence_sha256):
            raise ValueError("expected_sequence_sha256 must be a SHA256 hex digest")
    return {"accession": accession, "base_accession": base, "version": version,
            "namespace": namespace, "expected_taxon_id": expected_taxon_id,
            "expected_sequence_sha256": expected_sequence_sha256,
            "mapping_policy": MAPPING_POLICY}


def _new_record(identity):
    return {"schema_version": SCHEMA_VERSION, "request": identity, "state": "unavailable",
            "rationale": "", "entry": None, "entry_usable": False, "selected": None,
            "candidates": [], "responses": [], "retrieved_utc": None,
            "search": {"complete": False, "consumed_results": 0, "total_results": None},
            "source": {"name": "UniProtKB", "release": None, "release_date": None}}


def raw_entry(record: dict) -> dict | None:
    """Return the untouched selected raw entry only when this record permits its use."""
    return record.get("entry") if record.get("entry_usable") else None


def _atomic_json(destination: Path, value) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(prefix=f".{destination.name}.", suffix=".tmp",
                                     dir=destination.parent, delete=False) as stream:
        temporary = Path(stream.name)
        try:
            stream.write(_bytes(value))
            stream.flush()
            os.fsync(stream.fileno())
        except BaseException:
            stream.close()
            temporary.unlink(missing_ok=True)
            raise
    try:
        for attempt in range(5):
            try:
                os.replace(temporary, destination)
                break
            except PermissionError:
                if attempt == 4:
                    raise
                time.sleep(0.05 * (attempt + 1))
    finally:
        temporary.unlink(missing_ok=True)


def _publish(root: Path, record: dict) -> dict:
    record = dict(record)
    record_id = _sha(_bytes(record))
    record["record_id"] = record_id
    records = root / "records"
    path = records / f"{record_id}.json"
    if not path.exists():
        _atomic_json(path, record)
    request_key = _sha(_bytes(record["request"]))
    _atomic_json(root / "requests" / request_key / "latest.json",
                 {"schema_version": SCHEMA_VERSION, "record_id": record_id})
    return record


def _read_record(root: Path, identity, pin=None):
    try:
        if pin is None:
            request_key = _sha(_bytes(identity))
            index = json.loads((root / "requests" / request_key / "latest.json").read_text(encoding="utf-8"))
            if not isinstance(index, dict) or index.get("schema_version") != SCHEMA_VERSION:
                return None
            pin = index.get("record_id")
        if not isinstance(pin, str) or not _HEX.fullmatch(pin):
            return None
        record = json.loads((root / "records" / f"{pin}.json").read_text(encoding="utf-8"))
        if not isinstance(record, dict):
            return None
        payload = {k: v for k, v in record.items() if k != "record_id"}
        if (record.get("record_id") != pin or _sha(_bytes(payload)) != pin
                or record.get("schema_version") != SCHEMA_VERSION or record.get("request") != identity):
            return None
        if record.get("entry_usable"):
            entry = record.get("entry")
            if not _valid_entry(entry) or _candidate(entry, identity)["state"] != "acceptable":
                return None
        return record
    except (OSError, ValueError, TypeError):
        return None


def _legacy(cache_dir: Path, identity):
    path = cache_dir / f"{identity['accession']}.json"
    if not path.exists():
        return None
    record = _new_record(identity)
    record["source"] = {"name": "legacy_cache", "release": None, "release_date": None}
    try:
        data = path.read_bytes()
        record["legacy_sha256"] = _sha(data)
        entry = json.loads(data.decode("utf-8"))
        if entry is None or entry == {}:
            record.update(state="legacy_empty", rationale="historical empty value; availability is unverified")
        elif _valid_entry(entry):
            constrained = identity["expected_taxon_id"] is not None or identity["expected_sequence_sha256"] is not None
            record.update(state="legacy_unverified", entry=entry, entry_usable=not constrained,
                          rationale="historical raw entry; mapping, versions and retrieval time unverified"
                          + ("; expected identity constraints cannot be certified" if constrained else ""))
        else:
            record.update(state="legacy_invalid", rationale="historical value is not a valid raw entry")
    except (OSError, UnicodeError, ValueError):
        record.update(state="legacy_invalid", rationale="historical cache cannot be parsed")
    return record


def _valid_url(url: str, *, search: bool) -> bool:
    try:
        parsed = urllib.parse.urlsplit(url)
    except ValueError:
        return False
    return (parsed.scheme == "https" and parsed.netloc == _HOST and not parsed.fragment
            and (parsed.path == "/uniprotkb/search" if search else
                 bool(re.fullmatch(r"/uniprotkb/[A-Z0-9]+\.json", parsed.path))))


def _next_link(header: str):
    """Read returned next URLs without constructing cursors or trusting foreign links."""
    if not header:
        return None
    next_urls = []
    for part in re.split(r",(?=\s*<)", header):
        match = re.fullmatch(r'\s*<([^>]+)>\s*((?:;[^;]*)*)', part)
        if not match:
            raise ValueError("malformed pagination Link")
        rel = re.search(r';\s*rel\s*=\s*(?:"([^"]+)"|([^;\s]+))', match[2], re.I)
        if rel and "next" in (rel[1] or rel[2]).split():
            next_urls.append(match[1])
    if len(next_urls) > 1:
        raise ValueError("multiple pagination next links")
    if next_urls and not _valid_url(next_urls[0], search=True):
        raise ValueError("pagination next link has an unexpected HTTPS host/path")
    return next_urls[0] if next_urls else None


def _response(record, url, status, headers, body, *, final_url=None):
    headers = headers or {}
    meta = {"url": url, "final_url": final_url, "http_status": status, "sha256": _sha(body),
            "content_type": headers.get("Content-Type"),
            "release": headers.get("X-UniProt-Release"),
            "release_date": headers.get("X-UniProt-Release-Date"),
            "total_results": headers.get("X-Total-Results"), "link": headers.get("Link")}
    record["responses"].append(meta)
    for field, header in (("release", "X-UniProt-Release"), ("release_date", "X-UniProt-Release-Date")):
        value = headers.get(header)
        previous = record["source"][field]
        if previous and value and previous != value:
            return False
        if value:
            record["source"][field] = value
    return True


def _candidate(entry, identity):
    acc = identity["accession"]
    refs = entry.get("uniProtKBCrossReferences", [])
    matched = []
    state, rationale = "identity_unresolved", "requested exact accession is absent"
    if identity["namespace"] == "UniProtKB":
        if acc == entry["primaryAccession"] or acc in entry.get("secondaryAccessions", []):
            state, rationale = "acceptable", "requested primary/secondary UniProt accession verified"
    elif identity["namespace"] == "RefSeq":
        refs = [r for r in refs if r["database"] == "RefSeq"]
        matched = [r for r in refs if r["id"] == acc]
        related = [r for r in refs if r["id"].split(".")[0] == identity["base_accession"]]
        if matched:
            state, rationale = "acceptable", "exact protein cross-reference verified"
            if identity["version"] is None:
                rationale += "; requested version unspecified"
        elif related and identity["version"]:
            state = "version_mismatch" if any("." in r["id"] for r in related) else "version_unresolved"
            rationale = "protein cross-reference version differs or is unavailable"
        elif related and not identity["version"]:
            matched = related
            state, rationale = "acceptable", "unversioned requested accession matches protein xref base; version unspecified"
    else:
        matched = [r for r in refs if r["id"] == acc]
        state, rationale = "unsupported_namespace", "arbitrary cross-reference membership cannot certify protein identity"
    seq = entry.get("sequence", {}).get("value")
    sequence_sha256 = _sha(seq.encode("utf-8")) if seq else None
    taxon = entry.get("organism", {}).get("taxonId")
    if state == "acceptable":
        expected_taxon, expected_seq = identity["expected_taxon_id"], identity["expected_sequence_sha256"]
        if expected_taxon is not None and taxon is None:
            state, rationale = "organism_unresolved", "expected taxon cannot be checked"
        elif expected_taxon is not None and taxon != expected_taxon:
            state, rationale = "organism_mismatch", "returned taxon differs from expectation"
        elif expected_seq is not None and sequence_sha256 is None:
            state, rationale = "sequence_unresolved", "expected sequence cannot be checked"
        elif expected_seq is not None and sequence_sha256 != expected_seq:
            state, rationale = "sequence_mismatch", "returned sequence differs from expectation"
    return {"primary_accession": entry["primaryAccession"], "matched_cross_references": matched,
            "state": state, "rationale": rationale, "taxon_id": taxon,
            "sequence_sha256": sequence_sha256, "sequence_length": len(seq) if seq else None,
            "entry_version": entry.get("entryAudit", {}).get("entryVersion"),
            "sequence_version": entry.get("entryAudit", {}).get("sequenceVersion"),
            "entry_sha256": _sha(_bytes(entry))}


def _acquire(record, timeout):
    identity = record["request"]
    search = identity["namespace"] != "UniProtKB"
    query = (f"xref:RefSeq-{identity['base_accession']}" if identity["namespace"] == "RefSeq"
             else f"xref:{identity['accession']}")
    url = (_SEARCH + "?" + urllib.parse.urlencode({"query": query, "format": "json", "size": 500})
           if search else f"https://{_HOST}/uniprotkb/{identity['accession']}.json")
    record["query"] = query if search else None
    record["retrieved_utc"] = datetime.now(timezone.utc).isoformat()
    entries, seen = [], set()
    while url:
        if url in seen or len(seen) >= _MAX_PAGES:
            record.update(state="incomplete_search", rationale="pagination cycle or page bound reached")
            return record
        seen.add(url)
        try:
            with urllib.request.urlopen(url, timeout=timeout) as response:
                body = response.read()
                headers, status = response.headers, response.status
                final_url = response.geturl()
            stable_source = _response(record, url, status, headers, body, final_url=final_url)
            if not _valid_url(final_url, search=search):
                record.update(state="invalid_response", rationale="response redirected outside expected HTTPS host/path")
                return record
            if not stable_source:
                record.update(state="incomplete_search", rationale="source release changed during pagination")
                return record
        except urllib.error.HTTPError as exc:
            body = exc.read()
            _response(record, url, exc.code, exc.headers, body, final_url=exc.geturl())
            record.update(state=f"http_{exc.code}", rationale="HTTP failure; no biological absence inferred",
                          error_body=body.decode("utf-8", errors="replace"))
            return record
        except (TimeoutError, socket.timeout):
            record.update(state="timeout", rationale="request timed out")
            return record
        except (urllib.error.URLError, OSError) as exc:
            state = "timeout" if isinstance(getattr(exc, "reason", None), TimeoutError) else "transport_error"
            record.update(state=state, rationale="request transport failed", error_type=type(exc).__name__)
            return record
        if status != 200:
            record.update(state=f"http_{status}", rationale="non-success response; no biological absence inferred",
                          error_body=body.decode("utf-8", errors="replace"))
            return record
        try:
            data = json.loads(body.decode("utf-8"), parse_constant=_invalid_json_constant,
                              parse_float=_finite_json_float)
            # Validate every nested value before selection or content-addressed hashing.
            # JSON escapes can decode to lone surrogates that cannot represent UTF-8.
            _bytes(data)
        except (ValueError, UnicodeError, RecursionError):
            record.update(state="malformed_json", rationale="response is not JSON")
            return record
        page = data.get("results") if search and isinstance(data, dict) else ([data] if not search else None)
        if not isinstance(page, list) or any(not _valid_entry(entry) for entry in page):
            record.update(state="invalid_schema", rationale="response lacks complete valid entry/results structure")
            return record
        entries.extend(page)
        record["search"]["consumed_results"] = len(entries)
        try:
            total = headers.get("X-Total-Results")
            if total is not None:
                if not str(total).isdecimal():
                    raise ValueError("invalid total")
                total = int(total)
                previous = record["search"]["total_results"]
                if previous is not None and previous != total:
                    raise ValueError("changed total")
                record["search"]["total_results"] = total
            url = _next_link(headers.get("Link", "")) if search else None
        except ValueError as exc:
            record.update(state="incomplete_search", rationale=str(exc))
            return record
    total = record["search"]["total_results"]
    if search and total is not None and len(entries) != total:
        record.update(state="incomplete_search", rationale="consumed results differ from declared total")
        return record
    record["search"]["complete"] = True
    unique = {_sha(_bytes(entry)): entry for entry in entries}
    evaluated = [(_candidate(entry, identity), entry) for entry in unique.values()]
    evaluated.sort(key=lambda pair: (pair[0]["primary_accession"], pair[0]["entry_sha256"]))
    record["candidates"] = [candidate for candidate, _ in evaluated]
    accepted = [(candidate, entry) for candidate, entry in evaluated if candidate["state"] == "acceptable"]
    if len(accepted) == 1:
        selected, entry = accepted[0]
        record.update(state="hit", rationale=selected["rationale"], selected=selected, entry=entry, entry_usable=True)
    elif accepted:
        record.update(state="ambiguous", rationale="multiple acceptable distinct UniProt entries")
    elif identity["namespace"] == "cross_reference":
        record.update(state="unsupported_namespace", rationale="requested namespace cannot certify an exact protein identity")
    elif not evaluated and search:
        record.update(state="no_hit", rationale="valid complete exhausted empty search")
    else:
        states = {candidate["state"] for candidate, _ in evaluated}
        record.update(state=next(iter(states)) if len(states) == 1 else "identity_unresolved",
                      rationale="no acceptable exact identity; see candidate decisions")
    return record


def _invalid_json_constant(value):
    raise ValueError(f"non-JSON constant: {value}")


def _finite_json_float(value):
    number = float(value)
    if not math.isfinite(number):
        raise ValueError("JSON number exceeds finite float range")
    return number


def annotation_record(accession: str, *, use_cache: bool = True, timeout: int = 30,
                      cache_only: bool = False, refresh: bool = False,
                      expected_taxon_id: int | None = None, expected_sequence_sha256: str | None = None,
                      pinned_record_id: str | None = None, cache_dir: Path | None = None) -> dict:
    """NEW API: acquire/replay a provenance envelope, never pass it to an enzyme parser.

    Normal cache-first calls prefer new records then legacy; legacy is never silently
    refreshed. ``refresh=True`` explicitly acquires a new immutable attempt. Pins
    never fall back or request. ``cache_only`` never writes. ``use_cache=False``
    bypasses reads AND writes; combining it with a pin is invalid. Constrained legacy
    records retain raw data but cannot supply a usable entry. No freshness clock runs.
    """
    if cache_only and refresh:
        raise ValueError("cache_only and refresh are mutually exclusive")
    if pinned_record_id is not None and (not use_cache or refresh):
        raise ValueError("pins require cache reads and cannot be refreshed")
    identity = _request_identity(accession, expected_taxon_id, expected_sequence_sha256)
    cache_dir = Path(cache_dir) if cache_dir is not None else Path(
        os.environ.get("PREDICTOR_EFFECTOR_CACHE", resources.cache_path("effector")))
    root = cache_dir / CACHE_NAMESPACE
    if use_cache and not refresh:
        cached = _read_record(root, identity, pinned_record_id)
        if cached is not None:
            return cached
        if pinned_record_id is not None:
            record = _new_record(identity)
            record.update(state="pinned_record_unavailable", rationale="requested pin is absent, invalid or incompatible")
            return record
        legacy = _legacy(cache_dir, identity)
        if legacy is not None:
            return legacy
    record = _new_record(identity)
    if cache_only:
        record.update(state="cache_miss", rationale="cache-only lookup has no compatible record")
        return record
    record = _acquire(record, timeout)
    if use_cache:
        try:
            return _publish(root, record)
        except OSError as exc:
            # An acquisition can remain usable even if local publication fails.
            record["cache_publication"] = {"state": "failed", "error_type": type(exc).__name__}
    return record
