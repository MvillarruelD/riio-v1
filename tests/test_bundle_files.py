"""File publication and presentation contracts using only synthetic records."""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

from predictor.report import bundle_contract as contract
from predictor.report import outputs, render_report


@pytest.fixture
def bundle(tmp_path):
    root = tmp_path / "bundle"
    root.mkdir()
    (root / "result.txt").write_text("synthetic", encoding="utf-8")
    (root / "REPORT.html").write_text('<a href="result.txt">Data</a>', encoding="utf-8")
    manifest = {"schema_version": 2, "status": "complete", "report_generated": True,
                "required_artifacts": ["result.txt", "REPORT.html"], "artifacts": contract.inventory(root)}
    (root / "bundle_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    return root


def test_bundle_is_portable_and_checks_all_artifacts(bundle, tmp_path):
    relocated = tmp_path / "relocated"
    shutil.copytree(bundle, relocated)
    assert contract.validate_files(relocated) == []
    (relocated / "result.txt").write_text("altered", encoding="utf-8")
    assert any("checksum mismatch" in issue for issue in contract.validate_files(relocated))


@pytest.mark.parametrize(
    "path",
    [
        "",
        "../outside.txt",
        "a/../outside.txt",
        "/absolute",
        "C:/file",
        "a\\b",
        "./file",
        "a//b",
        "a/",
        "name.",
        "name ",
        "CON",
        "con.txt",
        "dir/PRN.dat",
        "CONIN$",
        "conin$.txt",
        "CONOUT$.log",
        "COM¹",
        "com².txt",
        "LPT³.dat",
        "CON .txt",
        "dir/AUX   .json",
    ],
)
def test_bundle_paths_are_contained(tmp_path, path):
    with pytest.raises(ValueError):
        contract.contained_file(tmp_path, path)


def test_python_310_fallback_rejects_extended_windows_devices(tmp_path, monkeypatch):
    monkeypatch.setattr(contract, "_OS_PATH_ISRESERVED", None)
    for path in ["CONIN$.txt", "CONOUT$", "COM¹.log", "LPT³", "NUL .txt"]:
        with pytest.raises(ValueError, match="reserved Windows basename"):
            contract.contained_file(tmp_path, path)


@pytest.mark.parametrize("ref", ["missing.png", "../escape.png", "%2e%2e/escape.png", "file:///tmp/x", "https://example.com/x.png"])
def test_report_cannot_depend_on_external_or_missing_images(bundle, ref):
    (bundle / "REPORT.html").write_text(f'<img src="{ref}">', encoding="utf-8")
    assert contract.report_asset_issues(bundle)


def test_external_citations_and_fragment_links_are_allowed(bundle):
    (bundle / "REPORT.html").write_text(
        '<a href="https://example.com/paper">HTTPS</a>'
        '<a href="http://example.com/paper">HTTP</a>'
        '<a href="mailto:author@example.com">Email</a>'
        '<a href="#files">Files</a>'
        '<img src="data:image/png;base64,AA==">',
        encoding="utf-8",
    )
    assert contract.report_asset_issues(bundle) == []


def _replace_report_with_valid_checksum(bundle, html):
    (bundle / "REPORT.html").write_text(html, encoding="utf-8")
    manifest_path = bundle / "bundle_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["artifacts"] = contract.inventory(bundle)
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")


def _set_manifest_version(bundle, version):
    manifest_path = bundle / "bundle_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if version == 1:
        manifest.pop("schema_version", None)
        manifest.pop("artifacts", None)
    else:
        manifest["schema_version"] = 2
        manifest["artifacts"] = contract.inventory(bundle)
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")


@pytest.mark.parametrize(
    "html",
    [
        "<style>body { background: url(https://example.invalid/x.png) }</style>",
        "<style>@import 'https://example.invalid/report.css';</style>",
        r"<style>body { background: u\72 l(https://example.invalid/x.png) }</style>",
        r"<style>@\69 mport 'https://example.invalid/report.css';</style>",
        "<style>body { background: u/**/rl(https://example.invalid/x.png) }</style>",
        '<div style="background: url(missing.png)"></div>',
        '<img srcset="https://example.invalid/x.png 2x">',
        '<object data="missing.pdf"></object>',
        '<base href="https://example.invalid/">',
        '<script src="local.js"></script>',
        '<div onclick="location.href=\'https://example.invalid/\'">go</div>',
        '<iframe srcdoc="<p>active</p>"></iframe>',
        '<style>body { background: image-set("remote.png" 1x) }</style>',
        '<style>body { background: -webkit-image-set("remote.png" 1x) }</style>',
        '<style>@font-face { src: src("remote.woff2") }</style>',
        '<style>body { width: expression(alert(1)) }</style>',
        '<div style="background:url(remote.png)" style="color:red"></div>',
        '<meta http-equiv="refresh" http-equiv="content-type" content="0;url=x">',
        '<svg><path clip-path="url(https://example.invalid/c.svg#x)"></path></svg>',
    ],
)
def test_unsupported_html_resource_surfaces_cannot_receive_v2_assurance(bundle, html):
    _replace_report_with_valid_checksum(bundle, html)
    audit = contract.audit_files(bundle)
    assert audit["valid"] is False
    assert audit["manifest_version"] == 2
    assert audit["assurance"] == "unavailable"
    assert audit["html_checks_performed"] is True
    assert audit["issues"]


@pytest.mark.parametrize(
    ("html", "expected_issue"),
    [
        (
            '<link rel="preload" as="image" href="https://example.invalid/x.png">',
            "unsupported element <link>",
        ),
        (
            '<link rel="preload" as="image" '
            'imagesrcset="https://example.invalid/x.png 2x">',
            "unsupported resource attribute link[imagesrcset]",
        ),
        (
            '<img src="data:image/png;base64,AA==" '
            'lowsrc="https://example.invalid/x.png">',
            "unsupported resource attribute img[lowsrc]",
        ),
        (
            '<svg><rect fill="url(https://example.invalid/a.svg#x)"></rect></svg>',
            "attribute rect[fill] contains CSS resource syntax",
        ),
        (
            '<svg><g filter="url(https://example.invalid/f.svg#x)"></g></svg>',
            "attribute g[filter] contains CSS resource syntax",
        ),
        (
            '<svg><path cursor="url(https://example.invalid/c.cur), auto"></path></svg>',
            "attribute path[cursor] contains CSS resource syntax",
        ),
    ],
)
def test_known_browser_resource_surfaces_fail_closed_with_valid_v2_metadata(
    bundle, html, expected_issue
):
    """Known browser-loading surfaces cannot retain checksum assurance."""
    _replace_report_with_valid_checksum(bundle, html)

    audit = contract.audit_files(bundle)

    assert audit["valid"] is False
    assert audit["manifest_version"] == 2
    assert audit["assurance"] == "unavailable"
    assert audit["html_checks_performed"] is True
    assert any(expected_issue in issue for issue in audit["issues"])


@pytest.mark.parametrize("version", [1, 2])
@pytest.mark.parametrize(
    "html",
    [
        "<!--><script src='https://example.invalid/x.js'></script>-->",
        "<!--><img src='https://example.invalid/x.png'>-->",
        "<![CDATA[><script src='https://example.invalid/x.js'></script>]]>",
    ],
)
def test_html_tokenization_breakouts_fail_closed_for_every_manifest_version(
    bundle, version, html
):
    _replace_report_with_valid_checksum(bundle, html)
    _set_manifest_version(bundle, version)

    audit = contract.audit_files(bundle)

    assert audit["valid"] is False
    assert audit["manifest_version"] == version
    assert audit["assurance"] == "unavailable"
    assert audit["html_checks_performed"] is True
    assert any("comment, CDATA, or declaration" in issue for issue in audit["issues"])


@pytest.mark.parametrize("version", [1, 2])
def test_simple_generated_doctype_remains_supported(bundle, version):
    _replace_report_with_valid_checksum(
        bundle,
        '<!doctype html><html><body><a href="result.txt">Data</a></body></html>',
    )
    _set_manifest_version(bundle, version)

    audit = contract.audit_files(bundle)

    assert audit["valid"] is True
    assert audit["manifest_version"] == version


def test_documented_harmless_attributes_remain_supported(bundle):
    _replace_report_with_valid_checksum(
        bundle,
        '<div class="panel" id="summary" title="Synthetic" lang="en" dir="ltr" '
        'role="region" tabindex="0" hidden aria-label="Synthetic panel" '
        'data-state="empty" style="color: red">Content</div>'
        '<a href="https://example.com/paper">Citation</a>'
        '<img src="data:image/png;base64,AA==" alt="Synthetic" loading="lazy">',
    )
    audit = contract.audit_files(bundle)
    assert audit["valid"] is True
    assert audit["assurance"] == "v2-complete-checksums"


def test_unicode_and_percent_encoded_local_assets_use_the_same_path_rules(bundle):
    asset = bundle / "figures" / "café β.png"
    asset.parent.mkdir()
    asset.write_bytes(b"synthetic image")
    (bundle / "REPORT.html").write_text(
        '<img src="figures/caf%C3%A9%20%CE%B2.png">', encoding="utf-8"
    )
    assert contract.report_asset_issues(bundle) == []


@pytest.mark.parametrize(
    "ref",
    [
        "%2e%2e/escape.png",
        "%2Fabsolute.png",
        "figures%2fimage.png",
        "figures%5Cimage.png",
        "bad%escape.png",
        "CON%2Epng",
        "name%2E",
    ],
)
def test_percent_encoded_references_cannot_bypass_portable_path_rules(bundle, ref):
    (bundle / "REPORT.html").write_text(f'<img src="{ref}">', encoding="utf-8")
    assert contract.report_asset_issues(bundle)


def test_manifest_casefold_collisions_fail_before_artifact_reads(bundle, monkeypatch):
    manifest = json.loads((bundle / "bundle_manifest.json").read_text(encoding="utf-8"))
    manifest["required_artifacts"] = ["REPORT.html", "Result.txt", "result.txt"]
    (bundle / "bundle_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")

    def fail_if_hashed(*args, **kwargs):
        raise AssertionError("malformed manifest must fail before artifact reads")

    monkeypatch.setattr(contract, "sha256", fail_if_hashed)
    issues = contract.validate_files(bundle)
    assert any("case-fold collision" in issue for issue in issues)


def test_manifest_nfc_collisions_fail_before_artifact_reads(bundle, monkeypatch):
    manifest = json.loads((bundle / "bundle_manifest.json").read_text(encoding="utf-8"))
    manifest["required_artifacts"] = ["REPORT.html", "café.txt", "cafe\u0301.txt"]
    (bundle / "bundle_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")

    def fail_if_hashed(*args, **kwargs):
        raise AssertionError("normalization collision must fail before artifact reads")

    monkeypatch.setattr(contract, "sha256", fail_if_hashed)
    assert any("NFC/case-fold collision" in issue for issue in contract.validate_files(bundle))


@pytest.mark.parametrize("path", ["a//b", "name.", "CON.txt", "../escape", "a\\b"])
def test_manifest_paths_use_the_portable_path_rule_before_artifact_reads(
    bundle, monkeypatch, path
):
    manifest = json.loads((bundle / "bundle_manifest.json").read_text(encoding="utf-8"))
    manifest["required_artifacts"] = [path]
    manifest["artifacts"] = {path: {"bytes": 0, "sha256": "0" * 64}}
    (bundle / "bundle_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")

    def fail_if_hashed(*args, **kwargs):
        raise AssertionError("nonportable manifest paths must fail before artifact reads")

    monkeypatch.setattr(contract, "sha256", fail_if_hashed)
    assert contract.validate_files(bundle)


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("schema_version", True, "schema version must be an integer"),
        ("schema_version", "2", "schema version must be an integer"),
        ("schema_version", 3, "unsupported bundle schema version"),
        ("status", [], "bundle status must be a string"),
        ("report_generated", "yes", "report_generated must be a boolean"),
        ("required_artifacts", {}, "required_artifacts must be a non-empty list"),
        ("artifacts", [], "artifacts must be an object"),
    ],
)
def test_manifest_shape_fails_before_artifact_reads(bundle, monkeypatch, field, value, message):
    manifest = json.loads((bundle / "bundle_manifest.json").read_text(encoding="utf-8"))
    manifest[field] = value
    (bundle / "bundle_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")

    def fail_if_hashed(*args, **kwargs):
        raise AssertionError("malformed manifest must fail before artifact reads")

    monkeypatch.setattr(contract, "sha256", fail_if_hashed)
    assert any(message in issue for issue in contract.validate_files(bundle))


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ({"status": "pending"}, "bundle is not complete"),
        ({"report_generated": False}, "report was not generated"),
        ({"report_generated": None}, "report_generated must be a boolean"),
    ],
)
def test_completion_fields_fail_before_artifact_reads(bundle, monkeypatch, mutation, message):
    manifest = json.loads((bundle / "bundle_manifest.json").read_text(encoding="utf-8"))
    manifest.update(mutation)
    (bundle / "bundle_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")

    def fail_if_hashed(*args, **kwargs):
        raise AssertionError("completion metadata must fail before artifact reads")

    monkeypatch.setattr(contract, "sha256", fail_if_hashed)
    assert any(message in issue for issue in contract.validate_files(bundle))


@pytest.mark.parametrize("constant", ["NaN", "Infinity", "-Infinity", "1e9999"])
def test_nonfinite_json_numbers_are_machine_readable_issues(bundle, monkeypatch, constant):
    raw = (bundle / "bundle_manifest.json").read_text(encoding="utf-8")
    raw = raw.replace('"schema_version": 2', f'"schema_version": {constant}')
    (bundle / "bundle_manifest.json").write_text(raw, encoding="utf-8")

    def fail_if_hashed(*args, **kwargs):
        raise AssertionError("invalid JSON constants must fail before artifact reads")

    monkeypatch.setattr(contract, "sha256", fail_if_hashed)
    audit = contract.audit_files(bundle)
    assert audit["valid"] is False
    assert audit["assurance"] == "unavailable"
    assert any("non-finite JSON number" in issue for issue in audit["issues"])


def test_excessive_json_nesting_is_rejected_before_artifact_reads(bundle, monkeypatch):
    manifest = json.loads((bundle / "bundle_manifest.json").read_text(encoding="utf-8"))
    nested = None
    for _ in range(contract._MAX_JSON_DEPTH + 1):
        nested = [nested]
    manifest["extra"] = nested
    (bundle / "bundle_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")

    def fail_if_hashed(*args, **kwargs):
        raise AssertionError("excessive JSON nesting must fail before artifact reads")

    monkeypatch.setattr(contract, "sha256", fail_if_hashed)
    assert any("maximum JSON nesting depth" in issue for issue in contract.validate_files(bundle))


def test_json_nesting_limit_allows_exactly_64_container_levels(bundle):
    manifest = json.loads((bundle / "bundle_manifest.json").read_text(encoding="utf-8"))
    nested = None
    for _ in range(contract._MAX_JSON_DEPTH):
        nested = [nested]
    manifest["extra"] = nested
    (bundle / "bundle_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    assert contract.validate_files(bundle) == []


def test_manifest_unicode_decode_failure_is_a_machine_readable_issue(bundle):
    (bundle / "bundle_manifest.json").write_bytes(b"\xff")
    audit = contract.audit_files(bundle)
    assert audit["valid"] is False
    assert audit["assurance"] == "unavailable"
    assert audit["issues"]


def test_manifest_parser_recursion_error_is_a_machine_readable_issue(
    bundle, monkeypatch
):
    def fail(*args, **kwargs):
        raise RecursionError("synthetic parser recursion")

    monkeypatch.setattr(contract.json, "loads", fail)
    audit = contract.audit_files(bundle)
    assert audit["valid"] is False
    assert audit["assurance"] == "unavailable"
    assert any("synthetic parser recursion" in issue for issue in audit["issues"])


def test_malformed_v2_metadata_does_not_claim_checksum_assurance(bundle):
    manifest = json.loads((bundle / "bundle_manifest.json").read_text(encoding="utf-8"))
    manifest["artifacts"]["result.txt"]["bytes"] = "synthetic"
    (bundle / "bundle_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    audit = contract.audit_files(bundle)
    assert audit["valid"] is False
    assert audit["manifest_version"] == 2
    assert audit["assurance"] == "unavailable"
    assert audit["hashes_checked"] == 0


def test_missing_report_generated_fails_before_artifact_reads(bundle, monkeypatch):
    manifest = json.loads((bundle / "bundle_manifest.json").read_text(encoding="utf-8"))
    manifest.pop("report_generated")
    (bundle / "bundle_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")

    def fail_if_hashed(*args, **kwargs):
        raise AssertionError("missing report flag must fail before artifact reads")

    monkeypatch.setattr(contract, "sha256", fail_if_hashed)
    assert any(
        "report_generated must be a boolean" in issue
        for issue in contract.validate_files(bundle)
    )


def test_unsupported_integer_manifest_version_is_reported_without_assurance(bundle):
    manifest = json.loads((bundle / "bundle_manifest.json").read_text(encoding="utf-8"))
    manifest["schema_version"] = 3
    (bundle / "bundle_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    audit = contract.audit_files(bundle)
    assert audit["valid"] is False
    assert audit["manifest_version"] == 3
    assert audit["assurance"] == "unavailable"
    assert audit["hashes_checked"] == 0
    assert any("unsupported bundle schema version: 3" in issue for issue in audit["issues"])


@pytest.mark.parametrize(
    "record",
    [
        None,
        [],
        {"bytes": True, "sha256": "0" * 64},
        {"bytes": -1, "sha256": "0" * 64},
        {"bytes": 1, "sha256": "not-a-sha256"},
        {"bytes": 1, "sha256": "A" * 64},
    ],
)
def test_malformed_artifact_records_fail_before_artifact_reads(bundle, monkeypatch, record):
    manifest = json.loads((bundle / "bundle_manifest.json").read_text(encoding="utf-8"))
    manifest["artifacts"]["result.txt"] = record
    (bundle / "bundle_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")

    def fail_if_hashed(*args, **kwargs):
        raise AssertionError("malformed artifact record must fail before artifact reads")

    monkeypatch.setattr(contract, "sha256", fail_if_hashed)
    assert contract.validate_files(bundle)


def test_duplicate_manifest_keys_are_rejected_before_artifact_reads(bundle, monkeypatch):
    (bundle / "bundle_manifest.json").write_text(
        '{"schema_version":2,"schema_version":1,"status":"complete",'
        '"required_artifacts":["result.txt"]}',
        encoding="utf-8",
    )

    def fail_if_hashed(*args, **kwargs):
        raise AssertionError("ambiguous manifest must fail before artifact reads")

    monkeypatch.setattr(contract, "sha256", fail_if_hashed)
    assert any("duplicate key" in issue for issue in contract.validate_files(bundle))


def test_legacy_and_v2_audits_report_distinct_assurance(bundle):
    v2 = contract.audit_files(bundle)
    assert v2 == {
        "valid": True,
        "scope": "files only",
        "assurance": "v2-complete-checksums",
        "manifest_version": 2,
        "hashes_checked": 2,
        "html_checks_performed": True,
        "issues": [],
    }

    manifest = json.loads((bundle / "bundle_manifest.json").read_text(encoding="utf-8"))
    manifest.pop("schema_version")
    manifest.pop("artifacts")
    (bundle / "bundle_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    legacy = contract.audit_files(bundle)
    assert legacy == {
        "valid": True,
        "scope": "files only",
        "assurance": "legacy-existence",
        "manifest_version": 1,
        "hashes_checked": 0,
        "html_checks_performed": True,
        "issues": [],
    }
    assert contract.validate_files(bundle) == []


def test_legacy_manifest_does_not_authenticate_an_optional_checksum_map(bundle):
    manifest = json.loads((bundle / "bundle_manifest.json").read_text(encoding="utf-8"))
    manifest.pop("schema_version")
    manifest["artifacts"]["result.txt"] = {"bytes": 0, "sha256": "0" * 64}
    (bundle / "bundle_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    audit = contract.audit_files(bundle)
    assert audit["valid"] is True
    assert audit["assurance"] == "legacy-existence"
    assert audit["hashes_checked"] == 0


@pytest.mark.parametrize("version", [1, 2])
def test_no_report_mode_rejects_true_report_flag_when_html_is_absent(bundle, version):
    (bundle / "REPORT.html").unlink()
    manifest_path = bundle / "bundle_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["required_artifacts"].remove("REPORT.html")
    manifest["report_generated"] = True
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    _set_manifest_version(bundle, version)

    audit = contract.audit_files(bundle, require_report=False)

    assert audit["valid"] is False
    assert audit["assurance"] == "unavailable"
    assert any(
        "report_generated is true but REPORT.html is missing" in issue
        for issue in audit["issues"]
    )


@pytest.mark.parametrize("version", [1, 2])
def test_no_report_mode_audits_present_html_and_rejects_false_report_flag(
    bundle, version
):
    _replace_report_with_valid_checksum(
        bundle,
        "<!--><img src='https://example.invalid/hidden.png'>-->",
    )
    manifest_path = bundle / "bundle_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["report_generated"] = False
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    _set_manifest_version(bundle, version)

    audit = contract.audit_files(bundle, require_report=False)

    assert audit["valid"] is False
    assert audit["assurance"] == "unavailable"
    assert audit["html_checks_performed"] is True
    assert any(
        "report_generated is false but REPORT.html is present" in issue
        for issue in audit["issues"]
    )
    assert any("comment, CDATA, or declaration" in issue for issue in audit["issues"])


@pytest.mark.parametrize("version", [1, 2])
def test_no_report_mode_accepts_declared_machine_readable_only_bundle(bundle, version):
    (bundle / "REPORT.html").unlink()
    manifest_path = bundle / "bundle_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["required_artifacts"].remove("REPORT.html")
    manifest["report_generated"] = False
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    _set_manifest_version(bundle, version)

    audit = contract.audit_files(bundle, require_report=False)

    assert audit["valid"] is True
    assert audit["manifest_version"] == version
    assert audit["html_checks_performed"] is False
    assert audit["assurance"] == (
        "legacy-existence" if version == 1 else "v2-complete-checksums"
    )


@pytest.mark.parametrize("version", [1, 2])
def test_control_manifest_must_use_exact_case_before_json_parse(
    bundle, version, monkeypatch
):
    _set_manifest_version(bundle, version)
    canonical = bundle / "bundle_manifest.json"
    temporary = bundle / "manifest-case-change.tmp"
    noncanonical = bundle / "BUNDLE_MANIFEST.JSON"
    canonical.rename(temporary)
    temporary.rename(noncanonical)

    def fail_if_parsed(*args, **kwargs):
        raise AssertionError("noncanonical control manifest must fail before JSON parse")

    monkeypatch.setattr(contract.json, "loads", fail_if_parsed)
    audit = contract.audit_files(bundle, require_report=False)

    assert audit["valid"] is False
    assert audit["manifest_version"] is None
    assert audit["hashes_checked"] == 0
    assert any("must be named exactly" in issue for issue in audit["issues"])


def test_inventory_rejects_noncanonical_control_manifest_spelling(bundle):
    canonical = bundle / "bundle_manifest.json"
    temporary = bundle / "manifest-case-change.tmp"
    noncanonical = bundle / "BUNDLE_MANIFEST.JSON"
    canonical.rename(temporary)
    temporary.rename(noncanonical)

    with pytest.raises(ValueError, match="must be named exactly"):
        contract.inventory(bundle)


def test_v2_requires_checksum_records_for_every_file(bundle):
    manifest = json.loads((bundle / "bundle_manifest.json").read_text(encoding="utf-8"))
    manifest["artifacts"].pop("result.txt")
    (bundle / "bundle_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    issues = contract.validate_files(bundle)
    assert any("no checksum record" in issue for issue in issues)


def test_linked_descendant_is_rejected_before_artifact_reads(bundle, monkeypatch):
    original_is_link = contract._is_link

    def synthetic_link(path):
        return path.name == "result.txt" or original_is_link(path)

    monkeypatch.setattr(contract, "_is_link", synthetic_link)

    def fail_if_hashed(*args, **kwargs):
        raise AssertionError("bundle links must fail before artifact reads")

    monkeypatch.setattr(contract, "sha256", fail_if_hashed)
    issues = contract.validate_files(bundle)
    assert any("linked path" in issue for issue in issues)


def test_artifact_read_failure_is_a_machine_readable_issue(bundle, monkeypatch):
    original_sha256 = contract.sha256

    def deny_one(path):
        if path.name == "result.txt":
            raise PermissionError("read denied")
        return original_sha256(path)

    monkeypatch.setattr(contract, "sha256", deny_one)
    audit = contract.audit_files(bundle)
    assert audit["valid"] is False
    assert audit["assurance"] == "unavailable"
    assert audit["hashes_checked"] == 1
    assert any("cannot inspect artifact result.txt" in issue for issue in audit["issues"])


def test_checksum_mismatch_does_not_claim_complete_assurance(bundle):
    (bundle / "result.txt").write_text("tampered", encoding="utf-8")
    audit = contract.audit_files(bundle)
    assert audit["valid"] is False
    assert audit["assurance"] == "unavailable"
    assert audit["hashes_checked"] == 2
    assert audit["html_checks_performed"] is True
    assert any("checksum mismatch" in issue for issue in audit["issues"])


def test_html_decode_failure_does_not_mark_html_check_complete(bundle):
    (bundle / "REPORT.html").write_bytes(b"\xff")
    manifest = json.loads((bundle / "bundle_manifest.json").read_text(encoding="utf-8"))
    manifest["artifacts"] = contract.inventory(bundle)
    (bundle / "bundle_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    audit = contract.audit_files(bundle)
    assert audit["valid"] is False
    assert audit["assurance"] == "unavailable"
    assert audit["html_checks_performed"] is False
    assert any("cannot inspect REPORT.html" in issue for issue in audit["issues"])


def test_html_path_recheck_failure_does_not_mark_html_check_complete(bundle, monkeypatch):
    original_entries = contract._bundle_entries
    calls = 0

    def fail_html_recheck(root):
        nonlocal calls
        calls += 1
        if calls == 2:
            return [], [], ["synthetic HTML path recheck failure"]
        return original_entries(root)

    monkeypatch.setattr(contract, "_bundle_entries", fail_html_recheck)
    audit = contract.audit_files(bundle)
    assert audit["valid"] is False
    assert audit["assurance"] == "unavailable"
    assert audit["html_checks_performed"] is False
    assert "synthetic HTML path recheck failure" in audit["issues"]


def test_unsupported_html_marks_check_performed_without_claiming_assurance(bundle):
    (bundle / "REPORT.html").write_text("<script>alert(1)</script>", encoding="utf-8")
    manifest = json.loads((bundle / "bundle_manifest.json").read_text(encoding="utf-8"))
    manifest["artifacts"] = contract.inventory(bundle)
    (bundle / "bundle_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    audit = contract.audit_files(bundle)
    assert audit["valid"] is False
    assert audit["assurance"] == "unavailable"
    assert audit["html_checks_performed"] is True
    assert any("unsupported active element" in issue for issue in audit["issues"])


def test_extra_file_is_detected_before_artifact_reads(bundle, monkeypatch):
    (bundle / "stale.txt").write_text("leftover", encoding="utf-8")

    def fail_if_hashed(*args, **kwargs):
        raise AssertionError("invalid filesystem inventory must fail before artifact reads")

    monkeypatch.setattr(contract, "sha256", fail_if_hashed)
    audit = contract.audit_files(bundle)
    assert audit["assurance"] == "unavailable"
    assert audit["hashes_checked"] == 0
    assert audit["html_checks_performed"] is False
    assert any("exactly" in issue for issue in audit["issues"])


def _record():
    dossier = {"tf_id": "Synthetic_report", "family": "test", "rescan": {"hits": []}, "regulon": []}
    primary = {"source": "test", "sequence": "", "score": 0, "provenance": {}, "status": "test"}
    return dossier, {"primary": primary, "ranked": [primary], "candidates": []}


def test_render_failure_does_not_publish_a_completed_bundle(tmp_path, monkeypatch):
    dossier, operators = _record()
    monkeypatch.setattr(outputs, "write_final_operators", lambda *args: 0)
    def fail(*args, **kwargs):
        raise RuntimeError("render failure")
    monkeypatch.setattr(render_report, "render", fail)
    with pytest.raises(RuntimeError, match="render failure"):
        outputs.write_stage_a(dossier, operators=operators, root=tmp_path)
    assert list(tmp_path.iterdir()) == []


def test_writer_publishes_a_valid_self_contained_synthetic_bundle(tmp_path, monkeypatch):
    dossier, operators = _record()
    monkeypatch.setattr(outputs, "write_final_operators", lambda *args: 0)
    monkeypatch.setattr(outputs, "runtime_provenance", lambda: {"source": "synthetic"})
    monkeypatch.setattr(render_report, "_gen_figures", lambda *args: {})
    result = outputs.write_stage_a(dossier, operators=operators, root=tmp_path)
    assert contract.validate_files(result) == []
    text = (result / "REPORT.html").read_text(encoding="utf-8")
    assert "lang='en'" in text and "Files &amp; provenance" in text
    assert "README.md" in text and "prefers-reduced-motion" in text


def test_writer_packages_exact_run_inputs_and_rewrites_cache_paths(tmp_path, monkeypatch):
    dossier, operators = _record()
    dossier.update({
        "tf_record": {"sequence": "MPEPTIDE"},
        "homolog_regions": {"pre_msa": ["AACCGG"], "expanded": ["AACCGG", "TTGGCC"]},
        "seed_selection": {"chosen_w": 6, "trials": [{"cluster_id": 0}],
                           "construction": {"top_k": 8}, "parameters": {"rescan_pvalue_thresh": 1e-4}},
    })
    cache = tmp_path / "outside-cache"
    cache.mkdir()
    cif = cache / "cached-model.cif"
    cif.write_text("data_synthetic\n", encoding="utf-8")
    msa = cache / "cached-alignment.a3m"
    msa.write_text(">query\nMPEPTIDE\n", encoding="utf-8")
    source_context = {
        "run_id": "run-123",
        "started_at_utc": "2026-09-09T12:00:00+00:00",
        "requested": {"name": "Synthetic_report", "family": "test"},
        "effective_config": {"allow_online": False, "allow_folds": True},
        "query": {"sequence": "MPEPTIDE", "input_kind": "raw_sequence"},
        "genome": {
            "accession": "SYNTHETIC_CONTIG", "sequence": "AACCGGTTAACC",
            "genes": [{"start": 1, "end": 8, "strand": "+", "name": "geneA",
                       "locus_tag": "LT1", "protein_id": "P1", "product": "synthetic"}],
            "contig_offsets": {"SYNTHETIC_CONTIG": 0}, "circular": False,
        },
        "windows": {"discovery_window": {"start": 0, "end": 6, "strand": "+", "sequence": "AACCGG"}},
        # `attempts` present, MSA with no `source`: the shape of 190 of the 238 canonical run records.
        # Without attempts the report's homolog ledger is skipped entirely, which is how a render bug
        # that killed every such bundle at the output stage passed this test (dd24d85).
        "homolog_collection": {"selected_route": "msa",
                               "attempts": [{"route": "ssn_cluster", "selected": True, "n_regions": 8}],
                               "msa": {"alignment_path": str(msa)}},
    }
    struct = {"folded": True, "backend": "synthetic", "dimer_cif": str(cif),
              "afdb_monomer_cif": None, "flags": []}
    monkeypatch.setattr(outputs, "write_final_operators", lambda *args: 0)
    monkeypatch.setattr(render_report, "_gen_figures", lambda *args: {})

    result = outputs.write_stage_a(
        dossier, operators=operators, struct=struct, source_context=source_context,
        root=tmp_path / "jobs",
    )

    expected = {
        "input/query.fasta", "input/run_parameters.json", "genome/genome.fna",
        "genome/genes.tsv", "genome/discovery_window.fasta",
        "homologs/source_alignment.a3m", "homologs/source_record.json",
        "motif/seed_selection.json", "run/run_record.json", "run/evidence_index.json",
        "structure/apo_oligomer.cif",
    }
    assert expected <= {p.relative_to(result).as_posix() for p in result.rglob("*") if p.is_file()}
    portable = json.loads((result / "dossier.json").read_text(encoding="utf-8"))
    assert portable["structure"]["dimer_cif"] == "structure/apo_oligomer.cif"
    run_record = json.loads((result / "run/run_record.json").read_text(encoding="utf-8"))
    assert run_record["query"]["length"] == 8
    assert run_record["genome"]["gene_count"] == 1
    assert run_record["stages"] and run_record["evidence_index"]
    serialized = "\n".join(p.read_text(encoding="utf-8", errors="ignore") for p in result.rglob("*") if p.is_file())
    assert str(cache) not in serialized
    html = (result / "REPORT.html").read_text(encoding="utf-8")
    assert "Run record &amp; traceability" in html
    assert "Pipeline stage ledger" in html
    assert "Evidence-to-artifact map" in html
    assert contract.validate_files(result) == []


def test_writer_rejects_unmapped_machine_local_paths(tmp_path, monkeypatch):
    dossier, operators = _record()
    dossier["unmapped_cache_record"] = {"path": str(tmp_path / "outside" / "cached.json")}
    monkeypatch.setattr(outputs, "write_final_operators", lambda *args: 0)
    with pytest.raises(ValueError, match="machine-local path"):
        outputs.write_stage_a(dossier, operators=operators, root=tmp_path / "jobs", render=False)
    assert not (tmp_path / "jobs" / "Synthetic_report").exists()


def test_report_requires_a_valid_dossier(tmp_path):
    (tmp_path / "dossier.json").write_text("[]", encoding="utf-8")
    with pytest.raises(ValueError, match="tf_id"):
        render_report.render(tmp_path)


def test_report_rejects_a_corrupt_existing_file_contract(tmp_path):
    (tmp_path / "bundle_manifest.json").write_text("{broken", encoding="utf-8")
    with pytest.raises(ValueError, match="valid object"):
        render_report.render(tmp_path)


def test_tabular_annotations_round_trip_without_extra_rows(tmp_path):
    import csv
    path = tmp_path / "annotations.tsv"
    value = 'an annotation\twith "quotes"\nand a newline'
    outputs._tsv(path, ["id", "annotation"], [[1, value]])
    with path.open(encoding="utf-8", newline="") as handle:
        assert list(csv.DictReader(handle, delimiter="\t")) == [{"id": "1", "annotation": value}]


def test_figure_export_closes_canvas_after_failure(tmp_path, monkeypatch):
    from predictor.report import figures
    fig = figures.plt.figure()
    def fail(*args, **kwargs):
        raise OSError("export failed")
    monkeypatch.setattr(fig, "savefig", fail)
    with pytest.raises(OSError, match="export failed"):
        figures._save(fig, Path(tmp_path), "synthetic")
    assert not figures.plt.fignum_exists(fig.number)


@pytest.fixture
def published_report(tmp_path, monkeypatch):
    from predictor import provenance
    from predictor.report import environment
    dossier, operators = _record()
    monkeypatch.setattr(outputs, "write_final_operators", lambda *args: 0)
    monkeypatch.setattr(outputs, "runtime_provenance", lambda: {"source": "original"})
    monkeypatch.setattr(provenance, "runtime_provenance", lambda: {"source": "synthetic report"})
    monkeypatch.setattr(environment, "snapshot", lambda: {"scope": "synthetic export", "marker": "original"})
    monkeypatch.setattr(render_report, "_gen_figures", lambda *args: {})
    return outputs.write_stage_a(dossier, operators=operators, root=tmp_path)


def test_new_report_links_the_manifest_and_environment(published_report):
    text = (published_report / "REPORT.html").read_text(encoding="utf-8")
    assert 'href="bundle_manifest.json"' in text
    assert 'href="report_environment.json"' in text
    assert "Browse all" in text
    assert contract.validate_files(published_report) == []


def test_failed_regeneration_preserves_every_original_file(published_report, monkeypatch):
    before = {p.relative_to(published_report): p.read_bytes() for p in published_report.rglob("*") if p.is_file()}
    def fail(data, figures):
        figures.mkdir()
        (figures / "partial.png").write_bytes(b"partial")
        raise RuntimeError("synthetic rendering failure")
    monkeypatch.setattr(render_report, "_gen_figures", fail)
    with pytest.raises(RuntimeError, match="synthetic rendering failure"):
        render_report.render(published_report)
    assert before == {p.relative_to(published_report): p.read_bytes() for p in published_report.rglob("*") if p.is_file()}
    assert not list(published_report.parent.glob(".report-stage-*"))
    assert contract.validate_files(published_report) == []


def test_regeneration_preserves_prediction_provenance_and_previous_folder(published_report, monkeypatch):
    from predictor.report import environment
    monkeypatch.setattr(environment, "snapshot", lambda: {"scope": "synthetic export", "marker": "new"})
    result = render_report.render(published_report)
    assert result == published_report / "REPORT.html"
    assert contract.validate_files(published_report) == []
    manifest = json.loads((published_report / "bundle_manifest.json").read_text(encoding="utf-8"))
    assert manifest["runtime"] == {"source": "original"}
    assert manifest["report_environment_history"][0]["marker"] == "original"
    assert len(list(published_report.parent.glob(".*.report-previous-*"))) == 1


def test_regeneration_migrates_a_verified_legacy_report(published_report):
    legacy_html = "<html><style>body{background:linear-gradient(#fff,#eee)}</style></html>"
    _replace_report_with_valid_checksum(published_report, legacy_html)
    assert any("unsupported CSS" in issue for issue in contract.validate_files(published_report))

    render_report.render(published_report)

    assert contract.validate_files(published_report) == []
    assert "linear-gradient" not in (published_report / "REPORT.html").read_text(encoding="utf-8")
    assert (published_report / "README.md").is_file()


def test_failed_report_promotion_restores_the_previous_folder(published_report, monkeypatch):
    original_replace = Path.replace
    def fail_promotion(path, target):
        if path.name.startswith(".report-stage-"):
            raise OSError("promotion denied")
        return original_replace(path, target)
    monkeypatch.setattr(Path, "replace", fail_promotion)
    with pytest.raises(OSError, match="promotion denied"):
        render_report.render(published_report)
    assert contract.validate_files(published_report) == []
    assert not list(published_report.parent.glob(".report-stage-*"))


def test_installed_file_audit_cli_returns_machine_readable_failure(tmp_path, capsys):
    assert contract.main([str(tmp_path)]) == 1
    audit = json.loads(capsys.readouterr().out)
    assert audit["valid"] is False
    assert audit["assurance"] == "unavailable"
    assert audit["manifest_version"] is None
    assert audit["hashes_checked"] == 0
    assert audit["html_checks_performed"] is False
    assert audit["issues"]


def test_installed_file_audit_cli_reports_v2_assurance(bundle, capsys):
    assert contract.main([str(bundle)]) == 0
    audit = json.loads(capsys.readouterr().out)
    assert audit["valid"] is True
    assert audit["assurance"] == "v2-complete-checksums"
    assert audit["manifest_version"] == 2
    assert audit["hashes_checked"] == 2
    assert audit["html_checks_performed"] is True
    assert audit["issues"] == []


@pytest.mark.parametrize("version", [1, 2])
def test_no_report_cli_never_skips_a_present_report(bundle, capsys, version):
    _replace_report_with_valid_checksum(
        bundle,
        "<![CDATA[><script src='https://example.invalid/hidden.js'></script>]]>",
    )
    _set_manifest_version(bundle, version)

    assert contract.main([str(bundle), "--no-report"]) == 1
    audit = json.loads(capsys.readouterr().out)
    assert audit["valid"] is False
    assert audit["html_checks_performed"] is True
    assert any("comment, CDATA, or declaration" in issue for issue in audit["issues"])


def test_interrupted_publication_removes_staging(tmp_path, monkeypatch):
    dossier, operators = _record()
    monkeypatch.setattr(outputs, "write_final_operators", lambda *args: 0)
    def interrupt(*args, **kwargs):
        raise KeyboardInterrupt()
    monkeypatch.setattr(render_report, "render", interrupt)
    with pytest.raises(KeyboardInterrupt):
        outputs.write_stage_a(dossier, operators=operators, root=tmp_path)
    assert list(tmp_path.iterdir()) == []


def test_machine_readable_export_has_a_valid_inventory(tmp_path, monkeypatch):
    dossier, operators = _record()
    monkeypatch.setattr(outputs, "write_final_operators", lambda *args: 0)
    monkeypatch.setattr(outputs, "runtime_provenance", lambda: {"source": "synthetic"})
    result = outputs.write_stage_a(dossier, operators=operators, root=tmp_path, render=False)
    assert contract.validate_files(result, require_report=False) == []
    assert not (result / "REPORT.html").exists()
    assert "--no-report" in (result / "README.md").read_text(encoding="utf-8")


def test_corrupt_optional_metadata_is_not_treated_as_absent(tmp_path):
    (tmp_path / "dossier.json").write_text(json.dumps(_record()[0]), encoding="utf-8")
    (tmp_path / "tf").mkdir()
    (tmp_path / "tf/tf_summary.json").write_text("{invalid", encoding="utf-8")
    with pytest.raises(ValueError, match="tf_summary.json must contain valid"):
        render_report.collect(tmp_path)


def test_optional_metadata_read_failure_propagates(tmp_path, monkeypatch):
    def fail(*args, **kwargs):
        raise PermissionError("read denied")
    monkeypatch.setattr(Path, "read_text", fail)
    with pytest.raises(PermissionError, match="read denied"):
        render_report._load(tmp_path / "metadata.json")
