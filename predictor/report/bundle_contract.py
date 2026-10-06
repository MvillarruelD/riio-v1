"""Filesystem-only checks for portable, complete report bundles.

This module deliberately does not interpret scientific records. It verifies paths,
checksums and a deliberately small local-resource surface in ``REPORT.html``.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import re
import unicodedata
from html.parser import HTMLParser
from pathlib import Path
from typing import Any
from urllib.parse import unquote, urlsplit


_WINDOWS_RESERVED = {"CON", "PRN", "AUX", "NUL", "CONIN$", "CONOUT$"}
_WINDOWS_RESERVED.update(f"COM{number}" for number in range(1, 10))
_WINDOWS_RESERVED.update(f"LPT{number}" for number in range(1, 10))
_WINDOWS_RESERVED.update(f"COM{number}" for number in "¹²³")
_WINDOWS_RESERVED.update(f"LPT{number}" for number in "¹²³")
_WINDOWS_INVALID = frozenset('<>:"\\|?*')
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_BAD_PERCENT_ESCAPE = re.compile(r"%(?![0-9a-fA-F]{2})")
_ENCODED_SEPARATOR = re.compile(r"%(?:2f|5c)", re.IGNORECASE)
_CSS_FUNCTION = re.compile(r"(?<![@\w-])([a-z_-][\w-]*)\s*\(", re.IGNORECASE)
_CSS_AT_RULE = re.compile(r"@([a-z][\w-]*)", re.IGNORECASE)
_SAFE_CSS_FUNCTIONS = {
    "blur",
    "clamp",
    "minmax",
    "nth-child",
    "repeat",
    "rgba",
    "var",
}
_SAFE_CSS_AT_RULES = {"media"}
_UNSUPPORTED_ACTIVE_TAGS = {"applet", "embed", "iframe", "object", "script"}
_SUPPORTED_HTML_TAGS = {
    "a",
    "aside",
    "b",
    "body",
    "br",
    "details",
    "div",
    "h1",
    "h2",
    "h3",
    "head",
    "header",
    "html",
    "i",
    "img",
    "li",
    "main",
    "meta",
    "nav",
    "p",
    "small",
    "span",
    "style",
    "strong",
    "summary",
    "table",
    "tbody",
    "td",
    "th",
    "thead",
    "title",
    "tr",
    "ul",
}
_GLOBAL_HTML_ATTRIBUTES = {
    "class",
    "dir",
    "hidden",
    "id",
    "lang",
    "role",
    "style",
    "tabindex",
    "title",
}
_TAG_HTML_ATTRIBUTES = {
    "a": {"href"},
    "details": {"open"},
    "img": {"alt", "height", "loading", "src", "width"},
    "meta": {"charset", "content", "name"},
    "td": {"colspan", "headers", "rowspan"},
    "th": {"colspan", "headers", "rowspan", "scope"},
}
_RESOURCE_ATTRIBUTES = {
    "action",
    "archive",
    "background",
    "cite",
    "classid",
    "codebase",
    "cursor",
    "data",
    "dynsrc",
    "fill",
    "filter",
    "formaction",
    "href",
    "imagesrcset",
    "longdesc",
    "lowsrc",
    "manifest",
    "marker",
    "marker-end",
    "marker-mid",
    "marker-start",
    "mask",
    "ping",
    "poster",
    "profile",
    "src",
    "srcset",
    "usemap",
    "xlink:href",
}
_CSS_RESOURCE_FUNCTION = re.compile(
    r"(?<![\w-])(?:url|image-set|-webkit-image-set|src|cross-fade|element|paint)\s*\(",
    re.IGNORECASE,
)
_CSS_RESOURCE_AT_RULE = re.compile(r"@(?:import|font-face|namespace)\b", re.IGNORECASE)
_MAX_JSON_DEPTH = 64
_OS_PATH_ISRESERVED = getattr(os.path, "isreserved", None)
_MANIFEST_NAME = "bundle_manifest.json"
_SIMPLE_HTML_DOCTYPE = re.compile(
    r"\A(?:\ufeff)?[\t\n\f\r ]*<!doctype[\t\n\f\r ]+html[\t\n\f\r ]*>",
    re.IGNORECASE,
)


def _fallback_windows_reserved(part: str) -> bool:
    """Match Windows device names on Python versions without ``isreserved``."""
    device_stem = part.split(".", 1)[0].rstrip(" .").upper()
    return device_stem in _WINDOWS_RESERVED


def _windows_reserved(part: str) -> bool:
    """Use the runtime check when present plus a stable Python-3.10 fallback."""
    runtime_reserved = (
        bool(_OS_PATH_ISRESERVED(part)) if _OS_PATH_ISRESERVED is not None else False
    )
    return runtime_reserved or _fallback_windows_reserved(part)


def _portable_parts(relative: str) -> tuple[str, ...]:
    """Return canonical POSIX path parts without consulting the filesystem."""
    if not isinstance(relative, str) or not relative:
        raise ValueError(f"invalid bundle path: {relative!r}")
    parts = relative.split("/")
    if any(not part for part in parts):
        raise ValueError(f"bundle path has a repeated or trailing separator: {relative!r}")
    if any(part in (".", "..") for part in parts):
        raise ValueError(f"bundle path escapes its root: {relative!r}")
    for part in parts:
        if part.endswith((".", " ")):
            raise ValueError(f"bundle path has a trailing dot or space: {relative!r}")
        if any(character in _WINDOWS_INVALID or ord(character) < 32 for character in part):
            raise ValueError(f"bundle path contains a nonportable character: {relative!r}")
        if _windows_reserved(part):
            raise ValueError(f"bundle path uses a reserved Windows basename: {relative!r}")
    return tuple(parts)


def _is_link(path: Path) -> bool:
    return path.is_symlink() or getattr(path, "is_junction", lambda: False)()


def contained_file(root: Path, relative: str) -> Path:
    """Resolve one canonical portable path, rejecting traversal and linked descendants."""
    parts = _portable_parts(relative)
    root = Path(root)
    if _is_link(root):
        raise ValueError("bundle root must not be a linked path")
    path = root
    for part in parts:
        path /= part
        if _is_link(path):
            raise ValueError(f"bundle path contains a linked descendant: {relative!r}")
    resolved_root = root.resolve()
    if not path.resolve().is_relative_to(resolved_root):
        raise ValueError(f"bundle path escapes its root: {relative!r}")
    return path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _portable_collision_key(relative: str) -> str:
    return "/".join(
        unicodedata.normalize("NFC", unicodedata.normalize("NFC", part).casefold())
        for part in relative.split("/")
    )


def _casefold_collision_issues(paths: list[str], *, context: str) -> list[str]:
    spellings: dict[str, str] = {}
    issues = []
    for relative in paths:
        folded = _portable_collision_key(relative)
        prior = spellings.setdefault(folded, relative)
        if prior != relative:
            issues.append(
                f"{context} has an NFC/case-fold collision: {prior!r} and {relative!r}"
            )
    return issues


def _bundle_entries(root: Path) -> tuple[list[str], list[str], list[str]]:
    """Inventory names without following a linked entry or opening any file body."""
    root = Path(root)
    if _is_link(root):
        return [], [], ["bundle root must not be a linked path"]
    if not root.is_dir():
        return [], [], [f"bundle directory does not exist: {root}"]

    files: list[str] = []
    all_entries: list[str] = []
    issues: list[str] = []

    def walk(directory: Path) -> None:
        with os.scandir(directory) as scan:
            entries = sorted(scan, key=lambda entry: entry.name)
        for entry in entries:
            path = Path(entry.path)
            relative = path.relative_to(root).as_posix()
            all_entries.append(relative)
            try:
                _portable_parts(relative)
            except ValueError as exc:
                issues.append(str(exc))
            if entry.is_symlink() or _is_link(path):
                issues.append(f"bundle contains a linked path: {relative}")
                continue
            if entry.is_dir(follow_symlinks=False):
                walk(path)
            elif entry.is_file(follow_symlinks=False):
                files.append(relative)
            else:
                issues.append(f"bundle contains a non-regular path: {relative}")

    try:
        walk(root)
    except OSError as exc:
        issues.append(f"cannot inspect bundle paths: {exc}")
    issues.extend(_casefold_collision_issues(all_entries, context="bundle"))
    return files, all_entries, sorted(set(issues))


def _control_manifest_issues(
    files: list[str], all_entries: list[str], *, required: bool
) -> list[str]:
    """Require one exact control-manifest spelling without opening its body."""
    manifest_key = _portable_collision_key(_MANIFEST_NAME)
    spellings = [
        relative
        for relative in all_entries
        if _portable_collision_key(relative) == manifest_key
    ]
    issues = [
        f"bundle control manifest must be named exactly {_MANIFEST_NAME!r}, "
        f"not {relative!r}"
        for relative in spellings
        if relative != _MANIFEST_NAME
    ]
    if _MANIFEST_NAME in all_entries and _MANIFEST_NAME not in files:
        issues.append(f"bundle control manifest is not a regular file: {_MANIFEST_NAME}")
    elif required and _MANIFEST_NAME not in files:
        issues.append(f"bundle control manifest is missing: {_MANIFEST_NAME}")
    return issues


def inventory(root: Path) -> dict[str, dict[str, int | str]]:
    """Hash every regular artifact except the manifest that contains these hashes."""
    root = Path(root)
    files, all_entries, issues = _bundle_entries(root)
    issues.extend(_control_manifest_issues(files, all_entries, required=False))
    if issues:
        raise ValueError(issues[0])
    result = {}
    for relative in files:
        if relative == _MANIFEST_NAME:
            continue
        checked = contained_file(root, relative)
        result[relative] = {"bytes": checked.stat().st_size, "sha256": sha256(checked)}
    return result


class _Resources(HTMLParser):
    """Collect the only supported resource attributes and flag broader HTML surfaces."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.references: list[tuple[str, str, str]] = []
        self.issues: list[str] = []
        self._style_depth = 0

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        tag = tag.casefold()
        attrs = [(name.casefold(), value) for name, value in attrs]
        names = [name for name, _ in attrs]
        duplicate_names = {name for name in names if names.count(name) > 1}
        for name in sorted(duplicate_names):
            self.issues.append(f"REPORT.html has duplicate attribute {tag}[{name}]")
        if tag == "style":
            self._style_depth += 1
        if tag not in _SUPPORTED_HTML_TAGS:
            self.issues.append(f"REPORT.html uses unsupported element <{tag}>")
        if tag in _UNSUPPORTED_ACTIVE_TAGS:
            self.issues.append(f"REPORT.html uses unsupported active element <{tag}>")
        if tag == "base":
            self.issues.append("REPORT.html uses an unsupported document base")
        if "srcset" in names:
            self.issues.append("REPORT.html uses unsupported srcset resources")
        if "srcdoc" in names:
            self.issues.append("REPORT.html uses unsupported srcdoc content")
        for name, value in attrs:
            if name.startswith("on"):
                self.issues.append(
                    f"REPORT.html uses unsupported event-handler attribute {tag}[{name}]"
                )
            allowed = (
                name in _GLOBAL_HTML_ATTRIBUTES
                or name in _TAG_HTML_ATTRIBUTES.get(tag, set())
                or name.startswith("aria-")
                or name.startswith("data-")
            )
            if not allowed:
                kind = "resource" if name in _RESOURCE_ATTRIBUTES else "unknown"
                self.issues.append(
                    f"REPORT.html uses unsupported {kind} attribute {tag}[{name}]"
                )
            if value and (
                _CSS_RESOURCE_FUNCTION.search(value) or _CSS_RESOURCE_AT_RULE.search(value)
            ):
                self.issues.append(
                    f"REPORT.html attribute {tag}[{name}] contains CSS resource syntax"
                )
            if name == "style" and value:
                self.issues.extend(_css_issues(value))
        if tag == "meta" and any(
            name == "http-equiv" and (value or "").casefold() == "refresh"
            for name, value in attrs
        ):
            self.issues.append("REPORT.html uses an unsupported meta refresh")

        for name, value in attrs:
            if not value:
                continue
            if (tag, name) in {("a", "href"), ("img", "src")}:
                self.references.append((tag, name, value))

    def handle_endtag(self, tag: str) -> None:
        if tag.casefold() == "style" and self._style_depth:
            self._style_depth -= 1

    def handle_data(self, data: str) -> None:
        if self._style_depth:
            self.issues.extend(_css_issues(data))


def _css_issues(css: str) -> list[str]:
    """Reject CSS outside the small function/at-rule subset emitted by reports."""
    issues = []
    if "\\" in css or "/*" in css or "*/" in css:
        issues.append("REPORT.html uses unsupported CSS escape or comment syntax")
    functions = {
        match.group(1).casefold()
        for match in _CSS_FUNCTION.finditer(css)
        if match.group(1).casefold() not in _SAFE_CSS_FUNCTIONS
    }
    if functions:
        issues.append(
            "REPORT.html uses unsupported CSS function(s): " + ", ".join(sorted(functions))
        )
    at_rules = {
        match.group(1).casefold()
        for match in _CSS_AT_RULE.finditer(css)
        if match.group(1).casefold() not in _SAFE_CSS_AT_RULES
    }
    if at_rules:
        issues.append(
            "REPORT.html uses unsupported CSS at-rule(s): " + ", ".join(sorted(at_rules))
        )
    return issues


def _html_lexical_issues(html: str) -> list[str]:
    """Reject browser-tokenization surfaces outside the generated HTML subset."""
    remainder = html
    doctype = _SIMPLE_HTML_DOCTYPE.match(remainder)
    if doctype is not None:
        remainder = remainder[doctype.end() :]
    issues = []
    if "<!" in remainder:
        issues.append(
            "REPORT.html uses unsupported HTML comment, CDATA, or declaration syntax"
        )
    if "<?" in html:
        issues.append("REPORT.html uses unsupported processing-instruction syntax")
    return issues


def _decoded_local_reference(reference: str, encoded_path: str) -> str:
    if _BAD_PERCENT_ESCAPE.search(encoded_path):
        raise ValueError(f"report reference has an invalid percent escape: {reference}")
    if _ENCODED_SEPARATOR.search(encoded_path):
        raise ValueError(f"report reference encodes a path separator: {reference}")
    try:
        decoded = unquote(encoded_path, encoding="utf-8", errors="strict")
    except UnicodeDecodeError as exc:
        raise ValueError(f"report reference is not valid UTF-8: {reference}") from exc
    _portable_parts(decoded)
    return decoded


def _report_asset_audit(root: Path) -> tuple[list[str], bool]:
    """Return HTML resource issues and whether decode/parse/check completed."""
    root = Path(root)
    files, _, path_issues = _bundle_entries(root)
    if path_issues:
        return path_issues, False
    file_set = set(files)
    folded_files = {_portable_collision_key(relative): relative for relative in files}
    report_path = contained_file(root, "REPORT.html")
    html = report_path.read_text(encoding="utf-8")
    lexical_issues = _html_lexical_issues(html)
    if lexical_issues:
        return lexical_issues, True
    parser = _Resources()
    parser.feed(html)
    parser.close()
    issues = list(parser.issues)
    for tag, attribute, reference in parser.references:
        try:
            url = urlsplit(reference)
        except ValueError as exc:
            issues.append(f"invalid report reference {reference!r}: {exc}")
            continue
        if url.scheme:
            if tag == "a" and attribute == "href":
                if url.scheme in ("https", "http") and url.netloc:
                    continue
                if url.scheme == "mailto" and url.path and not url.netloc:
                    continue
            if tag == "img" and url.scheme == "data" and reference.casefold().startswith("data:image/"):
                continue
            issues.append(f"report depends on an external or unsafe asset: {reference}")
            continue
        if url.netloc:
            issues.append(f"report depends on an external or unsafe asset: {reference}")
            continue
        if not url.path:
            continue
        try:
            relative = _decoded_local_reference(reference, url.path)
            contained_file(root, relative)
        except ValueError as exc:
            issues.append(str(exc))
            continue
        if relative not in file_set:
            differently_cased = folded_files.get(_portable_collision_key(relative))
            if differently_cased is not None:
                issues.append(
                    "report asset path spelling does not match bundle file: "
                    f"{reference} (file is {differently_cased})"
                )
            else:
                issues.append(f"report asset is missing: {reference}")
    return sorted(set(issues)), True


def report_asset_issues(root: Path) -> list[str]:
    """Check the constrained HTML resource surface without interpreting records."""
    return _report_asset_audit(root)[0]


class _JsonObject(dict[str, Any]):
    def __init__(self, pairs: list[tuple[str, Any]]) -> None:
        super().__init__(pairs)
        seen = set()
        self.duplicate_keys = []
        for key, _ in pairs:
            if key in seen:
                self.duplicate_keys.append(key)
            seen.add(key)


def _json_structure_issues(value: Any, location: str = "bundle manifest") -> list[str]:
    """Find duplicate keys and excessive nesting without recursive traversal."""
    issues = []
    pending = [(value, location, 0)]
    while pending:
        current, current_location, depth = pending.pop()
        if isinstance(current, (_JsonObject, list)) and depth > _MAX_JSON_DEPTH:
            issues.append(
                f"{current_location} exceeds maximum JSON nesting depth {_MAX_JSON_DEPTH}"
            )
            continue
        if isinstance(current, _JsonObject):
            issues.extend(
                f"{current_location} has a duplicate key: {key!r}"
                for key in current.duplicate_keys
            )
            pending.extend(
                (child, f"{current_location}.{key}", depth + 1)
                for key, child in current.items()
            )
        elif isinstance(current, list):
            pending.extend(
                (child, f"{current_location}[{index}]", depth + 1)
                for index, child in enumerate(current)
            )
    return issues


def _reject_json_constant(constant: str) -> None:
    raise ValueError(f"non-finite JSON number is not allowed: {constant}")


def _parse_json_float(value: str) -> float:
    parsed = float(value)
    if not math.isfinite(parsed):
        raise ValueError(f"non-finite JSON number is not allowed: {value}")
    return parsed


def _new_audit() -> dict[str, Any]:
    return {
        "valid": False,
        "scope": "files only",
        "assurance": "unavailable",
        "manifest_version": None,
        "hashes_checked": 0,
        "html_checks_performed": False,
        "issues": [],
    }


def audit_files(
    root: Path,
    *,
    require_report: bool = True,
    _audit_report_html: bool = True,
) -> dict[str, Any]:
    """Return a JSON-ready report of file-integrity assurance and defects.

    Manifest shape and artifact metadata are validated before any declared artifact
    is opened. The complete directory tree is then checked for linked paths before
    a target is read.
    """
    root = Path(root)
    audit = _new_audit()
    try:
        files, all_entries, path_issues = _bundle_entries(root)
    except (OSError, RecursionError, TypeError, UnicodeError, ValueError) as exc:
        audit["issues"] = [f"cannot inspect bundle paths: {exc}"]
        return audit
    path_issues.extend(_control_manifest_issues(files, all_entries, required=True))
    if path_issues:
        audit["issues"] = sorted(set(path_issues))
        return audit
    try:
        manifest_path = contained_file(root, _MANIFEST_NAME)
        manifest = json.loads(
            manifest_path.read_text(encoding="utf-8"),
            object_pairs_hook=_JsonObject,
            parse_constant=_reject_json_constant,
            parse_float=_parse_json_float,
        )
    except (json.JSONDecodeError, OSError, RecursionError, TypeError, UnicodeError, ValueError) as exc:
        audit["issues"] = [f"invalid bundle: {exc}"]
        return audit
    if not isinstance(manifest, dict):
        audit["issues"] = ["bundle manifest must be an object"]
        return audit

    issues = _json_structure_issues(manifest)
    version = manifest.get("schema_version", 1)
    supported_version: int | None = None
    if not isinstance(version, int) or isinstance(version, bool):
        issues.append("bundle schema version must be an integer")
    else:
        audit["manifest_version"] = version
        if version not in (1, 2):
            issues.append(f"unsupported bundle schema version: {version}")
        else:
            supported_version = version

    if not isinstance(manifest.get("status"), str):
        issues.append("bundle status must be a string")
    elif manifest["status"] != "complete":
        issues.append("bundle is not complete")
    report_generated = manifest.get("report_generated")
    if not isinstance(report_generated, bool):
        issues.append("report_generated must be a boolean")
    elif require_report and report_generated is not True:
        issues.append("report was not generated")

    required = manifest.get("required_artifacts")
    if not isinstance(required, list) or not required:
        issues.append("required_artifacts must be a non-empty list")
        required_paths: list[str] = []
    else:
        required_paths = []
        for index, relative in enumerate(required):
            if not isinstance(relative, str):
                issues.append(f"required_artifacts[{index}] must be a string")
                continue
            try:
                _portable_parts(relative)
            except ValueError as exc:
                issues.append(str(exc))
            required_paths.append(relative)
        if len(required_paths) != len(set(required_paths)):
            issues.append("required_artifacts contains duplicate paths")

    artifacts = manifest.get("artifacts", {})
    if not isinstance(artifacts, dict):
        issues.append("artifacts must be an object")
        artifact_paths: list[str] = []
        artifacts = {}
    else:
        artifact_paths = list(artifacts)
        for relative, record in artifacts.items():
            try:
                _portable_parts(relative)
            except ValueError as exc:
                issues.append(str(exc))
            if relative == _MANIFEST_NAME:
                issues.append(f"artifacts must not inventory {_MANIFEST_NAME} itself")
            if not isinstance(record, dict):
                issues.append(f"invalid artifact record: {relative}")
                continue
            byte_count = record.get("bytes")
            checksum = record.get("sha256")
            if not isinstance(byte_count, int) or isinstance(byte_count, bool) or byte_count < 0:
                issues.append(f"artifact record has invalid bytes: {relative}")
            if not isinstance(checksum, str) or _SHA256.fullmatch(checksum) is None:
                issues.append(f"artifact record has invalid sha256: {relative}")

    issues.extend(_casefold_collision_issues(required_paths, context="required_artifacts"))
    issues.extend(_casefold_collision_issues(artifact_paths, context="artifacts"))
    issues.extend(
        _casefold_collision_issues(required_paths + artifact_paths, context="manifest paths")
    )
    if version == 2:
        if not artifact_paths:
            issues.append("bundle checksum inventory is missing")
        missing_records = sorted(set(required_paths) - set(artifact_paths))
        issues.extend(f"required artifact has no checksum record: {path}" for path in missing_records)
    if issues:
        audit["issues"] = sorted(set(issues))
        return audit

    actual_files = set(files) - {_MANIFEST_NAME}
    for relative in required_paths:
        if relative not in actual_files:
            issues.append(f"required artifact is missing: {relative}")
    if version == 2 and actual_files != set(artifact_paths):
        issues.append("checksum inventory does not cover exactly the bundle files")
    if issues:
        audit["issues"] = sorted(set(issues))
        return audit

    report_present = "REPORT.html" in actual_files
    if report_generated is True and not report_present:
        issues.append("report_generated is true but REPORT.html is missing")
    elif report_generated is False and report_present:
        issues.append("report_generated is false but REPORT.html is present")

    if version == 2:
        for relative, record in artifacts.items():
            if relative not in actual_files:
                issues.append(f"artifact is missing: {relative}")
                continue
            path = contained_file(root, relative)
            try:
                size = path.stat().st_size
                checksum = sha256(path)
            except OSError as exc:
                issues.append(f"cannot inspect artifact {relative}: {exc}")
                continue
            audit["hashes_checked"] += 1
            if size != record["bytes"] or checksum != record["sha256"]:
                issues.append(f"artifact checksum mismatch: {relative}")
    if report_present and _audit_report_html:
        try:
            html_issues, html_checks_performed = _report_asset_audit(root)
            issues.extend(html_issues)
            audit["html_checks_performed"] = html_checks_performed
        except (OSError, RecursionError, TypeError, UnicodeError, ValueError) as exc:
            issues.append(f"cannot inspect REPORT.html: {exc}")
    audit["issues"] = sorted(set(issues))
    audit["valid"] = not audit["issues"]
    if audit["valid"]:
        audit["assurance"] = (
            "v2-complete-checksums" if supported_version == 2 else "legacy-existence"
        )
    return audit


def validate_files(
    root: Path,
    *,
    require_report: bool = True,
    _audit_report_html: bool = True,
) -> list[str]:
    """Compatibility wrapper returning only file-contract defects."""
    return audit_files(
        root,
        require_report=require_report,
        _audit_report_html=_audit_report_html,
    )["issues"]


def main(argv: list[str] | None = None) -> int:
    """Installed, dependency-light audit: python -m predictor.report.bundle_contract DIR."""
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path)
    parser.add_argument(
        "--no-report",
        action="store_true",
        help="allow a bundle whose manifest declares that REPORT.html was not generated",
    )
    args = parser.parse_args(argv)
    audit = audit_files(args.bundle, require_report=not args.no_report)
    print(json.dumps(audit, indent=2))
    return 0 if audit["valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
