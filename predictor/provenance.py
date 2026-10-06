"""Deterministic provenance for completed prediction bundles.

The package version alone cannot identify the code that produced a research result: editable
installs and unreleased fixes can share a version.  ``runtime_provenance`` therefore fingerprints
the installed Python sources and the packaged data manifest.  The value is independent of the
checkout path, timestamps, bytecode and generated caches, so it is suitable for comparing bundles
across machines.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from predictor import __version__


PACKAGE_ROOT = Path(__file__).resolve().parent
DATA_MANIFEST = PACKAGE_ROOT / "data" / "MANIFEST.json"


def _normalised(data: bytes) -> bytes:
    """Content with CRLF folded to LF.

    The repository stores text as LF (`.gitattributes: eol=lf`), but a Windows checkout can hold
    CRLF copies of the same commit. Hashing raw bytes made one commit fingerprint two ways
    depending on the checkout -- measured on the 2026-09-17 runs, where a worktree of the producing
    commit only matched after three CRLF files were copied in. The fingerprint is meant to identify
    code, not a checkout's line-ending state.
    """
    return data.replace(b"\r\n", b"\n")


def _digest_files(paths: list[Path], *, root: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted(paths, key=lambda item: item.relative_to(root).as_posix()):
        relative = path.relative_to(root).as_posix()
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        digest.update(_normalised(path.read_bytes()))
        digest.update(b"\0")
    return digest.hexdigest()


def git_provenance() -> dict[str, str | bool | None]:
    """The commit the package was loaded from, when it runs from a git checkout.

    Recorded beside the content fingerprint, never inside it: the fingerprint must stay comparable
    across machines and installs, while the commit is what a reader needs to check the code out.
    `git_dirty` covers only the package directory, because uncommitted notes elsewhere in the
    checkout do not change what ran. A wheel install, or a machine without git, reports
    "unavailable" rather than guessing.
    """
    import subprocess

    repo = PACKAGE_ROOT.parent
    if not (repo / ".git").exists():
        return {"git_commit": "unavailable", "git_dirty": None}
    try:
        commit = subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"], capture_output=True,
                                text=True, timeout=10, check=True).stdout.strip()
        status = subprocess.run(["git", "-C", str(repo), "status", "--porcelain", "--", PACKAGE_ROOT.name],
                                capture_output=True, text=True, timeout=10, check=True).stdout
    except (OSError, subprocess.SubprocessError):
        return {"git_commit": "unavailable", "git_dirty": None}
    return {"git_commit": commit, "git_dirty": bool(status.strip())}


def runtime_provenance() -> dict[str, str]:
    """Return the immutable identifiers needed to reproduce or compare a bundle."""
    source_files = [path for path in PACKAGE_ROOT.rglob("*.py") if "__pycache__" not in path.parts]
    source_sha256 = _digest_files(source_files, root=PACKAGE_ROOT)

    data_release = "unknown"
    data_manifest_sha256 = "unavailable"
    if DATA_MANIFEST.is_file():
        manifest_bytes = _normalised(DATA_MANIFEST.read_bytes())
        data_manifest_sha256 = hashlib.sha256(manifest_bytes).hexdigest()
        try:
            data_release = str(json.loads(manifest_bytes).get("database_release", "unknown"))
        except (TypeError, ValueError):
            data_release = "invalid"

    combined = hashlib.sha256()
    for value in (__version__, source_sha256, data_release, data_manifest_sha256):
        combined.update(value.encode("utf-8"))
        combined.update(b"\0")
    return {
        "package_version": __version__,
        "source_sha256": source_sha256,
        "data_release": data_release,
        "data_manifest_sha256": data_manifest_sha256,
        "runtime_sha256": combined.hexdigest(),
        **git_provenance(),
    }


if __name__ == "__main__":
    print(json.dumps(runtime_provenance(), indent=2, sort_keys=True))
