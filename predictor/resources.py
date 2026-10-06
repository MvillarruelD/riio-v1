"""Central paths for immutable data, writable caches, outputs, tools, and credentials.

The package is intentionally split into three namespaces:

* :data:`DATA_DIR` is immutable, versioned reference data shipped in the wheel.
* :data:`CACHE_DIR` is disposable machine-local state under the platform user cache.
* :data:`OUTPUT_DIR` is user-owned run output, defaulting to ``./results``.

Runtime code must never write beneath :data:`DATA_DIR`; doing so breaks a normal, non-editable
installation.  Database builders may read the shipped data, but write regenerated artifacts to the
cache (or to an explicit output path) before a maintainer deliberately vendors them.

Path overrides are locations only; none selects an algorithm:

``PREDICTOR_DATA_DIR``
    Entire packaged data tree (advanced deployments).
``PREDICTOR_REFS_DIR``
    Only the small family-reference subtree.
``PREDICTOR_CACHE_DIR``
    Writable cache root.
``PREDICTOR_OUTPUT_DIR``
    Run-output root.
``PREDICTOR_ESMFOLD_LIB`` / ``PREDICTOR_ENV_FILE``
    Optional folding adapter and credentials file.

Run ``python -m predictor.resources`` or ``tfop audit`` to inspect the resolved layout.
"""
from __future__ import annotations

import json
import os
import re
from pathlib import Path

from platformdirs import user_cache_path

PACKAGE_DIR = Path(__file__).resolve().parent
REPO = PACKAGE_DIR.parent

# Immutable package data.  Keeping these as real Paths (rather than Traversable resources) is safe because
# Python wheels are installed unpacked; the external engines also require filesystem paths.
DATA_DIR = Path(os.environ.get("PREDICTOR_DATA_DIR", PACKAGE_DIR / "data")).expanduser().resolve()
REFS_DIR = Path(os.environ.get("PREDICTOR_REFS_DIR", DATA_DIR / "refs")).expanduser().resolve()
SSN_DIR = DATA_DIR / "ssn" / "clusters"
SSN_DATABASE = DATA_DIR / "ssn" / "database"

# Writable state is never placed in the package or checkout.
CACHE_DIR = Path(
    os.environ.get("PREDICTOR_CACHE_DIR", user_cache_path("tf-operator-predictor", "tfop"))
).expanduser().resolve()
OUTPUT_DIR = Path(os.environ.get("PREDICTOR_OUTPUT_DIR", Path.cwd() / "results")).expanduser().resolve()


def cache_path(*parts: str, create: bool = False) -> Path:
    """Return a path below the writable user cache, optionally creating its directory."""
    path = CACHE_DIR.joinpath(*parts)
    if create:
        path.mkdir(parents=True, exist_ok=True)
    return path


def output_path(*parts: str, create: bool = False) -> Path:
    """Return a path below the user-owned output root, optionally creating its directory."""
    path = OUTPUT_DIR.joinpath(*parts)
    if create:
        path.mkdir(parents=True, exist_ok=True)
    return path


# --------------------------------------------------------------------------- external TOOL: ESMFold2 _lib
def esmfold_lib_dir() -> Path | None:
    """Directory that contains `_lib/` (the biohub ESMFold2 backend). Override `PREDICTOR_ESMFOLD_LIB`;
    otherwise use the optional adapter in a source checkout.  None means folding abstains cleanly."""
    env = os.environ.get("PREDICTOR_ESMFOLD_LIB")
    for cand in (Path(env).expanduser() if env else None, PACKAGE_DIR / "_engines" / "esmfold_lib", PACKAGE_DIR / "_engines" / "esmfold_lib", REPO / "env" / "tools" / "esmfold_lib"):
        if cand and (cand / "_lib").is_dir():
            return cand
    return None


# --------------------------------------------------------------------------- CREDENTIALS: the .env tokens
def env_file() -> Path | None:
    """The optional credentials file: explicit override, then a source-checkout ``.env``."""
    env = os.environ.get("PREDICTOR_ENV_FILE")
    for cand in (Path(env).expanduser() if env else None, REPO / ".env"):
        if cand and cand.is_file():
            return cand
    return None


_TOKEN_KEYS = ("BIOHUB_TOKEN", "NCBI_TOKEN", "NCBI_API_KEY", "NCBI_EMAIL")
#: Engine LOCATIONS, read from the same `.env`. Not credentials -- paths to tools that live outside
#: PATH. They belong here because the audit and the code that runs the tool must resolve them
#: identically: `tfop setup` reporting a tool MISSING while the pipeline would happily find it is a
#: false negative from the one report whose whole job is to say what is missing.
_PATH_KEYS = ("DATASETS_EXE", "BITACORA_DIR", "BITACORA_BIN", "BITACORA_WSL_DISTRO")


def load_tokens(*, override: bool = False) -> dict:
    """Populate os.environ with the API tokens from `env_file()` (so BOTH the ESMFold2 `_lib` token loader
    and `genome_resolver._ncbi_creds` see them regardless of where the .env lives). Existing env vars win
    unless `override`. Returns {key: present?} for the four token keys. Safe to call repeatedly + offline."""
    ef = env_file()
    if ef is not None:
        try:
            txt = ef.read_text(errors="ignore")
            for k in _TOKEN_KEYS + _PATH_KEYS:
                if override or not os.environ.get(k):
                    m = re.search(rf"^{k}\s*=\s*(\S+)", txt, re.MULTILINE)
                    if m and m.group(1):
                        os.environ[k] = m.group(1)
        except Exception:
            pass
    return {k: bool(os.environ.get(k)) for k in _TOKEN_KEYS}


def have_biohub_token() -> bool:
    load_tokens()
    tok = os.environ.get("BIOHUB_TOKEN", "").strip()
    return bool(tok) and tok != "your-token-here"


def database_info() -> dict:
    """Describe the immutable reference database shipped with this installation.

    This is deliberately manifest-based: it reports the released database bytes, not transient
    ``__pycache__`` files or other local state that may happen to sit below ``DATA_DIR``.
    """
    manifest_path = DATA_DIR / "MANIFEST.json"
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        files = manifest.get("files") or []
        total_bytes = sum(int(row.get("bytes") or 0) for row in files)
        return {
            "present": True,
            "release": str(manifest.get("database_release") or "unknown"),
            "file_count": int(manifest.get("file_count") or len(files)),
            "bytes": total_bytes,
            "mib": round(total_bytes / (1024 * 1024), 1),
            "manifest": str(manifest_path),
        }
    except (OSError, ValueError, TypeError):
        return {"present": False, "release": "unknown", "file_count": 0, "bytes": 0,
                "mib": 0.0, "manifest": str(manifest_path)}


# --------------------------------------------------------------------------- audit (WS7 verification)
def audit() -> dict:
    """Resolve every resource and report whether a normal installed package is operational."""
    tokens = load_tokens()
    lib = esmfold_lib_dir()
    ef = env_file()
    refs_ok = (REFS_DIR / "all_tfs.faa").exists() and (REFS_DIR / "curated_tfs.faa").exists()
    ssn_ok = SSN_DIR.is_dir() and SSN_DATABASE.is_dir()
    rows = [
        {"resource": "family references", "kind": "data", "path": str(REFS_DIR),
         "present": refs_ok, "source": "packaged" if refs_ok else "MISSING"},
        {"resource": "SSN cluster FASTAs", "kind": "data", "path": str(SSN_DIR),
         "present": SSN_DIR.is_dir(), "source": "packaged" if SSN_DIR.is_dir() else "MISSING"},
        {"resource": "SSN indexed database", "kind": "data", "path": str(SSN_DATABASE),
         "present": SSN_DATABASE.is_dir(), "source": "packaged" if SSN_DATABASE.is_dir() else "MISSING"},
        {"resource": "writable cache", "kind": "cache", "path": str(CACHE_DIR),
         "present": True, "source": "platform-user-cache"},
        {"resource": "run outputs", "kind": "output", "path": str(OUTPUT_DIR),
         "present": True, "source": "working-directory/override"},
        {"resource": "ESMFold2 _lib (folding tool)", "kind": "tool",
         "path": (str(lib) if lib else None), "present": lib is not None,
         "source": ("override/checkout" if lib else "absent")},
        {"resource": "API tokens (.env)", "kind": "credentials", "path": (str(ef) if ef else None),
         "present": ef is not None,
         "source": ("override/checkout" if ef else "absent")},
    ]
    return {"rows": rows, "tokens": tokens, "database": database_info(),
            "production_data_agnostic": refs_ok and ssn_ok,
            "data_dir": str(DATA_DIR), "cache_dir": str(CACHE_DIR), "output_dir": str(OUTPUT_DIR)}


def _print_audit() -> None:
    a = audit()
    db = a["database"]
    if db["present"] and a["production_data_agnostic"]:
        print(f"REFERENCE DATABASE: READY -- release {db['release']}, {db['file_count']} files, "
              f"{db['mib']:.1f} MiB installed")
        print("  Bundled with tf-operator-predictor; users do not download or build it separately.")
        print("  Query genomes/promoter context are fetched from NCBI as needed and cached locally.")
    else:
        print("REFERENCE DATABASE: INCOMPLETE -- reinstall the package, then run `tfop selftest`.")
    print()
    print(f"{'resource':<34} {'kind':<12} {'present':<8} {'source':<27} path")
    print("-" * 112)
    for r in a["rows"]:
        print(f"{r['resource']:<34} {r['kind']:<12} {str(r['present']):<8} {r['source']:<27} {r['path']}")
    print("-" * 112)
    print("tokens: " + ", ".join(f"{k}={'set' if v else '-'}" for k, v in a["tokens"].items()))
    print(f"PRODUCTION DATA-AGNOSTIC: {a['production_data_agnostic']}  "
          f"(immutable references packaged; caches and outputs external)")
    print("Optional large assets are NOT bundled: HomoDB (~13.7 GB) and MetalNet2 weights (~6 GB).")
    print("They add structural evidence only; a standard prediction does not require either one.")


def _demo() -> None:
    a = audit()
    # the two DATA resources must be vendored in-repo (never a 5.8 fallback)
    refs = next(r for r in a["rows"] if r["resource"].startswith("family references"))
    ssn = next(r for r in a["rows"] if r["resource"].startswith("SSN"))
    assert refs["present"] and refs["source"] == "packaged", refs
    assert ssn["present"] and ssn["source"] == "packaged", ssn
    assert a["production_data_agnostic"] is True, "packaged reference data are incomplete"
    # load_tokens is idempotent + offline-safe
    t1 = load_tokens(); t2 = load_tokens()
    assert set(t1) == set(_TOKEN_KEYS) and t1 == t2
    print("OK: packaged data resolve; writable caches and outputs are outside the package.")
    _print_audit()


if __name__ == "__main__":
    _demo()
