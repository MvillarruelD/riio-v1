"""
esmfold2_backend.py — Thin wrapper around the biohub ESMFold2 SDK.

Loads the API token from a project-local ``.env`` file (no env-var dance) and
exposes a single ``fold()`` entry point that takes a fully built
``StructurePredictionInput`` and returns an ``ESMFold2Result`` dataclass with
the artefacts our adapter needs (CIF text, pLDDT, PAE, pTM, ipTM, token chain
and residue ids).

Defaults are intentionally LOW (`num_loops=1`, `num_sampling_steps=16`) until
we have observed biohub's rate-limiting behaviour; callers can override per
job.
"""

from __future__ import annotations

import logging
import os
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

log = logging.getLogger(__name__)

# arsr_merr_families/scripts/_lib/esmfold2_backend.py → parents[2] = project root
PROJECT_ROOT = Path(__file__).resolve().parents[2]
DOTENV_PATH = PROJECT_ROOT / ".env"


# ── .env loading ─────────────────────────────────────────────────────────────

def _load_dotenv() -> None:
    """Populate os.environ from the project-local .env file.

    Uses python-dotenv when available; falls back to a tiny parser so the
    module can be imported even before the optional dep is installed.
    """
    if not DOTENV_PATH.exists():
        return
    try:
        from dotenv import load_dotenv as _load
        _load(DOTENV_PATH)
        return
    except ImportError:
        pass
    # Minimal fallback parser
    with open(DOTENV_PATH, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, _, val = line.partition("=")
            key = key.strip()
            val = val.strip().strip('"').strip("'")
            if key and key not in os.environ:
                os.environ[key] = val


def _load_token() -> str:
    """Return the biohub API token, or raise with a clear remediation message."""
    _load_dotenv()
    tok = os.environ.get("BIOHUB_TOKEN", "").strip()
    if not tok or tok == "your-token-here":
        raise RuntimeError(
            f"BIOHUB_TOKEN not found. Put it in {DOTENV_PATH} as:\n"
            f"  BIOHUB_TOKEN=<your-token>\n"
            f"(copy .env.example to .env and edit)."
        )
    return tok


# ── Result container ─────────────────────────────────────────────────────────

@dataclass
class ESMFold2Result:
    """Normalised view of one ESMFold2 fold call.

    Fields kept on the lean side — the full SDK response is preserved in
    ``raw_response`` for debugging without polluting the main contract.
    """
    cif_text: str
    plddt: Optional[np.ndarray]                    # per-residue 0–100
    pae: Optional[np.ndarray]                      # (n_tokens × n_tokens)
    ptm: Optional[float]
    iptm: Optional[float]
    token_chain_ids: List[str] = field(default_factory=list)
    token_res_ids: List[int] = field(default_factory=list)
    model_name: str = ""
    num_loops: int = 0
    num_sampling_steps: int = 0
    raw_response: Dict[str, Any] = field(default_factory=dict)


# ── Public fold() entry point ────────────────────────────────────────────────

def fold(
    structure_input,
    *,
    model: str = "esmfold2-fast-2026-05",
    url: str = "https://biohub.ai",
    num_loops: int = 1,
    num_sampling_steps: int = 16,
    include_pae: bool = True,
    timeout_s: int = 600,
    retries: int = 2,
) -> ESMFold2Result:
    """Call biohub's ESMFold2 and return a normalised result.

    Parameters
    ----------
    structure_input
        An ``esm.utils.structure.input_builder.StructurePredictionInput``.
    model, url, num_loops, num_sampling_steps, include_pae
        Forwarded to the SDK's FoldingConfig / client.
    timeout_s, retries
        Network-level resilience. One retry on transient 5xx; never on 4xx.
    """
    # Imports deferred so the module is importable without the SDK installed
    # (e.g. during static analysis or in CI that only needs the adapter).
    from esm.sdk import esmfold2_client
    from esm.sdk.api import FoldingConfig

    token = _load_token()
    client = esmfold2_client(model=model, url=url, token=token)
    config = FoldingConfig(
        num_loops=num_loops,
        num_sampling_steps=num_sampling_steps,
        include_pae=include_pae,
    )

    attempt = 0
    last_exc: Optional[Exception] = None
    while attempt <= retries:
        attempt += 1
        try:
            t0 = time.time()
            result = client.fold_all_atom(structure_input, config=config)

            # Validate: the SDK does NOT raise on 429 rate-limit responses — it
            # returns an empty result object with no structure/pLDDT.  Detect
            # this silent failure before logging "fold OK".
            _has_structure = (
                getattr(result, "complex", None) is not None
                or getattr(result, "structure", None) is not None
            )
            if not _has_structure:
                raise RuntimeError(
                    "ESMFold2 returned an empty result (no structure object); "
                    "possible 429 rate-limit or API error. "
                    "Raw repr: " + repr(result)[:300]
                )

            log.info(
                "ESMFold2 fold OK in %.1fs (attempt %d/%d)",
                time.time() - t0, attempt, retries + 1,
            )
            return _normalise(
                result,
                model_name=model,
                num_loops=num_loops,
                num_sampling_steps=num_sampling_steps,
            )
        except Exception as exc:  # noqa: BLE001 — SDK raises bare exceptions
            last_exc = exc
            msg = str(exc).lower()
            # 4xx-style errors + rate-limit: do not retry (won't help)
            non_retry_markers = (
                "400", "401", "403", "404",
                "unauthor", "forbid", "invalid",
                "429", "rate limit", "credit limit",
            )
            if any(m in msg for m in non_retry_markers):
                log.error("ESMFold2 non-retryable error: %s", exc)
                raise
            log.warning(
                "ESMFold2 attempt %d/%d failed: %s",
                attempt, retries + 1, exc,
            )
            if attempt > retries:
                break
            time.sleep(min(2 ** attempt, 30))
    raise RuntimeError(f"ESMFold2 failed after {retries + 1} attempts: {last_exc}")


# ── Normalisation ────────────────────────────────────────────────────────────

def _to_np(x) -> Optional[np.ndarray]:
    if x is None:
        return None
    try:
        arr = np.asarray(x, dtype=float)
        return arr
    except Exception:
        return None


def _first_not_none(*args):
    """Return the first argument that is not None.

    Safe for PyTorch tensors — avoids the ``bool(tensor)`` ambiguity that
    ``a or b`` triggers when ``a`` is a multi-element tensor.
    """
    for x in args:
        if x is not None:
            return x
    return None


def _normalise(
    raw_result,
    *,
    model_name: str,
    num_loops: int,
    num_sampling_steps: int,
) -> ESMFold2Result:
    """Best-effort extraction of CIF + confidences from the SDK result object.

    The SDK's public surface is still settling. We probe several plausible
    attribute names so the adapter doesn't break on minor API churn; anything
    we cannot find lands as ``None`` and the adapter logs it.
    """
    # mmCIF text — probe two attribute names; use _first_not_none so we never
    # invoke bool() on a tensor-valued attribute.
    cif_text = ""
    complex_obj = _first_not_none(
        getattr(raw_result, "complex", None),
        getattr(raw_result, "structure", None),
    )
    if complex_obj is not None:
        for attr in ("to_mmcif", "to_cif", "as_mmcif", "as_cif"):
            fn = getattr(complex_obj, attr, None)
            if callable(fn):
                try:
                    cif_text = fn()
                    break
                except Exception as exc:
                    log.debug("complex.%s() failed: %s", attr, exc)

    # Per-residue pLDDT — use _first_not_none; raw values may be tensors
    plddt = _first_not_none(
        getattr(raw_result, "plddt", None),
        getattr(raw_result, "per_residue_plddt", None),
        None if complex_obj is None else getattr(complex_obj, "plddt", None),
    )
    plddt_arr = _to_np(plddt)

    # PAE
    pae = _first_not_none(
        getattr(raw_result, "pae", None),
        getattr(raw_result, "predicted_aligned_error", None),
        None if complex_obj is None else getattr(complex_obj, "pae", None),
    )
    pae_arr = _to_np(pae)
    if pae_arr is not None and pae_arr.ndim != 2:
        log.warning("PAE has unexpected ndim %d; discarding", pae_arr.ndim)
        pae_arr = None

    # Global scores — convert carefully; tensors need item() before float()
    def _to_float(x) -> Optional[float]:
        if x is None:
            return None
        try:
            # 0-dim tensor or numpy scalar
            return float(x.item() if hasattr(x, "item") else x)
        except (TypeError, ValueError, RuntimeError):
            return None

    ptm = _to_float(getattr(raw_result, "ptm", None))
    iptm = _to_float(getattr(raw_result, "iptm", None))

    # Token chain / residue ids
    _tci_raw = _first_not_none(
        getattr(raw_result, "token_chain_ids", None),
        None if complex_obj is None else getattr(complex_obj, "token_chain_ids", None),
    )
    token_chain_ids = list(_tci_raw) if _tci_raw is not None else []

    _tri_raw = _first_not_none(
        getattr(raw_result, "token_res_ids", None),
        None if complex_obj is None else getattr(complex_obj, "token_res_ids", None),
    )
    token_res_ids = [
        int(x.item() if hasattr(x, "item") else x)
        for x in (_tri_raw if _tri_raw is not None else [])
    ]

    return ESMFold2Result(
        cif_text=cif_text,
        plddt=plddt_arr,
        pae=pae_arr,
        ptm=ptm,
        iptm=iptm,
        token_chain_ids=token_chain_ids,
        token_res_ids=token_res_ids,
        model_name=model_name,
        num_loops=num_loops,
        num_sampling_steps=num_sampling_steps,
        raw_response={"_sdk_object_repr": repr(raw_result)[:500]},
    )
