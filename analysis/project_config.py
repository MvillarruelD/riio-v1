"""Portable environment discovery shared by the production analysis drivers."""
from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent.parent
ANALYSIS_ROOT = PROJECT_ROOT / "analysis"


def predictor_root() -> Path:
    """Resolve the predictor checkout/install, preferring an explicit ``TFOP_REPO``."""
    explicit = os.environ.get("TFOP_REPO", "").strip()
    if explicit:
        root = Path(explicit).expanduser().resolve()
        if not (root / "predictor" / "__init__.py").is_file():
            raise RuntimeError(f"TFOP_REPO does not contain the predictor package: {root}")
        return root

    # Conventional checkout locations keep a source workspace usable before the editable install.
    # They are path-relative or home-relative, never tied to one account name.
    candidates = (
        PROJECT_ROOT,
        PROJECT_ROOT,
        PROJECT_ROOT.parent / "tf-operator-predictor",
        Path.home() / "dev" / "tf-operator-predictor",
    )
    for root in candidates:
        if (root / "predictor" / "__init__.py").is_file():
            return root.resolve()

    spec = importlib.util.find_spec("predictor")
    if spec is None or spec.origin is None:
        raise RuntimeError(
            "tf-operator-predictor is not importable; install it or set TFOP_REPO to its checkout"
        )
    return Path(spec.origin).resolve().parents[1]


PREDICTOR_ROOT = predictor_root()
PYTHON = os.environ.get("TFOP_PY", sys.executable)
