"""Small, Streamlit-independent helpers shared by GUI reruns and browser sessions."""
from __future__ import annotations

import re
import threading
from pathlib import Path


# Streamlit re-executes the page for every interaction. Keeping this lock in an imported module makes it
# process-wide instead of recreating it inside each page run. ``redirect_stdout`` needs that scope.
PIPELINE_RUN_LOCK = threading.Lock()


def safe_name(value: str, fallback: str = "query") -> str:
    """Return a conservative, bounded name suitable for a result directory."""
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", (value or "").strip()).strip("._")
    return (name or fallback)[:80]


def reserve_output_dir(root: Path, stem: str) -> Path:
    """Atomically reserve a new result directory, suffixing rather than overwriting."""
    root.mkdir(parents=True, exist_ok=True)
    for number in range(1, 10_000):
        name = stem if number == 1 else f"{stem}_{number}"
        candidate = root / name
        try:
            candidate.mkdir(exist_ok=False)
        except FileExistsError:
            continue
        return candidate
    raise RuntimeError(f"Could not allocate a new output folder for {stem!r}.")


def available_job_name(jobs_root: Path, stem: str) -> str:
    """Return a safe prediction name whose output directory does not yet exist."""
    base = safe_name(stem, "TF_prediction")
    for number in range(1, 10_000):
        name = base if number == 1 else f"{base}_{number}"
        if not (jobs_root / name).exists():
            return name
    raise RuntimeError(f"Could not allocate a new prediction name for {base!r}.")
