"""Portable canonical run configuration for the reproduction companion."""
from pathlib import Path
import os
try:
    import tomllib
except ImportError:
    import tomli as tomllib

def load_canonical():
    project = Path(__file__).resolve().parent.parent
    cfg = tomllib.loads((project / "analysis/canonical_runs.toml").read_text(encoding="utf-8"))
    base = Path(os.environ.get("RIIO_RUNS_DIR", project / "runs")).expanduser().resolve()
    for run in cfg["runs"].values():
        root = Path(run["root"])
        run["root"] = str(root if root.is_absolute() else base / root)
    return cfg
