"""Record the export environment without credentials, machine names or local paths."""
from __future__ import annotations

import platform
from datetime import datetime, timezone
from importlib.metadata import distributions


def snapshot() -> dict:
    """Installed versions at export time, distinct from the prediction's input provenance."""
    packages = sorted({(dist.metadata.get("Name", "unknown"), dist.version) for dist in distributions()})
    return {
        "schema_version": 1,
        "scope": "report export environment; not a record of external engine invocations",
        "captured_at_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "platform": {"system": platform.system(), "release": platform.release(), "machine": platform.machine()},
        "packages": [{"name": name, "version": version} for name, version in packages],
    }
