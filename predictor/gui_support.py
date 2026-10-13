"""Small, Streamlit-independent helpers shared by GUI reruns and browser sessions."""
from __future__ import annotations

import base64
import html
import io
import mimetypes
import os
import re
import threading
import zipfile
from pathlib import Path
from urllib.parse import unquote, urlsplit


# Streamlit re-executes the page for every interaction. Keeping this lock in an imported module makes it
# process-wide instead of recreating it inside each page run. ``redirect_stdout`` needs that scope.
PIPELINE_RUN_LOCK = threading.Lock()


def _linked(path: Path) -> bool:
    return path.is_symlink() or getattr(path, "is_junction", lambda: False)()


def _contained_asset(root: Path, reference: str) -> Path | None:
    """Resolve only local, non-linked report assets inside the selected folder."""
    try:
        parts = urlsplit(html.unescape(reference))
    except ValueError:
        return None
    if parts.scheme or parts.netloc or not parts.path:
        return None
    relative = Path(unquote(parts.path))
    if relative.is_absolute() or ".." in relative.parts:
        return None
    path = root / relative
    try:
        path.resolve().relative_to(root.resolve())
        if any(_linked(parent) for parent in (path, *path.parents) if parent != root and root in parent.parents):
            return None
        return path if path.is_file() else None
    except (OSError, ValueError):
        return None


def report_preview(text: str, root: Path) -> str:
    """Embed contained images and disable file links unavailable in the GUI iframe."""
    def image(match):
        path = _contained_asset(root, match.group(2))
        if path is None:
            return match.group(0)
        mime = mimetypes.guess_type(path.name)[0] or "application/octet-stream"
        encoded = base64.b64encode(path.read_bytes()).decode("ascii")
        return f'src="data:{mime};base64,{encoded}"'

    def anchor(match):
        reference = html.unescape(match.group(3))
        try:
            parts = urlsplit(reference)
        except ValueError:
            parts = None
        if reference.startswith("#") or (parts and parts.scheme in {"https", "http", "mailto"}):
            return match.group(0)
        return match.group(1) + 'aria-disabled="true" title="Open this file from the extracted full download"'

    text = re.sub(r"\bsrc\s*=\s*([\"'])(.*?)\1", image, text, flags=re.I)
    return re.sub(r"(<a\b[^>]*?)\bhref\s*=\s*([\"'])(.*?)\2", anchor, text, flags=re.I)


def context_archive(root: Path) -> tuple[bytes, list[str]]:
    """Snapshot a saved folder as a portable ZIP; skip linked or escaping descendants."""
    root = root.resolve(strict=True)
    if not root.is_dir():
        raise ValueError("The saved result folder is unavailable")
    data = io.BytesIO()
    skipped = []
    with zipfile.ZipFile(data, "w", zipfile.ZIP_DEFLATED) as archive:
        for directory, dirs, files in os.walk(root, followlinks=False):
            base = Path(directory)
            dirs.sort(key=str.casefold)
            for name in list(dirs):
                path = base / name
                if name.startswith(".") or _linked(path) or not path.resolve().is_relative_to(root):
                    dirs.remove(name)
                    skipped.append(path.relative_to(root).as_posix())
            for name in sorted(files, key=str.casefold):
                path = base / name
                relative = path.relative_to(root).as_posix()
                if _linked(path) or not path.resolve().is_relative_to(root):
                    skipped.append(relative)
                    continue
                archive.write(path, relative)
    return data.getvalue(), sorted(skipped)


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
