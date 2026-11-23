"""Small, Streamlit-independent helpers shared by GUI reruns and browser sessions."""
from __future__ import annotations

import base64
import csv
import html
import io
import json
import mimetypes
import os
import re
import shlex
import sys
import threading
import zipfile
from datetime import datetime
from functools import lru_cache
from pathlib import Path
from urllib.parse import unquote, urlsplit


# Streamlit re-executes the page for every interaction. Keeping this lock in an imported module makes it
# process-wide instead of recreating it inside each page run. ``redirect_stdout`` needs that scope.
PIPELINE_RUN_LOCK = threading.Lock()

#: The four inducer classes the GUI draws. ``unresolved`` is a result, never a silent fallback to organic.
CLASSES = ("metal", "redox", "organic", "unresolved")


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


# --------------------------------------------------------------------------- headline calls
def inducer_class(top: str | None, stored: str | None = None, *, has_stored: bool = False) -> str:
    """One of :data:`CLASSES`. A stored ``sensor_class`` wins; an empty inducer is ``unresolved``."""
    if not top:
        return "unresolved"
    if has_stored:
        return stored if stored in CLASSES[:3] else "unresolved"
    try:
        from predictor.effector import inducer_vocab
        found = inducer_vocab.inducer_class_of(top)
    except Exception:  # the legacy label regression is a convenience; its absence is not an error
        found = None
    return found if found in CLASSES[:3] else "unresolved"


def _float(value) -> float | None:
    try:
        return None if value in (None, "") else float(value)
    except (TypeError, ValueError):
        return None


def bundle_headline(bundle: Path) -> dict:
    """The headline calls of one finished bundle, read from ``dossier.json`` the way the report reads them."""
    bundle = Path(bundle)
    dossier = json.loads((bundle / "dossier.json").read_text(encoding="utf-8"))
    ind = dossier.get("inducers") or {}
    # A saved final null is a result: do not fall back to a stale effector field once the call was made.
    top = ind.get("top") if ("top" in ind or ind.get("finalized")) else (dossier.get("effector") or None)
    label = top or ""
    try:
        from predictor.report.figure_data import inducer_display
        label = inducer_display(ind) or label
    except Exception:  # display is presentation only
        pass
    primary = (dossier.get("operators") or {}).get("primary") or {}
    provenance = primary.get("provenance") or {}
    logo = dossier.get("genomic_logo")
    return {
        "family": dossier.get("family"),
        "family_support": dossier.get("family_support"),
        "genome": dossier.get("genome_accession"),
        "tf_locus": dossier.get("tf_locus"),
        "inducer": top,
        "inducer_label": label or None,
        "inducer_class": inducer_class(top, ind.get("sensor_class"), has_stored="sensor_class" in ind),
        "agreement": _float(ind.get("agreement")),
        "agreement_exact": bool(ind.get("finalized")),
        "ssn_cluster": dossier.get("ssn_cluster"),
        "operator": primary.get("sequence"),
        "operator_source": primary.get("source"),
        "operator_locality": provenance.get("locality"),
        "n_sites": len((dossier.get("rescan") or {}).get("hits") or []),
        "n_regulon_operons": len([o for o in (dossier.get("regulon") or [])
                                  if isinstance(o, dict) and o.get("first_gene")]),
        "motif": logo.get("consensus") if isinstance(logo, dict) else None,
    }


def result_dict(bundle: Path) -> dict:
    """Rebuild the ``run_novel`` result keys the GUI reads from a finished bundle on disk."""
    bundle = Path(bundle)
    dossier = json.loads((bundle / "dossier.json").read_text(encoding="utf-8"))
    head = bundle_headline(bundle)
    return {**dossier, "tf_id": dossier.get("tf_id") or bundle.name, "bundle": str(bundle),
            "ok": True, "required_stages_complete": True, "genome": head["genome"], "inducer": head["inducer"],
            "operators_primary": (dossier.get("operators") or {}).get("primary")}


# --------------------------------------------------------------------------- saved runs
def read_summary(out: Path) -> list[dict]:
    """Rows of ``summary.tsv`` with the inducer class normalised. Empty when the run has no summary yet."""
    path = Path(out) / "summary.tsv"
    if not path.is_file():
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        top = (row.get("inducer_top") or "").strip() or None
        explicit = row.get("inducer_class_provenance") == "explicit"
        row["inducer_top"] = top or ""
        row["class"] = (inducer_class(top, row.get("inducer_class") or None, has_stored=True)
                        if explicit else inducer_class(top))
        row["agreement_value"] = _float(row.get("agreement"))
    return rows


def read_manifest(out: Path) -> list[dict]:
    with (Path(out) / "manifest.csv").open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def is_done(status) -> bool:
    """`ok` from this invocation, or `done (resumed)` from an earlier one."""
    return str(status or "").startswith(("ok", "done"))


def is_failed(status) -> bool:
    return str(status or "").startswith(("failed", "timed out"))


def job_dir(out: Path, name: str) -> Path | None:
    """`jobs/<name>` of a run, or None: names come from a table on disk and must not leave the run."""
    if not name or name in (".", "..") or re.search(r"[/\\:]", name):
        return None
    jobs = Path(out) / "jobs"
    path = jobs / name
    try:
        return path if path.resolve().parent == jobs.resolve() and not _linked(path) else None
    except OSError:
        return None


@lru_cache(maxsize=4096)
def _headline_cached(bundle: str, stamp: int) -> dict | None:
    from predictor.genome_run import bundle_complete
    path = Path(bundle)
    try:
        return bundle_headline(path) if bundle_complete(path) else None
    except (OSError, ValueError, AttributeError):
        return None


def run_rows(out: Path, manifest: list[dict]) -> list[dict]:
    """Every candidate with its headline, as the run report reads it: the bundle wins over the table.

    summary.tsv is written when an invocation ends, so it lags an interrupted or still running one, and
    it stores the raw inducer, not the label the report shows (a shortlist must list its members).
    """
    summary = {row.get("run_name"): row for row in read_summary(out) if row.get("run_name")}
    listed = {row.get("run_name"): row for row in manifest if row.get("run_name")}
    rows = []
    for name in dict.fromkeys([*listed, *summary]):
        item = listed.get(name, {})
        row = {"run_name": name, "gene": item.get("gene", ""), "locus_tag": item.get("locus_tag", ""),
               "product": item.get("product", ""), "discovery_family": item.get("family", ""),
               "status": "not run", "class": "unresolved", "inducer_top": "", "agreement_value": None,
               **summary.get(name, {})}
        row["inducer_label"] = row["inducer_top"]
        bundle = job_dir(out, name)
        if bundle is not None and not is_failed(row["status"]):
            try:
                head = _headline_cached(str(bundle), (bundle / "bundle_manifest.json").stat().st_mtime_ns)
            except OSError:
                head = None
            if head:
                row.update({"family": head["family"], "inducer_top": head["inducer"] or "",
                            "inducer_label": head["inducer_label"] or "", "class": head["inducer_class"],
                            "agreement_value": head["agreement"], "primary_operator": head["operator"] or "",
                            "n_regulon_operons": head["n_regulon_operons"]})
                if not is_done(row["status"]):
                    row["status"] = "done (resumed)"
        rows.append(row)
    return rows


def _json_object(path: Path) -> dict:
    try:
        value = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return value if isinstance(value, dict) else {}


def run_kind(out: Path) -> str | None:
    """"batch" or "genome" as the run recorded it; None for runs older than run_input.json."""
    mode = _json_object(Path(out) / "run_input.json").get("mode")
    return None if not mode else "batch" if mode == "protein_list" else "genome"


def recorded_policy(out: Path) -> dict:
    """The options the first build fixed (the runner refuses any others in the same folder)."""
    return _json_object(Path(out) / "prediction_input.json").get("policy") or {}


def run_command(out: Path, *, windows: bool | None = None, **options) -> str | None:
    """The `tfop` command that continues this folder, or None when it cannot be rebuilt exactly.

    A wrong command is worse than none: the runner fingerprints its inputs and would refuse it.
    """
    out = Path(out)
    info = _json_object(out / "run_input.json")
    mode = info.get("mode")
    policy = info.get("discovery_policy") or {}
    if policy.get("census_only"):
        return None                         # a Pfam-only census has no `tfop genome` equivalent
    uploads = out / "uploads"
    roles = {item.get("role") for item in info.get("files") or [] if isinstance(item, dict)}

    def staged(*suffixes):
        if not uploads.is_dir():
            return None
        return next((str(path) for path in sorted(uploads.iterdir())
                     if path.stem == "input" and path.suffix.lower() in suffixes), None)

    common = dict(options, out=str(out), organism=info.get("organism"), windows=windows)
    if mode == "protein_list":
        source = staged(".faa", ".fa", ".fasta", ".tsv")
        return cli_command("batch", source=source, **common) if source and roles <= {"proteins"} else None
    common["census_fallback"] = policy.get("census_fallback")
    if mode == "genome_accession":
        return cli_command("genome", source=info["accession"], **common) if info.get("accession") else None
    inputs = {"local_genome": ("fasta", (".fna", ".fa", ".fasta")), "genbank": ("genbank", (".gb", ".gbk", ".gbff")),
              "proteome": ("proteome", (".faa", ".fa", ".fasta"))}
    if mode not in inputs:
        return None
    key, suffixes = inputs[mode]
    source = staged(*suffixes)
    gff = staged(".gff", ".gff3") if "gff" in roles else None
    if not source or ("gff" in roles and not gff):
        return None
    return cli_command("genome", gff=gff, **{key: source}, **common)


def list_saved_runs(root: Path, limit: int = 40, kind: str | None = None) -> list[dict]:
    """Multi-protein runs under ``root`` (newest first), so a finished run survives a browser refresh."""
    root = Path(root)
    found = []
    if not root.is_dir():
        return found
    for child in root.iterdir():
        try:
            if not child.is_dir() or _linked(child) or not (child / "manifest.csv").is_file():
                continue
            if kind and run_kind(child) not in (None, kind):
                continue
            with (child / "manifest.csv").open(encoding="utf-8", newline="") as handle:
                total = sum(1 for _ in csv.DictReader(handle))
            summary = read_summary(child)
            done = {row.get("run_name") for row in summary if is_done(row.get("status"))}
            failed = {row.get("run_name") for row in summary if is_failed(row.get("status"))}
            done |= {path.parent.name for path in child.glob("jobs/*/bundle_manifest.json")
                     if not path.parent.name.startswith(".") and path.parent.name not in failed}
            found.append({"name": child.name, "path": child, "total": total, "done": len(done),
                          "modified": datetime.fromtimestamp(child.stat().st_mtime)})
        except (OSError, csv.Error, UnicodeDecodeError):
            continue
    found.sort(key=lambda r: r["modified"], reverse=True)
    return found[:limit]


def list_saved_bundles(jobs_root: Path, limit: int = 40) -> list[dict]:
    """Finished single-TF bundles under ``jobs_root`` (newest first)."""
    jobs_root = Path(jobs_root)
    found = []
    if not jobs_root.is_dir():
        return found
    for child in jobs_root.iterdir():
        try:
            if child.is_dir() and not _linked(child) and (child / "dossier.json").is_file():
                found.append({"name": child.name, "path": child,
                              "modified": datetime.fromtimestamp(child.stat().st_mtime)})
        except OSError:
            continue
    found.sort(key=lambda r: r["modified"], reverse=True)
    return found[:limit]


# --------------------------------------------------------------------------- paths and terminal commands
def display_paths(path: Path) -> list[tuple[str, str]]:
    """Every spelling of ``path`` a person on this machine may need (WSL also has a Windows-side name)."""
    path = Path(path)
    shown = [("Folder", str(path))]
    distro = os.environ.get("WSL_DISTRO_NAME")
    text = path.as_posix()
    if distro and sys.platform.startswith("linux") and text.startswith("/"):
        parts = path.parts  # ('/', 'mnt', 'c', ...)
        if text.startswith("/mnt/") and len(parts) > 2 and len(parts[2]) == 1:
            shown.append(("In Windows", parts[2].upper() + ":\\" + "\\".join(parts[3:])))
        else:
            shown.append(("In Windows", "\\\\wsl.localhost\\" + distro + text.replace("/", "\\")))
    return shown


def _quote_windows(word: str) -> str:
    """One argument for PowerShell (5.1 and 7); also valid in cmd unless it holds $, a backtick or a quote."""
    if re.fullmatch(r"[A-Za-z0-9_+=:.,/\\-]+", word):
        return word
    if any(char in word for char in '$`"\u201c\u201d'):
        return "'" + word.replace("'", "''") + "'"     # literal string: nothing is expanded
    return '"' + word + '"'


def quote_command(parts: list[str], *, windows: bool | None = None) -> str:
    """One command line quoted for PowerShell on Windows and for POSIX shells elsewhere."""
    windows = (os.name == "nt") if windows is None else windows
    words = [str(p) for p in parts]
    return " ".join(map(_quote_windows, words)) if windows else shlex.join(words)


def cli_command(kind: str, *, windows: bool | None = None, **o) -> str:
    """The ``tfop`` command that reproduces a GUI run. Flags mirror ``predictor.run_cli``; none is invented."""
    parts = ["tfop"]

    def flag(on, name):
        if on:
            parts.append(name)

    if kind == "predict":
        parts += ["predict", o["source"]]
        if o.get("name"):
            parts += ["--name", o["name"]]
        if o.get("out"):
            parts += ["--out", o["out"]]
        for value, name in ((o.get("family"), "--family"), (o.get("organism"), "--organism"),
                            (o.get("effector"), "--effector")):
            if value:
                parts += [name, value]
        flag(o.get("fold"), "--fold")
        flag(o.get("offline"), "--offline")
        flag(not o.get("metalnet", True), "--no-metalnet")
        flag(not o.get("af3", True), "--no-af3")
    elif kind in {"genome", "batch", "scan"}:
        parts.append(kind)
        if o.get("source"):
            parts.append(o["source"])
        parts += ["--out", o["out"]]
        for key, name in (("fasta", "--fasta"), ("genbank", "--genbank"), ("proteome", "--proteome"),
                          ("gff", "--gff"), ("organism", "--organism")):
            if o.get(key):
                parts += [name, o[key]]
        if kind != "scan":
            if int(o.get("jobs") or 1) > 1:
                parts += ["--jobs", str(int(o["jobs"]))]
            if o.get("select"):
                parts += ["--select", ",".join(o["select"])]
            flag(o.get("fold"), "--fold")
            flag(o.get("offline"), "--offline")
            flag(not o.get("metalnet", True), "--no-metalnet")
            flag(not o.get("af3", True), "--no-af3")
            flag(o.get("discover_only"), "--discover-only")
            flag(kind == "genome" and o.get("census_fallback"), "--census-fallback")
    else:
        raise ValueError(f"unknown command kind {kind!r}")
    return quote_command(parts, windows=windows)


def pasted_protein_input(text: str) -> tuple[str, str]:
    """Normalise a pasted protein list to ``(suffix, content)`` for ``predict_many``.

    FASTA stays FASTA; anything else is read as accessions (one per line, ``#`` comments allowed) and
    becomes the accession TSV the batch runner already accepts.
    """
    text = (text or "").strip()
    if not text:
        raise ValueError("Paste FASTA records or one protein accession per line.")
    if text.startswith(">"):
        return ".faa", text + "\n"
    accessions: list[str] = []
    for line in text.splitlines():
        accessions.extend(line.split("#", 1)[0].replace(",", " ").split())
    if not accessions:
        raise ValueError("No accession was found in the pasted text.")
    for accession in accessions:
        if not re.fullmatch(r"[A-Za-z0-9_.\-]+", accession):
            raise ValueError(f"Not a protein accession: {accession!r}. Paste FASTA or one accession per line.")
    return ".tsv", "accession\n" + "\n".join(accessions) + "\n"
