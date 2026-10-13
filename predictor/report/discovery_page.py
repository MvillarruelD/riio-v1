"""The human-readable DISCOVERY_REPORT.html page, in the same style as the run report.

One row per candidate, grouped by the family that routed it, with its matched SSN-clade
profile(s), location and whether it was predicted; then one row per searched profile with links
to the raw BITACORA files behind it. Everything is read from stored files; discovery is never
rerun. The machine-readable candidate table stays in discovery_candidates.tsv.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import quote

from predictor.report.run_page import _clip, _e, _na

CSS = """
.profiles td{white-space:normal}
.files{columns:2;font-size:12px;margin:6px 0 0;padding-left:18px}
.files li{break-inside:avoid}
"""

GENERATED = {"DISCOVERY_REPORT.html", "discovery_candidates.tsv", "candidates.faa"}


def _href(relative: str) -> str:
    return quote(str(relative).replace("\\", "/"), safe="/")


def _inside(path: Path, root: Path) -> bool:
    try:
        return path.resolve().is_relative_to(root.resolve())
    except (OSError, ValueError):
        return False


def _folders(discovery: Path) -> list[Path]:
    return sorted(p for p in discovery.iterdir()
                  if p.is_dir() and not p.name.startswith(".") and _inside(p, discovery))


def _profiles(row: dict) -> list[str]:
    text = row.get("discovery_profile") or row.get("clade") or row.get("profile") or ""
    return [p.strip() for p in str(text).replace(",", "|").split("|") if p.strip()]


def _location(row: dict) -> str:
    contig, start, end = row.get("contig") or "", row.get("start") or "", row.get("end") or ""
    if not (start and end):
        return _na("not recorded")
    strand = row.get("strand") or ""
    return f'<span class="seq">{_e(contig)}:{_e(start)}–{_e(end)}</span> {_e(strand)}'


def _length(row: dict) -> str:
    if row.get("aa_len"):
        return _e(row["aa_len"])
    try:
        return _e(int(row["trim_end"]) - int(row["trim_start"]))
    except (KeyError, TypeError, ValueError):
        return _na()


def _score(row: dict) -> str:
    evalue, bits = row.get("evalue") or "", row.get("bits") or ""
    if not (evalue or bits):
        return _na()
    return _e(evalue) + (f' <span class="muted">{_e(bits)} bits</span>' if bits else "")


def _candidate_table(rows: list[dict], reports: dict[str, str]) -> str:
    # BITACORA records no per-candidate score, so the column appears only when some row has one.
    scored = any(r.get("evalue") or r.get("bits") for r in rows)
    labels = ("candidate", "product", "matched profile(s)", "location", "length (aa)",
              *(("best E-value",) if scored else ()), "prediction")
    head = "".join(f'<th scope="col" aria-sort="none"><button type="button">{label}</button></th>'
                   for label in labels)
    out = [f'<div class="scroll"><table class="hits"><thead><tr>{head}</tr></thead>']
    by_family: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        by_family[row.get("family") or "family not recorded"].append(row)
    for family in sorted(by_family):
        members = sorted(by_family[family], key=lambda r: (r.get("gene") or r.get("locus_tag") or r["run_name"]).lower())
        predicted = sum(r["run_name"] in reports for r in members)
        out.append(f'<tbody><tr class="group"><td colspan="{len(labels)}">{_e(family)}<span>{len(members)} '
                   f'candidate{"s" if len(members) != 1 else ""} · {predicted} predicted</span></td></tr>')
        for row in members:
            gene, locus = row.get("gene") or "", row.get("locus_tag") or ""
            protein = row.get("protein_id") or row["run_name"].split("__")[-1]
            name = gene or locus or protein
            sub = " · ".join(x for x in (locus if gene else "", protein if protein != name else "") if x)
            link = reports.get(row["run_name"])
            out.append("<tr>"
                       f'<td><span class="tf">{_e(name)}</span>'
                       + (f'<span class="locus">{_e(sub)}</span>' if sub else "") + "</td>"
                       f"<td>{_clip(row.get('product') or '', 48) or _na()}</td>"
                       f"<td>{_e(', '.join(_profiles(row))) or _na()}</td>"
                       f"<td>{_location(row)}</td><td>{_length(row)}</td>"
                       + (f'<td data-sort-value="{_e(row.get("evalue") or row.get("bits") or "")}">{_score(row)}</td>'
                          if scored else "")
                       + (f'<td><a href="{_e(link)}">open report</a></td>' if link
                          else f"<td>{_na('not predicted')}</td>") + "</tr>")
        out.append("</tbody>")
    out.append("</table></div>")
    return "".join(out)


def _profile_table(rows: list[dict], discovery: Path) -> str:
    """One row per searched profile: how many candidates it caught and its raw engine files."""
    caught = Counter(p for row in rows for p in _profiles(row))
    folders = _folders(discovery)
    # Saved profile labels are display values, never paths supplied to the filesystem.
    folder_by_name = {p.name: p for p in folders}
    names = sorted(set(caught) | {p.name for p in folders}, key=str.lower)
    if not names:
        return ""
    head = "".join(f'<th scope="col" aria-sort="none"><button type="button">{label}</button></th>'
                   for label in ("profile", "candidates", "raw files"))
    body = []
    for name in names:
        folder = folder_by_name.get(name)
        files = sorted(p for p in folder.rglob("*") if p.is_file() and _inside(p, discovery)) if folder else []
        empty = sum(not p.stat().st_size for p in files)
        listing = "".join(
            f'<li><a href="{_e(_href(p.relative_to(discovery).as_posix()))}">{_e(p.name)}</a>'
            + (' <span class="muted">(empty: no records)</span>' if not p.stat().st_size else "") + "</li>"
            for p in files)
        files_cell = (f'<details><summary>{len(files)} file{"s" if len(files) != 1 else ""}'
                      + (f", {empty} empty" if empty else "") + f'</summary><ul class="files">{listing}</ul></details>'
                      if files else _na("none retained"))
        body.append(f"<tr><td>{_e(name)}</td><td>{caught.get(name, 0) or _na('0')}</td><td>{files_cell}</td></tr>")
    return (f'<div class="scroll"><table class="hits profiles"><thead><tr>{head}</tr></thead><tbody>'
            + "".join(body) + "</tbody></table></div>")


def render(rows: list[dict], *, discovery: Path, reports: dict[str, str]) -> str:
    # Rows recovered from a raw engine table carry a protein id rather than a manifest run name.
    rows = [{**r, "run_name": r.get("run_name") or r.get("protein_id") or r.get("locus_tag") or ""}
            for r in rows]
    engines = Counter(r.get("discovery_engine") or "not recorded" for r in rows)
    families = Counter(r.get("family") or "family not recorded" for r in rows)
    profiles = {p for row in rows for p in _profiles(row)}
    searched = _folders(discovery)
    proteome = discovery / "proteins.faa"
    n_proteins = sum(line.startswith(">") for line in proteome.read_text(encoding="utf-8", errors="replace").splitlines()) \
        if proteome.is_file() and _inside(proteome, discovery) else 0
    engine_text = ", ".join(f"{k} ({v})" for k, v in sorted(engines.items()))
    lead = (f"Discovery recorded <b>{len(rows)}</b> candidate{'s' if len(rows) != 1 else ''} in "
            f"<b>{len(families)}</b> famil{'ies' if len(families) != 1 else 'y'}"
            + (f" from <b>{n_proteins:,}</b> proteins" if n_proteins else "")
            + f". Engine: {_e(engine_text)}.")
    tiles = [("candidates", len(rows)), ("families", len(families)),
             ("profiles searched", len(searched)) if searched else ("profiles matched", len(profiles)),
             ("predicted", sum(r["run_name"] in reports for r in rows))]
    body = [f'<p class="lead">{lead}</p>',
            '<div class="tiles">' + "".join(f'<div class="tile"><strong>{v}</strong><span>{_e(k)}</span></div>'
                                           for k, v in tiles) + "</div>",
            '<p class="downloads">'
            + ('<a href="../RUN_REPORT.html">Run report</a> ' if (discovery.parent / "RUN_REPORT.html").is_file()
               and _inside(discovery.parent / "RUN_REPORT.html", discovery.parent) else "")
            + '<a href="discovery_candidates.tsv">Candidate table (TSV)</a>'
            + (' <a href="candidates.faa">Candidate proteins (FASTA)</a>' if (discovery / "candidates.faa").is_file()
               and _inside(discovery / "candidates.faa", discovery) else "")
            + (' <a href="bitacora.log">BITACORA log</a>' if (discovery / "bitacora.log").is_file()
               and _inside(discovery / "bitacora.log", discovery) else "")
            + "</p>",
            "<h2>Candidates</h2>"]
    if rows:
        body.append(_candidate_table(rows, reports))
        body.append('<p class="notes"><b>Matched profile(s)</b>: the packaged SSN-clade profiles that caught '
                    "the protein; a protein caught by two profiles keeps both. <b>Family</b> is a routing hint: "
                    "the per-TF prediction makes its own family call. <b>Location</b> and <b>length</b> come "
                    "from the genome annotation; a <b>best E-value</b> column appears only when the engine "
                    "table recorded scores. Missing legacy evidence is shown as such; rebuilding this page never reruns "
                    "discovery.</p>")
    else:
        body.append('<p class="muted">No candidates were recorded.</p>')
    table = _profile_table(rows, discovery)
    if table:
        body += ["<h2>Profiles and raw evidence</h2>",
                 '<p class="notes">Each searched profile, the candidates it caught, and the raw BITACORA '
                 "files behind it. An empty file means that search retained no records.</p>", table]
    loose = sorted(p for p in discovery.iterdir() if p.is_file() and p.name not in GENERATED
                   and not p.name.startswith("proteins.faa.") and _inside(p, discovery))
    if loose:
        body.append("<details><summary>Other retained files</summary><ul class=\"files\">" + "".join(
            f'<li><a href="{_e(_href(p.name))}">{_e(p.name)}</a></li>' for p in loose) + "</ul></details>")
    return "".join(body)
