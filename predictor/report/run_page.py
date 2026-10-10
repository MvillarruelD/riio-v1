"""The human-readable RUN_REPORT.html page: one row per predicted TF, plain-language provenance.

Modelled on the survey's supplementary hit table: a TF per row, a handful of purposeful columns,
colour chips for class and ion, the primary operator in monospace, an own-region flag and the
regulated gene. Candidates that were discovered but never predicted are listed apart, so they are
never mistaken for abstentions. Everything is read from stored files; nothing is recomputed.
The complete per-candidate field set stays in Supplementary Table S1 (.tsv/.xlsx).
"""
from __future__ import annotations

import html
import json
from collections import Counter
from datetime import datetime, timezone

#: Shared with the manuscript figures (reports/figure_style.py) so an ion reads the same everywhere.
ION_COLOR = {"Zn2+": "#56B4E9", "Cu+": "#C96144", "Cd2+": "#8281B9", "Cd2+/Pb2+": "#4B4A8C",
             "As(III)": "#009E73", "Fe2+": "#9C5732", "Mn2+": "#CC79A7", "Ni2+": "#E69F00",
             "Co2+": "#F0C75E", "MoO42-": "#7A7A3A"}
ION_FALLBACK = "#9AA3AF"
#: Display spellings the inducer vocabulary emits for the same species as an ION_COLOR key.
ION_ALIAS = {"As3+": "As(III)", "arsenite": "As(III)", "Cu(I)": "Cu+", "Fe(II)": "Fe2+",
             "MoO4 2-": "MoO42-", "molybdate": "MoO42-"}
FAIL = ' class="fail"'  # row attribute for candidates that did not complete
CLASS_COLOR = {"metal": "#0F4D92", "redox": "#42949E", "organic": "#B9A17B", "unresolved": "#D8D8D8"}
CHANNEL = {"ssn_cluster": "SSN clade", "coordination": "coord.", "metalnet": "MetalNet2",
           "ligify": "Ligify", "ligify_db": "Ligify", "neighborhood": "neighbourhood",
           "neighbourhood": "neighbourhood"}
MODE = {"genome_accession": "Genome accession", "local_genome": "Local genome files",
        "genome_fasta": "Local genome files", "genbank": "Local GenBank file",
        "proteome": "Local proteome", "protein_list": "Supplied protein list",
        "single_protein": "Single protein"}

CSS = """
.lead{font-size:16px;max-width:980px}
.tiles{display:flex;gap:10px;flex-wrap:wrap;margin:18px 0 6px}
.tile{background:white;border:1px solid #dce5e7;border-radius:6px;padding:10px 16px;min-width:110px}
.tile strong{display:block;font-size:22px}.tile span{color:#637780;font-size:12px}
.downloads a{margin-right:14px}
table.hits{font-size:13px;border-top:2px solid #20313b}
table.hits th{background:white;font-size:12px;color:#20313b;border-bottom:1px solid #20313b;
vertical-align:bottom;padding:6px 8px}
table.hits td{padding:5px 8px;border-bottom:1px solid #eef2f3;overflow-wrap:normal;white-space:nowrap}
table.hits th{white-space:nowrap}table.hits td.wrap{white-space:normal;min-width:220px;max-width:340px}
table.hits tbody tr:nth-child(even){background:#f7f9f9}
table.hits tr.group td{background:white;font-weight:bold;font-size:14px;padding-top:14px;border-bottom:1px solid #b9c7cb}
table.hits tr.group span{font-weight:normal;color:#637780;font-size:12px;margin-left:8px}
.chip{display:inline-block;width:10px;height:10px;margin-right:6px;vertical-align:-1px;border-radius:1px}
.tf{font-style:italic;font-weight:bold}.locus{color:#637780;font-size:11px;display:block}
.seq{font-family:Consolas,'DejaVu Sans Mono',monospace;font-size:12px;white-space:nowrap}
.na{color:#9aa7ad}.muted{color:#637780}.yes{font-weight:bold}.gene{font-style:italic}
.fail td{background:#fdf0ec}
.notes{font-size:12px;color:#637780;max-width:1100px;line-height:1.5}
dl.prov{display:grid;grid-template-columns:200px 1fr;gap:8px 18px;background:white;padding:16px 20px;
border:1px solid #dce5e7;border-radius:6px;max-width:1100px}
dl.prov dt{font-weight:bold}dl.prov dd{margin:0}
.hash{font-family:Consolas,'DejaVu Sans Mono',monospace;font-size:12px;color:#637780}
.distributions{display:grid;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:18px;margin:18px 0}
.distributions dl{display:grid;grid-template-columns:1fr auto;gap:5px 14px;background:white;padding:14px}
.distributions dd{margin:0;font-weight:bold}
details{margin-top:12px}summary{cursor:pointer;color:#176b79}
"""


def _e(value) -> str:
    return html.escape("" if value is None else str(value), quote=True)


def _na(text="–") -> str:
    return f'<span class="na">{_e(text)}</span>'


def _chip(color: str) -> str:
    return f'<span class="chip" style="background:{color}"></span>'


def _gene_label(row: dict, dossier: dict) -> tuple[str, str]:
    """Gene name when the annotation has one, else the locus tag; the locus tag is the subtitle."""
    gene, locus = row.get("gene") or "", row.get("locus_tag") or ""
    name = gene or locus or dossier.get("tf_id") or row["run_name"]
    sub = locus if gene and locus else ("" if name == row["run_name"] else row["run_name"])
    return name, sub


def _clip(text: str, limit: int) -> str:
    """Escape, clipping long values to `limit` characters with the full value as a tooltip."""
    if len(text) <= limit:
        return _e(text)
    return f'<span title="{_e(text)}">{_e(text[:limit - 1].rstrip())}…</span>'


def _inducer_cell(row: dict, dossier: dict) -> str:
    inducers = dossier.get("inducers") or {}
    canonical = inducers.get("top") or row.get("inducer_top") or ""
    top = inducers.get("top_display") or canonical
    cls = row.get("inducer_class")
    if row.get("abstention"):
        return _na("no call (abstained)")
    color = next((ION_COLOR[ION_ALIAS.get(v, v)] for v in (canonical, top) if ION_ALIAS.get(v, v) in ION_COLOR), None)
    if cls == "metal" and color:
        return _chip(color) + _e(top)
    if cls == "metal":
        voters = sorted({CHANNEL.get(c.get("source"), c.get("source", "")) for c in inducers.get("calls") or []
                         if c.get("role") == "metal"})
        why = f' <span class="muted">({_e(", ".join(voters))})</span>' if voters else ""
        return ("unresolved ion" if "unresolved" in top.lower() or not top else _clip(top, 34)) + why
    return _clip(top, 44)


def _own_region(dossier: dict) -> str:
    """'yes (#k)' when one of the top 5 ranked operators lies in the regulator's own flanking region."""
    operators = dossier.get("operators") or {}
    ranked = operators.get("ranked") or [] if isinstance(operators, dict) else operators
    for rank, record in enumerate(ranked[:5], 1):
        if isinstance(record, dict) and (record.get("provenance") or {}).get("locality") == "proximal":
            return f'<span class="yes">yes (#{rank})</span>'
    return "no" if ranked else _na()


def _regulated(dossier: dict, primary_start) -> str:
    operators = dossier.get("operators") or {}
    primary = operators.get("primary") or {} if isinstance(operators, dict) else {}
    provenance = primary.get("provenance") or {}
    if primary_start is None or primary_start == "":
        return _na()
    for hit in dossier.get("per_hit_regulation") or []:
        operator = hit.get("operator") or {}
        if (operator.get("start") == primary_start
                and all(provenance.get(key) is None or operator.get(key) is None
                        or provenance[key] == operator[key] for key in ("end", "strand"))):
            gene, product = hit.get("regulated_gene"), hit.get("gene_product") or ""
            if not gene:
                break
            return (f'<span class="gene">{_e(gene)}</span>'
                    + (f' <span class="muted">{_clip(product, 34)}</span>' if product else ""))
    return _na()


def _head(labels) -> str:
    return "<tr>" + "".join(f'<th scope="col" aria-sort="none"><button type="button">{label}</button></th>'
                            for label in labels) + "</tr>"


def _predicted_table(predicted: list[dict]) -> str:
    head = _head(("TF", "class", "SSN clade", "inducer called", "primary operator (5'→3')",
                  "own-region<br>operator", "regulated gene", "report"))
    out = ['<div class="scroll"><table class="hits"><thead>', head, "</thead>"]
    by_family: dict[str, list[dict]] = {}
    for row in predicted:
        by_family.setdefault(row.get("family") or "unknown", []).append(row)
    for family in sorted(by_family):
        rows = sorted(by_family[family], key=lambda r: (_gene_label(r, r["_dossier"])[0].lower()))
        metal = sum(r.get("inducer_class") == "metal" for r in rows)
        out.append(f'<tbody><tr class="group"><td colspan="8">{_e(family)}<span>{len(rows)} predicted · '
                   f'{metal} metal call{"s" if metal != 1 else ""}</span></td></tr>')
        for row in rows:
            d = row["_dossier"]
            name, sub = _gene_label(row, d)
            cls = row.get("inducer_class") or "unresolved"
            primary = row.get("primary_operator")
            report = row.get("report_path")
            failed = not str(row.get("status", "")).startswith("complete")
            out.append(f'<tr{FAIL if failed else ""}>'
                       f'<td><span class="tf">{_e(name)}</span>'
                       + (f'<span class="locus">{_e(sub)}</span>' if sub else "") + "</td>"
                       f"<td>{_chip(CLASS_COLOR.get(cls, CLASS_COLOR['unresolved']))}{_e(cls)}</td>"
                       f"<td>{_e(row.get('ssn_cluster')) or _na()}</td>"
                       f'<td class="wrap">{_inducer_cell(row, d)}</td>'
                       f'<td class="seq">{_e(primary) if primary else _na("no operator")}</td>'
                       f"<td>{_own_region(d)}</td>"
                       f"<td>{_regulated(d, row.get('primary_start'))}</td>"
                       + "<td>" + (f'<a href="{_e(row["_href"])}">open</a>' if report else "")
                       + (f'<span class="locus">{_e(row.get("status"))}</span>' if failed else "")
                       + "</td></tr>")
        out.append("</tbody>")
    out.append("</table></div>")
    return "".join(out)


def _unpredicted_table(rows: list[dict]) -> str:
    head = _head(("locus tag", "gene", "product", "discovery family", "matched profile", "status"))
    body = []
    for row in sorted(rows, key=lambda r: (r.get("discovery_family") or "", r.get("locus_tag") or r["run_name"])):
        failed = str(row.get("status", "")).startswith(("failed", "timed out"))
        body.append(f'<tr{FAIL if failed else ""}>'
                    f"<td>{_e(row.get('locus_tag') or row['run_name'])}</td>"
                    f'<td class="gene">{_e(row.get("gene")) or _na()}</td>'
                    f"<td>{_clip(row.get('product') or '', 60) or _na()}</td>"
                    f"<td>{_e(row.get('discovery_family')) or _na()}</td>"
                    f"<td>{_e((row.get('discovery_profile') or '').replace('|', ', ')) or _na()}</td>"
                    f"<td>{_e(row.get('status'))}</td></tr>")
    return ('<div class="scroll"><table class="hits"><thead>' + head + "</thead><tbody>"
            + "".join(body) + "</tbody></table></div>")


def _runtime(row: dict) -> dict:
    """Bundle runtime fingerprint; legacy bundles stored it under `runtime_fingerprint`."""
    manifest = row["_bundle_manifest"]
    value = manifest.get("runtime") or manifest.get("runtime_fingerprint") or {}
    return value if isinstance(value, dict) else {"recorded": value}


def _when(stamp: str | None) -> str:
    try:
        value = datetime.fromisoformat(stamp[:-1] + "+00:00" if isinstance(stamp, str) and stamp.endswith("Z") else stamp)
        value = value.replace(tzinfo=timezone.utc) if value.tzinfo is None else value.astimezone(timezone.utc)
        return value.strftime("%Y-%m-%d %H:%M UTC")
    except (TypeError, ValueError):
        return ""


def _setting_line(records: list[dict]) -> str:
    """Effective settings, stated once when every bundle agrees and with counts when they do not."""
    labels = [("allow_online", "online retrieval"), ("allow_metalnet", "MetalNet2 metal-site channel"),
              ("allow_folds", "optional ESMFold2 fold"), ("emit_af3_jobs", "AlphaFold 3 job export")]
    parts = []
    for key, label in labels:
        values = Counter(bool((r.get("effective_config") or {}).get(key)) for r in records
                         if key in (r.get("effective_config") or {}))
        if not values:
            continue
        if len(values) == 1:
            parts.append(f"{label} {'on' if next(iter(values)) else 'off'}")
        else:
            parts.append(f"{label} on in {values[True]}, off in {values[False]}")
    return "; ".join(parts) or "not recorded in these bundles"


def _provenance(run_input: dict, rows: list[dict], predicted: list[dict]) -> str:
    items = []
    mode = (run_input or {}).get("mode") or ""
    source = (run_input or {}).get("accession") or (run_input or {}).get("organism") or ""
    files = (run_input or {}).get("files") or []
    what = MODE.get(mode, mode.replace("_", " ") or "")
    if source:
        what += f" <b>{_e(source)}</b>"
    if not what:
        # Runs made before v1.1 stored no run_input.json; the candidate manifest still names the genome.
        genomes = sorted({r.get("genome") or r.get("organism_acc") for r in rows} - {None, ""})
        what = (("Genome " + ", ".join(f"<b>{_e(g)}</b>" for g in genomes) + " (from the candidate "
                 "manifest; this run predates the stored run input)") if genomes else "not recorded")
    items.append(("Input", what))
    if files:
        listed = "".join(f'<div>{_e(f.get("role") or f.get("kind") or "file")}: {_e(f.get("name") or f.get("path"))}'
                         f' <span class="hash">sha256 {_e(str(f.get("sha256", ""))[:12])}</span></div>'
                         for f in files if isinstance(f, dict))
        items.append(("Input files", listed or _e(files)))
    engines = Counter(r.get("discovery_engine") or "not recorded" for r in rows)
    policy = (run_input or {}).get("discovery_policy") or {}
    engine_text = ", ".join(f"{_e(k)} ({v})" for k, v in sorted(engines.items()))
    if policy:
        engine_text += (" · Pfam census fallback " + ("allowed" if policy.get("census_fallback") else "not allowed")
                        + (" · census only" if policy.get("census_only") else ""))
    items.append(("Discovery", engine_text))
    runtimes = [_runtime(r) for r in predicted]
    versions = Counter((rt.get("package_version"), rt.get("git_commit"), rt.get("git_dirty"),
                        rt.get("data_release")) for rt in runtimes if rt)
    if versions:
        lines = []
        for (version, commit, dirty, data), n in versions.most_common():
            text = f"RIIO {_e(version)}"
            if commit:
                text += f' from source <span class="hash">{_e(commit[:10])}</span>'
                text += " with uncommitted changes" if dirty else ""
            text += f" · reference data {_e(data)}" if data else ""
            lines.append(text + (f" ({n} of {len(predicted)} TFs)" if len(versions) > 1 else ""))
        items.append(("Software", "<br>".join(lines)))
    records = [r["_run_record"] for r in predicted if r.get("_run_record")]
    if records:
        items.append(("Settings", _e(_setting_line(records))))
        starts = sorted(filter(None, (_when(r.get("started_at_utc")) for r in records)))
        ends = sorted(filter(None, (_when(r.get("published_at_utc")) for r in records)))
        if starts and ends:
            items.append(("Predicted", f"{starts[0]} to {ends[-1]}"))
        genomes = Counter(((r.get("genome") or {}).get("accession"), (r.get("genome") or {}).get("sha256"))
                          for r in records)
        items.append(("Genome scanned", "<br>".join(
            f"{_e(acc or 'not recorded')} <span class=\"hash\">sha256 {_e((sha or '')[:12])}</span>"
            + (f" ({n} TFs)" if len(genomes) > 1 else "") for (acc, sha), n in genomes.items())))
    dl = "".join(f"<dt>{_e(k)}</dt><dd>{v}</dd>" for k, v in items)
    # Each bundle keeps its full run record; embedding them all made a 39-TF page ~1 MB (97 % JSON).
    raw = {"run_input": run_input or "not recorded",
           "runtime_by_bundle": {r["run_name"]: _runtime(r) or "not recorded" for r in predicted}}
    records = "".join(f'<li><a href="{_e(r["_record_href"])}">{_e(r["run_name"])}</a></li>'
                      for r in predicted if r.get("_record_href"))
    return (f'<dl class="prov">{dl}</dl><details><summary>Machine-readable record (JSON)</summary><pre>'
            + _e(json.dumps(raw, indent=2, ensure_ascii=False)) + "</pre></details>"
            + (f"<details><summary>Per-TF run records ({sum(bool(r.get('_record_href')) for r in predicted)})</summary><ul>{records}</ul></details>"
               if records else ""))


def _distributions(predicted: list[dict]) -> str:
    sections = []
    for title, key in (("Saved family calls", "family"), ("Saved inducer classes", "inducer_class"),
                       ("Saved operator confidence", "operator_confidence")):
        recorded = [str(row[key]) for row in predicted if row.get(key) is not None and row.get(key) != ""]
        if recorded:
            counts = Counter(recorded)
            detail = '<dl>' + ''.join(f'<dt>{_e(label)}</dt><dd>{count}</dd>'
                                     for label, count in sorted(counts.items())) + '</dl>'
            if len(recorded) < len(predicted):
                detail += f'<p class="muted">Not recorded for {len(predicted) - len(recorded)} predicted TFs.</p>'
        else:
            detail = '<p class="muted">Unavailable: no saved values in these bundles.</p>'
        sections.append(f'<div><h3>{title}</h3>{detail}</div>')
    return '<div class="distributions">' + ''.join(sections) + '</div>'


def render(rows: list[dict], *, run_input: dict, links: list[str], af3_warnings) -> str:
    predicted = [r for r in rows if r["_predicted"]]
    others = [r for r in rows if not r["_predicted"]]
    failed = [r for r in rows if str(r.get("status", "")).startswith(("failed", "timed out"))]
    metal = sum(r.get("inducer_class") == "metal" for r in predicted)
    abstained = sum(bool(r.get("abstention")) for r in predicted)
    own = sum("yes" in _own_region(r["_dossier"]) for r in predicted)
    mode = MODE.get((run_input or {}).get("mode") or "", "")
    source = (run_input or {}).get("accession") or ""
    discovered = any((r.get("discovery_engine") or "") not in ("", "provided_protein") for r in rows)
    lead = (f"{'Discovery found' if discovered else 'The run received'} <b>{len(rows)}</b> candidate "
            f"transcription factor{'s' if len(rows) != 1 else ''}; <b>{len(predicted)}</b> "
            f"{'was' if len(predicted) == 1 else 'were'} predicted"
            + (f" and <b>{len(others)}</b> {'was' if len(others) == 1 else 'were'} not" if others else "")
            + f". {metal} of the predicted TFs received a metal inducer call"
            + (f"; {abstained} abstained" if abstained else "") + ".")
    tiles = [("candidates", len(rows)), ("predicted", len(predicted)), ("not predicted", len(others)),
             ("failed", len(failed)), ("metal calls", metal), ("own-region operator", own)]
    body = []
    if mode or source:
        body.append(f'<p class="muted">{_e(mode)} {_e(source)}</p>')
    body.append(f'<p class="lead">{lead}</p>')
    body.append('<div class="tiles">' + "".join(f'<div class="tile"><strong>{v}</strong><span>{_e(k)}</span></div>'
                                               for k, v in tiles) + "</div>")
    body.append('<p class="downloads">' + " ".join(links) + "</p>")
    body.append(_distributions(predicted))
    body.append("<h2>Predicted transcription factors</h2>")
    if predicted:
        body.append(_predicted_table(predicted))
        body.append(
            '<p class="notes">Use a column heading to sort within each family; headings support keyboard activation. '
            '<b>Class</b> and <b>inducer called</b>: the fused inducer call; an unresolved '
            "ion lists the channels that called metal. <b>Primary operator</b>: the saved rank-1 operator "
            "on its recorded strand. <b>Own-region operator</b>: one of the top 5 operators lies in the "
            "regulator's own flanking intergenic region (its rank in brackets). <b>Regulated gene</b>: the "
            "gene the primary operator is predicted to control, with its product. Rows shaded red did not "
            "complete. Every value is an unverified prediction; open a TF's report for its evidence.</p>")
    else:
        body.append('<p class="muted">No candidate has a prediction bundle in this run.</p>')
    if others:
        body.append(f"<h2>Discovered, not predicted ({len(others)})</h2>")
        body.append('<p class="notes">These candidates were found but have no prediction: they were not '
                    "selected (--select or the GUI subset), are still to run, or failed. They are not "
                    "abstentions. Rerunning the same command predicts them and reuses completed bundles.</p>")
        table = _unpredicted_table(others)
        body.append(table if len(others) <= 15 else
                    f"<details><summary>Show all {len(others)} candidates</summary>{table}</details>")
    if af3_warnings:
        body.append("<h2>AlphaFold 3 export warnings</h2><pre>" + _e(json.dumps(af3_warnings, indent=2)) + "</pre>")
    body.append("<h2>How this run was made</h2>")
    body.append(_provenance(run_input, rows, predicted))
    body.append('<p class="muted">Rebuilt from stored data; no prediction or discovery was rerun. '
                "Supplementary Table S1 holds every field for every candidate.</p>")
    return "".join(body)
