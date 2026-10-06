"""
render_report.py -- assemble the per-TF REPORT.html (the human-readable record) + its figures.

Reads a finished job folder (`results/jobs/<TF>/`: `dossier.json` + the subfolder JSONs), generates the
section figures into `figures/`, and writes a self-contained, sectioned `REPORT.html` with a verdict banner.
Re-runnable on a CACHED folder without touching the pipeline:

    python -m predictor.report.render_report --tf CzrA_Saureus

Sections: (1) the TF, (2) inducer and effector-associated genes, (3) operators and proposed operons, (4) the final
operator list with source + rationale, and (5) recorded files and provenance. Absent optional inputs
produce labeled placeholders; malformed files and unexpected rendering failures stop publication.

"""
from __future__ import annotations

import argparse
import copy
import html
import json
import shutil
import tempfile
import uuid
from pathlib import Path
from urllib.parse import quote

from predictor import resources
from predictor.report.bundle_contract import inventory, validate_files

from predictor.report import figures as F                           # noqa: E402
from predictor.report import figure_data as FD                       # noqa: E402
from predictor.report import operators as OPS                        # noqa: E402

JOBS = resources.output_path("jobs")
_IC = None


def _load(path: Path):
    """Allow absent optional files, but never hide malformed JSON or read failures."""
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except FileNotFoundError:
        return None
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"{path.name} must contain valid UTF-8 JSON: {exc}") from exc


def collect(job_dir: Path) -> dict:
    """Load the dossier and the TF summary written beside it."""
    job_dir = Path(job_dir)
    dossier = _load(job_dir / "dossier.json")
    if not isinstance(dossier, dict) or not dossier.get("tf_id"):
        raise ValueError("dossier.json must contain an object with a tf_id")
    return {
        "dossier": dossier,
        "job_dir": job_dir.resolve(),
        "tf_summary": _load(job_dir / "tf" / "tf_summary.json"),
        "run_record": _load(job_dir / "run" / "run_record.json") or {},
        "evidence_index": _load(job_dir / "run" / "evidence_index.json") or [],
        "report_environment": _load(job_dir / "report_environment.json") or {},
        "bundle_files": {p.relative_to(job_dir).as_posix() for p in job_dir.rglob("*") if p.is_file()},
    }


# --------------------------------------------------------------------------- figure generation
# NB: operator_pssm is still GENERATED (kept as a figure for other uses) but NOT embedded in the report;
# operator_context (genome neighborhood + promoter-zoom footprint) is the combined §3 figure.

def _gen_figures(data: dict, fig_dir: Path) -> dict:
    """Generate every report figure."""
    d = data["dossier"]
    figs: dict = {}
    # §1 apo fold -- only drawn when the optional structural branch actually ran (`--fold`). With it
    # off, `structure.folded` is False and there is no model to render; a placeholder card here reads
    # as a failed fold rather than a branch that was never asked for.
    if (d.get("structure") or {}).get("folded"):
        structure = copy.deepcopy(d.get("structure") or {})
        for field in ("dimer_cif", "afdb_monomer_cif"):
            value = structure.get(field)
            if value and not Path(value).is_absolute():
                structure[field] = str((fig_dir.parent / value).resolve())
        figs["apo_dimer"] = F.fig_apo_dimer(structure, fig_dir)
    # §2 logos (genomic + three sources) + operator track
    figs["logo_genomic"] = F.fig_operator_logo(d.get("genomic_logo") or {}, fig_dir,
                                               name="operator_logo_genomic",
                                               title="reported motif (kept rescan hits)")
    if d.get("seed_logo"):
        _sl = d["seed_logo"]
        figs["logo_seed"] = F.fig_operator_logo(
            _sl, fig_dir, name="operator_logo_seed",
            title=f"seed PWM (scans the genome), width {_sl.get('anchor_w')} bp")
    ms = d.get("multi_source_logos") or {}
    for key, lab in (("tf_only", "kept genome-rescan hits"), ("homolog", "homolog promoters"),
                     ("gathered", "homologs + MSA expansion")):
        figs[f"logo_{key}"] = F.fig_operator_logo(ms.get(key) or {}, fig_dir,
                                                  name=f"operator_logo_{key}", title=lab)
    figs["operator_track"] = F.fig_operator_track(FD.operator_track_rows(d), fig_dir)
    # operator PSSM: kept as a figure (other uses) but NOT embedded in the report
    figs["operator_pssm"] = F.fig_operator_pssm(d.get("genomic_logo") or {}, fig_dir)
    # §3 inducer + operator placement (TSS and CDS anchors) + the combined genome-context figure
    figs["inducer_evidence"] = F.fig_inducer_evidence(FD.inducer_evidence_rows(d), fig_dir)
    figs["operator_vs_tss"] = F.fig_operator_vs_tss(FD.tss_offsets(d), fig_dir,
                                                    cds_rows=FD.cds_offsets(d))
    figs["operator_context"] = F.fig_operator_context(FD.neighborhood_panels(d), fig_dir)
    # A failed render must never silently reuse a figure from an earlier rendering.
    return figs


# --------------------------------------------------------------------------- HTML
_CSS = """
:root{--ink:#272727;--muted:#60666d;--paper:#f4f6f8;--surface:#fcfcfd;--line:#d8dce1;
--accent:#1f3864;--accent-dark:#1f3864;--soft:#e9eff6;--amber:#a86b15;--red:#b34e35;
--operator:#2e8581;--shadow:rgba(31,56,100,.07)}
*{box-sizing:border-box}html{scroll-behavior:smooth;background:var(--paper)}
body{font-family:"Aptos","Segoe UI",Arial,sans-serif;color:var(--ink);max-width:1220px;margin:0 auto;
padding:28px clamp(16px,4vw,52px) 76px;line-height:1.58}main{background:var(--surface);border:1px solid var(--line);
border-radius:18px;padding:0 clamp(18px,4vw,54px) 58px;box-shadow:0 26px 70px var(--shadow)}
.skip{position:absolute;left:-999px;top:8px;background:var(--ink);color:white;padding:8px 12px;z-index:3}
.skip:focus{left:12px}.report-hero{padding:46px 0 28px;border-bottom:1px solid var(--line)}.report-hero-grid{display:grid;grid-template-columns:minmax(0,1.45fr) minmax(280px,.55fr);
gap:36px;align-items:end}.eyebrow{font-size:11px;letter-spacing:.16em;
text-transform:uppercase;color:var(--accent);font-weight:750}
h1{font-family:Charter,"Bitstream Charter",Georgia,serif;font-size:clamp(34px,5vw,62px);font-weight:650;
line-height:1.02;letter-spacing:-.045em;margin:12px 0 14px;
max-width:100%;overflow-wrap:anywhere}
h2{font-family:Charter,"Bitstream Charter",Georgia,serif;font-size:28px;font-weight:650;letter-spacing:-.025em;
border-top:1px solid var(--line);padding-top:24px;margin-top:58px;color:var(--ink);scroll-margin-top:72px}
h3{font-size:15px;color:var(--accent-dark);margin:26px 0 8px;letter-spacing:-.01em}
.sub{color:var(--muted);font-size:13px;margin-top:0;max-width:78ch;text-wrap:pretty}
.hero-primary{background:rgba(251,252,249,.86);border:1px solid var(--line);border-radius:12px;padding:17px 18px;
box-shadow:0 18px 44px rgba(35,63,50,.06)}.hero-primary b{display:block;font-size:10px;letter-spacing:.1em;
text-transform:uppercase;color:var(--muted);margin-bottom:7px}.hero-primary .sequence{font-size:15px;color:var(--ink);font-weight:700;
line-height:1.45;overflow-wrap:anywhere}.hero-primary small{display:block;margin-top:10px;color:var(--accent-dark);font-weight:650}
.run-stamp{margin-top:12px;padding-top:10px;border-top:1px solid var(--line);color:var(--muted);font-size:10px;line-height:1.55}
.report-nav{position:sticky;top:0;z-index:2;display:flex;gap:5px;overflow-x:auto;padding:9px 0;
background:rgba(251,252,248,.94);backdrop-filter:blur(12px);border-bottom:1px solid var(--line)}
.report-nav a{white-space:nowrap;text-decoration:none;color:var(--muted);font-size:12px;font-weight:650;
padding:7px 10px;border-radius:5px;transition:background .18s ease,color .18s ease}
.report-nav a:hover,.report-nav a:focus-visible{background:var(--soft);color:var(--accent-dark)}
.read-guide{border-left:3px solid var(--accent);padding:14px 17px;background:var(--soft);color:#31423a;
font-size:13px;margin:24px 0 28px}.summary-grid{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:16px;margin:16px 0 28px}
.summary-item{border-top:2px solid var(--ink);padding:10px 2px 0;min-width:0}.summary-item b{display:block;
font-size:10px;letter-spacing:.09em;text-transform:uppercase;color:var(--muted);margin-bottom:4px}.summary-item span{
display:block;font-size:14px;font-weight:650;overflow-wrap:anywhere;font-variant-numeric:tabular-nums}
.banner{display:flex;flex-wrap:wrap;gap:7px;margin:19px 0 4px}.chip{border-radius:6px;padding:7px 10px;
font-size:12px;color:var(--ink);min-width:126px;border:1px solid var(--line);border-left-width:3px;background:rgba(251,252,249,.8)}
.chip b{display:block;font-size:10px;letter-spacing:.04em;color:var(--muted);font-weight:650}.chip span{font-size:11px}
.green{border-left-color:var(--accent)}.amber{border-left-color:var(--amber)}.red{border-left-color:var(--red)}.na{border-left-color:#7b8580}
img{display:block;max-width:100%;height:auto;border:1px solid var(--line);border-radius:4px;background:var(--surface);
margin:10px 0 18px}
.table-wrap{overflow-x:auto;margin:10px 0 20px;border:1px solid var(--line);border-radius:9px;background:var(--surface)}
table{border-collapse:collapse;font-size:12px;width:100%;font-variant-numeric:tabular-nums}th,td{border-bottom:1px solid var(--line);
padding:8px 10px;text-align:left;vertical-align:top}th{background:#edf2ed;color:#3e4d45;font-size:10px;
letter-spacing:.045em;text-transform:uppercase;position:sticky;top:0}tr:last-child td{border-bottom:0}tr:nth-child(even){background:#faf9f4}
.mono{font-family:"Cascadia Mono","SFMono-Regular",Consolas,monospace;overflow-wrap:anywhere}
.ph{border:1px dashed #aab6ad;border-radius:8px;padding:15px;color:var(--muted);font-size:12px;
background:rgba(255,253,247,.62);margin:8px 0 18px}.kv{font-size:13px;display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:6px 20px}
.kv div{border-bottom:1px solid var(--line);padding:7px 0}.kv b{display:block;font-size:10px;letter-spacing:.06em;
text-transform:uppercase;color:var(--muted);margin-bottom:2px}.tag{display:inline-block;border-radius:4px;padding:2px 7px;
font-size:10px;color:#fff;font-weight:650}.natural{background:var(--accent)}.designed{background:#5d6e92}.pending{background:#7b8580}
@media(max-width:820px){.report-hero-grid{grid-template-columns:1fr;gap:18px}.summary-grid{grid-template-columns:repeat(2,minmax(0,1fr))}.kv{grid-template-columns:1fr}
body{padding:18px 16px 50px}.report-hero{padding-top:22px}.chip{min-width:108px}}
a{color:var(--accent-dark);text-underline-offset:3px}
a:focus-visible,summary:focus-visible{outline:2px solid var(--accent);outline-offset:3px}
td{overflow-wrap:anywhere}.file-list{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:8px 24px}
.file-list a{font-size:13px;overflow-wrap:anywhere}.provenance-note{font-size:12px;color:var(--muted)}
.audit-grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:11px;margin:14px 0 24px}
.audit-card{background:#f2f6f1;border:1px solid var(--line);border-radius:9px;padding:13px 14px;min-width:0}
.audit-card b{display:block;font-size:10px;letter-spacing:.08em;text-transform:uppercase;color:var(--muted);margin-bottom:5px}
.audit-card strong{display:block;font-size:17px;letter-spacing:-.02em;overflow-wrap:anywhere}.audit-card small{color:var(--muted)}
.status{display:inline-block;border-left:3px solid var(--line);padding-left:7px;font-weight:650}.status.recorded,.status.complete{border-left-color:var(--accent)}
.status.limited,.status.empty,.status.unresolved{border-left-color:var(--amber)}.status.not-run{border-left-color:#7b8580}
.run-note{background:#f6f2e8;border:1px solid #ded4bf;border-radius:8px;padding:13px 15px;color:#675533;font-size:12px;max-width:86ch}
.hash{font-size:10px;color:var(--muted);word-break:break-all}.section-lead{font-family:Charter,"Bitstream Charter",Georgia,serif;
font-size:17px;line-height:1.45;max-width:70ch;color:#34453c}.file-list a{padding:9px 0;border-bottom:1px solid var(--line);text-decoration:none}
.file-list a:hover{color:var(--accent);padding-left:3px}.compact td,.compact th{padding-top:6px;padding-bottom:6px}
@media(max-width:600px){.file-list{grid-template-columns:1fr}.report-nav{position:static}h2{margin-top:36px}}
@media(max-width:820px){.audit-grid{grid-template-columns:1fr}}
@media(prefers-reduced-motion:reduce){html{scroll-behavior:auto}*{transition:none!important}}
@media print{html,body,main{background:white}.report-nav,.skip{display:none}body{max-width:none;padding:0;color:#111}main{border:0;box-shadow:none;padding:0}
h1{font-size:24pt}h2{break-after:avoid}.table-wrap{overflow:visible}.summary-grid{grid-template-columns:repeat(4,minmax(0,1fr))}
th{position:static}thead{display:table-header-group}tr{break-inside:avoid}
img{break-inside:avoid;max-height:240mm;object-fit:contain}}
body{font-family:Arial,"Segoe UI",sans-serif}
main{border-radius:7px;box-shadow:0 12px 40px var(--shadow)}
h1,h2,.section-lead{font-family:Arial,"Segoe UI",sans-serif}
h1{font-size:clamp(32px,4vw,52px);font-weight:700}
h2{font-size:25px;font-weight:700}
h3{font-size:15px;margin-top:29px}
.hero-primary{background:var(--soft);border-radius:5px;box-shadow:none}
.report-nav{background:rgba(252,252,253,.97);backdrop-filter:none}
.read-guide{border:1px solid #c9d7e7;border-radius:5px;color:var(--ink)}
.summary-item{border-top-color:var(--accent)}
.chip{border-width:1px 1px 1px 4px;border-radius:4px;background:var(--surface)}
.chip.green{border-color:var(--line);border-left-color:var(--operator)}
.chip.amber{border-color:var(--line);border-left-color:var(--amber)}
.chip.red{border-color:var(--line);border-left-color:var(--red)}
.chip.na{border-color:var(--line);border-left-color:#9aa1a8}
.chip .mark{font-style:normal;font-weight:800;margin-right:4px}
.chip.green .mark{color:var(--operator)}.chip.amber .mark{color:var(--amber)}
.chip.red .mark{color:var(--red)}.chip.na .mark{color:var(--muted)}
.chip-key{font-size:11px;color:var(--muted);margin:6px 0 0}.chip-key i{font-style:normal;font-weight:800;margin:0 3px 0 9px}
.nb{display:inline-block;max-width:100%}
.sr-only{position:absolute;left:-9999px;width:1px;height:1px;overflow:hidden;white-space:nowrap}
th{background:#edf1f6;color:var(--accent)}tr:nth-child(even){background:#f7f9fb}
.status{border-left-width:1px}.status.recorded,.status.complete{border-left-color:var(--operator)}
.natural{background:var(--operator)}.audit-card{background:#f5f7fa;border-radius:5px}
.section-lead{font-size:15px;line-height:1.5;font-weight:400;color:var(--ink)}
.source-line{font-size:11px;color:var(--muted);margin:5px 0 21px;line-height:1.6}
.source-line b{letter-spacing:.06em;text-transform:uppercase;font-size:10px;margin-right:9px}
.source-line a{display:inline-block;margin-right:12px;white-space:nowrap}
.visual-stack{display:grid;grid-template-columns:1fr;gap:14px;margin:14px 0}
.visual-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:18px;align-items:start;margin:10px 0 14px}
.visual-figure{min-width:0;margin:0}
.visual-stack img,.visual-grid img{margin:0 0 5px;width:100%;object-fit:contain}
.figure-caption{font-size:12px;color:var(--muted);line-height:1.45}
.figure-caption b{color:var(--ink)}
.compact-figure{max-width:480px}.compact-figure img{max-height:230px;object-fit:contain}
.definition{font-size:12px;color:var(--muted);line-height:1.55;max-width:85ch;margin:6px 0 18px}
.definition strong{color:var(--ink)}
.evidence-detail summary{cursor:pointer;color:var(--accent);font-weight:650;white-space:nowrap}
.evidence-detail ul{margin:8px 0 0;padding-left:16px;max-width:60ch;line-height:1.5}
.evidence-detail li{overflow-wrap:anywhere}
@media(max-width:700px){.visual-grid{grid-template-columns:1fr}}
@media print{.source-line a{color:var(--ink)}main{box-shadow:none}}
"""


def _img(figs, key, rel="figures", *, alt="") -> str:
    p = figs.get(key)
    if p:
        return f'<img src="{html.escape(rel + "/" + Path(p).name, quote=True)}" alt="{html.escape(alt)}" loading="lazy">'
    return f'<div class="ph">— {html.escape(alt or key)} not available —</div>'


def _table(headers, rows) -> str:
    h = "".join(f"<th>{html.escape(str(x))}</th>" for x in headers)
    body = ""
    for r in rows:
        body += "<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>"
    return f'<div class="table-wrap"><table><thead><tr>{h}</tr></thead><tbody>{body}</tbody></table></div>'


def _esc(x) -> str:
    return html.escape("" if x is None else str(x))


def _artifact(path: str | None, label: str | None, files: set[str]) -> str:
    """Link only to an artifact present in this exact bundle."""
    label = label or path or "—"
    if path and path in files:
        return f'<a href="{quote(path, safe="/")}">{_esc(label)}</a>'
    return _esc(label)


def _sources(files: set[str], *paths: str) -> str:
    """Place links to the exact bundle records beside the result they explain."""
    links = [_artifact(path, path, files) for path in paths if path in files]
    return ('<p class="source-line"><b>Source files</b>' + ' · '.join(links) + '</p>') if links else ''


def _captioned_figure(figs: dict, key: str, title: str, detail: str = "") -> str:
    """A figure with a visible title, a one-line reading note and accessible alternative text."""
    if not figs.get(key):
        return ''
    note = f' {_esc(detail)}' if detail else ''
    return (f'<div class="visual-figure">{_img(figs, key, alt=title)}'
            f'<div class="figure-caption"><b>{_esc(title)}.</b>{note}</div></div>')


_CHIP_MARK = {"green": ("✓", "as expected"), "amber": ("!", "check"),
              "red": ("✗", "weak or absent"), "na": ("–", "not run")}


def _chip_html(chip: dict) -> str:
    """A verdict chip whose status is carried by a glyph and hidden text, not by colour alone."""
    mark, meaning = _CHIP_MARK.get(chip.get("status"), ("?", "unknown"))
    return (f'<div class="chip {_esc(chip.get("status"))}"><b><i class="mark" aria-hidden="true">{mark}</i>'
            f'{_esc(chip.get("label"))}<span class="sr-only"> ({meaning})</span></b>'
            f'<span>{_esc(chip.get("detail"))}</span></div>')


def _fmt_q(value) -> str:
    """q- and p-values at two significant figures; full precision stays in the TSV/JSON."""
    if isinstance(value, (int, float)) and value == value:
        return f"{value:.1e}" if value < 0.001 else f"{value:.3f}"
    return "—"


def _fmt_num(value, digits: int = 2) -> str:
    return f"{value:.{digits}f}" if isinstance(value, (int, float)) else _esc(value if value is not None else "—")


def _evidence_list(evidence: dict) -> str:
    """Render one source's recorded evidence as key: value lines rather than a JSON blob."""
    items = []
    for key, value in sorted((evidence or {}).items()):
        if value in (None, "", [], {}):
            continue
        if isinstance(value, dict):
            text = "; ".join(f"{k}: {v}" for k, v in value.items())
        elif isinstance(value, (list, tuple)):
            text = ", ".join(str(v) for v in value)
        else:
            text = str(value)
        items.append(f'<li><b>{_esc(key)}</b>: {_esc(text)}</li>')
    return f'<ul class="mono">{"".join(items)}</ul>' if items else '<ul><li>no fields recorded</li></ul>'


def _title_segments(tf: str) -> list[str]:
    """Split a run name after each "__" so the title wraps between its fields."""
    parts = str(tf).split("__")
    return [p + "__" for p in parts[:-1]] + [parts[-1]]


def _status_class(value) -> str:
    text = str(value or "").lower()
    if text in {"recorded", "complete"}:
        return text
    if "not run" in text or "available" in text:
        return "not-run"
    if "empty" in text or "no hit" in text:
        return "empty"
    if "limit" in text:
        return "limited"
    return "unresolved"


def build_html(data: dict, figs: dict) -> str:
    d = data["dossier"]
    run = data.get("run_record") or {}
    files = data.get("bundle_files") or set()
    tf = d.get("tf_id", "TF")
    fam = d.get("family", "?")
    # ---- banner
    chips = "".join(_chip_html(c) for c in FD.verdict(d))
    final_rows = FD.final_operator_rows(d)
    headline = next((r for r in final_rows if r.get("sequence") and r.get("status") == "ready"),
                    next((r for r in final_rows if r.get("sequence")), {}))
    primary_sequence = headline.get("sequence") or "not resolved"
    n_sites = len((d.get("rescan") or {}).get("hits") or [])
    run_id = run.get("run_id") or (d.get("run_record") or {}).get("run_id") or "legacy bundle"
    captured = run.get("started_at_utc") or (d.get("run_record") or {}).get("started_at_utc") or "not recorded"
    inducer_display = FD.inducer_display(d.get("inducers")) or "unresolved"
    parts = ['<header class="report-hero"><div class="report-hero-grid"><div>',
             '<span class="eyebrow">TF operator report · computational prediction</span>',
             # Break long run names at their "__" separators, never inside an accession. <wbr> is
             # outside the bundle contract's HTML subset, so each segment is an inline-block span.
             "<h1>" + "".join(f'<span class="nb">{_esc(seg)}</span>'
                              for seg in _title_segments(tf)) + "</h1>",
             f'<p class="sub">{_esc(fam)} · genome {_esc(d.get("genome_accession") or "unresolved")} · '
             f'inducer {_esc(inducer_display)} · '
             f'SSN {_esc(d.get("ssn_cluster") or "unresolved")}</p></div>',
             '<aside class="hero-primary"><b>Top-ranked operator</b>'
             f'<div class="sequence mono">{_esc(primary_sequence)}</div>'
             '<small>Computational prediction · validate experimentally</small>'
             f'<div class="run-stamp"><b>Run record</b><span class="mono">{_esc(run_id)}</span><br>'
             f'captured {_esc(captured)}</div></aside></div>',
             f'<div class="banner">{chips}</div>'
             '<p class="chip-key">Run checks, not validation:<i>✓</i>as expected<i>!</i>check'
             '<i>✗</i>weak or absent<i>–</i>not run</p></header>',
             '<nav class="report-nav" aria-label="Report sections">'
             '<a href="#decision">Summary</a><a href="#tf">TF</a><a href="#regulon">Inducer</a>'
             '<a href="#operators">Operators &amp; regulon</a>'
             '<a href="#final">Final operators</a><a href="#run">Run record</a>'
             '<a href="#files">Files</a></nav>']
    parts.append(
        '<div class="read-guide" id="decision">'
        '<b>How to read this report.</b> Section 2 shows how the inducer was called and which genome '
        'genes are annotated as handling it. Section 3 is independent of that call: a motif scan of the '
        'genome, its putative operators and the operons they would regulate. Section 4 lists the '
        'sequences to test first. The source-file links under each result open the exact record in '
        'this folder. Every result is a computational prediction, not experimental validation.</div>'
    )
    parts.append(
        '<div class="summary-grid">'
        f'<div class="summary-item"><b>Inducer</b><span>{_esc(inducer_display)}</span></div>'
        f'<div class="summary-item"><b>TF family</b><span>{_esc(fam)}</span></div>'
        f'<div class="summary-item"><b>Genome</b><span>{_esc(d.get("genome_accession") or "unresolved")}</span></div>'
        f'<div class="summary-item"><b>Putative sites</b><span>{n_sites}</span></div>'
        '</div>'
    )

    # ---- §1 TF
    locus = d.get("tf_locus") or []
    tfsum = data.get("tf_summary") or {}
    seq = (tfsum.get("sequence") or (d.get("tf_record") or {}).get("sequence") or "")
    parts.append('<h2 id="tf">1 · The transcription factor</h2>')
    parts.append('<p class="section-lead">The query protein, its genome location and the reference '
                 'used to identify its family. Identity is sequence similarity to the nearest recorded '
                 'reference; it does not measure operator or inducer accuracy. TF means transcription '
                 'factor; SSN means sequence-similarity network.</p>')
    tf_record = d.get("tf_record") or {}

    def _tf_field(key: str) -> str:
        value = tfsum.get(key) or tf_record.get(key)
        return _esc(value) if value not in (None, "") else "not recorded"

    identity = tfsum.get("nearest_identity") or tf_record.get("nearest_identity")
    identity_txt = (f"{100 * identity:.0f} %" if isinstance(identity, (int, float)) and identity <= 1
                    else _esc(identity) if identity is not None else "not recorded")
    organism = tfsum.get("organism") or tf_record.get("organism")
    # Genome-first runs record the assembly accession in this field; label it as what it is.
    organism_label = "assembly" if str(organism or "").startswith(("GCF_", "GCA_")) else "organism"
    parts.append('<div class="kv">'
                 f'<div><b>family</b>{_esc(fam)}</div>'
                 f'<div><b>{organism_label}</b>{_tf_field("organism")}</div>'
                 f'<div><b>genome / locus</b>{_esc(d.get("genome_accession"))} : '
                 f'{_esc(locus[0]) if locus else "?"}-{_esc(locus[1]) if len(locus)>1 else "?"} '
                 f'({_esc(locus[2]) if len(locus)>2 else "?"})</div>'
                 f'<div><b>query UniProt</b>{_tf_field("uniprot")}</div>'
                 f'<div><b>nearest-reference UniProt</b>{_tf_field("nearest_uniprot")}</div>'
                 f'<div><b>nearest reference</b>{_tf_field("nearest_tf")}</div>'
                 f'<div><b>identity to nearest reference</b>{identity_txt}</div>'
                 f'<div><b>length</b>{len(seq)} aa</div></div>')
    classification_flags = (tfsum.get("classification_flags") or
                            (d.get("tf_record") or {}).get("flags") or [])
    if classification_flags:
        parts.append('<h3>Family-classification evidence</h3><ul class="sub">'
                     + "".join(f'<li>{_esc(flag)}</li>' for flag in classification_flags) + '</ul>')
    if seq:
        parts.append(f'<h3>Protein sequence</h3><div class="mono" style="font-size:11px;word-break:break-all">{_esc(seq)}</div>')
    parts.append(_sources(files, "input/query.fasta", "tf/tf_summary.json", "genome/genes.tsv"))
    # A missing optional model should not display an empty confidence panel.
    s = d.get("structure") or {}
    if s.get("folded"):
        parts.append(f"<h3>Optional apo {_esc(s.get('state_name') or 'homo-oligomer')} fold</h3>")
        parts.append('<p class="definition">Apo means no effector is bound. This structure is an '
                     'optional quality view and is not used to rank operators or call the inducer. '
                     'pLDDT estimates local model confidence; AFDB RMSD compares matched Cα atoms '
                     'with the AlphaFold database model.</p>')
        parts.append('<div class="compact-figure">' +
                     _captioned_figure(figs, "apo_dimer", "Predicted apo protein structure",
                                       "Chains coloured separately.") + '</div>')
        parts.append('<div class="kv">'
                     f'<div><b>folded</b>{_esc(s.get("folded"))} ({_esc(s.get("backend"))})</div>'
                     f'<div><b>mean pLDDT</b>{_esc(s.get("plddt_mean"))}</div>'
                     f'<div><b>AFDB QC</b>RMSD {_esc(s.get("qc_rmsd"))} Å over '
                     f'{_esc(s.get("qc_n_matched"))} CA '
                     f'→ {"pass" if s.get("qc_pass") else "—"}</div></div>')
        parts.append(_sources(files, "structure/structure.json"))
    else:
        parts.append('<p class="sub">No structure model recorded. Optional structure modelling '
                     'was not run for this TF.</p>')

    # Build the independent operator path now, then place it after the inducer path below.
    operator_start = len(parts)
    parts.append('<h2 id="operators">3 · Operators &amp; proposed regulon</h2>')
    parts.append('<p class="section-lead">A seed position-weight matrix (PWM), built from candidate '
                 'operators, scans the genome. Its top hits are placed against predicted promoters. A '
                 'second motif, rebuilt from the best-aligned hits, rescans the genome and assigns '
                 'operons. Both steps are computational; no site here has been tested.</p>')
    cr = d.get("autoregulatory_recovery") or d.get("cognate_recovery") or {}
    tgt = d.get("target_operator") or {}
    glog = d.get("genomic_logo") or {}
    _auto = (f'autoregulatory operator RECOVERED ({_esc(cr.get("nearest_bp"))} bp from the self-promoter)'
             if cr.get("recovered") else
             (f'no operator at the self-promoter -- {_esc(tgt.get("n_sites_added"))} operator(s) found at '
              f'target genes (weak/non-autoregulator)' if tgt.get("n_sites_added")
              else f'no operator at the self-promoter ({_esc(cr.get("nearest_bp"))} bp to nearest hit)'))
    sel = FD.selectivity(d)
    sel_txt = (f' Operator selectivity: <b>{sel["hits_per_mb"]} hits/Mb</b> ({sel["selectivity"]} — '
               f'{_esc(sel["note"])}).' if sel else "")
    parts.append(f'<p class="sub">{len((d.get("rescan") or {}).get("hits") or [])} putative operators '
                 f'genome-wide. {_auto}. The reported motif kept '
                 f'{_esc(glog.get("n_kept"))}/{_esc(glog.get("n_hits"))} hits.{sel_txt}</p>')
    # Locality split. On the RegulonDB panel 12 of 13 primary-site recoveries were at the regulator's
    # OWN promoter (tier 1) and recovery given a distal top hit was 1/16 -- so a single recovery number
    # conflates "the motif found a site" with "the locality prior pointed at the self-promoter". State
    # which one this is; it is the difference between predicting operators and detecting autoregulation.
    loc = OPS.locality_summary(d)
    if loc["n"]:
        _pl = loc["primary_locality"]
        _note = ("the top site is at the regulator's OWN promoter, so this is an autoregulation call: "
                 "the locality prior placed that window first"
                 if _pl == "proximal" else
                 "the top site is DISTAL to the regulator, so the motif, not the locality prior, chose it"
                 if _pl == "distal" else
                 "the top site carries no locality tag (a pre-refactor bundle); it is counted as neither")
        parts.append(f'<p class="sub"><b>Locality:</b> primary site is <b>{_esc(_pl)}</b> &mdash; '
                     f'{_note}. Across all hits: {loc["proximal"]} proximal (tier&nbsp;1, the '
                     f'self-promoter), {loc["distal"]} distal, {loc["unknown"]} untagged. '
                     f'A distal prediction is the stronger claim and is the number to quote when the '
                     f'claim is &ldquo;we predict operators&rdquo;.</p>')
    parts.append('<h3>Genome scan and motif support</h3>')
    parts.append('<p class="definition"><strong>Putative sites</strong> are genome positions passing '
                 'the scan p-value cutoff. <strong>Hits/Mb</strong> divides their count by genome length; '
                 'a larger number means a less selective motif. <strong>Locality tier</strong> identifies '
                 'which genomic window produced a hit, not its experimental confidence.</p>')
    parts.append('<div class="visual-stack">' + _captioned_figure(
        figs, "operator_track", "Every putative operator along the genome",
        "Stem height is the site score (log-odds, bits); the red stem is the site at the regulator's "
        "own promoter.") + '</div>' if figs.get("operator_track")
        else _img(figs, "operator_track", alt="operator position track"))
    parts.append(_sources(files, "binding_sites.tsv", "motif/seed_pwm.json"))
    parts.append("<h3>Seed and reported motifs</h3>")
    # TWO matrices, and they are routinely conflated. The SEED PWM scans the genome and therefore
    # decides which hits exist at all; the GENOMIC logo is rebuilt afterwards from the hits the seed
    # found. Their support differs by an order of magnitude -- the seed is typically 1-4 sequences,
    # the genomic logo up to 40 -- so a statement about "the motif" is ambiguous unless it says which.
    _sl, _gl = d.get("seed_logo") or {}, d.get("genomic_logo") or {}
    if _sl or _gl:
        _rows = []
        if _sl:
            _rows.append(["<b>seed PWM</b> (scans the genome)",
                          f"w={_esc(_sl.get('anchor_w'))}, cluster {_esc(_sl.get('cluster_id'))}",
                          f"<b>{_esc(_sl.get('n_effective'))}</b> of {_esc(_sl.get('n_input'))} "
                          f"discovered candidate(s)",
                          f"<span class='mono'>{_esc(_sl.get('consensus'))}</span>"])
        if _gl:
            _rows.append(["<b>reported motif</b> (rebuilt from the hits)",
                          f"{_esc(_gl.get('n_hits'))} genome hit(s)",
                          f"<b>{_esc(_gl.get('n_kept'))}</b> kept by motif emergence",
                          f"<span class='mono'>{_esc(_gl.get('consensus'))}</span>"])
        parts.append(_table(["matrix", "built over", "n sequences", "consensus"], _rows))
        parts.append('<p class="sub">These are two different matrices. The seed PWM is what actually '
                     'scanned the genome, so every hit below exists because of it; the reported motif '
                     'is what those hits look like once aligned. A claim about "the motif" should say '
                     'which one, and quote its own n.</p>')
    seed_trials = (d.get("seed_selection") or {}).get("trials") or []
    if seed_trials:
        chosen_cluster = (d.get("seed_selection") or {}).get("chosen_cluster_id")
        trial_rows = []
        for trial in seed_trials:
            chosen = trial.get("cluster_id") == chosen_cluster
            trial_rows.append([
                '<b>chosen</b>' if chosen else "considered",
                _esc(trial.get("cluster_id")), _esc(trial.get("anchor_w")),
                _esc(trial.get("n_input")), _esc(trial.get("n_effective")),
                _esc(trial.get("n_independent")), _esc(trial.get("n_hits")),
                _esc(trial.get("n_keep")), _esc(trial.get("mean_ic")),
                f'<span class="mono">{_esc(trial.get("consensus"))}</span>',
            ])
        parts.append('<h3>Every seed-width trial</h3>')
        parts.append('<p class="definition">Each trial tests a candidate motif width before choosing '
                     'the seed PWM. <strong>Input</strong> counts candidate sequences; '
                     '<strong>effective</strong> counts the sequences retained after filtering; '
                     '<strong>independent</strong> discounts duplicate support. <strong>Mean IC</strong> '
                     'is average information content in bits per column. Hits and kept are the sites '
                     'found and retained by that trial.</p>')
        parts.append(_table(["decision", "cluster", "width", "input", "effective", "independent",
                             "hits", "kept", "mean IC", "consensus"], trial_rows))
        parts.append('<p class="sub">The complete candidate members, PWM prior, IC curve and ranked '
                     'hits for each trial are recorded in <span class="mono">motif/seed_selection.json</span>.</p>')
    def _n_of(logo: dict, *keys) -> str:
        n = next((logo.get(k) for k in keys if logo.get(k) is not None), None)
        return f"n = {n}" if n is not None else "n not recorded"

    # Wide motifs (often 30+ bp) stay legible only at full width, so the two primary logos stack.
    primary_logos = [
        _captioned_figure(figs, "logo_seed", "Seed PWM (scanned the genome)",
                          f"{_n_of(_sl, 'n_effective', 'n_seqs')} sequence(s); every hit in this "
                          "section exists because of this matrix."),
        _captioned_figure(figs, "logo_genomic", "Reported motif (rebuilt from the hits)",
                          f"{_n_of(_gl, 'n_kept', 'n_seqs')} best-aligned genome hits; this matrix "
                          "seeds the regulon rescan."),
    ]
    if any(primary_logos):
        parts.append('<div class="visual-stack">' + ''.join(primary_logos) + '</div>')
    ms = d.get("multi_source_logos") or {}
    # The TF-only pool is the same kept rescan hits the reported motif is built from; showing it again
    # when the consensus matches would present one result as two lines of evidence.
    tf_only_dup = bool(ms.get("tf_only", {}).get("consensus")
                       and ms["tf_only"].get("consensus") == _gl.get("consensus"))
    secondary_specs = [("logo_homolog", "Homolog promoters",
                        f"{_n_of(ms.get('homolog') or {}, 'n_seqs')}; motif aligned de novo across the "
                        "promoters of homologous regulators, independent of this genome's scan."),
                       ("logo_gathered", "Homologs plus MSA expansion",
                        f"{_n_of(ms.get('gathered') or {}, 'n_seqs')}; the same alignment over the "
                        "full gathered promoter set.")]
    if not tf_only_dup:
        secondary_specs.insert(0, ("logo_tf_only", "Kept genome-rescan hits",
                                   f"{_n_of(ms.get('tf_only') or {}, 'n_seqs')}; this genome only."))
    secondary_logos = [_captioned_figure(figs, key, title, note) for key, title, note in secondary_specs]
    if any(secondary_logos):
        parts.append('<h3>Supporting logos from other sequence pools</h3>')
        parts.append('<p class="definition">A motif that reappears in the promoters of homologous '
                     'regulators is conserved evidence that does not depend on this genome\'s scan. '
                     'Compare these with the reported motif above; disagreement means the scan and the '
                     'homologs point to different sites.'
                     + (' The kept-hit logo is omitted because it is identical to the reported motif.'
                        if tf_only_dup else '') + '</p>')
        parts.append('<div class="visual-grid">' + ''.join(secondary_logos) + '</div>')
    parts.append(_sources(files, "motif/seed_selection.json", "motif/three_logos.json",
                          "operators/aligned_operators_tf_only.fasta",
                          "operators/aligned_operators_homolog.fasta",
                          "operators/aligned_operators_gathered.fasta"))
    # how a site is scored (text only; the PSSM figure is generated but kept out of the report)
    parts.append("<h3>How an operator site is scored</h3>")
    parts.append('<p class="sub">Each site\'s <b>score</b> is the sum of the per-column '
                 '<b>log-odds</b> (aligned operators\' base frequencies vs the genome background, in bits) '
                 'along the bound bases. An exact <b>p-value</b> comes from the full null score distribution '
                 '(integer-lattice convolution of the per-column score PMFs — the FIMO/MOODS method), then '
                 'Benjamini–Hochberg <b>q-values</b>. The genome rescan multiplies in a <b>locality prior</b> '
                 '(tier 1 flanking the regulator = 1.0, tier 2 divergent promoters = 0.6, tier 3 other '
                 'intergenic = 0.35, tier 4 coding = 0.1); the ranking score is '
                 '<span class="mono">−log10(p) + log10(prior)</span>, so a strong distal site can still '
                 'outrank a weak nearby one while locality breaks ties.</p>')
    # operator table (richest: from per_hit_regulation) — now with the promoter elements occluded
    ph = d.get("per_hit_regulation") or []
    if ph:
        rows = []
        for r in ph[:30]:
            op = r.get("operator") or {}
            occ = ", ".join(r.get("occludes") or []) or "none"
            fp = (f'{_esc(r.get("op_start_to_tss"))}…{_esc(r.get("op_end_to_tss"))}'
                  if r.get("op_start_to_tss") is not None else _esc(r.get("op_to_tss_bp")))
            gene_bp = r.get("operator_to_gene_bp")
            rows.append([f'{_esc(op.get("start"))}-{_esc(op.get("end"))}', _esc(op.get("strand")),
                         f'{op.get("score"):.1f}' if isinstance(op.get("score"), (int, float)) else "",
                         _fmt_q(op.get("qvalue")),
                         f'<span class="mono">{_esc(op.get("seq"))}</span>',
                         _esc(r.get("regulated_gene") or "—"), f'<b>{_esc(occ)}</b>',
                         _esc(gene_bp) if gene_bp is not None else "—", fp or "—",
                         _esc(r.get("mode") or "—")])
        parts.append(f"<h3>Top {len(rows)} operators: where each sits, and what it covers</h3>")
        parts.append(_table(["locus", "strand", "score (bits)", "q-value", "sequence", "regulated gene",
                             "occludes", "upstream of CDS (bp)", "footprint vs TSS (bp)",
                             "geometry"], rows))
        parts.append('<p class="sub">A dash under <i>regulated gene</i> means no gene start lies within '
                     '400 bp downstream of the site. <b>Read the two distances differently.</b> '
                     '<i>Upstream of CDS</i> is the distance from the site centre to the start codon of '
                     'the assigned coding sequence (CDS), taken from the genome annotation, so it is '
                     'independent of the operator. '
                     '<i>Footprint vs TSS</i> gives the site\'s two edges relative to a σ70 '
                     'transcription start site (TSS; negative = upstream) predicted inside a window centred '
                     'on the operator and then chosen for maximal overlap with it — so its narrow '
                     'spread is partly a property of that selection rule, not only of biology. '
                     '<i>occludes</i> lists the promoter elements the operator physically covers; '
                     '<i>geometry</i> is a descriptive tag for the same footprint, computed '
                     'identically for every family. Neither asserts activation or repression.</p>')
    parts.append(_sources(files, "operators/all_operators.tsv", "regulation/per_hit_regulation.json"))

    operator_section = parts[operator_start:]
    del parts[operator_start:]

    # ---- §2 inducer evidence and effector-associated genes
    parts.append('<h2 id="regulon">2 · Inducer &amp; associated genes</h2>')
    parts.append('<p class="section-lead">Independent evidence channels call an inducer or abstain. '
                 'The final call and every source opinion are recorded here. Genome annotation is then '
                 'searched for genes associated with that effector; operator support for those genes '
                 'is checked afterward.</p>')
    parts.append("<h3>Inducer inference: supporting signals</h3>")
    _ind = d.get("inducers") or {}
    _ag, _gate = _ind.get("agreement"), _ind.get("coordination_gate")
    parts.append(
        f'<p class="sub"><b>Final call: {_esc(inducer_display)}</b>'
        + (f' · source agreement {_ag:.2f}' if isinstance(_ag, (int, float)) else '')
        + (f' · coordination gate {"passed" if _gate else "not passed"}' if _gate is not None else '')
        + '</p>')
    parts.append('<p class="definition"><strong>Call</strong> is the ligand named by a channel, if any. '
                 '<strong>Role</strong> is its metal, redox or non-metal classification. '
                 '<strong>Confidence</strong> is specific to that source. <strong>Abstained</strong> '
                 'means the channel named no ligand; coordination and MetalNet can abstain on the ion '
                 'while still recording a metal site, which the caveats below then cite. '
                 '<strong>Agreement</strong> summarizes concordant channels and is not a measured '
                 'probability of correctness.</p>')
    parts.append('<p class="definition">Source names: <strong>ssn_cluster</strong> reads the assigned '
                 'sequence-similarity clade and its curated anchors; <strong>coordination</strong> '
                 'checks sequence metal-binding signatures; <strong>metalnet</strong> predicts a '
                 'metal-binding site; <strong>neighborhood</strong> checks nearby gene annotations; '
                 '<strong>ligify_db</strong> matches a published sensor entry; '
                 '<strong>ligify</strong> derives an operon-chemistry candidate.</p>')
    # The inducer_evidence figure is still generated (bundle figures/), but not embedded: it restated
    # this table as an image of raw evidence dicts, too small to read at page width.
    inducer_calls = (d.get("inducers") or {}).get("calls") or []
    if inducer_calls:
        call_rows = []
        for call in inducer_calls:
            call_rows.append([
                f'<b>{_esc(call.get("source"))}</b>', _esc(call.get("ligand") or "abstained"),
                _fmt_num(call.get("confidence")), _esc(call.get("role") or "—"),
                f'<details class="evidence-detail"><summary>View evidence</summary>'
                f'{_evidence_list(call.get("evidence"))}</details>',
            ])
        parts.append(_table(["source", "call", "confidence", "role", "recorded evidence"], call_rows))
    parts.append(_sources(files, "ligand/inducer.json", "dossier.json"))
    # The reasoning notes. These were computed, stored in dossier.json and packaged by
    # `figure_data.inducer_evidence_rows()` -- and then rendered nowhere, so every caveat the fusion
    # produced was invisible to a human reader. They carry the parts a bare headline cannot: that a
    # metal site in this family is likely structural, that a cluster is assigned to several species,
    # that MetalNet dissented, and the standing warning that agreement=1.00 with one voter means
    # UNCONTESTED rather than corroborated.
    _notes = (d.get("inducers") or {}).get("notes") or []
    if _notes:
        parts.append("<h3>How that call was reached: caveats and dissent</h3>")
        parts.append("<ul class=\"sub\">"
                     + "".join(f"<li>{_esc(n)}</li>" for n in _notes)
                     + "</ul>")
    # if the cluster names several species, show the set alongside the single headline
    _ssn = next((c for c in ((d.get("inducers") or {}).get("calls") or [])
                 if c.get("source") == "ssn_cluster"), None)
    _cands = (_ssn or {}).get("candidates") or (d.get("inducers") or {}).get("candidates") or []
    _kind = ((d.get("inducers") or {}).get("mixture_kind") or (_ssn or {}).get("mixture_kind") or "")
    if len(_cands) > 1 and _kind == "unresolved":
        # A tie is a result. Do not let the single headline above read as a decided answer.
        parts.append(
            f'<p class="run-note">'
            f'<b>UNRESOLVED between {" and ".join(_esc(c) for c in _cands)}.</b> '
            f'The available signals do not separate them, so the headline above is one of a tied set rather '
            f'than a decided answer. Read the regulon and operator sections for <i>each</i> '
            f'candidate; treating the headline as settled would be precision this run does not have. '
            f'Why nothing separates them is in the caveats below.</p>')
    elif len(_cands) > 1:
        parts.append(f'<p class="sub"><b>Candidate set:</b> '
                     f'{", ".join(_esc(c) for c in _cands)} — the headline is one of these, and the '
                     f'others remain live hypotheses rather than rejected ones.</p>')
    # Display context and support fields already recorded in the dossier.
    _mc = d.get("metal_context") or {}
    if _mc:
        _site = _mc.get("has_metalnet_site")
        _eff = _mc.get("effective_500bp")
        _gene, _gap = _mc.get("nearest_metal_gene"), _mc.get("nearest_metal_gap_bp")
        _site_txt = ("a metal site" if _site else "no metal site" if _site is False
                     else "no MetalNet opinion")
        _nb_txt = (f"{_esc(_gene)} at {_esc(_gap)} bp" if _gene else
                   "no metal-homeostasis gene within 500 bp")
        parts.append(
            f'<p class="sub"><b>Recorded metal context:</b> '
            f'MetalNet finds {_site_txt}; nearest metal machinery is {_nb_txt}. '
            f'<b>&ldquo;Effective candidate&rdquo; (site AND neighbour &le;500 bp): '
            f'{"YES" if _eff else "no"}</b>.</p>')
        _sup = _mc.get("cluster_label_support") or {}
        if _sup:
            parts.append(
                f'<p class="sub"><b>Cluster label support:</b> {_esc(_sup.get("n_independent"))} '
                f'independent anchor(s) of {_esc(_sup.get("n_anchors"))} rows, confidence '
                f'<b>{_esc(_sup.get("tier"))}</b>'
                + (", and the anchors DISAGREE about the inducer — the label is a plurality, not a "
                   "consensus" if _sup.get("conflict") else "")
                + '. Independent means distinct sequences: identical database copies of one protein '
                  'count as one supporting observation, not many.</p>')
    parts.append('<h3>Effector-associated genes, and whether an operator sits at them</h3>')
    n_annotated = len(d.get("per_hit_regulation") or [])
    parts.append('<p class="definition">Genes are listed because their name or product matches the '
                 'called effector (for a metal: its transporters, efflux pumps and enzymes), ranked by '
                 'how specific that match is, at most 12. They are chosen from annotation alone, '
                 'before any operator evidence. <strong>Support status</strong> then checks the '
                 f'{n_annotated} top-scoring sites of the genome scan (section 3): '
                 '<i>operator hit passes q≤0.05</i> means one of them regulates the gene\'s operon at '
                 'that false-discovery rate; <i>ranked operator-supported candidate</i> means a site '
                 'regulates it but does not pass q≤0.05; <i>no operator support</i> means none of those '
                 'sites does. <strong>Predicted mode</strong> is the site\'s position relative to a '
                 'predicted promoter, not activation or repression.</p>')
    cg = d.get("confirmed_regulated_genes") or []
    if cg:
        rows = []
        for g in cg[:40]:
            rows.append([f'<span class="mono">{_esc(g.get("gene"))}</span>', _esc(g.get("product")),
                         _esc(g.get("support_label") or "—"), _esc(g.get("predicted_mode") or "—")])
        parts.append(_table(["gene", "product", "support status", "predicted mode"], rows))
    else:
        parts.append('<p class="ph">No effector-associated candidate genes were recorded for this TF.</p>')
    parts.append(_sources(files, "ligand/inducer.json", "regulation/per_hit_regulation.json",
                          "dossier.json"))

    # Legacy bundles may retain an operator-derived inducer source. Display it, but state its provenance.
    _reg = next((c for c in ((d.get("inducers") or {}).get("calls") or [])
                 if c.get("source") == "regulon"), None)
    _contrib = ((_reg or {}).get("evidence") or {}).get("contributors") or []
    if _contrib:
        _rows = [[f'<span class="mono">{_esc(c.get("gene"))}</span>', _esc(c.get("effector")),
                  _esc(c.get("detail")), _esc(c.get("evidence")),
                  (f'{c.get("weight"):.2f}' if isinstance(c.get("weight"), (int, float)) else "")]
                 for c in _contrib]
        parts.append("<h3>Co-regulated genes, and the metal each implies</h3>")
        parts.append(_table(["gene", "implies", "annotation", "support route", "weight"], _rows))
        parts.append(
            f'<p class="sub">These annotations are reached through operator hits and are shown as '
            f'legacy evidence. They do not revise the current inducer call. If the sites are real, '
            f'the associated genes would suggest {_esc((_reg or {}).get("ligand"))}.</p>')
    parts.extend(operator_section)
    parts.append("<h3>Where the operators sit relative to promoters and genes</h3>")
    parts.append('<p class="definition"><strong>Left:</strong> site centre relative to the predicted '
                 'transcription start site (TSS; negative = upstream). The TSS is searched in a window '
                 'centred on each site, so a narrow peak partly reflects that search and is not '
                 'independent evidence of regulation. <strong>Right:</strong> distance of the site '
                 'centre upstream of the annotated start codon of the assigned gene. This anchor does '
                 'not depend on the site, so it is the one to compare against chance. Bars count sites; '
                 'the line is a kernel density estimate scaled to the same counts.</p>')
    parts.append(_img(figs, "operator_vs_tss", alt="histograms of operator position relative to the "
                                                    "predicted TSS and to the annotated CDS start"))
    parts.append(_sources(files, "regulation/per_hit_regulation.json"))
    # combined genome context (detail): the neighborhood of the top natural operators + a promoter zoom
    # showing which element the autoregulatory operator occludes.
    parts.append("<h3>Genome context and promoter footprint</h3>")
    parts.append('<p class="sub">The regulator\'s own locus, then the lowest-q target operons. Genes are '
                 'strand-oriented arrows and the green box is the operator. Each inset magnifies that '
                 'promoter: the pink bar is the operator footprint, the dark boxes are the −35 and −10 '
                 'elements with the spacer between them, and the dashed line is the predicted TSS. The '
                 'inset title lists the elements the footprint covers.</p>')
    parts.append(_img(figs, "operator_context", alt="genome neighbourhoods of the top operators, each "
                                                     "with a promoter inset"))
    parts.append(_sources(files, "genome/genes.tsv", "regulation/per_hit_regulation.json",
                          "regulon.tsv"))
    parts.append('<h3>Proposed regulon operons</h3>')
    parts.append('<p class="definition">The reported motif rescans the intergenic genome (p &lt; 1e-4 '
                 'per site). Each site is assigned to the gene(s) whose start lies within 400 bp '
                 'downstream, and each first gene is extended into its operon through same-strand '
                 'neighbours at most 150 bp apart. One best site is kept per operon, and operons are '
                 'ranked by q-value. <strong>q-value</strong> is the genome-wide false-discovery '
                 'estimate for the site. No q-value cutoff is applied, so rows above 0.05 are '
                 'hypotheses, not significant calls. The reported motif was built from the '
                 f'{_esc(_gl.get("n_kept") or "top")} best seed-scan hits, so those same sites score '
                 'well against it by construction and their q-values here are optimistic. For '
                 'them, the seed-scan q-value in the top-operators table is the fairer estimate.</p>')
    regulon_records = [r for r in (d.get("regulon") or []) if isinstance(r, dict) and r.get("first_gene")]
    if regulon_records:
        regulon_rows = []
        for record in regulon_records[:50]:
            site = record.get("site") or [None, None]
            regulon_rows.append([
                f'<span class="mono">{_esc(record.get("first_gene"))}</span>',
                _esc(record.get("strand")), _esc(", ".join(record.get("genes") or [])),
                f'{_esc(site[0] if len(site) else None)}–{_esc(site[1] if len(site) > 1 else None)}',
                _fmt_q(record.get("qvalue")), _esc(record.get("tier")),
            ])
        stats = d.get("regulon_stats") or {}
        found = stats.get("n_operons_found")
        shown = len(regulon_records)
        n_sig = sum(1 for r in regulon_records
                    if isinstance(r.get("qvalue"), (int, float)) and r["qvalue"] <= 0.05)
        cap_note = (f'The {shown} lowest-q of {found} mapped operons are shown (report cap).'
                    if stats.get("truncated") and found is not None else
                    f'All {shown} mapped operon(s) are shown.')
        parts.append(f'<p class="sub">{cap_note} {n_sig} of those shown pass q ≤ 0.05.</p>')
        # "locality tier", not "tier": this report uses the word for two different things -- here the
        # rescan window a hit came from, and above the confidence of a clade's anchor support.
        parts.append(_table(["first gene", "strand", "operon", "operator locus", "q-value",
                             "locality tier"], regulon_rows))
        parts.append('<p class="sub">Locality tier is the kind of window the site was found in, not a '
                     'confidence: <b>1</b> the intergenic region flanking the regulator itself (its own '
                     'promoter), <b>2</b> a divergent promoter, <b>3</b> any other intergenic region, '
                     '<b>4</b> coding DNA. It raises a site\'s prior, so a strong hit in a distant tier '
                     'can still outrank a weak one nearby. '
                     'The TSV contains the operons retained for this report.</p>')
    else:
        parts.append('<p class="ph">No proposed operons were recorded for this TF.</p>')
    parts.append(_sources(files, "regulon.tsv", "dossier.json", "motif/seed_selection.json"))

    # ---- §4 final operators
    parts.append('<h2 id="final">4 · Ranked operator records</h2>')
    parts.append('<p class="section-lead">The prediction hand-off lists the preferred natural '
                 'operator and any designed or optional structure-related records. Status describes '
                 'whether a sequence is available; it is not an experimental binding result.</p>')
    # Scores of different kinds must not share an unlabelled column: a rescan site score (bits over
    # the whole site) and a logo's mean information content (bits per column) differ ~100-fold.
    _score_units = {"rescan_site_score": "site score, bits", "mean_ic_bits": "mean IC, bits/col"}
    rows, rationales = [], {}
    for r in final_rows:
        cat = r.get("category") or ""
        cls = ("natural" if "natural" in cat else ("pending" if r.get("status") in ("pending_af3", "absent")
               else "designed"))
        tag = f'<span class="tag {cls}">{_esc(cat)}</span>'
        score = r.get("score")
        unit = _score_units.get(r.get("score_type"))
        score_txt = (f'{score:.2f}' + (f' <span class="sub">({unit})</span>' if unit else '')
                     if isinstance(score, (int, float)) else "—")
        prov = r.get("provenance") or {}
        where = (f'{_esc(prov.get("start"))}-{_esc(prov.get("end"))} ({_esc(prov.get("strand"))}), '
                 f'{_esc(prov.get("locality") or "locality not recorded")}'
                 if prov.get("start") is not None else "—")
        rows.append([tag, f'<span class="mono">{_esc(r.get("sequence") or "—")}</span>', score_txt,
                     where, _esc(r.get("status"))])
        if r.get("rationale"):
            rationales.setdefault(cat, r["rationale"])
    parts.append(_table(["source", "sequence", "score", "genome site", "status"], rows))
    if rationales:
        parts.append('<p class="sub">' + ' '.join(f'<b>{_esc(cat)}:</b> {_esc(why)}'
                                                  for cat, why in rationales.items()) + '</p>')
    parts.append(_sources(files, "operators_ranked.json", "operators/FINAL_operators.tsv",
                          "dossier.json", "af3_jobs.json"))

    # ---- §5 run record: exact normalized inputs, effective parameters and evidence-to-file map
    parts.append('<h2 id="run">5 · Run record &amp; traceability</h2>')
    parts.append('<p class="section-lead">This section is the audit trail: what entered the run, what '
                 'each stage used, what it produced, and where the complete record lives inside this '
                 'TF directory.</p>')
    query = run.get("query") or {}
    genome = run.get("genome") or {}
    evidence = data.get("evidence_index") or run.get("evidence_index") or []
    config = run.get("effective_config") or {}
    parts.append('<div class="audit-grid">'
                 f'<div class="audit-card"><b>Normalized query</b><strong>{_esc(query.get("length") or len(seq))} aa</strong>'
                 f'<small class="hash">SHA-256 {_esc(query.get("sha256") or "not recorded")}</small></div>'
                 f'<div class="audit-card"><b>Scanned genome</b><strong>{_esc(genome.get("length") or (d.get("rescan") or {}).get("region_len") or 0):} bp</strong>'
                 f'<small>{_esc(genome.get("gene_count") or 0)} parsed gene features</small></div>'
                 f'<div class="audit-card"><b>Bundle coverage</b><strong>{len(files)} files</strong>'
                 f'<small>{len(d)} dossier fields · {len(evidence)} evidence groups</small></div>'
                 '</div>')

    requested = run.get("requested") or {}
    parameter_rows = []
    for key in sorted(set(requested) | set(config)):
        requested_value = requested.get(key, "—")
        effective_value = config.get(key, "—")
        parameter_rows.append([f'<span class="mono">{_esc(key)}</span>', _esc(requested_value),
                               _esc(effective_value)])
    if parameter_rows:
        parts.append('<h3>Requested inputs and effective configuration</h3>')
        parts.append(_table(["field", "requested", "effective / reachable"], parameter_rows))

    homolog_collection = run.get("homolog_collection") or {}
    attempts = homolog_collection.get("attempts") or []
    if attempts:
        attempt_rows = []
        for attempt in attempts:
            record = attempt.get("record") or {}
            assignment = attempt.get("cluster_assignment") or {}
            count = len(record.get("regions") or []) if record else attempt.get("n_regions")
            detail = assignment.get("cluster_id") or record.get("key") or "—"
            attempt_rows.append([
                _esc(attempt.get("route")), "yes" if attempt.get("selected") else "no",
                _esc(count or 0), _esc(detail),
                _esc(assignment.get("method") or record.get("sampling") or "—"),
            ])
        msa = homolog_collection.get("msa") or {}
        # Bound before the branch: the shortfall note below reads it whether or not the MSA ran, and a
        # binding inside `if` made every bundle whose MSA never ran die at the output stage.
        source_record = msa.get("source_record") or {}
        if msa.get("source"):
            attempt_rows.append([
                "MSA expansion", "yes" if msa.get("msa_added") else "no",
                _esc(msa.get("msa_added") or 0), _esc(msa.get("source")),
                f"{_esc(source_record.get('n_orthologs') or 0)} orthologs / "
                f"{_esc(source_record.get('n_mapped') or 0)} mapped",
            ])
        parts.append('<h3>Homolog-source ledger</h3>')
        parts.append(_table(["route", "selected", "regions added", "source / cluster", "method / mapping"],
                            attempt_rows))
        # A conservation set under the floor is the difference between a PWM resting on many sequences
        # and one resting on a handful. Saying the mapped count without saying it was short leaves the
        # reader to know the floor by heart.
        _sf = (d.get("homolog_regions") or {}).get("shortfall") or source_record.get("shortfall")
        if _sf:
            parts.append(
                f'<p class="sub"><b>Thin conservation set:</b> {_esc(_sf.get("n_mapped"))} mapped '
                f'ortholog(s) of {_esc(_sf.get("n_mappable"))} mappable, below the floor of '
                f'{_esc(_sf.get("floor"))} — {_esc(_sf.get("reason"))}. The homolog-derived logo and '
                f'any conservation weighting rest on that set, so treat them as weak here rather than '
                f'as evidence of a conserved site.</p>')

    stages = run.get("stages") or []
    if stages:
        stage_rows = []
        for stage in stages:
            status = stage.get("status") or "unresolved"
            stage_rows.append([
                f'<b>{_esc(stage.get("stage"))}</b>',
                f'<span class="status {_status_class(status)}">{_esc(status)}</span>',
                _esc(stage.get("used")), _esc(stage.get("result")),
                _artifact(stage.get("artifact"), stage.get("artifact"), files),
            ])
        parts.append('<h3>Pipeline stage ledger</h3>')
        parts.append(_table(["stage", "status", "information used", "recorded result", "artifact"], stage_rows))

    seed_selection = d.get("seed_selection") or {}
    construction = seed_selection.get("construction") or {}
    selection_params = seed_selection.get("parameters") or {}
    if construction or selection_params:
        exact_rows = []
        for group, values in (("seed construction", construction), ("seed selection", selection_params)):
            for key, value in values.items():
                if value is not None:
                    exact_rows.append([_esc(group), f'<span class="mono">{_esc(key)}</span>', _esc(value)])
        parts.append('<h3>Recorded motif and rescan parameters</h3>')
        parts.append(_table(["scope", "parameter", "value"], exact_rows))

    if evidence:
        evidence_rows = []
        for item in evidence:
            status = item.get("status") or "unresolved"
            anchor = item.get("report_anchor")
            label = (f'<a href="#{quote(str(anchor))}">{_esc(item.get("evidence"))}</a>'
                     if anchor else _esc(item.get("evidence")))
            evidence_rows.append([
                label, _esc(item.get("records")),
                f'<span class="status {_status_class(status)}">{_esc(status)}</span>',
                _artifact(item.get("artifact"), item.get("artifact"), files),
            ])
        parts.append('<h3>Evidence-to-artifact map</h3>')
        parts.append(_table(["information", "records / size", "capture", "complete artifact"], evidence_rows))

    caveats = d.get("caveats") or []
    if caveats:
        parts.append('<h3>Recorded interpretation limits</h3><ul class="sub">'
                     + "".join(f'<li>{_esc(note)}</li>' for note in caveats) + '</ul>')
    parts.append('<p class="run-note"><b>Bundle boundary.</b> The exact normalized query, scanned genome '
                 'sequence, parsed gene model, retained promoter sequences, optional source alignment, '
                 'effective settings, evidence tables, figures and decisions are stored here. External '
                 'binaries and the complete packaged reference database are identified by software/data '
                 'fingerprints; they are not duplicated into every TF folder. No report asset or recorded '
                 'structure/alignment path points to a cache outside this directory.</p>')

    # ---- §6 complete local file index
    downloads = [("dossier.json", "Full recorded result (JSON)"),
                 ("bundle_manifest.json", "File checksums and software provenance"),
                 ("report_environment.json", "Export environment and package versions"),
                 ("README.md", "Bundle guide and reproducibility limits"),
                 ("run/run_record.json", "Inputs, effective parameters and stage ledger"),
                 ("run/evidence_index.json", "Evidence-to-artifact index"),
                 ("input/query.fasta", "Normalized TF query used by the pipeline"),
                 ("input/run_parameters.json", "Requested and effective run parameters"),
                 ("genome/genome.fna", "Exact genome sequence scanned"),
                 ("genome/genes.tsv", "Parsed gene model used downstream"),
                 ("binding_sites.tsv", "Complete binding-site table"),
                 ("regulon.tsv", "Regulon table (the operons shown in this report)"),
                 ("homologs/homolog_promoters.fasta", "Recorded upstream sequences"),
                 ("homologs/promoter_records.tsv", "Promoter source identities and sequence hashes"),
                 ("homologs/source_record.json", "Homolog and MSA collection provenance"),
                 ("motif/three_logos.json", "Recorded logo matrices and alignments"),
                 ("operators/FINAL_operators.tsv", "Complete final operator table")]
    parts.append('<h2 id="files">6 · Files &amp; provenance</h2>')
    parts.append('<p class="section-lead">Open the complete recorded inputs and results. Links are '
                 'relative to this TF bundle, so they work when the folder is copied to another '
                 'computer. The manifest records file hashes for checking that copy.</p>')
    parts.append('<div class="file-list">')
    parts.extend(f'<a href="{path}">{label}</a>' for path, label in downloads if path in files)
    parts.append('</div><p class="provenance-note">Tables may be abbreviated; download the complete recorded exports below. '
                 'The bundle guide lists inputs and provenance that are not captured.</p>')
    parts.append(f'<details><summary>Browse all {len(files - {"REPORT.html"})} recorded files</summary><div class="file-list">')
    parts.extend(f'<a href="{quote(path, safe="/")}">{_esc(path)}</a>' for path in sorted(files)
                 if path != "REPORT.html")
    parts.append('</div></details>')

    return (f"<!doctype html><html lang='en'><head><meta charset='utf-8'>"
            f"<meta name='viewport' content='width=device-width,initial-scale=1'>"
            f"<meta name='description' content='TF operator prediction report for {_esc(tf)}'>"
            f"<title>{_esc(tf)} — TF operator report</title><style>{_CSS}</style></head>"
            f"<body><a class='skip' href='#decision'>Skip to result</a><main>{''.join(parts)}</main></body></html>")


# --------------------------------------------------------------------------- entry points
def _render_contents(job_dir: Path, manifest: dict | None, *, pending_manifest: bool = False) -> Path:
    """Render into an unpublished directory; callers own validation and publication."""
    manifest_path = job_dir / "bundle_manifest.json"
    fig_dir = job_dir / "figures"
    data = collect(job_dir)
    figs = _gen_figures(data, fig_dir)
    # Include newly generated figures in the file browser too.
    data["bundle_files"] = {p.relative_to(job_dir).as_posix() for p in job_dir.rglob("*") if p.is_file()}
    if pending_manifest:
        data["bundle_files"].add("bundle_manifest.json")
    out = job_dir / "REPORT.html"
    out.write_text(build_html(data, figs), encoding="utf-8")
    if manifest is not None:
        manifest["schema_version"] = 2
        manifest["artifacts"] = inventory(job_dir)
        manifest["required_artifacts"] = sorted(manifest["artifacts"])
        manifest["report_generated"] = True
        # Presentation regeneration does not change the original computation's runtime fingerprint.
        from predictor.provenance import runtime_provenance
        manifest["report_runtime"] = runtime_provenance()
        manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    n_ok = sum(1 for v in figs.values() if v)
    print(f"REPORT.html -> {out}  ({n_ok}/{len(figs)} figures)")
    return out


def render(job_dir, *, pending_manifest: bool = False) -> Path:
    """Regenerate a report without modifying the last valid bundle on failure.

    Newly created output staging directories have no manifest and are rendered in
    place. Existing bundles are validated, copied to staging and replaced only after
    the new report passes its file contract. Their original folders are preserved.
    """
    job_dir = Path(job_dir).resolve()
    manifest_path = job_dir / "bundle_manifest.json"
    if not manifest_path.exists():
        return _render_contents(job_dir, None, pending_manifest=pending_manifest)
    manifest_bytes = manifest_path.read_bytes()
    try:
        manifest = json.loads(manifest_bytes)
    except ValueError as exc:
        raise ValueError("bundle_manifest.json must contain a valid object") from exc
    if not isinstance(manifest, dict):
        raise ValueError("bundle_manifest.json must contain a valid object")
    # A verified legacy bundle can contain presentation markup that the current
    # renderer intentionally forbids. Check its complete inventory and hashes,
    # but defer the HTML policy audit until the replacement report exists.
    issues = validate_files(
        job_dir,
        require_report=False,
        _audit_report_html=False,
    )
    if issues:
        raise ValueError("cannot rerender an inconsistent bundle: " + "; ".join(issues))
    original_files = inventory(job_dir)  # Also rejects links before copytree can follow them.
    staging = Path(tempfile.mkdtemp(prefix=".report-stage-", dir=job_dir.parent))
    backup = job_dir.with_name(f".{job_dir.name}.report-previous-{uuid.uuid4().hex[:10]}")
    try:
        shutil.copytree(job_dir, staging, dirs_exist_ok=True)
        if inventory(staging) != original_files:
            raise ValueError("bundle changed while it was being copied for report regeneration")
        figures = staging / "figures"
        if figures.exists():
            shutil.rmtree(figures)
        from predictor.report.environment import snapshot
        previous_environment = _load(staging / "report_environment.json")
        if previous_environment is not None:
            manifest.setdefault("report_environment_history", []).append(previous_environment)
        (staging / "report_environment.json").write_text(json.dumps(snapshot(), indent=2), encoding="utf-8")
        from predictor.report.outputs import bundle_readme
        (staging / "README.md").write_text(bundle_readme(), encoding="utf-8")
        _render_contents(staging, manifest)
        issues = validate_files(staging)
        if issues:
            raise ValueError("regenerated report failed validation: " + "; ".join(issues))
        if inventory(job_dir) != original_files or manifest_path.read_bytes() != manifest_bytes:
            raise ValueError("source bundle changed during report regeneration")
        job_dir.replace(backup)
        try:
            staging.replace(job_dir)
        except BaseException:
            backup.replace(job_dir)
            raise
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    print(f"Previous report bundle preserved: {backup}")
    return job_dir / "REPORT.html"



def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--tf", required=True, help="job name under results/jobs/")
    ap.add_argument("--root", default=str(JOBS))
    a = ap.parse_args(argv)
    from predictor.report.outputs import job_path
    try:
        render(job_path(a.tf, root=Path(a.root)))
    except (OSError, ValueError, RuntimeError) as exc:
        ap.exit(1, f"Report regeneration failed: {exc}\n")


if __name__ == "__main__":
    main()
