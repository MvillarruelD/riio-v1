"""One-page web interface for TF prediction and genome/proteome census.

The UI is intentionally a thin client over the same public API and scanners used by ``tfop``. It owns
presentation and input guidance only; scientific behavior remains in ``predictor.api``. Every long run
writes the same folder layout as the command line, so a run can be started in either place and opened
in the other (see "Saved runs" and "Reproduce in a terminal").
"""
from __future__ import annotations

import contextlib
import html
import io
import queue
import threading
import time
from pathlib import Path

import streamlit as st

from predictor import api, resources
from predictor.annotate.genome_scan import (PIPELINE_FAMILIES, POTENTIAL_METALLOREGULATOR_FAMILIES,
                                            census_families)
from predictor.config import RunConfig
from predictor.gui_support import (CLASSES, PIPELINE_RUN_LOCK, available_job_name, bundle_headline, cli_command,
                                   context_archive, display_paths, is_done, is_failed, job_dir, list_saved_bundles,
                                   list_saved_runs, read_manifest, pasted_protein_input, recorded_policy,
                                   report_preview, reserve_output_dir, result_dict, run_command, run_kind,
                                   run_rows, safe_name)
from predictor.input_validation import read_protein


APP_TITLE = "TF operator predictor"  # Product-title placeholder: replace when the final project name is chosen.
GITHUB_URL = "https://github.com/MvillarruelD/riio-v1"

st.set_page_config(page_title=APP_TITLE, page_icon="◈", layout="wide", initial_sidebar_state="collapsed")

# Palette and type. Light only, on purpose: one deliberate theme (also set in .streamlit/config.toml and the
# launcher) instead of half a dark mode. Fonts are system stacks so the page works offline; the serif carries
# headings, the monospace carries identifiers and sequences, and figures are tabular.
_CSS = """
<style>
:root{
  --paper:#F4F1EA;--raised:#FBFAF6;--sunk:#EAE6DB;--ink:#1D1C19;--ink-2:#47443C;--muted:#675F55;
  --line:#D6D0C1;--line-2:#B9B2A0;--accent:#B4431F;--accent-ink:#8C3015;
  --metal:#2D5B8C;--redox:#A8690F;--organic:#5C7843;--unresolved:#8A8478;
  --serif:"Source Serif 4","Iowan Old Style","Palatino Linotype",Cambria,Georgia,"DejaVu Serif",serif;
  --mono:ui-monospace,"Cascadia Mono","SF Mono",Consolas,"DejaVu Sans Mono","Liberation Mono",monospace;
}
[data-testid="stAppViewContainer"]{background:var(--paper);color:var(--ink)}
[data-testid="stHeader"]{background:transparent}
[data-testid="stToolbar"]{visibility:hidden;height:0}
.block-container{max-width:1180px;padding-top:1.6rem;padding-bottom:4rem}
h1,h2,h3{font-family:var(--serif)!important;color:var(--ink);letter-spacing:-.012em;font-weight:600!important;text-wrap:balance}
h2{font-size:1.5rem!important;margin:0 0 .2rem!important;padding:0!important}
h3{font-size:1.12rem!important}
p,li,label{line-height:1.55}p{max-width:68ch}
a{color:var(--accent-ink)}
:focus-visible{outline:2px solid var(--accent)!important;outline-offset:2px}

/* masthead */
.mast{display:grid;grid-template-columns:minmax(0,1fr) auto;gap:2rem;align-items:end;border-bottom:2px solid var(--ink);padding-bottom:1rem}
.mast-kicker{font:500 .78rem var(--mono);color:var(--muted)}
.mast-title{font:600 clamp(1.9rem,3.2vw,2.6rem)/1.05 var(--serif);letter-spacing:-.02em;margin:.3rem 0 .45rem}
.mast-lead{margin:0;color:var(--ink-2);font-size:1rem;max-width:58ch}
.mast-meta{font:.76rem/1.7 var(--mono);color:var(--muted);text-align:right;white-space:nowrap}
.mast-meta b{color:var(--ink);font-weight:600}.mast-meta a{color:var(--ink)}
.skip-link{position:fixed;left:1rem;top:-4rem;z-index:100;background:var(--ink);color:var(--paper)!important;padding:.6rem .8rem;border-radius:3px}
.skip-link:focus{top:1rem}

/* section headings */
.sec{display:flex;align-items:baseline;gap:.85rem;border-top:1px solid var(--line-2);padding-top:.55rem;margin:2.4rem 0 .9rem}
.sec [data-testid="stHeaderActionElements"]{display:none}
.sec i{font:500 .78rem var(--mono);font-style:normal;color:var(--accent-ink)}
.sec h2{margin:0!important}.sec .sec-note{color:var(--muted);font-size:.86rem;margin-left:auto;text-align:right}
.lede{color:var(--ink-2);font-size:.93rem;margin:.1rem 0 1rem;max-width:64ch}
.fine{color:var(--muted);font-size:.8rem}.mono{font-family:var(--mono);font-size:.84em}

/* steps */
.steps{list-style:none;display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:.9rem;margin:1.1rem 0 1.5rem;padding:0}
.steps li{border-top:3px solid var(--line);padding:.5rem .1rem 0;color:var(--muted);font-size:.8rem;margin:0}
.steps li b{display:block;color:var(--ink);font:600 1rem var(--serif);margin:.05rem 0 .05rem}
.steps li.done{border-color:var(--ink)}.steps li.now{border-color:var(--accent)}
.steps li.now span{color:var(--accent-ink)}.steps li span{font:500 .72rem var(--mono)}

/* class marker, sheet, note */
.k{display:inline-flex;align-items:center;gap:.38rem;white-space:nowrap}
.k i{width:.62rem;height:.62rem;border-radius:1px;background:var(--unresolved);flex:none}
.k-metal i{background:var(--metal)}.k-redox i{background:var(--redox)}.k-organic i{background:var(--organic)}
.k-unresolved i{background:transparent;border:1.5px solid var(--unresolved)}
.sheet{display:grid;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:0;margin:1rem 0 1.1rem;border-top:2px solid var(--ink)}
.sheet>div{padding:.75rem 1.1rem .8rem 0;border-bottom:1px solid var(--line)}
.sheet dt{font-size:.74rem;color:var(--muted);margin:0 0 .2rem}
.sheet dd{margin:0;font:600 1.5rem/1.2 var(--serif);overflow-wrap:anywhere}
.sheet dd.mono{font:500 .92rem/1.4 var(--mono);letter-spacing:.005em}
.sheet small{display:block;margin-top:.3rem;color:var(--muted);font-size:.78rem;line-height:1.4}
.note{border-left:3px solid var(--accent);background:var(--sunk);padding:.65rem .95rem;margin:.9rem 0 1.2rem;font-size:.88rem;color:var(--ink-2);max-width:74ch}
.paths{font:.78rem/1.6 var(--mono);color:var(--muted);overflow-wrap:anywhere}.paths b{color:var(--ink-2);font-weight:600}

/* ledger table */
.ledger-wrap{overflow-x:auto;border-top:2px solid var(--ink);border-bottom:1px solid var(--line-2);margin:.4rem 0 1rem}
.ledger{width:100%;border-collapse:collapse;font-size:.86rem}
.ledger th{text-align:left;font:500 .74rem var(--mono);color:var(--muted);padding:.55rem .7rem .45rem 0;border-bottom:1px solid var(--line-2);white-space:nowrap}
.ledger td{padding:.5rem .7rem .5rem 0;border-bottom:1px solid var(--line);vertical-align:baseline}
.ledger tr:last-child td{border-bottom:0}.ledger tbody tr:hover td{background:rgba(234,230,219,.55)}
.ledger .num{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}
.ledger td small{display:block;color:var(--muted);font-size:.74rem;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;max-width:24ch}
.ledger td small .k{font-size:.76rem}
.ledger .id{font-family:var(--mono);font-size:.82rem;font-weight:600}.ledger tr.dim td{color:var(--muted)}
.ledger .op{font-family:var(--mono);font-size:.78rem;letter-spacing:.01em}
.meter{display:inline-block;width:3.2rem;height:.38rem;background:var(--sunk);vertical-align:middle;margin-right:.45rem;border-radius:1px}
.meter i{display:block;height:100%;background:var(--ink)}

/* widgets */
[data-testid="stForm"]{background:var(--raised);border:1px solid var(--line-2);border-radius:4px;padding:1.1rem 1.25rem 1.25rem}
[data-baseweb="input"],[data-baseweb="textarea"],[data-baseweb="select"]>div{background:var(--raised)!important;border-color:var(--line-2)!important;border-radius:3px!important}
[data-baseweb="input"]:focus-within,[data-baseweb="textarea"]:focus-within{border-color:var(--ink)!important}
[data-baseweb="tag"]{background:var(--sunk)!important;border:1px solid var(--line-2)!important;border-radius:2px!important;color:var(--ink)!important}
[data-baseweb="tag"] *{color:var(--ink)!important}
[data-testid="stFileUploaderDropzone"]{background:var(--raised);border:1px dashed var(--line-2);border-radius:3px}
[data-testid="stButton"] button,[data-testid="stFormSubmitButton"] button,[data-testid="stDownloadButton"] button{
  border-radius:3px;min-height:2.6rem;font-weight:600;touch-action:manipulation;border:1px solid var(--ink);
  background:transparent;color:var(--ink);transition:background-color .15s,color .15s,border-color .15s}
[data-testid="stButton"] button:hover,[data-testid="stDownloadButton"] button:hover,
[data-testid="stFormSubmitButton"] button:hover{background:var(--sunk);border-color:var(--ink);color:var(--ink)}
button[kind="primary"],button[kind="primaryFormSubmit"],[data-testid="stBaseButton-primary"],
[data-testid="stBaseButton-primaryFormSubmit"]{background:var(--ink)!important;color:var(--paper)!important;border-color:var(--ink)!important}
button[kind="primary"]:hover,button[kind="primaryFormSubmit"]:hover,[data-testid="stBaseButton-primary"]:hover,
[data-testid="stBaseButton-primaryFormSubmit"]:hover{background:var(--accent)!important;border-color:var(--accent)!important;color:#fff!important}
button:disabled,button[disabled]{background:var(--sunk)!important;border-color:var(--line)!important;color:#8b857a!important}
/* Segmented control: Streamlit 1.65 renders button[role=radio] inside stButtonGroup; 1.63 used stSegmentedControl. */
[data-testid="stSegmentedControl"],[data-testid="stButtonGroup"]{margin:.1rem 0 .5rem}
[data-testid="stSegmentedControl"] button,[data-testid="stButtonGroup"] button[role="radio"]{
  min-height:2.8rem;font-weight:600;background:transparent!important;color:var(--ink)!important;
  border:1px solid var(--line-2)!important;border-radius:3px!important}
[data-testid="stSegmentedControl"] button:hover,[data-testid="stButtonGroup"] button[role="radio"]:hover{background:var(--sunk)!important}
[data-testid="stSegmentedControl"] button[aria-pressed="true"],[data-testid="stSegmentedControl"] [role="radio"][aria-checked="true"],
[data-testid="stButtonGroup"] button[role="radio"][aria-checked="true"],[data-testid="stBaseButton-segmented_controlActive"]{
  background:var(--ink)!important;color:var(--paper)!important;border-color:var(--ink)!important}
[data-testid="stSegmentedControl"] button[aria-pressed="true"] *,[data-testid="stSegmentedControl"] [role="radio"][aria-checked="true"] *,
[data-testid="stButtonGroup"] button[role="radio"][aria-checked="true"] *,[data-testid="stBaseButton-segmented_controlActive"] *{color:var(--paper)!important}
[data-testid="stBaseButton-segmented_control"],[data-testid="stBaseButton-segmented_controlInactive"]{background:transparent!important;color:var(--ink)!important}
[data-testid="stExpander"]{border:1px solid var(--line-2);border-radius:3px;background:transparent}
[data-testid="stExpander"] summary:hover{background:var(--sunk)}
[data-testid="stDataFrame"]{border:1px solid var(--line-2);border-radius:3px;overflow:hidden}
[data-testid="stCode"] pre,[data-testid="stCodeBlock"] pre{background:var(--sunk)!important;border-radius:3px}
[data-testid="stAlert"]{border-radius:3px}
[data-testid="stStatusWidget"]{border-radius:3px}
@media(max-width:760px){
  .block-container{padding:1rem .9rem 3rem}.mast{grid-template-columns:1fr;gap:.8rem}.mast-meta{text-align:left;white-space:normal}
  .steps{grid-template-columns:1fr;gap:.2rem}.sec .sec-note{display:none}.sheet dd{font-size:1.3rem}
}
@media(prefers-reduced-motion:reduce){*,*::before,*::after{scroll-behavior:auto!important;animation:none!important;transition:none!important}}
</style>
"""
st.markdown(_CSS, unsafe_allow_html=True)

_PRODUCTION = "Production · BITACORA + SSN-clade profiles"
_QUICK = "Quick census · Pfam (not the production route)"
_COLORS = {"metal": "#2D5B8C", "redox": "#A8690F", "organic": "#5C7843", "unresolved": "#B3AD9F"}
_VEGA_THEME = {
    "view": {"stroke": None},
    "font": "Source Sans Pro, Segoe UI, system-ui, sans-serif",
    "axis": {"labelColor": "#47443C", "titleColor": "#675F55", "domainColor": "#B9B2A0", "tickColor": "#B9B2A0",
             "gridColor": "#E2DDD0", "labelFontSize": 12, "titleFontSize": 11, "titleFontWeight": "normal"},
    "legend": {"labelColor": "#47443C", "titleColor": "#675F55", "symbolType": "square", "orient": "top",
               "labelFontSize": 12, "titleFontSize": 11, "titleFontWeight": "normal"},
}


# ------------------------------------------------------------------------------------------- small helpers
def _e(value) -> str:
    # A line break inside an HTML block ends it, and the rest would be parsed as Markdown.
    text = "" if value is None else str(value).replace("\r", " ").replace("\n", " ")
    return html.escape(text, quote=True)


def _safe_name(value: str, fallback: str = "query") -> str:
    return safe_name(value, fallback)


def _available_output_path(stem: str) -> Path:
    """Return a scan output folder that never overwrites an earlier completed census."""
    return reserve_output_dir(resources.OUTPUT_DIR, stem)


def _available_job_name(stem: str) -> str:
    """Return a stable, filesystem-safe job name without making users resolve output collisions."""
    return available_job_name(resources.output_path("jobs"), stem)


def _inline_report_assets(text: str, bundle: Path) -> str:
    return report_preview(text, bundle)


def _klass(name: str) -> str:
    name = name if name in CLASSES else "unresolved"
    return f'<span class="k k-{name}"><i></i>{name}</span>'


def _section(number: str, title: str, note: str = "") -> None:
    tail = f'<small class="sec-note">{_e(note)}</small>' if note else ""
    st.markdown(f'<div class="sec"><i>{_e(number)}</i><h2>{_e(title)}</h2>{tail}</div>', unsafe_allow_html=True)


def _paths_block(path: Path) -> None:
    rows = "".join(f"<div><b>{_e(label)}</b> {_e(value)}</div>" for label, value in display_paths(path))
    st.markdown(f'<div class="paths">{rows}</div>', unsafe_allow_html=True)


@st.cache_data(ttl=600, show_spinner=False)
def _discovery_engine() -> tuple[bool, str]:
    """Whether production discovery (BITACORA) can run here; probing WSL is slow, so it is cached."""
    try:
        from predictor.discovery import bitacora
        return bitacora.engine_available()
    except Exception as exc:  # an unreadable engine is "not available", never a crashed page
        return False, f"{type(exc).__name__}: {exc}"


def _full_download(root: Path, kind: str) -> None:
    """Build potentially large archives only after an explicit user action."""
    key = f"context_archive_{kind}_{root.resolve()}"
    if st.button(f"Prepare full {kind} download", key=f"prepare_{key}"):
        try:
            with st.spinner("Preparing saved report files and assets…"):
                st.session_state[key] = context_archive(root)
        except (OSError, ValueError) as exc:
            st.error(f"Could not prepare the download: {exc}")
            st.session_state.pop(key, None)
    if key in st.session_state:
        data, skipped = st.session_state[key]
        st.download_button(f"Download full {kind} · ZIP", data,
                           file_name=f"{_safe_name(root.name)}_{kind.replace(' ', '_')}.zip",
                           mime="application/zip", key=f"download_{key}")
        st.caption("Extract the ZIP before opening its HTML reports. Prepare again to include newly saved files.")
        if skipped:
            st.warning("Linked or temporary paths were omitted: " + ", ".join(skipped))


class _LogWriter(io.TextIOBase):
    def __init__(self, sink: queue.Queue):
        self.sink = sink

    def write(self, text):
        if text:
            self.sink.put(str(text))
        return len(text or "")

    def flush(self):
        return None


def _run_with_log(fn):
    """Run a blocking pipeline call while keeping its elapsed time and latest progress visible."""
    sink: queue.Queue = queue.Queue()
    holder: dict = {}

    def target():
        try:
            # redirect_stdout is process-global, not thread-local. Serializing GUI runs prevents two
            # browser sessions from capturing each other's progress or corrupting their log displays.
            if not PIPELINE_RUN_LOCK.acquire(blocking=False):
                sink.put("Another analysis is running. This run will start when it finishes…\n")
                PIPELINE_RUN_LOCK.acquire()
            try:
                with contextlib.redirect_stdout(_LogWriter(sink)), contextlib.redirect_stderr(_LogWriter(sink)):
                    holder["result"] = fn()
            finally:
                PIPELINE_RUN_LOCK.release()
        except BaseException as exc:  # surface the worker exception in the UI thread
            holder["error"] = exc

    worker = threading.Thread(target=target, daemon=True)
    worker.start()
    started = time.monotonic()
    clock = st.empty()
    placeholder = st.empty()
    chunks: list[str] = []
    while worker.is_alive() or not sink.empty():
        while not sink.empty():
            chunks.append(sink.get_nowait())
        latest = "".join(chunks)[-9000:].strip()
        seconds = int(time.monotonic() - started)
        clock.caption(f"Running for {seconds // 60}:{seconds % 60:02d} · keep this tab open; if it closes, the "
                      "folder keeps what finished and the run can be reopened under Saved runs.")
        if latest:
            placeholder.code(latest, language=None)
        time.sleep(.15)
    worker.join()
    clock.empty()
    if "error" in holder:
        raise holder["error"]
    return holder.get("result"), "".join(chunks)


def _resolve_sequence(method: str, pasted: str, uploaded, accession: str, *, offline: bool) -> str:
    if method == "Paste protein sequence":
        source = pasted
    elif method == "Upload FASTA":
        source = uploaded.getvalue().decode("utf-8", errors="replace") if uploaded else ""
    else:
        acc = accession.strip()
        if not acc:
            raise ValueError("Enter a RefSeq/GenBank protein accession, for example NP_417751.1 (E. coli ZntR).")
        if offline:
            raise ValueError("Protein-accession lookup needs the network. Turn off offline mode or paste the sequence.")
        from predictor.annotate.protein_seqs import fetch_protein_sequences
        found = fetch_protein_sequences([acc], allow_ncbi=True, verbose=True)
        source = found.get(acc.split(".")[0], "")
        if not source:
            raise ValueError(f"No protein sequence was found for {acc!r}. Check the accession or paste FASTA.")
    read_protein(source)  # actionable validation before the long-running pipeline starts
    return source


# ------------------------------------------------------------------------------------------- masthead, nav
def _masthead(db: dict) -> None:
    st.markdown('<a class="skip-link" href="#analysis-workspace">Skip to analysis</a>', unsafe_allow_html=True)
    st.markdown(
        '<header class="mast"><div><div class="mast-kicker">RIIO · bacterial regulatory genomics</div>'
        f'<div class="mast-title">{_e(APP_TITLE)}</div>'
        '<p class="mast-lead">From a genome, a protein list or a single regulator to ranked operators, '
        'a possible inducer and a regulon.</p></div>'
        f'<div class="mast-meta"><div>reference data <b>v{_e(db.get("release"))}</b></div>'
        f'<div>runs on this computer</div><div><a href="{_e(GITHUB_URL)}" target="_blank" rel="noopener">source ↗</a></div>'
        '</div></header>', unsafe_allow_html=True)


def _workflow_selector() -> str:
    """Persistent primary navigation. Genome scanning is the first-run default."""
    st.markdown('<div id="analysis-workspace" style="height:1.1rem"></div>', unsafe_allow_html=True)
    choice = st.segmented_control(
        "Analysis type", ("Genome scan", "Protein list", "Single TF"), default="Genome scan",
        key="workflow_choice", label_visibility="collapsed", width="stretch",
    )
    if choice == "Protein list":
        st.markdown('<p class="lede">Predict every protein you supply, as given. Nothing is filtered out by discovery.</p>',
                    unsafe_allow_html=True)
        return "batch"
    if choice == "Single TF":
        st.markdown('<p class="lede">One accession or sequence to operators, a possible inducer and the regulon.</p>',
                    unsafe_allow_html=True)
        return "predict"
    st.markdown('<p class="lede">Find transcription-factor candidates in a genome, choose which to take forward, '
                'then build a full report for each.</p>', unsafe_allow_html=True)
    return "scan"


def _scope_panel() -> None:
    with st.expander("Scope: which families get a full report", expanded=False):
        st.markdown(f"The census screens {len(census_families())} bacterial TF families; complete "
                    "operator–inducer–regulon reports "
                    f"currently support {len(PIPELINE_FAMILIES)}.")
        st.write(" · ".join(PIPELINE_FAMILIES))
        st.caption("Blue marks families that include known metalloregulators. It is a family-level cue, "
                   "not evidence that an individual protein senses a metal.")


def _stepper(active: int, labels=("Find candidates", "Choose regulators", "Build reports")) -> None:
    items = []
    for index, label in enumerate(labels, 1):
        state = "done" if index < active else "now" if index == active else ""
        note = "done" if index < active else "current" if index == active else ""
        items.append(f'<li class="{state}"><span>{index:02d}{" · " + note if note else ""}</span><b>{_e(label)}</b></li>')
    st.markdown(f'<ol class="steps">{"".join(items)}</ol>', unsafe_allow_html=True)


# ------------------------------------------------------------------------------------------- saved work
def _under_allowed_roots(path: Path) -> bool:
    allowed = [resources.OUTPUT_DIR.resolve(), Path.cwd().resolve()]
    return any(path == root or root in path.parents for root in allowed)


def _saved_runs_picker(state_key: str) -> None:
    """Reopen a multi-protein run created here or by ``tfop genome`` / ``tfop batch``."""
    roots = [resources.OUTPUT_DIR, Path.cwd() / "runs"]
    kind = "batch" if state_key == "last_shared_batch" else "genome"
    runs, seen = [], set()
    for root in roots:
        for entry in list_saved_runs(root, kind=kind):
            if entry["path"].resolve() not in seen:
                seen.add(entry["path"].resolve())
                runs.append(entry)
    by_path = {str(r["path"]): r for r in runs}
    with st.expander(f"Saved runs ({len(runs)})" if runs else "Open a saved run", expanded=False):
        if runs:
            choice = st.selectbox(
                "Run", list(by_path), index=None, placeholder="Choose a saved run…", key=f"saved_pick_{state_key}",
                format_func=lambda p: f'{by_path[p]["name"]} · {by_path[p]["done"]}/{by_path[p]["total"]} predicted · '
                                      f'{by_path[p]["modified"]:%Y-%m-%d %H:%M}')
            if choice is not None and st.button("Open this run", key=f"saved_open_{state_key}"):
                st.session_state[state_key] = choice
                st.rerun()
        typed = st.text_input("Or a run folder path", key=f"saved_path_{state_key}",
                              placeholder=r"runs/ecoli   or   C:\data\runs\ecoli",
                              help="A folder holding manifest.csv, made by `tfop genome`, `tfop batch` or this page. "
                                   "Only folders inside the output folder or the working folder can be opened.")
        if typed.strip() and st.button("Open folder", key=f"saved_open_path_{state_key}"):
            candidate = Path(typed.strip().strip('"').strip("'")).expanduser()
            try:
                candidate = candidate.resolve(strict=True)
            except OSError:
                st.error("That folder does not exist.")
            else:
                if not (candidate / "manifest.csv").is_file():
                    st.error("No manifest.csv in that folder, so it is not a multi-protein run.")
                elif not _under_allowed_roots(candidate):
                    st.error("For safety, only folders inside the output folder or the working folder can be opened here.")
                elif run_kind(candidate) not in (None, kind):
                    st.error("That is a protein-list run; open it under Protein list." if kind == "genome"
                             else "That is a genome run; open it under Genome scan.")
                else:
                    st.session_state[state_key] = str(candidate)
                    st.rerun()


def _saved_bundles_picker() -> None:
    jobs = resources.output_path("jobs")
    bundles = {b["name"]: b for b in list_saved_bundles(jobs)}
    if not bundles:
        return
    with st.expander(f"Saved single-TF reports ({len(bundles)})", expanded=False):
        choice = st.selectbox("Report", list(bundles), index=None, placeholder="Choose a saved report…",
                              key="saved_bundle_pick",
                              format_func=lambda n: f'{n} · {bundles[n]["modified"]:%Y-%m-%d %H:%M}')
        if choice is not None and st.button("Open this report", key="saved_bundle_open"):
            try:
                st.session_state.last_result = api.Result.from_run(result_dict(bundles[choice]["path"]))
                st.session_state.last_command = None
                st.session_state.last_log = ""
                st.rerun()
            except (OSError, ValueError) as exc:
                st.error(f"Could not read that bundle: {exc}")


# ------------------------------------------------------------------------------------------- bundle view
def _sheet(head: dict) -> str:
    support = (head.get("family_support") or "not reported").replace("_", " ")
    agreement = head.get("agreement")
    basis = "exact agreement" if head.get("agreement_exact") else "legacy agreement"
    inducer = head.get("inducer_label") or "Unresolved"
    inducer_note = _klass(head["inducer_class"]) + (f" · {basis} {agreement:.2f}" if agreement is not None
                                                    and head.get("inducer") else "")
    operator = head.get("operator")
    op_note = " · ".join(x for x in (head.get("operator_locality"), head.get("operator_source")) if x)
    cells = (
        ("Family", f'<dd>{_e(head.get("family") or "—")}</dd><small>{_e(support)}</small>'),
        ("Possible inducer", f"<dd>{_e(inducer)}</dd><small>{inducer_note}</small>"),
        ("Primary operator", f'<dd class="mono">{_e(operator) if operator else "—"}</dd><small>{_e(op_note) or "&nbsp;"}</small>'),
        ("Regulon", f'<dd>{head.get("n_regulon_operons", 0)} operons</dd><small>{head.get("n_sites", 0)} putative sites</small>'),
    )
    body = "".join(f"<div><dt>{label}</dt>{value}</div>" for label, value in cells)
    return f'<dl class="sheet">{body}</dl>'


def _render_bundle(bundle: Path, key: str, *, log: str = "") -> None:
    """Headline calls, downloads and the full report for one finished bundle (single TF and batch alike)."""
    try:
        head = bundle_headline(bundle)
    except (OSError, ValueError) as exc:
        st.error(f"The saved bundle could not be read: {exc}")
        return
    st.markdown(_sheet(head), unsafe_allow_html=True)
    st.markdown('<div class="note"><b>PREDICTED_UNVERIFIED.</b> These are ranked computational calls for choosing '
                'experiments. A completed pipeline is not experimental validation.</div>', unsafe_allow_html=True)

    report = bundle / "REPORT.html"
    with st.expander("Download report & data", expanded=False):
        _full_download(bundle, "bundle")
        st.caption("Individual HTML downloads need their saved folder and report assets; the ZIP keeps them together.")
        downloads = st.columns(2)
        for index, (label, rel, mime) in enumerate((
            ("Report · HTML", "REPORT.html", "text/html"),
            ("Final operators · TSV", "operators/FINAL_operators.tsv", "text/tab-separated-values"),
            ("Ranked operators · JSON", "operators_ranked.json", "application/json"),
            ("Complete dossier · JSON", "dossier.json", "application/json"),
        )):
            path = bundle / rel
            if path.is_file():
                downloads[index % 2].download_button(label, path.read_bytes(), file_name=path.name, mime=mime,
                                                     width="stretch", key=f"dl_{key}_{rel}")
        _paths_block(bundle)

    if report.is_file():
        st.caption("Preview. File links are disabled here; use the full ZIP to open linked files.")
        st.iframe(_inline_report_assets(report.read_text(encoding="utf-8", errors="replace"), bundle),
                  width="stretch", height=1180)
    else:
        st.error("REPORT.html is missing from the completed bundle.")
    if log:
        with st.expander("Run log"):
            st.code(log[-12000:], language=None)


# ------------------------------------------------------------------------------------------- single TF
def _prediction_workspace() -> None:
    _section("01", "Analyze one regulator", "the report returns here when the run finishes")
    method = st.segmented_control(
        "Protein input", ("Protein accession", "Paste sequence", "Upload FASTA"),
        default="Protein accession", label_visibility="collapsed", width="stretch",
    )

    with st.form("predict_form", border=True):
        pasted, uploaded, protein_acc = "", None, ""
        if method == "Protein accession":
            protein_acc = st.text_input("Protein accession", placeholder="NP_417751.1",
                                        help="RefSeq/GenBank protein accessions are fetched from NCBI and cached.")
        elif method == "Paste sequence":
            pasted = st.text_area("Protein sequence", height=170,
                                  placeholder="Paste amino acids or a single-record FASTA…")
        else:
            uploaded = st.file_uploader("Single-protein FASTA", type=("fa", "faa", "fasta", "txt"))

        genome_acc = st.text_input(
            "Genome context", placeholder="Recommended · GCF_000005845.2 or NC_000913.3",
            help="Use the exact source assembly or replicon to prevent strain ambiguity.",
        )
        with st.expander("Advanced settings", expanded=False):
            run_name = st.text_input(
                "Output name (optional)", placeholder="Generated from the protein accession…",
                help="Completed results are never overwritten; a numeric suffix is added when needed.",
            )
            family_values = ["Auto-detect"] + sorted(PIPELINE_FAMILIES)
            family = st.selectbox("TF family", family_values,
                                  help="Auto-detect uses bundled references. Manual choices are limited to the "
                                       f"{len(PIPELINE_FAMILIES)} families supported by the full operator workflow.")
            effector = st.text_input("Known inducer / effector", placeholder="Optional, e.g. Zn2+")
            o1, o2 = st.columns(2)
            fold = o1.checkbox("Add apo-fold analysis", value=False,
                               help="Opt-in ESMFold2 enrichment; needs BIOHUB_TOKEN.")
            metalnet = o2.checkbox("Use MetalNet2 when available", value=True,
                                   help="If unavailable, the pipeline records an abstention and continues.")
            af3 = o1.checkbox("Write optional AF3 jobs", value=True,
                              help="Writes a manual hand-off file. Nothing in this pipeline reads it back.")
            offline = o2.checkbox("Use cached data only", value=False,
                                  help="Disables all network retrieval, including accession lookup.")

        submitted = st.form_submit_button("Predict operators, inducer & regulon", type="primary", width="stretch")

    if submitted:
        try:
            resolver_method = "Paste protein sequence" if method == "Paste sequence" else method
            source = _resolve_sequence(resolver_method, pasted, uploaded, protein_acc, offline=offline)
            cfg = RunConfig(allow_online=not offline, allow_folds=fold, allow_metalnet=metalnet, emit_af3_jobs=af3)
            default_name = protein_acc.strip() if method == "Protein accession" else "TF_prediction"
            job = _available_job_name(run_name.strip() or default_name)
            with st.status("Building the TF report…", expanded=True) as status:
                res, log = _run_with_log(lambda: api.predict(
                    source, name=job, family=None if family == "Auto-detect" else family,
                    organism=genome_acc.strip() or None, effector=effector.strip() or None,
                    fold=fold, cfg=cfg, verbose=True,
                ))
                if not res.ok:
                    raise RuntimeError(res.error or "The pipeline did not complete.")
                status.update(label="Report ready", state="complete", expanded=False)
            st.session_state.last_result = res
            st.session_state.last_log = log
            query = (protein_acc.strip() if method == "Protein accession"
                     else str(Path(res.bundle) / "input" / "query.fasta") if res.bundle else "protein.faa")
            # `tfop predict` refuses a finished job, so the reproduction goes to its own output folder.
            st.session_state.last_command = cli_command(
                "predict", source=query, name=job, out=str(resources.OUTPUT_DIR / "terminal"),
                family=None if family == "Auto-detect" else family, organism=genome_acc.strip() or None,
                effector=effector.strip() or None, fold=fold, offline=offline, metalnet=metalnet, af3=af3)
        except Exception as exc:
            st.error(f"Prediction stopped: {exc}")
            st.caption("Check the identifier or sequence, then open “Environment & data” below if the problem persists.")

    _saved_bundles_picker()
    result = st.session_state.get("last_result")
    if result is not None:
        _section("02", "Result", result.tf_id)
        bundle = Path(result.bundle) if result.bundle else None
        if not bundle:
            st.warning("The prediction returned without a written bundle.")
        else:
            _render_bundle(bundle, "single", log=st.session_state.get("last_log", ""))
            if st.session_state.get("last_command"):
                with st.expander("Reproduce in a terminal"):
                    st.code(st.session_state.last_command, language="bash")
                    st.caption("Writes to a separate output folder, so this saved result is left as it is.")


# ------------------------------------------------------------------------------------------- charts, ledgers
def _locus_label(row: dict) -> str:
    if row.get("contig") and row.get("start") not in (None, "") and row.get("end") not in (None, ""):
        try:
            return f'{row["contig"]}:{int(float(row["start"])):,}–{int(float(row["end"])):,} ({row.get("strand") or "?"})'
        except ValueError:
            return str(row["contig"])
    return "no coordinates"


def _short_name(row: dict) -> str:
    return row.get("gene") or row.get("locus_tag") or (row.get("run_name") or "?").split("__")[-1]


def _subtitle(row: dict) -> str:
    """Second line under a regulator: its annotated product, else its accession (never a repeat of the name)."""
    accession = (row.get("run_name") or "").split("__")[-1]
    text = (row.get("product") or "").strip() or (accession if accession != _short_name(row) else "")
    return f"<small>{_e(text)}</small>" if text else ""


def _chart(data: list[dict], spec: dict, height: int) -> None:
    st.vega_lite_chart(data, {**spec, "background": "transparent", "height": height, "config": _VEGA_THEME},
                       width="stretch", theme=None)


def _render_family_overview(rows: list[dict]) -> None:
    counts: dict[str, int] = {}
    for row in rows:
        family = row.get("family") or "Unresolved"
        counts[family] = counts.get(family, 0) + 1
    data = [{"family": family, "candidates": count,
             "kind": ("May include metalloregulators" if family in POTENTIAL_METALLOREGULATOR_FAMILIES else "Other family")}
            for family, count in counts.items()]
    _chart(data, {
        "mark": {"type": "bar", "cornerRadiusEnd": 1, "size": 15},
        "encoding": {
            "y": {"field": "family", "type": "nominal", "sort": "-x", "title": None},
            "x": {"field": "candidates", "type": "quantitative", "title": "candidate TFs", "axis": {"tickMinStep": 1}},
            "color": {"field": "kind", "type": "nominal", "title": None,
                      "scale": {"domain": ["May include metalloregulators", "Other family"], "range": ["#2D5B8C", "#B3AD9F"]}},
            "tooltip": [{"field": "family"}, {"field": "candidates"}, {"field": "kind", "title": "family type"}],
        },
    }, max(120, 30 * len(data) + 40))


def _render_inducer_overview(summary: list[dict]) -> None:
    counts: dict[tuple[str, str], int] = {}
    for row in summary:
        if not is_done(row.get("status")):
            continue
        key = (row.get("family") or row.get("discovery_family") or "Unresolved", row["class"])
        counts[key] = counts.get(key, 0) + 1
    if not counts:
        return
    data = [{"family": f, "inducer_class": c, "regulators": n, "class_order": CLASSES.index(c)}
            for (f, c), n in counts.items()]
    _chart(data, {
        "mark": {"type": "bar", "cornerRadiusEnd": 1, "size": 15},
        "encoding": {
            "y": {"field": "family", "type": "nominal", "sort": "-x", "title": None},
            "x": {"aggregate": "sum", "field": "regulators", "type": "quantitative", "title": "regulators", "axis": {"tickMinStep": 1}},
            "color": {"field": "inducer_class", "type": "nominal", "title": "inducer class",
                      "scale": {"domain": list(CLASSES), "range": [_COLORS[c] for c in CLASSES]}},
            "order": {"field": "class_order", "type": "quantitative"},
            "tooltip": [{"field": "family"}, {"field": "inducer_class", "title": "inducer class"}, {"field": "regulators"}],
        },
    }, max(120, 30 * len({d["family"] for d in data}) + 70))


def _candidate_ledger(rows: list[dict]) -> str:
    body = []
    for row in rows:
        body.append(
            f'<tr><td><span class="id">{_e(_short_name(row))}</span>{_subtitle(row)}</td>'
            f'<td>{_e(row.get("family"))}</td><td class="num">{_e(row.get("aa_len"))}</td>'
            f'<td class="op">{_e(_locus_label(row))}</td></tr>')
    return ('<div class="ledger-wrap"><table class="ledger"><thead><tr><th>regulator</th><th>family</th>'
            '<th class="num">aa</th><th>locus</th></tr></thead><tbody>' + "".join(body) + "</tbody></table></div>")


def _results_ledger(summary: list[dict]) -> str:
    body = []
    for row in summary:
        ok = is_done(row.get("status"))
        agreement = row.get("agreement_value")
        meter = (f'<span class="meter"><i style="width:{max(0, min(1, agreement)) * 100:.0f}%"></i></span>{agreement:.2f}'
                 if agreement is not None and row.get("inducer_top") else "—")
        operator = row.get("primary_operator") or ""
        shown = (operator[:24] + "…") if len(operator) > 24 else operator
        name = _short_name(row)
        status = "" if ok else f'<small>{_e(row.get("status") or "not run")}</small>'
        label = row.get("inducer_label") or row.get("inducer_top")
        inducer = ((_e(label) or "—") + "<small>" + _klass(row["class"]) + "</small>") if ok else "—"
        body.append(
            f'<tr class="{"" if ok else "dim"}"><td><span class="id">{_e(name)}</span>{status or _subtitle(row)}</td>'
            f'<td>{_e(row.get("family") or row.get("discovery_family"))}</td><td>{inducer}</td>'
            f'<td class="num">{meter if ok else "—"}</td>'
            f'<td class="op" title="{_e(operator)}">{_e(shown) or "—"}</td>'
            f'<td class="num">{_e(row.get("n_regulon_operons")) if ok and row.get("n_regulon_operons") != "" else "—"}</td></tr>')
    return ('<div class="ledger-wrap"><table class="ledger"><thead><tr><th>regulator</th><th>family</th><th>possible inducer</th>'
            '<th class="num">agreement</th><th>primary operator</th><th class="num">regulon</th></tr></thead><tbody>'
            + "".join(body) + "</tbody></table></div>")


def _scan_downloads(out: Path) -> None:
    files = (
        ("Summary · TSV", "summary.tsv", "text/tab-separated-values"),
        ("Candidate table · TSV", "candidate_tfs.tsv", "text/tab-separated-values"),
        ("TF loci · GFF3", "candidate_tfs.gff3", "text/plain"),
        ("Genome-browser track · BED", "candidate_tfs.bed", "text/plain"),
        ("Candidate proteins · FASTA", "candidate_tfs.faa", "text/plain"),
        ("Scan report · HTML", "SCAN_REPORT.html", "text/html"),
    )
    available = [(label, rel, mime) for label, rel, mime in files if (out / rel).is_file()]
    for start in range(0, len(available), 3):
        cols = st.columns(3)
        for col, (label, rel, mime) in zip(cols, available[start:start + 3]):
            path = out / rel
            col.download_button(label, path.read_bytes(), file_name=rel, mime=mime, width="stretch",
                                key=f"download_{out.name}_{rel}")


# ------------------------------------------------------------------------------------------- environment
def _environment_panel() -> None:
    with st.expander("Environment & data", expanded=False):
        st.markdown("The core reference database is already installed. This check covers executable tools and tokens.")
        if st.button("Check environment"):
            from predictor import setup_tools
            with st.spinner("Checking the same tool paths used by the pipeline…"):
                result = setup_tools.check()
            # `required` is True, False, or "genome" (needed only for genome discovery), as in `tfop setup`.
            needed = {True: "required", False: "optional"}
            rows = [{"engine": r["name"], "needed": needed.get(r["required"], "genome discovery only"),
                     "status": "ready" if r["present"] else "missing", "purpose": r["purpose"],
                     "how to install": "" if r["present"] else r["hint"]}
                    for r in result["rows"]]
            st.dataframe(rows, width="stretch", hide_index=True)
            missing = [r["name"] for r in result["rows"] if r["required"] is True and not r["present"]]
            discovery = [r["name"] for r in result["rows"] if r["required"] == "genome" and not r["present"]]
            if missing:
                st.warning("Install the missing required engines before a full prediction: " + ", ".join(missing))
            else:
                st.success("All required engines are available.")
            if discovery:
                st.info("Genome scans need " + ", ".join(discovery) + "; single-protein and protein-list "
                        "predictions do not.")
        st.caption("Run `tfop setup` in a terminal for the same table with install hints, or `tfop audit` for data and cache paths.")
        st.markdown(f'<div class="paths"><div><b>Results go to</b> {_e(resources.OUTPUT_DIR)}</div></div>',
                    unsafe_allow_html=True)


def _shared_exports(out: Path) -> None:
    _full_download(out, "run")
    st.caption("Individual HTML downloads require the saved run folder, linked job reports and assets.")
    for name, label in (("discovery/DISCOVERY_REPORT.html", "Discovery report"),
                        ("RUN_REPORT.html", "Run report"),
                        ("Supplementary_Table_S1.xlsx", "Supplementary table")):
        path = out / name
        if path.is_file():
            st.download_button(label, path.read_bytes(), file_name=path.name, key=f"shared_{path.name}")
    if (out / "af3").is_dir():
        _full_download(out / "af3", "AF3 jobs")
    if (out / "RUN_REPORT.html").is_file() and st.toggle("Preview run report"):
        st.caption("Preview: file links are disabled. Extract the full run ZIP to open linked reports and files.")
        st.iframe(_inline_report_assets((out / "RUN_REPORT.html").read_text(encoding="utf-8"), out), height=950)
    _paths_block(out)


# ------------------------------------------------------------------------------------------- genome / list
def _stage_discovery(protein_list: bool, source_mode: str, accession: str, upload, annotation, pasted: str,
                     name: str, engine: str, offline: bool) -> tuple[Path, dict]:
    # Validate everything that can be validated BEFORE a folder is reserved, so a mistyped form leaves nothing behind.
    pasted_input = None
    if source_mode == "Genome accession":
        if not accession.strip():
            raise ValueError("Enter a genome accession, for example GCF_000005845.2.")
    elif upload is None:
        if not protein_list:
            raise ValueError("Choose a local input file")
        pasted_input = pasted_protein_input(pasted)  # raises an actionable message when empty or malformed
    out = _available_output_path(_safe_name(name or accession or (upload.name if upload else "run"), "run"))
    kwargs = {"out_dir": out, "discover_only": True, "offline": offline}
    if source_mode == "Genome accession":
        kwargs["accession"] = accession.strip()
    else:
        staged = out / "uploads"
        staged.mkdir(parents=True, exist_ok=True)
        if pasted_input is not None:
            source = staged / ("input" + pasted_input[0])
            source.write_text(pasted_input[1], encoding="utf-8", newline="\n")
        else:
            source = staged / ("input" + Path(upload.name).suffix)
            source.write_bytes(upload.getvalue())
        if protein_list:
            kwargs["source"] = source
        else:
            kwargs[{"Local genome FASTA": "fasta", "GenBank": "genbank", "Proteome": "proteome"}[source_mode]] = source
        kwargs["organism"] = accession.strip() or None
        if annotation and source_mode == "Local genome FASTA":
            gff = staged / "input.gff"
            gff.write_bytes(annotation.getvalue())
            kwargs["gff"] = gff
    if not protein_list:
        kwargs.update(census_fallback=engine == _QUICK, census_only=engine == _QUICK)
    return out, {"kwargs": kwargs}


def _entry_form(protein_list: bool):
    """The discovery form. Returns the submitted values, or ``None`` before submission."""
    source_mode = "Protein list" if protein_list else st.radio(
        "Starting point", ("Genome accession", "Local genome FASTA", "GenBank", "Proteome"))
    engine_ok, engine_why = (True, "") if protein_list else _discovery_engine()
    if not protein_list and not engine_ok:
        st.warning("Production discovery (BITACORA) is not available on this computer, so the quick Pfam census is "
                   "selected. It is faster and coarser, and not the production route.")
        with st.expander("Why, and how to enable it"):
            st.write(engine_why)
            st.caption("Linux: `bash env/install_linux.sh --with-riio`. Windows: install WSL Ubuntu, then follow "
                       "docs/INSTALL_WINDOWS.md. `tfop setup` lists every engine.")
    with st.form("shared_entry_form"):
        accession = st.text_input("Genome accession (required)" if source_mode == "Genome accession"
                                  else "Genome accession (optional context)",
                                  placeholder="GCF_000005845.2" if source_mode == "Genome accession" else "")
        # Only the fields that apply to the chosen starting point: the radio above sits outside the form, so
        # changing it redraws immediately.
        kinds = {"Protein list": ("Protein file · FASTA or accession TSV", ("faa", "fa", "fasta", "tsv")),
                 "Local genome FASTA": ("Genome FASTA", ("fna", "fa", "fasta")),
                 "GenBank": ("GenBank file", ("gb", "gbk", "gbff")),
                 "Proteome": ("Proteome FASTA", ("faa", "fa", "fasta"))}
        upload = annotation = None
        if source_mode in kinds:
            upload = st.file_uploader(kinds[source_mode][0], type=kinds[source_mode][1])
        pasted = st.text_area("Or paste accessions / FASTA", height=120,
                              placeholder="NP_417751.1\nNP_415020.1\n…   or FASTA records") if protein_list else ""
        if source_mode == "Local genome FASTA":
            annotation = st.file_uploader("GFF3 annotation (optional)", type=("gff", "gff3"),
                                          help="Without it, genes are called with Pyrodigal.")
        name = st.text_input("Output name (optional)", placeholder="Defaults to the accession or file name")
        engine = _PRODUCTION if protein_list else st.radio(
            "Discovery engine", (_PRODUCTION, _QUICK), index=0 if engine_ok else 1,
            help="Production runs BITACORA with SSN-clade profiles. The quick census uses Pfam alone.")
        offline = st.checkbox("Use cached data only", value=False, help="No network retrieval at all.")
        submitted = st.form_submit_button("Read protein list" if protein_list else "Find candidate TFs",
                                          type="primary", width="stretch")
    if not submitted:
        return None
    return source_mode, accession, upload, annotation, pasted, name, engine, offline


def _run_panel(out: Path, state_key: str, *, protein_list: bool) -> None:
    try:
        rows = [row for row in read_manifest(out) if row.get("run_name")]
    except (OSError, ValueError) as exc:
        st.error(f"This run's manifest.csv cannot be read: {exc}")
        return
    summary = run_rows(out, rows)
    done = sum(1 for r in summary if is_done(r.get("status")))
    failed = sum(1 for r in summary if is_failed(r.get("status")))
    _section("02" if not done else "03", "Candidates" if not done else "Results",
             f'{out.name} · {len(rows)} candidates · {done} predicted' + (f" · {failed} failed" if failed else ""))
    _stepper(4 if done and done >= len(rows) else 3 if done else 2)
    flash = st.session_state.pop(f"{state_key}_flash", None)
    if flash:
        getattr(st, flash[0])(flash[1])

    if done:
        _render_inducer_overview(summary)
        st.markdown(_results_ledger(summary), unsafe_allow_html=True)
        tsv = out / "summary.tsv"
        if tsv.is_file():
            st.download_button("Summary · TSV", tsv.read_bytes(), file_name="summary.tsv",
                               mime="text/tab-separated-values", key=f"sum_{out.name}")
    else:
        _render_family_overview(rows)
        with st.expander(f"All {len(rows)} candidates", expanded=len(rows) <= 12):
            st.markdown(_candidate_ledger(rows), unsafe_allow_html=True)

    labels = {r["run_name"]: f'{_short_name(r)} · {r.get("family") or "Unresolved"}' for r in rows}
    if not done:  # nothing built yet: choosing and building is the next step
        st.subheader("Choose and build")
        _build_controls(out, rows, state_key, protein_list, labels)
    else:  # results first; building more is secondary
        completed = [r["run_name"] for r in summary if is_done(r.get("status"))
                     and (bundle := job_dir(out, r["run_name"])) and (bundle / "dossier.json").is_file()]
        if completed:
            st.subheader("Read one report")
            chosen = st.selectbox("Regulator", completed, index=None, placeholder="Choose a regulator…",
                                  key=f"report_{state_key}_{out.name}", format_func=lambda n: labels.get(n, n))
            if chosen:
                _render_bundle(out / "jobs" / chosen, f"{out.name}_{chosen}")
        with st.expander("Build more reports", expanded=False):
            _build_controls(out, rows, state_key, protein_list, labels)
    with st.expander("Exports and folder", expanded=False):
        _scan_downloads(out)
        _shared_exports(out)


def _build_controls(out: Path, rows: list[dict], state_key: str, protein_list: bool, labels: dict) -> None:
    """Family filter, optional per-regulator selection, run options and the build button."""
    from predictor import genome_run as GR
    families = sorted({r.get("family") or "Unresolved" for r in rows})
    chosen_families = st.multiselect("Families", families, default=families, key=f"fam_{state_key}_{out.name}")
    names = [r["run_name"] for r in rows if (r.get("family") or "Unresolved") in chosen_families]
    with st.expander(f"Fine-tune individual regulators ({len(names)})", expanded=False):
        selected = st.multiselect("Proteins to predict", names, default=names, format_func=lambda n: labels[n],
                                  key=f"pick_{state_key}_{out.name}_{'|'.join(chosen_families)}")
    st.markdown(f'<p class="fine"><b style="color:var(--ink)">{len(selected)} of {len(rows)}</b> regulators will be '
                'predicted. Unselected candidates stay in the manifest and reports.</p>', unsafe_allow_html=True)
    policy = recorded_policy(out)
    locked = bool(policy)
    tag = f"{state_key}_{out.name}_{'fixed' if locked else 'open'}"
    o1, o2, o3, o4 = st.columns(4)
    jobs = o1.number_input("Parallel proteins", min_value=1, max_value=16, value=1)
    fold = o2.checkbox("Add apo-fold analysis", value=bool(policy.get("fold", False)), disabled=locked,
                       key=f"fold_{tag}", help="Opt-in ESMFold2 enrichment; needs BIOHUB_TOKEN.")
    metalnet = o3.checkbox("Use MetalNet2 when available", value=bool(policy.get("metalnet", True)), disabled=locked,
                           key=f"metalnet_{tag}")
    af3 = o4.checkbox("Write AF3 jobs", value=bool(policy.get("af3", True)), disabled=locked, key=f"af3_{tag}",
                      help="A manual hand-off file; nothing reads it back.")
    offline = st.checkbox("Use cached data only", disabled=locked, key=f"offline_{tag}",
                          value=bool(policy.get("offline", st.session_state.get(f"{state_key}_offline", False))))
    if locked:
        st.caption("These options were fixed by the first build in this folder, so every report here is comparable. "
                   "For other options, start a new output folder.")
    if st.button("Build selected reports", disabled=not selected, type="primary"):
        try:
            # Every surface uses the shared runner; keep unselected discovery rows in the manifest.
            summary_path, _log = _run_with_log(lambda: GR.run_candidates(
                rows, out, jobs=int(jobs), offline=offline, fold=fold, metalnet=metalnet, af3=af3,
                selected_names=selected))
            from predictor.run_cli import _run_exit_status
            st.session_state[f"{state_key}_flash"] = (
                ("warning", "Some selected predictions failed. Review the run report and per-protein logs.")
                if _run_exit_status(summary_path) else ("success", "Selected reports finished."))
        except Exception as exc:  # a failed run is a message, not a traceback page
            st.error(str(exc) or type(exc).__name__)
        else:
            st.rerun()  # redraw so the results, charts and report picker reflect what was just written
    # `tfop batch` has no --select: for a subset there is no equivalent command, so none is shown.
    subset = 0 < len(selected) < len(rows)
    command = None if (protein_list and subset) else run_command(
        out, jobs=int(jobs), fold=fold, offline=offline, metalnet=metalnet, af3=af3,
        select=selected if subset else None)
    if command:
        with st.expander("Reproduce in a terminal"):
            st.code(command, language="bash")
            st.caption("Same inputs and folder layout: a run started here can be continued from the command line, "
                       "and the other way round.")


def _shared_workspace(*, protein_list=False) -> None:
    """Thin GUI over the same resumable runner and output layout as CLI/API."""
    state_key = "last_shared_batch" if protein_list else "last_shared_genome"
    _section("01", "Proteins" if protein_list else "Genome", "supply the input" if protein_list else "supply a genome")
    if not st.session_state.get(state_key):
        _stepper(1)
    submitted = _entry_form(protein_list)
    if submitted:
        source_mode, accession, upload, annotation, pasted, name, engine, offline = submitted
        try:
            out, plan = _stage_discovery(protein_list, source_mode, accession, upload, annotation, pasted,
                                         name, engine, offline)
            _run_with_log(lambda: (api.predict_many if protein_list else api.run_genome)(**plan["kwargs"]))
            st.session_state[state_key] = str(out)
            st.session_state[f"{state_key}_offline"] = offline
        except Exception as exc:  # parser and engine errors included
            st.error(str(exc) or type(exc).__name__)
    _saved_runs_picker(state_key)
    if st.session_state.get(state_key):
        out = Path(st.session_state[state_key])
        if (out / "manifest.csv").is_file():
            _run_panel(out, state_key, protein_list=protein_list)
        else:
            st.warning(f"{out} no longer holds a run. Choose another from Saved runs.")


def main() -> None:
    audit = resources.audit()
    _masthead(audit["database"])
    mode = _workflow_selector()
    if mode == "scan":
        _scope_panel()
    if mode in {"scan", "batch"}:
        _shared_workspace(protein_list=mode == "batch")
    else:
        _prediction_workspace()
    _environment_panel()
    st.markdown('<p class="fine" style="margin-top:2rem">Research software. Computational predictions are hypotheses '
                'and need experimental validation. Reference data are versioned, and completed results are never '
                'overwritten.</p>', unsafe_allow_html=True)


main()
