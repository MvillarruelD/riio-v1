"""One-page web interface for TF prediction and genome/proteome census.

The UI is intentionally a thin client over the same public API and scanners used by ``tfop``. It owns
presentation and input guidance only; scientific behavior remains in ``predictor.api``.
"""
from __future__ import annotations

import base64
import contextlib
import io
import mimetypes
import queue
import re
import tempfile
import threading
import time
from dataclasses import asdict
from pathlib import Path

import streamlit as st

from predictor import api, resources
from predictor.annotate.genome_scan import (PIPELINE_FAMILIES,
                                            POTENTIAL_METALLOREGULATOR_FAMILIES)
from predictor.config import RunConfig
from predictor.gui_support import (PIPELINE_RUN_LOCK, available_job_name,
                                   reserve_output_dir, safe_name)
from predictor.input_validation import known_families, read_protein


APP_TITLE = "TF operator predictor"  # Product-title placeholder: replace when the final project name is chosen.
GITHUB_URL = "https://github.com/MvillarruelD/riio-v1"


st.set_page_config(page_title=APP_TITLE, page_icon="◈", layout="wide",
                   initial_sidebar_state="collapsed")

_CSS = """
<style>
:root{
  --ink:#17231f;--muted:#647169;--paper:#f5f7f4;--surface:#ffffff;--line:#dbe2dd;
  --accent:#126652;--accent-dark:#0b4d3e;--soft:#e8f1ed;--warn:#8a5a20;--focus:#63a894;
}
[data-testid="stAppViewContainer"]{background:radial-gradient(circle at 88% 0,rgba(18,102,82,.08),transparent 30rem),var(--paper);color:var(--ink)}
[data-testid="stHeader"]{background:transparent}
[data-testid="stToolbar"]{visibility:hidden;height:0}
.block-container{max-width:1120px;padding-top:1.4rem;padding-bottom:4rem}
h1,h2,h3{color:var(--ink);letter-spacing:-.025em;text-wrap:balance}
h1{font-size:clamp(2.35rem,4vw,3.75rem)!important;line-height:1.02!important;max-width:18ch;margin:.7rem 0 .65rem!important}
h2{font-size:1.55rem!important;margin-top:1.7rem!important}
h3{font-size:1.08rem!important}
p,li,label{line-height:1.55} p{max-width:66ch}
[data-testid="stRadio"] label p,[data-testid="stCheckbox"] label p,
[data-testid="stFileUploader"] label p{color:var(--ink)!important}
.brand-row{display:flex;align-items:center;justify-content:space-between;gap:1rem;border-bottom:1px solid var(--line);
  padding-bottom:.72rem}.eyebrow{font-size:.69rem;letter-spacing:.14em;text-transform:uppercase;color:var(--accent);font-weight:750}
.project-link{color:var(--accent-dark)!important;text-decoration:underline;text-decoration-color:#96b5aa;text-underline-offset:3px;font-size:.82rem;font-weight:650;
  transition:text-decoration-color .2s ease,color .2s ease}.project-link:hover{color:var(--accent)!important;text-decoration-color:var(--accent)}
.lead{font-size:1.02rem;color:#4d5c54;max-width:60ch;margin-bottom:1rem;text-wrap:pretty}
.topline{display:flex;gap:.55rem;flex-wrap:wrap;color:var(--muted);font-size:.78rem;margin:.35rem 0 1.35rem}
.status-chip{display:inline-flex;align-items:center;gap:.38rem;background:rgba(255,255,255,.72);border:1px solid var(--line);
  border-radius:6px;padding:.34rem .55rem;font-variant-numeric:tabular-nums}.status-chip b{color:var(--ink)}
.dot{display:inline-block;width:.45rem;height:.45rem;border-radius:50%;background:var(--accent)}
.truth{background:var(--soft);border-left:3px solid var(--accent);padding:.85rem 1rem;margin:.8rem 0 1.4rem;
  color:#31423a;font-size:.9rem}
.mode-title{font-size:.73rem;font-weight:700;letter-spacing:.08em;text-transform:uppercase;color:var(--muted);
  margin:1.25rem 0 .3rem}
.scope-line{margin:.9rem 0 0;color:var(--muted);font-size:.8rem}.scope-line b{color:var(--ink)}
.workspace-head{display:flex;align-items:flex-end;justify-content:space-between;gap:1rem;margin-top:.75rem}
.workspace-head h2{margin-bottom:.15rem!important}.workspace-head p{margin:.15rem 0;color:var(--muted);font-size:.87rem}
[data-testid="stForm"]{background:rgba(255,255,255,.86);border:1px solid var(--line);border-radius:12px;
  padding:1rem 1.15rem 1.2rem;box-shadow:0 14px 36px rgba(31,56,45,.04)}
[data-testid="stTextInput"] input,[data-testid="stTextArea"] textarea{background:#fff;border-color:#cfd7cf}
[data-testid="stTextInput"] input:focus,[data-testid="stTextArea"] textarea:focus{border-color:var(--accent);
  box-shadow:0 0 0 1px var(--accent)}
[data-testid="stButton"] button,[data-testid="stFormSubmitButton"] button{border-radius:8px;transition:transform .18s ease,
  background-color .18s ease,border-color .18s ease;min-height:2.75rem;font-weight:650;touch-action:manipulation}
[data-testid="stButton"] button:hover,[data-testid="stFormSubmitButton"] button:hover{transform:translateY(-1px)}
[data-testid="stButton"] button:active,[data-testid="stFormSubmitButton"] button:active{transform:translateY(1px)}
button[kind="primary"]{background:var(--accent)!important;border-color:var(--accent)!important}
button[kind="primary"]:hover{background:var(--accent-dark)!important}
button:disabled,button[disabled]{background:#e3e8e4!important;border-color:#d5ddd7!important;color:#78847d!important;
  transform:none!important}
button:focus-visible,input:focus-visible,textarea:focus-visible,[role="radio"]:focus-visible{outline:3px solid rgba(99,168,148,.42)!important;outline-offset:2px}
[data-testid="stSegmentedControl"]{margin:.15rem 0 .4rem}
[data-testid="stSegmentedControl"] button,[data-testid="stSegmentedControl"] [role="radio"]{
  min-height:3.05rem;font-weight:650;background:#fff!important;color:var(--ink)!important;border-color:#cbd5ce!important}
[data-testid="stSegmentedControl"] button[aria-pressed="true"],
[data-testid="stSegmentedControl"] [role="radio"][aria-checked="true"]{
  background:var(--accent)!important;color:#fff!important;border-color:var(--accent)!important}
[data-testid="stSegmentedControl"] button[aria-pressed="true"] p,
[data-testid="stSegmentedControl"] [role="radio"][aria-checked="true"] p{color:#fff!important}
[data-testid="stBaseButton-segmented_controlActive"]{background:var(--accent)!important;color:#fff!important;
  border-color:var(--accent)!important}
[data-testid="stBaseButton-segmented_controlActive"] p{color:#fff!important}
[data-testid="stBaseButton-segmented_control"],[data-testid="stBaseButton-segmented_controlInactive"]{
  background:#fff!important;color:var(--ink)!important;border-color:#cbd5ce!important}
[data-testid="stBaseButton-segmented_control"] p,[data-testid="stBaseButton-segmented_controlInactive"] p{color:var(--ink)!important}
[data-testid="stMetric"]{background:transparent;border-top:2px solid var(--ink);padding-top:.7rem}
[data-testid="stMetricLabel"]{color:var(--muted);letter-spacing:.04em}
[data-testid="stMetricValue"]{font-variant-numeric:tabular-nums;font-size:1.35rem}
[data-testid="stExpander"]{border-color:var(--line);background:rgba(255,255,255,.58)}
[data-testid="stDataFrame"]{border:1px solid var(--line);border-radius:8px;overflow:hidden}
.result-head{margin-top:2.5rem;padding-top:1.2rem;border-top:1px solid var(--line)}
.status-good{color:var(--accent);font-weight:700}.status-warn{color:var(--warn);font-weight:700}
.flow{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:7px;margin:.75rem 0 1.1rem}
.flow-step{border-top:2px solid #c3ccc6;padding:.55rem .15rem;color:var(--muted);font-size:.72rem}
.flow-step b{display:block;color:var(--ink);font-size:.92rem;margin-top:.1rem}.flow-step.active{border-color:var(--accent)}
.flow-step.done{border-color:var(--accent-dark)}.section-kicker{font-size:.72rem;letter-spacing:.12em;text-transform:uppercase;
 color:var(--accent);font-weight:750;margin-top:1.8rem}.legend-note{font-size:.82rem;color:var(--muted);margin:.25rem 0 .8rem}
.legend-chip{display:inline-flex;align-items:center;gap:.35rem;margin-right:1rem}.legend-chip i{width:18px;height:8px;
 border-radius:2px;background:#79a998}.legend-chip.metal i{background:var(--accent-dark)}
.skip-link{position:fixed;left:1rem;top:-4rem;z-index:100;background:var(--ink);color:#fff!important;padding:.6rem .8rem;border-radius:6px}
.skip-link:focus{top:1rem}.technical-path{overflow-wrap:anywhere;font-family:ui-monospace,SFMono-Regular,Consolas,monospace;font-size:.75rem}
@media(max-width:700px){.block-container{padding:1rem .9rem 3rem}h1{font-size:2.35rem!important}.brand-row{align-items:center}
  .flow-step b{font-size:.78rem}.workspace-head{align-items:flex-start;flex-direction:column}
  .topline{gap:.4rem}.status-chip{font-size:.72rem}}
@media(prefers-reduced-motion:reduce){*,*::before,*::after{scroll-behavior:auto!important;animation:none!important;transition:none!important}}
</style>
"""
st.markdown(_CSS, unsafe_allow_html=True)


def _safe_name(value: str, fallback: str = "query") -> str:
    return safe_name(value, fallback)


def _available_output_path(stem: str) -> Path:
    """Return a scan output folder that never overwrites an earlier completed census."""
    return reserve_output_dir(resources.OUTPUT_DIR, stem)


_PRODUCTION = "Production · BITACORA + SSN-clade profiles"
_QUICK = "Quick census · Pfam (not the production route)"


@st.cache_data(ttl=600, show_spinner=False)
def _discovery_engine() -> tuple[bool, str]:
    """Whether production discovery (BITACORA) can run here; probing WSL is slow, so it is cached."""
    from predictor.discovery import bitacora
    return bitacora.engine_available()


def _available_job_name(stem: str) -> str:
    """Return a stable, filesystem-safe job name without making users resolve output collisions."""
    return available_job_name(resources.output_path("jobs"), stem)


def _inline_report_assets(text: str, bundle: Path) -> str:
    """Turn report-relative image paths into data URIs for the embedded same-page preview."""
    def repl(match):
        rel = match.group(1)
        path = bundle / rel
        if not path.is_file():
            return match.group(0)
        mime = mimetypes.guess_type(path.name)[0] or "application/octet-stream"
        encoded = base64.b64encode(path.read_bytes()).decode("ascii")
        return f'src="data:{mime};base64,{encoded}"'
    return re.sub(r'src="([^":]+)"', repl, text)


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
    """Run a blocking pipeline call while keeping its latest progress visible."""
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
    placeholder = st.empty()
    chunks: list[str] = []
    while worker.is_alive() or not sink.empty():
        while not sink.empty():
            chunks.append(sink.get_nowait())
        latest = "".join(chunks)[-9000:].strip()
        if latest:
            placeholder.code(latest, language=None)
        time.sleep(.15)
    worker.join()
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
            raise ValueError("Enter a RefSeq/GenBank protein accession, for example WP_000240047.1.")
        if offline:
            raise ValueError("Protein-accession lookup needs the network. Turn off offline mode or paste the sequence.")
        from predictor.annotate.protein_seqs import fetch_protein_sequences
        found = fetch_protein_sequences([acc], allow_ncbi=True, verbose=True)
        source = found.get(acc.split(".")[0], "")
        if not source:
            raise ValueError(f"No protein sequence was found for {acc!r}. Check the accession or paste FASTA.")
    read_protein(source)  # actionable validation before the long-running pipeline starts
    return source


def _workflow_selector() -> str:
    """Persistent primary navigation. Genome scanning is the first-run default."""
    st.markdown('<div class="mode-title">Start an analysis</div>', unsafe_allow_html=True)
    choice = st.segmented_control(
        "Analysis type", ("Genome scan", "Single TF"), default="Genome scan",
        key="workflow_choice", label_visibility="collapsed", width="stretch",
    )
    if choice == "Single TF":
        st.caption("Go directly from one protein accession or sequence to operators, possible inducer and regulon.")
        return "predict"
    st.caption("Find TF candidates in a genome, select one or more, then build complete reports.")
    return "scan"


def _scope_panel() -> None:
    """Keep the scientific scope visible without putting a 12-family paragraph in the main path."""
    st.markdown(
        '<p class="scope-line"><b>Scope:</b> the census screens 22 bacterial TF families; complete '
        'operator–inducer–regulon reports currently support 12.</p>', unsafe_allow_html=True,
    )
    with st.expander("View supported full-report families", expanded=False):
        st.write(" · ".join(PIPELINE_FAMILIES))
        st.caption(
            "Dark family bars flag families that include known metalloregulators. This is a family-level "
            "cue, not evidence that an individual protein senses a metal."
        )


def _render_prediction(res) -> None:
    st.markdown('<div class="result-head"><span class="eyebrow">Completed prediction</span></div>',
                unsafe_allow_html=True)
    d = res.raw
    cols = st.columns(3)
    cols[0].metric("Family", res.family or "—")
    cols[1].metric("Possible inducer", res.inducer or "—")
    cols[2].metric("Primary operator", res.primary_operator or "—")
    st.caption(f'Genome: {res.genome or "not resolved"} · Family support: '
               f'{(res.family_support or "not reported").replace("_", " ")}')
    st.markdown(
        '<div class="truth"><b>Interpretation.</b> This is a ranked <b>PREDICTED_UNVERIFIED</b> result. '
        'Use it to prioritize experiments; a green pipeline check is not experimental validation.</div>',
        unsafe_allow_html=True,
    )
    bundle = Path(res.bundle) if res.bundle else None
    if not bundle:
        st.warning("The prediction returned without a written bundle.")
        return

    report = bundle / "REPORT.html"
    st.subheader("Report")
    st.caption("Read the complete report below. Open exports only when you need to reuse the result.")
    with st.expander("Download report & data", expanded=False):
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
                                                     width="stretch")
        st.markdown(f'<div class="technical-path">Saved in {bundle}</div>', unsafe_allow_html=True)

    if report.is_file():
        embedded = _inline_report_assets(report.read_text(encoding="utf-8", errors="replace"), bundle)
        st.iframe(embedded, width="stretch", height=1180)
    else:
        st.error("REPORT.html is missing from the completed bundle.")

    with st.expander("Run log and machine-readable summary"):
        st.json({k: d.get(k) for k in ("family", "family_support", "genome", "tf_locus",
                                       "n_putative_sites", "motif_consensus", "motif_mean_ic",
                                       "inducer", "ssn_cluster", "designed_operator")})
        if st.session_state.get("last_log"):
            st.code(st.session_state.last_log[-12000:], language=None)


def _prediction_workspace() -> None:
    st.markdown('<div id="analysis-workspace" class="workspace-head"><div><h2>Analyze one TF</h2>'
                '<p>Start with the simplest identifier you have. The report returns here when the run finishes.</p>'
                '</div></div>', unsafe_allow_html=True)
    method = st.segmented_control(
        "Protein input",
        ("Protein accession", "Paste sequence", "Upload FASTA"),
        default="Protein accession",
        label_visibility="collapsed",
        width="stretch",
    )

    with st.form("predict_form", border=True):
        pasted, uploaded, protein_acc = "", None, ""
        if method == "Protein accession":
            protein_acc = st.text_input("Protein accession", placeholder="WP_000240047.1",
                                        help="RefSeq/GenBank protein accessions are fetched from NCBI and cached.")
        elif method == "Paste sequence":
            pasted = st.text_area("Protein sequence", height=170,
                                  placeholder="Paste amino acids or a single-record FASTA…")
        else:
            uploaded = st.file_uploader("Single-protein FASTA", type=("fa", "faa", "fasta", "txt"))

        genome_acc = st.text_input(
            "Genome context",
            placeholder="Recommended · GCF_000005845.2 or NC_000913.3",
            help="Use the exact source assembly or replicon to prevent strain ambiguity.",
        )
        with st.expander("Advanced settings", expanded=False):
            run_name = st.text_input(
                "Output name (optional)", placeholder="Generated from the protein accession…",
                help="Completed results are never overwritten; a numeric suffix is added when needed.",
            )
            family_values = ["Auto-detect"] + sorted(PIPELINE_FAMILIES)
            family = st.selectbox("TF family", family_values,
                                  help="Auto-detect uses bundled references. Manual choices are limited to the 12 families supported by the full operator workflow.")
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
            cfg = RunConfig(allow_online=not offline, allow_folds=fold, allow_metalnet=metalnet,
                            emit_af3_jobs=af3)
            default_name = protein_acc.strip() if method == "Protein accession" else "TF_prediction"
            requested_name = run_name.strip() or default_name
            with st.status("Building the TF report…", expanded=True) as status:
                res, log = _run_with_log(lambda: api.predict(
                    source, name=_available_job_name(requested_name),
                    family=None if family == "Auto-detect" else family,
                    organism=genome_acc.strip() or None, effector=effector.strip() or None,
                    fold=fold, cfg=cfg, verbose=True,
                ))
                if not res.ok:
                    raise RuntimeError(res.error or "The pipeline did not complete.")
                status.update(label="Report ready", state="complete", expanded=False)
            st.session_state.last_result = res
            st.session_state.last_log = log
        except Exception as exc:
            st.error(f"Prediction stopped: {exc}")
            st.caption("Check the identifier or sequence, then open “Environment & data” below if the problem persists.")

    if st.session_state.get("last_result") is not None:
        _render_prediction(st.session_state.last_result)


def _scan_flow(active: int) -> None:
    labels = (("1", "Supply a genome"), ("2", "Choose candidate TFs"),
              ("3", "Predict operators & inducers"))
    items = []
    for number, label in labels:
        index = int(number)
        state = "done" if index < active else "active" if index == active else ""
        items.append(f'<div class="flow-step {state}"><span>STEP {number}</span><b>{label}</b></div>')
    st.markdown(f'<div class="flow">{"".join(items)}</div>', unsafe_allow_html=True)


def _locus_label(row: dict) -> str:
    if row.get("contig") and row.get("start") is not None and row.get("end") is not None:
        return f'{row["contig"]}:{int(row["start"]):,}-{int(row["end"]):,} ({row.get("strand") or "?"})'
    return "Not available"


def _render_family_overview(rows: list[dict]) -> None:
    counts: dict[str, int] = {}
    for row in rows:
        family = row.get("family") or "Unresolved"
        counts[family] = counts.get(family, 0) + 1
    data = []
    for family, count in counts.items():
        category = ("Possible metalloregulator family"
                    if family in POTENTIAL_METALLOREGULATOR_FAMILIES
                    else "Other full-report family" if family in PIPELINE_FAMILIES else "Census only")
        data.append({"family": family, "candidates": count, "category": category})
    st.markdown('<div class="legend-note"><span class="legend-chip metal"><i></i>Possible metalloregulator family</span>'
                '<span class="legend-chip"><i></i>Other family</span></div>', unsafe_allow_html=True)
    st.caption("Darker bars mark families that include known metal-responsive regulators. The color is a family-level cue, not proof that an individual candidate senses a metal.")
    st.vega_lite_chart(data, {
        "height": max(150, 27 * len(data)),
        "mark": {"type": "bar", "cornerRadiusEnd": 3},
        "encoding": {
            "y": {"field": "family", "type": "nominal", "sort": "-x", "title": None},
            "x": {"field": "candidates", "type": "quantitative", "title": "Candidate TFs"},
            "color": {"field": "category", "type": "nominal", "legend": None,
                      "scale": {"domain": ["Possible metalloregulator family", "Other full-report family", "Census only"],
                                "range": ["#0e4f40", "#75a996", "#aeb9b3"]}},
            "tooltip": [{"field": "family", "type": "nominal"},
                        {"field": "candidates", "type": "quantitative"},
                        {"field": "category", "type": "nominal"}],
        },
    }, width="stretch")


def _render_inducer_overview(rows: list[dict], predictions: dict) -> None:
    counts: dict[tuple[str, str], int] = {}
    for row in rows:
        pred = predictions.get(row["protein_id"])
        if not pred:
            continue
        key = (row.get("family") or "Unresolved family", pred.get("inducer") or "Unresolved")
        counts[key] = counts.get(key, 0) + 1
    if not counts:
        return
    data = [{"family": family, "possible_inducer": inducer, "candidates": count}
            for (family, inducer), count in counts.items()]
    st.markdown('<div class="section-kicker">Prediction overview</div>', unsafe_allow_html=True)
    st.markdown("### Possible inducer distribution")
    st.caption("The family layout is reprised after prediction. Segments count computational inducer calls; unresolved calls remain visible.")
    st.vega_lite_chart(data, {
        "height": max(150, 30 * len({d["family"] for d in data})),
        "mark": {"type": "bar", "cornerRadiusEnd": 2},
        "encoding": {
            "y": {"field": "family", "type": "nominal", "sort": "-x", "title": None},
            "x": {"aggregate": "sum", "field": "candidates", "type": "quantitative", "title": "Predicted TFs"},
            "color": {"field": "possible_inducer", "type": "nominal", "title": "Possible inducer"},
            "tooltip": [{"field": "family", "type": "nominal"},
                        {"field": "possible_inducer", "type": "nominal"},
                        {"field": "candidates", "type": "quantitative"}],
        },
    }, width="stretch")


def _scan_downloads(out: Path) -> None:
    files = (
        ("TF loci · GFF3", "candidate_tfs.gff3", "text/plain"),
        ("Genome-browser track · BED", "candidate_tfs.bed", "text/plain"),
        ("Source annotation · GFF3", "source_annotation.gff3", "text/plain"),
        ("Source genome · FASTA", "source_genome.fna", "text/plain"),
        ("Candidate table · TSV", "candidate_tfs.tsv", "text/tab-separated-values"),
        ("Candidate proteins · FASTA", "candidate_tfs.faa", "text/plain"),
        ("Candidate data · JSON", "candidate_tfs.json", "application/json"),
        ("Scan report · HTML", "SCAN_REPORT.html", "text/html"),
    )
    available = [(label, rel, mime) for label, rel, mime in files if (out / rel).is_file()]
    for start in range(0, len(available), 4):
        cols = st.columns(4)
        for col, (label, rel, mime) in zip(cols, available[start:start + 4]):
            path = out / rel
            col.download_button(label, path.read_bytes(), file_name=rel, mime=mime,
                                width="stretch", key=f"download_{out.name}_{rel}")


def _render_batch_report(rows: list[dict], predictions: dict) -> None:
    completed = [r["protein_id"] for r in rows
                 if predictions.get(r["protein_id"], {}).get("bundle")]
    if not completed:
        return
    st.markdown("### Individual TF reports")
    chosen = st.selectbox("Report to preview", completed, key="batch_report_choice",
                          format_func=lambda pid: f'{pid} · {next(r.get("family") for r in rows if r["protein_id"] == pid)}')
    pred = predictions[chosen]
    bundle = Path(pred["bundle"])
    report = bundle / "REPORT.html"
    metrics = st.columns(3)
    metrics[0].metric("TF", chosen)
    metrics[1].metric("Possible inducer", pred.get("inducer") or "Unresolved")
    metrics[2].metric("Primary operator", pred.get("primary_operator") or "Unresolved")
    if report.is_file():
        st.download_button("Download this full report", report.read_bytes(), file_name=f"{_safe_name(chosen)}_REPORT.html",
                           mime="text/html")
        embedded = _inline_report_assets(report.read_text(encoding="utf-8", errors="replace"), bundle)
        st.iframe(embedded, width="stretch", height=1050)


def _scan_workspace() -> None:
    scan = st.session_state.get("last_scan")
    submitted = False
    if scan:
        head, action = st.columns([3, 1])
        with head:
            st.markdown('<div id="analysis-workspace" class="workspace-head"><div><h2>Genome scan</h2>'
                        f'<p>{scan.get("source", "Genome")} · candidate census complete</p></div></div>',
                        unsafe_allow_html=True)
        with action:
            if st.button("New genome scan", width="stretch"):
                st.session_state.pop("last_scan", None)
                st.session_state.pop("scan_selected_ids", None)
                st.rerun()
        _scan_flow(3 if st.session_state.get("scan_selected_ids") else 2)
    else:
        st.markdown('<div id="analysis-workspace" class="workspace-head"><div><h2>Scan a genome</h2>'
                    '<p>Use an NCBI accession for the shortest path, or upload local sequence files.</p>'
                    '</div></div>', unsafe_allow_html=True)
        _scan_flow(1)
        source_method = st.segmented_control(
            "Genome input", ("NCBI accession", "Upload FASTA"), default="NCBI accession",
            label_visibility="collapsed", width="stretch",
        )
        with st.form("scan_form", border=True):
            genome_acc, upload, annotation_upload = "", None, None
            if source_method == "NCBI accession":
                genome_acc = st.text_input(
                    "Genome accession", placeholder="GCF_000005845.2 or NC_000913.3",
                    help="NCBI Datasets fetches the genome and annotation once, then reuses the local cache.",
                )
            else:
                upload = st.file_uploader("Genome or proteome FASTA", type=("fna", "fa", "fasta", "faa"),
                                          key="scan_upload")
                annotation_upload = st.file_uploader(
                    "Genome annotation (optional GFF3)", type=("gff", "gff3"), key="scan_gff_upload",
                    help="Recommended for nucleotide genomes so source protein accessions, locus tags, genes and products are preserved.")
            engine_ok, engine_why = _discovery_engine()
            with st.expander("Advanced scan settings", expanded=not engine_ok):
                method = st.radio(
                    "Discovery method", (_PRODUCTION, _QUICK), index=0 if engine_ok else 1,
                    help="Production is the route the published runs used (`tfop genome`). The Pfam "
                         "census is faster and reaches nearly the same set, but it is a different method.",
                )
                if not engine_ok:
                    st.caption(f"Production discovery is unavailable here: {engine_why}. "
                               "See docs/ENGINES.md to install BITACORA.")
                families = st.multiselect(
                    "Limit the initial census", sorted(known_families()), placeholder="All TF families",
                    help="Leave empty to scan every family. Filter and select candidates after the scan.",
                )
                out_name = st.text_input(
                    "Output name (optional)", placeholder="Generated from the accession or filename…",
                    help="A numeric suffix is added automatically when an earlier scan uses the same name.",
                )
            submitted = st.form_submit_button("Find candidate TFs", type="primary", width="stretch")

    if submitted:
        if source_method == "NCBI accession" and not genome_acc.strip():
            st.error("Enter a genome accession, such as GCF_000005845.2.")
            return
        if source_method == "Upload FASTA" and upload is None:
            st.error("Choose a genome or proteome FASTA to scan.")
            return
        try:
            from predictor.annotate import genome_scan as GS

            production = method == _PRODUCTION

            def run_scan(source: Path, source_label: str, annotation: Path | None = None):
                suggested = out_name.strip() or genome_acc.strip() or Path(source_label).stem
                out = _available_output_path(f"scan_{_safe_name(suggested, 'genome_scan')}")
                label = ("Searching proteins with BITACORA and the SSN-clade profiles…" if production
                         else "Scanning proteins against the Pfam TF-family library…")
                with st.status(label, expanded=True) as status:
                    if production:
                        from predictor import genome_run as GR
                        found, _engine, kind = GR.discover(source, gff=annotation, work_dir=out / "bitacora")
                        cands, n_proteins = GR.as_census_candidates(found, source, annotation)
                        if families:
                            cands = [c for c in cands if c.family in families]
                    else:
                        cands, kind, n_proteins = GS.scan(source, gff=annotation,
                                                          families=set(families) or None)
                    GS.write_outputs(cands, out, kind=kind, n_proteins=n_proteins,
                                     source=source_label, source_path=source,
                                     annotation_path=annotation)
                    status.update(label=f"Census complete · {len(cands)} candidate TFs",
                                  state="complete", expanded=False)
                return cands, kind, n_proteins, out

            if source_method == "NCBI accession":
                from predictor.annotate import genome_mirror
                with st.spinner("Fetching the genome and annotation from NCBI (cached after the first run)…"):
                    source, annotation = genome_mirror.ensure_genome(genome_acc.strip(), allow_download=True)
                if source is None:
                    raise RuntimeError(
                        "The genome could not be fetched. Check the accession and confirm the NCBI "
                        "Datasets CLI is installed with `tfop audit`."
                    )
                cands, kind, n_proteins, out = run_scan(source, genome_acc.strip(), annotation)
                organism = genome_acc.strip()
            else:
                suffix = Path(upload.name).suffix or ".fasta"
                with tempfile.TemporaryDirectory(prefix="tfop_gui_scan_") as td:
                    source = Path(td) / f"input{suffix}"
                    source.write_bytes(upload.getvalue())
                    annotation = None
                    if annotation_upload is not None:
                        annotation = Path(td) / "input.gff3"
                        annotation.write_bytes(annotation_upload.getvalue())
                    cands, kind, n_proteins, out = run_scan(source, upload.name, annotation)
                organism = None
            st.session_state.last_scan = {
                "rows": [asdict(c) for c in cands], "out": str(out), "kind": kind,
                "n_proteins": n_proteins, "source": genome_acc.strip() if organism else upload.name,
                "organism": organism, "predictions": {}, "logs": {},
                "method": "production" if production else "pfam_census",
            }
            st.session_state.pop("scan_selected_ids", None)
            st.rerun()
        except Exception as exc:
            st.error(f"Scan stopped: {exc}")

    scan = st.session_state.get("last_scan")
    if not scan:
        return
    rows, out = scan["rows"], Path(scan["out"])
    st.markdown('<div class="result-head"><span class="eyebrow">Completed census</span></div>',
                unsafe_allow_html=True)
    located = sum(1 for row in rows if row.get("contig") and row.get("start") is not None)
    a, b, c, d = st.columns(4)
    a.metric("Input type", scan["kind"].replace("_", " "))
    b.metric("Proteins scanned", f"{scan['n_proteins']:,}")
    c.metric("Candidate TFs", f"{len(rows):,}")
    d.metric("Located on genome", f"{located:,}")
    if scan["kind"] == "proteome":
        st.info("This was a protein-only input, so genomic coordinates are unavailable. Upload the genome FASTA with its GFF3—or use an NCBI genome accession—to map each TF to a locus.")
    elif scan["kind"] == "genome":
        st.info("No source annotation was supplied. Loci below are predicted gene coordinates from Pyrodigal; the GFF3/BED outputs preserve that provenance.")
    else:
        st.success("Source annotation used: protein accessions and genomic loci were preserved from GFF3.")

    if rows:
        family_count = len({row.get("family") or "Unresolved" for row in rows})
        with st.expander(f"Family overview · {family_count} families", expanded=False):
            _render_family_overview(rows)

        st.markdown('<div class="section-kicker">Step 2</div>', unsafe_allow_html=True)
        st.markdown("### Choose TFs to analyze")
        st.caption("Filter by family, then select individual candidates or all supported candidates.")
        family_values = sorted({r.get("family") or "Unresolved" for r in rows})
        family_filter = st.multiselect("Show families", family_values, placeholder="All families",
                                       key=f"scan_family_filter_{out.name}")
        visible = [r for r in rows if not family_filter or (r.get("family") or "Unresolved") in family_filter]
        if "scan_selected_ids" not in st.session_state:
            st.session_state.scan_selected_ids = set()
        selected = set(st.session_state.scan_selected_ids)
        st.caption(f"{len(selected)} selected across all filters")
        controls = st.columns(3)
        if controls[0].button("Select shown", width="stretch"):
            st.session_state.scan_selected_ids = selected | {r["protein_id"] for r in visible}
            st.rerun()
        if controls[1].button("Select all supported", width="stretch"):
            st.session_state.scan_selected_ids = {r["protein_id"] for r in rows
                                                  if r.get("family") in PIPELINE_FAMILIES}
            st.rerun()
        if controls[2].button("Clear selection", width="stretch"):
            st.session_state.scan_selected_ids = set()
            st.rerun()

        editor_rows = []
        predictions = scan.get("predictions", {})
        for row in visible:
            pred = predictions.get(row["protein_id"], {})
            editor_rows.append({
                "run": row["protein_id"] in selected, "protein_id": row["protein_id"],
                "locus_tag": row.get("locus_tag") or "—", "locus": _locus_label(row),
                "family": row.get("family") or "Unresolved", "evalue": row["evalue"], "bits": row["bits"],
                "prediction": (pred.get("status") or "Ready")
                if row.get("family") in PIPELINE_FAMILIES else "Census only",
            })
        edited = st.data_editor(
            editor_rows, hide_index=True, width="stretch",
            disabled=["protein_id", "locus_tag", "locus", "family", "evalue", "bits", "prediction"],
            column_order=("run", "protein_id", "family", "evalue", "bits", "locus_tag", "locus", "prediction"),
            column_config={
                "run": st.column_config.CheckboxColumn("Analyze", help="Select this TF for the full operator, inducer and regulon workflow."),
                "protein_id": st.column_config.TextColumn("Protein"),
                "locus_tag": st.column_config.TextColumn("Locus tag", help="Stable gene identifier from the source GFF3 when available; predicted ORF id for an unannotated genome."),
                "locus": st.column_config.TextColumn("Genomic locus", help="Contig and 1-based inclusive coordinates. Hover a value to see the full text."),
                "family": st.column_config.TextColumn("Family"),
                "evalue": st.column_config.NumberColumn("E-value", format="%.2e", help="Expected number of equally strong HMM matches by chance. Smaller is more significant."),
                "bits": st.column_config.NumberColumn("Bits", format="%.1f", help="HMM alignment score in bits. Larger values indicate stronger family-match support."),
                "prediction": st.column_config.TextColumn("Status", help="Only the 12 supported families can enter the complete workflow."),
            },
            key=f"candidate_editor_{out.name}", height=min(520, 38 + 35 * len(editor_rows)),
        )
        edited_rows = edited.to_dict("records") if hasattr(edited, "to_dict") else edited
        visible_ids = {r["protein_id"] for r in visible}
        selected = (selected - visible_ids) | {r["protein_id"] for r in edited_rows if r.get("run")}
        st.session_state.scan_selected_ids = selected
        eligible = [r for r in rows if r["protein_id"] in selected and r.get("family") in PIPELINE_FAMILIES]
        pending = [r for r in eligible if predictions.get(r["protein_id"], {}).get("status") != "Complete"]
        unsupported = [r for r in rows if r["protein_id"] in selected and r.get("family") not in PIPELINE_FAMILIES]
        if unsupported:
            st.warning(f"{len(unsupported)} selected candidate(s) are census-only and will not be sent to the full pipeline. Their family and locus annotations remain downloadable.")

        st.markdown('<div class="section-kicker">Step 3</div>', unsafe_allow_html=True)
        st.markdown("### Build complete reports")
        st.caption("Each selected TF receives an operator, possible-inducer and regulon report in the same format.")
        with st.expander("Advanced run settings", expanded=False):
            o1, o2 = st.columns(2)
            batch_fold = o1.checkbox("Add apo-fold analysis", value=False, key="scan_batch_fold",
                                     help="Opt-in ESMFold2 enrichment; needs BIOHUB_TOKEN.")
            batch_metalnet = o2.checkbox("Use MetalNet2 when available", value=True, key="scan_batch_metalnet",
                                         help="If unavailable, the pipeline records an abstention and continues.")
            batch_af3 = o1.checkbox("Write optional AF3 jobs", value=True, key="scan_batch_af3",
                                    help="Writes a manual hand-off file. Nothing reads it back.")
            batch_offline = o2.checkbox("Use cached data only", value=False, key="scan_batch_offline",
                                        help="Disables network retrieval.")
        already_complete = len(eligible) - len(pending)
        if already_complete:
            st.caption(f"{already_complete} selected report{'s are' if already_complete != 1 else ' is'} already complete.")
        if len(pending) > 10:
            st.info(f"Large batch: {len(pending)} reports run sequentially and may take several hours. "
                    "You can narrow the family filter or select individual TFs first.")
        run_batch = st.button(f"Build {len(pending)} complete report{'s' if len(pending) != 1 else ''}",
                              type="primary", width="stretch", disabled=not pending)
        if run_batch:
            from predictor.annotate import genome_scan as GS
            cfg = RunConfig(allow_online=not batch_offline, allow_folds=batch_fold,
                            allow_metalnet=batch_metalnet, emit_af3_jobs=batch_af3)
            predictions = dict(scan.get("predictions", {}))
            logs = dict(scan.get("logs", {}))
            progress = st.progress(0, text="Preparing selected TFs…")
            with st.status(f"Building {len(pending)} complete TF report(s)…", expanded=True) as status:
                for index, row in enumerate(pending, 1):
                    pid = row["protein_id"]
                    progress.progress((index - 1) / len(pending),
                                      text=f"{index}/{len(pending)} · {pid} · {row.get('family')}")
                    try:
                        # Production discovery's family is a routing hint: as in `tfop genome`, the
                        # pipeline makes its own family call. A census family is passed through, as before.
                        hint = None if scan.get("method") == "production" else row.get("family")
                        res, log = _run_with_log(lambda row=row, pid=pid, hint=hint: api.predict(
                            row["sequence"], name=_available_job_name(f"{out.name}_{pid}"),
                            family=hint, organism=scan.get("organism"),
                            fold=batch_fold, cfg=cfg, verbose=True,
                        ))
                        logs[pid] = log
                        predictions[pid] = {
                            "inducer": res.inducer, "primary_operator": res.primary_operator,
                            "bundle": str(res.bundle) if res.bundle else None,
                            "status": "Complete" if res.ok else "Stopped",
                            "error": None if res.ok else (res.error or "The prediction did not complete."),
                        }
                    except Exception as exc:
                        predictions[pid] = {"inducer": None, "primary_operator": None,
                                            "bundle": None, "status": "Stopped", "error": str(exc)}
                progress.progress(1.0, text="Refreshing the combined scan report…")
                cands = [GS.Candidate(**{k: v for k, v in row.items()
                                         if k in GS.Candidate.__dataclass_fields__}) for row in rows]
                GS.write_outputs(cands, out, kind=scan["kind"], n_proteins=scan["n_proteins"],
                                 source=scan.get("source", ""), predictions=predictions)
                scan["predictions"], scan["logs"] = predictions, logs
                st.session_state.last_scan = scan
                completed = sum(1 for r in pending if predictions[r["protein_id"]]["status"] == "Complete")
                status.update(label=f"Reports ready · {completed}/{len(pending)} completed",
                              state="complete" if completed == len(pending) else "error", expanded=False)
            st.rerun()

        predictions = scan.get("predictions", {})
        _render_inducer_overview(rows, predictions)
        _render_batch_report(rows, predictions)

    with st.expander("Exports & reusable files", expanded=False):
        st.caption("GFF3 is the annotation overlay; BED opens in genome browsers. Predicted labels never overwrite curated source records.")
        _scan_downloads(out)
        st.markdown(f'<div class="technical-path">Saved in {out}</div>', unsafe_allow_html=True)
    report = out / "SCAN_REPORT.html"
    if report.is_file():
        preview = st.toggle("Preview combined scan report", value=False,
                            key=f"preview_scan_report_{out.name}")
        if preview:
            st.iframe(report, width="stretch", height=900)


def _environment_panel() -> None:
    with st.expander("Environment & data", expanded=False):
        st.markdown("The core reference database is already installed. This check covers executable tools and tokens.")
        if st.button("Check environment"):
            from predictor import setup_tools
            with st.spinner("Checking the same tool paths used by the pipeline…"):
                result = setup_tools.check()
            rows = [{"engine": r["name"], "needed": "required" if r["required"] else "optional",
                     "status": "ready" if r["present"] else "missing", "purpose": r["purpose"]}
                    for r in result["rows"]]
            st.dataframe(rows, width="stretch", hide_index=True)
            missing = [r["name"] for r in result["rows"] if r["required"] and not r["present"]]
            if missing:
                st.warning("Install the missing required engines before a full prediction: " + ", ".join(missing))
            else:
                st.success("All required engines are available.")


def main() -> None:
    audit = resources.audit()
    db = audit["database"]
    st.markdown('<a class="skip-link" href="#analysis-workspace">Skip to analysis</a>', unsafe_allow_html=True)
    st.markdown(
        f'<div class="brand-row"><span class="eyebrow">Bacterial regulatory genomics</span>'
        f'<a class="project-link" href="{GITHUB_URL}" target="_blank" rel="noopener">GitHub ↗</a></div>',
        unsafe_allow_html=True,
    )
    st.markdown(f'<h1>{APP_TITLE}</h1>', unsafe_allow_html=True)
    st.markdown(
        '<p class="lead">Scan a bacterial genome for transcription factors, then turn selected candidates '
        'into ranked operator, possible-inducer and regulon reports.</p>',
        unsafe_allow_html=True,
    )
    st.markdown(
        f'<div class="topline"><span class="status-chip"><i class="dot"></i><b>Database ready</b> '
        f'v{db["release"]}</span><span class="status-chip"><b>Runs locally</b></span>'
        '<span class="status-chip"><b>Predictions need validation</b></span></div>',
        unsafe_allow_html=True,
    )

    mode = _workflow_selector()
    _scope_panel()
    if mode == "scan":
        _scan_workspace()
    else:
        _prediction_workspace()
    _environment_panel()
    st.caption("Private research software · Computational predictions require experimental validation · "
               "Reference data are versioned · Completed results are never overwritten")


main()
