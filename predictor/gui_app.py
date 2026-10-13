"""One-page web interface for TF prediction and genome/proteome census.

The UI is intentionally a thin client over the same public API and scanners used by ``tfop``. It owns
presentation and input guidance only; scientific behavior remains in ``predictor.api``.
"""
from __future__ import annotations

import contextlib
import io
import queue
import threading
import time
from pathlib import Path

import streamlit as st

from predictor import api, resources
from predictor.annotate.genome_scan import (PIPELINE_FAMILIES,
                                            POTENTIAL_METALLOREGULATOR_FAMILIES)
from predictor.config import RunConfig
from predictor.gui_support import (PIPELINE_RUN_LOCK, available_job_name,
                                   context_archive, report_preview, reserve_output_dir, safe_name)
from predictor.input_validation import read_protein


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
    return report_preview(text, bundle)


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
        "Analysis type", ("Genome scan", "Protein list", "Single TF"), default="Genome scan",
        key="workflow_choice", label_visibility="collapsed", width="stretch",
    )
    if choice == "Protein list":
        return "batch"
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
        _full_download(bundle, "bundle")
        st.caption("Individual HTML downloads require their saved folder and report assets.")
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
        st.caption("Preview: file links are disabled. Use the full bundle download to open linked files.")
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
        _full_download(bundle, "bundle")
        st.caption("Individual HTML requires the saved folder. Preview file links are disabled; use the full bundle ZIP.")
        st.download_button("Download report HTML", report.read_bytes(), file_name=f"{_safe_name(chosen)}_REPORT.html",
                           mime="text/html")
        embedded = _inline_report_assets(report.read_text(encoding="utf-8", errors="replace"), bundle)
        st.iframe(embedded, width="stretch", height=1050)


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
    st.caption(f"Saved in {out}")


def _shared_workspace(*, protein_list=False) -> None:
    """Thin GUI over the same resumable runner and output layout as CLI/API."""
    import csv
    from predictor import genome_run as GR
    st.markdown("### Protein list" if protein_list else "### Genome discovery")
    source_mode = "Protein list" if protein_list else st.radio(
        "Starting point", ("Genome accession", "Local genome FASTA", "GenBank", "Proteome"))
    with st.form("shared_entry_form"):
        accession = st.text_input("Genome accession (required)" if source_mode == "Genome accession"
                                  else "Genome accession (optional context)")
        upload = st.file_uploader("Local sequence file", type=("fna", "fa", "fasta", "faa", "gb", "gbk", "gbff", "tsv"))
        annotation = None if protein_list else st.file_uploader("GFF3 annotation", type=("gff", "gff3"))
        name = st.text_input("Output name (optional)")
        engine = _PRODUCTION if protein_list else st.radio("Discovery engine", (_PRODUCTION, _QUICK))
        offline = st.checkbox("Use cached data only", value=False)
        submitted = st.form_submit_button("Read protein list" if protein_list else "Find candidate TFs")
    key = "last_shared_batch" if protein_list else "last_shared_genome"
    if submitted:
        try:
            out = _available_output_path(_safe_name(name or accession or (upload.name if upload else "run"), "run"))
            kwargs = {"out_dir": out, "discover_only": True, "offline": offline}
            if source_mode == "Genome accession":
                kwargs["accession"] = accession.strip() or None
            else:
                if upload is None:
                    raise ValueError("Choose a local input file")
                staged = out / "uploads"
                staged.mkdir(parents=True, exist_ok=True)
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
            _run_with_log(lambda: (api.predict_many if protein_list else api.run_genome)(**kwargs))
            st.session_state[key] = str(out)
        except (ValueError, RuntimeError, OSError) as exc:
            st.error(str(exc))
    if st.session_state.get(key):
        out = Path(st.session_state[key])
        with (out / "manifest.csv").open(encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        st.dataframe(rows, hide_index=True, width="stretch")
        selected = st.multiselect("Proteins to predict", [r["run_name"] for r in rows], default=[r["run_name"] for r in rows])
        jobs = st.number_input("Parallel proteins", min_value=1, max_value=16, value=1)
        fold = st.checkbox("Add apo-fold analysis", value=False)
        if st.button("Build selected reports", disabled=not selected):
            try:
                # Every surface uses the shared runner; keep unselected discovery rows in the manifest.
                summary, _log = _run_with_log(lambda: GR.run_candidates(rows, out, jobs=int(jobs), offline=offline,
                                                       fold=fold, selected_names=selected))
                from predictor.run_cli import _run_exit_status
                if _run_exit_status(summary):
                    st.warning("Some selected predictions failed. Review the run report and per-protein logs.")
                else:
                    st.success("Selected reports finished.")
            except (ValueError, RuntimeError, OSError) as exc:
                st.error(str(exc))
        _shared_exports(out)


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
    if mode in {"scan", "batch"}:
        _shared_workspace(protein_list=mode == "batch")
    else:
        _prediction_workspace()
    _environment_panel()
    st.caption("Research software · Computational predictions require experimental validation · "
               "Reference data are versioned · Completed results are never overwritten")


main()
