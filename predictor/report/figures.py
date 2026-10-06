"""
figures.py -- the individual figures for the per-TF report. Presentation only: each `fig_*` takes
already-loaded dicts/arrays + an `out_dir`, saves an editable SVG + a PNG, and returns the PNG path
(or None when optional inputs/tools are absent). Unexpected rendering failures propagate and prevent
publication; they are never replaced by an older figure.

Reuses: `report.logo` (sequence logos + the NMI palette/rcParams); PyMOL (pymol2) for the
apo-dimer cartoon; and scipy for the TSS KDE.

Run `python -m predictor.report.figures --self-test`.
"""
from __future__ import annotations

import sys
from functools import wraps
from pathlib import Path

import numpy as np

_REPO = Path(__file__).resolve().parents[2]

import matplotlib                                                    # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                     # noqa: E402
from matplotlib.patches import Rectangle, FancyArrow, ConnectionPatch  # noqa: E402
from matplotlib.lines import Line2D                                 # noqa: E402

from predictor.report.logo import draw_logo, NMI                             # noqa: E402
from predictor.motifs.information_content import per_column as _ic            # noqa: E402

_STATUS_COLOR = {"green": "#2E9E44", "amber": "#E2A12C", "red": "#D9544D", "na": "#A8A8A8"}


# --------------------------------------------------------------------------- infra
def _save(fig, out_dir: Path, name: str) -> Path:
    """Export editable SVG/PDF and a 300-dpi report PNG; always close the canvas."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    png = out_dir / f"{name}.png"
    try:
        with matplotlib.rc_context({"svg.fonttype": "none", "pdf.fonttype": 42}):
            for ext in ("svg", "pdf", "png"):
                fig.savefig(out_dir / f"{name}.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
    finally:
        plt.close(fig)
    return png


def _safe(fn):
    """Close new canvases and propagate unexpected failures instead of hiding them."""
    @wraps(fn)
    def wrapper(*a, **kw):
        before = set(plt.get_fignums())
        try:
            return fn(*a, **kw)
        except Exception as exc:
            raise RuntimeError(f"report figure {fn.__name__} failed: {exc}") from exc
        finally:
            for number in set(plt.get_fignums()) - before:
                plt.close(number)
    return wrapper


def _ic_of(pwm) -> np.ndarray:
    pwm = np.asarray(pwm, dtype=float)
    return _ic(pwm)


# --------------------------------------------------------------------------- §1 apo dimer
@_safe
def fig_apo_dimer(struct_d: dict, out_dir: Path) -> Path | None:
    """A minimal PyMOL cartoon of the apo homodimer (chains coloured), b-factor = pLDDT. Returns None when
    there is no on-disk dimer fold or PyMOL is unavailable (the report then shows the pLDDT/QC text only)."""
    cif = (struct_d or {}).get("dimer_cif")
    if not cif or not Path(cif).exists():
        return None
    try:
        import pymol2
    except Exception:
        return None
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    png = out_dir / "apo_dimer.png"
    with pymol2.PyMOL() as p:
        c = p.cmd
        c.load(cif, "dimer")
        c.hide("everything")
        c.show("cartoon")
        c.bg_color("white")
        c.set("ray_opaque_background", 1)
        c.set("cartoon_transparency", 0.0)
        c.color("teal", "chain A")
        c.color("wheat", "chain B")
        c.orient()
        c.ray(1400, 1000)
        c.png(str(png), dpi=200)
    return png if png.exists() else None


# --------------------------------------------------------------------------- §2 operator track
@_safe
def fig_operator_track(rows: list, out_dir: Path) -> Path | None:
    """Lollipop of every genome operator (x=dyad, height=score); the autoregulatory operator highlighted."""
    rows = [r for r in rows if r.get("dyad") is not None and r.get("score") is not None]
    fig, ax = plt.subplots(figsize=(7.4, 2.2))
    if not rows:
        ax.text(0.5, 0.5, "no operators found", ha="center", va="center", transform=ax.transAxes,
                color=NMI["neutral_dark"]); ax.axis("off")
        return _save(fig, out_dir, "operator_track")
    ys = [r["score"] for r in rows]
    for r in rows:
        auto = r.get("is_autoregulatory")
        col = NMI["delta_down"] if auto else NMI["baseline_mid"]
        x = r["dyad"] / 1e6
        ax.vlines(x, 0, r["score"], color=col, lw=(1.8 if auto else 0.8), zorder=3 if auto else 1)
        ax.plot(x, r["score"], "o", ms=(5 if auto else 3), color=col, zorder=4 if auto else 2)
    auto_rows = [r for r in rows if r.get("is_autoregulatory")]
    if auto_rows:
        a = auto_rows[0]
        ax.annotate(f"autoregulatory\n{a.get('regulated_gene') or ''}".strip(), (a["dyad"] / 1e6, a["score"]),
                    textcoords="offset points", xytext=(0, 8), ha="center", fontsize=6.5,
                    color=NMI["delta_down"], fontweight="bold")
    ax.set_xlabel("genome position (Mb)", fontsize=7.5)
    ax.set_ylabel("site score (bits)", fontsize=7.5)
    ax.set_ylim(0, max(ys) * 1.25)
    ax.set_title(f"{len(rows)} putative operators genome-wide (autoregulatory highlighted)",
                 fontsize=8, loc="left")
    return _save(fig, out_dir, "operator_track")


# --------------------------------------------------------------------------- §2 operator logos
@_safe
def fig_operator_logo(logo: dict, out_dir: Path, *, name: str, title: str) -> Path | None:
    """Sequence logo from a stored logo dict ({pwm, per_col_ic?, consensus,...}). Skips empty logos."""
    pwm = (logo or {}).get("pwm")
    if not pwm:
        return None
    pwm = np.asarray(pwm, dtype=float)
    if pwm.ndim != 2 or pwm.shape[0] != 4 or pwm.shape[1] == 0:
        return None
    ic = logo.get("per_col_ic")
    ic = np.asarray(ic, dtype=float) if ic else _ic_of(pwm)
    fig, ax = plt.subplots(figsize=(max(3.0, 0.28 * pwm.shape[1]), 1.7))
    draw_logo(ax, pwm, ic)
    ax.set_xlabel("position (bp)", fontsize=7.5)                    # override the shared default label
    cons = (logo or {}).get("consensus") or ""
    # Seed logos record `n_effective`, the reported logo `n_kept`; only the multi-source logos carry
    # `n_seqs`. Reading `n_seqs` alone printed "n=?" on the two logos the report leads with.
    n = next((logo.get(k) for k in ("n_seqs", "n_kept", "n_effective") if logo.get(k) is not None), None)
    n_txt = f", n = {n}" if n is not None else ""
    ax.set_title(f"{title}  ·  {cons}  ({float(ic.sum()):.1f} bits{n_txt})", fontsize=8, loc="left")
    return _save(fig, out_dir, name)


# --------------------------------------------------------------------------- §3 inducer evidence
@_safe
def fig_inducer_evidence(ev: dict, out_dir: Path) -> Path | None:
    """One row per inducer source with its ligand call + key evidence, plus the agreement headline."""
    rows = ev.get("rows") or []
    fig, ax = plt.subplots(figsize=(7.4, max(1.4, 0.5 + 0.42 * max(1, len(rows)))))
    ax.axis("off")
    top = ev.get("top"); raw_top = ev.get("raw_top", top); ag = ev.get("agreement")
    gate = ev.get("coordination_gate")
    head = f"inferred inducer:  {top or '—'}"
    if ag is not None:
        head += f"    (source agreement {ag:.2f}"
        head += f", coordination gate {'yes' if gate else 'no' if gate is False else '?'})"
    ax.text(0.0, 1.0, head, fontsize=9, fontweight="bold", va="top", color=NMI["baseline_dark"])
    if not rows:
        ax.text(0.0, 0.5, "no inducer inference (offline / no NCBI)", fontsize=8,
                color=NMI["neutral_dark"], va="top")
        return _save(fig, out_dir, "inducer_evidence")
    y = 0.80
    dy = 0.80 / max(1, len(rows))
    for r in rows:
        src = r.get("source", "?")
        lig = r.get("ligand")
        conf = r.get("confidence")
        ev_d = r.get("evidence") or {}
        ev_s = ", ".join(f"{k}={v}" for k, v in list(ev_d.items())[:4])
        col = NMI["delta_up"] if (lig and lig == raw_top) else NMI["neutral_dark"]
        ax.add_patch(Rectangle((0.0, y - dy * 0.46), 0.16, dy * 0.7, transform=ax.transAxes,
                               facecolor=NMI["bg_lilac"], edgecolor="none"))
        ax.text(0.01, y, src, fontsize=7.5, va="center", fontweight="bold")
        lig_s = str(lig or "—")
        lig_s = lig_s if len(lig_s) <= 22 else lig_s[:21] + "…"     # long ligand names ran into the evidence
        ax.text(0.18, y, lig_s + (f"  (conf {conf:.2f})" if isinstance(conf, (int, float)) else ""),
                fontsize=7.5, va="center", color=col)
        ev_s = ev_s if len(ev_s) <= 95 else ev_s[:94] + "…"
        ax.text(0.46, y, ev_s, fontsize=6.8, va="center", color=NMI["neutral_dark"])
        y -= dy
    return _save(fig, out_dir, "inducer_evidence")


# --------------------------------------------------------------------------- §3 operator vs TSS
@_safe
def _offset_panel(ax, vals: list, *, anchor: str, xlabel: str, empty: str) -> None:
    """One signed-distance histogram (+KDE) against a named anchor at x=0."""
    if not vals:
        ax.text(0.5, 0.5, empty, ha="center", va="center", transform=ax.transAxes,
                color=NMI["neutral_dark"]); ax.axis("off")
        return
    rng = (min(min(vals), 0) - 10, max(max(vals), 0) + 10)       # always include the anchor at x=0
    nb = max(8, min(40, len(vals)))
    _, edges, _ = ax.hist(vals, bins=nb, range=rng, color=NMI["baseline_soft"],
                          edgecolor=NMI["baseline_mid"], lw=0.5)
    if len(vals) >= 5 and len(set(vals)) > 1:
        try:
            from scipy.stats import gaussian_kde
            xs = np.linspace(rng[0], rng[1], 200)
            # Same axis, same units: density x n x bin width = expected sites per bin. A separate
            # twin axis autoscaled its own margin, so zero density was drawn above the baseline.
            ax.plot(xs, gaussian_kde(vals)(xs) * len(vals) * (edges[1] - edges[0]),
                    color=NMI["delta_down"], lw=1.3)
        except Exception:
            pass
    ax.set_ylim(bottom=0)
    ax.axvline(0, color=NMI["neutral_dark"], lw=1.0, ls="--")
    ax.text(0, 1.01, anchor, transform=ax.get_xaxis_transform(), fontsize=7, ha="center", va="bottom",
            color=NMI["neutral_dark"])
    ax.set_xlabel(xlabel, fontsize=7.5)
    ax.set_ylabel("sites", fontsize=7.5)


def fig_operator_vs_tss(offsets: list, out_dir: Path, *, cds_rows: list | None = None) -> Path | None:
    """Where the operators sit, against BOTH anchors: the predicted TSS and the assigned CDS start.

    The two panels are not interchangeable and the caption says so.

      * **TSS (left)** -- the sigma70 promoter is predicted inside a window CENTRED ON THE OPERATOR and
        then selected for maximal overlap with that same operator (`signals.promoter.classify_mode`), so
        the narrowness of this distribution is partly a property of the selection rule rather than of
        biology. It is the right panel for asking *which promoter elements does the operator cover*, and
        the wrong one for asking *is the placement better than chance*.
      * **CDS (right)** -- distance to the 5' end of the assigned gene, read from the genome annotation.
        Operator-INDEPENDENT, so it is the anchor a null model can be built against; that is exactly why
        the manuscript's F3 was rebuilt around it (`analysis/fig_operator_cds.py`).

    `cds_rows` is optional: without it this renders the single TSS panel it always did.
    """
    vals = [o["offset"] for o in offsets if o.get("offset") is not None]
    cvals = [o["offset"] for o in (cds_rows or []) if o.get("offset") is not None]
    if not cvals:
        fig, ax = plt.subplots(figsize=(5.2, 2.4))
        _offset_panel(ax, vals, anchor="TSS", empty="no TSS calls",
                      xlabel="operator centre − TSS (bp)   [<0 upstream · >0 downstream]")
        if vals:
            ax.set_title(f"operator position relative to the TSS (n={len(vals)})", fontsize=8, loc="left")
        return _save(fig, out_dir, "operator_vs_tss")

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.6))
    _offset_panel(axes[0], vals, anchor="TSS", empty="no TSS calls",
                  xlabel="site centre relative to TSS (bp; <0 upstream)")
    # `operator_to_gene_bp` is CDS start minus site centre on the gene's strand, i.e. positive UPSTREAM.
    # Labelled "operator centre − CDS start" it read with the opposite sign to the TSS panel beside it.
    _offset_panel(axes[1], cvals, anchor="start codon", empty="no CDS assignment",
                  xlabel="site centre upstream of CDS start (bp)")
    axes[0].set_title(f"vs predicted TSS (n = {len(vals)}); anchor depends on the site",
                      fontsize=7.5, loc="left", pad=12)
    axes[1].set_title(f"vs annotated CDS start (n = {len(cvals)}); independent anchor",
                      fontsize=7.5, loc="left", pad=12)
    fig.tight_layout()
    return _save(fig, out_dir, "operator_vs_tss")


# --------------------------------------------------------------------------- genome context (neighborhood + promoter-zoom callout)
def _draw_neighborhood(ax, p: dict, *, by: float = 0.62, op_label: bool = True) -> None:
    """One gene-neighborhood lane: flanking genes as strand-oriented arrows, the operator site boxed. `by`
    is the gene-line height (lowered when a zoom callout sits above); labels stagger below it."""
    genes, op = p["genes"], p.get("operator")
    xs = [g["start"] for g in genes] + [g["end"] for g in genes] + (list(op) if op else [])
    lo, hi = min(xs), max(xs)
    pad = max(200, int((hi - lo) * 0.06))
    lo, hi = lo - pad, hi + pad
    ax.set_xlim(lo, hi); ax.set_ylim(0, 1); ax.set_yticks([]); ax.set_xticks([])
    ax.plot([lo, hi], [by, by], color=NMI["neutral_dark"], lw=0.8, zorder=1)
    for gi, g in enumerate(sorted(genes, key=lambda g: g["start"])):
        x0, x1 = g["start"], g["end"]
        d0, dx = (x1, x0 - x1) if g.get("strand") == "-" else (x0, x1 - x0)
        ax.add_patch(FancyArrow(d0, by, dx, 0, width=0.11, head_width=0.22,
                                head_length=min((hi - lo) * 0.02, abs(dx) * 0.4 + 1),
                                length_includes_head=True, color=NMI["baseline_soft"], lw=0, zorder=2))
        ly = by - 0.11 - 0.10 * (gi % 2)              # stagger labels so dense neighbours don't overlap
        ax.text((x0 + x1) / 2, ly, g.get("name", ""), fontsize=5.6, ha="center", va="top",
                color=NMI["neutral_dark"])
    if op:
        ax.add_patch(Rectangle((op[0], by - 0.06), max((hi - lo) * 0.004, op[1] - op[0]), 0.12,
                               color=NMI["ours_large"], ec=NMI["delta_up"], lw=1.1, zorder=6))
        if op_label:
            ax.text((op[0] + op[1]) / 2, by + 0.12, "operator", fontsize=6, ha="center",
                    color=NMI["delta_up"])
    ttl = p.get("title", "")
    if p.get("subtitle"):
        ttl += "   ·   " + p["subtitle"]
    ax.set_title(ttl, fontsize=7.3, loc="left", color=NMI["baseline_dark"])
    ax.text(1.0, 1.02, f"{lo/1e6:.3f}–{hi/1e6:.3f} Mb", transform=ax.transAxes, ha="right", va="bottom",
            fontsize=5.4, color=NMI["neutral_dark"], style="italic")     # above the frame; no label clash
    for sp in ("top", "right", "left", "bottom"):
        ax.spines[sp].set_visible(False)


def _promoter_glyph(ax, z: dict, x0: float, y0: float, w: float, h: float) -> None:
    """Draw the compact promoter schematic (-35/spacer/-10 boxes, TSS, operator footprint) into the
    axes-fraction rectangle [x0,x0+w]x[y0,y0+h] -- the content of a zoom callout box."""
    el, op = z["elements"], z["op"]
    xs = [op[0], op[1], 0] + [v for e in el.values() for v in e]
    plo, phi = min(xs) - 3, max(xs) + 3
    padf = 0.07
    T = ax.transAxes

    def mx(v):
        return x0 + w * padf + (v - plo) / (phi - plo) * w * (1 - 2 * padf)

    ly = y0 + h * 0.44
    ax.plot([x0 + w * padf, x0 + w - w * padf], [ly, ly], transform=T, color=NMI["neutral_dark"],
            lw=0.6, zorder=8)
    ecol = {"-35": NMI["baseline_mid"], "-10": NMI["baseline_mid"], "spacer": NMI["baseline_soft"]}
    bh = h * 0.15
    for name in ("-35", "spacer", "-10"):
        if name not in el:
            continue
        a, b = mx(el[name][0]), mx(el[name][1])
        ax.add_patch(Rectangle((a, ly - bh / 2), max(0.001, b - a), bh, transform=T, color=ecol[name],
                               ec="none", zorder=9))
        # zorder above the inset's white box (6); at the default 3 these labels were never visible.
        ax.text((a + b) / 2, ly - bh * 1.2 - 0.004, name, transform=T, fontsize=5.5, ha="center", va="top",
                color=NMI["neutral_dark"], zorder=11)
    xt = mx(0)
    ax.plot([xt, xt], [ly - h * 0.22, ly + h * 0.24], transform=T, color=NMI["delta_up"], lw=0.9,
            ls="--", zorder=9)
    ax.text(xt, ly + h * 0.25, "TSS", transform=T, fontsize=5.5, ha="center", va="bottom",
            color=NMI["delta_up"], zorder=11)
    xoa, xob = mx(op[0]), mx(op[1])
    ax.add_patch(Rectangle((xoa, ly - bh * 1.15), max(0.002, xob - xoa), bh * 2.3, transform=T,
                           facecolor=NMI["ours_large"], ec=NMI["delta_down"], lw=0.7, alpha=0.5, zorder=10))


def _add_zoom_callout(ax, p: dict) -> None:
    """Place a SMALL promoter-zoom box ABOVE the neighborhood, over the operator, with leader lines from the
    operator up to the box -- a magnifier callout showing which promoter section the zoom is of."""
    z = p.get("zoom") or {}
    op = p.get("operator")
    if not z.get("elements") or not op:
        return
    lo, hi = ax.get_xlim()
    op_frac = min(max(((op[0] + op[1]) / 2 - lo) / (hi - lo), 0.0), 1.0)
    zw, zh, zy0 = 0.38, 0.30, 0.57                    # smaller than the lane; sits above it
    zx0 = min(max(op_frac - zw / 2, 0.02), 0.98 - zw)
    ax.add_patch(Rectangle((zx0, zy0), zw, zh, transform=ax.transAxes, facecolor="white",
                           ec=NMI["neutral_dark"], lw=0.7, zorder=6))
    _promoter_glyph(ax, z, zx0, zy0, zw, zh)
    occ = ", ".join(z.get("occludes") or []) or "—"
    ax.text(zx0 + zw / 2, zy0 + zh + 0.015, f"promoter zoom · occludes {occ}", transform=ax.transAxes,
            fontsize=5.2, ha="center", va="bottom", color=NMI["baseline_dark"])
    for gx, fx in ((op[0], zx0), (op[1], zx0 + zw)):  # leader lines: operator -> zoom box bottom corners
        con = ConnectionPatch(xyA=(gx, 0.28), coordsA="data", axesA=ax,
                              xyB=(fx, zy0), coordsB="axes fraction", axesB=ax,
                              color=NMI["neutral_mid"], lw=0.5, zorder=5)
        ax.add_artist(con)


@_safe
def fig_gene_neighborhoods(panels: list, out_dir: Path) -> Path | None:
    """Standalone gene-neighborhood figure (one lane per panel). Retained for reuse; the report embeds the
    combined `fig_operator_context` instead."""
    panels = [p for p in (panels or []) if p.get("genes")]
    if not panels:
        return None
    fig, axes = plt.subplots(len(panels), 1, figsize=(7.6, 1.3 * len(panels) + 0.3), squeeze=False)
    for ax, p in zip(axes[:, 0], panels):
        _draw_neighborhood(ax, p)
    fig.tight_layout()
    return _save(fig, out_dir, "gene_neighborhoods")


@_safe
def fig_operator_context(panels: list, out_dir: Path) -> Path | None:
    """The genome context of the top natural operators: one gene-neighborhood lane per operator, each with a
    SMALL promoter-zoom callout above it (operator footprint over the -35/spacer/-10 boxes vs TSS), joined by
    leader lines showing which promoter section the zoom is of. Placed as a §3 detail."""
    panels = [p for p in (panels or []) if p.get("genes")]
    if not panels:
        return None
    # panels WITH a resolvable promoter zoom get a taller lane (neighborhood in the lower band + callout
    # above); panels WITHOUT one render as a compact normal lane (no empty upper space).
    has_zoom = [bool((p.get("zoom") or {}).get("elements")) for p in panels]
    hr = [2.0 if z else 1.0 for z in has_zoom]
    h = sum(2.15 if z else 1.2 for z in has_zoom) + 0.3
    fig, axes = plt.subplots(len(panels), 1, figsize=(7.6, h), squeeze=False,
                             gridspec_kw={"height_ratios": hr})
    for ax, p, z in zip(axes[:, 0], panels, has_zoom):
        if z:
            _draw_neighborhood(ax, p, by=0.22, op_label=False)   # neighborhood low; zoom callout above
            _add_zoom_callout(ax, p)
        else:
            _draw_neighborhood(ax, p, by=0.5)                    # normal centred lane, no empty band
    fig.subplots_adjust(hspace=0.55)
    return _save(fig, out_dir, "operator_context")


# --------------------------------------------------------------------------- §2 operator PSSM / scoring
@_safe
def fig_operator_pssm(logo: dict, out_dir: Path) -> Path | None:
    """How an operator site is SCORED: the log-odds PSSM (4xW, ACGT vs a uniform background, in bits) as a
    heatmap with the consensus path boxed, over a per-column information (bits) strip. A site's score is the
    sum of the boxed log-odds along its bases; the p-value comes from the exact null score distribution."""
    pwm = (logo or {}).get("pwm")
    if not pwm:
        return None
    P = np.asarray(pwm, dtype=float)
    if P.shape[0] != 4:
        P = P.T
    W = P.shape[1]
    log_odds = np.log2(np.clip(P, 1e-6, None) / 0.25)        # bits vs uniform 0.25 background
    ic = _ic_of(P)
    cons = "".join("ACGT"[i] for i in P.argmax(axis=0))
    fig, (axic, axh) = plt.subplots(2, 1, figsize=(max(4.6, 0.34 * W + 1.2), 3.1),
                                    gridspec_kw={"height_ratios": [1, 3]}, sharex=True)
    axic.bar(range(W), ic, color=NMI["baseline_mid"], width=0.9)
    axic.set_ylabel("bits", fontsize=7)
    axic.set_title("per-column information (bits) — how sharp each position is", fontsize=7.5, loc="left")
    axic.tick_params(labelsize=6)
    for sp in ("top", "right"):
        axic.spines[sp].set_visible(False)
    im = axh.imshow(log_odds, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    axh.set_yticks(range(4)); axh.set_yticklabels(list("ACGT"), fontsize=7)
    axh.set_xticks(range(W)); axh.set_xticklabels(list(cons), fontsize=6.3)
    axh.set_xlabel("operator position (consensus base labelled)  —  site score = Σ log-odds along the "
                   "bound bases", fontsize=6.8)
    for j, ch in enumerate(cons):                             # box the consensus (max-scoring) path
        axh.add_patch(Rectangle((j - 0.5, "ACGT".index(ch) - 0.5), 1, 1, fill=False,
                                ec=NMI["delta_up"], lw=1.2))
    cb = fig.colorbar(im, ax=axh, fraction=0.025, pad=0.01)
    cb.set_label("log-odds (bits)", fontsize=6.5); cb.ax.tick_params(labelsize=6)
    fig.tight_layout()
    return _save(fig, out_dir, "operator_pssm")


# --------------------------------------------------------------------------- §3 operator footprint vs TSS
@_safe
def fig_operator_footprint(offsets: list, out_dir: Path) -> Path | None:
    """The FULL operator footprint (not just its centre) relative to the TSS: one horizontal bar per operator
    from its start to end edge, x=0 at the TSS, coloured by regulation mode, annotated with which promoter
    elements it occludes. Shows, e.g., that an ArsR operator physically covers the -10 box and the TSS."""
    rows = [o for o in (offsets or []) if o.get("span")]
    fig, ax = plt.subplots(figsize=(6.4, max(1.7, 0.30 * max(1, len(rows)) + 0.7)))
    if not rows:
        ax.text(0.5, 0.5, "no operator footprint (needs TSS calls)", ha="center", va="center",
                transform=ax.transAxes, color=NMI["neutral_dark"]); ax.axis("off")
        return _save(fig, out_dir, "operator_footprint")
    rows = sorted(rows, key=lambda o: -(o.get("score") or 0))[:20]
    modecol = {"elongated-spacer-overlap": NMI["baseline_mid"],
               "promoter-core-overlap": NMI["delta_down"]}
    for i, o in enumerate(rows):
        s, e = o["span"]
        y = len(rows) - i
        col = modecol.get(o.get("mode"), NMI["neutral_dark"])
        ax.plot([s, e], [y, y], color=col, lw=4.2, solid_capstyle="butt", zorder=3)
        occ = ",".join(o.get("occludes") or [])
        lab = (o.get("gene") or "") + (f"  [{occ}]" if occ else "")
        ax.text(e + (max(1, (max(r["span"][1] for r in rows) - min(r["span"][0] for r in rows)) * 0.01)),
                y, lab, fontsize=6, va="center", color=NMI["neutral_dark"])
    ax.axvline(0, color=NMI["delta_up"], lw=1.2, ls="--", zorder=2)
    ax.text(0, len(rows) + 0.6, "TSS", fontsize=7, color=NMI["delta_up"], ha="center")
    ax.set_ylim(0.3, len(rows) + 1.2); ax.set_yticks([])
    ax.set_xlabel("operator footprint relative to TSS (bp)   [<0 upstream · >0 downstream]", fontsize=7.5)
    ax.set_title(f"operator footprints vs TSS (n={len(rows)}; each bar = the whole operator span)",
                 fontsize=8, loc="left")
    ax.legend(handles=[Line2D([0], [0], color=NMI["baseline_mid"], lw=4, label="elongated spacer"),
                       Line2D([0], [0], color=NMI["delta_down"], lw=4, label="promoter core")],
              fontsize=6, loc="lower right", frameon=False)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    return _save(fig, out_dir, "operator_footprint")


# --------------------------------------------------------------------------- self-test
def _demo() -> None:
    import tempfile
    td = Path(tempfile.mkdtemp(prefix="figtest_"))
    pwm = np.array([[0.7, 0.1, 0.1, 0.85], [0.1, 0.1, 0.1, 0.05],
                    [0.1, 0.7, 0.1, 0.05], [0.1, 0.1, 0.7, 0.05]])
    assert fig_operator_logo({"pwm": pwm.tolist(), "consensus": "ACGT", "n_seqs": 8}, td,
                             name="logo_test", title="test") is not None
    assert fig_operator_logo({"pwm": None}, td, name="empty", title="x") is None       # empty -> None
    rows = [{"dyad": 100, "score": 9.0, "is_autoregulatory": True, "regulated_gene": "geneA"},
            {"dyad": 900, "score": 4.0, "is_autoregulatory": False, "regulated_gene": None}]
    assert fig_operator_track(rows, td) is not None
    assert fig_operator_track([], td) is not None                                       # empty -> placeholder fig
    ev = {"top": "Zn2+", "agreement": 0.33, "coordination_gate": True,
          "rows": [{"source": "ssn_cluster", "ligand": "Zn2+", "confidence": 1.0,
                    "evidence": {"cluster": "ArsR_c2"}}]}
    assert fig_inducer_evidence(ev, td) is not None
    offs = [{"offset": v} for v in (-19, -22, -25, -18, -30, -12, -20, -24)]
    assert fig_operator_vs_tss(offs, td) is not None
    assert fig_operator_vs_tss([], td) is not None                                      # placeholder
    # NEW: gene neighborhoods, operator PSSM, operator footprint
    panels = [{"title": "TF — autoregulatory locus", "subtitle": "own promoter",
               "genes": [{"name": "tfR", "start": 500, "end": 800, "strand": "+"},
                         {"name": "targ", "start": 200, "end": 480, "strand": "-"}],
               "operator": [485, 505]},
              {"title": "target: copA operon", "subtitle": "q=1e-9 · copA",
               "genes": [{"name": "copA", "start": 1200, "end": 3600, "strand": "+"}], "operator": [1150, 1180]}]
    assert fig_gene_neighborhoods(panels, td) is not None
    assert fig_gene_neighborhoods([], td) is None                                       # nothing -> None
    # combined context = neighborhood lanes, each with its own promoter-zoom callout above
    zoom = {"op": [-16, 4], "elements": {"-35": [-35, -29], "spacer": [-29, -12], "-10": [-12, -6]},
            "occludes": ["-10", "TSS"], "mode": "promoter-core-overlap", "spacer_len": 17}
    pz = [dict(panels[0], zoom=zoom), dict(panels[1], zoom={})]                          # one has a zoom, one doesn't
    assert fig_operator_context(pz, td) is not None
    assert fig_operator_context([], td) is None                                         # nothing -> None
    assert fig_operator_pssm({"pwm": pwm.tolist(), "consensus": "ACGT"}, td) is not None
    assert fig_operator_pssm({"pwm": None}, td) is None
    fp = [{"span": [-29, -9], "occludes": ["-10", "TSS"], "mode": "promoter-core-overlap",
           "gene": "czrB", "score": 17.8},
          {"span": [-2, 18], "occludes": ["spacer"], "mode": "elongated-spacer-overlap",
           "gene": "cueR", "score": 12.0}]
    assert fig_operator_footprint(fp, td) is not None
    assert fig_operator_footprint([{"offset": -19}], td) is not None                    # no span -> placeholder
    # the CDS panel is additive: with no CDS distances the figure is the TSS panel it always was
    assert fig_operator_vs_tss([{"offset": -19}], td, cds_rows=[{"offset": -64}]) is not None
    # absent-input degradation -> None
    assert fig_apo_dimer({"dimer_cif": None}, td) is None
    print(f"OK: figures render (logo/track/neighborhoods/pssm/inducer/tss+cds/footprint) + degrade "
          f"to None on absent inputs. -> {td}")


if __name__ == "__main__":
    if "--self-test" in sys.argv or True:
        _demo()
