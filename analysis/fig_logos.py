#!/usr/bin/env python
"""F4 — operator sequence logos at three levels, in SEPARATE figures:
   F4a per family, F4b per SSN cluster, F4c per TF.
Degenerate logos are handled honestly: fixed 0-2 bit axis + a total-information-content (IC) readout and
a colour chip per logo, so low conservation is QUANTIFIED, not hidden. A family logo is expected to be
degenerate (a family pools heterogeneous cluster operators) -- that is itself the finding."""
from __future__ import annotations
import csv, json, sys
from collections import defaultdict
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig_style import save, require, focus_family, FAM, BASE, INK, NEUTRAL_M, NEUTRAL_L
import matplotlib.pyplot as plt


# One definition of which candidates this run is about; see run_paths.MANIFEST for why there
# is no local fallback.
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent))
import run_paths as RP  # noqa: E402


def _manifest_path():
    """The manifest for THIS run. Defined once, in run_paths."""
    return RP.MANIFEST

PROJ = Path(__file__).resolve().parent.parent
# Resolved from the run tag, never from a literal path: a figure drawn from an archived
# run into a current report is the expensive kind of error, because the report still
# looks finished.
JOBS = RP.JOBS
FIG = RP.FIGURES
FAMS = ("ArsR", "MerR", "Fur")


#: REMOVED: the family must come from the bundle, not the run name. `fig_style.focus_family`
#: reads `dossier["family"]`; the run name encodes what DISCOVERY guessed from a partial
#: profile match, which the full-sequence classifier can and does contradict.


def _collect():
    """{name: (family, cluster, gene, pwm)} for the focus TFs that have a genomic_logo PWM."""
    man = {r["run_name"]: r for r in csv.DictReader(_manifest_path().open(encoding="utf-8"))}
    out = {}
    for name, m in man.items():
        if not (JOBS / name / "dossier.json").exists():
            continue
        j = json.loads((JOBS / name / "dossier.json").read_text(encoding="utf-8"))
        fam = focus_family(j.get("family"))
        if fam is None:
            continue
        pwm = (j.get("genomic_logo") or {}).get("pwm")
        if pwm:
            out[name] = (fam, j.get("ssn_cluster"), m.get("gene") or name.split("__")[-1], np.asarray(pwm, float))
    return out


def _col_ic(p):
    p = p / p.sum(0, keepdims=True).clip(1e-9)
    return np.log2(4) + (p * np.log2(p.clip(1e-9))).sum(0)


def _peak_align(pwms, width):
    half = width // 2
    acc = np.zeros((4, width)); n = np.zeros(width)
    for p in pwms:
        peak = int(np.argmax(_col_ic(p)))
        for j in range(p.shape[1]):
            col = half + (j - peak)
            if 0 <= col < width:
                acc[:, col] += p[:, j] / p[:, j].sum().clip(1e-9); n[col] += 1
    n[n == 0] = 1
    return acc / n


def _draw_logo(ax, pwm, title, fam):
    """IC-scaled letter stack (bits). Returns total IC."""
    try:
        import logomaker, pandas as pd
        p = pwm / pwm.sum(0, keepdims=True).clip(1e-9)
        ic = _col_ic(p)
        mat = (p * ic).T
        df = pd.DataFrame(mat, columns=["A", "C", "G", "T"])
        logomaker.Logo(df, ax=ax, color_scheme={b: BASE[b] for b in "ACGT"}, show_spines=False)
        total_ic = float(ic.sum())
    except Exception:
        ax.text(0.5, 0.5, "logo n/a", ha="center", va="center", fontsize=6, transform=ax.transAxes)
        total_ic = 0.0
    ax.set_ylim(0, 2); ax.set_xticks([]); ax.set_yticks([0, 1, 2])
    ax.tick_params(axis="y", labelsize=5, length=2)
    ax.set_ylabel("bits", fontsize=5.5, labelpad=1)
    # title (family-coloured) + IC readout + degeneracy chip
    ax.set_title(title, fontsize=7, color=FAM.get(fam, INK), loc="left", pad=2)
    chip = "#2E9E44" if total_ic >= 8 else ("#E28E2C" if total_ic >= 4 else NEUTRAL_L)
    ax.text(1.0, 1.14, f"IC {total_ic:.0f} bits", transform=ax.transAxes, ha="right", va="center",
            fontsize=5.6, color=INK)
    ax.scatter([1.0], [1.14], transform=ax.transAxes, s=26, color=chip, clip_on=False, zorder=5,
               edgecolor="white", linewidth=0.5)
    ax.text(1.055, 1.14, "", transform=ax.transAxes)
    return total_ic


def fig_family(data):
    # `np.median([])` is NaN, and `int(NaN)` raised "cannot convert float NaN to integer" -- a
    # traceback that named numpy rather than the absent bundles that actually caused it.
    require(data, "operator PWMs from this run's bundles",
            "the batch stage (run_phase1.py)",
            "Check TFOP_RUN_TAG and that jobs/ holds this run's bundles.")
    fig, axes = plt.subplots(len(FAMS), 1, figsize=(4.6, 3.2))
    for ax, fam in zip(axes, FAMS):
        pwms = [pw for (f, c, g, pw) in data.values() if f == fam]
        require(pwms, f"operator PWMs for family {fam}", "the batch stage (run_phase1.py)")
        width = int(np.median([p.shape[1] for p in pwms]))
        _draw_logo(ax, _peak_align(pwms, width), f"{fam}   (n={len(pwms)} TFs, pooled)", fam)
    fig.suptitle("Operator logo per FAMILY", fontsize=9, fontweight="bold", x=0.02, ha="left", y=1.0)
    fig.text(0.02, -0.03, "A family pools heterogeneous cluster operators, so its logo is expected to be "
             "low-IC (degenerate) — the signal lives at the cluster level (F4b). Chip: green ≥8 bits, "
             "amber 4–8, grey <4.", fontsize=5.4, color=NEUTRAL_M, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save(fig, FIG / "F4a_operator_logo_per_family")
    print("  F4a_operator_logo_per_family.svg/.png")


def fig_cluster(data):
    clu = defaultdict(list); clu_fam = {}
    for (f, c, g, pw) in data.values():
        if c:
            clu[c].append(pw); clu_fam[c] = f
    items = [(c, ps) for c, ps in clu.items() if len(ps) >= 2]
    # sort: family, then descending pooled IC (informative first)
    def pic(ps):
        w = int(np.median([p.shape[1] for p in ps])); return float(_col_ic(_peak_align(ps, w)).sum())
    items.sort(key=lambda cp: (FAMS.index(clu_fam[cp[0]]), -pic(cp[1])))
    rows = len(items)
    if not rows:
        # A single-genome run can have no SSN cluster carrying two regulators, and `subplots(0, 1)`
        # raised. Draw a labelled placeholder (as fig_cluster_motifs does for its F4b) so figure QA
        # finds a well-formed file and the report says why the panel is empty.
        fig = plt.figure(figsize=(4.8, 1.6))
        fig.text(0.02, 0.62, "Operator logo per SSN CLUSTER — not applicable to this run",
                 fontsize=9, fontweight="bold", ha="left")
        fig.text(0.02, 0.30, "No SSN cluster holds two or more of this run's regulators, so there is "
                 "no cluster to pool.", fontsize=6.5, color=NEUTRAL_M, ha="left")
        save(fig, FIG / "F4b_operator_logo_per_cluster")
        print("  F4b_operator_logo_per_cluster SKIPPED (no cluster with >=2 TFs); placeholder drawn")
        return
    fig, axes = plt.subplots(rows, 1, figsize=(4.8, 0.82 * rows))
    axes = np.atleast_1d(axes)
    for ax, (c, ps) in zip(axes, items):
        width = int(np.median([p.shape[1] for p in ps]))
        _draw_logo(ax, _peak_align(ps, width), f"{c}   (n={len(ps)})", clu_fam[c])
    fig.suptitle("Operator logo per SSN CLUSTER (peak-IC aligned, ≥2 TFs)", fontsize=9, fontweight="bold",
                 x=0.02, ha="left", y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    save(fig, FIG / "F4b_operator_logo_per_cluster")
    print("  F4b_operator_logo_per_cluster.svg/.png")


def fig_tf(data):
    # grouped by family then cluster; 3 columns of compact per-TF logos
    items = sorted(data.items(), key=lambda kv: (FAMS.index(kv[1][0]), kv[1][1] or "zz", kv[1][2]))
    n = len(items); ncol = 3; nrow = -(-n // ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(7.2, 0.66 * nrow))
    axes = np.atleast_2d(axes)
    for k, (name, (f, c, g, pw)) in enumerate(items):
        ax = axes[k // ncol, k % ncol]
        _draw_logo(ax, pw, f"{g[:15]} · {c or '-'}", f)
        ax.title.set_fontsize(5.4)
    for k in range(n, nrow * ncol):
        axes[k // ncol, k % ncol].set_axis_off()
    fig.suptitle("Operator logo per TF (grouped by family · cluster; chip = per-TF IC)", fontsize=9,
                 fontweight="bold", x=0.02, ha="left", y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.99))
    save(fig, FIG / "F4c_operator_logo_per_tf")
    print("  F4c_operator_logo_per_tf.svg/.png")


if __name__ == "__main__":
    FIG.mkdir(parents=True, exist_ok=True)
    data = _collect()
    which = sys.argv[1:] or ["family", "cluster", "tf"]
    if "family" in which: fig_family(data)
    if "cluster" in which: fig_cluster(data)
    if "tf" in which: fig_tf(data)
