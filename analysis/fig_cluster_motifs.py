#!/usr/bin/env python
"""F4b (rebuilt) -- do the regulators of an SSN cluster share an operator motif?

WHY THIS REPLACES THE OLD F4b.  The previous panel ranked SSN clusters by the total information content
(IC) of a pooled logo built by centring each member PWM on its OWN maximum-IC column and averaging.  Both
steps inflate IC, and averaging does so in inverse proportion to the number of members, so the published
ranking was very nearly a plot of 1/n: Spearman rho(n_regulators, IC) = -0.78 (-0.84 excluding the MerR_9
outlier).  Every 2-member cluster scored 10-38 bits; the 8-member cluster scored 3.  The showcase
"MerR_9/ZntR = 38 bits" is two orthologous ZntR proteins from two Vibrio genomes -- near-identical copies
of one regulator, not a crisp operator.

The rebuilt figure fixes three things:
  (1) ALIGNMENT.  Cluster members are aligned to each other by maximising mean per-column correlation over
      all offsets and both orientations (a Tomtom-style motif alignment), not each on its own peak column.
  (2) A STATISTIC THAT IS COMPARABLE ACROSS CLUSTER SIZES.  The primary readout is the mean pairwise motif
      similarity WITHIN a cluster, compared against a null of between-cluster pairs.  Unlike pooled IC,
      this does not shrink as n grows, and it answers the actual biological question.
  (3) NO IC RANKING.  Clusters are ordered by family then name.  IC is still shown, with n beside it and an
      explicit note that IC is not comparable across n; the IC-vs-n relationship is plotted so the reader
      can see the artefact rather than be misled by it.

Run with the `data` env (needs logomaker + scipy):
    $TFOP_PY fig_cluster_motifs.py
"""
from __future__ import annotations
import csv, json, sys
from collections import defaultdict
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig_style import save, require, focus_family, FAM, BASE, INK, NEUTRAL_M, NEUTRAL_L
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import spearmanr, mannwhitneyu


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
OUT = RP.RESULTS
FAMS = ("ArsR", "MerR", "Fur")
MIN_OVERLAP = 8          # columns required to score an alignment offset
LN2 = np.log(2)


#: REMOVED: the family must come from the bundle, not the run name. `fig_style.focus_family`
#: reads `dossier["family"]`; the run name encodes what DISCOVERY guessed from a partial
#: profile match, which the full-sequence classifier can and does contradict.


def _norm(p):
    p = np.asarray(p, float)
    if p.shape[0] != 4:
        p = p.T
    return p / p.sum(0, keepdims=True).clip(1e-9)


def _rc(p):
    """Reverse complement of a 4xW PWM (rows A,C,G,T)."""
    return p[::-1, ::-1]


def _col_ic(p):
    return 2.0 + (p * np.log2(p.clip(1e-9))).sum(0)


def _corrected_ic(p, n_sites):
    """Total IC with the Schneider et al. (1986) small-sample correction e(n) = (s-1)/(2 ln2 n), s=4."""
    e = 3.0 / (2 * LN2 * max(1, n_sites))
    return float(np.clip(_col_ic(p) - e, 0, None).sum())


def similarity(P, Q):
    """Tomtom-style: max over offsets and orientations of the mean per-column Pearson r."""
    best = -1.0
    for Qo in (Q, _rc(Q)):
        w1, w2 = P.shape[1], Qo.shape[1]
        for off in range(-w2 + MIN_OVERLAP, w1 - MIN_OVERLAP + 1):
            a0, b0 = max(0, off), max(0, -off)
            ov = min(w1 - a0, w2 - b0)
            if ov < MIN_OVERLAP:
                continue
            A, B = P[:, a0:a0 + ov], Qo[:, b0:b0 + ov]
            A = A - A.mean(0, keepdims=True)
            B = B - B.mean(0, keepdims=True)
            den = np.sqrt((A ** 2).sum(0) * (B ** 2).sum(0)).clip(1e-9)
            best = max(best, float(((A * B).sum(0) / den).mean()))
    return best


def align_to_reference(pwms, width=None):
    """Align every PWM to the member that is on average most similar to the others, then average."""
    n = len(pwms)
    S = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            S[i, j] = S[j, i] = similarity(pwms[i], pwms[j])
    ref = int(np.argmax(S.sum(1)))
    R = pwms[ref]
    W = width or R.shape[1]
    acc = np.zeros((4, W)); cnt = np.zeros(W)
    for k, P in enumerate(pwms):
        bestscore, bestoff, bestP = -1.0, 0, P
        for Po in (P, _rc(P)):
            for off in range(-Po.shape[1] + MIN_OVERLAP, W - MIN_OVERLAP + 1):
                a0, b0 = max(0, off), max(0, -off)
                ov = min(W - a0, Po.shape[1] - b0)
                if ov < MIN_OVERLAP:
                    continue
                A, B = R[:, a0:a0 + ov], Po[:, b0:b0 + ov]
                A = A - A.mean(0, keepdims=True); B = B - B.mean(0, keepdims=True)
                den = np.sqrt((A ** 2).sum(0) * (B ** 2).sum(0)).clip(1e-9)
                s = float(((A * B).sum(0) / den).mean())
                if s > bestscore:
                    bestscore, bestoff, bestP = s, off, Po
        a0, b0 = max(0, bestoff), max(0, -bestoff)
        ov = min(W - a0, bestP.shape[1] - b0)
        acc[:, a0:a0 + ov] += bestP[:, b0:b0 + ov]; cnt[a0:a0 + ov] += 1
    cnt[cnt == 0] = 1
    return _norm(acc / cnt), S


def collect():
    man = {r["run_name"]: r for r in csv.DictReader(_manifest_path().open(encoding="utf-8"))}
    out = {}
    for name, m in man.items():
        if not (JOBS / name / "dossier.json").exists():
            continue
        j = json.loads((JOBS / name / "dossier.json").read_text(encoding="utf-8"))
        fam = focus_family(j.get("family"))
        if fam is None:
            continue
        gl = j.get("genomic_logo") or {}
        if gl.get("pwm"):
            out[name] = dict(fam=fam, clu=j.get("ssn_cluster"),
                             pwm=_norm(gl["pwm"]), n_sites=int(gl.get("n_kept") or 0))
    return out


def draw_logo(ax, pwm, n_sites):
    import logomaker, pandas as pd
    ic = _col_ic(pwm)
    df = pd.DataFrame((pwm * ic).T, columns=list("ACGT"))
    logomaker.Logo(df, ax=ax, color_scheme={b: BASE[b] for b in "ACGT"}, show_spines=False)
    ax.set_ylim(0, 2); ax.set_xticks([]); ax.set_yticks([0, 2])
    ax.tick_params(axis="y", labelsize=4.6, length=1.5)
    ax.set_ylabel("bits", fontsize=4.8, labelpad=1)


def _skip_not_applicable(n_found: int, need: int) -> None:
    """Emit an explicit 'not applicable' F4b when the run has too few multi-regulator clusters.

    The within- vs between-cluster comparison needs several clusters each carrying two or more
    regulators with a genomic logo; a single-organism validation run typically has too few. Rather
    than hard-fail and block the report, draw a labelled placeholder and record the reason, so
    figure QA still finds a well-formed F4b in all four formats and the report states plainly why
    the panel is empty. The four-genome production run has >= `need` clusters and draws the real
    figure by the path below.
    """
    fig = plt.figure(figsize=(7.2, 4.0))
    fig.suptitle("Do the regulators of an SSN cluster share an operator motif?", fontsize=9,
                 fontweight="bold", x=0.02, ha="left", y=0.97)
    ax = fig.add_axes([0.06, 0.10, 0.88, 0.72]); ax.axis("off")
    ax.text(0.5, 0.5,
            "Not applicable to this run.\n\n"
            f"This panel compares operator-motif similarity within an SSN cluster against a\n"
            f"between-cluster null, which needs at least {need} clusters each carrying two or more\n"
            f"regulators with a genomic logo. This run has {n_found}. A single-organism validation\n"
            "run typically has too few; the four-genome production run does not.",
            ha="center", va="center", fontsize=7.5, color=INK, linespacing=1.6)
    save(fig, FIG / "F4b_cluster_motif_sharing")

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "F4b_cluster_motif_stats.json").write_text(json.dumps(
        {"skipped": True, "reason": "insufficient multi-regulator clusters",
         "clusters_found": n_found, "clusters_required": need}, indent=2), encoding="utf-8")
    with (OUT / "cluster_motif_similarity.csv").open("w", newline="", encoding="utf-8") as fh:
        csv.writer(fh).writerow(
            ["cluster", "family", "n_regulators", "n_sites", "mean_within_similarity", "total_IC_bits"])
    print(f"  F4b SKIPPED (not applicable): {n_found} multi-regulator cluster(s), need >={need}. "
          "Drew a labelled placeholder so the report can build.")


def main():
    data = collect()
    clu = defaultdict(list)
    for v in data.values():
        if v["clu"]:
            clu[v["clu"]].append(v)
    items = {c: vs for c, vs in clu.items() if len(vs) >= 2}
    fam_of = {c: vs[0]["fam"] for c, vs in items.items()}
    order = sorted(items, key=lambda c: (FAMS.index(fam_of[c]), c))

    # This panel compares motifs WITHIN a cluster against a between-cluster null, so it needs several
    # clusters carrying two or more regulators each. Below that the right-hand column's row slice
    # collapses and matplotlib raises "GridSpec slice would result in no space allocated for
    # subplot" -- an error that says nothing about the real cause, which is simply too little data.
    # Partway through a batch that is the NORMAL state, so say so and skip rather than crash.
    MIN_CLUSTERS = 6
    require(data, "genomic logos in this run's bundles",
            "the batch stage (run_phase1.py)",
            "No bundle carried a genomic_logo -- check TFOP_RUN_TAG and that jobs/ holds this "
            "run's bundles rather than an archived set.")
    # A cross-cluster comparison genuinely cannot be drawn from fewer than MIN_CLUSTERS multi-
    # regulator clusters (a single-organism validation run usually has one). Draw an explicit
    # 'not applicable' placeholder and exit 0 rather than hard-failing and blocking the report.
    # The four-genome production run has >= MIN_CLUSTERS and falls through to the real figure.
    if len(order) < MIN_CLUSTERS:
        _skip_not_applicable(len(order), MIN_CLUSTERS)
        return

    # ---- per-cluster aligned consensus + within-cluster similarity
    rows, within, pooled = {}, {}, {}
    for c in order:
        pwms = [v["pwm"] for v in items[c]]
        avg, S = align_to_reference(pwms)
        iu = np.triu_indices(len(pwms), 1)
        within[c] = float(S[iu].mean())
        pooled[c] = avg
        rows[c] = dict(family=fam_of[c], n_regulators=len(pwms),
                       n_sites=int(sum(v["n_sites"] for v in items[c])),
                       mean_within_similarity=within[c],
                       total_IC_bits=_corrected_ic(avg, sum(v["n_sites"] for v in items[c])))

    # ---- between-cluster null (same family where possible, else any)
    between = []
    keys = list(items)
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            for va in items[a][:3]:
                for vb in items[b][:3]:
                    between.append(similarity(va["pwm"], vb["pwm"]))
    between = np.asarray(between)
    w = np.asarray([within[c] for c in order])
    u, p_sim = mannwhitneyu(w, between, alternative="greater")

    # ---- IC vs n artefact
    ns = np.asarray([rows[c]["n_regulators"] for c in order], float)
    ics = np.asarray([rows[c]["total_IC_bits"] for c in order], float)
    rho, p_rho = spearmanr(ns, ics)

    # =============================================================== figure
    nrow = len(order)
    fig = plt.figure(figsize=(7.2, 0.52 * nrow + 2.5))
    gs = gridspec.GridSpec(nrow + 1, 2, width_ratios=[1.55, 1.0], hspace=0.85, wspace=0.32,
                           left=0.09, right=0.98, top=0.93, bottom=0.07)
    for k, c in enumerate(order):
        ax = fig.add_subplot(gs[k, 0])
        draw_logo(ax, pooled[c], rows[c]["n_sites"])
        ax.set_title(f"{c}   n={rows[c]['n_regulators']} regulators · within-cluster r = {within[c]:.2f} · "
                     f"IC {rows[c]['total_IC_bits']:.0f} bits",
                     fontsize=5.6, color=FAM[fam_of[c]], loc="left", pad=1.5)

    # right column: two diagnostics
    # leave a blank band of rows between the two diagnostics so the 2-line titles never collide
    gap = 3
    hi = max(2, (nrow - gap) // 2)
    axA = fig.add_subplot(gs[0:hi, 1])
    axA.hist(between, bins=24, color=NEUTRAL_L, edgecolor="white", linewidth=0.3, density=True,
             label=f"between-cluster pairs (n={between.size})")
    for c in order:
        axA.axvline(within[c], color=FAM[fam_of[c]], lw=0.9, alpha=0.85)
    axA.set_xlabel("motif similarity (mean column r)", fontsize=6.2)
    axA.set_ylabel("density", fontsize=6.2)
    axA.set_title(f"Within-cluster similarity vs between-cluster null\n"
                  f"median within {np.median(w):.2f} vs between {np.median(between):.2f}, P = {p_sim:.1e}",
                  fontsize=6.4, color=INK, loc="left", pad=3)
    axA.tick_params(labelsize=5.6)
    axA.legend(fontsize=5.0, frameon=False, loc="upper left")

    axB = fig.add_subplot(gs[hi + gap:nrow, 1])
    for c in order:
        axB.scatter(rows[c]["n_regulators"], rows[c]["total_IC_bits"], s=16, color=FAM[fam_of[c]],
                    alpha=0.9, edgecolor="white", linewidth=0.4, zorder=3)
    top = max(rows[c]["total_IC_bits"] for c in order)
    out_c = max(order, key=lambda c: rows[c]["total_IC_bits"])
    # ha="right" so the label grows LEFTWARDS, into the panel. Left-aligned from x=3.4 it ran off
    # the rightmost panel and off the figure, and since `save()` exports with bbox_inches="tight"
    # that overflow became figure width -- this annotation, not the title, was what held F4b at
    # 185 mm against the 183 mm limit.
    axB.annotate(f"{out_c}\n(2 orthologues of one regulator)",
                 xy=(rows[out_c]["n_regulators"], top), xytext=(3.4, top - 1.5),
                 fontsize=5.0, color=NEUTRAL_M, ha="right",
                 arrowprops=dict(arrowstyle="-", lw=0.5, color=NEUTRAL_M))
    axB.set_xlabel("regulators contributing to the cluster logo", fontsize=6.2)
    axB.set_ylabel("total IC of pooled logo (bits)", fontsize=6.2)
    # Wrapped onto four short lines. `loc="left"` anchors the title at the axes' left edge, so a
    # long line runs RIGHT past the panel and, on the rightmost panel, past the figure -- which is
    # what pushed this figure over the 183 mm publication limit.
    # Three lines was not enough: the statistic shared a line with the claim, making that line the
    # longest in the figure and leaving the export 2.2 mm over on 2026-09-04. The claim and the
    # statistic are now on separate lines, so the longest line is prose of a fixed length rather
    # than one that grows with the formatted numbers.
    axB.set_title("Why IC must NOT be ranked\n"
                  f"Spearman ρ = {rho:+.2f} (P = {p_rho:.3f})\n"
                  "pooled IC falls as more regulators\n"
                  "are averaged — it measures n, not conservation",
                  fontsize=6.4, color=INK, loc="left", pad=3)
    axB.tick_params(labelsize=5.6)

    fig.suptitle("Do the regulators of an SSN cluster share an operator motif?", fontsize=9,
                 fontweight="bold", x=0.02, ha="left", y=0.985)
    save(fig, FIG / "F4b_cluster_motif_sharing")

    OUT.mkdir(parents=True, exist_ok=True)
    summary = {"clusters": rows,
               "between_cluster_null": {"n_pairs": int(between.size), "median": float(np.median(between))},
               "within_vs_between": {"median_within": float(np.median(w)), "U": float(u), "P": float(p_sim)},
               "IC_vs_n": {"spearman_rho": float(rho), "P": float(p_rho)}}
    (OUT / "F4b_cluster_motif_stats.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    with (OUT / "cluster_motif_similarity.csv").open("w", newline="", encoding="utf-8") as fh:
        wtr = csv.writer(fh)
        wtr.writerow(["cluster", "family", "n_regulators", "n_sites", "mean_within_similarity", "total_IC_bits"])
        for c in order:
            r = rows[c]
            wtr.writerow([c, r["family"], r["n_regulators"], r["n_sites"],
                          f"{r['mean_within_similarity']:.3f}", f"{r['total_IC_bits']:.1f}"])
    print("  F4b_cluster_motif_sharing.svg/.png")
    # NB print these separately -- both dicts carry a "P" key and merging them silently hides one
    print("within_vs_between:", json.dumps(summary["within_vs_between"]))
    print("between_null:     ", json.dumps(summary["between_cluster_null"]))
    print("IC_vs_n:          ", json.dumps(summary["IC_vs_n"]))
    for c in order:
        print(f"  {c:14s} n={rows[c]['n_regulators']:2d}  r={within[c]:+.2f}  IC={rows[c]['total_IC_bits']:5.1f}")


if __name__ == "__main__":
    main()
