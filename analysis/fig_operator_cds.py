#!/usr/bin/env python
"""F3 (rebuilt) -- operator position relative to the CDS start codon, against a genomic null.

WHY THIS REPLACES THE OLD F3.  The previous figure plotted `op_to_tss_bp`, the distance from the operator
to a sigma70 TSS.  That TSS was predicted inside a window CENTRED ON THE OPERATOR and then selected for
maximal overlap with that same operator (`promoter.classify_mode`), so the distance was largely determined
by the selection rule, not by biology -- the whole distribution spanned only -38..+27 bp.

The rebuilt figure uses an operator-INDEPENDENT anchor: the 5' end of the coding sequence (start codon) of
the assigned gene, taken straight from the paper's own GFF3.  Distance is `operator_to_gene_bp`, i.e.
(gene 5' end - operator centre), strand-aware, as assigned by the pipeline's `_regulated_gene` rule
(nearest gene 5' end downstream of the operator, within 400 bp).

Two things the old figure lacked and a reviewer will demand:
  (1) A NULL.  Random genomic positions are pushed through the IDENTICAL assignment rule, so the observed
      distribution is compared against "where would any position fall relative to the nearest CDS start?"
      given the gene density of these genomes.  Two nulls are computed: unrestricted (any position) and
      the conservative intergenic-restricted one (positions outside annotated CDS).
  (2) NO PSEUDOREPLICATION in the inferential claim.  The old n=1500 was 25 regulators x their top 60 hits.
      Here the test statistic is the PER-REGULATOR median (n = 25/24/8 regulators), compared against
      per-regulator medians drawn from the null with matched hit counts.

Run with the `data` env (needs scipy):
    $TFOP_PY fig_operator_cds.py
"""
from __future__ import annotations
import csv, json, os, sys
from collections import defaultdict
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig_style import save, require, focus_family, FAM, INK, NEUTRAL_M, NEUTRAL_L
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu


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
MIRROR = Path(os.environ.get("PREDICTOR_MIRROR", Path.home() / ".predictor" / "genome_mirror"))
FAMS = ("ArsR", "MerR", "Fur")
MAX_DIST = 400          # the pipeline's _regulated_gene window
N_NULL_DRAWS = 400_000  # random positions per genome
SEED = 20260720


#: REMOVED: the family must come from the bundle, not the run name. `fig_style.focus_family`
#: reads `dossier["family"]`; the run name encodes what DISCOVERY guessed from a partial
#: profile match, which the full-sequence classifier can and does contradict.


# --------------------------------------------------------------------------- genome annotation
def load_cds(acc):
    """{contig: (starts, ends, strands)} for CDS features, 1-based inclusive, from the paper's own GFF."""
    gff = MIRROR / f"{acc}.gff"
    if not gff.exists():
        raise SystemExit(f"missing GFF for {acc}: {gff}")
    per = defaultdict(list)
    with gff.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8 or f[2] != "CDS":
                continue
            per[f[0]].append((int(f[3]), int(f[4]), f[6]))
    out = {}
    for ctg, rows in per.items():
        rows.sort()
        out[ctg] = (np.array([r[0] for r in rows]), np.array([r[1] for r in rows]),
                    np.array([r[2] for r in rows]))
    return out


def contig_lengths(acc):
    """{contig: length} from the GFF ##sequence-region headers."""
    out = {}
    with (MIRROR / f"{acc}.gff").open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("##sequence-region"):
                p = line.split()
                if len(p) >= 4:
                    out[p[1]] = int(p[3])
            elif not line.startswith("#"):
                break
    return out


def assign_distance(pos, starts, ends, strands):
    """Vectorised `_regulated_gene`: nearest gene 5' end downstream of `pos` within MAX_DIST.
    For '+' genes the 5' end is `start` and d = start - pos; for '-' genes it is `end` and d = pos - end."""
    plus = strands == "+"
    p5 = np.sort(starts[plus])                 # '+' 5' ends, ascending
    m5 = np.sort(ends[~plus])                  # '-' 5' ends, ascending
    pos = np.asarray(pos)
    best = np.full(pos.shape, MAX_DIST + 1, dtype=np.int64)
    if p5.size:                                # smallest start >= pos
        i = np.searchsorted(p5, pos, side="left")
        ok = i < p5.size
        d = np.where(ok, p5[np.clip(i, 0, p5.size - 1)] - pos, MAX_DIST + 1)
        best = np.minimum(best, np.where(ok & (d >= 0), d, MAX_DIST + 1))
    if m5.size:                                # largest end <= pos
        j = np.searchsorted(m5, pos, side="right") - 1
        ok = j >= 0
        d = np.where(ok, pos - m5[np.clip(j, 0, m5.size - 1)], MAX_DIST + 1)
        best = np.minimum(best, np.where(ok & (d >= 0), d, MAX_DIST + 1))
    return best[best <= MAX_DIST]


def cds_mask(ctg_len, starts, ends):
    m = np.zeros(ctg_len + 2, dtype=bool)
    for s, e in zip(starts, ends):
        m[s:e + 1] = True
    return m


def build_null(accs, rng):
    """Distances for random positions under the identical assignment rule.
    Returns (unrestricted, intergenic_restricted) arrays pooled over the given assemblies."""
    unres, inter = [], []
    for acc in accs:
        cds = load_cds(acc)
        lens = contig_lengths(acc)
        for ctg, (st, en, sd) in cds.items():
            L = lens.get(ctg) or int(en.max()) + 1000
            n = max(1000, int(N_NULL_DRAWS * L / max(1, sum(lens.values() or [L]))))
            pos = rng.integers(1, L + 1, size=n)
            unres.append(assign_distance(pos, st, en, sd))
            mask = cds_mask(L, st, en)
            ig = pos[~mask[pos]]
            if ig.size:
                inter.append(assign_distance(ig, st, en, sd))
    return (np.concatenate(unres) if unres else np.array([]),
            np.concatenate(inter) if inter else np.array([]))


# --------------------------------------------------------------------------- observed
def collect_observed():
    man = {r["run_name"]: r for r in csv.DictReader(_manifest_path().open(encoding="utf-8"))}
    obs_hits, obs_med, accs = defaultdict(list), defaultdict(list), set()
    for name, m in man.items():
        if not (JOBS / name / "dossier.json").exists():
            continue
        j = json.loads((JOBS / name / "dossier.json").read_text(encoding="utf-8"))
        if focus_family(j.get("family")) is None:
            continue
        accs.add(m["organism_acc"])
        d = [h["operator_to_gene_bp"] for h in (j.get("per_hit_regulation") or [])
             if h.get("operator_to_gene_bp") is not None]
        if not d:
            continue
        fam = focus_family(j.get("family"))
        obs_hits[fam] += d
        obs_med[fam].append(float(np.median(d)))
    return obs_hits, obs_med, sorted(accs)


def null_per_regulator(null_pool, counts, rng):
    """Per-regulator medians drawn from the null with matched hit counts (kills pseudoreplication)."""
    return [float(np.median(rng.choice(null_pool, size=max(1, c), replace=True))) for c in counts]


def main():
    rng = np.random.default_rng(SEED)
    obs_hits, obs_med, accs = collect_observed()
    # Without this, a run whose bundles are absent or archived drew three empty histograms, wrote
    # the statistics file, printed a summary and exited 0 -- a complete-looking F3 about nothing.
    require(any(len(v) for v in obs_hits.values()),
            "operator-to-CDS distances from this run's bundles",
            "the batch stage (run_phase1.py), then aggregate",
            "Check TFOP_RUN_TAG points at a run whose jobs/ still holds its bundles.")
    require(accs, "genome assemblies referenced by the manifest", "the manifest and the batch stage")
    print(f"assemblies: {accs}")
    null_all, null_ig = build_null(accs, rng)
    require(null_ig.size or null_all.size, "the genomic null distribution",
            "this script, from the run's GFF3 annotations",
            "The null could not be built, so there is nothing to compare the operators against.")
    print(f"null draws: unrestricted {null_all.size:,}  intergenic {null_ig.size:,}")

    stats = {}
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.0), height_ratios=[2.2, 1.0])
    bins = np.arange(0, MAX_DIST + 1, 10)
    for k, fam in enumerate(FAMS):
        ax = axes[0, k]
        d = np.asarray(obs_hits[fam], float)
        # CDS body: the gene starts at 0 and runs to the right (operators sit upstream => negative axis)
        ax.axvspan(0, 60, color=NEUTRAL_L, alpha=0.75, zorder=0, lw=0)
        # x in DATA coords, y in AXES fraction -- never place an artist at a data y before the scale is set
        ax.text(30, 1.01, "CDS", ha="center", va="bottom", fontsize=5.4, color=NEUTRAL_M,
                transform=ax.get_xaxis_transform())
        ax.hist(-d, bins=-bins[::-1], color=FAM[fam], alpha=0.9, edgecolor="white",
                linewidth=0.3, zorder=2, density=True, label="operator hits")
        for pool, ls, lab in ((null_ig, "-", "null (intergenic)"), (null_all, ":", "null (any position)")):
            if pool.size:
                h, e = np.histogram(pool, bins=bins, density=True)
                ax.plot(-(e[:-1] + 5), h, ls, color=INK, lw=0.9, zorder=3, label=lab)
        om, nm = float(np.median(d)), float(np.median(null_ig))
        med = np.asarray(obs_med[fam])
        nullmed = np.asarray(null_per_regulator(null_ig, [len(obs_hits[fam]) // max(1, len(med))] * len(med), rng))
        u, p = mannwhitneyu(med, nullmed, alternative="less")
        stats[fam] = {"n_regulators": len(med), "n_hits": int(d.size),
                      "median_obs_bp": om, "median_null_intergenic_bp": nm,
                      "median_of_per_regulator_medians_bp": float(np.median(med)),
                      "mannwhitney_U": float(u), "p_one_sided_closer_than_null": float(p)}
        ax.set_title(f"{fam}  (n={len(med)} regulators)", color=FAM[fam], fontweight="bold", fontsize=8)
        ax.text(0.03, 0.96, f"median {om:.0f} bp\nnull {nm:.0f} bp\nP = {p:.1e}",
                transform=ax.transAxes, va="top", fontsize=5.6)
        ax.set_xlim(-MAX_DIST, 60)
        ax.set_xlabel("operator centre → CDS start (bp)", fontsize=6.3)
        ax.tick_params(labelsize=5.8)
        if k == 0:
            ax.set_ylabel("density", fontsize=6.5)
            ax.legend(fontsize=4.9, frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.62))

        # ---- lower strip: ONE POINT PER REGULATOR (no pseudoreplication)
        axs = axes[1, k]
        y = rng.normal(0, 0.055, size=len(med))
        axs.scatter(-med, y, s=9, color=FAM[fam], alpha=0.85, edgecolor="white", linewidth=0.3, zorder=3)
        axs.axvline(-np.median(nullmed), color=INK, ls="-", lw=0.9, zorder=2)
        axs.text(-np.median(nullmed), 0.30, " null", fontsize=5.2, color=INK, va="center")
        axs.set_xlim(-MAX_DIST, 60); axs.set_ylim(-0.32, 0.32)
        axs.set_yticks([]); axs.tick_params(labelsize=5.8)
        axs.set_xlabel("per-regulator median (bp)", fontsize=6.0)
        for sp in ("left", "right", "top"):
            axs.spines[sp].set_visible(False)

    fig.suptitle("Operator position relative to the CDS start codon, against a genomic null",
                 fontsize=9, fontweight="bold", x=0.02, ha="left", y=1.02)
    # WRAPPED, deliberately. As one line this caption was ~250 characters -- about 10 in at 5.2 pt --
    # and `bbox_inches="tight"` widened the saved page to 204.8 mm, over the 183 mm limit. The text
    # is unchanged; only the line breaks are new.
    fig.text(0.02, 0.008,
             "Anchor is the annotated CDS 5′ end (paper's own GFF3), independent of the operator.\n"
             "Null = random positions through the identical gene-assignment rule. Lower strips show\n"
             "one point per regulator; P from Mann–Whitney U on per-regulator medians.",
             fontsize=5.2, color=NEUTRAL_M, ha="left", linespacing=1.5)
    fig.tight_layout(rect=(0, 0.05, 1, 0.97))
    save(fig, FIG / "F3_operator_cds_distance")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "F3_operator_cds_stats.json").write_text(json.dumps(stats, indent=2), encoding="utf-8")
    print("  F3_operator_cds_distance.svg/.png")
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
