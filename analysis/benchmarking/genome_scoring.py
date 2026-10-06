"""genome_scoring.py -- the operator and inducer scoring shared by both genome benchmarks.

One definition, used identically for E. coli and Salmonella:

* **Ranked operators.** Every genome-rescan hit of the TF (all pass p <= 1e-4), sorted by PWM
  score. Rank 1 is the primary operator the pipeline recommends.
* **Recovered at top-k.** At least one of the k highest-ranked hits overlaps a truth span within
  `SLOP` bp. Reported for k = 1, 5, 10 and for every hit ("all").
* **Chance expectation.** For each TF and each k, the probability that k sites of the same widths,
  placed uniformly at random in the genome's intergenic DNA (the rescan's own search space), would
  overlap a truth span. Summed over TFs this is the number expected by chance, and the observed
  count is tested against it with an exact Poisson-binomial tail. A lenient metric is only
  informative when it beats this expectation, so the two are always reported together.

All coordinates are the bundle's concatenated "scan" coordinates (`genome/genes.tsv: scan_start,
scan_end`), which is what the rescan hits use on multi-replicon genomes.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

SLOP = 25
TOP_K = (1, 5, 10, "all")
N_NULL = 2000
SEED = 20260924


def load_genes(bundle: Path) -> list[dict]:
    with (bundle / "genome" / "genes.tsv").open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        r["s"], r["e"] = int(r["scan_start"]), int(r["scan_end"])
    return sorted(rows, key=lambda r: r["s"])


def intergenic(genes: list[dict]) -> np.ndarray:
    """[(start, end)] of the gaps between consecutive genes, in scan coordinates."""
    gaps, reach = [], None
    for g in genes:
        if reach is not None and g["s"] > reach + 1:
            gaps.append((reach + 1, g["s"] - 1))
        reach = g["e"] if reach is None else max(reach, g["e"])
    return np.array(gaps, dtype=np.int64)


def ranked_hits(dossier: dict) -> list[dict]:
    hits = (dossier.get("rescan") or {}).get("hits") or []
    return sorted(hits, key=lambda h: (-(h.get("score") or 0), h.get("start") or 0))


def overlaps(a0: int, a1: int, spans: list[tuple[int, int]], slop: int = SLOP) -> bool:
    return any(a0 - slop <= b1 and b0 - slop <= a1 for b0, b1 in spans)


def recovered_at(hits: list[dict], spans: list[tuple[int, int]]) -> dict:
    """{k: bool} -- does any of the top-k hits overlap a truth span?"""
    out = {}
    for k in TOP_K:
        sub = hits if k == "all" else hits[:k]
        out[k] = bool(spans) and any(overlaps(h["start"], h["end"], spans) for h in sub)
    return out


def best_rank(hits: list[dict], spans: list[tuple[int, int]]) -> int | None:
    for i, h in enumerate(hits, 1):
        if overlaps(h["start"], h["end"], spans):
            return i
    return None


def null_probability(k: int, width: int, spans: list[tuple[int, int]], gaps: np.ndarray,
                     rng: np.random.Generator, n: int = N_NULL) -> float:
    """P(at least one of k random intergenic sites of `width` bp overlaps a truth span)."""
    if not spans or k == 0 or len(gaps) == 0:
        return 0.0
    lengths = gaps[:, 1] - gaps[:, 0] + 1
    cum = np.cumsum(lengths)
    draws = rng.integers(0, cum[-1], size=(n, k))
    idx = np.searchsorted(cum, draws, side="right")
    starts = gaps[idx, 0] + (draws - np.concatenate([[0], cum[:-1]])[idx])
    ends = starts + width - 1
    hit = np.zeros((n, k), dtype=bool)
    for b0, b1 in spans:
        hit |= (starts - SLOP <= b1) & (b0 - SLOP <= ends)
    return float(hit.any(axis=1).mean())


def poisson_binomial_tail(ps: list[float], observed: int) -> float:
    """P(X >= observed) for X = sum of independent Bernoulli(p_i)."""
    dist = np.array([1.0])
    for p in ps:
        dist = np.convolve(dist, [1 - p, p])
    return float(dist[observed:].sum()) if observed < len(dist) else 0.0


def operator_block(rows: list[dict]) -> dict:
    """Summarise rows carrying `rec_<k>` (bool) and `null_<k>` (prob) into observed vs chance."""
    out = {"n": len(rows)}
    for k in TOP_K:
        obs = sum(bool(r[f"rec_{k}"]) for r in rows)
        ps = [r[f"null_{k}"] for r in rows]
        out[str(k)] = {"observed": obs, "expected_by_chance": round(sum(ps), 2),
                       "p_vs_chance": float(f"{poisson_binomial_tail(ps, obs):.3g}") if rows else None}
    out["median_hits_per_tf"] = int(np.median([r["n_hits"] for r in rows])) if rows else 0
    return out


def score_operator(bundle: Path, dossier: dict, spans: list[tuple[int, int]],
                   gaps: np.ndarray, rng: np.random.Generator) -> dict:
    hits = ranked_hits(dossier)
    rec = recovered_at(hits, spans)
    width = int(np.median([h["end"] - h["start"] + 1 for h in hits])) if hits else 20
    row = {"n_hits": len(hits), "n_truth_spans": len(spans), "best_rank": best_rank(hits, spans),
           "primary_locality": ("proximal" if hits and hits[0].get("generator") == "rescan_tier1"
                                else ("distal" if hits else ""))}
    for k in TOP_K:
        kk = len(hits) if k == "all" else min(k, len(hits))
        row[f"rec_{k}"] = rec[k]
        row[f"null_{k}"] = null_probability(kk, width, spans, gaps, rng)
    return row


def autoregulated(dossier: dict, k: int = 5) -> bool:
    """A top-k operator lies in the regulator's own flanking intergenic region (rescan tier 1)."""
    return any(h.get("generator") == "rescan_tier1" for h in ranked_hits(dossier)[:k])


def load_dossier(bundle: Path) -> dict:
    return json.loads((bundle / "dossier.json").read_text(encoding="utf-8"))
