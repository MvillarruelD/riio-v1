"""
pwm_scan.py -- FIMO-equivalent PWM scanner, pure Python (numpy/scipy).

Replaces MEME-suite **FIMO** so the pipeline needs no FIMO binary (no WSL / Docker).
Given a motif (counts or probabilities) and DNA sequence(s) it:

  1. builds a log-odds PSSM against a background model,
  2. scans BOTH strands,
  3. assigns an EXACT p-value to every score via an integer-lattice convolution of the
     per-column score distributions -- the same approach FIMO and MOODS use,
  4. converts reported p-values to Benjamini-Hochberg q-values,
  5. returns hits ranked by p-value.

This is the engine behind `signals/motif_rescan.py` (genome rescan, Phase 1.6A) and the
family-consensus motif floor (Phase 1). Coordinates are 0-based half-open on the + strand of
the input sequence; `signals/coords.py` lifts them to the canonical genome frame.

Run `python pwm_scan.py` for a self-test (plants a site and recovers it).
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

ALPH = "ACGT"
IDX = {c: i for i, c in enumerate(ALPH)}
_COMP = np.array([3, 2, 1, 0])          # A<->T, C<->G  (row permutation for reverse-complement)


@dataclass
class Hit:
    seq_id: str
    start: int          # 0-based, + strand
    end: int            # exclusive
    strand: str         # '+' or '-'
    score: float        # log-odds (bits when log_base=2)
    pvalue: float
    qvalue: float = float("nan")
    matched: str = ""   # + strand substring under the motif


# ----------------------------------------------------------------------------- matrices
def pwm_from_counts(counts, pseudocount: float = 0.5) -> np.ndarray:
    """counts (4xW or Wx4, A,C,G,T order) -> probability matrix (4xW)."""
    counts = np.asarray(counts, dtype=float)
    if counts.shape[0] != 4 and counts.shape[1] == 4:
        counts = counts.T
    counts = counts + pseudocount
    return counts / counts.sum(axis=0, keepdims=True)


def background_from_seqs(seqs) -> np.ndarray:
    codes = np.concatenate([_encode(s) for s in seqs])
    codes = codes[codes >= 0]
    bg = np.array([(codes == a).mean() for a in range(4)], dtype=float)
    bg = np.where(bg <= 0, 1e-9, bg)
    return bg / bg.sum()


def log_odds(prob, background=None, log_base: float = 2.0) -> np.ndarray:
    prob = np.asarray(prob, dtype=float)
    if prob.shape[0] != 4:
        prob = prob.T
    bg = np.full(4, 0.25) if background is None else np.asarray(background, float)
    bg = bg / bg.sum()
    return np.log(prob / bg[:, None]) / np.log(log_base)


def _rc_pssm(pssm: np.ndarray) -> np.ndarray:
    return pssm[_COMP][:, ::-1]


# ----------------------------------------------------------------------------- encoding
def _encode(seq: str) -> np.ndarray:
    """DNA string -> int8 array; non-ACGT -> -1."""
    b = np.frombuffer(seq.upper().encode("ascii", "replace"), dtype=np.uint8)
    out = np.full(b.shape[0], -1, dtype=np.int8)
    for c, i in IDX.items():
        out[b == ord(c)] = i
    return out


def _scan_scores(pssm: np.ndarray, codes: np.ndarray) -> np.ndarray:
    """Score every length-W window; windows containing N -> -inf."""
    _, W = pssm.shape
    n = codes.shape[0] - W + 1
    if n <= 0:
        return np.empty(0)
    scores = np.zeros(n)
    valid = np.ones(n, dtype=bool)
    for j in range(W):
        c = codes[j:j + n]
        bad = c < 0
        scores += np.where(bad, 0.0, pssm[np.where(bad, 0, c), j])
        valid &= ~bad
    scores[~valid] = -np.inf
    return scores


# --------------------------------------------------------------- exact null distribution
def _null_distribution(pssm: np.ndarray, background: np.ndarray, scale: float):
    """Integer-lattice convolution of per-column score PMFs under `background`
    (the FIMO / MOODS method). Returns (offset, pmf): raw_score = (offset + i) / scale."""
    cols = np.rint(pssm * scale).astype(int)        # (4, W)
    dist = None
    offset = 0
    for j in range(cols.shape[1]):
        col = cols[:, j]
        cmin = int(col.min())
        cd = np.zeros(int(col.max()) - cmin + 1)
        for a in range(4):
            cd[col[a] - cmin] += background[a]
        if dist is None:
            dist, offset = cd, cmin
        else:
            dist = np.convolve(dist, cd)
            offset += cmin
    dist = dist / dist.sum()
    return offset, dist


def _pvalue_fn(offset: int, dist: np.ndarray, scale: float):
    """Survival function P(score >= s) over the null lattice."""
    sf = np.cumsum(dist[::-1])[::-1]                  # sf[k] = sum_{>=k}
    n = sf.shape[0]

    def pv(raw_score: float) -> float:
        k = int(np.ceil(raw_score * scale)) - offset
        if k < 0:
            return 1.0
        if k >= n:
            return float(dist[-1]) if n else 0.0
        return float(sf[k])
    return pv


def bh_qvalues(pvals) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    n = p.shape[0]
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0.0, 1.0)
    return out


# ----------------------------------------------------------------------------- public API
def scan(sequences, counts=None, prob=None, background=None,
         pvalue_thresh: float = 1e-3, both_strands: bool = True,
         pseudocount: float = 0.5, scale: float = 50.0, log_base: float = 2.0):
    """Scan `sequences` (str | {id: str} | [str]) with a motif. Returns ranked list[Hit].

    Provide either `counts` (raw alignment counts) or `prob` (probabilities). Background
    defaults to the composition of the input sequences (FIMO's default behaviour).
    """
    if isinstance(sequences, str):
        sequences = {"seq": sequences}
    elif not isinstance(sequences, dict):
        sequences = {f"seq{i}": s for i, s in enumerate(sequences)}

    if prob is None:
        if counts is None:
            raise ValueError("provide counts= or prob=")
        prob = pwm_from_counts(counts, pseudocount)
    if background is None:
        background = background_from_seqs(sequences.values())
    background = np.asarray(background, float)
    background = background / background.sum()

    pssm = log_odds(prob, background, log_base)
    offset, dist = _null_distribution(pssm, background, scale)
    pv = _pvalue_fn(offset, dist, scale)
    W = pssm.shape[1]

    strands = [("+", pssm)] + ([("-", _rc_pssm(pssm))] if both_strands else [])
    hits: list[Hit] = []
    for sid, seq in sequences.items():
        codes = _encode(seq)
        for strand, mat in strands:
            for i, s in enumerate(_scan_scores(mat, codes)):
                if not np.isfinite(s):
                    continue
                p = pv(s)
                if p <= pvalue_thresh:
                    hits.append(Hit(sid, i, i + W, strand, float(s), p,
                                    matched=seq[i:i + W]))
    if hits:
        for h, q in zip(hits, bh_qvalues([h.pvalue for h in hits])):
            h.qvalue = q
    hits.sort(key=lambda h: (h.pvalue, -h.score))
    return hits


# ----------------------------------------------------------------------------- self-test
def _demo() -> None:
    import random
    sites = ["TTGACA", "TTGACA", "TTGTCA", "TTGACA", "TAGACA"]   # toy -35-like motif
    W = len(sites[0])
    counts = np.zeros((4, W))
    for s in sites:
        for j, ch in enumerate(s):
            counts[IDX[ch], j] += 1

    random.seed(0)
    bg = "".join(random.choice("ACGT") for _ in range(500))
    planted = bg[:200] + "TTGACA" + bg[200:]                     # site at 0-based 200
    hits = scan({"chr": planted}, counts=counts, pvalue_thresh=1e-3)

    print(f"reported {len(hits)} hit(s); top 5:")
    for h in hits[:5]:
        print(f"  {h.seq_id}:{h.start}-{h.end}({h.strand})  "
              f"score={h.score:.2f}  p={h.pvalue:.2e}  q={h.qvalue:.2e}  {h.matched}")
    assert any(h.start == 200 and h.strand == "+" for h in hits), "planted site not recovered"
    # reverse-complement of the planted site should also be detectable on the - strand
    rc = "TGTCAA"
    planted2 = bg[:300] + rc + bg[300:]
    hits2 = scan({"chr": planted2}, counts=counts, pvalue_thresh=1e-3)
    assert any(h.start == 300 and h.strand == "-" for h in hits2), "minus-strand site missed"
    print("OK: + strand site @200 and - strand site @300 both recovered.")


if __name__ == "__main__":
    _demo()
