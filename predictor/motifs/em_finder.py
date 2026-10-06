"""
em_finder.py -- de-novo motif discovery by EM, pure Python (numpy).

Replaces MEME-suite **MEME** for the de-novo floor (Phase 1) so the pipeline needs no MEME
binary (no WSL / Docker). It fits a single ungapped motif of fixed width by Expectation-
Maximization under a one-occurrence-per-sequence (OOPS) model, scanning BOTH strands, with an
optional **palindrome constraint** (symmetrize the PWM with its reverse complement each
M-step) -- which is exactly right for the dyad-symmetric operators this project targets and is
something stock MEME does not do natively.

Output is a probability matrix (4xW) consumable by `pwm_scan.scan(prob=...)`, plus the per-
sequence best site. For multiple widths, call `discover()` over a width range and keep the best
log-likelihood. EM is for modest sequence sets (e.g. homolog inter-operon regions), not genome
scale -- genome scanning is `pwm_scan`.

Run `python em_finder.py` for a self-test (plants palindromic sites and recovers the motif).
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

try:                                    # works both as `python em_finder.py` and as a package
    from .pwm_scan import _encode, IDX, _COMP, log_odds  # noqa: F401 - IDX re-exported for callers
except ImportError:
    from .pwm_scan import _encode, _COMP, log_odds


@dataclass
class Motif:
    prob: np.ndarray            # (4, W) probability matrix
    width: int
    loglik: float
    sites: list                 # [(seq_index, start, strand)]
    background: np.ndarray


def _window_scores(lo: np.ndarray, codes: np.ndarray) -> np.ndarray:
    """Log-odds score of every length-W window; N positions contribute 0."""
    _, W = lo.shape
    n = codes.shape[0] - W + 1
    if n <= 0:
        return np.empty(0)
    out = np.zeros(n)
    for j in range(W):
        c = codes[j:j + n]
        bad = c < 0
        out += np.where(bad, 0.0, lo[np.where(bad, 0, c), j])
    return out


def _symmetrize(prob: np.ndarray) -> np.ndarray:
    """Average a PWM with its reverse complement -> perfectly palindromic motif."""
    return 0.5 * (prob + prob[_COMP][:, ::-1])


def fit(sequences, width: int, *, palindrome: bool = True, n_iter: int = 200,
        pseudocount: float = 0.5, tol: float = 1e-4, seed: int = 0) -> Motif:
    """Fit one width-`width` motif by OOPS EM over both strands."""
    rng = np.random.default_rng(seed)
    codes = [_encode(s) for s in sequences]
    bg = np.concatenate(codes)
    bg = bg[bg >= 0]
    background = np.array([(bg == a).mean() for a in range(4)], dtype=float)
    background = np.where(background <= 0, 1e-9, background)
    background /= background.sum()

    # ---- init from a random valid window of a random sequence
    usable = [k for k, c in enumerate(codes) if c.shape[0] >= width]
    if not usable:
        raise ValueError("all sequences shorter than width")
    s0 = codes[rng.choice(usable)]
    p0 = int(rng.integers(0, s0.shape[0] - width + 1))
    counts = np.full((4, width), pseudocount)
    for j in range(width):
        b = s0[p0 + j]
        if b >= 0:
            counts[b, j] += 1.0
    prob = counts / counts.sum(axis=0, keepdims=True)
    if palindrome:
        prob = _symmetrize(prob)

    prev = -np.inf
    sites: list = []
    for _ in range(n_iter):
        lo = log_odds(prob, background, log_base=np.e)   # natural log for EM math
        lo_rc = lo[_COMP][:, ::-1]
        expected = np.full((4, width), pseudocount)
        loglik = 0.0
        sites = []
        for si, c in enumerate(codes):
            n = c.shape[0] - width + 1
            if n <= 0:
                continue
            sp = _window_scores(lo, c)
            sm = _window_scores(lo_rc, c)
            allsc = np.concatenate([sp, sm])            # 2n responsibilities (OOPS)
            m = allsc.max()
            w = np.exp(allsc - m)
            tot = w.sum()
            loglik += m + np.log(tot)
            r = w / tot
            rp, rm = r[:n], r[n:]

            best = int(np.argmax(allsc))
            sites.append((si, best % n, "+" if best < n else "-"))

            for j in range(width):
                base = c[j:j + n]
                ok = base >= 0
                # + strand: window k, offset j -> column j
                np.add.at(expected[:, j], base[ok], rp[ok])
                # - strand: same windows, complemented base into mirror column
                np.add.at(expected[:, width - 1 - j], _COMP[base[ok]], rm[ok])

        prob = expected / expected.sum(axis=0, keepdims=True)
        if palindrome:
            prob = _symmetrize(prob)
        if np.isfinite(prev) and abs(loglik - prev) < tol * max(1.0, abs(prev)):
            break
        prev = loglik

    return Motif(prob=prob, width=width, loglik=prev, sites=sites, background=background)


def discover(sequences, widths=range(12, 40), *, palindrome: bool = True,
             restarts: int = 4, **kw) -> Motif:
    """Try a range of widths and EM restarts; return the highest-loglik motif (length-normalized so
    widths are comparable). The default width range spans real operator sizes for both archetypes --
    ArsR/SmtB ~14-32 bp and MerR ~28-39 bp -- without reaching the 40-47 bp double operators. Widths
    longer than the shortest input sequence are skipped (no valid EM window)."""
    best = None
    best_key = -np.inf
    lens = [len(s) for s in sequences]
    n_seq = len(lens)
    for W in widths:
        if sum(1 for L in lens if L >= W) < min(3, n_seq):   # need a few sequences to host this width
            continue
        for r in range(restarts):
            try:
                m = fit(sequences, W, palindrome=palindrome, seed=r, **kw)
            except ValueError:
                break
            key = m.loglik / W                      # normalize for width comparison
            if key > best_key:
                best, best_key = m, key
    return best


# ----------------------------------------------------------------------------- self-test
def _demo() -> None:
    import random
    from .pwm_scan import scan

    random.seed(1)
    core = "TTGACCTAGGTCAA"            # a palindrome (rc == itself)
    assert core == core.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    seqs = []
    for _ in range(12):
        flank1 = "".join(random.choice("ACGT") for _ in range(60))
        flank2 = "".join(random.choice("ACGT") for _ in range(60))
        # implant with a couple of point mutations to make it non-trivial
        site = list(core)
        for _ in range(2):
            site[random.randrange(len(site))] = random.choice("ACGT")
        seqs.append(flank1 + "".join(site) + flank2)

    m = fit(seqs, width=len(core), palindrome=True, seed=0)
    consensus = "".join("ACGT"[i] for i in m.prob.argmax(axis=0))
    print(f"width={m.width}  loglik={m.loglik:.1f}")
    print(f"true core : {core}")
    print(f"recovered : {consensus}")
    # the recovered consensus should match the planted core at most positions
    matches = sum(a == b for a, b in zip(core, consensus))
    print(f"consensus match: {matches}/{len(core)}")
    assert matches >= len(core) - 3, "EM failed to recover the planted palindrome"

    # the recovered PWM should score the implant sites strongly via pwm_scan
    hits = scan({f"s{i}": s for i, s in enumerate(seqs)}, prob=m.prob, pvalue_thresh=1e-3)
    found = {h.seq_id for h in hits}
    print(f"pwm_scan re-detected the motif in {len(found)}/{len(seqs)} sequences.")
    assert len(found) >= len(seqs) - 2, "recovered PWM does not rediscover its own sites"
    print("OK: de-novo palindrome recovered and round-trips through pwm_scan.")


if __name__ == "__main__":
    _demo()
