"""
operator_logo.py -- the logo built from the OPERATOR SEQUENCES FOUND IN THE GENOME (distinct from the
DeepPBS structural-readout logo), plus a motif-emergence criterion for how many genome hits to keep.

Two products:
  * `build_operator_logo(seqs)` -- a position-frequency logo from the matched operator sequences the
    genome rescan returned. They arrive already in the consensus-PWM frame (the rescan used the
    contact-trimmed DeepPBS consensus), so this logo is implicitly **trimmed by the TF-DNA contacts**;
    we then trim the low-information flanks so only the motif remains.
  * `select_by_motif(ranked_seqs)` -- the data-driven answer to "how many Snowprint/rescan hits should we
    keep?". Real operator hits reinforce a motif (high per-column information); adding false hits dilutes
    it. We add hits top-down and keep the largest prefix whose mean per-column information stays within
    `keep_frac` of the running peak -- i.e. keep hits while the motif holds, stop when it washes out.

Pure / offline. Run `python operator_logo.py` for a self-test (motif emerges from planted hits, dilutes
when random hits are appended, and the chosen cut sits at the motif boundary).
"""
from __future__ import annotations

import collections
from dataclasses import dataclass, field

import numpy as np

from predictor.motifs.information_content import per_column as _ic

_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}
_RC = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}


def _rc(s: str) -> str:
    return "".join(_RC.get(c, "N") for c in reversed(s))


def _consensus(seqs, W):
    c = np.zeros((4, W))
    for s in seqs:
        for j, ch in enumerate(s):
            if ch in _IDX:
                c[_IDX[ch], j] += 1
    return "".join("ACGT"[i] for i in c.argmax(axis=0))


def _match_frac(s, cons):
    return sum(a == b for a, b in zip(s, cons)) / max(1, len(cons))


def _orient_to_consensus(seqs, *, margin: float = 0.1):
    """Flip any sequence whose reverse-complement matches the running consensus CLEARLY better (by
    `margin`). A no-op when the inputs are already in one orientation (and for palindromes, where the two
    orientations tie). Guards against mixed-strand inputs washing out the logo -- see coords.from_genome_hit.
    """
    cons = _consensus(seqs, len(seqs[0]))
    out = []
    for s in seqs:
        r = _rc(s)
        out.append(r if _match_frac(r, cons) > _match_frac(s, cons) + margin else s)
    return out


@dataclass(frozen=True)
class PreparedSequences:
    """The exact equal-width sequence stack used to build a count matrix.

    ``input_sequences`` is the normalized, non-empty input and ``sequences`` is the modal-width,
    orientation-normalized subset that actually contributes.  Keeping both makes the long-standing
    modal-length filtering explicit without changing its numerical behaviour.
    """

    input_sequences: tuple[str, ...] = ()
    sequences: tuple[str, ...] = ()
    modal_w: int | None = None

    @property
    def n_input(self) -> int:
        return len(self.input_sequences)

    @property
    def n_effective(self) -> int:
        return len(self.sequences)


def prepare_sequences(seqs, *, orient: bool = True) -> PreparedSequences:
    """Normalize and select the exact modal-length sequences used by :func:`counts_from_seqs`.

    Modal ties retain ``Counter.most_common``'s first-observed behaviour, matching the previous inline
    implementation exactly.  Orientation is likewise applied only after modal-length filtering.
    """
    normalized = tuple(s.upper() for s in seqs if s)
    if not normalized:
        return PreparedSequences()
    modal_w = collections.Counter(len(s) for s in normalized).most_common(1)[0][0]
    use = [s for s in normalized if len(s) == modal_w]
    if orient and len(use) > 1:
        use = _orient_to_consensus(use)
    return PreparedSequences(
        input_sequences=normalized,
        sequences=tuple(use),
        modal_w=modal_w,
    )


def counts_from_seqs(seqs, *, orient: bool = True) -> np.ndarray | None:
    """4xW count matrix from equal-length sequences (uses the modal length; others are dropped). When
    `orient` (default), each sequence is flipped to the consensus orientation first, so mixed-strand
    inputs still stack in register (defensive -- callers should already orient via coords.from_genome_hit).
    """
    return counts_from_prepared(prepare_sequences(seqs, orient=orient))


def counts_from_prepared(prepared: PreparedSequences) -> np.ndarray | None:
    """Build counts from an already-prepared stack, with no second filtering/orientation pass."""
    if not prepared.sequences or prepared.modal_w is None:
        return None
    c = np.zeros((4, prepared.modal_w))
    for s in prepared.sequences:
        for j, ch in enumerate(s):
            if ch in _IDX:
                c[_IDX[ch], j] += 1
    return c


@dataclass
class OperatorLogo:
    pwm: np.ndarray | None = None                        # 4xW position-frequency, trimmed to the motif
    per_col_ic: np.ndarray | None = None
    total_ic: float = 0.0
    consensus: str = ""
    n_seqs: int = 0
    trim_span: tuple = (0, 0)                             # (lo, hi) columns kept relative to the input frame


def build_operator_logo(seqs, *, pseudocount: float = 0.25, trim_ic: float = 0.35,
                        min_w: int = 6) -> OperatorLogo:
    """Position-frequency logo from genomic operator sequences, trimmed to the contiguous informative
    core (columns with IC > `trim_ic`) -- "trim until a motif appears"."""
    counts = counts_from_seqs(seqs)
    if counts is None:
        return OperatorLogo()
    n = int(counts.sum(axis=0).max())
    pwm = (counts + pseudocount) / (counts + pseudocount).sum(axis=0, keepdims=True)
    ic = _ic(pwm)
    cols = np.where(ic > trim_ic)[0]
    if len(cols) >= min_w:
        lo, hi = int(cols.min()), int(cols.max()) + 1
    else:
        lo, hi = 0, pwm.shape[1]
    pwm_t = pwm[:, lo:hi]
    ic_t = _ic(pwm_t)
    cons = "".join("ACGT"[i] for i in pwm_t.argmax(axis=0))
    return OperatorLogo(pwm=pwm_t, per_col_ic=ic_t, total_ic=float(ic_t.sum()), consensus=cons,
                        n_seqs=n, trim_span=(lo, hi))


# --------------------------------------------------------------------------- motif-AWARE logo (WS3)
# ONE width search range for every family. It spans the operator widths seen across all of them: a compact
# dyad with a short spacer at the low end (~12 bp) up to a long dyad with a ~19-bp spacer at the high end
# (~41 bp). Two families used to carry narrower bespoke ranges, which meant the EM was allowed to consider
# widths for one TF that it was forbidden from considering for another. The EM already selects the width by
# likelihood, so letting every family search the same range is both simpler and less presumptuous.
OPERATOR_WIDTHS = range(12, 42)


def _trimmed_logo(pwm: np.ndarray, *, n_seqs: int, trim_ic: float, min_w: int) -> OperatorLogo:
    ic = _ic(pwm)
    cols = np.where(ic > trim_ic)[0]
    lo, hi = (int(cols.min()), int(cols.max()) + 1) if len(cols) >= min_w else (0, pwm.shape[1])
    pwm_t = pwm[:, lo:hi]
    ic_t = _ic(pwm_t)
    cons = "".join("ACGT"[i] for i in pwm_t.argmax(axis=0))
    return OperatorLogo(pwm=pwm_t, per_col_ic=ic_t, total_ic=float(ic_t.sum()), consensus=cons,
                        n_seqs=n_seqs, trim_span=(lo, hi))


def build_aligned_logo(seqs, *, widths=None, palindrome: bool = True,
                       trim_ic: float = 0.35, min_w: int = 6, restarts: int = 4,
                       return_motif: bool = False):
    """Motif-AWARE operator logo: FIND + ALIGN the dyad-symmetric operator inside the input sequences --
    which may be long, unaligned promoter regions, NOT pre-framed operator instances -- then build a trimmed
    logo from the aligned PWM. Replaces the naive fixed-width column stack of `build_operator_logo` (which
    assumes inputs already share a frame and washes out on raw promoters).

    Alignment is `motifs.em_finder` (OOPS EM, BOTH strands, palindrome-symmetrized M-step) -- the
    offset+strand+palindrome-aware alignment the operators need -- over `OPERATOR_WIDTHS` for every TF."""
    seqs = [s.upper() for s in seqs if s and len(s) >= min_w]
    if len(seqs) < 2:
        lg = build_operator_logo(seqs)
        return (lg, None) if return_motif else lg
    from predictor.motifs import em_finder as emf
    m = emf.discover(seqs, widths=(widths if widths is not None else OPERATOR_WIDTHS),
                     palindrome=palindrome, restarts=restarts)
    if m is None:
        lg = build_operator_logo(seqs)
        return (lg, None) if return_motif else lg
    lg = _trimmed_logo(m.prob, n_seqs=len(seqs), trim_ic=trim_ic, min_w=min_w)
    return (lg, m) if return_motif else lg


@dataclass
class MultiSourceLogos:
    """The three operator logos the redesign calls for, each from a different sequence pool."""
    tf_only: OperatorLogo = field(default_factory=OperatorLogo)       # the TF's own genome rescan hits
    homolog: OperatorLogo = field(default_factory=OperatorLogo)       # operators within homolog promoters
    gathered: OperatorLogo = field(default_factory=OperatorLogo)      # operators within the full MSA/SSN set


def aligned_instances(regions, motif, *, min_w: int = 6) -> list[str]:
    """Reconstruct the ALIGNED operator instances an EM `motif` found inside `regions` -- the actual
    sequence strings stacked to build the logo (the spec's 'FASTA of the aligned operators'). Re-applies
    the SAME filter `build_aligned_logo` used so `motif.sites` (seq_index into the filtered list) line up,
    then slices each site to `motif.width`, reverse-complementing '-'-strand sites into the logo frame."""
    if motif is None:
        return []
    seqs = [s.upper() for s in regions if s and len(s) >= min_w]
    W = int(getattr(motif, "width", 0))
    out = []
    for (si, start, strand) in getattr(motif, "sites", []):
        if not (0 <= si < len(seqs)):
            continue
        inst = seqs[si][start:start + W]
        if len(inst) != W:
            continue
        out.append(_rc(inst) if strand == "-" else inst)
    return out


def three_logos(tf_hits, homolog_regions, gathered_regions, *, return_aligned: bool = False):
    """Build the three-source operator logos. `tf_hits` are in-frame rescan operator instances (naive stack);
    `homolog_regions` and `gathered_regions` are long promoter windows -> motif-aware EM alignment.

    With `return_aligned=True`, also return the aligned operator INSTANCES behind each logo
    (dict tf_only/homolog/gathered -> list[str]) so the caller can export them as FASTA -- the EM alignment
    is reconstructed from each logo's `Motif.sites` (which is otherwise discarded)."""
    tf_seqs = [s for s in tf_hits if s] if tf_hits else []
    if return_aligned:
        if homolog_regions:
            h_lg, h_m = build_aligned_logo(homolog_regions, return_motif=True)
        else:
            h_lg, h_m = OperatorLogo(), None
        if gathered_regions:
            g_lg, g_m = build_aligned_logo(gathered_regions, return_motif=True)
        else:
            g_lg, g_m = OperatorLogo(), None
        logos = MultiSourceLogos(
            tf_only=build_operator_logo(tf_seqs) if tf_seqs else OperatorLogo(),
            homolog=h_lg, gathered=g_lg)
        aligned = {"tf_only": tf_seqs,
                   "homolog": aligned_instances(homolog_regions, h_m) if homolog_regions else [],
                   "gathered": aligned_instances(gathered_regions, g_m) if gathered_regions else []}
        return logos, aligned
    return MultiSourceLogos(
        tf_only=build_operator_logo(tf_seqs) if tf_seqs else OperatorLogo(),
        homolog=build_aligned_logo(homolog_regions) if homolog_regions else OperatorLogo(),
        gathered=build_aligned_logo(gathered_regions) if gathered_regions else OperatorLogo())


@dataclass
class MotifSelection:
    n_keep: int = 0
    logo: OperatorLogo = field(default_factory=OperatorLogo)
    ic_curve: list = field(default_factory=list)         # mean per-column IC as each hit is added
    peak_n: int = 0


def select_by_motif(ranked_seqs, *, min_hits: int = 3, max_hits: int = 40, keep_frac: float = 0.75,
                    trim_ic: float = 0.35) -> MotifSelection:
    """Choose how many genome hits to keep by motif emergence. `ranked_seqs` are the matched operator
    sequences, best first. Add them top-down; track mean per-column IC of the trimmed logo; keep the
    LARGEST prefix whose mean IC stays >= `keep_frac` x the running peak (the motif still holds). Returns
    the chosen count + the logo built from exactly those hits.

    **`max_hits` BINDS in practice -- it is not a generous safety bound.** Measured over the 17 panel
    bundles (`project_archive/benchmarking_regulondb/max_hits_probe.py`), `n_keep` sits at the 40 ceiling in
    **17/17**: at `keep_frac=0.75` the walk never finds a natural stopping point inside 400 hits, so
    this constant, not the emergence rule, is what decides how many hits build the motif.

    Raising it degrades the logo monotonically on average -- mean per-column IC changes -0.040 bits at
    80, -0.097 at 160 and -0.157 at 400, where every one of the 17 is worse. So the deep tail really is
    mostly noise for motif building.

    Two caveats before anyone retunes this on that evidence:

    * **IC is a SELF-FIT score.** A logo built from fewer sequences is sharper by construction, so
      "IC falls as hits are added" partly just restates that. Judging this constant by IC repeats the
      mistake that retired the mean-IC seed selector (handoff §9.1). The honest criterion is operator
      RECOVERY, which needs the scorer, not this function.
    * **`n_keep` is not cosmetic.** `pipeline` builds the regulon from `ranked[:n_keep]`, so changing
      it changes which operons are proposed -- not just the picture in the report.

    Nothing is thrown away either way: every hit stays in `rescan.hits` and `binding_sites.tsv`. This
    chooses only which prefix defines the motif.
    """
    seqs = [s for s in ranked_seqs if s]
    if len(seqs) < min_hits:
        lg = build_operator_logo(seqs, trim_ic=trim_ic)
        return MotifSelection(n_keep=len(seqs), logo=lg, ic_curve=[lg.total_ic / max(1, lg.pwm.shape[1])]
                              if lg.pwm is not None else [], peak_n=len(seqs))
    curve, mean_ics = [], []
    cap = min(max_hits, len(seqs))
    for n in range(1, cap + 1):
        lg = build_operator_logo(seqs[:n], trim_ic=trim_ic)
        mic = (lg.total_ic / lg.pwm.shape[1]) if (lg.pwm is not None and lg.pwm.shape[1]) else 0.0
        mean_ics.append(mic)
        curve.append(round(mic, 3))
    # peak (ignoring the tiny-n regime) then walk out while the motif holds
    start = min_hits - 1
    peak_i = start + int(np.argmax(mean_ics[start:])) if len(mean_ics) > start else 0
    peak = mean_ics[peak_i]
    n_keep = peak_i + 1
    for n in range(peak_i + 1, cap):
        if mean_ics[n] >= keep_frac * peak:
            n_keep = n + 1
        else:
            break
    n_keep = max(min_hits, n_keep)
    logo = build_operator_logo(seqs[:n_keep], trim_ic=trim_ic)
    return MotifSelection(n_keep=n_keep, logo=logo, ic_curve=curve, peak_n=peak_i + 1)


# --------------------------------------------------------------------------- self-test
def _demo() -> None:
    import random
    rng = random.Random(0)
    motif = "TTGACATGAACATATTCAGTCAA"                    # a dyad-ish operator
    W = len(motif)

    def jitter(s, p=0.12):                               # a real-operator instance with a few mismatches
        return "".join(c if rng.random() > p else rng.choice("ACGT") for c in s)

    true = [jitter(motif) for _ in range(8)]             # 8 genuine operator hits
    noise = ["".join(rng.choice("ACGT") for _ in range(W)) for _ in range(12)]   # 12 false hits

    # logo from the true hits alone: a clear motif
    lg = build_operator_logo(true)
    print(f"true-only logo: consensus={lg.consensus} W={lg.pwm.shape[1]} totalIC={lg.total_ic:.1f} "
          f"meanIC={lg.total_ic / lg.pwm.shape[1]:.2f}")
    assert lg.total_ic / lg.pwm.shape[1] > 0.8, "a real operator set should give an informative logo"

    # orientation defense: RC half the hits (as a buggy mixed-strand input would) -> the oriented logo
    # must still recover the motif, whereas stacking them as-is washes it out
    mixed = true[:4] + [_rc(s) for s in true[4:]]
    ic_oriented = build_operator_logo(mixed).total_ic / max(1, build_operator_logo(mixed).pwm.shape[1])
    c_raw = counts_from_seqs(mixed, orient=False)
    pwm_raw = (c_raw + 0.25) / (c_raw + 0.25).sum(axis=0, keepdims=True)
    ic_raw = float(_ic(pwm_raw).mean())
    print(f"mixed-strand: oriented meanIC={ic_oriented:.2f} vs un-oriented meanIC={ic_raw:.2f}")
    assert ic_oriented > ic_raw + 0.2, "orientation defense did not restore the washed-out logo"

    # ranked list: 8 true then 12 noise -> the selector should keep ~8 and stop before the noise washes it out
    sel = select_by_motif(true + noise, min_hits=3, keep_frac=0.75)
    print(f"select_by_motif: n_keep={sel.n_keep} (peak@{sel.peak_n})  consensus={sel.logo.consensus}")
    print(f"  IC curve (mean bits/col): {sel.ic_curve}")
    assert 5 <= sel.n_keep <= 11, f"should keep roughly the true hits, kept {sel.n_keep}"
    # the kept logo still resembles the motif (majority of columns match)
    match = sum(a == b for a, b in zip(sel.logo.consensus, motif[:len(sel.logo.consensus)]))
    assert match >= 0.6 * len(sel.logo.consensus), "kept logo lost the motif"

    # all-noise -> low IC, selector keeps the minimum and the logo is flat
    sel0 = select_by_motif(noise, min_hits=3)
    flat_ic = sel0.logo.total_ic / max(1, sel0.logo.pwm.shape[1]) if sel0.logo.pwm is not None else 0
    print(f"all-noise: n_keep={sel0.n_keep} meanIC={flat_ic:.2f}")
    assert flat_ic < 0.6, "random sequences should not form an informative motif"

    # WS3: motif-AWARE logo on LONG UNALIGNED promoter regions (varying offsets) -- naive stacking fails here
    pal = "TTGACATGCATGTCAA"                              # 16-bp palindrome (rc == itself)
    rng2 = random.Random(7)
    promoters = []
    for _ in range(8):
        n1, n2 = rng2.randint(80, 160), rng2.randint(80, 160)
        f1 = "".join(rng2.choice("ACGT") for _ in range(n1))
        f2 = "".join(rng2.choice("ACGT") for _ in range(n2))
        site = "".join(c if rng2.random() > 0.1 else rng2.choice("ACGT") for c in pal)
        promoters.append(f1 + site + f2)
    naive = build_operator_logo(promoters)
    aligned = build_aligned_logo(promoters)
    nmic = naive.total_ic / max(1, naive.pwm.shape[1]) if naive.pwm is not None else 0.0
    amic = aligned.total_ic / max(1, aligned.pwm.shape[1]) if aligned.pwm is not None else 0.0
    print(f"WS3 unaligned promoters: naive meanIC={nmic:.2f} vs EM-aligned meanIC={amic:.2f} "
          f"consensus={aligned.consensus}")
    assert amic > 0.8 and amic > nmic + 0.3, "motif-aware EM alignment must beat naive stacking on raw promoters"

    # WS3: three_logos builds all three pools (tf-only stack + EM-aligned homolog/gathered)
    tl = three_logos(true, promoters, promoters)
    assert tl.tf_only.pwm is not None and tl.homolog.pwm is not None and tl.gathered.pwm is not None
    print(f"WS3 three_logos: tf_only IC={tl.tf_only.total_ic:.1f}  homolog IC={tl.homolog.total_ic:.1f}  "
          f"gathered IC={tl.gathered.total_ic:.1f}")

    print("OK: genomic-operator logo emerges from real hits, trims to the motif, and the motif-emergence "
          "selector keeps the operator-supported hits and stops when noise dilutes the signal.")


if __name__ == "__main__":
    _demo()
