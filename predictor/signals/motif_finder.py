"""
motif_finder.py -- our best de-novo operator finder (Phase 1, the structure-free floor).

Designed from Snowprint's algorithm and its documented failure modes (see WORKPLAN "How Snowprint
works"). Given a promoter / inter-operon region it generates operator candidates from three pluggable
sources and ranks them, optionally boosted by **phylogenetic conservation** across homolog regions:

  inverted_repeat  -- dyad-symmetric palindromes (Snowprint's seed step), imperfect/discontinuous,
                      with FAMILY-AWARE half-site & spacer priors (ArsR/SmtB 4-8/0-6; MerR 5-9/15-21,
                      the 19-bp signature) -- priors Snowprint lacks.
  direct_repeat    -- direct repeats (Snowprint's stated gap).
  family_seeded    -- scan with the curated family consensus PWM (known_operators) via pwm_scan.

Conservation (when homolog inter-operon regions are supplied): a candidate that recurs across homologs
is up-weighted -- the phylogenetic-footprinting signal, without needing an MSA tool. Output candidates
carry a PWM + per-base importance, consumable by motif_rescan / the closed loop.

Run `python motif_finder.py` for a self-test (plants a MerR-style dyad + homologs and recovers it).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from predictor.motifs.pwm_scan import scan, pwm_from_counts
from predictor.motifs.information_content import per_column as _ic
from predictor.annotate import known_operators

_COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}

# Half-site / spacer priors (bp). FAMILY_PRIORS is an explicit experimental ablation; it is generous
# enough to capture full-length operators -- ArsR/SmtB
# now spans the ~6-8-6 dyad (half 6, spacer 8, ~20 bp) AND wider ~12-2-12-ish sites up to ~30 bp (the
# old spacer cap of 6 could NOT find a spacer-8 operator), and MerR keeps its ~19-bp-spacer signature.
# The ArsR upper bound (half 9 + spacer 12 + half 9 = 30 bp) stays below 40-47 bp "double operators", so
# those are captured as two contiguous single dyads (the slider emits both), not merged. DEFAULT_PRIOR is
# the production family-agnostic search. The 2026-07-29 held-out-window test recovered 8/8 known sites in
# its top eight versus 7/8 with family geometry (and 5/8 versus 3/8 at rank one), so family geometry must
# not exclude candidates in production.
# Fur: the "Fur box" is classically a 19 bp AT-rich inverted repeat with a 7-1-7 arrangement (two 7 bp
# half-sites separated by 1 bp), so the prior brackets that (half 6-9, spacer 0-3 => 12-21 bp sites).
# NB: Fur proteins can oligomerise along the operator and footprint considerably larger regions than the
# 19 bp core; those are captured as adjacent dyads rather than one wide site, matching how the ArsR prior
# handles its 40-47 bp "double operators". Irr is the known exception -- it reads ICE elements, not a Fur
# box -- so an Irr-cluster query is expected to fit this prior poorly.
FAMILY_PRIORS = {
    "ArsR/SmtB": dict(half=(4, 9), spacer=(0, 12)),
    "MerR": dict(half=(5, 9), spacer=(15, 21)),
    "Fur": dict(half=(6, 9), spacer=(0, 3)),
    # 2026 drop: operator geometries from the literature (see analysis/ssn_final_state/SSN_MASTER.operator).
    # CsoR/RcnR: G/C-rich pseudopalindrome (TtCsoR CCCCACCCCA..TGGGGTGGGG ~22 bp; RcnR TACT-G6N-AGTA).
    "CsoR/FrmR": dict(half=(4, 10), spacer=(0, 6)),
    # NikR: E. coli 28 bp GTATGA-N16-TCATAC -- a 6 bp half-site with a long ~16 bp spacer.
    "NikR": dict(half=(5, 8), spacer=(12, 18)),
    # DtxR/IdeR/MntR: 19 bp interrupted palindrome ("iron box"), tandem CCTAA repeats, two dimers.
    "DtxR/MntR": dict(half=(6, 9), spacer=(0, 3)),
    # MarR/SlyA: winged-HTH inverted repeats (marO; AdcR ADC-box TTAAC-NN-GTTAA ~12 bp core).
    "MarR/SlyA": dict(half=(4, 8), spacer=(0, 6)),
    # CopY/BlaI/MecI: TACAnnTGTA "cop/bla box" inverted repeat (~5 bp half, short spacer), often tandem.
    "CopY": dict(half=(4, 7), spacer=(0, 4)),
    # Rrf2/IscR/NsrR: long inverted repeats (IscR type-1 ~24 bp; ScNsrR two 11 bp IRs + 1 bp).
    "Rrf2": dict(half=(6, 11), spacer=(0, 3)),
}
DEFAULT_PRIOR = dict(half=(4, 9), spacer=(0, 21))


def _rc(s: str) -> str:
    return "".join(_COMP.get(c, "N") for c in reversed(s))


def _complexity(s: str) -> float:
    """Low-complexity penalty in [0.2, 1.0] for a half-site. Homopolymers (AAAA) and period-2 tandem
    repeats (ATAT, GCGC) form spurious *perfect* palindromes that would otherwise outrank real operators;
    they get penalised, while a diverse half-site (e.g. TTGACA) keeps the full weight. The spacer is NOT
    scored (the IR score is computed on the contacted half-sites only), so an AT-rich spacer is unaffected."""
    s = (s or "").upper()
    if len(s) < 2:
        return 1.0
    distinct = len({c for c in s if c in _IDX})
    div = min(1.0, distinct / 3.0)                                   # 1 base->0.33, 2->0.67, 3+->1.0
    per2 = sum(1 for k in range(len(s) - 2) if s[k] == s[k + 2]) / max(1, len(s) - 2)
    rep = 1.0 - 0.6 * per2                                           # full period-2 repeat -> 0.4
    return max(0.2, div * rep)


def _counts(seq: str) -> np.ndarray:
    W = len(seq)
    c = np.zeros((4, W))
    for j, ch in enumerate(seq.upper()):
        if ch in _IDX:
            c[_IDX[ch], j] += 1
        else:
            c[:, j] += 0.25
    return c


@dataclass
class OperatorCandidate:
    center: float
    start: int
    end: int
    seq: str
    score: float                       # sequence_score (repeat quality), pre-conservation
    kind: str                          # 'IR' | 'DR' | 'family'
    half: int | None = None
    spacer: int | None = None
    half_sites: tuple | None = None
    conservation: float = 0.0
    final: float = 0.0
    pwm: np.ndarray | None = None
    per_base_importance: np.ndarray | None = None


# --------------------------------------------------------------------------- candidate generators
#: How a dyad's quality becomes a number. `"raw"` is production and the historical behaviour.
#:
#: `raw` = 2m - h, which for a PERFECT dyad (m = h) is simply **h**: the score is the half-site
#: length, rises with width, and is bounded only by the prior's h ceiling. Measured consequence on
#: the Salmonella `fur` promoter: the top candidate is a perfect 9 bp half-site pair scoring 9.0 and
#: the real Fur box a 6-5-6 scoring 6.0 -- the true site loses on width alone. The snowstream
#: `tetr_137` panel scored this ranking at AUC 0.250 against known operators, i.e. anti-correlated.
#:
#: `surprise` corrects for chance (-log10 of the binomial tail) but is still monotone in h, because a
#: perfect 9-mer dyad genuinely IS rarer than a perfect 6-mer one. It is included to SHOW that, not
#: because it is expected to win.
#:
#: `density` multiplies by the fraction of the site that is actually half-site, 2h/width, so a dyad
#: that spends 12 bp on spacer is discounted against a compact one. This is the only variant that can
#: reorder the Fur case, and it is family-agnostic and composition-blind -- note that COMPOSITION
#: re-ranking is a different idea and has been measured to fail twice (`seed_clusters` docstring).
SCORE_MODES = ("raw", "surprise", "density", "surprise_density")


@lru_cache(maxsize=None)
def _binom_surprise(m: int, h: int, p: float = 0.25) -> float:
    """-log10 P(X >= m) for X ~ Binomial(h, p). Chance-corrected dyad quality."""
    tail = sum(math.comb(h, k) * (p ** k) * ((1 - p) ** (h - k)) for k in range(m, h + 1))
    return -math.log10(max(tail, 1e-300))


def _dyad_score(m: int, h: int, width: int, mode: str) -> float:
    if mode == "raw":
        return float(2 * m - h)
    if mode == "surprise":
        return _binom_surprise(m, h)
    density = (2.0 * h) / max(1, width)
    if mode == "density":
        return float(2 * m - h) * density
    if mode == "surprise_density":
        return _binom_surprise(m, h) * density
    raise ValueError(f"unknown score mode {mode!r}; expected one of {SCORE_MODES}")


def find_inverted_repeats(seq, *, half_range, spacer_range, min_frac=0.75, score_mode="raw"):
    seq = seq.upper()
    L = len(seq)
    out = []
    for h in range(half_range[0], half_range[1] + 1):
        for i in range(0, L - 2 * h + 1):
            left = seq[i:i + h]
            if "N" in left:
                continue
            for s in range(spacer_range[0], spacer_range[1] + 1):
                j = i + h + s
                if j + h > L:
                    break
                right = seq[j:j + h]
                if "N" in right:
                    continue
                m = sum(a == b for a, b in zip(left, _rc(right)))
                if m / h >= min_frac:
                    out.append(OperatorCandidate(
                        center=(i + j + h) / 2.0, start=i, end=j + h, seq=seq[i:j + h],
                        score=_dyad_score(m, h, 2 * h + s, score_mode),
                        kind="IR", half=h, spacer=s,
                        half_sites=((i, i + h), (j, j + h))))
    return out


def find_direct_repeats(seq, *, half_range, spacer_range, min_frac=0.8, score_mode="raw"):
    seq = seq.upper()
    L = len(seq)
    out = []
    for h in range(half_range[0], half_range[1] + 1):
        for i in range(0, L - 2 * h + 1):
            left = seq[i:i + h]
            if "N" in left:
                continue
            for s in range(spacer_range[0], spacer_range[1] + 1):
                j = i + h + s
                if j + h > L:
                    break
                right = seq[j:j + h]
                if "N" in right:
                    continue
                m = sum(a == b for a, b in zip(left, right))
                if m / h >= min_frac:
                    out.append(OperatorCandidate(
                        center=(i + j + h) / 2.0, start=i, end=j + h, seq=seq[i:j + h],
                        score=_dyad_score(m, h, 2 * h + s, score_mode) * 0.9,  # IR preferred
                        kind="DR", half=h, spacer=s,
                        half_sites=((i, i + h), (j, j + h))))
    return out


def family_seeded(seq, family, *, exclude_tf=None, strict_exclusion=False):
    if not family:
        return []
    pwm, _n = known_operators.family_consensus(
        family, exclude_tf=exclude_tf, strict_exclusion=strict_exclusion,
    )
    if pwm is None:
        return []
    out = []
    for h in scan(seq, prob=pwm, pvalue_thresh=0.05, both_strands=True):
        out.append(OperatorCandidate(
            center=(h.start + h.end) / 2.0, start=h.start, end=h.end,
            seq=_rc(h.matched) if h.strand == "-" else h.matched,
            score=-math.log10(max(h.pvalue, 1e-300)), kind="family"))
    return out


def pwm_seeded(seq, seed_pwm, *, kind="seeded", pvalue_thresh=0.05):
    """Scan the region with an externally-supplied PWM (e.g. the DeepPBS structure-readout, or any
    learned motif) and emit candidates. This is how a STRUCTURE-derived base-specificity profile seeds
    the de-novo operator search -- the structurally-read positions point the search at the operator."""
    if seed_pwm is None:
        return []
    out = []
    for h in scan(seq, prob=seed_pwm, pvalue_thresh=pvalue_thresh, both_strands=True):
        out.append(OperatorCandidate(
            center=(h.start + h.end) / 2.0, start=h.start, end=h.end,
            seq=_rc(h.matched) if h.strand == "-" else h.matched,
            score=-math.log10(max(h.pvalue, 1e-300)), kind=kind))
    return out


# --------------------------------------------------------------------------- conservation
#: Positional drift tolerated between homologs when deciding whether their hits are "the same site".
#: An operator's spacing to its gene is constrained by the promoter architecture it sits in, so genuine
#: orthologous sites land at a similar distance from the gene start; +/-`_POS_TOL` bp is generous enough
#: to absorb indels in the intergenic region without letting an unrelated hit 200 bp away count as the
#: same site. Regions are compared by DISTANCE FROM THE GENE-PROXIMAL (3') END, which is invariant to the
#: 5' clipping that happens near a contig edge.
_POS_TOL = 50

#: Two promoters sharing more than this fraction of their k-mers are treated as the same observation.
_DUP_JACCARD = 0.90
_DUP_K = 8


def _kmers(seq, k=_DUP_K):
    s = (seq or "").upper()
    return {s[i:i + k] for i in range(max(0, len(s) - k + 1))}


def independence_weights(regions):
    """Weight each region by 1 / (size of its near-duplicate group).

    Homolog promoters are sourced from cluster members' genomes, and public genome sets are dominated by
    heavily re-sequenced organisms: a dozen near-identical strains of one species are a dozen copies of
    ONE evolutionary observation. Counting them as a dozen independent confirmations is what makes a
    plain "fraction of regions with a hit" score reward the most re-sequenced clade rather than the most
    conserved site. Grouping by k-mer Jaccard is order- and length-tolerant, so it also catches the same
    locus fetched at slightly different window bounds."""
    n = len(regions)
    sets = [_kmers(r) for r in regions]
    group = list(range(n))                      # union-find over near-duplicates

    def find(i):
        while group[i] != i:
            group[i] = group[group[i]]
            i = group[i]
        return i

    for i in range(n):
        if not sets[i]:
            continue
        for j in range(i + 1, n):
            if not sets[j]:
                continue
            inter = len(sets[i] & sets[j])
            if not inter:
                continue
            if inter / len(sets[i] | sets[j]) >= _DUP_JACCARD:
                group[find(i)] = find(j)
    sizes = {}
    for i in range(n):
        sizes[find(i)] = sizes.get(find(i), 0) + 1
    return [1.0 / sizes[find(i)] for i in range(n)]


def conservation(cand_seq, homolog_regions, *, alpha=0.05):
    """Phylogenetic-footprinting score: how much INDEPENDENT evidence puts this candidate at the SAME
    place in the homolog promoters. Returns 0..1 (the weighted fraction of independent homolog groups
    whose best significant hit is positionally coherent with the others).

    Three things the earlier "fraction of regions containing any significant hit" score got wrong, all
    of which inflate it for non-operators:

      * **Position was ignored.** A hit anywhere in a ~380 bp window counted. Regions are now compared by
        distance from the gene-proximal end -- invariant to 5' clipping at a contig edge -- and only hits
        that agree on a position (within `_POS_TOL`) count as the same conserved site. This requires the
        transcriptional-orientation convention (`genome_resolver.orient_promoter`); before it, plus- and
        minus-strand promoters were mirrored and no positional comparison was meaningful.
      * **Redundancy counted as evidence.** Near-identical strains are collapsed by
        `independence_weights`, so a site conserved across one over-sequenced species no longer outscores
        one conserved across genuinely diverse genomes.
      * **Only the count of regions mattered**, so a candidate hitting many regions at scattered,
        unrelated offsets scored the same as a real, positionally fixed operator.

    Each region is still tested at its own **Bonferroni-corrected** threshold `alpha / len(region)`, so a
    soft single-sequence PWM cannot hit any long region by chance. Scans stay both-strands: a palindromic
    operator reads the same either way, so requiring strand agreement would penalise exactly the dyads
    this finder is built to recover.

    MEASURED: this scores differently from the old version per candidate, but on the RegulonDB panel it
    produces BYTE-IDENTICAL genome results -- same sites, same true sites, same top-decile precision, and
    the seed PWM changed for 0 of 7 TFs. `seed_seqs` ranks by `score * (1 + conservation)`, so with
    conservation in [0,1] the repeat-quality `score` dominates and no reweighting inside that range
    reorders the candidates; forcing conservation to dominate instead cost 28% of true sites for a
    top-decile gain confined to the SSN-anchored TFs. It is kept for the bias it removes (re-sequenced
    strains counted as independent evidence), which that panel cannot detect because E. coli SSN clusters
    are diverse -- NOT for a demonstrated improvement. Do not report it as one.
    See the conservation ablation, recorded in the analysis repo before commit 5815559.
    """
    if not homolog_regions or len(cand_seq) < 6:
        return 0.0
    pwm = pwm_from_counts(_counts(cand_seq))
    # One pooled scan for the whole homolog set: rebuilding the exact PSSM null per region per candidate
    # (~60 candidates x 25 homologs) made this take minutes, and a pooled background is the correct
    # comparison basis across homologs anyway.
    sequences = {f"h{i}": r for i, r in enumerate(homolog_regions)}
    thresholds = {sid: alpha / max(1, len(seq)) for sid, seq in sequences.items()}
    hits = scan(sequences, prob=pwm, pvalue_thresh=max(thresholds.values()), both_strands=True)

    # best significant hit per region, as distance from the gene-proximal (3') end
    best: dict = {}
    for h in hits:
        if h.pvalue > thresholds.get(h.seq_id, 0.0):
            continue
        prev = best.get(h.seq_id)
        if prev is None or h.pvalue < prev[0]:
            region_len = len(sequences[h.seq_id])
            centre = (h.start + h.end) / 2.0
            best[h.seq_id] = (h.pvalue, region_len - centre)
    if not best:
        return 0.0

    weights = independence_weights(list(homolog_regions))
    w_of = {f"h{i}": w for i, w in enumerate(weights)}
    total = sum(weights) or 1.0

    # Positional coherence: keep the ±_POS_TOL window carrying the most independent weight. Sweeping the
    # observed distances is enough -- an optimal window can always be anchored on one of them.
    scored = [(d, w_of[sid]) for sid, (_p, d) in best.items()]
    best_w = 0.0
    for centre, _w in scored:
        acc = sum(w for d, w in scored if abs(d - centre) <= _POS_TOL)
        best_w = max(best_w, acc)
    return best_w / total


def _dedupe(cands, min_sep=4):
    kept = []
    for c in sorted(cands, key=lambda x: -x.score):
        if all(abs(c.center - k.center) > min_sep for k in kept):
            kept.append(c)
    return kept


# --------------------------------------------------------------------------- public API
def predict(region_seq, *, family=None, homolog_regions=None, top_k=8, seed_pwm=None,
            seed_kind="seeded", family_informed=False, use_family_seed=False,
            exclude_family_tf=None, score_mode="raw"):
    """Return ranked OperatorCandidates for `region_seq`. The IR/DR width search uses the shared
    family-agnostic `DEFAULT_PRIOR` by default. Pass `family_informed=True` only for an explicit
    family-geometry ablation; family geometry must not exclude candidates in production.
    `use_family_seed=True` explicitly mixes in the curated family PWM. Production keeps it separate:
    the leave-protein-out panel covered only 11/17 TFs and recovered 7 true sites among 1,413 calls.
    `homolog_regions` enables the phylogenetic-conservation boost; `seed_pwm` adds structure-seeded
    candidates. All genome rescan hits are kept downstream regardless of `top_k`."""
    prior = (FAMILY_PRIORS.get(family, DEFAULT_PRIOR) if family_informed else DEFAULT_PRIOR)
    cands = (find_inverted_repeats(region_seq, half_range=prior["half"], spacer_range=prior["spacer"],
                                   score_mode=score_mode)
             + find_direct_repeats(region_seq, half_range=prior["half"], spacer_range=prior["spacer"],
                                   score_mode=score_mode)
             + (family_seeded(region_seq, family, exclude_tf=exclude_family_tf)
                if use_family_seed else [])
             + pwm_seeded(region_seq, seed_pwm, kind=seed_kind))
    cands = _dedupe(cands)
    for c in cands:
        c.conservation = conservation(c.seq, homolog_regions) if homolog_regions else 0.0
        c.final = c.score * (1.0 + c.conservation)        # conserved candidates win
        c.pwm = pwm_from_counts(_counts(c.seq))
        c.per_base_importance = _ic(c.pwm)                 # per-column information content (bits)
    cands.sort(key=lambda x: -x.final)
    return cands[:top_k]


def at_content(seq: str) -> float:
    """Fraction of A/T among the ACGT bases. Ambiguity codes are ignored rather than counted as GC."""
    s = (seq or "").upper()
    n = sum(1 for c in s if c in "ACGT")
    return (sum(1 for c in s if c in "AT") / n) if n else 0.0


def seed_clusters(cands, *, tol=3, max_clusters=None):
    """Partition the candidates into width-coherent seed sets -- ALL of them, best-ranked first.

    This replaces the single-anchor gamble that `seed_seqs` documents. The old step picked one anchor
    width and deleted every candidate outside +/- `tol` of it, so a wrong winner did not rank the right
    answer badly, it removed it (TtgR, finding 3 of the TetR-137 benchmark, recorded in the analysis
    repo before commit 5815559). Re-ranking the anchor does not fix that, it only moves who loses:
    measured on the RegulonDB panel (`project_archive/benchmarking_regulondb/probe_seed_selection.tsv`),
    ranking the anchor by composition rescued
    TtgR and ModE and *evicted CueR's true operator*, because CueR's Cu(I) operator is 19 bp at AT 0.368
    while the winning candidate was 33 bp.

    So no candidate is discarded here at all. Each width cluster becomes its own seed set, the caller
    rescans the genome with each, and the DATA -- which cluster's genome-wide hits form a coherent
    motif -- decides, rather than a dyad score measured to have AUC 0.250 against known operators.

    Clusters are ordered by their best member's `final`, i.e. the long-standing ordering, so cluster 0
    is exactly what the old code would have chosen. Returns `[(seqs, anchor_width), ...]`.

    Ordering deliberately does NOT use a composition prior. That was tried and **measured worse**: over
    the 72 RegulonDB sites the finders cover, scoring the anchor by `final x AT` retains 33 and by
    `final x composition_prior` 39, against 37 for plain `final` -- and on the metal families
    specifically, 15 and 20 against **23**. Composition is a real family property (Fur 0.767 mean
    operator AT vs MerR 0.554) but it does not rank anchors better, so it is not wired in. The ordering
    is in any case immaterial to the shipped result: `pipeline.run_novel` rescans EVERY cluster and
    picks by genome-wide motif coherence, so this only decides cluster 0 and any `max_clusters` cut.

    `max_clusters` defaults to None -- **do not set it lightly**. Any cap reintroduces the very bug this
    function exists to remove, just further down the ranking: with a cap of 3, CueR's true 19 bp operator
    lands in cluster 3 and is deleted again (measured, `probe_seed_selection.py`). It is exposed only
    because each cluster costs the caller one genome rescan, so a caller under a hard time budget can
    trade correctness for speed KNOWINGLY. `predict(top_k=8)` bounds this at 8 clusters; 3-5 is typical."""
    cs = [c for c in cands if getattr(c, "seq", None)]
    if not cs:
        return []
    dyads = [c for c in cs if c.kind in ("IR", "DR")]
    base = dyads or cs
    conserved = [c for c in base if getattr(c, "conservation", 0.0) > 0]

    pool = sorted(conserved or base, key=lambda c: -getattr(c, "final", c.score))
    out, remaining = [], list(pool)
    while remaining and (max_clusters is None or len(out) < max_clusters):
        anchor = remaining[0]
        w0 = len(anchor.seq)
        members = [c for c in remaining if abs(len(c.seq) - w0) <= tol]
        out.append(([c.seq for c in members], w0))
        remaining = [c for c in remaining if abs(len(c.seq) - w0) > tol]
    return out


def seed_seqs(cands, *, tol=3):
    """Width-coherent seed sequences for building the genome-rescan PWM. With the generous width search the
    candidates span many widths AND mix generator types; `counts_from_seqs` uses the MODAL width, so a mixed
    pool can collapse the seed PWM onto a short, non-specific (AT-rich) motif that rescans to noise. This
    anchors the seed on REAL palindromes/repeats: it ranks IR/DR candidates (the de-novo dyads -- NOT the
    family/PWM-seeded single hits, whose -log10(p) scores are AT-biased and on a different scale), prefers
    conserved ones, takes the top candidate's width, and keeps every candidate within `tol` bp of it.
    Returns sequences (best first). Selects the SEED only; every genome rescan hit is still kept downstream.

    THE ANCHOR IS THE POINT. Which candidate wins does not merely set the ranking, it sets `w0` and
    therefore which candidates survive the +/- `tol` bp filter -- so a wrong winner does not rank the right
    answer badly, it deletes it from the seed set. Measured on TtgR (finding 3 of the TetR-137
    benchmark, recorded in the analysis repo before commit 5815559): the true 28 bp operator was
    generated byte-identical at rank 3 of 8 and evicted by a 23 bp anchor set
    by a GC-rich decoy, after which the seed PWM was built from three decoys, the rescan returned 110 hits
    none on the site, and all five exported AF3 jobs pointed at the wrong DNA. The dyad score behind that
    ordering has AUC 0.250 (0.557 length-normalised, i.e. no signal) against 203 known operators.

    **`pipeline.run_novel` no longer uses this** -- it uses `seed_clusters`, which keeps every width
    instead of choosing between them, because re-ranking the anchor was measured only to move who gets
    evicted (it rescues TtgR and ModE and loses CueR). This remains as the single-best seed set for
    callers that genuinely want one, and returns exactly `seed_clusters(...)[0]`.

    Loosening `tol` is NOT the alternative: the tolerance exists because a mixed-width pool collapses
    onto a short AT-rich motif that rescans to noise (`counts_from_seqs` keeps only the modal width, so
    a mixed pool silently discards most of itself)."""
    clusters = seed_clusters(cands, tol=tol, max_clusters=1)
    return clusters[0][0] if clusters else []


# --------------------------------------------------------------------------- self-test
def _demo():
    import random
    random.seed(7)

    half = "TTGACATG"                          # 8-bp half-site (distinctive)
    op = half + "".join(random.choice("ACGT") for _ in range(19)) + _rc(half)   # MerR-style 19-bp dyad
    center0 = 100 + len(op) / 2.0
    assert len(op) == 35

    def bg(n):
        return "".join(random.choice("ACGT") for _ in range(n))

    region = bg(100) + op + bg(120)            # operator planted at 100..135
    cands = predict(region, family="MerR")
    top = cands[0]
    print(f"top de-novo candidate: center={top.center:.1f} kind={top.kind} "
          f"half={top.half} spacer={top.spacer} score={top.score:.1f}  {top.seq}")
    assert abs(top.center - center0) <= 4, "did not center on the planted dyad"
    assert top.kind == "IR" and top.spacer == 19, "should detect the 19-bp inverted-repeat dyad"

    # conservation: same operator across homolog regions should boost the right candidate
    homologs = []
    for _ in range(6):
        o = list(op)
        o[random.randrange(6, 25)] = random.choice("ACGT")     # mutate the spacer, keep half-sites
        homologs.append(bg(40) + "".join(o) + bg(40))
    cands2 = predict(region, family="MerR", homolog_regions=homologs)
    t2 = cands2[0]
    print(f"with conservation: top center={t2.center:.1f} conservation={t2.conservation:.2f} "
          f"final={t2.final:.1f}")
    assert abs(t2.center - center0) <= 4 and t2.conservation >= 0.8, "conservation boost failed"

    # family-agnostic call still finds it (default spacer prior spans 19)
    c3 = predict(region)
    assert any(abs(c.center - center0) <= 4 for c in c3[:5]), "default-prior finder missed the dyad"
    print("OK: de-novo IR dyad recovered (family + default), conservation boosts the conserved site.")


if __name__ == "__main__":
    _demo()
