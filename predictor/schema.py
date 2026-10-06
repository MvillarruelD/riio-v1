"""
schema.py -- the CANONICAL coordinate frame + operator types every signal is lifted into (Phase 0.1).

The four-plus evidence types emit coordinates in *different* frames:

  * motif_finder  -- offset within the extracted inter-operon region (region-relative, + strand),
  * motif_rescan  -- already genome coordinates (GenomeHit), but its own dataclass,
  * regulon       -- the operator site that drives an operon (genome coords, RegulonOperon),
  * literature    -- curated genome coordinates (the oracle ceiling),
  * structure/DeepPBS (later) -- offset over the input complex's DNA.

Nothing is comparable until they share ONE frame. This module defines it:

    (genome_accession, strand, dyad_center, half_site_spans, offset_to_tss)

with the operator **dyad center** as the single zero. `signals/coords.py` provides a
`to_genome_coords()` adapter per method that returns a `CanonicalOperator` in this frame, plus a
Tomtom-style `align_pwms()` (motif PCC is undefined without an offset/orientation search). Fusion,
voting and the benchmark all consume `CanonicalOperator` / `OperatorPrediction`, never a method's
native type.

Run `python schema.py` for a self-test (construction, half-site/spacer invariants, TSS offset sign).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

# confidence / provenance vocabulary
CURATED = "CURATED"                       # experimentally backed (co-crystal / wet-lab); truth-eligible
DERIVED = "DERIVED"                        # transferred from a curated neighbour by homology
PREDICTED_UNVERIFIED = "PREDICTED_UNVERIFIED"   # de-novo, not yet structurally/experimentally checked
ORACLE = "ORACLE"                          # read from the answer key; excluded from fusion/voting
_CONF_KINDS = {CURATED, DERIVED, PREDICTED_UNVERIFIED, ORACLE}


@dataclass(frozen=True)
class Provenance:
    """Where a prediction came from and how much to trust it."""
    method: str                            # 'motif_finder' | 'motif_rescan' | 'regulon' | ...
    confidence_kind: str = PREDICTED_UNVERIFIED
    oracle: bool = False                   # literature=True -> never fused/voted
    notes: tuple = ()

    def __post_init__(self):
        if self.confidence_kind not in _CONF_KINDS:
            raise ValueError(f"confidence_kind must be one of {_CONF_KINDS}, got {self.confidence_kind!r}")


@dataclass
class HalfSite:
    """One half of the dyad, in canonical genome coordinates (0-based, half-open, + strand)."""
    start: int
    end: int
    strand: str = "+"

    @property
    def width(self) -> int:
        return self.end - self.start


@dataclass
class CanonicalOperator:
    """An operator candidate in the canonical genome frame. The dyad center is the zero.

    All integer coordinates are genome 0-based half-open on the + strand. `strand` records which
    strand the site was *read* on (palindromic operators are symmetric; it matters for direct repeats
    and for placing the TSS). Scores are filled progressively -- `None` means "not yet evaluated".
    """
    genome_accession: str
    dyad_center: int                       # canonical zero (genome coord)
    start: int                             # full operator span, genome coords
    end: int
    strand: str = "+"
    half_sites: tuple | None = None        # (HalfSite, HalfSite) or None
    spacer_len: int | None = None          # MerR ~19 is the mechanistic signature
    seq: str = ""
    pwm: np.ndarray | None = None          # 4xL (A,C,G,T), probability matrix
    per_base_importance: np.ndarray | None = None   # length-L (conservation or contact weight)
    offset_to_tss: int | None = None       # signed: dyad_center -> TSS (downstream +, upstream -)
    # progressive scores (see WORKPLAN Appendix A.1)
    sequence_score: float | None = None    # IR/DR quality x conservation (structure-free)
    structure_score: float | None = None   # operator-likeness from a fold; None until verified
    combined_score: float | None = None
    score: float | None = None             # the method's own headline score; the genome rescan ranks operators by it
    pvalue: float | None = None             # scan-level null probability, when the method supplies one
    qvalue: float | None = None             # multiple-testing-adjusted value, when supplied
    qvalue_kind: str | None = None          # e.g. exact_bh | conservative_bh_upper_bound
    confidence: float | None = None        # calibrated [0,1] when available
    generator: str = ""                    # which candidate generator emitted it (IR/DR/family/...)
    provenance: Provenance | None = None
    # SSN-cluster + template-retrieval provenance (populated by the homolog-selection / retrieval stages)
    ssn_cluster: str | None = None         # SSN isofunctional cluster id, e.g. "MerR_9" (None if unassigned)
    homolog_selection: str = ""            # how the homolog set was chosen: "ssn" | "blast" | ""
    #: Structural-template retrieval provenance. No production route fills this: the retrievers it was
    #: built for (Foldseek, Folddisco) were measured and removed. Kept as an additive field.
    template_provenance: list = field(default_factory=list)  # [{retriever, target, iptm|tm, rmsd}]

    @property
    def width(self) -> int:
        return self.end - self.start

    def __post_init__(self):
        if self.end < self.start:
            raise ValueError(f"operator end {self.end} < start {self.start}")
        if not (self.start <= self.dyad_center <= self.end):
            raise ValueError(f"dyad_center {self.dyad_center} outside span [{self.start}, {self.end}]")
        if self.pwm is not None:
            p = np.asarray(self.pwm, dtype=float)
            if p.ndim != 2 or 4 not in p.shape:
                raise ValueError(f"pwm must be 4xL, got shape {p.shape}")


@dataclass
class OperatorPrediction:
    """A method's ranked output, all candidates already lifted to the canonical frame."""
    candidates: list                       # list[CanonicalOperator], best first
    method: str
    genome_accession: str | None = None
    failure_mode: str | None = None        # Snowprint taxonomy (WORKPLAN A.7)
    provenance: Provenance | None = None
    notes: list = field(default_factory=list)

    @property
    def best(self):
        return self.candidates[0] if self.candidates else None

    @property
    def oracle(self) -> bool:
        return bool(self.provenance and self.provenance.oracle)


@dataclass
class InducerCall:
    """One source's opinion about the cognate inducer, kept side-by-side with the others.

    The sources that emit one in production are the SSN cluster, Ligify (the vendored database and the
    live operon-chemistry route), MetalNet2 and the metal-coordination signature. `source` is not
    validated against that list: an abstaining source is still recorded, and an experiment may add one.
    """
    source: str                            # 'ssn_cluster'|'ligify_db'|'ligify'|'metalnet'|'coordination'
    ligand: str | None                     # 'Zn2+' | 'L-arabinose' | None (source made no call)
    confidence: float = 0.0                # [0,1] -- clamped in __post_init__
    role: str = ""                         # 'metal' | 'non-metal' | 'redox' | 'multidrug' | ''
    evidence: dict = field(default_factory=dict)   # source-specific audit (cluster_id, enzyme, rhea, ...)
    #: Ordered components when this source names MORE THAN ONE species -- e.g. an SSN cluster whose
    #: anchors co-sense `Ni2+/Co2+`. Empty means "not a mixture", never "unknown". `ligand` keeps its
    #: meaning (the single best candidate) so every existing consumer is unaffected.
    candidates: tuple = ()
    #: "" | "co_sensed" | "unresolved". `unresolved` means the evidence genuinely does not separate the
    #: candidates -- not that the source failed. Reporting one of them as the answer would be a
    #: fabricated precision; reporting both, and why nothing separates them, is the result.
    mixture_kind: str = ""

    def __post_init__(self):
        # Enforce the documented [0,1] range at the type. Ligify derives its confidence from a rank that
        # can be NEGATIVE, which put out-of-range values in this field and made them awkward to reason
        # about downstream. Clamping does NOT change who votes: `from_calls` enfranchises
        # `confidence > 0`, so a clamped 0.0 abstains exactly as a negative did -- which is the intended
        # meaning ("this source has no usable evidence") and the guard that stops a random operon
        # metabolite becoming a metalloregulator's inducer.
        if self.confidence is not None:
            self.confidence = min(1.0, max(0.0, float(self.confidence)))


@dataclass
class InducerConsensus:
    """Unified inducer inference that RETAINS every source's call. Per the design: agreement across
    sources strengthens the call, but disagreement leaves all candidates valid (so `calls` is never
    discarded). `top` is the headline pick; `coordination_gate` records whether the metal-coordination
    signature was confirmed (None = not checked) -- enforcement of the gate is the caller's policy."""
    top: str | None = None
    agreement: float = 0.0                 # fraction of *calling* sources that name `top`
    coordination_gate: bool | None = None
    calls: list = field(default_factory=list)      # list[InducerCall], all sources retained
    notes: tuple = ()
    #: The fused candidate set when any source named a mixture. `top` is still the single headline, so
    #: this is additive: a reader who wants one answer keeps getting one.
    candidates: tuple = ()
    #: "" | "co_sensed" | "unresolved" -- see InducerCall.mixture_kind.
    mixture_kind: str = ""
    #: Human-facing label set by effector.inducer._finalize. It normally equals `top`; for unresolved
    #: class placeholders it can expose the named candidate shortlist. DISPLAY ONLY: `top` remains the
    #: sole input to inducer_class_of and to every derived class/count.
    top_display: str = ""

    @staticmethod
    def from_calls(calls, *, coordination_gate: bool | None = None) -> "InducerConsensus":
        """Default fusion: headline = highest-confidence named call; agreement = share of *voting*
        sources that name it. Sources that abstained (ligand=None) are kept but don't vote.

        A named call with NON-POSITIVE confidence does not vote either: such a source is asserting "I have
        no usable evidence", and letting it become the headline merely because it was the only one to speak
        produces a confident-looking wrong answer. That was a real failure mode -- for an UNLABELLED SSN
        cluster the only naming source is Ligify, which can report a negative rank-derived confidence
        (e.g. -0.2), so a random operon metabolite (D-threo-isocitrate) became the inducer of a
        metalloregulator at agreement 1.0. With no positive voter we return top=None, which routes to the
        honest class-level fallback in `effector.inducer._finalize` ("non-metal effector (undetermined)" /
        "divalent metal (ion unresolved)") instead of inventing a ligand.
        """
        named = [c for c in calls if c.ligand]
        voting = [c for c in named if c.confidence > 0]
        if not voting:
            return InducerConsensus(calls=list(calls), coordination_gate=coordination_gate)
        top = max(voting, key=lambda c: c.confidence).ligand
        # A source that names a MIXTURE containing the headline agrees with it. Comparing raw strings
        # counted `Ni2+/Co2+` as disagreeing with `Ni2+`, which is the same defect class as counting
        # `Ni` and `Ni2+` as two inducers -- a flattened label can never equal one of its own members.
        def _backs(c) -> bool:
            return c.ligand == top or (bool(c.candidates) and top in c.candidates)

        agreement = sum(1 for c in voting if _backs(c)) / len(voting)
        # Union of the candidates, headline first, order-stable. Only sources that actually VOTE
        # contribute, plus any explicit mixture components. A source asserting "I have no usable
        # evidence" (confidence <= 0) must not put a ligand into the candidate set: on IdeR that swept
        # in a down-weighted Ligify guess of `nitrate` alongside the two real alternatives, which
        # makes a considered tie look like a shrug.
        #
        # Deduplicate on the CANONICAL token rather than the raw string. Sources spell one species
        # differently -- the SSN cluster table writes arsenite `As3+` where the anchor KB writes
        # `As(III)` -- so a raw dedup lets the same ion appear twice. Seen in the v7 run: E. coli ArsR
        # came out ('As3+', 'antimonite', 'As(III)'), reading as a three-way ambiguity when it is a
        # two-way one. The RAW spelling is what gets kept, since `top` is raw and the report prints
        # these beside it; only the "have I got this one already?" identity is canonicalised.
        from predictor.effector import inducer_vocab as _vocab

        cands: list = [top]
        seen = {_vocab.normalise(top) or top}
        for c in calls:
            explicit = tuple(c.candidates or ())
            implicit = (c.ligand,) if (c.ligand and c.confidence > 0) else ()
            for x in explicit or implicit:
                if not x:
                    continue
                key = _vocab.normalise(x) or x
                if key not in seen:
                    seen.add(key)
                    cands.append(x)
        kinds = {c.mixture_kind for c in calls if c.mixture_kind}
        return InducerConsensus(top=top, agreement=agreement, coordination_gate=coordination_gate,
                                calls=list(calls),
                                candidates=tuple(cands) if len(cands) > 1 else (),
                                mixture_kind=("unresolved" if "unresolved" in kinds
                                              else ("co_sensed" if kinds else "")))

    @property
    def metal(self) -> bool:
        return any(c.ligand == self.top and c.role == "metal" for c in self.calls)


# --------------------------------------------------------------------------- self-test
def _demo() -> None:
    pwm = np.full((4, 10), 0.25)
    op = CanonicalOperator(
        genome_accession="NC_000964.3", dyad_center=1050, start=1035, end=1070, strand="+",
        half_sites=(HalfSite(1035, 1043, "+"), HalfSite(1062, 1070, "+")),
        spacer_len=19, seq="A" * 35, pwm=pwm, offset_to_tss=-12,
        provenance=Provenance("motif_finder", PREDICTED_UNVERIFIED))
    print(f"operator span={op.start}-{op.end} (w={op.width}) dyad={op.dyad_center} "
          f"spacer={op.spacer_len} tss_offset={op.offset_to_tss} prov={op.provenance.method}")
    assert op.width == 35 and op.half_sites[0].width == 8
    assert op.offset_to_tss == -12, "TSS upstream of dyad should be negative"

    pred = OperatorPrediction([op], method="motif_finder", genome_accession="NC_000964.3",
                              provenance=Provenance("motif_finder"))
    assert pred.best is op and not pred.oracle

    lit = OperatorPrediction([], method="literature",
                             provenance=Provenance("literature", CURATED, oracle=True))
    assert lit.oracle and lit.best is None, "literature must be flagged oracle"

    # SSN/template provenance fields are additive (default empty) and carry through construction
    op2 = CanonicalOperator("NC_x", dyad_center=50, start=40, end=60, ssn_cluster="MerR_9",
                            homolog_selection="ssn",
                            template_provenance=[{"retriever": "example", "target": "1abc", "tm": 0.7}])
    assert op2.ssn_cluster == "MerR_9" and op2.homolog_selection == "ssn"
    assert op2.template_provenance[0]["retriever"] == "example"

    # unified inducer consensus: agreement strengthens, abstainers retained, disagreement keeps all valid
    calls = [InducerCall("ssn_cluster", "Zn2+", 0.8, "metal", {"cluster_id": "MerR_9"}),
             InducerCall("coordination", "Zn2+", 0.6, "metal"),
             InducerCall("ligify", None, 0.0)]                  # abstains
    cons = InducerConsensus.from_calls(calls, coordination_gate=True)
    assert cons.top == "Zn2+" and abs(cons.agreement - 1.0) < 1e-9 and cons.metal
    assert len(cons.calls) == 3 and cons.coordination_gate is True, "all sources retained + gate recorded"
    split = InducerConsensus.from_calls([InducerCall("ssn_cluster", "Ni2+", 0.7, "metal"),
                                         InducerCall("ligify", "L-lactate", 0.9, "non-metal")])
    assert split.top == "L-lactate" and abs(split.agreement - 0.5) < 1e-9, "disagreement -> 0.5 agreement"

    # invariants reject malformed inputs
    for bad in (
        lambda: CanonicalOperator("X", dyad_center=5, start=10, end=20),     # dyad outside span
        lambda: CanonicalOperator("X", dyad_center=15, start=20, end=10),    # end < start
        lambda: Provenance("m", confidence_kind="MAYBE"),                    # bad confidence kind
    ):
        try:
            bad()
        except ValueError:
            pass
        else:
            raise AssertionError("invariant not enforced")

    print("OK: canonical operator/prediction types + provenance vocabulary + invariants verified.")


if __name__ == "__main__":
    _demo()
