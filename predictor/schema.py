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
    #: Empty defaults preserve legacy readers; new producers state identity and class separately.
    identity_status: str = ""
    sensor_class: str | None = None
    class_status: str = ""
    class_evidence: dict = field(default_factory=dict)

    def __post_init__(self):
        self.candidates = tuple(self.candidates or ())
        # Enforce the documented [0,1] range at the type. Ligify derives its confidence from a rank that
        # can be NEGATIVE, which put out-of-range values in this field and made them awkward to reason
        # about downstream. Clamping does NOT change who votes: `from_calls` enfranchises
        # `confidence > 0`, so a clamped 0.0 abstains exactly as a negative did -- which is the intended
        # meaning ("this source has no usable evidence") and the guard that stops a random operon
        # metabolite becoming a metalloregulator's inducer.
        if self.confidence is not None:
            self.confidence = min(1.0, max(0.0, float(self.confidence)))


#: EXP072 source precedence for the HEADLINE: clade (SSN) > published database > local operon chemistry.
#: A clade label rests on curated anchors; a published Ligify record is a reviewed tool's prediction for
#: this protein (or, at >= 99 % identity, its near-identical homologue); the local route is our own
#: operon-context guess, whose confidence is a Ligify rank and not a probability. Confidence therefore
#: only breaks ties WITHIN a tier -- a local tie at 0.77 must not outrank a clade call at 0.50 (the
#: E. coli FrmR case). On the MG1655 panel the marginal rates point the same way (truth in the source's
#: candidate set: SSN 13/13, published 4/11, local 4/17 -- EXP072 prereg section 5.1), but the panel
#: holds a single head-to-head disagreement, so the order is the owner's stated policy and not a fitted
#: result. A source not listed sits with the local route.
_SOURCE_TIER = {"ssn_cluster": 0, "coordination": 0, "ligify_db": 1, "ligify_db_homolog": 1}
_LOWEST_TIER = 2


def source_tier(source: str) -> int:
    """Headline precedence of a source; smaller wins."""
    return _SOURCE_TIER.get(source or "", _LOWEST_TIER)


def precedence_key(call) -> tuple:
    """Sort key for choosing among exact-eligible calls: source tier first, then confidence."""
    return (source_tier(call.source), -(call.confidence or 0.0))


@dataclass
class InducerConsensus:
    """Unified inducer inference that RETAINS every raw source's call, including rejected proposals.
    Final policy determines which evidence and alternatives remain accepted. `top` is the headline;
    `coordination_gate` records whether the metal-coordination
    signature was confirmed (None = not checked) -- enforcement of the gate is the caller's policy."""
    top: str | None = None
    agreement: float = 0.0                 # fraction of positive calling sources giving exact support
    coordination_gate: bool | None = None
    calls: list = field(default_factory=list)      # list[InducerCall], all sources retained
    notes: tuple = ()
    #: The fused candidate set when any source named a mixture. `top` is still the single headline, so
    #: this is additive: a reader who wants one answer keeps getting one.
    candidates: tuple = ()
    #: "" | "co_sensed" | "unresolved" -- see InducerCall.mixture_kind.
    mixture_kind: str = ""
    #: Human-facing label; use sensor_class for new classification, never this display string.
    top_display: str = ""
    #: Finalization is explicit so a serialized null/empty set cannot revive a legacy proposal.
    finalized: bool = False
    #: Set membership is compatibility, not exact identification or independent corroboration.
    compatibility: float = 0.0
    support_sources: tuple = ()
    compatible_sources: tuple = ()
    #: Policy-accepted evidence; `calls` above always retains the original source output.
    accepted_calls: list = field(default_factory=list)
    #: Portable record acquired during the original local chemistry inference, never re-fetched.
    chemistry_snapshot: dict | None = None
    #: Candidate-only acceptances by existing finalization policy (not additional confidence votes).
    policy_candidates: tuple = ()
    policy_mixture_kind: str = ""
    identity_status: str = ""
    sensor_class: str | None = None
    class_status: str = ""
    class_evidence: dict = field(default_factory=dict)

    @staticmethod
    def from_calls(calls, *, coordination_gate: bool | None = None) -> "InducerConsensus":
        """Default fusion: headline = the named call of the highest-precedence SOURCE (`source_tier`),
        highest confidence within it (EXP072; before that, highest confidence alone); agreement = exact
        support. A lower-tier call is not discarded -- it still feeds `candidates` and `compatibility`.
        Positive candidate-only sets may contribute alternatives and compatibility, never a named
        headline. Abstentions and nonpositive source proposals are retained without operative support.

        A named call with NON-POSITIVE confidence does not vote either: such a source is asserting "I have
        no usable evidence", and letting it become the headline merely because it was the only one to speak
        produces a confident-looking wrong answer. That was a real failure mode -- for an UNLABELLED SSN
        cluster the only naming source is Ligify, which can report a negative rank-derived confidence
        (e.g. -0.2), so a random operon metabolite (D-threo-isocitrate) became the inducer of a
        metalloregulator at agreement 1.0. With no positive voter we return top=None, which routes to the
        honest class-level fallback in `effector.inducer._finalize` ("non-metal effector (undetermined)" /
        "divalent metal (ion unresolved)") instead of inventing a ligand.
        """
        from predictor.effector.chemical_identity import exact_eligible
        calls = list(calls)
        voting = [c for c in calls if exact_eligible(c)]
        cons = InducerConsensus(top=min(voting, key=precedence_key).ligand if voting else None,
                                calls=calls, coordination_gate=coordination_gate)
        cons.refresh([c for c in calls if (c.confidence or 0) > 0])
        return cons

    def refresh(self, accepted_calls, *, finalized: bool = False) -> "InducerConsensus":
        """Rebuild all headline-dependent fields after policy, retaining the raw source calls.

        Agreement uses positive calling sources as its denominator, including rejected dissent.
        An unresolved source is only compatible with a member, even if its raw first proposal happens
        to equal the headline. Effective coordination proposals can be recorded separately without
        overwriting an original abstention. New identity keys preserve charge and transport roles.
        """
        from predictor.effector.chemical_identity import (
            identity_key, chemistry, display_label, exact_eligible, unique_names, SENSOR_CLASSES,
        )

        positive = [c for c in accepted_calls if (c.confidence or 0) > 0]
        self.accepted_calls = positive

        def key(x):
            return identity_key(x) or x

        top_key = key(self.top) if self.top else None
        exact, compatible = [], []
        cands, seen = [], set()
        for c in positive:
            components = tuple(c.candidates or ()) or ((c.ligand,) if c.ligand else ())
            keys = {key(x) for x in components if x}
            matches = bool(top_key and c.ligand and key(c.ligand) == top_key)
            if top_key and (matches or top_key in keys):
                if not exact_eligible(c):
                    compatible.append(c.source)
                elif matches or c.mixture_kind == "co_sensed":
                    exact.append(c.source)
                else:
                    compatible.append(c.source)
            for x in components:
                if x and key(x) not in seen:
                    seen.add(key(x))
                    cands.append(x)
        for x in self.policy_candidates:
            if x and key(x) not in seen:
                seen.add(key(x))
                cands.append(x)
        # Put a named final headline first; class placeholders do not manufacture candidates.
        if top_key in seen:
            cands = [self.top] + [x for x in cands if key(x) != top_key]
        denominator = {c.source for c in self.calls
                       if (c.confidence or 0) > 0 and (c.ligand or c.candidates)}
        denominator.update(c.source for c in positive if c.ligand or c.candidates)
        self.support_sources = tuple(dict.fromkeys(exact))
        self.compatible_sources = tuple(dict.fromkeys(compatible))
        self.agreement = len(self.support_sources) / len(denominator) if denominator else 0.0
        self.compatibility = len(self.compatible_sources) / len(denominator) if denominator else 0.0
        # A singleton represented by the named headline needs no separate set. A singleton that
        # differs from a class/null/other headline is retained evidence, without implying a mixture.
        self.candidates = tuple(cands) if len(cands) > 1 or (cands and key(cands[0]) != top_key) else ()
        kinds = {c.mixture_kind for c in positive if c.mixture_kind and c.candidates}
        if self.policy_candidates and self.policy_mixture_kind:
            kinds.add(self.policy_mixture_kind)
        self.mixture_kind = ("unresolved" if "unresolved" in kinds else
                             "co_sensed" if kinds else "") if len(cands) > 1 else ""
        self.identity_status = ("unresolved" if not self.top and len(cands) > 1 else
                                chemistry(cands[0])[1] if not self.top and cands else chemistry(self.top)[1])
        if not self.top and cands and self.identity_status == "resolved":
            self.identity_status = "unresolved"
        if not self.top and len(cands) == 1:
            declared = {c.identity_status for c in positive if c.identity_status
                        and any(key(x) == key(cands[0]) for x in (c.candidates or (c.ligand,)))}
            if len(declared) == 1 and not declared.intersection({"resolved", "co_sensed"}):
                self.identity_status = next(iter(declared))
        if self.top and chemistry(self.top)[1] == "class_only" and cands:
            self.identity_status = "unresolved"
        if any(c.ligand == self.top and c.mixture_kind == "co_sensed" and exact_eligible(c)
               and len(unique_names(c.candidates)) > 1 for c in positive):
            self.identity_status = "co_sensed"
        self.top_display = display_label(self.top, self.candidates)
        # Class support belongs to the sensed signal, independently of exact identity.
        class_records = []
        uncertain = []
        policy_input = self.class_evidence.get("policy_input")
        if self.class_evidence.get("policy") and self.sensor_class in SENSOR_CLASSES:
            class_records.append({**self.class_evidence["policy"], "class": self.sensor_class})
        elif self.class_evidence.get("policy") and self.sensor_class is not None:
            policy_input = self.sensor_class
            self.sensor_class = None
        if policy_input is not None:
            uncertain.append({"source": "policy", "reported_class": policy_input, "status": "unknown"})
        for c in positive:
            cls = c.sensor_class
            basis = c.class_evidence
            invalid_class = cls is not None and cls not in SENSOR_CLASSES
            if invalid_class or c.class_status not in {"", "supported"}:
                uncertain.append({**basis, "source": c.source, "reported_class": cls,
                                  "status": "conflicting" if c.class_status == "conflicting" else "unknown"})
                continue
            if not cls and c.source not in {"operon_chem", "ligify", "ligify_db", "ligify_db_homolog"}:
                cls = {"metal": "metal", "redox": "redox", "non-metal": "organic",
                       "multidrug": "organic"}.get(c.role)
                basis = {"basis": "explicit_source_role"} if cls else {}
            if c.source not in {"operon_chem", "ligify"}:
                from predictor.effector.chemical_identity import collection_class
                names = tuple(c.candidates or ()) or ((c.ligand,) if c.ligand else ())
                collection_cls, collection_status, records = collection_class(names)
                if collection_status == "conflicting":
                    for candidate_cls in dict.fromkeys(record[0] for record in records if record[0]):
                        class_records.append({"source": c.source, "class": candidate_cls,
                                              "basis": "conflicting_signal_candidates"})
                if not cls:
                    cls = collection_cls
                    basis = {"basis": "explicit_signal_vocabulary", "chemistry": records} if cls else {}
            if cls:
                class_records.append({**basis, "source": c.source, "class": cls})
        classes = {r["class"] for r in class_records}
        policy = self.class_evidence.get("policy") if self.sensor_class in SENSOR_CLASSES else None
        if not policy:
            self.sensor_class = next(iter(classes)) if len(classes) == 1 else None
        self.class_status = ("conflicting" if len(classes) > 1 or any(r["status"] == "conflicting" for r in uncertain)
                             else "unknown" if uncertain else "supported" if classes else "unknown")
        self.class_evidence = {"records": class_records, **({"policy": policy} if policy else {}),
                               **({"source_uncertainty": uncertain} if uncertain else {}),
                               **({"policy_input": policy_input} if policy_input is not None else {})}
        self.finalized = finalized
        return self

    @property
    def metal(self) -> bool:
        if self.class_status:
            return self.sensor_class == "metal"
        return bool(self.top) and any(c.ligand == self.top and c.role == "metal" and (c.confidence or 0) > 0
                   for c in (self.accepted_calls if self.finalized else self.calls))


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
                                         InducerCall("ligify", "citrate", 0.9, "non-metal")])
    # EXP072: the clade outranks a louder local operon-chemistry call; the local name stays a candidate
    assert split.top == "Ni2+" and abs(split.agreement - 0.5) < 1e-9, "disagreement -> 0.5 agreement"
    assert "citrate" in split.candidates, "a lower-precedence source is kept, not discarded"
    same_tier = InducerConsensus.from_calls([InducerCall("ligify", "citrate", 0.9, "non-metal"),
                                             InducerCall("operon_chem", "salicylate", 0.4, "non-metal")])
    assert same_tier.top == "citrate", "confidence decides inside one tier"

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
