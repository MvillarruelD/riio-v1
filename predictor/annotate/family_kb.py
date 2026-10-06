"""
family_kb.py -- read the vendored per-cluster ANCHOR SETS.

An SSN cluster used to carry one hand-picked `representative` whose inducer became the cluster's
label. That is why handoff §5.7 has to warn that 7 of the 17 benchmark regulators *are* the anchor of
the cluster they get assigned to, and why a label carried no notion of agreement, conflict or
support: with one anchor there is nothing to agree with.

An anchor set is every protein in the cluster with independently known biology, built in
`analysis/family_kb` and exported to `predictor/data/refs/family_kb/<TAG>.json`. Each anchor keeps its
accession, organism, inducer, evidence level and `source_route`, so a caller can weigh agreement and
a benchmark can EXCLUDE anchors by provenance -- which the expanded set makes mandatory rather than
optional, since more anchors means more benchmark TFs become anchored.

Two routes contribute, and they are kept distinguishable on purpose:
  `kb`      a census protein carrying sourced biology (RegulonDB, PRODORIC, DBTBS, curated, ...)
  `figure`  an anchor caption the survey's authors drew on their own network figure, resolved to an
            accession and placed by exact sequence -- independent of our curation, so agreement
            between the two is corroboration rather than restatement.

This module is READ-ONLY and never rebuilds. Regenerating the vendored files is a research-workspace
job (`analysis/family_kb/vendor_family_kb.py`), which carries the empty-rebuild guard. If the files
are absent the module degrades to "no anchors known" rather than failing, exactly as the SSN layer
does when a family has no cluster DB.

Run `python -m predictor.annotate.family_kb` for a self-test.
"""
from __future__ import annotations

import functools
import json
from dataclasses import dataclass, field
from pathlib import Path

from predictor import resources

_REPO = Path(__file__).resolve().parents[2]
KB_DIR = resources.REFS_DIR / "family_kb"

#: label-confidence tiers, matching the vocabulary used throughout this project's reports
TIER_HIGH = "high"
TIER_TENTATIVE = "tentative"
TIER_MIXED = "mixed-anchor"
TIER_UNKNOWN = "unknown"
#: written by tools/apply_site_typing.py (2026-09-17) for a clade whose curated anchors all describe
#: the STRUCTURAL metal site (e.g. Fur's Zn(Cys)4 cage), so none states the sensed metal and the
#: effector comes from the clade label alone
TIER_STRUCTURAL_ONLY = "structural-only"
TIERS = (TIER_HIGH, TIER_TENTATIVE, TIER_MIXED, TIER_UNKNOWN, TIER_STRUCTURAL_ONLY)


@dataclass(frozen=True)
class Anchor:
    accession: str
    #: hash of the protein sequence. Load-bearing for holdouts: 64 census rows carry E. coli CueR's
    #: exact sequence under different accessions, so dropping one accession is not a holdout.
    seq_md5: str | None = None
    #: EVERY accession carrying this exact sequence, `accession` first. The vendored KB stores one
    #: record per distinct sequence -- 2,857 rows collapse to 164 sequences, and MerR alone holds 745
    #: rows for 23 proteins -- so the other accessions have nowhere else to live. Every
    #: accession-keyed lookup below must search this, not `accession` alone, or holding out a
    #: non-representative accession silently drops nothing and the holdout is theatre.
    accessions: tuple[str, ...] = ()
    name: str | None = None
    organism: str | None = None
    inducer: str | None = None
    evidence_level: str | None = None
    #: How the anchor's biology was established: "literature" (a curator read a paper), "omics"
    #: (a high-throughput record such as RegulonDB ChIP-seq), or "computational".
    #:
    #: This was `evidence_tier = 1 | 2 | 3`. It is words now because the project used the word "tier"
    #: for four unrelated things -- this, a clade's label support, an operator's locality rank, and
    #: the tier-1/tier-2 contract that decides what the paper may state as a result. A bare "tier 1"
    #: in a table could mean any of them, and the numeric spelling collided with the one that
    #: matters most to a reader.
    #:
    #: NOT interchangeable with `evidence_level` (experimental / curated / inferred): a RegulonDB
    #: ChIP-seq record is `experimental` at level but "omics" here, so reading the level as a proxy
    #: silently admits omics anchors -- precisely what a literature-only inducer rule exists to
    #: exclude. None means the vendored blob predates this field, never "unknown by nature".
    evidence_grade: str | None = None
    source_route: str = ""
    n_operators: int = 0
    n_regulon_genes: int = 0

    def has_accession(self, acc_upper: str) -> bool:
        """Does this anchor carry `acc_upper` (already upper-cased) under any of its accessions?"""
        return (self.accession.upper() == acc_upper
                or any(x.upper() == acc_upper for x in self.accessions))

    @property
    def is_experimental(self) -> bool:
        return self.evidence_level == "experimental"

    @property
    def is_literature(self) -> bool:
        """A curator read a paper. The only grade permitted to STATE an inducer."""
        return self.evidence_grade == "literature"


@dataclass(frozen=True)
class AnchorSet:
    """Every anchor of one cluster, plus what they collectively say."""
    cluster_id: str
    n_anchors: int = 0
    #: distinct SEQUENCES, not rows. The census holds 64 identical TrEMBL copies of E. coli CueR;
    #: counting rows would report overwhelming support for a single piece of evidence.
    n_independent_anchors: int = 0
    n_distinct_inducers: int = 0
    inducer_consensus: str | None = None
    conflict: bool = False
    confidence_tier: str = TIER_UNKNOWN
    anchors: tuple[Anchor, ...] = field(default_factory=tuple)

    @property
    def tier(self) -> str:
        """Coarse tier, for callers that do not want to parse the full string."""
        return self.confidence_tier.split(" (")[0] if self.confidence_tier else TIER_UNKNOWN

    @property
    def distinct_inducers(self) -> tuple[str, ...]:
        """Canonical inducers stated by these anchors, spelling merged and abstentions removed.

        `n_distinct_inducers` and `conflict` above come from the vendored blob and count RAW strings.
        This is the corrected view: `NikR_c1` states `Ni` and `Ni2+`, which is one inducer written two
        ways, and two anchors elsewhere carry sentences meaning "none is known" rather than a ligand.
        """
        from predictor.effector import inducer_vocab
        return tuple(inducer_vocab.distinct(a.inducer for a in self.anchors if a.inducer))

    @property
    def real_conflict(self) -> bool:
        """Do these anchors disagree about the inducer, after spelling is normalised?

        Prefer this to `conflict` when the answer changes behaviour. Measured over the vendored KB,
        the two differ for exactly one cluster -- but that cluster is having its confidence scaled
        down for a typographic difference.
        """
        return len(self.distinct_inducers) > 1

    @property
    def agreement(self) -> float:
        """Share of inducer-stating anchors backing the consensus, in [0,1].

        0.0 when nothing states an inducer -- a cluster with anchors but no inducer opinion is not
        the same as a cluster whose anchors agree, and a caller must be able to tell them apart.
        """
        stated = [a for a in self.anchors if a.inducer]
        if not stated or not self.n_distinct_inducers:
            return 0.0
        return 1.0 / self.n_distinct_inducers if self.n_distinct_inducers > 1 else 1.0

    def excluding(self, accessions) -> "AnchorSet":
        """This set without the named accessions -- leave-N-anchors-out, for benchmarks.

        The holdout is always at SEQUENCE level, and cannot be otherwise. Excluding `P0A9G4` while
        63 identical TrEMBL copies of E. coli CueR still anchor `MerR_7` is not a holdout: the label
        survives on exactly the evidence that was meant to be withheld. Since the vendored KB stores
        one record per distinct sequence with all its accessions, naming any one of them removes the
        protein outright, and the leaky variant is no longer representable. There used to be a
        `by_sequence=False` switch for it; it was removed rather than left as a no-op, because a
        parameter that silently stopped meaning anything is worse than no parameter.

        Recomputes the tier from what remains, so holding out an anchor genuinely weakens the label
        instead of leaving a stale confidence behind it.
        """
        drop = {a.upper() for a in accessions}
        kept = tuple(a for a in self.anchors if not any(a.has_accession(d) for d in drop))
        # Canonical spellings, not raw strings: `Ni` and `Ni2+` are one inducer, and an abstention
        # sentence ("none identified (...)") is not one at all. See `effector.inducer_vocab`.
        from predictor.effector import inducer_vocab
        inducers = set(inducer_vocab.distinct(a.inducer for a in kept if a.inducer))
        # independent support is distinct sequences among what remains
        n_ind = len({a.seq_md5 or a.accession for a in kept})
        if not kept:
            tier = TIER_UNKNOWN
        elif len(inducers) > 1:
            tier = f"{TIER_MIXED} ({len(inducers)} inducers)"
        elif n_ind >= 2:
            tier = f"{TIER_HIGH} ({n_ind} consistent anchors)"
        else:
            tier = f"{TIER_TENTATIVE} (single anchor)"
        return AnchorSet(self.cluster_id, len(kept), n_ind, len(inducers),
                         (next(iter(inducers)) if len(inducers) == 1 else None),
                         len(inducers) > 1, tier, kept)


#: the numeric spelling this field used to carry, for blobs written before the rename
_GRADE_FROM_TIER = {1: "literature", 2: "omics", 3: "computational"}


def _grade(anchor: dict) -> str | None:
    """`evidence_grade`, accepting the older numeric `evidence_tier` spelling."""
    g = anchor.get("evidence_grade")
    if g:
        return str(g)
    t = anchor.get("evidence_tier")
    if t in (None, ""):
        return None
    try:
        return _GRADE_FROM_TIER.get(int(t))
    except (TypeError, ValueError):
        return None


@functools.lru_cache(maxsize=1)
def _load() -> dict[str, AnchorSet]:
    """{cluster_id: AnchorSet} across every vendored family. Absent files -> empty, never an error."""
    out: dict[str, AnchorSet] = {}
    if not KB_DIR.is_dir():
        return out
    for path in sorted(KB_DIR.glob("*.json")):
        try:
            blob = json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            continue
        for cid, c in (blob.get("clusters") or {}).items():
            out[cid] = AnchorSet(
                cluster_id=cid,
                n_anchors=int(c.get("n_anchors") or 0),
                n_independent_anchors=int(c.get("n_independent_anchors") or 0),
                n_distinct_inducers=int(c.get("n_distinct_inducers") or 0),
                inducer_consensus=c.get("inducer_consensus"),
                conflict=bool(c.get("conflict")),
                confidence_tier=c.get("confidence_tier") or TIER_UNKNOWN,
                anchors=tuple(Anchor(**{k: a.get(k) for k in
                                        ("accession", "seq_md5", "name", "organism", "inducer",
                                         "evidence_level", "source_route")},
                                     # absent in blobs written before sequence-deduplication, where
                                     # each accession still had its own row
                                     accessions=tuple(a.get("accessions")
                                                      or (a.get("accession"),)),
                                     evidence_grade=_grade(a),
                                     n_operators=int(a.get("n_operators") or 0),
                                     n_regulon_genes=int(a.get("n_regulon_genes") or 0))
                              for a in (c.get("anchors") or []) if a.get("accession")),
            )
    return out


def anchor_set(cluster_id: str) -> AnchorSet | None:
    """The anchor set for `cluster_id`, or None if the cluster has none vendored."""
    return _load().get(cluster_id)


def anchors_for(cluster_id: str) -> tuple[Anchor, ...]:
    s = anchor_set(cluster_id)
    return s.anchors if s else ()


def is_anchor(accession: str, cluster_id: str | None = None) -> bool:
    """Is this accession an anchor -- of `cluster_id`, or of any cluster?

    The question a benchmark must ask before counting a result as independent evidence.
    """
    acc = (accession or "").upper()
    if not acc:
        return False
    sets = [anchor_set(cluster_id)] if cluster_id else list(_load().values())
    return any(s and any(a.has_accession(acc) for a in s.anchors) for s in sets)


def clusters_anchored_by(accession: str) -> list[str]:
    """Every cluster this accession anchors."""
    acc = (accession or "").upper()
    return sorted(cid for cid, s in _load().items()
                  if any(a.has_accession(acc) for a in s.anchors))


def available() -> bool:
    return bool(_load())


def _demo() -> None:
    sets = _load()
    if not sets:
        print("no vendored anchor sets found at "
              f"{KB_DIR.relative_to(_REPO)} -- module degrades to 'no anchors known'")
        return
    n_anchors = sum(len(s.anchors) for s in sets.values())
    print(f"{len(sets)} clusters, {n_anchors} anchors vendored")

    multi = [s for s in sets.values() if s.n_independent_anchors > 1]
    conflict = [s for s in sets.values() if s.conflict]
    print(f"  {len(multi)} clusters with >1 independent anchor; {len(conflict)} with a conflict")
    assert sets, "vendored anchor sets must not be empty"
    assert n_anchors > 100, f"suspiciously few anchors ({n_anchors}) -- check the vendored files"

    # a known multi-anchor cluster, and leave-one-out must weaken it
    ex = next((s for s in sets.values() if s.n_independent_anchors > 1 and s.anchors), None)
    if ex:
        first = ex.anchors[0].accession
        red = ex.excluding([first])
        print(f"  {ex.cluster_id}: {ex.n_anchors} anchors, tier={ex.tier!r}; excluding {first} -> "
              f"{red.n_anchors} remain, tier={red.tier!r}")
        assert red.n_anchors == ex.n_anchors - 1, "holdout must drop exactly the held-out protein"
        assert is_anchor(first, ex.cluster_id), "the excluded accession must be a known anchor"
        assert not any(a.has_accession(first.upper()) for a in red.anchors),             "excluded anchor still present under one of its other accessions"

    # a conflicted cluster must report it rather than average it away
    if conflict:
        c = conflict[0]
        print(f"  conflict example {c.cluster_id}: {c.n_distinct_inducers} inducers, "
              f"tier={c.confidence_tier!r}, agreement={c.agreement:.2f}")
        assert c.agreement < 1.0, "a conflicted cluster must not report full agreement"
    print("OK: anchor sets load, leave-N-out works, conflicts are visible.")


if __name__ == "__main__":
    _demo()
