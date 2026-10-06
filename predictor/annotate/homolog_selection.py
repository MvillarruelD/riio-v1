"""
homolog_selection.py -- SSN-cluster-restricted homolog selection (Phase 1, Step 1).

The homolog-expansion stage today widens the set by sequence similarity (MSA close-orthologs / MMseqs).
In a divergent family that neighbourhood mixes paralogs that sense *different* metals; since the operator
dyad geometry is conserved but the half-site bases differ by inducer, mixing them blurs the motif. This
module restricts the homolog set to the query's **SSN cluster** (isofunctional group), so the operator is
built from operator-coherent homologs.

Two ways a query relates to the SSN:
  * by **sequence** (robust, the default): MMseqs the query against all SSN members and take the cluster of
    its best hits. This works even when the query's own MSA orthologs are NOT SSN members -- which is the
    common case (e.g. NmtR's Mycobacterium orthologs are absent from the AK22-centric ArsR clusters).
  * by **accession** (cheap, when applicable): if the query's orthologs ARE SSN members, filter them by the
    acc->cluster map directly (`restrict_orthologs`).

Cluster assignment reuses `homologs.homolog_sweep` (MMseqs2); members are read via `ssn_clusters`. The
selected member accessions feed the existing UniProt->genome->promoter path (`msa_homologs`) unchanged.

Run `python -m predictor.annotate.homolog_selection --assign <fasta>` to assign a TF to its SSN cluster.
"""
from __future__ import annotations

import argparse
import functools
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from predictor import resources

_REPO = Path(__file__).resolve().parents[2]
from predictor.annotate import homologs, ssn_clusters as ssn          # noqa: E402

_SSN_DB_DIR = resources.SSN_DATABASE
_SEP = "__"                                                   # member header = "<cluster_id>__<acc>"

# The only families with a curated SSN (cluster FASTAs + inducer labels). SSN-cluster inference is RESERVED
# to these -- for any other family there is no cluster DB, and searching one would let a divergent TF match
# a foreign cluster and inherit its inducer (a cross-family bleed). Callers gate on this set.
# It tracks ssn.FAMILIES, so the 2026 drop (all 12 families now have cluster FASTAs + member DBs) is
# searchable here automatically. The guard below still rejects any family with NO cluster DB on disk, which
# otherwise would MMseqs against a missing/foreign DB and inherit a spurious cluster.
_FAMILIES = set(ssn.FAMILIES)


#: Minimum top-hit identity for an assignment to be made at all. MEASURED, not chosen: a two-
#: population sweep (`analysis/family_kb/calibrate_floors.py`, 400 clade members vs 400 family
#: proteins belonging to NO labelled clade) gives
#:
#:      min_identity   members assigned / correct      proteins-with-no-clade wrongly assigned
#:      0.00 (before)      100.0 %   /  99.2 %                     77.5 %
#:      0.50                99.8 %   /  99.2 %                     43.2 %
#:      0.60                99.2 %   /  99.2 %                      8.2 %
#:      0.70                96.0 %   /  99.2 %                      3.8 %
#:
#: 0.60 is the knee: false assignment falls 9.4x for 0.8 points of member recall, and accuracy on
#: the members that still assign is unchanged. This matters because 57 % of the ArsR family and 69 %
#: of MerR sit outside every labelled clade, and each spurious assignment hands that protein the
#: clade's INDUCER.
MIN_IDENTITY = 0.60

#: Minimum size of the cluster an assignment may land in. Default 0 = OFF, and that is a measured
#: result rather than an omission: over the same sweep a size floor never reduced false assignment
#: (a remote query votes toward the LARGEST clade, not a small one) and only cost member recall --
#: 100 % -> 95.5 % at 100 sequences, -> 90.5 % at 500. The parameter is kept because it is the right
#: knob for a caller who wants to ignore tiny clades, but it is not part of the abstention rule.
MIN_CLUSTER_SIZE = 0


@dataclass
class ClusterAssignment:
    cluster_id: str | None                 # winning SSN cluster, or None if unassigned
    support: float                         # winner's share of total hit bit-score [0,1]
    n_hits: int                            # SSN members hit
    top_identity: float                    # best-hit % identity (sanity on assignment quality)
    votes: dict = field(default_factory=dict)   # cluster_id -> summed bits (audit)
    info: "ssn.ClusterInfo | None" = None  # curated inducer/role/gate for the winner
    method: str = "mmseqs"                 # "direct" | "high_identity" | "mmseqs"
    #: WHY there is or is not a cluster. Callers must be able to tell "this protein is not in the
    #: SSN" from "this family has no SSN" from "we were not confident enough" -- previously all
    #: three were an indistinguishable `cluster_id=None`.
    status: str = "assigned"

    @property
    def assigned(self) -> bool:
        return self.cluster_id is not None


#: every value `ClusterAssignment.status` can take
STATUS_ASSIGNED = "assigned"                       # direct / high-identity / vote, above the floors
STATUS_NO_SSN = "no_ssn_for_family"                # the family has no SSN in this build
STATUS_NO_HIT = "unassigned_no_hit"                # nothing in the family DB matched at all
STATUS_LOW_IDENTITY = "unassigned_low_identity"    # best hit below MIN_IDENTITY
STATUS_LOW_SUPPORT = "unassigned_low_support"      # vote too split to name a winner
STATUS_SMALL_CLUSTER = "unassigned_small_cluster"  # winner smaller than min_cluster_size


# --------------------------------------------------------------------------- SSN member DB
def build_ssn_db(family: str, *, force: bool = False) -> Path:
    """Concatenate a family's cluster FASTAs into one MMseqs target file, header `<cluster_id>__<acc>`.

    Raises ValueError for a family with no SSN members. Writing the empty file instead left a 0-byte
    FASTA on disk that every later call happily returned, and MMseqs then failed several layers away with
    "The input files have no entry" -- an error that says nothing about the actual cause (a family such as
    LacI/GalR simply has no SSN). Failing here names the problem where it happens."""
    packaged = _SSN_DB_DIR / f"ssn_members_{ssn.family_tag(family)}.fasta"
    if packaged.exists() and packaged.stat().st_size > 0 and not force:
        return packaged
    # Wheels are immutable. A forced or missing-index rebuild belongs in the user cache.
    out = resources.cache_path("ssn", packaged.name)
    if out.exists() and out.stat().st_size > 0 and not force:
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with out.open("w", encoding="ascii", errors="ignore") as fh:
        for c in ssn.load_clusters(family):
            for acc, seq in ssn.iter_sequences(c):
                if seq:
                    fh.write(f">{c.cluster_id}{_SEP}{acc}\n{seq}\n")
                    n += 1
    if n == 0:
        out.unlink(missing_ok=True)
        raise ValueError(f"no SSN members for family {family!r} -- it has no SSN in this build "
                         f"(the twelve are {', '.join(ssn.FAMILIES)})")
    (out.with_suffix(".count")).write_text(str(n))
    return out


def _target_cluster(target: str) -> str:
    return target.split(_SEP, 1)[0]


# --------------------------------------------------------------------------- assignment
@functools.lru_cache(maxsize=None)
def _cluster_sizes() -> dict:
    """{cluster_id: n_seqs}. Cached: counting is one pass over every clade FASTA."""
    return {c.cluster_id: c.n_seqs for c in ssn.load_clusters(with_counts=True)}


def assign_cluster(seq: str, family: str, *, uniprot_acc: str | None = None,
                   top_n: int = 50, sensitivity: float = 7.5,
                   min_support: float = 0.5, min_identity: float = MIN_IDENTITY,
                   min_cluster_size: int = MIN_CLUSTER_SIZE) -> ClusterAssignment:
    """Assign `seq` to an SSN cluster, or ABSTAIN and say why.

    Three-tier lookup (cheapest first):
      1. Direct accession lookup: if `uniprot_acc` is known and present in the SSN member map →
         O(1), no MMseqs. This is the common case for natural TFs, which are almost certainly
         already in the SSN or will be mapped via their nearest curated neighbour.
      2. High-identity MMseqs hit: if top hit has pident >= 0.90, the assignment is unambiguous
         without a full vote (the same or nearly-identical sequence in the SSN).
      3. Bit-score voting over top_n hits: the fallback for more divergent sequences.

    Tiers 2 and 3 are then gated by `min_identity` (see MIN_IDENTITY: measured, and the reason this
    function can now abstain at all) and optionally `min_cluster_size`. Tier 1 is exempt -- a direct
    accession hit means the protein IS a member, which is not a similarity judgement.

    `ClusterAssignment.status` always says which of those happened; `cluster_id=None` no longer
    conflates "not in the SSN", "this family has no SSN" and "too divergent to call"."""
    # --- tier 1: direct accession lookup (natural TFs already in the SSN) -----------------------
    if uniprot_acc:
        amap = ssn.build_accession_map()
        cid = amap.get(uniprot_acc)
        if cid:
            return ClusterAssignment(cid, 1.0, 1, 1.0, {cid: 1.0},
                                     ssn.info_for_cluster(cid), method="direct")

    # --- tier 2 & 3: MMseqs search against the SSN member FASTA --------------------------------
    # These tiers search a FAMILY SSN DB. All twelve `ssn.FAMILIES` have carried one since the 2026 drop
    # (MerR, ArsR/SmtB, Fur, CopY, CsoR/FrmR, DtxR/MntR, GntR, LysR-type, MarR/SlyA, NikR, Rrf2, TetR/AcrR);
    # before that only MerR and ArsR/SmtB did and `build_ssn_db` fell back to the ArsR members, so a
    # divergent TF could spuriously match an ArsR cluster and inherit its metal inducer. The guard below is
    # what prevented that and still bounds the search to families that actually have an SSN; tier 1 (the
    # direct accession lookup above) already handled any TF literally in the SSN, regardless of family.
    if family not in _FAMILIES:
        return ClusterAssignment(None, 0.0, 0, 0.0, status=STATUS_NO_SSN)
    db = build_ssn_db(family)
    sweep = homologs.homolog_sweep(seq, db, sensitivity=sensitivity, evalue=1e-3,
                                   max_seqs=max(top_n * 4, 200), cluster=False)
    hits = sweep.hits[:top_n]
    if not hits:
        return ClusterAssignment(None, 0.0, 0, 0.0, status=STATUS_NO_HIT)

    top_identity = max(h.pident for h in hits)
    votes: dict = defaultdict(float)
    for h in hits:
        votes[_target_cluster(h.target)] += h.bits
    audit = {k: round(v, 1) for k, v in sorted(votes.items(), key=lambda kv: -kv[1])}

    # the identity floor gates BOTH similarity tiers: a query this far from everything in the family
    # is not a member of any clade, and naming one would hand it that clade's inducer
    if top_identity < min_identity:
        return ClusterAssignment(None, 0.0, len(hits), top_identity, audit,
                                 None, method="mmseqs", status=STATUS_LOW_IDENTITY)

    sizes = _cluster_sizes()

    # tier 2: single high-identity hit — essentially the same sequence, no vote needed
    if hits[0].pident >= 0.90:
        cid = _target_cluster(hits[0].target)
        if sizes.get(cid, 0) < min_cluster_size:
            return ClusterAssignment(None, hits[0].pident, len(hits), top_identity, audit,
                                     None, method="high_identity", status=STATUS_SMALL_CLUSTER)
        return ClusterAssignment(cid, hits[0].pident, len(hits), top_identity,
                                 {_target_cluster(h.target): round(h.bits, 1) for h in hits[:5]},
                                 ssn.info_for_cluster(cid), method="high_identity")

    # tier 3: bit-score weighted vote over all hits
    total = sum(votes.values()) or 1.0
    cid, best = max(votes.items(), key=lambda kv: kv[1])
    support = best / total
    if support < min_support:
        return ClusterAssignment(None, support, len(hits), top_identity, audit,
                                 None, method="mmseqs", status=STATUS_LOW_SUPPORT)
    if sizes.get(cid, 0) < min_cluster_size:
        return ClusterAssignment(None, support, len(hits), top_identity, audit,
                                 None, method="mmseqs", status=STATUS_SMALL_CLUSTER)
    return ClusterAssignment(cid, support, len(hits), top_identity, audit,
                             ssn.info_for_cluster(cid), method="mmseqs")


# --------------------------------------------------------------------------- member selection
def select_cluster_members(cluster_id: str, *, n: int = 60, pool: int = 4000, seed: int = 11,
                           dedup_sequences: bool = True) -> list:
    """Representative member accessions of a cluster, one per distinct SEQUENCE.

    These feed the UniProt->genome->promoter path, and the promoters feed motif discovery -- so
    redundancy here is not free. Deduplicating only identical ACCESSIONS (the previous behaviour)
    left the same protein in the sample many times over: measured across the twelve families, 52
    clusters drew duplicate sequences and `CopY_c4` sampled 60 members that were only **27 distinct
    proteins**. A motif built from 33 redundant promoters looks far more conserved than the evidence
    warrants, and handoff §6 already records that deeper/oversampled homolog sets hurt rather than
    help.

    Sequence deduplication happens BEFORE sampling, so `n` means n distinct proteins rather than n
    rows. Pass `dedup_sequences=False` for the old behaviour.
    """
    import hashlib
    import random
    info = ssn.info_for_cluster(cluster_id)
    if info is None:
        return []
    if not dedup_sequences:
        accs = list(dict.fromkeys(acc for acc, _ in ssn.iter_sequences(info, limit=pool)))
    else:
        seen_seq: set[str] = set()
        accs = []
        for acc, seq in ssn.iter_sequences(info, limit=pool):
            if not seq:
                continue
            h = hashlib.md5(seq.strip().upper().encode()).hexdigest()
            if h in seen_seq:
                continue
            seen_seq.add(h)
            accs.append(acc)
        accs = list(dict.fromkeys(accs))
    rng = random.Random(seed)
    return rng.sample(accs, n) if len(accs) > n else accs


# --------------------------------------------------------------------------- accession-path restriction
def restrict_orthologs(orthologs: list, cluster_id: str, *, amap: dict | None = None) -> list:
    """Keep only orthologs whose accession is a member of `cluster_id` (pure; for the case where the
    query's MSA orthologs ARE SSN members). `orthologs` is a list of dicts with an 'acc' key."""
    amap = amap if amap is not None else ssn.build_accession_map()
    return [o for o in orthologs if amap.get(o.get("acc")) == cluster_id]


# --------------------------------------------------------------------------- self-test / CLI
def _demo() -> None:
    # offline: a synthetic query identical to a real cluster member must assign to that member's cluster
    info = ssn.info_for_cluster("MerR_9")                     # ecZntR (Zn)
    acc, seq = next(ssn.iter_sequences(info))
    a = assign_cluster(seq, ssn.MERR, top_n=30)
    print(f"self-query ({acc}) -> {a.cluster_id} support={a.support:.2f} id={a.top_identity:.2f} "
          f"hits={a.n_hits}")
    assert a.cluster_id == "MerR_9", f"a member of MerR_9 must self-assign to MerR_9, got {a.cluster_id}"
    assert a.top_identity > 0.95, "self-hit identity should be ~1.0"

    # accession-path restriction is exact
    orth = [{"acc": acc}, {"acc": "NOT_A_MEMBER"}]
    kept = restrict_orthologs(orth, "MerR_9")
    assert len(kept) == 1 and kept[0]["acc"] == acc
    print("OK: MMseqs cluster assignment + accession-path restriction verified.")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--assign", help="FASTA/seq to assign to an SSN cluster")
    ap.add_argument("--family", default=None, choices=list(ssn.FAMILIES),
                    help="the TF's family (required with --assign)")
    ap.add_argument("--build-db", action="store_true")
    a = ap.parse_args(argv)
    if a.assign and not a.family:
        ap.error("--assign needs --family (the SSN member DB is per family)")
    if a.build_db:
        for fam in ssn.FAMILIES:
            p = build_ssn_db(fam, force=True)
            print(f"built {p.relative_to(_REPO)} ({(p.with_suffix('.count')).read_text()} members)")
        return
    if a.assign:
        seq = a.assign
        if Path(a.assign).exists():
            seq = "".join(l.strip() for l in Path(a.assign).read_text().splitlines()
                          if not l.startswith(">"))
        res = assign_cluster(seq, a.family)
        print(f"assigned cluster: {res.cluster_id}  support={res.support:.2f}  "
              f"top_identity={res.top_identity:.2f}  hits={res.n_hits}")
        if res.info:
            print(f"  inducer={res.info.inducer} role={res.info.role} gate={res.info.gate or '-'}")
        print(f"  votes: {res.votes}")
        return
    _demo()


if __name__ == "__main__":
    main()
