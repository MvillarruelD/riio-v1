"""
operators.py -- the multi-source operator-candidate table (WS5 of the redesign).

Design runs AFTER structure, and the operator a user can act on may come from several independent sources.
This module assembles them into ONE ranked artifact (`operators.json`) so the provenance + score of every
candidate is explicit rather than buried in the dossier. Sources (each tagged with `source` + `provenance`):

  * conservation_natural  -- REAL sites the genome-wide rescan found (the cognate operator + its paralogues);
                             score = the rescan site score; provenance = genome locus + TF-neighborhood flag.
  * conservation_logo     -- the consensus designed from the genomic operator logo (argmax of the PWM);
                             score = the logo mean information content (bits/col).

Two further sources, `structural` and `deeppbs_iterative`, were listed here until 2026-09-02. Both were
Stage-B products of `finalize_af3`, whose engines (DeepPBS, FoldX) had already been deleted from the
package -- so both rows were emitted `pending_af3` on EVERY run and could never become anything else. A
table row that is structurally incapable of being filled is not "complete by design", it is noise in the
one table a user is told to act on. Removed with Stage B; see CHANGELOG 2026-09-02.

`build_operators(dossier, struct=None, family=...)` is pure (no I/O) and testable; the caller writes the JSON.
Run `python -m predictor.report.operators` for a self-test (offline, synthetic dossier)."""
from __future__ import annotations

from dataclasses import dataclass, field, asdict


@dataclass
class OperatorCandidate:
    source: str                       # conservation_natural | conservation_logo
    sequence: str | None              # the operator DNA (None when the source produced nothing)
    score: float | None               # source-specific score (see score_type)
    score_type: str                   # "rescan_site_score" | "mean_ic_bits"
    provenance: dict = field(default_factory=dict)
    status: str = "ready"             # "ready" | "unreliable" | "absent"


#: `motif_rescan` tags every hit with the locality tier of the window it came from. Tier 1 is the
#: intergenic region flanking the regulator itself, i.e. its own promoter.
PROXIMAL_TIER = 1


def locality_of(generator: str | None) -> str:
    """'proximal' | 'distal' | 'unknown' from a `rescan_tierN` generator tag.

    This distinction is the difference between two very different claims. Measured on the RegulonDB
    panel (handoff §13.2), **12 of 13 primary-site recoveries were tier-1 proximal** and recovery given
    a distal top hit was 1/16 -- so the operator axis is behaving largely as an autoregulation
    detector. A proximal recovery says the locality prior put the regulator's own promoter first; a
    distal one says the motif found a site the prior did not hand it. Reporting a single number hides
    which of the two happened, so the report must not.

    `unknown` is never folded into either: pre-refactor bundles carry no generator tag, and counting
    those as distal would credit them with predictions they never made.
    """
    if not generator:
        return "unknown"
    try:
        tier = int(str(generator).rsplit("tier", 1)[-1])
    except (ValueError, IndexError):
        return "unknown"
    return "proximal" if tier <= PROXIMAL_TIER else "distal"


def locality_summary(dossier: dict) -> dict:
    """Counts of proximal / distal / unknown over the genome-rescan hits of one bundle."""
    hits = (dossier.get("rescan") or {}).get("hits") or []
    out = {"proximal": 0, "distal": 0, "unknown": 0, "n": len(hits)}
    for h in hits:
        out[locality_of(h.get("generator"))] += 1
    top = max(hits, key=lambda h: (h.get("score") or 0), default=None)
    out["primary_locality"] = locality_of(top.get("generator")) if top else "unknown"
    return out


def _natural_candidates(dossier: dict, *, top: int = 5) -> list[OperatorCandidate]:
    """Top genome-rescan sites as actionable natural operators (the cognate operator + its paralogues)."""
    hits = (dossier.get("rescan") or {}).get("hits") or []
    neigh = dossier.get("neighborhood") or {}
    win = neigh.get("window")
    out = []
    for h in sorted(hits, key=lambda h: -(h.get("score") or 0))[:top]:
        seq = h.get("seq") or ""
        if not seq:
            continue
        in_n = bool(win and h.get("dyad") is not None and win[0] <= h["dyad"] <= win[1])
        out.append(OperatorCandidate(
            source="conservation_natural", sequence=seq.upper(),
            score=(float(h["score"]) if h.get("score") is not None else None),
            score_type="rescan_site_score",
            provenance={"accession": h.get("accession"), "start": h.get("start"), "end": h.get("end"),
                        "strand": h.get("strand"), "dyad": h.get("dyad"),
                        "in_tf_neighborhood": in_n, "generator": h.get("generator"),
                        # `in_tf_neighborhood` is NOT a synonym: it uses a looser window, and current
                        # bundles carry tier-3 hits that are still flagged in-neighborhood. The tier is
                        # what the locality claim rests on.
                        "locality": locality_of(h.get("generator"))}))
    return out


def _logo_candidate(dossier: dict) -> OperatorCandidate | None:
    """The consensus designed from the genomic operator logo (argmax of the PWM)."""
    glog = dossier.get("genomic_logo") or {}
    cons = glog.get("consensus")
    des = (dossier.get("designed_operator") or {})
    seq = des.get("sequence") or cons
    if not seq:
        return None
    total_ic = glog.get("total_ic")
    ic_per_col = (total_ic / len(cons)) if (total_ic and cons) else None
    return OperatorCandidate(
        source="conservation_logo", sequence=str(seq).upper(),
        score=(round(float(ic_per_col), 3) if ic_per_col is not None else None),
        score_type="mean_ic_bits",
        provenance={"method": des.get("method") or "genomic-logo argmax",
                    "note": des.get("note"), "n_kept": glog.get("n_kept"), "total_ic": total_ic})


def build_operators(dossier: dict, *, struct: dict | None = None, family: str | None = None,
                    top_natural: int = 5) -> dict:
    """Assemble the multi-source operator table. Returns dict(candidates=[...], ranked=[...], primary=...).
    `ranked` is the actionable cross-source ranking (natural sites first by score, then the logo consensus);
    `primary` is the single best recommendation. Pure -- the caller serialises it."""
    fam = family or dossier.get("family")
    cands: list[OperatorCandidate] = []
    cands += _natural_candidates(dossier, top=top_natural)
    lg = _logo_candidate(dossier)
    if lg is not None:
        cands.append(lg)

    # actionable ranking: ready natural sites by score, then the logo consensus as the designed fallback
    ready_nat = [c for c in cands if c.source == "conservation_natural" and c.status == "ready"]
    ready_nat.sort(key=lambda c: -(c.score if c.score is not None else -1e9))
    ranked = ready_nat + ([lg] if lg is not None else [])
    primary = None
    if ranked:
        primary = asdict(ranked[0])
    return {"candidates": [asdict(c) for c in cands],
            # Ranked rows are the programmatic starting point and must retain the locus/method evidence
            # already present on candidates. Dropping it made every non-primary row impossible to trace.
            "ranked": [asdict(c) for c in ranked],
            "primary": primary,
            "family": fam,
            "design_after_structure": True,
            "structure_available": bool(struct and struct.get("folded"))}


def _demo() -> None:
    doss = {
        "family": "MerR",
        "rescan": {"hits": [
            {"accession": "NC_1", "start": 100, "end": 140, "strand": "+", "dyad": 120, "score": 9.1,
             "seq": "ACGTTGACCTTAGGTCAACGT", "generator": "rescan"},
            {"accession": "NC_1", "start": 900, "end": 940, "strand": "-", "dyad": 920, "score": 4.2,
             "seq": "TTGACCAATTGGTCAA", "generator": "rescan"}]},
        "neighborhood": {"window": [0, 500]},
        "genomic_logo": {"consensus": "TTGACCTTAGGTCAA", "total_ic": 18.0, "n_kept": 6, "n_hits": 9},
        "designed_operator": {"sequence": "TTGACCTTAGGTCAA", "method": "conservation consensus"},
    }
    struct = {"folded": True, "plddt_mean": 88.0, "qc_rmsd": 1.9, "qc_pass": True}
    out = build_operators(doss, struct=struct, family="MerR")
    srcs = [c["source"] for c in out["candidates"]]
    assert "conservation_natural" in srcs and "conservation_logo" in srcs
    assert out["primary"]["source"] == "conservation_natural"
    assert out["primary"]["sequence"] == "ACGTTGACCTTAGGTCAACGT", out["primary"]
    # the top natural site (score 9.1) outranks the weaker one (4.2) and the logo
    assert out["ranked"][0]["score"] == 9.1 and out["ranked"][-1]["source"] == "conservation_logo"
    # EVERY candidate must be producible. `structural` and `deeppbs_iterative` were emitted here on every
    # run and could never be filled, because the engines behind them no longer existed.
    assert "structural" not in srcs and "deeppbs_iterative" not in srcs, srcs
    assert all(c["status"] != "pending_af3" for c in out["candidates"]), out["candidates"]
    # the table must not change shape with the family -- no family-name verdicts
    for fam in ("ArsR/SmtB", "MerR", "TetR/AcrR", None):
        o = build_operators({**doss, "family": fam}, struct=struct, family=fam)
        assert [c["source"] for c in o["candidates"]] == srcs, (fam, o["candidates"])
    print(f"OK: multi-source operators -- sources={srcs}; primary={out['primary']['source']} "
          f"({out['primary']['sequence']}); no unfillable rows; shape identical across families.")


if __name__ == "__main__":
    _demo()
