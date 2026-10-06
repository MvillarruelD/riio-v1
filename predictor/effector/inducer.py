"""
inducer.py -- unified inducer inference that fuses every source but keeps each retrievable (Phase 3).

Several independent signals say what the cognate inducer is:
  * **ssn_cluster** -- the SSN isofunctional cluster's consensus inducer (annotate.homolog_selection +
    annotate.ssn_clusters), with the metal label GATED on the coordination signature;
  * **coordination** -- structure.metal_site: does the metal-coordination site exist? (gates the metal call);
  * **metalnet** -- structure.metalnet: MetalNet2's family-AGNOSTIC metal-site detector, consulted where
    our own chemistry has no branch. The coordination gate covers six of the twelve families, so for the
    other six it is silent on every member (82 of the 150 run-5 candidates); MetalNet closes that gap.
    It contributes a metal / non-metal ROLE only -- `ligand` is deliberately None. Its metal-TYPE
    prediction is kept in `evidence` as advisory and never names the headline ion;
  * **neighborhood** -- annotate.neighborhood: is there metal-homeostasis machinery (a transporter,
    efflux pump, reductase or siderophore system) within 500 bp of this regulator? The survey's own
    second axis. Like MetalNet it contributes a metal ROLE only and `ligand` is None -- a copA
    neighbour says copper is handled nearby, not that this TF senses copper. It is the WEAKEST metal
    evidence carried here and sits at the bottom of the precedence ladder: it may establish a metal
    class only where our own chemistry is silent AND MetalNet found no site, and never against a
    negative family gate. Measured worth on the survey's 150: it is the only metal evidence for
    exactly ONE candidate, and corroborates 12 (analysis/benchmarking/discovery/);
  * **ligify** -- effector.ligify operon-chemistry (optional; needs genome context);
    (enzyme-reaction-only) misses and GO can't resolve;

Per the design: agreement across sources strengthens the call, but disagreement leaves all candidates
valid -- so every source's `InducerCall` is retained in the returned `InducerConsensus` (schema). The
gate is applied as policy here (not in the schema): an SSN *metal* label whose coordination signature is
absent is demoted to advisory.

Run `python -m predictor.effector.inducer --tf CueR_Ecoli --family MerR`.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from predictor import resources

_REPO = Path(__file__).resolve().parents[2]
from predictor.schema import InducerCall, InducerConsensus                        # noqa: E402
from predictor.annotate import homolog_selection as hsel     # noqa: E402
from predictor.structure import metal_site                                        # noqa: E402

_METAL_ROLES = {"metal"}
#: A Ligify (operon-chemistry) call whose gate weight falls below this is SUPPRESSED: it is noted as
#: down-weighted and, when the coordination gate says metal, cannot hold the headline.
LIGIFY_SUPPRESSED_WEIGHT = 0.5


# --------------------------------------------------------------------------- member-level inducer
def _fam_tag(family: str | None) -> str | None:
    """The on-disk tag for a family's curated-binder file, via the single source of truth
    (`ssn_clusters.family_tag`). This used to be a hand-written MerR/ArsR lookup returning None for
    everything else, so a curated-binder set generated for any other family could never be loaded no
    matter that it existed. Now every family resolves; families with no file on disk simply read empty."""
    if not family:
        return None
    from predictor.annotate import ssn_clusters as _ssn
    return _ssn.family_tag(family)


def _pct_identity(q: str, t: str) -> float:
    """Global %identity of two protein sequences, normalised by the shorter length (so a query that is a
    sub/super-sequence of a curated binder still scores ~1.0). Biopython global alignment, BLOSUM62."""
    try:
        from Bio import Align
        from Bio.Align import substitution_matrices
    except Exception:
        # crude fallback: identical-prefix fraction (only correct for already-aligned equal seqs)
        n = sum(1 for a, b in zip(q, t) if a == b)
        return n / max(1, min(len(q), len(t)))
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aln = aligner.align(q.upper(), t.upper())[0]
    n_id = sum(1 for a, b in zip(str(aln[0]), str(aln[1])) if a == b and a != "-")
    return n_id / max(1, min(len(q), len(t)))


_BINDER_CACHE: dict = {}


# --------------------------------------------------------------------------- manual ion anchors
# The curated-binder ions are mined from UniProt COFACTOR / metal-BINDING annotations. A handful of
# ESTABLISHED metal sensors carry no such annotation upstream, so mining correctly yields ion=None and the
# member-level override below cannot fire -- the query then falls back to its SSN cluster's plurality
# label, which is exactly wrong for a heterogeneous cluster.
#
# NB: ion=None is the CORRECT answer for most un-ionised binders (YgaV/HlyU/BigR = thiol/RSS, GlnR/TnrA =
# nitrogen, CarH = B12/light, mta/BluR = multidrug). Only add an accession here when the protein is an
# established METAL sensor whose ion is missing upstream -- with the citation. Keep this list short; it is
# a patch for source-data gaps, not a place to encode guesses.
_MANUAL_ION = {
    # KmtR (M. tuberculosis Rv0827c): the established Ni(II)/Co(II) sensor, NmtR's counterpart. UniProt
    # O53838 has no COFACTOR and no metal BINDING feature -- the only trace of the metal is a reference
    # TITLE ("Mycobacterial cells have dual nickel-cobalt sensors...", Campbell et al. 2007, JBC
    # 282:32298), which is not machine-derivable. Without this, KmtR bins to the heterogeneous ArsR_c2
    # (CzrA/SmtB/ZiaR=Zn + NmtR/KmtR=Ni + CadC=Cd) and inherits its Zn2+ plurality label.
    "O53838": {"ion": "Ni2+", "why": "Campbell 2007 JBC 282:32298; UniProt lacks COFACTOR/BINDING"},
}


def _curated_binders(family: str) -> list:
    tag = _fam_tag(family)
    if tag is None:
        return []
    if tag not in _BINDER_CACHE:
        p = resources.SSN_DATABASE / f"curated_binders_{tag}.json"
        try:
            rows = json.loads(p.read_text(encoding="utf-8")) if p.exists() else []
        except Exception:
            rows = []
        for b in rows:                                     # fill source-data gaps, never override a mined ion
            m = _MANUAL_ION.get(b.get("acc"))
            if m and not b.get("ion"):
                b["ion"], b["ion_source"] = m["ion"], f"manual_anchor ({m['why']})"
            # A COFACTOR cluster is not a sensed inducer. The ions here are mined from UniProt
            # metal-BINDING features, which do not distinguish "this protein senses X" from "this
            # protein contains X" -- so every SoxR ortholog arrives carrying ion="Fe" from its
            # [2Fe-2S] cluster. Measured: 7 such rows in the MerR set, and leave-one-sequence-out
            # showed the route then names Fe2+ for SoxR at 0.96 identity, a false metal for a
            # superoxide sensor (`project_archive/benchmarking_regulondb/curated_binder_loo.py`). Production masks
            # it today only because the redox short-circuit in `_finalize` fires first, which is luck,
            # not a safeguard. Same principle the Rrf2 gate already states: an [Fe-S] cluster is real
            # chemistry and is not a metal-ion call.
            detail = (b.get("ion_detail") or "").lower()
            if b.get("ion") and any(k in detail for k in
                                    ("fe-s", "[2fe", "[4fe", "[3fe", "iron-sulfur", "cluster")):
                b["ion_not_inducer"] = b.pop("ion")
                b["ion_source"] = f"cofactor cluster ({b.get('ion_detail')}) -- not a sensed inducer"
        _BINDER_CACHE[tag] = rows
    return _BINDER_CACHE[tag]


def curated_member_ion(seq: str, family: str, *, min_identity: float = 0.5) -> dict | None:
    """Member-level inducer: among the family's CURATED ion-bearing binders, return the one most similar
    to the query (ion, identity, acc, name) when it clears `min_identity`. Resolves the lumped/mis-labelled
    cluster failure mode -- e.g. E. coli CueR (Cu+) bins to the Cd/Hg-labelled MerR_1, but its own curated
    ortholog names Cu+ directly. Cluster-agnostic on purpose: the query's nearest curated ion-binder is a
    more reliable inducer signal than the cluster's plurality label. Returns None for a genuinely novel TF
    not close to any ion-bearing curated binder (-> falls back to the cluster label)."""
    # site_type=="structural" excludes a constitutive structural-site metal (e.g. the Fur
    # Zn(Cys4) cage) from the member-level override: it is not the sensed effector, and letting
    # it win would mislabel Mn/Ni/Fe sensors as Zn (see data/refs/family_kb site-typing).
    binders = [b for b in _curated_binders(family)
               if b.get("ion") and b.get("seq") and b.get("site_type") != "structural"]
    if not binders:
        return None
    best = None
    for b in binders:
        idn = _pct_identity(seq, b["seq"])
        if best is None or idn > best["identity"]:
            best = {"ion": b["ion"], "identity": idn, "acc": b.get("acc"), "name": b.get("name")}
    if best and best["identity"] >= min_identity:
        return best
    return None


def infer_inducer(seq: str, family: str, *, genome_ctx=None, tf_gene=None, allow_ncbi: bool = False,
                  ligify_fn=None, metalnet=None, metalnet_fn=None, protein_acc: str | None = None,
                  tf_id=None, allow_colabfold: bool = False, fold_cif: str | None = None,
                  op_cif: str | None = None, verbose: bool = False) -> InducerConsensus:
    """Fuse inducer sources into a gated, source-retaining consensus.

    Sources (each retained side-by-side): SSN cluster consensus, the metal-coordination gate, MetalNet2's
    family-agnostic metal-site detector, Ligify operon chemistry (needs `genome_ctx` + `tf_gene`; the
    network fetch is gated on `allow_ncbi`)Always proposes >=1 hypothesis: when no
    source names a ligand, a class-level call is read off the gates (so the dossier never goes silent).

    `metalnet` accepts an already-computed `structure.metalnet.predict_site` result (or a bool); pass
    `metalnet_fn` to inject one, or leave both None to have it resolved here. The engine is optional --
    when it is not installed the source abstains and every verdict is exactly what it was before."""
    calls: list[InducerCall] = []

    # --- source 1: SSN cluster consensus inducer (sequence-assigned) ---
    try:
        assign = hsel.assign_cluster(seq, family)
        if assign.cluster_id and assign.info:
            ligand, role = assign.info.inducer, assign.info.role
            ev = {"cluster": assign.cluster_id, "identity": round(assign.top_identity, 2),
                  "representative": assign.info.representative}
            # MEMBER-LEVEL override: the query's nearest CURATED ion-bearing ortholog names the ion
            # directly. SSN clusters are often lumped/polyspecific and their plurality label can be wrong
            # for a specific member (the canonical CueR/Cu+ -> Cd/Hg-labelled MerR_1 failure). When a
            # confident member match disagrees with the cluster label, trust the member.
            unresolved_with = None
            mem = curated_member_ion(seq, family)
            if mem:
                ev["curated_member"] = {"acc": mem["acc"], "name": mem["name"], "ion": mem["ion"],
                                        "identity": round(mem["identity"], 2)}
                # The override exists to fix a LUMPED cluster -- where the query's true ion is not
                # the cluster's label at all (CueR/Cu+ landing in a Cd/Hg-labelled MerR cluster). It
                # must NOT arbitrate between two ions that are both anchors of the very cluster the
                # query was assigned to: there it is choosing between co-sensed alternatives on a
                # hair's-breadth identity margin, and the cluster consensus is the better default.
                #
                # Measured on M. tuberculosis IdeR, the archetypal iron sensor: its nearest curated
                # binder is an IdeR ORTHOLOG whose census label is "Cobalt" (from the in vitro
                # promiscuity the Chem. Rev. review describes -- DtxR/IdeR is activated by Co, Fe, Mn
                # and Ni in vitro). Both are tier-1 anchors of DtxR_IdeR, whose own label is Fe2+, and
                # the override flipped the headline to Co2+ at identity 1.00. The alternative is not
                # lost: it stays in `candidates`.
                _aset = assign.info.anchor_set
                _cluster_ions = set(_aset.distinct_inducers) if _aset else set()
                if mem["ion"] != ligand and mem["ion"] in _cluster_ions and ligand in _cluster_ions:
                    # A genuine TIE, and it is reported as one. Both ions are tier-1 anchors of the
                    # cluster the query was assigned to, and identity to the nearest binder is not a
                    # basis for choosing between them -- for IdeR the "winner" was an ortholog carrying
                    # an in vitro Cobalt label, ahead of the cluster's own correct Fe2+, on a
                    # two-residue length difference.
                    #
                    # Collapsing that to one ion would be fabricated precision. The call is marked
                    # UNRESOLVED, both ions ride in `candidates`, and the downstream evidence is
                    # given the chance to separate them on co-regulated handler chemistry -- which is
                    # evidence about this organism rather than about sequence similarity.
                    unresolved_with = mem["ion"]
                    ev["unresolved_between"] = sorted({ligand, mem["ion"]})
                    ev["why_unresolved"] = (
                        "both are tier-1 anchors of this cluster and the nearest-binder identity "
                        "margin is not evidence about which is sensed; DtxR-family sites in "
                        "particular are activated in vitro by Co, Fe, Mn and Ni alike")
                    if verbose:
                        print(f"  ssn_cluster: {ligand} vs {mem['ion']} UNRESOLVED -- both are "
                              f"anchors of {assign.cluster_id}; reporting both")
                elif mem["ion"] != ligand:
                    ev["cluster_label"] = ligand                  # keep the raw cluster label for audit
                    ligand, role = mem["ion"], (role or "metal")
                    if verbose:
                        print(f"  ssn_cluster: member override {ev['cluster_label']} -> {ligand} "
                              f"(nearest curated binder {mem['acc']} id={mem['identity']:.2f})")
            # ANCHOR SUPPORT. `assign.support` measures how cleanly the query landed in the cluster;
            # it says nothing about whether the cluster's LABEL is trustworthy. That label used to
            # rest on a single representative, so there was nothing else to report. It now rests on
            # an anchor set, and 9 of 71 clusters have anchors that DISAGREE about the inducer
            # (`ArsR_c2` carries Zn(II), Ni(II) and Pb(II) under a flat `Zn2+`). Confidence is
            # therefore assignment support scaled by label agreement, and the anchor evidence is
            # recorded so a reader can see which of the two is weak.
            aset = assign.info.anchor_set
            confidence = round(assign.support, 2)
            if aset is not None:
                ev["anchors"] = {
                    "n": aset.n_anchors, "n_independent": aset.n_independent_anchors,
                    "tier": aset.confidence_tier, "conflict": aset.conflict,
                    "n_distinct_inducers": aset.n_distinct_inducers,
                    "consensus": aset.inducer_consensus,
                }
                if aset.conflict:
                    # Penalise a disagreement about the CLASS, not about the ion. When every
                    # conflicting anchor names a metal, the class is certain and only the ion is
                    # open -- and that uncertainty is already expressed, exactly and without loss, by
                    # the candidate set below. Halving the confidence as well double-counts it, and
                    # the cost is real: on M. tuberculosis IdeR the cluster call (support 0.99, two
                    # tier-1 metal anchors) was cut to 0.50 and lost the headline to a regulon vote of
                    # 0.632 built on two genes -- `arsA`/`arsB1` -- from an operator scan whose
                    # measured precision is ~2%. The pipeline reported As3+ for the archetypal iron
                    # sensor.
                    from predictor.effector import inducer_vocab as _v
                    _classes = {_v.inducer_class_of(x) for x in aset.distinct_inducers}
                    _classes.discard(None)
                    if len(_classes) > 1:
                        confidence = round(assign.support * aset.agreement, 2)
                        if verbose:
                            print(f"  ssn_cluster: {assign.cluster_id} anchors disagree about the "
                                  f"CLASS ({', '.join(sorted(_classes))}); confidence "
                                  f"{assign.support:.2f} -> {confidence:.2f}")
                    elif verbose:
                        print(f"  ssn_cluster: {assign.cluster_id} anchors name several "
                              f"{next(iter(_classes), 'unknown')} inducers; confidence kept at "
                              f"{confidence:.2f} and the ion uncertainty carried as a candidate set")
            # FLOW FIX (2026-09-15): the flat cluster label (`ClusterInfo.inducer`) is blank for some
            # clusters that nonetheless carry a clean, high-tier CURATED-ANCHOR consensus -- ArsR_c3 has
            # no plurality label but two independent anchors that both say As(III). The anchor set is the
            # better-sourced evidence (cited members, not a flat table cell; see ClusterInfo.anchor_set),
            # so a silent label must not silence the source: when there is no label the anchors supply the
            # call, at the confidence they actually support. A single agreed inducer becomes the ligand and
            # its class the role; an ion-level split whose anchors still share ONE class fixes the class and
            # leaves the ion to `candidates` below; a cross-class split asserts nothing and rides entirely
            # in `candidates`. This only fills a blank -- it never overrides a stated label.
            if not ligand and aset is not None:
                from predictor.effector import inducer_vocab as _vfix
                _role_of = {"metal": "metal", "redox": "redox", "organic": "non-metal"}
                _di = aset.distinct_inducers
                if len(_di) == 1 and not aset.real_conflict:
                    ligand = _di[0]
                    _nr = _role_of.get(_vfix.inducer_class_of(ligand))
                    if _nr and role in (None, "", "unknown"):
                        role = _nr
                    ev["ligand_from"] = "anchor_consensus"
                elif len(_di) > 1:
                    _cls = {_vfix.inducer_class_of(x) for x in _di}
                    _cls.discard(None)
                    if len(_cls) == 1:                       # ion open, but the class is unanimous
                        _nr = _role_of.get(next(iter(_cls)))
                        if _nr and role in (None, "", "unknown"):
                            role = _nr
                        ev["ligand_from"] = "anchor_class_consensus"
                    # cross-class disagreement: assert nothing here; the shortlist rides in `candidates`.

            # A cluster label may name MORE THAN ONE species -- the SSN table stores co-sensing
            # flattened into one string (`Ni2+/Co2+`, `Cd2+/Pb2+`, `Cu+/persulfide`), and the anchor
            # set can carry several distinct inducers for a genuinely mixed clade. Emit the components
            # so the fusion can credit a source whose mixture CONTAINS the headline, instead of
            # counting a flattened label as disagreeing with one of its own members.
            from predictor.effector import inducer_vocab as _vocab
            candidates = _vocab.parse_mixture(ligand)
            if aset is not None and aset.real_conflict:
                # the anchors themselves disagree: the cluster is assigned to several species, and
                # that set IS the claim
                candidates = tuple(dict.fromkeys(tuple(aset.distinct_inducers) + tuple(candidates)))
            if len(candidates) > 1:
                ev["cluster_assigned_to"] = list(candidates)
                if verbose:
                    print(f"  ssn_cluster: {assign.cluster_id} is assigned to "
                          f"{', '.join(candidates)} -- carried as a candidate set, not one label")
            if unresolved_with and unresolved_with not in candidates:
                candidates = tuple(candidates) + (unresolved_with,)
            calls.append(InducerCall(
                source="ssn_cluster", ligand=ligand, confidence=confidence,
                role=role, evidence=ev,
                candidates=candidates if len(candidates) > 1 else (),
                mixture_kind=("unresolved" if unresolved_with else
                              ("co_sensed" if len(candidates) > 1 else ""))))
        else:
            # An abstention is a result, and it has a reason. `homolog_selection` computes and names
            # that reason -- low identity, no hit, no SSN for the family -- and it used to be printed
            # and then discarded, so the run log knew why 88 of 140 candidates had no clade and the
            # bundles did not. Emit it as a call that VOTES FOR NOTHING (ligand=None, confidence=0.0,
            # the same shape the coordination gate uses when it is silent), so the reason travels
            # with the prediction without touching the fusion.
            calls.append(InducerCall(
                source="ssn_cluster", ligand=None, confidence=0.0, role="",
                evidence={"status": assign.status,
                          "top_identity": round(getattr(assign, "top_identity", 0.0) or 0.0, 3),
                          "support": round(getattr(assign, "support", 0.0) or 0.0, 3),
                          "family": family or ""}))
            if verbose:
                print(f"  ssn_cluster: no cluster ({assign.status}, "
                      f"best identity {getattr(assign, 'top_identity', 0.0):.2f})")
    except Exception as e:
        if verbose:
            print(f"  ssn_cluster source skipped ({type(e).__name__}: {e})")

    # --- source 2: coordination gate (confirms / denies the metal class; may also name an ion) ---
    #
    # This source contributed a metal/non-metal ROLE and nothing else until 2026-09-02. It now also
    # carries the `metal_motifs` catalogue's read: which documented coordination signatures fire, and --
    # for the four measured to concentrate in a single SSN clade -- which ion they argue for.
    #
    # It may name an ion ONLY where the SSN clade did not. That is a deliberate ceiling, not caution
    # for its own sake: a clade label rests on curated anchors with citations, a motif rests on a regex,
    # and the clade is the better evidence wherever it exists. Where the two disagree `_finalize` keeps
    # the clade as the headline and reports both; where there is NO clade the motif is the only
    # protein-level evidence about which metal, and on the last full run that was 88 of 140 candidates.
    g = metal_site.coordination_gate(seq, family, fold_cif=fold_cif, op_cif=op_cif)
    gate_metal = g.get("has_site")
    _ions = dict(g.get("implied_ions") or {})
    # A holo fold, when supplied, types the metal sites and contributes ONLY the regulatory-site
    # metal (never the structural Zn(Cys4) cage). It joins the sequence-motif implied ions; the clade
    # ceiling still applies in `_finalize` (a structure ion is promoted only where no clade spoke).
    for _ion, _why in (g.get("structure_implied_ions") or {}).items():
        _ions.setdefault(_ion, []).extend(_why)
    _coord_lig = sorted(_ions)[0] if len(_ions) == 1 else None
    calls.append(InducerCall(
        source="coordination",
        # `ligand` stays None here; `_finalize` promotes it only when no clade spoke. Naming it at this
        # point would let `InducerConsensus.from_calls` pick it as the headline on confidence alone and
        # silently outrank the clade, which is the one thing this axis must not do.
        ligand=None, confidence=0.0,
        role=("metal" if gate_metal else ("non-metal" if gate_metal is False else "")),
        candidates=tuple(sorted(_ions)),
        mixture_kind=("unresolved" if len(_ions) > 1 else ""),
        evidence={**{k: g.get(k) for k in ("gate", "chemistry", "n_cys", "n_his") if k in g},
                  "motifs": list(g.get("motifs") or []),
                  "implied_ions": {k: list(v) for k, v in _ions.items()},
                  # the single ion the motifs argue for, if they agree on one; `_finalize` reads this
                  "motif_ion": _coord_lig,
                  # carried so `_finalize` can apply family-scoped policy on the re-fuse path too --
                  # a re-fusing caller has no sequence and no family of its own
                  "family": family or ""}))

    # --- source 2b: MetalNet2 -- the family-AGNOSTIC metal-site detector ---
    # Runs where our chemistry does not: the coordination gate has branches for six of the twelve
    # families, so for GntR/TetR/MarR/DtxR/CopY/LysR it is silent on every member. MetalNet reads
    # co-evolving CHED constellations out of the MSA and is therefore evidence about THIS protein, not
    # about its family -- the thing those 82 candidates had none of.
    #
    # It votes on the metal/non-metal ROLE only. `ligand` is None ON PURPOSE: MetalNet's second model
    # does emit a metal type (ZN, CU, FE, SF4, ...), but that call is not trusted to name an inducer
    # here -- it is trained on structural sites of all kinds, and a Ca/Mg/structural-Zn site is not the
    # sensed effector of a regulator. The type is carried in `evidence` for the reader, and
    # `_finalize`'s class-level fallback says "divalent metal (ion unresolved)" rather than guessing.
    mn = metalnet
    if mn is None:
        try:
            if metalnet_fn is not None:
                mn = metalnet_fn(seq)
            else:
                from predictor.structure import metalnet as _mn
                mn = _mn.predict_site(seq, tf_id=tf_id, family=family,
                                      allow_colabfold=allow_colabfold, verbose=verbose)
        except Exception as e:
            mn = None
            if verbose:
                print(f"  metalnet source skipped ({type(e).__name__}: {e})")
    mn_gate = metal_site.metalnet_verdict(mn)
    if mn_gate is not None or isinstance(mn, dict):
        # order matters: the report figure shows the first four evidence keys, so lead with the verdict
        # and what it rests on, not with the advisory metal type
        ev = {k: mn.get(k) for k in ("has_site", "n_site_residues", "site_residues", "msa_source",
                                     "n_site_pairs", "max_prob", "metal_type_top", "metal_type_votes",
                                     "n_msa_used", "reason")
              if isinstance(mn, dict) and k in mn}
        ev["metal_type_is_advisory"] = True
        calls.append(InducerCall(
            source="metalnet", ligand=None, confidence=0.0,
            role=("metal" if mn_gate else ("non-metal" if mn_gate is False else "")), evidence=ev))
        if verbose:
            from predictor.structure import metalnet as _mn
            print(f"  {_mn.summary(mn) if isinstance(mn, dict) else f'metalnet: has_site={mn_gate}'}")

    # --- source 3: genomic neighbourhood (needs a resolved genome context + the TF gene) ---
    # A metalloregulator usually sits beside the transporter or reductase it controls, and the survey
    # (Rondon et al.) used exactly that as its second axis. This is the WEAKEST metal evidence we carry
    # -- co-location is circumstantial, and plenty of regulators sit next to a handler for a metal they
    # do not sense -- so it names no ion (`ligand=None`, like MetalNet) and, per the precedence in
    # `_finalize`, may only establish a metal CLASS where both our own coordination chemistry and
    # MetalNet are silent about a site. Measured on the RegulonDB panel: fires on 3/9 true metal
    # sensors, 0 false metals on the MarR/AcrR/GntR organic controls
    # (`analysis/benchmarking/discovery/neighborhood_panel.py`).
    if genome_ctx is not None and tf_gene is not None:
        try:
            from predictor.annotate import neighborhood as _nb
            nb = _nb.analyse(genome_ctx, tf_gene)
            nb_call = _nb.inducer_call(nb)
            if nb_call is not None:
                calls.append(nb_call)
            if verbose:
                if nb_call is not None:
                    _n = nb.nearest_metal
                    print(f"  neighborhood: {_n.gene or _n.product[:30]} at {_n.gap} bp "
                          f"({_n.category}) -> metal-class evidence (names no ion)")
                else:
                    print(f"  neighborhood: no metal-homeostasis gene within {_nb.NEAR_BP} bp")
        except Exception as e:
            if verbose:
                print(f"  neighborhood source skipped ({type(e).__name__}: {e})")

    # --- source 4: Ligify operon chemistry (needs a resolved genome context + the TF gene) ---
    # NB: effectors_from_context(ctx, tf_gene) -- NOT predict_effector(seq); the network UniProt fetch is
    # gated on allow_ncbi (offline runs still get SSN + coordination). `ligify_fn` injects a test stub.
    # Ligify infers an effector from the TF's own operon chemistry, which for a metal sensor is usually a
    # metabolic enzyme in a neighbouring operon rather than the sensed metal -- so unweighted it contributes
    # most of the spurious non-metal / organic-ligand calls. Rather than switching the source on or off,
    # its call is WEIGHTED by whether this TF looks like a metal or reactive-species sensor: `ligify_weight`
    # returns ~0.1 when the coordination gate fired or a reactive-Cys/[Fe-S] site is present, and ~1.0 for a
    # TF with no such site. The same rule is applied to every family; nothing is hardcoded per family.
    _lw, _lw_reason = metal_site.ligify_weight(seq, family, gate=g, metalnet=mn)

    # 4a. The PUBLISHED Ligify prediction, when this exact protein is one of their 3,164. Checked first
    # and OFFLINE -- it needs no genome context and no network, so it reaches queries the derivation
    # below cannot (an offline run, or a pasted sequence with no resolvable locus). Reported under its
    # own source name so a reader can tell a published answer from one we recomputed. It is weighted by
    # exactly the same metal/redox rule: their operon-chemistry call is no more entitled to name the
    # inducer of a metalloregulator than ours is.
    _ldb_rec = None
    try:
        from predictor.effector import ligify_db as _ldb
        # The published database is keyed by RefSeq accession, with the protein sequence as the
        # second key. `protein_acc` is optional and the SEQUENCE is what always works -- a caller that
        # passes the wrong kind of accession (a UniProt id, say) simply misses, it never mismatches,
        # because there is no fuzzy matching. Measured on the last manifest: 13 of 140 by accession,
        # 10 of 140 by sequence alone, so the accession is worth passing when the caller has it.
        _ldb_rec = _ldb.lookup(seq=seq, acc=protein_acc)
    except Exception:                                       # noqa: BLE001
        _ldb_rec = None
    if _ldb_rec is not None and _ldb_rec.ligand:
        calls.append(InducerCall(
            source="ligify_db", ligand=_ldb_rec.ligand,
            confidence=round(max(0.0, _lw * float(_ldb_rec.rank or 0) / 100.0), 2),
            role="", candidates=tuple(_ldb_rec.ligands),
            mixture_kind=("co_sensed" if len(_ldb_rec.ligands) > 1 else ""),
            evidence={"refseq": _ldb_rec.refseq, "matched_by": _ldb_rec.matched_by,
                      "their_rank": _ldb_rec.rank, "their_metrics": _ldb_rec.rank_metrics,
                      "rhea": _ldb_rec.rhea_id, "equation": _ldb_rec.equation,
                      "annotation": _ldb_rec.annotation,
                      "weight": _lw, "weight_reason": _lw_reason,
                      "source_db": "groov-bio/ligify-ui (MIT)"}))
        if verbose:
            print(f"  ligify_db: published prediction '{_ldb_rec.ligand}' "
                  f"(their rank {_ldb_rec.rank}, matched by {_ldb_rec.matched_by})")

    # 4b. The local derivation, for everything else -- which is the great majority.
    # NB: effectors_from_context(ctx, tf_gene) -- NOT predict_effector(seq); the network UniProt fetch
    # is gated on allow_ncbi, so an offline run still gets SSN + coordination (+ ligify_db above).
    if ligify_fn is not None or (genome_ctx is not None and tf_gene is not None and allow_ncbi):
        try:
            from predictor.effector import ligify
            cands = (ligify_fn(genome_ctx, tf_gene) if ligify_fn is not None
                     else ligify.effectors_from_context(genome_ctx, tf_gene))
            top = cands[0] if cands else None
            if top is not None:
                calls.append(InducerCall(
                    source="ligify", ligand=getattr(top, "ligand", None),
                    # `confidence` is documented [0,1] (schema.InducerCall) but Ligify can return a
                    # NEGATIVE rank, which put an out-of-range value in the field. Clamp to the contract.
                    # This does NOT change who votes: `from_calls` enfranchises `confidence > 0`, so a
                    # clamped 0.0 abstains exactly as the negative did -- which is the intended meaning
                    # ("no usable evidence"), and the guard against a random operon metabolite becoming a
                    # metalloregulator's inducer. The GntR flip is fixed by the organic-family guard in
                    # `_finalize`, not here.
                    confidence=round(max(0.0, _lw * float(getattr(top, "rank", 50)) / 100.0), 2),
                    role="",
                    evidence={"enzyme": getattr(top, "enzyme", None), "rhea": getattr(top, "rhea", None),
                              "equation": getattr(top, "equation", None),
                              "weight": _lw, "weight_reason": _lw_reason}))
                if verbose:
                    print(f"  ligify: '{getattr(top, 'ligand', None)}' x weight {_lw:.2f} ({_lw_reason})")
        except Exception as e:
            if verbose:
                print(f"  ligify source skipped ({type(e).__name__}: {e})")

    # --- source 4 (removed): structural classifier -------------------------------------------------
    # Folddisco structural retrieval was a fourth source. It is gone, and not as a judgement call:
    # measured over the last run's 139 replayable candidates it ANSWERED ZERO and changed no
    # headline. A source that never speaks cannot be defended to a referee, and keeping it "for
    # completeness" cost a WSL bridge, a Foldseek index and two modules.

    sensor_cls, _cls_why = metal_site.inducer_class(seq, family, g, mn)
    return _finalize(calls, gate_metal, gate_name=g.get("gate"), gate_chem=g.get("chemistry"),
                     sensor_class=sensor_cls)


def _finalize(calls, gate_metal, *, gate_name=None, gate_chem=None, sensor_class=None) -> InducerConsensus:
    """Fuse the per-source calls into a consensus and apply the gate policy + class-level fallback.
    Used by `infer_inducer`; kept separate so the gate notes and the never-silent fallback live in
    one place rather than being duplicated at the call site."""
    cons = InducerConsensus.from_calls(calls, coordination_gate=gate_metal)

    # Recover the redox verdict when the caller did not supply sensor_class. A caller that re-fuses
    # without recomputing it (it needs the sequence) would otherwise see a reactive-species sensor's
    # headline revert from the redox string to the SSN
    # cluster's cofactor-ion label -- the E. coli SoxR regression: [2Fe-2S] cluster labelled "Fe" instead
    # of "reactive oxygen/nitrogen species (redox)". A source that voted role="redox" (the SSN
    # [2Fe-2S]/reactive-Cys cluster) IS that signal, so honour it. Metal/organic roles are untouched.
    if sensor_class is None and any((c.role or "") == "redox" for c in calls):
        sensor_class = "redox"
    # Same recovery for the MetalNet route: a caller that re-finalizes without recomputing sensor_class
    # would otherwise drop a metal class that only MetalNet established, and the organic-family guard
    # would then fire on a TF that does have a metal site.
    if sensor_class in (None, "organic") and gate_metal is None and any(
            c.source == "metalnet" and c.role == "metal" for c in calls):
        sensor_class = "metal"
    # The neighbourhood is the weakest metal evidence we carry, so it sits at the bottom of the same
    # precedence ladder: it may establish a metal class only where our own family chemistry is SILENT
    # (never against a negative gate -- that gate exists to catch the non-metal members of metal
    # families) and where MetalNet has not already settled it either way. A MetalNet negative does not
    # veto it, because a MetalNet negative is documented weak evidence: measured sensitivity ~70 %, and
    # it misses genuine sensors including ZntR, ArsR and ModE (METALNET_INTEGRATION.md 5b).
    if sensor_class in (None, "organic") and gate_metal is None and not any(
            c.source == "metalnet" and c.role == "metal" for c in calls) and any(
            c.source == "neighborhood" and c.role == "metal" for c in calls):
        sensor_class = "metal"

    # reactive-species (redox) sensors carry a specific [2Fe-2S]/reactive-Cys signature (SoxR, NsrR-like)
    # that a RegulonDB benchmark showed is separable from true metal sensors -- their effector is an
    # oxidative/nitrosative signal, not a metal ion or an organic ligand, so it takes the headline directly.
    if sensor_class == "redox":
        cons.top = "reactive oxygen/nitrogen species (redox)"
        cons.notes = ("[2Fe-2S]/reactive-Cys signature -> reactive-species (redox) sensor: the effector is "
                      "an oxidative/nitrosative signal (superoxide, H2O2 or NO), not a metal ion or organic "
                      "ligand",)
        return cons

    # --- gate policy: a metal label without the coordination signature is advisory, not confident ---
    notes = []
    top_call = next((c for c in cons.calls if c.ligand == cons.top), None)
    # The MetalNet verdict is recovered from the retained calls, so a re-fuse
    # (which has no sequence and cannot recompute it) applies exactly the same policy as the first pass.
    coord_call = next((c for c in cons.calls if c.source == "coordination"), None)
    fam_name = ((coord_call.evidence or {}).get("family") or "") if coord_call else ""
    struct_caveat = metal_site.structural_site_caveat(fam_name)
    mn_call = next((c for c in cons.calls if c.source == "metalnet"), None)
    mn_gate = None if mn_call is None else (True if mn_call.role == "metal"
                                            else False if mn_call.role == "non-metal" else None)
    # The structural verdict the guards below act on: our own family chemistry when it has an opinion,
    # MetalNet only where it is silent. Same precedence as `metal_site.inducer_class` rule 4b.
    struct_gate = gate_metal if gate_metal is not None else mn_gate
    if top_call and top_call.role in _METAL_ROLES and gate_metal is False:
        notes.append(f"metal label '{cons.top}' NOT confirmed by coordination gate -> advisory")
    if top_call and top_call.role in _METAL_ROLES and gate_metal is True:
        notes.append(f"metal label confirmed by {gate_name} ({gate_chem or 'thiol'})")
    if mn_call is not None:
        ev = mn_call.evidence or {}
        if gate_metal is None and mn_gate is True and struct_caveat:
            # State the finding AND the reason it is not being read as a sensor call. Silence here
            # would be worse than either: the site is real, MetalNet found it, and a reader comparing
            # our organic headline against the survey's metal-site count deserves to see why they
            # differ rather than to wonder whether we missed it.
            notes.append(
                f"MetalNet finds a metal site ({ev.get('n_site_residues')} co-evolving CHED residues"
                + (f": {', '.join(ev.get('site_residues') or [])}" if ev.get("site_residues") else "")
                + f") -- but {struct_caveat}. Reported as a "
                  f"structural/cofactor site; it does NOT make this a metal sensor, and it does not "
                  f"suppress the organic evidence")
        elif gate_metal is None and mn_gate is True:
            mt = ev.get("metal_type_top")
            notes.append(
                f"no coordination chemistry is defined for this family, but MetalNet finds a metal site "
                f"({ev.get('n_site_residues')} co-evolving CHED residues"
                + (f": {', '.join(ev.get('site_residues') or [])}" if ev.get("site_residues") else "")
                + ") -> metal sensor class from MetalNet"
                + (f"; its metal-type suggestion '{mt}' is ADVISORY and does not set the ion" if mt else ""))
        elif gate_metal is None and mn_gate is False:
            notes.append("neither the coordination gate (no chemistry for this family) nor MetalNet "
                         "found a metal site")
        elif gate_metal is True and mn_gate is True:
            notes.append(f"MetalNet independently confirms the metal site ({gate_name})")
        elif gate_metal is True and mn_gate is False:
            notes.append(f"MetalNet found NO metal site, disagreeing with the family gate {gate_name}; "
                         f"the family-specific chemistry keeps the call -- treat the ion as less certain")
        elif gate_metal is False and mn_gate is True:
            notes.append("MetalNet finds a metal site although the family coordination gate is negative; "
                         "the family-specific gate keeps the call (it exists to catch the non-metal "
                         "members of metal families) -- MetalNet dissent recorded, not acted on")
        elif mn_gate is None:
            notes.append(f"MetalNet gave no opinion ({ev.get('reason') or 'not available'})")
    nb_call = next((c for c in cons.calls if c.source == "neighborhood"), None)
    if nb_call is not None:
        nev = nb_call.evidence or {}
        established = (gate_metal is None and mn_gate is not True)
        notes.append(
            f"genomic neighbourhood: {nev.get('nearest_metal_gene')} at {nev.get('gap_bp')} bp "
            f"({nev.get('category')})"
            + (" -> metal sensor class from co-location, no structural site was found by either method; "
               "co-location names no ion and is the weakest evidence here"
               if established else " (corroborating; the structural verdict already decided the class)"))
    # A cluster assigned to several species: say so, say which, and say that the headline is the
    # fusion of that set with the other sources rather than a pick from inside it.
    ssn_call = next((c for c in cons.calls if c.source == "ssn_cluster"), None)
    if ssn_call is not None and len(ssn_call.candidates) > 1:
        others = [c.source for c in cons.calls
                  if c.source != "ssn_cluster" and (c.ligand or c.role)]
        notes.append(
            f"SSN cluster {(ssn_call.evidence or {}).get('cluster', '?')} is assigned to "
            f"{', '.join(ssn_call.candidates)} -- a genuine multi-species cluster, not a failure to "
            f"decide. The headline '{cons.top}' is the fusion of that set with "
            f"{', '.join(others) if others else 'no other source'}; the full set is kept in "
            f"`candidates` and any member remains a live hypothesis")
    # An UNRESOLVED tie is a result, not a failure to decide. Say which ions are tied and why
    # nothing separates them. This used to record whether the regulon broke the tie; under the tier
    # contract it cannot, because the regulon is reached THROUGH the operator hits and a tier-2
    # product may not settle a tier-1 claim.
    if ssn_call is not None and ssn_call.mixture_kind == "unresolved":
        _tied = (ssn_call.evidence or {}).get("unresolved_between") or list(ssn_call.candidates)
        _why = (ssn_call.evidence or {}).get("why_unresolved") or ""
        notes.append(
            f"UNRESOLVED between {' and '.join(_tied)}: {_why}. No protein-level source separates "
            f"them. Both remain live hypotheses; reporting one as the answer would be precision we "
            f"do not have")
    # --- the coordination axis: corroborate a clade, contradict it in the open, or speak where there
    #     is none. Precedence rule 9, added 2026-09-02.
    #
    # The ceiling is deliberate and it is asymmetric: a clade label rests on curated anchors with
    # citations, a motif rests on a regex measured to concentrate in one clade. Where both exist the
    # clade is the better evidence and keeps the headline -- agreement raises confidence, disagreement
    # is REPORTED rather than resolved, because a sequence motif is not entitled to overturn curated
    # biology. Where there is no clade the motif is the only protein-level evidence about which metal,
    # and on the last full production run that was 88 of 140 candidates, all of them otherwise reduced
    # to "divalent metal (ion unresolved)".
    motif_ion = (coord_call.evidence or {}).get("motif_ion") if coord_call else None
    motif_ions = (coord_call.evidence or {}).get("implied_ions") or {} if coord_call else {}
    motif_names = (coord_call.evidence or {}).get("motifs") or [] if coord_call else []
    clade_named = ssn_call is not None and bool(ssn_call.ligand)
    if motif_ions:
        _srcs = sorted({m for v in motif_ions.values() for m in v})
        if clade_named and ssn_call.ligand in motif_ions:
            notes.append(
                f"coordination signature {', '.join(_srcs)} independently implies "
                f"{ssn_call.ligand}, agreeing with SSN cluster "
                f"{(ssn_call.evidence or {}).get('cluster', '?')} -- two axes that do not read each "
                f"other (curated clade anchors vs. the protein's own ligand set)")
        elif clade_named:
            cons.candidates = tuple(dict.fromkeys((*cons.candidates, *sorted(motif_ions))))
            cons.mixture_kind = "unresolved"
            notes.append(
                f"DISAGREEMENT: SSN cluster {(ssn_call.evidence or {}).get('cluster', '?')} says "
                f"'{ssn_call.ligand}', the coordination signature ({', '.join(_srcs)}) implies "
                f"{' or '.join(sorted(motif_ions))}. The clade keeps the headline -- it rests on "
                f"curated anchors and the motif on a sequence pattern -- but both are live and the "
                f"full set is in `candidates`")
        elif cons.top is None and motif_ion and sensor_class != "organic":
            # No clade, no other source named anything, and the motifs agree on ONE ion.
            cons.top = motif_ion
            coord_call.ligand = motif_ion
            coord_call.confidence = 0.35
            cons.agreement = 1.0
            notes.append(
                f"no SSN clade for this protein, so the headline '{motif_ion}' rests on the "
                f"coordination signature alone ({', '.join(motif_ions[motif_ion])}). This is weaker "
                f"than a clade call and is not corroborated by any independent source; treat the ion "
                f"as the hypothesis and the metal CLASS as the result")
        elif cons.top is None and len(motif_ions) > 1:
            cons.candidates = tuple(sorted(motif_ions))
            cons.mixture_kind = "unresolved"
            notes.append(
                f"no SSN clade, and the coordination signature ({', '.join(_srcs)}) is consistent with "
                f"{' or '.join(sorted(motif_ions))} without separating them -- the motif is shared by "
                f"those sensors. Reported as a class with the candidate set, not as one ion")
    elif motif_names:
        notes.append(f"coordination signature(s) {', '.join(motif_names)} present; none of them is "
                     f"clade-diagnostic, so they support the metal class without naming an ion")

    # SSN gives metal CLASS reliably but not always the exact ion (e.g. CueR/Cu -> CadR/Cd cluster)
    if top_call and top_call.source == "ssn_cluster" and top_call.role in _METAL_ROLES:
        notes.append("ion identity from cluster only; corroborate with the coordination "
                     "gate or MetalNet")
    # `agreement` is the share of VOTING sources naming `top` -- with a single voter it is 1.0 by
    # construction and means "uncontested", NOT "corroborated". Say so, or a reader mistakes it for
    # confidence (the ArsR_c9 -> organic-metabolite failure read as agreement=1.0).
    n_voting = sum(1 for c in cons.calls if c.ligand and c.confidence > 0)
    if cons.top is not None and n_voting <= 1:
        notes.append(f"single-source call ({top_call.source if top_call else '?'}): agreement=1.00 means "
                     f"UNCONTESTED, not corroborated -- no independent source named an inducer")
    # transparency: say when the operon-chemistry (Ligify) call was suppressed by the metal/redox gate
    lig_call = next((c for c in cons.calls if c.source == "ligify"), None)
    if lig_call and float(lig_call.evidence.get("weight", 1.0)) < LIGIFY_SUPPRESSED_WEIGHT:
        notes.append(f"Ligify '{lig_call.ligand}' down-weighted x{lig_call.evidence.get('weight')}: "
                     f"{lig_call.evidence.get('weight_reason')}")
    # metal gate positive but the only ligand named is a suppressed Ligify organic -> the metal wins the
    # headline (a suppressed operon-chemistry guess must not be the top call for a metal sensor)
    if (struct_gate is True and top_call is not None and top_call.source == "ligify"
            and float(top_call.evidence.get("weight", 1.0)) < LIGIFY_SUPPRESSED_WEIGHT):
        metal_call = next((c for c in cons.calls if c.role in _METAL_ROLES and c.ligand
                           and c.source != "ligify"), None)
        cons.top = metal_call.ligand if metal_call else "divalent metal (ion unresolved)"
        notes.append(f"metal gate positive -> headline is {cons.top}, not the suppressed Ligify "
                     f"'{top_call.ligand}'")
        top_call = next((c for c in cons.calls if c.ligand == cons.top), metal_call)

    # metalloregulator-family guard: a known metal-sensor family (sensor_class == "metal") must not have its
    # headline set by a NON-metal ORGANIC vote from a weak single source (the regulon substrate-map or Ligify
    # operon chemistry). A dispersed His/Asp/carboxylate site (e.g. MntR's Mn(II), DtxR/MntR) leaves the
    # sequence gate silent -> gate_metal is None -> the gate-positive branch above never fires, and without
    # this the headline is hijacked (the E. coli MntR regression: "multidrug" from one efflux-pump regulon
    # vote wins over the Mn(II) sensor). Demote to the class-level metal call unless a metal-role source
    # named an ion.
    if (sensor_class == "metal" and top_call is not None and top_call.role not in _METAL_ROLES
            and top_call.source == "ligify"):
        metal_call = next((c for c in cons.calls if c.role in _METAL_ROLES and c.ligand
                           and c.source != "ligify"), None)
        demoted = cons.top
        cons.top = metal_call.ligand if metal_call else "divalent metal (ion unresolved)"
        notes.append(f"metalloregulator family -> headline is {cons.top}, not the non-metal '{demoted}' "
                     f"from {top_call.source} (metal sensor class; coordination gate silent)")
        top_call = next((c for c in cons.calls if c.ligand == cons.top), metal_call)

    # ORGANIC-FAMILY GUARD -- the exact mirror of the metalloregulator guard above. A TF with NO metal
    # evidence must not have its headline set to a METAL by a weak single source. The regulon substrate map
    # votes on what the CO-REGULATED genes transport, and an organic sensor's regulon can easily contain a
    # transporter that happens to move a metal; that is corroboration at best, never the naming source.
    # E. coli GntR is the case this was written for (online A/B, 2026-07-29): gluconate/LacI-GalR, no
    # coordination site, correct "organic" headline -- until retained p/q admitted enough operons to give it
    # a regulon, whose Zn2+ vote (conf 0.393) became the ONLY positive voter (Ligify sat at -0.10) and took
    # the headline at agreement 1.00. The existing gate policy could not stop it: it only appends a NOTE, and
    # it is written `gate_metal is False`, so it never fires when the gate is SILENT (None) -- which is the
    # case for every family the coordination gate has no branch for, LacI/GalR included. Hence `is not True`.
    # Deliberately NOT excluding ssn_cluster: a cluster's metal label is real evidence about the protein,
    # unlike a substrate vote about its neighbours.
    # `struct_gate is not True` rather than `gate_metal is not True` since 2026-08: with MetalNet wired
    # in, the guard's trigger is "NEITHER method found a site" (handoff §7.1). A TF whose family has no
    # coordination chemistry but for which MetalNet DOES find a site is no longer a TF with no metal
    # evidence, so the guard must stand aside there -- MetalNet is direct evidence about the protein, of
    # the same kind the guard already exempts `ssn_cluster` for, not a substrate vote about its neighbours.
    if (sensor_class == "organic" and struct_gate is not True
            and top_call is not None and top_call.role in _METAL_ROLES
            and top_call.source == "ligify"):
        # prefer another source that actually voted (positive confidence) for a non-metal ligand; a
        # non-positive call is a source saying "no usable evidence" and must not be promoted either
        alt = next((c for c in cons.calls if c.ligand and c.confidence > 0
                    and c.role not in _METAL_ROLES and c is not top_call), None)
        demoted = cons.top
        cons.top = alt.ligand if alt else None          # None -> class-level organic fallback below
        notes.append(f"organic-ligand sensor class with no coordination signature -> '{demoted}' from "
                     f"{top_call.source} is corroboration, not the headline (a substrate vote about the "
                     f"co-regulated genes cannot introduce a metal for a TF with no metal site)")
        top_call = alt

    # A POSITIVE coordination gate also outranks a NON-METAL SSN CLUSTER label. The guard above deliberately
    # excludes `ssn_cluster` (a cluster label is normally the best inducer evidence there is), but that left a
    # hole exactly where the SSN is known to fail: the metal boundary INSIDE a family, where a query lands in
    # a neighbouring clade with the wrong chemistry. E. coli RcnR is the case -- it assigns to `CsoR_ecFrmR`
    # (the rcnR-rich FrmR clade) and took its "formaldehyde" label as the headline, while BOTH independent
    # signals said metal: its CsoR-His/Cys coordination site was present AND the regulon substrate-map voted
    # Co2+. A positive, family-specific coordination site is direct structural evidence about THIS protein; a
    # cluster label is a plurality property of a group it was assigned to by similarity. When they conflict
    # and the gate has actually fired, the gate wins -- the same strength-of-evidence ordering
    # metal_site.inducer_class applies. Requires gate_metal is True, not merely sensor_class == "metal": when
    # the gate is SILENT (a dispersed site, or no branch for the family) the cluster label is still the better
    # evidence and must be left alone.
    if (gate_metal is True and top_call is not None and top_call.source == "ssn_cluster"
            and top_call.role not in _METAL_ROLES):
        metal_call = next((c for c in cons.calls if c.role in _METAL_ROLES and c.ligand
                           and c.source not in ("ligify", "ssn_cluster")), None)
        demoted = cons.top
        cons.top = metal_call.ligand if metal_call else "divalent metal (ion unresolved)"
        notes.append(
            f"coordination gate POSITIVE ({gate_name or 'gate'}{', ' + gate_chem if gate_chem else ''}) "
            f"-> headline is {cons.top}, not the non-metal '{demoted}' inherited from SSN cluster "
            f"{top_call.evidence.get('cluster', '?')}; the cluster label disagrees with this protein's own "
            f"metal site (metal-boundary misassignment -- corroborate the ion)")
        top_call = next((c for c in cons.calls if c.ligand == cons.top), metal_call)

    # --- class-level fallback: never go silent (user policy "propose at least one probable inducer") ---
    if cons.top is None:
        if struct_gate is True or sensor_class == "metal":
            cons.top = "divalent metal (ion unresolved)"
            notes.append("no source named an ion; metal sensor class (gate, MetalNet or metalloregulator "
                         "family) -> class-level metal call (low confidence -- corroborate with "
                         "ligify / cluster annotation)")
        elif gate_metal is False or sensor_class == "organic":
            cons.top = "non-metal effector (undetermined)"
            notes.append("organic-ligand sensor family / no metal-coordination signature -> non-metal "
                         "effector class (identity undetermined; e.g. the TetR/AcrR lipophilic-ligand "
                         "sensors, whose Ligify operon-chemistry candidate is advisory)")
        else:
            notes.append("inducer undetermined (no cluster label, no operon chemistry, gate inconclusive)")
    cons.notes = tuple(notes)
    # Human-facing label only. `cons.top` remains the sole class/count input.
    from predictor.effector import inducer_vocab as _vocab
    cons.top_display = _vocab.display_label(cons.top, cons.candidates)
    return cons


def _fmt(cons: InducerConsensus) -> str:
    lines = [f"  TOP: {cons.top}  agreement={cons.agreement:.2f}  gate={cons.coordination_gate}"]
    for c in cons.calls:
        lines.append(f"    [{c.source:12s}] ligand={c.ligand} role={c.role or '-'} conf={c.confidence} "
                     f"{c.evidence}")
    for n in cons.notes:
        lines.append(f"  note: {n}")
    return "\n".join(lines)


def _self_test():
    from predictor.annotate.af3_msa import parse_a3m, _op_msa

    def qseq(tf):
        a = _op_msa(tf)
        return "".join(c for c in next(parse_a3m(a))["seq"] if c.isalpha()).upper() if a else None

    # The MetalNet source is stubbed to ABSTAIN here on purpose. A real call loads ESM-2 650M and the
    # AutoGluon bag (~2 min per sequence), which would blow the 120 s per-module budget `tfop selftest`
    # allows and would make this self-test depend on an optional engine. The engine itself is covered by
    # `structure.metalnet --self-test`; the fusion policy is covered below with injected verdicts.
    def _no_metalnet(_seq):
        return None

    for tf, fam in [("CueR_Ecoli", "MerR"), ("CzrA_Saureus", "ArsR/SmtB"), ("NmtR_Mtb", "ArsR/SmtB")]:
        s = qseq(tf)
        if not s:
            continue
        cons = infer_inducer(s, fam, metalnet_fn=_no_metalnet)
        print(f"\n{tf} ({fam}):")
        print(_fmt(cons))
        assert any(c.source == "ssn_cluster" for c in cons.calls)
        assert any(c.source == "coordination" for c in cons.calls)

    # --- metal-boundary guard: a POSITIVE gate outranks a non-metal SSN cluster label (E. coli RcnR) ------
    # RcnR assigns to the rcnR-rich CsoR_ecFrmR clade and inherited its "formaldehyde" label as the headline
    # even though its own CsoR-His/Cys site was present. Verified online: the as-run bundle said
    # formaldehyde. The gate cannot NAME the ion by itself, so the correct outcome is the class-level
    # call -- naming Co2+ needed the regulon, which is tier-2 and no longer feeds this decision.
    rcnr_calls = [
        InducerCall(source="ssn_cluster", ligand="formaldehyde", confidence=1.0, role="non-metal",
                    evidence={"cluster": "CsoR_ecFrmR"}),
        InducerCall(source="coordination", ligand=None, confidence=0.0, role="metal",
                    evidence={"gate": "CsoR-His/Cys", "chemistry": "His/Cys (Cu(I) or Ni/Co)"}),
    ]
    got = _finalize(rcnr_calls, True, gate_name="CsoR-His/Cys", sensor_class="metal")
    assert got.top == "divalent metal (ion unresolved)",         f"positive gate must outrank the non-metal cluster label, got {got.top!r}"
    assert any("metal-boundary" in n for n in got.notes), f"the override must be explained: {got.notes}"

    # ... but a SILENT gate must NOT trigger it: there the cluster label is still the better evidence.
    silent = _finalize([c for c in rcnr_calls if c.source != "coordination"], None, sensor_class=None)
    assert silent.top == "formaldehyde", f"silent gate must leave the cluster label alone, got {silent.top!r}"

    # ... and a NEGATIVE gate must leave it alone too (E. coli FrmR: gate says no metal site, cluster says
    # formaldehyde -- they agree, and the answer is formaldehyde).
    frmr = _finalize(
        [InducerCall(source="ssn_cluster", ligand="formaldehyde", confidence=1.0, role="non-metal",
                     evidence={"cluster": "CsoR_ecFrmR"}),
         InducerCall(source="coordination", ligand=None, confidence=0.0, role="non-metal",
                     evidence={"gate": "CsoR-His/Cys"})], False, sensor_class="organic")
    assert frmr.top == "formaldehyde", f"FrmR must keep formaldehyde, got {frmr.top!r}"
    print("OK: positive coordination gate outranks a non-metal SSN cluster label (metal boundary); "
          "silent and negative gates leave it alone.")

    # --- MetalNet2 as the fifth source (2026-08) -------------------------------------------------------
    # The gap it closes: for the six families the coordination gate has no chemistry for, `gate_metal` is
    # None, so every guard written against `gate_metal is False` was inert and the organic-family guard
    # fired on absence of evidence rather than on evidence of absence.
    def _mn_call(role, **ev):
        return InducerCall(source="metalnet", ligand=None, confidence=0.0, role=role,
                           evidence={"metal_type_is_advisory": True, **ev})

    # the organic-family guard: a lone weak single-source vote must not introduce a metal for a TF
    # with no metal evidence. The regression that named this was a regulon substrate vote; the
    # regulon no longer feeds the inducer call, so it is expressed with the weak source that remains.
    gntr_calls = [
        InducerCall(source="coordination", ligand=None, confidence=0.0, role="", evidence={"gate": ""}),
        InducerCall(source="ligify", ligand="Zn2+", confidence=0.393, role="metal",
                    evidence={"weight": 1.0}),
    ]
    guarded = _finalize(list(gntr_calls), None, sensor_class="organic")
    assert guarded.top == "non-metal effector (undetermined)", \
        f"organic-family guard must still demote a lone substrate vote, got {guarded.top!r}"
    # ...but once MetalNet finds an actual site on the protein, the guard's premise is gone. Its trigger
    # is now "NEITHER method fired", so the metal call stands.
    with_mn = _finalize(list(gntr_calls) + [_mn_call("metal", has_site=True, n_site_residues=4,
                                                     site_residues=["C12", "C15", "H40", "C44"],
                                                     metal_type_top="ZN")],
                        None, sensor_class="metal")
    assert with_mn.top == "Zn2+", f"a MetalNet site must let the metal call stand, got {with_mn.top!r}"
    assert any("MetalNet finds a metal site" in n for n in with_mn.notes), \
        f"the MetalNet override must be explained: {with_mn.notes}"
    # MetalNet must never NAME the ion: its type model is trained on structural sites of every kind, so
    # a Ca/Mg/structural-Zn call is not the sensed effector. With no other source the honest answer is
    # the CLASS, not MetalNet's "ZN".
    alone = _finalize([InducerCall(source="coordination", ligand=None, confidence=0.0, role="",
                                   evidence={"gate": ""}),
                       _mn_call("metal", has_site=True, n_site_residues=4, metal_type_top="ZN")],
                      None, sensor_class="metal")
    assert alone.top == "divalent metal (ion unresolved)", \
        f"MetalNet's metal type must stay advisory, got {alone.top!r}"
    assert not any("ZN" == (c.ligand or "") for c in alone.calls), "the metalnet source must not name a ligand"
    # An engine that is absent or abstains must leave the run exactly as it was.
    abstain = _finalize(list(gntr_calls) + [_mn_call("", has_site=None, reason="engine unavailable")],
                        None, sensor_class="organic")
    assert abstain.top == guarded.top, "an abstaining MetalNet must change nothing"
    # A family gate that RAN keeps the verdict; the disagreement is recorded, not acted on.
    conflict = _finalize([InducerCall(source="ssn_cluster", ligand="formaldehyde", confidence=1.0,
                                      role="non-metal", evidence={"cluster": "CsoR_ecFrmR"}),
                          InducerCall(source="coordination", ligand=None, confidence=0.0, role="non-metal",
                                      evidence={"gate": "CsoR-His/Cys"}),
                          _mn_call("metal", has_site=True, n_site_residues=3)], False,
                        gate_name="CsoR-His/Cys", sensor_class="organic")
    assert conflict.top == "formaldehyde", f"the family gate must win a conflict, got {conflict.top!r}"
    assert any("dissent recorded, not acted on" in n for n in conflict.notes), conflict.notes
    print("OK: MetalNet is the fifth source -- it gates the metal CLASS where our chemistry is silent, "
          "never names the ion, and never overrides a family gate that ran.")

    print("\nOK: unified inducer fuses SSN + coordination + MetalNet "
          "(per-source retained, gate applied).")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--tf"); ap.add_argument("--family", default="MerR")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args(argv)
    if a.self_test or not a.tf:
        _self_test()
        return
    from predictor.annotate.af3_msa import parse_a3m, _op_msa
    s = "".join(c for c in next(parse_a3m(_op_msa(a.tf)))["seq"] if c.isalpha()).upper()
    print(_fmt(infer_inducer(s, a.family, verbose=True)))


if __name__ == "__main__":
    main()
