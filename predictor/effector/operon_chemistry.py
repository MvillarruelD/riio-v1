"""
operon_chemistry.py -- predict a TF's cognate effector from its DIVERGENT OPERON's chemistry.

Operon-context scoring after d'Oelsnitz et al. (Ligify, ACS Synth. Biol. 2024 / Ligify 2.0 2026).
The SCORING FORMULA is Ligify's (`effector/rank.py`, a faithful port). The NEIGHBOURHOOD is ours and
is deliberately NOT Ligify's flat box: EXP065 measured that the old `effector/ligify.py` read a strand-
blind +/-10 kb / 10-genes-a-side window (median 20 genes), which spent ~90 formula points on window size
before any evidence was weighed and credited enzymes that could not share a promoter with the TF. This
module replaces that window with a promoter-gated operon walk, so the route is named for the evidence
type it produces (`source="operon_chem"`), not for a method it only half-implemented. The published
database lookup keeps the name Ligify -- that one really is their answer (`effector/ligify_db.py`).

End to end:  TF sequence  ->  source genome + TF locus  (genome_resolver, local mirror first)
             ->  the TF's DIVERGENT PARTNER  (opposite-strand neighbour sharing the promoter, gap
                 <= MAX_PROMOTER_BP -- Ligify's `getOperon` promoter rule), AND, since EXP071, the
                 regulator's OWN co-transcribed run as a second labelled route; neither -> no call
             ->  the regulated OPERON  (co-directional walk from the partner: shared `regulon.operon_of`
                 with a wider per-gap stop, an outer span cap reproducing Ligify's 8 kb reach, and
                 non-coding genes terminating the unit)
             ->  for each operon enzyme: its reaction metabolites  (effector.ligand, cached UniProt)
             ->  score each (TF, ligand) by operon context  (effector.rank, Ligify's formula)
             ->  ranked candidate effectors, aggregated per ligand.

Ranks describe enzyme context, not chemical specificity or calibrated compound probabilities. Reaction
participants tied on rank and distance remain unresolved in the inducer caller. Name order is display only.

Why the two bounds, and not a single number (EXP066):
  * MAX_PROMOTER_BP (1000) is Ligify's own divergent-promoter rule: a candidate with no divergent
    partner within 1 kb is dropped.
  * OPERON_CHEM_MAX_GAP_BP is the per-gap STOP of the co-directional walk. It is the only knob that
    sets where the operon ends. The 300 bp default is where the operon-structure sweep (EXP066)
    saturates. Replicon identity is checked explicitly for the TF, partner and operon members rather
    than inferred from this distance. `regulon`'s own route
    keeps its 150 bp default -- the shared `operon_of` is reused but its default is untouched, so the
    regulon predictions are unchanged.
  * OPERON_CHEM_MAX_SPAN_BP (8000) is the walk's REACH, the number `getOperon` actually bounds on
    (`THEIR_WALK_BP` in the EXP065 port): no operon gene may lie more than 8 kb from the partner. The
    per-gap stop almost always fires first; the span cap is the guarantee that a gene-dense region
    cannot run the walk away.

Run `python operon_chemistry.py` for a self-test (synthetic divergent operon; live path SKIPs).
"""
from __future__ import annotations

from dataclasses import dataclass, field

from predictor.effector.ligand import uniprot_entry, is_enzyme, candidate_ligands, _synthetic_entry
from predictor.effector.rank import calculate_rank
from predictor.annotate import genome_resolver, genome_db, context as ctxmod
from predictor.annotate.neighborhood import _gap
from predictor.signals import regulon
from predictor.signals.motif_rescan import Gene, GenomeContext

__all__ = [
    "EffectorCandidate", "MAX_PROMOTER_BP", "OPERON_CHEM_MAX_GAP_BP", "OPERON_CHEM_MAX_SPAN_BP",
    "divergent_partner", "operon_from_partner", "effectors_from_context", "predict_effector",
    "TIE_BREAK_SIGNALS", "tie_break_evidence", "tie_break_key", "operative_tie",
]

#: Largest TF->divergent-partner intergenic gap that still counts as a shared promoter. Ligify's
#: `getOperon` drops any candidate with no divergent intergenic region <= 1000 bp; this is that rule.
MAX_PROMOTER_BP = 1000

#: Per-gap STOP of the co-directional operon walk (EXP066). Chosen on operon structure, NOT fitted to
#: the inducer endpoint: within-operon gaps are enriched < ~50 bp and depleted beyond a few hundred
#: (Salgado et al. 2000, PNAS 97:6652; Moreno-Hagelsieb & Collado-Vides 2002, Bioinformatics 18:S329),
#: and the EXP066 sweep over the canonical bundles saturates by 300 bp. Explicit replicon checks
#: independently prevent the walk from accepting genes across a concatenation boundary.
#: Deliberately WIDER than `regulon.DEFAULT_OPERON_MAX_GAP` (150): the regulon route feeds FDR-screened
#: operator sites and wants precision, while this route is a heavily down-weighted, gate-suppressed
#: corroborating vote that can trade a little precision for enzyme recall. The two routes share the
#: walker, not the stop.
OPERON_CHEM_MAX_GAP_BP = 300

#: Outer REACH of the walk, measured from the partner's start -- Ligify's `getOperon` 8 kb bound
#: (`THEIR_WALK_BP`, EXP065). No operon gene may lie farther than this from the partner, so a dense
#: same-strand run cannot produce a runaway operon even though the per-gap stop is permissive.
OPERON_CHEM_MAX_SPAN_BP = 8000


@dataclass
class EffectorCandidate:
    ligand: str
    rank: int
    color: str
    chebi: list = field(default_factory=list)
    equation: str = ""
    rhea: str = ""
    enzyme: str = ""            # supporting enzyme gene name / protein_id
    enzyme_product: str = ""
    distance: int = 0           # enzyme-regulator gene distance in the operon (TF at index end)
    reaction_chebi: list = field(default_factory=list)
    occurrences: list = field(default_factory=list)
    identity_status: str = ""
    headline_eligible: bool = True
    route: str = "divergent"    # EXP071: "divergent" | "cotranscribed" (the winning occurrence's route)
    #: EXP072: within-tie ordering key (smaller first) and the signals it was built from. The Ligify
    #: rank is per ENZYME, so this is the only thing that can separate one enzyme's participants.
    tiebreak: tuple = ()
    tiebreak_evidence: dict = field(default_factory=dict)


# --------------------------------------------------------------------------- within-tie ordering
#: EXP072 (prereg section 3.2). Signals that CAN order candidates sharing one (rank, distance), most
#: specific first; the key is compared lexicographically. Each is read off the candidate's own reaction
#: occurrences in the regulator's unit(s):
#:   import      the species is carried out->in by a transporter of the unit;
#:   substrate   it is on the substrate side of a chemical (non-transport) reaction, using UniProt's
#:               `physiologicalReactions` direction where one is annotated and the written order otherwise;
#:   recurrence  the number of distinct unit enzymes whose reactions involve it.
#: The evidence (`tie_break_evidence`) is always computed and recorded on every candidate.
ALL_TIE_BREAK_SIGNALS = ("import", "substrate", "recurrence")
#: The signals actually applied. EMPTY ON PURPOSE -- a recorded negative, not an unfinished feature.
#: Measured cumulatively on the 12 compound-specific organic cases of the MG1655 panel (EXP072 results,
#: plain and self-held-out alike): `import` removed the true effector from FrlR's tie and named nothing
#: correctly; `substrate` named LldR's lactate but dropped UxuR's fructuronate (13 of 19 narrowed ties
#: rested on the written equation order alone, and UxuR's is written against its physiological
#: direction); `recurrence` dropped GlcC's glycolate. Truth inside the operative tie fell from 4 to 1
#: of 12 for one correct name. The preregistered rule ships the longest prefix whose every arm passed,
#: and the first arm failed. With an empty tuple every participant of the best enzyme ties (EXP071).
TIE_BREAK_SIGNALS = ()


def _effective_side(occurrence: dict) -> str:
    """'left' = substrate side once a sourced right-to-left direction is applied."""
    side = occurrence.get("side") or ""
    if occurrence.get("physiological_direction") == "right-to-left":
        return {"left": "right", "right": "left"}.get(side, side)
    return side


def tie_break_evidence(occurrences) -> dict:
    """The three EXP072 signals for one candidate identity, with the basis of the direction used."""
    imported = substrate = False
    basis = set()
    for o in occurrences or ():
        side = _effective_side(o)
        basis.add("physiological" if o.get("physiological_direction") else "as_written")
        compartment = o.get("compartment") or ""
        if compartment:
            imported = imported or (compartment == "out" and side == "left") or (compartment == "in" and side == "right")
        elif side == "left":
            substrate = True
    enzymes = {o.get("enzyme") for o in occurrences or () if o.get("enzyme")}
    return {"import": int(imported), "substrate": int(substrate), "recurrence": len(enzymes),
            "direction_basis": sorted(basis)}


def tie_break_key(evidence: dict, signals=None) -> tuple:
    """Lexicographic key over the enabled signals; smaller sorts first (more evidence wins)."""
    signals = TIE_BREAK_SIGNALS if signals is None else signals
    return tuple(-int((evidence or {}).get(s, 0)) for s in signals)


def operative_tie(cands) -> tuple:
    """(top, rank_tie, operative) for the inducer caller.

    `rank_tie` is every candidate sharing the best (rank, distance) -- one enzyme's participants, as
    before EXP072. `operative` is the headline-eligible part of it that also shares the best tie-break
    key. A generic participant (`AH2`) is never a member, so it cannot win the tie and silence the call.
    Candidates without a key (an injected stub) all tie, which is the EXP071 behaviour."""
    cands = list(cands or ())
    if not cands:
        return None, [], []

    def rd(c):
        return (-getattr(c, "rank", 50), getattr(c, "distance", 0))
    top = min(cands, key=rd)
    rank_tie = [c for c in cands if rd(c) == rd(top)]
    eligible = [c for c in rank_tie if getattr(c, "ligand", None) and getattr(c, "headline_eligible", True)]
    if not eligible:
        return top, rank_tie, []

    def tb(c):
        return tuple(getattr(c, "tiebreak", ()) or ())
    best = min(map(tb, eligible))
    return top, rank_tie, [c for c in eligible if tb(c) == best]


# --------------------------------------------------------------------------- neighbourhood
def _tf_index(genes_sorted: list, tf_gene: Gene) -> int | None:
    for k, g in enumerate(genes_sorted):
        if g is tf_gene:
            return k
    for k, g in enumerate(genes_sorted):            # fall back to coordinates (a copied Gene)
        if g.start == tf_gene.start and g.end == tf_gene.end:
            return k
    return None


def _gene_replicon(ctx: GenomeContext, gene: Gene) -> str | None:
    """Resolve the entire half-open gene interval, or reject unavailable/invalid provenance.

    Older single-record synthetic contexts may omit `locate`; their accession is sufficient only
    when no contig map is present. A mapped context must provide its coordinate resolver.
    """
    if not 0 <= gene.start < gene.end <= len(ctx.sequence):
        return None
    locate = getattr(ctx, "locate", None)
    offsets = getattr(ctx, "contig_offsets", None)
    if not callable(locate):
        return None if offsets else (getattr(ctx, "accession", None) or None)
    start_rep, start_pos = locate(gene.start)
    end_rep, end_pos = locate(gene.end - 1)
    if not start_rep or start_rep != end_rep or start_pos < 0 or end_pos < start_pos:
        return None
    if offsets:
        offset = offsets.get(start_rep)
        # `locate` falls back to accession/global coordinates before the first mapped offset.
        # An accession that is also a map key must not make that fallback valid provenance.
        if (offset is None or gene.start < offset
                or start_pos != gene.start - offset or end_pos != gene.end - 1 - offset):
            return None
    return start_rep


def divergent_partner(genes_sorted: list, tf_gene: Gene, *, max_promoter_bp: int = MAX_PROMOTER_BP,
                      ctx: GenomeContext | None = None):
    """The TF's divergent partner: the opposite-strand gene on the TF's promoter side, within
    `max_promoter_bp`. TF on '+' diverges from its LEFT neighbour (which must be '-'); TF on '-'
    diverges from its RIGHT neighbour (which must be '+'). Returns a Gene or None (the hard promoter
    gate -- no partner, or too far, means no operon-chemistry call). With `ctx`, both entire intervals
    must resolve to the same replicon."""
    tf_rep = _gene_replicon(ctx, tf_gene) if ctx is not None else None
    if ctx is not None and tf_rep is None:
        return None
    i = _tf_index(genes_sorted, tf_gene)
    if i is None:
        return None
    if tf_gene.strand == "+":
        j = i - 1
        if j < 0:
            return None
        nb = genes_sorted[j]
        side_ok = nb.strand == "-"
    elif tf_gene.strand == "-":
        j = i + 1
        if j >= len(genes_sorted):
            return None
        nb = genes_sorted[j]
        side_ok = nb.strand == "+"
    else:
        return None
    if ctx is not None and _gene_replicon(ctx, nb) != tf_rep:
        return None
    if side_ok and _gap(tf_gene.start, tf_gene.end, nb.start, nb.end) <= max_promoter_bp:
        return nb
    return None


def operon_from_partner(genes_sorted: list, partner: Gene, *, max_gap: int = OPERON_CHEM_MAX_GAP_BP,
                        max_span: int = OPERON_CHEM_MAX_SPAN_BP,
                        ctx: GenomeContext | None = None) -> list:
    """The partner's operon: the shared co-directional walk (`regulon.operon_of`) with a wider per-gap
    stop, truncated by the 8 kb reach and by the first non-coding gene. Returns Genes in transcription
    order starting at the partner, or [] when the partner itself is non-coding. With `ctx`, terminate
    before the first member whose entire interval does not resolve to the partner's replicon."""
    partner_rep = _gene_replicon(ctx, partner) if ctx is not None else None
    if ctx is not None and partner_rep is None:
        return []
    if not partner.protein_id:
        return []                                   # a tRNA/ncRNA head carries no operon chemistry
    walk = regulon.operon_of(genes_sorted, partner, max_gap=max_gap)
    out = []
    for g in walk:
        if ctx is not None and _gene_replicon(ctx, g) != partner_rep:
            break
        if g is not partner and not g.protein_id:
            break                                   # non-coding gene ends the transcription unit
        if abs(g.start - partner.start) > max_span:
            break                                   # beyond Ligify's getOperon reach
        out.append(g)
    return out


# --------------------------------------------------------------------------- core (testable offline)
def effectors_from_context(ctx: GenomeContext, tf_gene: Gene, *, entry_fetcher=uniprot_entry,
                           max_promoter_bp: int = MAX_PROMOTER_BP,
                           max_gap: int = OPERON_CHEM_MAX_GAP_BP,
                           max_span: int = OPERON_CHEM_MAX_SPAN_BP) -> list:
    """Score candidate effectors from a resolved genome context, using the TF's divergent operon as the
    enzyme set. `entry_fetcher(protein_id)->entry|None` is injectable so the operon/enzyme/rank logic is
    testable without network. Empty list when the promoter gate finds no divergent partner."""
    genes = sorted(ctx.genes, key=lambda g: g.start)
    if _tf_index(genes, tf_gene) is None:
        return []
    units = []
    partner = divergent_partner(genes, tf_gene, max_promoter_bp=max_promoter_bp, ctx=ctx)
    if partner is not None:
        operon = operon_from_partner(genes, partner, max_gap=max_gap, max_span=max_span, ctx=ctx)
        if operon:
            units.append(("divergent", operon))
    # EXP071: the regulator's OWN transcription run (frlABCD-R, lldP-R-D, arsR-BC, nikABCDE-R). This is
    # the arrangement Ligify's formula was written for -- the regulator inside the operon -- and the
    # divergent-only rule never read it (20 of 39 MG1655 panel regulators sit in such a run).
    run = cotranscribed_run(genes, tf_gene, max_gap=max_gap, max_span=max_span, ctx=ctx)
    if run:
        units.append(("cotranscribed", run))
    if not units:
        return []

    by_ligand, entries = {}, {}
    from predictor.effector.chemical_identity import identity_key, chemistry
    for route, members in units:
        # Ligify's formula scores an enzyme by its GENE DISTANCE to the regulator inside the operon. A
        # divergent TF is not in its target operon, so it is placed back into the ordered run: in genome
        # order it sits next to the partner (distance 1), the next operon gene is distance 2, and so on.
        # A co-transcribed TF is already in its run. Either way n includes the regulator, exactly as
        # Ligify's `getOperon` set does.
        scoring = sorted(members + [tf_gene], key=lambda g: g.start)
        pos = {id(g): k for k, g in enumerate(scoring)}
        reg_index = pos[id(tf_gene)]
        descriptions = [f"{g.product} {g.name}" for g in scoring]
        for g in members:
            if not g.protein_id:
                continue
            if g.protein_id not in entries:            # one acquisition per protein across both routes
                entries[g.protein_id] = entry_fetcher(g.protein_id)
            entry = entries[g.protein_id]
            if not entry or not is_enzyme(entry):
                continue
            enz_index = pos[id(g)]
            sc = calculate_rank(scoring, reg_index=reg_index, enz_index=enz_index, descriptions=descriptions)
            _add_candidates(by_ligand, entry, g, sc, abs(enz_index - reg_index), route,
                            identity_key, chemistry)
    for cand in by_ligand.values():                    # EXP072: order inside one enzyme's participant set
        cand.tiebreak_evidence = tie_break_evidence(cand.occurrences)
        cand.tiebreak = tie_break_key(cand.tiebreak_evidence)
    return sorted(by_ligand.values(), key=lambda c: (-c.rank, c.distance, c.tiebreak, c.ligand))


def cotranscribed_run(genes_sorted: list, tf_gene: Gene, *, max_gap: int = OPERON_CHEM_MAX_GAP_BP,
                      max_span: int = OPERON_CHEM_MAX_SPAN_BP, ctx: GenomeContext | None = None) -> list:
    """EXP071: genes co-directional with the TF and contiguous with it (every gap <= `max_gap`), walked
    BOTH ways from the TF, excluding the TF. Stops at a non-coding gene, a different replicon, or a gene
    starting more than `max_span` from the TF. Returned in genome order."""
    tf_rep = _gene_replicon(ctx, tf_gene) if ctx is not None else None
    if ctx is not None and tf_rep is None:
        return []
    i = _tf_index(genes_sorted, tf_gene)
    if i is None or tf_gene.strand not in "+-":
        return []
    out = []
    for step in (-1, 1):
        j, prev = i + step, tf_gene
        while 0 <= j < len(genes_sorted):
            g = genes_sorted[j]
            gap = (g.start - prev.end) if step == 1 else (prev.start - g.end)
            if (g.strand != tf_gene.strand or gap > max_gap or not g.protein_id
                    or abs(g.start - tf_gene.start) > max_span
                    or (ctx is not None and _gene_replicon(ctx, g) != tf_rep)):
                break
            out.append(g)
            prev, j = g, j + step
    return sorted(out, key=lambda g: g.start)


def _add_candidates(by_ligand, entry, g, sc, distance, route, identity_key, chemistry):
    for lig in candidate_ligands(entry):
        cand = EffectorCandidate(
            ligand=lig["name"], rank=sc["rank"], color=sc["color"], chebi=lig["chebi"],
            equation=lig["equation"], rhea=lig["rhea"],
            enzyme=g.name or g.protein_id, enzyme_product=g.product,
            distance=distance, reaction_chebi=lig["reaction_chebi"],
            occurrences=[{**lig, "enzyme": g.name or g.protein_id, "route": route}],
            identity_status=chemistry(lig["name"])[1], headline_eligible=lig["headline_eligible"],
            route=route)
        key = identity_key(lig["name"])
        prev = by_ligand.get(key)
        if prev is None:
            by_ligand[key] = cand
        elif (cand.rank, -cand.distance) > (prev.rank, -prev.distance):
            cand.occurrences = prev.occurrences + cand.occurrences
            by_ligand[key] = cand
        else:
            prev.occurrences.extend(cand.occurrences)


# --------------------------------------------------------------------------- public API
def _resolve_context(seq, *, db, allow_ncbi):
    """Resolve the TF's genome + locus, then load its context OFFLINE from the mirror if present,
    else fetch a window via EFetch. Returns (GenomeContext, tf_gene) or (None, None)."""
    cands = genome_resolver.resolve(seq, db=db, top_n=5)
    if not cands:
        return None, None
    c = cands[0]
    ctx = ctxmod.genome_context_from_mirror(c.accession)        # offline if mirrored
    if ctx is not None:
        tf = ctxmod.find_tf_gene(ctx, c.tf_start, c.tf_end)
        return ctx, tf
    if not allow_ncbi:
        return None, None
    pad = 12000
    s1, e1 = max(1, c.tf_start + 1 - pad), c.tf_end + pad
    gb = genome_resolver.fetch_region(c.accession, s1, e1, rettype="gb")
    ctx = ctxmod.from_genbank(gb)
    tf = ctxmod.find_tf_gene(ctx, c.tf_start - (s1 - 1), c.tf_end - (s1 - 1))
    return ctx, tf


def predict_effector(seq, *, db=None, allow_ncbi: bool = False, top: int = 10,
                     accession: str | None = None, use_published_db: bool = True) -> list:
    """TF sequence -> ranked candidate effectors from its divergent operon. Offline against the mirror
    where possible; UniProt enzyme reactions are cached.

    Looks the protein up in the published Ligify database first (`use_published_db`, default on) and
    falls through to the local operon-chemistry derivation when it is not there."""
    if use_published_db:
        try:
            from predictor.effector import ligify_db as _ldb
            rec = _ldb.lookup(seq=seq if isinstance(seq, str) else None, acc=accession)
        except Exception:                               # noqa: BLE001
            rec = None
        if rec is not None and rec.ligands:
            rank = rec.rank if rec.rank is not None else 0
            color = ("#02a602" if rank >= 70 else "#d4d400" if rank >= 50
                     else "#d48302" if rank >= 30 else "#f50b02")
            return [EffectorCandidate(ligand=lig, rank=rank, color=color, chebi=[],
                                      equation=rec.equation, rhea=rec.rhea_id,
                                      enzyme=rec.annotation, enzyme_product="", distance=0)
                    for lig in rec.ligands][:top]
    db = db or genome_db.default_db()
    ctx, tf = _resolve_context(seq, db=db, allow_ncbi=allow_ncbi)
    if ctx is None or tf is None:
        return []
    return effectors_from_context(ctx, tf)[:top]


# --------------------------------------------------------------------------- self-test
def _demo() -> None:
    # Divergent layout: araC(-) regulator at 100-400, promoter, araA/araB/araD(+) operon transcribed
    # rightward AWAY from araC; a non-coding tRNA at 2450 that must TERMINATE the walk; a distant gene.
    genes = [
        Gene(100, 400, "-", "araC", protein_id="WP_araC", product="AraC family transcriptional regulator"),
        Gene(450, 1450, "+", "araA", protein_id="WP_araA", product="L-arabinose isomerase"),
        Gene(1470, 2400, "+", "araB", protein_id="WP_araB", product="ribulokinase"),
        Gene(2450, 2520, "+", "trnX", protein_id="", product="tRNA-Gly"),      # ends the operon
        Gene(2600, 3600, "+", "araD", protein_id="WP_araD", product="L-ribulose-5-phosphate epimerase"),
    ]
    ctx = GenomeContext("SYN_ARA.1", "N" * 9000, genes)
    tf = genes[0]                                        # araC on '-' -> partner is the RIGHT neighbour

    part = divergent_partner(sorted(genes, key=lambda g: g.start), tf)
    assert part is not None and part.name == "araA", f"divergent partner should be araA, got {part}"
    operon = operon_from_partner(sorted(genes, key=lambda g: g.start), part)
    assert [g.name for g in operon] == ["araA", "araB"], f"tRNA must terminate the walk: {[g.name for g in operon]}"

    entries = {"WP_araA": _synthetic_entry()}            # only araA is an enzyme here
    cands = effectors_from_context(ctx, tf, entry_fetcher=lambda pid: entries.get(pid))
    print(f"{len(cands)} candidate effector(s):")
    for c in cands[:6]:
        print(f"  {c.ligand:<28} rank={c.rank} ({c.color})  via {c.enzyme} d={c.distance}  [{c.rhea}]")
    names = [c.ligand for c in cands]
    assert names, "no effectors predicted"
    assert "beta-L-arabinopyranose" in names, "the arabinose substrate should be a candidate"
    assert "ATP" not in names and "ADP" not in names, "currency leaked into effectors"
    top = cands[0]
    # araA is adjacent to araC across the promoter (distance 1); scoring list is [araC,araA,araB] (n=3)
    assert top.distance == 1, f"partner enzyme should be distance 1, got {top.distance}"
    assert "araD" not in [c.enzyme for c in cands], "gene past the terminating tRNA wrongly used"

    # EXP071: a TF with a CO-directional neighbour has no divergent partner, but that neighbour is in
    # the regulator's OWN transcription run, which is read as the second, labelled route.
    conv = [Gene(100, 400, "+", "tfA", protein_id="WP_tfA", product="regulator"),
            Gene(450, 1450, "+", "enzB", protein_id="WP_enzB", product="dehydrogenase")]
    ctx2 = GenomeContext("SYN_CONV.1", "N" * 3000, conv)
    own = effectors_from_context(ctx2, conv[0], entry_fetcher=lambda pid: _synthetic_entry())
    assert own and all(c.route == "cotranscribed" for c in own), \
        "a co-transcribed enzyme must be read through the regulator's own run"
    # ... and a TF with neither a divergent partner nor a co-directional neighbour in reach makes no call
    lone = [Gene(100, 400, "+", "tfA", protein_id="WP_tfA", product="regulator"),
            Gene(2000, 3000, "+", "enzB", protein_id="WP_enzB", product="dehydrogenase")]
    ctx3 = GenomeContext("SYN_LONE.1", "N" * 4000, lone)
    assert effectors_from_context(ctx3, lone[0], entry_fetcher=lambda pid: _synthetic_entry()) == [], \
        "no divergent partner and no co-transcribed run -> no call"
    # EXP072: the within-tie signals are recorded on every candidate and applied to none (shipped off)
    assert TIE_BREAK_SIGNALS == () and all(c.tiebreak == () and c.tiebreak_evidence for c in cands)
    print("OK: divergent operon walked, tRNA terminates it, currency excluded; the regulator's own run is "
          "a second route; no partner and no run closes the gate; tie-break evidence recorded, not applied.")

    try:
        out = predict_effector("MKKLTVSDLA", allow_ncbi=False)
        assert isinstance(out, list)
        print(f"OK: predict_effector() wired (returned {len(out)} candidates from local resolve).")
    except Exception as ex:
        print(f"SKIP (no local genome DB / resolve): {type(ex).__name__}: {ex}")


if __name__ == "__main__":
    _demo()
