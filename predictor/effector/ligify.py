"""
ligify.py -- predict a TF's cognate effector from genome context (local reimplementation of Ligify).

End to end:  TF sequence  ->  source genome + TF locus  (genome_resolver, local mirror first)
             ->  genome neighbourhood / operon  (context, OFFLINE from the .gff mirror when available)
             ->  for each neighbouring enzyme: its reaction metabolites  (effector.ligand, cached UniProt)
             ->  score each (TF, ligand) by operon context  (effector.rank, Ligify's formula)
             ->  ranked candidate effectors, aggregated per ligand.

This is the INPUT side of the biosensor loop (ligand -> TF -> operator/regulon). It is deliberately a
"small tool" copied locally (per project policy) so it runs offline against the dereplicated mirror and
only touches UniProt (a database, cached) for enzyme reactions. Honest expectation: Ligify recovered
31/100 validated biosensors -- treat outputs as ranked hypotheses, corroborate with the structural
effector signal (AlphaFill metal/ligand transplant) and the curated family priors.

Run `python ligify.py` for a self-test (synthetic operon + injected enzyme entries; live path SKIPs).
"""
from __future__ import annotations

from dataclasses import dataclass, field

from predictor.effector.ligand import uniprot_entry, is_enzyme, candidate_ligands, _synthetic_entry
from predictor.effector.rank import calculate_rank
from predictor.annotate import genome_resolver, genome_db, context as ctxmod
from predictor.signals.motif_rescan import Gene, GenomeContext


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
    distance: int = 0           # enzyme-regulator gene distance in the neighbourhood


# --------------------------------------------------------------------------- neighbourhood
def neighbourhood(ctx: GenomeContext, tf_gene: Gene, *, max_kb: float = 10.0,
                  max_genes: int = 10) -> list:
    """Genes flanking the TF within +/-max_kb (and at most max_genes either side), ordered by position
    -- Ligify's operon-context window. Returns the ordered Gene list (includes the TF)."""
    genes = sorted(ctx.genes, key=lambda g: g.start)
    if tf_gene not in genes:
        return [tf_gene]
    i = genes.index(tf_gene)
    span = max_kb * 1000
    lo = i
    while lo - 1 >= 0 and i - (lo - 1) <= max_genes and tf_gene.start - genes[lo - 1].end <= span:
        lo -= 1
    hi = i
    while hi + 1 < len(genes) and (hi + 1) - i <= max_genes and genes[hi + 1].start - tf_gene.end <= span:
        hi += 1
    return genes[lo:hi + 1]


# --------------------------------------------------------------------------- core (testable offline)
def effectors_from_context(ctx: GenomeContext, tf_gene: Gene, *, entry_fetcher=uniprot_entry,
                           max_kb: float = 10.0, max_genes: int = 10) -> list:
    """Score candidate effectors from a resolved genome context. `entry_fetcher(protein_id)->entry|None`
    is injectable so the operon/enzyme/rank logic is testable without network."""
    hood = neighbourhood(ctx, tf_gene, max_kb=max_kb, max_genes=max_genes)
    reg_index = hood.index(tf_gene)
    descriptions = [f"{g.product} {g.name}" for g in hood]

    by_ligand = {}
    for enz_index, g in enumerate(hood):
        if g is tf_gene or not g.protein_id:
            continue
        entry = entry_fetcher(g.protein_id)
        if not entry or not is_enzyme(entry):
            continue
        sc = calculate_rank(hood, reg_index=reg_index, enz_index=enz_index, descriptions=descriptions)
        for lig in candidate_ligands(entry):
            cand = EffectorCandidate(
                ligand=lig["name"], rank=sc["rank"], color=sc["color"], chebi=lig["chebi"],
                equation=lig["equation"], rhea=lig["rhea"],
                enzyme=g.name or g.protein_id, enzyme_product=g.product,
                distance=abs(enz_index - reg_index))
            prev = by_ligand.get(lig["name"])
            if prev is None or cand.rank > prev.rank:
                by_ligand[lig["name"]] = cand          # keep the best-supported occurrence per ligand
    return sorted(by_ligand.values(), key=lambda c: (-c.rank, c.distance, c.ligand))


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


def _from_published_db(seq, acc) -> list:
    """The PUBLISHED Ligify prediction for this exact protein, if it is one of their 3,164.

    Checked before re-deriving, because when the same protein is in their released database their
    answer IS the answer: we would be recomputing an approximation of a result they have already
    published, with the same formula, from the same kind of evidence. Returning theirs makes the
    provenance honest (`source="ligify_db"` rather than `"ligify"`) and lets a reader see agreement or
    disagreement with the published tool instead of a silently different number.

    Coverage is modest -- 13 of 140 candidates on the last production manifest -- and concentrated in
    the organic-ligand families where our own inducer call is weakest. Nothing fuzzy: a near neighbour
    is a different protein and gets the local derivation like everything else.
    """
    try:
        from predictor.effector import ligify_db as _ldb
    except Exception:                                   # noqa: BLE001
        return []
    rec = _ldb.lookup(seq=seq if isinstance(seq, str) else None, acc=acc)
    if rec is None or not rec.ligands:
        return []
    rank = rec.rank if rec.rank is not None else 0
    return [EffectorCandidate(
        ligand=lig, rank=rank, color=_rank_color(rank),
        chebi=[], equation=rec.equation, rhea=rec.rhea_id,
        enzyme=rec.annotation, enzyme_product="", distance=0)
        for lig in rec.ligands]


def _rank_color(rank: int) -> str:
    """Ligify's own colour bands, so a DB record and a locally derived one render identically."""
    return ("#02a602" if rank >= 70 else "#d4d400" if rank >= 50
            else "#d48302" if rank >= 30 else "#f50b02")


def predict_effector(seq, *, db=None, allow_ncbi: bool = False, max_kb: float = 10.0,
                     max_genes: int = 10, top: int = 10, accession: str | None = None,
                     use_published_db: bool = True) -> list:
    """TF sequence -> ranked candidate effectors. Offline against the mirror where possible; UniProt
    enzyme reactions are cached. Returns list[EffectorCandidate] (best rank first).

    Looks the protein up in the published Ligify database first (`use_published_db`, default on) and
    falls through to the local genome-context derivation when it is not there -- which is the great
    majority of queries. Pass `accession` when the caller knows it; otherwise the sequence hash is the
    join key."""
    if use_published_db:
        published = _from_published_db(seq, accession)
        if published:
            return published[:top]
    db = db or genome_db.default_db()
    ctx, tf = _resolve_context(seq, db=db, allow_ncbi=allow_ncbi)
    if ctx is None or tf is None:
        return []
    return effectors_from_context(ctx, tf, max_kb=max_kb, max_genes=max_genes)[:top]


# --------------------------------------------------------------------------- self-test
def _demo() -> None:
    # operon: divergent regulator (araC-like) next to two metabolic enzymes + a transporter
    genes = [
        Gene(100, 400, "-", "araC", protein_id="WP_araC", product="AraC family transcriptional regulator"),
        Gene(450, 1450, "+", "araA", protein_id="WP_araA", product="L-arabinose isomerase"),
        Gene(1470, 2400, "+", "araB", protein_id="WP_araB", product="ribulokinase"),
        Gene(2450, 3600, "+", "araD", protein_id="WP_araD", product="L-ribulose-5-phosphate epimerase"),
        Gene(8000, 8600, "+", "far", protein_id="WP_far", product="hypothetical protein"),  # out of window
    ]
    ctx = GenomeContext("SYN_ARA.1", "N" * 9000, genes)
    tf = genes[0]

    # inject enzyme entries: araA is the arabinose isomerase; araB/araD have no entry (non-enzyme here)
    entries = {"WP_araA": _synthetic_entry()}
    cands = effectors_from_context(ctx, tf, entry_fetcher=lambda pid: entries.get(pid), max_kb=10)
    print(f"{len(cands)} candidate effector(s):")
    for c in cands[:6]:
        print(f"  {c.ligand:<28} rank={c.rank} ({c.color})  via {c.enzyme} d={c.distance}  "
              f"[{c.rhea}]  {c.equation}")

    names = [c.ligand for c in cands]
    assert names, "no effectors predicted"
    assert "beta-L-arabinopyranose" in names, "the arabinose substrate should be a candidate"
    assert "ATP" not in names and "ADP" not in names, "currency leaked into effectors"
    top = cands[0]
    # araA is adjacent to araC (distance 1) in a window with one regulator -> high rank
    assert top.distance == 1 and top.rank >= 70, f"adjacent-enzyme rank too low: {top.rank}"
    assert "far" not in [c.enzyme for c in cands], "out-of-window gene wrongly used"
    print("OK: effector candidates ranked from operon context; currency + out-of-window excluded.")

    # live path is best-effort (needs a built genome DB + network); just confirm it returns a list
    try:
        out = predict_effector("MKKLTVSDLA", allow_ncbi=False)
        assert isinstance(out, list)
        print(f"OK: predict_effector() wired (returned {len(out)} candidates from local resolve).")
    except Exception as ex:
        print(f"SKIP (no local genome DB / resolve): {type(ex).__name__}: {ex}")


if __name__ == "__main__":
    _demo()
