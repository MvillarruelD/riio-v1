"""
homolog_regions.py -- collect homolog inter-operon (promoter) regions (Snowprint's comparative step).

For a TF sequence, the local genome DB already yields homolog loci across many genomes
(`genome_resolver.resolve` returns the source genome at high identity and homologs at lower identity).
This module extracts each homolog's **promoter-side window** (offline, via `blastdbcmd -entry`), so the
de-novo finder can score **phylogenetic conservation** of an operator candidate across homologs --
the signal that separates a real operator from a chance palindrome (and rescues Snowprint failure
mode #3: too few homologs is now visible as a short region list).

These region sequences feed `signals/motif_finder.predict(..., homolog_regions=...)`.

Run `python homolog_regions.py` for a self-test (collects AdhR homolog promoter regions).
"""
from __future__ import annotations

from pathlib import Path

from predictor.annotate import genome_resolver, genome_db
from predictor.annotate.genome_resolver import _blast, _run


def _extract(db: str, accession: str, start1: int, stop1: int) -> str:
    """blastdbcmd a 1-based inclusive sub-range of a genome accession; '' on failure."""
    try:
        res = _run([_blast("blastdbcmd"), "-db", db, "-entry", accession,
                    "-range", f"{max(1, start1)}-{stop1}", "-outfmt", "%s"])
    except Exception:
        return ""
    return "".join(res.stdout.split()).upper()


def regions_from_candidates(cands, db, *, upstream: int = genome_resolver.PROMOTER_UPSTREAM,
                            downstream: int = genome_resolver.PROMOTER_DOWNSTREAM) -> list:
    """Extract the promoter-side window of each resolved homolog candidate (strand-aware, LOCAL DB).
    Returned in transcriptional orientation, the shared convention (`genome_resolver.orient_promoter`)."""
    regions = []
    for c in cands:
        s1, e1 = genome_resolver.promoter_window_bounds(c.tf_start, c.tf_end, c.tf_strand,
                                                        upstream=upstream, downstream=downstream)
        seq = _extract(db, c.accession, s1, e1)
        if len(seq) >= 30:
            regions.append(genome_resolver.orient_promoter(seq, c.tf_strand))
    return regions


def regions_from_candidates_ncbi(cands, *, upstream: int = 350, downstream: int = 30) -> list:
    """ONLINE twin of `regions_from_candidates`: pull each locus's promoter window via NCBI EFetch
    (for IPG-derived candidates whose accession is not in the local mirror)."""
    regions = []
    for c in cands:
        seq = genome_resolver.fetch_upstream_region(c, upstream=upstream, downstream=downstream)
        if seq:
            regions.append(seq)
    return regions


def collect_homolog_regions(seq, *, n: int = 10, db=None, min_identity: float = 0.25,
                            upstream: int = 350, downstream: int = 30,
                            allow_ncbi: bool = False, min_regions: int = 4,
                            ncbi_top_n: int = 25) -> list:
    """sequence -> list of homolog promoter-region sequences (incl. the source genome).

    Default path is the LOCAL dereplicated mirror (offline, reproducible). When the local mirror yields
    too few homolog regions to make the conservation signal usable (`< min_regions` -- the *sparse*
    stratum from `homologs.py`) and `allow_ncbi=True`, augment with NCBI IPG/EFetch-derived regions
    (the online fallback). NCBI regions are de-duplicated against the local ones by sequence."""
    db = db or genome_db.default_db()
    regions = []
    if db is not None:
        cands = genome_resolver.resolve(seq, db=db, top_n=n, min_identity=min_identity)
        regions = regions_from_candidates(cands, db, upstream=upstream, downstream=downstream)

    if allow_ncbi and len(regions) < min_regions:
        try:
            nc = genome_resolver.resolve_ncbi(seq, top_n=ncbi_top_n, max_loci=ncbi_top_n)
            extra = regions_from_candidates_ncbi(nc, upstream=upstream, downstream=downstream)
            have = set(regions)
            regions += [r for r in extra if r not in have]
        except Exception:
            pass            # network/credentials unavailable -> keep the local (possibly sparse) set
    return regions


def _demo() -> None:
    import yaml
    root = Path(__file__).resolve().parents[2].parent / "5.8 Promoters" / "arsr_merr_families"
    seq = yaml.safe_load(open(root / "transcription_factors" / "merr" / "AdhR_Bsubtilis.yaml",
                              encoding="utf-8"))["sequence"]["protein"]
    if genome_db.default_db() is None:
        genome_db.build_local_db()
    regions = collect_homolog_regions(seq, n=10)
    print(f"collected {len(regions)} homolog promoter regions (lengths: "
          f"{[len(r) for r in regions]})")
    assert len(regions) >= 1, "no homolog regions collected"
    assert all(len(r) >= 30 for r in regions)

    # the conservation signal should now be usable by motif_finder
    from predictor.signals.motif_finder import predict
    cands = predict(regions[0], family="MerR", homolog_regions=regions, top_k=5)
    if cands:
        top = cands[0]
        print(f"top operator in source promoter: {top.kind} spacer={top.spacer} "
              f"conservation={top.conservation:.2f}  {top.seq}")
    print("OK: homolog promoter regions collected and usable for conservation scoring.")


if __name__ == "__main__":
    _demo()
