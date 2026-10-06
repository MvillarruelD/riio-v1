"""
auto_locus.py -- fully-automatic genome + TF-locus resolution from sequence alone (no hardcoded accessions).

The validated runs hardcode `mirror_acc` / `tf_locus` for the awkward cases (plasmid TFs, multi-contig
assemblies, strain-specific operons). Those hardcodes exist because the LOCAL genome BLAST DB either lacks
the genome or concatenates contigs (so coordinates don't line up). NCBI's Identical Protein Groups (IPG)
removes the need for all of them: given the TF protein, IPG returns every genome carrying that EXACT protein
with its authoritative (nucleotide accession, start, stop, strand) -- including plasmids (returns the plasmid
accession) and WGS contigs (returns the specific contig in ITS OWN coordinate frame, so no concatenation
mismatch). This module:

  1. remote blastp the query -> top protein accession(s)   (the protein's own NCBI record)
  2. IPG those accessions -> candidate genomic loci         (reuses genome_resolver._parse_ipg)
  3. rank: prefer a complete RefSeq chromosome/plasmid (NC_ / NZ_CP) over a WGS contig, then by being a
     plasmid when the protein is plasmid-borne -> return the best (accession, start, end, strand, organism)
     and the alternatives.

So `auto_locus(seq)` gives what `mirror_acc` + `tf_locus` were hardcoded for. Cached to results/auto_locus/.
Network + NCBI (slow remote blastp); offline self-test parses a canned IPG table.

  python -m predictor.annotate.auto_locus --tf CadC_Saureus_pI258   # demo on a plasmid TF
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

from predictor import resources

_REPO = Path(__file__).resolve().parents[2]
from predictor.annotate import genome_resolver as gr                          # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

CACHE = resources.cache_path("auto_locus")


def _acc_rank(accession: str) -> tuple:
    """Sort key (higher = better source genome). Prefer RefSeq complete records, then GenBank complete,
    then WGS contigs; plasmids rank with chromosomes (the operon may be plasmid-borne)."""
    acc = accession.upper()
    refseq = acc.startswith(("NC_", "NZ_CP", "AC_"))
    complete = bool(re.match(r"(NC_|NZ_CP|CP|AE|AL|AP|BA|U\d)", acc))
    wgs = bool(re.search(r"\d{6,}", acc)) and not complete
    return (refseq, complete, not wgs)


def top_protein_hits(seq: str, *, n: int = 4, hitlist: int = 25) -> list:
    """Remote blastp the query -> top protein accessions (best identity first). The top hit is usually the
    TF itself; we keep a few in case the exact protein is not the #1 record."""
    from io import StringIO
    from Bio.Blast import NCBIWWW, NCBIXML
    seq = gr._query_seq(seq)
    with NCBIWWW.qblast("blastp", "nr", seq, hitlist_size=hitlist) as h:
        rec = NCBIXML.read(StringIO(h.read()))
    accs = []
    for aln in rec.alignments[:n]:
        m = re.search(r"\b([A-Z]{1,3}_?\d{5,})\.?\d*\b", aln.accession or aln.hit_id or "")
        if m:
            accs.append(aln.accession or m.group(1))
    return accs


def auto_locus(seq: str, *, key: str | None = None, use_cache: bool = True, n_hits: int = 4,
               verbose: bool = True) -> dict | None:
    """Authoritative genome + locus for the TF from its sequence (via blastp -> IPG). Returns
    dict(accession, start, end, strand, organism, n_genomes, alternatives) or None. start/end are
    0-based half-open + strand coords (same frame as GenomeCandidate)."""
    if key and use_cache:
        cp = CACHE / f"{key}.json"
        if cp.exists():
            d = json.loads(cp.read_text(encoding="utf-8"))
            if verbose:
                print(f"  [cache] auto_locus {key}: {d['accession']}:{d['start']}-{d['end']} ({d['strand']})")
            return d
    try:
        accs = top_protein_hits(seq, n=n_hits)
    except Exception as e:
        if verbose:
            print(f"  auto_locus blastp failed: {type(e).__name__}: {e}")
        return None
    if verbose:
        print(f"  auto_locus: top protein hits = {accs}")
    cands = []
    Entrez = gr._entrez()
    for acc in accs:
        try:
            with Entrez.efetch(db="protein", id=acc, rettype="ipg", retmode="text") as h:
                data = h.read()
            txt = data.decode("utf-8", "replace") if isinstance(data, (bytes, bytearray)) else data
            cands += gr._parse_ipg(txt)[:300]              # cap per hit (a WP_ IPG can list 1000s of genomes)
        except Exception:
            continue
        # ACCUMULATE across ALL top hits rather than breaking on the first: the #1 blastp hit is often a
        # strain-specific GenBank protein whose IPG lists only ITS genome (e.g. a Shigella contig), which
        # would beat the canonical RefSeq chromosome. Gathering every hit's loci lets the RefSeq-complete
        # genome (NC_/NZ_CP) win the rank below. Stop early only once a RefSeq-complete candidate is in hand.
        if any(_acc_rank(c.accession)[0] for c in cands):
            break
    if not cands:
        if verbose:
            print("  auto_locus: no IPG loci")
        return None
    # dedupe loci, rank by source-genome quality
    seen, uniq = set(), []
    for c in cands:
        k = (c.accession, c.tf_start, c.tf_strand)
        if k not in seen:
            seen.add(k)
            uniq.append(c)
    uniq.sort(key=lambda c: _acc_rank(c.accession), reverse=True)
    best = uniq[0]
    out = {"accession": best.accession, "start": best.tf_start, "end": best.tf_end,
           "strand": best.tf_strand, "organism": best.organism, "n_genomes": len(uniq),
           "alternatives": [{"accession": c.accession, "start": c.tf_start, "end": c.tf_end,
                             "strand": c.tf_strand, "organism": c.organism} for c in uniq[1:6]]}
    if key and use_cache:
        CACHE.mkdir(parents=True, exist_ok=True)
        (CACHE / f"{key}.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    if verbose:
        print(f"  auto_locus: {out['accession']}:{out['start']}-{out['end']} ({out['strand']}) "
              f"{out['organism']}  [{out['n_genomes']} genomes carry this protein]")
    return out


def _self_test():
    # offline: ranking + IPG parse
    canned = ("Id\tSource\tNucleotide Accession\tStart\tStop\tStrand\tProtein\tProtein Name\tOrganism\n"
              "1\tRefSeq\tNC_013319.1\t10145\t10513\t-\tWP_x\tCadC\tStaphylococcus aureus\n"
              "1\tINSDC\tABC12345.1\t999\t1300\t+\tABC\tCadC\tWGS strain\n")
    loci = gr._parse_ipg(canned)
    assert loci and loci[0].accession == "NC_013319.1"
    loci.sort(key=lambda c: _acc_rank(c.accession), reverse=True)
    assert loci[0].accession == "NC_013319.1", "RefSeq complete plasmid must rank first"
    assert _acc_rank("NC_013319.1") > _acc_rank("ABUW01000005.1"), "RefSeq > WGS"
    print("OK: auto_locus ranking + IPG parse (offline).")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", default=None)
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args(argv)
    if a.self_test:
        _self_test(); return
    if a.fasta:
        seq = "".join(l.strip() for l in Path(a.fasta).read_text().splitlines() if not l.startswith(">"))
        res = auto_locus(seq, key=Path(a.fasta).stem)
    else:
        ap.error("pass --fasta or --self-test")
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
