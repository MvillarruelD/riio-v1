"""
protein_seqs.py -- fetch CDS protein sequences for RefSeq/GenBank accessions (WP_/NP_/YP_/AP_...), the
sequence supply that lets the ESM-C regulon-member annotator FIRE on real runs.

The offline GFF `GenomeContext` (annotate.context) carries each gene's `protein_id` and free-text
`product`, but NOT its translation. So a regulon member annotated "hypothetical protein" reaches the
ESM-C fallback (since removed) with no sequence to embed. This module closes
that gap: given the protein accessions of the members the substrate map could not place, it pulls their
sequences from NCBI (EFetch `db=protein`, batched), caches them to disk, and returns `{base_acc: seq}`.

Design (same discipline as effector.ligand's UniProt cache):
  * **cache-first** -- one `.faa` per accession under `$PREDICTOR_PROTSEQ_CACHE`
    (default `results/protein_seqs/`); a cached accession is never re-fetched.
  * **gated** -- the network fetch only runs with `allow_ncbi=True`; offline it returns whatever is cached
    (so a production run with no network simply embeds fewer members -- never errors).
  * **batched + version-agnostic** -- many ids per EFetch call; keyed by version-stripped accession so a
    GFF `WP_x` query resolves the EFetch `WP_x.1` record.

Run `python -m predictor.annotate.protein_seqs WP_000148532.1 P0A9...   (live)
    python -m predictor.annotate.protein_seqs --self-test                (offline cache logic)`.
"""
from __future__ import annotations

import argparse
import os
import re
import time
from pathlib import Path

from predictor import resources

_REPO = Path(__file__).resolve().parents[2]
_CACHE = Path(os.environ.get("PREDICTOR_PROTSEQ_CACHE", resources.cache_path("protein_seqs")))


def _base(acc: str) -> str:
    """Version-stripped accession used as the stable cache + lookup key (WP_123.1 -> WP_123)."""
    return (acc or "").strip().split(".")[0]


def _cache_path(acc: str) -> Path:
    return _CACHE / f"{re.sub(r'[^A-Za-z0-9_.-]', '_', _base(acc))}.faa"


def _read_cache(acc: str) -> str | None:
    p = _cache_path(acc)
    if not p.exists():
        return None
    try:
        txt = p.read_text(encoding="ascii", errors="ignore")
    except Exception:
        return None
    seq = "".join(l.strip() for l in txt.splitlines() if l and not l.startswith(">"))
    return seq or ""          # "" = cached-but-empty (a confirmed miss: don't refetch)


def _write_cache(acc: str, seq: str) -> None:
    _CACHE.mkdir(parents=True, exist_ok=True)
    _cache_path(acc).write_text(f">{_base(acc)}\n{seq}\n", encoding="ascii", errors="ignore")


def _parse_fasta(text: str) -> dict[str, str]:
    """Multi-FASTA text -> {base_accession: sequence}. The header's first token is the accession."""
    out, acc, buf = {}, None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if acc is not None:
                out[_base(acc)] = "".join(buf)
            acc = line[1:].split()[0] if len(line) > 1 else None
            buf = []
        elif acc is not None:
            buf.append(line.strip())
    if acc is not None:
        out[_base(acc)] = "".join(buf)
    return out


def fetch_protein_sequences(protein_ids, *, allow_ncbi: bool = False, batch: int = 50,
                            verbose: bool = False) -> dict[str, str]:
    """`{base_acc: seq}` for the given protein accessions. Cache-first; only the still-missing accessions
    are fetched (and only when `allow_ncbi`). Empty/failed accessions are cached as '' so they are not
    retried. Non-protein-looking ids (no WP_/NP_/.. and not a UniProt acc) are skipped."""
    ids = [i for i in {(_base(p)) for p in (protein_ids or []) if p} if i]
    out: dict[str, str] = {}
    missing: list[str] = []
    for acc in ids:
        cached = _read_cache(acc)
        if cached is None:
            missing.append(acc)
        elif cached:
            out[acc] = cached
    if missing and allow_ncbi:
        try:
            try:
                from predictor.annotate import genome_resolver as gr
            except ImportError:
                from predictor.annotate import genome_resolver as gr
            Entrez = gr._entrez()
        except Exception as e:
            if verbose:
                print(f"  EFetch unavailable ({type(e).__name__}: {e}); cache-only")
            return out
        for i in range(0, len(missing), batch):
            chunk = missing[i:i + batch]
            try:
                with Entrez.efetch(db="protein", id=",".join(chunk), rettype="fasta",
                                   retmode="text") as h:
                    parsed = _parse_fasta(h.read())
            except Exception as e:
                if verbose:
                    print(f"  EFetch failed for {len(chunk)} id(s) ({type(e).__name__}: {e})")
                parsed = {}
            for acc in chunk:
                seq = parsed.get(acc, "")
                _write_cache(acc, seq)               # cache hits AND confirmed misses
                if seq:
                    out[acc] = seq
            time.sleep(0.12)                          # be polite to NCBI (no key) -- ~8 req/s ceiling
        if verbose:
            print(f"  fetched {sum(1 for a in missing if out.get(a))}/{len(missing)} protein sequence(s)")
    return out


def make_fetcher(*, allow_ncbi: bool = False, verbose: bool = False):
    """A 1-arg `seq_fetcher(ids) -> {base_acc: seq}` closure for injection into
    `effector.regulon_enrichment.enrich_genes` / `effector.inducer.augment_with_regulon`."""
    def _f(ids):
        return fetch_protein_sequences(ids, allow_ncbi=allow_ncbi, verbose=verbose)
    return _f


# --------------------------------------------------------------------------- self-test
def _self_test() -> None:
    import tempfile
    global _CACHE
    saved = _CACHE
    _CACHE = Path(tempfile.mkdtemp(prefix="protseq_test_"))
    try:
        # cache logic: pre-seed one accession, confirm cache-first returns it with NO network
        _write_cache("WP_TEST1.1", "MKKLTVSDLA")
        got = fetch_protein_sequences(["WP_TEST1.1"], allow_ncbi=False)
        assert got == {"WP_TEST1": "MKKLTVSDLA"}, f"cache-first failed: {got}"
        # version-agnostic lookup: querying the base accession hits the same cache
        assert fetch_protein_sequences(["WP_TEST1"], allow_ncbi=False) == {"WP_TEST1": "MKKLTVSDLA"}
        # offline miss returns {} (no exception) and does NOT create a bogus cache entry
        assert fetch_protein_sequences(["WP_NOPE9"], allow_ncbi=False) == {}
        assert not _cache_path("WP_NOPE9").exists(), "offline miss must not write a cache file"
        # confirmed-empty cache is treated as a miss, not refetched
        _write_cache("WP_EMPTY1", "")
        assert fetch_protein_sequences(["WP_EMPTY1"], allow_ncbi=False) == {}
        # multi-fasta parser
        p = _parse_fasta(">WP_A.1 zinc transporter\nMKKL\nTVS\n>WP_B.2 hypothetical\nMQRS\n")
        assert p == {"WP_A": "MKKLTVS", "WP_B": "MQRS"}, f"fasta parse: {p}"
        print("  cache-first OK; version-agnostic OK; offline miss {} OK; empty-cache miss OK; fasta-parse OK")
        print("OK: protein-sequence fetcher is cache-first, gated, batched, version-agnostic.")
    finally:
        import shutil
        shutil.rmtree(_CACHE, ignore_errors=True)
        _CACHE = saved


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ids", nargs="*", help="protein accessions to fetch (live)")
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--offline", action="store_true", help="cache-only (no EFetch)")
    a = ap.parse_args(argv)
    if a.self_test or not a.ids:
        _self_test()
        return
    res = fetch_protein_sequences(a.ids, allow_ncbi=not a.offline, verbose=True)
    for acc, seq in res.items():
        print(f"{acc}\t{len(seq)} aa\t{seq[:40]}{'...' if len(seq) > 40 else ''}")


if __name__ == "__main__":
    main()
