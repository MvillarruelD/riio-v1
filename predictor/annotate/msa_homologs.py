"""
msa_homologs.py -- turn a deep MSA (AF3 server, or ColabFold without it) into homolog PROMOTER regions.

The conservation-first discovery signal is only as good as the homolog set behind it. Today it draws on a
~92-sequence LOCAL genome BLAST DB; a deep MSA carries ~10^4 homologs (CzrA: 15,627) with UniRef/UniProt
IDs + taxonomy. This module bridges the gap the project lacked: MSA homolog proteins -> their genomes ->
their promoter windows, which feed `signals/motif_finder.predict(..., homolog_regions=...)` exactly like the
local/NCBI regions do (same + strand sub-range convention, directly comparable).

Pipeline: parse a3m -> CLOSE orthologs only (>=~45% id; operator conservation tracks close orthology, distant
paralogs carry different operators -- the af3_msa insight) -> UniProt REST maps each accession to a RefSeq/
EMBL protein -> NCBI IPG gives the genomic locus/loci of that protein (reuses genome_resolver._parse_ipg) ->
EFetch the promoter window (reuses genome_resolver.fetch_upstream_region). Everything is cached to
`results/msa_homologs/<key>.json` so a TF is mapped once.

MSA source: the AF3 server already wrote `results/af3_decision/*_op_af3decision/msas/*unpaired*.a3m` (reused
via `af3_msa._op_msa`). For a NOVEL TF with no AF3 fold, `colabfold_a3m(seq)` fetches the SAME MMseqs2 MSA
from the public ColabFold server (https://api.colabfold.com) -- the exact pipeline the AlphaFold2 ColabFold
notebook uses, no AF3 account, no folding.

  python -m predictor.annotate.msa_homologs --a3m <file.a3m> --key TF     # collect + cache regions
  python -m predictor.annotate.msa_homologs --colabfold tf.fasta --key TF # fetch MSA then collect
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

from predictor import resources

_REPO = Path(__file__).resolve().parents[2]
from predictor.annotate.af3_msa import parse_a3m, _op_msa                       # noqa: E402
from predictor.annotate import genome_resolver as gr                            # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

CACHE = resources.cache_path("msa_homologs")
_UA = {"User-Agent": "predictor-msa-homologs/1.0"}

#: Tag for the promoter-window convention baked into every cached region set. Regions are cached as raw
#: sequence, so a change in HOW they are extracted silently invalidates every cache without changing its
#: key. Bump this whenever `genome_resolver.orient_promoter` / `promoter_window_bounds` change meaning.
#: "t1" = transcriptional orientation (minus-strand windows reverse-complemented).
_REGION_CONVENTION = "t1"


def _get(url, timeout=60):
    req = urllib.request.Request(url, headers=_UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


#: Below this many MAPPED orthologs the pooled logo rests on too few sequences to mean much; far
#: above it, low-identity noise starts to outweigh the added information. The floor's job is not to
#: force a retry -- with the UniRef member a 6000-sequence a3m clears it easily -- but to turn a
#: silent zero into a stated one, recorded as `shortfall` so a thin conservation axis is visible in
#: the bundle rather than inferable only from a log line.
MIN_MAPPABLE_ORTHOLOGS = 20

#: The UniProt accession format (UniProtKB 'Accession number' spec). A whitelist, because the
#: alternative -- splitting any header on its first underscore -- silently invents accessions.
_UNIPROT_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$")


def _acc_of(uniref: str) -> str:
    """UniRef90_G5JME4 -> G5JME4 ; bare UniProt accession -> itself ; anything else -> ''.

    Returning '' for an unrecognised header is the point. ColabFold's environmental members carry
    SRA-style headers like `SRR6516225_2681333`; the previous split-on-first-underscore returned
    `2681333`, the sequence INDEX, which is not an accession of any kind and which UniProt will
    never resolve. Verified live on one candidate: 0 of 12 such "accessions" returned.

    That alone did not cause the empty conservation axis -- a correct parse yields `SRR6516225`,
    an SRA run with no UniProt entry and no assembled genome, so still unmappable (the real fix is
    asking ColabFold for the UniRef member; see expand_homolog_regions). What it did cause was
    `close_orthologs` reporting orthologs the pipeline could never use, so a systematic source
    problem read in the logs as a per-candidate biological one.
    """
    head = uniref.split()[0] if uniref else ""
    if head.startswith("UniRef"):
        head = head.split("_", 1)[1] if "_" in head else ""
    return head if _UNIPROT_RE.match(head) else ""


def close_orthologs(a3m_path, *, min_id: float = 0.45, max_id: float = 0.97, cap: int = 80) -> list:
    """Close orthologs (mappable UniProt accessions) from an a3m, dedup by accession, best-first.
    Skips UniParc (UPI*) entries -- they rarely carry a clean genome cross-reference."""
    seen, out = set(), []
    for r in parse_a3m(a3m_path, max_seqs=20000):
        if r["identity"] >= 1.0:
            continue
        if not (min_id <= r["identity"] <= max_id):
            continue
        acc = _acc_of(r["uniref"])
        if not acc or acc.startswith("UPI") or acc in seen:
            continue
        seen.add(acc)
        out.append({"acc": acc, "identity": round(r["identity"], 3), "tax": r["taxname"]})
        if len(out) >= cap:
            break
    out.sort(key=lambda r: -r["identity"])
    return out


def uniprot_to_protein(accs: list) -> dict:
    """Batch UniProt REST: {uniprot_acc -> a RefSeq (preferred) or EMBL protein accession | None}.
    RefSeq WP_* proteins are preferred because their IPG maps to many genomes (more homolog regions)."""
    out = {}
    for i in range(0, len(accs), 80):
        chunk = accs[i:i + 80]
        url = ("https://rest.uniprot.org/uniprotkb/accessions?accessions=" + ",".join(chunk) +
               "&fields=accession,xref_refseq,xref_embl&format=json")
        try:
            d = json.loads(_get(url).decode("utf-8", "replace"))
        except Exception:
            continue
        for e in d.get("results", []):
            acc = e.get("primaryAccession")
            xr = e.get("uniProtKBCrossReferences", [])
            refseq = [x["id"] for x in xr if x.get("database") == "RefSeq"]
            embl = []
            for x in xr:
                if x.get("database") == "EMBL":
                    props = {p["key"]: p["value"] for p in x.get("properties", [])}
                    pid = props.get("ProteinId")
                    if pid and pid != "-":
                        embl.append(pid)
            out[acc] = refseq[0] if refseq else (embl[0] if embl else None)
        time.sleep(0.2)
    return out


def _read_text(h) -> str:
    """efetch handles return bytes for some rettypes (e.g. 'ipg') and str for others -> always str."""
    data = h.read()
    return data.decode("utf-8", "replace") if isinstance(data, (bytes, bytearray)) else data


def loci_from_protein(prot_id: str, *, cap: int = 3) -> list:
    """NCBI IPG -> genomic loci of an identical protein (reuses genome_resolver._parse_ipg).
    The IPG efetch handle yields BYTES; decode before parsing (the shared _parse_ipg splits on a str)."""
    Entrez = gr._entrez()
    try:
        with Entrez.efetch(db="protein", id=prot_id, rettype="ipg", retmode="text") as h:
            return gr._parse_ipg(_read_text(h))[:cap]
    except Exception:
        return []


def _upstream(cand, *, upstream: int = gr.PROMOTER_UPSTREAM,
              downstream: int = gr.PROMOTER_DOWNSTREAM) -> str:
    """Promoter-side window of a locus, byte-safe, in transcriptional orientation -- the shared
    convention of the local/NCBI paths (`genome_resolver.orient_promoter`)."""
    s1, e1 = gr.promoter_window_bounds(cand.tf_start, cand.tf_end, cand.tf_strand,
                                       upstream=upstream, downstream=downstream)
    if e1 <= s1:
        return ""
    try:
        Entrez = gr._entrez()
        with Entrez.efetch(db="nuccore", id=cand.accession, rettype="fasta", retmode="text",
                           seq_start=str(s1), seq_stop=str(e1)) as h:
            fa = _read_text(h)
    except Exception:
        return ""
    seq = "".join(l.strip() for l in fa.splitlines() if not l.startswith(">")).upper()
    if len(seq) < 30:
        return ""
    return gr.orient_promoter(seq, cand.tf_strand)


def collect_msa_regions(a3m_path, *, key: str | None = None, min_id: float = 0.45, max_id: float = 0.97,
                        cap: int = 80, max_regions: int = 60, upstream: int = 350, downstream: int = 30,
                        loci_per_protein: int = 2,
                        use_cache: bool = True, verbose: bool = True) -> dict:
    """a3m -> homolog promoter regions (+ provenance), cached. Returns dict(regions, n_orthologs,
    n_mapped, n_loci, sources).

    The cache key names the SOURCE ALIGNMENT, not just the TF. It used to be
    `key or Path(a3m_path).stem`, and every caller passes `key=tf_id`, so the a3m never reached the
    key at all -- two different alignments for one TF shared one cache entry. That made the UniRef
    switch inert on any warm machine: the new a3m was fetched, then its regions were served from the
    environmental run's cached zero. Observed 2026-09-04 as
    `[cache] MaviHomi__GntR__WP_003873121.1_t1: 0 MSA regions` on a run whose own log line above it
    read `ColabFold MSA -> ..._colabfold_uniref.a3m (7063 sequences)`.
    """
    # The a3m stem usually already begins with the TF id, so concatenating both repeated it
    # ("Mavi__GntR__WP_x__Mavi__GntR__WP_x_colabfold_uniref_t1"). Prefer the stem, and only prepend
    # `key` when the stem does not already carry it -- the invariant that matters is that the key
    # names the ALIGNMENT, not merely the TF.
    _stem = Path(a3m_path).stem
    _base = _stem if (not key or _stem.startswith(key)) else f"{key}__{_stem}"
    ck = f"{_base}_{_REGION_CONVENTION}"
    cpath = CACHE / f"{ck}.json"
    if use_cache and cpath.exists():
        d = json.loads(cpath.read_text(encoding="utf-8"))
        if verbose:
            print(f"  [cache] {ck}: {len(d['regions'])} MSA regions")
        return d

    orth = close_orthologs(a3m_path, min_id=min_id, max_id=max_id, cap=cap)
    mapping = uniprot_to_protein([o["acc"] for o in orth]) if orth else {}
    regions, sources, loci_seen = [], [], set()
    for o in orth:
        pid = mapping.get(o["acc"])
        if not pid:
            continue
        for c in loci_from_protein(pid, cap=loci_per_protein):
            lk = (c.accession, c.tf_start, c.tf_strand)
            if lk in loci_seen:
                continue
            loci_seen.add(lk)
            seq = _upstream(c, upstream=upstream, downstream=downstream)
            if seq:
                regions.append(seq)
                sources.append({"uniprot": o["acc"], "protein": pid, "nuc": c.accession,
                                "identity": o["identity"], "organism": c.organism})
            if len(regions) >= max_regions:
                break
        if len(regions) >= max_regions:
            break
    n_mapped = sum(bool(v) for v in mapping.values())
    out = {"key": ck, "regions": regions, "n_orthologs": len(orth), "n_mapped": n_mapped,
           "n_loci": len(loci_seen), "sources": sources}
    if n_mapped < MIN_MAPPABLE_ORTHOLOGS:
        out["shortfall"] = {
            "n_mappable": len(orth), "n_mapped": n_mapped, "floor": MIN_MAPPABLE_ORTHOLOGS,
            "reason": ("no mappable accessions in the a3m -- is this the uniref member?"
                       if not orth else
                       "UniProt returned no RefSeq/EMBL cross-reference for these accessions"),
        }
    if use_cache:
        CACHE.mkdir(parents=True, exist_ok=True)
        cpath.write_text(json.dumps(out, indent=2), encoding="utf-8")
    if verbose:
        # "mappable orthologs", not "close orthologs": since _acc_of whitelists UniProt accessions
        # this count is what the pipeline can actually follow up, not what the a3m happens to hold.
        print(f"  {ck}: {len(orth)} mappable orthologs -> {n_mapped} mapped -> "
              f"{len(regions)} promoter regions")
        if "shortfall" in out:
            print(f"  WARNING: {n_mapped} mapped < floor {MIN_MAPPABLE_ORTHOLOGS} -- the "
                  f"conservation axis will be weak or empty ({out['shortfall']['reason']})")
    return out


# --------------------------------------------------------------------------- SSN-cluster membership -> promoters (WS2)
def _cluster_by_id(cluster_id):
    """Locate an SSN cluster from its id, across ALL twelve SSN families.

    This used to read `fam = "MerR" if id.startswith("MerR") else "ArsR/SmtB"` -- written when the SSN had
    only those two families. After the 2026 drop extended it to twelve, every other cluster (Fur_ecFur,
    Rrf2_ecIscR, CsoR_ecRcnR, NikR_c1, DtxR_ecMntR, LysR_c2, ...) was looked up in the ArsR/SmtB list,
    never found, and returned None -- silently disabling the SSN-cluster homolog route for ten of the
    twelve families. A cluster id is `<family_tag>_<...>`, so resolve the tag properly."""
    from predictor.annotate import ssn_clusters as _ssn
    cid = str(cluster_id)
    tag = cid.split("_", 1)[0]
    fams = [f for f in _ssn.FAMILIES if _ssn.family_tag(f) == tag] or list(_ssn.FAMILIES)
    for fam in fams:
        try:
            for c in _ssn.load_clusters(fam):
                if c.cluster_id == cid:
                    return c
        except Exception:
            continue
    return None


def _seq_identity(a: str, b: str) -> float:
    from Bio.Align import PairwiseAligner, substitution_matrices
    if not a or not b:
        return 0.0
    al = PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score = -11, -1
    try:
        aln = al.align(a, b)[0]
    except Exception:
        return 0.0
    m = sum(a[qb + k] == b[tb + k] for (qb, qe), (tb, te) in zip(*aln.aligned) for k in range(qe - qb))
    return m / max(1, min(len(a), len(b)))


def _sample_members(cluster_info, query_seq, *, n: int, sampling: str, scan_cap: int = 200) -> list:
    """Pick `n` member accessions from a cluster by `sampling`:
       top_similar    -> rank scanned members by identity to `query_seq`, best n (needs query_seq);
       representative -> evenly spaced across the scanned window (DIVERSITY proxy, the footprint-friendly one);
       all            -> first n in file order (baseline)."""
    from predictor.annotate import ssn_clusters as _ssn
    members = list(_ssn.iter_sequences(cluster_info, limit=scan_cap))   # [(acc, seq), ...]
    if not members:
        return []
    if sampling == "top_similar" and query_seq:
        scored = sorted(((acc, _seq_identity(query_seq, ms)) for acc, ms in members), key=lambda x: -x[1])
        return [a for a, _ in scored[:n]]
    accs = [a for a, _ in members]
    if sampling == "representative" and n < len(accs):
        idx = sorted(set(round(i * (len(accs) - 1) / (n - 1)) for i in range(n))) if n > 1 else [0]
        return [accs[i] for i in idx]
    return accs[:n]


def _mapped_members(accs: list, want: int, *, verbose: bool = False) -> list:
    """Keep the members that actually resolve to a protein accession, in their given order, up to `want`.

    Most cluster members are unreviewed TrEMBL entries with no usable genome cross-reference, so sampling
    exactly `want` and mapping them wastes most of the budget -- measured, 25 sampled members yield only
    ~9 mapped. Mapping is a BATCHED UniProt call (80 accessions/request), so testing several times more
    candidates costs one or two extra requests and turns the sample into `want` members that genuinely map.
    Order is preserved, so the caller's sampling strategy still decides WHICH members are kept.

    Members are sampled 1:1 with `want`, never oversampled. Oversampling to raise the region yield was
    measured and rejected: at 6x the region count nearly triples (Fur 17 -> 50, IscR 18 -> 50) while
    recovered true sites FALL 23 -> 15 and top-decile precision moves 0.141 -> 0.151 (about one site,
    inside the noise). This mirrors the depth sweep. Region COUNT is not the binding constraint; region
    HOMOGENEITY is -- extra promoters that are not really homologous operators average into a blurrier
    PWM. A future attempt should filter regions for homogeneity, not collect more of them."""
    mapping = uniprot_to_protein(accs) if accs else {}
    keep = [a for a in accs if mapping.get(a)][:want]
    if verbose:
        print(f"  member selection: {len(accs)} candidates -> {sum(1 for a in accs if mapping.get(a))} "
              f"map to a protein -> keeping {len(keep)}")
    return keep, mapping


def cluster_promoter_regions(query_seq: str | None, ssn_cluster: str, *,
                             n_cluster_genomes: int = 25, cluster_sampling: str = "representative",
                             upstream: int = 350, downstream: int = 30, loci_per_protein: int = 2,
                             max_regions: int = 60, scan_cap: int = 200, use_cache: bool = True,
                             verbose: bool = True) -> dict:
    """Homolog PROMOTER regions sourced DIRECTLY from the TF's SSN-cluster members (the isofunctional set),
    bypassing the query MSA -- the "use the sequences from the SSN, find their genomes, run the same
    autoregulatory-promoter palindrome search" path. Samples `n_cluster_genomes` members by
    `cluster_sampling` -> UniProt -> IPG locus -> autoregulatory promoter window (same + strand convention
    as the local/NCBI/MSA paths). Feeds the motif finder's `homolog_regions`.

    The pipeline always calls this with the defaults; depth 25 and representative sampling are the measured
    optimum (deeper blurs the PWM -- see the depth sweep). The parameters exist for the offline sweep
    scripts that re-derive that result, not as a runtime choice.
    Online (UniProt REST + NCBI IPG/EFetch); cached per (cluster, sampling, n)."""
    info = _cluster_by_id(ssn_cluster)
    if info is None:
        if verbose:
            print(f"  cluster_promoter_regions: unknown cluster {ssn_cluster}")
        return {"regions": [], "cluster": ssn_cluster, "n_members": 0, "n_mapped": 0, "n_loci": 0, "sources": []}
    ck = f"cluster_{ssn_cluster}_{cluster_sampling}_n{n_cluster_genomes}_{_REGION_CONVENTION}"
    cpath = CACHE / f"{ck}.json"
    if use_cache and cpath.exists():
        d = json.loads(cpath.read_text(encoding="utf-8"))
        if verbose:
            print(f"  [cache] {ck}: {len(d['regions'])} cluster-promoter regions")
        return d
    cands = _sample_members(info, query_seq, n=n_cluster_genomes,
                            sampling=cluster_sampling, scan_cap=scan_cap)
    accs, mapping = _mapped_members(cands, n_cluster_genomes, verbose=verbose)
    regions, sources, loci_seen = [], [], set()
    for acc in accs:
        pid = mapping.get(acc)
        if not pid:
            continue
        for c in loci_from_protein(pid, cap=loci_per_protein):
            lk = (c.accession, c.tf_start, c.tf_strand)
            if lk in loci_seen:
                continue
            loci_seen.add(lk)
            s = _upstream(c, upstream=upstream, downstream=downstream)
            if s:
                regions.append(s)
                sources.append({"uniprot": acc, "protein": pid, "nuc": c.accession, "organism": c.organism})
            if len(regions) >= max_regions:
                break
        if len(regions) >= max_regions:
            break
    out = {"key": ck, "cluster": ssn_cluster, "sampling": cluster_sampling, "n_requested": n_cluster_genomes,
           "regions": regions, "n_members": len(accs), "n_mapped": sum(bool(v) for v in mapping.values()),
           "n_loci": len(loci_seen), "sources": sources}
    if use_cache:
        CACHE.mkdir(parents=True, exist_ok=True)
        cpath.write_text(json.dumps(out, indent=2), encoding="utf-8")
    if verbose:
        print(f"  {ck}: {len(accs)} members -> {out['n_mapped']} mapped -> {len(regions)} promoter regions")
    return out


# --------------------------------------------------------------------------- ColabFold MSA (no AF3 server)
def colabfold_a3m(seq: str, out_path, *, mode: str = "env", timeout: int = 1200, verbose: bool = True,
                  prefer: str | None = None) -> Path | None:
    """Fetch the MMseqs2 MSA for `seq` from the public ColabFold server -> write an a3m. This is the exact
    MSA the AlphaFold2 ColabFold notebook builds (UniRef30 + environmental), with NO AF3 account and NO
    folding. Submit -> poll -> download the result tar -> extract the largest a3m. None on failure.

    `prefer` selects a specific member of the result tar by filename substring instead of taking the
    largest -- `prefer="uniref"` picks `uniref.a3m`, which is what `structure.metalnet` needs because
    MetalNet2 was trained and benchmarked on the UniRef alignment specifically (`model/train/cmd.sh`:
    `msa_type=uni`), not on the merged environmental one. Falls back to the largest if absent."""
    out_path = Path(out_path)
    base = "https://api.colabfold.com"
    # the MMseqs2 server expects a FASTA-formatted query (`>101\n<seq>`), not a bare sequence
    q = f">101\n{seq.strip().replace(chr(10), '')}\n"
    data = urllib.parse.urlencode({"q": q, "mode": mode}).encode()
    try:
        # The public server rate-limits ticket submission with HTTP 429. A single sequence rarely trips
        # it; a batch (structure.metalnet over a 150-candidate run) trips it repeatedly, and an
        # unretried 429 silently costs that TF its alignment -- which reads downstream as "no opinion"
        # rather than as a transient network failure. Back off and retry instead.
        req = urllib.request.Request(base + "/ticket/msa", data=data, headers=_UA)
        tk = None
        for attempt in range(5):
            try:
                with urllib.request.urlopen(req, timeout=120) as r:
                    tk = json.loads(r.read().decode())
                break
            except urllib.error.HTTPError as e:
                if e.code != 429 or attempt == 4:
                    raise
                wait = 20 * (attempt + 1)
                if verbose:
                    print(f"  ColabFold rate-limited (429); retrying in {wait}s "
                          f"[{attempt + 1}/4]")
                time.sleep(wait)
        tid = (tk or {}).get("id")
        if not tid:
            if verbose:
                print(f"  ColabFold: no ticket id ({tk})")
            return None
        if verbose:
            print(f"  ColabFold ticket {tid} ({tk.get('status')}); polling...")
        t0 = time.time()
        while time.time() - t0 < timeout:
            with urllib.request.urlopen(base + f"/ticket/{tid}", timeout=60) as r:
                st = json.loads(r.read().decode()).get("status")
            if st == "COMPLETE":
                break
            if st in ("ERROR", "UNKNOWN", "MAINTENANCE"):
                if verbose:
                    print(f"  ColabFold status={st}")
                return None
            time.sleep(8)
        else:
            if verbose:
                print("  ColabFold timed out")
            return None
        raw = _get(base + f"/result/download/{tid}", timeout=300)
        import io
        import tarfile
        out_path.parent.mkdir(parents=True, exist_ok=True)
        best, best_n, wanted = None, -1, None
        with tarfile.open(fileobj=io.BytesIO(raw), mode="r:gz") as tf:
            for m in tf.getmembers():
                if m.name.endswith(".a3m"):
                    txt = tf.extractfile(m).read().decode("utf-8", "replace")
                    n = txt.count(">")
                    if prefer and prefer in Path(m.name).name:
                        wanted = txt
                    if n > best_n:
                        best, best_n = txt, n
        best = wanted if wanted is not None else best
        if best is None:
            return None
        # the server can emit NUL bytes inside the alignment; they are never meaningful and, left in,
        # they shift the columns of every consumer that reads the rows positionally
        out_path.write_text(best.replace("\x00", ""), encoding="utf-8")
        if verbose:
            print(f"  ColabFold MSA -> {out_path} ({best_n} sequences)")
        return out_path
    except Exception as e:
        if verbose:
            print(f"  ColabFold failed: {type(e).__name__}: {e}")
        return None


def expand_homolog_regions(local_regions: list, *, tf_id: str | None = None, seq: str | None = None,
                           min_local: int = 8, allow_colabfold: bool = False, max_regions: int = 60,
                           verbose: bool = True) -> tuple:
    """Augment a LOCAL homolog-region set with MSA-derived regions ONLY when local is sparse (the
    Snowprint failure-mode-#3 regime, where the experiment showed a real motif gain -- NmtR 4->64 regions,
    IC 0.845->0.909). When local is already rich (>= min_local), returns it unchanged (no network).

    Returns (regions, info) where info records whether MSA was used and how many regions it added. MSA
    source: an AF3 op a3m on disk (free), else a ColabFold MSA if `allow_colabfold` (the AF3-free route)."""
    info = {"local": len(local_regions), "msa_added": 0, "source": None,
            "alignment_path": None, "source_record": {}}
    if len(local_regions) >= min_local:
        return local_regions, info
    a3m = _op_msa(tf_id) if tf_id else None
    info["source"] = "af3_msa" if a3m else None
    if a3m is None and seq and allow_colabfold:
        # prefer="uniref", and a DIFFERENT cache filename so an environmental a3m left on disk by an
        # earlier run is not silently reused. The ColabFold tar carries both a UniRef-only alignment
        # and a merged UniRef+environmental one; taking the largest member took the merged one, which
        # for one measured candidate was 2019 SRR + 266 ERR + 24 UniRef. SRA metagenome sequences
        # have no UniProt entry and no assembled genome, so no promoter can be fetched from them:
        # 59 of the 86 candidates that reached this step in run A mapped 0 orthologs.
        # structure/metalnet.py has always asked for the UniRef member; this path needs it for a
        # different reason -- mappable accessions, not a trained-on alignment -- and never said so.
        a3m = colabfold_a3m(seq, CACHE / f"{tf_id or 'query'}_colabfold_uniref.a3m",
                            prefer="uniref", verbose=verbose)
        info["source"] = "colabfold" if a3m else None
    if a3m is None:
        if verbose:
            print(f"  MSA expansion: local sparse ({len(local_regions)}) but no MSA source "
                  f"(af3 a3m absent; allow_colabfold={allow_colabfold})")
        return local_regions, info
    info["alignment_path"] = str(a3m)
    try:
        msa_record = collect_msa_regions(a3m, key=tf_id, max_regions=max_regions,
                                         verbose=verbose)
        msa_regions = msa_record["regions"]
        info["source_record"] = {k: v for k, v in msa_record.items() if k != "regions"}
        info["source_records"] = [
            {"sequence": region, **(source if isinstance(source, dict) else {})}
            for region, source in zip(msa_regions, msa_record.get("sources") or [])
        ]
    except Exception as e:
        if verbose:
            print(f"  MSA expansion failed: {type(e).__name__}: {e}")
        return local_regions, info
    have = set(local_regions)
    added = [r for r in msa_regions if r not in have]
    info["msa_added"] = len(added)
    if verbose:
        print(f"  MSA expansion: local {len(local_regions)} + {len(added)} MSA regions "
              f"(source={info['source']})")
    return local_regions + added, info


def regions_for_tf(tf_id: str, *, seq: str | None = None, allow_colabfold: bool = False, **kw) -> list:
    """High-level: homolog promoter regions for a TF from its MSA. Prefers an AF3 op MSA on disk; else
    (novel TF) fetches a ColabFold MSA if `allow_colabfold`. Returns a list of region sequences ('' filtered)."""
    a3m = _op_msa(tf_id)
    if a3m is None and seq and allow_colabfold:
        a3m = colabfold_a3m(seq, CACHE / f"{tf_id}_colabfold_uniref.a3m", prefer="uniref")
    if a3m is None:
        return []
    return collect_msa_regions(a3m, key=tf_id, **kw)["regions"]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--a3m", default=None)
    ap.add_argument("--colabfold", default=None, help="FASTA -> fetch ColabFold MSA first")
    ap.add_argument("--key", default=None)
    ap.add_argument("--no-cache", action="store_true")
    a = ap.parse_args(argv)
    a3m = a.a3m
    if a.colabfold:
        seq = "".join(l.strip() for l in Path(a.colabfold).read_text().splitlines() if not l.startswith(">"))
        a3m = colabfold_a3m(seq, CACHE / f"{a.key or 'query'}_colabfold_uniref.a3m",
                            prefer="uniref")
    if not a3m:
        ap.error("pass --a3m or --colabfold")
    out = collect_msa_regions(a3m, key=a.key, use_cache=not a.no_cache)
    print(json.dumps({k: v for k, v in out.items() if k != "regions"}, indent=2))
    print(f"region lengths: {[len(r) for r in out['regions'][:20]]}{' ...' if len(out['regions'])>20 else ''}")


if __name__ == "__main__":
    main()
