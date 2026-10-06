"""
genome_resolver.py -- discover the source genome from a TF sequence ALONE (Phase 0.2 front step).

The pipeline's only input is a TF sequence (amino-acid or DNA); no genome accession is supplied. We
BLAST the sequence against a genomic database, keep the top-scoring genomes, and first-pass locate the
TF CDS within each -- returning ranked `GenomeCandidate`s carrying the TF locus. Mirrors the
operators-pipeline approach (blast -> top genomes -> find the TF) but as the genuinely-missing front
step (`fetch_ncbi_cds.js` started one step too late, assuming the accession was already known).

Two backends:
  local : tblastn (protein) / blastn (DNA) against a local nucleotide DB (makeblastdb). Offline and
          deterministic -- the engine; tested here.
  ncbi  : remote blastp -> Identical Protein Groups (IPG) -> genomic loci, + EFetch to pull the region.
          The production "sequence only" path (network + NCBI API key from arsr_merr_families/.env).

The resolved (accession, TF locus, genome record) feeds `annotate/context.py` (inter-operon
extraction) and `signals/motif_rescan.py`.

Run `python genome_resolver.py` for a self-test (offline local backend; live NCBI fetch if reachable).
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

def _blast_dir() -> Path | None:
    """Directory holding the BLAST+ binaries: `$BLAST_BIN` if set, else wherever `blastn` sits on PATH.
    None when BLAST+ is not installed -- callers degrade rather than shelling out to a missing binary.
    (Previously defaulted to one developer's ncbi-blast install, which exists on no other machine.)"""
    env = os.environ.get("BLAST_BIN")
    if env:
        return Path(env)
    import shutil
    exe = shutil.which("blastn") or shutil.which("blastn.exe")
    return Path(exe).parent if exe else None


def _env_file():
    """The .env carrying the NCBI/biohub tokens, resolved via the central resource registry (repo-local
    `.env` / `PREDICTOR_ENV_FILE` override, falling back to the legacy 5.8 location). Imported lazily so
    this module keeps working standalone if resources.py is unavailable."""
    try:
        from .. import resources
    except Exception:
        try:
            from predictor import resources
        except Exception:
            return Path(__file__).resolve().parents[2].parent / "5.8 Promoters" / "arsr_merr_families" / ".env"
    return resources.env_file() or (Path(__file__).resolve().parents[2].parent /
                                    "5.8 Promoters" / "arsr_merr_families" / ".env")

_DNA = set("ACGTNacgtn")
_OUTFMT = "6 sseqid pident length qlen sstart send evalue bitscore sstrand"


@dataclass
class GenomeCandidate:
    accession: str
    tf_start: int            # 0-based half-open, + strand coords
    tf_end: int
    tf_strand: str           # '+' or '-'
    identity: float          # fraction in [0,1]
    bitscore: float
    coverage: float          # fraction of the query aligned
    organism: str = ""
    source: str = ""


# --------------------------------------------------------------------------- helpers
def is_dna(seq: str) -> bool:
    s = seq.strip().replace("\n", "")
    if not s:
        return False
    return sum(c in _DNA for c in s) / len(s) > 0.9


def _blast(name: str) -> str:
    d = _blast_dir()
    if d is not None:
        exe = d / (name + (".exe" if os.name == "nt" else ""))
        if exe.exists():
            return str(exe)
    found = shutil.which(name)
    if found:
        return found
    raise FileNotFoundError(
        f"BLAST+ tool '{name}' not found. Install BLAST+ and put it on PATH, or set $BLAST_BIN to the "
        f"directory holding the binaries (searched: {d or 'PATH only'}).")


def _run(cmd, cwd=None):
    res = subprocess.run([str(c) for c in cmd], cwd=cwd and str(cwd),
                         capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"{Path(str(cmd[0])).name} failed (rc={res.returncode}):\n"
                           f"{(res.stderr or res.stdout)[-1500:]}")
    return res


def _write_fasta(path: Path, records) -> None:
    with path.open("w", encoding="ascii") as fh:
        for rid, seq in records:
            fh.write(f">{rid}\n{seq}\n")


def _clean_acc(sid: str) -> str:
    """Normalize a BLAST sseqid (e.g. 'ref|NC_000964.3|' from a -parse_seqids DB) to the accession."""
    if "|" not in sid:
        return sid
    parts = [p for p in sid.split("|") if p]
    for p in parts:
        if "." in p and any(ch.isdigit() for ch in p):
            return p
    return parts[-1] if parts else sid


def _query_seq(query) -> str:
    if isinstance(query, (str, os.PathLike)) and Path(str(query)).exists():
        txt = Path(query).read_text(encoding="ascii")
        return "".join(l.strip() for l in txt.splitlines() if not l.startswith(">"))
    return str(query).strip()


# --------------------------------------------------------------------------- local backend
def _search_db(seq, db_prefix, *, top_n, min_identity, evalue, threads, workdir):
    """Run tblastn/blastn of `seq` against a prebuilt nucleotide DB prefix; parse to candidates."""
    qf = workdir / "q.fasta"
    _write_fasta(qf, [("query", seq)])
    prog = "blastn" if is_dna(seq) else "tblastn"
    out = workdir / "hits.tsv"
    cmd = [_blast(prog), "-query", qf, "-db", db_prefix, "-evalue", evalue,
           "-num_threads", threads, "-outfmt", _OUTFMT, "-out", out]
    if prog == "blastn":
        cmd += ["-task", "blastn"]
    _run(cmd, cwd=workdir)

    best = {}                                                    # sseqid -> best HSP fields
    for line in out.read_text().splitlines():
        if not line.strip():
            continue
        sid, pid, length, qlen, ss, se, ev, bits, sstrand = line.split("\t")
        bits = float(bits)
        if sid not in best or bits > best[sid][0]:
            best[sid] = (bits, float(pid) / 100.0, int(length), int(qlen),
                         int(ss), int(se), sstrand)
    cands = []
    for sid, (bits, ident, length, qlen, ss, se, sstrand) in best.items():
        start0, end0 = min(ss, se) - 1, max(ss, se)
        strand = "+" if sstrand.lower().startswith("plus") else "-"
        cov = length / qlen if qlen else 0.0
        if ident >= min_identity:
            cands.append(GenomeCandidate(_clean_acc(sid), start0, end0, strand, ident, bits, cov,
                                         source=f"local-{prog}"))
    cands.sort(key=lambda c: c.bitscore, reverse=True)
    return cands[:top_n]


def resolve_local(query, genome_fasta=None, *, db=None, top_n: int = 10,
                  min_identity: float = 0.30, evalue: float = 1e-5,
                  threads: int | None = None, workdir: str | None = None):
    """tblastn (protein) / blastn (DNA) the query against a local nucleotide DB; return ranked
    GenomeCandidates (best HSP per subject) with the located TF locus.

    Pass `db=` (a prebuilt makeblastdb prefix, e.g. from genome_db) to skip DB construction, or
    `genome_fasta=` to build a one-shot DB in a temp dir.
    """
    if db is None and genome_fasta is None:
        raise ValueError("provide db= (prebuilt prefix) or genome_fasta=")
    seq = _query_seq(query)
    threads = threads or max(1, (os.cpu_count() or 2) // 2)
    owns = workdir is None
    tmp = Path(workdir or tempfile.mkdtemp(prefix="genres_"))
    tmp.mkdir(parents=True, exist_ok=True)
    try:
        if db is None:
            db = str(tmp / "db")
            _run([_blast("makeblastdb"), "-in", genome_fasta, "-dbtype", "nucl", "-out", db], cwd=tmp)
        return _search_db(seq, str(db), top_n=top_n, min_identity=min_identity,
                          evalue=evalue, threads=threads, workdir=tmp)
    finally:
        if owns:
            shutil.rmtree(tmp, ignore_errors=True)


def resolve(query, *, db=None, genome_fasta=None, allow_ncbi: bool = False, **kw):
    """Dispatcher: prefer an explicit/prebuilt local DB (fast, reproducible), else a one-shot DB
    from `genome_fasta`, else the default local genome DB, else (if allowed) remote NCBI."""
    if db is None and genome_fasta is None:
        try:
            from .genome_db import default_db
        except ImportError:
            from .genome_db import default_db
        db = default_db()
    if db is not None:
        return resolve_local(query, db=db, **kw)
    if genome_fasta is not None:
        return resolve_local(query, genome_fasta=genome_fasta, **kw)
    if allow_ncbi:
        return resolve_ncbi(query, **kw)
    raise RuntimeError("no local genome DB found; build one with genome_db.build_local_db() "
                       "or pass allow_ncbi=True")


# --------------------------------------------------------------------------- NCBI backend
def _ncbi_creds():
    # accept either NCBI_API_KEY or NCBI_TOKEN (the name the user dropped into the .env), from env or file
    email = os.environ.get("NCBI_EMAIL", "")
    key = os.environ.get("NCBI_API_KEY", "") or os.environ.get("NCBI_TOKEN", "")
    _ef = _env_file()
    if (not email or not key) and _ef and _ef.exists():
        txt = _ef.read_text(errors="ignore")
        email = email or (re.search(r"NCBI_EMAIL\s*=\s*(\S+)", txt) or [None, ""])[1]
        key = key or (re.search(r"NCBI_(?:API_KEY|TOKEN)\s*=\s*(\S+)", txt) or [None, ""])[1]
    return email or "anonymous@example.com", key or None


def _entrez():
    from Bio import Entrez
    email, key = _ncbi_creds()
    Entrez.email = email
    if key:
        Entrez.api_key = key
    return Entrez


def fetch_region(accession: str, start: int | None = None, end: int | None = None,
                 *, rettype: str = "gb") -> str:
    """EFetch a nuccore record (or sub-region) as GenBank/FASTA text. `start`/`end` are 1-based
    inclusive NCBI coordinates when given. Consumed by `annotate/context.py`."""
    Entrez = _entrez()
    kw = dict(db="nuccore", id=accession, rettype=rettype, retmode="text")
    if start is not None and end is not None:
        kw.update(seq_start=str(start), seq_stop=str(end))
    with Entrez.efetch(**kw) as h:
        return h.read()


# --------------------------------------------------------------------------- promoter window convention
#: THE promoter-window convention, shared by every homolog-region extractor (local blastdbcmd, NCBI
#: EFetch, and the SSN-cluster/IPG path). Keep these three in lockstep -- a region set built under two
#: conventions cannot be compared position-by-position.
PROMOTER_UPSTREAM = 350
PROMOTER_DOWNSTREAM = 30


def promoter_window_bounds(tf_start: int, tf_end: int, tf_strand: str, *,
                           upstream: int = PROMOTER_UPSTREAM,
                           downstream: int = PROMOTER_DOWNSTREAM) -> tuple[int, int]:
    """1-based inclusive (start, stop) of a locus's promoter-side window on the genome's + strand."""
    if tf_strand == "+":
        s1, e1 = tf_start + 1 - upstream, tf_start + downstream
    else:
        s1, e1 = tf_end - downstream, tf_end + upstream
    return max(1, s1), e1


def orient_promoter(seq: str, tf_strand: str) -> str:
    """Put a + strand promoter window into TRANSCRIPTIONAL orientation: 5' -> 3' toward the gene.

    Every extractor slices the genome's + strand, so a MINUS-strand gene's promoter comes back mirrored
    AND complemented relative to a plus-strand one. The old convention left it that way and described the
    regions as "directly comparable for the conservation step" -- they were only comparable to a score
    that ignores position and strand, which is exactly what the old `motif_finder.conservation` did (any
    significant hit, either strand, anywhere in the window). Reverse-complementing the minus-strand window
    puts every region in one frame, with the gene start at a fixed offset (~`upstream` from the 5' end),
    so aligned positional conservation becomes meaningful. Scans stay both-strands, so a palindromic
    operator is unaffected either way."""
    if not seq:
        return ""
    return _revcomp(seq) if tf_strand == "-" else seq


def fetch_upstream_region(cand: GenomeCandidate, *, upstream: int = PROMOTER_UPSTREAM,
                          downstream: int = PROMOTER_DOWNSTREAM) -> str:
    """EFetch the promoter-side (5') window of a homolog locus -- the online twin of
    `homolog_regions._extract` (blastdbcmd) for accessions that are NOT in the local mirror.

    Returned in transcriptional orientation (`orient_promoter`), identical in convention to the local and
    SSN-cluster paths, so regions from all three sources are comparable position-by-position. Returns ''
    on failure or a too-short window."""
    s1, e1 = promoter_window_bounds(cand.tf_start, cand.tf_end, cand.tf_strand,
                                    upstream=upstream, downstream=downstream)
    if e1 <= s1:
        return ""
    try:
        fa = fetch_region(cand.accession, s1, e1, rettype="fasta")
    except Exception:
        return ""
    seq = "".join(l.strip() for l in fa.splitlines() if not l.startswith(">")).upper()
    if len(seq) < 30:
        return ""
    return orient_promoter(seq, cand.tf_strand)


def resolve_ncbi(query, *, top_n: int = 10, db: str = "nr", hitlist: int = 50,
                 max_loci: int | None = None):
    """Production 'sequence only' path: remote blastp -> top protein accessions -> IPG -> genomic
    loci. Network + NCBI key. (Remote BLAST is slow; not exercised by the offline self-test.)

    `top_n` caps how many protein hits are expanded via IPG; `max_loci` caps the returned genomic loci
    (defaults to top_n). For the Snowprint conservation fallback, pass a larger `top_n`/`max_loci` so
    enough homolog promoter regions are collected to make the conservation signal usable."""
    from io import StringIO
    from Bio.Blast import NCBIWWW, NCBIXML
    Entrez = _entrez()
    seq = _query_seq(query)
    program = "blastx" if is_dna(seq) else "blastp"

    with NCBIWWW.qblast(program, db, seq, hitlist_size=hitlist) as h:
        rec = NCBIXML.read(StringIO(h.read()))
    accs = []
    for aln in rec.alignments[:top_n]:
        m = re.search(r"\b([A-Z]{1,3}_?\d{5,})\.?\d*\b", aln.accession or aln.hit_id or "")
        if m:
            accs_id = aln.accession or m.group(1)
            accs.append(accs_id)

    # Per-accession failures were swallowed silently. That is what let a TypeError on EVERY
    # accession look like a genome with no homologs: 25 exceptions, an empty list, and no output.
    # Failures are still tolerated -- one bad accession must not lose the other 24 -- but they are
    # now counted, and an ALL-failed batch says so, because "every lookup raised" and "this protein
    # has no identical-protein group" are not the same answer and used to be indistinguishable.
    cands = []
    errors: dict[str, int] = {}
    for acc in accs:
        try:
            with Entrez.efetch(db="protein", id=acc, rettype="ipg", retmode="text") as h:
                cands += _parse_ipg(h.read())
        except Exception as e:
            k = type(e).__name__
            errors[k] = errors.get(k, 0) + 1
            continue
    if accs and len(accs) == sum(errors.values()):
        import warnings
        warnings.warn(
            f"IPG lookup failed for ALL {len(accs)} blastp hits ({errors}); "
            f"resolve_ncbi is returning no loci. This is a defect, not an empty result.",
            RuntimeWarning, stacklevel=2)
    # dedupe by (accession, locus): one protein maps to many genomes via IPG, each a distinct homolog
    # promoter for the conservation step -- so key on the locus, not just the contig accession.
    seen = {}
    for c in cands:
        c.source = "ncbi-ipg"
        key = (c.accession, c.tf_start, c.tf_strand)
        if key not in seen:
            seen[key] = c
    return list(seen.values())[:(max_loci or top_n)]


def _parse_ipg(text):
    """Parse an Identical Protein Groups TSV into GenomeCandidates (genomic loci of the protein).

    Accepts BYTES OR STR. `Entrez.efetch(rettype="ipg")` returns bytes even under
    `retmode="text"`, and splitting bytes on a str separator raises
    `TypeError: a bytes-like object is required, not 'str'`.

    This is not hypothetical, and it is the reason the remote homolog fallback never worked.
    `resolve_ncbi` passed the handle straight through (`_parse_ipg(h.read())`) while
    `msa_homologs.loci_from_protein` decoded first; only the second worked. Since the caller wraps
    each accession in `except Exception: continue`, ALL of them raised, `cands` came back empty, and
    the fallback returned zero homolog regions on all 84 invocations of run A -- 30 min to 6 h of
    blastp discarded silently each time, and the direct cause of that run's eight timeouts. Nothing
    surfaced, because an empty result is indistinguishable from "the local mirror was enough".

    Normalising here rather than at the one call site is deliberate: the type an efetch handle
    yields is not something every caller should have to remember, and the version that remembered
    and the version that forgot sat two modules apart for months.
    """
    if isinstance(text, (bytes, bytearray)):
        text = text.decode("utf-8", "replace")
    out, header = [], None
    for line in text.splitlines():
        cols = line.split("\t")
        if header is None:
            header = {name: i for i, name in enumerate(cols)}
            continue
        try:
            acc = cols[header.get("Nucleotide Accession", 2)]
            start = int(cols[header.get("Start", 3)])
            stop = int(cols[header.get("Stop", 4)])
            strand = cols[header.get("Strand", 5)]
            org = cols[header.get("Organism", 8)] if len(cols) > 8 else ""
        except (IndexError, ValueError):
            continue
        if not acc:
            continue
        out.append(GenomeCandidate(acc, min(start, stop) - 1, max(start, stop),
                                   "+" if strand == "+" else "-", 1.0, 0.0, 1.0, org))
    return out


# --------------------------------------------------------------------------- self-test
_CODON = {"A": "GCT", "R": "CGT", "N": "AAT", "D": "GAT", "C": "TGT", "Q": "CAA", "E": "GAA",
          "G": "GGT", "H": "CAT", "I": "ATT", "L": "CTT", "K": "AAA", "M": "ATG", "F": "TTT",
          "P": "CCT", "S": "TCT", "T": "ACT", "W": "TGG", "Y": "TAT", "V": "GTT", "*": "TAA"}


def _revtrans(prot: str) -> str:
    return "".join(_CODON.get(a, "NNN") for a in prot)


def _revcomp(dna: str) -> str:
    return dna.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _demo() -> None:
    import random
    random.seed(11)

    prot = ("MSEKQDLTVKDLAKETGLSVHTLRYYERIGLLPEPDRSEGNYRLYTQAHLERLAFIKRAKRLGFSL"
            "EEIAELLALWDDRHRASADVKAIAQAHLAE")
    cds = _revtrans(prot)

    def bg(n):
        return "".join(random.choice("ACGT") for _ in range(n))

    # contig_A: TF CDS on + strand at offset 500;  contig_B: same CDS on - strand at 800;  contig_C: decoy
    cA = bg(500) + cds + bg(400)
    cB = bg(800) + _revcomp(cds) + bg(300)
    cC = bg(1500)

    tmp = Path(tempfile.mkdtemp(prefix="genresdemo_"))
    try:
        gfa = tmp / "genomes.fasta"
        _write_fasta(gfa, [("contig_A", cA), ("contig_B", cB), ("contig_C_decoy", cC)])
        cands = resolve_local(prot, gfa, workdir=str(tmp / "work"))
        print(f"resolved {len(cands)} candidate genome(s):")
        for c in cands:
            print(f"  {c.accession:<16} id={c.identity:.0%} cov={c.coverage:.2f} "
                  f"bits={c.bitscore:.0f} locus={c.tf_start}-{c.tf_end}({c.tf_strand}) [{c.source}]")

        top = cands[0]
        assert top.accession == "contig_A", "top genome should be the + strand contig"
        assert top.tf_strand == "+" and abs(top.tf_start - 500) <= 3, "TF locus/strand wrong on +"
        b = next(c for c in cands if c.accession == "contig_B")
        assert b.tf_strand == "-" and abs(b.tf_start - 800) <= 3, "TF locus/strand wrong on -"
        assert all(c.accession != "contig_C_decoy" for c in cands), "decoy should not resolve"
        print("OK: TF located on the correct genomes/strands; decoy excluded.")

        # optional live NCBI fetch (fast) -- proves the production fetch path
        try:
            gb = fetch_region("NC_000913.3", 1, 800, rettype="gb")
            assert "LOCUS" in gb and "ORIGIN" in gb
            print(f"OK (network): EFetch returned {len(gb)} bytes of GenBank for NC_000913.3:1-800.")
        except Exception as e:
            print(f"SKIP (network): NCBI fetch not exercised ({type(e).__name__}: {e}).")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    _demo()
