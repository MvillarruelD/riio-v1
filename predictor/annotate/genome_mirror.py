"""
genome_mirror.py -- build / refresh a DEREPLICATED representative genome mirror (FASTA + GFF) for the
local homolog+context layer, via the **NCBI Datasets** command-line tool.

Why a dereplicated mirror rather than "all of NCBI":
  * "All" prokaryotic assemblies are millions of genomes / multiple TB, mostly redundant, and stale the
    day they are downloaded. Snowprint itself never mirrors -- it queries NCBI per request.
  * For *phylogenetic footprinting* (the conservation signal that separates a real operator from a
    chance palindrome), dereplication IS the quality control: thousands of near-identical E. coli
    genomes inflate apparent conservation with redundancy, not evolutionary signal. So we mirror a
    *representative* set -- GTDB species reps (~110k, R220) or RefSeq "reference/representative"
    prokaryotes -- which is tens of GB as genome+annotation, not TB.
  * The operator lives in the **intergenic region of the homolog's operon**, so the gene COORDINATES
    (GFF) are the real requirement, not just the protein. A protein-only DB (or ESM Atlas) cannot give
    the DNA context. Every mirrored genome is therefore stored as `<acc>.fna` + `<acc>.gff`.

The mirror feeds `genome_db.build_local_db()` (BLAST DB + `blastdbcmd` region extraction). NCBI
IPG/EFetch stays as the *online fallback* for the sparse-homolog regime (see `homolog_regions.py`).

Layout: `$PREDICTOR_GENOME_MIRROR` (default `~/.predictor/genome_mirror`, space-free -- BLAST+ 2.17
mis-splits DB paths at a space). The NCBI Datasets CLI is located via `$DATASETS_EXE` -> PATH; install
it from https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/ (the
single-binary `datasets`), or `conda install -c conda-forge ncbi-datasets-cli`.

  add-taxon TAXON   download the RefSeq reference set for a taxon -> mirror/  (e.g. add-taxon bacteria)
  build-set FILE    download genomes for a newline-delimited accession list (e.g. a GTDB rep list)
  info              summarize the mirror (n genomes, how many carry GFF, size on disk)
  rebuild-db        (re)build the local BLAST DB including the mirror (delegates to genome_db)
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

from predictor.annotate.genome_db import MIRROR_DIR, mirror_genome_fastas, build_local_db

# datasets packs each genome as ncbi_dataset/data/<accession>/{<acc>_<asm>_genomic.fna, genomic.gff}
_INCLUDE = "genome,gff3"
_INSTALL_HINT = (
    "NCBI Datasets CLI not found. Install the single-binary `datasets` from "
    "https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/ "
    "(or `conda install -c conda-forge ncbi-datasets-cli`), or set $DATASETS_EXE to its path."
)


# --------------------------------------------------------------------------- CLI plumbing
def find_datasets() -> str:
    """Locate the NCBI Datasets entry point: $DATASETS_EXE -> `datasets` on PATH."""
    env = os.environ.get("DATASETS_EXE")
    if env and Path(env).exists():
        return env
    found = shutil.which("datasets")
    if found:
        return found
    raise FileNotFoundError(_INSTALL_HINT)


def _run_datasets(datasets: str, args, *, cwd: Path) -> subprocess.CompletedProcess:
    res = subprocess.run([datasets, *map(str, args)], cwd=str(cwd),
                         capture_output=True, text=True)
    if res.returncode != 0:
        tail = (res.stderr or res.stdout or "")[-2000:]
        raise RuntimeError(f"datasets {args[0]} {args[1] if len(args) > 1 else ''} "
                           f"failed (rc={res.returncode}):\n{tail}")
    return res


# --------------------------------------------------------------------------- zip ingestion
def _pick_fna(acc_dir: Path) -> Path | None:
    """The genomic FASTA inside an accession dir (largest *.fna/_genomic.fna)."""
    fnas = sorted(acc_dir.glob("*.fna")) + sorted(acc_dir.glob("*_genomic.fna"))
    fnas = sorted(set(fnas), key=lambda p: p.stat().st_size, reverse=True)
    return fnas[0] if fnas else None


def ingest_zip(zip_path, *, mirror_dir: Path = MIRROR_DIR) -> list:
    """Unpack a `datasets download genome` zip into the mirror as `<acc>.fna` (+ `<acc>.gff`).

    Returns the list of accessions ingested. Genomes without a FASTA are skipped; a missing GFF is
    tolerated (logged via the return of `mirror_info`) but means that genome can only contribute
    homolog hits, not gene-aware context extraction.
    """
    mirror_dir = Path(mirror_dir)
    mirror_dir.mkdir(parents=True, exist_ok=True)
    ingested = []
    with tempfile.TemporaryDirectory(prefix="ncbi_zip_") as td:
        tdp = Path(td)
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(tdp)
        data_root = tdp / "ncbi_dataset" / "data"
        if not data_root.exists():
            raise FileNotFoundError(f"unexpected datasets zip layout (no ncbi_dataset/data in {zip_path})")
        for acc_dir in sorted(p for p in data_root.iterdir() if p.is_dir()):
            acc = acc_dir.name
            fna = _pick_fna(acc_dir)
            if fna is None:
                continue
            shutil.copyfile(fna, mirror_dir / f"{acc}.fna")
            gff = acc_dir / "genomic.gff"
            if gff.exists():
                shutil.copyfile(gff, mirror_dir / f"{acc}.gff")
            ingested.append(acc)
    return ingested


# --------------------------------------------------------------------------- downloads
def download_taxon(taxon: str, *, reference: bool = True, assembly_level: str | None = None,
                   limit: int | None = None, dry_run: bool = False, mirror_dir: Path = MIRROR_DIR):
    """Download a taxon's RefSeq set (default: `--reference`, one representative per species) with GFF.

    `assembly_level` (e.g. "complete") and `limit` narrow the pull. With `dry_run`, only the command is
    returned (no network) -- useful to inspect / log the exact `datasets` invocation.
    """
    ds = find_datasets()
    args = ["download", "genome", "taxon", taxon, "--include", _INCLUDE,
            "--assembly-source", "RefSeq"]
    if reference:
        args.append("--reference")
    if assembly_level:
        args += ["--assembly-level", assembly_level]
    if limit:
        args += ["--limit", limit]
    if dry_run:
        return {"cmd": [ds, *map(str, args)], "dry_run": True}
    with tempfile.TemporaryDirectory(prefix="ncbi_dl_") as td:
        zp = Path(td) / "genomes.zip"
        _run_datasets(ds, [*args, "--filename", zp], cwd=Path(td))
        accs = ingest_zip(zp, mirror_dir=mirror_dir)
    return {"taxon": taxon, "ingested": len(accs), "accessions": accs[:20], "mirror": str(mirror_dir)}


def build_set(accession_file, *, batch: int = 2000, dry_run: bool = False,
              mirror_dir: Path = MIRROR_DIR):
    """Download genomes for a newline-delimited accession list (e.g. a GTDB species-rep list) with GFF.

    Large lists are pulled in batches (datasets accepts an --inputfile of accessions). `dry_run` returns
    the planned batches/commands without touching the network.
    """
    ds = find_datasets()
    accs = [l.strip() for l in Path(accession_file).read_text(encoding="ascii").splitlines()
            if l.strip() and not l.startswith("#")]
    if not accs:
        raise ValueError(f"no accessions in {accession_file}")
    batches = [accs[i:i + batch] for i in range(0, len(accs), batch)]
    if dry_run:
        return {"n_accessions": len(accs), "n_batches": len(batches), "datasets": ds}
    ingested = []
    for bi, chunk in enumerate(batches):
        with tempfile.TemporaryDirectory(prefix="ncbi_dl_") as td:
            tdp = Path(td)
            lst = tdp / "accs.txt"
            lst.write_text("\n".join(chunk), encoding="ascii")
            zp = tdp / "genomes.zip"
            _run_datasets(ds, ["download", "genome", "accession", "--inputfile", lst,
                               "--include", _INCLUDE, "--filename", zp], cwd=tdp)
            ingested += ingest_zip(zp, mirror_dir=mirror_dir)
    return {"n_accessions": len(accs), "ingested": len(ingested), "mirror": str(mirror_dir)}


# --------------------------------------------------------------------------- per-accession, ON DEMAND
def _assembly_for_nuccore(nuccore_acc: str, *, retries: int = 4) -> str | None:
    """Map a nucleotide accession (NC_/NZ_/CP...) to its genome ASSEMBLY accession (GCF_/GCA_) via NCBI
    elink. NCBI Datasets downloads by assembly, but our resolver returns nucleotide accessions, so this
    bridges the two. EUtils is flaky (transient EOF/500), so each step is retried with backoff.
    Best-effort: returns None only after exhausting retries."""
    import time
    from Bio import Entrez
    Entrez.email = os.environ.get("NCBI_EMAIL", "predictor@example.com")
    if os.environ.get("NCBI_API_KEY"):
        Entrez.api_key = os.environ["NCBI_API_KEY"]

    def _retry(fn):
        last = None
        for i in range(retries):
            try:
                return fn()
            except Exception as e:                        # transient eutils errors -> back off and retry
                last = e
                time.sleep(1.5 * (i + 1))
        raise last

    try:
        h = _retry(lambda: Entrez.read(Entrez.esearch(db="nuccore", term=nuccore_acc)))
        if not h["IdList"]:
            return None
        link = _retry(lambda: Entrez.read(Entrez.elink(dbfrom="nuccore", db="assembly",
                                                       id=h["IdList"][0])))
        ls = link[0].get("LinkSetDb") or []
        if not ls:
            return None
        asm_uid = ls[0]["Link"][0]["Id"]
        summ = _retry(lambda: Entrez.read(Entrez.esummary(db="assembly", id=asm_uid)))
        doc = summ["DocumentSummarySet"]["DocumentSummary"][0]
        return doc.get("AssemblyAccession")
    except Exception:
        return None


def fetch_accession(accession: str, *, dry_run: bool = False, mirror_dir: Path = MIRROR_DIR) -> dict:
    """Download ONE genome on demand (FASTA + GFF) into the mirror -- the per-genome acquisition the
    pipeline triggers when it has resolved a TF's source genome and needs gene-aware / genome-wide context.
    This is deliberately per-accession (a few MB for a bacterium), NOT a bulk mirror and NOT a local nt
    BLAST DB. Accepts an assembly accession (GCF_/GCA_) directly, or a nucleotide accession (mapped via
    `_assembly_for_nuccore`). Needs the NCBI Datasets CLI (`find_datasets`)."""
    ds = find_datasets()                                  # raises a clear install hint if absent
    asm = accession
    if not accession.upper().startswith(("GCF_", "GCA_")):
        asm = _assembly_for_nuccore(accession)
        if asm is None:
            raise RuntimeError(f"could not map {accession} to a genome assembly (network/elink); "
                               f"pass a GCF_/GCA_ accession explicitly")
    args = ["download", "genome", "accession", asm, "--include", _INCLUDE]
    if dry_run:
        return {"cmd": [ds, *map(str, args)], "assembly": asm, "from": accession, "dry_run": True}
    with tempfile.TemporaryDirectory(prefix="ncbi_one_") as td:
        zp = Path(td) / "genome.zip"
        _run_datasets(ds, [*args, "--filename", zp], cwd=Path(td))
        accs = ingest_zip(zp, mirror_dir=mirror_dir)
    return {"requested": accession, "assembly": asm, "ingested": accs, "mirror": str(mirror_dir)}


def ensure_genome(accession: str, *, allow_download: bool = True, mirror_dir: Path = MIRROR_DIR):
    """Return (fna_path, gff_path|None) for an accession, mirror-first and acquiring on demand.

    Resolution order: (1) already mirrored -> return; (2) `allow_download` + NCBI Datasets CLI ->
    `fetch_accession` (per-genome). When the CLI is absent the caller should fall back to the EFetch
    *window* path (`genome_resolver.fetch_region`), which needs no whole-genome download and is what the
    structure-free operator discovery already uses. Returns (None, None) if not acquired."""
    mirror_dir = Path(mirror_dir)
    # the mirror keys on whatever accession was ingested; match by stem prefix for GCF/GCA or exact for nuccore
    for fna in list(mirror_dir.glob(f"{accession}*.fna")) + ([mirror_dir / f"{accession}.fna"]
                                                             if (mirror_dir / f"{accession}.fna").exists() else []):
        gff = fna.with_suffix(".gff")
        return fna, (gff if gff.exists() else None)
    if not allow_download:
        return None, None
    try:
        info = fetch_accession(accession, mirror_dir=mirror_dir)
    except (FileNotFoundError, RuntimeError):
        return None, None
    acc = (info.get("ingested") or [None])[0]
    if not acc:
        return None, None
    fna = mirror_dir / f"{acc}.fna"
    gff = fna.with_suffix(".gff")
    return (fna if fna.exists() else None), (gff if gff.exists() else None)


# --------------------------------------------------------------------------- info / rebuild
def mirror_info(mirror_dir: Path = MIRROR_DIR) -> dict:
    mirror_dir = Path(mirror_dir)
    fnas = mirror_genome_fastas() if mirror_dir == MIRROR_DIR else sorted(
        str(p) for p in mirror_dir.glob("*.fna"))
    n_gff = sum(1 for f in fnas if Path(f).with_suffix(".gff").exists())
    size = sum(Path(f).stat().st_size for f in fnas) + sum(
        Path(f).with_suffix(".gff").stat().st_size for f in fnas
        if Path(f).with_suffix(".gff").exists())
    return {"mirror": str(mirror_dir), "n_genomes": len(fnas), "n_with_gff": n_gff,
            "size_mb": round(size / 1e6, 1)}


def rebuild_db():
    """(Re)build the local BLAST DB including the mirror (curated + mirror + extra)."""
    return build_local_db()


# --------------------------------------------------------------------------- CLI / self-test
def _main(argv):
    cmd = argv[0] if argv else "info"
    if cmd == "add-taxon":
        if len(argv) < 2:
            print("usage: add-taxon TAXON [--all] [--level complete] [--limit N] [--dry-run]")
            return
        kw = dict(reference="--all" not in argv, dry_run="--dry-run" in argv)
        if "--level" in argv:
            kw["assembly_level"] = argv[argv.index("--level") + 1]
        if "--limit" in argv:
            kw["limit"] = int(argv[argv.index("--limit") + 1])
        print(json.dumps(download_taxon(argv[1], **kw), indent=2))
    elif cmd == "build-set":
        print(json.dumps(build_set(argv[1], dry_run="--dry-run" in argv), indent=2))
    elif cmd == "fetch":
        if len(argv) < 2:
            print("usage: fetch <GCF_/GCA_ or NC_/NZ_ accession> [--dry-run]")
            return
        print(json.dumps(fetch_accession(argv[1], dry_run="--dry-run" in argv), indent=2))
    elif cmd == "info":
        print(json.dumps(mirror_info(), indent=2))
    elif cmd == "rebuild-db":
        prefix, n = rebuild_db()
        print(f"rebuilt {prefix} from {n} sequences (curated + mirror + extra)")
    else:
        print(__doc__)


def _demo() -> None:
    """Offline self-test: synthesize a `datasets`-shaped zip, ingest it, verify FASTA+GFF pairing.

    Exercises the ingestion + mirror bookkeeping (the network-bound download is SKIPped unless the
    `datasets` CLI is installed, in which case a real --dry-run command build is checked)."""
    tmp = Path(tempfile.mkdtemp(prefix="mirrordemo_"))
    try:
        # build a fake datasets zip: 2 accessions, one with a GFF, one without
        data = tmp / "pkg" / "ncbi_dataset" / "data"
        a1 = data / "GCF_000005845.2"
        a2 = data / "GCF_000009999.1"
        a1.mkdir(parents=True)
        a2.mkdir(parents=True)
        (a1 / "GCF_000005845.2_ASM584v2_genomic.fna").write_text(">seqA\nACGTACGTACGT\n", encoding="ascii")
        (a1 / "genomic.gff").write_text("##gff-version 3\nseqA\t.\tgene\t1\t9\t.\t+\t.\tID=g1\n",
                                        encoding="ascii")
        (a2 / "GCF_000009999.1_ASM999v1_genomic.fna").write_text(">seqB\nTTTTGGGGCCCC\n", encoding="ascii")
        # a2 deliberately has no genomic.gff
        zp = tmp / "genomes.zip"
        with zipfile.ZipFile(zp, "w") as zf:
            for f in (tmp / "pkg").rglob("*"):
                if f.is_file():
                    zf.write(f, f.relative_to(tmp / "pkg"))

        mirror = tmp / "mirror"
        accs = ingest_zip(zp, mirror_dir=mirror)
        print(f"ingested {len(accs)} accessions: {accs}")
        assert set(accs) == {"GCF_000005845.2", "GCF_000009999.1"}, "accession ingestion wrong"
        assert (mirror / "GCF_000005845.2.fna").exists(), "fna not placed"
        assert (mirror / "GCF_000005845.2.gff").exists(), "gff not paired"
        assert not (mirror / "GCF_000009999.1.gff").exists(), "non-existent gff should not appear"

        info = mirror_info(mirror)
        print(f"mirror_info: {info}")
        assert info["n_genomes"] == 2 and info["n_with_gff"] == 1, "mirror_info miscount"

        # ensure_genome is mirror-first: the just-ingested accession resolves with no download
        fna, gff = ensure_genome("GCF_000005845.2", mirror_dir=mirror, allow_download=False)
        print(f"ensure_genome (mirror hit): fna={None if not fna else fna.name} gff={None if not gff else gff.name}")
        assert fna and fna.name == "GCF_000005845.2.fna" and gff, "ensure_genome should hit the mirror"

        # command construction is correct even without the CLI installed (dry-run path)
        try:
            cmd = download_taxon("bacteria", dry_run=True)
            assert "--reference" in cmd["cmd"] and "gff3" in " ".join(cmd["cmd"]), "taxon cmd malformed"
            print(f"datasets dry-run cmd: {' '.join(cmd['cmd'])}")
            one = fetch_accession("GCF_000005845.2", dry_run=True)      # per-accession on-demand command
            assert one["assembly"] == "GCF_000005845.2" and "accession" in one["cmd"], "fetch cmd malformed"
            print(f"per-accession fetch dry-run: {' '.join(one['cmd'])}")
        except FileNotFoundError:
            print("SKIP: NCBI Datasets CLI not installed (download path not exercised; ingestion OK).")
        print("OK: datasets zip ingested into FASTA+GFF mirror; ensure_genome + per-accession fetch wired.")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        _main(sys.argv[1:])
    else:
        _demo()
