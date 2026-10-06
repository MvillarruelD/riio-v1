"""
genome_db.py -- build / maintain a LOCAL genome BLAST DB for fast, offline, reproducible genome
resolution (vs slow & flaky remote blastp per TF).

The DB is built from three sources, concatenated and dereplicated by accession:
  * the curated pipeline genomes (`arsr_merr_families/genomes/*/*_genomic.fna`),
  * the **dereplicated representative mirror** (`genome_mirror.py` -> `~/.predictor/genome_mirror/*.fna`,
    GTDB/RefSeq species reps with GFF coordinates -- the bulk local-homolog+context layer), and
  * ad-hoc genomes appended by accession via EFetch (`extra/`).
`genome_resolver.resolve()` uses this DB by default and only falls back to remote NCBI IPG/EFetch when a
sequence has too few confident hits locally (the sparse-stratum fallback).

  build       python genome_db.py build      # (re)build ~/.predictor/genome_db/arsr_merr_genomes
  info        python genome_db.py info
  add ACC...  python genome_db.py add NC_000913.3   # fetch + append a genome, then rebuild
"""
from __future__ import annotations

import glob
import os
import sys
from pathlib import Path

from predictor import resources

from predictor.annotate.genome_resolver import _blast, _run, fetch_region

_REPO = Path(__file__).resolve().parents[2]
_CURATED_GENOMES = _REPO.parent / "5.8 Promoters" / "arsr_merr_families" / "genomes"
# BLAST DBs must live in a SPACE-FREE path: BLAST+ 2.17 mis-splits the absolute DB path at the space
# in "5.11 Predictor" during its internal DB self-check. The DB is a build artifact, not source, so it
# lives under the user profile (override with $PREDICTOR_GENOME_DB).
_DB_DIR = Path(os.environ.get("PREDICTOR_GENOME_DB", resources.cache_path("genome_db")))
_EXTRA_DIR = _DB_DIR / "extra"
DEFAULT_DB = _DB_DIR / "arsr_merr_genomes"

# the dereplicated representative mirror (built by genome_mirror.py). Same space-free-path reasoning as
# the DB dir; override with $PREDICTOR_GENOME_MIRROR. genome_mirror imports this constant (one-way dep).
MIRROR_DIR = Path(os.environ.get("PREDICTOR_GENOME_MIRROR",
                                 resources.cache_path("genome_mirror")))


def curated_genome_fastas():
    return sorted(glob.glob(str(_CURATED_GENOMES / "*" / "*_genomic.fna")))


def mirror_genome_fastas():
    """Representative-mirror FASTAs (one .fna per accession, with a sibling .gff for context)."""
    return sorted(glob.glob(str(MIRROR_DIR / "*.fna")))


def extra_genome_fastas():
    return sorted(glob.glob(str(_EXTRA_DIR / "*.fna")))


def default_db():
    """Path prefix of the local DB if it exists (checks for a BLAST index file), else None."""
    for ext in (".nin", ".nal", ".ndb"):
        if Path(str(DEFAULT_DB) + ext).exists():
            return str(DEFAULT_DB)
    return None


def build_local_db(out_prefix=DEFAULT_DB, *, extra=True, title="arsr_merr local genomes"):
    """Concatenate curated (+ extra) genomic FASTAs and run makeblastdb (nucleotide)."""
    out_prefix = Path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fastas = curated_genome_fastas() + mirror_genome_fastas() + (extra_genome_fastas() if extra else [])
    if not fastas:
        raise FileNotFoundError(
            f"no genome FASTAs under {_CURATED_GENOMES}, {MIRROR_DIR} or {_EXTRA_DIR}")

    combined = out_prefix.parent / "_combined.fna"
    seen, n, skip = set(), 0, False
    with combined.open("w", encoding="ascii", errors="replace") as out:
        for fa in fastas:
            for line in Path(fa).read_text(encoding="ascii", errors="replace").splitlines():
                if line.startswith(">"):
                    acc = line[1:].split()[0]
                    skip = acc in seen                   # skip duplicate accessions across files
                    if not skip:
                        seen.add(acc)
                        n += 1
                if not skip:
                    out.write(line + "\n")

    # run from the output dir with relative names so a space in the repo path can't break BLAST args.
    # -parse_seqids enables `blastdbcmd -entry <accession>` (offline region extraction for homolog
    # inter-operon collection). It re-opens the DB by absolute path, so the DB MUST live in a
    # space-free dir (it does: ~/.predictor/genome_db) or BLAST+ 2.17 mis-splits the path at a space.
    _run([_blast("makeblastdb"), "-in", combined.name, "-dbtype", "nucl",
          "-out", out_prefix.name, "-parse_seqids", "-title", title.replace(" ", "_")],
         cwd=out_prefix.parent)
    combined.unlink(missing_ok=True)
    return str(out_prefix), n


def add_genome(accession: str):
    """Fetch a genome (FASTA) by NCBI accession into env/genome_db/extra/ and rebuild the DB."""
    _EXTRA_DIR.mkdir(parents=True, exist_ok=True)
    fa = _EXTRA_DIR / f"{accession}.fna"
    if not fa.exists():
        fa.write_text(fetch_region(accession, rettype="fasta"), encoding="ascii", errors="replace")
    return build_local_db()


def db_info(prefix=None):
    prefix = prefix or default_db()
    if not prefix:
        return "no local DB built yet"
    return _run([_blast("blastdbcmd"), "-db", prefix, "-info"]).stdout


def _main(argv):
    cmd = argv[0] if argv else "build"
    if cmd == "build":
        print(f"sources: {len(curated_genome_fastas())} curated + "
              f"{len(mirror_genome_fastas())} mirror + {len(extra_genome_fastas())} extra FASTA files")
        prefix, n = build_local_db()
        print(f"built {prefix} from {n} sequences")
        print(db_info(prefix))
    elif cmd == "info":
        print(db_info())
    elif cmd == "add":
        for acc in argv[1:]:
            prefix, n = add_genome(acc)
            print(f"added {acc}; DB now {n} sequences at {prefix}")
    else:
        print(__doc__)


if __name__ == "__main__":
    _main(sys.argv[1:])
