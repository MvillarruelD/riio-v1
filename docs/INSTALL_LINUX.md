# Native Ubuntu installation

Supported target: Ubuntu 22.04/24.04. Install Git and a conda frontend (micromamba, mamba, or conda).
Clone the repository into a Linux filesystem; for WSL use a path under `~/`, not `/mnt/c`, for best I/O behavior.

```bash
git clone https://github.com/MvillarruelD/riio-v1.git
cd riio-v1
bash env/install_linux.sh --with-riio
conda activate riio-engines
tfop setup
tfop audit
tfop selftest
```

With micromamba use `micromamba activate riio-engines`. The installer also supports `--miniforge`
when no conda frontend is available; it downloads the Miniforge installer into your home directory.
No sudo is needed for the engine environment. The conda specification includes Python 3.12,
HMMER, BLAST+, MMseqs2, NCBI Datasets, Perl and Git. BITACORA is cloned at
`715950426b70b029279343dcc541115dfee20352`, then the packaged coverage patches are applied.
`tfop setup` validates the installed patch contents.

For plots, install `fonts-liberation` through the Ubuntu package manager when permitted.
Matplotlib uses available sans-serif fallbacks. Graphical desktop services are unnecessary for reports;
`MPLBACKEND=Agg` selects a headless backend. The browser GUI runs locally with `tfop-gui`.

The source installer writes `.env` in the checkout. Installed packages resolve settings using
`PREDICTOR_ENV_FILE`, then checkout settings, then the per-user configuration directory reported by
`tfop audit`. Environment variables take precedence over settings files.

```bash
tfop genome --fasta tests/data/mini_genome.fna --gff tests/data/mini_genome.gff \
  --out runs/discovery-smoke --discover-only --offline
tfop batch proteins.faa --organism GCF_000005845.2 --out runs/proteins
```

The first command exercises real discovery without network after installation. Complete prediction
also needs homolog/promoter context; `--offline` does not manufacture missing caches.

Troubleshooting: rerun the installer to repair missing tools/patches, inspect `tfop setup`, confirm
`BITACORA_BIN` points to the environment's native `bin` directory, and inspect `logs/` for per-TF failures.
Profile-build fallback is recorded as `[pfam_fallback]`; it should be investigated before comparing discovery
manifests. MetalNet2 and folding are optional, separate installs; see [ENGINES](ENGINES.md).

## Verify real engines and complete predictions

The GitHub Actions `e2e-linux` job runs on every push and pull request, alongside Ubuntu and Windows
unit tests. It installs native engines, pins/patches BITACORA, runs `tfop setup` and `tfop selftest`,
then invokes the same helper available in a source checkout:

```bash
python tools/linux_e2e.py --out /tmp/riio-e2e-fresh --timeout 1200
```

Use a fresh output directory. The helper first runs real BITACORA on the bundled K-12 fragment with
network disabled. It then downloads the unmodified complete K-12 assembly `GCF_000005845.2`, discovers
its candidates and identifies zntR by the exact fixture protein sequence. Complete single-protein,
protein-list, local-genome and genome-accession predictions run with this full context. Genome routes
select only that known positive while retaining every discovered candidate in the reports.

The full assembly matters: the small discovery fixture excludes distant regulatory targets and is
insufficient for an offline operator prediction. Complete integration therefore **requires network**
for NCBI genome/homolog-promoter retrieval. Network failures fail the job; they are not converted into
successful skips. Optional folding/MetalNet2 are disabled, while real discovery, homolog searches,
operator prediction, reports and AF3 serialization remain enabled.

Every completed bundle is checked by `tools/validate_bundle.py`; HTML/XLSX links, retained raw tables,
and nonempty apo/operator jobs in both AF3 formats are validated. Command deadlines bound child
process groups. `integration_result.json`, command logs and per-TF logs record failures, and CI
uploads the output directory even when a command fails. The job has a 90-minute overall deadline.
