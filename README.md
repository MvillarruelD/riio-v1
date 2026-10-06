# RIIO v1

**Regulator, Inducer & Operator inference**, named with a nod to Río de la Plata.

RIIO takes a bacterial transcription-factor protein or genome and produces operator, motif,
inducer and regulon predictions, with evidence and provenance in a portable `REPORT.html` bundle.
This repository contains the active RIO software and the production analysis companion used by the
project. Version 1.0.0 preserves the `tfop` commands and `predictor` Python API.

## Download and install

Use Python **3.11 or 3.12** for the complete repository workflow. The predictor package also supports
Python 3.10. Linux is the simplest platform for installing the external bioinformatics tools.

```bash
git clone https://github.com/MvillarruelD/riio-v1.git
cd riio-v1
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install ".[genome,gui]"
```

On Windows PowerShell, activate with `.venv\Scripts\Activate.ps1` instead. If activation is disabled,
use `.venv\Scripts\python.exe -m pip` and `.venv\Scripts\tfop.exe` directly.

Alternatively, download the wheel from [Releases](https://github.com/MvillarruelD/riio-v1/releases)
and install it with `python -m pip install "./tf_operator_predictor-1.0.0-py3-none-any.whl[genome,gui]"`.
The wheel contains the predictor and its reference database. The source archive/Git checkout also
contains the analysis companion and canonical candidate inputs.

An installation directly from the frozen release is also available:

```bash
python -m pip install "tf-operator-predictor[genome,gui] @ git+https://github.com/MvillarruelD/riio-v1.git@v1.0.0"
```

## Install the search tools

Python installation does not install MMseqs2, NCBI BLAST+, or the NCBI Datasets executable.
With conda or mamba:

```bash
conda create -n riio-tools -c conda-forge -c bioconda mmseqs2 blast ncbi-datasets-cli
conda activate riio-tools
```

Run the RIIO commands with these tools available on `PATH`, or configure their executable paths as
described in [engine setup](docs/ENGINES.md). Windows MMseqs2 needs its upstream Windows distribution;
BITACORA needs Linux or WSL. Each tool retains its own license.

```bash
tfop --version
tfop audit
tfop setup
tfop selftest --timeout 120
```

`audit` checks the packaged database. `setup` identifies missing tools. A skipped optional engine
means that evidence source is unavailable; installing only Python is not a complete engine setup.

## Run the pipeline

Launch the graphical workspace:

```bash
tfop-gui
```

Or use the command line:

```bash
# Protein accession, amino-acid sequence, or a single-record protein FASTA
tfop predict WP_000240047.1 --name MyTF --organism GCF_000005845.2 --out results/MyTF

# Candidate census from a genome or proteome
tfop scan proteins.faa --out results/census

# Genome-first workflow; inspect all options before starting a large batch
tfop genome --help
tfop genome GCF_000005845.2 --discover-only --out results/genome
```

Open the resulting `REPORT.html` in a browser. Reports retain machine-readable predictions,
evidence, input/export links, and checksums. The core reference database is included automatically;
query genomes and online annotations are fetched as needed. Network access is needed for that route.

Credentials are optional for installation. Copy `.env.example` to `.env` and fill in the services
you use. Set `PREDICTOR_ENV_FILE` to its absolute path when running an installed wheel from another
directory. Never commit `.env`. NCBI requests work without a key at the public rate limit.

MetalNet2 and BITACORA are separately installed integrations. The canonical paper runs enabled
MetalNet2 and used BITACORA discovery. Follow [the pinned engine setup](docs/ENGINES.md) to reproduce
that configuration. Folding is optional and off by default.

## Reproduce the project workflow

The [production guide](docs/PRODUCTION.md) documents the three canonical candidate manifests,
238 accompanying protein FASTAs, portable output roots, batch/resume commands, and downstream
analysis. It distinguishes this release from the paper's earlier frozen predictor commit.

```bash
python -m pip install -e ".[dev,genome,gui]" -r analysis/requirements.txt
python -m pytest
python -m pytest -m released_data
python tools/build_data_manifest.py --check
python tools/release_check.py --selftest-timeout 120
```

The software checks do not rerun the full multi-hour online prediction batch. See
[the release audit](docs/RELEASE_AUDIT.md) for exactly what was checked.

## Repository contents

| Path | Purpose |
|---|---|
| `predictor/` | Active CLI, GUI, API, inference modules and packaged reference database |
| `predictor/_engines/` | Owned engine runners and required adapter helpers |
| `analysis/` | Production stage drivers, canonical manifests, input FASTAs and contract tests |
| `analysis/benchmarking/` | Active canonical scoring helpers and required reference inputs |
| `reports/` | Current report/figure generation programs |
| `tests/`, `tools/` | Software tests, data integrity checks and distribution verification |
| `docs/`, `LICENSES/` | Installation, reproducibility, source provenance and upstream notices |

Generated runs, private credentials, environments, binaries/model weights, historical experiments
and manuscript documents are excluded. The optional composite-figure stage skips manuscript artwork
when it is absent. Core predictions and standalone figures do not require that artwork.

## Credit and citation

See [AUTHORS.md](AUTHORS.md), [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md), and
[CITATION.cff](CITATION.cff). The final software author list and paper DOI will be updated by the
project authors. Cite the exact software version used; no archived DOI has been assigned yet.

Original RIIO code is released under [MIT](LICENSE). The included BITACORA-derived Perl patch is
GPL-3.0-only, separately identified in `LICENSES/`; upstream data and engines retain their terms.
Predictions are computational hypotheses for experimental prioritization.
