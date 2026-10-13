# RIIO v1.1

**Regulator, Inducer & Operator inference**.

RIIO takes a bacterial transcription-factor protein or genome and produces operator, motif,
inducer and regulon predictions, with evidence and provenance in a portable `REPORT.html` bundle.
This repository contains the active RIIO software and the production analysis companion used by the
project. Version 1.1.0 preserves the `tfop` commands and `predictor` Python API.

## Quick start

RIIO runs on Linux and on Windows 10/11. One installer per system sets up RIIO, its search engines and
its settings. The step-by-step guides assume no previous experience with the command line tools involved.

**Linux, or Windows through WSL (recommended):** [installation and first steps](docs/INSTALL_LINUX.md)

```bash
git clone https://github.com/MvillarruelD/riio-v1.git
cd riio-v1
bash env/install_linux.sh --miniforge --with-riio
conda activate riio-engines        # or: source ~/miniforge3/bin/activate riio-engines
```

**Windows, natively in PowerShell:** [installation and first steps](docs/INSTALL_WINDOWS.md), which also
covers WSL. Install Python 3.12 and Git first.

```powershell
git clone https://github.com/MvillarruelD/riio-v1.git
cd riio-v1
powershell -ExecutionPolicy Bypass -File env\install_windows.ps1
.\.venv\Scripts\Activate.ps1
```

Then check the installation and predict the operators, inducer and regulon of *E. coli* ZntR, whose
sequence ships with RIIO (3 to 5 minutes; needs internet access):

```bash
tfop setup
tfop selftest
tfop predict tests/data/mini_tf.faa --name ZntR --organism GCF_000005845.2 --out runs/ZntR
```

Open `runs/ZntR/jobs/ZntR/REPORT.html` in a browser. Predictions are computational hypotheses for
experimental prioritization.

## Workflows

The CLI, Python API and browser GUI use the same resumable per-protein runner. Start the GUI with
`tfop-gui`, then open <http://localhost:8501>.

| Starting input | CLI | GUI | Needs |
|---|---|---|---|
| One protein sequence, FASTA or accession | `tfop predict protein.faa --name MyTF --organism GCF_...` | Single TF | MMseqs2, BLAST+, NCBI Datasets, internet |
| Multi-record protein FASTA or accession TSV | `tfop batch proteins.faa --organism GCF_... --out runs/list` | Protein list | as above |
| Genome accession | `tfop genome GCF_000005845.2 --out runs/accession` | Genome scan / accession | as above, plus BITACORA |
| Local nucleotide genome and annotation | `tfop genome --fasta genome.fna --gff genome.gff --out runs/local` | Genome scan / local files | as above, plus BITACORA |
| Quick Pfam census of a proteome or genome | `tfop scan proteome.faa` | Genome scan / Quick census | nothing else (offline) |

`--organism` is the NCBI assembly accession of the protein's source genome. Local genome context is
copied into the run and reaches each selected protein prediction. GenBank and proteome inputs are also
supported. Complete predictions from a protein list need an organism or resolvable source genome for
genomic context. Changed sequences, local files or prediction settings require a new output directory.
`tfop <command> --help` explains every option; see also [entry modes and API examples](docs/ENTRY_POINTS.md).

Runs include per-protein `jobs/*/REPORT.html`, `discovery/DISCOVERY_REPORT.html`, `RUN_REPORT.html`
and `Supplementary_Table_S1.xlsx`/`.tsv`. [Report rebuilding](docs/REPORTS.md) reads existing results
with `tfop report runs/local`. [Optional AF3 exports](docs/AF3.md) use
`tfop af3 runs/local --format both`; they do not affect predictions. Keep the complete run tree
when sharing so relative HTML and spreadsheet links remain usable.

## Other ways to install

- **Your own Python environment** (Python 3.10 or newer): install from the RIIO folder in editable
  mode, `python -m pip install -e ".[genome,gui]"`, so that RIIO reads the folder's settings file. The
  search engines are separate programs: `tfop setup` lists any that are missing, with install
  instructions for your system, and [engine setup](docs/ENGINES.md) describes them all.
- **Conda on Linux or macOS:** in the RIIO folder, `conda env create -f environment.yml` and
  `conda activate riio` install RIIO with MMseqs2, BLAST+, NCBI Datasets, HMMER and Perl. Genome
  discovery also needs the pinned BITACORA checkout: `ENV_NAME=riio bash env/install_linux.sh` adds it.
- **Wheel from [Releases](https://github.com/MvillarruelD/riio-v1/releases):**
  `python -m pip install "./tf_operator_predictor-1.1.0-py3-none-any.whl[genome,gui]"`. The wheel
  contains the predictor and its reference database; the source archive/Git checkout also contains the
  analysis companion, the installers and the canonical candidate inputs. An installed wheel keeps its
  settings in a per-user file, which `tfop setup` names.
- **From the release tag**, once v1.1.0 is tagged:
  `python -m pip install "tf-operator-predictor[genome,gui] @ git+https://github.com/MvillarruelD/riio-v1.git@v1.1.0"`.

## Settings and optional engines

The installers record the engine locations in `.env` in the RIIO folder; `tfop setup` shows the settings
file in use and `.env.example` documents every key. Credentials are optional: NCBI requests work
without a key at the public rate limit, and an NCBI API key (`NCBI_TOKEN`) makes large runs faster.
Never commit `.env`.

MetalNet2 and BITACORA are separately installed integrations. The canonical paper runs enabled
MetalNet2 and used BITACORA discovery. Follow [the pinned engine setup](docs/ENGINES.md) to reproduce
that configuration. Folding is optional and off by default. Each external tool retains its own license.

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
[the v1.1 preparation audit](docs/RELEASE_AUDIT_V1_1.md) for exactly what was checked.

## Repository contents

| Path | Purpose |
|---|---|
| `predictor/` | Active CLI, GUI, API, inference modules and packaged reference database |
| `predictor/_engines/` | Owned engine runners and required adapter helpers |
| `env/` | Linux and Windows installers and the Linux engine environment |
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
