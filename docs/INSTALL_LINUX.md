# Installing and using RIIO on Linux

This guide assumes no previous experience with conda or Python packaging. Copy each command into a
terminal, press Enter, and wait for it to finish before typing the next one. The same steps work inside
WSL Ubuntu on Windows (see [the Windows guide](INSTALL_WINDOWS.md)).

**You need:** a 64-bit Linux computer (tested on Ubuntu 22.04, 24.04 and 26.04), an internet connection
and about 3 GB of free disk space. Administrator (`sudo`) rights are needed only if `git` or `curl` is
missing. Everything else installs inside your home folder.

## 1. Check the two prerequisite tools

```bash
git --version
curl --version
```

If either answers `command not found`, install both (this asks for your password):

```bash
sudo apt update && sudo apt install -y git curl
```

## 2. Download RIIO

```bash
cd ~
git clone https://github.com/MvillarruelD/riio-v1.git
cd riio-v1
```

Keep RIIO in your Linux home folder. In WSL, do not use a folder under `/mnt/c`: it works, but every
file operation is much slower there.

## 3. Run the installer

```bash
bash env/install_linux.sh --miniforge --with-riio
```

It takes about 5 to 15 minutes, mostly downloads. It needs no `sudo` and installs, in your home folder:

- Miniforge, a small conda distribution, into `~/miniforge3` (skipped if you already have conda, mamba
  or micromamba; `--miniforge` then does nothing);
- a conda environment named `riio-engines` with Python 3.12, MMseqs2, NCBI BLAST+, HMMER, NCBI Datasets
  and Perl;
- BITACORA, the genome-discovery engine, pinned to the exact commit the publication used, in
  `env/tools/bitacora`, with RIIO's patches applied;
- RIIO itself, into the same environment;
- a settings file, `.env` in the RIIO folder, recording where each tool lives.

It ends with `tfop setup`, a readiness table in which every engine should say `OK`. The optional
MetalNet2 row may say `MISSING`; RIIO then skips that one evidence source (see [ENGINES](ENGINES.md)).
If the installer stops with an error, read its last lines, fix the cause and run it again: every step
can be repeated safely.

## 4. Activate the environment

Do this **in every new terminal** before using RIIO:

```bash
conda activate riio-engines
```

If the terminal answers `conda: command not found`, or asks you to run `conda init`, use this instead:

```bash
source ~/miniforge3/bin/activate riio-engines
```

To make the short form work from now on, run `~/miniforge3/bin/conda init bash` once and open a new
terminal. With micromamba, use `micromamba activate riio-engines`. The prompt now starts with
`(riio-engines)`; `cd ~/riio-v1` returns to the RIIO folder.

## 5. Check the installation

```bash
tfop selftest
```

Every module checks itself, in one to three minutes. A complete installation ends with a line such as
`48 ok, 0 skipped (missing engine/network), 0 timed out (>300s), 0 failed`. A `skip` means a module
needs an engine or network access that is not available; `tfop setup` shows which.

## 6. Make a first prediction

RIIO includes the *E. coli* ZntR protein as a test input. Predict its operators, inducer and regulon
in the context of the *E. coli* K-12 genome:

```bash
tfop predict tests/data/mini_tf.faa --name ZntR --organism GCF_000005845.2 --out runs/ZntR
```

This takes about 3 to 5 minutes and needs internet access: RIIO downloads the genome and homolog
promoters from NCBI and caches them for later runs. The command ends with a list of the files it wrote.
Open the report in a browser:

```bash
xdg-open runs/ZntR/jobs/ZntR/REPORT.html
```

In WSL, open it with Windows instead:

```bash
explorer.exe "$(wslpath -w runs/ZntR/jobs/ZntR/REPORT.html)"
```

`operators/FINAL_operators.tsv` holds the final operator list and `operators_ranked.json` the same in a
machine-readable form. All results are computational predictions to be tested experimentally.

## 7. Analyse your own proteins and genomes

| You have | Command |
|---|---|
| One protein: a sequence, a one-record FASTA, or an NCBI/UniProt accession | `tfop predict MyTF.faa --name MyTF --organism GCF_xxxxxxxxx.x --out runs/MyTF` |
| Several proteins: a multi-record FASTA, or a TSV with an `accession` column | `tfop batch proteins.faa --organism GCF_xxxxxxxxx.x --out runs/list` |
| A genome accession | `tfop genome GCF_000005845.2 --out runs/ecoli` |
| A genome on disk: FASTA + GFF3, or GenBank | `tfop genome --fasta genome.fna --gff genome.gff --out runs/local` |

`--organism` is the NCBI assembly accession (`GCF_...` or `GCA_...`) of the genome the protein comes
from; it is shown on the organism's NCBI Assembly or Datasets page. Without it, RIIO tries to find the
genome from the protein sequence. Add `--discover-only` to `tfop genome` to list the candidate
regulators first (`runs/ecoli/manifest.csv` and `runs/ecoli/discovery/DISCOVERY_REPORT.html`), then
predict a subset with `--select`. Every command explains its options with `--help`, for example
`tfop genome --help`. More entry modes are described in [ENTRY_POINTS](ENTRY_POINTS.md).

Each run folder holds one `jobs/<name>/REPORT.html` per protein, plus `RUN_REPORT.html` and
`Supplementary_Table_S1.xlsx` for multi-protein runs. Keep a run folder whole when you copy or share it:
its reports link to each other with relative paths. Use a new `--out` folder whenever you change the
input or the options; RIIO refuses to mix different inputs in one run folder.

## 8. Use the graphical interface (optional)

```bash
tfop-gui
```

Then open <http://localhost:8501> in a web browser (from WSL, the Windows browser works). The interface
offers the same single-protein, protein-list and genome-scan workflows. Press `Ctrl+C` in the terminal
to stop it. It accepts connections from this computer only; to serve it to other machines on purpose,
start it with `tfop-gui --server.address=0.0.0.0`.

## Optional: an NCBI API key

RIIO works without one, but NCBI limits anonymous requests. A free key, created in your NCBI account
settings, makes large runs faster and more reliable. Add it to the `.env` file in the RIIO folder:

```bash
nano .env
```

Add the lines `NCBI_TOKEN=your-key` and `NCBI_EMAIL=you@example.org`, save with `Ctrl+O` then Enter,
and leave with `Ctrl+X`. `tfop setup` then shows `NCBI_TOKEN=set`. The file `.env.example` describes
every setting.

## Updating and uninstalling

To update, run `git pull` in the RIIO folder, then the installer again. To uninstall, remove the
environment with `conda env remove -n riio-engines`, then delete the RIIO folder and the download
cache `~/.cache/tf-operator-predictor`. Miniforge stays in `~/miniforge3` for other projects; delete
that folder too if you do not need it.

## Troubleshooting

| Message or symptom | What to do |
|---|---|
| `conda: command not found` | Use `source ~/miniforge3/bin/activate riio-engines` (step 4). |
| `tfop: command not found` | The environment is not active in this terminal: repeat step 4. |
| `python: command not found`, or `ensurepip is not available` | Only the manual route below uses the system Python; the installer brings its own. |
| A `MISSING` row in `tfop setup` | Run the installer again; it repairs a partial installation. |
| `INPUT ERROR: [Errno 36] File name too long` | Fixed in RIIO 1.1 (proteins longer than 255 amino acids); update RIIO. |
| `could not classify the query as a supported TF family` | The protein is not a transcription factor of a family RIIO supports. |
| Network errors, or HTTP 429 from NCBI | Check the connection, wait a minute and run the same command again; an NCBI API key helps. |
| `belongs to different input; choose a new output directory` | The input or options changed since that run: use a new `--out` folder. |
| Anything else | Run `tfop setup` and `tfop audit`, and read the per-protein logs in `runs/<name>/logs/`. |

## Manual installation (advanced)

If you manage Python environments yourself, install the engines and RIIO separately.
`bash env/install_linux.sh` without `--with-riio` installs only the engines and BITACORA. Then, in a
Python 3.10 or newer environment (on Ubuntu, `sudo apt install python3-venv` first):

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[genome,gui]"
```

Install from the RIIO folder in editable mode (`-e`) as shown: that is how RIIO finds the folder's
`.env`. A package installed any other way (for example the wheel from the Releases page) reads the
per-user settings file `~/.config/tf-operator-predictor/.env`, or the file named by
`PREDICTOR_ENV_FILE`. `tfop setup` prints which settings file it read; [ENGINES](ENGINES.md) lists every
setting. For plots, `fonts-liberation` gives the intended typeface; `MPLBACKEND=Agg` runs without a display.

## For developers: verify real engines and complete predictions

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

Discovery with the installed engines can also be checked offline, on the bundled genome fragment:

```bash
tfop genome --fasta tests/data/mini_genome.fna --gff tests/data/mini_genome.gff \
  --out runs/discovery-smoke --discover-only --offline
```

A profile-build fallback is recorded as `[pfam_fallback]`; investigate it before comparing discovery
manifests. MetalNet2 and folding are optional, separate installs; see [ENGINES](ENGINES.md).
