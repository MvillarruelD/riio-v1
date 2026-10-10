# Distributing tf-operator-predictor

Release and installation policy for the private 0.3 series.

## How scientific software is usually distributed

For a Python bioinformatics tool like this, five patterns dominate, roughly in order of how self-contained
they are:

1. **PyPI + `pip`** — `pip install toolname`. The default for Python tools. Reference data is either bundled
   as *package data* (small) or downloaded on first run (large). Console entry points make it "run as a
   program". Its cousin, **`pip install git+https://…`**, installs the same way straight from a repo and is
   the norm for unpublished / pre-print / private research code.
2. **`pipx`** — `pipx install toolname` (or `pipx install "git+…"`). Installs a *CLI application* into its own
   isolated environment and puts its commands on your PATH. This is the modern recommended way to hand
   someone a Python program they just *run*, with nothing to activate.
3. **Bioconda / conda** — `conda install -c bioconda toolname`. The de-facto standard for published
   bioinformatics tools, because conda also packages the **non-Python engines** (MMseqs2, BLAST+, HMMER,
   Foldseek are all on Bioconda). One command pulls the tool *and* its whole toolchain.
4. **Containers (Docker / Apptainer-Singularity)** — `docker run …` / `apptainer run …`. The reproducibility
   gold standard for multi-engine pipelines: the image pins Python, every engine, and every version. This is
   how nf-core, Galaxy tools, and most "complex pipeline" papers ship.
5. **Standalone executables** (PyInstaller, `shiv`, `zipapp`) — a single file that runs without a Python
   install. Rare for bioinformatics: fragile with compiled wheels (numpy, pyhmmer) and useless for the
   external engines. Not recommended here.

## What fits this tool

`tf-operator-predictor` is a Python CLI that shells out to optional external engines (MMseqs2, BLAST+,
Foldseek, NCBI datasets and FoldX) and ships versioned immutable reference data inside `predictor/data/`.
It is currently **private / unpublished, shared with named collaborators**.

Given that, the fit is:

- **Now (private, unpublished): install a wheel directly from GitHub with `pip` or `pipx`.** An editable
  checkout remains available to developers, but is not required to locate reference data.
- **When published: Bioconda + a Docker image.** Bioconda handles the external engines (the hard part of
  bioinformatics distribution); the container gives bit-reproducible runs. Both need the tool public and a
  recipe / Dockerfile.
- **PyPI** becomes appropriate after the project is public and the licence permits redistribution.

## Wheel contract

- `predictor/data/` is read-only, versioned release data. `MANIFEST.json` records hashes, sizes and
  record counts; JSON schemas document the main resources.
- Installing the wheel installs the complete core reference database. Release 0.3.0 is about 34 MiB
  compressed and 141 MiB installed; there is no separate user-facing database download or build step.
- Query genomes are fetched and cached on demand. HomoDB and MetalNet2 weights remain optional external
  assets and must not be confused with the core runtime database.
- Generated databases and network/engine caches use the platform user cache directory reported by
  `tfop audit` (override with `PREDICTOR_CACHE_DIR`).
- User-facing results default to `./results` (override with `PREDICTOR_OUTPUT_DIR`).
- Source material, historical derivation artifacts and local verification evidence are archived
  outside the active checkout and remain excluded from wheels and production lookup.
- The wheel must pass `tfop audit`, `tfop selftest` and a representative `tfop scan` from a working
  directory outside the checkout before release.
- The optional `gui` extra installs the one-page Streamlit workspace and the `tfop-gui` entry point; it
  calls the same public API as the CLI and embeds the same generated `REPORT.html`.

Release 0.3.0 was verified by installing its wheel without editable mode into a clean Python 3.12
environment. The installed CLI found all 71 SSN clusters and 285,672 accession mappings, scanned the
packaged curated TF FASTA, and completed the self-test with no failures.

## Release checklist

For the narrower generic checks that do not read shipped reference records, use
[TESTING.md](TESTING.md). Those checks are intentionally distinct from this full
release and data-validation workflow.

The shipped-data integration assertions are a separate, explicit gate because the
default pytest configuration excludes them:

```bash
python -m pytest -m released_data
```

The artifact, isolated-install and runtime release workflow is:

```bash
python tools/release_check.py
```

`tools/release_check.py` checks the packaged-data manifest, compiles, lints, runs the default
non-`released_data` suite and the separate `released_data` integration gate,
builds **both** a wheel and a source distribution, creates a clean virtual environment, installs the
wheel **non-editably**, changes to a directory outside the checkout, and only then runs
`tfop --version`, `tfop audit`, `tfop selftest` and a packaged-data `tfop scan`. Finally it verifies
that every file listed in `predictor/data/MANIFEST.json` exists in the *installed* package and that
nothing was written inside site-packages. About 7 minutes.

Running from outside the checkout is the whole point. Inside it, the source tree satisfies imports
and relative paths that a real installation would not, so a check that passes there proves nothing —
which is exactly how two packaged PDFs stayed listed in the release manifest while being absent from
every wheel.

It needs the `build` frontend, which is in the dev extra:

```bash
python -m pip install -e ".[dev]"
```

The individual steps, if you need to run one on its own:

```bash
python tools/build_data_manifest.py --check
python -m ruff check predictor tests tools
python -m pytest
python -m pytest -m released_data
python -m build --wheel --sdist
```

Then install the wheel in a clean environment outside the repository and run:

```bash
tfop --version
tfop audit
tfop selftest --timeout 30 -v
tfop scan /path/to/curated_tfs.faa
```
