# External engines and settings

The Python package includes the reference data and the pyhmmer classifier. The search engines below are
separate programs. The installers set up every required one: `env/install_linux.sh` on Linux and in WSL
([guide](INSTALL_LINUX.md)), `env/install_windows.ps1` on Windows ([guide](INSTALL_WINDOWS.md)).
`tfop setup` shows what is ready, which settings file it read, and how to install anything missing.

| Engine | Purpose | Needed for | Resolution / setting | Linux / WSL | Windows |
|---|---|---|---|---|---|
| MMseqs2 | Homolog/SSN collection | every prediction | `MMSEQS`, vendored `env/tools/mmseqs`, PATH | conda (bioconda) | official `mmseqs-win64.zip` |
| BLAST+ | Locus/genome and BITACORA protein searches | every prediction | `BLAST_BIN`, PATH | conda (bioconda) | official NCBI `x64-win64` archive |
| NCBI Datasets | Genome accession acquisition | genome downloads | `DATASETS_EXE`, PATH | conda-forge `ncbi-datasets-cli` | official `datasets.exe` |
| pyhmmer | Packaged Pfam family classifier | every prediction | Python dependency | pip | pip |
| HMMER | BITACORA HMM search/clade models | genome discovery | PATH or `BITACORA_BIN` | conda (bioconda) | inside WSL only |
| Perl / patched BITACORA | Production candidate discovery | genome discovery | `BITACORA_DIR`, `BITACORA_BIN`, `BITACORA_WSL_DISTRO` | pinned installer checkout | inside WSL (`install_windows.ps1 -WithWsl`) |
| MetalNet2 | Optional metal-site evidence | optional | `METALNET_PYTHON`, `METALNET2_DIR` | separate environment (below) | separate environment (below) |
| ESMFold adapter | Optional folding | `--fold` only | `PREDICTOR_ESMFOLD_LIB`, `BIOHUB_TOKEN` | packaged adapter + token | packaged adapter + token |

Bioconda publishes no Windows builds, so on Windows the official upstream builds are used instead.
Native Linux uses native subprocesses, with `BITACORA_BIN` prepended to PATH.
Only Windows invokes WSL. The installers pin BITACORA to
`715950426b70b029279343dcc541115dfee20352` and copy `predictor/_engines/bitacora_patches/*.pl`.
`tfop setup` checks those files against the installed Scripts directory.

```bash
python -c "from predictor.discovery import bitacora; print(bitacora.ensure_patched())"
tfop setup
```

The repair command deliberately updates installed engine scripts; `tfop setup` itself only inspects.
Clade-model builds use private scratch directories. When a clade model cannot be built, warnings and
`[pfam_fallback]` profile tags expose the family-model fallback; investigate before comparing candidates.

## Settings file

Settings files accept unquoted or quoted values, spaces in quotes, and inline comments outside quotes.
They are read as UTF-8 (a byte-order mark is ignored). Read order: explicit `PREDICTOR_ENV_FILE`, the
checkout's `.env` (source and editable installs), then the per-user file
(`~/.config/tf-operator-predictor/.env` on Linux, `%LOCALAPPDATA%\tfop\tf-operator-predictor\.env` on
Windows); existing environment variables win. `.env.example` documents each key. Supported keys:

- `MMSEQS`, `BLAST_BIN`, `DATASETS_EXE`, `BITACORA_DIR`, `BITACORA_BIN`, `BITACORA_PY`, `BITACORA_WSL_DISTRO`;
- `METALNET_PYTHON`, `METALNET2_DIR`;
- `BIOHUB_TOKEN`, `NCBI_TOKEN`, `NCBI_API_KEY`, `NCBI_EMAIL`.

Path and output overrides include `PREDICTOR_CACHE_DIR`, `PREDICTOR_OUTPUT_DIR`,
`PREDICTOR_GENOME_MIRROR` and `PREDICTOR_ENV_FILE`; inspect their resolution with `tfop audit`.
Tokens are optional for public NCBI requests but may affect service rate limits; never share `.env`.
Folding and MetalNet2 absence remain explicit abstentions: the evidence source abstains and every
other verdict is unchanged.

`PREDICTOR_ALLOW_LEGACY_WORKSPACE=1` is for maintainers rebuilding the source references from the
original development corpora; it is off by default and ordinary runs use the packaged reference bytes.
It does not supply a missing deployed database.

## MetalNet2 (optional)

MetalNet2 (<https://github.com/wangchulab/MetalNet2>) predicts metal-binding sites from sequence and
alignment. The publication's canonical runs enabled it; without it, that one inducer evidence source
abstains. It needs its own Python 3.9 environment, because its classifiers are AutoGluon 0.8.0 models
pickled with numpy 1.25, pandas 1.5.3 and scikit-learn 1.2.2, plus about 6 GB of model weights.
The authors recorded and verified the steps below on Windows, where they reproduce upstream's shipped
example exactly. The same steps apply on Linux, but this release has not verified them there.

1. Check out the pinned commit into `env/tools/metalnet2` in the RIIO folder:

   ```bash
   git clone https://github.com/wangchulab/MetalNet2.git env/tools/metalnet2
   git -C env/tools/metalnet2 checkout 73841055b720c25289422e5480afe6f4b19a4dad
   ```

2. Create the environment and install the versions the models were written with (pin `spacy`: without
   it, pip tries to build a version that needs Python 3.10 and the whole install aborts):

   ```bash
   conda create -y -n metalnet python=3.9.16
   conda run -n metalnet python -m pip install "torch==1.13.1" "numpy==1.25.0" "pandas==1.5.3" \
       "scikit-learn==1.2.2" "scipy==1.10.1" "spacy==3.5.3" "fastai==2.7.12" "lightgbm==3.3.5" \
       "catboost==1.2" "xgboost==1.7.5" "autogluon.tabular==0.8.0" "fair-esm==2.0.0" "biopython==1.81" \
       "absl-py" "more-itertools" "networkx==3.1" "graphviz==0.20.1" "pyyaml" "gdown"
   ```

3. Download the two model bundles (about 1.9 GB) and unzip `model_1.zip` into
   `env/tools/metalnet2/model/train/` and `model_2.zip` into `env/tools/metalnet2/extra/train/`:

   ```bash
   conda run -n metalnet python -c "import gdown; gdown.download(id='1Qihh1iqIosVKZO0-OH800qj0eXAvCfMu', output='model_1.zip')"
   conda run -n metalnet python -c "import gdown; gdown.download(id='1GSOpS0AFeraaRUTWoFQS-hiVNh1p8Azr', output='model_2.zip')"
   ```

4. Download the four ESM checkpoints (about 4 GB) into the torch hub cache
   (`~/.cache/torch/hub/checkpoints`, or `%USERPROFILE%\.cache\torch\hub\checkpoints` on Windows) with
   `curl -L -O`: `esm_msa1b_t12_100M_UR50S.pt` and `esm2_t33_650M_UR50D.pt` from
   `https://dl.fbaipublicfiles.com/fair-esm/models/`, and `esm_msa1b_t12_100M_UR50S-contact-regression.pt`
   and `esm2_t33_650M_UR50D-contact-regression.pt` from `https://dl.fbaipublicfiles.com/fair-esm/regression/`.
   On Windows, torch's own downloader fails inside this Python 3.9 environment with an SSL error, so use
   `curl.exe`.

5. Record the interpreter in the settings file and check it: set `METALNET_PYTHON` to the environment's
   `python` (`.../envs/metalnet/bin/python` on Linux, `...\envs\metalnet\python.exe` on Windows; `conda
   env list` shows the folder), set `METALNET2_DIR` if the checkout is elsewhere, then run `tfop setup`
   (the MetalNet2 row reads `OK`) and `python -m predictor.structure.metalnet --self-test`.

The engine loads its models once per process, which dominates its runtime; per-sequence results are
cached. Alignments come, in order, from a cached a3m, an AF3 server a3m on disk, ColabFold (UniRef,
network) or a local MMseqs2 alignment against the packaged family members (offline, shallower;
recorded as `msa_source="local_family"`).
