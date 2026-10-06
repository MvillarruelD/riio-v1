# External engine setup

Run `tfop setup` after installing RIIO. It reports the resolver used by the pipeline. Use `MMSEQS`
and `DATASETS_EXE` for tools outside `PATH`; BLAST+ tools should be on `PATH`. `pyhmmer` is installed
with the predictor, and does not substitute for the HMMER executables needed by BITACORA.

## BITACORA

The canonical source pin is `715950426b70b029279343dcc541115dfee20352`.

```bash
git -c core.autocrlf=false clone https://github.com/molevol-ub/bitacora bitacora
git -C bitacora checkout 715950426b70b029279343dcc541115dfee20352
conda create -n bitacora -c conda-forge -c bioconda hmmer "blast>=2.12"
export BITACORA_DIR="$PWD/bitacora"
python -m predictor.discovery.bitacora --check
```

Run with HMMER/BLAST on `PATH` on Linux. Windows uses WSL: set `BITACORA_WSL_DISTRO` and
`BITACORA_BIN` to the distribution and directory containing the WSL executables. Configure these
through real environment variables or the credentials/location file. Preserve LF line endings.

RIIO carries the original pipeline's GPL-3.0-only Perl patch. It makes coverage fractions configurable
and is copied into the user's BITACORA checkout by the existing adapter. The pipeline's established
coverage defaults are unchanged. See the packaged patch README and full upstream license.
The engine itself is installed separately. The runner and patch now survive a wheel installation.

## MetalNet2

The canonical source pin is `73841055b720c25289422e5480afe6f4b19a4dad`.

```bash
git clone https://github.com/wangchulab/MetalNet2 metalnet2
git -C metalnet2 checkout 73841055b720c25289422e5480afe6f4b19a4dad
conda create -n metalnet python=3.9.16
conda activate metalnet
python -m pip install "torch==1.13.1" "numpy==1.25.0" "pandas==1.5.3" "scikit-learn==1.2.2" \
  "scipy==1.10.1" "spacy==3.5.3" "fastai==2.7.12" "lightgbm==3.3.5" "catboost==1.2" \
  "xgboost==1.7.5" "autogluon.tabular==0.8.0" "fair-esm==2.0.0" "biopython==1.81" \
  absl-py more-itertools "networkx==3.1" "graphviz==0.20.1" pyyaml
```

Follow the pinned upstream repository's model-download instructions: both AutoGluon model stacks
and the ESM checkpoints are required. These multi-GB weights are not in RIIO downloads. Preserve
their upstream directory layout (`model/train/` and `extra/train/`). Set `METALNET2_DIR` to the
absolute engine checkout and `METALNET_PYTHON` to the dedicated environment's Python executable.
Return to the RIIO environment and run:

```bash
python -m predictor.structure.metalnet --status
tfop setup
```

The engine uses an older isolated dependency stack because its serialized models require it.
The owned subprocess runner is packaged; availability of that runner is distinct from availability
of the external models. Missing MetalNet2 is reported as unavailable evidence.

## Optional folding

The active folding adapter is packaged, but folding remains off by default. Biohub credentials and
service access are external; configure `BIOHUB_TOKEN` only if using `--fold`. No optional structure
retrieval experiments, FoldX, HomoDB, or their installation programs are included.
