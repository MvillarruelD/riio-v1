#!/usr/bin/env bash
# install_linux.sh -- install RIIO's external engines on Linux (tested on Ubuntu 22.04/24.04/26.04).
#
# What it does, in order (every step is idempotent, so re-running it repairs a partial install):
#   1. finds micromamba, mamba or conda (or installs Miniforge into ~/miniforge3 with --miniforge)
#   2. creates/updates the conda env `riio-engines` from env/engines-linux.yml
#      (HMMER, BLAST+, MMseqs2, NCBI Datasets, Perl)
#   3. clones BITACORA at the pinned commit into env/tools/bitacora (or $BITACORA_DIR)
#   4. writes the engine locations into the settings file (.env) the pipeline reads
#   5. optionally (--with-riio) pip-installs RIIO itself, with the genome + GUI extras, into that env
#   6. runs `tfop setup` to show what is ready
#
# It never needs sudo. MetalNet2 (optional, ~2 GB of model weights) is not installed here; see
# docs/ENGINES.md.
#
# Usage:
#   bash env/install_linux.sh                 # engines only (you install the Python package yourself)
#   bash env/install_linux.sh --with-riio     # engines + RIIO installed in the same env (one env)
#   bash env/install_linux.sh --miniforge     # also install Miniforge if no conda is found
#   ENV_NAME=myenv bash env/install_linux.sh  # a different conda env name
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_NAME="${ENV_NAME:-riio-engines}"
BITACORA_URL="https://github.com/molevol-ub/bitacora"
# The commit every published run used. Change it only together with a re-validation of discovery.
BITACORA_COMMIT="715950426b70b029279343dcc541115dfee20352"
BITACORA_DIR="${BITACORA_DIR:-$REPO/env/tools/bitacora}"
WANT_RIIO=0
WANT_MINIFORGE=0
for arg in "$@"; do
  case "$arg" in
    --with-riio|--venv) WANT_RIIO=1 ;;
    --miniforge) WANT_MINIFORGE=1 ;;
    -h|--help) sed -n '2,22p' "$0"; exit 0 ;;
    *) echo "unknown option: $arg" >&2; exit 2 ;;
  esac
done

say() { printf '\n==> %s\n' "$*"; }

# ---------------------------------------------------------------- 1. a conda front-end
find_conda() {
  for c in micromamba mamba conda; do
    if command -v "$c" >/dev/null 2>&1; then echo "$c"; return; fi
  done
  for c in "$HOME/miniforge3/bin/mamba" "$HOME/miniforge3/bin/conda" "$HOME/mambaforge/bin/mamba" \
           "$HOME/miniconda3/bin/conda" "$HOME/anaconda3/bin/conda"; do
    if [ -x "$c" ]; then echo "$c"; return; fi
  done
}
CONDA="$(find_conda || true)"
if [ -z "$CONDA" ]; then
  if [ "$WANT_MINIFORGE" = 1 ]; then
    say "Installing Miniforge into ~/miniforge3"
    arch="$(uname -m)"; os="$(uname -s)"
    installer="$(mktemp -d)/miniforge.sh"
    curl -fsSL -o "$installer" \
      "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-${os}-${arch}.sh"
    bash "$installer" -b -p "$HOME/miniforge3"
    CONDA="$HOME/miniforge3/bin/mamba"
  else
    echo "No micromamba/mamba/conda found. Install Miniforge (https://conda-forge.org/download/)" >&2
    echo "or re-run with --miniforge to have this script install it into ~/miniforge3." >&2
    exit 1
  fi
fi
say "Using $CONDA"

# ---------------------------------------------------------------- 2. the engine env
if "$CONDA" env list 2>/dev/null | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  say "Updating conda env $ENV_NAME"
  "$CONDA" env update -y -n "$ENV_NAME" -f "$REPO/env/engines-linux.yml" 2>/dev/null \
    || "$CONDA" env update -n "$ENV_NAME" -f "$REPO/env/engines-linux.yml"
else
  say "Creating conda env $ENV_NAME (HMMER, BLAST+, MMseqs2, NCBI Datasets, Perl)"
  "$CONDA" env create -y -n "$ENV_NAME" -f "$REPO/env/engines-linux.yml" 2>/dev/null \
    || "$CONDA" env create -n "$ENV_NAME" -f "$REPO/env/engines-linux.yml"
fi
ENV_PREFIX="$("$CONDA" env list | awk -v n="$ENV_NAME" '$1==n {print $NF}')"
[ -n "$ENV_PREFIX" ] && [ -d "$ENV_PREFIX/bin" ] || { echo "could not locate env $ENV_NAME" >&2; exit 1; }
BIN="$ENV_PREFIX/bin"
for tool in hmmsearch hmmalign hmmbuild blastp makeblastdb mmseqs datasets perl git; do
  [ -x "$BIN/$tool" ] || { echo "env $ENV_NAME is missing $tool" >&2; exit 1; }
done

# ---------------------------------------------------------------- 3. BITACORA at the pinned commit
if [ ! -d "$BITACORA_DIR/.git" ]; then
  say "Cloning BITACORA into $BITACORA_DIR"
  mkdir -p "$(dirname "$BITACORA_DIR")"
  "$BIN/git" clone --quiet "$BITACORA_URL" "$BITACORA_DIR"
fi
say "Pinning BITACORA to ${BITACORA_COMMIT:0:7}"
# The tracked coverage patch (predictor/_engines/bitacora_patches/) is copied over Scripts/ before
# every search, so a modified working tree there is expected; reset it before checking out.
"$BIN/git" -C "$BITACORA_DIR" checkout --quiet -- . 2>/dev/null || true
"$BIN/git" -C "$BITACORA_DIR" fetch --quiet origin || true
"$BIN/git" -C "$BITACORA_DIR" checkout --quiet "$BITACORA_COMMIT"
chmod +x "$BITACORA_DIR"/*.sh "$BITACORA_DIR"/Scripts/*.pl 2>/dev/null || true
for patch in "$REPO"/predictor/_engines/bitacora_patches/*.pl; do
  cp "$patch" "$BITACORA_DIR/Scripts/$(basename "$patch")"
done

# ---------------------------------------------------------------- 4. settings (.env)
# A source checkout reads <repo>/.env; an installed wheel reads ~/.config/tfop/.env. Write the one
# that applies, keeping every key that is not ours.
if [ -f "$REPO/pyproject.toml" ] && [ -d "$REPO/predictor" ]; then
  ENV_FILE="$REPO/.env"
else
  ENV_FILE="${XDG_CONFIG_HOME:-$HOME/.config}/tfop/.env"
fi
mkdir -p "$(dirname "$ENV_FILE")"
touch "$ENV_FILE"
set_key() {
  local key="$1" val="$2" tmp
  tmp="$(mktemp)"
  grep -v -E "^[[:space:]]*(export[[:space:]]+)?${key}[[:space:]]*=" "$ENV_FILE" > "$tmp" || true
  printf '%s="%s"\n' "$key" "$val" >> "$tmp"
  mv "$tmp" "$ENV_FILE"
}
say "Writing engine locations to $ENV_FILE"
set_key BITACORA_DIR "$BITACORA_DIR"
set_key BITACORA_BIN "$BIN"
set_key BLAST_BIN "$BIN"
set_key MMSEQS "$BIN/mmseqs"
set_key DATASETS_EXE "$BIN/datasets"

# ---------------------------------------------------------------- 5. optional Python env
# The env carries its own Python 3.12, so this needs neither sudo nor the python3-venv package.
if [ "$WANT_RIIO" = 1 ]; then
  say "Installing RIIO (genome + GUI extras) into $ENV_NAME"
  "$BIN/python" -m pip install --quiet --upgrade pip
  "$BIN/python" -m pip install --quiet -e "$REPO[genome,gui]"
fi

# ---------------------------------------------------------------- 6. readiness report
TFOP=""
if [ -x "$BIN/tfop" ]; then TFOP="$BIN/tfop"; elif command -v tfop >/dev/null; then TFOP="tfop"; fi
if [ -n "$TFOP" ]; then
  # Apply the tracked scripts, then verify discovery readiness explicitly. `tfop setup`
  # distinguishes genome-only tools from required per-TF tools, so its exit status alone
  # does not certify the production discovery route.
  if [ "$WANT_RIIO" = 1 ]; then
    "$BIN/python" - <<'PY'
from predictor.discovery import bitacora
bitacora.ensure_patched()
for check in (bitacora.engine_available, bitacora.patch_status):
    ok, why = check()
    if not ok:
        raise SystemExit(why)
PY
  fi
  say "tfop setup"
  "$TFOP" setup
else
  say "Engines installed. Now install RIIO itself (see README), then run: tfop setup"
fi
if [ "$WANT_RIIO" = 1 ]; then
  say "Done. Activate the env with:  conda activate $ENV_NAME   (micromamba: micromamba activate $ENV_NAME)"
fi
