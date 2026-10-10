# External engines and settings

The Python wheel includes the reference data and pyhmmer classifier. External executables remain separate.

| Engine | Purpose | Resolution / setting | Install |
|---|---|---|---|
| MMseqs2 | Homolog/SSN collection | `MMSEQS`, vendored engine, PATH | Linux bioconda; Windows official archive |
| BLAST+ | Locus/genome and BITACORA protein searches | `BLAST_BIN`, PATH | bioconda / official NCBI binary |
| NCBI Datasets | Genome accession acquisition | `DATASETS_EXE`, PATH | conda-forge `ncbi-datasets-cli` |
| pyhmmer | Packaged Pfam family classifier | Python dependency | pip dependency |
| HMMER | BITACORA HMM search/clade models | native PATH or `BITACORA_BIN`; WSL on Windows | bioconda `hmmer` |
| Perl / patched BITACORA | Production candidate discovery | `BITACORA_DIR`, `BITACORA_BIN`, `BITACORA_WSL_DISTRO` | pinned installer checkout |
| MetalNet2 | Optional metal-site evidence | `METALNET_PYTHON`, `METALNET2_DIR` | separate compatible model environment |
| ESMFold adapter | Optional folding | `PREDICTOR_ESMFOLD_LIB`, `BIOHUB_TOKEN` | separate adapter/token |

Native Linux uses native subprocesses, with `BITACORA_BIN` prepended to PATH.
Only Windows invokes WSL. The installer pins BITACORA to
`715950426b70b029279343dcc541115dfee20352` and copies `predictor/_engines/bitacora_patches/*.pl`.
`tfop setup` checks those files against the installed Scripts directory.

```bash
python -c "from predictor.discovery import bitacora; print(bitacora.ensure_patched())"
tfop setup
```

The repair command deliberately updates installed engine scripts; `tfop setup` itself only inspects.
Clade-model builds use private scratch directories. When a clade model cannot be built, warnings and
`[pfam_fallback]` profile tags expose the family-model fallback; investigate before comparing candidates.

Settings files accept unquoted or quoted values, spaces in quotes, and inline comments outside quotes.
Read order: explicit `PREDICTOR_ENV_FILE`, checkout `.env`, user configuration file; existing environment
variables win. Supported keys include:

- `MMSEQS`, `BLAST_BIN`, `DATASETS_EXE`, `BITACORA_DIR`, `BITACORA_BIN`, `BITACORA_PY`, `BITACORA_WSL_DISTRO`;
- `METALNET_PYTHON`, `METALNET2_DIR`;
- `BIOHUB_TOKEN`, `NCBI_TOKEN`, `NCBI_API_KEY`, `NCBI_EMAIL`.

Path and output overrides include `PREDICTOR_CACHE_DIR`, `PREDICTOR_OUTPUT_DIR`,
`PREDICTOR_GENOME_MIRROR` and `PREDICTOR_ENV_FILE`; inspect their resolution with `tfop audit`.
Tokens are optional for public NCBI requests but may affect service rate limits; never share `.env`.
MetalNet2 requires its model-compatible Python environment and weights; the maintainer record in
`env/SETUP.md` describes that legacy installation. Folding and MetalNet2 absence remain explicit abstentions.

Unpublished neighbouring `5.8 Promoters` corpora are disabled by default. Only maintainers rebuilding
source references should set `PREDICTOR_ALLOW_LEGACY_WORKSPACE=1`; ordinary runs use packaged reference bytes.
This flag does not supply a missing deployed database. See the per-OS installation guides for native paths.
