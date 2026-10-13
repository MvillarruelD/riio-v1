# Changelog


## 1.1.0 — 2026-10-09

- Add native Linux engine installation and explicit Windows/WSL instructions.
- Share sequence, protein-list, genome-accession and local-genome execution across CLI/API/GUI.
- Preserve local genome/GFF/contig/locus context and reject incompatible input/settings on resume.
- Generate portable discovery/run HTML reports and linked supplementary XLSX/TSV tables.
- Export optional server/local AF3 job sets separately from prediction and report rebuilding.
- Retain the frozen Figure 4 recipe and pinned canonical reference data. Software/export checks
  and prediction-neutrality comparisons remain distinct validation steps.
- Fix `tfop predict` on Linux with Python 3.10-3.13 for proteins longer than 255 amino acids (also the
  README example): a raw sequence was probed as a file name, which Linux rejects with an error.
- Add `env/install_windows.ps1`, a native Windows installer for RIIO, MMseqs2, BLAST+ and NCBI Datasets
  that needs no administrator rights; `-WithWsl` adds BITACORA genome discovery through WSL.
- Rewrite the Linux and Windows installation guides for new users, with a first prediction on the
  bundled ZntR sequence, and replace the README example accession, which was not a TF.
- Give `tfop setup` installation hints that work on each OS and print the settings file it reads;
  read settings files as UTF-8 (byte-order mark ignored) so non-ASCII user folders survive on Windows.
- Serve `tfop-gui` on this computer only and without usage statistics; extra Streamlit flags pass
  through. Document every CLI option and correct the per-user settings path on Linux.
- Make `environment.yml` install this checkout instead of the 1.0.0 tag.

## 1.0.0 — 2026-10-06

First clean RIIO distribution, exported from RIO `f80c90a` and the analysis source commit recorded in
`docs/SOURCE_MANIFEST.json`. The paper's canonical prediction baseline is RIO `e9ba665`.

- Retain the supported CLI, GUI, API, active production drivers and required reference inputs.
- Package the owned BITACORA/MetalNet runners, BITACORA patch, and active folding-adapter helpers.
- Resolve the consolidated predictor locally and make canonical output roots portable.
- Include full upstream license notices, contribution evidence and provisional software citation.
- Exclude experimental code, historical runs, credentials, external engines and model weights.
- Preserve core reference database bytes and inference algorithms; release verification is documented.
