# RIIO v1.1 release preparation audit

These checks cover the local source/export snapshots listed below. A passing software
check does not establish experimental biological accuracy or historical prediction neutrality.

## Source provenance

- RIO: `c20b2a318434ccbcaa9f7181a00c6142d1dad55b`
- jjr-ssn-analysis: `18d0baedcf15964559c7f5afb8cf7c10a7aee0d5`
- publication_checkout: `2810de652ab70d150e63c283e7cd3b77b5eb3d8e`

Complete file checksums are in [SOURCE_MANIFEST.json](SOURCE_MANIFEST.json).
The [v1.0 release audit](RELEASE_AUDIT.md) and `EXPORT_AUDIT.json` are historical records.

## Completed checks

- Predictor source suites: Ruff, generic and released-data pytest passed on Windows Python 3.14.3 for the code of predictor 12cde97; 138 focused CLI, GUI, portability and packaging tests passed at c20b2a3.
- Analysis release-tooling tests: analysis/tests/test_v11_implementation.py: 18 passed, including the source-archive contract that now requires env/install_windows.ps1.
- Final staging generic tests (Windows Python 3.12.13, README developer install): 869 passed, 8 skipped, 51 deselected.
- Final staging released-data tests: 48 passed, 3 skipped (two need MMseqs2 on PATH, one needs a Git work tree).
- Ruff and packaged-data manifest: ruff check predictor tests tools analysis reports and tools/build_data_manifest.py --check both passed.
- CI install mimic: with only .[dev,gui] installed, as the v1.1 workflow does: 867 passed, 10 skipped on the preceding staging, which differs only in docs/INSTALL_LINUX.md, docs/INSTALL_WINDOWS.md, env/install_windows.ps1 and the selftest timeout default in predictor/run_cli.py.
- Isolated release check: tools/release_check.py --skip-tests --selftest-timeout 120 passed in 7.0 min: wheel (36.6 MB) and sdist, clean non-editable install, 479 manifest files, installed CLI selftest 45 ok / 3 engine skips / 0 timed out, outputs outside the package, package left unmodified.
- Source archive contract: the 1,043-member sdist includes both installers, the Linux engine environment, Streamlit configuration, .env.example, environment.yml and the smoke fixtures.
- Fresh Linux installation (WSL Ubuntu 26.04): env/install_linux.sh --miniforge --with-riio; tfop setup with every engine including BITACORA; tfop selftest 48 ok; bundled ZntR (202 s) and ModE NP_415282.1 (262 aa, which crashed before this fix) predictions completed; the former README accession no longer crashes; tfop-gui listened on 127.0.0.1 only; the environment.yml route was ready; run on the preceding staging, whose Linux installer is byte-identical.
- Fresh Windows installation (Python 3.12.13, minimal PATH): env/install_windows.ps1 completed in 9 min (venv, MMseqs2 with hard-linked BusyBox, MD5-verified BLAST+ 2.17.0, NCBI Datasets, .env); tfop setup ready; bundled ZntR and ModE predictions completed; tfop-gui answered on 127.0.0.1 only. Its 120 s selftest timed out signals.operator_logo under concurrent load (79 s when run alone); the default is now 300 s.
- Windows genome discovery through WSL: install_windows.ps1 -WithWsl installed the WSL engine environment and the pinned BITACORA in env/tools/bitacora and recorded BITACORA_DIR, BITACORA_WSL_DISTRO and BITACORA_BIN; tfop setup showed BITACORA OK, and native-Windows discovery on the bundled K-12 fragment returned zntR (NP_417751.1) from profile ecZntR|PbrR.

## Outstanding gates

- Exact historical 12-TF prediction-neutrality comparison.
- Remote Linux/Windows CI after an approved push.
- Manual AF3 server upload and acceptance.
- PDF conversion of the new review DOCX (converter unavailable).

No push, tag, release or main-branch merge is represented by this audit.
