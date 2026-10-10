# Windows installation

Two supported interfaces are available: a native Windows Python application with BITACORA through WSL,
or the complete native-Linux application inside WSL Ubuntu. WSL-only installation follows
[INSTALL_LINUX](INSTALL_LINUX.md); keep that checkout under the Linux home directory.

For native Python, clone a checkout and create an environment:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install ".[genome,gui]"
tfop audit
tfop setup
```

Install BLAST+ and NCBI Datasets from their official distributions or conda-forge. Download the
MMseqs2 Windows distribution from its [official releases](https://github.com/soedinglab/MMseqs2/releases)
and keep its complete directory (including workflow shell utilities). Set `MMSEQS` to its executable
or `.bat` launcher, `BLAST_BIN` to the BLAST binary directory, and `DATASETS_EXE` to `datasets.exe`.
`tfop setup` reports the exact paths selected. Do not assume a conda Windows executable is under `bin`;
it may be under `Library/bin` or another directory shown by the installation.

For BITACORA, install WSL Ubuntu with Perl, HMMER and BLAST in its conda environment. The native Windows
adapter requires the BITACORA checkout in a Windows-accessible path, for example
`C:/riio/env/tools/bitacora`; WSL sees that directory as `/mnt/c/riio/env/tools/bitacora`.
Run the Linux installer from the corresponding checkout with `BITACORA_DIR` set to that WSL path.
Then set the native application's `.env` locations using Windows paths for native components:

```dotenv
BITACORA_DIR="C:/riio/env/tools/bitacora"
BITACORA_WSL_DISTRO=Ubuntu
BITACORA_BIN=$HOME/miniforge3/envs/riio-engines/bin
BLAST_BIN="C:/Tools/ncbi-blast/bin"
MMSEQS="C:/Tools/mmseqs/mmseqs.bat"
DATASETS_EXE="C:/Tools/datasets.exe"
```

Replace examples with actual locations. `BITACORA_BIN` is the Linux bin directory inside WSL,
while `BLAST_BIN`, `MMSEQS`, `DATASETS_EXE` and `BITACORA_DIR` are native Windows paths.
The installer writes Linux settings; use a separate `PREDICTOR_ENV_FILE` for native Windows if both
applications share a checkout. A `\\wsl$` engine path is not the drive-to-WSL bridge supported here.
The simplest way to avoid the mixed path arrangement is to run the whole application in WSL.

Tracked patches are applied before discovery. For an already-installed checkout, repair them with:

```powershell
python -c "from predictor.discovery import bitacora; print(bitacora.ensure_patched())"
tfop setup
tfop selftest
tfop-gui
```

If activation is restricted by PowerShell policy, invoke `.venv/Scripts/python.exe` and
`.venv/Scripts/tfop.exe` directly instead of changing system policy. Logs and output directories
remain ordinary files; preserve directory structure when sharing reports.
