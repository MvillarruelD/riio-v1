# Installing and using RIIO on Windows

This guide assumes no previous experience with WSL, conda or Python packaging. RIIO runs on Windows 10
(version 2004 or newer) and Windows 11, by either of two routes:

| | Route A: inside WSL (recommended) | Route B: native Windows |
|---|---|---|
| What works | Everything, exactly as on Linux | Everything; genome discovery runs its Linux part in WSL |
| You install | WSL (needs administrator rights once) | Python 3.12 and Git (no administrator rights) |
| One-step installer | `bash env/install_linux.sh --miniforge --with-riio` | `powershell -ExecutionPolicy Bypass -File env\install_windows.ps1` |
| Disk space | about 3 GB | about 2 GB (another 2.5 GB in WSL for genome discovery) |

WSL (Windows Subsystem for Linux) is Microsoft's built-in Linux environment. Route A uses the same
Linux installer that the project's end-to-end integration test exercises, and is the simplest to keep
working. Choose Route B if WSL is not allowed on your computer, or if you only need protein-level
predictions.

## Route A: RIIO inside WSL

### A1. Install WSL (once per computer)

Open PowerShell **as administrator** (Start menu, type `PowerShell`, right-click it, choose
*Run as administrator*) and run:

```powershell
wsl --install
```

Restart Windows when asked. Then open **Ubuntu** from the Start menu. The first time, it asks you to
choose a Linux user name and password. The password is not shown while you type; that is normal.
Remember it: Ubuntu asks for it when installing system packages.

If installation fails with an error mentioning virtualization (for example `0x80370102`), enable
virtualization (Intel VT-x or AMD-V/SVM) in your computer's UEFI/BIOS settings and try again.

### A2. Install and use RIIO

In the Ubuntu window, follow [the Linux guide](INSTALL_LINUX.md) from its step 1. Every command there
works unchanged. Keep RIIO in your Ubuntu home folder (`~/riio-v1`), not under `/mnt/c`.

### A3. Working between Windows and Ubuntu

- **Open a report in your Windows browser:**
  `explorer.exe "$(wslpath -w runs/ZntR/jobs/ZntR/REPORT.html)"`
- **Browse your RIIO folder in Windows Explorer:** run `explorer.exe .` in the RIIO folder, or type
  `\\wsl.localhost\Ubuntu\home\<your Linux user name>\riio-v1` in the Explorer address bar.
- **Use a file from Windows:** your `C:` drive appears in Ubuntu as `/mnt/c`, for example
  `cp /mnt/c/Users/<you>/Downloads/MyTF.faa ~/riio-v1/`.
- **Graphical interface:** run `tfop-gui` in Ubuntu, then open <http://localhost:8501> in your Windows
  browser.

## Route B: RIIO directly on Windows

### B1. Install Python and Git (once per computer)

- **Python 3.12:** download the *Windows installer (64-bit)* from
  <https://www.python.org/downloads/windows/> and run it with the default options. It also installs
  the `py` launcher, which the RIIO installer uses, so ticking *Add python.exe to PATH* is optional.
- **Git:** download it from <https://git-scm.com/download/win> and install it with the default options.

Alternatively, in PowerShell: `winget install -e --id Python.Python.3.12` and
`winget install -e --id Git.Git`. Close PowerShell afterwards and open a **new** window, so it finds
the new programs.

### B2. Download RIIO

Open PowerShell (not as administrator) and run:

```powershell
cd $HOME
git clone https://github.com/MvillarruelD/riio-v1.git
cd riio-v1
```

Use a short folder outside OneDrive, such as `C:\Users\<you>\riio-v1` or `C:\riio-v1`. Windows limits
full paths to 260 characters, and OneDrive synchronization can lock files while RIIO writes them.

### B3. Run the installer

```powershell
powershell -ExecutionPolicy Bypass -File env\install_windows.ps1
```

It takes about 10 minutes and needs no administrator rights. `-ExecutionPolicy Bypass` allows this one
script to run without changing any Windows setting. It:

- creates a Python environment in the `.venv` folder and installs RIIO into it;
- downloads the official Windows builds of MMseqs2, NCBI BLAST+ (checksum-verified) and NCBI Datasets
  into `env\tools` (about 200 MB);
- writes their locations into the settings file `.env` in the RIIO folder;
- runs `tfop setup` and `tfop selftest` (about 3 minutes).

Every engine row of `tfop setup` should say `OK`, except BITACORA, which stays `MISSING` until you add
genome discovery (B6), and the optional MetalNet2 (see [ENGINES](ENGINES.md)). If the script stops with
an error, read the message, fix the cause and run it again: every step can be repeated safely.

### B4. Activate RIIO in each new PowerShell window

```powershell
cd $HOME\riio-v1
.\.venv\Scripts\Activate.ps1
```

The prompt now starts with `(.venv)`. If PowerShell answers that *running scripts is disabled on this
system*, allow scripts that you create or that come from trusted sources, for your Windows account only,
then try again:

```powershell
Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
```

If you prefer not to change that setting, skip activation and type `.venv\Scripts\tfop` wherever this
guide says `tfop`.

### B5. Make a first prediction

RIIO includes the *E. coli* ZntR protein as a test input:

```powershell
tfop predict tests\data\mini_tf.faa --name ZntR --organism GCF_000005845.2 --out runs\ZntR
start runs\ZntR\jobs\ZntR\REPORT.html
```

The prediction takes about 3 to 5 minutes and needs internet access; `start` opens the report in your
browser. [Step 7 of the Linux guide](INSTALL_LINUX.md#7-analyse-your-own-proteins-and-genomes) explains
how to run your own proteins and genomes; the commands are identical (Windows also accepts `/` in paths).
For the graphical interface, run `tfop-gui` and open <http://localhost:8501>.

### B6. Add genome discovery (optional)

`tfop genome` and the GUI genome scan use BITACORA, which needs Linux tools. Install WSL as in step A1,
open Ubuntu once to create your Linux user, then run the installer again with `-WithWsl`:

```powershell
powershell -ExecutionPolicy Bypass -File env\install_windows.ps1 -WithWsl
```

It installs Miniforge, the Linux search tools and the pinned BITACORA inside WSL (5 to 15 minutes), and
records them in the RIIO folder's `.env`. `tfop setup` then shows BITACORA as `OK`, and
`tfop genome GCF_000005845.2 --out runs\ecoli --discover-only` lists the candidate regulators of
*E. coli* K-12. If your WSL distribution is not called `Ubuntu` (see `wsl -l -v`), add
`-WslDistro <name>`. Discovery is slower than in Route A because its files cross between Windows and WSL.

## Optional: an NCBI API key

RIIO works without one, but NCBI limits anonymous requests. A free key, created in your NCBI account
settings, makes large runs faster and more reliable. Open the settings file in Notepad (Route A: use
`nano .env` in Ubuntu instead):

```powershell
notepad .env
```

Add the lines `NCBI_TOKEN=your-key` and `NCBI_EMAIL=you@example.org` and save. `tfop setup` then shows
`NCBI_TOKEN=set`. The file `.env.example` describes every setting.

## Updating and uninstalling

To update, run `git pull` in the RIIO folder, then the installer of your route again.

To uninstall Route B, delete the RIIO folder and the cache folder `%LOCALAPPDATA%\tfop`. If you used
`-WithWsl`, also remove its environment inside WSL:
`wsl -d Ubuntu -- ~/miniforge3/bin/conda env remove -y -n riio-engines`. Route A is uninstalled as on
Linux; WSL itself can be removed in *Settings > Apps*.

## Troubleshooting

| Message or symptom | What to do |
|---|---|
| *running scripts is disabled on this system* | See step B4, or use `.venv\Scripts\tfop`. |
| `python` opens the Microsoft Store, or *Python was not found* | Install Python from python.org (B1); the RIIO installer uses the `py` launcher. |
| *no supported Python (3.10 to 3.14) was found* | Install Python 3.12 (B1), then open a new PowerShell window. |
| `git` or `tfop` *is not recognized* | Open a new window after installing Git; activate RIIO (B4) before using `tfop`. |
| A download fails | Check the connection or proxy and run the installer again. `tfop setup` lists manual download links. |
| Antivirus software quarantines `busybox.exe` or `mmseqs.exe` | They come from the official MMseqs2 release; restore them, then run the installer again. |
| Errors about long paths, or files that cannot be found | Move the RIIO folder to a short path such as `C:\riio-v1` and run the installer again. |
| `WslRegisterDistribution failed` or `0x80370102` | Enable virtualization in UEFI/BIOS (A1). |
| The GUI page does not open | Use <http://localhost:8501>, and keep the window running `tfop-gui` open. |
| BITACORA `MISSING` after `-WithWsl` | Run `tfop setup`; its *why:* line says what WSL lacks. Run the installer with `-WithWsl` again. |
| Anything else | Run `tfop setup` and `tfop audit`, and read the per-protein logs in `runs\<name>\logs\`. |

## Manual configuration (advanced)

Each installer records the engine locations in `.env` in the RIIO folder; `tfop setup` prints the
settings file it read and [ENGINES](ENGINES.md) lists every setting. Install RIIO from its folder in
editable mode (`.venv\Scripts\python -m pip install -e ".[genome,gui]"`), as the installer does, for it
to read that file. A package installed any other way, such as the wheel from the Releases page, reads
`%LOCALAPPDATA%\tfop\tf-operator-predictor\.env` instead, or the file named by `PREDICTOR_ENV_FILE`.

For a hand-made genome-discovery bridge, keep BITACORA in a folder that both systems can reach (for
example `env\tools\bitacora` in the RIIO folder) and set three keys in the Windows `.env`:
`BITACORA_DIR` (its Windows path), `BITACORA_WSL_DISTRO` (the distribution name) and `BITACORA_BIN`
(the Linux folder holding `hmmsearch`, `blastp` and `perl` inside WSL). `\\wsl$` paths are not supported
for `BITACORA_DIR`. The tracked BITACORA patches are applied before each discovery; repair an existing
checkout with:

```powershell
python -c "from predictor.discovery import bitacora; print(bitacora.ensure_patched())"
tfop setup
```

Logs and output folders are ordinary files; keep each run folder whole when sharing its reports.
