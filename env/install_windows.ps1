<#
.SYNOPSIS
    Install RIIO natively on Windows: Python environment, search engines and settings, in one step.

.DESCRIPTION
    Run it from the RIIO folder in a normal (not administrator) PowerShell window:

        powershell -ExecutionPolicy Bypass -File env\install_windows.ps1

    What it does, in order (every step is idempotent, so re-running it repairs a partial install):
      1. finds Python 3.12 (else 3.11, 3.13, 3.10, 3.14) through the py launcher or PATH
      2. creates the virtual environment .venv in the RIIO folder and installs RIIO into it
         (editable, with the genome and GUI extras), so `tfop` reads the folder's .env
      3. downloads the official Windows builds of MMseqs2, NCBI BLAST+ (checksum-verified) and
         NCBI Datasets into env\tools\
      4. writes their locations into the settings file .env (other lines there are kept)
      5. with -WithWsl: installs BITACORA and its Linux tools inside WSL and connects them
      6. runs `tfop setup` and `tfop selftest`, then prints how to use RIIO from a new window

    It needs no administrator rights and changes nothing outside the RIIO folder except the Python
    package cache (and, with -WithWsl, a conda environment inside WSL). "-ExecutionPolicy Bypass"
    applies to this one run only.

    Genome discovery (`tfop genome`, the GUI genome scan) runs BITACORA, which needs Linux tools
    (HMMER, Perl). -WithWsl sets it up inside WSL Ubuntu; see docs\INSTALL_WINDOWS.md. Everything
    else runs natively.

.PARAMETER Python
    Full path of the Python interpreter to use instead of searching for one.

.PARAMETER SkipEngines
    Install only the Python environment (configure MMseqs2, BLAST+ and NCBI Datasets yourself).

.PARAMETER SkipSelftest
    Skip the module self-tests at the end (about 3 minutes).

.PARAMETER Force
    Download the engines again even if they are already in env\tools\.

.PARAMETER WithWsl
    Also install BITACORA for genome discovery, through WSL (needs a working WSL Linux distribution).

.PARAMETER WslDistro
    The WSL distribution to use with -WithWsl (default: Ubuntu; list yours with: wsl -l -v).

.PARAMETER WslEnvName
    The conda environment created inside WSL with -WithWsl (default: riio-engines).

.EXAMPLE
    powershell -ExecutionPolicy Bypass -File env\install_windows.ps1

.EXAMPLE
    powershell -ExecutionPolicy Bypass -File env\install_windows.ps1 -WithWsl

.EXAMPLE
    powershell -ExecutionPolicy Bypass -File env\install_windows.ps1 -Python "C:\Python312\python.exe"
#>
[CmdletBinding()]
param(
    [string]$Python = "",
    [switch]$SkipEngines,
    [switch]$SkipSelftest,
    [switch]$Force,
    [switch]$WithWsl,
    [string]$WslDistro = "Ubuntu",
    [string]$WslEnvName = "riio-engines"
)

$ErrorActionPreference = "Stop"
# Windows PowerShell 5.1 redraws a progress bar per received block, which slows downloads many-fold.
$ProgressPreference = "SilentlyContinue"
# Older .NET defaults can omit TLS 1.2, which GitHub and NCBI require.
[Net.ServicePointManager]::SecurityProtocol = [Net.ServicePointManager]::SecurityProtocol -bor [Net.SecurityProtocolType]::Tls12

$Repo = Split-Path -Parent $PSScriptRoot
$Tools = Join-Path $Repo "env\tools"
$Downloads = Join-Path $Tools "downloads"
$Venv = Join-Path $Repo ".venv"
$EnvFile = Join-Path $Repo ".env"
$MmseqsUrl = "https://github.com/soedinglab/MMseqs2/releases/latest/download/mmseqs-win64.zip"
$BlastBase = "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
$DatasetsUrl = "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/win64/datasets.exe"

function Say([string]$Message) {
    Write-Host ""
    Write-Host "==> $Message" -ForegroundColor Cyan
}

function Fail([string]$Message) {
    Write-Host ""
    Write-Host "ERROR: $Message" -ForegroundColor Red
    exit 1
}

function Invoke-Native([string]$Exe, [string[]]$Arguments) {
    # Run a program, show its output, return its exit code. Windows PowerShell 5.1 turns a program's
    # stderr lines into error records when output is redirected, which "Stop" would make fatal even
    # for a successful run, so only the exit code decides here.
    $ErrorActionPreference = "Continue"
    & $Exe @Arguments | Out-Host
    return $LASTEXITCODE
}

function Invoke-Checked([string]$What, [string]$Exe, [string[]]$Arguments) {
    $code = Invoke-Native $Exe $Arguments
    if ($code -ne 0) { Fail "$What failed (exit code $code)." }
}

function Get-NativeText([string]$Exe, [string[]]$Arguments) {
    # A program's standard output as one string, or $null when it fails.
    $ErrorActionPreference = "Continue"
    $out = & $Exe @Arguments 2>$null
    if ($LASTEXITCODE -ne 0) { return $null }
    return ($out | Out-String)
}

function ConvertTo-WslPath([string]$WindowsPath) {
    # C:\Users\me\riio-v1 -> /mnt/c/Users/me/riio-v1 (WSL's default automount), as the BITACORA
    # adapter itself maps paths. Computed here so a non-ASCII folder name is never re-encoded.
    $full = [System.IO.Path]::GetFullPath($WindowsPath)
    return "/mnt/" + $full.Substring(0, 1).ToLower() + $full.Substring(2).Replace("\", "/")
}

function Save-Url([string]$Url, [string]$OutFile) {
    Write-Host "    downloading $Url"
    try {
        Invoke-WebRequest -UseBasicParsing -Uri $Url -OutFile $OutFile
    } catch {
        Fail "could not download $Url ($($_.Exception.Message)). Check the internet connection and re-run."
    }
}

function Get-PythonVersion([string[]]$Command) {
    # "major.minor" of a candidate interpreter, or $null. The Microsoft Store "python" placeholder
    # prints an advert and fails, which counts as absent.
    $ErrorActionPreference = "Continue"
    $exe = $Command[0]
    $rest = @()
    if ($Command.Length -gt 1) { $rest = @($Command[1..($Command.Length - 1)]) }
    try {
        $out = & $exe @rest -c "import sys; print('%d.%d' % sys.version_info[:2])" 2>$null
        if ($LASTEXITCODE -eq 0 -and "$out" -match '^\d+\.\d+$') { return "$out" }
    } catch { }
    return $null
}

function Set-Setting([string]$Key, [string]$Value) {
    # Replace one KEY="value" line in .env and keep every other line. Written as UTF-8 without a BOM,
    # so a non-ASCII user folder survives and the first key stays readable.
    $lines = @()
    if (Test-Path -LiteralPath $EnvFile) { $lines = @(Get-Content -LiteralPath $EnvFile -Encoding UTF8) }
    $pattern = '^\s*(export\s+)?' + [regex]::Escape($Key) + '\s*='
    $kept = @($lines | Where-Object { $_ -notmatch $pattern })
    $kept += ('{0}="{1}"' -f $Key, $Value)
    [System.IO.File]::WriteAllLines($EnvFile, [string[]]$kept, (New-Object System.Text.UTF8Encoding($false)))
}

# ---------------------------------------------------------------- 0. where we are
if (-not (Test-Path (Join-Path $Repo "pyproject.toml")) -or -not (Test-Path (Join-Path $Repo "predictor"))) {
    Fail "run this script from the RIIO folder: it must sit in env\ next to pyproject.toml."
}
if ($Repo.Length -gt 120) {
    Write-Warning ("The RIIO folder path is $($Repo.Length) characters long. Windows limits full paths to 260 " +
                   "characters, so deep result folders may fail. A short location such as C:\riio-v1 is safer.")
}
if ($Repo -match "OneDrive") {
    Write-Warning ("The RIIO folder is inside OneDrive. Syncing can lock files in .venv and runs\; a folder " +
                   "outside OneDrive (for example C:\riio-v1) avoids that.")
}
Say "RIIO folder: $Repo"

# ---------------------------------------------------------------- 1. Python
$PythonCmd = $null
if ($Python) {
    $PythonCmd = @($Python)
    $version = Get-PythonVersion $PythonCmd
    if (-not $version) { Fail "the interpreter $Python does not run." }
} else {
    $candidates = @()
    if (Get-Command py -ErrorAction SilentlyContinue) {
        foreach ($v in @("3.12", "3.11", "3.13", "3.10", "3.14")) { $candidates += , @("py", "-$v") }
    }
    if (Get-Command python -ErrorAction SilentlyContinue) { $candidates += , @("python") }
    foreach ($candidate in $candidates) {
        $version = Get-PythonVersion $candidate
        if ($version -and @("3.10", "3.11", "3.12", "3.13", "3.14") -contains $version) {
            $PythonCmd = $candidate
            break
        }
    }
    if (-not $PythonCmd) {
        Fail ("no supported Python (3.10 to 3.14) was found. Install Python 3.12 from " +
              "https://www.python.org/downloads/windows/ (or run: winget install -e --id Python.Python.3.12), " +
              "open a NEW PowerShell window and run this script again.")
    }
}
Say "Using Python $version ($($PythonCmd -join ' '))"

# ---------------------------------------------------------------- 2. virtual environment + RIIO
$VenvPython = Join-Path $Venv "Scripts\python.exe"
if (-not (Test-Path $VenvPython)) {
    Say "Creating the virtual environment .venv"
    $rest = @()
    if ($PythonCmd.Length -gt 1) { $rest = @($PythonCmd[1..($PythonCmd.Length - 1)]) }
    Invoke-Checked "creating .venv" $PythonCmd[0] ($rest + @("-m", "venv", $Venv))
}
Say "Installing RIIO into .venv (a few minutes the first time)"
Invoke-Checked "upgrading pip" $VenvPython @("-m", "pip", "install", "--quiet", "--upgrade", "pip")
Invoke-Checked "installing RIIO" $VenvPython @("-m", "pip", "install", "--quiet", "-e", ($Repo + "[genome,gui]"))
$Tfop = Join-Path $Venv "Scripts\tfop.exe"

# ---------------------------------------------------------------- 3. engines
if (-not $SkipEngines) {
    New-Item -ItemType Directory -Force $Downloads | Out-Null

    # MMseqs2: homolog search. Bioconda has no Windows build; the official one ships BusyBox for its
    # workflow scripts. Its mmseqs.bat links the BusyBox tools as SYMBOLIC links, which need an
    # administrator prompt; HARD links (bin\bash.exe, ...) behave the same and need no rights.
    $MmseqsDir = Join-Path $Tools "mmseqs"
    $MmseqsExe = Join-Path $MmseqsDir "bin\mmseqs.exe"
    if ($Force -or -not (Test-Path $MmseqsExe)) {
        Say "Installing MMseqs2 (official Windows build) into env\tools\mmseqs"
        $zip = Join-Path $Downloads "mmseqs-win64.zip"
        Save-Url $MmseqsUrl $zip
        if (Test-Path $MmseqsDir) { Remove-Item -Recurse -Force $MmseqsDir }
        Expand-Archive -Path $zip -DestinationPath $Tools -Force
    }
    $mmBin = Join-Path $MmseqsDir "bin"
    $hasBash = @(Get-ChildItem -LiteralPath $mmBin -Force -ErrorAction SilentlyContinue |
                 Where-Object { $_.Name -eq "bash" -or $_.Name -eq "bash.exe" }).Count -gt 0
    if (-not $hasBash) {
        Invoke-Checked "installing the BusyBox helper tools" (Join-Path $mmBin "busybox.exe") @("--install", $mmBin)
    }
    Invoke-Checked "checking MMseqs2" $MmseqsExe @("version")

    # NCBI BLAST+: the newest release, verified against NCBI's published MD5 checksum.
    $blast = Get-ChildItem -LiteralPath $Tools -Directory -Filter "ncbi-blast-*" -ErrorAction SilentlyContinue |
             Sort-Object Name | Select-Object -Last 1
    if ($Force -or -not $blast -or -not (Test-Path (Join-Path $blast.FullName "bin\blastp.exe"))) {
        Say "Installing NCBI BLAST+ into env\tools"
        $listing = Join-Path $Downloads "blast-latest.html"
        Save-Url $BlastBase $listing
        $match = [regex]::Match((Get-Content -Raw -LiteralPath $listing), 'ncbi-blast-[0-9.]+\+-x64-win64\.tar\.gz(?=")')
        if (-not $match.Success) { Fail "no Windows BLAST+ archive is listed at $BlastBase" }
        $archive = Join-Path $Downloads $match.Value
        Save-Url ($BlastBase + $match.Value) $archive
        Save-Url ($BlastBase + $match.Value + ".md5") ($archive + ".md5")
        $expected = ((Get-Content -Raw -LiteralPath ($archive + ".md5")).Trim() -split '\s+')[0].ToLower()
        $actual = (Get-FileHash -Algorithm MD5 -LiteralPath $archive).Hash.ToLower()
        if ($expected -ne $actual) { Fail "BLAST+ download is corrupt (MD5 $actual, expected $expected); re-run." }
        # Windows' own tar (bsdtar). Git's GNU tar would read "C:" as a remote host.
        $tar = Join-Path $env:SystemRoot "System32\tar.exe"
        if (-not (Test-Path $tar)) { Fail "tar.exe is missing; Windows 10 (1803) or newer is required." }
        Invoke-Checked "extracting BLAST+" $tar @("-xzf", $archive, "-C", $Tools)
        $blast = Get-ChildItem -LiteralPath $Tools -Directory -Filter "ncbi-blast-*" |
                 Sort-Object Name | Select-Object -Last 1
    }
    $BlastBin = Join-Path $blast.FullName "bin"
    Invoke-Checked "checking BLAST+" (Join-Path $BlastBin "blastp.exe") @("-version")

    # NCBI Datasets: genome downloads for `tfop predict --organism` and `tfop genome <accession>`.
    $DatasetsExe = Join-Path $Tools "ncbi-datasets\datasets.exe"
    if ($Force -or -not (Test-Path $DatasetsExe)) {
        Say "Installing NCBI Datasets into env\tools\ncbi-datasets"
        New-Item -ItemType Directory -Force (Split-Path -Parent $DatasetsExe) | Out-Null
        Save-Url $DatasetsUrl $DatasetsExe
    }
    Invoke-Checked "checking NCBI Datasets" $DatasetsExe @("--version")

    # ------------------------------------------------------------ 4. settings file
    Say "Writing the engine locations to $EnvFile"
    Set-Setting "MMSEQS" $MmseqsExe
    Set-Setting "BLAST_BIN" $BlastBin
    Set-Setting "DATASETS_EXE" $DatasetsExe
    Remove-Item -Recurse -Force $Downloads -ErrorAction SilentlyContinue
}

# ---------------------------------------------------------------- 5. BITACORA through WSL (optional)
if ($WithWsl) {
    Say "Installing BITACORA and its Linux tools inside WSL ($WslDistro); 5 to 15 minutes"
    $wsl = Join-Path $env:SystemRoot "System32\wsl.exe"
    if (-not (Test-Path $wsl)) {
        Fail ("WSL is not installed. In an administrator PowerShell run  wsl --install -d Ubuntu-24.04 , restart " +
              "Windows, open Ubuntu once to create your Linux user, then run this script again with -WithWsl.")
    }
    if ((Invoke-Native $wsl @("-d", $WslDistro, "--exec", "true")) -ne 0) {
        Fail ("the WSL distribution '$WslDistro' does not start. List yours with  wsl -l -v  and pass its " +
              "name with -WslDistro, or install one with  wsl --install -d Ubuntu-24.04 .")
    }
    # The Linux installer works on this same folder: BITACORA lands in env\tools\bitacora, where the
    # Windows side finds it. Its Linux paths go to a scratch file, never into the Windows .env.
    $wslSettings = "/tmp/riio-wsl-engines-" + [guid]::NewGuid().ToString("N") + ".env"
    Invoke-Checked "the WSL engine installer" $wsl @(
        "-d", $WslDistro, "--exec", "env", "ENV_NAME=$WslEnvName", "bash",
        ((ConvertTo-WslPath $Repo) + "/env/install_linux.sh"), "--miniforge", "--settings-file", $wslSettings)
    $text = Get-NativeText $wsl @("-d", $WslDistro, "--exec", "cat", $wslSettings)
    $null = Get-NativeText $wsl @("-d", $WslDistro, "--exec", "rm", "-f", $wslSettings)
    $match = [regex]::Match("$text", '(?m)^BITACORA_BIN="([^"]+)"')
    if (-not $match.Success) { Fail "the WSL engine installer did not report where its tools are." }
    Say "Connecting BITACORA (WSL) to this installation in $EnvFile"
    Set-Setting "BITACORA_DIR" (Join-Path $Tools "bitacora")
    Set-Setting "BITACORA_WSL_DISTRO" $WslDistro
    Set-Setting "BITACORA_BIN" $match.Groups[1].Value
}

# ---------------------------------------------------------------- 6. readiness
Say "tfop setup"
$ready = ((Invoke-Native $Tfop @("setup")) -eq 0)
if (-not $SkipSelftest) {
    # A fresh environment builds matplotlib's font cache on first use; 300 s per module leaves room for
    # that on a slow computer, where the default 120 s can time out a module that is working.
    Say "tfop selftest (each module checks itself; about 3 minutes)"
    if ((Invoke-Native $Tfop @("selftest", "--timeout", "300")) -ne 0) { $ready = $false }
}

Say "Done."
Write-Host "In every NEW PowerShell window, go to the RIIO folder and activate its environment:"
Write-Host "    cd `"$Repo`""
Write-Host "    .\.venv\Scripts\Activate.ps1"
Write-Host "If PowerShell answers that running scripts is disabled, allow it for your account once:"
Write-Host "    Set-ExecutionPolicy -Scope CurrentUser RemoteSigned"
Write-Host "(or skip activation and type .venv\Scripts\tfop instead of tfop)."
Write-Host ""
Write-Host "First prediction (about 5 minutes; needs internet):"
Write-Host "    tfop predict tests\data\mini_tf.faa --name ZntR --organism GCF_000005845.2 --out runs\ZntR"
Write-Host "Graphical interface:  tfop-gui   (then open http://localhost:8501 in your browser)"
if (-not $WithWsl) {
    Write-Host "Genome discovery (tfop genome) needs BITACORA through WSL: re-run this script with -WithWsl"
    Write-Host "once WSL is installed (docs\INSTALL_WINDOWS.md)."
}
if (-not $ready) {
    Write-Host ""
    Write-Host "Some checks above did not pass; read their messages, fix, and re-run this script." -ForegroundColor Yellow
    exit 1
}
