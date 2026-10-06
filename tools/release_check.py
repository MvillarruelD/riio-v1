#!/usr/bin/env python
"""release_check.py -- prove a release is installable and complete, from a clean environment.

    python tools/release_check.py
    python tools/release_check.py --keep-temp      # leave the venv and wheel for inspection
    python tools/release_check.py --skip-tests     # geometry only; for iterating on packaging

What it is for. Every check below has failed at least once in a way the ordinary test suite could
not see, because the suite runs from the CHECKOUT: the checkout has the reference data at the repo
root, an editable install puts the source tree on `sys.path`, and the working directory is the one
place every relative path happens to resolve. A wheel installed elsewhere has none of that. So this
script builds the artifact, installs it into a fresh virtual environment, changes to a directory
outside the checkout, and only then asks whether the software works.

The steps, in order:

  1. packaged-data manifest is current
  2. sources compile
  3. Ruff
  4. pytest
  5. build a wheel AND a source distribution
  6. create a clean virtual environment
  7. install the wheel non-editably
  8. change to a clean directory OUTSIDE the checkout
  9. `tfop --version`, `tfop audit`, `tfop selftest`, and a packaged-data `tfop scan`
 10. every manifest-listed file exists in the installed package
 11. caches and results are written OUTSIDE site-packages
 12. clean up (unless --keep-temp)

One environmental trap is designed around explicitly. This machine has an unrelated `struct.py`
sitting in the Windows temp root, and Python puts a script's own directory on `sys.path` -- so
running anything from that directory shadows the standard library's `struct` and produces failures
that look like they come from the code under test. Every temporary path here is therefore a fresh
SUBDIRECTORY, never the temp root itself.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
#: A short protein used only to prove the installed CLI can resolve its packaged data and run.
SMOKE_PROTEIN = (
    "MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRGRRQPIPKARRPEGRTWAQPGYPWPLYGNEG"
)


class Failure(Exception):
    """A check that did not pass. Carries the message the operator needs, and nothing else."""


def run(cmd, *, cwd=None, env=None, timeout=1800, check=True, label=None):
    """Run a command, echoing it, and return the completed process."""
    shown = " ".join(str(c) for c in cmd)
    print(f"    $ {shown[:160]}")
    result = subprocess.run([str(c) for c in cmd], cwd=str(cwd) if cwd else None,
                            env=env, capture_output=True, text=True, timeout=timeout)
    if check and result.returncode != 0:
        tail = (result.stdout or "")[-1500:] + (result.stderr or "")[-1500:]
        raise Failure(f"{label or shown} failed (exit {result.returncode})\n{tail}")
    return result


def venv_python(venv: Path) -> Path:
    """The interpreter inside a virtual environment, on either layout."""
    for rel in ("Scripts/python.exe", "bin/python"):
        candidate = venv / rel
        if candidate.exists():
            return candidate
    raise Failure(f"no interpreter inside {venv}")


def venv_script(venv: Path, name: str) -> Path:
    for rel in (f"Scripts/{name}.exe", f"bin/{name}"):
        candidate = venv / rel
        if candidate.exists():
            return candidate
    raise Failure(f"console script {name!r} was not installed into {venv}")


# --------------------------------------------------------------------------- checks
def check_data_manifest(py: Path) -> None:
    run([py, REPO / "tools" / "build_data_manifest.py", "--check"], cwd=REPO,
        label="packaged-data manifest check")


def check_compiles(py: Path) -> None:
    run([py, "-m", "compileall", "-q", "predictor", "tests", "tools"], cwd=REPO,
        label="compileall")


def check_lint(py: Path) -> None:
    run([py, "-m", "ruff", "check", "predictor", "tests", "tools"], cwd=REPO, label="ruff")


def check_tests(py: Path) -> None:
    run([py, "-m", "pytest", "-q"], cwd=REPO, label="generic pytest")
    run([py, "-m", "pytest", "-q", "-m", "released_data"], cwd=REPO,
        label="released-data pytest")


def build_artifacts(py: Path, outdir: Path) -> tuple[Path, Path]:
    """Build both a wheel and an sdist. Both must build; a release ships both."""
    outdir.mkdir(parents=True, exist_ok=True)
    try:
        run([py, "-m", "build", "--wheel", "--sdist", "--outdir", outdir], cwd=REPO,
            label="python -m build")
    except Failure as exc:
        if "No module named build" not in str(exc):
            raise
        # The documented command needs the `build` frontend from the dev extra. Fall back so a
        # missing frontend does not masquerade as a packaging defect.
        print("    ! `build` frontend absent; falling back to `pip wheel` (install .[dev] to fix)")
        run([py, "-m", "pip", "wheel", "--no-deps", "-w", outdir, REPO], cwd=REPO,
            label="pip wheel")
    wheels = sorted(outdir.glob("*.whl"))
    if not wheels:
        raise Failure(f"no wheel was produced in {outdir}")
    sdists = sorted(outdir.glob("*.tar.gz"))
    return wheels[0], (sdists[0] if sdists else None)


def make_clean_venv(py: Path, venv: Path) -> Path:
    run([py, "-m", "venv", venv], label="create virtual environment")
    vpy = venv_python(venv)
    run([vpy, "-m", "pip", "install", "--upgrade", "pip", "--quiet"], label="upgrade pip")
    return vpy


def install_wheel(vpy: Path, wheel: Path) -> None:
    # NON-editable, and from the wheel rather than the source tree: an editable install would put
    # the checkout on sys.path and hide exactly the packaging defects this script exists to find.
    run([vpy, "-m", "pip", "install", "--quiet", wheel], label="install wheel")


def installed_package_dir(vpy: Path, workdir: Path) -> Path:
    """Where the WHEEL put the package -- resolved from outside the checkout.

    `python -c` puts the current directory on `sys.path`, so running this from the repository
    imported the checkout and reported it as "the installed package". Every check downstream then
    inspected the source tree instead of site-packages and passed for the wrong reason -- which is
    the exact class of defect this whole script exists to catch, so it had to not do it itself.
    """
    env = dict(os.environ)
    env.pop("PYTHONPATH", None)
    out = run([vpy, "-c", "import predictor, pathlib; print(pathlib.Path(predictor.__file__).parent)"],
              cwd=workdir, env=env, label="locate installed package")
    package = Path(out.stdout.strip().splitlines()[-1])
    if REPO in package.parents or package == REPO / "predictor":
        raise Failure(
            f"the installed package resolved to the CHECKOUT ({package}), not to site-packages.\n"
            "    The wheel install did not take effect, or something put the source tree on the path."
        )
    return package


def check_cli(venv: Path, workdir: Path, *, timeout: int) -> None:
    """Exercise the installed console script from OUTSIDE the checkout."""
    tfop = venv_script(venv, "tfop")
    env = dict(os.environ)
    env.pop("PYTHONPATH", None)          # nothing may reach back into the checkout
    env["PREDICTOR_OUTPUT_DIR"] = str(workdir / "results")
    env["XDG_CACHE_HOME"] = str(workdir / "cache")

    run([tfop, "--version"], cwd=workdir, env=env, label="tfop --version")
    run([tfop, "audit"], cwd=workdir, env=env, label="tfop audit")
    run([tfop, "selftest", "--timeout", str(timeout)], cwd=workdir, env=env,
        timeout=timeout * 20, label="tfop selftest")

    # A scan over a tiny protein FASTA: it reads the packaged Pfam HMM and reference data, which is
    # the part a wheel most easily omits.
    fasta = workdir / "smoke.faa"
    fasta.write_text(f">smoke_tf\n{SMOKE_PROTEIN}\n", encoding="utf-8")
    run([tfop, "scan", str(fasta), "--out", str(workdir / "scan")], cwd=workdir, env=env,
        timeout=900, label="tfop scan (packaged data)")
    report = workdir / "scan" / "scan_manifest.json"
    if not report.is_file():
        raise Failure(f"tfop scan produced no manifest at {report}")


def check_manifest_files_installed(vpy: Path, package_dir: Path) -> None:
    """Every file the data manifest lists must exist in the INSTALLED package."""
    manifest = json.loads((REPO / "predictor" / "data" / "MANIFEST.json").read_text(encoding="utf-8"))
    data_dir = package_dir / "data"
    missing = [entry["path"] for entry in manifest["files"]
               if not (data_dir / entry["path"]).is_file()]
    if missing:
        preview = "\n      ".join(missing[:15])
        raise Failure(f"{len(missing)} manifest file(s) absent after installation:\n      {preview}"
                      + (f"\n      ... and {len(missing) - 15} more" if len(missing) > 15 else "")
                      + "\n    Add the pattern to [tool.setuptools.package-data] in pyproject.toml.")
    print(f"    all {len(manifest['files'])} manifest files present in the installed package")
    required = [
        "_engines/bitacora_runner.py", "_engines/metalnet_runner.py",
        "_engines/bitacora_patches/get_blastp_parsed_newv2.pl",
        "_engines/esmfold_lib/_lib/esmfold2_backend.py",
        "discovery/vendor/LICENSE", "discovery/vendor/CITATION.cff",
    ]
    absent = [name for name in required if not (package_dir / name).is_file()]
    if absent:
        raise Failure(f"engine adapters or upstream notices absent from installed wheel: {absent}")
    print("    engine adapters, patch, and upstream notices present")


def check_nothing_written_into_site_packages(package_dir: Path, before: dict) -> None:
    """An installed package is read-only. Caches and results belong outside it."""
    after = {p: p.stat().st_mtime for p in package_dir.rglob("*") if p.is_file()}
    added = sorted(set(after) - set(before))
    # `.pyc` files are written by the interpreter itself on import; they are not the runtime
    # writing data, and excluding them keeps this check about what the code does.
    added = [p for p in added if p.suffix != ".pyc" and "__pycache__" not in p.parts]
    if added:
        preview = "\n      ".join(str(p.relative_to(package_dir)) for p in added[:10])
        raise Failure(f"runtime wrote {len(added)} file(s) inside the installed package:\n"
                      f"      {preview}\n"
                      "    Caches belong under the platform user cache, results under "
                      "PREDICTOR_OUTPUT_DIR.")
    print("    nothing was written inside site-packages")


def check_outputs_landed_outside(workdir: Path) -> None:
    produced = [p for p in (workdir / "results", workdir / "cache", workdir / "scan")
                if p.exists()]
    if not produced:
        raise Failure("no results or cache directory was created outside the package; the run may "
                      "have written somewhere unexpected")
    print(f"    outputs landed outside the package: "
          f"{', '.join(p.name for p in produced)}")


# --------------------------------------------------------------------------- driver
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--keep-temp", action="store_true",
                    help="do not delete the temporary venv, wheel and working directory")
    ap.add_argument("--skip-tests", action="store_true",
                    help="skip compileall/ruff/pytest; check packaging geometry only")
    ap.add_argument("--selftest-timeout", type=int, default=120,
                    help="per-module timeout passed to `tfop selftest` (default 120s)")
    a = ap.parse_args(argv)

    py = Path(sys.executable)
    # A fresh SUBDIRECTORY, never the temp root: the root can contain modules that shadow the
    # standard library for anything run from it. See the module docstring.
    root = Path(tempfile.mkdtemp(prefix="tfop_release_"))
    venv = root / "venv"
    dist = root / "dist"
    workdir = root / "work"
    workdir.mkdir(parents=True)

    started = time.time()
    steps: list[tuple[str, str]] = []
    print(f"release check for {REPO}")
    print(f"  temporary root: {root}\n")

    def step(name, fn, *args, **kwargs):
        print(f"[{name}]")
        t0 = time.time()
        fn(*args, **kwargs)
        steps.append((name, f"ok ({time.time() - t0:.0f}s)"))
        print()

    try:
        step("data manifest", check_data_manifest, py)
        if not a.skip_tests:
            step("compile", check_compiles, py)
            step("lint", check_lint, py)
            step("tests", check_tests, py)
        else:
            steps.append(("compile/lint/tests", "SKIPPED (--skip-tests)"))

        print("[build]")
        wheel, sdist = build_artifacts(py, dist)
        size_mb = wheel.stat().st_size / 1e6
        print(f"    wheel  {wheel.name}  ({size_mb:.1f} MB)")
        print(f"    sdist  {sdist.name if sdist else 'NOT BUILT'}")
        if sdist is None:
            raise Failure("no source distribution was produced; a release ships both")
        steps.append(("build", f"ok (wheel {size_mb:.1f} MB + sdist)"))
        print()

        print("[clean venv]")
        t0 = time.time()
        vpy = make_clean_venv(py, venv)
        steps.append(("clean venv", f"ok ({time.time() - t0:.0f}s)"))
        print()

        step("install wheel", install_wheel, vpy, wheel)

        package_dir = installed_package_dir(vpy, workdir)
        print(f"[installed] {package_dir}\n")
        before = {p: p.stat().st_mtime for p in package_dir.rglob("*") if p.is_file()}

        step("manifest files installed", check_manifest_files_installed, vpy, package_dir)
        step("installed CLI (outside the checkout)", check_cli, venv, workdir,
             timeout=a.selftest_timeout)
        step("outputs outside the package", check_outputs_landed_outside, workdir)
        step("package stayed read-only", check_nothing_written_into_site_packages,
             package_dir, before)

    except Failure as exc:
        print(f"\nRELEASE CHECK FAILED\n\n{exc}\n")
        if a.keep_temp:
            print(f"temporary tree kept at {root}")
        else:
            shutil.rmtree(root, ignore_errors=True)
        return 1
    except subprocess.TimeoutExpired as exc:
        print(f"\nRELEASE CHECK FAILED: timed out after {exc.timeout}s running {exc.cmd}\n")
        if not a.keep_temp:
            shutil.rmtree(root, ignore_errors=True)
        return 1

    print("=" * 72)
    for name, outcome in steps:
        print(f"  {name:<40} {outcome}")
    print("=" * 72)
    print(f"RELEASE CHECK PASSED in {(time.time() - started) / 60:.1f} min")
    if a.keep_temp:
        print(f"temporary tree kept at {root}")
    else:
        shutil.rmtree(root, ignore_errors=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
