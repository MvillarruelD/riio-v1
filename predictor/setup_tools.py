"""
setup_tools.py -- `tfop setup`: verify the external engines the pipeline shells out to and print a
readiness table + install hints.

The orchestrator is light Python; the heavy work is done by external tools. This command does NOT
auto-install them (they have their own installers) -- it checks the SAME resolution logic the pipeline
itself uses (not a naive PATH lookup: MMseqs2 and HMMER have vendored/pip paths that a plain PATH check
would miss), reports what's ready, and tells you how to get what's missing. Engines the code can't find
degrade gracefully: the source that needs one abstains and says so, rather than the run failing.
"""
from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path


def _which_check(*names, env_var: str | None = None):
    """Is one of `names` runnable? `env_var` names an override the CONSUMER already honours.

    Passing it matters: the audit and the code that actually runs the tool must agree. NCBI Datasets
    is located by `genome_mirror` via `$DATASETS_EXE` first, so a check that only consulted PATH
    reported MISSING over a working install -- a false negative from the one report whose job is to
    tell you what is missing.
    """
    def _check():
        if env_var:
            override = os.environ.get(env_var)
            if override and Path(override).exists():
                return True, override
        for n in names:
            p = shutil.which(n)
            if p:
                return True, p
        return False, None
    return _check


def _mmseqs_check():
    # Reuse the pipeline's OWN resolver (env var -> vendored env/tools/mmseqs -> PATH), not a naive PATH
    # check: bioconda ships no win-64 build at all, so on Windows this is normally a vendored binary, not PATH.
    try:
        from predictor.annotate.homologs import find_mmseqs
        return True, find_mmseqs()
    except Exception:
        return False, None


def _pyhmmer_check():
    # The code classifies TF families via the `pyhmmer` PYTHON PACKAGE (HMMER bundled in the wheel), not a
    # standalone `hmmscan` binary on PATH -- pyhmmer ships wheels for Windows/macOS/Linux alike.
    try:
        import pyhmmer
        return True, getattr(pyhmmer, "__version__", "installed")
    except Exception:
        return False, None


def _metalnet_check():
    # MetalNet2's classifiers are AutoGluon 0.8.0 pickles that need Python 3.9 / numpy 1.25 / pandas
    # 1.5.3 / sklearn 1.2.2, so they cannot be loaded by this interpreter and the engine lives in its own
    # conda env. Ask the adapter, which resolves the env, the vendored checkout AND the model weights --
    # a PATH check would mean nothing here.
    try:
        from predictor.structure import metalnet as _mn
        st = _mn.status()
        return bool(st.get("ready")), (st.get("python") if st.get("ready") else None) or st.get("detail")
    except Exception:
        return False, None


def _bitacora_check():
    # Ask the adapter: it resolves $BITACORA_DIR, $BITACORA_BIN and (on Windows) WSL exactly as the
    # genome route does, so this row cannot disagree with `tfop genome`.
    try:
        from predictor.discovery import bitacora as _b
        ok, why = _b.engine_available()
        if ok:
            ok, why = _b.patch_status()
        return ok, (str(_b.ENGINE_DIR) if ok else why)
    except Exception as exc:
        return False, f"{type(exc).__name__}: {exc}"


def _blast_check(tool: str):
    def _check():
        base = os.environ.get("BLAST_BIN")
        if base:
            for name in (tool, tool + ".exe"):
                cand = Path(base) / name
                if cand.exists():
                    return True, str(cand)
        p = shutil.which(tool)
        return (True, p) if p else (False, None)
    return _check


def _datasets_hint() -> str:
    exe = "datasets.exe" if sys.platform.startswith("win") else "datasets"
    return ("conda install -n ncbi-datasets -c conda-forge ncbi-datasets-cli, then set $DATASETS_EXE to "
            f"<env>/bin/{exe} (env/install_linux.sh does this on Linux) -- a separate env keeps the "
            "solver away from the working one")


def _bitacora_hint() -> str:
    if sys.platform.startswith("win"):
        return ("BITACORA needs HMMER, which has no Windows build: install it inside WSL Ubuntu -- see "
                "docs/INSTALL_WINDOWS.md")
    return "bash env/install_linux.sh  (or see docs/INSTALL_LINUX.md, section BITACORA)"


def _mmseqs_hint() -> str:
    if sys.platform.startswith("win"):
        return ("no bioconda win-64 build exists; download mmseqs-win64.zip from "
                "https://github.com/soedinglab/MMseqs2/releases and extract to env/tools/mmseqs/ "
                "(see docs/ENGINES.md), or set $MMSEQS to the exe path")
    return "conda install -c bioconda mmseqs2"


# (display name, checker() -> (present, detail), required?, purpose, install hint)
_TOOLS = [
    ("MMseqs2", _mmseqs_check, True, "homolog search / SSN clustering", _mmseqs_hint()),
    ("BLAST+ blastp", _blast_check("blastp"), True, "TF-locus / genome resolution",
     "conda install -c bioconda blast"),
    ("BLAST+ makeblastdb", _blast_check("makeblastdb"), True, "per-genome BLAST DB",
     "conda install -c bioconda blast"),
    ("pyhmmer (Pfam/HMMER)", _pyhmmer_check, True, "Pfam family classification", "pip install pyhmmer"),
    ("NCBI datasets", _which_check("datasets", env_var="DATASETS_EXE"), True, "genome acquisition",
     _datasets_hint()),
    ("BITACORA", _bitacora_check, "genome", "candidate discovery (`tfop genome`, GUI Genome tab)",
     _bitacora_hint()),
    ("MetalNet2", _metalnet_check, False, "family-agnostic metal-site gate (inducer source 5)",
     "needs the `metalnet` conda env + the two AutoGluon model bundles -- see docs/ENGINES.md "
     "(section MetalNet2); without it the source abstains and every other verdict is unchanged"),
]


def check() -> dict:
    """Resolve every external engine + the tokens. Pure inspection; safe to run anywhere/offline."""
    from predictor import resources
    # Load `.env` FIRST: it carries the engine locations as well as the tokens, and a check that ran
    # before it would report a tool missing that the pipeline itself would find.
    resources.load_tokens()
    rows = []
    for name, checker, required, purpose, hint in _TOOLS:
        present, detail = checker()
        rows.append({"name": name, "required": required, "present": present, "path": detail,
                     "purpose": purpose, "hint": hint})
    return {"rows": rows,
            # DeepPBS and FoldX were listed here until 2026-09-02. Nothing in the package calls them
            # any more -- their only consumer was the Stage-B `tfop finalize`, removed with them -- and
            # reporting an engine as "optional" implies some code path would use it if present.
            "optional": {"ESMFold2 adapter": resources.esmfold_lib_dir() is not None},
            "tokens": resources.load_tokens()}


def main(argv=None) -> int:
    r = check()
    print("=" * 78)
    print("tfop setup -- external engine readiness")
    print("=" * 78)
    print("Core reference database: bundled with the Python package (no separate download).")
    print("This check covers executable tools used to fetch genomes and run searches.\n")
    print(f"{'engine':<22} {'req':<4} {'status':<9} purpose")
    print("-" * 78)
    missing, genome_missing = [], []
    for row in r["rows"]:
        status = "OK" if row["present"] else "MISSING"
        req = {True: "yes", False: "opt"}.get(row["required"], "gen")
        print(f"{row['name']:<22} {req:<4} {status:<9} {row['purpose']}")
        if row["present"] and row.get("path"):
            print(f"{'':<22} {'':<4} {'':<9} -> {row['path']}")
        if row["required"] is True and not row["present"]:
            missing.append(row)
        elif row["required"] == "genome" and not row["present"]:
            genome_missing.append(row)
    print("-" * 78)
    print("req: yes = every run; gen = whole-genome discovery only; opt = optional evidence source")
    print("optional:  " + ", ".join(f"{k}={'OK' if v else '-'}" for k, v in r["optional"].items()))
    print("tokens:    " + ", ".join(f"{k}={'set' if v else '-'}" for k, v in r["tokens"].items())
          + "   (set them in .env; see .env.example)")
    for row in genome_missing:
        print(f"\nNOTE: {row['name']} is not ready, so `tfop genome` and the GUI Genome tab cannot "
              f"discover candidates (single-protein and protein-list runs are unaffected).\n"
              f"  why: {row['path']}\n  fix: {row['hint']}")
    if missing:
        print("\nMissing REQUIRED engines -- install, then re-run `tfop setup`:")
        seen = set()
        for row in missing:
            if row["hint"] not in seen:
                print(f"  - {row['name']:<22} {row['hint']}")
                seen.add(row["hint"])
        return 1
    print("\nREADY for standard predictions. Run `tfop audit` to see the bundled database, cache, "
          "output, and credential locations.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
