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


def _mmseqs_hint() -> str:
    if sys.platform.startswith("win"):
        return ("no bioconda win-64 build exists; download mmseqs-win64.zip from "
                "https://github.com/soedinglab/MMseqs2/releases and extract to env/tools/mmseqs/ "
                "(see docs/ENGINES.md), or set $MMSEQS to the exe path")
    return "conda install -c bioconda mmseqs2"


# (display name, checker() -> (present, detail), required?, purpose, install hint)
_TOOLS = [
    ("MMseqs2", _mmseqs_check, True, "homolog search / SSN clustering", _mmseqs_hint()),
    ("BLAST+ blastp", _which_check("blastp"), True, "TF-locus / genome resolution",
     "conda install -c bioconda blast"),
    ("BLAST+ makeblastdb", _which_check("makeblastdb"), True, "per-genome BLAST DB",
     "conda install -c bioconda blast"),
    ("pyhmmer (Pfam/HMMER)", _pyhmmer_check, True, "Pfam family classification", "pip install pyhmmer"),
    ("NCBI datasets", _which_check("datasets", env_var="DATASETS_EXE"), True, "genome acquisition",
     "conda install -n ncbi-datasets -c conda-forge ncbi-datasets-cli, then set $DATASETS_EXE to "
     "<env>/bin/datasets.exe -- a separate env keeps the solver away from the working `data` one"),
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
    missing = []
    for row in r["rows"]:
        status = "OK" if row["present"] else "MISSING"
        print(f"{row['name']:<22} {('yes' if row['required'] else 'opt'):<4} {status:<9} {row['purpose']}")
        if row["present"] and row.get("path"):
            print(f"{'':<22} {'':<4} {'':<9} -> {row['path']}")
        if row["required"] and not row["present"]:
            missing.append(row)
    print("-" * 78)
    print("optional:  " + ", ".join(f"{k}={'OK' if v else '-'}" for k, v in r["optional"].items()))
    print("tokens:    " + ", ".join(f"{k}={'set' if v else '-'}" for k, v in r["tokens"].items())
          + "   (set them in .env; see .env.example)")
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
