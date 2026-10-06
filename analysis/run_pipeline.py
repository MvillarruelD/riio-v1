#!/usr/bin/env python
"""run_pipeline.py -- the single driver: candidates in, report out.

This replaces the hand-run sequence of `run_phase1.py` -> `aggregate_phase1.py` ->
`regenerate_regulons.py` -> the figure scripts -> `build_report_v5_docx.py` -> `make_share_bundle.py`.
It does not reimplement any of them; it runs them in order, with the two things the manual chain never
had.

**1. It refuses to start when there is NOTHING LEFT TO DO.**

HANDOFF 5.2: `run_phase1.py` and `regenerate_regulons.py` both skip any candidate that already has
output. Pointed at a COMPLETE run `jobs/`, a "rerun" therefore finishes in seconds, reports
success, and changes nothing. That has already cost this project a cycle.

The test is "nothing pending", not "something is present". Those are different, and using the second
was a bug: it refused every legitimate resume, and since the only way past it was `--allow-stale`,
reaching for that flag then waved the real no-op through as well -- the guard rejected the safe case
and passed the dangerous one. A partially-complete run is a RESUME and needs no flag; `run_phase1`
skips what is finished and works on the rest, so interrupting a long batch to change `--jobs` or
`--timeout` is normal and cheap.

**2. It refuses to run stages whose inputs are stale.**

Each stage declares what it consumes. A figure or report stage will not run over bundles that predate
the current predictor commit unless `--allow-stale` says so explicitly, because a report built from
half-old bundles is worse than no report: it looks complete.

    python analysis/run_pipeline.py --check                 # what would run, and what blocks it
    python analysis/run_pipeline.py --tag run7 --archive-as retry1  # archive this run's bundles
    python analysis/run_pipeline.py --all --jobs 4          # the batch, then everything downstream
    python analysis/run_pipeline.py --from regulons         # resume at a stage
    python analysis/run_pipeline.py --only batch --dry-run
    python analysis/run_pipeline.py --tag t --only batch --all   # resumes; no flag needed
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import manifest as MF
from project_config import PREDICTOR_ROOT, PYTHON

PROJ = Path(__file__).resolve().parent.parent
ANALYSIS = PROJ / "analysis"
REPO = PREDICTOR_ROOT
PY = PYTHON
#: Which candidate list the run is about. Mutable so `--manifest` can point the WHOLE driver -- the
#: archive guard included -- at a genome-first manifest. If only the batch honoured it, `--archive-as`
#: would compute "the run's bundles" from the wrong list and leave the previous run's bundles in
#: place for the batch to skip.
#:
#: The default is derived from `TFOP_RUN_TAG`, and there is deliberately NO fallback list. The old
#: fallback was `run_manifest.csv`, the survey's frozen 150 -- and because 132 of the 140 production
#: candidates share a name with one of those 150, pointing at it does not fail. `--check` under an
#: exported tag reported "132 complete, 18 pending" for a finished 140-candidate run, which reads as
#: an interrupted batch and is not one. `--tag` overrides this in main().
_ENV_TAG = os.environ.get("TFOP_RUN_TAG", "").strip()
MANIFEST = ANALYSIS / f"run_manifest_{_ENV_TAG}.csv" if _ENV_TAG else None
# A placeholder used only before argument resolution. `main()` replaces it from `run_paths.env()`;
# no production operation is allowed to proceed without a run tag.
JOBS = Path(os.environ.get("TFOP_JOBS", ANALYSIS / f"run_{_ENV_TAG or '__unset__'}" / "jobs"))

STAGES = [
    ("discover", "run_discovery.py", "genomes -> candidates -> manifest (BITACORA)"),
    ("batch", "run_phase1.py", "the per-TF predictions -- NOVEL TFs run 30-120 min EACH; raise --timeout above 7200 for Vibrio"),
    ("aggregate", "aggregate_phase1.py", "bundles -> one master table"),
    ("regulons", "regenerate_regulons.py", "per-regulator regulons from each bundle"),
    ("figures", None, "the run-derivable figures only (see FIGURE_SCRIPTS)"),
    ("report", "build_report_v5_docx.py", "the .docx for the authors"),
    ("share", "make_share_bundle.py", "the share-ready folder"),
]
#: Figures THIS run can actually produce -- i.e. whose inputs are bundles, the aggregate table or the
#: regulon summary, all written by the batch/aggregate/regulons stages.
#:
#: This is the whole figure set. The scripts that fanned out to one-off statistics files no
#: pipeline stage writes -- `operator_null.json`, the `af3_*` structural tables, the `motif_*`
#: benchmark sweeps -- have been removed rather than left to fail, because on a dirty disk they
#: do not fail: they quietly draw an old number into a new report.
#: `classify_figures.py` re-derives this split; run it if the figure set changes.
FIGURE_SCRIPTS = ("fig_pipeline.py",        # F1  schematic
                  "fig_corroboration.py",   # F2  per-family metal calls vs the survey
                  "fig_evidence_ledger.py",  # F6  evidence tier per family + axis coverage
                  "fig_operator_cds.py",    # F3  operator -> CDS start vs genomic null
                  "fig_logos.py",           # F4a/F4c family and per-TF operator logos
                  "fig_cluster_motifs.py",  # F4b within- vs between-cluster motif sharing
                  "fig_four_genomes.py",    # F22
                  # ORDER MATTERS for the last three, and it is a chain, not a preference:
                  #   fig5_complement WRITES run_paths.COMPLEMENT_TSV
                  #   fig5_panelE     READS it to draw F23 (and exits if it is absent)
                  #   fig5_compose    READS F23 to build the F24 composite
                  # Listed panelE-before-complement, the stage fails on a clean run and "succeeds"
                  # on a dirty one by reading the previous run's leftovers.
                  "fig5_complement.py",     # writes the per-regulator table
                  "fig5_panelE.py",         # F23
                  "fig5_compose.py")        # F24 composite



def manifest_names() -> list[str]:
    if MANIFEST is None or not MANIFEST.exists():
        return []                      # no tag yet, or discovery has not run: both normal states
    with MANIFEST.open(encoding="utf-8") as fh:
        return [r["run_name"] for r in csv.DictReader(fh)]


def manifest_issues() -> list[str]:
    """Structural/input errors that would make this manifest unsafe to execute.

    Delegated to `analysis/manifest.py` so the driver, `run_phase1.py` and `aggregate_phase1.py`
    apply ONE definition. When each had its own, running a stage script directly skipped the
    driver's checks entirely.
    """
    if MANIFEST is None:
        return ["no manifest selected"]
    return MF.validate(MANIFEST, PROJ)


def report_manifest_validation() -> bool:
    if MANIFEST is None:
        print("manifest validation   : FAILED (no manifest selected)")
        return False
    return MF.report(MANIFEST, PROJ)


def candidate_bundles() -> list[Path]:
    """Bundles that belong to the run, i.e. the manifest's -- never the benchmark ones."""
    if not JOBS.is_dir():
        return []
    names = set(manifest_names())
    return [p for p in JOBS.iterdir() if p.is_dir() and p.name in names]


def is_done(name: str) -> bool:
    """Complete in the sense `run_phase1.bundle_done` means it: the bundle manifest is complete.

    A bundle directory or one headline JSON is not sufficient. The predictor publishes bundles
    atomically and records the required artifacts in `bundle_manifest.json`; the orchestration uses
    that contract rather than inferring completion from directory presence.
    """
    bundle = JOBS / name
    try:
        complete = json.loads((bundle / "bundle_manifest.json").read_text(
            encoding="utf-8")).get("status") == "complete"
    except (OSError, ValueError):
        return False
    return complete and all((bundle / rel).is_file() for rel in (
        "dossier.json", "operators_ranked.json", "binding_sites.tsv", "regulon.tsv"
    ))


_CURRENT_RUNTIME: dict[str, str] | None = None


def current_runtime() -> dict[str, str]:
    """Fingerprint the predictor exactly as the batch interpreter imports it."""
    global _CURRENT_RUNTIME
    if _CURRENT_RUNTIME is not None:
        return _CURRENT_RUNTIME
    code = ("import json; from predictor.provenance import runtime_provenance; "
            "print(json.dumps(runtime_provenance()))")
    env = dict(os.environ)
    # Running from the predictor checkout makes this work for both editable and ordinary source
    # installs, and matches the working directory used by the batch itself.
    result = subprocess.run([PY, "-c", code], cwd=str(REPO), env=env, capture_output=True,
                            text=True, timeout=60)
    if result.returncode != 0:
        raise RuntimeError(f"cannot fingerprint predictor runtime: {(result.stderr or '').strip()[:300]}")
    try:
        _CURRENT_RUNTIME = json.loads(result.stdout.strip().splitlines()[-1])
    except (IndexError, json.JSONDecodeError) as exc:
        raise RuntimeError("predictor runtime returned invalid provenance") from exc
    return _CURRENT_RUNTIME


def bundle_runtime(name: str) -> dict[str, str] | None:
    try:
        payload = json.loads((JOBS / name / "bundle_manifest.json").read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    runtime = payload.get("runtime")
    return runtime if isinstance(runtime, dict) else None


def stale_bundles() -> list[tuple[str, str]]:
    """Completed manifest bundles whose code/data fingerprint differs from this runtime."""
    expected = current_runtime().get("runtime_sha256")
    stale: list[tuple[str, str]] = []
    for name in manifest_names():
        if not is_done(name):
            continue
        recorded = (bundle_runtime(name) or {}).get("runtime_sha256")
        if recorded != expected:
            stale.append((name, recorded or "unrecorded"))
    return stale


def validate_bundle_inputs(*, allow_stale: bool) -> bool:
    """Fail closed before any stage that consumes completed prediction bundles."""
    done, pending = run_progress()
    if pending:
        print(f"    refusing downstream work: {len(pending)}/{len(done) + len(pending)} bundle(s) "
              "are incomplete")
        return False
    stale = stale_bundles()
    if stale and not allow_stale:
        print(f"    refusing downstream work: {len(stale)} bundle(s) were produced by another "
              "code/data runtime")
        print(f"    current runtime: {current_runtime()['runtime_sha256'][:16]}")
        for name, fingerprint in stale[:5]:
            print(f"      {name}: {fingerprint[:16]}")
        print("    archive and rerun those bundles, or pass --allow-stale for explicit historical inspection")
        return False
    if stale:
        print(f"    WARNING: --allow-stale accepted {len(stale)} bundle(s) from another runtime")
    return True


def run_progress() -> tuple[list, list]:
    """(done, pending) over the manifest, by the batch's own definition of done."""
    names = manifest_names()
    done = [n for n in names if is_done(n)]
    return done, [n for n in names if n not in set(done)]


def foreign_bundles() -> list[Path]:
    """Directories in this run's job root that are not named in its manifest.

    Production never reads these. Their presence usually means an explicit `TFOP_JOBS` override was
    pointed at a legacy shared directory; benchmark outputs belong under `analysis/benchmarking/`.
    """
    if not JOBS.is_dir():
        return []
    names = set(manifest_names())
    return [p for p in JOBS.iterdir() if p.is_dir() and p.name not in names]


def predictor_commit() -> str:
    try:
        out = subprocess.run(["git", "-C", str(REPO), "rev-parse", "--short", "HEAD"],
                             capture_output=True, text=True, timeout=30)
        return out.stdout.strip() or "unknown"
    except Exception:
        return "unknown"


def archive(tag: str) -> int:
    """Move the run's candidate bundles aside. Benchmarks stay: the scorers read them."""
    dest = JOBS.parent / f"{JOBS.name}_{tag}_archive"
    if dest.exists():
        print(f"  ! {dest} already exists -- pick another tag rather than merging two runs into one")
        return 1
    cands = candidate_bundles()
    if not cands:
        print(f"  nothing to archive: no manifest bundles in {JOBS}")
        return 0
    dest.mkdir(parents=True)
    for p in cands:
        shutil.move(str(p), str(dest / p.name))
    print(f"  archived {len(cands)} candidate bundle(s) -> {dest}")
    print(f"  left {len(foreign_bundles())} unrelated bundle(s) in place")
    return 0


#: Analysis outputs archived on 2026-08-20 because they hold RUN-3 / RUN-4 data. Thirty figure and
#: report scripts still read them, and surfacing that is the point of the check below: a run-3
#: statistic inside a new report is the expensive kind of silent error, because the report still looks
#: complete. `ARCHIVE_MANIFEST.md` at the archive root says where each one went.
#: The seven scripts of the figure/report/share/regulon stages were re-pointed at `run_paths.py` on
#: 2026-08-20, so those stages no longer read archived data. What remains on this list is the ~22
#: one-off analysis scripts that still name a historical run; they are not part of a production run,
#: and the check below exists so that is a visible fact rather than a surprise.
ARCHIVED_INPUTS = ("run3_clean", "run4_regulon", "tetr_bench", "talk_group_meeting")


def stale_input_scripts() -> dict:
    """{script: [archived dirs it reads]} -- the stages that would otherwise use superseded data."""
    import re

    pat = re.compile("|".join(ARCHIVED_INPUTS))
    out = {}
    for py in sorted(list(ANALYSIS.glob("*.py")) + list(ANALYSIS.glob("*/*.py"))):
        # These two hold the archived names as DATA -- run_pipeline's own guard list, and the
        # checker whose entire job is to detect them. Scanning them reports the detector as the
        # defect.
        if py.name in ("run_pipeline.py", "check_run_paths.py"):
            continue
        try:
            hits = sorted(set(pat.findall(py.read_text(encoding="utf-8", errors="replace"))))
        except Exception:
            continue
        if hits:
            out[py.name] = hits
    return out


#: Windows refuses a path longer than this unless long-path support is on, and the failure surfaces
#: as `FileNotFoundError: [WinError 206]` from deep inside a stage -- after the expensive work.
_WINDOWS_MAX_PATH = 260
#: The deepest path a bundle creates, relative to `jobs/`: the atomic staging directory (the run name
#: plus a random suffix) and the AF3 hand-off tree beneath it. Measured from a real run.
_DEEPEST_BUNDLE_SUFFIX = r"\.{name}.stage-abcdefgh\af3_validation\af3_results\seed-1_sample-0\model.cif"


def path_headroom() -> tuple[int, int, str]:
    """(longest path this run will create, headroom before the OS limit, the candidate driving it).

    Worth checking because the cost is asymmetric: the check takes no time, and the failure it
    prevents arrives 65 minutes in, after the prediction has succeeded and while the bundle is being
    published. A run root only ~70 characters deeper than this repository's own layout is enough to
    trigger it, which a clone into a nested directory easily is.
    """
    names = manifest_names() or ["x" * 31]
    longest = max(names, key=len)
    total = len(str(JOBS)) + len(_DEEPEST_BUNDLE_SUFFIX.format(name=longest))
    return total, _WINDOWS_MAX_PATH - total, longest


def check() -> int:
    """Report whether a genuine run is possible, and what would silently no-op."""
    cands, foreign = candidate_bundles(), foreign_bundles()
    print(f"predictor commit      : {predictor_commit()}")
    try:
        print(f"runtime fingerprint   : {current_runtime()['runtime_sha256']}")
    except RuntimeError as exc:
        print(f"runtime fingerprint   : ERROR -- {exc}")
        return 2
    print(f"jobs                  : {JOBS}")
    longest_path, headroom, driver = path_headroom()
    if sys.platform == "win32" and headroom < 0:
        print(f"path length           : FAIL -- {longest_path} chars, {-headroom} over the "
              f"{_WINDOWS_MAX_PATH}-character Windows limit")
        print(f"                        driven by {driver}; the bundle's AF3 hand-off tree is the "
              "deepest part")
        print("                        Use a shorter run root (TFOP_RUN_ROOT), or enable Windows "
              "long-path support.")
        print("                        A run started now would predict successfully and then fail "
              "while PUBLISHING the bundle.")
        return 2
    print(f"path length           : {longest_path} chars, {headroom} under the "
          f"{_WINDOWS_MAX_PATH}-character Windows limit"
          f"{'  (tight)' if headroom < 40 else ''}")
    print(f"job directories       : {len(cands)} candidate + {len(foreign)} unrelated")
    if foreign:
        print("                        unrelated bundles are ignored; keep benchmarks outside run jobs")
    if MANIFEST is None:
        print("manifest              : NONE -- set TFOP_RUN_TAG or pass --tag")
        print("                        (without it, progress below is meaningless)")
    else:
        print(f"manifest              : {MANIFEST.name} -- {len(manifest_names())} candidates")
    if not report_manifest_validation():
        return 2
    stale = stale_input_scripts()
    if stale:
        print(f"\nARCHIVED INPUTS: {len(stale)} analysis script(s) read run-3/run-4 data that was")
        print("archived on 2026-08-20. They will FAIL rather than quietly draw a figure from the old")
        print("run -- which is the intended behaviour: a run-3 statistic inside a new report is the")
        print("expensive kind of silent error, because the report still looks complete.")
        print("Re-point them at the new run's outputs before the figures/report stages.")
        for name, dirs in list(stale.items())[:6]:
            print(f"    {name:<30} reads {', '.join(dirs)}")
        if len(stale) > 6:
            print(f"    ... and {len(stale) - 6} more")

    done, pending = run_progress()
    stale_bundles_found = stale_bundles()
    partial = len(cands) - len(done)          # a directory but no operators_ranked.json
    print(f"progress              : {len(done)} complete, {len(pending)} pending"
          f"{f' ({partial} started but unfinished, will be redone)' if partial > 0 else ''}")
    print(f"runtime provenance    : {len(done) - len(stale_bundles_found)} current, "
          f"{len(stale_bundles_found)} stale/unrecorded")

    if not done:
        print("\nREADY: no candidate has output yet, so the batch will do real work.")
        return 0
    if pending:
        print(f"\nRESUME: {len(pending)} candidate(s) still to run. Start the batch normally -- "
              f"run_phase1\n        skips the {len(done)} finished ones. No flag is needed, and "
              f"changing --jobs or\n        --timeout between resumes is fine.")
        return 0
    print(f"\nBLOCKED: all {len(done)} manifest candidates already have output.")
    print("Both run_phase1.py and regenerate_regulons.py SKIP any candidate that already has output,")
    print("so a run started now would finish quickly, report success, and change nothing.")
    print("\nArchive them first:")
    print("    python analysis/run_pipeline.py --archive-as <tag>")
    print("\n(or, by hand -- note this must move ONLY the manifest's bundles, never the benchmarks)")
    return 1


#: Set by main() so every stage writes into the SAME run directory. Without this the figure and
#: report stages fall back to `analysis/run_current/` while the batch fills a tagged manifest, and
#: the two halves of one run end up in different places.
RUN_TAG = ""


def run_paths_env(tag: str) -> dict:
    """`run_paths.env()` evaluated for `tag`.

    Imported in a subprocess rather than in-process because `run_paths` resolves its paths at import
    time from the environment: importing it here once would freeze them at whatever the tag was on
    the first call, and a second run in the same process would silently reuse the first run's
    directories.
    """
    code = ("import sys,json;"
            "sys.path.insert(0, r'%s');"
            "import run_paths;print(json.dumps(run_paths.env()))" % str(ANALYSIS))
    try:
        child_env = dict(os.environ, TFOP_RUN_TAG=tag)
        if MANIFEST is not None:
            child_env["AN_MANIFEST"] = str(MANIFEST)
        out = subprocess.run([PY, "-c", code], capture_output=True, text=True, timeout=60,
                             env=child_env)
        if out.returncode == 0 and out.stdout.strip():
            import json
            return json.loads(out.stdout.strip().splitlines()[-1])
        print(f"    (run_paths.env failed: {(out.stderr or '').strip()[:200]})")
    except Exception as exc:
        print(f"    (run_paths.env failed: {exc})")
    raise RuntimeError(f"could not resolve paths for run tag {tag!r}")


def write_figure_manifest(produced: list, failed: list, missing: list | None = None) -> None:
    """Record which figures this run has, and which it deliberately does not.

    Without this, a reader of the report cannot tell an absent figure from one that was never
    attempted -- and "F7 is missing" reads as an oversight when it actually means "the AF3 structural
    analysis was not part of this run".
    """
    import json

    env = run_paths_env(RUN_TAG) if RUN_TAG else {}
    out_dir = Path(env.get("AN_OUT", ANALYSIS / "results"))
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "figures_manifest.json").write_text(json.dumps({
            "run_tag": RUN_TAG,
            "predictor_commit": predictor_commit(),
            "produced": produced,
            "failed": failed,
            # Split out of `failed` so a reader can tell "the script ran and errored" from "the
            # script is not in the tree at all". Both fail the stage; they are fixed differently.
            "missing_scripts": list(missing or []),
            "configured": list(FIGURE_SCRIPTS),
            "note": ("`configured` is the complete production figure set. A figure absent from "
                     "`produced` was NOT quietly skipped -- it is listed in `failed`, the stage "
                     "returned non-zero, and no report was built from a partial set. Figures "
                     "needing statistics that no pipeline stage writes are not configured here at "
                     "all, and must never be filled in from an older run's files."),
        }, indent=2), encoding="utf-8")
        print(f"    wrote {out_dir / 'figures_manifest.json'}")
    except Exception as exc:
        print(f"    (could not write the figure manifest: {exc})")


def run_script(name: str, extra: list[str], *, dry: bool) -> int:
    cmd = [PY, str(ANALYSIS / name), *extra]
    print(f"    $ {' '.join(cmd[1:])}")
    if dry:
        return 0
    env = dict(os.environ, KMP_DUPLICATE_LIB_OK="TRUE")
    if RUN_TAG:
        # Not just TFOP_RUN_TAG: the ~24 one-off figure scripts each default to `analysis/run3_clean`
        # and are re-pointed by AN_RUN / AN_FIG / AN_OUT / AN_VALID / TFOP_JOBS instead. Setting the
        # tag alone moves the seven report-chain scripts and leaves the rest reading run-3 data.
        env.update(run_paths_env(RUN_TAG))
    t0 = time.time()
    rc = subprocess.run(cmd, cwd=str(PROJ), env=env).returncode
    print(f"    -> rc={rc} in {time.time() - t0:.0f}s")
    return rc


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="report readiness and exit")
    ap.add_argument("--manifest", default=None,
                    help="candidate list under analysis/ (default: run_manifest_<tag>.csv). Pass "
                         "this only to override the tag's own list.")
    ap.add_argument("--tag", help="discovery tag; with --only/--from discover, names the manifest "
                                  "that stage writes (run_manifest_<tag>.csv)")
    ap.add_argument("--archive-as", metavar="TAG",
                    help="move this run's candidate bundles to jobs_<TAG>_archive and exit")
    ap.add_argument("--all", action="store_true", help="pass --all to the batch (all 150)")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--timeout", type=int, default=7200)
    # Structure OFF is the default because the canonical run is structure OFF (every 0917 bundle
    # records allow_folds=False). The old default passed --fold, so the runbook's command silently
    # produced a run that could not be compared with the canonical one.
    fold_opt = ap.add_mutually_exclusive_group()
    fold_opt.add_argument("--fold", dest="fold", action="store_true", default=False,
                          help="also fold each candidate (apo-dimer; needs BIOHUB_TOKEN). NOT the "
                               "canonical configuration -- a folded run is not comparable to it")
    fold_opt.add_argument("--no-fold", dest="fold", action="store_false",
                          help="no structure prediction (the default; accepted for old commands)")
    ap.add_argument("--offline", action="store_true",
                    help="run the batch with no network sources (passes --offline to every "
                         "prediction). For a machine without network, an NCBI outage, or a fast "
                         "end-to-end check that does not wait on the BLAST queue.")
    ap.add_argument("--only", choices=[s[0] for s in STAGES])
    ap.add_argument("--from", dest="from_stage", choices=[s[0] for s in STAGES])
    ap.add_argument("--dry-run", action="store_true",
                    help="PLAN the run: resolve the tag, validate the manifest and every FASTA, "
                         "confirm each configured figure script exists, and print every stage and "
                         "command -- without executing any of them. It deliberately does NOT "
                         "require the bundles a real run would have produced by then, so it works "
                         "on a run that has not started; it reports pending counts instead. Exit 0 "
                         "means the plan is executable, not that the products exist.")
    ap.add_argument("--allow-stale", action="store_true",
                    help="run downstream stages over bundles older than the current predictor "
                         "commit. NOT needed to resume a batch -- a partially-complete run resumes "
                         "on its own.")
    a = ap.parse_args(argv)

    global JOBS, MANIFEST, RUN_TAG
    # Tag precedence: --tag, then a --manifest named run_manifest_<tag>.csv, then TFOP_RUN_TAG.
    # An explicit --manifest always wins for the PATH, but it must not silently un-set the tag.
    RUN_TAG = a.tag or (
        a.manifest[len("run_manifest_"):-len(".csv")]
        if (a.manifest or "").startswith("run_manifest_") and a.manifest.endswith(".csv")
        else _ENV_TAG)
    if a.manifest:
        MANIFEST = ANALYSIS / a.manifest
    elif RUN_TAG:
        MANIFEST = ANALYSIS / f"run_manifest_{RUN_TAG}.csv"
    else:
        MANIFEST = None

    if RUN_TAG:
        resolved = run_paths_env(RUN_TAG)
        JOBS = Path(resolved["TFOP_JOBS"])

    if a.archive_as:
        if not RUN_TAG:
            print("--archive-as requires --tag, a tagged --manifest, or TFOP_RUN_TAG")
            return 2
        return archive(a.archive_as)
    if a.check:
        return check()

    names = [s[0] for s in STAGES]
    todo = [a.only] if a.only else names[names.index(a.from_stage):] if a.from_stage else names
    # A full dry-run starting at discovery cannot validate the manifest yet: discovery is the
    # stage that would create it, and dry-run deliberately executes no stages.  Defer only this
    # future-manifest check.  A real run, or a batch-only dry-run, still fails closed below.
    manifest_planned_by_discovery = (
        a.dry_run
        and "discover" in todo
        and "batch" in todo
        and todo.index("discover") < todo.index("batch")
        and MANIFEST is not None
        and not MANIFEST.exists()
    )

    # The guard, and what it is actually for. Landmine 5.2 is that a batch over an ALREADY-COMPLETE
    # output directory finishes in seconds, reports success and changes nothing. The danger is
    # therefore "nothing left to do", not "something is already there" -- and the old test was the
    # latter, so it refused every legitimate resume while still passing the no-op through whenever
    # anyone reached for --allow-stale to get past it. That is exactly backwards.
    #
    # A partially-complete run is a RESUME: `run_phase1` skips the finished candidates and does real
    # work on the rest. It needs no flag, and interrupting a long batch to change --jobs or --timeout
    # is a normal thing to want to do.
    if "batch" in todo:
        done, pending = run_progress()
        stale = stale_bundles() if done else []
        if stale:
            print("=== refusing to resume over stale bundles ===\n")
            print(f"{len(stale)} completed bundle(s) were made by a different or unrecorded runtime.")
            print("Archive this run's bundles before rerunning so the batch cannot mix implementations:")
            print("    python analysis/run_pipeline.py --archive-as <tag>")
            return 1
        if done and not pending:
            print("=== refusing to start ===\n")
            print(f"All {len(done)} manifest candidates already have output, so the batch would skip")
            print("every one of them, finish in seconds and report success without changing anything.")
            print("\nArchive them first, then re-run:")
            print("    python analysis/run_pipeline.py --archive-as <tag>")
            # Deliberately NOT bypassable by --allow-stale. That flag is about running downstream
            # stages over older bundles; letting it also wave through a no-op batch is what made the
            # old guard worse than useless, since the flag was the documented way past it.
            return 1
        if done:
            print(f"[resume] {len(done)} of {len(done) + len(pending)} candidates already complete; "
                  f"{len(pending)} to run\n")

    # The figure/report/share stages write into analysis/run_<tag>/; create it once here rather than
    # leaving each script to discover it is missing halfway through a long run.
    if RUN_TAG:
        try:
            resolved = run_paths_env(RUN_TAG)
            for key in ("TFOP_JOBS", "TFOP_RUN_LOGS", "AN_FIG", "AN_OUT", "TFOP_RUN_REGULON"):
                Path(resolved[key]).mkdir(parents=True, exist_ok=True)
        except Exception as exc:
            print(f"    (could not pre-create the run directories: {exc})")

    print(f"predictor commit {predictor_commit()} · stages: {', '.join(todo)}"
          f"{' (dry run)' if a.dry_run else ''}\n")
    for stage, script, desc in STAGES:
        if stage not in todo:
            continue
        print(f"[{stage}] {desc}")
        # A dry run is a PLANNING run: it resolves the tag, validates the manifest and the static
        # inputs, and prints every stage and command -- but it must not require the bundles that
        # only the (simulated) batch would have produced. Demanding them made `--dry-run` on a new
        # run fail at the first downstream stage, which told the operator nothing about the plan.
        if stage not in {"discover", "batch"}:
            if a.dry_run:
                if manifest_planned_by_discovery:
                    print("    (plan only: candidate count is unknown until discovery; "
                          "a real run requires every bundle complete)")
                else:
                    done, pending = run_progress()
                    print(f"    (plan only: {len(done)} complete, {len(pending)} pending bundle(s); "
                          "a real run requires 0 pending)")
            elif not validate_bundle_inputs(allow_stale=a.allow_stale):
                return 1
        if stage == "discover":
            if not a.tag:
                print("    ! --tag is required to run discovery (it names the manifest)")
                return 2
            rc = run_script(script, ["--tag", a.tag], dry=a.dry_run)
        elif stage == "batch":
            if MANIFEST is None:
                print("  no candidate list: pass --tag or export TFOP_RUN_TAG")
                return 2
            if manifest_planned_by_discovery:
                print(f"    (plan only: {MANIFEST.name} will be created by discovery; "
                      "manifest validation is deferred)")
            elif not report_manifest_validation():
                return 2
            extra = ["--jobs", str(a.jobs), "--timeout", str(a.timeout),
                     "--manifest", MANIFEST.name]
            if a.offline:
                extra.append("--offline")
            if a.all:
                extra.append("--all")
            extra.append("--fold" if a.fold else "--no-fold")
            rc = run_script(script, extra, dry=a.dry_run)
        elif stage == "figures":
            # Run every panel so the manifest is complete, then fail the stage if any panel or the
            # publication-export QA fails. A partial figure set must never flow into the report.
            produced, failed, missing = [], [], []
            # Stamped BEFORE the first figure runs and handed to figure_qa as --since. Figure
            # scripts write side-car stats into results/, and without this cut-off the QA counts
            # those as the run's INPUTS -- making every figure generated before the last side-car
            # "stale" and the gate unpassable. See figure_qa._newest_input_mtime.
            figures_started = time.time()
            for f in FIGURE_SCRIPTS:
                # A configured production figure whose script is absent is a DEFECT, not an
                # omission. Skipping it left the report one panel short while the stage still
                # exited 0 -- the silent-partial-success failure this driver exists to prevent.
                if not (ANALYSIS / f).exists():
                    print(f"    ! {f}: configured figure script is MISSING from {ANALYSIS}")
                    missing.append(f)
                    failed.append(f)
                    continue
                (produced if run_script(f, [], dry=a.dry_run) == 0 else failed).append(f)
            rc = 1 if failed or not produced else 0
            print(f"\n    figures produced {len(produced)}, failed {len(failed)}"
                  f"{f' ({len(missing)} script(s) missing)' if missing else ''}")
            if failed:
                print(f"    FAILED: {', '.join(failed)}")
            if not a.dry_run:
                write_figure_manifest(produced, failed, missing)
                if rc == 0:
                    rc = run_script("figure_qa.py", ["--since", f"{figures_started:.3f}"],
                                    dry=False)
        else:
            rc = run_script(script, [], dry=a.dry_run)
        if rc:
            print(f"\n!! stage '{stage}' returned {rc} -- stopping so the failure is not carried "
                  f"into the report")
            return rc
        print()
    print("done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
