#!/usr/bin/env python
"""
run_phase1.py -- drive the Phase-1 batch: the 49 validated ArsR/MerR candidates from the Rondon SSN
paper's four genomes through tf-operator-predictor.

Resumable: a run whose bundle manifest is complete is skipped, so the batch can be re-launched after
an interrupt without losing work. Individual failures are logged; the process returns non-zero after
all scheduled candidates finish so downstream stages cannot consume a partial run.

  TFOP_RUN_TAG=<tag> python analysis/run_phase1.py --all
  python analysis/run_phase1.py --jobs 3        # parallelism (default 2)
  python analysis/run_phase1.py --no-fold       # skip the ESMFold2 apo-dimer (much faster)
  python analysis/run_phase1.py --dry-run
"""
from __future__ import annotations
import argparse, json, os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from project_config import PREDICTOR_ROOT, PYTHON

PROJ = Path(__file__).resolve().parent.parent
REPO = PREDICTOR_ROOT
PY = PYTHON
sys.path.insert(0, str(REPO))
from predictor.provenance import runtime_provenance  # noqa: E402

RUNTIME_SHA256 = runtime_provenance()["runtime_sha256"]
sys_path = str(PROJ / "analysis")
if sys_path not in sys.path:
    sys.path.insert(0, sys_path)
import manifest as MF  # noqa: E402
import run_paths as RP  # noqa: E402  (the run tag must be resolved before any output is selected)

JOBS = RP.JOBS
LOGS = RP.LOGS
LOGS.mkdir(parents=True, exist_ok=True)
PROGRESS = LOGS / "phase1_progress.tsv"


def bundle_status(name: str) -> str:
    """Return missing, incomplete, stale, or complete for one output bundle."""
    bundle = JOBS / name
    manifest = bundle / "bundle_manifest.json"
    if not bundle.exists():
        return "missing"
    try:
        payload = json.loads(manifest.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return "incomplete"
    complete = payload.get("status") == "complete" and all((bundle / rel).is_file() for rel in (
        "dossier.json", "operators_ranked.json", "binding_sites.tsv", "regulon.tsv"
    ))
    if not complete:
        return "incomplete"
    recorded = (payload.get("runtime") or {}).get("runtime_sha256")
    return "complete" if recorded == RUNTIME_SHA256 else "stale"


def bundle_done(name: str) -> bool:
    return bundle_status(name) == "complete"


def run_one(r: dict, *, fold: str | None, timeout: int, offline: bool = False) -> dict:
    name = r["run_name"]
    t0 = time.time()
    if bundle_done(name):
        return {"run": name, "status": "skip(done)", "secs": 0}
    if bundle_status(name) == "stale":
        return {"run": name, "status": "STALE(archive-first)", "secs": 0}
    # `--mode production` was dropped: the pipeline no longer has modes, and a production run
    # (online + agnostic of pre-gathered local data) is now the only behaviour there is.
    cmd = [PY, "-m", "predictor.run_cli", "predict", str(PROJ / r["fasta"]),
           "--name", name, "--organism", r["organism_acc"]]
    if r["family_flag"] != "(auto)":
        cmd += ["--family", r["family_flag"]]
    if fold:
        cmd += [fold]
    if offline:
        # The predictor's own control. Without a way to reach it from here, a batch on a machine
        # with no network -- or during an NCBI outage -- has to bypass the driver entirely, which
        # loses the manifest gate, the staleness guard and the progress record.
        cmd += ["--offline"]
    env = dict(os.environ, KMP_DUPLICATE_LIB_OK="TRUE",
               PREDICTOR_OUTPUT_DIR=str(JOBS.parent), TFOP_JOBS=str(JOBS))
    log = LOGS / f"{name}.log"
    try:
        with log.open("w", encoding="utf-8", errors="replace") as fh:
            p = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, cwd=str(REPO),
                               env=env, timeout=timeout)
        rc = p.returncode
    except subprocess.TimeoutExpired:
        return {"run": name, "status": f"TIMEOUT({timeout}s)", "secs": round(time.time() - t0)}
    except Exception as e:
        return {"run": name, "status": f"ERROR:{type(e).__name__}", "secs": round(time.time() - t0)}
    ok = bundle_done(name)
    return {"run": name, "status": ("ok" if ok else f"FAIL(rc={rc},no-bundle)"), "secs": round(time.time() - t0)}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true",
                    help="deprecated compatibility flag; the full manifest is now the default")
    ap.add_argument("--family", help="run only rows whose family == this (e.g. Fur)")
    ap.add_argument("--focus", action="store_true", help="the 57 focus-family TFs (ArsR, MerR, Fur)")
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--timeout", type=int, default=5400, help="per-run seconds (default 90 min)")
    ap.add_argument("--no-fold", action="store_true")
    ap.add_argument("--fold", action="store_true")
    ap.add_argument("--offline", action="store_true",
                    help="pass --offline to every prediction: no network sources. Produces a real "
                         "bundle from packaged data alone, so it is also the fast way to exercise "
                         "the chain end to end without waiting on the NCBI BLAST queue.")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--manifest", default=None,
                    help="manifest under analysis/ (default: run_manifest_<TFOP_RUN_TAG>.csv)")
    a = ap.parse_args(argv)
    fold = "--no-fold" if a.no_fold else ("--fold" if a.fold else None)

    mpath = Path(a.manifest) if a.manifest else RP.MANIFEST
    if not mpath.is_absolute():
        mpath = PROJ / "analysis" / mpath
    # Validated through the SHARED module, not just read. Running this script directly is the
    # normal way to resume a long batch by hand, and it used to skip every check `run_pipeline.py`
    # performs -- so a bad row surfaced hours later as one failed candidate among many.
    rows = MF.load_or_exit(mpath, PROJ)
    print(f"manifest: {mpath.name} ({len(rows)} rows, validated)")
    if a.focus:
        todo = [r for r in rows if r["family"] in ("ArsR", "MerR", "Fur")]
    elif a.family:
        todo = [r for r in rows if r["family"] == a.family]
    else:
        todo = rows
    pending = [r for r in todo if not bundle_done(r["run_name"])]
    stale = [r["run_name"] for r in todo if bundle_status(r["run_name"]) == "stale"]
    print(f"phase 1: {len(todo)} runs selected | {len(todo)-len(pending)} already done | {len(pending)} pending")
    print(f"  repo={REPO}  bundles={JOBS}  workers={a.jobs}  "
          f"fold={fold or 'default(auto)'}  timeout={a.timeout}s"
          f"{'  OFFLINE' if a.offline else ''}")
    if stale:
        print(f"refusing to mix runtimes: {len(stale)} completed bundle(s) are stale or have no "
              "runtime provenance")
        for name in stale[:5]:
            print(f"  {name}")
        print("archive the run's bundles with run_pipeline.py --archive-as <tag>, then resume")
        return 2
    if a.dry_run:
        for r in pending[:6]:
            print("   would run:", r["run_name"], "|", r["family_flag"], "|", r["organism_acc"])
        return 0
    if not PROGRESS.exists():
        PROGRESS.write_text("run\tstatus\tsecs\n", encoding="utf-8")

    t0 = time.time(); done = 0; failures = []
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = {ex.submit(run_one, r, fold=fold, timeout=a.timeout, offline=a.offline): r
                for r in pending}
        for f in as_completed(futs):
            res = f.result(); done += 1
            with PROGRESS.open("a", encoding="utf-8") as fh:
                fh.write(f"{res['run']}\t{res['status']}\t{res['secs']}\n")
            if res["status"] not in ("ok", "skip(done)"):
                failures.append(res)
            el = time.time() - t0
            rate = el / max(1, done)
            print(f"[{done}/{len(pending)}] {res['status']:22s} {res['secs']:5d}s  {res['run']}"
                  f"   (elapsed {el/60:.0f}m, ~{rate*(len(pending)-done)/60:.0f}m left)", flush=True)
    print(f"\nphase 1 batch finished in {(time.time()-t0)/60:.0f} min")
    if failures:
        print(f"FAILED: {len(failures)} candidate(s) did not produce a complete bundle; "
              f"see {LOGS} and resume this tag after correcting them.")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
