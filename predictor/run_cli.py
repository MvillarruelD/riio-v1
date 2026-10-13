#!/usr/bin/env python
"""
run_cli.py -- the friendly top-level CLI for the TF -> operator pipeline (installed as `tfop`).

One command from a TF protein sequence to a complete per-TF output folder (operators, structure, motif,
regulation, inducer, AF3 hand-off). Wraps `predictor.pipeline.run_novel` with sensible defaults and prints
a labelled map of every file it produced.

USAGE
  # first time: check the engines, the bundled data, then every module's self-test:
  tfop setup
  tfop audit
  tfop selftest

  # predict from a sequence, a one-record FASTA, or a protein accession (online by default):
  tfop predict mytf.fasta --name MyTF
  tfop predict NP_417751.1 --name ZntR --organism GCF_000005845.2 --out runs/ZntR
  tfop predict MARSADKQ... --name MyTF --family ArsR/SmtB --organism GCF_000013425.1

  # several proteins (multi-record FASTA, or a TSV with an `accession` column):
  tfop batch proteins.faa --organism GCF_000005845.2 --out runs/list

  # OR go genome-first (the production route): BITACORA discovery -> one bundle per candidate:
  tfop genome GCF_000005845.2 --out runs/ecoli --jobs 4    # -> runs/ecoli/summary.tsv + jobs/<TF>/
  tfop genome GCF_000005845.2 --out runs/ecoli --discover-only   # stop after manifest.csv
  tfop genome --fasta genome.fna --gff genome.gff --out runs/local   # a local genome + annotation

  # a quick Pfam-HMM census of a proteome or genome (not the production discovery route):
  tfop scan proteome.faa                                   # -> candidate TFs by family (fast)
  tfop scan genome.fna --families ArsR/SmtB,MerR,Fur       # nucleotide -> gene-called, metalloregulators

  # optionally add an apo-dimer fold (off by default; needs BIOHUB_TOKEN):
  tfop predict mytf.fasta --name MyTF --fold

  # after a run: rebuild its HTML/XLSX reports, or export the optional AlphaFold 3 jobs:
  tfop report runs/ecoli
  tfop af3 runs/ecoli --format both

(You can equivalently run `python -m predictor.run_cli <cmd>`.)
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

# the per-file guide printed after a run (filename -> what it holds)
_OUTPUT_MAP = [
    ("REPORT.html", "OPEN THIS: the human-readable report (verdict banner + every figure + table, sectioned)"),
    ("figures/", "all standalone figures (apo dimer, logos, operator track, inducer evidence, TSS, ...)"),
    ("operators/FINAL_operators.tsv", "the final operator list: source (natural/designed) + why each was chosen"),
    ("operators_ranked.json", "ranked operator recommendations + the single `primary` pick"),
    ("operators/all_operators.tsv", "every genome operator + regulated gene, mode, TSS, surrounding genes"),
    ("operators/aligned_operators_*.fasta", "the aligned operator sequences behind each logo (TF/homolog/gathered)"),
    ("tf/tf.fasta + tf_summary.json", "the TF sequence + identity (family, organism, locus, UniProt, cluster)"),
    ("dossier.json", "the full machine-readable record (every field; substrate for any view)"),
    ("bundle_manifest.json", "bundle completeness, required artifacts, and report warnings"),
    ("binding_sites.tsv", "every genome-wide putative operator with coordinates + calibrated p/q-values"),
    ("regulon.tsv", "predicted operons the TF regulates (operator-PWM genome scan -> operons)"),
    ("structure/structure.json", "apo-dimer fold (ESMFold2) + AFDB-monomer integrity QC (RMSD/pLDDT)"),
    ("ligand/inducer.json", "inferred inducer/effector + regulated-gene promoters"),
    ("genome/genome.json", "genome accession, TF locus, neighborhood, promoter occlusion"),
    ("regulation/per_hit_regulation.json", "per-hit TSS/occlusion plus ranked-vs-significant gene support"),
    ("af3_jobs.json", "OPTIONAL: TF-dimer + operator jobs you may fold on alphafoldserver.com by hand; "
                      "nothing in the pipeline reads the result back"),
]


_PROTEIN_ACCESSION = re.compile(r"^(?:[A-Z]{1,4}_\d+(?:\.\d+)?|[A-Z0-9]{6,10})$")


def _resolve_prediction_source(value: str, *, offline: bool = False) -> str:
    """Accept raw protein, one-record FASTA/path, or a protein accession at the CLI boundary."""
    from predictor.input_validation import read_protein
    try:
        read_protein(value)
        return value
    except ValueError as input_error:
        accession = str(value or "").strip()
        if not _PROTEIN_ACCESSION.fullmatch(accession) or not any(c.isdigit() for c in accession):
            raise input_error
        if offline:
            raise ValueError("protein-accession lookup needs the network; paste the sequence for --offline")
        from predictor.annotate.protein_seqs import fetch_protein_sequences
        found = fetch_protein_sequences([accession], allow_ncbi=True, verbose=True)
        sequence = found.get(accession.split(".")[0], "")
        if not sequence:
            raise ValueError(f"no protein sequence was found for accession {accession!r}")
        return sequence


def _print_output_map(bundle: Path) -> None:
    import glob as _glob
    print("\n" + "=" * 78)
    print("PREDICTION COMPLETE")
    print(f"OUTPUT FOLDER:  {bundle}")
    print("=" * 78)
    for fname, desc in _OUTPUT_MAP:
        token = fname.split(" ")[0].rstrip("/")          # first path token (drop "+ ...", trailing /)
        exists = bool(_glob.glob(str(bundle / token))) or (bundle / token).exists()
        mark = "+" if exists else " "
        print(f"  [{mark}] {fname:<40} {desc}")
    print("=" * 78)
    print(f"START HERE (human):       {bundle / 'REPORT.html'}")
    print(f"USE / TEST THESE:         {bundle / 'operators' / 'FINAL_operators.tsv'}")
    print(f"START HERE (programmatic): {bundle / 'operators_ranked.json'}")
    print("Interpretation: operators and regulons are ranked PREDICTED_UNVERIFIED hypotheses, "
          "not experimentally confirmed sites.")
    print("Legend: [+] written; [ ] optional/not produced. Missing optional files do not make the run fail.")


def cmd_predict(a) -> int:
    import os
    if a.out:
        os.environ["PREDICTOR_OUTPUT_DIR"] = str(Path(a.out).expanduser().resolve())
    from predictor.config import RunConfig
    from predictor import api
    try:
        source = _resolve_prediction_source(a.seq, offline=a.offline)
    except (ValueError, OSError) as exc:
        print(f"INPUT ERROR: {exc}")
        return 2
    # Folding is opt-in. Its old environment-dependent default made two users get different bundles
    # from the same input merely because one happened to have a token configured.
    fold = False if a.no_fold else (True if a.fold else None)
    cfg = api._resolve_folds(RunConfig(allow_metalnet=not a.no_metalnet,
                                       emit_af3_jobs=not a.no_af3,
                                       allow_online=not a.offline), fold)
    print(f"=== predict: {a.name}  [fold={cfg.allow_folds} metalnet={cfg.allow_metalnet} "
          f"af3={cfg.emit_af3_jobs} online={cfg.allow_online}] ===")
    try:
        res = api.predict(source, name=a.name, family=a.family, organism=a.organism, genome_acc=a.genome_acc,
                          effector=a.effector, cfg=cfg, fold=cfg.allow_folds,
                          local_genome=a.local_genome, local_gff=a.local_gff, local_contig=a.local_contig,
                          tf_locus=(int(a.tf_locus[0]), int(a.tf_locus[1]), a.tf_locus[2]) if a.tf_locus else None)
    except (ValueError, OSError) as exc:
        print(f"INPUT ERROR: {exc}", file=sys.stderr)
        return 2
    if not res.ok:
        stage = res.raw.get("failed_stage")
        print(f"\nERROR{f' [{stage}]' if stage else ''}: {res.error}")
        print("Next: run `tfop setup` to check engines and `tfop audit` to check packaged data/paths.")
        return 1
    d = res.raw
    print("\n--- result ---")
    for k in ("family", "family_support", "genome", "tf_locus", "n_putative_sites", "motif_consensus", "motif_mean_ic",
              "autoregulatory_recovered", "nearest_bp",
              "structure_folded", "structure_plddt", "structure_qc_rmsd",
              "inducer", "ssn_cluster", "designed_operator"):
        if d.get(k) is not None:
            print(f"  {k:24}: {d[k]}")
    if res.primary_operator:
        print(f"  {'primary_operator':18}: {res.primary_operator}  [{res.primary_source}]")
    if res.bundle:
        _print_output_map(res.bundle)
    return 0


def cmd_scan(a) -> int:
    from predictor.annotate import genome_scan as GS
    from predictor import api
    fams = {f.strip() for f in a.families.split(",")} if a.families else None
    print(f"=== scan: {a.input}  [families={a.families or 'all 19'}] ===")
    cands, kind, n = GS.scan(a.input, families=fams, cpus=a.cpus)
    out = Path(a.out) if a.out else Path("results") / f"scan_{Path(a.input).stem}"
    GS.write_outputs(cands, out, kind=kind, n_proteins=n, source=a.input)
    print(f"\n[{kind}] {n} proteins -> {len(cands)} candidate TF(s)   (written to {out})")
    for fam, k in GS.family_counts(cands).most_common():
        print(f"  {str(fam):18s} {k}")
    print(f"\nreport: {out/'SCAN_REPORT.html'}   data: candidate_tfs.json/.tsv/.faa + scan_manifest.json")
    if not a.predict:
        print("Avenue A only. Re-run with --predict (and --organism) to add each metalloregulator "
              "candidate's inducer + operators.")
        return 0
    # avenue B: run the per-TF operator pipeline on the candidates whose family it supports.
    #   --select a,b   -> only these protein_ids (individual selection)
    #   --families F   -> already narrowed the census (so "by family" too)
    todo = [c for c in cands if c.family in GS.PIPELINE_FAMILIES]
    if a.select:
        want = {s.strip() for s in a.select.split(",") if s.strip()}
        # NARROW the supported list; never rebuild from `cands`. The census recognises more families
        # than the operator pipeline supports, so selecting out of `cands` would hand an unsupported
        # family to `api.predict` -- which fails per candidate, late, after the batch has started.
        todo = [c for c in todo if c.protein_id in want]
        missing = want - {c.protein_id for c in todo}
        if missing:
            known = {c.protein_id: c.family for c in cands}
            for pid in sorted(missing):
                why = (f"family {known[pid]} is not one the operator pipeline supports"
                       if pid in known else "not found in this scan")
                print(f"  [skip] {pid}: {why}")
    if a.limit:
        todo = todo[:a.limit]
    print(f"\n--- avenue B: predicting operators for {len(todo)} candidate(s) ---")
    rows, preds = [], {}
    for i, c in enumerate(todo, 1):
        print(f"[{i}/{len(todo)}] {c.protein_id} ({c.family}) ...", flush=True)
        try:
            res = api.predict(c.sequence, name=c.protein_id, family=c.family, organism=a.organism,
                              genome_acc=a.genome_acc, fold=False, verbose=False)
            rows.append([c.protein_id, c.family, res.inducer, res.primary_operator,
                         "ok" if res.ok else (res.error or "err")])
            preds[c.protein_id] = {"inducer": res.inducer, "primary_operator": res.primary_operator,
                                   "bundle": str(res.bundle) if res.bundle else None}
        except Exception as e:                                  # one bad candidate must not sink the batch
            rows.append([c.protein_id, c.family, None, None, f"{type(e).__name__}"])
    with (out / "candidate_predictions.tsv").open("w", encoding="utf-8") as fh:
        fh.write("protein_id\tfamily\tinducer\tprimary_operator\tstatus\n")
        for r in rows:
            fh.write("\t".join(str(x if x is not None else "") for x in r) + "\n")
    # regenerate the JSON + HTML report + manifest with the predictions merged in
    GS.write_outputs(cands, out, kind=kind, n_proteins=n, source=a.input, predictions=preds)
    print(f"\n{'candidate':26s} {'family':12s} {'inducer':10s} primary operator")
    for pid, fam, ind, op, st in rows:
        print(f"  {pid[:24]:26s} {str(fam):12s} {str(ind):10s} {op or ('-' if st == 'ok' else st)}")
    print(f"\nreport: {out/'SCAN_REPORT.html'}   predictions: {out/'candidate_predictions.tsv'}  "
          f"(full per-candidate bundles under results/jobs/)")
    return 0


def cmd_audit(_a) -> int:
    from predictor import resources
    resources._print_audit()
    return 0


def cmd_genome(a) -> int:
    """Whole genome -> candidate TFs (BITACORA) -> one bundle per candidate -> summary.tsv."""
    from predictor import genome_run
    try:
        out = genome_run.run_genome(a.accession, Path(a.out), jobs=a.jobs, census_fallback=a.census_fallback,
                                    offline=a.offline, fold=a.fold, timeout=a.timeout,
                                    discover_only=a.discover_only, verbose=a.verbose, fasta=a.fasta, gff=a.gff,
                                    genbank=a.genbank, proteome=a.proteome, organism=a.organism,
                                    metalnet=not a.no_metalnet, af3=not a.no_af3,
                                    selected_names=a.select.split(",") if a.select else None)
    except (RuntimeError, ValueError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(f"\n-> {out}")
    return _run_exit_status(out)


def cmd_batch(a) -> int:
    from predictor.api import predict_many
    try:
        path = predict_many(a.input, a.out, organism=a.organism, jobs=a.jobs, offline=a.offline,
                           fold=a.fold, timeout=a.timeout, discover_only=a.discover_only,
                           metalnet=not a.no_metalnet, af3=not a.no_af3,
                           genome_fasta=a.genome_fasta, genome_gff=a.genome_gff, genome_genbank=a.genome_genbank)
        print(path)
    except (ValueError, RuntimeError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return _run_exit_status(path)


def _run_exit_status(path) -> int:
    """An emitted summary remains inspectable even when a candidate failed."""
    import csv
    path = Path(path)
    if path.name != "summary.tsv":
        return 0
    with path.open(encoding="utf-8") as handle:
        return int(any(r.get("status", "").startswith(("failed", "timed out"))
                       for r in csv.DictReader(handle, delimiter="\t")))


def cmd_report(a) -> int:
    from predictor.report.run_report import write_discovery_report, write_run_report
    try:
        path = Path(a.run_dir)
        if not path.is_dir() or not any((path / marker).exists() for marker in ("manifest.csv", "summary.tsv", "jobs")):
            raise ValueError("report requires an existing run directory with manifest, summary, or jobs")
        write_discovery_report(path)
        print(write_run_report(path))
    except (ValueError, RuntimeError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 0


def cmd_af3(a) -> int:
    from predictor.structure import af3_jobset
    try:
        path = Path(a.path)
        if not path.is_dir():
            raise ValueError("AF3 export requires an existing run or bundle directory")
        export = af3_jobset.export_bundle if (path / "dossier.json").exists() else af3_jobset.export_run
        print(export(path, top_n=a.top_n, format=a.format, batch_size=a.batch_size))
    except (ValueError, RuntimeError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 0


def cmd_params(a) -> int:
    """Every run parameter, read live from the code that uses it (predictor.params)."""
    from predictor import params
    if a.markdown:
        print(params.markdown())
        return 0
    stage = None
    for r in params.table():
        if r["stage"] != stage:
            stage = r["stage"]
            print(f"\n[{stage}]")
        print(f"  {r['name']:<34} {str(r['value']):<32} {r['where']}")
    return 0


def cmd_selftest(a) -> int:
    """Verify this install by running every module's built-in self-test in a subprocess.

    This is the user-facing counterpart to the pytest suite (which needs the dev extra and a checkout):
    after `pip install`, `tfop selftest` answers "is my install actually working?". Modules are run
    isolated so one crash cannot take the run down, and a module that needs an absent external engine or
    the network is reported as SKIP rather than failing the install check.
    """
    import subprocess
    import time
    from pathlib import Path

    root = Path(__file__).resolve().parent
    # Two conventions coexist: most modules expose a `--self-test` flag, a few just run their `_demo()`
    # on a bare `python module.py`. Detect either, and fall back to a flagless invocation below.
    # `run_cli` and `gui` are entry points, not testable modules -- and run_cli matches the discovery
    # heuristic only because this very function mentions "--self-test".
    excluded = {"run_cli", "gui", "gui_app", "setup_tools"}
    mods = sorted(
        "predictor." + ".".join(p.relative_to(root).with_suffix("").parts)
        for p in root.rglob("*.py")
        if "__pycache__" not in p.parts and p.name != "__init__.py" and p.stem not in excluded
        and any(k in p.read_text(encoding="utf-8", errors="replace")
                for k in ("_self_test", "--self-test", "def _demo"))
    )
    if a.only:
        want = [s.strip() for s in a.only.split(",") if s.strip()]
        mods = [m for m in mods if any(w in m for w in want)]

    ok = fail = skip = timed_out = 0
    failures = []
    print(f"running {len(mods)} module self-test(s) with a {a.timeout}s timeout each\n")
    for m in mods:
        t0 = time.time()
        label = m.replace("predictor.", "")
        try:
            r = subprocess.run([sys.executable, "-m", m, "--self-test"], cwd=root.parent,
                               capture_output=True, text=True, timeout=a.timeout)
            rc, out = r.returncode, (r.stderr or r.stdout or "")
            if rc != 0 and "unrecognized arguments" in out:
                r = subprocess.run([sys.executable, "-m", m], cwd=root.parent,
                                   capture_output=True, text=True, timeout=a.timeout)
                rc, out = r.returncode, (r.stderr or r.stdout or "")
        except subprocess.TimeoutExpired:
            rc, out = None, "timed out"
        dt = time.time() - t0
        if rc == 0:
            ok += 1
            status = "ok"
        elif rc is None:
            # A module that did not finish is neither a pass nor evidence of a missing engine: report
            # it as its own status so a slow or hung test cannot hide among the environmental skips.
            timed_out += 1
            status = "timeout"
        elif _is_environmental(out):
            skip += 1
            status = "skip"
        else:
            fail += 1
            status = "FAIL"
            failures.append((label, out.strip().splitlines()[-1][:160] if out.strip() else "no output"))
        if status != "ok" or a.verbose:
            print(f"  {status:7s} {dt:5.1f}s  {label}")

    print(f"\n{ok} ok, {skip} skipped (missing engine/network), {timed_out} timed out "
          f"(>{a.timeout}s), {fail} failed")
    if timed_out:
        print(f"A timeout is not a pass: rerun those modules with a longer --timeout "
              f"(e.g. `tfop selftest --only <module> --timeout {a.timeout * 4}`).")
    for label, why in failures:
        print(f"  FAIL {label}: {why}")
    if fail:
        print("\nA failure usually means a missing reference file or a broken install. "
              "Run `tfop audit` to see how each external resource resolves.")
    return 1 if fail else 0


def _is_environmental(output: str) -> bool:
    """True when a self-test failed for a reason that is about the ENVIRONMENT, not the code: an absent
    external engine, no network, or a missing optional dependency. Those are skips, not failures."""
    markers = ("No module named", "not found", "FileNotFoundError", "ConnectionError", "URLError",
               "HTTPError", "timed out", "Temporary failure in name resolution", "mmseqs", "blast",
               "foldseek", "datasets", "BIOHUB_TOKEN", "NCBI", "CalledProcessError", "WinError 2")
    return any(m.lower() in output.lower() for m in markers)


def cmd_setup(_a) -> int:
    from predictor import setup_tools
    return setup_tools.main()


def main(argv=None) -> int:
    # Output redirected on Windows (`tfop predict ... > log.txt`) is encoded with the ANSI code page;
    # escape a character it cannot represent (e.g. in a gene product name) instead of crashing the run.
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(errors="backslashreplace")
        except (AttributeError, ValueError):
            pass
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=("First time? Run:  tfop audit  ->  tfop setup  ->  tfop selftest\n"
                "The core reference database is bundled with the package; there is no separate data install."),
    )
    from predictor import __version__
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("predict", help="TF protein -> operators, motif, inducer, regulon, and HTML report")
    p.add_argument("seq", help="TF protein sequence, one-record FASTA path, or protein accession (WP_/NP_/UniProt)")
    p.add_argument("--name", default="query", help="a name for the run (the output folder)")
    p.add_argument("--out", default=None,
                   help="output root (default: ./results; bundle is written under <out>/jobs/<name>)")
    p.add_argument("--family", default=None, help="TF family (else auto-classified)")
    p.add_argument("--organism", default=None, help="RefSeq/GCF accession of the source genome (optional; "
                                                    "it can also be found from the sequence via NCBI)")
    p.add_argument("--genome-acc", dest="genome_acc", default=None,
                   help="genome assembly accession to use; takes precedence over --organism")
    p.add_argument("--effector", default=None, help="known inducer ion/ligand (else inferred)")
    p.add_argument("--fold", action="store_true",
                   help="add the optional apo-dimer ESMFold2 fold (OFF by default; needs BIOHUB_TOKEN)")
    p.add_argument("--no-fold", dest="no_fold", action="store_true",
                   help=argparse.SUPPRESS)  # retained for backwards-compatible scripts; OFF is the default
    # COST levers. Each turns an expensive optional step off; none selects a different algorithm --
    # the path is identical and the switched-off source abstains in the open, exactly as it does on a
    # machine where that engine is not installed.
    p.add_argument("--no-metalnet", dest="no_metalnet", action="store_true",
                   help="skip the MetalNet2 metal-site gate (~2 min/sequence); it abstains instead")
    p.add_argument("--no-af3", dest="no_af3", action="store_true",
                   help="do not write the optional AlphaFold-3 hand-off JSONs")
    p.add_argument("--offline", action="store_true",
                   help="no network: no NCBI, no ColabFold, no AFDB (uses the local genome mirror)")
    p.add_argument("--local-genome", help=argparse.SUPPRESS)
    p.add_argument("--local-gff", help=argparse.SUPPRESS)
    p.add_argument("--local-contig", help=argparse.SUPPRESS)
    p.add_argument("--tf-locus", nargs=3, help=argparse.SUPPRESS)
    p.set_defaults(func=cmd_predict)

    psc = sub.add_parser("scan", help="quick Pfam-HMM census of a genome/proteome; --predict adds operators "
                                      "(for production discovery use `tfop genome`)")
    psc.add_argument("input", help="a proteome FASTA (protein) OR a genome FASTA (nucleotide -> gene-called)")
    psc.add_argument("--families", default=None,
                     help="comma list of families to keep, e.g. 'ArsR/SmtB,MerR,Fur' (default: all 19)")
    psc.add_argument("--out", default=None, help="output dir (default: results/scan_<input>)")
    psc.add_argument("--cpus", type=int, default=0, help="hmmsearch threads (0 = all cores)")
    psc.add_argument("--predict", action="store_true",
                     help="avenue B: run the operator pipeline on each candidate (heavier)")
    psc.add_argument("--select", default=None,
                     help="with --predict, run only these candidate protein_ids (comma list) -- individual "
                          "selection; combine with --families for by-family")
    psc.add_argument("--organism", default=None,
                     help="genome accession for --predict (homolog + genome resolution)")
    psc.add_argument("--genome-acc", dest="genome_acc", default=None,
                     help="with --predict, genome assembly accession; takes precedence over --organism")
    psc.add_argument("--limit", type=int, default=0, help="with --predict, cap the number of candidates")
    psc.set_defaults(func=cmd_scan)

    pa = sub.add_parser("audit", help="show bundled database release/size plus cache and output paths")
    pa.set_defaults(func=cmd_audit)

    ps = sub.add_parser("setup", help="check required engines and explain how to install anything missing")
    ps.set_defaults(func=cmd_setup)

    pp = sub.add_parser("params", help="print every run parameter, read live from the code that uses it")
    pp.add_argument("--markdown", action="store_true", help="emit the Methods-ready markdown table")
    pp.set_defaults(func=cmd_params)

    pg = sub.add_parser("genome", help="whole genome -> candidate TFs (BITACORA) -> a bundle per candidate "
                                       "-> summary.tsv (the production route)")
    pg.add_argument("accession", nargs="?", help="assembly accession, e.g. GCF_000005845.2 (mirrored or downloaded)")
    pg.add_argument("--out", required=True, help="run directory (manifest, FASTAs, jobs/, logs/, summary.tsv)")
    pg.add_argument("--jobs", type=int, default=1, help="candidates predicted in parallel (default 1)")
    pg.add_argument("--timeout", type=int, default=None, help="per-candidate timeout in seconds")
    pg.add_argument("--census-fallback", action="store_true",
                    help="use the Pfam census if BITACORA is unavailable (recorded as a different engine)")
    pg.add_argument("--discover-only", action="store_true", help="stop after writing manifest.csv")
    pg.add_argument("--offline", action="store_true", help="no network (genome must already be mirrored)")
    pg.add_argument("--fold", action="store_true", help="also fold each candidate (off by default)")
    pg.add_argument("-v", "--verbose", action="store_true",
                    help="print candidate-discovery (BITACORA / gene-calling) progress details")
    inputs = pg.add_mutually_exclusive_group()
    inputs.add_argument("--fasta", help="local nucleotide genome FASTA")
    inputs.add_argument("--genbank", help="local GenBank genome (one or more records)")
    inputs.add_argument("--proteome", help="local protein FASTA for discovery")
    pg.add_argument("--gff", help="GFF3 annotation for --fasta")
    pg.add_argument("--organism", help="optional source genome accession")
    pg.add_argument("--select", help="comma-separated manifest run_names to predict (all candidates remain in reports)")
    pg.add_argument("--no-metalnet", action="store_true",
                    help="skip the optional MetalNet2 metal-site gate; that evidence source abstains")
    pg.add_argument("--no-af3", action="store_true", help="do not write the optional AlphaFold-3 hand-off JSONs")
    pg.set_defaults(func=cmd_genome)

    pb = sub.add_parser("batch", help="predict every protein in multi-FASTA or accession TSV")
    pb.add_argument("input", help="multi-record protein FASTA, or a TSV whose header has an 'accession' column")
    pb.add_argument("--out", required=True, help="run directory (manifest, FASTAs, jobs/, logs/, summary.tsv)")
    pb.add_argument("--organism", help="source genome accession (e.g. GCF_000005845.2) for genomic context")
    context_inputs = pb.add_mutually_exclusive_group()
    context_inputs.add_argument("--genome-fasta", help="local source genome for supplied proteins")
    context_inputs.add_argument("--genome-genbank", help="local GenBank source genome")
    pb.add_argument("--genome-gff", help="annotation for --genome-fasta")
    pb.add_argument("--jobs", type=int, default=1, help="proteins predicted in parallel (default 1)")
    pb.add_argument("--timeout", type=int, help="per-protein timeout in seconds")
    pb.add_argument("--offline", action="store_true", help="no network; accession lookups are refused")
    pb.add_argument("--fold", action="store_true", help="also fold each protein (off by default)")
    pb.add_argument("--discover-only", action="store_true", help="stop after writing manifest.csv")
    pb.add_argument("--no-metalnet", action="store_true",
                    help="skip the optional MetalNet2 metal-site gate; that evidence source abstains")
    pb.add_argument("--no-af3", action="store_true", help="do not write the optional AlphaFold-3 hand-off JSONs")
    pb.set_defaults(func=cmd_batch)

    pr = sub.add_parser("report", help="rebuild run reports without re-predicting")
    pr.add_argument("run_dir", help="an existing run directory (holds manifest.csv, summary.tsv or jobs/)")
    pr.set_defaults(func=cmd_report)

    paf = sub.add_parser("af3", help="rebuild AF Server/local jobs outside existing bundles")
    paf.add_argument("path", help="a run directory (exports every jobs/* bundle) or one TF bundle")
    paf.add_argument("--top-n", type=int, default=3, help="DNA jobs for up to this many top operators (default 3)")
    paf.add_argument("--format", choices=["server", "local", "both"], default="both",
                     help="AF Server JSON, local AlphaFold 3 JSON, or both (default both)")
    paf.add_argument("--batch-size", type=int, default=30, help="jobs per AF Server batch file (default 30)")
    paf.set_defaults(func=cmd_af3)

    pst = sub.add_parser("selftest", help="verify this install by running every module's built-in self-test")
    pst.add_argument("--only", default=None,
                     help="comma-separated substrings; run only matching modules (e.g. 'metal_site,inducer')")
    # 300 s, not 120: signals.operator_logo's EM alignment alone takes ~80 s on an idle desktop, so a
    # slower or busy computer timed out a working module (seen in fresh Windows installs).
    pst.add_argument("--timeout", type=int, default=300, help="per-module timeout in seconds (default 300)")
    pst.add_argument("-v", "--verbose", action="store_true", help="list passing modules too")
    pst.set_defaults(func=cmd_selftest)

    a = ap.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    raise SystemExit(main())
