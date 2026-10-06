"""genome_run.py -- a whole genome in, a directory of prediction bundles and one summary table out.

The genome-first path, packaged so the CLI (`tfop genome`), the GUI and the analysis drivers all run the
same code:

  1. `ensure_genome`       the genome and its annotation, from the local mirror or NCBI Datasets
  2. `discover`            candidate transcription factors: BITACORA with the SSN-clade seed profiles
                           (the production route). The Pfam census is used only when asked for, and is
                           then recorded as a different engine on every row.
  3. `write_manifest`      one row per candidate, full-length FASTA per candidate, same schema as the
                           analysis run manifests
  4. `predict_all`         one `tfop predict` per candidate, each in its own process (a crash or timeout
                           cannot take the run down); resumes past bundles that are already complete
  5. `summarise`           summary.tsv: one row per candidate with the headline calls

Nothing here changes a prediction: step 4 is exactly `tfop predict` per candidate.
"""
from __future__ import annotations

import csv
import json
import os
import re
import subprocess
import sys
from pathlib import Path

MANIFEST_COLS = ["run_name", "fasta", "genome", "organism_acc", "family", "family_flag", "validated",
                 "contig", "start", "end", "strand", "aa_len", "locus_tag", "gene", "product"]
#: Provenance. Kept apart from the columns above: a table that cannot tell a BITACORA row from a
#: census row will eventually be asked to compare the two.
PROV_COLS = ["discovery_engine", "discovery_profile", "protein_source", "trim_start", "trim_end"]


def safe(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s or "").strip("_") or "x"


def _unescape(s: str) -> str:
    return re.sub(r"%([0-9A-Fa-f]{2})", lambda m: chr(int(m.group(1), 16)), s).replace(",", " ")


def ensure_genome(accession: str, *, allow_download: bool = True) -> tuple[Path, Path | None]:
    from predictor.annotate import genome_mirror
    fna, gff = genome_mirror.ensure_genome(accession, allow_download=allow_download)
    if fna is None:
        raise RuntimeError(f"genome {accession} is not in the local mirror and could not be downloaded "
                           f"(needs the NCBI Datasets CLI; see `tfop setup`)")
    return Path(fna), (Path(gff) if gff else None)


def gff_index(gff: Path | None) -> dict:
    """{protein_id: {contig, start, end, strand, locus_tag, gene, product}} from the genome's own GFF.

    BITACORA runs in protein mode and reports no coordinates; they come from the same annotation the
    proteome was read from, which keeps start/end on the genome's own frame.
    """
    out: dict[str, dict] = {}
    if not gff or not Path(gff).exists():
        return out
    for line in Path(gff).read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        pid = attrs.get("protein_id", "")
        if not pid:
            continue
        rec = out.setdefault(pid, {"contig": f[0], "start": int(f[3]), "end": int(f[4]), "strand": f[6],
                                   "locus_tag": attrs.get("locus_tag", ""), "gene": attrs.get("gene", ""),
                                   "product": _unescape(attrs.get("product", ""))})
        rec["start"], rec["end"] = min(rec["start"], int(f[3])), max(rec["end"], int(f[4]))  # split CDS
    return out


def discover(fna: Path, *, gff: Path | None = None, work_dir: Path | None = None,
             census_fallback: bool = False, verbose: bool = False) -> tuple[list, str, str]:
    """(candidates, engine, protein_source) for one genome.

    BITACORA is the production engine. Without it this raises, unless `census_fallback` is set; the
    census reaches nearly the same set but is a different method, so it is never substituted silently.
    `gff` is only needed when the annotation is not a sibling `.gff` of `fna` (an upload). `work_dir`
    must be given for any input whose file name is not unique, because BITACORA's default work
    directory is keyed on it and caches the proteome there.
    """
    from predictor.discovery import bitacora

    ok, why = bitacora.engine_available()
    if ok:
        res = bitacora.find_candidates(fna, gff, out_dir=work_dir, verbose=verbose)
        if res.status == "ok":
            return list(res.candidates), "bitacora", res.protein_source
        if not census_fallback:
            raise RuntimeError(f"BITACORA failed on {fna.name}: {res.reason[:400]} -- pass "
                               f"census_fallback to use the Pfam census, recorded as a different engine")
    elif not census_fallback:
        raise RuntimeError(f"the BITACORA engine is not available ({why}); install it (docs/ENGINES.md) or "
                           f"allow the Pfam census fallback, which is recorded as a different engine")
    from predictor.annotate import genome_scan
    recs, kind = genome_scan.proteins_from_input(str(fna), gff=gff, verbose=verbose)
    cands = [c for c in genome_scan.scan_proteins(recs) if c.family in genome_scan.PIPELINE_FAMILIES]
    return cands, "pfam_census", kind


def as_census_candidates(cands: list, fna: Path, gff: Path | None = None) -> tuple[list, int]:
    """Production candidates as `genome_scan.Candidate` rows, located on the genome, for the GUI table.

    BITACORA runs in protein mode and reports no coordinates; they are read back from the same
    proteome reader it searched (the genome's own GFF, or called genes when there is none). Returns
    (rows, number of proteins searched). The Pfam columns stay empty: these were not found by Pfam.
    """
    from predictor.annotate import genome_scan as GS

    prots, _kind, meta = GS.proteins_from_input(str(fna), gff=gff, verbose=False, return_metadata=True)
    rows = []
    for c in cands:
        ev, m = getattr(c, "evidence", {}) or {}, meta.get(c.protein_id, {})

        def _num(k, ev=ev):
            try:
                return float(ev[k])
            except (KeyError, TypeError, ValueError):
                return float("nan")
        rows.append(GS.Candidate(
            protein_id=c.protein_id, family=c.family or None, pfam_acc=None, pfam_name=None,
            evalue=_num("evalue"), bits=_num("bits"), length=len(c.sequence), sequence=c.sequence,
            contig=m.get("contig"), start=m.get("start"), end=m.get("end"), strand=m.get("strand"),
            locus_tag=m.get("locus_tag"), gene=m.get("gene"), product=m.get("product"),
            coordinate_source=m.get("coordinate_source")))
    return rows, len(prots)


def manifest_rows(cands: list, *, short: str, organism: str, accession: str, engine: str,
                  protein_source: str, coords: dict, fasta_dir: Path, rel_to: Path | None = None) -> list[dict]:
    """Write one full-length FASTA per familied candidate and return its manifest rows."""
    fasta_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for c in cands:
        fam = c.family or ""
        if not fam:
            continue                           # unfamilied hits are census noise, not candidates
        pid, seq = c.protein_id, c.sequence
        fasta = fasta_dir / f"{safe(fam)}__{safe(pid)}.fasta"
        fasta.write_text(f">{pid} {fam}\n{seq}\n", encoding="utf-8", newline="\n")
        g = coords.get(pid, {})
        ts, te = getattr(c, "trim_start", None), getattr(c, "trim_end", None)
        rows.append({
            "run_name": f"{short}__{safe(fam)}__{safe(pid)}",
            "fasta": (fasta.relative_to(rel_to) if rel_to else fasta).as_posix(),
            "genome": organism, "organism_acc": accession, "family": fam,
            # `(auto)`: the family above is discovery's routing hint; the pipeline's own family call
            # must stay free to disagree with it.
            "family_flag": "(auto)", "validated": "no",
            "contig": g.get("contig", ""), "start": g.get("start", ""), "end": g.get("end", ""),
            "strand": g.get("strand", ""), "aa_len": len(seq), "locus_tag": g.get("locus_tag", ""),
            "gene": g.get("gene", ""), "product": g.get("product", ""),
            "discovery_engine": engine, "discovery_profile": getattr(c, "profile", ""),
            "protein_source": protein_source,
            "trim_start": "" if ts is None else ts, "trim_end": "" if te is None else te,
        })
    return rows


def write_manifest(rows: list[dict], path: Path) -> None:
    with Path(path).open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=MANIFEST_COLS + PROV_COLS)
        w.writeheader()
        w.writerows(rows)


def bundle_complete(bundle: Path) -> bool:
    try:
        status = json.loads((bundle / "bundle_manifest.json").read_text(encoding="utf-8")).get("status")
    except (OSError, ValueError):
        return False
    return status == "complete" and (bundle / "dossier.json").is_file()


def predict_one(row: dict, *, out_dir: Path, offline: bool = False, fold: bool = False,
                timeout: int | None = None) -> tuple[str, str]:
    """`tfop predict` for one manifest row, in its own process. Returns (run_name, status)."""
    name = row["run_name"]
    if bundle_complete(out_dir / "jobs" / name):
        return name, "done (resumed)"
    cmd = [sys.executable, "-m", "predictor.run_cli", "predict", str(out_dir / row["fasta"]),
           "--name", name, "--organism", row["organism_acc"]]
    if offline:
        cmd.append("--offline")
    if fold:
        cmd.append("--fold")
    env = dict(os.environ, PREDICTOR_OUTPUT_DIR=str(out_dir))
    (out_dir / "logs").mkdir(parents=True, exist_ok=True)
    with (out_dir / "logs" / f"{name}.log").open("w", encoding="utf-8") as log:
        try:
            rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, env=env, timeout=timeout).returncode
        except subprocess.TimeoutExpired:
            return name, "timed out"
    return name, ("ok" if rc == 0 and bundle_complete(out_dir / "jobs" / name) else f"failed (exit {rc})")


def predict_all(rows: list[dict], *, out_dir: Path, jobs: int = 1, **kw) -> dict[str, str]:
    from concurrent.futures import ThreadPoolExecutor

    from concurrent.futures import as_completed

    status: dict[str, str] = {}
    with ThreadPoolExecutor(max_workers=max(1, jobs)) as pool:
        # Report in completion order: `pool.map` yields in submission order, so one slow first
        # candidate hid every later completion behind it.
        futures = [pool.submit(predict_one, r, out_dir=out_dir, **kw) for r in rows]
        for i, fut in enumerate(as_completed(futures), 1):
            name, st = fut.result()
            status[name] = st
            print(f"  [{i}/{len(rows)}] {name}: {st}", flush=True)
    # Callers read this in manifest order; completion order is only for the progress lines.
    return {r["run_name"]: status[r["run_name"]] for r in rows if r["run_name"] in status}


SUMMARY_COLS = ["run_name", "locus_tag", "gene", "product", "discovery_family", "family", "ssn_cluster",
                "inducer_top", "inducer_class", "agreement", "primary_operator", "primary_start",
                "primary_end", "primary_locality", "n_regulon_operons", "status"]


def summarise(rows: list[dict], status: dict[str, str], *, out_dir: Path) -> Path:
    """summary.tsv: one row per candidate with the headline calls, read back from each bundle."""
    from predictor.effector import inducer_vocab as V

    out = []
    for r in rows:
        name = r["run_name"]
        rec = {"run_name": name, "locus_tag": r["locus_tag"], "gene": r["gene"], "product": r["product"],
               "discovery_family": r["family"], "status": status.get(name, "not run")}
        bundle = out_dir / "jobs" / name
        if bundle_complete(bundle):
            d = json.loads((bundle / "dossier.json").read_text(encoding="utf-8"))
            ind = d.get("inducers") or {}
            prim = (d.get("operators") or {}).get("primary") or {}
            pv = prim.get("provenance") or {}
            rec.update({"family": d.get("family"), "ssn_cluster": d.get("ssn_cluster"),
                        "inducer_top": ind.get("top"), "inducer_class": V.inducer_class_of(ind.get("top")),
                        "agreement": ind.get("agreement"), "primary_operator": prim.get("sequence"),
                        "primary_start": pv.get("start"), "primary_end": pv.get("end"),
                        "primary_locality": pv.get("locality"),
                        "n_regulon_operons": len(d.get("regulon") or [])})
        out.append(rec)
    path = out_dir / "summary.tsv"
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=SUMMARY_COLS, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(out)
    return path


def run_genome(accession: str, out_dir: Path, *, jobs: int = 1, census_fallback: bool = False,
               offline: bool = False, fold: bool = False, timeout: int | None = None,
               discover_only: bool = False, verbose: bool = False) -> Path:
    """The whole genome-first run. Returns the manifest path (discover_only) or summary.tsv."""
    # The NCBI Datasets CLI path ($DATASETS_EXE) and the API tokens live in the predictor's .env. `tfop
    # predict` loads them, but on a cold genome cache this function needs them first, to download.
    from predictor import resources
    resources.load_tokens()
    out_dir = Path(out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = out_dir / "manifest.csv"
    if manifest.exists():
        with manifest.open(encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        # Resuming a directory written for ANOTHER genome would silently predict that genome's
        # candidates and report them under this accession.
        recorded = {r.get("organism_acc") for r in rows} - {None, ""}
        if recorded and recorded != {accession}:
            raise RuntimeError(f"{out_dir} holds a run for {', '.join(sorted(recorded))}, not {accession}; "
                               "use a different --out, or remove that directory to start over")
        print(f"resuming: {len(rows)} candidate(s) from {manifest}")
    else:
        fna, gff = ensure_genome(accession, allow_download=not offline)
        cands, engine, psource = discover(fna, census_fallback=census_fallback, verbose=verbose)
        rows = manifest_rows(cands, short=safe(accession), organism=accession, accession=accession,
                             engine=engine, protein_source=psource, coords=gff_index(gff),
                             fasta_dir=out_dir / "fastas", rel_to=out_dir)
        if not rows:
            raise RuntimeError(f"no candidate transcription factors found in {accession}")
        write_manifest(rows, manifest)
        print(f"{len(rows)} candidate(s) found with {engine} ({psource}) -> {manifest}")
    if discover_only:
        return manifest
    status = predict_all(rows, out_dir=out_dir, jobs=jobs, offline=offline, fold=fold, timeout=timeout)
    return summarise(rows, status, out_dir=out_dir)
