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
import hashlib
import shutil
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


def _candidate_name(name: str) -> None:
    from predictor.report.bundle_contract import _portable_parts
    if name != safe(name) or len(_portable_parts(name)) != 1:
        raise ValueError("candidate names must be safe path components")
    # Leave room for the FASTA/log suffix on filesystems with 255-byte components.
    if len(name.encode("utf-8")) > 240:
        raise ValueError("candidate name is too long for portable filenames (maximum 240 bytes)")


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
             census_fallback: bool = False, verbose: bool = False, census_only: bool = False) -> tuple[list, str, str]:
    """(candidates, engine, protein_source) for one genome.

    BITACORA is the production engine. Without it this raises, unless `census_fallback` is set; the
    census reaches nearly the same set but is a different method, so it is never substituted silently.
    `gff` is only needed when the annotation is not a sibling `.gff` of `fna` (an upload). `work_dir`
    must be given for any input whose file name is not unique, because BITACORA's default work
    directory is keyed on it and caches the proteome there.
    """
    from predictor.discovery import bitacora

    ok, why = (False, "census explicitly selected") if census_only else bitacora.engine_available()
    if ok:
        res = bitacora.find_candidates(fna, gff, out_dir=work_dir, verbose=verbose)
        if res.status == "ok":
            return list(res.candidates), "bitacora", res.protein_source
        if not census_fallback:
            raise RuntimeError(f"BITACORA failed on {fna.name}: {res.reason[:400]} -- pass "
                               f"census_fallback to use the Pfam census, recorded as a different engine")
    elif not census_fallback and not census_only:
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
    tokens = [f"{safe(c.family or '')}__{safe(c.protein_id)}" for c in cands
              if c.family or engine == "provided_protein"]
    for token in tokens:
        _candidate_name(token)
        _candidate_name(f"{short}__{token}")
    if len(tokens) != len({token.casefold() for token in tokens}):
        raise ValueError("duplicate or colliding protein identifier")
    fasta_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    used = set()
    for c in cands:
        fam = c.family or ""
        if not fam and engine != "provided_protein":
            continue                           # unfamilied hits are census noise, not candidates
        pid, seq = c.protein_id, c.sequence
        token = f"{safe(fam)}__{safe(pid)}"
        _candidate_name(token)
        if token.casefold() in used:
            raise ValueError(f"duplicate or colliding protein identifier: {pid!r}")
        used.add(token.casefold())
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
                timeout: int | None = None, metalnet: bool = True, af3: bool = True) -> tuple[str, str]:
    """`tfop predict` for one manifest row, in its own process. Returns (run_name, status)."""
    name = row["run_name"]
    if bundle_complete(out_dir / "jobs" / name):
        record_path = out_dir / "jobs" / name / "run" / "run_record.json"
        if record_path.is_file():
            record = json.loads(record_path.read_text(encoding="utf-8"))
            requested = {"allow_online": not offline, "allow_folds": bool(fold),
                         "allow_metalnet": bool(metalnet), "emit_af3_jobs": bool(af3)}
            if any(record.get("effective_config", {}).get(key) != value for key, value in requested.items()):
                return name, "failed (resume configuration changed; choose a new output directory)"
            from predictor.input_validation import read_protein
            sequence = read_protein(str(out_dir / row["fasta"]))
            if record.get("query", {}).get("sha256") != hashlib.sha256(sequence.encode()).hexdigest():
                return name, "failed (resume protein input changed; choose a new output directory)"
        return name, "done (resumed)"
    cmd = [sys.executable, "-m", "predictor.run_cli", "predict", str(out_dir / row["fasta"]),
           "--name", name]
    if row.get("organism_acc"):
        cmd += ["--organism", row["organism_acc"]]
    local = out_dir / "inputs" / "genome.fna"
    annotation = out_dir / "inputs" / "genome.gff"
    if local.is_file() and annotation.is_file():
        cmd += ["--local-genome", str(local), "--local-gff", str(annotation)]
        if row.get("contig"):
            cmd += ["--local-contig", row["contig"]]
        if row.get("start") and row.get("end"):
            cmd += ["--tf-locus", str(int(row["start"]) - 1), str(row["end"]), row.get("strand") or "+"]
    if not metalnet:
        cmd.append("--no-metalnet")
    if not af3:
        cmd.append("--no-af3")
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
        futures = {pool.submit(predict_one, r, out_dir=out_dir, **kw): r["run_name"] for r in rows}
        for i, fut in enumerate(as_completed(futures), 1):
            try:
                name, st = fut.result()
            except Exception as exc:
                name, st = futures[fut], f"failed ({type(exc).__name__}: {exc})"
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


def input_identity(mode: str, *, files=(), accession=None, organism=None) -> dict:
    """Portable identity that prevents resuming unrelated input in the same output directory."""
    return {"mode": mode, "accession": accession, "organism": organism,
            "files": [{"role": role, "sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest()}
                      for role, path in files if path]}


def check_run_identity(out_dir: Path, identity: dict) -> None:
    path = out_dir / "run_input.json"
    if path.exists() and json.loads(path.read_text(encoding="utf-8")) != identity:
        raise ValueError(f"{out_dir} belongs to different input; choose a new output directory")
    if not path.exists() and (out_dir / "manifest.csv").exists() and identity["mode"] != "genome_accession":
        raise ValueError("existing legacy run has no input fingerprint; rebuild reports or use a new directory")
    path.write_text(json.dumps(identity, indent=2) + "\n", encoding="utf-8")


def finish_run(out_dir: Path, *, af3: bool = True) -> None:
    from predictor.report.run_report import write_discovery_report, write_run_report
    write_discovery_report(out_dir)
    if af3:
        from predictor.structure import af3_jobset
        try:
            af3_jobset.export_run(out_dir)
            (out_dir / "af3" / "export_warnings.json").unlink(missing_ok=True)
        except ValueError as exc:
            directory = out_dir / "af3"
            directory.mkdir(exist_ok=True)
            (directory / "export_warnings.json").write_text(json.dumps([str(exc)], indent=2), encoding="utf-8")
    write_run_report(out_dir)


def run_candidates(rows: list[dict], out_dir: Path, *, jobs: int = 1, discover_only: bool = False,
                   offline: bool = False, fold: bool = False, timeout=None, metalnet=True, af3=True,
                   selected_names=None) -> Path:
    """Shared resumable per-TF runner used by every multi-protein entry surface."""
    out_dir = Path(out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    if any(not isinstance(row, dict) or not row.get("run_name") or not row.get("fasta") for row in rows):
        raise ValueError("candidate manifest requires nonempty run_name and fasta fields on every row")
    names = [r["run_name"] for r in rows]
    if len(names) != len({n.casefold() for n in names}):
        raise ValueError("candidate names must be unique safe path components")
    from predictor.report.bundle_contract import contained_file
    for name in names:
        _candidate_name(name)
        contained_file(out_dir, f"jobs/{name}")
        contained_file(out_dir, f"logs/{name}.log")
    for row in rows:
        contained_file(out_dir, row["fasta"])
    if selected_names is not None and set(selected_names) - set(names):
        raise ValueError("selection contains unknown candidate names")
    if discover_only:
        write_manifest(rows, out_dir / "manifest.csv")
        from predictor.report.run_report import write_discovery_report, write_run_report
        write_discovery_report(out_dir)
        write_run_report(out_dir)
        return out_dir / "manifest.csv"
    chosen = rows if selected_names is None else [r for r in rows if r["run_name"] in selected_names]
    identity = {"candidates": [{"row": {key: "" if value is None else str(value) for key, value in row.items()},
                                "sha256": hashlib.sha256(contained_file(out_dir, row["fasta"]).read_bytes()).hexdigest()}
                               for row in rows],
                "local_context": [{"file": name, "sha256": hashlib.sha256(contained_file(out_dir, "inputs/" + name).read_bytes()).hexdigest()}
                                  for name in ("genome.fna", "genome.gff") if (out_dir / "inputs" / name).is_file()],
                "policy": {"offline": bool(offline), "fold": bool(fold), "metalnet": bool(metalnet), "af3": bool(af3)}}
    identity_path = out_dir / "prediction_input.json"
    if identity_path.is_file() and json.loads(identity_path.read_text(encoding="utf-8")) != identity:
        raise ValueError("prediction inputs or configuration changed; choose a new output directory")
    identity_path.write_text(json.dumps(identity, indent=2) + "\n", encoding="utf-8")
    write_manifest(rows, out_dir / "manifest.csv")
    status = predict_all(chosen, out_dir=out_dir, jobs=jobs, offline=offline, fold=fold,
                         timeout=timeout, metalnet=metalnet, af3=af3)
    summary = summarise(rows, status, out_dir=out_dir)
    finish_run(out_dir, af3=af3)
    return summary


def prepare_local(fasta, gff, genbank, proteome, out_dir: Path):
    """Copy local inputs and preserve supplied or called CDS annotations for the prediction core."""
    from predictor.annotate import genome_scan as GS
    dest = out_dir / "inputs"
    dest.mkdir(parents=True, exist_ok=True)
    if genbank:
        from Bio import SeqIO
        from predictor.annotate.context import _GENE_FEATURES
        fna, annotation = dest / "genome.fna", dest / "genome.gff"
        records = list(SeqIO.parse(str(genbank), "genbank"))
        if not records:
            raise ValueError("GenBank input contains no records")
        ids = [r.id for r in records]
        if len(ids) != len(set(ids)):
            raise ValueError("GenBank contig identifiers must be unique")
        SeqIO.write(records, str(fna), "fasta")
        lines = ["##gff-version 3"]
        from urllib.parse import quote
        translations = []
        for record in records:
            if str(record.annotations.get("topology", "")).lower() == "circular":
                lines.append("\t".join([record.id, "GenBank", "region", "1", str(len(record)), ".", "+", ".", "Is_circular=true"]))
            for k, feature in enumerate(record.features):
                if feature.type not in _GENE_FEATURES:
                    continue
                if feature.location is None:
                    raise ValueError(f"GenBank {feature.type} feature has no location")
                parts = list(getattr(feature.location, "parts", [feature.location]))
                if len(parts) > 1 and any(int(p.start) == 0 for p in parts) and any(int(p.end) == len(record) for p in parts):
                    raise ValueError("origin-wrapped GenBank feature cannot be represented safely by the linear schema")
                attrs = {key: feature.qualifiers[key][0] for key in ("protein_id", "locus_tag", "gene", "product") if feature.qualifiers.get(key)}
                offset = 0
                if feature.type == "CDS":
                    attrs.setdefault("protein_id", attrs.get("locus_tag", f"{record.id}_CDS{k}"))
                    offset = int((feature.qualifiers.get("codon_start") or ["1"])[0]) - 1
                    if offset not in (0, 1, 2):
                        raise ValueError("GenBank codon_start must be 1, 2, or 3")
                    table = int((feature.qualifiers.get("transl_table") or ["11"])[0])
                    translation = (feature.qualifiers.get("translation") or [None])[0]
                    if translation is None:
                        nucleotides = feature.extract(record.seq)[offset:]
                        nucleotides = nucleotides[:len(nucleotides) - len(nucleotides) % 3]
                        translation = str(nucleotides.translate(table=table, to_stop=True))
                    translations.append((attrs["protein_id"], translation))
                    attrs["transl_table"] = str(table)
                attrs["ID"] = f"{record.id}_{feature.type}{k}"
                consumed = 0
                for part in parts:
                    phase = str(offset if consumed == 0 else (3 - (consumed - offset) % 3) % 3) if feature.type == "CDS" else "."
                    lines.append("\t".join([record.id, "GenBank", feature.type, str(int(part.start)+1), str(int(part.end)), ".", "-" if part.strand == -1 else "+", phase, ";".join(f"{a}={quote(str(b))}" for a,b in attrs.items())]))
                    consumed += len(part)
        if len(translations) != len({pid for pid, _ in translations}):
            raise ValueError("GenBank protein identifiers must be unique")
        annotation.write_text("\n".join(lines)+"\n", encoding="utf-8")
        if translations:
            (dest / "genome.proteins.faa").write_text("".join(f">{pid}\n{seq}\n" for pid, seq in translations), encoding="utf-8")
        return fna, annotation
    source = Path(proteome or fasta)
    fna = dest / ("proteome.faa" if proteome else "genome.fna")
    shutil.copyfile(source, fna)
    annotation = None
    sibling = next((p for p in (source.with_suffix(".gff"), source.with_suffix(".gff3")) if p.exists()), None)
    supplied = gff or (sibling if not proteome else None)
    if supplied:
        shutil.copyfile(supplied, dest / "source_annotation.gff")
        annotation = dest / "genome.gff"
        shutil.copyfile(supplied, annotation)
    contigs = [pid for pid, _seq in GS.read_fasta(fna)]
    if len(contigs) != len(set(contigs)):
        raise ValueError("local FASTA identifiers must be unique")
    proteins, kind, metadata = GS.proteins_from_input(fna, gff=annotation, return_metadata=True)
    if supplied and kind == "genome":
        raise ValueError("supplied annotation has no usable CDS; original preserved in inputs/source_annotation.gff")
    identifiers = [pid for pid, _seq in proteins]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("local protein identifiers must be unique")
    if proteome and kind != "proteome":
        raise ValueError("--proteome requires protein FASTA, not nucleotide sequence")
    if kind == "genome":
        from urllib.parse import quote
        annotation = dest / "genome.gff"
        lines = ["##gff-version 3"]
        for pid, _seq in proteins:
            m = metadata[pid]
            lines.append("\t".join([m["contig"], "Pyrodigal", "CDS", str(m["start"]), str(m["end"]), ".", m["strand"], "0", f"ID={quote(pid)};protein_id={quote(pid)};locus_tag={quote(pid)}"]))
        annotation.write_text("\n".join(lines)+"\n", encoding="utf-8")
    elif kind == "proteome" and not proteome:
        target = dest / "proteome.faa"
        fna.rename(target)
        fna = target
    return fna, annotation


def run_genome(accession: str | None = None, out_dir: Path | None = None, *, jobs: int = 1, census_fallback: bool = False,
               offline: bool = False, fold: bool = False, timeout: int | None = None,
               discover_only: bool = False, verbose: bool = False, fasta=None, gff=None,
               genbank=None, proteome=None, organism=None, metalnet=True, af3=True, census_only=False,
               selected_names=None) -> Path:
    """Discover and predict an accession, local genome/GFF, GenBank, or proteome."""
    from predictor import resources
    resources.load_tokens()
    if out_dir is None:
        raise ValueError("out_dir is required")
    if sum(bool(v) for v in (accession, fasta, genbank, proteome)) != 1:
        raise ValueError("choose exactly one accession, fasta, genbank, or proteome input")
    if gff and not fasta:
        raise ValueError("gff requires a local genome FASTA")
    out_dir = Path(out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    mode = "genome_accession" if accession else "genbank" if genbank else "proteome" if proteome else "local_genome"
    files = [("fasta", fasta), ("gff", gff), ("genbank", genbank), ("proteome", proteome)]
    if fasta and not gff:
        sibling = next((p for p in (Path(fasta).with_suffix(".gff"), Path(fasta).with_suffix(".gff3")) if p.exists()), None)
        if sibling:
            files.append(("gff", sibling))
    identity = input_identity(mode, files=files, accession=accession, organism=organism)
    identity["discovery_policy"] = {"census_only": bool(census_only), "census_fallback": bool(census_fallback)}
    manifest = out_dir / "manifest.csv"
    if manifest.exists() and accession:
        with manifest.open(encoding="utf-8") as fh:
            previous = list(csv.DictReader(fh))
        recorded = {r.get("organism_acc") for r in previous} - {None, ""}
        if recorded and recorded != {accession}:
            raise RuntimeError(f"{out_dir} holds a run for {', '.join(sorted(recorded))}, not {accession}")
    check_run_identity(out_dir, identity)
    if manifest.exists():
        with manifest.open(encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
    else:
        if accession:
            fna, annotation = ensure_genome(accession, allow_download=not offline)
        else:
            fna, annotation = prepare_local(fasta, gff, genbank, proteome, out_dir)
        cands, engine, psource = discover(fna, gff=annotation, work_dir=out_dir / "discovery", census_fallback=census_fallback, verbose=verbose, census_only=census_only)
        from dataclasses import asdict, is_dataclass
        discovery_dir = out_dir / "discovery"
        discovery_dir.mkdir(parents=True, exist_ok=True)
        evidence = [asdict(c) if is_dataclass(c) else dict(vars(c)) for c in cands]
        (discovery_dir / "candidates.json").write_text(json.dumps(evidence, indent=2, default=str), encoding="utf-8")
        from predictor.annotate import genome_scan as GS
        _proteins, _kind, coords = GS.proteins_from_input(fna, gff=annotation, return_metadata=True)
        label = accession or organism or Path(genbank or proteome or fasta).stem
        rows = manifest_rows(cands, short=safe(label), organism=organism or label, accession=accession or organism or "",
                             engine=engine, protein_source=psource, coords=coords,
                             fasta_dir=out_dir / "fastas", rel_to=out_dir)
        if not rows:
            raise RuntimeError("no candidate transcription factors found")
    return run_candidates(rows, out_dir, jobs=jobs, discover_only=discover_only, offline=offline,
                          fold=fold, timeout=timeout, metalnet=metalnet, af3=af3, selected_names=selected_names)
