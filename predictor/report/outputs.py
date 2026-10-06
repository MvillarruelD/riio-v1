"""Publish recorded per-TF outputs as a portable report folder.

`write_stage_a` stages machine-readable records, optional HTML and figures, and a
checksum manifest. Successful publication establishes the file contract, not
completion of external computation or capture of every original input.
Historical design rationale is retained in the review archive and changelog.
"""
from __future__ import annotations

import copy
import csv
import hashlib
import json
import shutil
import tempfile
import uuid
from collections import defaultdict, deque
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath, PureWindowsPath

from predictor import resources
from predictor.provenance import runtime_provenance
from predictor.report.bundle_contract import inventory, report_asset_issues

JOBS = resources.output_path("jobs")


def bundle_readme() -> str:
    """Return the portable bundle guide shared by new and migrated exports."""
    return (
        "# Individual report bundle\n\n"
        "If HTML rendering was enabled, open REPORT.html in a browser. Keep this entire directory together when sharing it; "
        "its figures and downloads use relative links. No web service is needed to view it.\n\n"
        "## Recorded files\n\n"
        "dossier.json is the complete recorded result. The input/, genome/, homologs/, motif/, "
        "operators/, ligand/, regulation/, structure/ and run/ directories contain the normalized "
        "inputs, evidence, decisions and convenience exports needed to retrace the recorded run. "
        "An empty export means that the corresponding source was unavailable or abstained.\n\n"
        "bundle_manifest.json identifies the software runtime and records each artifact's size "
        "and SHA-256 digest. The manifest excludes its own digest. Preserve the original bundle "
        "before editing it. Report regeneration preserves the previous folder automatically.\n\n"
        "report_environment.json records Python, platform and installed package versions at "
        "export time. It intentionally excludes credentials, local paths and command arguments. "
        "run/run_record.json separately records the effective pipeline arguments, stage ledger, "
        "input hashes and evidence counts without credentials or machine-specific paths.\n\n"
        "Verify the folder using: python -m predictor.report.bundle_contract PATH_TO_BUNDLE "
        "(add --no-report for a machine-readable-only export).\n\n"
        "## Reproducibility limits\n\n"
        "This folder is a self-contained audit record, not a vendored installation of every external "
        "program or the complete packaged reference database. The exact query, scanned genome, gene "
        "model, retained homolog promoters, source alignment when available, effective parameters, "
        "software/data fingerprints and every recorded result are included. Re-executing network or "
        "licensed engines still requires those engines and the identified package data release.\n"
    )


def _tsv(path: Path, columns, rows) -> None:
    """Quote embedded tabs/newlines so annotations cannot corrupt table boundaries."""
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sequence_record(sequence: str | None) -> dict:
    sequence = str(sequence or "").upper()
    return {
        "length": len(sequence),
        "sha256": _sha256_bytes(sequence.encode("ascii")) if sequence else None,
    }


def _external_path_references(value, trail: str = "record") -> list[str]:
    """Return machine-local absolute paths hidden in a record, excluding web identifiers."""
    issues = []
    if isinstance(value, dict):
        for key, child in value.items():
            issues.extend(_external_path_references(child, f"{trail}.{key}"))
    elif isinstance(value, (list, tuple)):
        for index, child in enumerate(value):
            issues.extend(_external_path_references(child, f"{trail}[{index}]"))
    elif isinstance(value, str) and "://" not in value and not value.startswith("mailto:"):
        if PureWindowsPath(value).is_absolute() or PurePosixPath(value).is_absolute():
            issues.append(f"{trail} points outside the bundle: {value}")
    return issues


def _copy_recorded_file(source, destination: Path) -> dict | None:
    """Copy a used cache/file input into the bundle and return portable metadata."""
    if not source:
        return None
    try:
        src = Path(source)
        if not src.is_file():
            return None
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, destination)
        payload = destination.read_bytes()
        return {
            "path": destination.as_posix(),
            "source_name": src.name,
            "bytes": len(payload),
            "sha256": _sha256_bytes(payload),
        }
    except (OSError, TypeError, ValueError):
        return None


def _portable_structure(struct: dict | None, bundle: Path) -> dict:
    """Copy structure inputs out of caches and rewrite their paths relative to the bundle."""
    payload = copy.deepcopy(struct or {})
    targets = {
        "dimer_cif": "structure/apo_oligomer.cif",
        "afdb_monomer_cif": "structure/afdb_monomer.cif",
    }
    copied = {}
    for field, relative in targets.items():
        original = payload.get(field)
        record = _copy_recorded_file(original, bundle / relative)
        if record:
            payload[field] = relative
            record["path"] = relative
            copied[field] = record
        else:
            payload[field] = None
            if original:
                payload.setdefault("flags", []).append(f"{field}_not_packaged")
    if copied:
        payload["packaged_files"] = copied
    return payload


def _write_genome_snapshot(source_context: dict, destination: Path) -> dict:
    genome = copy.deepcopy(source_context.get("genome") or {})
    sequence = str(genome.pop("sequence", "") or "").upper()
    genes = list(genome.pop("genes", []) or [])
    destination.mkdir(parents=True, exist_ok=True)
    accession = genome.get("accession") or "genome"
    _fasta([(str(accession), sequence)], destination / "genome.fna")
    offsets = genome.get("contig_offsets") or {}
    ordered_offsets = sorted(offsets.items(), key=lambda item: item[1])

    def locate(position):
        selected, offset = accession, 0
        for contig, candidate in ordered_offsets:
            if candidate <= int(position or 0):
                selected, offset = contig, candidate
        return selected, int(position or 0) - offset

    gene_rows = []
    for gene in genes:
        contig, local_start = locate(gene.get("start"))
        _, local_end = locate(gene.get("end"))
        gene_rows.append([
            contig, local_start, local_end,
            gene.get("start"), gene.get("end"), gene.get("strand"), gene.get("name", ""),
            gene.get("locus_tag", ""), gene.get("protein_id", ""), gene.get("product", ""),
        ])
    _tsv(
        destination / "genes.tsv",
        ["replicon", "replicon_start", "replicon_end", "scan_start", "scan_end", "strand",
         "name", "locus_tag", "protein_id", "product"],
        gene_rows,
    )
    (destination / "contig_offsets.json").write_text(
        json.dumps(offsets, indent=2), encoding="utf-8"
    )
    return {
        **genome,
        **_sequence_record(sequence),
        "gene_count": len(genes),
        "files": {
            "sequence": "genome/genome.fna",
            "genes": "genome/genes.tsv",
            "contig_offsets": "genome/contig_offsets.json",
        },
    }


def _stage_ledger(dossier: dict, run_record: dict) -> list[dict]:
    """A compact human-readable index over the complete machine-readable record."""
    rescan = dossier.get("rescan") or {}
    seed = dossier.get("seed_selection") or {}
    homologs = dossier.get("homolog_regions") or {}
    inducers = dossier.get("inducers") or {}
    structure = dossier.get("structure") or {}
    regulon = [r for r in (dossier.get("regulon") or []) if isinstance(r, dict)]
    return [
        {"stage": "input", "status": "recorded", "used": "normalized TF amino-acid sequence",
         "result": f"{run_record['query']['length']} aa", "artifact": "input/query.fasta"},
        {"stage": "family", "status": "complete" if dossier.get("family") else "unresolved",
         "used": "requested family or family-classification evidence",
         "result": dossier.get("family") or "unresolved", "artifact": "tf/tf_summary.json"},
        {"stage": "genome", "status": "complete" if dossier.get("genome_accession") else "unresolved",
         "used": "exact genome sequence and parsed gene model",
         "result": f"{rescan.get('region_len') or 0:,} bp; {run_record['genome'].get('gene_count', 0):,} genes",
         "artifact": "genome/genome.fna"},
        {"stage": "TF locus", "status": "complete" if dossier.get("tf_locus") else "unresolved",
         "used": "query-vs-genome locus resolution",
         "result": str(dossier.get("tf_locus") or "unresolved"), "artifact": "genome/genome.json"},
        {"stage": "homolog promoters", "status": "complete" if homologs.get("expanded") else "limited",
         "used": "SSN/BLAST promoters with optional MSA expansion",
         "result": f"{len(homologs.get('expanded') or [])} retained regions",
         "artifact": "homologs/homolog_promoters.fasta"},
        {"stage": "motif seeds", "status": "complete" if seed.get("trials") else "unresolved",
         "used": "palindrome candidates and homolog conservation",
         "result": f"{len(seed.get('trials') or [])} width trial(s); chosen w={seed.get('chosen_w')}",
         "artifact": "motif/seed_selection.json"},
        {"stage": "genome rescan", "status": "complete" if rescan.get("hits") else "no hits",
         "used": "chosen seed PWM and locality prior",
         "result": f"{len(rescan.get('hits') or [])} putative sites",
         "artifact": "binding_sites.tsv"},
        {"stage": "inducer", "status": "complete" if inducers.get("top") else "unresolved",
         "used": f"{len(inducers.get('calls') or [])} recorded evidence source(s)",
         "result": inducers.get("top") or "unresolved", "artifact": "ligand/inducer.json"},
        {"stage": "regulon", "status": "complete" if regulon else "empty",
         "used": "motif-coherent genome hits and gene annotation",
         "result": f"{len(regulon)} proposed operon record(s)", "artifact": "regulon.tsv"},
        {"stage": "structure", "status": "complete" if structure.get("folded") else "not run/available",
         "used": "optional apo-oligomer model and AFDB QC",
         "result": structure.get("backend") or "no structure model", "artifact": "structure/structure.json"},
        {"stage": "operator ranking", "status": "complete" if (dossier.get("operators") or {}).get("primary") else "unresolved",
         "used": "recorded conservation-first candidates",
         "result": f"{len((dossier.get('operators') or {}).get('candidates') or [])} candidate(s)",
         "artifact": "operators/FINAL_operators.tsv"},
    ]


def _evidence_index(dossier: dict, run_record: dict) -> list[dict]:
    """Map every evidence family to its complete artifact and HTML section."""
    entries = [
        ("Normalized TF input", run_record["query"].get("length"), "input/query.fasta", "tf"),
        ("Effective run parameters", len(run_record.get("effective_config") or {}), "input/run_parameters.json", "run"),
        ("Scanned genome", run_record["genome"].get("length"), "genome/genome.fna", "run"),
        ("Parsed gene annotations", run_record["genome"].get("gene_count"), "genome/genes.tsv", "run"),
        ("Homolog promoter inputs", len((dossier.get("homolog_regions") or {}).get("expanded") or []),
         "homologs/homolog_promoters.fasta", "operators"),
        ("Seed-width trials", len((dossier.get("seed_selection") or {}).get("trials") or []),
         "motif/seed_selection.json", "operators"),
        ("Genome-wide binding sites", len((dossier.get("rescan") or {}).get("hits") or []),
         "binding_sites.tsv", "operators"),
        ("Per-hit promoter records", len(dossier.get("per_hit_regulation") or []),
         "regulation/per_hit_regulation.json", "operators"),
        ("Inducer evidence calls", len((dossier.get("inducers") or {}).get("calls") or []),
         "ligand/inducer.json", "regulon"),
        ("Regulon records", len(dossier.get("regulon") or []), "regulon.tsv", "regulon"),
        ("Ranked operator candidates", len((dossier.get("operators") or {}).get("candidates") or []),
         "operators/operators.json", "final"),
        ("AF3 hand-off jobs", len(dossier.get("af3_jobs") or []), "af3_jobs.json", "final"),
        ("Complete dossier fields", len(dossier), "dossier.json", "run"),
    ]
    return [
        {"evidence": label, "records": int(count or 0), "artifact": artifact,
         "report_anchor": anchor, "status": "recorded" if count else "empty/not available"}
        for label, count, artifact, anchor in entries
    ]


def _write_traceability(dossier: dict, bundle: Path, source_context: dict | None) -> dict:
    """Write normalized inputs and the audit ledger; return the portable run record."""
    source_context = copy.deepcopy(source_context or {})
    input_dir = bundle / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    query_sequence = str(
        (source_context.get("query") or {}).get("sequence")
        or (dossier.get("tf_record") or {}).get("sequence")
        or ""
    ).upper()
    _fasta([(str(dossier.get("tf_id") or "TF"), query_sequence)], input_dir / "query.fasta")
    parameters = {
        "requested": source_context.get("requested") or {},
        "effective_config": source_context.get("effective_config") or {},
        "input_kind": (source_context.get("query") or {}).get("input_kind") or "normalized_sequence",
    }
    (input_dir / "run_parameters.json").write_text(json.dumps(parameters, indent=2), encoding="utf-8")

    genome_record = _write_genome_snapshot(source_context, bundle / "genome")
    windows = copy.deepcopy(source_context.get("windows") or {})
    for name, record in windows.items():
        if not isinstance(record, dict):
            continue
        sequence = record.pop("sequence", "") or ""
        record.update(_sequence_record(sequence))
        relative = f"genome/{name}.fasta"
        _fasta([(f"{dossier.get('genome_accession') or 'genome'}:{record.get('start')}-{record.get('end')}", sequence)],
               bundle / relative)
        record["artifact"] = relative

    homolog_record = copy.deepcopy(source_context.get("homolog_collection") or {})
    (bundle / "homologs").mkdir(parents=True, exist_ok=True)
    alignment_path = ((homolog_record.get("msa") or {}).pop("alignment_path", None)
                      if isinstance(homolog_record.get("msa"), dict) else None)
    alignment = _copy_recorded_file(alignment_path, bundle / "homologs" / "source_alignment.a3m")
    if alignment:
        alignment["path"] = "homologs/source_alignment.a3m"
        homolog_record.setdefault("msa", {})["packaged_alignment"] = alignment
    metadata_by_sequence = defaultdict(deque)
    for attempt in homolog_record.get("attempts") or []:
        record = attempt.get("record") or {}
        for sequence, metadata in zip(record.get("regions") or [], record.get("sources") or []):
            metadata_by_sequence[str(sequence)].append({"route": attempt.get("route"), **(metadata or {})})
    for record in (homolog_record.get("msa") or {}).get("source_records") or []:
        metadata_by_sequence[str(record.get("sequence") or "")].append({"route": "msa_expansion", **record})
    promoter_rows = []
    for index, sequence in enumerate((dossier.get("homolog_regions") or {}).get("expanded") or [], 1):
        metadata = metadata_by_sequence[str(sequence)].popleft() if metadata_by_sequence[str(sequence)] else {}
        promoter_rows.append([
            f"homolog_promoter_{index}", metadata.get("route", "unmapped_source"),
            metadata.get("uniprot", ""), metadata.get("protein", ""), metadata.get("nuc", ""),
            metadata.get("organism", ""), metadata.get("identity", ""),
            len(sequence), _sequence_record(sequence).get("sha256"),
        ])
    _tsv(
        bundle / "homologs" / "promoter_records.tsv",
        ["export_id", "source_route", "uniprot", "protein", "nucleotide_accession", "organism",
         "identity", "length", "sequence_sha256"],
        promoter_rows,
    )
    (bundle / "homologs" / "source_record.json").write_text(
        json.dumps(homolog_record, indent=2), encoding="utf-8"
    )

    motif_dir = bundle / "motif"
    motif_dir.mkdir(parents=True, exist_ok=True)
    (motif_dir / "seed_selection.json").write_text(
        json.dumps(dossier.get("seed_selection") or {}, indent=2), encoding="utf-8"
    )
    (motif_dir / "literature_seed_candidates.json").write_text(
        json.dumps(dossier.get("literature_seed_candidates") or [], indent=2), encoding="utf-8"
    )
    (motif_dir / "seed_pwm.json").write_text(
        json.dumps(dossier.get("seed_logo") or {}, indent=2), encoding="utf-8"
    )

    run_record = {
        "schema_version": 1,
        "run_id": source_context.get("run_id") or uuid.uuid4().hex,
        "started_at_utc": source_context.get("started_at_utc"),
        "published_at_utc": datetime.now(timezone.utc).isoformat(),
        "query": {**_sequence_record(query_sequence), "artifact": "input/query.fasta"},
        "requested": parameters["requested"],
        "effective_config": parameters["effective_config"],
        "genome": genome_record,
        "windows": windows,
        "homolog_collection": homolog_record,
        "software": runtime_provenance(),
        "dossier_fields": sorted(dossier),
        "notes": [
            "Paths in this record are relative to the individual TF report directory.",
            "Credentials, host names and original cache locations are intentionally excluded.",
        ],
    }
    run_record["stages"] = _stage_ledger(dossier, run_record)
    evidence = _evidence_index(dossier, run_record)
    run_record["evidence_index"] = evidence
    run_dir = bundle / "run"
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "run_record.json").write_text(json.dumps(run_record, indent=2), encoding="utf-8")
    (run_dir / "evidence_index.json").write_text(json.dumps(evidence, indent=2), encoding="utf-8")
    data_manifest = Path(__file__).resolve().parents[1] / "data" / "MANIFEST.json"
    if data_manifest.is_file():
        shutil.copy2(data_manifest, run_dir / "reference_data_manifest.json")
    return run_record

def job_path(tf_id: str, *, root: Path = JOBS) -> Path:
    """Return a contained job path without touching the filesystem."""
    root = Path(root).resolve()
    name = str(tf_id).strip()
    if not name or name in {".", ".."} or Path(name).is_absolute() or Path(name).name != name:
        raise ValueError("job name must be one non-empty path component")
    if any(ord(char) < 32 or char in '<>:"/\\|?*' for char in name):
        raise ValueError(f"job name contains a filesystem-unsafe character: {name!r}")
    if len(name) > 100 or name.endswith((" ", ".")):
        raise ValueError("job name must be at most 100 characters and cannot end in a space or period")
    if name.upper().split(".", 1)[0] in {"CON", "PRN", "AUX", "NUL", *(f"COM{i}" for i in range(1, 10)),
                                         *(f"LPT{i}" for i in range(1, 10))}:
        raise ValueError(f"job name is reserved by Windows: {name!r}")
    d = (root / name).resolve()
    if d.parent != root:
        raise ValueError("job path escapes the configured output root")
    return d


def write_binding_sites(dossier: dict, path: Path) -> int:
    """Every putative binding site with coordinates and calibrated significance -> ranked TSV."""
    hits = (dossier.get("rescan") or {}).get("hits") or []
    neigh = dossier.get("neighborhood") or {}
    win = neigh.get("window")
    cols = ["rank", "accession", "start", "end", "strand", "dyad_center", "score", "pvalue",
            "qvalue", "qvalue_kind", "generator", "in_TF_neighborhood", "sequence"]
    rows = []
    for i, h in enumerate(sorted(hits, key=lambda h: -(h.get("score") or 0)), 1):
        in_n = bool(win and win[0] <= h.get("dyad", -1) <= win[1])
        rows.append([i, h.get("accession", ""), h.get("start"), h.get("end"), h.get("strand", "+"),
                     h.get("dyad"), (f"{h['score']:.3f}" if h.get("score") is not None else ""),
                     (f"{h['pvalue']:.3e}" if h.get("pvalue") is not None else ""),
                     (f"{h['qvalue']:.3e}" if h.get("qvalue") is not None else ""),
                     h.get("qvalue_kind", ""), h.get("generator", ""), in_n, h.get("seq", "")])
    _tsv(path, cols, rows)
    return len(rows)


def write_regulon(dossier: dict, path: Path) -> int:
    """The proposed regulon -> TSV: one row per operon the TF likely drives (operator-PWM genome scan ->
    operons), with the operator site, q-value and locality tier. Empty file (header only) when no regulon."""
    operons = [o for o in (dossier.get("regulon") or []) if isinstance(o, dict) and "first_gene" in o]
    cols = ["rank", "first_gene", "strand", "operon_genes", "site_start", "site_end", "qvalue", "tier"]
    rows = []
    for i, o in enumerate(sorted(operons, key=lambda o: (o.get("qvalue") if o.get("qvalue") is not None
                                                          else 1.0)), 1):
        site = o.get("site") or [None, None]
        rows.append([i, o.get("first_gene", ""), o.get("strand", ""),
                     ";".join(o.get("genes") or []), site[0], site[1],
                     (f"{o['qvalue']:.3g}" if o.get("qvalue") is not None else ""), o.get("tier", "")])
    _tsv(path, cols, rows)
    return len(rows)


# --------------------------------------------------------------------------- richer data exports (report)
def _fasta(records, path: Path) -> int:
    """records = list[(header, seq)] -> FASTA; returns count."""
    lines = []
    for hdr, seq in records:
        if seq:
            sequence = str(seq).replace("\n", "").replace("\r", "")
            wrapped = "\n".join(sequence[i:i + 80] for i in range(0, len(sequence), 80))
            lines.append(f">{hdr}\n{wrapped}")
    path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")
    return len(lines)


def write_tf_fasta(dossier: dict, tf_dir: Path) -> dict:
    """tf/tf.fasta + tf/tf_summary.json from the threaded TF record (best-effort; skips the FASTA if the
    sequence wasn't captured)."""
    tf_dir.mkdir(parents=True, exist_ok=True)
    rec = dossier.get("tf_record") or {}
    seq = rec.get("sequence") or ""
    summary = {"tf_id": dossier.get("tf_id"), "family": dossier.get("family"),
               "family_support": dossier.get("family_support"),
               "organism": rec.get("organism"), "genome_accession": dossier.get("genome_accession"),
               "tf_locus": dossier.get("tf_locus"), "uniprot": rec.get("uniprot"),
               "nearest_uniprot": rec.get("nearest_uniprot"),
               "nearest_tf": rec.get("nearest_tf"), "pfam": rec.get("pfam"),
               "nearest_identity": rec.get("nearest_identity"),
               "nearest_evalue": rec.get("nearest_evalue"),
               "classification_flags": rec.get("flags") or [],
               "ssn_cluster": dossier.get("ssn_cluster"), "length_aa": len(seq) if seq else None}
    (tf_dir / "tf_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    if seq:
        _fasta([(f"{dossier.get('tf_id','TF')} {dossier.get('family','')}".strip(), seq)],
               tf_dir / "tf.fasta")
    return summary


def write_homologs(dossier: dict, hom_dir: Path) -> int:
    """homologs/homolog_promoters.fasta from the threaded homolog promoter regions (pre-MSA set)."""
    hom_dir.mkdir(parents=True, exist_ok=True)
    hr = dossier.get("homolog_regions") or {}
    regs = hr.get("pre_msa") or hr.get("expanded") or []
    recs = [(f"homolog_promoter_{i+1}", s) for i, s in enumerate(regs) if s]
    return _fasta(recs, hom_dir / "homolog_promoters.fasta")


def write_aligned_operators(dossier: dict, op_dir: Path) -> dict:
    """operators/aligned_operators_{tf_only,homolog,gathered}.fasta -- the aligned operator instances behind
    each logo (threaded from operator_logo.three_logos(return_aligned=True))."""
    op_dir.mkdir(parents=True, exist_ok=True)
    ms = dossier.get("multi_source_logos") or {}
    counts = {}
    for key in ("tf_only", "homolog", "gathered"):
        seqs = (ms.get(key) or {}).get("aligned_seqs") or []
        recs = [(f"{key}_op_{i+1}", s) for i, s in enumerate(seqs) if s]
        counts[key] = _fasta(recs, op_dir / f"aligned_operators_{key}.fasta")
    return counts


def write_all_operators(dossier: dict, path: Path) -> int:
    """operators/all_operators.tsv -- every genome operator merged with its per-hit regulation (gene, mode,
    TSS) + the surrounding genes; the enriched, usable operator list."""
    hits = (dossier.get("rescan") or {}).get("hits") or []
    neigh = dossier.get("neighborhood") or {}
    win = neigh.get("window")
    genes = neigh.get("genes") or []
    by_dyad = {}
    auto = None
    for r in (dossier.get("per_hit_regulation") or []):
        op = r.get("operator") or {}
        if op.get("dyad") is not None:
            by_dyad[op["dyad"]] = r
    # autoregulatory dyad
    center = neigh.get("center")
    if center is not None and hits:
        auto = min((h for h in hits if h.get("dyad") is not None),
                   key=lambda h: abs(h["dyad"] - center), default=None)
    auto_dyad = auto["dyad"] if auto else None
    cols = ["rank", "accession", "start", "end", "strand", "dyad", "score", "in_TF_neighborhood",
            "is_autoregulatory", "regulated_gene", "gene_product", "operator_to_gene_bp", "mode",
            "tss", "op_to_tss_bp", "sequence", "surrounding_genes"]
    rows = []
    for i, h in enumerate(sorted(hits, key=lambda h: -(h.get("score") or 0)), 1):
        dy = h.get("dyad")
        r = by_dyad.get(dy, {})
        sur = ";".join(g.get("name", "") for g in genes
                       if abs(((g.get("start", 0) + g.get("end", 0)) // 2) - (dy or 0)) <= 3000)[:120]
        rows.append([i, h.get("accession", ""), h.get("start"), h.get("end"), h.get("strand", "+"),
                     dy, (f"{h['score']:.3f}" if h.get("score") is not None else ""),
                     bool(win and dy is not None and win[0] <= dy <= win[1]),
                     (dy is not None and dy == auto_dyad),
                     r.get("regulated_gene", ""), r.get("gene_product", ""),
                     r.get("operator_to_gene_bp", ""), r.get("mode", ""), r.get("tss", ""),
                     r.get("op_to_tss_bp", ""), h.get("seq", ""), sur])
    _tsv(path, cols, rows)
    return len(rows)


def write_final_operators(dossier: dict, op_dir: Path) -> int:
    """operators/FINAL_operators.{tsv,fasta} -- the headline list: source category + rationale + status."""
    from predictor.report import figure_data as _fd
    op_dir.mkdir(parents=True, exist_ok=True)
    rows = _fd.final_operator_rows(dossier)
    cols = ["source", "category", "sequence", "score", "score_type", "status", "rationale"]
    trows = [[r.get("source", ""), r.get("category", ""), r.get("sequence") or "",
              (f'{r["score"]:.3f}' if isinstance(r.get("score"), (int, float)) else ""),
              r.get("score_type", ""), r.get("status", ""), r.get("rationale", "")] for r in rows]
    _tsv(op_dir / "FINAL_operators.tsv", cols, trows)
    _fasta([(f'{r.get("source")}|{r.get("category")}|{r.get("status")}', r.get("sequence"))
            for r in rows if r.get("sequence")], op_dir / "FINAL_operators.fasta")
    return len(rows)


# --------------------------------------------------------------------------- Stage A (pre-AF3) folder
def write_operators_tsv(operators: dict, path: Path) -> int:
    """The multi-source operator table -> TSV (one row per candidate, ranked sources first)."""
    cands = operators.get("candidates") or []
    cols = ["source", "status", "score", "score_type", "sequence", "provenance"]
    rows = []
    for c in cands:
        rows.append([c.get("source", ""), c.get("status", ""),
                     (f"{c['score']:.3f}" if c.get("score") is not None else ""),
                     c.get("score_type", ""), c.get("sequence") or "",
                     json.dumps(c.get("provenance") or {}, separators=(",", ":"))])
    _tsv(path, cols, rows)
    return len(rows)


def write_stage_a(dossier: dict, *, struct: dict | None = None, operators: dict | None = None,
                  source_context: dict | None = None, root: Path = JOBS,
                  render: bool = True) -> Path:
    """The self-contained per-TF folder with structured subfolders. Writes the structure / operators /
    ligand / genome sections alongside the dossier and the optional AF3 hand-off. This bundle is the
    COMPLETE result: nothing later fills it in.

        results/jobs/<TF>/
          ├── dossier.json  operators_ranked.json  binding_sites.tsv  regulon.tsv  af3_jobs.json
          ├── structure/structure.json          (apo dimer fold + AFDB-monomer QC; only with --fold)
          ├── operators/operators.json + .tsv    (multi-source candidates; WS5)
          ├── ligand/inducer.json                (inducer consensus + regulated-gene promoters)
          ├── genome/genome.json                 (accession, TF locus, neighborhood, promoter occlusion)
          ├── motif/three_logos.json             (TF-only / homolog / gathered logos; WS3)
          └── regulation/per_hit_regulation.json (per-hit TSS/occlusion + gene support; WS4)
    """
    tf = dossier["tf_id"]
    root = Path(root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    final = job_path(tf, root=root)
    if (final / "operators_ranked.json").exists():
        raise FileExistsError(
            f"completed job already exists: {final}; use a new --name or archive the existing bundle"
        )
    if not operators or not operators.get("primary"):
        raise ValueError("Stage-A bundle requires a primary operator")
    d = Path(tempfile.mkdtemp(prefix=f".{tf}.stage-", dir=root))
    try:
        return _write_stage_a_contents(dossier, d=d, final=final, struct=struct, operators=operators,
                                       source_context=source_context, render=render)
    except BaseException:
        shutil.rmtree(d, ignore_errors=True)
        raise


def _write_stage_a_contents(dossier: dict, *, d: Path, final: Path, struct: dict | None,
                            operators: dict | None, source_context: dict | None,
                            render: bool) -> Path:
    """Write a complete bundle into staging, then atomically promote it."""
    portable = copy.deepcopy(dossier)
    if operators is not None:
        portable["operators"] = copy.deepcopy(operators)
    tf = portable["tf_id"]

    # structure/ (WS1)
    (d / "structure").mkdir(parents=True, exist_ok=True)
    portable_structure = _portable_structure(struct or portable.get("structure"), d)
    portable["structure"] = portable_structure
    (d / "structure" / "structure.json").write_text(
        json.dumps(portable_structure, indent=2), encoding="utf-8"
    )

    # The normalized inputs and run ledger are written before the dossier so every path inside the
    # published JSON is relative to this directory and every cache-backed file has already been copied.
    run_record = _write_traceability(portable, d, source_context)
    msa_expansion = (portable.get("conservation") or {}).get("msa_expansion")
    if isinstance(msa_expansion, dict) and "alignment_path" in msa_expansion:
        packaged = ((run_record.get("homolog_collection") or {}).get("msa") or {}).get("packaged_alignment")
        msa_expansion["alignment_path"] = "homologs/source_alignment.a3m" if packaged else None
    portable["run_record"] = {
        "run_id": run_record["run_id"],
        "started_at_utc": run_record.get("started_at_utc"),
        "published_at_utc": run_record["published_at_utc"],
        "artifact": "run/run_record.json",
    }
    external_paths = _external_path_references(portable) + _external_path_references(run_record)
    if external_paths:
        raise ValueError("bundle record retains a machine-local path: " + external_paths[0])
    (d / "dossier.json").write_text(json.dumps(portable, indent=2), encoding="utf-8")
    n = write_binding_sites(portable, d / "binding_sites.tsv")
    nr = write_regulon(portable, d / "regulon.tsv")
    (d / "af3_jobs.json").write_text(json.dumps(portable.get("af3_jobs") or [], indent=2), encoding="utf-8")

    # operators/ (WS5) + the top-level ranked operator list
    no = 0
    (d / "operators").mkdir(parents=True, exist_ok=True)
    if operators is not None:
        (d / "operators" / "operators.json").write_text(json.dumps(operators, indent=2), encoding="utf-8")
        no = write_operators_tsv(operators, d / "operators" / "operators_ranked.tsv")
        (d / "operators_ranked.json").write_text(
            json.dumps({"ranked": operators.get("ranked"), "primary": operators.get("primary")}, indent=2),
            encoding="utf-8")

    # ligand/ (inducer + regulated-gene promoters)
    (d / "ligand").mkdir(parents=True, exist_ok=True)
    (d / "ligand" / "inducer.json").write_text(json.dumps({
        "inducers": portable.get("inducers") or {}, "effector": portable.get("effector"),
        "ion": portable.get("ion"), "ligand": portable.get("ligand"),
        "ligand_genes": portable.get("ligand_genes") or [],
        "ligand_keywords": portable.get("ligand_keywords") or []}, indent=2), encoding="utf-8")

    # genome/ (accession + locus + neighborhood + mode)
    (d / "genome").mkdir(parents=True, exist_ok=True)
    (d / "genome" / "genome.json").write_text(json.dumps({
        "genome_accession": portable.get("genome_accession"), "tf_locus": portable.get("tf_locus"),
        "genome_topology": portable.get("genome_topology") or {},
        "neighborhood": portable.get("neighborhood") or {}, "operator_window": portable.get("operator_window"),
        "promoter_occlusion": portable.get("promoter_occlusion") or {},
        "autoregulatory_recovery": (portable.get("autoregulatory_recovery")
                                    or portable.get("cognate_recovery") or {}),
        "target_operator": portable.get("target_operator") or {},
        "run_input": run_record.get("genome") or {}}, indent=2), encoding="utf-8")

    # motif/ + regulation/ (WS3/WS4)
    (d / "motif").mkdir(parents=True, exist_ok=True)
    (d / "motif" / "three_logos.json").write_text(
        json.dumps(portable.get("multi_source_logos") or {}, indent=2), encoding="utf-8")
    (d / "regulation").mkdir(parents=True, exist_ok=True)
    (d / "regulation" / "per_hit_regulation.json").write_text(json.dumps({
        "per_hit": portable.get("per_hit_regulation") or [],
        "confirmed_genes": portable.get("confirmed_regulated_genes") or []}, indent=2), encoding="utf-8")

    # richer report exports (TF fasta+summary, homolog promoters, aligned-operator FASTAs, enriched + final tables)
    write_tf_fasta(portable, d / "tf")
    write_homologs(portable, d / "homologs")
    write_aligned_operators(portable, d / "operators")
    write_all_operators(portable, d / "operators" / "all_operators.tsv")
    write_final_operators(portable, d / "operators")
    from predictor.report.environment import snapshot
    (d / "report_environment.json").write_text(json.dumps(snapshot(), indent=2), encoding="utf-8")
    (d / "README.md").write_text(bundle_readme(), encoding="utf-8")

    # A requested report is a required deliverable. Failed rendering must not publish success.
    if render:
        from predictor.report import render_report as _rr
        _rr.render(d, pending_manifest=True)

    required = ["dossier.json", "operators_ranked.json", "binding_sites.tsv", "regulon.tsv"]
    if render:
        required.append("REPORT.html")
    missing = [rel for rel in required if not (d / rel).is_file()]
    if missing:
        raise RuntimeError(f"Stage-A staging bundle is missing required artifact(s): {', '.join(missing)}")
    manifest = {
        "schema_version": 2,
        "status": "complete",
        "tf_id": tf,
        "required_artifacts": required,
        "runtime": runtime_provenance(),
        "report_generated": (d / "REPORT.html").is_file(),
        "warnings": [],
        "artifacts": inventory(d),
    }
    (d / "bundle_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    if render:
        issues = report_asset_issues(d)
        if issues:
            raise RuntimeError("report assets are incomplete: " + "; ".join(issues))

    backup = None
    if final.exists():
        backup = final.parent / f".{final.name}.incomplete-{uuid.uuid4().hex[:10]}"
        final.replace(backup)
    try:
        d.replace(final)
    except BaseException:
        if backup is not None and backup.exists() and not final.exists():
            backup.replace(final)
        raise
    print(f"Stage-A bundle -> {final}  ({n} binding sites, {nr} regulon operons, {no} operator candidates, "
          f"{len(portable.get('af3_jobs') or [])} AF3 jobs)"
          + (f"; previous incomplete bundle preserved at {backup}" if backup else ""))
    return final
