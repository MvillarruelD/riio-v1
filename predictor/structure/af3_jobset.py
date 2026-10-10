"""af3_jobset.py -- the AlphaFold-3 job SET for one TF: apo, operator-bound, metal-bound, and both.

`af3_export` turns operator picks into AF Server jobs. This module decides WHICH jobs a TF gets and
writes them for both places AlphaFold 3 runs:

  kind            what is folded                                   why
  --------------  -----------------------------------------------  ---------------------------------
  apo             the oligomer alone                               the reference conformation
  dna_r1..dna_rN  the oligomer + operator rank k (60-bp duplex)    which operator the protein reads
  metal           the oligomer + its cognate metal                 does the predicted site form?
  metal_dna_r1    the oligomer + metal + first accepted operator    an exploratory loaded complex

`metal` and `metal_dna_r1` exist only for TFs whose inducer is a metal or metalloid.

**Two formats.**
* *AF Server* (alphafoldserver.com, `dialect: alphafoldserver`) accepts exactly ten ions
  (`af3_export.AF3_IONS`). A cognate metal outside them is replaced by the closest ion chosen for
  COORDINATION, and the job is flagged `analog` in the manifest. Fe2+ -> MN (Mn(II), the redox-inert
  Fe(II) stand-in used in Fur and DtxR/IdeR crystallography) is the case that matters most. A metal
  with no defensible stand-in (molybdate) gets no server metal job, and the manifest says so.
* *Local AlphaFold 3* (`dialect: alphafold3`) takes any PDB Chemical Component Dictionary code, so it
  gets verified CCD species: FE2, CU1, NI, CD, HG, PB, AST (free arsenite), MOO (molybdate).
  No local antimonite recipe is asserted: SB is a bare Sb(III) atom, not the oxyanion.

**Counts.** The protein is folded at its family's oligomeric state (`structure.oligomer`, plus NikR as
a tetramer here), and ions are counted per protomer from `FAMILY_SITES` -- one site per protomer by
default; DtxR/IdeR/MntR uses two as an explicit export assumption. Structural Zn and additional
Fur sites are not inferred from family membership.

**Operators** come from the bundle's ranked list in its own order, deduplicated by locus AND by
sequence: a repeated element (four identical ArsR sites in E. coli) would otherwise spend three of the
three DNA slots on one duplex.

Rebuild for a finished bundle or a whole run, without re-predicting:

    tfop af3 RUN_DIR_or_BUNDLE [--top-n 3] [--format both]
    python -m predictor.structure.af3_jobset        # offline self-test
"""
from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path

from predictor.structure import af3_export as X

#: Default number of operator-bound jobs per TF.
TOP_N = 3
#: Duplex length, the same as the pipeline's long-standing hand-off (short duplexes break docking).
DNA_LENGTH = 60
#: Packaging default per upload file; this does not guarantee any account quota.
SERVER_BATCH = 30
#: Local AlphaFold 3 input version written. Version 1 is read by every AF3 release.
LOCAL_VERSION = 1


@dataclass(frozen=True)
class MetalSpec:
    server_ion: str | None     # one of af3_export.AF3_IONS, or None: no defensible server stand-in
    local_ccd: str | None      # PDB CCD code for local AF3, or None
    match: str                 # exact | close | weak | same element, wrong oxidation state | none
    why: str


SOURCE_SERVER = "https://github.com/google-deepmind/alphafold/blob/main/server/README.md"
SOURCE_LOCAL = "https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md"
CCD_SOURCE = "https://www.rcsb.org/ligand/"
# Export recipes, not inducer predictions. AST is the free arsenite oxyanion;
# it does not assert the chemistry of protein-bound As(Cys)3.
METALS = {
    "Zn2+": MetalSpec("ZN", "ZN", "exact", "zinc(II)"),
    "Fe2+": MetalSpec("MN", "FE2", "close", "Mn(II) coordination surrogate; not iron"),
    "Mn2+": MetalSpec("MN", "MN", "exact", "manganese(II)"),
    "Ni2+": MetalSpec("CO", "NI", "close", "Co(II) coordination surrogate; not nickel"),
    "Co2+": MetalSpec("CO", "CO", "exact", "cobalt(II)"),
    "Cu+": MetalSpec("CU", "CU1", "same element, wrong oxidation state", "server CU is Cu(II)"),
    "Cd2+": MetalSpec("ZN", "CD", "close", "Zn(II) surrogate; not cadmium"),
    "Pb2+": MetalSpec("ZN", "PB", "weak", "Zn(II) placeholder does not reproduce lead chemistry"),
    "Hg2+": MetalSpec("ZN", "HG", "weak", "Zn(II) placeholder does not reproduce mercury chemistry"),
    "As(III)": MetalSpec("ZN", "AST", "weak", "server Zn placeholder; local AST is free AsO3, not thiolate-bound arsenic"),
    "antimonite": MetalSpec("ZN", None, "weak", "server Zn placeholder; no verified antimonite CCD (SB is bare Sb(III), not oxyanion)"),
    "MoO42-": MetalSpec(None, "MOO", "none", "no defensible allowed server analogue"),
}
# Preserve Figure 4's existing byte output by leaving its archived recipe alone.
DNA_BOUND_FORM = {"ArsR/SmtB": "apo", "CsoR/FrmR": "apo", "MerR": "apo",
                  "Fur": "holo", "DtxR/MntR": "holo", "NikR": "holo"}
# Site numbers are export assumptions. Family membership alone does not establish
# structural Zn or a second regulatory site, so neither is silently inserted.
FAMILY_SITES = {"dtxr": 2, "ider": 2, "mntr": 2}


@dataclass(frozen=True)
class Job:
    name: str
    kind: str
    protein: str
    copies: int
    dna: str | None = None
    ions: tuple = ()                 # (cognate token, MetalSpec, count)
    operator_rank: int | None = None
    locus: str = ""
    expected_dna_form: str = "unknown"
    assumptions: str = ""


@dataclass
class JobSet:
    jobs: list[Job] = field(default_factory=list)
    diagnostics: list[dict] = field(default_factory=list)


def _value(row, key, default=None):
    return row.get(key, default) if isinstance(row, dict) else getattr(row, key, default)


def ranked_operators(payload):
    """Read canonical ranked objects or legacy flat operator arrays, in stored order."""
    if isinstance(payload, dict):
        payload = payload.get("ranked", payload.get("operators", payload.get("candidates", [])))
    out = []
    for row in payload or []:
        if not isinstance(row, dict):
            out.append(row)
            continue
        provenance = row.get("provenance") or {}
        out.append({**provenance, **row,
                    "genome_accession": row.get("genome_accession") or provenance.get("accession") or row.get("accession", ""),
                    "dyad_center": row.get("dyad_center", provenance.get("dyad")),
                    "seq": row.get("seq") or row.get("sequence", "")})
    return out


def metal_spec(label: str) -> MetalSpec | None:
    from predictor.effector.inducer_vocab import normalise
    return METALS.get(normalise(label))


def build_job_set(seq, ranked_ops, effector, family, genome_seq=None, top_n=TOP_N,
                  *, name_prefix="TF", copies=None, sites_per_protomer=None,
                  genome_accession=None) -> JobSet:
    """Build independent export jobs in ranked order, without modifying predictions.

    A missing/ambiguous 60 bp construct skips only its DNA job. No bases are
    padded or replaced. Counts are explicit modelling assumptions and overridable.
    """
    from predictor.effector.inducer_vocab import is_metal, parse_mixture
    from predictor.structure.oligomer import oligomeric_state, is_default
    if type(top_n) is not int or top_n < 0:
        raise ValueError("top_n must be a nonnegative integer")
    result = JobSet()
    protein = str(seq or "").upper()
    if not protein or any(c not in "ACDEFGHIKLMNPQRSTVWY" for c in protein):
        result.diagnostics.append({"kind": "all", "status": "skipped", "reason": "missing or invalid protein sequence"})
        return result
    fam = str(family or "")
    n = copies if copies is not None else (4 if "nikr" in fam.lower() else oligomeric_state(fam))
    sites = sites_per_protomer if sites_per_protomer is not None else next(
        (v for k, v in FAMILY_SITES.items() if k in fam.lower()), 1)
    if type(n) is not int or n < 1 or type(sites) is not int or sites < 1:
        raise ValueError("copies and sites_per_protomer must be positive integers")
    assumptions = (f"oligomer copies={n} ({'override' if copies is not None else 'family export assumption'}); "
                   f"sites/protomer={sites} ({'override' if sites_per_protomer is not None else 'export assumption'}); "
                   "no structural ions inferred")
    if is_default(fam) and "nikr" not in fam.lower():
        assumptions += "; unknown-family dimer default"
    form = DNA_BOUND_FORM.get(fam, "unknown")
    def add(kind, dna=None, ions=(), rank=None, locus="", suffix=""):
        result.jobs.append(Job(X.sanitize_job_name(f"{name_prefix}__{kind}{suffix}"), kind,
                               protein, n, dna, ions, rank, locus, form, assumptions))
    add("apo")
    loci, sequences, accepted = set(), set(), []
    for rank, op in enumerate(ranked_operators(ranked_ops), 1):
        if len(accepted) >= top_n:
            break
        acc = str(_value(op, "genome_accession", "") or _value(op, "accession", ""))
        try:
            start = int(_value(op, "start", 0) or 0)
            end = int(_value(op, "end", 0) or 0)
            dyad = _value(op, "dyad_center")
            dyad = int(dyad if dyad is not None else (start + end) // 2)
        except (TypeError, ValueError, OverflowError):
            result.diagnostics.append({"kind": "dna", "operator_rank": rank, "status": "skipped",
                                       "reason": "invalid stored operator coordinates"})
            continue
        valid_locus = bool(acc) and 0 <= start < end and start <= dyad <= end
        locus_key = (acc, round(dyad / 10)) if valid_locus else None
        locus = f"{acc}:{start}-{end}" if valid_locus else "not recorded"
        dna = _value(op, "dna_top")
        g = genome_seq.get(acc) if isinstance(genome_seq, dict) else genome_seq
        if genome_accession and acc and acc != genome_accession and not isinstance(genome_seq, dict):
            g = None
        if not dna and g and valid_locus and end <= len(g):
            dna = X.center_construct(str(g), dyad, DNA_LENGTH)
        if not dna:
            candidate = str(_value(op, "seq", "") or "")
            dna = candidate if len(candidate) == DNA_LENGTH else None
        dna = str(dna or "").upper()
        if len(dna) != DNA_LENGTH or any(b not in "ACGT" for b in dna):
            result.diagnostics.append({"kind": "dna", "operator_rank": rank, "locus": locus,
                                       "status": "skipped", "reason": "exact unambiguous 60 bp construct unavailable"})
            continue
        duplex = min(dna, X.revcomp(dna))
        if (locus_key is not None and locus_key in loci) or duplex in sequences:
            result.diagnostics.append({"kind": "dna", "operator_rank": rank, "locus": locus,
                                       "status": "skipped", "reason": "duplicate locus or duplex sequence"})
            continue
        if locus_key is not None:
            loci.add(locus_key)
        sequences.add(duplex)
        accepted.append((rank, dna, locus))
        add(f"dna_r{len(accepted)}", dna, rank=rank, locus=locus)
    if not accepted:
        result.diagnostics.append({"kind": "dna", "status": "skipped", "reason": "no valid ranked operators available"})
    if not is_metal(effector):
        result.diagnostics.append({"kind": "metal", "status": "skipped", "reason": "nonmetal, missing, or mixed-class inducer"})
        return result
    species = parse_mixture(effector)
    for token in species:
        spec = METALS.get(token)
        if spec is None:
            result.diagnostics.append({"kind": "metal", "cognate": token, "status": "skipped",
                                       "reason": "no verified CCD/analogue recipe; inducer unchanged"})
            continue
        suffix = "__" + X.sanitize_job_name(token) if len(species) > 1 else ""
        ions = ((token, spec, n * sites),)
        add("metal", ions=ions, suffix=suffix)
        if accepted:
            rank, dna, locus = accepted[0]
            add("metal_dna_r1", dna, ions, rank, locus, suffix)
        else:
            result.diagnostics.append({"kind": "metal_dna_r1", "cognate": token, "status": "skipped",
                                       "reason": "no valid operator construct"})
    return result


def to_server_json(job: Job) -> dict:
    validate_job(job)
    if any(spec.server_ion is None for _, spec, _ in job.ions):
        raise ValueError("no allowed server analogue")
    payload = X.af3_job(job.name, job.protein, job.dna, copies=job.copies, seeds=())
    payload["modelSeeds"] = ["1"]  # server documentation: uint32 strings
    for _, spec, count in job.ions:
        payload["sequences"].append({"ion": {"ion": spec.server_ion, "count": count}})
    X.validate_jobs([payload])
    return payload


def _chain_id(index):
    result = ""
    while True:
        result = chr(65 + index % 26) + result
        index = index // 26 - 1
        if index < 0:
            return result


def to_local_json(job: Job) -> dict:
    validate_job(job)
    ids = [_chain_id(i) for i in range(job.copies)]
    sequences = [{"protein": {"id": ids, "sequence": job.protein}}]
    cursor = job.copies
    if job.dna is not None:
        for dna in (job.dna, X.revcomp(job.dna)):
            sequences.append({"dna": {"id": _chain_id(cursor), "sequence": dna}})
            cursor += 1
    for _, spec, count in job.ions:
        if spec.local_ccd is None:
            raise ValueError("no verified local CCD")
        sequences.append({"ligand": {"id": [_chain_id(i) for i in range(cursor, cursor + count)],
                                       "ccdCodes": [spec.local_ccd]}})
        cursor += count
    payload = {"name": job.name, "modelSeeds": [1], "sequences": sequences,
               "dialect": "alphafold3", "version": LOCAL_VERSION}
    validate_local(payload)
    return payload


def validate_job(job):
    """Enforce the four recipe kinds before dialect serialization."""
    wants_dna = job.kind.startswith("dna_r") or job.kind == "metal_dna_r1"
    wants_metal = job.kind in ("metal", "metal_dna_r1")
    if job.kind != "apo" and not wants_dna and not wants_metal:
        raise ValueError("unknown AF3 job kind")
    if wants_dna != (job.dna is not None) or wants_metal != bool(job.ions):
        raise ValueError("AF3 job entities disagree with kind")
    if job.dna is not None and (len(job.dna) != DNA_LENGTH or any(c not in "ACGT" for c in job.dna)):
        raise ValueError("job needs an exact 60 bp duplex")
    if any(type(count) is not int or count < 1 for _, _, count in job.ions):
        raise ValueError("ion counts must be positive integers")


def validate_local(payload):
    """Validate the emitted local subset (not the complete upstream AF3 schema)."""
    import re
    if payload.get("dialect") != "alphafold3" or payload.get("version") != LOCAL_VERSION:
        raise ValueError("invalid local dialect/version")
    seeds = payload.get("modelSeeds")
    if not seeds or any(type(x) is not int or x < 0 for x in seeds):
        raise ValueError("local modelSeeds must be nonempty nonnegative integers")
    seen, proteins, dnas = set(), [], []
    for entity in payload.get("sequences", []):
        if len(entity) != 1:
            raise ValueError("one entity type required")
        kind, item = next(iter(entity.items()))
        if kind not in ("protein", "dna", "ligand"):
            raise ValueError("unsupported emitted local entity")
        ids = item.get("id")
        ids = ids if isinstance(ids, list) else [ids]
        if not ids:
            raise ValueError("entity needs ids")
        for identifier in ids:
            if not isinstance(identifier, str) or not re.fullmatch("[A-Z]+", identifier) or identifier in seen:
                raise ValueError("unique uppercase local entity ids required")
            seen.add(identifier)
        if kind == "protein":
            proteins.append(item)
            if not item.get("sequence") or any(c not in "ACDEFGHIKLMNPQRSTVWY" for c in item["sequence"]):
                raise ValueError("invalid protein")
        elif kind == "dna":
            dnas.append(item.get("sequence", ""))
        elif not item.get("ccdCodes") or "smiles" in item:
            raise ValueError("CCD ligand required")
    if len(proteins) != 1 or len(dnas) not in (0, 2):
        raise ValueError("one protein entity and zero/two DNA entities required")
    if dnas and (not dnas[0] or any(c not in "ACGT" for c in "".join(dnas)) or X.revcomp(dnas[0]) != dnas[1]):
        raise ValueError("invalid duplex")


def _read_fasta(path):
    records = {}
    if not path.exists():
        return records
    key = None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            key = line[1:].split()[0]
            records[key] = []
        elif key:
            records[key].append(line.strip().upper())
    return {key: "".join(parts) for key, parts in records.items()}


def _bundle_jobset(bundle, top_n):
    dossier = json.loads((bundle / "dossier.json").read_text(encoding="utf-8"))
    proteins = _read_fasta(bundle / "input" / "query.fasta")
    seq = next(iter(proteins.values()), "") or (dossier.get("tf_record") or {}).get("sequence", "")
    if not seq:
        legacy = dossier.get("af3_jobs") or []
        legacy_path = bundle / "af3_jobs.json"
        if not legacy and legacy_path.exists():
            legacy = json.loads(legacy_path.read_text(encoding="utf-8"))
        seq = next((s["proteinChain"]["sequence"] for j in legacy for s in j.get("sequences", [])
                    if "proteinChain" in s), "")
    op_path = bundle / "operators_ranked.json"
    ops = json.loads(op_path.read_text(encoding="utf-8")) if op_path.exists() else []
    ops = ranked_operators(ops)
    genomes = _read_fasta(bundle / "genome" / "genome.fna")
    acc = str(dossier.get("genome_accession") or "")
    if not genomes and acc:
        # Legacy bundles predate portable genome inputs. Resolve local mirror
        # only; never download or rerun prediction during a rebuild.
        from predictor.annotate.genome_mirror import ensure_genome
        fna, _ = ensure_genome(acc, allow_download=False)
        if fna:
            genomes = _read_fasta(Path(fna))
    if len(genomes) == 1 and acc and acc not in genomes:
        genomes[acc] = next(iter(genomes.values()))
    return build_job_set(seq, ops, dossier.get("effector"), dossier.get("family"), genomes,
                         top_n, name_prefix=str(dossier.get("tf_id") or bundle.name))


def _export(bundles, destination, top_n, format, batch_size):
    if format not in ("server", "local", "both"):
        raise ValueError("format must be server, local, or both")
    if type(batch_size) is not int or batch_size < 1:
        raise ValueError("batch_size must be positive")
    destination = Path(destination).resolve()
    if any((parent / "dossier.json").is_file() for parent in (destination, *destination.parents)):
        raise ValueError("AF3 sidecars must be outside existing bundles")
    for bundle in bundles:
        if destination == bundle.resolve() or bundle.resolve() in destination.parents:
            raise ValueError("AF3 sidecars must be outside existing bundles")
    destination.mkdir(parents=True, exist_ok=True)
    old_files = set()
    old_manifest = destination / "af3_manifest.tsv"
    if old_manifest.exists():
        with old_manifest.open(encoding="utf-8", newline="") as stream:
            old_files = {row.get("file", "") for row in csv.DictReader(stream, delimiter="\t")}
    rows, server_jobs, local_count, warnings = [], [], 0, []
    names = set()
    for bundle in bundles:
        try:
            result = _bundle_jobset(bundle, top_n)
        except (OSError, ValueError, TypeError, KeyError) as exc:
            result = JobSet(diagnostics=[{"kind": "all", "status": "skipped", "reason": str(exc)}])
        for diagnostic in result.diagnostics:
            rows.append({"tf": bundle.name, "format": "all", **diagnostic})
            warnings.append(f"{bundle.name}: {diagnostic['reason']}")
        for job in result.jobs:
            for dialect in ("server", "local"):
                if format not in (dialect, "both"):
                    continue
                row = {"tf": bundle.name, "kind": job.kind, "format": dialect, "job": job.name,
                       "operator_rank": job.operator_rank, "locus": job.locus,
                       "expected_dna_form": job.expected_dna_form, "assumptions": job.assumptions,
                       "cognate": ";".join(t for t, _, _ in job.ions),
                       "ion": ";".join(str(s.server_ion if dialect == "server" else s.local_ccd) for _, s, _ in job.ions),
                       "ion_count": ";".join(str(c) for _, _, c in job.ions),
                       "analog": any(s.match != "exact" for _, s, _ in job.ions) if dialect == "server" else False,
                       "chemistry_note": "; ".join(s.why for _, s, _ in job.ions),
                       "ccd_source": ";".join(CCD_SOURCE + s.local_ccd for _, s, _ in job.ions if s.local_ccd)}
                try:
                    key = (dialect, job.name.casefold())
                    if key in names:
                        raise ValueError("duplicate job name across TFs")
                    payload = to_server_json(job) if dialect == "server" else to_local_json(job)
                    names.add(key)
                    if dialect == "server":
                        server_jobs.append(payload)
                    else:
                        from predictor.report.bundle_contract import contained_file
                        target = contained_file(destination, f"local/{job.name}.json")
                        target.parent.mkdir(parents=True, exist_ok=True)
                        target.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
                        row["file"] = f"local/{job.name}.json"
                        local_count += 1
                    row.update(status="written", reason="")
                except ValueError as exc:
                    row.update(status="skipped", reason=str(exc))
                    warnings.append(f"{bundle.name}/{job.kind}/{dialect}: {exc}")
                rows.append(row)
    for offset in range(0, len(server_jobs), batch_size):
        filename = f"af3_server_batch_{offset // batch_size + 1:02d}.json"
        from predictor.report.bundle_contract import contained_file
        target = contained_file(destination, f"server/{filename}")
        target.parent.mkdir(parents=True, exist_ok=True)
        chunk = server_jobs[offset:offset + batch_size]
        target.write_text(json.dumps(chunk, indent=2) + "\n", encoding="utf-8")
        chunk_names = {j["name"] for j in chunk}
        for row in rows:
            if row.get("format") == "server" and row.get("status") == "written" and row.get("job") in chunk_names:
                row["file"] = f"server/{filename}"
    columns = ["tf", "kind", "format", "status", "reason", "job", "file", "cognate", "ion", "ion_count",
               "analog", "operator_rank", "locus", "expected_dna_form", "assumptions", "chemistry_note", "ccd_source"]
    manifest_path = destination / "af3_manifest.tsv"
    with manifest_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    # Retire only files named in our prior manifest, preserving unrelated files.
    # Resolve and verify targets before unlinking; existing bundles never qualify.
    new_files = {row.get("file", "") for row in rows}
    for relative in old_files - new_files:
        target = (destination / relative).resolve()
        if relative and target.parent in (destination / "server", destination / "local") and target.suffix == ".json":
            target.unlink(missing_ok=True)
    summary = {"output_dir": str(destination), "manifest": str(manifest_path), "server_jobs": len(server_jobs),
               "local_jobs": local_count, "bundles": len(bundles), "warnings": warnings,
               "schema_sources": [SOURCE_SERVER, SOURCE_LOCAL], "policy": "export assumptions; predictions unchanged"}
    (destination / "export_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def export_bundle(bundle, output_dir=None, *, top_n=TOP_N, format="both", batch_size=SERVER_BATCH):
    """Read one bundle; write only sibling/custom sidecars, never bundle contents."""
    bundle = Path(bundle).resolve()
    destination = Path(output_dir) if output_dir else bundle.with_name(bundle.name + "_af3")
    return _export([bundle], destination, top_n, format, batch_size)


def export_run(run_dir, top_n=TOP_N, format="both", batch_size=SERVER_BATCH):
    """Export a run's jobs/* bundles, or one bundle into a sibling directory."""
    run_dir = Path(run_dir).resolve()
    if (run_dir / "dossier.json").exists():
        return export_bundle(run_dir, top_n=top_n, format=format, batch_size=batch_size)
    bundles = sorted(p.parent for p in (run_dir / "jobs").glob("*/dossier.json"))
    if not bundles:
        bundles = sorted(p.parent for p in run_dir.glob("*/dossier.json"))
    if not bundles:
        raise ValueError(f"no completed TF bundles found in {run_dir}")
    return _export(bundles, run_dir / "af3", top_n, format, batch_size)
