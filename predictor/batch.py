"""Explicit protein-list input: every supplied protein reaches the unchanged per-TF core."""
from pathlib import Path
from types import SimpleNamespace
import csv


def run_batch(source, out_dir, *, organism=None, jobs=1, offline=False, fold=False,
              timeout=None, discover_only=False, metalnet=True, af3=True,
              genome_fasta=None, genome_gff=None, genome_genbank=None):
    from predictor import genome_run as gr
    from predictor.input_validation import read_protein
    from predictor.run_cli import _resolve_prediction_source
    source, out_dir = Path(source), Path(out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    if genome_fasta and genome_genbank:
        raise ValueError("choose one local genome FASTA or GenBank context")
    if genome_gff and not genome_fasta:
        raise ValueError("genome_gff requires genome_fasta")
    effective_gff = genome_gff
    if genome_fasta and not effective_gff:
        effective_gff = next((p for p in (Path(genome_fasta).with_suffix(".gff"), Path(genome_fasta).with_suffix(".gff3")) if p.is_file()), None)
    identity = gr.input_identity("protein_list", files=[("proteins", source), ("genome_fasta", genome_fasta),
                                ("genome_gff", effective_gff), ("genome_genbank", genome_genbank)], organism=organism)
    gr.check_run_identity(out_dir, identity)
    manifest = out_dir / "manifest.csv"
    if manifest.exists():
        with manifest.open(encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
    else:
        coords = {}
        if genome_fasta or genome_genbank:
            _local, annotation = gr.prepare_local(genome_fasta, genome_gff, genome_genbank, None, out_dir)
            coords = gr.gff_index(annotation)
        text = source.read_text(encoding="utf-8-sig")
        candidates = []
        if text.lstrip().startswith(">"):
            from Bio import SeqIO
            for record in SeqIO.parse(str(source), "fasta"):
                candidates.append(SimpleNamespace(protein_id=record.id, sequence=read_protein(str(record.seq)), family=""))
        else:
            table = csv.DictReader(text.splitlines(), delimiter="\t")
            if not table.fieldnames or "accession" not in table.fieldnames:
                raise ValueError("accession TSV requires an 'accession' header")
            for row in table:
                accession = (row.get("accession") or "").strip()
                if not accession:
                    raise ValueError("accession TSV contains an empty accession")
                candidates.append(SimpleNamespace(protein_id=accession, sequence=read_protein(_resolve_prediction_source(accession, offline=offline)), family=""))
        if not candidates:
            raise ValueError("protein list contains no proteins")
        rows = gr.manifest_rows(candidates, short=gr.safe(source.stem), organism=organism or "",
                                accession=organism or "", engine="provided_protein", protein_source="user_protein_list",
                                coords=coords, fasta_dir=out_dir / "fastas", rel_to=out_dir)
    return gr.run_candidates(rows, out_dir, jobs=jobs, offline=offline, fold=fold, timeout=timeout,
                             discover_only=discover_only, metalnet=metalnet, af3=af3)
