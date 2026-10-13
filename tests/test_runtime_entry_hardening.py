"""Software-boundary regression tests: inputs, saved-run reuse and sidecar failures."""
import csv
import json
from types import SimpleNamespace

import pytest

from predictor import genome_run as gr
from predictor.batch import run_batch
from predictor.structure.af3_jobset import build_job_set, export_bundle, export_run


PROTEIN = "MKKLTVSDLAREHGYEVVETPEQVTAL"


def rows(directory, identifier="P"):
    candidate = SimpleNamespace(protein_id=identifier, sequence=PROTEIN, family="")
    return gr.manifest_rows([candidate], short="supplied", organism="", accession="",
                            engine="provided_protein", protein_source="user", coords={},
                            fasta_dir=directory / "fastas", rel_to=directory)


@pytest.mark.parametrize("name", ["CON", "con.txt", "trailing.", "a" * 241, "../outside", "a/b"])
def test_runner_rejects_nonportable_or_long_names_before_writes(tmp_path, name):
    with pytest.raises(ValueError):
        gr.run_candidates([{"run_name": name, "fasta": "input.faa"}], tmp_path, discover_only=True)
    assert not (tmp_path / "manifest.csv").exists()


@pytest.mark.parametrize("source", ["../outside.faa", "C:/outside.faa", "fastas/../../outside.faa"])
def test_manifest_fasta_cannot_escape_run(tmp_path, source):
    with pytest.raises(ValueError):
        gr.run_candidates([{"run_name": "safe", "fasta": source}], tmp_path, discover_only=True)
    assert not (tmp_path / "manifest.csv").exists()


def test_manifest_missing_required_fields_gets_actionable_input_error(tmp_path):
    with pytest.raises(ValueError, match="run_name and fasta"):
        gr.run_candidates([{"fasta": "input.faa"}], tmp_path, discover_only=True)


def test_incomplete_accession_table_row_gets_input_error(tmp_path):
    source = tmp_path / "source.tsv"
    source.write_text("label\taccession\nmissing-accession\n")
    with pytest.raises(ValueError, match="empty accession"):
        run_batch(source, tmp_path / "run", discover_only=True)


def test_runner_fingerprint_survives_csv_resume_and_rejects_changed_policy_or_fasta(tmp_path, monkeypatch):
    records = rows(tmp_path)
    calls = []
    monkeypatch.setattr(gr, "predict_all", lambda entries, **kw: calls.append(entries) or {})
    monkeypatch.setattr(gr, "finish_run", lambda *args, **kw: None)
    gr.run_candidates(records, tmp_path, offline=True, af3=False)
    with (tmp_path / "manifest.csv").open() as stream:
        saved = list(csv.DictReader(stream))
    gr.run_candidates(saved, tmp_path, offline=True, af3=False)
    with pytest.raises(ValueError, match="inputs or configuration changed"):
        gr.run_candidates(saved, tmp_path, offline=False, af3=False)
    (tmp_path / records[0]["fasta"]).write_text(f">P\n{PROTEIN}M\n")
    with pytest.raises(ValueError, match="inputs or configuration changed"):
        gr.run_candidates(saved, tmp_path, offline=True, af3=False)
    assert len(calls) == 2


def test_resume_checks_existing_bundle_record_config_and_input(tmp_path):
    records = rows(tmp_path)
    directory = tmp_path / "jobs" / records[0]["run_name"]
    (directory / "run").mkdir(parents=True)
    (directory / "bundle_manifest.json").write_text('{"status":"complete"}')
    (directory / "dossier.json").write_text("{}")
    record = {"effective_config": {"allow_online": True, "allow_folds": False,
                                   "allow_metalnet": True, "emit_af3_jobs": True},
              "query": {"sha256": "wrong"}}
    (directory / "run" / "run_record.json").write_text(json.dumps(record))
    assert "configuration changed" in gr.predict_one(records[0], out_dir=tmp_path, offline=True)[1]
    assert "protein input changed" in gr.predict_one(records[0], out_dir=tmp_path)[1]


def test_all_failed_run_keeps_global_report_and_af3_failure(tmp_path):
    records = rows(tmp_path)
    gr.write_manifest(records, tmp_path / "manifest.csv")
    gr.summarise(records, {records[0]["run_name"]: "failed (exit 1)"}, out_dir=tmp_path)
    gr.finish_run(tmp_path)
    report = (tmp_path / "RUN_REPORT.html").read_text()
    assert "failed (exit 1)" in report
    assert "no completed TF bundles" in report
    assert "af3/export_warnings.json" in report


def test_report_keeps_failed_resume_status_even_when_old_bundle_is_complete(tmp_path):
    from predictor.report.run_report import write_run_report
    records = rows(tmp_path)
    gr.write_manifest(records, tmp_path / "manifest.csv")
    directory = tmp_path / "jobs" / records[0]["run_name"]
    directory.mkdir(parents=True)
    (directory / "dossier.json").write_text('{"family":"Fur"}')
    (directory / "bundle_manifest.json").write_text('{"status":"complete"}')
    (directory / "REPORT.html").write_text("<html>old bundle</html>")
    gr.summarise(records, {records[0]["run_name"]: "failed (resume configuration changed)"}, out_dir=tmp_path)
    report = write_run_report(tmp_path)
    assert "failed (resume configuration changed)" in report.read_text()


def test_batch_resuming_checks_automatic_sibling_annotation(tmp_path):
    proteins = tmp_path / "provided.faa"
    proteins.write_text(f">P\n{PROTEIN}\n")
    genome = tmp_path / "genome.fna"
    genome.write_text(">c\n" + "ATG" + "AAA" * 25 + "TAA" + "\n")
    annotation = genome.with_suffix(".gff")
    annotation.write_text("c\tsource\tCDS\t1\t81\t.\t+\t0\tprotein_id=P\n")
    run_batch(proteins, tmp_path / "run", genome_fasta=genome, discover_only=True)
    annotation.write_text(annotation.read_text() + "# changed annotation\n")
    with pytest.raises(ValueError, match="different input"):
        run_batch(proteins, tmp_path / "run", genome_fasta=genome, discover_only=True)


def test_genbank_joined_cds_frame_table_and_noncoding_context(tmp_path):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
    from Bio.SeqRecord import SeqRecord
    from predictor.annotate import context, genome_scan
    sequence = "C" * 30 + "AATGT" + "N" * 6 + "GA" + "AAA" * 20 + "TAA" + "C" * 100
    record = SeqRecord(Seq(sequence), id="contig")
    record.annotations["molecule_type"] = "DNA"
    joined = CompoundLocation([SimpleLocation(30, 35, strand=1), SimpleLocation(41, 106, strand=1)])
    record.features = [SeqFeature(joined, type="CDS", qualifiers={"protein_id": ["joined"],
                       "codon_start": ["2"], "transl_table": ["4"]}),
                       SeqFeature(SimpleLocation(140, 165, strand=1), type="tRNA", qualifiers={"gene": ["trnX"]})]
    source = tmp_path / "source.gb"
    SeqIO.write([record], source, "genbank")
    fna, gff = gr.prepare_local(None, None, source, None, tmp_path / "run")
    proteins, kind = genome_scan.proteins_from_input(fna, gff=gff, verbose=False)
    assert kind == "annotated_genome"
    assert proteins == [("joined", "MW" + "K" * 20)]
    assert [line.split("\t")[7] for line in gff.read_text().splitlines() if "\tCDS\t" in line] == ["1", "2"]
    assert any(gene.name == "trnX" for gene in context.from_gff(gff, fna).genes)


def test_af3_invalid_saved_coordinates_skip_row_and_preserve_apo():
    result = build_job_set(PROTEIN, [{"start": "malformed", "end": 20}], None, "Fur")
    assert [job.kind for job in result.jobs] == ["apo"]
    assert any("stored operator coordinates" in row["reason"] for row in result.diagnostics)


def test_af3_sidecars_cannot_write_another_existing_bundle(tmp_path):
    first, second = tmp_path / "first", tmp_path / "second"
    first.mkdir()
    second.mkdir()
    for directory in (first, second):
        (directory / "dossier.json").write_text("{}")
    with pytest.raises(ValueError, match="outside existing bundles"):
        export_bundle(first, second / "nested")
    assert not (second / "nested").exists()


def test_case_colliding_af3_names_are_not_overwritten(tmp_path):
    for directory_name, name in (("one", "TF"), ("two", "tf")):
        directory = tmp_path / "jobs" / directory_name
        directory.mkdir(parents=True)
        (directory / "dossier.json").write_text(json.dumps({
            "tf_id": name, "tf_record": {"sequence": PROTEIN}, "family": "Fur"}))
    result = export_run(tmp_path)
    assert result["local_jobs"] == result["server_jobs"] == 1
    assert any("duplicate job name" in warning for warning in result["warnings"])


def test_cli_missing_local_input_returns_error_without_traceback(monkeypatch, capsys):
    from predictor import api, run_cli
    monkeypatch.setattr(api, "predict", lambda *a, **kw: (_ for _ in ()).throw(FileNotFoundError("missing local genome")))
    assert run_cli.main(["predict", PROTEIN, "--offline", "--local-genome", "missing.fna", "--local-gff", "missing.gff"]) == 2
    assert "missing local genome" in capsys.readouterr().err
