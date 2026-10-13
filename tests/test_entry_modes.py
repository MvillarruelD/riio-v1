"""Entry-mode contracts: exact supplied proteins, local coordinates, and safe resume."""
import csv
import sys
import types
from pathlib import Path

import pytest

from predictor import genome_run as gr
from predictor.batch import run_batch


@pytest.fixture
def reports(monkeypatch):
    calls = []
    module = types.ModuleType("predictor.report.run_report")
    module.write_discovery_report = lambda path: calls.append(("discovery", path))
    module.write_run_report = lambda path: calls.append(("run", path))
    monkeypatch.setitem(sys.modules, module.__name__, module)
    return calls


def test_batch_preserves_every_protein_and_order_without_discovery(tmp_path, reports):
    source = tmp_path / "supplied.faa"
    source.write_text(">second\nMKKLTVSDLAREHGYEVVETPEQVTAL\n>first\nMEQRLAGVSVTEQALNQYVRALQEK\n")
    path = run_batch(source, tmp_path / "run", offline=True, discover_only=True)
    with path.open() as handle:
        rows = list(csv.DictReader(handle))
    assert [Path(r["fasta"]).stem for r in rows] == ["x__second", "x__first"]
    assert all(r["discovery_engine"] == "provided_protein" and r["family"] == "" for r in rows)
    assert [call[0] for call in reports] == ["discovery", "run"]


def test_batch_rejects_duplicate_and_sanitization_collisions(tmp_path, reports):
    source = tmp_path / "supplied.faa"
    source.write_text(">a/b\nMKKLTVSDLAREHGYEVVETPEQVTAL\n>a_b\nMEQRLAGVSVTEQALNQYVRALQEK\n")
    with pytest.raises(ValueError, match="colliding"):
        run_batch(source, tmp_path / "run", discover_only=True)


@pytest.mark.parametrize("identifiers", [("a", "a"), ("Upper", "upper"), ("a/b", "a_b")])
def test_batch_collision_preflight_preserves_existing_fastas(tmp_path, reports, identifiers):
    source = tmp_path / "supplied.faa"
    source.write_text("".join(f">{name}\nMKKLTVSDLAREHGYEVVETPEQVTAL\n" for name in identifiers))
    fastas = tmp_path / "run" / "fastas"
    fastas.mkdir(parents=True)
    previous = fastas / "previous.fasta"
    previous.write_text("preserved\n")
    with pytest.raises(ValueError, match="colliding"):
        run_batch(source, tmp_path / "run", discover_only=True)
    assert list(fastas.iterdir()) == [previous]
    assert previous.read_text() == "preserved\n"


def test_resume_rejects_changed_protein_input(tmp_path, reports):
    source = tmp_path / "supplied.faa"
    source.write_text(">a\nMKKLTVSDLAREHGYEVVETPEQVTAL\n")
    run_batch(source, tmp_path / "run", discover_only=True)
    source.write_text(">b\nMEQRLAGVSVTEQALNQYVRALQEK\n")
    with pytest.raises(ValueError, match="different input"):
        run_batch(source, tmp_path / "run", discover_only=True)


def test_accession_table_uses_supplied_order(tmp_path, reports, monkeypatch):
    source = tmp_path / "supplied.tsv"
    source.write_text("accession\nWP_000240047.1\nWP_000240048.1\n")
    seen = []
    def resolve(accession, **kw):
        seen.append(accession)
        return "MKKLTVSDLAREHGYEVVETPEQVTAL"
    monkeypatch.setattr("predictor.run_cli._resolve_prediction_source", resolve)
    run_batch(source, tmp_path / "run", discover_only=True)
    assert seen == ["WP_000240047.1", "WP_000240048.1"]


def test_local_prediction_subprocess_receives_exact_contig_and_half_open_locus(tmp_path, monkeypatch):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    (inputs / "genome.fna").write_text(">second\nATG\n")
    (inputs / "genome.gff").write_text("##gff-version 3\n")
    seen = []
    monkeypatch.setattr(gr.subprocess, "run", lambda cmd, **kw: seen.append(cmd) or types.SimpleNamespace(returncode=1))
    gr.predict_one({"run_name": "my_tf", "fasta": "fastas/query.fasta", "contig": "second",
                    "start": "11", "end": "100", "strand": "-", "organism_acc": ""}, out_dir=tmp_path)
    command = seen[0]
    assert command[command.index("--tf-locus")+1:] == ["10", "100", "-"]
    assert command[command.index("--local-contig")+1] == "second"
    assert "--family" not in command and "--organism" not in command


def test_multirecord_genbank_preserves_contigs_and_cds_ids(tmp_path):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, SimpleLocation
    records = []
    for name in ("chromosome", "plasmid"):
        record = SeqRecord(Seq("ATG" + "AAA"*25 + "TAA"), id=name)
        record.annotations["molecule_type"] = "DNA"
        record.features = [SeqFeature(SimpleLocation(0, len(record), strand=1), type="CDS",
                           qualifiers={"protein_id": [name+"_protein"], "locus_tag": [name+"_locus"]})]
        records.append(record)
    source = tmp_path / "input.gb"
    SeqIO.write(records, source, "genbank")
    fna, gff = gr.prepare_local(None, None, source, None, tmp_path / "run")
    assert [r.id for r in SeqIO.parse(fna, "fasta")] == ["chromosome", "plasmid"]
    assert set(gr.gff_index(gff)) == {"chromosome_protein", "plasmid_protein"}


def test_local_genome_runner_keeps_discovery_coordinates(tmp_path, reports, monkeypatch):
    source = tmp_path / "g.fna"
    source.write_text(">contig\n" + "ATG" + "AAA"*25 + "TAA" + "\n")
    gff = tmp_path / "g.gff"
    gff.write_text("contig\tRefSeq\tCDS\t1\t81\t.\t+\t0\tID=cds;protein_id=P;gene=testR\n")
    candidate = types.SimpleNamespace(protein_id="P", sequence="M"+"K"*25, family="MerR", profile="test")
    monkeypatch.setattr(gr, "discover", lambda *a, **kw: ([candidate], "bitacora", "annotated_genome"))
    path = gr.run_genome(None, tmp_path / "run", fasta=source, gff=gff, discover_only=True)
    with path.open() as handle:
        row = next(csv.DictReader(handle))
    assert (row["contig"], row["start"], row["end"]) == ("contig", "1", "81")
    assert (tmp_path / "run" / "inputs" / "genome.gff").is_file()


def test_cli_exposes_new_modes_and_local_options(monkeypatch):
    from predictor import run_cli
    calls = []
    monkeypatch.setattr(run_cli, "cmd_genome", lambda a: calls.append(a) or 0)
    assert run_cli.main(["genome", "--fasta", "g.fna", "--gff", "g.gff", "--out", "run"]) == 0
    assert calls[0].accession is None and calls[0].fasta == "g.fna"
