"""genome_run: discovery rows, resume, and the summary table -- offline, no engine, no genome."""
from __future__ import annotations

import csv
import json
from dataclasses import dataclass

from predictor import genome_run as gr


@dataclass
class _Cand:
    protein_id: str
    sequence: str
    family: str = "MerR"
    profile: str = "ecZntR"
    trim_start: int | None = 0
    trim_end: int | None = 20


def _rows(tmp_path):
    cands = [_Cand("WP_1.1", "M" * 30), _Cand("WP_2.1", "K" * 25, family=""),   # no family -> skipped
             _Cand("WP_3.1", "A" * 40, family="ArsR/SmtB", trim_start=None, trim_end=None)]
    coords = {"WP_1.1": {"contig": "c1", "start": 10, "end": 100, "strand": "+", "locus_tag": "b0001",
                         "gene": "zntR", "product": "Zn regulator"}}
    return gr.manifest_rows(cands, short="G", organism="G", accession="G", engine="bitacora",
                            protein_source="annotated_genome", coords=coords,
                            fasta_dir=tmp_path / "fastas", rel_to=tmp_path)


def test_manifest_rows_write_full_length_fastas_and_keep_provenance(tmp_path):
    rows = _rows(tmp_path)
    assert [r["run_name"] for r in rows] == ["G__MerR__WP_1.1", "G__ArsR_SmtB__WP_3.1"]
    assert rows[0]["locus_tag"] == "b0001" and rows[0]["family_flag"] == "(auto)"
    assert rows[0]["discovery_engine"] == "bitacora" and rows[1]["trim_start"] == ""
    text = (tmp_path / rows[0]["fasta"]).read_text()
    assert text.splitlines()[1] == "M" * 30, "the FASTA must hold the full-length sequence"
    gr.write_manifest(rows, tmp_path / "manifest.csv")
    with (tmp_path / "manifest.csv").open(encoding="utf-8") as fh:
        assert list(csv.DictReader(fh).fieldnames) == gr.MANIFEST_COLS + gr.PROV_COLS


def test_complete_bundles_are_resumed_not_rerun(tmp_path, monkeypatch):
    rows = _rows(tmp_path)
    b = tmp_path / "jobs" / rows[0]["run_name"]
    b.mkdir(parents=True)
    (b / "bundle_manifest.json").write_text(json.dumps({"status": "complete"}))
    (b / "dossier.json").write_text(json.dumps({"family": "MerR", "inducers": {"top": "Zn2+", "agreement": 1.0},
                                                "operators": {"primary": {"sequence": "ACGT",
                                                              "provenance": {"start": 5, "end": 8,
                                                                             "locality": "proximal"}}},
                                                "regulon": [{}, {}]}))
    ran = []
    monkeypatch.setattr(gr.subprocess, "run", lambda *a, **k: ran.append(a) or type("R", (), {"returncode": 1})())
    assert gr.predict_one(rows[0], out_dir=tmp_path) == (rows[0]["run_name"], "done (resumed)")
    assert not ran, "a complete bundle must not be predicted again"
    assert gr.predict_one(rows[1], out_dir=tmp_path)[1].startswith("failed")


def test_summary_reads_headline_calls_back_from_the_bundle(tmp_path):
    rows = _rows(tmp_path)
    b = tmp_path / "jobs" / rows[0]["run_name"]
    b.mkdir(parents=True)
    (b / "bundle_manifest.json").write_text(json.dumps({"status": "complete"}))
    (b / "dossier.json").write_text(json.dumps({"family": "MerR", "ssn_cluster": "MerR_9",
                                                "inducers": {"top": "Zn2+", "agreement": 1.0},
                                                "operators": {"primary": {"sequence": "ACGT",
                                                              "provenance": {"start": 5, "end": 8,
                                                                             "locality": "proximal"}}},
                                                "regulon": [{}, {}]}))
    path = gr.summarise(rows, {rows[0]["run_name"]: "ok"}, out_dir=tmp_path)
    with path.open(encoding="utf-8") as fh:
        got = list(csv.DictReader(fh, delimiter="\t"))
    assert got[0]["inducer_top"] == "Zn2+" and got[0]["inducer_class"] == "metal"
    assert got[0]["primary_locality"] == "proximal" and got[0]["n_regulon_operons"] == "2"
    assert got[1]["status"] == "not run" and got[1]["inducer_top"] == ""


def test_production_candidates_are_located_from_the_genomes_own_annotation(tmp_path):
    """BITACORA reports no coordinates; the GUI table reads them back from the annotation it searched."""
    cds = "ATG" + "AAA" * 8 + "TAA"                       # M K8 stop
    (tmp_path / "g.fna").write_text(">c1\n" + "CC" + cds + "GGGG\n")
    (tmp_path / "g.gff").write_text("##gff-version 3\n" + "\t".join(
        ["c1", "RefSeq", "CDS", "3", str(2 + len(cds)), ".", "+", "0",
         "ID=cds1;protein_id=WP_9.1;locus_tag=b9999;gene=tstR;product=test%2C regulator"]) + "\n")
    cand = type("C", (), {"protein_id": "WP_9.1", "sequence": "MKKKKKKKK", "family": "MerR",
                          "evidence": {"evalue": "1e-20", "bits": "80.5"}})()
    rows, n = gr.as_census_candidates([cand], tmp_path / "g.fna", tmp_path / "g.gff")
    assert n == 1
    r = rows[0]
    assert (r.contig, r.start, r.end, r.strand, r.locus_tag) == ("c1", 3, 2 + len(cds), "+", "b9999")
    assert r.evalue == 1e-20 and r.bits == 80.5 and r.pfam_acc is None and r.family == "MerR"


def test_resuming_a_directory_written_for_another_genome_is_refused(tmp_path):
    import pytest
    gr.write_manifest([{"run_name": "A__ArsR__X", "fasta": "fastas/x.fasta", "organism_acc": "GCF_000000001.1"}],
                      tmp_path / "manifest.csv")
    with pytest.raises(RuntimeError, match="holds a run for GCF_000000001.1"):
        gr.run_genome("GCF_000005845.2", tmp_path, discover_only=True)
    # the same accession resumes as before
    assert gr.run_genome("GCF_000000001.1", tmp_path, discover_only=True) == tmp_path / "manifest.csv"


def test_predict_all_reports_in_completion_order_but_returns_manifest_order(tmp_path, monkeypatch, capsys):
    import time
    delays = {"slow": 0.3, "fast": 0.0}

    def fake(row, *, out_dir, **kw):
        time.sleep(delays[row["run_name"]])
        return row["run_name"], "ok"

    monkeypatch.setattr(gr, "predict_one", fake)
    status = gr.predict_all([{"run_name": "slow"}, {"run_name": "fast"}], out_dir=tmp_path, jobs=2)
    lines = [ln for ln in capsys.readouterr().out.splitlines() if ln.strip().startswith("[")]
    assert "fast" in lines[0] and "slow" in lines[1]
    assert list(status) == ["slow", "fast"]
