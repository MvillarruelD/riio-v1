"""Regressions from the 2026-11-22 end-to-end review: one test (or group) per defect found.

Offline and engine-free: every network or engine boundary is replaced.
"""
from __future__ import annotations

import csv
import io
import json
import os
import random
import shlex
import socket
import subprocess
import types
import urllib.error

import pytest

PROTEIN_A = "MKKLTVSDLAREHGYEVVETPEQVTAL"
PROTEIN_B = "MEQRLAGVSVTEQALNQYVRALQEKWW"


def _dna(n, seed=0):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def _rc(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _cds(contig, start, end, strand, attrs, phase="0"):
    return f"{contig}\tT\tCDS\t{start}\t{end}\t.\t{strand}\t{phase}\t{attrs}"


# ----------------------------------------------------------------------------- caches: identity
def test_auto_locus_cache_is_the_protein_not_the_job_name(tmp_path, monkeypatch):
    from predictor.annotate import auto_locus, genome_resolver
    monkeypatch.setattr(auto_locus, "CACHE", tmp_path)
    monkeypatch.setattr(auto_locus, "top_protein_hits", lambda seq, **kw: [seq])

    class Entrez:
        @staticmethod
        def efetch(id, **kw):
            genome = "NC_000001.1" if id == PROTEIN_A else "NC_000002.1"
            return io.StringIO("Id\tSource\tNucleotide Accession\tStart\tStop\tStrand\tProtein\tName\tOrganism\n"
                               f"1\tRefSeq\t{genome}\t101\t400\t+\tP\tTF\torg\n")

    monkeypatch.setattr(genome_resolver, "_entrez", lambda: Entrez)
    assert auto_locus.auto_locus(PROTEIN_A, key="query", verbose=False)["accession"] == "NC_000001.1"
    assert auto_locus.auto_locus(PROTEIN_B, key="query", verbose=False)["accession"] == "NC_000002.1"


def test_alignment_file_is_shared_only_by_the_same_sequence(tmp_path, monkeypatch):
    from predictor.annotate import msa_homologs
    monkeypatch.setattr(msa_homologs, "CACHE", tmp_path)
    assert msa_homologs.alignment_path("query", PROTEIN_A) != msa_homologs.alignment_path("query", PROTEIN_B)
    legacy = tmp_path / "query_colabfold_uniref.a3m"
    legacy.write_text(f">101\n{PROTEIN_A}\n>hit\n{PROTEIN_B}\n")
    assert msa_homologs.alignment_path("query", PROTEIN_A) == legacy, "existing caches stay valid"
    assert msa_homologs.alignment_path("query", PROTEIN_B) != legacy


def test_fold_cache_is_the_sequence_not_the_job_name(tmp_path, monkeypatch):
    from predictor.signals import structure_backend
    tags = []
    adapter = types.SimpleNamespace(write_af3_layout=lambda result, tf_id, **kw: tags.append(tf_id) or
                                    types.SimpleNamespace(cif_path="x", folder="x", plddt_mean=1, n_tokens=1))
    inputs = types.SimpleNamespace(build_apo_homomer_input=lambda *a, **k: types.SimpleNamespace(
        structure_input=None, chain_ids=["A", "B"]))
    backend = types.SimpleNamespace(fold=lambda *a, **k: types.SimpleNamespace(ptm=1, iptm=1))
    monkeypatch.setattr(structure_backend, "_lib", lambda: (inputs, backend, adapter, None))
    monkeypatch.setattr(structure_backend, "WORK", tmp_path)
    for seq in (PROTEIN_A, PROTEIN_B):
        structure_backend.fold_apo_homomer("query", seq)
    assert tags[0] != tags[1]


def test_bitacora_default_work_folder_is_the_genome_content(tmp_path):
    from predictor.discovery import bitacora
    keys = []
    for name in ("one", "two"):
        path = tmp_path / name / "genome.fna"
        path.parent.mkdir()
        path.write_text(f">c\n{_dna(60, len(keys))}\n")
        keys.append(bitacora._work_key(path))
    assert keys[0] != keys[1]


# ----------------------------------------------------------------------------- caches: failures
def test_failed_protein_fetch_is_retried_and_old_empty_entries_heal(tmp_path, monkeypatch):
    from predictor.annotate import genome_resolver, protein_seqs
    monkeypatch.setattr(protein_seqs, "_CACHE", tmp_path)

    class Down:
        @staticmethod
        def efetch(**kw):
            raise OSError("network down")

    class Up:
        @staticmethod
        def efetch(**kw):
            return io.StringIO(f">WP_000001.1 x\n{PROTEIN_A}\n")

    monkeypatch.setattr(genome_resolver, "_entrez", lambda: Down)
    assert protein_seqs.fetch_protein_sequences(["WP_000001.1"], allow_ncbi=True) == {}
    assert not list(tmp_path.glob("*.faa"))
    (tmp_path / "WP_000001.faa").write_text(">WP_000001\n\n")  # written by earlier versions
    monkeypatch.setattr(genome_resolver, "_entrez", lambda: Up)
    assert protein_seqs.fetch_protein_sequences(["WP_000001.1"], allow_ncbi=True) == {"WP_000001": PROTEIN_A}


def test_uniprot_accessions_resolve_from_ncbi_fasta_headers():
    from predictor.annotate import protein_seqs
    assert protein_seqs._parse_fasta(f">sp|P0A9E0.1|ARAC_ECOLI x\n{PROTEIN_A}\n")["P0A9E0"] == PROTEIN_A


@pytest.mark.parametrize("error,cached", [(urllib.error.URLError("down"), False),
                                          (urllib.error.HTTPError("u", 404, "no", {}, None), True)])
def test_clade_promoters_are_cached_only_when_complete(tmp_path, monkeypatch, error, cached):
    from predictor.annotate import msa_homologs
    monkeypatch.setattr(msa_homologs, "CACHE", tmp_path)
    monkeypatch.setattr(msa_homologs, "_cluster_by_id", lambda cid: object())
    monkeypatch.setattr(msa_homologs, "_sample_members", lambda *a, **k: ["P12345"])
    monkeypatch.setattr(msa_homologs.time, "sleep", lambda s: None)

    def fail(*a, **k):
        raise error

    monkeypatch.setattr(msa_homologs, "_get", fail)
    record = msa_homologs.cluster_promoter_regions(PROTEIN_A, "Fur_x", verbose=False)
    assert record["regions"] == []
    assert bool(list(tmp_path.glob("*.json"))) is cached
    assert ("incomplete" in record) is not cached


@pytest.mark.parametrize("error,cached", [(urllib.error.URLError("down"), False),
                                          (urllib.error.HTTPError("u", 404, "no", {}, None), True)])
def test_afdb_outage_is_not_remembered_as_absence(tmp_path, monkeypatch, error, cached):
    from predictor.structure import afdb
    monkeypatch.setattr(afdb, "_CACHE", tmp_path)

    def fail(*a, **k):
        raise error

    monkeypatch.setattr(afdb.urllib.request, "urlopen", fail)
    assert afdb.lookup("P0ACS5") is None
    assert bool(list(tmp_path.glob("*.json"))) is cached


def test_metalnet_engine_failure_is_not_cached(tmp_path, monkeypatch):
    from predictor.structure import metalnet
    monkeypatch.setattr(metalnet, "CACHE", tmp_path)
    monkeypatch.setattr(metalnet, "status", lambda: {"ready": True, "python": "python"})
    monkeypatch.setattr(metalnet, "a3m_for", lambda *a, **k: (tmp_path / "x.a3m", "test"))

    def engine(cmd, **kw):
        out = cmd[cmd.index("--out") + 1]
        sid = f"mn_{metalnet._key(PROTEIN_A)}"
        with open(out, "w") as handle:
            json.dump({"results": {sid: {"ok": False, "reason": "out of memory"}}}, handle)
        return types.SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(metalnet.subprocess, "run", engine)
    assert metalnet.predict_site(PROTEIN_A)["has_site"] is None
    assert not list(tmp_path.glob("*.json"))


def test_atomic_publication_leaves_the_old_file_on_failure(tmp_path):
    from predictor import resources
    target = tmp_path / "cache.json"
    target.write_text("old")
    with pytest.raises(RuntimeError):
        with resources.atomic_path(target) as temporary:
            temporary.write_text("half")
            raise RuntimeError
    assert target.read_text() == "old" and [p.name for p in tmp_path.iterdir()] == ["cache.json"]


# ----------------------------------------------------------------------------- genome identity + annotation
def test_mirror_never_matches_a_partial_accession(tmp_path):
    from predictor.annotate import genome_mirror
    for stem in ("GCF_000005845.2", "GCF_000006745.1_ASM674v1"):
        (tmp_path / f"{stem}.fna").write_text(">c\nACGT\n")

    def find(accession):
        return genome_mirror.ensure_genome(accession, allow_download=False, mirror_dir=tmp_path)[0]

    for wrong in ("GCF_0", "GCF_00000584", "GCF_000005845.", "*", "../GCF_000005845.2"):
        assert find(wrong) is None, wrong
    assert find("GCF_000005845.2").stem == "GCF_000005845.2"
    assert find("GCF_000005845").stem == "GCF_000005845.2"
    assert find("GCF_000006745.1").stem == "GCF_000006745.1_ASM674v1"


def test_identical_protein_copies_are_not_joined(tmp_path):
    from predictor import genome_run
    from predictor.annotate import genome_scan
    cds = "ATG" + "AAAGAACTG" * 9 + "TAA"
    genome = "CC" * 50 + _rc(cds) + "GG" * 500 + cds + "TT" * 50
    second = 100 + len(cds) + 1000
    fna, gff = tmp_path / "g.fna", tmp_path / "g.gff"
    fna.write_text(f">c\n{genome}\n")
    gff.write_text("\n".join([
        _cds("c", 101, 100 + len(cds), "-", "ID=cds-WP_1;protein_id=WP_1.1"),
        _cds("c", second + 1, second + len(cds), "+", "ID=cds-WP_1-2;protein_id=WP_1.1")]) + "\n")
    proteins, meta = genome_scan.proteins_from_gff(fna, gff, return_metadata=True)
    assert proteins == [("WP_1.1", "M" + "KEL" * 9)]
    assert (meta["WP_1.1"]["start"], meta["WP_1.1"]["end"]) == (101, 100 + len(cds))
    locus = genome_run.gff_index(gff)["WP_1.1"]
    assert (locus["start"], locus["end"]) == (101, 100 + len(cds))


def test_segments_of_one_feature_are_still_joined(tmp_path):
    from predictor.annotate import genome_scan
    cds = "ATG" + "AAAGAACTG" * 9 + "TAA"
    fna, gff = tmp_path / "g.fna", tmp_path / "g.gff"
    fna.write_text(f">c\n{cds[:30]}{'C' * 20}{cds[30:]}\n")
    gff.write_text("\n".join([_cds("c", 1, 30, "+", "ID=cds-1;protein_id=P1"),
                              _cds("c", 51, 20 + len(cds), "+", "ID=cds-1;protein_id=P1")]) + "\n")
    assert genome_scan.proteins_from_gff(fna, gff) == [("P1", "M" + "KEL" * 9)]


@pytest.mark.parametrize("strand", ["+", "-"])
def test_partial_cds_respects_gff_phase(tmp_path, strand):
    from predictor.annotate import genome_scan
    cds = "GG" + "AAAGAACTG" * 9
    fna, gff = tmp_path / "g.fna", tmp_path / "g.gff"
    fna.write_text(f">c\n{cds if strand == '+' else _rc(cds)}\n")
    gff.write_text(_cds("c", 1, len(cds), strand, "ID=x;protein_id=P1", phase="2") + "\n")
    assert genome_scan.proteins_from_gff(fna, gff) == [("P1", "KEL" * 9)]


def test_product_names_are_decoded_as_utf8(tmp_path):
    from predictor import genome_run
    gff = tmp_path / "g.gff"
    gff.write_text(_cds("c", 1, 90, "+", "ID=x;protein_id=P1;product=r%C3%A9gulateur") + "\n")
    assert genome_run.gff_index(gff)["P1"]["product"] == "r\u00e9gulateur"


def _genbank(path, second_translation):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, SimpleLocation
    from Bio.SeqRecord import SeqRecord
    record = SeqRecord(Seq(_dna(600)), id="chr1")
    record.annotations["molecule_type"] = "DNA"
    for start, translation in ((0, PROTEIN_A), (300, second_translation)):
        record.features.append(SeqFeature(SimpleLocation(start, start + 90, strand=1), type="CDS",
                                          qualifiers={"protein_id": ["WP_1.1"], "translation": [translation]}))
    SeqIO.write(record, str(path), "genbank")


def test_refseq_genbank_with_paralogous_copies_is_accepted(tmp_path):
    from predictor import genome_run
    _genbank(tmp_path / "ok.gbff", PROTEIN_A)
    fna, _ = genome_run.prepare_local(None, None, tmp_path / "ok.gbff", None, tmp_path / "ok")
    assert fna.with_suffix(".proteins.faa").read_text().count(">") == 1
    _genbank(tmp_path / "bad.gbff", PROTEIN_B)
    with pytest.raises(ValueError, match="must be unique"):
        genome_run.prepare_local(None, None, tmp_path / "bad.gbff", None, tmp_path / "bad")


def test_soft_masked_genomes_read_like_unmasked_ones():
    from predictor.annotate import context
    gff = "##gff-version 3\n" + _cds("c", 11, 40, "+", "ID=x;protein_id=P1")
    assert context.from_gff(gff, ">c\nacgtACGTacgt\n").sequence == "ACGTACGTACGT"


def test_circularity_comes_from_the_replicon_that_was_loaded():
    from predictor.annotate import context
    gff = "##gff-version 3\nchr\tT\tregion\t1\t60\t.\t+\t.\tIs_circular=true\n"
    fasta = f">chr\n{_dna(60)}\n>linear\n{_dna(40)}\n"
    assert context.from_gff(gff, fasta, accession="linear").circular is False
    assert context.from_gff(gff, fasta, accession="chr").circular is True
    assert context.from_gff(gff, fasta).circular is True


def test_local_assembly_is_scanned_whole_like_an_accession(tmp_path, monkeypatch):
    from predictor import pipeline
    from predictor.config import RunConfig
    fna, gff = tmp_path / "g.fna", tmp_path / "g.gff"
    fna.write_text(f">chr\n{_dna(3000, 1)}\n>plasmid\n{_dna(1500, 2)}\n")
    gff.write_text("##gff-version 3\n" + _cds("plasmid", 601, 900, "+", "ID=x;protein_id=P1") + "\n")

    class Reached(Exception):
        pass

    def stop(ctx, start, end, strand):
        raise Reached(ctx, start, end)

    monkeypatch.setattr(pipeline.context, "locate_operator_window", stop)
    with pytest.raises(Reached) as caught:
        pipeline.run_novel(PROTEIN_A, name="local_tf", family="Fur", local_genome=fna, local_gff=gff,
                           local_contig="plasmid", tf_locus=(600, 900, "+"),
                           cfg=RunConfig(allow_online=False, allow_folds=False, write=False, verbose=False))
    ctx, start, end = caught.value.args
    offset = ctx.contig_offsets["plasmid"]
    assert set(ctx.contig_offsets) == {"chr", "plasmid"} and (start, end) == (600 + offset, 900 + offset)
    assert ctx.locate(start) == ("plasmid", 600)


def test_bitacora_hit_ids_resolve_to_the_longest_identifier(tmp_path):
    from predictor._engines import bitacora_runner
    (tmp_path / "X_proteins_trimmed.fasta").write_text(f">contig_12_dom\n{PROTEIN_B}\n")
    rows = bitacora_runner.collect(tmp_path, {"contig_1": PROTEIN_A, "contig_12": PROTEIN_B})
    assert [row["protein_id"] for row in rows] == ["contig_12"]


def _engine_call(tmp_path, monkeypatch, folder):
    from predictor._engines import bitacora_runner
    seen = []

    def engine(cmd, **kw):
        seen.append((cmd, kw.get("cwd")))
        if os.name != "nt":                    # upstream: `makeblastdb -in <proteins>` writes beside it
            with open(cmd[cmd.index("-p") + 1] + ".phr", "w") as index:
                index.write("index")
        return types.SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(bitacora_runner.subprocess, "run", engine)
    out = tmp_path / folder
    out.mkdir()
    bitacora_runner.run_bitacora(engine=tmp_path, profiles=tmp_path, proteins=tmp_path / "p.faa", genome=None,
                                 gff=None, out=out, name="my genome;x", threads=1, evalue="1e-3", timeout=5)
    return out, *seen[0]


@pytest.mark.parametrize("folder", ["my run", "R&D (2) $(touch pwned)", "a;b"])
def test_engine_never_receives_a_path_it_would_split_or_execute(tmp_path, monkeypatch, folder):
    from predictor._engines import bitacora_runner
    out, cmd, cwd = _engine_call(tmp_path, monkeypatch, folder)
    if os.name == "nt":
        launched = cmd[-1].rsplit(" && ", 1)[-1]
        assert launched.startswith('bash "$stage/engine/runBITACORA_command_line.sh"')
        assert folder not in launched and "my_genome_x" in launched
        assert "rm -r" not in cmd[-1]
    else:
        assert all(bitacora_runner.engine_safe(word) for word in cmd) and bitacora_runner.engine_safe(cwd)
        assert "my_genome_x" in cmd
        assert (tmp_path / "p.faa.phr").is_file(), "the index lands beside the real proteome"
        assert not os.path.lexists(os.path.dirname(cwd)), "the link folder is gone"
    assert out.is_dir(), "removing the links must leave their targets alone"


def test_engine_gets_safe_paths_unchanged(tmp_path, monkeypatch):
    from predictor._engines import bitacora_runner
    if not bitacora_runner.engine_safe(bitacora_runner.win2wsl(tmp_path)):
        pytest.skip("the temporary folder itself needs links")
    out, cmd, cwd = _engine_call(tmp_path, monkeypatch, "plain")
    assert "stage" not in " ".join(cmd)
    if os.name != "nt":
        assert cwd == str(out)


# ----------------------------------------------------------------------------- speedups are exact
def test_flank_index_matches_the_linear_scan():
    from predictor.signals import motif_rescan
    rng = random.Random(4)
    for _ in range(300):
        genes = [motif_rescan.Gene(s, s + rng.randint(1, 30), rng.choice("+-"))
                 for s in (rng.randint(0, 120) for _ in range(rng.randint(0, 14)))]
        index = motif_rescan._FlankIndex(genes)
        a = rng.randint(0, 160)
        gap = (a, a + rng.randint(0, 15))
        left, right = motif_rescan._flanking(genes, gap)
        assert index(gap)[0] is left and index(gap)[1] is right


def test_vectorised_pvalues_match_the_scalar_definition():
    import numpy as np
    from predictor.motifs import pwm_scan
    rng = np.random.default_rng(3)
    for width in (6, 17, 33):
        prob = rng.dirichlet([.4] * 4, size=width).T
        bg = np.array([.3, .2, .2, .3])
        pssm = pwm_scan.log_odds(prob, bg)
        offset, dist, sf = pwm_scan._cached_null(pssm.tobytes(), width, bg.tobytes(), 50.0)
        scalar = pwm_scan._pvalue_fn(offset, dist, 50.0)
        scores = np.concatenate([rng.uniform(pssm.min(0).sum() - 1, pssm.max(0).sum() + 1, 500),
                                 [pssm.min(0).sum(), pssm.max(0).sum()]])
        assert pwm_scan._pvalues(scores, offset, dist, sf, 50.0).tolist() == [scalar(s) for s in scores]


# ----------------------------------------------------------------------------- input, API, runner
def test_deep_run_folders_can_be_read(tmp_path):
    from predictor import input_validation
    folder = tmp_path
    while len(str(folder)) < 250:
        folder = folder / ("d" * 30)
    try:
        folder.mkdir(parents=True)
        (folder / "query.fasta").write_text(f">q\n{PROTEIN_A}\n")
    except OSError:
        pytest.skip("long paths are disabled on this system")
    assert input_validation.read_protein(str(folder / "query.fasta")) == PROTEIN_A


@pytest.mark.parametrize("value", ["protein.faa", "runs/x/query.fasta", r"C:\data\my tf.fasta"])
def test_a_mistyped_path_says_file_not_found(value):
    from predictor import input_validation
    with pytest.raises(ValueError, match="file not found"):
        input_validation.read_protein(value)


def test_predict_does_not_modify_the_callers_config(monkeypatch):
    from predictor import api, pipeline
    from predictor.config import RunConfig
    monkeypatch.setattr(pipeline, "run_novel", lambda *a, **k: {"ok": False, "error": "stop"})
    cfg = RunConfig()
    api.predict(PROTEIN_A, cfg=cfg, verbose=False)
    assert cfg.allow_folds is True and cfg.verbose is True


def _run(tmp_path, names=("G__Fur__one", "G__Fur__two"), complete=("G__Fur__one",)):
    from predictor import genome_run
    rows = []
    for name in names:
        row = {key: "" for key in genome_run.MANIFEST_COLS + genome_run.PROV_COLS}
        row.update(run_name=name, fasta=f"fastas/{name}.fasta", family="Fur")
        rows.append(row)
        if name in complete:
            bundle = tmp_path / "jobs" / name
            bundle.mkdir(parents=True)
            (bundle / "bundle_manifest.json").write_text('{"status": "complete"}')
            (bundle / "dossier.json").write_text(json.dumps({
                "family": "Fur", "operators": {"primary": {"sequence": "ACGT"}}, "regulon": [],
                "inducers": {"top": None, "finalized": True, "sensor_class": None,
                             "identity_status": "unresolved", "candidates": ["Zn2+", "Cd2+"]}}))
    genome_run.write_manifest(rows, tmp_path / "manifest.csv")
    return rows


def test_a_second_partial_selection_keeps_earlier_results(tmp_path):
    from predictor import genome_run
    rows = _run(tmp_path)
    genome_run.summarise(rows, {"G__Fur__two": "failed (exit 1)"}, out_dir=tmp_path)
    with (tmp_path / "summary.tsv").open() as handle:
        status = {r["run_name"]: r["status"] for r in csv.DictReader(handle, delimiter="\t")}
    assert status == {"G__Fur__one": "done (resumed)", "G__Fur__two": "failed (exit 1)"}


def test_control_characters_do_not_lose_the_workbook(tmp_path):
    from predictor.report import run_report
    run_report._write_xlsx(tmp_path / "t.xlsx", [{"product": "bad\x07name"}], ["product", "report_path"])
    assert (tmp_path / "t.xlsx").is_file()


@pytest.mark.parametrize("message,skipped", [("AssertionError: motif not found", False),
                                             ("NameError: name 'blast_dir' is not defined", False),
                                             ("ModuleNotFoundError: No module named 'predictor.gone'", False),
                                             ("FileNotFoundError: BLAST+ tool 'blastn' not found", True)])
def test_selftest_does_not_hide_code_defects_as_missing_engines(message, skipped):
    from predictor import run_cli
    assert run_cli._is_environmental(f"Traceback...\n{message}\n") is skipped


def test_a_selftest_timeout_is_not_a_pass(monkeypatch):
    from predictor import run_cli

    def slow(cmd, **kw):
        raise subprocess.TimeoutExpired(cmd, 1)

    monkeypatch.setattr(subprocess, "run", slow)
    assert run_cli.cmd_selftest(types.SimpleNamespace(only="input_validation,config", timeout=1, verbose=False)) == 1


def test_scan_rejects_an_unknown_family(tmp_path, capsys):
    from predictor import run_cli
    source = tmp_path / "p.faa"
    source.write_text(f">a\n{PROTEIN_A}\n")
    assert run_cli.main(["scan", str(source), "--families", "Furr", "--out", str(tmp_path / "o")]) == 2
    assert "unknown family" in capsys.readouterr().err


# ----------------------------------------------------------------------------- GUI
def test_gui_counts_resumed_and_unsummarised_bundles(tmp_path):
    from predictor import genome_run, gui_support
    run = tmp_path / "run"
    run.mkdir()
    rows = _run(run)
    assert gui_support.run_rows(run, rows)[0]["status"] == "done (resumed)", "an interrupted run has no summary"
    genome_run.summarise(rows, {"G__Fur__one": "done (resumed)"}, out_dir=run)
    shown = gui_support.run_rows(run, rows)
    assert [gui_support.is_done(r["status"]) for r in shown] == [True, False]
    assert not any(gui_support.is_failed(r["status"]) for r in shown)
    assert gui_support.list_saved_runs(tmp_path)[0]["done"] == 1


def test_gui_ledger_lists_the_members_of_a_shortlist(tmp_path):
    from predictor import gui_support
    row = gui_support.run_rows(tmp_path, _run(tmp_path))[0]
    assert "Zn2+" in row["inducer_label"] and "Cd2+" in row["inducer_label"]
    assert row["class"] == "unresolved"


@pytest.mark.parametrize("name", ["..", "../../secrets", "a\\b", "C:x"])
def test_gui_never_follows_a_table_name_out_of_the_run(tmp_path, name):
    from predictor import gui_support
    assert gui_support.job_dir(tmp_path, name) is None


def test_gui_prints_no_command_it_cannot_rebuild(tmp_path):
    from predictor import gui_support
    assert gui_support.run_command(tmp_path) is None
    for info in ({"mode": "protein_list", "files": [{"role": "proteins"}]},
                 {"mode": "local_genome", "files": [{"role": "fasta"}]},
                 {"mode": "genome_accession", "accession": None}):
        (tmp_path / "run_input.json").write_text(json.dumps(info))
        assert gui_support.run_command(tmp_path) is None, info


@pytest.mark.parametrize("protein_list", [True, False])
def test_gui_printed_command_passes_the_runners_fingerprint(tmp_path, monkeypatch, protein_list):
    """Discover as the page does, then run the printed command through the real parser and identity check."""
    from predictor import api, genome_run, gui_support, run_cli
    from predictor.report import run_report
    monkeypatch.setattr(run_report, "write_discovery_report", lambda out: None)
    monkeypatch.setattr(run_report, "write_run_report", lambda out: None)
    out = tmp_path / "my run"
    (out / "uploads").mkdir(parents=True)
    source = out / "uploads" / ("input.faa" if protein_list else "input.fna")
    source.write_text(f">a\n{PROTEIN_A}\n")
    common = dict(out_dir=out, discover_only=True, organism="GCF_000005845.2")
    if protein_list:
        api.predict_many(source=source, **common)
    else:
        source.write_text(f">c\n{_dna(300)}\n")
        monkeypatch.setattr(genome_run, "prepare_local", lambda *a: (source, None))
        monkeypatch.setattr(genome_run, "discover", lambda *a, **k: (
            [types.SimpleNamespace(protein_id="a", sequence=PROTEIN_A, family="Fur")], "bitacora", "genome"))
        from predictor.annotate import genome_scan
        monkeypatch.setattr(genome_scan, "proteins_from_input", lambda *a, **k: ([], "genome", {}))
        api.run_genome(fasta=source, census_fallback=False, census_only=False, **common)
    command = gui_support.run_command(out, windows=False)
    reached = []
    monkeypatch.setattr(genome_run, "run_candidates", lambda *a, **k: reached.append(a) or out / "manifest.csv")
    assert run_cli.main(shlex.split(command)[1:]) == 0 and reached, command


@pytest.mark.skipif(os.name != "nt", reason="PowerShell quoting is a Windows question")
@pytest.mark.parametrize("word", [r"C:\R&D\run (2)", "a,b", "Cu+/Ag+", "$HOME", "x;y", "@list", "O'Brien", "100%"])
def test_windows_quoting_survives_powershell(tmp_path, word):
    import sys
    from predictor import gui_support
    script = tmp_path / "argv.py"
    script.write_text("import sys\nprint(sys.argv[1])\n")
    line = "& " + " ".join("'" + part.replace("'", "''") + "'" for part in (sys.executable, str(script)))
    result = subprocess.run(["powershell", "-NoProfile", "-Command",
                             f"{line} {gui_support.quote_command([word], windows=True)}"],
                            capture_output=True, text=True, timeout=60)
    assert result.stdout.strip() == word, result.stderr


def test_launcher_sees_a_listener_on_every_interface():
    from predictor import gui
    with socket.socket() as server:
        server.bind(("0.0.0.0", 0))
        server.listen(1)
        assert gui._port_free(server.getsockname()[1]) is False


def test_launcher_accepts_port_equals(monkeypatch):
    from predictor import gui
    seen = []
    monkeypatch.setattr(gui.subprocess, "call", lambda cmd: seen.append(cmd) or 0)
    assert gui.main(["--port=8600", "--no-browser"]) == 0 and "--server.port=8600" in seen[0]
    assert gui.main(["--port", "70000", "--no-browser"]) == 2


def test_gui_page_locks_recorded_options_and_shows_resumed_results(tmp_path, monkeypatch):
    pytest.importorskip("streamlit")
    from pathlib import Path

    import streamlit as st
    from streamlit.testing.v1 import AppTest

    from predictor import genome_run, resources
    from predictor.discovery import bitacora
    run = tmp_path / "run"
    run.mkdir()
    rows = _run(run)
    genome_run.summarise(rows, {"G__Fur__one": "done (resumed)"}, out_dir=run)
    (run / "prediction_input.json").write_text(json.dumps(
        {"policy": {"offline": True, "fold": False, "metalnet": False, "af3": True}}))
    monkeypatch.setattr(resources, "OUTPUT_DIR", tmp_path)
    monkeypatch.setattr(resources, "audit", lambda: {"database": {"release": "test"}})
    monkeypatch.setattr(bitacora, "engine_available", lambda: (True, ""))
    st.cache_data.clear()
    app = AppTest.from_file(str(Path(resources.__file__).with_name("gui_app.py")), default_timeout=60)
    app.session_state["last_shared_genome"] = str(run)
    app.run()
    assert not app.exception
    page = "\n".join(element.value for element in app.markdown)
    assert "1 predicted" in page and "Zn2+" in page and "failed" not in page
    metalnet = next(box for box in app.checkbox if box.label == "Use MetalNet2 when available")
    assert metalnet.value is False and metalnet.disabled
