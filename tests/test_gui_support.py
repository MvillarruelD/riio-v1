"""Unit tests for GUI state that must remain safe across Streamlit reruns."""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import io
import zipfile

import pytest

from predictor.gui_support import (available_job_name, context_archive, report_preview,
                                   reserve_output_dir, safe_name)


def test_safe_name_is_bounded_and_portable():
    assert safe_name("  Zn sensor / strain:1  ") == "Zn_sensor_strain_1"
    assert safe_name("...", "TF_prediction") == "TF_prediction"
    assert len(safe_name("x" * 200)) == 80


def test_available_job_name_suffixes_existing_results(tmp_path):
    (tmp_path / "my_TF").mkdir()
    (tmp_path / "my_TF_2").mkdir()
    assert available_job_name(tmp_path, "my TF") == "my_TF_3"


def test_output_directory_reservation_is_atomic(tmp_path):
    with ThreadPoolExecutor(max_workers=6) as pool:
        paths = list(pool.map(lambda _: reserve_output_dir(tmp_path, "scan_genome"), range(12)))

    assert len(set(paths)) == 12
    assert all(path.is_dir() for path in paths)
    assert {path.name for path in paths} == {"scan_genome", *(f"scan_genome_{i}" for i in range(2, 13))}


def test_full_context_archive_preserves_relative_report_links_and_excludes_drafts(tmp_path):
    for name, data in {"RUN_REPORT.html": b'<a href="jobs/TF/REPORT.html">TF</a>',
                       "jobs/TF/REPORT.html": b'<img src="figures/logo.png">',
                       "jobs/TF/figures/logo.png": b"saved image", "af3/job.json": b"{}",
                       "discovery/DISCOVERY_REPORT.html": b"discovery",
                       "jobs/.TF.stage-interrupted/dossier.json": b"draft"}.items():
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
    before = {p: p.read_bytes() for p in tmp_path.rglob("*") if p.is_file()}
    data, skipped = context_archive(tmp_path)
    with zipfile.ZipFile(io.BytesIO(data)) as archive:
        assert set(archive.namelist()) == {
            "RUN_REPORT.html", "jobs/TF/REPORT.html", "jobs/TF/figures/logo.png",
            "af3/job.json", "discovery/DISCOVERY_REPORT.html"}
        assert archive.read("jobs/TF/figures/logo.png") == b"saved image"
    assert skipped == ["jobs/.TF.stage-interrupted"]
    assert all(path.read_bytes() == value for path, value in before.items())


def test_archive_and_preview_do_not_follow_external_links(tmp_path):
    root = tmp_path / "result"
    root.mkdir()
    secret = tmp_path / "outside.png"
    secret.write_bytes(b"outside must not be read")
    try:
        (root / "linked.png").symlink_to(secret)
        (root / "linked_dir").symlink_to(tmp_path, target_is_directory=True)
    except OSError:
        pytest.skip("symlink creation unavailable")
    data, skipped = context_archive(root)
    with zipfile.ZipFile(io.BytesIO(data)) as archive:
        assert not archive.namelist()
    assert skipped == ["linked.png", "linked_dir"]
    preview = report_preview('<img src="linked.png"><img src="../outside.png">', root)
    assert "base64" not in preview


def test_preview_embeds_only_contained_images_and_disables_file_anchors(tmp_path):
    (tmp_path / "figures").mkdir()
    (tmp_path / "figures" / "my logo.png").write_bytes(b"saved image")
    text = ('<img src="figures/my%20logo.png"><a href="jobs/TF/REPORT.html">TF</a>'
            "<a href='operators.tsv'>Data</a><a href='#section'>Section</a>"
            '<a href="https://example.org">Reference</a><a href="file:///private">File</a>')
    preview = report_preview(text, tmp_path)
    assert 'src="data:image/png;base64,' in preview
    assert preview.count('aria-disabled="true"') == 3
    assert 'href="#section"' not in preview  # original single-quoted fragment remains intact
    assert "href='#section'" in preview
    assert 'href="https://example.org"' in preview
    assert 'href="jobs/' not in preview and "href='operators" not in preview


# ----------------------------------------------------------------------------- headline calls and saved runs
import json  # noqa: E402
import os  # noqa: E402
import shlex  # noqa: E402
from pathlib import Path  # noqa: E402

from predictor import gui_support as support  # noqa: E402


def _bundle(path, inducers, **extra):
    path.mkdir(parents=True)
    dossier = {"family": "Fur", "family_support": "family_validated", "genome_accession": "NC_000913.3",
               "inducers": inducers, "operators": {"primary": {"sequence": "TTACTTTGTTATAAGACGTGA",
                                                               "source": "conservation_natural",
                                                               "provenance": {"locality": "distal"}}},
               "rescan": {"hits": [{}, {}, {}]}, "regulon": [{"first_gene": "a"}, {"first_gene": "b"}, "x"], **extra}
    (path / "dossier.json").write_text(json.dumps(dossier), encoding="utf-8")
    return path


def test_headline_reads_the_stored_class_and_exact_agreement(tmp_path):
    head = support.bundle_headline(_bundle(tmp_path / "b", {
        "top": "Fe2+", "agreement": 1.0, "sensor_class": "metal", "finalized": True, "top_display": "Fe2+"}))
    assert (head["inducer"], head["inducer_class"], head["agreement"], head["agreement_exact"]) == \
        ("Fe2+", "metal", 1.0, True)
    assert head["inducer_label"] == "Fe2+"
    assert (head["n_sites"], head["n_regulon_operons"], head["operator_locality"]) == (3, 2, "distal")


def test_a_saved_final_null_is_unresolved_not_organic(tmp_path):
    """The rerun stored 31 bundles with a null top inducer; they must never be drawn as organic."""
    head = support.bundle_headline(_bundle(tmp_path / "b", {"top": None, "finalized": True, "sensor_class": None,
                                                           "top_display": ""}, effector="Zn2+"))
    assert head["inducer"] is None and head["inducer_class"] == "unresolved"
    assert head["inducer_label"] is None, "a stale effector field must not replace a saved final null"


def test_summary_rows_normalise_class_without_inventing_one(tmp_path):
    (tmp_path / "summary.tsv").write_text(
        "run_name\tinducer_top\tinducer_class\tinducer_class_provenance\tagreement\tstatus\n"
        "a\tFe2+\tmetal\texplicit\t1.0\tok\n"
        "b\t\t\texplicit\t0.0\tok\n"
        "c\tfoo\torganic\texplicit\t\tok\n"
        "d\tfoo\tbogus\texplicit\t0.5\tok\n"
        "e\t\t\t\t\tnot run\n", encoding="utf-8")
    rows = support.read_summary(tmp_path)
    assert [r["class"] for r in rows] == ["metal", "unresolved", "organic", "unresolved", "unresolved"]
    assert rows[0]["agreement_value"] == 1.0 and rows[2]["agreement_value"] is None
    assert support.read_summary(tmp_path / "missing") == []


def test_saved_runs_are_listed_newest_first_and_ignore_other_folders(tmp_path):
    for name, stamp in (("old", 1_700_000_000), ("new", 1_800_000_000)):
        run = tmp_path / name
        run.mkdir()
        (run / "manifest.csv").write_text("run_name,family\nA,Fur\nB,Fur\n", encoding="utf-8")
        (run / "summary.tsv").write_text("run_name\tstatus\nA\tok\nB\tfailed\n", encoding="utf-8")
        os.utime(run, (stamp, stamp))
    (tmp_path / "not_a_run").mkdir()
    found = support.list_saved_runs(tmp_path)
    assert [r["name"] for r in found] == ["new", "old"]
    assert (found[0]["total"], found[0]["done"]) == (2, 1)
    assert support.list_saved_runs(tmp_path / "missing") == []


def test_saved_bundle_rebuilds_the_result_keys_the_page_reads(tmp_path):
    bundle = _bundle(tmp_path / "jobs" / "zntR", {"top": "Zn2+", "agreement": 0.5, "finalized": True,
                                                  "sensor_class": "metal"}, tf_id="zntR")
    raw = support.result_dict(bundle)
    assert raw["inducer"] == "Zn2+" and raw["operators_primary"]["sequence"].startswith("TTACTT")
    assert raw["bundle"] == str(bundle) and raw["genome"] == "NC_000913.3"
    assert [b["name"] for b in support.list_saved_bundles(tmp_path / "jobs")] == ["zntR"]


# ----------------------------------------------------------------------------- terminal commands
def test_commands_are_quoted_for_each_shell():
    spaced = "C:\\Users\\Ana Lopez\\runs\\ecoli"
    windows = support.cli_command("genome", windows=True, source="GCF_000005845.2", out=spaced, jobs=4,
                                  select=["a", "b"], offline=True, metalnet=False)
    assert windows == ('tfop genome GCF_000005845.2 --out "C:\\Users\\Ana Lopez\\runs\\ecoli" --jobs 4 '
                       '--select a,b --offline --no-metalnet')
    posix = support.cli_command("genome", windows=False, source="GCF_000005845.2", out="/home/ana/my runs/ecoli")
    assert posix == "tfop genome GCF_000005845.2 --out '/home/ana/my runs/ecoli'"


@pytest.mark.parametrize("kind,options", [
    ("genome", dict(source="GCF_000005845.2", out="o", jobs=2, select=["x"], fold=True, offline=True,
                    metalnet=False, af3=False, discover_only=True, census_fallback=True)),
    ("genome", dict(fasta="g.fna", gff="g.gff", organism="NC_1", out="o")),
    ("genome", dict(genbank="g.gbk", out="o")),
    ("genome", dict(proteome="p.faa", out="o")),
    ("batch", dict(source="p.faa", out="o", organism="GCF_1", jobs=3, fold=True)),
    ("predict", dict(source="WP_1.1", name="X", family="Fur", organism="GCF_1", effector="Zn2+",
                     fold=True, offline=True, metalnet=False, af3=False)),
])
def test_printed_commands_are_accepted_by_the_real_tfop_parser(monkeypatch, kind, options):
    """A flag the GUI prints must exist in ``tfop``: parse what it generates (nothing is executed)."""
    from predictor import run_cli
    ran = []
    for name in ("cmd_predict", "cmd_genome", "cmd_batch"):
        monkeypatch.setattr(run_cli, name, lambda a, _n=name: ran.append(_n) or 0)
    argv = shlex.split(support.cli_command(kind, windows=False, **options))[1:]
    assert run_cli.main(argv) == 0 and len(ran) == 1, argv


def test_unknown_command_kind_is_rejected():
    with pytest.raises(ValueError):
        support.cli_command("nope", out="x")


def test_pasted_input_becomes_fasta_or_an_accession_table():
    assert support.pasted_protein_input(">A\nMKV\n") == (".faa", ">A\nMKV\n")
    suffix, text = support.pasted_protein_input("WP_000240047.1, NP_417751.1  # zinc\n\nQ9X0Y1\n")
    assert suffix == ".tsv" and text == "accession\nWP_000240047.1\nNP_417751.1\nQ9X0Y1\n"
    for bad in ("", "   ", "# only a comment", "MKV*;rm"):
        with pytest.raises(ValueError):
            support.pasted_protein_input(bad)


@pytest.mark.skipif(os.name == "nt", reason="the WSL path spelling only exists on a Linux host")
def test_wsl_paths_get_a_windows_spelling(monkeypatch):
    monkeypatch.setattr(support.sys, "platform", "linux")
    monkeypatch.setenv("WSL_DISTRO_NAME", "Ubuntu")
    assert support.display_paths(Path("/home/ana/results/run1")) == [
        ("Folder", "/home/ana/results/run1"), ("In Windows", "\\\\wsl.localhost\\Ubuntu\\home\\ana\\results\\run1")]
    assert support.display_paths(Path("/mnt/c/Users/ana/run1"))[1] == ("In Windows", "C:\\Users\\ana\\run1")
    monkeypatch.delenv("WSL_DISTRO_NAME")
    assert len(support.display_paths(Path("/home/ana/run1"))) == 1
