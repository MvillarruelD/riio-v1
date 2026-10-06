"""Publication failure/recovery tests; no scientific records or databases are read."""
import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest


@pytest.fixture
def share(tmp_path, monkeypatch):
    analysis = tmp_path / "analysis"
    analysis.mkdir()
    rp = SimpleNamespace(TAG="synthetic", ANALYSIS=analysis, JOBS=tmp_path / "jobs",
                         MANIFEST=analysis / "run_manifest_synthetic.csv",
                         REPORT_DOCX=tmp_path / "report.docx", COMPLEMENT_TSV=tmp_path / "results.tsv",
                         REGULON_SUMMARY=tmp_path / "regulon.tsv", FIGURES=tmp_path / "figures",
                         RESULTS=tmp_path / "results")
    monkeypatch.setitem(sys.modules, "run_paths", rp)
    spec = importlib.util.spec_from_file_location("share_under_test", Path(__file__).parents[1] / "make_share_bundle.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(module, "PROJ", tmp_path)
    monkeypatch.setattr(module, "OUT", tmp_path / "SHARE_FOR_AUTHORS_synthetic")
    return module


def test_failed_assembly_leaves_prior_share_intact(share, monkeypatch):
    share.OUT.mkdir()
    (share.OUT / "old.txt").write_text("original", encoding="utf-8")
    def fail(destination, **kwargs):
        (destination / "partial.txt").write_text("partial", encoding="utf-8")
        raise OSError("copy failed")
    monkeypatch.setattr(share, "assemble", fail)
    assert share.main([]) == 1
    assert (share.OUT / "old.txt").read_text(encoding="utf-8") == "original"
    assert not list(share.PROJ.glob(".share-*"))


def test_success_preserves_old_release_without_mixing(share, monkeypatch):
    share.OUT.mkdir()
    (share.OUT / "old.txt").write_text("original", encoding="utf-8")
    def assemble(destination, **kwargs):
        (destination / "new.txt").write_text("new release", encoding="utf-8")
    monkeypatch.setattr(share, "assemble", assemble)
    share.publish()
    assert not (share.OUT / "old.txt").exists()
    assert (share.OUT / "new.txt").is_file()
    backups = list(share.PROJ.glob(".SHARE_FOR_AUTHORS_*.previous-*"))
    assert len(backups) == 1 and (backups[0] / "old.txt").is_file()


def test_summary_share_copies_inputs_and_keeps_benchmark_panels_separate(share):
    (share.PROJ / "input.fasta").write_text(">synthetic\nplaceholder\n", encoding="utf-8")
    share.RP.MANIFEST.write_text("run_name,fasta\nexample,input.fasta\n", encoding="utf-8")
    for path in (share.RP.REPORT_DOCX, share.RP.REGULON_SUMMARY):
        path.write_text("synthetic artifact", encoding="utf-8")
    share.RP.COMPLEMENT_TSV.write_text(
        "inducer_class\tssn_ligand\tgate\nmetal\tZn2+\tTrue\norganic\t\tFalse\n",
        encoding="utf-8")
    share.RP.FIGURES.mkdir()
    for name in ("F1.svg", "F1.pdf", "F21_example.png"):
        (share.RP.FIGURES / name).write_text("synthetic figure", encoding="utf-8")
    share.publish(no_bundles=True)
    assert (share.OUT / "03_data/inputs/example.fasta").is_file()
    assert "inputs/example.fasta" in (share.OUT / "03_data/candidate_manifest.csv").read_text(encoding="utf-8")
    assert (share.OUT / "02_figures/F1.pdf").is_file()
    assert (share.OUT / "06_benchmark_context/F21_example.png").is_file()
    manifest = json.loads((share.OUT / "share_manifest.json").read_text(encoding="utf-8"))
    assert manifest["individual_bundles_included"] is False
    assert "03_data/inputs/example.fasta" in manifest["artifacts"]


def test_missing_required_artifact_is_an_error(share):
    with pytest.raises(FileNotFoundError, match="missing report"):
        share.copy(share.PROJ / "absent", share.PROJ / "target", "report")


def test_full_share_contains_a_portable_verified_individual_bundle(share):
    from predictor.report.bundle_contract import inventory, validate_files
    (share.PROJ / "input.fasta").write_text(">synthetic\nplaceholder\n", encoding="utf-8")
    share.RP.MANIFEST.write_text("run_name,fasta\nexample,input.fasta\n", encoding="utf-8")
    for path in (share.RP.REPORT_DOCX, share.RP.REGULON_SUMMARY):
        path.write_text("synthetic artifact", encoding="utf-8")
    share.RP.COMPLEMENT_TSV.write_text(
        "inducer_class\tssn_ligand\tgate\nmetal\tZn2+\tTrue\norganic\t\tFalse\n",
        encoding="utf-8")
    share.RP.FIGURES.mkdir()
    (share.RP.FIGURES / "F1.svg").write_text("synthetic figure", encoding="utf-8")
    bundle = share.JOBS / "example"
    bundle.mkdir(parents=True)
    (bundle / "record.txt").write_text("synthetic record", encoding="utf-8")
    (bundle / "REPORT.html").write_text('<a href="record.txt">Record</a>', encoding="utf-8")
    (bundle / "bundle_manifest.json").write_text(json.dumps({
        "schema_version": 2, "status": "complete", "report_generated": True,
        "required_artifacts": ["record.txt", "REPORT.html"], "artifacts": inventory(bundle),
    }), encoding="utf-8")
    share.publish()
    assert validate_files(share.OUT / "05_per_regulator/example") == []
    package = json.loads((share.OUT / "share_manifest.json").read_text(encoding="utf-8"))
    assert "05_per_regulator/example/bundle_manifest.json" in package["artifacts"]


def test_failed_promotion_restores_previous_release(share, monkeypatch):
    share.OUT.mkdir()
    (share.OUT / "original.txt").write_text("original", encoding="utf-8")
    monkeypatch.setattr(share, "assemble", lambda *args, **kwargs: None)
    original_replace = Path.replace
    def fail_promotion(path, target):
        if path.name.startswith(".share-"):
            raise OSError("promotion failure")
        return original_replace(path, target)
    monkeypatch.setattr(Path, "replace", fail_promotion)
    with pytest.raises(OSError, match="promotion failure"):
        share.publish()
    assert (share.OUT / "original.txt").is_file()
    assert not list(share.PROJ.glob(".share-*"))
