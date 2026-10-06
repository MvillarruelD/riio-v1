"""Orchestration contracts for `run_pipeline.py`.

None of these run a prediction. The batch is hours long and needs external engines, so bundles are
synthesised on disk and the stage scripts are stubbed; what is under test is the DRIVER's decisions:
when it refuses to start, when it resumes, what it archives, and what it refuses to carry into a
report.

The driver keeps its resolved run in module globals (`MANIFEST`, `JOBS`, `RUN_TAG`) and reaches the
predictor through subprocesses. `pipeline_at()` rebinds those globals to a temporary run and stubs
the three subprocess boundaries, so each test exercises the real decision logic without a predictor,
a network, or a run directory in the repository.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

import run_pipeline as RPL

REQUIRED_ARTIFACTS = ("dossier.json", "operators_ranked.json", "binding_sites.tsv", "regulon.tsv")
CURRENT_FINGERPRINT = "current" * 8


def write_manifest(path: Path, rows: list[dict], fasta: Path) -> Path:
    """A manifest with `rows`, every row pointing at `fasta`.

    The path is written RELATIVE, because that is what the validator requires and what makes a
    manifest distributable -- an absolute path pins it to the machine that generated it.
    """
    header = "run_name,fasta,genome,organism_acc,family,family_flag\n"
    lines = [header]
    for r in rows:
        lines.append(",".join([
            r["run_name"], r.get("fasta", fasta.name), r.get("genome", "Genome sp."),
            r.get("organism_acc", "GCF_000005845.2"), r.get("family", "MerR"),
            r.get("family_flag", "(auto)"),
        ]) + "\n")
    path.write_text("".join(lines), encoding="utf-8")
    return path


def make_bundle(jobs: Path, name: str, *, status: str = "complete",
                fingerprint: str = CURRENT_FINGERPRINT,
                artifacts: tuple[str, ...] = REQUIRED_ARTIFACTS,
                manifest: bool = True) -> Path:
    """A bundle on disk, exactly as the predictor publishes one."""
    bundle = jobs / name
    bundle.mkdir(parents=True, exist_ok=True)
    for rel in artifacts:
        (bundle / rel).write_text("{}" if rel.endswith(".json") else "x\n", encoding="utf-8")
    if manifest:
        (bundle / "bundle_manifest.json").write_text(json.dumps({
            "status": status, "runtime": {"runtime_sha256": fingerprint},
        }), encoding="utf-8")
    return bundle


@pytest.fixture
def pipeline_at(tmp_path, monkeypatch, protein_fasta):
    """Point the driver at a temporary run and stub every subprocess boundary."""
    def _setup(names: list[str], *, tag: str = "testtag", rows: list[dict] | None = None):
        fasta = protein_fasta()
        run_root = tmp_path / f"run_{tag}"
        jobs = run_root / "jobs"
        jobs.mkdir(parents=True)
        manifest_path = tmp_path / f"run_manifest_{tag}.csv"
        write_manifest(manifest_path, rows or [{"run_name": n} for n in names], fasta)

        monkeypatch.setattr(RPL, "MANIFEST", manifest_path)
        monkeypatch.setattr(RPL, "JOBS", jobs)
        monkeypatch.setattr(RPL, "RUN_TAG", tag)
        monkeypatch.setattr(RPL, "PROJ", tmp_path)
        monkeypatch.setattr(RPL, "ANALYSIS", tmp_path)
        # Subprocess boundaries: the predictor fingerprint, the git commit, and run_paths resolution.
        monkeypatch.setattr(RPL, "current_runtime",
                            lambda: {"runtime_sha256": CURRENT_FINGERPRINT})
        monkeypatch.setattr(RPL, "predictor_commit", lambda: "deadbee")
        monkeypatch.setattr(RPL, "run_paths_env", lambda _tag: {
            "TFOP_JOBS": str(jobs), "TFOP_RUN_LOGS": str(run_root / "logs"),
            "AN_FIG": str(run_root / "figures"), "AN_OUT": str(run_root / "results"),
            "TFOP_RUN_REGULON": str(run_root / "regulon"),
        })
        return jobs, manifest_path, run_root
    return _setup


# --------------------------------------------------------------------------- progress and state
def test_new_run_is_ready_and_reports_everything_pending(pipeline_at, capsys):
    pipeline_at(["a", "b", "c"])
    assert RPL.check() == 0
    out = capsys.readouterr().out
    assert "0 complete, 3 pending" in out
    assert "READY" in out


def test_partially_complete_run_is_a_resume_and_needs_no_flag(pipeline_at, capsys):
    jobs, _, _ = pipeline_at(["a", "b", "c"])
    make_bundle(jobs, "a")
    assert RPL.check() == 0
    out = capsys.readouterr().out
    assert "1 complete, 2 pending" in out
    assert "RESUME" in out


def test_completed_run_is_refused_because_a_rerun_would_change_nothing(pipeline_at, capsys):
    jobs, _, _ = pipeline_at(["a", "b"])
    for n in ("a", "b"):
        make_bundle(jobs, n)
    assert RPL.check() == 1
    assert "BLOCKED" in capsys.readouterr().out


def test_a_bundle_without_its_manifest_is_not_done(pipeline_at):
    jobs, _, _ = pipeline_at(["a"])
    make_bundle(jobs, "a", manifest=False)
    assert RPL.is_done("a") is False
    assert RPL.run_progress() == ([], ["a"])


def test_a_bundle_missing_a_required_artifact_is_not_done(pipeline_at):
    jobs, _, _ = pipeline_at(["a"])
    make_bundle(jobs, "a", artifacts=("dossier.json", "operators_ranked.json"))
    assert RPL.is_done("a") is False


def test_a_bundle_marked_incomplete_is_not_done(pipeline_at):
    jobs, _, _ = pipeline_at(["a"])
    make_bundle(jobs, "a", status="partial")
    assert RPL.is_done("a") is False


# --------------------------------------------------------------------------- staleness
def test_a_bundle_from_another_runtime_is_stale(pipeline_at):
    jobs, _, _ = pipeline_at(["a", "b"])
    make_bundle(jobs, "a")
    make_bundle(jobs, "b", fingerprint="older" * 10)
    assert [n for n, _ in RPL.stale_bundles()] == ["b"]


def test_downstream_work_is_refused_over_stale_bundles(pipeline_at, capsys):
    jobs, _, _ = pipeline_at(["a"])
    make_bundle(jobs, "a", fingerprint="older" * 10)
    assert RPL.validate_bundle_inputs(allow_stale=False) is False
    assert "another code/data runtime" in capsys.readouterr().out


def test_allow_stale_permits_downstream_inspection_but_says_so(pipeline_at, capsys):
    """`--allow-stale` is for historical inspection only, and must announce what it accepted."""
    jobs, _, _ = pipeline_at(["a"])
    make_bundle(jobs, "a", fingerprint="older" * 10)
    assert RPL.validate_bundle_inputs(allow_stale=True) is True
    assert "WARNING" in capsys.readouterr().out


def test_downstream_work_is_refused_while_any_bundle_is_pending(pipeline_at, capsys):
    jobs, _, _ = pipeline_at(["a", "b"])
    make_bundle(jobs, "a")
    assert RPL.validate_bundle_inputs(allow_stale=True) is False
    assert "incomplete" in capsys.readouterr().out


def test_allow_stale_cannot_wave_through_a_no_op_batch(pipeline_at, capsys):
    """The flag covers downstream stages; letting it also pass a complete batch is the old bug."""
    jobs, manifest, _ = pipeline_at(["a"])
    make_bundle(jobs, "a")
    rc = RPL.main(["--tag", "testtag", "--manifest", manifest.name,
                   "--only", "batch", "--allow-stale", "--dry-run"])
    assert rc == 1
    assert "refusing to start" in capsys.readouterr().out


# --------------------------------------------------------------------------- archiving
def test_archive_moves_only_the_manifests_own_bundles(pipeline_at, capsys):
    jobs, _, _ = pipeline_at(["a", "b"])
    make_bundle(jobs, "a")
    make_bundle(jobs, "b")
    make_bundle(jobs, "benchmark_thing")          # not in the manifest
    assert RPL.archive("t1") == 0
    dest = jobs.parent / f"{jobs.name}_t1_archive"
    assert sorted(p.name for p in dest.iterdir()) == ["a", "b"]
    assert (jobs / "benchmark_thing").is_dir(), "a foreign bundle must be left in place"


def test_foreign_bundles_are_reported_but_never_counted_as_progress(pipeline_at):
    jobs, _, _ = pipeline_at(["a"])
    make_bundle(jobs, "benchmark_thing")
    assert [p.name for p in RPL.foreign_bundles()] == ["benchmark_thing"]
    assert RPL.candidate_bundles() == []
    assert RPL.run_progress() == ([], ["a"])


def test_archive_refuses_to_merge_two_runs_into_one_directory(pipeline_at, capsys):
    jobs, _, _ = pipeline_at(["a"])
    make_bundle(jobs, "a")
    assert RPL.archive("t1") == 0
    make_bundle(jobs, "a")
    assert RPL.archive("t1") == 1
    assert "already exists" in capsys.readouterr().out


# --------------------------------------------------------------------------- manifest gating
def test_a_duplicate_row_blocks_the_run(pipeline_at, capsys):
    pipeline_at(["a"], rows=[{"run_name": "a"}, {"run_name": "a"}])
    assert RPL.check() == 2
    assert "duplicate output target" in capsys.readouterr().out


def test_run_names_differing_only_in_case_are_one_output_target(pipeline_at, capsys):
    pipeline_at(["a"], rows=[{"run_name": "Cand_1"}, {"run_name": "cand_1"}])
    assert RPL.check() == 2
    assert "duplicate output target" in capsys.readouterr().out


def test_a_missing_fasta_blocks_the_run(pipeline_at, capsys):
    pipeline_at(["a"], rows=[{"run_name": "a", "fasta": "analysis/does_not_exist.fasta"}])
    assert RPL.check() == 2
    assert "FASTA missing" in capsys.readouterr().out


def test_a_nucleotide_fasta_is_rejected_by_the_production_validator(tmp_path, pipeline_at, capsys):
    bad = tmp_path / "nuc.fasta"
    bad.write_text(">n\n" + "ATGC" * 40 + "\n", encoding="utf-8")
    pipeline_at(["a"], rows=[{"run_name": "a", "fasta": bad.name}])
    assert RPL.check() == 2
    assert "nucleotide" in capsys.readouterr().out


def test_a_multi_record_fasta_is_rejected(tmp_path, pipeline_at, capsys):
    bad = tmp_path / "two.fasta"
    bad.write_text(">a\nMSTNPKPQRKTKRNTNRRPQ\n>b\nMPEPTIDESEQUENCEAAAA\n", encoding="utf-8")
    pipeline_at(["a"], rows=[{"run_name": "a", "fasta": bad.name}])
    assert RPL.check() == 2
    assert "exactly one FASTA record" in capsys.readouterr().out


# --------------------------------------------------------------------------- stage sequencing
def test_a_failing_stage_stops_every_downstream_stage(pipeline_at, monkeypatch, capsys):
    jobs, manifest, _ = pipeline_at(["a"])
    make_bundle(jobs, "a")
    ran: list[str] = []

    def fake_run_script(name, extra, *, dry):
        ran.append(name)
        return 1 if name == "aggregate_phase1.py" else 0

    monkeypatch.setattr(RPL, "run_script", fake_run_script)
    rc = RPL.main(["--tag", "testtag", "--manifest", manifest.name, "--from", "aggregate"])
    assert rc == 1
    assert ran == ["aggregate_phase1.py"], "downstream stages must not run after a failure"
    assert "stopping so the failure is not carried" in capsys.readouterr().out


@pytest.mark.parametrize("flags, expected", [([], "--no-fold"), (["--no-fold"], "--no-fold"),
                                             (["--fold"], "--fold")])
def test_the_batch_folds_only_when_asked(pipeline_at, monkeypatch, flags, expected):
    """The canonical run is structure OFF, so a batch started without a flag must not fold."""
    _, manifest, _ = pipeline_at(["a"])
    seen: list[list[str]] = []

    def fake_run_script(name, extra, *, dry):
        seen.append(extra)
        return 0

    monkeypatch.setattr(RPL, "run_script", fake_run_script)
    assert RPL.main(["--tag", "testtag", "--manifest", manifest.name, "--only", "batch", *flags]) == 0
    (extra,) = seen
    assert expected in extra and ({"--fold", "--no-fold"} - {expected}).isdisjoint(extra)


def test_fold_and_no_fold_cannot_both_be_given(pipeline_at):
    _, manifest, _ = pipeline_at(["a"])
    with pytest.raises(SystemExit):
        RPL.main(["--tag", "testtag", "--manifest", manifest.name, "--fold", "--no-fold"])


def test_a_missing_figure_script_fails_the_stage_and_is_recorded(pipeline_at, monkeypatch, capsys):
    """A configured figure whose script is absent must fail -- never be skipped with exit 0."""
    jobs, manifest, run_root = pipeline_at(["a"])
    make_bundle(jobs, "a")
    present = "fig_present.py"
    monkeypatch.setattr(RPL, "FIGURE_SCRIPTS", (present, "fig_absent.py"))
    (RPL.ANALYSIS / present).write_text("# stub\n", encoding="utf-8")
    monkeypatch.setattr(RPL, "run_script", lambda name, extra, *, dry: 0)

    rc = RPL.main(["--tag", "testtag", "--manifest", manifest.name, "--only", "figures"])
    assert rc == 1
    out = capsys.readouterr().out
    assert "fig_absent.py" in out and "MISSING" in out

    recorded = json.loads((run_root / "results" / "figures_manifest.json").read_text(encoding="utf-8"))
    assert recorded["missing_scripts"] == ["fig_absent.py"]
    assert "fig_absent.py" in recorded["failed"]
    assert recorded["configured"] == [present, "fig_absent.py"]


def test_figure_qa_runs_only_after_every_panel_succeeded(pipeline_at, monkeypatch):
    """A partial figure set must never reach QA, and so must never reach the report."""
    jobs, manifest, _ = pipeline_at(["a"])
    make_bundle(jobs, "a")
    monkeypatch.setattr(RPL, "FIGURE_SCRIPTS", ("fig_one.py", "fig_two.py"))
    for f in ("fig_one.py", "fig_two.py"):
        (RPL.ANALYSIS / f).write_text("# stub\n", encoding="utf-8")
    ran: list[str] = []

    def fake_run_script(name, extra, *, dry):
        ran.append(name)
        return 1 if name == "fig_two.py" else 0

    monkeypatch.setattr(RPL, "run_script", fake_run_script)
    assert RPL.main(["--tag", "testtag", "--manifest", manifest.name, "--only", "figures"]) == 1
    assert "figure_qa.py" not in ran


# --------------------------------------------------------------------------- dry run
def test_dry_run_plans_a_run_that_has_not_started(pipeline_at, monkeypatch, capsys):
    """The planning contract: no stage executes, and downstream stages do not demand bundles."""
    _, manifest, _ = pipeline_at(["a", "b"])
    # The temporary ANALYSIS holds no figure scripts, and a missing one is (correctly) a failure --
    # so give the plan a figure set that exists. What is under test here is the PLANNING contract.
    monkeypatch.setattr(RPL, "FIGURE_SCRIPTS", ("fig_stub.py",))
    (RPL.ANALYSIS / "fig_stub.py").write_text("# stub\n", encoding="utf-8")
    executed: list[str] = []
    real_run_script = RPL.run_script

    def spy(name, extra, *, dry):
        assert dry is True, "a dry run must not execute a stage"
        executed.append(name)
        return real_run_script(name, extra, dry=True)

    monkeypatch.setattr(RPL, "run_script", spy)
    rc = RPL.main(["--tag", "testtag", "--manifest", manifest.name, "--dry-run"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "plan only: 0 complete, 2 pending" in out
    assert "run_phase1.py" in executed and "aggregate_phase1.py" in executed


def test_dry_run_from_discovery_defers_validation_of_the_future_manifest(
        pipeline_at, monkeypatch, capsys):
    """A new tag is plannable before discovery has generated its manifest."""
    _, manifest, _ = pipeline_at(["a"])
    manifest.unlink()
    monkeypatch.setattr(RPL, "FIGURE_SCRIPTS", ("fig_stub.py",))
    (RPL.ANALYSIS / "fig_stub.py").write_text("# stub\n", encoding="utf-8")

    assert RPL.main(["--tag", "testtag", "--dry-run"]) == 0
    out = capsys.readouterr().out
    assert "will be created by discovery" in out
    assert "manifest validation is deferred" in out
    assert "candidate count is unknown until discovery" in out


def test_batch_only_dry_run_still_rejects_a_missing_manifest(
        pipeline_at, monkeypatch, capsys):
    """Planning cannot invent a batch input unless discovery is part of the same plan."""
    _, manifest, _ = pipeline_at(["a"])
    manifest.unlink()

    assert RPL.main(["--tag", "testtag", "--only", "batch", "--dry-run"]) == 2
    assert "manifest validation   : FAILED" in capsys.readouterr().out


# --------------------------------------------------------------------------- run-tag isolation
def test_two_tags_never_share_a_job_directory(pipeline_at):
    jobs_a, _, root_a = pipeline_at(["x"], tag="alpha")
    make_bundle(jobs_a, "x")
    jobs_b, _, root_b = pipeline_at(["x"], tag="beta")
    assert jobs_a != jobs_b
    assert RPL.run_progress() == ([], ["x"]), "a new tag must not see another tag's bundles"


# --------------------------------------------------------------------------- Windows path limit
def test_path_headroom_is_measured_from_the_deepest_bundle_path(pipeline_at):
    """The number reported must be the path that actually gets created, not the jobs root."""
    jobs, _, _ = pipeline_at(["short"])
    total, headroom, driver = RPL.path_headroom()
    assert driver == "short"
    assert total > len(str(jobs)), "headroom must account for the bundle tree, not just the root"
    assert headroom == RPL._WINDOWS_MAX_PATH - total


def test_the_longest_run_name_drives_the_measurement(pipeline_at):
    long_name = "G" * 60
    pipeline_at([], rows=[{"run_name": "a"}, {"run_name": long_name}])
    _total, _headroom, driver = RPL.path_headroom()
    assert driver == long_name


def test_a_run_root_that_exceeds_the_windows_limit_is_refused(pipeline_at, monkeypatch, capsys):
    """This is the failure a real one-candidate run hit: the prediction succeeded after 65 minutes
    and then died with WinError 206 while PUBLISHING the bundle. One second of checking replaces it.
    """
    monkeypatch.setattr(RPL.sys, "platform", "win32")
    pipeline_at(["cand"])
    monkeypatch.setattr(RPL, "JOBS", Path("C:/" + "d" * 240 + "/jobs"))
    assert RPL.check() == 2
    out = capsys.readouterr().out
    assert "over the 260-character Windows limit" in out
    assert "PUBLISHING the bundle" in out


def test_a_comfortable_run_root_passes_and_reports_its_headroom(pipeline_at, capsys):
    pipeline_at(["cand"])
    assert RPL.check() == 0
    assert "under the 260-character Windows limit" in capsys.readouterr().out


def test_the_limit_is_not_enforced_off_windows(pipeline_at, monkeypatch, capsys):
    """POSIX has no 260-character ceiling; the number is still reported, but it never blocks."""
    monkeypatch.setattr(RPL.sys, "platform", "linux")
    pipeline_at(["cand"])
    # POSIX permits long total paths, but each individual component must stay below NAME_MAX.
    monkeypatch.setattr(RPL, "JOBS", Path("/" + "d" * 150 + "/" + "d" * 150 + "/jobs"))
    assert RPL.check() == 0
