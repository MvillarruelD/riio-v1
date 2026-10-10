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
