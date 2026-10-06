"""Unit tests for GUI state that must remain safe across Streamlit reruns."""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor

from predictor.gui_support import available_job_name, reserve_output_dir, safe_name


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
