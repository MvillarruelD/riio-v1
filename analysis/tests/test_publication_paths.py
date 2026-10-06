"""The consolidated release must not select an unrelated developer checkout."""
from pathlib import Path

from canonical_config import load_canonical
from project_config import PROJECT_ROOT, predictor_root


def test_local_predictor_takes_precedence(monkeypatch):
    monkeypatch.delenv("TFOP_REPO", raising=False)
    assert predictor_root() == PROJECT_ROOT


def test_canonical_roots_follow_user_output_directory(monkeypatch, tmp_path):
    monkeypatch.setenv("RIIO_RUNS_DIR", str(tmp_path))
    cfg = load_canonical()
    assert len(cfg["runs"]) == 3
    assert sum(run["n_candidates"] for run in cfg["runs"].values()) == 238
    for run in cfg["runs"].values():
        assert Path(run["root"]).parent == tmp_path
