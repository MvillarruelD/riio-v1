"""Display-only inducer labels must not alter the analysis classification contract."""
from __future__ import annotations

import json
import os
import sys

from project_config import PREDICTOR_ROOT

os.environ.setdefault("TFOP_RUN_TAG", "inducerdisplaytest")
sys.path.insert(0, str(PREDICTOR_ROOT))

import fig5_complement  # noqa: E402
from inducer_evidence import display_inducer  # noqa: E402


def _legacy_inducers():
    return {
        "top": "divalent metal (ion unresolved)",
        "candidates": ["Co2+", "formaldehyde"],
        "agreement": 0.5,
    }


def test_display_prefers_persisted_value_and_backfills_legacy_bundle():
    legacy = _legacy_inducers()
    assert display_inducer(legacy).startswith("inconclusive: Co2+ or formaldehyde")
    assert display_inducer({**legacy, "top_display": "persisted label"}) == "persisted label"


def test_fig5_uses_display_text_but_classifies_and_counts_from_top(tmp_path, monkeypatch):
    manifest = tmp_path / "manifest.csv"
    manifest.write_text(
        "run_name,family,organism_acc,gene,product\nrcnr,CsoR,GCF_000195955.2,rcnR,sensor\n",
        encoding="utf-8",
    )
    jobs = tmp_path / "jobs"
    bundle = jobs / "rcnr"
    bundle.mkdir(parents=True)
    (bundle / "dossier.json").write_text(
        json.dumps({"inducers": _legacy_inducers()}), encoding="utf-8"
    )
    monkeypatch.setattr(fig5_complement.RP, "MANIFEST", manifest)
    monkeypatch.setattr(fig5_complement, "JOBS", jobs)

    row = fig5_complement.collect()[0]
    assert row["inducer"].startswith("inconclusive: Co2+ or formaldehyde")
    assert row["inducer_top"] == "divalent metal (ion unresolved)"
    assert row["inducer_class"] == "metal"
