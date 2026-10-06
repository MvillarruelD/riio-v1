"""A conservation set under the floor has to SAY so, in the dossier and in the report.

`msa_homologs.collect_msa_regions` records `shortfall` when fewer than MIN_MAPPABLE_ORTHOLOGS
orthologs map. It reached `run_record.json` and stopped there -- but the report, the figures and the
analysis scorers read `dossier.json`, so a thin conservation axis was invisible where anyone looks for
it. Measured on the 2026-09-17 canonical runs: 9 of 238 candidates were short and no dossier said so.
"""
from __future__ import annotations

from predictor.annotate import msa_homologs as mh
from predictor.report.render_report import build_html

#: The real shape, from SentSL1344__MerR__WP_000122327.1 of the canonical run.
SHORTFALL = {"n_mappable": 36, "n_mapped": 12, "floor": 20,
             "reason": "UniProt returned no RefSeq/EMBL cross-reference for these accessions"}


def _data(dossier_extra: dict) -> dict:
    return {
        "dossier": {"tf_id": "T", "family": "MerR", "homolog_regions": {"pre_msa": [], "expanded": [],
                                                                       **dossier_extra}},
        "run_record": {"homolog_collection": {
            "attempts": [{"route": "blast_similarity", "selected": True}],
            "msa": {"msa_added": 3, "source": "colabfold",
                    "source_record": {"n_orthologs": 36, "n_mapped": 12}}}},
        "bundle_files": set(),
    }


def test_the_floor_is_the_value_the_handoff_specified():
    assert mh.MIN_MAPPABLE_ORTHOLOGS == 20


def test_report_states_the_shortfall_and_what_it_costs():
    html = build_html(_data({"shortfall": SHORTFALL}), {})
    assert "Thin conservation set" in html, "a short conservation set must be stated, not implied"
    assert "12" in html and "36" in html and "20" in html, "counts and the floor must all appear"
    assert "no RefSeq/EMBL cross-reference" in html, "the recorded reason must reach the reader"


def test_a_met_floor_says_nothing_about_a_shortfall():
    """Absence is the signal that the set was adequate; the phrase must not appear unprompted."""
    assert "Thin conservation set" not in build_html(_data({}), {})


def test_report_renders_when_the_msa_never_ran():
    """The common case, and the one dd24d85 broke: homolog attempts exist but no MSA expansion ran.

    190 of 238 canonical candidates are like this. The shortfall note read a variable that was only
    bound inside the MSA branch, so every one of them died with UnboundLocalError at the output stage
    -- after the whole pipeline had run. The fixtures above all carried an MSA record, so they never
    reached that path; this one does.
    """
    data = _data({})
    data["run_record"]["homolog_collection"] = {
        "attempts": [{"route": "ssn_cluster", "selected": True}], "msa": {"source": None}}
    html = build_html(data, {})
    assert "Homolog-source ledger" in html and "Thin conservation set" not in html
