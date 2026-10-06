"""Defects found by exercising a complete production bundle rather than isolated helpers."""
from __future__ import annotations

import csv

from predictor.effector import substrate_map
from predictor.report import operators, outputs


def test_actp_is_not_mislabeled_as_a_copper_atpase():
    hit = substrate_map.classify_gene("acetate/glycolate:cation symporter", "actP")
    assert hit is None
    assert "actp" not in substrate_map.keywords_for_effector("Cu+")


def test_ranked_operators_retain_locus_provenance():
    dossier = {
        "family": "MerR",
        "rescan": {"hits": [{
            "accession": "NC_1.1", "start": 100, "end": 116, "strand": "+", "dyad": 108,
            "score": 12.5, "seq": "GACAGCGTTGCCAGCT", "generator": "rescan_tier1",
        }]},
        "neighborhood": {"window": [50, 150]},
    }

    result = operators.build_operators(dossier)

    assert result["ranked"][0]["provenance"]["start"] == 100
    assert result["ranked"][0]["provenance"]["locality"] == "proximal"
    assert result["primary"] == result["ranked"][0]


def test_binding_site_export_keeps_calibrated_significance(tmp_path):
    dossier = {
        "rescan": {"hits": [{
            "accession": "NC_1.1", "start": 100, "end": 116, "strand": "+", "dyad": 108,
            "score": 12.5, "pvalue": 1.25e-8, "qvalue": 3.75e-4,
            "qvalue_kind": "conservative_bh_upper_bound", "generator": "rescan_tier1",
            "seq": "GACAGCGTTGCCAGCT",
        }]},
        "neighborhood": {"window": [50, 150]},
    }
    path = tmp_path / "binding_sites.tsv"

    assert outputs.write_binding_sites(dossier, path) == 1
    with path.open(encoding="utf-8") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))

    assert row["pvalue"] == "1.250e-08"
    assert row["qvalue"] == "3.750e-04"
    assert row["qvalue_kind"] == "conservative_bh_upper_bound"
    assert row["generator"] == "rescan_tier1"
