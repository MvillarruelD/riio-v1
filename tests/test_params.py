"""The parameter registry resolves, and the values the canonical runs used are pinned.

A failing pin is not a bug to silence: it means a prediction-affecting parameter changed. Update the
pin in the same commit, say why, and expect the benchmark and survey numbers to move.
"""
from __future__ import annotations

import inspect

from predictor import params

#: values of the 2026-09-17 canonical runs (predictor 9fe0c32)
PINNED = {
    "BITACORA E-value": "1e-5",
    "HSP coverage, query": "0.50",
    "HSP coverage, subject": "0.50",
    "family E-value": 1e-3,
    "clade minimum identity": 0.60,
    "clade minimum support": 0.5,
    "curated-member identity": 0.5,
    "Ligify suppression weight": 0.5,
    "seed motifs kept (top_k)": 8,
    "seed width tolerance (bp)": 3,
    "seed selection method": "cluster0",
    "seed minimum hits": 3,
    "seed maximum hits": 40,
    "seed keep fraction": 0.75,
    "seed trim IC (bits)": 0.35,
    "rescan scope": "intergenic",
    "rescan p-value": 1e-4,
    "locality priors by tier": {1: 1.0, 2: 0.6, 3: 0.35, 4: 0.1},
    "regulon site p-value": 1e-4,
    "operon maximum gap (bp)": 150,
    "site-to-gene distance (bp)": 400,
    "operon display cap": 25,
    "minimum SSN promoter set": 4,
    "autoregulation distance (bp)": 80,
}


def test_every_parameter_resolves():
    rows = params.table()
    assert len(rows) == len(params.PARAMS)
    assert all(r["value"] is not inspect.Parameter.empty for r in rows)
    assert len({r["name"] for r in rows}) == len(rows), "parameter names must be unique"


def test_canonical_values_are_pinned():
    got = {r["name"]: r["value"] for r in params.table()}
    drift = {k: (v, got.get(k)) for k, v in PINNED.items() if got.get(k) != v}
    assert not drift, f"prediction parameters changed (pinned, now): {drift}"


def test_markdown_table_lists_every_row():
    md = params.markdown().splitlines()
    assert len(md) == len(params.PARAMS) + 2
