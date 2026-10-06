"""Evidence that is computed must reach the call and the dossier (docs/DATA_FLOW.md).

The defect these lock: on 2026-09-15 the SSN source recorded a clean curated-anchor consensus in its
evidence while emitting no ligand, because it read only the flat cluster label -- so two ArsR_c3
regulators were reported as "non-metal (undetermined)" with the answer (As(III)) sitting in their own
evidence dict. Fixed in 418f5d6; untested until now.
"""
from __future__ import annotations

from dataclasses import asdict

import pytest

from predictor.annotate import homolog_selection as hsel
from predictor.annotate import ssn_clusters as ssn
from predictor.effector import inducer as ind
from predictor.schema import InducerCall, InducerConsensus


def _assigned_to(cluster_id: str):
    info = next(c for c in ssn.load_clusters() if c.cluster_id == cluster_id)
    return hsel.ClusterAssignment(cluster_id=cluster_id, support=0.95, n_hits=10, top_identity=0.8,
                                  info=info, method="mmseqs")


@pytest.mark.released_data
def test_blank_cluster_label_is_filled_by_a_unanimous_anchor_set(monkeypatch):
    assert _assigned_to("ArsR_c3").info.inducer in (None, ""), "precondition: ArsR_c3 has no flat label"
    monkeypatch.setattr(hsel, "assign_cluster", lambda seq, family, **k: _assigned_to("ArsR_c3"))
    monkeypatch.setattr(ind, "curated_member_ion", lambda *a, **k: None)
    cons = ind.infer_inducer("M" + "A" * 100, "ArsR/SmtB", metalnet=False)
    ssn_call = next(c for c in cons.calls if c.source == "ssn_cluster")
    assert ssn_call.ligand == "As(III)", ssn_call
    assert ssn_call.evidence.get("ligand_from") == "anchor_consensus", ssn_call.evidence


def test_every_call_and_its_evidence_survive_serialisation():
    calls = [InducerCall("ssn_cluster", "Zn2+", 0.9, "metal", evidence={"cluster": "Fur_ecZur", "x": 1}),
             InducerCall("ligify", None, 0.0, None, evidence={"reason": "no operon enzymes"}),
             InducerCall("coordination", "divalent metal", 0.6, "metal", evidence={"motif": "CXXC"})]
    d = asdict(InducerConsensus.from_calls(calls, coordination_gate=True))
    assert [c["source"] for c in d["calls"]] == [c.source for c in calls], "an abstaining source was dropped"
    for sent, got in zip(calls, d["calls"]):
        assert got["evidence"] == sent.evidence, f"{sent.source}: evidence lost in the dossier"
