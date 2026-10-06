"""Regression tests for correctness issues found in the 2026-07 whole-project review."""

from types import SimpleNamespace as NS

import numpy as np


def test_confident_hmm_call_does_not_require_blast(monkeypatch):
    from predictor.annotate import family_db

    monkeypatch.setattr(family_db, "hmm_classify",
                        lambda _p: ("PF00001", "TF_domain", "MerR", 1e-30, 100.0))
    monkeypatch.setattr(family_db, "blastp_classify",
                        lambda *_a, **_k: (_ for _ in ()).throw(FileNotFoundError("blastp")))
    monkeypatch.setattr(family_db, "_structure_axis", lambda *_a, **_k: None)
    call = family_db.classify_family("M" * 120)
    assert call.family == "MerR" and call.method == "hmm"
    assert call.nearest_identity is None
    assert any("blastp corroboration unavailable" in flag for flag in call.flags)


def test_tf_record_keeps_identity_distinct_from_family_confidence(monkeypatch):
    from predictor.annotate import family_db, tf_record

    call = family_db.FamilyCall(
        family="MerR", method="hmm", confidence=0.93, nearest_tf="ref",
        nearest_identity=0.41, nearest_evalue=1e-12,
    )
    monkeypatch.setattr(family_db, "classify_family", lambda *_a, **_k: call)
    monkeypatch.setattr(family_db, "_ALL_JSON", NS(exists=lambda: False))
    rec = tf_record.classify("M" * 120)
    assert rec.nearest_identity == 0.41
    assert rec.nearest_evalue == 1e-12
    assert rec.nearest_identity != call.confidence


def test_nearest_reference_uniprot_is_not_query_identity(monkeypatch):
    from predictor.annotate import family_db, tf_record

    call = family_db.FamilyCall(family="MerR", method="hmm", nearest_tf="ReferenceTF")
    monkeypatch.setattr(family_db, "classify_family", lambda *_a, **_k: call)
    monkeypatch.setattr(
        family_db, "_ALL_JSON",
        NS(exists=lambda: True, read_text=lambda **_k: '{"ReferenceTF":{"uniprot":"P12345"}}'),
    )
    rec = tf_record.classify("M" * 120)
    assert rec.uniprot is None
    assert rec.nearest_uniprot == "P12345"


def test_apo_builder_honors_requested_stoichiometry(monkeypatch):
    from predictor._engines.esmfold_lib._lib import esmfold2_inputs

    class Protein:
        def __init__(self, **kwargs):
            self.kwargs = kwargs

    class Structure:
        def __init__(self, **kwargs):
            self.kwargs = kwargs

    monkeypatch.setattr(
        esmfold2_inputs, "_sdk",
        lambda: NS(ProteinInput=Protein, StructurePredictionInput=Structure),
    )
    built = esmfold2_inputs.build_apo_homomer_input("CsoR", "M" * 100, copies=4)
    assert built.chain_ids == ["A", "B", "C", "D"]
    assert built.structure_input.kwargs["sequences"][0].kwargs["id"] == built.chain_ids


def test_family_seed_can_be_ablated(monkeypatch):
    from predictor.signals import motif_finder

    monkeypatch.setattr(
        motif_finder, "family_seeded",
        lambda *_a, **_k: (_ for _ in ()).throw(AssertionError("family seed was called")),
    )
    motif_finder.predict("ACGT" * 30, family="MerR", use_family_seed=False)


def test_shared_geometry_is_the_default(monkeypatch):
    from predictor.signals import motif_finder

    seen = []
    original = motif_finder.find_inverted_repeats

    def capture(seq, **kwargs):
        seen.append(kwargs)
        return original(seq, **kwargs)

    monkeypatch.setattr(motif_finder, "find_inverted_repeats", capture)
    motif_finder.predict("ACGT" * 30, family="MerR", use_family_seed=False)
    assert seen[0]["spacer_range"] == motif_finder.DEFAULT_PRIOR["spacer"]


def test_family_operator_pool_can_leave_query_tf_out(monkeypatch):
    from predictor.annotate import family_db, known_operators

    monkeypatch.setattr(
        known_operators, "load_index",
        lambda: {
            "by_tf": {"CueR_Ecoli": ["ACGTACGT"], "ZntR_Ecoli": ["TGCATGCA"]},
            "by_family": {"MerR": ["ACGTACGT", "TGCATGCA"]},
        },
    )
    monkeypatch.setattr(
        family_db, "_ALL_JSON",
        NS(
            exists=lambda: True,
            read_text=lambda **_k: (
                '{"CueR_Ecoli":{"family":"MerR","uniprot":"P0A9G4"},'
                '"ZntR_Ecoli":{"family":"MerR","uniprot":"P0ACS5"}}'
            ),
        ),
    )
    assert known_operators.family_operator_sequences("MerR", exclude_tf="CueR_bench") == ["TGCATGCA"]


def test_repeat_search_includes_last_valid_window():
    from predictor.signals.motif_finder import find_direct_repeats, find_inverted_repeats

    # With L == 2h, the only valid zero-spacer window starts at zero.
    assert find_direct_repeats("ACGTACGT", half_range=(4, 4), spacer_range=(0, 0), min_frac=1.0)
    assert find_inverted_repeats("ACGTACGT", half_range=(4, 4), spacer_range=(0, 0), min_frac=1.0)


def test_seeded_minus_strand_sequence_is_in_reading_orientation(monkeypatch):
    from predictor.signals import motif_finder

    hit = NS(start=0, end=6, matched="AAATTC", strand="-", pvalue=1e-6)
    monkeypatch.setattr(motif_finder, "scan", lambda *_a, **_k: [hit])
    got = motif_finder.pwm_seeded("AAATTC", np.full((4, 6), 0.25))
    assert got[0].seq == "GAATTT"


def test_regulator_can_be_resolved_by_protein_id():
    from predictor.signals.motif_rescan import Gene, _resolve_regulator

    gene = Gene(10, 100, "+", "locus_tag", protein_id="WP_123")
    assert _resolve_regulator([gene], "WP_123") is gene


def test_genbank_parser_includes_ncrna_features():
    from predictor.annotate.context import from_genbank

    gb = """LOCUS       SYN                      100 bp    DNA     linear   BCT 01-JAN-2000
DEFINITION  synthetic.
ACCESSION   SYN
VERSION     SYN.1
FEATURES             Location/Qualifiers
     source          1..100
     tRNA            10..30
                     /gene="trnX"
     CDS             40..90
                     /gene="cdsA"
                     /protein_id="WP_1"
ORIGIN
        1 aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa
       61 aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa
//
"""
    ctx = from_genbank(gb)
    assert {g.name for g in ctx.genes} == {"trnX", "cdsA"}


def test_canonical_rescan_operator_retains_p_and_q_values():
    from predictor.signals.coords import from_genome_hit
    from predictor.signals.motif_rescan import GenomeHit

    hit = GenomeHit("SYN", 10, 16, "+", 13, 5.0, 1e-6, 2e-4, 1, 1.0, 5.0, "AAAAAA")
    op = from_genome_hit(hit)
    assert op.pvalue == 1e-6
    assert op.qvalue == 2e-4


def test_organic_ligand_is_not_misread_as_cobalt():
    from predictor.report.ligand_regulon import _element_key

    assert _element_key("Co2+") == "CO"
    assert _element_key("cobalt") == "CO"
    assert _element_key("carbon monoxide") is None


def test_ranked_operator_support_is_not_called_confirmation():
    from predictor.report.ligand_regulon import confirm_regulated_genes

    rec = {
        "operator": {"score": 9.0, "qvalue": 0.8},
        "regulated_operon": ["target"],
        "regulated_gene": "target",
        "mode": "promoter-core-overlap",
    }
    out = confirm_regulated_genes([rec], [{"gene": "target"}])
    assert out[0]["operator_supported"] is True
    assert out[0]["confirmed"] is False
    assert "ranked" in out[0]["support_label"]


def test_conservation_does_not_mutate_qvalue():
    from predictor.signals.motif_rescan import Gene, GenomeContext, GenomeHit
    from predictor.signals.regulon import reconstruct_regulon

    ctx = GenomeContext("SYN", "A" * 500, [Gene(100, 200, "+", "target")])
    hit = GenomeHit("SYN", 70, 80, "+", 75, 8.0, 1e-5, 0.02, 1, 1.0, 8.0, "A" * 10)
    result = reconstruct_regulon(ctx, hits=[hit], conserved_targets={"target"})
    assert result.operons[0].qvalue == 0.02
    assert result.operons[0].conserved_target is True
