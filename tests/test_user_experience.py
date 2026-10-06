"""Public onboarding, CLI input, and report-understanding contracts."""
from __future__ import annotations

import pytest

from predictor import resources, run_cli
from predictor.annotate import genome_scan
from predictor.report import render_report


@pytest.mark.released_data
def test_audit_explains_that_the_core_database_is_bundled(capsys):
    info = resources.database_info()
    assert info["present"] is True
    assert info["file_count"] > 400
    assert info["mib"] > 100

    resources._print_audit()
    out = capsys.readouterr().out
    assert "REFERENCE DATABASE: READY" in out
    assert "do not download or build it separately" in out
    assert "Query genomes/promoter context are fetched from NCBI as needed" in out


def test_cli_accepts_raw_protein_without_network():
    seq = "MKKLTVSDLAREHGYEVVETPEQVTAL"
    assert run_cli._resolve_prediction_source(seq, offline=True) == seq


def test_cli_resolves_a_protein_accession(monkeypatch):
    expected = "MKKLTVSDLAREHGYEVVETPEQVTAL"
    monkeypatch.setattr(
        "predictor.annotate.protein_seqs.fetch_protein_sequences",
        lambda ids, allow_ncbi, verbose: {"WP_000240047": expected},
    )
    assert run_cli._resolve_prediction_source("WP_000240047.1") == expected


def test_report_leads_with_a_decision_summary_and_validation_warning():
    dossier = {
        "tf_id": "Synthetic_report",
        "family": "Unspecified",
        "genome_accession": "SYNTHETIC_CONTIG",
        "tf_locus": [100, 430, "+"],
        "inducers": {"top": "Synthetic label", "agreement": 0.67},
        "rescan": {"hits": [{"score": 9.2, "seq": "SYNTHETIC_SEQUENCE"}]},
        "operators": {"candidates": [{
            "source": "conservation_natural", "sequence": "SYNTHETIC_SEQUENCE",
            "score": 9.2, "status": "ready",
        }]},
    }
    html = render_report.build_html({"dossier": dossier, "tf_summary": {}}, {})
    assert "summary-grid" in html
    assert "SYNTHETIC_SEQUENCE" in html
    assert "not experimental validation" in html
    assert "computational prediction" in html
    assert "Current full-report family scope" not in html
    assert "predicted evidence" not in html
    assert 'class="report-nav"' in html
    assert 'id="final"' in html
    assert "viewport" in html
    # The report must not advertise a structural stage that no longer exists, nor tell the reader to
    # run a command that raises ModuleNotFoundError. Removed 2026-09-02 with Stage B.
    for gone in ("DeepPBS", "FoldX", "tfop finalize", "Structural validation", 'id="structure"'):
        assert gone not in html, gone
    # This dossier has no `structure` block at all, i.e. the default no-fold run. The fold panel must be
    # replaced by a statement, not rendered empty with "folded None / mean pLDDT None".
    assert "No structure model recorded" in html
    assert "mean pLDDT" not in html


def test_report_prefers_inducer_display_and_backfills_legacy_dossiers():
    base = {
        "tf_id": "RcnR_like", "family": "CsoR/FrmR", "genome_accession": "SYNTHETIC",
        "operators": {"candidates": []},
    }
    explicit = render_report.build_html({"dossier": {**base, "inducers": {
        "top": "divalent metal (ion unresolved)",
        "top_display": "explicit display label",
        "candidates": ["Co2+", "formaldehyde"],
        "agreement": 0.5,
    }}, "tf_summary": {}}, {})
    assert "explicit display label" in explicit

    legacy = render_report.build_html({"dossier": {**base, "inducers": {
        "top": "divalent metal (ion unresolved)",
        "candidates": ["Co2+", "formaldehyde"],
        "agreement": 0.5,
    }}, "tf_summary": {}}, {})
    assert "inconclusive: Co2+ or formaldehyde" in legacy


def test_report_shows_the_fold_panel_only_when_a_fold_happened():
    base = {
        "tf_id": "Synthetic_report", "family": "Unspecified", "genome_accession": "SYNTHETIC_CONTIG",
        "operators": {"candidates": []},
    }
    folded = render_report.build_html(
        {"dossier": {**base, "structure": {"folded": True, "backend": "esmfold2",
                                           "plddt_mean": 88.4, "state_name": "homodimer"}},
         "tf_summary": {}}, {})
    assert "mean pLDDT" in folded and "88.4" in folded
    assert "No structure model recorded" not in folded

    unfolded = render_report.build_html(
        {"dossier": {**base, "structure": {"folded": False, "backend": None, "plddt_mean": None}},
         "tf_summary": {}}, {})
    assert "No structure model recorded" in unfolded
    assert "mean pLDDT" not in unfolded and "AFDB QC" not in unfolded


def test_report_places_each_gene_table_by_its_source_and_links_bundle_records():
    dossier = {
        "tf_id": "Synthetic_report", "family": "ArsR/SmtB", "genome_accession": "SYNTHETIC",
        "inducers": {"top": "Zn2+", "calls": [{"source": "ssn_cluster", "ligand": "Zn2+",
                                           "confidence": 1.0, "role": "metal", "evidence": {}}]},
        "confirmed_regulated_genes": [{"gene": "zntA", "product": "zinc exporter",
                                       "support_label": "no operator support"}],
        "regulon": [{"first_gene": "otherA", "strand": "+", "genes": ["otherA", "otherB"],
                     "site": [100, 120], "qvalue": 0.2, "tier": 3}],
        "regulon_stats": {"n_operons_found": 43, "truncated": True},
        "operators": {"candidates": []},
    }
    files = {"ligand/inducer.json", "regulation/per_hit_regulation.json", "regulon.tsv",
             "dossier.json", "binding_sites.tsv"}
    report = render_report.build_html({"dossier": dossier, "bundle_files": files}, {})
    genes_heading = "Effector-associated genes, and whether an operator sits at them"
    assert report.index('id="regulon"') < report.index(genes_heading)
    assert report.index(genes_heading) < report.index('id="operators"')
    assert report.index('id="operators"') < report.index("Proposed regulon operons")
    assert 'href="ligand/inducer.json"' in report
    assert 'href="regulation/per_hit_regulation.json"' in report
    assert 'href="regulon.tsv"' in report
    assert "The 1 lowest-q of 43 mapped operons are shown (report cap)." in report
    assert "0 of those shown pass q ≤ 0.05" in report
    assert "No q-value cutoff is applied" in report
    assert "1 shown of 43 mapped operons (report cap)" in report
    # all three support labels are defined, not only the two a reader meets first
    for label in ("operator hit passes q≤0.05", "ranked operator-supported candidate",
                  "no operator support"):
        assert label in report, label
    # the cap is reported from the data, never as a hard-coded number
    assert "25-operon" not in report


def test_report_explains_empty_gene_and_operon_results():
    report = render_report.build_html({"dossier": {"tf_id": "Empty_report"},
                                       "bundle_files": {"dossier.json"}}, {})
    assert "Effector-associated genes" in report
    assert "No effector-associated candidate genes were recorded" in report
    assert "Proposed regulon operons" in report
    assert "No proposed operons were recorded" in report
    assert 'class="visual-stack"' not in report and 'class="visual-grid"' not in report


def test_report_formats_scores_with_units_and_statuses_without_colour_alone():
    dossier = {
        "tf_id": "Synthetic_report", "family": "ArsR/SmtB", "genome_accession": "SYNTHETIC",
        "genomic_logo": {"consensus": "ACGT", "total_ic": 1.984, "width": 4, "n_kept": 13},
        "regulon": [{"first_gene": "geneA", "strand": "+", "genes": ["geneA"], "site": [1, 20],
                     "qvalue": 1.6365858681664355e-10, "tier": 1}],
        "operators": {"candidates": [
            {"source": "conservation_natural", "sequence": "AAAA", "score": 52.12,
             "score_type": "rescan_site_score", "status": "ready",
             "provenance": {"start": 10, "end": 45, "strand": "+", "locality": "proximal"}},
            {"source": "conservation_logo", "sequence": "ACGT", "score": 0.496,
             "score_type": "mean_ic_bits", "status": "ready", "provenance": {}}]},
    }
    report = render_report.build_html({"dossier": dossier, "bundle_files": set()}, {})
    assert "1.6e-10" in report and "1.6365858681664355e-10" not in report
    assert "site score, bits" in report and "mean IC, bits/col" in report
    # the rationale is stated once per category, and no longer claims conservation support
    assert report.count("A genome site from the PWM rescan") == 1
    assert "supported by phylogenetic conservation" not in report
    # 0.496 bits/col is below the 0.5 threshold and must not print as "0.50"
    assert "reported motif 0.49 bits/col" in report
    # status is carried by a glyph and hidden text as well as colour
    assert 'class="mark"' in report and "(weak or absent)" in report


def test_scan_report_matches_the_decision_first_visual_system(tmp_path):
    candidate = genome_scan.Candidate(
        protein_id="WP_000240047.1", family="ArsR/SmtB", pfam_acc="PF01022",
        pfam_name="HATPase_c", evalue=1e-20, bits=88.0, length=110,
        contig="NC_007795.1", start=100, end=430, strand="+", locus_tag="czrA",
        coordinate_source="source GFF3",
    )
    report = tmp_path / "SCAN_REPORT.html"
    genome_scan.write_scan_report([candidate], report, kind="annotated_genome",
                                  n_proteins=4300, source="GCF_000005845.2")
    html = report.read_text(encoding="utf-8")

    assert "Genome TF census" in html
    assert "Scope and interpretation" in html
    assert "The census screens 22 families" in html
    assert "github.com/MvillarruelD/riio-v1" in html
    assert "not experimental" in html
    assert "table-wrap" in html
    assert "viewport" in html
    assert "1.00e-20" in html
    assert "Expected number of equally strong HMM matches" in html
    assert "HMM alignment score in bits" in html
    assert "NC_007795.1:100-430 (+)" in html
    assert "family can include metalloregulators" in html


def test_scan_report_reprises_family_view_for_inducers(tmp_path):
    candidate = genome_scan.Candidate(
        protein_id="WP_1", family="ArsR/SmtB", pfam_acc="PF01022", pfam_name="ArsR",
        evalue=2e-18, bits=72.4, length=105,
    )
    report = tmp_path / "SCAN_REPORT.html"
    genome_scan.write_scan_report(
        [candidate], report, predictions={"WP_1": {
            "inducer": "Zn2+", "primary_operator": "ATATGCAT", "bundle": None,
            "status": "Complete",
        }},
    )
    html = report.read_text(encoding="utf-8")
    assert "Possible inducer distribution" in html
    assert "Zn2+" in html
    assert "candidate_predictions.tsv" in html


def test_report_title_wraps_at_run_name_fields_within_the_html_contract():
    from predictor.report.bundle_contract import _Resources
    report = render_report.build_html({"dossier": {"tf_id": "GCF_000005845.2__ArsR__NP_417153.1"},
                                       "bundle_files": set()}, {})
    assert '<span class="nb">GCF_000005845.2__</span><span class="nb">ArsR__</span>' in report
    parser = _Resources()
    parser.feed(report)
    assert parser.issues == []
