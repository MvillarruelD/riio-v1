"""`tfop scan` must use a genome's own annotation when it has one.

Pyrodigal invents protein ids and picks its own start codons, so a re-called genome cannot be joined
to anything keyed on accession. Measured on the survey's four genomes, the accession match rate
against the published 150 candidates is 149/150 from the annotation and 0/150 from a re-call
(`analysis/discovery_bench/DISCOVERY_DECISION_20260819.md`). The regulon, the neighbourhood flags and
every join to the collaborator's tables are keyed on accession, so this is a correctness property,
not a preference.
"""
from __future__ import annotations

import json

import pytest

from predictor.annotate import genome_scan as GS

# Two CDS on one contig: one forward, one reverse. Padded so the genome is long enough to be
# recognisably DNA and the coordinates are not degenerate.
FWD = "ATG" + "GCTAGCAAAGAACTGATT" + "TAA"          # 24 nt -> 7 aa
REV = "ATG" + "CATCACCATCACGAATTC" + "TGA"          # 24 nt -> 7 aa
PAD = "ACGT" * 30


def _revcomp(s: str) -> str:
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def _write_genome(tmp_path, with_gff: bool):
    contig = PAD + FWD + PAD + _revcomp(REV) + PAD
    fwd_start = len(PAD) + 1
    fwd_end = fwd_start + len(FWD) - 1
    rev_start = len(PAD) + len(FWD) + len(PAD) + 1
    rev_end = rev_start + len(REV) - 1

    fna = tmp_path / "tinygenome.fna"
    fna.write_text(f">contig1 test\n{contig}\n", encoding="utf-8")

    if with_gff:
        gff = tmp_path / "tinygenome.gff"
        rows = [
            f"contig1\ttest\tCDS\t{fwd_start}\t{fwd_end}\t.\t+\t0\tID=cds-1;protein_id=WP_000000001.1",
            f"contig1\ttest\tCDS\t{rev_start}\t{rev_end}\t.\t-\t0\tID=cds-2;protein_id=WP_000000002.1",
        ]
        gff.write_text("##gff-version 3\n" + "\n".join(rows) + "\n", encoding="utf-8")
    return fna


def test_annotated_genome_route_preserves_accessions(tmp_path):
    fna = _write_genome(tmp_path, with_gff=True)
    recs, kind = GS.proteins_from_input(str(fna), verbose=False)

    assert kind == "annotated_genome", "a genome with a sibling .gff must not be gene-called"
    ids = {pid for pid, _ in recs}
    assert ids == {"WP_000000001.1", "WP_000000002.1"}, (
        "the annotation's own protein accessions are the join key and must survive verbatim"
    )


def test_annotated_route_translates_both_strands(tmp_path):
    fna = _write_genome(tmp_path, with_gff=True)
    recs, _kind = GS.proteins_from_input(str(fna), verbose=False)
    by_id = dict(recs)

    # to_stop=True, so the terminal stop codon is not part of the protein
    assert by_id["WP_000000001.1"] == "MASKELI"
    assert by_id["WP_000000002.1"] == "MHHHHEF", "a minus-strand CDS must be reverse-complemented first"


def test_annotated_route_can_return_locus_metadata(tmp_path):
    fna = _write_genome(tmp_path, with_gff=True)
    recs, kind, metadata = GS.proteins_from_input(str(fna), verbose=False, return_metadata=True)

    assert kind == "annotated_genome"
    assert len(recs) == 2
    assert metadata["WP_000000001.1"]["contig"] == "contig1"
    assert metadata["WP_000000001.1"]["start"] == len(PAD) + 1
    assert metadata["WP_000000002.1"]["strand"] == "-"
    assert metadata["WP_000000001.1"]["coordinate_source"] == "source GFF3"


def test_unannotated_genome_still_falls_back_to_a_gene_call(tmp_path):
    pytest.importorskip("pyrodigal")
    fna = _write_genome(tmp_path, with_gff=False)
    _recs, kind = GS.proteins_from_input(str(fna), verbose=False)

    assert kind == "genome", "without an annotation the gene call is the correct fallback"


def test_a_malformed_annotation_degrades_instead_of_raising(tmp_path):
    pytest.importorskip("pyrodigal")
    fna = _write_genome(tmp_path, with_gff=True)
    (tmp_path / "tinygenome.gff").write_text("##gff-version 3\nnot\ta\tvalid\trow\n", encoding="utf-8")

    _recs, kind = GS.proteins_from_input(str(fna), verbose=False)
    assert kind == "genome", "a broken GFF must fall back, not sink the scan"


def test_multi_segment_cds_is_joined_not_truncated(tmp_path):
    """One protein_id spread over several GFF lines yields the whole protein, not its first segment."""
    contig = PAD + FWD + PAD
    start = len(PAD) + 1
    fna = tmp_path / "split.fna"
    fna.write_text(f">contig1\n{contig}\n", encoding="utf-8")
    gff = tmp_path / "split.gff"
    gff.write_text(
        "##gff-version 3\n"
        f"contig1\tt\tCDS\t{start}\t{start + 11}\t.\t+\t0\tID=a;protein_id=WP_1.1\n"
        f"contig1\tt\tCDS\t{start + 12}\t{start + len(FWD) - 1}\t.\t+\t0\tID=a;protein_id=WP_1.1\n",
        encoding="utf-8",
    )
    recs, kind = GS.proteins_from_input(str(fna), verbose=False)

    assert kind == "annotated_genome"
    assert dict(recs)["WP_1.1"] == "MASKELI", "segments must be joined before translation"


def test_scan_outputs_include_gff3_and_bed_locus_overlays(tmp_path):
    source = _write_genome(tmp_path, with_gff=True)
    candidate = GS.Candidate(
        protein_id="WP_000000001.1", family="ArsR/SmtB", pfam_acc="PF01022",
        pfam_name="HATPase_c", evalue=1e-20, bits=88.0, length=7, sequence="MASKELI",
        contig="contig1", start=121, end=144, strand="+", locus_tag="gene_1",
        coordinate_source="source GFF3",
    )
    out = tmp_path / "scan"
    GS.write_outputs([candidate], out, kind="annotated_genome", n_proteins=2,
                     source="GCF_test", source_path=source,
                     annotation_path=source.with_suffix(".gff"))

    gff = (out / "candidate_tfs.gff3").read_text(encoding="utf-8")
    bed = (out / "candidate_tfs.bed").read_text(encoding="utf-8")
    tsv = (out / "candidate_tfs.tsv").read_text(encoding="utf-8")
    assert "contig1\ttfop\tgene\t121\t144" in gff
    assert "tf_family=ArsR%2FSmtB" in gff
    assert "contig1\t120\t144\tgene_1" in bed, "BED coordinates must be 0-based, half-open"
    assert "locus_tag\tgene\tcontig\tstart\tend\tstrand" in tsv
    assert (out / "source_genome.fna").is_file()
    assert (out / "source_annotation.gff3").is_file()


def test_regenerated_scan_manifest_keeps_copied_source_files(tmp_path):
    """Adding predictions must not make existing source inputs disappear from the manifest."""
    source = _write_genome(tmp_path, with_gff=True)
    candidate = GS.Candidate(
        protein_id="WP_000000001.1", family="ArsR/SmtB", pfam_acc="PF01022",
        pfam_name="HATPase_c", evalue=1e-20, bits=88.0, length=7, sequence="MASKELI",
        contig="contig1", start=121, end=144, strand="+", locus_tag="gene_1",
        coordinate_source="source GFF3",
    )
    out = tmp_path / "scan"
    GS.write_outputs([candidate], out, kind="annotated_genome", n_proteins=2,
                     source="GCF_test", source_path=source,
                     annotation_path=source.with_suffix(".gff"))

    GS.write_outputs(
        [candidate], out, kind="annotated_genome", n_proteins=2, source="GCF_test",
        predictions={candidate.protein_id: {"status": "complete", "inducer": "Zn(II)"}},
    )

    files = json.loads((out / "scan_manifest.json").read_text(encoding="utf-8"))["files"]
    assert files["source_sequence"] == "source_genome.fna"
    assert files["source_annotation"] == "source_annotation.gff3"
