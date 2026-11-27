"""An alternative start codon is translated as methionine when it initiates translation.

In the bacterial code (NCBI translation table 11) GTG, TTG, CTG, ATT, ATC and ATA initiate besides ATG,
and the initiator is read as Met. A codon-by-codon translation gives Val for a GTG start and Leu for a
TTG start: until this fix the GFF route did that, and 54 of the 238 canonical query proteins began with
V or L where the NCBI protein record has M (for example E. coli MarR, NP_416047.4, VKSTSDLF... against
MKSTSDLF...). A CDS whose 5' end is missing keeps its literal first codon, which is not a start then.
"""
from __future__ import annotations

from predictor.annotate import genome_scan as GS

GTG_CDS = "GTG" + "AAAGCTCTGATT" + "TAA"      # V K A L I -> M K A L I
TTG_CDS = "TTG" + "CATCACGAATTC" + "TGA"      # L H H E F -> M H H E F
PAD = "ACGT" * 30


def _revcomp(s: str) -> str:
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def test_translate_cds_reads_alternative_starts_as_met():
    assert GS.translate_cds(GTG_CDS) == "MKALI"
    assert GS.translate_cds(TTG_CDS) == "MHHEF"
    assert GS.translate_cds("ATG" + "AAAGCT" + "TAA") == "MKA"


def test_translate_cds_keeps_the_first_codon_of_a_five_prime_partial_cds():
    assert GS.translate_cds(GTG_CDS, five_prime_complete=False) == "VKALI"
    # a first codon that is not a start codon is never rewritten
    assert GS.translate_cds("GCT" + "AAA" + "TAA") == "AK"


def _genome(tmp_path, extra_attrs: str = ""):
    contig = PAD + GTG_CDS + PAD + _revcomp(TTG_CDS) + PAD
    f0 = len(PAD) + 1
    f1 = f0 + len(GTG_CDS) - 1
    r0 = f1 + len(PAD) + 1
    r1 = r0 + len(TTG_CDS) - 1
    fna = tmp_path / "g.fna"
    fna.write_text(f">contig1\n{contig}\n", encoding="utf-8")
    gff = tmp_path / "g.gff"
    gff.write_text("##gff-version 3\n"
                   f"contig1\tt\tCDS\t{f0}\t{f1}\t.\t+\t0\tID=cds-1;protein_id=P1{extra_attrs}\n"
                   f"contig1\tt\tCDS\t{r0}\t{r1}\t.\t-\t0\tID=cds-2;protein_id=P2\n", encoding="utf-8")
    return fna, gff


def test_gff_route_translates_gtg_and_ttg_starts_as_met_on_both_strands(tmp_path):
    fna, gff = _genome(tmp_path)
    prots = dict(GS.proteins_from_gff(fna, gff))
    assert prots == {"P1": "MKALI", "P2": "MHHEF"}


def test_gff_route_keeps_the_first_codon_when_ncbi_marks_the_five_prime_end_missing(tmp_path):
    fna, gff = _genome(tmp_path, extra_attrs=";partial=true;start_range=.,121")
    prots = dict(GS.proteins_from_gff(fna, gff))
    assert prots["P1"] == "VKALI"
    assert prots["P2"] == "MHHEF"


def test_genbank_route_translates_an_untranslated_gtg_cds_as_met(tmp_path):
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO

    from predictor import genome_run

    contig = PAD + GTG_CDS + PAD
    rec = SeqRecord(Seq(contig), id="contig1", name="contig1", description="test",
                    annotations={"molecule_type": "DNA", "topology": "linear"})
    loc = FeatureLocation(len(PAD), len(PAD) + len(GTG_CDS), strand=1)
    rec.features.append(SeqFeature(loc, type="CDS", qualifiers={"protein_id": ["P1"], "transl_table": ["11"]}))
    gbk = tmp_path / "g.gbk"
    SeqIO.write([rec], str(gbk), "genbank")

    genome_run.prepare_local(None, None, gbk, None, tmp_path / "out")
    sidecar = (tmp_path / "out" / "inputs" / "genome.proteins.faa").read_text(encoding="utf-8")
    assert ">P1\nMKALI\n" in sidecar
