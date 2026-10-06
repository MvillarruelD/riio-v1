"""The IPG parser must accept what `Entrez.efetch` actually returns.

REGRESSION. `Entrez.efetch(db="protein", rettype="ipg", retmode="text")` yields BYTES despite the
"text" retmode. `_parse_ipg` split on a str separator, so bytes raised
`TypeError: a bytes-like object is required, not 'str'`.

Two callers existed. `msa_homologs.loci_from_protein` decoded first and worked;
`genome_resolver.resolve_ncbi` passed the handle straight through and did not. Because resolve_ncbi
wraps each accession in `except Exception: continue`, all 25 raised, the candidate list came back
empty, and the remote homolog fallback returned ZERO regions on all 84 invocations of run A --
30 min to 6 h of blastp thrown away each time, and the direct cause of that run's eight timeouts.
Nothing surfaced: an empty result is indistinguishable from "this protein has no homologs".

These tests are offline and use a captured IPG fragment, so they cost nothing and cannot be skipped
by a network outage.
"""
from __future__ import annotations

import warnings

from predictor.annotate import genome_resolver as gr

#: A real IPG response shape: header, then one row per genomic locus of the identical protein.
IPG_TSV = (
    "Id\tSource\tNucleotide Accession\tStart\tStop\tStrand\tProtein\tProtein Name\tOrganism\n"
    "123\tRefSeq\tNZ_ACFI01000060.1\t26209\t26698\t-\tWP_003872472.1\tArsR family\tM. avium\n"
    "123\tRefSeq\tNZ_AGAP01000430.1\t8665\t9154\t+\tWP_003872472.1\tArsR family\tM. avium\n"
)


def test_parse_ipg_accepts_str():
    got = gr._parse_ipg(IPG_TSV)
    assert len(got) == 2
    assert got[0].accession == "NZ_ACFI01000060.1"
    assert got[0].tf_strand == "-"


def test_parse_ipg_accepts_bytes_the_way_efetch_returns_them():
    """The whole defect in one assertion: this raised TypeError and silently cost run A its
    homolog fallback on every no-clade candidate."""
    got = gr._parse_ipg(IPG_TSV.encode("utf-8"))
    assert len(got) == 2, "bytes must parse identically to str"


def test_bytes_and_str_agree_exactly():
    a = gr._parse_ipg(IPG_TSV)
    b = gr._parse_ipg(IPG_TSV.encode("utf-8"))
    assert [(c.accession, c.tf_start, c.tf_end, c.tf_strand) for c in a] == \
           [(c.accession, c.tf_start, c.tf_end, c.tf_strand) for c in b]


def test_undecodable_bytes_do_not_raise():
    """A malformed byte must degrade to a replacement character, not take out the whole lookup --
    one unparseable accession should never cost the other twenty-four."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        got = gr._parse_ipg(IPG_TSV.encode("utf-8") + b"\xff\xfe garbage\n")
    assert len(got) == 2
