"""
sequence_checks.py — protein sequence sanity checks before folding a TF monomer.

The ceiling is a COST AND ABSURDITY backstop, not a curation check.

It began as one: when this library only handled ArsR/SmtB (~100–130 aa) and MerR
(~120–145 aa) monomers fetched from UniProt, anything over 150 aa was almost
always the wrong entry — a full-length operon protein instead of the regulator —
and refusing to fold it was the right call.

That premise no longer holds. The pipeline now spans twelve families and takes
its sequences from BITACORA over annotated genomes (`protein_source=
annotated_genome`), not from UniProt lookups, so length is no longer a proxy for
a misfetch. Several families are genuinely long: MocR/GabR-subfamily GntR
regulators carry a C-terminal aminotransferase domain and run 470–500 aa,
LysR-type regulators ~300 aa, TetR/AcrR ~200 aa. Measured on the 155-candidate
production manifest, a 150 aa ceiling rejected 89 of them (57 %) — 30 GntR and
21 TetR among them — every one a real regulator.

The ceiling is therefore set at 500 aa: above the longest genuine candidate in
either production manifest (490 aa) and still low enough to refuse something
absurd before it costs a fold from the daily quota. Curation of implausible
entries is handled upstream, by the family-median length check the production
pre-registration specifies — not here. A wrong protein of ordinary length was
never caught by this guard anyway.

Known documented exceptions are listed in EXEMPT_TFS; they predate the raise and
are kept because they document real architectures.
"""

from __future__ import annotations

from typing import Iterable

MAX_MONOMER_LENGTH_AA = 500

# Known documented exceptions to the 150 aa ceiling:
# CoaR_Synechocystis  — ~365 aa MerR fused to a C-terminal precorrin isomerase
#                       domain (Rutherford 1999 JBC). Fusion is real biology.
# GolS_Styphimurium   — 154 aa; gold-sensor MerR; slightly above ceiling but
#                       manually verified as A0A636FGX4 / WP_058107416.1
#                       (session-6 curation). Not a misfetch.
# BmrR_Bsubtilis      — 281 aa; UniProt P39075 is the genuine BmrR MerR-family
#                       transcriptional activator (NOT BmrA ABC transporter).
#                       Architecture: N-HTH (1-90) + linker (91-120) + C-terminal
#                       drug-binding/dimerization BRC domain (121-279).
#                       PDB 1BOW/1EXI/3Q1M confirm full-length structure.
#                       Zheleznova 1999 Cell; Zheleznova 2000 Nature.
# BluR_Ecoli         — 243 aa; b1162/ycgE, a MerR-family regulator with a
#                       C-terminal extension beyond the canonical ~130-aa core
#                       (Tulin 2023). Real length, not a misfetch.
# CarH_Mxanthus      — 302 aa; B12 (AdoCbl)-dependent MerR-family photoreceptor
#                       with a C-terminal cobalamin-binding domain (Ortiz-Guerrero
#                       2011; Jost 2015). Real length, not a misfetch.
EXEMPT_TFS = frozenset({"CoaR_Synechocystis", "GolS_Styphimurium",
                        "BmrR_Bsubtilis", "BluR_Ecoli", "CarH_Mxanthus"})


class SequenceCheckError(ValueError):
    """Raised when a TF protein sequence fails sanity checks."""


def validate_tf_sequence(
    tf_id: str,
    sequence: str,
    *,
    max_len: int = MAX_MONOMER_LENGTH_AA,
    exempt: Iterable[str] = EXEMPT_TFS,
) -> None:
    """Raise SequenceCheckError if the monomer is implausibly long.

    Parameters
    ----------
    tf_id : str
        Canonical TF id (used to check the exemption list).
    sequence : str
        Amino-acid sequence string.
    max_len : int
        Ceiling above which we refuse to fold (default 500 aa).
    exempt : Iterable[str]
        TF ids that are allowed to exceed max_len.
    """
    if not sequence:
        raise SequenceCheckError(f"{tf_id}: empty protein sequence")

    n = len(sequence)
    if n <= max_len:
        return
    if tf_id in set(exempt):
        return
    raise SequenceCheckError(
        f"{tf_id}: monomer length {n} aa > {max_len} aa ceiling; "
        f"this is an absurdity backstop, so verify the sequence is a regulator "
        f"and not a full-length operon protein before folding"
    )
