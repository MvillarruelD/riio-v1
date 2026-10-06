"""Current per-TF data assembled by :func:`predictor.pipeline.run_novel`.

The dataclass contains only fields populated by the live conservation-first pipeline. Operator ranking,
design, structure, provenance and traceability fields that are computed later are attached to the plain
dictionary in ``run_novel`` before publication.
"""
from __future__ import annotations

from dataclasses import dataclass, field

@dataclass
class TFDossier:
    tf_id: str
    family: str
    genome_accession: str | None = None
    tf_locus: tuple | None = None
    effector: str | None = None
    ion: str | None = None
    operator_window: str = ""
    rescan: dict = field(default_factory=dict)            # {region_len, plant_center, n_hits, hits:[...]}
    neighborhood: dict = field(default_factory=dict)      # real-genome view: {center, window, genes, operator, tss}
    genomic_logo: dict = field(default_factory=dict)      # logo from operators FOUND IN THE GENOME (motif-emergence trimmed)
    #: The SEED PWM -- the matrix built from the discovered candidates, which SCANS the genome and so
    #: decides which hits exist at all. Distinct from `genomic_logo`, which is rebuilt afterwards from
    #: the hits it found. The two are routinely conflated and have very different support: the seed is
    #: typically 1-4 sequences, the genomic logo up to 40. Both are reported, each with its own n.
    seed_logo: dict = field(default_factory=dict)
    ligand: list = field(default_factory=list)            # inferred cognate ligand(s) (Ligify)
    ligand_genes: list = field(default_factory=list)      # ligand-associated genes -> extra candidate regions
    ligand_keywords: list = field(default_factory=list)
    regulon: list = field(default_factory=list)           # operons the TF likely drives (operator PWM scan)
    #: {n_operons_found, n_operons_reported, truncated, n_sites}. `regulon` above is capped at
    #: `signals.regulon.DEFAULT_MAX_OPERONS`, a cap nearly every candidate reaches, so its length reads
    #: as a measured regulon size when it is the cap. This says which of the two a reader is looking at.
    regulon_stats: dict = field(default_factory=dict)
    af3_jobs: list = field(default_factory=list)
    #: Metal-context evidence that is computed but has no other home: the survey's own "effective
    #: candidate" definition (a predicted metal SITE plus a metal-homeostasis gene within 500 bp), the
    #: nearest metal neighbour, and how well-anchored the assigned SSN cluster's label is. Without
    #: this the report cannot reproduce the published analysis it is meant to be compared against.
    metal_context: dict = field(default_factory=dict)
    caveats: list = field(default_factory=list)
    provenance: str = "PREDICTED_UNVERIFIED"
    # SSN-cluster routing + unified inducer inference.
    ssn_cluster: str | None = None                       # SSN isofunctional cluster, e.g. "MerR_9"
    inducers: dict = field(default_factory=dict)         # serialized InducerConsensus: {top, agreement, coordination_gate, calls:[...]}
