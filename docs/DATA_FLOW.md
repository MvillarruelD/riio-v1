# Data flow: what each stage computes, and where it lands in the bundle

One prediction runs `pipeline.run_novel` → `_run_novel_impl`. Every stage below writes into the
dossier (`dossier.json`); nothing is fed back from a later stage into an earlier one.

| # | stage | main function | dossier keys it produces |
|---|---|---|---|
| 1 | family | `annotate.tf_record.classify` → `annotate.family_db` | `family`, `family_support`, `tf_record` |
| 2 | genome + TF locus | `annotate.genome_mirror` / `genome_resolver.resolve`, `annotate.context` | `genome_accession`, `genome_topology`, `tf_locus`, `operator_window`, `neighborhood` |
| 3 | SSN clade + homologue promoters | `annotate.homolog_selection.assign_cluster`, `annotate.msa_homologs`, `annotate.homolog_regions` | `ssn_cluster`, `homolog_regions` (with `shortfall` when under the ortholog floor), `conservation` |
| 4 | inducer / metal | `effector.inducer.infer_inducer` (fusion in `_finalize`, `schema.InducerConsensus`) | `inducers` (every source call with its full evidence), `effector`, `ion`, `metal_context` |
| 5 | seed motifs | `signals.motif_finder.predict` → `signals.seed_selection` | `seed_selection`, `seed_logo`, `literature_seed_candidates` |
| 6 | genome rescan | `signals.motif_rescan.rescan_genome` | `rescan` (every hit: coordinates, score, p, locality tier) |
| 7 | operator motif + ranking | `signals.operator_logo.select_by_motif`, `report.operators.build_operators` | `genomic_logo`, `multi_source_logos`, `operators` (`primary` + ranked candidates) |
| 8 | autoregulation, promoter | `pipeline` (distance to the TF's own window), `signals.promoter.classify_mode` | `autoregulatory_recovery`, `promoter_occlusion` |
| 9 | regulon | `report.ligand_regulon.regulon` → `signals.regulon.reconstruct_regulon` | `regulon`, `regulon_stats`, `per_hit_regulation`, `confirmed_regulated_genes`, `ligand`, `ligand_genes`, `ligand_keywords` |
| 10 | structure (off by default) | `structure.early_structure.run_structure_step` | `structure` |
| 11 | hand-off + record | `structure.af3_export.to_af3_jobs`, `report.outputs.write_stage_a` | `af3_jobs`, `designed_operator`, `run_record`, `caveats`, `provenance` |

Every parameter these stages use is listed, with its value and location, by `tfop params`.

## Citations in the source

Comments in this package cite two things that are not part of it:

- **`handoff §N`** — a section of the project handoff document, which recorded each measurement as it
  was made. It was deleted by commit 5815559 of the analysis repository ("Freeze the workspace for
  publication") and is readable in that repository's history. The section numbers are kept because they
  are how a note is found there.
- **`analysis/...` and `project_archive/...` paths** — files in the analysis repository, which is
  separate from this package and is not distributed with it.

Both are provenance for a stated measurement, never something the code reads. A comment that cites one
also states the number it is citing, so the reasoning stands without the source at hand.

## Rules

1. **Evidence a source computes is either used by its call or reported as unused.** The 2026-09-15
   defect was the opposite: the SSN source recorded a unanimous curated-anchor consensus in its
   evidence and emitted no ligand, because it read only the flat clade label.
   `tests/test_evidence_flow.py` locks the fix and checks that every source call, including
   abstaining ones, reaches the dossier with its evidence intact.
2. **Abstention is explicit.** A source with nothing to say is kept in `inducers.calls` with
   `ligand=None`; it never disappears.
3. **No stage reads a later stage's output.** The regulon is reported alongside the inducer call
   but is not one of its inputs.

## Fields that carry no information in the canonical runs

Measured with `tools/audit_dossier_fields.py` over the 239 dossiers of the 2026-09-17 runs
(structure prediction off). A field is listed when it was empty, or held one value, in every dossier.

| field(s) | observed | why |
|---|---|---|
| `structure.*` (14 fields) | empty or `folded: false` everywhere | folding was off for every run; the block is filled only when `--fold` is used |
| `seed_selection.policy_name`, `.chosen_cluster_id`, `.chosen_why`, `.eligible_cluster_ids`, `.fallback_used`, `.fallback_reason`, `.parameters`, `.construction` | one value | the single production policy restated per bundle |
| `seed_logo.cluster_id` | always `0` | follows from the policy (cluster 0) |
| `designed_operator.method`, `.note` | one value | one design method |
| `operators.design_after_structure`, `.structure_available` | one value | fixed by the structure-off configuration |
| `metal_context.definition` | one value | a definition string repeated per bundle |
| `provenance` | always `PREDICTED_UNVERIFIED` | every output is an unverified prediction |
| `tf_record.is_dna` | always `false` | inputs were proteins |
| `tf_record.uniprot` | always empty | candidates came from genome discovery, not UniProt |
| `genome_topology.origin_wrapped_features` | always empty | no gene spanned a replicon origin in these genomes |

None of these affects a prediction. Whether to keep them (as per-bundle self-description) or drop
them from the output is an open decision, recorded here rather than made silently.
