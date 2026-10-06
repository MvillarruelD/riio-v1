# Packaged reference data

This directory is the immutable data layer used by `tf-operator-predictor`. It is included in source
distributions and wheels; runtime code must never write here. Machine-local derived files belong under
the cache reported by `tfop audit`, while prediction bundles belong under `PREDICTOR_OUTPUT_DIR` (or
`./results` by default).

## What an installed user downloads

The package installer downloads this database together with the Python code. There is no separate
database command, build step, or first-run reference download. In release 0.3.0 the complete wheel is
about 34 MiB compressed and this immutable data tree is about 141 MiB after installation.

Query-specific genomes and homolog promoter regions are different: the pipeline retrieves them from
NCBI on demand and stores reusable copies in the writable cache. It does not ask each user to download a
global genome mirror.

Two very large experimental resources are intentionally outside the package: HomoDB (~13.7 GB) and the
MetalNet2 model stack (~6 GB). They add optional structural evidence and are not prerequisites for a
standard prediction. `tfop audit` reports the core database; `tfop setup` reports the tools and extras.

## Layout

- `refs/` — the Pfam HMM library and family map, broad and curated TF FASTA/metadata pairs, curated
  operator index, anchor knowledge base, and precomputed family PWMs.
- `ssn/clusters/` — 71 cluster FASTAs spanning the 12 supported SSN families.
- `ssn/database/` — runtime indexes derived from the cluster FASTAs: accession-to-cluster lookup,
  family member search FASTAs, cluster inducer table, and curated binder evidence.
- `schemas/` — JSON Schemas for the three databases most useful to downstream work.
- `MANIFEST.json` — file sizes, SHA-256 digests, record counts, and database release metadata. Regenerate
  it with `python tools/build_data_manifest.py` after an intentional data update.

## Stability and provenance

`MANIFEST.json` is the database release boundary. A production result should record the package version;
the manifest digest then identifies the exact reference bytes behind it. Family labels use the canonical
strings exposed by `predictor.annotate.ssn_clusters.FAMILIES`.

The JSON files are ordinary UTF-8 and the FASTAs are standard text, so they can be used independently of
the pipeline. Programmatic users should prefer the read-only loaders in `predictor.annotate.family_kb`,
`known_operators`, and `ssn_clusters`, which preserve normalization and overlap-resolution rules.

## Updating the database

Database construction and benchmarking are deliberately outside the runtime path. Build and validate a
candidate release in the research workspace, run the full test suite and database integrity checks, copy
only the accepted immutable artifacts here, regenerate `MANIFEST.json`, and review the data diff before
release. Never promote HTTP caches, holdout FASTAs, logs, structure predictions, or benchmark outputs.
