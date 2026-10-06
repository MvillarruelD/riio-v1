# RIIO 1.0.0 software release audit

Prepared 2026-10-06. This review concerns software distribution and portability; no biological
inference rule, candidate sequence, or reference record was redesigned or recurated.

## Source and scope

- RIO predictor source: `f80c90a` (full commit in `SOURCE_MANIFEST.json`).
- Analysis source: `6363de9f27343e0e39cd4a7cd24f495d10eaa809`.
- Canonical paper predictor baseline: `e9ba665`, tag `v1.0.0-rc`, package/reference release 0.3.0.
- New software version: 1.0.0. Immutable database release remains 0.3.0 with original bytes/hashes.
- Three canonical manifests and all 238 referenced single-protein FASTAs are included.

The clean export retains all runtime-reachable predictor modules and the production analysis/report
dependency closure. Four reference-only discovery scripts, thirteen unused structural helpers,
exploratory/maintenance avenues, raw working materials, generated outputs, engine binaries/weights,
and private environment files are excluded. File provenance is recorded in `SOURCE_MANIFEST.json`.

## Release changes

Owned BITACORA and MetalNet runners, required folding helpers, and the inherited BITACORA patch now
reside in the installable Python package. Installed-wheel verification explicitly checks them.
`METALNET2_DIR` is a location override; no model parameter or inference rule changed. The consolidated
analysis checkout takes precedence over unrelated local predictor checkouts. Canonical run roots are
resolved relative to `runs/` or `RIIO_RUNS_DIR`. Missing discovery diagnostic-cache tables stop scoring.

Benchmark family constants and exact paper accession membership are retained without shipping an
unused database builder or full supplementary sequence duplicates. Source FASTA hashes document the
membership transformation. The report dependency `cnsplots` is now declared. GUI/report GitHub links,
version metadata, citation and installation instructions point to this repository.

The original MIT code notice and separately licensed GPL-3.0-only BITACORA patch are distinguished.
Ligify and metalloregulator-mining full notices are preserved in repository and distribution artifacts.
Manuscript names are acknowledged; final software authorship and paper DOI remain pending author input.

## Verification evidence

- Original predictor and analysis generic suites passed before export.
- Clean Python 3.12 installation with predictor, GUI, genome and analysis dependencies succeeded.
- Clean exported source/analysis suite passed; three skips concern figure-specific test applicability.
- Packaged-data suite: 47 passed initially; four failures were resolved and their four tests passed.
  Two required an explicitly configured MMseqs2 executable, omitted from the clean export by design.
  Two assumed software and immutable database releases must share a version; those checks now validate
  package-version coherence and the independent reference release separately.
- Data manifest is current; original immutable reference bytes are retained. Git attributes preserve
  reference line endings so a fresh clone cannot invalidate recorded hashes.
- Compilation and Ruff checks passed. Report and benchmark imports resolve without the unused external
  database-builder source mirror. Monorepo and relocated-output tests passed.
- `CITATION.cff` validates against the official citation-file-format schema.
- Publication scan found no private-key/GitHub/OpenAI credential patterns or private/generated paths;
  all exported files are below GitHub's 100 MB per-file limit.

Installed-wheel and remote CI results are recorded below when complete.

## Reproducibility limits

The online, multi-hour canonical prediction batch was not rerun. Optional engine/service availability
and current network database responses are not certified by software tests. Full offline recreation
of the paper requires the recorded generated bundles, genome/annotation caches and discovery tables,
which are not code assets and are not distributed here. Optional manuscript artwork is excluded.
The current release follows later RIO report/genome orchestration changes and is not claimed to be
byte-identical to the executable used for the frozen paper batch.
