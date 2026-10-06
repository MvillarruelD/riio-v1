# Publication repository v1: prompt and execution plan

## Rewritten request

Audit the software that implements the current RIO pipeline, then create a clean GitHub repository for its first publication release. Include only the active predictor, its required reference assets, the production analysis drivers, installation and usage documentation, meaningful tests, and release tooling. Exclude abandoned experiments, generated results, credentials, machine-local environments, and historical working materials. Preserve the existing inference behavior and biological content. Document source commits, external dependencies, limitations of reproducibility, licenses, upstream attribution, and verified contributor roles. Give the project an Argentina-inspired name that covers regulator, inducer, and operator inference, and publish a versioned release after the software checks pass. Do not invent paper authorship, permissions, or test results.

## Phase 0 — establish the release boundary

- Read both repository READMEs, `PRODUCTION_RUNBOOK.md`, `analysis/run_pipeline.py`, `analysis/canonical_runs.toml`, predictor `pyproject.toml`, `LICENSE`, `docs/TESTING.md`, and `tools/release_check.py`.
- Use documented `tfop predict`, `tfop scan`, `tfop genome`, `tfop audit`, `tfop setup`, and `tfop selftest` interfaces. Preserve `predictor` imports and `tfop` commands to avoid a behavior-changing rename.
- Trace the production stage graph and required input files; identify upstream code/data notices and verified Git contributors.
- Record the paper baseline separately from current source: canonical runs used predictor `e9ba665` (`v1.0.0-rc`); current source is a later commit.
- Verification: clean source trees, exact commits, explicit inclusion manifest; no inferred contributor roles or licenses.

## Phase 1 — final software review

- Run existing predictor and analysis contract suites and packaged-reference checks. Review source-to-wheel resource resolution, external engine adapters, and paths used by the retained production scripts.
- Restrict changes to release packaging, portability, attribution, and demonstrable software defects. Reference the existing release-check pattern rather than inventing another prediction workflow.
- Verification: record actual results and distinguish offline software checks from a full prediction rerun; preserve canonical reference bytes.

## Phase 2 — construct a clean, independent repository

- Export an explicit allowlist into a new checkout outside the research workspace. Copy only tracked predictor runtime sources/reference files, relevant tests, active engine adapters, and the production analysis dependency closure.
- Include the three canonical candidate manifests and only their referenced protein inputs. Exclude caches, `.env`, generated output, binaries/model weights, experimental studies, database construction utilities, and research history.
- Add installation/CLI/GUI instructions, portable production commands, source provenance, contribution guidance, citation metadata, third-party notices, and a release audit.
- Make analysis resolve the predictor in the new repository before other local checkouts. Preserve upstream notices in both repository and installable artifacts.
- Verification: allowlist and hash audit; no accidental private files; only required reference assets; run retained tests inside the clean checkout.

## Phase 3 — verify installable artifacts

- Follow `tools/release_check.py`: build wheel and sdist, install the wheel non-editably into a fresh environment, run CLI/data smoke checks outside the source tree, and verify every manifest asset survived packaging.
- Add CI for source tests, released-data checks, analysis contracts, and isolated wheel installation.
- Verification: wheel and sdist exist; packaged-reference hashes match source; caches/results remain outside the installed package; report optional-engine skips explicitly.

## Phase 4 — GitHub and v1 release

- Create a new GitHub repository under the authenticated account with a fresh history and documented source provenance. Provisional name: RIIO (Regulator, Inducer & Operator inference), a nod to Río de la Plata; repository `riio-v1`.
- User decision: publish `MvillarruelD/riio-v1` publicly, with MIT for original code and preserved upstream licenses. Record the supplied manuscript author names; final software author roles and ORCIDs remain for the authors to update.
- Push reviewed files, check remote CI, tag `v1.0.0`, and attach wheel/sdist and audit information to the release once required checks pass.
- Verification: remote commit matches local commit; visibility/license accurately described; release assets downloadable; paper citation instructions do not claim a DOI that does not exist.

## Scope and evidence limits

No biological redesign, database recuration, manuscript rewriting, or hours-long production prediction rerun is part of this software release. Existing canonical outputs are not copied into the code repository. Optional third-party engines and online services must be documented; an offline smoke test cannot certify their availability or reproduce changing network responses.
