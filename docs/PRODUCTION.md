# Production and canonical reproduction

Use a Git checkout/source archive for the analysis companion. Install its dependencies:

```bash
python -m pip install -e ".[dev,genome,gui]" -r analysis/requirements.txt
```

The three supplied canonical manifests contain 39 E. coli, 45 Salmonella, and 154 survey candidates.
Their 238 FASTAs are included at the original relative paths. No candidate or biological rule was
changed for this export. Copy a manifest to `analysis/run_manifest_<newtag>.csv` to start a new run.

```bash
cp analysis/run_manifest_survey_rerun20260924.csv analysis/run_manifest_my_run.csv
export TFOP_RUN_ROOT="$PWD/runs/run_my_run"
python analysis/manifest.py run_manifest_my_run.csv
python analysis/run_pipeline.py --tag my_run --check
python analysis/run_pipeline.py --tag my_run --only batch --jobs 4 --timeout 21600
python analysis/run_pipeline.py --tag my_run --from aggregate
```

PowerShell equivalents for copying and setting the run root:

```powershell
Copy-Item analysis/run_manifest_survey_rerun20260924.csv analysis/run_manifest_my_run.csv
$env:TFOP_RUN_ROOT = Join-Path $PWD 'runs/run_my_run'
python analysis/manifest.py run_manifest_my_run.csv
python analysis/run_pipeline.py --tag my_run --check
python analysis/run_pipeline.py --tag my_run --only batch --jobs 4 --timeout 21600
python analysis/run_pipeline.py --tag my_run --from aggregate
```

The driver propagates failures and refuses stale/incomplete inputs. Reissue the batch command to
resume a partial run. Use a new tag for a new code/data runtime. Keep folding off for the canonical
configuration. Do not change source during a batch. Genome acquisition, external engine availability,
and online annotation responses affect reproducibility; software tests do not reproduce a full batch.

Run `python analysis/run_discovery.py --help` for fresh discovery. The frozen canonical candidate
set includes a documented manual exclusion and cannot be assumed identical to a newly discovered set.
The supplied manifests preserve that earlier choice. General users can use `tfop genome` directly.

## Canonical outputs and downstream reports

`analysis/canonical_runs.toml` records the original paper-run provenance and retains original tags.
Its run roots are relative to `runs/` by default. Set `RIIO_RUNS_DIR` to another directory containing
the three run folders. Update the entries/provenance deliberately when promoting new canonical runs.
This v1 release is a distribution of the current software, not the exact original paper executable:
the canonical batch used `e9ba665`, while current RIO source `f80c90a` includes later presentation
and genome orchestration changes. Reference bytes and inference modules are preserved.

Once complete canonical bundles and the original discovery-cache tables are available:

```bash
python analysis/benchmarking/ecoli_kb/discovery_recall.py
python analysis/benchmarking/ecoli_kb/score_discovered.py
python analysis/benchmarking/salmonella/score_discovered.py
python reports/compute_numbers.py
python reports/build_benchmark_figures.py
python reports/build_survey_figures.py
python reports/build_reports.py
python reports/build_legends.py
```

Discovery miss diagnostics require BITACORA `*/hmmer/*.tblout` files. Set `TFOP_BITACORA_CACHE` to
their parent cache directory. The program stops if those inputs are absent rather than classifying
missing diagnostic evidence. Generated canonical bundles and the discovery cache are not code assets
and are not copied into this repository. There is no claim of completely offline paper reproduction.

The scoring helpers use the supplied benchmark tables and an exact accession-membership export
from the paper's per-family FASTA headers. Its source hashes and transformation are recorded alongside
`paper_family_members.json`. The obsolete database-building program and full external RegulonDB mirror
are unnecessary for scoring and excluded. No family assignment logic was changed.

The manuscript composite figure is optional and skips when manuscript artwork is unavailable.
Standalone prediction figures, report code and software contracts remain included.
