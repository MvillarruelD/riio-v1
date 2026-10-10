# Generic software checks

The following command group is the reviewed generic boundary for packaging and local
report-file behavior. It uses explicit test nodes and does not run release validation,
reference integrity suites, external engines, network work, or predictions:

```powershell
# Run from the predictor repository root with its environment activated.
python -m pytest `
  tests/test_packaging.py::test_every_module_imports `
  tests/test_packaging.py::test_no_developer_home_paths_in_shipped_code `
  tests/test_packaging.py::test_no_committed_credentials `
  tests/test_packaging.py::TestProjectMetadata `
  tests/test_packaging.py::TestReleaseCoherence::test_non_runtime_trees_are_outside_the_package `
  tests/test_packaging.py::TestReleaseCoherence::test_research_archive_is_never_imported `
  tests/test_packaging.py::TestReleaseCoherence::test_packages_find_includes_only_the_predictor_package `
  tests/test_packaging.py::TestReleaseCoherence::test_the_build_frontend_is_a_declared_dev_dependency `
  tests/test_bundle_files.py `
  tests/test_user_experience.py::test_report_leads_with_a_decision_summary_and_validation_warning `
  tests/test_user_experience.py::test_report_shows_the_fold_panel_only_when_a_fold_happened
```

Tests marked `released_data` are excluded by default. The marker covers all known
tests that open shipped package records, references, data manifests or generated
reference measurements. Mixed modules use exact function markers; the fully
record-backed Ligify module uses a module marker. Use `-m released_data` only when
those reads are intended.

`python -m pytest` is the broad suite excluding known marked record reads, not the narrow command above.
`python -m pytest -m released_data` is its separate shipped-data integration gate.
`tools/release_check.py` runs both pytest gates separately and then performs its full
artifact/install checks. None of those broader commands is implied by this generic
software check.

The project supports Python 3.10. Tests import the standard-library `tomllib` on
Python 3.11+ and fall back to the conditional `tomli` development dependency on
Python 3.10. A metadata assertion alone does not demonstrate that the suite ran on
3.10; record that interpreter cell as untested when it is unavailable.
