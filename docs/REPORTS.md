# Discovery and global run reports

`tfop genome`, `tfop batch`, and GUI multi-protein runs write reports outside the per-TF bundles.

| File | Content |
|---|---|
| `discovery/DISCOVERY_REPORT.html` | One row per candidate grouped by discovery family (gene, product, matched SSN-clade profiles, location, length, link to its per-TF report when predicted); one row per searched profile with its raw BITACORA files, empty files labelled |
| `discovery/discovery_candidates.tsv` / `candidates.faa` | Reusable discovery exports |
| `RUN_REPORT.html` | One row per predicted TF, sortable within each family using keyboard-accessible headings; discovered-but-unpredicted candidates listed apart; saved family/class counts and available operator-confidence counts; provenance with the raw record behind a toggle and links to each bundle's `run/run_record.json` |
| `Supplementary_Table_S1.xlsx` / `.tsv` | All manifest candidates with relative report paths; XLSX hyperlinks for existing reports |
| `summary.tsv` | Headline prediction calls read back from each complete bundle |

Actual BITACORA search outputs stay under the run's `discovery/` directory. Existing historical runs may
lack those raw tables or evidence; rebuilt reports state available evidence and do not invent it.
Protein-list runs label the supplied proteins as `provided_protein`; no discovery filtering is claimed.

```bash
tfop report runs/ecoli
tfop af3 runs/ecoli --format both
tfop report runs/ecoli
```

The final report rebuild adds links to the AF3 manifest produced by the intervening export.
Commands require an existing run directory. They read manifests, summaries, dossiers and ranked operators
without prediction, and do not write inside `jobs/*`. Missing per-TF HTML remains an unlinked row.
Keep `RUN_REPORT.html`, supplementary files, `discovery/`, `af3/` and `jobs/` together when sharing;
relative links depend on that structure. Opening the report as a local file preserves navigation better
than an embedded preview with no filesystem base URL.

The GUI offers **Prepare full run download** and **Prepare full bundle download** actions that build
ZIP snapshots only when requested. Extract them before opening HTML so linked reports and images
retain their relative paths. Individual HTML downloads require the saved folder. GUI previews embed
contained images and disable file links; external references and section links remain available.
Downloads omit linked paths and hidden temporary directories, reporting omissions. Interrupted
hidden publication drafts are also excluded from stored-run candidate counts.

For canonical runs, hash/validate existing bundles before and after rebuilding, and record the source
prediction commit separately from the export commit. A successful report rebuild proves export integrity;
it does not independently verify prediction neutrality or the biological calls.
