# Entry points and resume behavior

All interfaces use the same per-TF prediction core. A discovery family is a routing hint;
the core makes its own family call. No entry mode modifies operator or inducer ranking.

```bash
# Raw sequence, one-record FASTA, or protein accession (accessions require network)
tfop predict protein.faa --name MyTF --organism GCF_000005845.2
# Predict every provided record; no BITACORA/Pfam filtering
tfop batch proteins.faa --organism GCF_000005845.2 --out runs/list --jobs 2
# Production discovery from accession or local input
tfop genome GCF_000005845.2 --out runs/accession
tfop genome --fasta genome.fna --gff genome.gff --out runs/local
tfop genome --genbank genome.gbff --out runs/genbank
tfop genome --proteome proteins.faa --out runs/proteome
# Cheap census, explicitly a different discovery method
tfop scan proteins.faa --out runs/census
```

An accession list is tab-separated with an `accession` header, one identifier per row.
Multi-FASTA and TSV lists preserve record order; duplicate identifiers and sanitized filename
collisions are errors. Provide amino-acid sequences, not untranslated nucleotide input.
A supplied list predicts every record, including proteins outside the validated family subset;
capability/abstention labels remain visible in the per-TF evidence.

Nucleotide FASTA preserves a supplied/sibling GFF3 when available, otherwise calls genes with Pyrodigal.
GenBank inputs support multiple records, preserve CDS identifiers, annotated translations and contigs,
and reject origin-wrapped CDS that the linear coordinate schema cannot represent safely.
Local genome files are copied into `inputs/`; per-TF subprocesses receive that exact context and locus.
An organism label alone does not replace local nucleotide context. Proteomes need an accession or
network-resolved source genome for complete downstream operator prediction.

Resume checks compare the stored input fingerprint with the supplied sequence/context files,
accession and discovery policy. Candidate names must be unique portable path components; manifest
FASTA paths must stay inside the run. Before reusing a completed bundle with a stored run record,
the runner checks the exact protein hash and effective online/folding/MetalNet2/AF3 settings. A mismatch
becomes a visible failed resume rather than silently reusing another input or configuration. Use a new
`--out` directory when changing these inputs/settings. Legacy local/list runs without fingerprints
must be rebuilt as sidecar reports or restarted in a fresh directory; report/AF3 rebuilding never
re-predicts or rewrites existing bundles.

Production discovery is BITACORA with packaged SSN-clade profiles. `--census-fallback` allows the Pfam
census when production discovery is unavailable, and records the method on every row. GUI Quick census
explicitly selects the Pfam route. `--discover-only` writes the manifest and rebuildable reports without
predicting candidates. `--offline` prohibits network retrieval; it does not fill missing genome/homolog caches.

```python
from predictor.api import predict, predict_many, discover, run_genome

result = predict("protein.faa", name="MyTF", organism="GCF_000005845.2")
summary = predict_many("proteins.faa", "runs/list", organism="GCF_000005845.2")
summary = run_genome(out_dir="runs/local", fasta="genome.fna", gff="genome.gff")
candidates, engine, protein_source = discover("genome.fna", gff="genome.gff")
```

Runs contain `manifest.csv`, `run_input.json`, `fastas/`, `jobs/`, `logs/`, `discovery/`, reports,
supplementary tables and AF3 exports. The input fingerprint prevents resuming changed local files or
lists in the same output directory. Complete bundles resume without prediction; failed candidates can
be retried. Use a new output directory when inputs change. Legacy accession manifests remain resumable
when the requested genome matches; legacy local/list runs without fingerprints require a fresh directory.
The reporting commands support older layouts without requiring prediction.

The GUI Genome scan and Protein list views call the same runner and allow a subset selection while
keeping the complete discovery manifest/table. Failed and unselected rows remain visible.
