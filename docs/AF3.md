# AlphaFold 3 job sets

AF3 exports are a manual structure-modeling hand-off; they do not feed operator or inducer predictions.

```bash
tfop af3 runs/ecoli --format both --top-n 3 --batch-size 30
tfop af3 runs/ecoli/jobs/MyTF --format local
```

Run export writes `af3/server/af3_server_batch_NN.json`, `af3/local/<job>.json`,
`af3/af3_manifest.tsv` and export warnings. Bundle rebuild writes a sibling `<bundle>_af3` directory,
never into the existing bundle. The manifest records TF, kind, operator rank/locus, cognate chemistry,
server analogue, expected DNA form and exported paths. The 30-job chunk size is a configurable batching
convention, not a promise about current service quotas.

Each supported TF gets an apo oligomer, DNA jobs for up to the top three unique 60-bp duplexes,
a cognate-metal job when chemistry is supported, and a metal plus top-one DNA job. Operators retain
ranked-file order and are deduplicated by locus and duplex sequence. No operator means no DNA jobs;
unavailable chemistry is recorded instead of emitting an invalid job. Protein copy counts use the
family oligomer model (including supported tetramers); family/site assumptions are export metadata,
not experimentally established stoichiometry for each protein.

AF Server JSON uses the allowlisted server ion identifiers; local AF3 JSON uses verified CCD species.
Examples: Fe(II) uses server MN with an analogue flag and local FE2; nickel uses CO/NI; Cu(I) uses
server CU/local CU1 and explicitly records the oxidation-state difference. Zinc placeholders for
cadmium, lead, mercury or arsenic are marked approximations. Local AST denotes free arsenite, not a
protein-coordinated arsenic complex. Antimonite has no verified local CCD mapping, so local metal jobs
are omitted with a warning. Molybdate has local MOO but no defensible allowed server analogue.
Organic/redox/persulfide inducers produce no metal jobs. Fur structural zinc is not inserted solely
from a family label; additional metal-site claims require evidence.

Server jobs and local jobs use their distinct dialects. Local `ccdCodes` represent ligands; input files
are checked internally before writing. Internal contract validation does not replace the current AF3
installation/server upload schema check. Actual upload validation requires the target installation or
user's AF Server account. Review analogue flags before interpreting a structure as cognate chemistry.
