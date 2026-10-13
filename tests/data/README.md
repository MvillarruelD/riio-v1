# Linux engine fixture

Unmodified nucleotide fragment from the locally cached RefSeq assembly
GCF_000005845.2 (Escherichia coli K-12 MG1655), contig NC_000913.3,
1-based coordinates 3437205–3440630, containing zntR and neighbouring CDS.
The GFF preserves wholly contained RefSeq CDS with coordinates lifted to the fragment;
no synthetic or engineered sequence is introduced. mini_tf.faa is the source zntR translation.

Source genomic FASTA SHA256: `53bb6a51b6e92139ced1e38f74b7938781027c52200922ff03718c2237d23bb4`.
Fragment FASTA SHA256: `1fc7a4a9de35281a5a215a41efe4b7a315f9ca47c769e3cdeee31d76843ca51e`.

This fixture certifies real BITACORA discovery against bundled seed profiles.
It does not provide an offline homolog-promoter database for complete operator prediction.
It also excludes distant zntR regulatory targets. `tools/linux_e2e.py` verifies this fragment's
discovery first, then uses the unmodified complete assembly and online homolog-promoter context
for complete predictions. It selects the same zntR protein by exact sequence in each genome route.
