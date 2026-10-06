"""discovery -- genome in, candidate transcription factors out.

Two routes to the same list, and which one ran is always recorded:

* `annotate.genome_scan` -- our Pfam-HMM census (89 models across 22 families, at curated gathering
  thresholds). Offline,
  ~5 s per genome, no external engine beyond the bundled `pyhmmer`.
* `discovery.bitacora` -- the survey's own tool (Rondon et al.), run by us rather than consumed as a
  foreign artifact. Needs the BITACORA checkout plus BLAST+ and HMMER; absent, it abstains.

Measured on the survey's four genomes, matched by protein ACCESSION: the census recovers **149/150** of
BITACORA's published candidates (the 2026-08-19 discovery decision, recorded in the analysis repo
before commit 5815559). So the census is a validated route, not a guess -- but the two are not interchangeable in provenance, and a
table built from one must never be labelled as the other.

`vendor/` holds the collaborator's own result-processing scripts verbatim, under their MIT licence.
They are reference, not our code path: do not edit them, and compute nothing from them that belongs in
a file of ours.
"""
