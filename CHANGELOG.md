# Changelog

## 1.2.2 — 2026-11-26

A reproducibility fix for the homologue-promoter top-up. Results are unchanged wherever the ColabFold
fetch succeeds. 1.2.1 was not published as a release of its own; this release carries its fixes too.

- **The ColabFold MSA expansion reads its cache first, retries transient failures and records a
  failure.** A regulator with fewer than 8 local homologue promoters gets about 60 more from a
  ColabFold alignment of its sequence. That fetch was the one step that made a rerun of the canonical
  candidates non-reproducible: across two runs of the 238, 6 fetches failed (read timeouts and refused
  connections), each leaving 4-6 promoters instead of about 65, and 3 of the 6 moved the primary
  operator.
  - A cached alignment of the same sequence is used before the server is asked.
  - Refused or dropped connections, read timeouts, HTTP 429/5xx and a busy submission status are
    retried at submission, polling and download. Other 4xx answers are not retried.
  - A needed expansion that fails is recorded as `homolog_regions.expansion_failure` in the dossier,
    in `run_record.json` and as a note in `REPORT.html`. Before, it looked like an expansion that was
    never needed.
- **Canonical companion.** The supplied manifests are still those of the 2026-11-10 run. A complete
  rerun of its 238 candidates with 1.2.1 left 235 primary operators unchanged; the other three are
  survey regulators whose ColabFold fetch failed in the 2026-11-10 run. A run with this release
  therefore differs from that run on the fields the 1.2.1 entry names and, when ColabFold answers, at
  those three primaries.

## 1.2.1 — 2026-11-23

Three defects of 1.2.0, fixed. Query protein sequences change at residue 1 for regulators with an
alternative start codon, and the operator anchor moves for + strand regulators; a run made with this
version is not byte-comparable with a 1.2.0 run on those fields.

- **Alternative start codons are translated as methionine.** In the bacterial code (NCBI translation
  table 11) GTG, TTG, CTG, ATT, ATC and ATA can initiate, and the initiator codon is read as Met. The GFF
  route translated codon by codon, giving Val for a GTG start and Leu for a TTG start: 54 of the 238
  canonical query proteins began with V or L where the NCBI protein record has M. A new
  `genome_scan.translate_cds` is used by the GFF route, the GenBank route (CDSs without `/translation`)
  and the DNA-query helpers. A CDS marked as missing its 5' end, or starting out of frame, keeps its
  literal first codon.
- **The + strand operator anchor uses the real intergenic gap.** `locate_operator_window` took the
  regulator itself as its upstream neighbour, so every + strand anchor was the 160-bp fallback. It now
  mirrors the - strand rule. Over the 238 canonical runs, 61 anchors move and `nearest_bp` changes in
  60; no autoregulation call and no operator changes.
- **The regulator is passed to the genome rescan as its own gene, not by name.** Where a gene name
  repeats in a genome, the rescan could anchor on another copy. A string regulator is now resolved by
  accession or locus tag first, and by name only when the name is unique (a repeated name raises).

## 1.2.0 — 2026-11-23

Inducer calls change in this release. Operator and regulon calls change only where the regulator was
previously resolved to an overlapping neighbouring gene.

- **Inducer evidence (changes predictions).** The derived effector route reads the regulator's own
  operon (divergent and co-transcribed), with each source labelled. Sources are ranked by kind (SSN
  clade, then published database, then local operon chemistry) before confidence. Published homologues
  at 99 % identity or more may vote for their organic ligands; lower-identity ones are reported as
  context only. Bare element symbols resolve to the explicit ion, and a clade set that spans metal and
  organic classes is split by the protein's own coordination gate. Chemical-identity uncertainty and an
  explicit sensor class are kept and shown in the report.
- **Inducer fusion fixes (changes predictions).** A negative coordination gate is a non-metal call when
  no accepted source carries a metal or redox class. An assigned non-metal SSN clade outranks the bare
  family prior. A negative gate drops the metal member of a mixed clade set. Metal ions named by operon
  chemistry lose their transport-compartment tag (`Cu(+)(in)` is reported as `Cu+`).
- **ArsR/SmtB metal gate (changes predictions).** The α5 channel requires the catalogued
  `D.H.{10}H.{2}[HE]` site instead of a bare `H..H` pair, and the `ArsR_c7` seed label is corrected to
  the persulfide clade (YgaV/HlyU/BigR).
- **Regulator choice.** Where neighbouring genes overlap, the regulator is the gene with the best
  overlap with the query, not the upstream one.
- **GUI.** A redesigned interface started with `tfop gui` (first free port from 8501, `--port`,
  `--no-browser`). Saved runs reopen in the page, and each run shows the equivalent `tfop` command for
  cmd/PowerShell or POSIX shells. Results views read the stored report fields.
- **Data integrity.**
  - Caches are keyed by sequence, not job name. Network failures are no longer cached as answers.
    Cache files are written atomically, so parallel jobs never read half a file.
  - A partial accession no longer matches another mirrored genome.
  - Copies of one protein sharing a RefSeq `WP_` id are no longer concatenated, and CDS phase is
    honoured.
  - Local multi-replicon genomes are scanned whole.
  - RefSeq GenBank files with duplicated identical proteins are accepted.
  - UniProt `sp|ACC|` headers resolve.
- **Runner and CLI.**
  - BITACORA works with paths containing spaces or shell characters, and never executes them.
  - Deep run folders work.
  - `tfop selftest` reports code errors and timeouts as failures, not as a missing engine.
  - `tfop scan` validates its families and exits non-zero on failure.
- **Speed.** The genome rescan no longer recomputes null distributions and flanking genes for every
  intergenic gap. Output is identical and a typical prediction is about 40 % faster.
- **Canonical companion.** The supplied manifests and FASTAs are now those of the 2026-11-10 run
  (`analysis/canonical_runs.toml`), the run this release reproduces; the 2026-09-24 lists are removed.
- The GUI's example accession is ZntR (`NP_417751.1`); the previous one was a histidine kinase.
- **Known defects, not fixed in this release because fixing them changes recorded values:**
  - GFF-route proteins keep V/L at GTG/TTG start codons.
  - On the + strand, the operator anchor takes the regulator itself as its upstream neighbour.
  - A regulator named by a repeated gene name can resolve to another copy.

## 1.1.0 — 2026-10-09

- Add native Linux engine installation and explicit Windows/WSL instructions.
- Share sequence, protein-list, genome-accession and local-genome execution across CLI/API/GUI.
- Preserve local genome/GFF/contig/locus context and reject incompatible input/settings on resume.
- Generate portable discovery/run HTML reports and linked supplementary XLSX/TSV tables.
- Export optional server/local AF3 job sets separately from prediction and report rebuilding.
- Retain the frozen Figure 4 recipe and pinned canonical reference data. Software/export checks
  and prediction-neutrality comparisons remain distinct validation steps.
- Fix `tfop predict` on Linux with Python 3.10-3.13 for proteins longer than 255 amino acids (also the
  README example): a raw sequence was probed as a file name, which Linux rejects with an error.
- Add `env/install_windows.ps1`, a native Windows installer for RIIO, MMseqs2, BLAST+ and NCBI Datasets
  that needs no administrator rights; `-WithWsl` adds BITACORA genome discovery through WSL.
- Rewrite the Linux and Windows installation guides for new users, with a first prediction on the
  bundled ZntR sequence, and replace the README example accession, which was not a TF.
- Give `tfop setup` installation hints that work on each OS and print the settings file it reads;
  read settings files as UTF-8 (byte-order mark ignored) so non-ASCII user folders survive on Windows.
- Serve `tfop-gui` on this computer only and without usage statistics; extra Streamlit flags pass
  through. Document every CLI option and correct the per-user settings path on Linux.
- Make `environment.yml` install this checkout instead of the 1.0.0 tag.

## 1.0.0 — 2026-10-06

First clean RIIO distribution, exported from RIO `f80c90a` and the analysis source commit recorded in
`docs/SOURCE_MANIFEST.json`. The paper's canonical prediction baseline is RIO `e9ba665`.

- Retain the supported CLI, GUI, API, active production drivers and required reference inputs.
- Package the owned BITACORA/MetalNet runners, BITACORA patch, and active folding-adapter helpers.
- Resolve the consolidated predictor locally and make canonical output roots portable.
- Include full upstream license notices, contribution evidence and provisional software citation.
- Exclude experimental code, historical runs, credentials, external engines and model weights.
- Preserve core reference database bytes and inference algorithms; release verification is documented.
