# Operator-motif discovery: how it works, what each parameter does, and why the defaults are what they are

This note documents the operator-discovery layer of the pipeline and the controlled study that fixed its
defaults. It answers three questions a user or reviewer will ask: *how* is the operator motif built, *what*
happens if you change how it is built, and *why* is the shipped behaviour the one that ships. The study was
run on the 57 ArsR/SmtB, MerR and Fur metalloregulators of the Rondón/Antelo/Giedroc/Capdevila SSN across
four genomes (two *Mycobacterium*, two *Vibrio*); the full analysis, figures and the AlphaFold-3 hand-off
live in the companion research workspace.

---

## 1. How operator discovery works

For each TF the pipeline builds one position weight matrix (PWM) and scans the TF's own genome with it.

1. **Candidate operator.** `signals/motif_finder.predict` proposes candidate operators from the TF's own
   promoter, family-aware (half-site + spacer geometry for MerR, short inverted repeat for ArsR/SmtB, the
   Fur box for Fur).
2. **Phylogenetic footprint.** The candidate is located in the promoters of the TF's **SSN-cluster
   homologues** (their own genomes, not the study genomes) at a per-region Bonferroni threshold; the sites
   that are recovered are aligned and counted into the matrix. A functional operator is conserved across
   orthologous promoters while the surrounding sequence drifts — that conservation is the signal.
3. **Symmetrisation.** Because these regulators bind as dimers, the matrix is dyad-symmetrised.
4. **Genome-wide rescan.** `motifs/pwm_scan` scans all intergenic DNA on both strands. The score threshold
   comes from an **exact integer-lattice null** (the FIMO/MOODS method), which yields a p-value per site and
   Benjamini–Hochberg q-values across the scan. A locality prior re-ranks hits by tier — tier 1 flanks the
   regulator (its own promoter), tier 2 is divergent, tier 3 is everywhere else — but never gates them.
5. **Regulon.** The top hits are assigned to genes; the regulon is an *output* of the scan, not an input.

The autoregulatory operator (the site in the TF's own promoter) and the regulon then fall out as findings:
is there a strong operator at the regulator's own gene, and which other genes does the same motif hit?

---

## 2. What each parameter does (the controlled study)

Every lever below was swept independently and scored against the controls of §3. The consistent result is
that levers which improve the matrix *as a motif* — a deeper footprint, a sharper or more symmetric logo,
higher information content — do **not** improve the *predictions*, because motif quality is not what limits
them.

| Lever | What it changes | Effect on the motif | Effect on prediction |
|---|---|---|---|
| **Homologue sampling** — diverse (`representative`) → closest (`top_similar`) | draws the matrix from the query's nearest cluster members instead of a taxonomically diverse spread | deeper footprint (median 4→6 promoters), higher enrichment (1.4→1.9), **more cluster-specific** (same-cluster classifier AUC 0.55→0.67, bootstrap P=0.97) | **none** — corroboration and false-positive load unchanged |
| **Number of homologues** — fixed 10–40, or adaptive "collect until X carry a site" | how many promoters enter the matrix | information content climbs with the count | **none** — the median regulator has only ~5 homologues that carry a detectable operator, so an adaptive rule cannot exceed availability |
| **Palindromicity** λ (0 → 1) | blends the matrix with its reverse complement, from none to full dyad symmetry | sharper dyad symmetry, fewer genome-wide hits | false-positive *rate* flat; λ = 0 best for enrichment |
| **≥ 2-genera footprint filter** | keeps a footprint only if its sites span ≥ 2 genera | negligible | none |
| **EM refinement** (conservation-guarded) | one expectation-maximisation pass over the genome hits, accepted only if it still recovers the homologues | information content 0.25→0.31 bits/col | corroboration 27%→27% — unchanged |
| **Literature-operator seed** on/off | seeds the search with curated family operators | byte-identical output | none — the seed is inert |
| **Effector-guided boost** — *remove* | stops up-weighting the promoters of effector-inferred metal genes during the scan | — | **the one consequential change** — see below |

### The effector-guided boost is now OFF by default

`pipeline.run_novel` used to feed the promoters of genes inferred *from the effector call* back into the
scan, both as a locality boost and as the substrate for a second motif (`_target_promoter_discovery`). That
lets the effector call raise the score of operators at exactly the genes a metal-gene-based validation then
looks for, so an operator hit at a metal gene is no longer independent evidence — it can be a consequence of
having up-weighted that promoter. In the study this inflated regulon precision (0.23 → 0.65) **and** the
permutation null by the same mechanism, so enrichment over chance *fell*. **Both routes have been removed
from the pipeline** -- there is no flag to re-enable them. Operator discovery depends only on the
regulator's own promoter and its SSN-cluster homologues, so the operator-conditioned test is genuinely
independent on every run. Candidate handler genes are still reported; they simply never enter ranking.

### Bottom line

The shipped motif is at or near the ceiling of what the SSN alone supports. The binding constraint is not
the matrix mathematics but the availability of conserved sites — only ~40 % of a cluster's homologues retain
a recognisable operator, and the median regulator recovers ~5. Real improvement therefore comes from
**closer homologues or additional genomes**, not from a cleverer motif model. `top_similar` sampling does
produce measurably better-quality, more cluster-specific matrices, so it is worth exposing for tasks that
consume motif *distinctiveness* (telling clusters apart) — but it must not be sold as better *prediction*.

---

## 3. The controls (what keeps the comparison honest)

- **The no-boost baseline.** All comparisons run with the effector-guided boost off, so a hit at a
  cognate-metal gene is independent of the metal call.
- **Enrichment over the pairing null, not the raw rate.** A motif that simply hits more genes corroborates
  more often for a trivial reason; the scale-free metric divides observed corroboration by corroboration
  under a within-genome permutation of the metal↔operator pairing, rewarding specificity over promiscuity.
- **The paired per-regulator test + bootstrap.** Set-level curves can cross by chance on ~25 regulators.
  Each regulator contributes one number per motif — the rank of the first cognate-metal gene among its
  score-ranked hits, i.e. the regulon size at which it first corroborates — and motifs are compared as
  matched pairs (signed-rank + sign test + a bootstrap over regulators). No rebuilt matrix beat the shipped
  one (best config P = 0.73 for "genuinely better", against the 0.95 required).
- **The cluster-discrimination bootstrap.** Whether rebuilt motifs separate clusters better is decided by
  the area under a same-cluster classifier over all regulator pairs, bootstrapped over regulators.
- **The synthetic-operator AF3 controls.** A PWM cannot test whether its own consensus is the best binder
  (the consensus wins by construction). Degenerate-flipped, core-disrupted and scrambled variants of the
  consensus are folded alongside it, each in the same natural flank; the shipped-vs-rebuilt consensus fold
  is the decisive test. AlphaFold-3 confirms the operators as real binding sites (a scrambled operator folds
  worst) but does **not** resolve single-base changes, and the rebuilt consensus is not a better binder.
- **ipTM calibration.** Where the synthetic jobs lacked confidence summaries, ipTM was recomputed from the
  predicted-aligned-error matrix and calibrated against the complexes that carry both (Pearson r = 0.98,
  RMSE = 0.03), so the recomputed values are trustworthy.

---

## 4. Structural validation (AlphaFold-3 hand-off)

Structure prediction shares no machinery with the sequence pipeline, so it is independent corroboration:

- **Operator complexes.** Across 27 folded complexes the regulators come to rest on the predicted operators
  (contacted base pairs fall inside the operator span ~1.5× more than chance), and the contact profiles are
  bimodal — a dimer engaging the two half-sites of an inverted repeat.
- **Holo dimers (metal sites).** The Fur superfamily seats all four requested ions in fully coordinated
  sites, reproducing the two-site chemistry — a Cys4 structural zinc plus a His/Asp/Glu regulatory site —
  without being told the sites exist; ArsR/SmtB seats thiolate-rich sites; MerR, the structurally weakest
  family, barely resolves one.
- **The MerR DBD-motif correction.** The published MerR recognition motif's fourth position is a typo
  (tyrosine); the sequences carry phenylalanine/leucine/methionine there, and the corrected motif sits in
  the top decile of DNA-contacting segments.

MerR is the weakest family on six independent measures (lowest corroboration; fails its own permutation
null; no positional enrichment; lowest within-cluster motif similarity; contacts not concentrated on the
operator; no complex reaching the display-quality interface threshold), two of which come from structure and
so are not restatements of the sequence result. Treat MerR operator calls with corresponding caution.

---

## 5. Practical guidance

- Operator hits are independent evidence you can use to corroborate the effector call -- the pipeline has
  no effector-guided mode to turn on, so this holds unconditionally.
- Quote the corroboration rate *with its false-discovery threshold*: it runs from ~16 % at q ≤ 1e-6 to a
  plateau of ~59 % once every returned site qualifies. The threshold-free statement is the enrichment over
  the permutation null, which is what the biological claim rests on.
- To improve operator recovery, add genomes or draw closer homologues into the cluster; do not expect a
  different matrix recipe to help.
- The per-regulator operator sequences (autoregulatory site, top regulon-target site, consensus, degenerate
  IUPAC) are the natural probes for electrophoretic mobility shift assays, which can resolve the single-base
  operator differences that co-folding cannot.
