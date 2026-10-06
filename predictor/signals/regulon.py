"""
regulon.py -- reconstruct a TF's REGULON (its full target-gene set) from a discovered operator motif.

The operator scanner / motif_finder gives the *autoregulatory* operator (the site 5' of the TF). This
module takes that operator as a PWM and extends it to the whole regulon -- exactly RegPredict's
"regulon reconstruction for a known motif" workflow, and the structure-free output side of the biosensor
loop (ligand -> TF -> *regulon*):

  1. scan the genome's intergenic DNA with the operator PWM   (motif_rescan -> exact-p hits, tiered),
  2. assign each significant site to the operon(s) it drives   (genes transcribed AWAY from the site),
  3. extend each first gene into its operon                    (co-directional, short intergenic gaps),
  4. emit the predicted regulon                                (operons + member genes + the site).

It consumes a `GenomeContext` -- which now comes OFFLINE from the dereplicated mirror
(`context.genome_context_from_mirror`, .fna + .gff) so a genome-wide scan needs no EFetch. The conserved
homolog-neighborhood signal (genes recurrently adjacent to TF homologs) is an independent reinforcement
hook (`conserved_targets=`); validate the output against RegPrecise / RegulonDB where the organism is
covered.

Run `python regulon.py` for a self-test (regulator + two target operons; both recovered).
"""
from __future__ import annotations

from dataclasses import dataclass, field

from predictor.signals.motif_rescan import Gene, GenomeContext, rescan_genome
from predictor.motifs.pwm_scan import pwm_from_counts


@dataclass
class RegulonOperon:
    first_gene: str             # the operon's leading gene (the one the promoter directly drives)
    strand: str
    genes: list                 # member gene names, in transcription order
    protein_ids: list           # CDS protein accessions (parallel to genes), for effector inference
    site_start: int             # the regulating operator site (genome 0-based half-open)
    site_end: int
    site_strand: str
    score: float
    pvalue: float
    qvalue: float
    tier: int
    conserved_target: bool = False


@dataclass
class Regulon:
    accession: str
    operons: list = field(default_factory=list)      # list[RegulonOperon], best-evidence first
    n_sites: int = 0
    #: How many operons passed the site filter BEFORE `max_operons` truncated the list. Recorded
    #: because the cap is reached by nearly every candidate, so `len(operons)` alone reads as a
    #: measured regulon size when it is the cap being reported. 0 for an untruncated empty result.
    n_operons_found: int = 0

    @property
    def truncated(self) -> bool:
        return self.n_operons_found > len(self.operons)

    @property
    def target_genes(self):
        seen, out = set(), []
        for op in self.operons:
            for g in op.genes:
                if g and g not in seen:
                    seen.add(g)
                    out.append(g)
        return out


# --------------------------------------------------------------------------- operon assembly
#: Largest intergenic gap still treated as "same operon". Set from the genome-organisation literature,
#: NOT fitted to the benchmark: adjacent gene pairs inside an operon are strongly enriched below ~50 bp
#: and depleted beyond a few hundred (Salgado et al. 2000, PNAS 97:6652; Moreno-Hagelsieb & Collado-Vides
#: 2002, Bioinformatics 18:S329), so operon predictors conventionally cut around 100-150 bp. The previous
#: 60 bp was tight enough to truncate real operons at their first wide gap. The hard stop is kept on
#: purpose -- a large gap genuinely marks a new transcription unit; what was wrong was the threshold, and
#: (see annotate.context) an ncRNA-blind gene model that measured gaps straight across tRNA/rRNA genes.
DEFAULT_OPERON_MAX_GAP = 150

#: How many operons to keep, ranked by operator score. This is a REPORTING bound -- enough to cover the
#: regulons this tool is used on without dumping thousands of rows into a dossier -- and deliberately NOT
#: a tuned parameter. Selection is by rank rather than significance because a q-gate provably cannot
#: discriminate here: with an honest genome-wide FDR (`motif_rescan._genome_wide_qvalues`) a hit needs
#: p <= ~5e-9 to reach q <= 0.05 while these PWMs deliver p ~ 1e-4, so the gate admits either everything
#: or nothing.
#:
#: What the RegulonDB sweep actually showed (project_archive/benchmarking_regulondb/): the score ORDER does carry signal
#: -- sites overlapping a true target are ~4x enriched in the top decile -- but only ~2% of genome-wide
#: sites lie near a true target at all, and 12 of 17 panel TFs get no true site whatsoever. So no value of
#: this cap rescues the operator route; tightening it only trades recall for precision along a curve.
#: A cap of 3 posts the best MACRO F1 (0.108) purely as an artifact of several 3-gene regulans in the
#: panel: pooled across genes it collapses to micro-recall 0.021 and gives Fur 0 of its 282 genes. It is
#: set at 25 for that reason -- large enough not to truncate a real regulon, and not chosen to flatter a
#: benchmark. Operator-route precision is bounded by motif specificity, not by post-processing.
DEFAULT_MAX_OPERONS = 25


def operon_of(genes_sorted, first_gene, *, max_gap: int = DEFAULT_OPERON_MAX_GAP) -> list:
    """Extend a first gene into its operon: consecutive co-directional genes with intergenic gap
    <= max_gap, walked in the direction of transcription. Returns Genes in transcription order."""
    idx = genes_sorted.index(first_gene)
    out = [first_gene]
    if first_gene.strand == "+":
        i = idx
        while i + 1 < len(genes_sorted):
            nxt = genes_sorted[i + 1]
            if nxt.strand == "+" and nxt.start - genes_sorted[i].end <= max_gap:
                out.append(nxt)
                i += 1
            else:
                break
    else:                                            # '-' operon transcribes toward lower coordinates
        i = idx
        while i - 1 >= 0:
            prv = genes_sorted[i - 1]
            if prv.strand == "-" and genes_sorted[i].start - prv.end <= max_gap:
                out.append(prv)
                i -= 1
            else:
                break
    return out


def regulated_first_genes(ctx: GenomeContext, center: int, *, max_dist: int = 400,
                          look_through: int = 2) -> list:
    """The gene(s) a site at `center` directly drives: the nearest genes transcribed AWAY from the site
    (a divergent promoter drives one on each strand).

    `look_through` is how many genes deep to search on each side before giving up. Considering only the
    single nearest gene -- and discarding it outright when its strand pointed the wrong way -- meant one
    convergent neighbour hid the real target completely. A small ORF or an antisense gene sitting between
    the site and the operon it drives is common enough that the first gene is often not the right one.
    Everything remains bounded by `max_dist`, so looking deeper cannot reach past the promoter window."""
    genes = ctx.genes
    out = []
    right = sorted((g for g in genes if g.start >= center), key=lambda g: g.start)[:look_through]
    for g in right:
        if g.start - center > max_dist:
            break
        if g.strand == "+":
            out.append(g)
            break                                    # the first co-oriented gene is the operon head
    left = sorted((g for g in genes if g.end <= center), key=lambda g: -g.end)[:look_through]
    for g in left:
        if center - g.end > max_dist:
            break
        if g.strand == "-":
            out.append(g)
            break
    return out


# --------------------------------------------------------------------------- public API
def reconstruct_regulon(ctx: GenomeContext, *, pwm=None, counts=None, scope: str = "intergenic",
                        qvalue_thresh: float = 1.0, pvalue_thresh: float = 1e-4,
                        max_promoter_dist: int = 400,
                        operon_max_gap: int = DEFAULT_OPERON_MAX_GAP,
                        max_operons: int | None = DEFAULT_MAX_OPERONS,
                        conserved_targets=None, regulator=None, hits=None) -> Regulon:
    """Operator PWM -> predicted regulon. Significant intergenic sites are mapped to the operons they
    drive and de-duplicated by first gene; the best-evidence site per operon is kept.

    `conserved_targets`: optional set of gene names recurrently found next to TF homologs (the
    Snowprint comparative signal). This evidence is annotated separately; it never modifies a p/q-value.
    Combining heterogeneous evidence by multiplying an FDR value is not statistically valid."""
    if hits is None:
        if pwm is None and counts is None:
            raise ValueError("provide pwm=, counts= or hits=")
        if pwm is None:
            pwm = pwm_from_counts(counts)
        # `hits=` lets a caller supply an already-computed genome scan. The scan dominates the runtime,
        # while everything after it (promoter assignment, operon assembly, ranking) is cheap -- so reusing
        # one scan across many parameter settings turns a parameter sweep from hours into seconds.
        hits = rescan_genome(ctx, pwm=pwm, regulator=regulator, scope=scope,
                             pvalue_thresh=pvalue_thresh)
    genes_sorted = sorted(ctx.genes, key=lambda g: g.start)
    conserved = set(conserved_targets or [])

    best = {}                                        # first-gene key -> RegulonOperon (best evidence)
    n_sites = 0
    for h in hits:
        firsts = regulated_first_genes(ctx, h.dyad_center, max_dist=max_promoter_dist)
        if not firsts:
            continue
        n_sites += 1
        for fg in firsts:
            # a missing/NaN q-value must fail CLOSED. It used to fall back to 0.0, i.e. "perfectly
            # significant", so any hit whose FDR could not be computed sailed through the gate.
            q = h.qvalue if h.qvalue == h.qvalue else 1.0
            if q > qvalue_thresh:
                continue
            members = operon_of(genes_sorted, fg, max_gap=operon_max_gap)
            op = RegulonOperon(
                first_gene=fg.name, strand=fg.strand,
                genes=[g.name for g in members], protein_ids=[g.protein_id for g in members],
                site_start=h.start, site_end=h.end, site_strand=h.strand,
                score=h.score, pvalue=h.pvalue, qvalue=q, tier=h.tier,
                conserved_target=fg.name in conserved)
            key = (fg.start, fg.strand)
            if key not in best or op.score > best[key].score:
                best[key] = op

    operons = sorted(best.values(), key=lambda o: (o.qvalue, not o.conserved_target, -o.score))
    n_found = len(operons)
    if max_operons is not None:
        operons = operons[:max_operons]
    return Regulon(ctx.accession, operons, n_sites, n_operons_found=n_found)


# --------------------------------------------------------------------------- self-test
def _demo() -> None:
    import random
    import numpy as np

    motif = "TTGACCTAGGTCAA"                          # the operator (matches motif_rescan demo)
    W = len(motif)
    idx = {"A": 0, "C": 1, "G": 2, "T": 3}
    counts = np.zeros((4, W))
    for j, ch in enumerate(motif):
        counts[idx[ch], j] += 1

    random.seed(7)
    L = 4000
    seq = list("".join(random.choice("ACGT") for _ in range(L)))

    def plant(at):
        seq[at:at + W] = list(motif)

    # layout: autoregulated TF (tfR) at 500-800(+) with operator at 470; a 2-gene + operon (targA/targB)
    # at 1500/1620 driven by an operator at 1460; a single - gene (targC) at 2600-2400 driven at 2640.
    plant(470)       # drives tfR (+)
    plant(1460)      # drives targA (+) -> operon targA,targB
    plant(2645)      # drives targC (-)
    seq = "".join(seq)

    genes = [
        Gene(500, 800, "+", "tfR", protein_id="WP_tfR"),
        Gene(1500, 1600, "+", "targA", protein_id="WP_A", product="dehydrogenase"),
        Gene(1620, 1750, "+", "targB", protein_id="WP_B", product="transporter"),
        Gene(2350, 2620, "-", "targC", protein_id="WP_C", product="hydrolase"),
        Gene(3300, 3600, "+", "unrelated", protein_id="WP_U"),    # >600 bp from any site -> excluded
    ]
    ctx = GenomeContext("SYN_REG.1", seq, genes)

    # The 3 planted operators are exact matches; a chance partial match in 4 kb of random sequence is a
    # real scanner FP, and the two are separable by q-value. NOTE the absolute q scale moved when the FDR
    # became genome-wide (motif_rescan._genome_wide_qvalues): the planted sites now sit at ~7e-6 rather
    # than ~1e-9, because the correction finally counts every window tested instead of only the handful
    # inside one intergenic interval. The SEPARATION is what matters and it is unchanged -- planted ~7e-6
    # vs decoy ~1.4e-3, a ~200x gap -- so the cutoff below sits between them rather than at 1e-7.
    reg = reconstruct_regulon(ctx, counts=counts, qvalue_thresh=1e-4, regulator="tfR")
    print(f"regulon of {reg.accession}: {len(reg.operons)} operon(s), {reg.n_sites} site(s)")
    for op in reg.operons:
        print(f"  {op.first_gene}({op.strand}) operon={op.genes} "
              f"site={op.site_start}-{op.site_end} q={op.qvalue:.1e} tier={op.tier}")
    print(f"target genes: {reg.target_genes}")

    first_genes = {op.first_gene for op in reg.operons}
    assert {"tfR", "targA", "targC"} == first_genes, f"unexpected operon set: {first_genes}"
    targA_op = next(op for op in reg.operons if op.first_gene == "targA")
    assert targA_op.genes == ["targA", "targB"], f"targA operon not assembled: {targA_op.genes}"
    assert targA_op.protein_ids == ["WP_A", "WP_B"], "protein_ids not carried to regulon"
    assert "unrelated" not in reg.target_genes, "weak chance site should fall below the q cutoff"

    # at a lenient cutoff the chance site reappears but ranks strictly below all real operators
    loose = reconstruct_regulon(ctx, counts=counts, qvalue_thresh=0.1, regulator="tfR")
    real_q = max(op.qvalue for op in loose.operons if op.first_gene in {"tfR", "targA", "targC"})
    decoy = next((op for op in loose.operons if op.first_gene == "unrelated"), None)
    assert decoy is None or decoy.qvalue > real_q, "decoy should rank below all real operators"
    print("OK: regulon reconstructed (autoreg TF + 2-gene operon + divergent - gene); "
          "chance FP excluded by q-cutoff and out-ranked by all real sites.")


if __name__ == "__main__":
    _demo()
