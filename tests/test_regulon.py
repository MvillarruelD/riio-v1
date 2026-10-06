"""Regulon reconstruction: operon assembly, promoter assignment, duplicate-site merge, honest FDR.

These pin behaviour that a RegulonDB benchmark showed was wrong, and in two cases pin a *statistical*
property that is easy to silently undo.
"""
import numpy as np

from predictor.signals import regulon as RG
from predictor.signals.motif_rescan import Gene, GenomeContext, rescan_genome, _merge_adjacent
from predictor.motifs.pwm_scan import pwm_from_counts

MOTIF = "TTGACCTAGGTCAA"


def _pwm(motif=MOTIF):
    idx = {"A": 0, "C": 1, "G": 2, "T": 3}
    counts = np.zeros((4, len(motif)))
    for j, ch in enumerate(motif):
        counts[idx[ch], j] += 1
    return pwm_from_counts(counts)


def _genome(planted, genes, length=4000, seed=7):
    import random
    random.seed(seed)
    seq = list("".join(random.choice("ACGT") for _ in range(length)))
    for at in planted:
        seq[at:at + len(MOTIF)] = list(MOTIF)
    return GenomeContext("SYN.1", "".join(seq), genes)


class TestOperonAssembly:
    """An operator upstream of a polycistronic unit regulates every gene in it. The walk used a 60 bp
    hard cutoff, which truncated real operons at their first wide gap."""

    def test_expands_through_a_gap_wider_than_the_old_60bp_cutoff(self):
        genes = [Gene(1000, 1100, "+", "a"), Gene(1220, 1300, "+", "b")]   # 120 bp gap
        got = [g.name for g in RG.operon_of(genes, genes[0])]
        assert got == ["a", "b"], f"120 bp gap should stay in-operon at the {RG.DEFAULT_OPERON_MAX_GAP} bp default"

    def test_still_stops_at_a_genuinely_large_gap(self):
        genes = [Gene(1000, 1100, "+", "a"), Gene(2000, 2100, "+", "b")]   # 900 bp
        assert [g.name for g in RG.operon_of(genes, genes[0])] == ["a"]

    def test_stops_at_a_strand_change(self):
        genes = [Gene(1000, 1100, "+", "a"), Gene(1120, 1200, "-", "b")]
        assert [g.name for g in RG.operon_of(genes, genes[0])] == ["a"]

    def test_minus_strand_walks_toward_lower_coordinates(self):
        genes = [Gene(1000, 1100, "-", "up"), Gene(1150, 1250, "-", "head")]
        got = [g.name for g in RG.operon_of(genes, genes[1])]
        assert got == ["head", "up"]


class TestPromoterAssignment:
    """Considering only the single nearest gene per side -- and discarding it when its strand pointed the
    wrong way -- let one convergent neighbour hide the real target."""

    def test_looks_through_a_wrong_strand_neighbour(self):
        # site at 1000; nearest gene to the right is '-' (wrong way), the real target sits just past it
        genes = [Gene(1010, 1080, "-", "decoy"), Gene(1100, 1200, "+", "target")]
        ctx = GenomeContext("S", "A" * 3000, genes)
        assert "target" in [g.name for g in RG.regulated_first_genes(ctx, 1000)]

    def test_respects_the_distance_bound(self):
        genes = [Gene(2500, 2600, "+", "far")]
        ctx = GenomeContext("S", "A" * 4000, genes)
        assert RG.regulated_first_genes(ctx, 1000, max_dist=400) == []

    def test_divergent_promoter_drives_both_sides(self):
        genes = [Gene(900, 980, "-", "left"), Gene(1020, 1100, "+", "right")]
        ctx = GenomeContext("S", "A" * 3000, genes)
        names = {g.name for g in RG.regulated_first_genes(ctx, 1000)}
        assert names == {"left", "right"}


class TestDuplicateSiteMerge:
    """One real palindromic operator surfaces as several hits (both strands + shifted windows). Unmerged
    they trebled the site count and distorted the FDR denominator -- visible in the benchmark as true
    sites appearing in adjacent-rank pairs."""

    def test_both_strand_readings_collapse_to_one_site(self):
        ctx = _genome([470], [Gene(600, 900, "+", "g")])
        hits = rescan_genome(ctx, pwm=_pwm(), scope="intergenic", pvalue_thresh=1e-4)
        near = [h for h in hits if abs(h.dyad_center - 477) <= 10]
        assert len(near) == 1, f"the planted palindrome should be ONE site, got {len(near)}"

    def test_merge_keeps_the_best_scoring_representative(self):
        from predictor.signals.motif_rescan import GenomeHit
        mk = lambda start, score: GenomeHit("A", start, start + 14, "+", start + 7,
                                            score, 1e-5, 1e-5, 3, 1.0, 5.0, "")
        kept = _merge_adjacent([mk(100, 5.0), mk(103, 9.0), mk(500, 4.0)], 14)
        assert sorted(h.score for h in kept) == [4.0, 9.0]


class TestHonestFDR:
    """The q-value must reflect the multiple-testing burden of the WHOLE genome scan.

    Before, BH ran inside each intergenic interval, so a lone hit got q == p and any downstream q-gate was
    structurally inert. A later bug ran BH over only the hits that survived the p-filter, which assigned
    every survivor q == p and passed them all. The denominator has to be the windows actually tested.
    """

    def test_qvalue_reflects_the_genome_wide_burden(self):
        ctx = _genome([470, 1460, 2645],
                      [Gene(600, 900, "+", "a"), Gene(1600, 1800, "+", "b"), Gene(2300, 2600, "-", "c")])
        hits = rescan_genome(ctx, pwm=_pwm(), scope="intergenic", pvalue_thresh=1e-4)
        assert hits
        for h in hits:
            assert h.qvalue >= h.pvalue, "a q-value can never be smaller than its p-value"
            # with ~thousands of windows tested, an exact match at p~4e-9 must be inflated well above p
            assert h.qvalue > h.pvalue * 10, (
                f"q={h.qvalue:.2e} barely exceeds p={h.pvalue:.2e} -- the FDR denominator looks like it "
                f"counts only the surviving hits again, not the windows tested")

    def test_missing_qvalue_fails_closed(self):
        """A hit whose FDR could not be computed must be rejected, not treated as perfectly significant."""
        from predictor.signals.motif_rescan import GenomeHit
        nan_hit = GenomeHit("A", 1000, 1014, "+", 1007, 9.0, 1e-9, float("nan"), 1, 1.0, 9.0, "")
        genes = [Gene(1100, 1300, "+", "target")]
        ctx = GenomeContext("A", "A" * 3000, genes)
        rg = RG.reconstruct_regulon(ctx, hits=[nan_hit], qvalue_thresh=0.05)
        assert rg.operons == [], "a NaN q-value must fail closed"


class TestRankSelection:
    """Selection is by rank, because a q-gate provably cannot discriminate at these motif strengths."""

    def test_max_operons_caps_the_output(self):
        ctx = _genome([470, 1460, 2645],
                      [Gene(600, 900, "+", "a"), Gene(1600, 1800, "+", "b"), Gene(2300, 2600, "-", "c")])
        rg = RG.reconstruct_regulon(ctx, pwm=_pwm(), max_operons=1)
        assert len(rg.operons) == 1

    def test_uncapped_when_none(self):
        ctx = _genome([470, 1460, 2645],
                      [Gene(600, 900, "+", "a"), Gene(1600, 1800, "+", "b"), Gene(2300, 2600, "-", "c")])
        rg = RG.reconstruct_regulon(ctx, pwm=_pwm(), max_operons=None)
        assert len(rg.operons) >= 3

    def test_precomputed_hits_match_an_internal_scan(self):
        """`hits=` must be equivalent to letting reconstruct_regulon scan itself -- it is what makes a
        parameter sweep cheap, so a divergence would silently invalidate every swept result."""
        ctx = _genome([470, 1460], [Gene(600, 900, "+", "a"), Gene(1600, 1800, "+", "b")])
        pwm = _pwm()
        a = RG.reconstruct_regulon(ctx, pwm=pwm)
        b = RG.reconstruct_regulon(ctx, hits=rescan_genome(ctx, pwm=pwm, scope="intergenic",
                                                           pvalue_thresh=1e-4))
        assert [o.first_gene for o in a.operons] == [o.first_gene for o in b.operons]


class TestGeneModel:
    """Non-coding genes must be in the model: an operon walk blind to a tRNA measures the gap across it."""

    def test_ncrna_features_are_parsed(self, tmp_path):
        from predictor.annotate import context as CTX
        gff = tmp_path / "t.gff"
        gff.write_text(
            "##gff-version 3\n"
            "c1\t.\tregion\t1\t3000\t.\t+\t.\tID=r\n"
            "c1\t.\tCDS\t100\t400\t.\t+\t0\tID=c1;gene=cdsA\n"
            "c1\t.\ttRNA\t450\t530\t.\t+\t.\tID=t1;gene=trnW\n"
            "c1\t.\tCDS\t600\t900\t.\t+\t0\tID=c2;gene=cdsB\n", encoding="utf-8")
        fna = tmp_path / "t.fna"
        fna.write_text(">c1\n" + "ACGT" * 750 + "\n", encoding="utf-8")
        ctx = CTX.from_gff(str(gff), str(fna))
        assert "trnW" in {g.name for g in ctx.genes}, "tRNA genes must be in the gene model"


class TestDeepHomologRoute:
    """The SSN-cluster homolog route must be reachable for EVERY family that has an SSN.

    Two hardcoded two-family assumptions, written when the SSN held only MerR and ArsR/SmtB, survived the
    2026 drop that took it to twelve and silently disabled the deep-homolog route for the other ten:
    `pipeline._ssn_families` and `msa_homologs._cluster_by_id`, the latter resolving any non-MerR cluster
    id to ArsR/SmtB so `Fur_ecFur`, `Rrf2_ecIscR`, `NikR_c1` ... were looked up in the wrong family list
    and returned None. Those families then fell back to a ~12-sequence BLAST set. Enabling the route
    raised true-operator recovery 12 -> 21 sites and top-decile precision 0.114 -> 0.136 on the RegulonDB
    panel, so a silent regression here would cost real accuracy.
    """

    def test_cluster_lookup_resolves_every_family(self):
        from predictor.annotate import ssn_clusters as ssn
        from predictor.annotate.msa_homologs import _cluster_by_id
        missing = []
        for fam in ssn.FAMILIES:
            clusters = ssn.load_clusters(fam)
            if not clusters:
                continue
            cid = clusters[0].cluster_id
            if _cluster_by_id(cid) is None:
                missing.append(f"{fam} -> {cid}")
        assert not missing, f"cluster ids unresolvable (deep-homolog route disabled): {missing}"

    def test_pipeline_gate_covers_every_ssn_family(self):
        """The pipeline's gate must be driven off ssn.FAMILIES, never a hardcoded subset of families."""
        import inspect
        import re
        from predictor import pipeline
        from predictor.annotate import ssn_clusters as ssn
        src = inspect.getsource(pipeline)
        gate = [ln for ln in src.splitlines() if "in set(" in ln and "FAMILIES" in ln]
        assert gate, "the SSN-cluster gate should test membership of ssn.FAMILIES"
        # no family name may appear as a literal anywhere in the gate's function -- that is how the
        # two-family freeze got in last time
        named = {f for f in ssn.FAMILIES if re.search(rf'["\']{re.escape(f)}["\']', src)}
        assert not named, f"pipeline hardcodes family names: {sorted(named)}"
        assert len(ssn.FAMILIES) == 12

    def test_no_module_hardcodes_a_two_family_split(self):
        """No module may branch on one family name with another as the else -- the exact shape of the
        freeze that disabled ten families in msa_homologs and cluster_fold."""
        import ast
        from pathlib import Path
        import predictor
        from predictor.annotate import ssn_clusters as ssn
        root = Path(predictor.__file__).parent
        fams = set(ssn.FAMILIES)
        offenders = []
        for p in root.rglob("*.py"):
            tree = ast.parse(p.read_text(encoding="utf-8", errors="replace"))
            for node in ast.walk(tree):
                # `"<family>" if <cond> else "<family>"` -- real code only; the AST never sees the
                # docstrings and comments in which the old bug is deliberately recorded.
                if not isinstance(node, ast.IfExp):
                    continue
                ends = [node.body, node.orelse]
                if all(isinstance(e, ast.Constant) and e.value in fams for e in ends):
                    offenders.append(f"{p.relative_to(root)}:{node.lineno}")
        assert not offenders, f"two-family hardcoded split reintroduced: {offenders}"
