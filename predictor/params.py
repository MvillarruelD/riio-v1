"""params.py -- one table of every run parameter, read live from the code that uses it.

Nothing here holds a value. Each row names where a parameter lives (a module constant, a dataclass
default, a function-signature default or the production seed policy) and reads it at call time, so the
table cannot drift from the code. `tfop params` prints it; `tfop params --markdown` prints the Methods
table.

Add a row whenever a new threshold enters the prediction path.
"""
from __future__ import annotations

import dataclasses
import importlib
import inspect
from dataclasses import dataclass


@dataclass(frozen=True)
class Param:
    stage: str
    name: str
    where: str          # "module:SYMBOL", "module:Class.field", "module:func(arg)" or "module:POLICY.attr"
    meaning: str

    def value(self):
        mod_name, ref = self.where.split(":", 1)
        mod = importlib.import_module(f"predictor.{mod_name}")
        if "(" in ref:                                            # function-signature default
            func, arg = ref.rstrip(")").split("(")
            return inspect.signature(getattr(mod, func)).parameters[arg].default
        head, _, attr = ref.partition(".")
        obj = getattr(mod, head)
        if not attr:
            return obj
        if dataclasses.is_dataclass(obj) and isinstance(obj, type):      # dataclass field default
            return next(f.default for f in dataclasses.fields(obj) if f.name == attr)
        return getattr(obj, attr)                                         # attribute of an instance


PARAMS: tuple[Param, ...] = (
    # --- discovery -----------------------------------------------------------------------------------
    Param("discovery", "BITACORA E-value", "discovery.bitacora:SEARCH_EVALUE",
          "clade-profile search cut (hmmsearch + blastp)"),
    Param("discovery", "HSP coverage, query", "discovery.bitacora:COVERAGE_QCOV", "minimum query coverage"),
    Param("discovery", "HSP coverage, subject", "discovery.bitacora:COVERAGE_SCOV", "minimum subject coverage"),
    # --- family and clade ----------------------------------------------------------------------------
    Param("family", "family E-value", "annotate.family_db:FAMILY_EVALUE",
          "Pfam HMM hit / nearest-reference blastp E-value that supports a family call"),
    Param("clade", "clade minimum identity", "annotate.homolog_selection:MIN_IDENTITY",
          "best-hit identity needed to assign an SSN clade"),
    Param("clade", "clade minimum support", "annotate.homolog_selection:assign_cluster(min_support)",
          "share of voting homologues that must agree on the clade"),
    # --- inducer --------------------------------------------------------------------------------------
    Param("inducer", "curated-member identity", "effector.inducer:curated_member_ion(min_identity)",
          "identity to a curated metal binder needed to transfer its ion"),
    Param("inducer", "Ligify suppression weight", "effector.inducer:LIGIFY_SUPPRESSED_WEIGHT",
          "below this gate weight an operon-chemistry call is suppressed"),
    Param("inducer", "coordination distance (A)", "structure.metal_site:_COORD_CUT",
          "donor-to-metal distance counted as coordination (structure only)"),
    Param("inducer", "MetalNet2 minimum MSA depth", "structure.metalnet:MIN_MSA_DEPTH",
          "shallower alignments abstain"),
    Param("inducer", "neighbourhood flank (genes)", "annotate.neighborhood:FLANK_GENES",
          "genes each side scanned for metal-handling neighbours"),
    Param("inducer", "neighbourhood 'near' (bp)", "annotate.neighborhood:NEAR_BP", "distance counted as near"),
    Param("inducer", "neighbourhood 'very near' (bp)", "annotate.neighborhood:VERY_NEAR_BP",
          "distance counted as very near"),
    # --- operator -------------------------------------------------------------------------------------
    Param("operator", "promoter window upstream (bp)", "annotate.genome_resolver:PROMOTER_UPSTREAM",
          "upstream promoter DNA taken per homologue"),
    Param("operator", "promoter window downstream (bp)", "annotate.genome_resolver:PROMOTER_DOWNSTREAM",
          "downstream promoter DNA taken per homologue"),
    Param("operator", "minimum SSN promoter set", "pipeline:MIN_SSN_PROMOTER_REGIONS",
          "fewest clade homologue promoters that count as a usable set"),
    Param("operator", "BLAST fallback promoters", "pipeline:BLAST_HOMOLOG_REGIONS",
          "homologue promoters requested when the clade set is too small"),
    Param("operator", "minimum mappable orthologues", "annotate.msa_homologs:MIN_MAPPABLE_ORTHOLOGS",
          "UniRef orthologues needed from the MSA route"),
    Param("operator", "seed motifs kept (top_k)", "signals.seed_selection:SeedConstruction.top_k",
          "de-novo seed motifs carried forward"),
    Param("operator", "seed width tolerance (bp)", "signals.seed_selection:SeedConstruction.tol",
          "widths merged into one seed cluster"),
    Param("operator", "PWM pseudocount per base", "signals.seed_selection:SeedConstruction.pseudocount_per_base",
          "added to every count column"),
    Param("operator", "seed selection policy", "signals.seed_selection:DEFAULT_POLICY.policy_name",
          "the one production policy"),
    Param("operator", "seed selection method", "signals.seed_selection:DEFAULT_POLICY.selection_method",
          "cluster0 = largest-support anchor width"),
    Param("operator", "seed minimum hits", "signals.seed_selection:DEFAULT_POLICY.min_hits", ""),
    Param("operator", "seed maximum hits", "signals.seed_selection:DEFAULT_POLICY.max_hits", ""),
    Param("operator", "seed keep fraction", "signals.seed_selection:DEFAULT_POLICY.keep_frac", ""),
    Param("operator", "seed trim IC (bits)", "signals.seed_selection:DEFAULT_POLICY.trim_ic",
          "flanking columns below this are trimmed"),
    Param("operator", "rescan scope", "signals.seed_selection:DEFAULT_POLICY.rescan_scope",
          "genome DNA searched for operator sites"),
    Param("operator", "rescan p-value", "signals.seed_selection:DEFAULT_POLICY.rescan_pvalue_thresh",
          "per-site cut for the genome rescan (not FDR-controlled)"),
    Param("operator", "locality priors by tier", "signals.motif_rescan:DEFAULT_PRIORS",
          "1 = own flanking intergenic, 2 = divergent, 3 = other intergenic, 4 = coding"),
    Param("operator", "operator candidates reported", "report.operators:build_operators(top_natural)",
          "ranked natural sites listed after the primary"),
    Param("operator", "autoregulation distance (bp)", "pipeline:AUTOREG_MAX_DISTANCE_BP",
          "a hit this close to the TF's own operator window counts as autoregulation"),
    Param("operator", "promoter canonical spacer (bp)", "signals.promoter:_CANONICAL_SPACER", "-35/-10 spacer"),
    Param("operator", "promoter long spacer (bp)", "signals.promoter:_LONG_SPACER",
          "spacer at or above which an operator-spanning promoter is flagged"),
    # --- regulon --------------------------------------------------------------------------------------
    Param("regulon", "regulon site p-value", "pipeline:REGULON_PVALUE", "per-site cut for regulon sites"),
    Param("regulon", "operon maximum gap (bp)", "signals.regulon:DEFAULT_OPERON_MAX_GAP",
          "largest intergenic gap inside one operon"),
    Param("regulon", "site-to-gene distance (bp)", "signals.regulon:regulated_first_genes(max_dist)",
          "furthest first gene a site can regulate"),
    Param("regulon", "operon display cap", "signals.regulon:DEFAULT_MAX_OPERONS",
          "operons reported per TF; every regulon count is a function of this cap"),
    Param("regulon", "per-hit annotations", "pipeline:PER_HIT_REGULATION_MAX_HITS",
          "rescan hits annotated with a regulated gene in the report"),
    # --- hand-off -------------------------------------------------------------------------------------
    Param("hand-off", "AF3 loci exported", "pipeline:AF3_TOP_LOCI", "distinct loci written as AlphaFold-3 jobs"),
    Param("hand-off", "AF3 operator length (bp)", "pipeline:AF3_OPERATOR_LENGTH", "genomic DNA per AF3 job"),
)


def table() -> list[dict]:
    return [{"stage": p.stage, "name": p.name, "value": p.value(), "where": p.where, "meaning": p.meaning}
            for p in PARAMS]


def markdown() -> str:
    rows = ["| stage | parameter | value | defined in | meaning |", "|---|---|---|---|---|"]
    for r in table():
        rows.append(f"| {r['stage']} | {r['name']} | `{r['value']}` | `{r['where']}` | {r['meaning']} |")
    return "\n".join(rows)


if __name__ == "__main__":
    import sys
    if "--self-test" in sys.argv:
        t = table()
        assert len(t) == len(PARAMS) and all(r["value"] is not inspect.Parameter.empty for r in t), t
        print(f"OK: {len(t)} parameters resolve from their defining modules")
    else:
        print(markdown())
