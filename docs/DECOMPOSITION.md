# How to fully decompose `run_novel` (the pipeline control plane)

`predictor/pipeline.py::run_novel` is the orchestrator: one ~570-line function that runs ~18 numbered
stages, threading ~50 local variables between them. It reads well top-to-bottom but is hard to iterate on
(reorder/insert a stage) and impossible to unit-test a stage in isolation. This document specifies the
**target design** and a **safe, incremental migration** so the science is never silently broken.

> Status: **step 0 (lighter-touch helper extraction) is done** — the self-contained try/except "islands"
> (inducer inference, target-promoter discovery, regulon enrichment) are already module-level `_stage_*`
> helpers. The remaining steps below (a real `PipelineContext` + extracting the stages that share many
> locals) are deferred until the science stabilizes, because they need a **live** golden test to verify.

## Target design

A `PipelineContext` dataclass holds all shared state; each stage is a function `_stage_x(ctx) -> None`
that reads and writes `ctx` fields; `run_novel` becomes a short, readable driver:

```python
def run_novel(seq, *, name="query", ..., cfg=None) -> dict:
    ctx = PipelineContext.from_inputs(seq, name=name, ..., cfg=cfg)
    for stage in (_stage_config, _stage_family, _stage_structure, _stage_genome, _stage_locus,
                  _stage_conservation, _stage_inducer, _stage_ligand_promoters, _stage_target_discovery,
                  _stage_rescan, _stage_motif, _stage_regulon, _stage_regulon_enrichment,
                  _stage_logos_and_perhit, _stage_mode, _stage_design, _stage_assemble, _stage_operators):
        stage(ctx)
        if ctx.error:                 # stages set ctx.error instead of `return {...}`; the driver bails
            return ctx.result()
    if ctx.write:
        _stage_write(ctx)             # bundle + AF3 hand-off
    return ctx.result()
```

### `PipelineContext` fields (grouped by the stage that first sets them)

```python
@dataclass
class PipelineContext:
    # -- inputs (from_inputs) --
    seq: str; name: str; family: str | None; organism: str | None; genome_acc: str | None
    tf_locus: tuple | None; effector: str | None; cfg: RunConfig
    # convenience mirrors of cfg (stage_config): allow_ncbi, allow_colabfold, allow_folds, write, verbose
    # -- control --
    error: str | None = None
    # -- stage_family --      rec, family
    # -- stage_structure --   struct
    # -- stage_genome --      acc, ctx (GenomeContext), auto
    # -- stage_locus --       tf_start, tf_end, tf_strand, op_center, gap_lo, gap_hi
    # -- stage_conservation - reg_start, reg_seq, homs, homs_pre, msa_info, ssn_cluster, mf, sf_seqs,
    #                         pwm, top, op_seq_start, op_seq_end
    # -- stage_inducer --     inducer_cons, ssn_cl   (may update effector)
    # -- stage_ligand_promoters -- tf_gene_obj, ligand, lgenes, lig_kws, boost_regions
    # -- stage_target_discovery -- pwm_target, target_consensus
    # -- stage_rescan --      reg, hits, rescan_hits, n_target_added, nearest, autoreg_recovered
    # -- stage_motif --       ranked, sel, glog, mean_ic
    # -- stage_regulon --     kept, op_counts, regulon      (may update effector)
    # -- stage_logos_and_perhit -- multi_logos, per_hit_reg, operator-supported genes
    # -- stage_mode --        mode
    # -- stage_design --      designed, design_note, design_method
    # -- stage_assemble --    ion, af3_jobs, doss_d
    # -- stage_operators --   operators
    # -- stage_write --       bundle
```

Notes / gotchas that make this non-trivial:
- **`effector` is mutated by three stages** (inducer, target-discovery source, regulon-enrichment). On a
  `PipelineContext` this becomes `ctx.effector` reassignment — trivial once it is a field, but it is the
  main reason a naive "extract to a function returning a tuple" is error-prone (you must return + reassign).
- **Early exits.** Stages `genome`, `locus`, `conservation` currently `return {...error...}`. In the
  target design they set `ctx.error = "..."` and the driver bails. `ctx.result()` must reproduce the exact
  error dicts (`{"tf_id","family","path":"novel", ...}`).
- **`ctx.result()`** must reproduce today's final return dict byte-for-byte (same keys, same rounding).
- **Silent-degradation trap.** Most stages are wrapped in `try/except Exception`. If a refactor drops a
  variable a stage needs, the exception is swallowed and the feature **silently stops working** — a passing
  offline self-test will NOT catch it. Hence the live golden test below is mandatory for this step.

## Migration (one stage at a time)

1. Introduce `PipelineContext` with `from_inputs` + `result()`. Keep `run_novel` as-is but build a `ctx`
   at the top and, at the end, assert `ctx.result() == <today's dict>` on the golden TF (temporary guard).
2. Move stages **bottom-up** (assemble/operators/write first — they only read state), one per commit. After
   each: run the live golden test (below) and confirm the dossier is identical.
3. Convert the early-exit stages last; replace their `return {...}` with `ctx.error = ...` and verify the
   error paths with a deliberately bad input (e.g. an unlocatable locus).
4. Delete the temporary `assert` guard once every stage is migrated.

## Verification — the live golden test (required)

The offline self-test (`python -m predictor.pipeline --self-test`) only covers helpers, NOT `run_novel`.
To verify a decomposition step, diff a **real run** before and after:

```bash
# prerequisites (once): tokens in .env (BIOHUB_TOKEN, NCBI_TOKEN) and engines on PATH
tfop setup && tfop audit                    # both green
# pick a small, well-behaved TF (e.g. an ArsR/SmtB) and a known genome:
tfop predict CzrA.fasta --name GOLD --family "ArsR/SmtB" --organism <GCF_...> --no-fold
cp results/jobs/GOLD/dossier.json /tmp/dossier.before.json
# ... apply one decomposition step ...
tfop predict CzrA.fasta --name GOLD --family "ArsR/SmtB" --organism <GCF_...> --no-fold
diff <(jq -S . /tmp/dossier.before.json) <(jq -S . results/jobs/GOLD/dossier.json)   # must be empty
```

`--no-fold` keeps it deterministic and free (no biohub quota). An empty `diff` across a couple of TFs
(one ArsR, one MerR) is the acceptance criterion for each step.
