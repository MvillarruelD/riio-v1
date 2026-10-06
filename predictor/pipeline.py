"""
pipeline.py -- THE operator-discovery pipeline (the "final" entry).

One sequence in -> (1) every putative operator/promoter of the TF in its genome, (2) the binding motif,
(3) the best DNA the TF could bind (a design candidate) -- plus promoter position, mode of regulation,
the regulon, the inducer, and an AF3 hand-off bundle.

EVERY FAMILY TAKES THE SAME PATH. There is one route -- conservation-first discovery -- and it is applied
identically to all of `annotate.ssn_clusters.FAMILIES` and to any family outside them:

    1. classify the family, fold the apo dimer if a backend is reachable
    2. acquire the genome, locate the TF gene, take its autoregulatory window
    3. gather homolog promoters: the TF's SSN cluster when its family has one, else BLAST similarity;
       deepen with an MSA expansion
    4. discover the motif in that window, then rescan the whole genome with it
    5. build the operator logo, the regulon, the inducer consensus, the mode of regulation

Earlier revisions routed MerR and ArsR/SmtB differently from everything else -- a structural arm for MerR,
conservation-only for ArsR, neutral text for the rest -- and carried per-family claims about which
in-silico signals could be trusted. That routing is gone: no family gets a discovery step, a scoring rule,
a width prior or a mechanism label that another family does not.

Run (from an environment with the package installed):
    python -m predictor.pipeline --self-test                            # offline smoke test
    tfop predict mytf.fasta --name MyTF --organism GCF_000195955.2      # the supported entry point
"""
from __future__ import annotations

import argparse
import os
import sys
import uuid
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

_PRED = Path(__file__).resolve().parent            # .../predictor
_REPO = _PRED.parent


# Run parameters that were bare literals inside run_novel. Named here so predictor.params can report
# them and a reader can find them; the values are exactly the ones the canonical runs used.
TF_LOCUS_TOP_N = 3                  # genome hits considered when locating the TF's own gene
MIN_SSN_PROMOTER_REGIONS = 4        # fewest SSN-clade homolog promoters that count as a usable set
BLAST_HOMOLOG_REGIONS = 12          # homolog promoters requested from the BLAST-similarity fallback
AUTOREG_MAX_DISTANCE_BP = 80        # a rescan hit this close to the TF's own operator window = autoregulation
REGULON_PVALUE = 1e-4               # per-site p-value for sites that seed the regulon
PER_HIT_REGULATION_MAX_HITS = 60    # rescan hits annotated with a regulated gene in the report
AF3_TOP_LOCI = 5                    # distinct loci exported as AlphaFold-3 jobs
AF3_OPERATOR_LENGTH = 60            # bp of genomic DNA per exported AF3 job
if __package__ in (None, ""):                      # support `python predictor/pipeline.py`
    sys.path.insert(0, str(_REPO))

from predictor.annotate import context, genome_resolver, genome_mirror              # noqa: E402
from predictor.annotate import homolog_regions as hr                                # noqa: E402
from predictor.annotate import tf_record                                            # noqa: E402
from predictor.signals import motif_finder, operator_logo as ol, motif_rescan, coords, promoter  # noqa: E402
from predictor.signals import seed_selection as ss                                  # noqa: E402
from predictor.structure import af3_export                                          # noqa: E402
from predictor.report.dossier import TFDossier                                      # noqa: E402
from predictor.report import outputs                                                # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

# --------------------------------------------------------------------------- helpers
def _read_seq(s: str) -> str:
    from predictor.input_validation import read_protein
    return read_protein(s)


def _failure(name, family, stage, message, **extra):
    """One machine-readable failure shape for every mandatory production stage."""
    return {
        "tf_id": name,
        "family": family,
        "path": "novel",
        "ok": False,
        "required_stages_complete": False,
        "failed_stage": stage,
        "error": str(message),
        **extra,
    }


def _designed_operator(logo) -> str | None:
    """The idealised strongest operator = the most-informative base at every column of the genomic motif."""
    if logo is None or logo.pwm is None or logo.pwm.shape[1] == 0:
        return None
    return "".join("ACGT"[i] for i in np.asarray(logo.pwm).argmax(axis=0))


# This package ships one prediction path: run_novel() below. There is no separate route for
# characterised TFs; benchmarks call the same path.


# --------------------------------------------------------------------------- run_novel stage helpers
# These are the self-contained "islands" of run_novel (each already best-effort/try-except), pulled out as
# named functions so the orchestrator is shorter and each is independently testable. Behaviour is unchanged
# -- bodies are identical to the inline blocks; they take explicit inputs and return their outputs (note the
# `effector` fill-in propagates back via the return value). See docs/DECOMPOSITION.md for the full plan.
def _infer_inducer(seq, family, ctx, tf_start, effector, *, allow_ncbi, verbose,
                   allow_metalnet=True, protein_acc=None,
                   tf_id=None, allow_colabfold=False, fold_cif=None, op_cif=None):
    """Stage 4.5: unified inducer inference (SSN cluster + coordination gate + MetalNet + Ligify
    [+ structural]). Returns (inducer_cons, ssn_cl, effector); effector may be filled in from a NAMED
    source -- note MetalNet is deliberately NOT among those sources: it establishes that a metal site
    exists, not which metal, so it must never name the analyte.

    `fold_cif`/`op_cif` (optional): a holo fold with bound metals (and, optionally, its operator
    complex) lets the coordination gate type the metal sites and contribute the REGULATORY-site metal
    for two-site regulators. Absent (the default, and the structure-off batch), the gate is unchanged."""
    inducer_cons, ssn_cl, failure = None, None, None
    try:
        from predictor.effector import inducer as inducer_mod
        tf_gene = next((gg for gg in ctx.genes if gg.start <= tf_start < gg.end), None)
        # `allow_metalnet=False` injects an abstaining stub -- exactly what a machine without the
        # engine produces, so the source reports "no opinion" in the open rather than the path changing.
        _mn_fn = (lambda _s: None) if not allow_metalnet else None
        # `protein_acc` is the TF's own RefSeq id when the locus resolved -- the join key for the
        # published Ligify database. Never substituted from a neighbour or a nearest reference: a wrong
        # accession here would attribute another protein's published biosensor ligand to this one.
        _acc = protein_acc or (getattr(tf_gene, "protein_id", None) if tf_gene is not None else None)
        inducer_cons = inducer_mod.infer_inducer(seq, family, genome_ctx=ctx, tf_gene=tf_gene, allow_ncbi=allow_ncbi,
                                                 tf_id=tf_id, allow_colabfold=allow_colabfold,
                                                 protein_acc=_acc, fold_cif=fold_cif, op_cif=op_cif,
                                                 metalnet_fn=_mn_fn, verbose=verbose)
        ssn_cl = next((c.evidence.get("cluster") for c in inducer_cons.calls
                       if c.source == "ssn_cluster"), None)
        if effector is None:        # derive the analyte from a NAMED source (not the class-level fallback)
            named = next((c.ligand for c in inducer_cons.calls
                          if c.ligand and c.source in ("ssn_cluster", "ligify_db", "ligify")), None)
            if named:
                effector = named
        if verbose:
            print(f"inducer: top={inducer_cons.top} agreement={inducer_cons.agreement:.2f} "
                  f"gate={inducer_cons.coordination_gate} cluster={ssn_cl} "
                  f"(sources: {', '.join(c.source for c in inducer_cons.calls)})")
    except Exception as e:
        failure = f"{type(e).__name__}: {e}"
        if verbose:
            print(f"  inducer inference skipped ({type(e).__name__}: {e})")
    return inducer_cons, ssn_cl, effector, failure


# Stage 6.55 (removed): the reconstructed regulon used to be re-fused into the inducer call as a
# seventh source, naming the ion from its members' transporter substrate specificity. It is reached
# THROUGH the operator hits, whose site-level precision is about 2 %, so it is a tier-2 product and
# may not raise the confidence of a tier-1 claim. Measured over 139 candidates, removing it changes
# 16 headline calls and moves the metal fraction by ZERO -- every change within-class. The regulon is
# still reconstructed and still reported, as a hypothesis with a score.


def _construct_seed_clusters(reg_seq, *, family, homolog_regions, construction):
    """Route private benchmark construction values through the one production generator path."""
    policy = ss.SeedPolicy(
        **{
            **ss.DEFAULT_POLICY.__dict__,
            "pseudocount_per_base": construction.pseudocount_per_base,
        }
    )
    candidates = motif_finder.predict(
        reg_seq,
        family=family,
        homolog_regions=homolog_regions or None,
        top_k=construction.top_k,
        family_informed=False,
        use_family_seed=False,
        score_mode=construction.score_mode,
    )
    clustered_sequences = motif_finder.seed_clusters(
        candidates,
        tol=construction.tol,
        max_clusters=None,
    )
    clusters = ss.prepare_seed_clusters(candidates, clustered_sequences, policy=policy)
    return candidates, clustered_sequences, clusters, policy


def _seed_selection_report(decision, trials, construction):
    report = decision.to_dict(trials)
    report["construction"] = construction.to_dict()
    return report


# --------------------------------------------------------------------------- novel-TF path (0 folds)
def _run_novel_impl(seq: str, *, name: str = "query", family: str | None = None,
                    organism: str | None = None, genome_acc: str | None = None, tf_locus=None,
                    effector: str | None = None, cfg=None,
                    seed_construction: ss.SeedConstruction) -> dict:
    """The cold-start path for a TF that was never folded/curated. Conservation-first discovery, identical
    for every family. Emits operators + motif + mode + design candidate + the job bundle.

    `cfg` (RunConfig) carries the environment capabilities of this run -- whether a fold backend, a
    ColabFold MSA and the network are reachable -- not a choice of algorithm. One config is built here when
    the caller does not supply one."""
    original_input = seq
    requested = {
        "name": name,
        "family": family,
        "organism": organism,
        "genome_accession": genome_acc,
        "tf_locus": list(tf_locus) if tf_locus else None,
        "effector": effector,
    }
    input_kind = "normalized_sequence"
    try:
        candidate_path = Path(str(original_input))
        if len(str(original_input)) < 1024 and candidate_path.is_file():
            input_kind = "fasta_file"
    except (OSError, ValueError):
        pass
    run_started = datetime.now(timezone.utc).isoformat()
    run_id = uuid.uuid4().hex
    try:
        seq = _read_seq(seq)
    except ValueError as exc:
        return _failure(name, family, "input", exc)
    try:
        outputs.job_path(name)
    except ValueError as exc:
        return _failure(name, family, "output", exc)

    # ---- 0. run config: what this environment can reach (network / fold backend / ColabFold) ---------
    try:
        from predictor import resources as _res          # WS7: load API tokens from the repo-local .env
        _res.load_tokens()
    except Exception:
        pass
    from predictor.config import RunConfig
    if cfg is None:
        cfg = RunConfig()
    allow_ncbi, allow_colabfold = cfg.allow_online, cfg.allow_colabfold
    write, verbose = cfg.write, cfg.verbose
    if write and (outputs.job_path(name) / "operators_ranked.json").exists():
        return _failure(
            name, family, "output",
            f"completed job already exists: {outputs.job_path(name)}; use a new --name or archive it first",
        )

    # ---- 1. family ----------------------------------------------------------------------------------
    rec = None
    if not family:
        rec = tf_record.classify(seq)
        family = rec.family
        if not family:
            evidence = "; ".join(getattr(rec, "flags", []) or [])
            return _failure(name, None, "family", "could not classify the query as a supported TF family",
                            family_evidence=evidence)
    else:
        try:
            from predictor.input_validation import validate_family
            family = validate_family(family)
        except ValueError as exc:
            return _failure(name, family, "family", exc)
    from predictor.input_validation import family_support as _family_support
    support_level = _family_support(family)
    if verbose:
        print(f"=== pipeline (novel): {name} [{family or 'family?'}] ===")

    # ---- 2. EARLY structure step (WS1): apo dimer fold (quota-aware) + AFDB-monomer integrity QC ----
    # Folds the apo homodimer up front (when cfg.allow_folds + a backend token), cross-checks it against
    # the free AFDB monomer, and makes the structure available to the ligand step + operator design. With
    # no backend / no folds it degrades to a request descriptor + AFDB-only QC (never hard-fails).
    struct = None
    if cfg.allow_folds or cfg.allow_online:
        try:
            from predictor.structure import early_structure as _es
            if rec is None:
                rec = tf_record.classify(seq)
            # Sequence-only inputs have no query UniProt accession. Never substitute the nearest
            # reference: that would turn query fold-integrity QC into an unlabelled homolog comparison.
            struct = _es.run_structure_step(seq, family, uniprot=getattr(rec, "uniprot", None),
                                            name=name, cfg=cfg, verbose=verbose)
        except Exception as e:
            if verbose:
                print(f"  structure step skipped ({type(e).__name__}: {e})")

    # ---- 2. genome: on-demand fetch -> sequence-resolution ------------------------------------------
    # The run is agnostic of any pre-gathered local data: the named genome is fetched (and cached) from
    # NCBI rather than read out of a mirror that only exists on a development machine.
    acc = genome_acc or organism
    ctx = None
    if acc:
        # `ensure_genome` is mirror-FIRST: it returns an already-cached genome before it considers
        # downloading, and `allow_download` gates only the download. Gating the whole call on
        # `allow_ncbi` therefore made an offline run refuse a genome sitting in its own cache --
        # "could not acquire the explicitly requested genome" for a file already on disk. Only the
        # fetch is online; reading the cache never was.
        fna, _ = genome_mirror.ensure_genome(acc, allow_download=bool(allow_ncbi))
        if fna:
            ctx = context.genome_context_from_mirror(Path(fna).stem)
    # EXPLICIT --organism / --genome-acc: honor it ABSOLUTELY. When the user names a nuccore accession
    # (NC_/NZ_/CP_...), fetch THAT EXACT genome via EFetch -- never silently search a different genome.
    # This is the safeguard for the failure mode where ensure_genome/Datasets is unavailable and the
    # sequence-based auto-locus would otherwise place the TF in a homolog's genome (e.g. Shigella vs the
    # E. coli K-12 the user asked for). Only the Datasets/assembly path (GCF_/GCA_) still needs auto-locus.
    if ctx is None and acc and allow_ncbi and not str(acc).upper().startswith(("GCF_", "GCA_")):
        try:
            # gbwithparts (not plain gb): RefSeq complete genomes are served as CON records whose ORIGIN
            # sequence is referenced, not inlined -- plain "gb" yields a feature-only record and
            # context.from_genbank then raises UndefinedSequenceError. gbwithparts expands it to full sequence.
            gb = genome_resolver.fetch_region(acc, rettype="gbwithparts")
            ctx = context.from_genbank(gb)
            if verbose:
                print(f"  genome: honored --organism {acc} via EFetch ({len(ctx.genes)} genes)")
        except Exception as e:
            if verbose:
                print(f"  could NOT fetch the requested genome {acc} ({type(e).__name__}: {e}); "
                      f"falling back to sequence-based auto-locus")

    # AUTHORITATIVE fallback (fully automatic, no hardcodes): when no genome was named/acquired -- the
    # plasmid / multi-contig / wrong-strain cases that used to need a hardcoded mirror_acc+tf_locus -- ask
    # NCBI IPG for the exact source contig + coords and work IN THAT CONTIG'S frame (avoids the concatenation
    # mismatch). Gated by allow_ncbi (remote blastp). Sets `auto` so the locus step reuses the same answer.
    auto = None
    if ctx is None and acc:
        msg = (f"could not acquire the explicitly requested genome {acc}; refusing to substitute a "
               "different sequence-resolved genome")
        if verbose:
            print(msg)
        return _failure(name, family, "genome", msg)
    if ctx is None and allow_ncbi:
        try:
            from predictor.annotate import auto_locus as al
            auto = al.auto_locus(seq, key=name, verbose=verbose)
            if auto:
                gb = genome_resolver.fetch_region(auto["accession"], rettype="gbwithparts")
                ctx = context.from_genbank(gb)
                # surface the alternative source genomes so the user can pin one with --organism if the
                # auto-pick is not the intended strain (the operator search runs in ctx.accession's frame).
                if verbose and auto.get("alternatives"):
                    alts = ", ".join(f"{a['accession']} ({a['organism']})" for a in auto["alternatives"][:4])
                    print(f"  auto_locus picked {auto['accession']} ({auto.get('organism')}); "
                          f"alternatives carrying this protein: {alts}")
                    print("  -> re-run with --organism <accession> to force a specific source genome")
        except Exception as e:
            if verbose:
                print(f"  auto_locus fallback failed ({type(e).__name__}: {e})")
            auto = None
    if ctx is None:
        msg = f"could not acquire a genome for {name} (organism={organism}, acc={genome_acc}); " \
              f"pass --organism <RefSeq/GCF accession> and check network access"
        if verbose:
            print(msg)
        return _failure(name, family, "genome", msg)
    if verbose:
        print(f"genome: {ctx.accession} ({len(ctx.sequence):,} bp, {len(ctx.genes)} genes)")

    # ---- 3. TF locus + autoregulatory operator window ----------------------------------------------
    if tf_locus:
        tf_start, tf_end, tf_strand = tf_locus
    elif auto and auto.get("accession") == ctx.accession:
        tf_start, tf_end, tf_strand = auto["start"], auto["end"], auto["strand"]   # authoritative IPG locus
    else:
        # Resolve the locus AGAINST THE FETCHED GENOME ONLY (a one-shot DB built from ctx), never a global
        # pre-gathered mirror DB -- so the answer depends only on the genome this run actually acquired.
        import shutil as _sh
        import tempfile as _tf
        _td = _tf.mkdtemp(prefix="locus_")
        try:
            _fa = Path(_td) / f"{ctx.accession}.fna"
            _fa.write_text(f">{ctx.accession}\n{ctx.sequence}\n", encoding="utf-8")
            loc = genome_resolver.resolve(seq, genome_fasta=str(_fa), top_n=TF_LOCUS_TOP_N)
        finally:
            _sh.rmtree(_td, ignore_errors=True)
        if not loc:
            return _failure(name, family, "locus",
                            "could not locate the TF gene in the genome (tblastn miss); pass --tf-locus")
        tf_start, tf_end, tf_strand = loc[0].tf_start, loc[0].tf_end, loc[0].tf_strand
    if getattr(ctx, "circular", False) and (tf_start < 400 or tf_end > len(ctx.sequence) - 400):
        return _failure(
            name, family, "circular_origin",
            "the TF promoter may cross the circular origin, which the linear operator coordinate schema "
            "cannot represent safely",
            genome=ctx.accession,
        )
    op_center, gap_lo, gap_hi = context.locate_operator_window(ctx, tf_start, tf_end, tf_strand)
    if verbose:
        print(f"TF locus {tf_start}-{tf_end} ({tf_strand}); autoregulatory operator anchor ~{op_center} "
              f"(gap {gap_lo}-{gap_hi})")

    # ---- 4. conservation-first discovery (Snowprint's evidence type + family palindrome priors) -----
    if tf_strand == "-":
        reg_start = gap_lo
        reg_seq = ctx.sequence[reg_start: min(len(ctx.sequence), reg_start + 400)]
    else:
        reg_start = max(0, gap_hi - 400)
        reg_seq = ctx.sequence[reg_start: gap_hi]
    # Homolog autoregulatory-promoter set (the conservation signal). Every family follows the same two
    # steps, tried in COST order:
    #   * SSN-cluster promoters -- isofunctional cluster members, for any family that has an SSN;
    #   * BLAST-similarity homologs (`collect_homolog_regions`) -- for families with no SSN, and as the
    #     fallback when the cluster yield is thin (< 4 regions).
    # The SSN path is tried FIRST because the BLAST path costs a remote NCBI blastp plus up to 12 promoter
    # EFetches (multi-minute); running it unconditionally and then discarding it was pure wasted work.
    # Cluster assignment: direct accession lookup (O(1) when the TF's UniProt is in the SSN map), then
    # high-identity MMseqs hit (>= 0.90), then voting fallback.
    homs, _ssn_ok = [], False
    homolog_collection = {"selected_route": None, "attempts": []}
    from predictor.annotate import ssn_clusters as _ssn_mod
    if cfg.allow_online and family in set(_ssn_mod.FAMILIES):
        try:
            from predictor.annotate import msa_homologs as _mh
            from predictor.annotate import homolog_selection as _hsel
            # `rec.uniprot` belongs to the nearest curated neighbour, not necessarily this query.
            # Passing it here turns the O(1) accession lookup into a circular cluster assignment.
            _asg = _hsel.assign_cluster(seq, family)
            _cid = _asg.cluster_id
            if _cid:
                _creg_record = _mh.cluster_promoter_regions(seq, _cid, verbose=verbose)
                _creg = _creg_record.get("regions", [])
                homolog_collection["attempts"].append({
                    "route": "ssn_cluster",
                    "selected": len(_creg) >= MIN_SSN_PROMOTER_REGIONS,
                    "cluster_assignment": {
                        "cluster_id": _cid,
                        "method": _asg.method,
                        "support": float(_asg.support or 0.0),
                    },
                    "record": _creg_record,
                })
                if len(_creg) >= MIN_SSN_PROMOTER_REGIONS:
                    homs, _ssn_ok = _creg, True
                    homolog_collection["selected_route"] = "ssn_cluster"
                    if verbose:
                        print(f"  SSN cluster: {_cid} [{_asg.method}, support {_asg.support:.2f}] -> "
                              f"{len(_creg)} promoter regions (redundant BLAST set skipped)")
                elif verbose:
                    print(f"  SSN cluster {_cid} yield thin ({len(_creg)}); falling back to BLAST homologs")
        except Exception as e:
            if verbose:
                print(f"  SSN-cluster homologs skipped ({type(e).__name__}: {e})")
    if not _ssn_ok:                                    # fallback: BLAST-similarity homolog promoters
        # The local genome mirror runs; the remote NCBI blastp -> IPG -> EFetch fallback is OFF.
        # It works (its bytes/str parsing defect was fixed 2026-09-05, tests/test_ipg_parse.py), but it
        # costs 30 min to 6 h per candidate against ~5 min for the UniRef MSA path that now supplies
        # this axis, so re-enabling it is a cost decision for whoever plans a run (CHANGELOG
        # 2026-09-05). It is hard-wired False, not a toggle: an environment switch here would select
        # which algorithm builds the conservation set, which the single-path rule forbids
        # (tests/test_production_hardening.py::test_no_environment_variable_selects_an_algorithm).
        try:
            homs = hr.collect_homolog_regions(seq, n=BLAST_HOMOLOG_REGIONS, allow_ncbi=False)
            homolog_collection["selected_route"] = "blast_similarity"
            homolog_collection["attempts"].append({
                "route": "blast_similarity", "selected": True,
                "allow_remote_ncbi": False, "n_regions": len(homs),
            })
        except Exception as e:
            if verbose:
                print(f"  homolog collection failed ({e}); proceeding palindrome-only")
            homs = []
            homolog_collection["attempts"].append({
                "route": "blast_similarity", "selected": False,
                "allow_remote_ncbi": False, "error_type": type(e).__name__,
            })
    homs_pre = list(homs)    # WS3: homolog promoter set BEFORE MSA expansion (feeds the 'homolog' logo)
    # MSA homolog expansion (augment-when-sparse): deepens the conservation set via an MSA when the
    # SSN/BLAST homolog set is thin (Snowprint failure mode #3). AF3-free MSA via ColabFold.
    msa_info = {"local": len(homs), "msa_added": 0, "source": None}
    if allow_ncbi or allow_colabfold:
        try:
            from predictor.annotate import msa_homologs as mh
            homs, msa_info = mh.expand_homolog_regions(homs, tf_id=name, seq=seq,
                                                       allow_colabfold=allow_colabfold, verbose=verbose)
        except Exception as e:
            if verbose:
                print(f"  MSA expansion skipped ({type(e).__name__}: {e})")
    homolog_collection["msa"] = msa_info
    homolog_collection["pre_msa_count"] = len(homs_pre)
    homolog_collection["retained_count"] = len(homs)
    # top_k=8. Generation is NOT the bottleneck: measured over the 101 RegulonDB sites located in the
    # local genome (the generation diagnosis, recorded in the analysis repo before commit 5815559),
    # the finders already produce a candidate covering the true site for **90 of 101** -- 0 sites
    # fall outside the geometry prior and
    # 0 fail to be a detectable dyad. Run 6 therefore raised top_k to 24 to stop truncating the answer
    # away, and it worked at RETENTION: seed-pool site coverage 39% -> 63% overall, 44% -> 73% on the
    # metal families.
    #
    # It did not survive selection, and that is why this is back at 8. At top_k=24 the extra candidates
    # move the modal length, so a DIFFERENT candidate builds the PWM: IscR and MntR are lost while only
    # ArsR is gained (7/17 -> 6/17 complete-panel primary-site recovery), and the held-out result is
    # identical either way (3/10 non-anchored under both leave-one-TF-out and leave-one-family-out).
    # Extra retention that no truth-blind selector can convert is cost without benefit -- and it costs
    # ~1.4x the rescans, since cost scales with the number of width clusters (measured mean 3.97 at
    # top_k=8 vs 5.55 at 24). Full working: the 2026-08-19 seed decision, recorded in the analysis
    # repo before commit 5815559.
    mf, seed_sets, seed_clusters, seed_policy = _construct_seed_clusters(
        reg_seq,
        family=family,
        homolog_regions=homs,
        construction=seed_construction,
    )
    # Literature-derived family PWMs are a separate candidate route, never fused into the de-novo
    # seed or counted as independent validation. Leave the named query out wherever TF-level
    # provenance permits; incomplete source mappings simply make this route unavailable.
    literature_candidates = []
    try:
        literature_candidates = motif_finder.family_seeded(
            reg_seq, family, exclude_tf=name,
        )[:8]
    except Exception as e:
        if verbose:
            print(f"  literature-seed route skipped ({type(e).__name__}: {e})")
    if not seed_sets:
        return _failure(name, family, "motif_seed", "no motif_finder candidates in the autoregulatory window",
                        genome=ctx.accession)
    if not seed_clusters:
        return _failure(name, family, "motif_seed", "no motif_finder candidates in the autoregulatory window",
                        genome=ctx.accession)
    sf_seqs = seed_clusters[0].prepared.sequences
    if verbose:
        print(f"seed widths: {[cluster.anchor_w for cluster in seed_clusters]} "
              f"({[cluster.n_effective for cluster in seed_clusters]} effective seed(s) each) -- "
              f"each is rescanned; the genome decides, not the anchor")
    top = mf[0]
    # Genome-absolute coordinates of the top operator (width = what motif_finder actually found)
    op_seq_start = top.start + reg_start
    op_seq_end = top.end + reg_start

    # ---- 4.5 unified inducer inference (SSN cluster + coordination gate + Ligify [+ structural]) ----
    inducer_cons, ssn_cl, effector, inducer_failure = _infer_inducer(
        seq, family, ctx, tf_start, effector,
        allow_ncbi=allow_ncbi, verbose=verbose,
        allow_metalnet=getattr(cfg, "allow_metalnet", True),
        tf_id=name, allow_colabfold=allow_colabfold,
    )
    if inducer_cons is None:
        return _failure(name, family, "inducer", inducer_failure or "inducer consensus was not produced",
                        genome=ctx.accession)

    # ---- 4.6 candidate regulated genes (reported, never fed back into discovery) --------------------
    # Use the inferred inducer to flag likely-regulated genes (transport/efflux/reductase/the metal name).
    # These are REPORTED alongside the operators and used to confirm hits after the fact; they never enter
    # operator discovery. `lgenes` are inferred FROM THE EFFECTOR CALL, so feeding their promoters back in
    # -- as a locality boost on the rescan, or as the substrate for a second motif -- would let the effector
    # call raise the score of operators at precisely the genes the downstream validation then looks for. An
    # operator hit at a metal gene has to be a finding, not a consequence of having up-weighted that
    # promoter. Annotation-light: ligand_genes uses GFF product names (offline).
    from predictor.report import ligand_regulon as lr
    tf_gene_obj = next((gg for gg in ctx.genes if gg.start <= tf_start < gg.end), None)
    ligand = lr.infer_ligand(ctx, tf_gene_obj) if (tf_gene_obj is not None and allow_ncbi) else []
    lgenes, lig_kws = lr.ligand_genes(ctx, effector=effector, ligand_names=[l["ligand"] for l in ligand])
    if verbose:
        print(f"candidate handler genes: {len(lgenes)} (kw={lig_kws[:6]}); "
              f"ligand={ligand[0]['ligand'] if ligand else '(none)'} "
              f"-- reported only, not used to rank operators")

    # ---- 5. genome-WIDE rescan -> ALL putative operators -------------------------------------------
    reg = next((g.name or g.protein_id for g in ctx.genes if g.start <= tf_start <= g.end), None)
    # One rescan PER SEED WIDTH, then the genome picks the winner. The old code rescanned once with the
    # single surviving width, so a mis-chosen anchor was unrecoverable -- every downstream artifact, up
    # to the exported AF3 jobs, pointed at the wrong DNA. Each width now gets a fair genome-wide test.
    # Only the winner's hits go downstream (the logo, the regulon and `counts_from_seqs` all assume a
    # single width), but every cluster's outcome is recorded so the choice is auditable.
    #
    # Selection is a replayable pure policy, and the promoted one is `cluster0`: take the trial built
    # from the largest-support anchor width. It was selected by the held-out harness in 139 of 139
    # applications; the run-6 `mean_ic` rule it replaces was selected in none, and scores 1/17 against
    # this policy's 7/17 on the complete panel. Mean IC is a SELF-FIT score -- a cluster's PWM creates
    # the hits its own mean IC is then computed from -- so a one-sequence matrix is maximally sharp,
    # matches only near-exact genomic copies, and wins a criterion that measures nothing external.
    # The mean-IC, support and shrinkage rules stay in `seed_selection` as named benchmark policies.
    seed_trials = []
    for cluster in seed_clusters:
        _h = motif_rescan.rescan_genome(
            ctx,
            pwm=cluster.pwm,
            regulator=reg,
            scope=seed_policy.rescan_scope,
            pvalue_thresh=seed_policy.rescan_pvalue_thresh,
        )
        _rh = [coords.from_genome_hit(x) for x in _h]
        seed_trials.append(ss.build_trial(cluster, _rh, policy=seed_policy))
    decision = ss.choose_trial(seed_trials, policy=seed_policy)
    best = decision.chosen
    seed_trial_report = _seed_selection_report(decision, seed_trials, seed_construction)
    if verbose and len(seed_trials) > 1:
        for trial in seed_trials:
            audit = trial.audit
            print(f"  seed w={trial.anchor_w:>3}: {audit['n_effective']:>2} effective/"
                  f"{audit['n_input']:>2} input seed(s) {audit['n_hits']:>4} hits  "
                  f"meanIC={audit['mean_ic']:.3f}  {audit['consensus']}"
                  f"{'   <== chosen' if trial is best else ''}")
        print(f"  -> {seed_trial_report['chosen_why']}")
    rescan_hits = list(best.hits)
    if not rescan_hits:
        return _failure(name, family, "genome_rescan", "the selected motif produced no genome-wide hits",
                        genome=ctx.accession)
    # Report the exact modal-length sequences that built the chosen PWM, not the raw mixed-width pool.
    sf_seqs = best.cluster.prepared.sequences
    dists = [abs(h.dyad_center - op_center) for h in rescan_hits]
    nearest = min(dists) if dists else None
    # AUTOREGULATORY-operator recovery: is there a strong operator at the TF's OWN promoter? (a FINDING --
    # a "no" means the TF does not strongly autoregulate, e.g. CueR, NOT that discovery failed.)
    autoreg_recovered = bool(nearest is not None and nearest <= AUTOREG_MAX_DISTANCE_BP)

    # ---- 6. the binding motif (data-driven hit selection) ------------------------------------------
    ranked = [getattr(h, "seq", "") for h in sorted(rescan_hits, key=lambda h: -(h.score or 0))
              if getattr(h, "seq", "")]
    sel = ol.select_by_motif(ranked, min_hits=3)
    glog = sel.logo
    mean_ic = (glog.total_ic / glog.pwm.shape[1]) if (glog.pwm is not None and glog.pwm.shape[1]) else 0.0

    # ---- 6.5 regulon proposal (operator PWM -> operons the TF likely drives) ------------------------
    kept = ranked[:sel.n_keep]
    op_counts = ol.counts_from_seqs(kept) if kept else None
    regulon_stats: dict = {}
    regulon = (lr.regulon(ctx, op_counts, regulator=reg, pvalue_thresh=REGULON_PVALUE,
                          stats=regulon_stats) if op_counts is not None else [])
    if verbose:
        # The list is capped at DEFAULT_MAX_OPERONS, which nearly every candidate reaches, so say when
        # the reported number IS the cap rather than a count of what was found.
        cap_note = f" (capped; {regulon_stats['n_operons_found']} passed)" if regulon_stats.get("truncated") else ""
        print(f"regulon: {len([o for o in regulon if 'first_gene' in o])} operon(s) proposed{cap_note} "
              f"(from {len(kept)} motif-coherent sites)")

    # ---- 6.55 (removed) --------------------------------------------------------------------------
    # The reconstructed regulon used to be re-fused into the inducer call as a seventh source, naming
    # the ion from its members' transporter substrate specificity. It no longer is, and this note
    # stands in its place because the absence is a decision, not an oversight.
    #
    # The regulon is reached THROUGH the operator hits, whose site-level precision is about 2 %. It is
    # a tier-2 product, and a tier-2 product may not raise the confidence of a tier-1 claim. Measured
    # over the 139 replayable candidates of the last run, removing it changes 16 headline calls and
    # moves the metal fraction by ZERO: every change is within-class, typically a named ion falling
    # back to "divalent metal (ion unresolved)". So the cost is ion specificity, not metal recall,
    # and that is the honest trade.
    #
    # The regulon is still reconstructed and still reported -- as a hypothesis with a score.

    # ---- 6.6 WS3 three motif-aware logos + WS4 per-hit regulation (additive; best-effort) -----------
    multi_logos = {}
    try:
        ml, aligned = ol.three_logos(kept, homs_pre, homs, return_aligned=True)

        def _lg(o, key):
            base = ({"consensus": o.consensus, "total_ic": round(float(o.total_ic), 2), "n_seqs": o.n_seqs,
                     "width": int(o.pwm.shape[1]), "pwm": o.pwm.tolist()} if o.pwm is not None else {})
            base["aligned_seqs"] = aligned.get(key, [])      # the aligned operator instances behind the logo
            return base
        multi_logos = {"tf_only": _lg(ml.tf_only, "tf_only"), "homolog": _lg(ml.homolog, "homolog"),
                       "gathered": _lg(ml.gathered, "gathered")}
        if verbose:
            print(f"three logos: tf_only IC={ml.tf_only.total_ic:.1f} ({ml.tf_only.consensus}); "
                  f"homolog IC={ml.homolog.total_ic:.1f}; gathered IC={ml.gathered.total_ic:.1f}")
    except Exception as e:
        if verbose:
            print(f"  WS3 three_logos skipped ({type(e).__name__}: {e})")
    per_hit_reg, confirmed_genes = [], []
    try:
        # classify more hits (was 25) so the higher-q target operons -- whose operators can score below the
        # top 25 -- also get a promoter/TSS/footprint record (feeds the per-operator zoom callouts).
        if verbose:
            print(f"promoter positioning: classifying up to {min(60, len(rescan_hits))} operator "
                  "neighbourhoods...", flush=True)
        per_hit_reg = lr.per_hit_regulation(ctx, rescan_hits, max_hits=PER_HIT_REGULATION_MAX_HITS)
        confirmed_genes = lr.confirm_regulated_genes(per_hit_reg, lgenes)
        if verbose:
            nsup = sum(1 for c in confirmed_genes if c.get("operator_supported"))
            nconf = sum(1 for c in confirmed_genes if c.get("confirmed"))
            print(f"per-hit regulation: {len(per_hit_reg)} hit(s) classified; "
                  f"{nsup}/{len(confirmed_genes)} inducer-inferred gene(s) have a ranked operator hit; "
                  f"{nconf} pass q<=0.05")
    except Exception as e:
        if verbose:
            print(f"  WS4 per-hit regulation skipped ({type(e).__name__}: {e})")

    # ---- 7. promoter position + promoter occlusion ---------------------------------------------------
    # Identical for every family: `classify_mode` measures where the operator sits relative to the -35/-10
    # boxes and the TSS and returns (a) `occludes`, the promoter elements the operator physically covers,
    # and (b) a descriptive geometry tag for the same footprint ('promoter-core-overlap',
    # 'elongated-spacer-overlap', 'no-promoter-overlap'). It never asserts a mechanism -- which is why the
    # dossier key is `promoter_occlusion` and not the old `mode_of_regulation`.
    mode = {}
    try:
        mc = promoter.classify_mode(reg_seq, top.start, top.end)
        mode = {"mode": mc.mode, "spacer_len": mc.spacer_len, "confidence": mc.confidence,
                "evidence": mc.evidence, "tss": mc.promoter.get("tss"), "flags": mc.flags,
                "occludes": mc.occludes, "overlaps": mc.overlaps, "footprint": mc.footprint}
    except Exception as e:
        if verbose:
            print(f"  mode classification skipped ({type(e).__name__})")

    # ---- 8. best DNA (design candidate) ------------------------------------------------------------
    designed = _designed_operator(glog)
    design_method = "conservation consensus (genomic operator logo)"
    design_note = ("conservation consensus from the genomic operator logo. No in-silico affinity model is "
                   "applied, so treat this as a conservation-based starting point and confirm it "
                   "experimentally.")

    if verbose:
        print(f"conservation: {len(homs)} homolog region(s); {len(sf_seqs)} motif seed(s); "
              f"top kind={top.kind} half={top.half} spacer={top.spacer} cons={top.conservation:.2f}")
        print(f"genome rescan: {len(rescan_hits)} putative site(s) from the w={best.anchor_w} seed; "
              f"nearest to autoregulatory operator = {nearest} bp")
        print(f"binding motif: consensus={glog.consensus} meanIC={mean_ic:.2f} (kept {sel.n_keep}/{len(ranked)})")
        print(f"AUTOREGULATORY OPERATOR at own promoter: {'YES' if autoreg_recovered else 'no'} "
              f"({nearest} bp from self-promoter anchor)")
        print(f"designed operator (best DNA): {designed}  [{design_method}]")

    # ---- 9. assemble + write the job bundle ---------------------------------------------------------
    win = [op_center - 1800, op_center + 1800]
    rescan_payload = {"region_len": len(ctx.sequence), "plant_center": op_center, "n_hits": len(rescan_hits),
                      "hits": [{"accession": h.genome_accession, "start": h.start, "end": h.end,
                                "strand": h.strand, "dyad": h.dyad_center,
                                "score": (float(h.score) if h.score is not None else None),
                                "pvalue": (float(h.pvalue) if h.pvalue is not None else None),
                                "qvalue": (float(h.qvalue) if h.qvalue is not None else None),
                                "qvalue_kind": getattr(h, "qvalue_kind", None),
                                "seq": getattr(h, "seq", "") or "", "generator": h.generator}
                               for h in rescan_hits]}
    genomic_logo = ({"pwm": glog.pwm.tolist(), "per_col_ic": glog.per_col_ic.tolist(),
                     "consensus": glog.consensus, "total_ic": glog.total_ic,
                     "n_kept": sel.n_keep, "n_hits": len(ranked)} if glog.pwm is not None else {})
    # The seed PWM in the same shape as `genomic_logo`, so the report can draw both with one renderer
    # and state each one's support. `best` is the chosen trial; its cluster carries the matrix that
    # actually scanned the genome.
    seed_logo = {}
    try:
        import numpy as _np
        _spwm = getattr(best.cluster, "pwm", None)
        if _spwm is not None:
            _p = _np.asarray(_spwm, dtype=float)
            _ic = (2.0 + (_p * _np.log2(_np.clip(_p, 1e-9, None))).sum(axis=0))
            _idx = "ACGT"
            seed_logo = {
                "pwm": _p.tolist(), "per_col_ic": _ic.tolist(),
                "consensus": "".join(_idx[i] for i in _p.argmax(axis=0)),
                "total_ic": float(_ic.sum()),
                "anchor_w": best.anchor_w, "cluster_id": best.cluster_id,
                "n_input": best.cluster.prepared.n_input,
                "n_effective": best.cluster.prepared.n_effective,
                "n_independent": getattr(best.cluster, "n_independent", None),
            }
    except Exception as e:
        if verbose:
            print(f"  seed logo skipped ({type(e).__name__}: {e})")

    neigh = {"center": op_center, "window": win,
             "genes": [{"name": g.name or g.protein_id, "start": g.start, "end": g.end, "strand": g.strand}
                       for g in ctx.genes if g.end > win[0] and g.start < win[1]],
             "operator": [op_seq_start, op_seq_end]}

    # AF3 hand-off: the top putative operators as TF-dimer + dsDNA jobs (validate the best DNA). Dedupe by
    # LOCUS first: a near-palindromic operator (e.g. NmtR) is independently detected on BOTH strands at
    # (near-)identical coordinates -- the same physical dsDNA duplex -- so without this the same site could
    # consume 2 of the 5 AF3 job slots (and produce duplicate job names). Keep the best-scoring hit per
    # locus bin (same round(dyad/10) tolerance signals.regulon uses for site dedup), then take the top 5
    # distinct loci.
    ion = af3_export.map_effector_to_ion(effector)
    best_by_locus: dict = {}
    for h in rescan_hits:
        key = (h.genome_accession, round((h.dyad_center or 0) / 10))
        if key not in best_by_locus or (h.score or 0) > (best_by_locus[key].score or 0):
            best_by_locus[key] = h
    picks = sorted(best_by_locus.values(), key=lambda h: -(h.score or 0))[:AF3_TOP_LOCI] or []
    af3_jobs = (af3_export.to_af3_jobs(picks, seq, genome_seq=ctx.sequence, length=AF3_OPERATOR_LENGTH, ion=ion,
                                       name_prefix=name)
                if (picks and getattr(cfg, "emit_af3_jobs", True)) else [])
    try:
        af3_export.validate_jobs(af3_jobs)
    except Exception:
        af3_jobs = []

    # Metal context: the survey's second axis and the anchor support behind the cluster label. Both
    # are already decided by this point and were simply not being written anywhere a reader looks.
    metal_context = {}
    try:
        from predictor.annotate import neighborhood as _nbmod
        if inducer_cons is not None:
            metal_context = dict(_nbmod.survey_effective(inducer_cons))
            _ssnc = next((c for c in inducer_cons.calls if c.source == "ssn_cluster"), None)
            _anch = ((_ssnc.evidence or {}).get("anchors") if _ssnc else None) or {}
            if _anch:
                metal_context["cluster_label_support"] = {
                    "n_anchors": _anch.get("n"), "n_independent": _anch.get("n_independent"),
                    "tier": _anch.get("tier"), "conflict": _anch.get("conflict"),
                    "consensus": _anch.get("consensus"),
                }
            if _ssnc is not None and len(_ssnc.candidates) > 1:
                metal_context["cluster_assigned_to"] = list(_ssnc.candidates)
    except Exception as e:
        if verbose:
            print(f"  metal context skipped ({type(e).__name__}: {e})")

    doss = TFDossier(
        tf_id=name, family=family, genome_accession=ctx.accession,
        tf_locus=[tf_start, tf_end, tf_strand], effector=effector, ion=ion,
        operator_window=ctx.sequence[op_seq_start:op_seq_end].upper(),
        rescan=rescan_payload, neighborhood=neigh, genomic_logo=genomic_logo, af3_jobs=af3_jobs,
        metal_context=metal_context, seed_logo=seed_logo,
        ssn_cluster=ssn_cl,
        inducers=(asdict(inducer_cons) if inducer_cons is not None else {}),
        ligand=ligand, ligand_genes=lgenes, ligand_keywords=lig_kws, regulon=regulon,
        regulon_stats=regulon_stats,
        caveats=[
            "Novel-TF cold start: conservation-first operator discovery (homolog phylogenetic "
            "footprinting + palindrome priors). An apo structure may be generated for QC/inducer "
            "evidence when a configured backend is reachable, but it does not score operator sites.",
            "No in-silico structural energy/affinity model is applied: discovery rests on conservation "
            "alone, and any structural arm added later is exploratory.",
            "Operators are found by a genome-wide rescan of the TF's own motif; the regulon lists the "
            "operons that PWM drives. Candidate regulated genes are reported separately and never used to "
            "rank operators, so an operator at a handler gene stays independent evidence. Regulated-gene "
            "flagging depends on genome annotation.",
        ])
    # attach extra design + recovery fields the dossier schema doesn't have a slot for
    doss_d = asdict(doss)
    doss_d["family_support"] = support_level
    doss_d["genome_topology"] = {
        "circular": bool(getattr(ctx, "circular", False)),
        "origin_spanning_sites_scanned": False if getattr(ctx, "circular", False) else None,
        "origin_wrapped_features": list(getattr(ctx, "origin_wrapped_features", ()) or ()),
    }
    if getattr(ctx, "circular", False):
        doss_d["caveats"].append(
            "Circular replicon: the origin-spanning intergenic interval is deliberately omitted because "
            "wrapped operators cannot be represented by the current linear coordinate schema."
        )
    doss_d["designed_operator"] = {"sequence": designed, "method": design_method, "note": design_note}
    # autoregulation as a FINDING: is a strong operator present at the TF's OWN promoter? ("no" = weak/non-
    # autoregulator, e.g. CueR, not a discovery failure).
    doss_d["autoregulatory_recovery"] = {"recovered": autoreg_recovered, "nearest_bp": nearest}
    doss_d["promoter_occlusion"] = mode
    doss_d["conservation"] = {"n_homologs": len(homs), "n_seeds": len(sf_seqs),
                              "top_conservation": round(float(top.conservation), 2),
                              "msa_expansion": msa_info}
    # every seed width that was tried, its genome-wide outcome, and which one the data chose -- so a
    # reader can see whether the motif rested on a contested call or an uncontested one
    doss_d["seed_selection"] = seed_trial_report
    doss_d["multi_source_logos"] = multi_logos                   # WS3: tf_only / homolog / gathered
    doss_d["per_hit_regulation"] = per_hit_reg                   # WS4: occlusion + TSS/CDS per genome hit
    doss_d["confirmed_regulated_genes"] = confirmed_genes        # compatibility key; includes support labels
    struct_d = struct.as_dict() if struct is not None else None  # WS1: apo dimer fold + AFDB QC
    doss_d["structure"] = struct_d
    # report exports: the TF record (sequence/organism/nearest neighbour) + the homolog promoter regions
    tf_record_payload = asdict(rec) if rec is not None else {}
    tf_record_payload.update({"sequence": seq, "family": family, "organism": organism})
    doss_d["tf_record"] = tf_record_payload
    # `shortfall` says the conservation set came in under MIN_MAPPABLE_ORTHOLOGS, and why. It is
    # recorded by msa_homologs.collect_msa_regions and reaches run_record.json through
    # homolog_collection, but the dossier is what the report, the figures and the analysis scorers
    # read -- a thin conservation axis that appears only in the provenance file is not visible where
    # anyone looks at it. Absent when the floor was met, or when the MSA path never ran.
    doss_d["homolog_regions"] = {"pre_msa": homs_pre, "expanded": homs}
    _msa_shortfall = ((msa_info.get("source_record") or {}) or {}).get("shortfall")
    if _msa_shortfall:
        doss_d["homolog_regions"]["shortfall"] = _msa_shortfall
    doss_d["literature_seed_candidates"] = [
        {
            "sequence": c.seq, "start": reg_start + c.start, "end": reg_start + c.end,
            "score": round(float(c.score), 4), "source": "curated-family-pwm",
            "independent_validation": False,
        }
        for c in literature_candidates
    ]

    # ---- WS5: multi-source operator candidates (design AFTER structure) -----------------------------
    operators = None
    try:
        from predictor.report import operators as _ops
        operators = _ops.build_operators(doss_d, struct=struct_d, family=family)
        doss_d["operators"] = operators
        if verbose and operators.get("primary"):
            p = operators["primary"]
            print(f"operators: {len(operators['candidates'])} candidate(s) across sources; "
                  f"primary={p['source']} ({p['sequence']})")
    except Exception as e:
        if verbose:
            print(f"  WS5 operators table skipped ({type(e).__name__}: {e})")
        return _failure(name, family, "operators", f"{type(e).__name__}: {e}", genome=ctx.accession)
    if not operators or not operators.get("primary"):
        return _failure(name, family, "operators", "no primary operator could be assembled",
                        genome=ctx.accession)

    bundle = None
    if write:
        # WS6 Stage-A: the self-contained per-TF folder (structure/ operators/ ligand/ genome/ motif/ regulation/)
        try:
            source_context = {
                "run_id": run_id,
                "started_at_utc": run_started,
                "requested": requested,
                "effective_config": cfg.as_dict() if hasattr(cfg, "as_dict") else {},
                "query": {"sequence": seq, "input_kind": input_kind},
                "genome": {
                    "accession": ctx.accession,
                    "sequence": ctx.sequence,
                    "genes": [asdict(g) for g in ctx.genes],
                    "circular": bool(getattr(ctx, "circular", False)),
                    "origin_wrapped_features": list(getattr(ctx, "origin_wrapped_features", ()) or ()),
                    "contig_offsets": dict(getattr(ctx, "contig_offsets", {}) or {}),
                },
                "windows": {
                    "discovery_window": {
                        "start": reg_start,
                        "end": reg_start + len(reg_seq),
                        "strand": tf_strand,
                        "sequence": reg_seq,
                    },
                    "top_candidate_window": {
                        "start": op_seq_start,
                        "end": op_seq_end,
                        "strand": tf_strand,
                        "sequence": ctx.sequence[op_seq_start:op_seq_end].upper(),
                    },
                },
                "homologs": {
                    "pre_msa": homs_pre,
                    "expanded": homs,
                },
                "homolog_collection": homolog_collection,
            }
            bundle = str(outputs.write_stage_a(
                doss_d, struct=struct_d, operators=operators, source_context=source_context
            ))
        except Exception as exc:
            return _failure(name, family, "output", f"{type(exc).__name__}: {exc}", genome=ctx.accession)

    return {"tf_id": name, "family": family, "family_support": support_level,
            "path": "novel", "route": "conservation-first", "ok": True,
            "required_stages_complete": True,
            "genome": ctx.accession, "tf_locus": [tf_start, tf_end, tf_strand],
            "n_putative_sites": len(rescan_hits), "motif_consensus": glog.consensus,
            "motif_mean_ic": round(float(mean_ic), 2), "n_homologs": len(homs),
            "msa_expansion": msa_info,
            "autoregulatory_recovered": autoreg_recovered, "nearest_bp": nearest,
            "designed_operator": designed, "design_method": design_method, "design_note": design_note,
            "mode": mode.get("mode"), "operator_window": doss_d["operator_window"],
            "inducer": (inducer_cons.top if inducer_cons is not None else None),
            "ssn_cluster": ssn_cl, "ligand": (ligand[0]["ligand"] if ligand else None),
            "n_regulated_gene_promoters": len(lgenes),
            "n_regulon_operons": len([o for o in regulon if isinstance(o, dict) and "first_gene" in o]),
            "structure_folded": (struct.folded if struct is not None else False),
            "structure_plddt": (struct.plddt_mean if struct is not None else None),
            "structure_qc_rmsd": (struct.qc_rmsd if struct is not None else None),
            "operators_primary": (operators.get("primary") if operators else None),
            "bundle": bundle}


# --------------------------------------------------------------------------- public + benchmark seams
def run_novel(seq: str, *, name: str = "query", family: str | None = None,
              organism: str | None = None, genome_acc: str | None = None, tf_locus=None,
              effector: str | None = None, cfg=None) -> dict:
    """Run the single production path with the frozen production seed construction."""
    return _run_novel_impl(
        seq,
        name=name,
        family=family,
        organism=organism,
        genome_acc=genome_acc,
        tf_locus=tf_locus,
        effector=effector,
        cfg=cfg,
        seed_construction=ss.DEFAULT_CONSTRUCTION,
    )


def _run_seed_construction_benchmark(
    seq: str,
    *,
    construction: ss.SeedConstruction,
    name: str,
    family: str,
    organism: str,
    cfg=None,
) -> dict:
    """Internal Phase-B seam; it changes construction inputs, never the pipeline path."""
    if construction.source_run_id is None:
        raise ValueError("benchmark construction requires complete immutable source provenance")
    return _run_novel_impl(
        seq,
        name=name,
        family=family,
        organism=organism,
        cfg=cfg,
        seed_construction=construction,
    )


# --------------------------------------------------------------------------- dispatcher
def predict(tf: str, *, name: str | None = None, **kw) -> dict:
    """The single entry. `tf` = a TF protein sequence or a FASTA path -> the production novel-TF path."""
    return run_novel(tf, name=name or "query", **kw)


def _print_summary(s: dict) -> None:
    print("\n--- summary ---")
    for k in ("tf_id", "family", "path", "route", "genome", "n_putative_sites", "motif_consensus",
              "motif_mean_ic", "autoregulatory_recovered", "nearest_bp", "designed_operator",
              "structure_folded", "structure_plddt", "structure_qc_rmsd", "mode", "bundle"):
        if k in s and s[k] is not None:
            print(f"  {k:18}: {s[k]}")


# --------------------------------------------------------------------------- self-test (offline)
def _self_test() -> None:
    # 1) sequence reader: a raw sequence (possibly multiline) is whitespace-joined; a FASTA *file* has its
    #    '>' header stripped
    _test_protein = "MSTNPKPQRKTKRNTNRRPQDVKFPGG"
    assert _read_seq(_test_protein) == _test_protein
    assert _read_seq(f"  {_test_protein[:12]}\n{_test_protein[12:]}  ") == _test_protein
    import tempfile as _tmp
    _fh = _tmp.NamedTemporaryFile("w", suffix=".fasta", delete=False)
    _fh.write(f">tf\n{_test_protein[:12]}\n{_test_protein[12:]}\n"); _fh.close()
    try:
        assert _read_seq(_fh.name) == _test_protein
    finally:
        os.unlink(_fh.name)
    # 2) designed operator = the argmax base at every column of the genomic logo
    class _L:
        pwm = np.array([[1.0, 0, 0, 0], [0, 1.0, 0, 0], [0, 0, 1.0, 0], [0, 0, 0, 1.0]])
    assert _designed_operator(_L()) == "ACGT"
    print("OK: pipeline production helpers self-test passed (offline).")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seq", help="a TF protein sequence or FASTA path (novel-TF path)")
    ap.add_argument("--name", default="query")
    ap.add_argument("--family", default=None, help="TF family (else auto-classified)")
    ap.add_argument("--organism", default=None, help="RefSeq/GCF accession of the source genome")
    ap.add_argument("--genome-acc", default=None)
    ap.add_argument("--tf-locus", default=None, help="start,end,strand to override tblastn locus resolution")
    ap.add_argument("--effector", default=None)
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args(argv)

    if a.self_test:
        _self_test(); return
    if a.seq:
        locus = tuple(x.strip() for x in a.tf_locus.split(",")) if a.tf_locus else None
        if locus:
            locus = (int(locus[0]), int(locus[1]), locus[2])
        s = run_novel(a.seq, name=a.name, family=a.family, organism=a.organism, genome_acc=a.genome_acc,
                      tf_locus=locus, effector=a.effector)
        _print_summary(s); return
    ap.print_help()


if __name__ == "__main__":
    main()
