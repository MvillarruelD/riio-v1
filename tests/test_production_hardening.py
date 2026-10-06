"""Production-boundary and shareable-bundle safety contracts."""
from __future__ import annotations

import json
import re
from pathlib import Path
from types import SimpleNamespace

import pytest

from predictor import pipeline
import predictor
from predictor.api import Result
from predictor.annotate import context
from predictor.input_validation import family_support, read_protein, validate_family
from predictor.report import outputs
from predictor.signals.motif_rescan import GenomeHit, _genome_wide_qvalues, intergenic_intervals


PROTEIN = "MSTNPKPQRKTKRNTNRRPQDVKFPGG"


def test_single_protein_input_accepts_raw_and_one_fasta(tmp_path):
    assert read_protein(PROTEIN.lower()) == PROTEIN
    fasta = tmp_path / "tf.faa"
    fasta.write_text(f">tf\n{PROTEIN[:12]}\n{PROTEIN[12:]}*\n", encoding="utf-8")
    assert read_protein(str(fasta)) == PROTEIN


@pytest.mark.parametrize(
    "value, message",
    [
        (">a\nMSTNPKPQRKTKRNTNRRPQ\n>b\nMPEPTIDESEQUENCEAAAAA\n", "exactly one FASTA record"),
        ("ATGCGTACGTACGTACGTACGTACGTACGT", "appears to be nucleotide"),
        ("MSTNPKPQRK!KRNTNRRPQ", "invalid residue"),
        ("MPEPTIDE", "implausibly short"),
    ],
)
def test_invalid_public_inputs_abstain(value, message):
    with pytest.raises(ValueError, match=message):
        read_protein(value)


def test_manual_family_must_be_canonical():
    assert validate_family("MerR") == "MerR"
    with pytest.raises(ValueError, match="unknown TF family"):
        validate_family("MerRR-typo")
    assert family_support("MerR") == "family_validated"
    assert family_support("Fur") == "partially_validated"
    assert family_support("TetR/AcrR") == "ssn_supported_unvalidated"
    assert family_support("AraC/XylS") == "generic_unvalidated"


def test_unclassified_query_returns_a_required_stage_failure(monkeypatch):
    monkeypatch.setattr(pipeline.tf_record, "classify", lambda _seq: SimpleNamespace(family=None, flags=["none"]))
    got = pipeline.run_novel(PROTEIN, cfg=SimpleNamespace(allow_online=False, allow_colabfold=False,
                                                         allow_folds=False, allow_metalnet=False,
                                                         emit_af3_jobs=False, write=False, verbose=False))
    assert got["failed_stage"] == "family"
    assert not got["required_stages_complete"]
    assert not Result.from_run(got).ok


def test_job_names_are_contained(tmp_path):
    assert outputs.job_path("safe run", root=tmp_path).parent == tmp_path.resolve()
    for name in ("../escape", "a/b", str(tmp_path / "absolute"), "bad:name", "CON", "trailing."):
        with pytest.raises(ValueError):
            outputs.job_path(name, root=tmp_path)


def _minimal_bundle(tf_id="TF1"):
    dossier = {
        "tf_id": tf_id,
        "family": "MerR",
        "rescan": {"hits": []},
        "regulon": [],
        "af3_jobs": [],
        "tf_record": {"sequence": PROTEIN},
    }
    primary = {"source": "conservation_logo", "sequence": "ACGTACGT", "score": 1.0,
               "score_type": "mean_ic_bits", "provenance": {}}
    operators = {"candidates": [], "ranked": [primary], "primary": primary}
    return dossier, operators


def _synthetic_runtime():
    return {
        "package_version": "test",
        "data_release": "synthetic",
        "source_sha256": "0" * 64,
        "manifest_sha256": "1" * 64,
        "runtime_sha256": "2" * 64,
    }


def test_stage_a_is_atomic_and_preserves_an_incomplete_predecessor(tmp_path, monkeypatch):
    monkeypatch.setattr(outputs, "runtime_provenance", _synthetic_runtime)
    incomplete = tmp_path / "TF1"
    incomplete.mkdir()
    (incomplete / "old.partial").write_text("recover me", encoding="utf-8")
    dossier, operators = _minimal_bundle()
    final = outputs.write_stage_a(dossier, operators=operators, root=tmp_path, render=False)
    assert final == incomplete
    assert not (final / "old.partial").exists()
    manifest = json.loads((final / "bundle_manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "complete"
    assert manifest["runtime"] == _synthetic_runtime()
    assert re.fullmatch(r"[0-9a-f]{64}", manifest["runtime"]["runtime_sha256"])
    backups = list(tmp_path.glob(".TF1.incomplete-*"))
    assert len(backups) == 1 and (backups[0] / "old.partial").read_text() == "recover me"


def test_stage_a_refuses_to_mix_with_a_completed_bundle(tmp_path, monkeypatch):
    monkeypatch.setattr(outputs, "runtime_provenance", _synthetic_runtime)
    dossier, operators = _minimal_bundle()
    outputs.write_stage_a(dossier, operators=operators, root=tmp_path, render=False)
    with pytest.raises(FileExistsError, match="completed job already exists"):
        outputs.write_stage_a(dossier, operators=operators, root=tmp_path, render=False)


def test_circular_gff_is_recorded_and_terminal_gap_is_not_linearized(tmp_path):
    fasta = tmp_path / "c.fna"
    fasta.write_text(">c1\n" + "A" * 1000 + "\n", encoding="utf-8")
    gff = tmp_path / "c.gff"
    gff.write_text(
        "##gff-version 3\n"
        "c1\t.\tregion\t1\t1000\t.\t+\t.\tID=r;Is_circular=true\n"
        "c1\t.\tCDS\t100\t200\t.\t+\t0\tID=a;gene=a\n"
        "c1\t.\tCDS\t700\t800\t.\t+\t0\tID=b;gene=b\n",
        encoding="utf-8",
    )
    ctx = context.from_gff(gff, fasta)
    assert ctx.circular
    assert intergenic_intervals(ctx.genes, len(ctx.sequence), circular=True) == [(200, 699)]


def test_filtered_bh_values_are_explicitly_upper_bounds():
    hit = GenomeHit("A", 1, 11, "+", 6, 10.0, 1e-6, 1e-6, 1, 1.0, 6.0, "ACGTACGTAC")
    adjusted = _genome_wide_qvalues([hit], 1_000_000)[0]
    assert adjusted.qvalue == 1.0
    assert adjusted.qvalue_kind == "conservative_bh_upper_bound"


def test_result_requires_mandatory_stage_completion():
    assert not Result.from_run({"tf_id": "x", "required_stages_complete": False}).ok
    assert Result.from_run({"tf_id": "x", "required_stages_complete": True}).ok


def test_package_versions_stay_synchronized():
    pyproject = (Path(predictor.__file__).resolve().parents[1] / "pyproject.toml").read_text(encoding="utf-8")
    match = re.search(r'^version\s*=\s*"([^"]+)"', pyproject, flags=re.MULTILINE)
    assert match and match.group(1) == predictor.__version__


class TestPipelineStagesActuallyRun:
    """A stage that silently returns None is worse than one that raises.

    `pipeline._infer_inducer` once lost its entire body to an over-eager edit: it kept its signature,
    caught nothing, and returned `(None, None, effector, None)` for every protein. The full test
    suite passed and `tfop selftest` passed, because neither runs the novel-TF pipeline -- the
    failure only appeared as "inducer consensus was not produced" on 45 of 45 candidates in a batch.

    These are cheap structural assertions that would have caught it. They are not a substitute for an
    end-to-end run, but they fail in seconds rather than after a genome fetch.
    """

    def test_infer_inducer_calls_the_inducer(self):
        import inspect
        from predictor import pipeline
        src = inspect.getsource(pipeline._infer_inducer)
        assert "infer_inducer(" in src, "the inducer stage does not call infer_inducer"
        assert "return" in src

    def test_infer_inducer_returns_a_consensus_for_a_real_sequence(self):
        from types import SimpleNamespace as NS
        from predictor import pipeline
        # SmtB/ArsR-like; any real sequence will do -- the point is that a consensus comes back
        seq = ("MTDNNQALKDAGLKVTLPRLKILEVLQEPDNHHVSAEDLYKRLIDMGEEIGLATVYRVLNQFDDAGIVTRHNFEGGKS"
               "VFELTQQHHHDHLICLDCGKVIEFSDDSIEARQREIAAKHGIRLTNHSLYLYGHCAEGDCREDEHAHEGK")
        ctx = NS(accession="SYN.1", sequence="N" * 3000,
                 genes=[NS(start=100, end=550, strand="+", name="tf", protein_id="WP_X",
                           product="transcriptional regulator")])
        cons, _ssn, _eff, failure = pipeline._infer_inducer(
            seq, "Fur", ctx, 100, None, allow_ncbi=False, verbose=False, allow_metalnet=False)
        assert cons is not None, f"no inducer consensus was produced (failure={failure!r})"
        assert cons.top, "the consensus has no headline"
        assert any(c.source == "coordination" for c in cons.calls), \
            "the coordination gate did not contribute -- the stage is not wired to the sources"


# --------------------------------------------------------------------------- replicon concatenation
def test_replicon_pad_exceeds_every_span_that_could_cross_it():
    """A multi-replicon assembly is concatenated into ONE coordinate space, separated by
    `REPLICON_PAD` N's. That pad is the only thing stopping the concatenation from inventing
    relationships across a boundary, and it works solely because it is larger than every span the
    pipeline reaches across. Raising any of these limits past the pad silently re-opens boundary
    leakage, in a way no output would look wrong -- so assert the relation here, where a change to
    any one of four modules trips it.
    """
    from predictor.annotate.context import REPLICON_PAD
    from predictor.signals.operator_logo import OPERATOR_WIDTHS
    from predictor.signals.regulon import DEFAULT_OPERON_MAX_GAP

    assert REPLICON_PAD > max(OPERATOR_WIDTHS), "a motif could be called across a replicon boundary"
    assert REPLICON_PAD > DEFAULT_OPERON_MAX_GAP, "an operon could be walked across a boundary"

    import inspect
    from predictor.report import ligand_regulon
    from predictor.signals import regulon as regmod
    for fn in (ligand_regulon._regulated_gene, regmod.regulated_first_genes):
        md = inspect.signature(fn).parameters["max_dist"].default
        assert REPLICON_PAD > md, (
            f"{fn.__qualname__} reaches {md} bp, which crosses a {REPLICON_PAD} bp pad: an operator "
            f"could be assigned a gene on another replicon")


def test_wrap_guard_is_scoped_to_the_feature_s_own_replicon():
    """The origin-wrap guard drops a merged feature spanning >80 % of its replicon. Measured against
    the whole concatenation instead, a feature wrapping the origin of a small plasmid spans a
    fraction of a percent and sails through the test written to catch it."""
    from predictor.annotate.context import from_gff, REPLICON_PAD

    big, small = 20000, 2000
    fasta = f">chr\n{'ACGT' * (big // 4)}\n>plasmid\n{'ACGT' * (small // 4)}\n"
    # a plasmid feature merged across its own origin: one ID, two parts, near each end
    gff = "\n".join([
        "chr\t.\tCDS\t100\t400\t.\t+\t0\tID=g1;gene=normal",
        "plasmid\t.\tCDS\t50\t200\t.\t+\t0\tID=w1;gene=wrapped",
        f"plasmid\t.\tCDS\t{small - 150}\t{small - 10}\t.\t+\t0\tID=w1;gene=wrapped",
    ])
    ctx = from_gff(gff, fasta)
    assert "wrapped" in ctx.origin_wrapped_features, (
        "a feature spanning ~95 % of its 2 kb plasmid must be dropped even though it is only ~9 % "
        "of the 22 kb concatenation")
    assert {g.name for g in ctx.genes} == {"normal"}
    assert ctx.contig_offsets["plasmid"] == big + REPLICON_PAD


def test_locate_resolves_a_concatenated_coordinate_to_its_replicon():
    from predictor.annotate.context import from_gff, REPLICON_PAD

    big, small = 20000, 2000
    fasta = f">chr\n{'ACGT' * (big // 4)}\n>plasmid\n{'ACGT' * (small // 4)}\n"
    one_line_gff = "chr\t.\tCDS\t100\t400\t.\t+\t0\tID=g1;gene=normal\n"
    ctx = from_gff(one_line_gff, fasta)
    assert ctx.locate(10) == ("chr", 10)
    assert ctx.locate(big + REPLICON_PAD + 7) == ("plasmid", 7)
    # a single-record genome needs no special case at any call site
    solo = from_gff(one_line_gff, f">chr\n{'ACGT' * 100}\n")
    assert solo.locate(42) == (solo.accession, 42)


def test_no_environment_variable_selects_an_algorithm():
    """`config.RunConfig` states the single-path rule -- "every lever is about COST or REACHABILITY
    ... none of them selects an algorithm, a family branch or a scoring rule" -- but asserts it only
    over its own fields. A dataclass default elsewhere reading os.environ slips straight past that,
    which is how TFOP_SCORE_MODE came to pick the dyad scoring rule for a whole run, bound once at
    import, with the shell as the only record of it.

    Env vars for COST and REACHABILITY are fine (where a binary lives, an API key, a search
    threshold handed to an external tool). What is not fine is one that changes which algorithm the
    predictor runs.
    """
    import re
    from pathlib import Path

    pkg = Path(__file__).resolve().parents[1] / "predictor"
    #: The baseline, every one of them reachability or cost: where a binary lives, which credential
    #: to present, where a cache goes, what an EXTERNAL tool is allowed to spend. None of them
    #: changes which algorithm this package runs. Listed one by one on purpose -- a new name has to
    #: be classified by a person, which is the whole point of the guard.
    ALLOWED = {
        # where things are
        "PREDICTOR_DATA_DIR", "PREDICTOR_REFS_DIR", "PREDICTOR_CACHE_DIR", "PREDICTOR_OUTPUT_DIR",
        "PREDICTOR_ENV_FILE", "PREDICTOR_ESMFOLD_LIB",
        "PREDICTOR_GENOME_DB", "PREDICTOR_GENOME_MIRROR", "PREDICTOR_PROTSEQ_CACHE",
        "PREDICTOR_EFFECTOR_CACHE", "PREDICTOR_STRUCT_CACHE",
        "BITACORA_DIR", "BITACORA_BIN", "BITACORA_PY", "BITACORA_PROFILES", "BITACORA_WSL_DISTRO",
            "BLAST_BIN", "MMSEQS", "DATASETS_EXE", "DEEPPBS_PYTHON", "METALNET_PYTHON", "METALNET2_DIR",
        # who we say we are
        "BIOHUB_TOKEN", "ESM_FORGE_TOKEN", "FORGE_TOKEN", "NCBI_TOKEN", "NCBI_EMAIL",
        "NCBI_API_KEY",
        # what an EXTERNAL search may spend -- BITACORA's own thresholds, not ours
        "BITACORA_EVALUE", "BITACORA_QCOV", "BITACORA_SCOV",
        # which remote/local ESM weights to call. The greyest two on this list: they name an
        # external resource, but a different checkpoint would give different embeddings. They are
        # here because the choice is about what is REACHABLE, not about what we compute from it --
        # if that ever stops being true they belong in an explicit argument like score_mode.
        "ESM_FORGE_MODEL", "ESM_LOCAL_MODEL",
    }
    # catch the idioms, not one spelling: os.environ[...], os.environ.get(...), os.getenv(...),
    # and the same through an alias or a `from os import environ, getenv`.
    pat = re.compile(r"""(?:environ\s*(?:\.get\s*)?[\(\[]|getenv\s*\()\s*["']([A-Z][A-Z0-9_]*)["']""")
    found = {}
    for f in pkg.rglob("*.py"):
        if "__pycache__" in f.parts:
            continue
        for name in pat.findall(f.read_text(encoding="utf-8", errors="replace")):
            found.setdefault(name, set()).add(f.relative_to(pkg).as_posix())

    unexpected = {k: sorted(v) for k, v in found.items()
                  if k not in ALLOWED and not k.startswith(("BITACORA_HMM_", "HTTP", "HTTPS", "NO_"))}
    assert not unexpected, (
        "new environment variable(s) read inside predictor/. If it selects an algorithm, a family "
        "branch or a scoring rule, it breaks the single-path rule and must become an explicit "
        f"argument instead; if it is genuinely cost/reachability, add it to ALLOWED here: {unexpected}")


# --------------------------------------------------------------------------- tfop scan --select
def _cand(pid, family):
    from predictor.annotate.genome_scan import Candidate
    return Candidate(protein_id=pid, family=family, pfam_acc=None, pfam_name=None,
                     evalue=1e-30, bits=100.0, length=120, sequence=PROTEIN)


def _scan_args(tmp_path, **over):
    base = dict(input="genome.faa", families=None, cpus=1, out=str(tmp_path / "scan"),
                predict=True, select=None, limit=None, organism=None, genome_acc=None)
    base.update(over)
    return SimpleNamespace(**base)


def test_scan_select_cannot_route_an_unsupported_family_into_the_pipeline(tmp_path, monkeypatch):
    """`--select` narrows the SUPPORTED list; it never re-widens it to the whole census.

    The census recognises more families than the operator pipeline supports. Rebuilding the todo
    list from every candidate let `--select` hand e.g. an AraC/XylS protein to `api.predict`, which
    fails per candidate, late, after the batch has started.
    """
    from predictor import run_cli
    from predictor.annotate import genome_scan as GS
    from predictor import api

    supported, unsupported = _cand("supported_1", "MerR"), _cand("unsupported_1", "AraC/XylS")
    assert unsupported.family not in GS.PIPELINE_FAMILIES, "fixture must use an unsupported family"

    monkeypatch.setattr(GS, "scan", lambda *_a, **_k: ([supported, unsupported], "proteins", 2))
    monkeypatch.setattr(GS, "write_outputs", lambda *_a, **_k: None)
    predicted = []

    def _fake_predict(_seq, name=None, **_k):
        predicted.append(name)
        return SimpleNamespace(inducer=None, primary_operator=None, ok=True, error=None, bundle=None)

    monkeypatch.setattr(api, "predict", _fake_predict)
    (tmp_path / "scan").mkdir()

    rc = run_cli.cmd_scan(_scan_args(tmp_path, select="supported_1,unsupported_1"))
    assert rc == 0
    assert predicted == ["supported_1"], f"unsupported family reached the pipeline: {predicted}"


def test_scan_select_reports_ids_it_drops(tmp_path, monkeypatch, capsys):
    """A requested id that is dropped must be named, with the reason -- never silently ignored."""
    from predictor import run_cli
    from predictor.annotate import genome_scan as GS
    from predictor import api

    cands = [_cand("supported_1", "MerR"), _cand("unsupported_1", "AraC/XylS")]
    monkeypatch.setattr(GS, "scan", lambda *_a, **_k: (cands, "proteins", 2))
    monkeypatch.setattr(GS, "write_outputs", lambda *_a, **_k: None)
    monkeypatch.setattr(api, "predict", lambda _s, name=None, **_k: SimpleNamespace(
        inducer=None, primary_operator=None, ok=True, error=None, bundle=None))
    (tmp_path / "scan").mkdir()

    run_cli.cmd_scan(_scan_args(tmp_path, select="unsupported_1,never_scanned"))
    out = capsys.readouterr().out
    assert "unsupported_1" in out and "AraC/XylS" in out
    assert "never_scanned" in out and "not found in this scan" in out


# --------------------------------------------------------------------------- offline genome cache
def test_ensure_genome_returns_a_cached_genome_without_download(tmp_path):
    """The contract the offline path depends on: mirror FIRST, download gated separately.

    `ensure_genome` resolves an already-mirrored accession before it ever considers downloading, so
    `allow_download=False` still returns a cached genome. Its caller in `pipeline` used to gate the
    whole call on `allow_ncbi`, which made an offline run refuse a genome sitting in its own cache:
    "could not acquire the explicitly requested genome GCF_000195955.2" for a file already on disk.
    """
    from predictor.annotate import genome_mirror

    acc = "GCF_000195955.2"
    (tmp_path / f"{acc}.fna").write_text(">chr\nACGT\n", encoding="utf-8")
    (tmp_path / f"{acc}.gff").write_text("##gff-version 3\n", encoding="utf-8")

    fna, gff = genome_mirror.ensure_genome(acc, allow_download=False, mirror_dir=tmp_path)
    assert fna is not None and fna.name == f"{acc}.fna"
    assert gff is not None and gff.name == f"{acc}.gff"


def test_ensure_genome_reports_absence_rather_than_downloading_when_offline(tmp_path):
    from predictor.annotate import genome_mirror
    fna, gff = genome_mirror.ensure_genome("GCF_999999999.9", allow_download=False,
                                           mirror_dir=tmp_path)
    assert (fna, gff) == (None, None)


def test_the_pipeline_reads_the_genome_mirror_even_with_no_network(monkeypatch):
    """Offline must still CONSULT the cache; only the fetch is online."""
    from predictor import pipeline as P

    calls = []

    def fake_ensure(acc, *, allow_download=True, **kw):
        calls.append({"acc": acc, "allow_download": allow_download})
        return None, None                      # nothing cached -> pipeline fails later, as it should

    monkeypatch.setattr(P.genome_mirror, "ensure_genome", fake_ensure)
    monkeypatch.setattr(P.tf_record, "classify",
                        lambda _seq: SimpleNamespace(family="MerR", flags=[], uniprot=None))
    P.run_novel(PROTEIN, name="q", family="MerR", organism="GCF_000195955.2",
                cfg=SimpleNamespace(allow_online=False, allow_colabfold=False, allow_folds=False,
                                    allow_metalnet=False, emit_af3_jobs=False, write=False,
                                    verbose=False))
    assert calls, "the genome mirror was never consulted offline"
    assert calls[0]["allow_download"] is False, "offline must not attempt a download"
