"""AF3 export policies must not invent DNA or mutate existing scientific bundles."""
import csv
import hashlib
import json
import os
import shutil
from pathlib import Path

import pytest

from predictor.structure import af3_export as legacy
from predictor.structure.af3_jobset import (
    build_job_set, export_bundle, export_run, to_local_json, to_server_json, validate_local,
)

PROTEIN = "MACKDHILNPQRSTVWYEFG"
DNA = "ACGT" * 15


def op(dna=DNA, start=100, accession="NC_TEST"):
    return dict(dna_top=dna, start=start, end=start + 20, dyad_center=start + 10,
                genome_accession=accession)


def test_four_kinds_dual_dialects_and_fe_neutrality():
    result = build_job_set(PROTEIN, [op(), op("A" * 60, 200), op("C" * 60, 300)], "Fe2+", "Fur")
    assert [j.kind for j in result.jobs] == ["apo", "dna_r1", "dna_r2", "dna_r3", "metal", "metal_dna_r1"]
    for job in result.jobs:
        legacy.validate_jobs([to_server_json(job)])
        validate_local(to_local_json(job))
    metal = result.jobs[-2]
    assert to_server_json(metal)["sequences"][-1]["ion"] == {"ion": "MN", "count": 2}
    assert to_local_json(metal)["sequences"][-1]["ligand"]["ccdCodes"] == ["FE2"]
    assert legacy.map_effector_to_ion("Fe2+") == "FE"  # legacy dossier policy unchanged
    assert "no structural ions inferred" in metal.assumptions


def test_missing_dna_keeps_apo_metal_and_no_fabrication():
    result = build_job_set(PROTEIN, [{"seq": "ACGT", "start": 2, "end": 6}], "As(III)", "ArsR/SmtB")
    assert [j.kind for j in result.jobs] == ["apo", "metal"]
    assert to_local_json(result.jobs[-1])["sequences"][-1]["ligand"]["ccdCodes"] == ["AST"]
    assert "free AsO3" in result.jobs[-1].ions[0][1].why
    assert any("60 bp" in r["reason"] for r in result.diagnostics)


def test_duplex_locus_dedup_preserves_rank_and_tetramer():
    seq = "A" * 20 + "C" * 40
    ops = [op(seq), op(legacy.revcomp(seq), 200), op("C" * 60, 100), op("G" * 20 + "A" * 40, 400)]
    result = build_job_set(PROTEIN, ops, "Ni2+", "NikR")
    dnas = [j for j in result.jobs if j.kind.startswith("dna_")]
    assert [j.operator_rank for j in dnas] == [1, 4]
    assert result.jobs[0].copies == 4
    assert result.jobs[-1].ions[0][2] == 4
    assert build_job_set(PROTEIN, [], None, "CsoR/FrmR").jobs[0].copies == 4


def test_unresolved_and_polyspecific_metals_are_not_invented():
    assert [j.kind for j in build_job_set(PROTEIN, [], "divalent metal (ion unresolved)", "ArsR").jobs] == ["apo"]
    result = build_job_set(PROTEIN, [], "Cd2+/Pb2+", "ArsR")
    metals = [j for j in result.jobs if j.kind == "metal"]
    assert len(metals) == 2 and len({j.name for j in metals}) == 2
    assert {j.ions[0][0] for j in metals} == {"Cd2+", "Pb2+"}
    mixed = build_job_set(PROTEIN, [], "Cu+/persulfide", "CsoR")
    assert [j.kind for j in mixed.jobs] == ["apo"]


def test_molybdate_and_antimonite_abstention():
    molybdate = build_job_set(PROTEIN, [], "molybdate", "ModE").jobs[-1]
    with pytest.raises(ValueError, match="server analogue"):
        to_server_json(molybdate)
    assert to_local_json(molybdate)["sequences"][-1]["ligand"]["ccdCodes"] == ["MOO"]
    sb = build_job_set(PROTEIN, [], "antimonite", "ArsR").jobs[-1]
    with pytest.raises(ValueError, match="local CCD"):
        to_local_json(sb)


def test_local_validation_rejects_duplicate_ids():
    payload = to_local_json(build_job_set(PROTEIN, [op()], None, "ArsR").jobs[-1])
    payload["sequences"][-1]["dna"]["id"] = "A"
    with pytest.raises(ValueError, match="unique"):
        validate_local(payload)


def make_bundle(path, effector="Fe2+"):
    path.mkdir(parents=True)
    (path / "input").mkdir()
    (path / "input" / "query.fasta").write_text(f">TF\n{PROTEIN}\n")
    (path / "dossier.json").write_text(json.dumps(dict(tf_id=path.name, effector=effector, family="Fur")))
    (path / "operators_ranked.json").write_text(json.dumps([op()]))
    (path / "bundle_manifest.json").write_text('{"unchanged": true}')


def snapshot(path):
    return {str(p.relative_to(path)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in path.rglob("*") if p.is_file()}


def test_sidecar_run_rebuild_leaves_bundle_bytes_untouched(tmp_path):
    run = tmp_path / "run"
    bundle = run / "jobs" / "TF"
    make_bundle(bundle)
    before = snapshot(bundle)
    summary = export_run(run, batch_size=2)
    assert summary["server_jobs"] == summary["local_jobs"] == 4
    assert snapshot(bundle) == before
    batches = list((run / "af3" / "server").glob("*.json"))
    assert len(batches) == 2
    for path in batches:
        legacy.validate_jobs(json.loads(path.read_text()))
    with (run / "af3" / "af3_manifest.tsv").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert all((run / "af3" / row["file"]).is_file() for row in rows if row["status"] == "written")
    with pytest.raises(ValueError, match="outside"):
        export_bundle(bundle, bundle / "af3")


def test_legacy_protein_fallback_and_missing_input_diagnostics(tmp_path):
    bundle = tmp_path / "old"
    make_bundle(bundle, "molybdate")
    (bundle / "input" / "query.fasta").unlink()
    (bundle / "af3_jobs.json").write_text(json.dumps([legacy.af3_job("old", PROTEIN, DNA)]))
    result = export_bundle(bundle)
    assert result["server_jobs"] == 2 and result["local_jobs"] == 4
    assert any("server analogue" in w for w in result["warnings"])
    (bundle / "af3_jobs.json").unlink()
    result = export_bundle(bundle, tmp_path / "no_protein")
    assert result["server_jobs"] == result["local_jobs"] == 0
    assert any("protein sequence" in w for w in result["warnings"])


def test_contig_mismatch_cannot_slice_wrong_genome():
    result = build_job_set(PROTEIN, [{"start": 40, "end": 60, "genome_accession": "OTHER"}],
                           None, "ArsR", "ACGT" * 100, genome_accession="NC_TEST")
    assert [j.kind for j in result.jobs] == ["apo"]
    missing = build_job_set(PROTEIN, [{"sequence": "ACGT"}], None, "ArsR", "ACGT" * 100)
    assert [j.kind for j in missing.jobs] == ["apo"]
    explicit = build_job_set(PROTEIN, [{"dna_top": "A" * 60}, {"dna_top": "C" * 60}], None, "ArsR")
    assert [j.kind for j in explicit.jobs] == ["apo", "dna_r1", "dna_r2"]


def test_canonical_ranked_object_bundle_format(tmp_path):
    bundle = tmp_path / "canonical"
    make_bundle(bundle)
    (bundle / "genome").mkdir()
    genome = "ACGT" * 200
    (bundle / "genome" / "genome.fna").write_text(f">NC_TEST\n{genome}\n")
    ranked = {"ranked": [{"sequence": "ACGT", "source": "conservation_natural",
                          "provenance": {"accession": "NC_TEST", "start": 98, "end": 102, "dyad": 100}}]}
    (bundle / "operators_ranked.json").write_text(json.dumps(ranked))
    result = export_bundle(bundle)
    assert result["server_jobs"] == result["local_jobs"] == 4
    jobs = json.loads(next((tmp_path / "canonical_af3" / "server").glob("*.json")).read_text())
    dna_job = next(j for j in jobs if j["name"].endswith("__dna_r1"))
    assert dna_job["sequences"][1]["dnaSequence"]["sequence"] == genome[70:130]


@pytest.mark.skipif(not os.environ.get("RIIO_CANONICAL_RUN"), reason="optional canonical run integration fixture")
def test_real_canonical_export_readonly(tmp_path):
    run = Path(os.environ["RIIO_CANONICAL_RUN"])
    source = next((run / "jobs").glob("*/dossier.json")).parent
    target = tmp_path / source.name
    target.mkdir()
    for relative in ("dossier.json", "operators_ranked.json", "af3_jobs.json", "input/query.fasta", "genome/genome.fna"):
        original = source / relative
        if original.is_file():
            copied = target / relative
            copied.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(original, copied)
    before = snapshot(target)
    summary = export_bundle(target)
    assert snapshot(target) == before
    assert summary["server_jobs"] >= 2 and summary["local_jobs"] >= 2
    rows = list(csv.DictReader(Path(summary["manifest"]).open(encoding="utf-8"), delimiter="\t"))
    assert any(r["kind"].startswith("dna_r") and r["status"] == "written" for r in rows)
