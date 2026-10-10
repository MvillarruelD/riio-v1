#!/usr/bin/env python
"""make_af3_figure4_jobs.py -- the AlphaFold 3 hand-off that completes Figure 4.

Figure 4 shows six regulators, each with two reserved panels: the holo dimer with its called ion
bound, and a zoom on the metal-binding site. Folding is deliberately not part of the pipeline
(structure prediction was off in the canonical run), so this writes the jobs the user uploads by
hand. Per regulator, two alphafoldserver jobs:

  * **apo + operator DNA** -- the homodimer on the 60 bp centred on its primary operator, no metal;
  * **holo + metal**       -- the homodimer with its called effector metal, no DNA.

Nothing is recomputed. The six regulators and their primary operators come from
`reports/numbers.json` ("context", written by reports/compute_numbers.py); the protein sequence and
the 60 bp duplex come from each bundle's own `af3_jobs.json`, whose rank-1 entry is the same
operator -- asserted per regulator rather than assumed. The predictor is frozen at v1.0.0-rc, so
this edits nothing there and imports its ion vocabulary and helpers.

    python analysis/make_af3_figure4_jobs.py

Writes `analysis/af3_figure4/`: the job JSON, a manifest TSV and a README stating every substitution.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from canonical_config import load_canonical  # noqa: E402

PROJ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJ / "analysis"))
sys.path.insert(0, str(PROJ / "analysis" / "benchmarking"))
from project_config import PREDICTOR_ROOT  # noqa: E402

sys.path.insert(0, str(PREDICTOR_ROOT))
import genome_scoring as GS  # noqa: E402
from predictor.structure.af3_export import (AF3_IONS, map_effector_to_ion,  # noqa: E402
                                            revcomp, sanitize_job_name)

OUT_DIR = PROJ / "analysis" / "af3_figure4"
#: one effector site per protomer, so a homodimer carries two ions
IONS_PER_DIMER = 2
#: the duplex length the pipeline's own hand-off uses (short duplexes break dimer docking)
DNA_LEN = 60

# Exact historical recipe; shared only to preserve the paper hand-off.
from predictor.structure.af3_legacy_policy import DNA_BOUND_FORM, SURROGATE  # noqa: E402


def af3_ion(called: str) -> tuple[str, str, str]:
    """(server ion, quality, why) for a called effector metal, cross-checked against the predictor."""
    if called in SURROGATE:
        ion, quality, why = SURROGATE[called]
    elif (direct := str(called).strip().upper().replace(" ", "")) in AF3_IONS:
        ion, quality, why = direct, "exact", "the server supports this ion directly"
    else:
        raise SystemExit(f"no surrogate recorded for the called metal {called!r} -- add one to "
                         f"SURROGATE with its coordination rationale")
    if ion not in AF3_IONS:
        raise SystemExit(f"{called!r} -> {ion!r} is not one of the server's ions {sorted(AF3_IONS)}")
    # the predictor has its own map; where it resolves, say so rather than diverging silently
    theirs = map_effector_to_ion(called)
    note = "" if theirs in (None, ion) else f"the pipeline's own hand-off uses {theirs} here"
    return ion, quality, "; ".join(x for x in (why, note) if x)


def job(name: str, tf_seq: str, *, dna_top: str | None = None, ion: str | None = None) -> dict:
    """One alphafoldserver job: the homodimer, plus either a DNA duplex or an ion."""
    seqs: list[dict] = [{"proteinChain": {"sequence": tf_seq.upper(), "count": 2}}]
    if dna_top:
        seqs.append({"dnaSequence": {"sequence": dna_top.upper(), "count": 1}})
        seqs.append({"dnaSequence": {"sequence": revcomp(dna_top.upper()), "count": 1}})
    if ion:
        seqs.append({"ion": {"ion": ion, "count": IONS_PER_DIMER}})
    return {"name": sanitize_job_name(name), "modelSeeds": [1], "sequences": seqs,
            "dialect": "alphafoldserver", "version": 1}


def validate(jobs: list[dict]) -> None:
    """Raise on anything the server would reject. Unlike the predictor's validator this accepts both
    shapes: a DNA job (two complementary chains, no ion) and a metal job (one ion, no DNA)."""
    problems = []
    seen: set[str] = set()
    for j in jobs:
        nm = j.get("name", "<unnamed>")
        if nm in seen:
            problems.append(f"{nm}: duplicate job name")
        seen.add(nm)
        if j.get("dialect") != "alphafoldserver" or j.get("version") != 1:
            problems.append(f"{nm}: bad dialect/version")
        seqs = j.get("sequences", [])
        for s in seqs:
            for key, ent in s.items():
                if "count" not in ent:
                    problems.append(f"{nm}: {key} entity is missing required 'count'")
        prots = [s["proteinChain"] for s in seqs if "proteinChain" in s]
        dnas = [s["dnaSequence"] for s in seqs if "dnaSequence" in s]
        ions = [s["ion"] for s in seqs if "ion" in s]
        if len(prots) != 1 or prots[0].get("count", 1) != 2:
            problems.append(f"{nm}: expected one homodimeric proteinChain (count 2)")
        elif any(c not in "ACDEFGHIKLMNPQRSTVWY" for c in prots[0]["sequence"]):
            problems.append(f"{nm}: protein has non-standard residues")
        if dnas and ions:
            problems.append(f"{nm}: this hand-off keeps DNA and metal in separate jobs")
        if not dnas and not ions:
            problems.append(f"{nm}: neither DNA nor ion")
        if dnas:
            if len(dnas) != 2:
                problems.append(f"{nm}: expected two dnaSequence chains, got {len(dnas)}")
            else:
                top, bot = dnas[0]["sequence"].upper(), dnas[1]["sequence"].upper()
                if any(b not in "ACGT" for b in top + bot):
                    problems.append(f"{nm}: DNA has non-ACGT characters")
                elif revcomp(top) != bot:
                    problems.append(f"{nm}: the two DNA chains are not reverse-complements")
        for io in ions:
            if io.get("ion") not in AF3_IONS:
                problems.append(f"{nm}: unsupported ion {io.get('ion')!r}")
    if problems:
        raise SystemExit("invalid AF3 jobs:\n  " + "\n  ".join(problems))


def main() -> int:
    cfg = load_canonical()
    root = Path(cfg["runs"]["survey"]["root"])
    N = json.loads((PROJ / "reports" / "numbers.json").read_text(encoding="utf-8"))
    selected = N["context"]["selected"]

    jobs, rows = [], []
    for s in selected:
        bundle = root / "jobs" / s["run_name"]
        hit = GS.ranked_hits(GS.load_dossier(bundle))[0]
        pipeline_jobs = json.loads((bundle / "af3_jobs.json").read_text(encoding="utf-8"))
        first = pipeline_jobs[0]
        # the bundle's rank-1 hand-off must be the operator Figure 4 draws, or the duplex below
        # would belong to a different site
        if f"{hit['start']}-{hit['end']}" not in first["name"]:
            raise SystemExit(f"{s['run_name']}: the bundle's first AF3 job is not the rank-1 hit")
        tf_seq = next(x["proteinChain"]["sequence"] for x in first["sequences"] if "proteinChain" in x)
        dna_top = next(x["dnaSequence"]["sequence"] for x in first["sequences"] if "dnaSequence" in x)
        op = s["operator"]["seq"].upper()
        if op not in dna_top.upper() and revcomp(op) not in dna_top.upper():
            raise SystemExit(f"{s['run_name']}: the operator is not inside its own {DNA_LEN} bp construct")

        ion, quality, why = af3_ion(s["ion"])
        tag = s["run_name"].split("__")[0]
        gene = s["tf_gene"]
        apo, holo = f"{gene}__{tag}__apo_DNA", f"{gene}__{tag}__holo_{ion}"
        jobs.append(job(apo, tf_seq, dna_top=dna_top))
        jobs.append(job(holo, tf_seq, ion=ion))
        common = {"tf": gene, "genome": s["genome"], "family": s["family"], "clade": s["clade"],
                  "run_name": s["run_name"], "protein_len": len(tf_seq)}
        rows.append({**common, "job": sanitize_job_name(apo), "kind": "apo + operator DNA",
                     "ion": "", "called_metal": "", "substitution": "",
                     "dna_bp": len(dna_top), "operator": op,
                     "operator_locus": f"{s['replicon']}:{s['operator']['start']}-{s['operator']['end']}",
                     "note": f"DNA-binding form for this family: {DNA_BOUND_FORM.get(s['family'], '?')}"})
        rows.append({**common, "job": sanitize_job_name(holo), "kind": "holo + metal",
                     "ion": f"{ion} x{IONS_PER_DIMER}", "called_metal": s["ion"],
                     "substitution": quality, "dna_bp": 0, "operator": "", "operator_locus": "",
                     "note": why})

    validate(jobs)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "af3_figure4_jobs.json").write_text(json.dumps(jobs, indent=1) + "\n", encoding="utf-8", newline="\n")
    cols = ["tf", "genome", "family", "clade", "job", "kind", "called_metal", "ion", "substitution",
            "dna_bp", "operator", "operator_locus", "protein_len", "run_name", "note"]
    with (OUT_DIR / "af3_figure4_manifest.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    p = N["provenance"]
    readme = [
        "# AlphaFold 3 jobs for Figure 4",
        "",
        f"Generated by `analysis/make_af3_figure4_jobs.py` from the canonical run "
        f"({', '.join(r['tag'] for r in p['runs'].values())}, predictor {p['predictor_commit']}).",
        "Upload `af3_figure4_jobs.json` to alphafoldserver.com; `af3_figure4_manifest.tsv` says what "
        "each job is.",
        "",
        f"{len(jobs)} jobs, two per regulator:",
        "",
        "* **apo + operator DNA** -- the homodimer with the 60 bp centred on its primary operator "
        "(both strands, as the server wants a duplex), no metal.",
        f"* **holo + metal** -- the homodimer with its called effector metal, {IONS_PER_DIMER} ions "
        "(one per protomer), no DNA.",
        "",
        "## Ion substitutions",
        "",
        "The server accepts exactly ten ions -- MG ZN CL CA NA MN K FE (Fe3+) CU (Cu2+) CO -- so a "
        "called metal outside that set needs a surrogate chosen for coordination:",
        "",
        "| called | job ion | match | why |",
        "|---|---|---|---|",
    ]
    for called, (ion, quality, _) in sorted(SURROGATE.items()):
        used = [r for r in rows if r["called_metal"] == called]
        if used:
            readme.append(f"| {called} | {ion} | {quality} | {af3_ion(called)[2]} |")
    readme += [
        "",
        "## Which form binds DNA",
        "",
        "The apo + DNA job is the physiological complex only for the families whose apo form holds "
        "the operator (ArsR/SmtB, CsoR/FrmR, MerR). For Fur and DtxR/MntR the metal-loaded form is "
        "the DNA-binding one, so there the apo + DNA job is an empty-site control and the "
        "physiological complex would be a holo + DNA job (not generated here):",
        "",
        "| regulator | family | DNA-binding form |",
        "|---|---|---|",
    ]
    for r in rows:
        if r["kind"].startswith("apo"):
            readme.append(f"| {r['tf']} | {r['family']} | {DNA_BOUND_FORM.get(r['family'], '?')} |")
    readme += ["", "Every protein chain is the full-length sequence from the run manifest, not "
               "BITACORA's trimmed record (IdeR is 229 aa here against a 154 aa trim).", ""]
    (OUT_DIR / "README.md").write_text("\n".join(readme), encoding="utf-8", newline="\n")

    print(f"wrote {len(jobs)} jobs for {len(selected)} regulators -> {OUT_DIR}")
    for r in rows:
        extra = f"  [{r['called_metal']} -> {r['ion']}, {r['substitution']}]" if r["ion"] else \
                f"  [{r['dna_bp']} bp]"
        print(f"  {r['job']:<44} {r['kind']}{extra}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
