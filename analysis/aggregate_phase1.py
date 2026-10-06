#!/usr/bin/env python
"""
aggregate_phase1.py -- collapse the per-TF bundles into one master table.

Field semantics follow the skill's outputs_and_representations.md, NOT guesswork:
  * `autoregulatory_recovery.recovered == false` is a LEGITIMATE finding (weak/non-autoregulator),
    not a pipeline failure -- reported as-is, never as an error.
  * every row keeps its `provenance` tag (PREDICTED_UNVERIFIED / DERIVED / CURATED) verbatim.
  * selectivity = genome-wide putative sites per Mb; >=120/Mb flags a degenerate (low-information)
    operator whose hits are probably noise.
"""
from __future__ import annotations
import csv, json


# One definition of which candidates this run is about; see run_paths.MANIFEST for why there
# is no local fallback.
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent))
import manifest as MF  # noqa: E402
import run_paths as RP  # noqa: E402

PROJ = _Path(__file__).resolve().parent.parent


def _manifest_path():
    """The manifest for THIS run. Defined once, in run_paths."""
    return RP.MANIFEST

# Resolved from the run tag, never from a literal path: a figure drawn from an archived
# run into a current report is the expensive kind of error, because the report still
# looks finished.
JOBS = RP.JOBS
OUT = RP.RESULTS

def _bundle_complete(name: str) -> bool:
    bundle = JOBS / name
    try:
        status = json.loads((bundle / "bundle_manifest.json").read_text(
            encoding="utf-8")).get("status")
    except (OSError, ValueError):
        return False
    return status == "complete" and all((bundle / rel).is_file() for rel in (
        "dossier.json", "operators_ranked.json", "binding_sites.tsv", "regulon.tsv",
    ))



def main() -> int:
    """Collapse this run's bundles into one master table. Returns a process exit code."""
    # Validated, not merely read: this stage is routinely run on its own to rebuild the master
    # table, and a bare DictReader here meant a direct invocation skipped every manifest check.
    man = {r["run_name"]: r for r in MF.load_or_exit(_manifest_path(), PROJ)}
    missing = [name for name in man if not _bundle_complete(name)]
    if missing:
        preview = ", ".join(missing[:5])
        raise SystemExit(
            f"refusing to aggregate an incomplete run: {len(missing)}/{len(man)} bundle(s) are missing "
            f"required artifacts ({preview}{', ...' if len(missing) > 5 else ''})"
        )
    rows = []
    for name, m in man.items():
        d = JOBS / name / "dossier.json"
        if not d.exists():
            continue
        j = json.loads(d.read_text(encoding="utf-8"))
        ops = j.get("operators") or {}
        prim = (ops.get("primary") or {})
        ind = j.get("inducers") or {}
        st = j.get("structure") or {}
        # The pre-2026-09-02 key `mode_of_regulation` is no longer read: every canonical bundle
        # carries only `promoter_occlusion`, and the predictor dropped the fallback too.
        mor = j.get("promoter_occlusion") or {}
        aut = j.get("autoregulatory_recovery") or {}
        logo = j.get("genomic_logo") or {}
        cons = j.get("conservation") or {}
        n_sites = (j.get("rescan") or {}).get("n_hits") or len(logo.get("hits") or [])
        glen = 4.4e6
        calls = {c.get("source"): c.get("ligand") for c in (ind.get("calls") or [])}
        rows.append({
            "run_name": name, "genome": m["genome"], "locus_tag": m.get("locus_tag", ""),
            "gene": m.get("gene", ""), "product": (m.get("product", "") or "")[:60],
            "family": j.get("family"), "ssn_cluster": j.get("ssn_cluster"),
            "primary_operator": prim.get("sequence"), "operator_source": prim.get("source"),
            "operator_score": round(prim.get("score", 0) or 0, 2), "score_type": prim.get("score_type"),
            "n_operator_candidates": len(ops.get("candidates") or ops.get("ranked") or []),
            "motif_consensus": logo.get("consensus"),
            "motif_mean_ic": round((logo.get("total_ic") or 0) / max(1, len(logo.get("per_col_ic") or [1])), 3),
            "inducer_top": ind.get("top"), "inducer_agreement": ind.get("agreement"),
            "coordination_gate": ind.get("coordination_gate"),
            "inducer_from_cluster": calls.get("ssn_cluster"), "inducer_from_coordination": calls.get("coordination"),
            "inducer_from_ligify": calls.get("ligify"), "ion": j.get("ion"),
            "mode": mor.get("mode") or j.get("mode"), "mode_confidence": mor.get("confidence"),
            "autoreg_recovered": aut.get("recovered"), "autoreg_nearest_bp": aut.get("nearest_bp"),
            "regulon_operons": len(j.get("regulon") or []),
            "confirmed_regulated_genes": sum(1 for g in (j.get("confirmed_regulated_genes") or []) if g.get("confirmed")),
            "n_putative_sites": n_sites,
            "selectivity_per_Mb": round(n_sites / (glen / 1e6), 1) if n_sites else None,
            "n_homologs": cons.get("n_homologs"), "structure_route": j.get("structure_route"),
            "folded": st.get("folded"), "plddt": round(st.get("plddt_mean") or 0, 1) if st.get("plddt_mean") else None,
            "fold_qc_rmsd": st.get("qc_rmsd"), "fold_qc_pass": st.get("qc_pass", st.get("qc")),
            "af3_jobs": len(j.get("af3_jobs") or []), "provenance": j.get("provenance"),
            "n_caveats": len(j.get("caveats") or []),
        })

    rows.sort(key=lambda r: (r["genome"], r["family"], r["run_name"]))
    # Created here rather than at import: an import-time mkdir materialises the run's output directory
    # merely because something imported this module, which is how an empty, plausible-looking results
    # directory appears for a run that never happened.
    RP.ensure()
    if not rows:
        raise SystemExit("refusing to write an empty master table: no bundle yielded a "
                         "dossier.json. Run the batch stage first.")
    p = OUT / "operator_predictions.csv"
    with p.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f"wrote {p}  ({len(rows)} TFs)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
