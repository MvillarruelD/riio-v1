#!/usr/bin/env python
"""Apply the improved regulon reconstruction to the paper's four genomes, offline.

The regulon step is the only thing that changed (genome-wide FDR, duplicate-site merge, real operon
assembly, promoter look-through, ncRNA-aware gene model). Everything upstream of it -- motif finding,
homolog gathering, the inducer call, structure -- is untouched, so a full online re-run would recompute
identical answers for those stages at ~20-110 min per TF. Instead this rebuilds the genome from the local
mirror, reuses each bundle's own operator PWM, and regenerates just the regulon. Minutes, no network.

Writes `regulon_v2.json` into every bundle (leaving the original `dossier.json` untouched) plus a summary
TSV across all 150 candidates.

    $TFOP_PY regenerate_regulons.py
    $TFOP_PY regenerate_regulons.py --genome Mtub
"""
import argparse
import csv
import json
import sys
import time
from pathlib import Path

from project_config import PREDICTOR_ROOT

PROJ = Path(__file__).resolve().parent.parent
REPO = PREDICTOR_ROOT
sys.path.insert(0, str(REPO))

from predictor import resources                           # noqa: E402
from predictor.annotate import context as CTX             # noqa: E402
from predictor.annotate.genome_db import MIRROR_DIR       # noqa: E402
from predictor.signals import regulon as RG               # noqa: E402
from predictor.signals.motif_rescan import rescan_genome  # noqa: E402

# the report chain's output locations, so a figure is never drawn from an archived run
sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_paths as RP


def _manifest_path():
    """The manifest for THIS run. Defined once, in run_paths."""
    return RP.MANIFEST

MIRROR = MIRROR_DIR
JOBS = RP.JOBS
MANIFEST = _manifest_path()


def _refresh_bundle_manifest(bundle: Path) -> None:
    """Keep a bundle's manifest checksum inventory exact after adding regulon_v2.json.

    The predictor writes bundle_manifest.json at prediction time; this stage adds regulon_v2.json
    afterwards. Under the schema-v2 bundle contract the manifest's `artifacts` inventory must cover
    exactly the bundle files, so a file added downstream has to be folded in or the bundle fails
    validation (and share publication). Recompute the inventory with the predictor's own helper so
    analysis and predictor agree on the checksum rule; the manifest excludes itself from inventory,
    so rewriting it here does not invalidate what it records.
    """
    mpath = bundle / "bundle_manifest.json"
    if not mpath.is_file():
        return
    try:
        from predictor.report import bundle_contract
    except Exception:
        return
    manifest = json.loads(mpath.read_text(encoding="utf-8"))
    if "artifacts" not in manifest:
        return
    manifest["artifacts"] = bundle_contract.inventory(bundle)
    mpath.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
OUT = RP.REGULON_SUMMARY

_CTX: dict = {}
_SCAN: dict = {}          # (genome accession, PWM bytes) -> genome scan; see the loop below


def genome_ctx(acc: str):
    """GenomeContext rebuilt from the mirrored GFF/FNA, cached per accession (4 genomes, 150 TFs)."""
    if acc not in _CTX:
        gff, fna = MIRROR / f"{acc}.gff", MIRROR / f"{acc}.fna"
        if not gff.is_file():
            _CTX[acc] = None
        else:
            _CTX[acc] = CTX.from_gff(str(gff), str(fna) if fna.is_file() else None)
    return _CTX[acc]


#: manifest `family` -> the packaged family-consensus PWM.
#: Rrf2 and CopY are absent on purpose: the vendored operator index carries no operators for them
#: (checked -- the index spans LuxR/GntR/TetR/LacI/LysR/AraC/Fur/MarR/ArsR/MerR/ArgR/IclR/TrpR/DtxR/
#: NiaR/VanR/CsoR/SsuR/HxlR), so no family consensus can be built and those TFs are reported unavailable
#: rather than given a fabricated motif.
_FAMILY_PWM = {
    "ArsR": "ArsR_SmtB", "MerR": "MerR", "Fur": "Fur", "GntR": "GntR", "TetR": "TetR_AcrR",
    "MarR": "MarR_SlyA", "LysR": "LysR_type__LTTR_", "CsoR": "CsoR_FrmR", "DtxR": "DtxR_MntR",
    "LacI": "LacI_GalR", "AraC": "AraC_XylS",
}
_PWM_DIR = resources.REFS_DIR / "family_pwms"


def load_pwm(bundle: Path, family: str = ""):
    """The best available motif for this TF, and where it came from.

    `per_tf_motif` -- the motif the pipeline built for THIS protein from its own gathered operator set.
    `family_prior` -- the family consensus PWM, used when the TF has no completed bundle. Strictly weaker
    (it cannot capture a TF's own specificity) and labelled as such everywhere it is reported.
    """
    import numpy as np
    p = bundle / "motif" / "three_logos.json"
    if p.is_file():
        d = json.loads(p.read_text(encoding="utf-8"))
        for key in ("gathered", "homolog", "tf_only"):
            blk = d.get(key) or {}
            if blk.get("pwm"):
                m = np.asarray(blk["pwm"], dtype=float)
                return (m if m.shape[0] == 4 else m.T), f"per_tf_motif:{key}"
    stem = _FAMILY_PWM.get((family or "").split("/")[0].strip())
    if stem:
        f = _PWM_DIR / f"{stem}.npy"
        if f.is_file():
            m = np.load(f)
            return (m if m.shape[0] == 4 else m.T), f"family_prior:{stem}"
    return None, None


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", help="only run names starting with this tag (Mtub/Mavi/Vcho/Vvul)")
    ap.add_argument("--resume", action="store_true",
                    help="skip TFs that already have regulon_v2.json (monotonic progress across calls)")
    ap.add_argument("--limit", type=int, default=0,
                    help="process at most N TFs this invocation, then write and exit")
    a = ap.parse_args(argv)

    rows = list(csv.DictReader(MANIFEST.open(encoding="utf-8")))
    if a.genome:
        rows = [r for r in rows if r["run_name"].startswith(a.genome)]
    OUT.parent.mkdir(parents=True, exist_ok=True)

    # Long detached runs get killed in this environment, so the script is built to be driven in chunks:
    # --resume skips finished TFs and --limit bounds one invocation, making repeated short calls add up.
    if a.resume:
        rows = [r for r in rows
                if not (JOBS / r["run_name"] / "regulon_v2.json").is_file()
                and not (OUT.parent / "family_prior" / f"{r['run_name']}.json").is_file()]
        print(f"resume: {len(rows)} TF(s) still to do")
    if a.limit:
        rows = rows[:a.limit]

    out, t0, missing, skipped = [], time.time(), 0, []
    for i, r in enumerate(rows, 1):
        name, acc = r["run_name"], r["organism_acc"]
        bundle = JOBS / name
        ctx = genome_ctx(acc)
        pwm, src = load_pwm(bundle, r.get("family", ""))
        if ctx is None or pwm is None:
            missing += 1
            skipped.append((name, r.get("family", "?"),
                            "no genome" if ctx is None else "no per-TF motif and no family PWM"))
            continue
        try:
            # A genome scan is the whole cost here (~40-60 s). Two TFs that share a family-prior PWM AND a
            # genome scan the identical motif against the identical sequence, so the result is the same --
            # cache on (genome, motif bytes) and the family-prior tier collapses from ~69 scans to ~20.
            key = (acc, pwm.tobytes())
            hits = _SCAN.get(key)
            if hits is None:
                hits = rescan_genome(ctx, pwm=pwm, scope="intergenic", pvalue_thresh=1e-4)
                _SCAN[key] = hits
            rg = RG.reconstruct_regulon(ctx, hits=hits)
        except Exception as e:
            print(f"  {name}: {type(e).__name__}: {e}")
            continue
        payload = {
            "run_name": name, "genome": acc, "pwm_source": src,
            "n_sites_scanned": len(hits), "n_operons": len(rg.operons),
            "target_genes": rg.target_genes,
            "operons": [{"first_gene": o.first_gene, "strand": o.strand, "genes": o.genes,
                         "site": [o.site_start, o.site_end], "score": o.score,
                         "qvalue": o.qvalue, "tier": o.tier} for o in rg.operons],
            "generated_by": "regenerate_regulons.py (offline; genome-wide FDR + operon assembly)",
        }
        # Per-TF-motif results go next to their bundle. Family-prior TFs have NO bundle (they were never
        # run), and fabricating an empty bundle directory for them would corrupt the "how many TFs have
        # actually been run" count that drives the batch's own resume logic -- so they are written to a
        # clearly separate tree instead.
        if bundle.is_dir():
            dest = bundle / "regulon_v2.json"
        else:
            dest = OUT.parent / "family_prior" / f"{name}.json"
            dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(json.dumps(payload, indent=1), encoding="utf-8")
        if bundle.is_dir():
            # regulon_v2.json is a new bundle artifact added after the predictor wrote the manifest.
            # Fold it into the manifest's checksum inventory so the schema-v2 bundle contract still
            # sees the inventory covering exactly the bundle files (else share publication rejects it).
            _refresh_bundle_manifest(bundle)
        out.append(dict(run_name=name, genome=acc, family=r.get("family_flag", ""),
                        pwm_source=src, n_sites=len(hits), n_operons=len(rg.operons),
                        n_genes=len(rg.target_genes),
                        top_operon=(rg.operons[0].first_gene if rg.operons else ""),
                        genes=";".join(rg.target_genes[:30])))
        if i % 20 == 0:
            print(f"  {i}/{len(rows)}  ({time.time()-t0:.0f}s)", flush=True)

    # merge with any previous chunk rather than overwriting it -- chunked runs must accumulate
    if out:
        merged = {}
        if OUT.is_file():
            for old in csv.DictReader(OUT.open(encoding="utf-8"), delimiter="\t"):
                merged[old["run_name"]] = old
        for r in out:
            merged[r["run_name"]] = {k: str(v) for k, v in r.items()}
        allrows = [merged[k] for k in sorted(merged)]
        with OUT.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t")
            w.writeheader(); w.writerows(allrows)
        print(f"  summary now holds {len(allrows)} TF(s)")
    by_genome: dict = {}
    for r in out:
        by_genome.setdefault(r["genome"], []).append(r)
    print(f"\nregenerated {len(out)}/{len(rows)} regulons in {time.time()-t0:.0f}s")
    n_own = sum(1 for r in out if r["pwm_source"].startswith("per_tf"))
    print(f"  motif source: {n_own} per-TF motif, {len(out)-n_own} family prior")
    for g, rs in sorted(by_genome.items()):
        ng = sum(r["n_genes"] for r in rs)
        print(f"  {g:20s} {len(rs):>4d} TFs   {ng:>6d} predicted regulon genes   "
              f"median operons/TF = {sorted(r['n_operons'] for r in rs)[len(rs)//2]}")
    if skipped:
        print(f"\n  {len(skipped)} TF(s) with NO usable motif (reported unavailable, not fabricated):")
        for n, fam, why in skipped:
            print(f"    {n:44s} family={fam:6s} {why}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
