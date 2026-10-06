#!/usr/bin/env python
"""Complement manuscript Figure 5 (four-genome candidate regulators) with per-REGULATOR inducer calls.

WHAT THE MANUSCRIPT'S FIGURE 5 DOES (v1.3, caption at para [97])
    A  regulators per family per genome (HMM/BITACORA, E <= 1e-5)
    B  total candidates per genome
    C  genomic organisation of three selected candidates
    D  "Potential coverage of metals by genome. For each metal ..., each square represents a FAMILY
       present in that genome capable of sensing that ion; the colour indicates whether the FAMILY
       possesses a predicted CHDE coordination site."

Panel D is therefore **family-level and potential**: it answers "could this genome sense this ion,
given the families it carries?". The manuscript states the corresponding limitation itself --
"inferring their specificity for a given metal requires integrating additional information" -- and, for
the genomic-context route, that Fur-like master regulators exert "distributed regulation across multiple
loci", which is why their context analysis found 0% metal-gene enrichment for Fur.

WHAT THIS ADDS
    A per-REGULATOR, evidence-fused inducer call for each of the same 150 candidates, from an
    independent pipeline: SSN-cluster consensus + a sequence coordination gate + Ligify operon
    chemistry + the substrate classes of an operator-derived regulon. That converts "this genome has a
    family that could sense Zn" into "this specific regulator is predicted to sense Zn, on this
    evidence" -- and, because the regulon is reconstructed genome-wide from the operator PWM rather than
    from adjacency, it is not blind to the distributed regulators their context analysis misses.

THREE INDEPENDENT AXES OF COMPARISON
  1. coordination gate  vs  MetalNet CHDE (41/150, ~27%): two independent structural predictions of the
     same property. Agreement corroborates; disagreement localises where each is weak.
  2. realised vs potential metal coverage: panel D counts a family square per genome; we count actual
     per-regulator metal calls. Expect realised << potential -- the interesting, testable difference.
  3. the three panel-C examples (CmtR/ctpG, M. avium ArsR+arsM, V. vulnificus CueR): direct agreement
     checks on named cases the manuscript already validates.

    py fig5_complement.py            # writes the TSV + prints the comparison
"""
import argparse
import collections
import csv
import json
import sys
from pathlib import Path


def _manifest_path():
    """The manifest for THIS run. Defined once, in run_paths."""
    return RP.MANIFEST

HERE = Path(__file__).resolve().parent
PROJ = HERE.parent
sys.path.insert(0, str(HERE))
import run_paths as RP  # noqa: E402
sys.path.insert(0, str(RP.REPO))
from inducer_evidence import display_inducer  # noqa: E402
JOBS = RP.JOBS   # honour TFOP_JOBS; this was a dev-machine literal
OUT = RP.RESULTS

#: manifest organism accession -> the genome as named in the manuscript
GENOME = {
    "GCF_000195955.2": "M. tuberculosis H37Rv",
    "GCF_004345205.2": "M. avium subsp. hominissuis",
    "GCF_008369605.1": "V. cholerae RFB16",
    "GCF_002224265.1": "V. vulnificus NBRC15645",
}

#: MetalNet CHDE coordination-site counts reported in the manuscript (para [86]), per family.
#: `None` = the manuscript gives no per-family figure (ArsR is quoted only as "about 20%").
MS_CHDE = {"Fur": (8, 8), "DtxR": (4, 4), "CsoR": (4, 6), "GntR": (9, 32), "MarR": (1, 17),
           "TetR": (0, 24), "CopY": (0, 3), "LysR": (0, 2), "ArsR": None, "MerR": None, "Rrf2": None}
MS_CHDE_TOTAL = (41, 150)

#: the three worked examples in panel C
MS_EXAMPLES = {
    "M. tuberculosis H37Rv": ("CmtR (ArsR) adjacent to ctpG, a metal-transporting ATPase", "cmtR"),
    "M. avium subsp. hominissuis": ("ArsR/SmtB adjacent to a heavy-metal ATPase and arsM", "arsR"),
    "V. vulnificus NBRC15645": ("CueR-type associated with an iron-uptake locus", "cueR"),
}

# Element symbols are matched ONLY as an ion token (`Zn2+`, `Cu+`, `Cd2+/Pb2+`) and metal names ONLY as
# whole words. A bare substring test is wrong here and silently inflates the metal count: "steroid-CoA"
# contains "co" (cobalt), "nitrogen" contains "ni" (nickel), "gluconate" contains "co" -- all three were
# classed as metals before this was fixed, which would have overstated the realised metal repertoire.
_ION_RE = __import__("re").compile(
    r"\b(Zn|Cu|Fe|Mn|Ni|Co|Cd|Hg|As|Pb|Ag|Au|Mo|Cr|Sb|Se|Te|W|V)\d?[+-]")
_METAL_NAME_RE = __import__("re").compile(
    r"\b(zinc|copper|iron|ferrous|ferric|manganese|nickel|cobalt|cadmium|mercur\w*|arsen\w*|"
    r"lead|silver|gold|molybd\w*|chrom\w*|antimon\w*|selen\w*|tellur\w*|tungsten|vanad\w*|"
    r"divalent metal|metalloid|heavy metal)\b")
_REDOX_RE = __import__("re").compile(
    r"(redox|reactive oxygen|reactive nitrogen|peroxide|rooh|superoxide|\brss\b|persulfide|"
    r"sulfane|nitrosat\w*|\bno\b|oxidative)")


def inducer_class(top: str) -> str:
    """metal / redox / organic / undetermined for a headline inducer string.

    Delegates to the pipeline's own vocabulary rather than re-deriving the classification here. The
    previous local version used its own regexes and read a bare element symbol as organic: `Fe`
    carries no charge, so the ion pattern missed it and the element-NAME pattern only matched
    spelled-out names, leaving one genuine iron call filed as an organic one. That is what a second
    copy of a classification always eventually does -- it drifts from the first, and the two disagree
    by a count nobody can explain (here, 51 metal calls against the checker's 52).
    """
    t = (top or "").strip()
    if not t:
        return "undetermined"
    import sys as _sys
    _sys.path.insert(0, str(RP.REPO))
    from predictor.effector import inducer_vocab as _v
    cls = _v.inducer_class_of(t)
    if cls:
        return cls
    # the vocabulary abstains on a cross-class mixture and on an empty/abstention label
    return "undetermined" if ("undetermined" in t.lower() or "unresolved" in t.lower()) else "organic"


def candidate_species(cands) -> list:
    """The candidate set as one species per entry, deduplicated on the canonical token.

    The consensus deduplicates on the raw string, so a source that stores co-sensing as a FLATTENED
    label puts that label into the set as a member beside its own components: one MerR came out
    `Cd2+/Pb2+ · Cd2+ · Pb2+ · Cu(+)(in) · Cu+`, five entries for three species. `parse_mixture`
    knows which slashes separate species and which belong to one name (`persulfide/RSS`), so the
    split is the vocabulary's decision and not a `str.split` here.
    """
    import sys as _sys
    _sys.path.insert(0, str(RP.REPO))
    from predictor.effector import inducer_vocab as _v
    out = []
    for c in cands or ():
        for part in _v.parse_mixture(c):
            if part not in out:
                out.append(part)
    return out


def collect():
    rows = []
    for r in csv.DictReader(_manifest_path().open(encoding="utf-8")):
        f = JOBS / r["run_name"] / "dossier.json"
        if not f.is_file():
            continue
        d = json.loads(f.read_text(encoding="utf-8"))
        ind = d.get("inducers") or {}
        calls = {c["source"]: c for c in ind.get("calls", [])}
        coord = calls.get("coordination") or {}
        ev = coord.get("evidence") or {}
        st = d.get("structure") or {}
        ops = d.get("operators") or {}
        prim = (ops.get("primary") or {}) if isinstance(ops, dict) else {}
        top = ind.get("top") or ""
        rows.append(dict(
            run_name=r["run_name"], family=r["family"], genome=GENOME.get(r["organism_acc"], "?"),
            gene=r.get("gene") or "", product=(r.get("product") or "")[:60],
            # `inducer` is presentation-only; keep the class-bearing token beside it so every count
            # remains tied to the exact stored call even for legacy bundles backfilled on the fly.
            inducer=display_inducer(ind), inducer_top=top, inducer_class=inducer_class(top),
            # The fused candidate set, carried BESIDE the headline and never instead of it. The
            # consensus retains every source's ligand and flags a genuine mixture. Quantitative
            # analysis reads `inducer_top`; this display column and the candidate column expose the
            # ambiguity without changing those calculations.
            candidates="/".join(candidate_species(ind.get("candidates"))),
            mixture_kind=ind.get("mixture_kind") or "",
            agreement=round(float(ind.get("agreement") or 0), 3),
            gate=("" if ind.get("coordination_gate") is None else str(ind.get("coordination_gate"))),
            gate_chem=ev.get("gate") or "", gate_chemistry=ev.get("chemistry") or "",
            ssn_cluster=(calls.get("ssn_cluster", {}).get("evidence") or {}).get("cluster") or "",
            ssn_ligand=calls.get("ssn_cluster", {}).get("ligand") or "",
            regulon_ligand=calls.get("regulon", {}).get("ligand") or "",
            folded=bool(st.get("folded")), plddt=st.get("plddt_mean"),
            primary_operator=prim.get("sequence") or "",
            n_operons=len([o for o in (d.get("regulon") or []) if isinstance(o, dict)]),
        ))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quiet", action="store_true")
    ap.parse_args()   # --quiet is accepted for symmetry with the other figure scripts

    rows = collect()
    if not rows:
        print("no dossiers yet")
        return 1
    OUT.mkdir(parents=True, exist_ok=True)
    tsv = RP.COMPLEMENT_TSV
    with tsv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)

    n = len(rows)
    print(f"per-regulator calls available: {n}/150\n")

    # --- axis 1: our coordination gate vs the manuscript's MetalNet CHDE -------------------------
    print("AXIS 1 -- coordination gate (ours) vs MetalNet CHDE (manuscript), per family")
    print(f"{'family':8s}{'n':>4s}{'gate+':>7s}{'gate-':>7s}{'silent':>8s}{'ours %':>8s}   {'MetalNet CHDE':>14s}")
    by_fam = collections.defaultdict(list)
    for r in rows:
        by_fam[r["family"]].append(r)
    g_pos = g_tot = 0
    for fam in sorted(by_fam):
        rs = by_fam[fam]
        pos = sum(1 for r in rs if r["gate"] == "True")
        neg = sum(1 for r in rs if r["gate"] == "False")
        sil = sum(1 for r in rs if r["gate"] == "")
        g_pos += pos; g_tot += len(rs)
        ms = MS_CHDE.get(fam)
        ms_s = f"{ms[0]}/{ms[1]} ({100*ms[0]//ms[1]}%)" if ms else "not reported"
        print(f"{fam:8s}{len(rs):>4d}{pos:>7d}{neg:>7d}{sil:>8d}{100*pos/len(rs):>7.0f}%   {ms_s:>14s}")
    print(f"{'TOTAL':8s}{g_tot:>4d}{g_pos:>7d}{'':>7s}{'':>8s}{100*g_pos/g_tot:>7.0f}%   "
          f"{MS_CHDE_TOTAL[0]}/{MS_CHDE_TOTAL[1]} ({100*MS_CHDE_TOTAL[0]//MS_CHDE_TOTAL[1]}%)")

    # --- axis 2: realised vs potential metal coverage -------------------------------------------
    print("\nAXIS 2 -- inducer CLASS per genome (realised, per regulator)")
    print(f"{'genome':30s}{'n':>4s}{'metal':>7s}{'redox':>7s}{'organic':>9s}{'undet':>7s}")
    by_gen = collections.defaultdict(list)
    for r in rows:
        by_gen[r["genome"]].append(r)
    for g in sorted(by_gen):
        rs = by_gen[g]
        c = collections.Counter(r["inducer_class"] for r in rs)
        print(f"{g:30s}{len(rs):>4d}{c['metal']:>7d}{c['redox']:>7d}{c['organic']:>9d}"
              f"{c['undetermined']:>7d}")

    print("\nAXIS 2b -- specific ions actually called (realised metal repertoire)")
    for g in sorted(by_gen):
        ions = collections.Counter(r["inducer_top"] for r in by_gen[g]
                                   if r["inducer_class"] == "metal")
        if ions:
            print(f"  {g:30s} " + ", ".join(f"{k}({v})" for k, v in ions.most_common()))

    # --- axis 3: the three worked examples ------------------------------------------------------
    print("\nAXIS 3 -- the manuscript's panel-C worked examples")
    for g, (desc, needle) in MS_EXAMPLES.items():
        cands = [r for r in by_gen.get(g, [])
                 if needle in (r["gene"] or "").lower() or needle in r["run_name"].lower()]
        print(f"  {g}\n    manuscript: {desc}")
        if not cands:
            print("    ours      : (regulator not yet run / gene name not matched)")
        for r in cands[:3]:
            print(f"    ours      : {r['run_name'][:40]} -> {r['inducer']} "
                  f"[{r['inducer_class']}] gate={r['gate'] or 'silent'} "
                  f"operons={r['n_operons']}")

    print(f"\nwrote {tsv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
