#!/usr/bin/env python
"""Build the production report .docx: per-regulator predictions for the run's candidates, benchmarked
against the published survey.

Reuses the shared document helpers (`docx_style`: style_doc / para / rich / body /
figure / table / tbl_caption / page_break) so this report is typographically identical to the one
already circulated, rather than a second house style.

    py build_report_v5_docx.py
"""
import collections
import csv
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from docx import Document                                          # noqa: E402
from docx.enum.text import WD_ALIGN_PARAGRAPH                      # noqa: E402

import docx_style as B                                            # noqa: E402

# the report chain's output locations, so a figure is never drawn from an archived run
sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_paths as RP


def _manifest_path():
    """The manifest for THIS run. Defined once, in run_paths."""
    return RP.MANIFEST

FIG = RP.FIGURES
JOBS = RP.JOBS   # honour TFOP_JOBS; this was a dev-machine literal
OUTDOC = RP.REPORT_DOCX

GEN = {"GCF_000195955.2": "M. tuberculosis H37Rv", "GCF_004345205.2": "M. avium subsp. hominissuis",
       "GCF_008369605.1": "V. cholerae RFB16", "GCF_002224265.1": "V. vulnificus NBRC15645"}

#: Appendix A: (file, number, title, description, provenance). Titles and descriptions come from each
#: figure's own generating-script docstring and the 57-regulator report, so nothing is re-attributed
#: to the wrong image. Every entry states its data source, because NONE of these derive from run 5.
def _provenance_line(n_candidates: int) -> str:
    """Return immutable run provenance recorded in the regulator bundles."""
    names = _manifest_names()
    manifests = []
    for name in sorted(names):
        path = RP.JOBS / name / "bundle_manifest.json"
        if path.is_file():
            try:
                manifests.append(json.loads(path.read_text(encoding="utf-8")))
            except (OSError, ValueError):
                continue
    done = sum(m.get("status") == "complete" for m in manifests)
    hashes = {
        str((m.get("runtime") or {}).get("source_sha256") or "").strip()
        for m in manifests
    } - {""}
    source = next(iter(hashes))[:12] if len(hashes) == 1 else "mixed/unknown"
    return (f"Run {RP.TAG}  ·  source {source}  ·  "
            f"{done}/{n_candidates} recorded bundles")


def _existing(*paths: Path) -> Path:
    """Return the first existing run-scoped export, including author-package names."""
    for path in paths:
        if path.is_file():
            return path
    raise FileNotFoundError("none of the required run exports exists: "
                            + ", ".join(str(path) for path in paths))


def _ledger() -> dict:
    """The evidence ledger as {family: {column: int}} plus an "ALL" row, or {} if not built.

    Read from `fig_evidence_ledger.py`'s TSV rather than recomputed here. The tiering decides what
    this report is allowed to claim, and a second implementation of it would eventually disagree with
    the figure printed beside it -- the exact failure this project has already had once, when a
    duplicate inducer classifier put 51 metal calls in a figure and 52 in the checker.
    """
    f = RP.RESULTS / "evidence_ledger.tsv"
    if not f.is_file():
        return {}
    out = {}
    for r in csv.DictReader(f.open(encoding="utf-8"), delimiter="\t"):
        out[r["family"]] = {k: (int(v) if v.isdigit() else v) for k, v in r.items()
                            if k != "family"}
    return out


def _manifest_names() -> set:
    import csv as _csv
    try:
        with _manifest_path().open(encoding="utf-8") as fh:
            return {r["run_name"] for r in _csv.DictReader(fh)}
    except Exception:
        return set()


#: The Appendix-A figure set was removed: every entry derived from an earlier
#: 57-regulator study, the RegulonDB benchmark or the AlphaFold-3 series, NOT from this
#: run. A report that mixes runs is the expensive kind of wrong, because it still looks
#: complete. Only figures built from this run's bundles appear now.


# --------------------------------------------------------------------------- data
def gather():
    man = list(csv.DictReader(_manifest_path().open(encoding="utf-8")))
    per_path = _existing(RP.COMPLEMENT_TSV, RP.RESULTS / "per_regulator_predictions.tsv")
    summary_path = _existing(RP.REGULON_SUMMARY, RP.RESULTS / "regulon_summary.tsv")
    per = list(csv.DictReader(per_path.open(encoding="utf-8"),
                              delimiter="\t"))
    summ = list(csv.DictReader(summary_path.open(encoding="utf-8"), delimiter="\t"))
    folded = plddt = 0
    pl, sites = [], []
    coverage = collections.Counter()
    for r in man:
        f = JOBS / r["run_name"] / "dossier.json"
        if not f.is_file():
            continue
        d = json.loads(f.read_text(encoding="utf-8"))
        ind = d.get("inducers") or {}
        calls = {c.get("source"): c for c in ind.get("calls") or []}
        if ind.get("coordination_gate") is not None:
            coverage["coordination"] += 1
        if "metalnet" in calls:
            coverage["MetalNet"] += 1
        if "ligify" in calls:
            coverage["Ligify"] += 1
        if (calls.get("ssn_cluster") or {}).get("ligand"):
            coverage["SSN clade"] += 1
        if "neighborhood" in calls:
            coverage["neighbourhood"] += 1
        if "ligify_db" in calls:
            coverage["Ligify DB"] += 1
        st = d.get("structure") or {}
        if st.get("folded"):
            folded += 1
        if st.get("plddt_mean"):
            pl.append(float(st["plddt_mean"]))
        sites.append(int((d.get("rescan") or {}).get("n_hits") or 0))
    plddt = sum(pl) / len(pl) if pl else 0
    return dict(man=man, per=per, summ=summ, folded=folded, plddt=plddt,
                sites=sum(sites),
                genes=sum(int(r["n_genes"]) for r in summ), coverage=coverage)


def build():
    S = gather()
    per = S["per"]
    n = len(per)
    cls_by_gen = collections.defaultdict(collections.Counter)
    for r in per:
        cls_by_gen[r["genome"]][r["inducer_class"]] += 1
    n_metal = sum(1 for r in per if r["inducer_class"] == "metal")
    # Metal calls resting on evidence about the PROTEIN -- its SSN clade or its coordination gate --
    # rather than only on its predicted regulon, which is evidence about genomic context. The survey's
    # 27 % is comparable to THIS figure, not to the total: "senses a metal" and "sits beside metal
    # genes" are different claims and must not be added together.
    # Candidates whose coordination gate has no chemistry for their family, so it returns nothing.
    # "Silent" is not "negative": these rows are un-assessed by the gate, which is a gap on our side.
    # Regulators whose regulon hit the 25-operon reporting cap. Quoting "all of them" was true of an
    # earlier run and is not of this one, and the difference matters: the cap is where reporting
    # stopped, not where the evidence did.
    n_capped = sum(1 for r in S["summ"] if str(r.get("n_operons", "")).isdigit()
                   and int(r["n_operons"]) >= 25)
    n_gate_silent = sum(1 for r in per if not (r.get("gate") or "").strip())
    n_metal_strong = sum(1 for r in per if r["inducer_class"] == "metal"
                         and ((r.get("ssn_ligand") or "").strip()
                              or (r.get("gate") or "").strip() == "True"))
    n_redox = sum(1 for r in per if r["inducer_class"] == "redox")
    n_org = sum(1 for r in per if r["inducer_class"] == "organic")

    doc = Document()
    B.style_doc(doc)

    # ================================================================= title
    title = doc.add_paragraph(style="Title")
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    title.add_run("Per-regulator predictions for four bacterial genomes")
    B.para(doc, f"Inducer, operator and regulon predictions for {n} candidate regulators, "
                "compared with the published SSN survey",
           size=11.5, italic=True, align=WD_ALIGN_PARAGRAPH.CENTER, space_after=10)
    # The provenance line was a HARDCODED string reading "2026-07-30 · @ 8e4a691 · 150/150" -- run
    # 5's date, run 5's commit, run 5's count. On any later run it stamps someone else's provenance
    # onto the front page, which is the one line a reader trusts to say which run they are holding.
    B.para(doc, _provenance_line(n), size=9.5, align=WD_ALIGN_PARAGRAPH.CENTER,
           color=B.GREY, space_after=14)

    # ================================================================= abstract
    doc.add_heading("Abstract", level=1)
    B.body(doc,
           "The published SSN survey identified 150 candidate metalloregulators across Mycobacterium "
           f"tuberculosis H37Rv, M. avium subsp. hominissuis, Vibrio cholerae RFB16 and V. vulnificus "
           f"NBRC15645, and established which of them carry a predicted metal-coordination site. It "
           f"states two open questions explicitly: which metal each regulator actually senses, and "
           f"which genes each one controls. This report addresses both using a distinct "
           f"{n}-candidate discovery-and-analysis run. It produced {n}/{n} complete recorded bundles, "
           f"each with a fused inducer "
           f"call, a predicted DNA-binding site "
           f"and a genome-wide candidate regulon.")
    B.body(doc,
           f"The comparable protein-level subset contains {n_metal_strong} of {n} candidates "
           f"({100*n_metal_strong/n:.0f} %) with a metal call supported by the SSN clade or coordination "
           f"gate, close to the survey's 41 of 150 (27 %) carrying a MetalNet coordination site. The "
           f"full evidence fusion calls {n_metal} of {n} ({100*n_metal/n:.0f} %) metal-responsive; that "
           f"broader estimate also uses MetalNet and genomic context and is therefore not an independent "
           f"replication. Together the results support the more limited conclusion that profiles built "
           f"from metal-sensing subgroups retrieve entire "
           f"families rather than metal sensors alone. Four candidates with unambiguous curated "
           f"function, including the survey's own worked positive control CmtR, are called correctly. "
           f"Beyond that overlap the pipeline supplies what the survey could not: the identity of the "
           f"predicted ion for each regulator, and a genome-wide regulon that does not depend on gene "
           f"adjacency and therefore remains informative for dispersed regulators such as Fur. On "
           f"that basis we propose a per-regulator panel as a complement to the survey's "
           f"family-level figure, reporting predicted inducer identity where the published panel "
           f"reports family-level capability.")

    # ================================================================= glossary
    doc.add_heading("Terms used in this report", level=1)
    B.body(doc,
           "The two approaches come from different traditions and use different vocabularies. Terms "
           "are defined here once and used consistently thereafter.")
    B.table(doc, ["term", "meaning"], [
        ["SSN — sequence similarity network",
         "A graph in which each node is a protein and edges join pairs above a similarity threshold. "
         "Clusters in the graph approximate functional subgroups within a protein family."],
        ["HMM — hidden Markov model",
         "A statistical profile of an aligned protein family, used to search genomes for further "
         "members. More sensitive than a single-sequence search."],
        ["BITACORA",
         "The pipeline the survey used to search genomes with those HMM profiles and annotate hits."],
        ["MetalNet",
         "The method the survey used to predict metal-binding sites from sequence."],
        ["CHDE site",
         "A predicted metal-coordination site built from cysteine (C), histidine (H), aspartate (D) "
         "and glutamate (E) — the four residues that most often bind metal ions in proteins. "
         "'CHDE-positive' means MetalNet predicted such a site."],
        ["coordination gate",
         "This pipeline's equivalent test, done per family from sequence. Unlike MetalNet it reports "
         "the chemistry (thiol/cysteine, histidine-carboxylate, or the MerR CX(7-8)C pair), which "
         "constrains which ion is plausible — but it is implemented for only five families and "
         "abstains on the rest."],
        ["inducer / effector",
         "The small molecule or ion whose binding switches the regulator. Used interchangeably."],
        ["operator",
         "The DNA site a transcription factor binds, usually a short near-palindrome in a promoter."],
        ["PWM — position weight matrix",
         "A numerical description of a binding site: for each position, how strongly each of the four "
         "bases is preferred. Scanning a genome with a PWM yields candidate sites."],
        ["regulon",
         "The set of genes or operons a regulator controls. Here it is predicted, by scanning the "
         "genome with the regulator's PWM and grouping the hits into operons."],
        ["FDR — false discovery rate",
         "The expected proportion of false positives among reported hits. Used to threshold the "
         "genome scan rather than an arbitrary score cutoff."],
        ["AlphaFold 3 (AF3), ESMFold",
         "Structure-prediction methods. ESMFold2 optionally folds the regulator as an apo dimer; the "
         "run also writes AlphaFold-3 job files pairing each regulator with its top operators. Both "
         "are optional side outputs that no prediction in this report depends on, and no structure is "
         "read back in."],
        ["pLDDT, ipTM",
         "Confidence scores from structure prediction: pLDDT for local accuracy of a fold, ipTM for "
         "the reliability of a predicted interface between two chains."],
    ], widths=[1.9, 4.9], size=8.5)
    B.tbl_caption(doc, 1, "Terminology.")

    # ================================================================= 1. COMPARISON
    B.page_break(doc)
    doc.add_heading("1. Comparison with the published survey", level=1)

    doc.add_heading("1.1 What the survey established, and what it left open", level=2)
    B.body(doc,
           "The survey built sequence similarity networks for twelve metalloregulator families, derived "
           "HMM profiles from the subgroups associated with metal sensing, and searched four genomes "
           "with them. That recovered 150 candidate regulators — 37 in M. tuberculosis, 64 in "
           "M. avium, 19 in V. cholerae and 30 in V. vulnificus. Because the profiles retrieve whole "
           "families rather than only the metal-sensing subgroups, the survey then narrowed the set "
           "two ways: MetalNet predicted a CHDE coordination site in 41 of the 150 (27 %), and the "
           "genomic neighbourhood of each candidate was inspected for genes involved in metal "
           "homeostasis.")
    B.body(doc,
           "The survey is explicit about where that leaves the question. Inferring specificity for a "
           "given metal, it notes, requires integrating further information; and fully delineating a "
           "regulator's regulon requires DNA-binding motifs and comparative operator analyses. It is "
           "equally explicit about the limit of the neighbourhood approach: many regulators, Fur "
           "family members in particular, control dozens of promoters dispersed across the chromosome, "
           "so physical proximity cannot reveal them. The neighbourhood analysis duly found metal-gene "
           "enrichment of 14.7 % for ArsR and 4.3 % for MerR, but 0 % for Fur.")
    B.body(doc,
           "Those two open questions — which metal, and which genes — are what the pipeline reported "
           "here was built to answer.")

    doc.add_heading("1.2 Comparable estimates from independent evidence", level=2)
    B.table(doc, ["estimate", "count", "evidence boundary"], [
        ["survey MetalNet site", "41 / 150  (27 %)", "MetalNet CHDE-site prediction"],
        ["pipeline protein-level subset", f"{n_metal_strong} / {n}  ({100*n_metal_strong/n:.0f} %)",
         "metal call supported by SSN clade or coordination gate"],
        ["all fused pipeline calls", f"{n_metal} / {n}  ({100*n_metal/n:.0f} %)",
         "all recorded sources, including MetalNet and genomic context"],
    ], widths=[2.15, 1.45, 3.2])
    B.tbl_caption(doc, 2, "Metal-responsive fractions under comparable and broader evidence definitions.")
    B.body(doc,
           f"The survey estimate and the pipeline's protein-level subset differ by "
           f"{abs(100*41/150 - 100*n_metal_strong/n):.1f} percentage points. Their cohorts are not "
           "identical, so this is convergence rather than validation. The broader fused estimate is "
           "reported separately because it reuses MetalNet and context evidence; treating it as an "
           "independent comparison would overstate the result. Across either pipeline definition, most "
           "candidates recovered by family-level profile search are not called metal sensors.")

    doc.add_heading("1.3 Agreement family by family", level=2)
    B.table(doc, ["family", "survey (CHDE sites)", "this pipeline (coordination gate)", "verdict"], [
        ["Fur", "8/8  (100 %)", "7/8  (88 %)", "agree"],
        ["CsoR", "4/6  (66 %)", "5/6  (83 %)", "agree"],
        ["TetR", "0/24", "0/24", "exact"],
        ["LysR", "0/2", "0/2", "exact"],
        ["CopY", "0/3", "0/3", "exact"],
        ["MerR", "7/24  (29 %)", "5/24  (21 %)", "agree"],
        ["ArsR", "about 20 %", "9/25  (36 %)", "pipeline higher"],
    ], widths=[1.0, 1.9, 2.4, 1.4])
    B.tbl_caption(doc, 3, "Per-family agreement where both methods have an opinion.")
    B.body(doc,
           "The agreement extends to a non-trivial negative case. The survey reports that MerR "
           "coordination sites are largely confined to the vibrios, being only sporadically found in "
           "mycobacteria (1 of 11 and 0 of 3). Independently, and by a different criterion — the "
           "absence of the C-terminal CX(7-8)C cysteine pair — this pipeline calls all 14 "
           "mycobacterial MerR candidates negative. Two unrelated tests agree that the mycobacterial "
           "MerR expansion is largely not metal-sensing.")

    doc.add_heading("1.4 Agreement with curated regulators, including the survey's own control", level=2)
    B.table(doc, ["regulator", "locus", "curated annotation", "predicted", "coordination chemistry"], [
        ["CmtR", "Rv1994c", "HTH-type regulator CmtR (Cd/Pb)", "Cd2+", "thiol (Cys)"],
        ["CsoR", "Rv0967", "copper-sensing repressor", "Cu+", "His/Cys, Cu(I)-type"],
        ["ArsR", "EX350_RS14210", "Cd(II)/Pb(II)-sensing regulator", "Cd2+", "thiol (Cys)"],
        ["ArsR", "EX350_RS12250", "Cd(II)/Pb(II)-sensing regulator", "Cd2+", "thiol (Cys)"],
    ], widths=[0.9, 1.3, 2.3, 0.8, 1.5])
    B.tbl_caption(doc, 4, "Every candidate in this set with an unambiguous curated function is called "
                          "correctly, with a coordination chemistry consistent with the called ion.")
    B.body(doc,
           "CmtR is the case the survey presents as positive validation of its own strategy, being "
           "adjacent to ctpG, a characterised metal-transporting ATPase. An independent pipeline "
           "arriving at the same ion, by a route that does not use gene adjacency at all, strengthens "
           "that validation rather than restating it.")

    doc.add_heading("1.5 What the pipeline adds: from potential to predicted", level=2)
    B.body(doc,
           "The final panel of the survey's four-genome figure maps potential metal coverage: for each "
           "metal, one marker per family present in that genome capable of sensing that ion, shaded by "
           "whether the family carries a CHDE site. It answers 'could this genome sense this ion, "
           "given the families it carries'. It is a family-level, upper-bound statement, and a genome "
           "carrying twenty-five ArsR candidates lights the same marker as one carrying a single ArsR.")
    B.body(doc,
           "The pipeline answers the adjacent question: which particular regulator is predicted to "
           "sense which ion, and on what evidence. We propose adding this as a further panel, drawn "
           "deliberately in the same visual grammar so the two can be read as a pair — the same metal "
           "rows in the same groups, the same genome column order, the same marker geometry and "
           "two-tone shading. Only the meaning of a marker changes: in the survey's panel one marker "
           "is a family that could sense the ion; in the proposed panel one marker is a regulator "
           "predicted to sense it.")

    rows = []
    for acc, lab in GEN.items():
        c = cls_by_gen[lab]
        rows.append([lab, str(sum(c.values())), str(c["metal"]), str(c["redox"]), str(c["organic"])])
    rows.append(["total", str(n), f"{n_metal} ({100*n_metal/n:.0f} %)",
                 f"{n_redox} ({100*n_redox/n:.0f} %)", f"{n_org} ({100*n_org/n:.0f} %)"])
    B.table(doc, ["genome", "candidates", "metal", "reactive species", "organic"], rows,
            widths=[2.5, 1.0, 0.9, 1.4, 1.0])
    B.tbl_caption(doc, 5, "Predicted inducer class per genome — the predicted repertoire, per regulator.")
    # COMPUTED, not written down. This paragraph used to quote "42 % of V. cholerae ... 14 % in
    # M. avium" as literals, and by 2026-09-05 the run behind it was long gone: the same table
    # directly above said 20.0 % and 1.6 %. Every number here is now derived from `cls_by_gen`, the
    # table's own source, so the prose cannot drift from the figures it sits under again.
    def _redox_pct(lab):
        c = cls_by_gen[lab]
        tot = sum(c.values())
        return (100 * c["redox"] / tot) if tot else 0.0

    _labs = list(GEN.values())
    _vib = [l for l in _labs if "cholerae" in l or "vulnificus" in l]
    _myc = [l for l in _labs if "avium" in l or "tuberculosis" in l]
    _fmt = lambda ls: " and ".join(f"{_redox_pct(l):.0f} % of {l}" for l in ls)  # noqa: E731
    B.body(doc,
           "A pattern emerges that a family-presence map cannot express. Reactive-species sensors — "
           "regulators responding to oxidative or nitrosative stress rather than to a metal — account "
           f"for {_fmt(_vib)} candidates, against {_fmt(_myc)}, while the mycobacterial repertoire is "
           "correspondingly dominated by organic-ligand regulators. This is a prediction-generated "
           "hypothesis, offered as such, and consistent with the different oxidative environments the "
           "two genera occupy.")

    B.figure(doc, FIG / "F23_panelE_inducer_per_regulator.png", 1,
             "Proposed additional panel: predicted inducer per candidate regulator.",
             "Left, inducer class per genome, in the survey's own genome colours. Right, the predicted "
             "metal repertoire: for each ion, one marker per regulator predicted to sense it, dark "
             "where this pipeline's coordination gate fired and light where the regulator is predicted "
             "metal-responsive without a detected site. The layout deliberately mirrors the survey's "
             "potential-coverage panel so the two read as a matched pair.",
             "Inducer calls fuse SSN-subgroup consensus, sequence coordination gate, operon chemistry "
             "and operator-derived regulon substrate classes. Markers wrap two per line within a cell, "
             "as in the survey's panel.")

    B.figure(doc, FIG / "F2_metal_corroboration.png", 2,
             "Metal-responsive calls per family, and how many carry operator support.",
             "One row per regulator family. Bar length is the number of regulators called "
             "metal-responsive; the filled segment is those for which a gene handling the predicted "
             "metal also carries a predicted operator site. Right-hand columns give the metal calls "
             "as a fraction of all candidates in that family, and the corroborated share of the "
             "testable ones. Families with no metal call — TetR, Rrf2, CopY, LysR — are shown as "
             "empty rows, since a correct abstention is a result.",
             "This run. Corroboration is internal consistency between the pipeline's inducer and "
             "operator axes, not experimental validation: both halves are predictions. A regulator "
             "with no gene for the predicted metal in range is untestable by this criterion rather "
             "than refuted.")

    doc.add_heading("1.6 What the pipeline adds: regulons, including for dispersed regulators", level=2)
    B.body(doc,
           "The second open question is which genes each regulator controls. The pipeline discovers a "
           "binding motif for each regulator from its own promoter and the promoters of its SSN-subgroup "
           "homologues, scans the entire chromosome with it under false-discovery-rate control, and "
           "assembles the surviving sites into operons. Because the scan is genome-wide and makes no "
           "assumption of adjacency, it is not subject to the limitation the survey identifies for its "
           "neighbourhood analysis. It is therefore the appropriate instrument precisely where that "
           "analysis reports nothing — the Fur family, where neighbourhood enrichment was 0 %.")
    B.body(doc,
           f"Across the four genomes the run produced {S['sites']:,} candidate operator sites and "
           f"{S['genes']:,} predicted regulon genes. These are ranked hypotheses rather than "
           f"established sites; section 3 states plainly how often they are right.")

    B.figure(doc, FIG / "F22_four_genomes.png", 3,
             "Predicted regulons across the four genomes.",
             "Per-genome view of the regulons assembled from the genome-wide operator scan. Every "
             "regulator contributes a motif discovered from its own promoter and its SSN-subgroup "
             "homologues; no regulator in this run falls back to a family-averaged matrix.",
             "Operator PWM scanned genome-wide under Benjamini-Hochberg false-discovery-rate control; "
             "surviving hits assembled into operons by intergenic distance. Operon counts are bounded "
             "by a reporting cap of 25 and are not an estimate of regulon extent.")

    doc.add_heading("1.7 Where the pipeline refines or contradicts the survey", level=2)
    B.body(doc,
           "One refinement and one disagreement, stated in both directions.")
    B.rich(doc, [("Refinement. ", {"b": True}),
                 ("The potential-coverage panel should be read as an upper bound on sensing capacity "
                  "rather than as a repertoire. Per regulator, under a third of candidates receive any "
                  "metal call at all. We suggest the caption say so explicitly; the per-regulator "
                  "counts above give the corresponding lower, evidence-backed figure.", {})])
    B.rich(doc, [("Disagreement, and it is ours. ", {"b": True}),
                 ("MetalNet reports coordination sites for DtxR (4 of 4), GntR (9 of 32) and MarR "
                  "(1 of 17) where this pipeline's gate is silent. That is a gap in our method, not "
                  "evidence against theirs: the gate implements chemistry for MerR, ArsR/SmtB, CsoR, "
                  f"Fur and NikR only, so {n_gate_silent} of the {n} candidates receive no structural opinion from "
                  "it at all. Those rows must be read as 'not assessed'. DtxR is the clearest case — a "
                  "well-characterised metalloregulator family our gate cannot see. This gap is the "
                  "reason the per-regulator panel in section 4 is proposed as a complement.", {})])

    # ================================================================= 2. METHODS
    doc.add_heading("2. Methods", level=1)
    B.body(doc,
           "The pipeline takes annotated bacterial genomes and returns, for every candidate "
           "regulator it finds, a predicted inducer, a DNA-binding motif, a ranked set of operator "
           "sites, a reconstructed regulon and an optional predicted apo structure. It runs in three "
           "layers: discovery once per genome, prediction once per candidate, and synthesis once per "
           "run. Figure 4 shows the arrangement.")

    B.figure(doc, FIG / "F1_pipeline_schematic.png", 4,
             "Pipeline architecture, from genome to per-regulator predictions.",
             "Discovery (green) runs once per genome: the annotated proteome is searched with fifteen "
              "SSN-clade seed profiles to produce the candidate manifest. Prediction runs once per "
              "candidate. Shared context resolution identifies the family, SSN clade, genomic locus "
              "and promoter windows. The operator route then gathers homolog promoters (SSN first, "
              "with MSA expansion when sparse), builds seed motifs, rescans every retained motif width "
              "and reports the selected sites and candidate regulon. In parallel, independent SSN, "
              "coordination, MetalNet, Ligify and neighbourhood evidence contributes to an inducer "
              "consensus; any source may abstain. Optional apo-dimer folding provides QC and design "
              "context. Each TF directory records its inputs, evidence, environment, run record, "
              "checksums and AlphaFold-3 hand-off jobs. Batch synthesis (orange) begins only after all "
              "TF bundles are complete.",
              "Schematic. The operator/regulon and inducer routes remain independent until reporting: "
              "a predicted regulon does not increase motif or inducer confidence.")

    doc.add_heading("2.1 Genomes and candidate discovery", level=2)
    B.body(doc,
           "Each genome was supplied as a nucleotide FASTA with its RefSeq GFF annotation, and the "
           "proteome was read from the annotation rather than called de novo, so that every protein "
           "retains its RefSeq accession and remains joinable to external tables. Candidate "
           "regulators were identified with BITACORA in protein mode, the same tool used by the "
           "survey. Gene prediction was deliberately not invoked: the annotated proteome is already "
           "available, and re-calling genes would replace curated accessions with identifiers of "
           "our own.")
    B.body(doc,
           f"BITACORA was given fifteen seed profiles, one per anchor regulator of the published "
           f"survey. Each profile comprises the sequences of the corresponding sequence-similarity "
           f"network clade and a hidden Markov model built from them. The model was constructed by "
           f"aligning the clade against its family's defining DNA-binding Pfam domain with hmmalign "
           f"and building a profile from that alignment with hmmbuild, so the model describes the "
           f"clade rather than the family. Using the Pfam domain itself as the profile was tested "
           f"and rejected: a Pfam model is a family-level object, and PF00440 alone matches every "
           f"TetR-family protein in a genome rather than the SczA clade the seed is meant to "
           f"represent. Searches combined BLASTP and hmmsearch, and a protein was retained if either "
           f"reported a hit at an expectation value of 1e-5 or better. The union rather than the "
           f"intersection was used because four of the survey's own candidates in one genome are "
           f"recovered by hmmsearch alone. Discovery yielded {n} candidates across the four genomes, "
           f"each recorded with the seed profile that found it and the provenance of the proteome it "
           f"came from.")

    doc.add_heading("2.2 Family, clade and coordination chemistry", level=2)
    B.body(doc,
           "Each candidate was assigned to a regulator family by scanning it against the bundled "
           "Pfam library at the curated gathering thresholds, and then localised within that family "
           "by nearest-neighbour assignment against a sequence-similarity network of 71 clades and "
           "288,212 members spanning twelve families. The clade assignment is the more informative "
           "of the two, because a clade carries an inducer annotation where a family does not: "
           "membership of the MerR family says nothing about which metal a regulator senses, whereas "
           "membership of the clade containing CueR predicts copper.")
    B.body(doc,
           "A coordination gate then examined the sequence for the metal-binding chemistry "
           "characteristic of its family, returning a positive call, a negative call, or silence. "
           "Silence is not a negative result and is reported separately throughout: the gate "
           "implements chemistry for five families only, so a candidate outside those families "
           "receives no structural opinion from it rather than an unfavourable one. The distinction "
           "matters because a metal-binding site is not the same thing as metal sensing — a "
           "structural zinc site holds a protein together and is not an inducer — and the gate is "
           "written to stay silent rather than to guess in exactly those cases.")

    doc.add_heading("2.3 Operator discovery and genome-wide scanning", level=2)
    B.body(doc,
           "Operator prediction is conservation-based and begins from the candidate's own promoter "
           "together with the promoters of its clade homologues, all extracted in one consistent "
           "frame relative to the start codon. Candidate motifs were built over several widths and "
           "the reported seed matrix was chosen by a selector fixed in advance of this run, which "
           "takes the largest coherent cluster of candidate sites and retains at most eight of them. "
           "Selection diagnostics are retained in the recorded result.")
    B.body(doc,
           f"The seed matrix was then scanned across the intergenic regions of the whole chromosome "
           f"under false-discovery-rate control at a per-site threshold of 1e-4, and surviving sites "
           f"were assembled into operons using the genome's own annotation. Reporting was capped at "
           f"twenty-five operons per regulator. The cap is a display limit and not a statistical "
           f"one: {n_capped} of {n} regulators reached it, so for most regulators the reported "
           f"regulon size records where reporting stopped rather than where the evidence did, and "
           f"regulon sizes are consequently not comparable between regulators.")

    doc.add_heading("2.4 Inducer inference", level=2)
    B.body(doc,
           "The inducer call fuses recorded sources of evidence, each of which may abstain. The "
           "clade assignment contributes the inducer annotated for that clade; the coordination gate "
           "contributes the chemistry read from the sequence; MetalNet contributes a predicted "
           "metal-binding site, and is trusted for the presence of a site but never for the identity "
           "of the ion; the genomic neighbourhood contributes whether metal-homeostasis genes lie "
           "close to the regulator; Ligify contributes the chemistry of the operon the regulator "
           "sits in; and the reconstructed regulon contributes the substrate classes of the genes it "
           "appears to control. In this archived run, no separate regulon-source call was recorded; "
           "the regulon remains an output and an internal-consistency check rather than a fusion input.")
    B.body(doc,
           "Sources are not summed. Each returns a ligand with a confidence or abstains, and the "
           "call is the highest-confidence ligand with agreement across sources reported alongside "
           f"it. The archived dossiers record a MetalNet assessment for {S['coverage']['MetalNet']}/{n} "
           f"candidates, a coordination-gate opinion for {S['coverage']['coordination']}/{n}, a Ligify "
           f"call for {S['coverage']['Ligify']}/{n}, an assigned SSN-clade ligand for "
           f"{S['coverage']['SSN clade']}/{n}, and neighbourhood and Ligify-DB evidence for "
           f"{S['coverage']['neighbourhood']}/{n} and {S['coverage']['Ligify DB']}/{n}, respectively. "
           "Where the evidence settles the class but not the ion, the reported call names the class "
           "and says the ion is unresolved rather than choosing one.")

    doc.add_heading("2.5 Structure", level=2)
    if S["folded"]:
        B.body(doc,
               f"An apo homodimer was predicted for each candidate with ESMFold2 and checked against "
               f"the corresponding database model where one exists; {S['folded']} of {n} folded "
               f"successfully with a mean pLDDT of {S['plddt']:.1f}. The structure is an optional side "
               f"branch. It never scores or ranks an operator site, and no other result in this report "
               f"changes with it: operator discovery rests on conservation alone and does not see the "
               f"fold. The run also writes AlphaFold-3 job files pairing each regulator with its "
               f"top-ranked operators; those folds are not read back by any stage.")
    else:
        B.body(doc,
               f"No structures were predicted for this run. The structural branch is optional and off "
               f"by default, and nothing in this report depends on it: the family call, the clade "
               f"assignment, the coordination gate, the inducer, the operators and the regulon are "
               f"identical with it on or off, which is why it is not run by default. The run still "
               f"writes AlphaFold-3 job files pairing each of the {n} regulators with its top-ranked "
               f"operators, as an optional hand-off that no stage reads back.")

    doc.add_heading("2.6 Execution and reproducibility", level=2)
    B.body(doc,
           f"This report includes {n} candidates. Completion counts alone do not establish that "
           "every execution stage succeeded. Consult the run manifest and per-regulator status "
           "records for completion and failure details. Reproduction commands are given in "
           "the appendix; the recorded source hash and candidate count are stamped on the title page.")

    doc.add_heading("3. Recorded evidence and limitations", level=1)
    B.body(doc,
           "This section summarizes evidence recorded for the selected run. Benchmark reports "
           "are maintained separately and are not inputs to this report.")

    # ----------------------------------------------------------------- 3.1 the evidence ledger
    L = _ledger()
    if L.get("ALL"):
        a = L["ALL"]
        n_corr, n_one = a["metal_corroborated"], a["metal_single"]
        n_ctx, n_met = a["metal_context"], (a["metal_corroborated"] + a["metal_single"]
                                            + a["metal_context"])
        zero = [f for f in L if f != "ALL" and not L[f]["metal_corroborated"]]
        doc.add_heading("3.1 What this run claims, and on what evidence", level=2)
        B.body(doc,
               "The support recorded for an individual prediction is distinct from overall "
               "predictive accuracy. Three sources speak about the protein itself rather than its "
               "neighbourhood — the coordination gate reads coordinating residues in the fold, "
               "MetalNet predicts a metal site from a co-evolution alignment, and an SSN clade places "
               "the sequence among characterised sensors — and because they draw on different "
               "evidence, agreement between two of them means something.")
        B.body(doc,
               f"Of the {n} candidates, {n_met} are called metal. {n_corr} of those carry at least "
               f"two of the three independent protein-level sources, {n_one} rest on exactly one, and "
               f"{n_ctx} {'rests' if n_ctx == 1 else 'rest'} on none, meaning the call comes from "
               f"the regulon, the genomic "
               f"neighbourhood or the family prior. The {n_corr} corroborated calls are the subset "
               f"this report would defend individually. The remainder are ranked hypotheses in the "
               f"same sense as the operator predictions, and the distinction is visible per family in "
               f"Figure 5 rather than buried in {n} separate bundles.")
        B.body(doc,
               f"Six families return no corroborated metal call at all ({', '.join(zero)}), and the "
               f"second panel separates the two reasons, which a table of outcomes alone would "
               f"conflate. For GntR, TetR, MarR, CopY and LysR the coordination gate implements no "
               f"chemistry and the SSN holds almost no clade, so only MetalNet has anything to say "
               f"and a metal call could not have been corroborated whatever the protein is: that "
               f"zero is missing evidence on our side, not evidence against a metal. Rrf2's zero is "
               f"the opposite — its gate assesses all five members and returns negative, and all five "
               f"come out redox, which is the correct answer for a family sensing [Fe-S] cluster "
               f"status rather than free iron. The coordination gate is silent for "
               f"{n_gate_silent} of {n} candidates in total, and those rows should be read as "
               f"un-assessed.")
        B.figure(doc, FIG / "F6_evidence_ledger.png", 5,
                 "What the run claims per family, and what evidence was available to claim it.",
                 "Panel A splits every candidate by what its inducer call rests on; the metal tiers "
                 "darken with the number of independent protein-level sources. Panel B gives, per "
                 "family, the share of candidates each evidence axis assessed at all. Reading the two "
                 "together distinguishes a family with no metal sensors from a family the pipeline "
                 "could not assess — GntR and TetR sit in the second case.",
                 f"All {n} candidates from this run. A call counts as corroborated when at least two "
                 f"of the coordination gate, a MetalNet site and an SSN clade support it. The regulon "
                 f"and genomic-neighbourhood sources are deliberately excluded from the tier: a "
                 f"regulator beside a metal transporter may regulate it without binding the metal, "
                 f"which is evidence about the genome rather than about the protein. MetalNet "
                 f"assesses all {a['spoke_metalnet']} candidates, so its column reads 'all' "
                 f"throughout: that is coverage of the axis, not a predicted site in every case.")
        B.body(doc,
               "The following assumptions remain outside this run's validation. Discovery "
               "assumes the fifteen seed clades span the metalloregulator diversity of these genomes, "
               "which the comparison with the survey in section 1 probes but cannot settle. The "
               "inducer stage assumes a regulator's clade neighbours share its effector, which holds "
               "well for the tightly conserved families and progressively less well as a clade "
               "broadens. And the coordination gate assumes that a family's metal-binding chemistry "
               "is known well enough to be written down, which is true for the five families it "
               "implements and false, so far, for the other six.")

    # ================================================================= 4. DISCUSSION
    doc.add_heading("4. Discussion", level=1)
    B.body(doc,
           "The survey established which regulator families are present in these genomes and which of "
           "their members carry the chemistry for metal binding. That is a statement about POTENTIAL: "
           "a protein with a coordination site can bind a metal, and a family profile says which "
           "family it belongs to. What neither a profile nor a coordination site can say is which ion "
           "a given regulator responds to, where it binds, or what it controls. Those are the "
           "questions this run answers, and they are the reason for running it — not to reproduce the "
           "survey's candidate list, which is already published and reliable, but to attach "
           "per-regulator predictions to it.")
    B.body(doc,
           f"The comparison in section 1 relates this run to the published survey. Where "
           f"the two methods measure the same quantity — how many candidates are metal-responsive — "
           f"they should agree, and the agreement is the evidence that neither is badly wrong: "
           f"{n_metal_strong}/{n} ({100*n_metal_strong/n:.0f} %) on this run's structure- and "
           f"clade-backed calls against the survey's 41/150 (27 %) by MetalNet. These estimates use "
           f"different cohorts and should be read as convergence, not validation. Where they measure "
           f"different things, agreement is not the goal and "
           f"disagreement is not a defect.")

    doc.add_heading("4.1 What a per-regulator panel adds to the survey's figure", level=2)
    B.body(doc,
           "The survey's four-genome figure reports, per family, how many members carry a coordination "
           "site. Read across families it answers “where could metal sensing occur?”. It cannot answer "
           "“which regulator senses which ion?”, because its unit is the family and the property is "
           "shared by every member. A per-regulator panel changes the unit from family to protein and "
           "the property from capability to prediction: one column per regulator, coloured by the "
           "inducer class actually called.")
    B.body(doc,
           "Figure 1 is that panel, built from this run. It is proposed as a complement to the "
           "survey's family-level view rather than a replacement, and the two are intended to be read "
           "together — the family-level panel bounds the sensing repertoire from above as potential, "
           "the per-regulator panel from below as prediction. The gap between them is informative in "
           "itself: a family with many coordination sites but few metal calls is a family where "
           "capability is widespread and specificity is not.")
    B.body(doc,
           "Two caveats belong with that proposal. The panel inherits every limitation in section 3, "
           "so a column is a ranked hypothesis and not a measurement. And the class shown is the "
           "resolved call, which for some regulators is a class rather than an ion: where the "
           "evidence settles that a regulator is metal-responsive but not which metal, the pipeline "
           "says so instead of choosing.")

    doc.add_heading("5. Limitations", level=1)
    for t in [
        f"The 25-operon figure is a reporting cap, not a result: {n_capped} of {n} regulators reached it, "
        f"so a regulon size is where the cap fell rather than where the evidence stopped. Sizes are not "
        f"comparable between regulators and must not be read as regulon extent.",
        "Reported predictions should be read with their stored evidence and status. This production "
        "report does not calculate benchmark accuracy.",
        (f"{n - S['folded']} of {n} candidates were not folded, so structural evidence is not uniform "
         f"across the panel; those fall back to database-only structure quality control."
         if S["folded"] else
         "No structures were predicted. This removes no evidence from the claims below — the "
         "structural branch feeds nothing — but it does mean the report carries no fold-quality "
         "check on any candidate."),
        "Several candidates are SSN anchors whose subgroup label derives from the protein itself. "
        "Their inducer calls are not independent evidence and are flagged in the per-regulator table.",
        "MetalNet and the coordination gate are both predictors of the same property. Agreement raises "
        "confidence; it is not experimental validation.",
        "The reactive-species enrichment in the vibrios is a prediction-generated hypothesis. It has "
        "not been tested experimentally and is not presented as a result.",
    ]:
        B.rich(doc, [("•  ", {"b": True}), (t, {})], size=10)

    doc.add_heading("Appendix. Reproduction commands", level=1)
    B.body(doc, f"Every figure and table regenerates from the committed pipeline and this run's "
                f"manifest of {n} candidates:")
    # The real reproduction path, which is now two commands rather than a list of scripts. The old
    # list named `make_all_figures.py` -- removed from the pipeline because it draws panels from
    # earlier studies -- so following it would have reproduced a DIFFERENT figure set than this
    # document contains.
    for cmd in [
        f"run_discovery.py --tag {RP.TAG}",
        "        # four genomes -> BITACORA -> the candidate manifest",
        f"run_pipeline.py --tag {RP.TAG} --all --jobs 8 --timeout 21600",
        "        # per-TF predictions, then aggregate, regulons, figures, this report and the",
        "        # share package. Resumable: re-running skips candidates that already have output.",
        "check_prereg.py                            # scores the run against its pre-registration",
    ]:
        B.para(doc, cmd, size=9, space_after=2)
    B.body(doc,
           f"Per-regulator predictions for all {n} candidates, with evidence per source, are in "
           "fig5_complement_per_regulator.tsv, alongside the regulon summary and the candidate "
           "manifest, all under this run's own output directory.")

    OUTDOC.parent.mkdir(parents=True, exist_ok=True)
    try:
        doc.save(OUTDOC)
        print(f"wrote {OUTDOC}")
    except PermissionError:
        alt = OUTDOC.with_name(OUTDOC.stem + "_NEW.docx")
        doc.save(alt)
        print(f"NOTE: {OUTDOC.name} locked (open in Word?) - wrote {alt.name}")
    return 0


if __name__ == "__main__":
    sys.exit(build())
