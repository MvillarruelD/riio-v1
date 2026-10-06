"""build_legends.py -- the figure list of the run document, and LEGENDS.md.

`figure_entries()` is the single list of figures the document carries: label, title, file stem, main
or supplementary, and legend. build_docx.py lays out exactly these entries, so the document cannot
carry a figure this list does not describe. Each legend comes from its `legend_*` function below;
until one is written it returns a visible placeholder. LEGEND_BRIEF.md says what each must contain.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
PLACEHOLDER = "[Legend to be written.]"

#: (file stem, label, title, main?) in document order
FIGURES = [
    ("fig1_pipeline", "Figure 1", "The genome-to-bundle prediction pipeline", True),
    ("fig2_benchmark", "Figure 2", "Benchmark on E. coli K-12 and S. Typhimurium SL1344", True),
    ("fig3_survey", "Figure 3", "Four-genome survey", True),
    ("fig4_context_structures", "Figure 4", "Operator-bound metalloregulators in genomic context", True),
    ("figS1_per_regulator", "Figure S1", "Every discovered regulator, axis by axis", False),
    ("figS2_own_promoter", "Figure S2", "Own-promoter recovery", False),
    ("figS3_pwm_quality", "Figure S3", "Operator and PWM quality", False),
    ("figS4_ssn_clade", "Figure S4", "The SSN clade and the metal call", False),
    ("figS5_survey_metal_table", "Figure S5", "The survey's metal calls", False),
    # S5bis: every survey candidate, one page per genome group (build_survey_table.BIS_PAGES)
    ("figS5bis_survey_table_1", "Figure S5bis (1/3)", "Every survey candidate: M. tuberculosis", False),
    ("figS5bis_survey_table_2", "Figure S5bis (2/3)", "Every survey candidate: M. avium", False),
    ("figS5bis_survey_table_3", "Figure S5bis (3/3)",
     "Every survey candidate: V. cholerae and V. vulnificus", False),
]


def provenance_line(N: dict) -> str:
    p = N["provenance"]
    tag = f" (tag {p['predictor_tag']})" if p.get("predictor_tag") else ""
    runs = ", ".join(r["tag"] for r in p["runs"].values())
    return (f"Runs {runs}; predictor commit {p['predictor_commit']}{tag}; "
            f"structure prediction {'on' if p.get('structure') else 'off'}.")


# ---------------------------------------------------------------- legends
# One function per figure, each taking numbers.json as `N` and returning the legend text. Every count,
# fraction or P value in a legend must be read from N, never typed. `**bold**` and `*italic*` spans are
# rendered by build_docx. What each legend must say is specified in LEGEND_BRIEF.md.
_WORDS = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight",
          9: "nine", 10: "ten", 11: "eleven", 12: "twelve"}


def w(n: int) -> str:
    """A small count as a word, the way Nature sets it; anything larger stays a numeral."""
    return _WORDS.get(n, str(n))


def pval(p: float) -> str:
    """1.72e-15 -> '1.7 x 10^-15'; anything down to 0.001 stays decimal."""
    if p >= 0.001:
        return f"{p:.3g}"
    mant, exp = f"{p:.1e}".split("e")
    return f"{mant} × 10^{int(exp)}"


def survey_genomes(N: dict) -> list[str]:
    """The survey genomes, in the order every figure uses."""
    from figure_style import SURVEY_ORDER
    per = N["survey"]["per_genome"]
    return [g for g in SURVEY_ORDER if g in per]


def _species(names) -> str:
    """Italicised genome names as a list: 'a', 'a and b', 'a, b and c'."""
    it = [f"*{g}*" for g in names]
    return it[0] if len(it) == 1 else ", ".join(it[:-1]) + " and " + it[-1]


def _channel_names(N: dict) -> list[str]:
    """The inducer channels, named as Figure S4 names them."""
    from build_supplementary_figures import CHANNEL_LABEL, CHANNELS
    return [CHANNEL_LABEL.get(c, c) for c in CHANNELS]


def _top_k(N: dict) -> int:
    """The k of the own-region operator rule (top-k operators)."""
    return N["survey"]["all"]["autoregulation"]["top_k"]


def expected(x: float) -> str:
    """A chance expectation: two decimals, or '< 0.01' rather than a misleading '0.00'."""
    return "< 0.01" if x < 0.005 else f"{x:.2f}"


def bench_genomes(N: dict) -> list[str]:
    """The benchmark genome keys that numbers.json carries, in figure order."""
    from figure_style import BENCH_GENOMES
    return [k for k, _ in BENCH_GENOMES if k in N["benchmark_genome_composition"]]


def n_all_genomes(N: dict) -> int:
    return len(survey_genomes(N)) + len(bench_genomes(N))


def caveat(N: dict) -> str:
    """The benchmark independence caveat, set for a legend: species italic, no file path."""
    text = re.sub(r"\s*\([^()]*\.py\)", "", N["benchmark"]["independence_caveat"])
    return text.replace("E. coli", "*E. coli*").replace(">=", "≥")


def legend_fig1(N: dict) -> str:
    chans = _channel_names(N)
    struct = "on" if N["provenance"].get("structure") else "off"
    return (
        "Schematic; no run data are plotted. Each step is drawn as the data it produces; boxed names are "
        "the established software that runs it (filled: the step's main engine). "
        "**a**, Candidate discovery, once per genome, in four numbered steps. The proteome is every "
        "protein of the genome's own annotation (NCBI Datasets). Selected metal-sensing clades of the "
        "sequence-similarity network (SSN) each give one profile HMM, and BITACORA searches the proteome "
        "with those HMMs (HMMER) and the clades' sequences (BLAST+) at E ≤ 1e-5. Candidates are the "
        "proteins caught (dark bars). The profiles resolve clades, not whole Pfam families, so a family "
        "member outside every selected clade never becomes a candidate. "
        "**b**, Per-TF prediction, one independent run per candidate. The candidate's locus (regulator "
        "dark, flanking intergenic DNA shaded) is assigned a family (pyhmmer, BLAST+) and an SSN clade "
        "(MMseqs2 homologue hits). The operator path (teal) takes the promoter regions of orthologues "
        "(teal windows; MMseqs2), aligns them to find a conserved inverted repeat (two facing half-sites) "
        "and builds a PWM from it, rescans the genome's intergenic DNA with that PWM (dashed line: the "
        "p ≤ 1e-4 cut), and ranks the kept sites; #1 is the primary operator. The regulon and the mode "
        "of regulation (PromoterCalculator) follow from the ranked sites. "
        f"The inducer path (blue) runs {w(len(chans))} "
        "evidence channels: the SSN clade's curated label, the Cys/His coordination motif, a predicted "
        "metal site (MetalNet2), metal-handling genes in the neighbourhood, and operon chemistry "
        "(Ligify). Each calls on its own or abstains (a dashed chip); the fused call keeps every call "
        "with its evidence behind one class and ion. Neither path reads the other's output. They meet "
        "only afterwards, drawn as the regulator's own locus: does an operator (teal) sit at the TF's "
        "own promoter, upstream of the regulator (blue), for a call that names an ion? "
        f"The structure add-on (dashed) is optional and was {struct} in this run: an "
        "ESMFold2 apo dimer, checked against the AlphaFold DB monomer, which adds a fold and changes no "
        "call. Each run ends in a per-TF bundle: motif and ranked operators, regulon, the inducer call "
        "with each channel's evidence, and the run record, plus AlphaFold 3 job files as a hand-off the "
        "pipeline never runs. Nothing computed later is fed back into the motif or the inducer call; "
        "every output is an unverified prediction."
    )


def legend_fig2(N: dict) -> str:
    b = N["benchmark"]
    e, s = b["ecoli"]["scores"], b["salmonella"]["scores"]
    cen = b["metal_sensor_accounting"]["ecoli_census_totals"]
    nfam = len(N["families"])
    k = N["benchmark_genome_composition"]["ecoli"]["autoregulation"]["top_k"]
    n_bench = len(bench_genomes(N))
    # discovered TFs no census family accounts for: named Pfam families, and unmatched candidates
    named = sorted({f for g in bench_genomes(N) for f in b[g]["discovered_outside_census"]
                    if f != "no UniProt match"})
    no_up = len(b["salmonella"]["discovered_outside_census"].get("no UniProt match", []))
    tan = "; ".join(([", ".join(named)] if named else [])
                    + ([f"{no_up} SL1344 candidates with no UniProt entry"] if no_up else [])) or "none"
    conf = [(r["truth_class"], r["pred_class"]) for g in bench_genomes(N)
            for r in b[g]["per_tf"] if r.get("truth_class")]
    diag = sum(1 for t, p in conf if t == p)
    eo, so = e["operator"]["all"]["1"], s["operator"]["all"]["1"]
    return (
        f"Every TF discovered in the {w(n_bench)} benchmark genomes, scored on each axis for which truth "
        "exists; n/N is printed on every mark. "
        f"**a**, Discovered/known per studied family. Census: *E. coli* K-12, {cen['tfs']} "
        f"transcription factors (RegulonDB {b['ecoli']['sources']['release']} plus UniProt proteins "
        f"carrying a studied family's signature Pfam), {cen['in_study_families']} in the {w(nfam)} "
        f"studied families; SL1344, {s['discovery']['in_family_total']} UniProtKB proteins carrying such a "
        "Pfam. In both, a protein the paper's own family sequence set lists counts in that family, which "
        "places ModE in LysR although UniProt cross-references only its molybdate-binding domain. "
        f"Tan: discovered TFs outside the census families ({tan}). "
        "**b**, Inducer class of the final fused call. "
        "**c**, Per genome, on its own y axis: metal calls (final call a metal) over discovered TFs, "
        "ion-resolved calls (naming one species, stacked by ion) over metal calls, and ion-resolved calls "
        f"with one of the top {k} operators in the regulator's own flanking intergenic region. "
        "**d**, Accuracy per axis over the TFs with truth on that axis. The primary operator (rank-1 "
        "rescan site) is correct when it overlaps a truth site within 25 bp (*E. coli*: classical "
        f"binding at RegulonDB confidence Strong or Confirmed): {eo['observed']}/{e['operator']['all']['n']} "
        f"in *E. coli* ({expected(eo['expected_by_chance'])} expected by chance, "
        f"P = {pval(eo['p_vs_chance'])}) and {so['observed']}/{s['operator']['all']['n']} in SL1344 "
        f"({expected(so['expected_by_chance'])} expected, P = {pval(so['p_vs_chance'])}); exact "
        "Poisson-binomial test against sites of the same widths placed uniformly in intergenic DNA. "
        f"**e**, Class confusion, both genomes pooled ({diag}/{len(conf)} correct). "
        "**f**, Every known metal sensor: found, called metal, exact ion (amber: the ion is in the "
        "shortlist but not first); *E. coli* sensors outside the studied families carry their Pfam "
        "family. " + caveat(N)
    )


def legend_fig3(N: dict) -> str:
    sv, per = N["survey"]["all"], N["survey"]["per_genome"]
    genomes = survey_genomes(N)
    nfam = len(N["families"])
    k = per[genomes[0]]["autoregulation"]["top_k"]
    f, op = sv["metal_funnel"], sv["operator_supports_ion_call"]
    counts = ", ".join(f"*{g}* {per[g]['n']}" for g in genomes)
    outside = ", ".join(sv["families_outside_study"])
    return (
        f"The {sv['n']} candidates discovered in the {w(len(genomes))} survey genomes ({counts}). No operator or "
        "inducer truth set exists for these genomes, so no panel here is an accuracy claim: each shows "
        "what the pipeline returned. "
        f"**a**, Candidates per family, by the pipeline's own family call: the {w(nfam)} studied families "
        "in the fixed study order every figure shares (an order of presentation, not a biological one), "
        f"then, in tan, the other families that candidates fall in ({outside}); a blank cell means no "
        "candidate, and each row label carries its genome's candidate count. "
        "**b**, Inducer class of the final fused call, as a percentage of each genome's candidates, with "
        "the count inside each segment. "
        "**c**, As Figure 2c, one panel per genome on its own y axis: metal calls over candidates, "
        "ion-resolved calls over metal calls (stacked and labelled by ion), and ion-resolved calls whose "
        f"regulator has one of its top {k} operators in its own flanking intergenic region. Pooled over "
        f"the {w(len(genomes))} genomes these are {f['metal_calls']}/{sv['n']} metal calls, "
        f"{f['ion_resolved']}/{f['metal_calls']} of them ion-resolved, and "
        f"{op['n']}/{op['of_ion_resolved']} with operator support. Because the rescan's locality prior "
        "favours the regulator's own region, that last column is agreement between the inducer path and "
        "the operator path, not independent proof."
    )


def legend_fig4(N: dict) -> str:
    return (
        "CueR, Fur (fcrX) and DtxR/IdeR (EX350_RS17345) are the three displayed predictions whose "
        "AlphaFold 3 dimers engage their predicted operators. **Genomic context:** annotated genes "
        "are directional arrows drawn to scale; each regulator's arrow and metal badge match its "
        "protein colour, and the badge names the called effector. The green mark is the rank-1 "
        "operator; q is the operator "
        "call's adjusted value. Coordinates and a 1-kb scale are shown. "
        "**Operator complex:** predicted dimer bound to a 60-bp duplex. Each predicted operator "
        "was extended to 60 bp with its bordering genomic sequence from the same region; the "
        "operator is shown in pale green and marked by a green bracket, and gold marks contacted "
        "nucleotides. The Fur and DtxR DNA complexes include modeled metal; the CueR complex is "
        "apo. **Holo dimer and metal site:** each regulator has a distinct protein palette, with "
        "its two protomers in dark and light tints and modeled ions as grey spheres. The rightmost view enlarges "
        "one site; coordinating side chains are sticks, with dark dotted links for donor atoms "
        "within 3.0 Å and lighter dashed links for 3.0–4.3 Å. AlphaFold 3 used Cu(II) in place "
        "of the called Cu(I) for CueR and Mn(II) in place of the called Fe(II) for Fur and DtxR; "
        "Fur also contains structural Zn(II). The structures and contacts are predictions, "
        "not experimental evidence of binding or metal specificity."
    )


def legend_figS1(N: dict) -> str:
    from build_benchmark_figures import has_truth_beyond_family
    b = N["benchmark"]
    eco, sal = b["ecoli"]["per_tf"], b["salmonella"]["per_tf"]
    n_e = sum(1 for r in eco if has_truth_beyond_family(r))
    n_s = sum(1 for r in sal if has_truth_beyond_family(r))
    n_osman = b["metal_sensor_accounting"]["salmonella_recovery"]["sensors"]
    return (
        "One row per discovered benchmark TF that carries truth beyond its family — "
        f"{n_e} of the {len(eco)} discovered in *E. coli* K-12 and {n_s} of the {len(sal)} in SL1344 "
        "— with one mark per axis: family, inducer class, metal ion and primary operator. A mark is "
        "correct, partial (the right class, with the ion in the shortlist but not first), incorrect, or "
        "not scored where that TF has no truth on that axis, as in the key. Names in bold are metal "
        "sensors; a trailing asterisk marks a TF whose own sequence anchors a shipped SSN clade, which is "
        "curated support and not an independent prediction. Bottom right, the remaining TFs, which carry "
        f"a family call and nothing else ({len(eco) - n_e} in *E. coli* K-12, {len(sal) - n_s} in "
        "SL1344), listed by name with their family tally; SL1344 candidates with no census entry "
        "have no family truth either and are listed separately as not scored. For SL1344, class, ion and "
        f"operator truth exist only for the {n_osman} Osman 2019 sensors."
    )


def legend_figS2(N: dict) -> str:
    per = N["survey"]["per_genome"]
    genomes = survey_genomes(N)
    k = per[genomes[0]]["autoregulation"]["top_k"]
    return (
        f"Per genome, the percentage of candidates for which one of the regulator's top {k} operators "
        "lies in its own flanking intergenic region; n/N above each bar is those candidates over every "
        f"candidate discovered in that genome. The {w(len(genomes))} survey genomes are blue and the "
        f"{w(len(bench_genomes(N)))} benchmark genomes, right of the dashed line, grey. The rescan's locality prior ranks the "
        "regulator's own region first, so this rate is in part a property of the ranking rather than a "
        "measurement of autoregulation; it is not an accuracy claim."
    )


def legend_figS3(N: dict) -> str:
    import statistics

    from build_supplementary_figures import genome_rows
    rows = genome_rows(N)
    per_tf = [t for _, c in rows for t in c["pwm_quality"]["per_tf"]]
    tot = sum(c["pwm_quality"]["n"] for _, c in rows)
    one = sum(c["pwm_quality"]["n_effective_counts"].get("1", 0) for _, c in rows)
    med = statistics.median(t["palindromicity"] for t in per_tf)
    n_fam = len({t["family"] for t in per_tf if t["family"]})
    return (
        f"The {tot} seed PWMs of all {w(n_all_genomes(N))} genomes, pooled or side by side; the survey genomes have no "
        "truth set, so no panel is an accuracy claim. "
        "**a**, Median palindromicity per family — the Pearson correlation between a PWM and its "
        "reverse complement, 1 a perfect inverted repeat and 0 no dyad symmetry — over the "
        f"{n_fam} families in which a PWM was built, with n beside each bar. Blue marks the families whose "
        "studied clades are metal-sensing by construction, a property of the study design rather than a "
        f"result; the dashed line is the overall median, r = {med:.2f}. "
        "**b**, The effective number of sequences behind each seed PWM, as a percentage of each genome's "
        f"PWMs. Pooled over the {w(n_all_genomes(N))} genomes, {one} of {tot} PWMs rest on a single sequence. "
        "**c**, Palindromicity against effective sequences, one point per PWM (jittered horizontally), "
        "the bar the median and n above each column."
    )


def legend_figS4(N: dict) -> str:
    from build_supplementary_figures import STATUS, genome_rows
    rows = genome_rows(N)
    chans = _channel_names(N)
    st = {key: sum(c["ssn_contribution"]["status"][key] for _, c in rows) for key, _, _ in STATUS}
    by = sum(c["ssn_contribution"]["ion_named_by_clade"] for _, c in rows)
    other = sum(c["ssn_contribution"]["ion_named_otherwise"] for _, c in rows)
    return (
        f"Every metal call of the {w(n_all_genomes(N))} genomes ({sum(st.values())} in all); a metal call is a TF whose "
        "final fused inducer call is a metal. "
        "**a**, What the SSN clade said for each of a genome's metal calls: the clade agreed with at "
        f"least one other channel ({st['clade_and_other']} pooled), carried the call alone "
        f"({st['clade_only']}), stayed silent while the other channels carried it ({st['other_only']}), "
        f"or itself called non-metal ({st['clade_disagrees']}). The total is printed above each bar. "
        f"**b**, Of the ion-resolved calls, whether the ion was named by the clade ({by} pooled) or by "
        f"another channel ({other}). "
        f"**c**, Which channels called metal, pooled over the {w(n_all_genomes(N))} genomes, as an UpSet-style "
        "matrix: each "
        "bar counts the metal calls made by exactly the combination of channels marked beneath it, out of "
        f"the {w(len(chans))} — {', '.join(chans[:-1])} and {chans[-1]} — and a blue bar marks "
        "the combinations that include the clade."
    )


def _survey_table_columns(N: dict, k: int) -> str:
    """The column definitions Figures S5 and S5bis (1/3) share. The table's own footnote states the
    regulated-gene and own-region rules in full, so the legend names them and points there."""
    has_shortlist = any(r["shortlist"] for r in N["survey_table"]["rows"] if r["metal_call"])
    unresolved = ("followed by the clade's candidate ions where it has a shortlist and otherwise by the "
                  "channels that called metal" if has_shortlist else
                  "followed by the channels that called metal (no clade shortlist exists for any such call "
                  "in this run)")
    return (
        "Columns: the TF (gene name in italic, otherwise the tail of its locus tag); the family and SSN "
        "clade the pipeline assigned (– where no clade matched); the inducer called, as a coloured square "
        f"and its ion, or 'unresolved' for a metal call naming no single ion, {unresolved}; the primary "
        f"operator, 5′→3′ on its own strand; whether one of the top {k} operators lies in the regulator's "
        "own flanking intergenic region, with its rank, scored only for calls that name one ion "
        "(· otherwise); and the regulated gene, with its product in grey where the gene has no name. The "
        "rules for the regulated gene and the own-region operator are stated under the table."
    )


def legend_figS5(N: dict) -> str:
    T = N["survey_table"]
    return (
        f"The {T['n_metal_calls']} metal calls among the {T['n']} survey candidates, one block per genome "
        f"in the order every figure uses ({_species(survey_genomes(N))}) and, within a block, by the "
        "study's family order and then by name; each block header repeats its genome's metal calls and "
        "candidates. " + _survey_table_columns(N, T["top_k"]) +
        " The survey genomes have no truth set, so nothing in this table is an accuracy claim."
    )


def legend_figS5bis(N: dict, part: int) -> str:
    """One legend per S5bis page; `part` is 1, 2 or 3 (build_survey_table.BIS_PAGES)."""
    from build_survey_table import BIS_PAGES
    T = N["survey_table"]
    genomes = survey_genomes(N)
    here = [g for g in BIS_PAGES[part - 1] if g in genomes]
    n_here = [sum(r["genome"] == g for r in T["rows"]) for g in here]
    counts = (f"{n_here[0]} candidates" if len(here) == 1 else
              ", ".join(str(n) for n in n_here[:-1]) + f" and {n_here[-1]} candidates")
    head = (f"Every one of the {T['n']} survey candidates, not only the metal calls, split over "
            f"{w(len(BIS_PAGES))} pages by genome; this page carries {_species(here)} ({counts}), each "
            "block headed by its genome's metal calls and candidates and ordered by the study's family "
            "order, then by name. ")
    tail = "The survey genomes have no truth set, so nothing in this table is an accuracy claim."
    if part != 1:
        return (head + f"Columns, marks and rules as in Figure S5bis (1/{len(BIS_PAGES)}), including the "
                "inducer-class square and the rule for the regulated gene. " + tail)
    return (head + "One column is added to Figure S5: the inducer class of the final call, as a coloured "
            "square. 'non-metal, no ligand' is the pipeline's 'non-metal effector (undetermined)' call, "
            "which the inducer vocabulary classes as organic. "
            + _survey_table_columns(N, T["top_k"]) + " " + tail)


#: file stem -> legend function
LEGENDS = {
    "fig1_pipeline": legend_fig1,
    "fig2_benchmark": legend_fig2,
    "fig3_survey": legend_fig3,
    "fig4_context_structures": legend_fig4,
    "figS1_per_regulator": legend_figS1,
    "figS2_own_promoter": legend_figS2,
    "figS3_pwm_quality": legend_figS3,
    "figS4_ssn_clade": legend_figS4,
    "figS5_survey_metal_table": legend_figS5,
    "figS5bis_survey_table_1": lambda N: legend_figS5bis(N, 1),
    "figS5bis_survey_table_2": lambda N: legend_figS5bis(N, 2),
    "figS5bis_survey_table_3": lambda N: legend_figS5bis(N, 3),
}


def figure_entries(N: dict | None = None) -> list[dict]:
    if N is None:
        N = json.loads((HERE / "numbers.json").read_text(encoding="utf-8"))
    return [{"key": k, "label": lab, "title": t, "main": m, "legend": LEGENDS[k](N)}
            for k, lab, t, m in FIGURES]


def main() -> int:
    N = json.loads((HERE / "numbers.json").read_text(encoding="utf-8"))
    L = ["# Figure legends", "", f"Generated from `numbers.json`. {provenance_line(N)}", ""]
    for e in figure_entries(N):
        L += [f"## {e['label']} | {e['title']}", "", f"`figures/{e['key']}.*`", "", e["legend"], ""]
    (HERE / "LEGENDS.md").write_text("\n".join(L), encoding="utf-8")
    print("  wrote LEGENDS.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
