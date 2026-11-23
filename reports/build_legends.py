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
    ("fig1_pipeline", "Figure 1", "RIIO combines sequence and genomic context to predict operators and inducers", True),
    ("fig2_benchmark", "Figure 2", "Benchmark on E. coli K-12 and S. Typhimurium SL1344", True),
    ("fig3_survey", "Figure 3", "Four-genome survey", True),
    # Figure 4 is the SDP figure (build_fig_sdp_structures.py; Fur-family row V. cholerae Zur, the user's
    # choice of 2026-10-15). It replaced fig4_context_structures (CueR/Fur/IdeR, no SDPs).
    ("fig_sdp_representatives", "Figure 4",
     "Discovered metalloregulators coloured by specificity-determining positions", True),
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
#: rows per page of the all-candidates figure; build_fig_tried_candidates.ROWS_PER_PAGE, kept here so
#: this module does not import the builder (test_sdp_structures asserts the two agree)
ALL_CANDIDATES_PER_PAGE = 5


def all_candidates_figures() -> list[tuple[str, str, str, bool]]:
    """The supplementary pages that show every modelled regulator, read from the figure's own layout
    (reports/structures/tried/regulators_sdp.json) so the labels cannot drift from the pages."""
    path = HERE / "structures" / "tried" / "regulators_sdp.json"
    if not path.exists():
        return []
    rows = [(k, v) for k, v in json.loads(path.read_text(encoding="utf-8")).items() if not k.startswith("_")]
    pages = [rows[i:i + ALL_CANDIDATES_PER_PAGE] for i in range(0, len(rows), ALL_CANDIDATES_PER_PAGE)]
    n = len(pages) + 1                       # page 1 is the summary table
    out = [("fig_all_candidates_p1", f"Figure S6 (1/{n})",
            "Every regulator modelled for Figure 4, at a glance", False)]
    for i, page in enumerate(pages, 2):
        groups: list[str] = []
        for _, r in page:
            g = r["group"].split(" (")[0].split(",")[0]   # the block's name, without its gloss
            if g not in groups:
                groups.append(g)
        out.append((f"fig_all_candidates_p{i}", f"Figure S6 ({i}/{n})",
                    "Every regulator modelled for Figure 4: " + _and(groups), False))
    return out


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
    return (
        "Conceptual overview; drawings illustrate states rather than measured data. The TF (a CueR "
        "dimer), its sequence and its SSN clade (MerR cluster 7) are drawn from E. coli CueR as an "
        "example; motifs, scan scores and loci are illustrative. **a**, Four input modes are "
        "available through the CLI, "
        "Python API or browser GUI. Supplied proteins bypass discovery without a family filter. "
        "Genome routes obtain a proteome from annotation or gene calls, discover candidates with "
        "BITACORA against curated sequence-similarity-network (SSN) clade profiles "
        "(HMMER/BLAST+; default E ≤ 1e-5), and select TFs for prediction. Discovery may stop here. "
        "**b**, Both routes converge on TF sequence, family, SSN clade and available locus/flanking "
        "DNA. The teal operator path retrieves orthologue promoters, derives an inverted-repeat "
        "position-weight matrix (PWM), scans intergenic DNA (default per-site p ≤ 1e-4), and ranks "
        "candidate operators with regulon and mode annotations. The amber inducer path fuses "
        "curated SSN-clade labels, Cys/His motifs, optional MetalNet2 sequence/MSA site evidence, "
        "nearby-gene functions and operon-chemistry ligand evidence. Channels may abstain; ions are named "
        "only if resolved. Inducer-derived gene labels annotate operator hits without boosting "
        "their ranking. An apo fold and fold QC are optional (off by default). **c**, Both outputs "
        "support a regulatory hypothesis linking the TF, candidate operators/targets and putative "
        "inducer. Solid arrows indicate data flow; dashed links indicate unvalidated biological "
        "relationships. Per-TF bundles retain predictions and evidence; run and discovery reports "
        "and Supplementary Table S1 are rebuilt from saved data. AlphaFold 3 apo, DNA and "
        "metal-context input jobs are exported for server/local submission, never executed. "
        "Single-protein calls enter prediction directly; list/genome runs check input and policy "
        "fingerprints to reuse bundles and retry failures. Outputs require experimental validation."
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
        "**b**, Inducer class of the final fused call. " + site_only_sentence(
            [N["benchmark_genome_composition"][g] for g in bench_genomes(N)]) +
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
        "**f**, Every known metal sensor: found, called metal, exact ion (amber: the call is a shortlist "
        "of ions that holds the known one without ranking it first; the shortlist's members are written "
        "in the note); *E. coli* sensors outside the studied families carry their Pfam "
        "family. " + caveat(N)
    )


def _strains(N: dict) -> str:
    """One sentence naming the survey strains (numbers.json survey.strains), or '' if absent."""
    st = N["survey"].get("strains") or {}
    if not st:
        return ""
    names = [f"*{st[g]['organism']}* ({st[g]['accession']})" for g in survey_genomes(N) if g in st]
    note = (" The *V. cholerae* genome is strain N16961, not the RFB16 strain of the published survey; its "
            "candidates were discovered afresh on N16961." if any("N16961" in v["organism"] for v in st.values()) else "")
    return "Genomes: " + _and(names) + "." + note + " "


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
        f"{_strains(N)}"
        f"**a**, Candidates per family, by the pipeline's own family call: the {w(nfam)} studied families "
        "in the fixed study order every figure shares (an order of presentation, not a biological one), "
        f"then, in tan, the other families that candidates fall in ({outside}); a blank cell means no "
        "candidate, and each row label carries its genome's candidate count. "
        "**b**, Inducer class of the final fused call, as a percentage of each genome's candidates, with "
        "the count inside each segment. " + site_only_sentence([sv]) +
        "**c**, As Figure 2c, one panel per genome on its own y axis: metal calls over candidates, "
        "ion-resolved calls over metal calls (stacked and labelled by ion), and ion-resolved calls whose "
        f"regulator has one of its top {k} operators in its own flanking intergenic region. Pooled over "
        f"the {w(len(genomes))} genomes these are {f['metal_calls']}/{sv['n']} metal calls, "
        f"{f['ion_resolved']}/{f['metal_calls']} of them ion-resolved, and "
        f"{op['n']}/{op['of_ion_resolved']} with operator support. Because the rescan's locality prior "
        "favours the regulator's own region, that last column is agreement between the inducer path and "
        "the operator path, not independent proof."
    )


# ---------------------------------------------------------------- Figure 4: SDP representatives
# numbers.json "sdp" rows: the three drawn (figure_row) and V. cholerae Fur, the Fur-family candidate
# the drawn Zur row was chosen over, kept so the document can state that comparison.
#: the protein names the text uses for each row (labels, not results)
SDP_NAME = {"smtB": "SmtB", "RS09920": "EX350_RS09920", "RS12295": "EX350_RS12295", "fcrX": "Fur",
            "zur": "Zur", "RS08730": "Zur"}
_AA3 = {"A": "Ala", "C": "Cys", "D": "Asp", "E": "Glu", "F": "Phe", "G": "Gly", "H": "His", "I": "Ile",
        "K": "Lys", "L": "Leu", "M": "Met", "N": "Asn", "P": "Pro", "Q": "Gln", "R": "Arg", "S": "Ser",
        "T": "Thr", "V": "Val", "W": "Trp", "Y": "Tyr"}
_ROMAN = {"1": "I", "2": "II", "3": "III"}


def _ox(ion: str) -> str:
    """'Cd2+' -> 'Cd(II)'."""
    m = re.fullmatch(r"([A-Z][a-z]?)(\d)\+", ion)
    return f"{m[1]}({_ROMAN[m[2]]})" if m else ion


def _residues(t: dict) -> str:
    """'D116:F H118:F' -> 'Asp116/His118' (a test's residue set, in sequence order)."""
    return "/".join(_AA3[x[0]] + x[1:x.index(":")] for x in t["residues"].split())


def _site(r: dict, name: str) -> str:
    """A registry site name as text: 'alpha5' -> 'α5 site (Asp116/...)', "Cys77'/... (CadR-type)" ->
    'Cys77′/... site (CadR-type)'; names that already say 'site' stay as they are."""
    s = name.replace("alpha", "α").replace("'", "′")
    if "site" not in s:
        s = s.replace(" (", " site (", 1) if " (" in s else s + " site"
    if not re.search(r"\d\d", s):
        s += f" ({_residues(r['tests']['H1 functional'])})"
    return s


def _pooled(N: dict, hyp: str, cls: str, docked: bool = False) -> dict:
    subset = f"pDNA >= {N['sdp']['cuts']['docked_iptm']}" if docked else "all"
    return next(p for p in N["sdp"]["pooled"]
                if (p["hypothesis"], p["sdp_class"], p["subset"]) == (hyp, cls, subset))


def _frac(t: dict) -> str:
    return f"{t['in_class']}/{t['n']}"


def _and(items: list[str]) -> str:
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + " and " + items[-1]


def _sdp_rows(N: dict) -> tuple[list[str], str, list[str]]:
    """(the drawn rows in manuscript order, the drawn Fur-family row, the Fur-family rows not drawn)."""
    R = N["sdp"]["rows"]
    rows = [k for k, r in R.items() if r["figure_row"]]
    fur = next(k for k in rows if R[k]["family"] == "Fur")
    return rows, fur, [k for k, r in R.items() if not r["figure_row"] and r["family"] == "Fur"]


def _uncharacterised(R: dict, rows: list[str]) -> bool:
    """True when no drawn row has experimental characterisation of its own (the registry's record)."""
    return all(R[k].get("characterised") is False for k in rows)


def _char_note(r: dict) -> str:
    return r["characterisation_note"].split(" (UniProt and PubMed")[0]


def _who(R: dict, k: str) -> str:
    r = R[k]
    return f"*{r['genome']}* {SDP_NAME[k]} ({r['family'].split('/')[0]} cluster {r['sdp_cluster']})"


def legend_fig4_sdp(N: dict) -> str:
    S = N["sdp"]
    R, cut = S["rows"], S["cuts"]
    rows, fur, _ = _sdp_rows(N)
    op = R[fur]["operator"]
    h1 = _pooled(N, "H1", "functional")

    if op["own_region"]:
        where = (f"{SDP_NAME[fur]}'s is its rank-{op['rank']} site, in its own promoter region"
                 + (f", shared with *{op['gene']}*." if op["gene"] else "."))
    elif op["rank"] == 1:
        where = (f"{SDP_NAME[fur]}'s is its top-ranked site, in the *{op['gene']}* promoter; no site in its "
                 "own promoter region is among its top five.")
    else:
        where = (f"{SDP_NAME[fur]}'s is the rank-{op['rank']} site, in the *{op['gene']}* promoter: no model "
                 "docked on the rank-1 site.")
    def zoom(k: str) -> str:
        sites = [_site(R[k], s) for s, role in R[k]["sites"].items()
                 if role in ("regulatory", "structural", "secondary") and (k != "smtB" or role == "regulatory")
                 and not s.startswith("Cys site")]
        crystal = ", as in this protein's crystal structure (PDB 2W57)" if "MN" in R[k]["standin"] else ""
        return f"{SDP_NAME[k]}, {_and(sites)}{crystal}"

    standins = [f"{_ox(el.title() + '2+')} for {_ox(ion)} in {SDP_NAME[k]}"
                for k in rows for el, ion in R[k]["standin"].items()]
    structural = any("H2 structural" in R[k]["tests"] for k in rows)
    return (
        f"One discovered regulator per family, coloured by the specificity-determining positions (SDPs) of "
        f"its SSN cluster: {_and([_who(R, k) for k in rows])}. "
        "**Context:** the green box is the pipeline's predicted operator used for the DNA complex (rank, "
        f"q value); no curated site was used. {where} "
        "**AlphaFold 3 panels,** left to right: DNA-bound complex (60-bp duplex; operator bracketed, "
        "contacted nucleotides gold) with its protein–DNA chain-pair ipTM, red below "
        f"{cut['docked_iptm']}; DNA interface, side chains within {cut['dna_contact_A']} Å of the DNA "
        "labelled; holo dimer in the complex's pose; metal-site zooms "
        f"({'; '.join(zoom(k) for k in rows)}). Metal–ligand bonds ≤ {cut['coordination_A']} Å"
        + (f"; AlphaFold 3 used {_and(standins)}" if standins else "") + ". "
        "**Colours:** as in the manuscript's Figs 4C and 5C, functional SDPs red and structural SDPs "
        "blue by relative score; other aligned residues near-white, unaligned grey. "
        "**Bottom:** SDP class of the regulatory-metal ligands, their second shell (within "
        f"{cut['second_shell_A']:g} Å of the ion), "
        + ("the structural Zn ligands, " if structural else "")
        + f"DNA base contacts (hatched below ipTM {cut['docked_iptm']}) and all residues. Over "
        f"{w(len(h1['regulators']))} modelled regulators, {h1['in_class']} of {h1['n']} regulatory-metal "
        f"ligands are functional SDPs ({h1['expected']:.0f} expected; circular-shift P = "
        f"{pval(h1['p_rotation'])}); as AlphaFold 3 and the SDPs share evolutionary signal, this is "
        "co-localisation, not independent validation."
    )


def legend_all_candidates(N: dict, key: str) -> str:
    """One legend per page of the all-candidates supplement: the summary table, then the row pages."""
    S, dk = N["sdp"], N["sdp"]["cuts"]["docked_iptm"]
    cfg = json.loads((HERE / "structures" / "tried" / "regulators_sdp.json").read_text(encoding="utf-8"))
    rows = [(k, v) for k, v in cfg.items() if not k.startswith("_")]
    docked = sum(v["protein_dna_iptm"] >= dk for _, v in rows)
    h1 = _pooled(N, "H1", "functional")
    common = (
        f"Every one of the {len(rows)} regulators modelled for Figure 4 that was folded on an operator the "
        "pipeline predicted, grouped by the transcription factor it is. Three further regulators "
        "(M. tuberculosis NmtR and the two ZntRs) were folded only on a curated box and are left out, "
        "because the figure's operators must be predictions. Colours are those of Figure 4: functional "
        "SDPs red, structural SDPs blue, intensity the relative score; the DtxR and CsoR families have no "
        "SDP table and are drawn in neutral grey.")
    if key.endswith("_p1"):
        return (
            "**The summary table.** " + common + " For each regulator the table gives its SSN cluster and "
            "identity to that cluster's consensus, the called inducer, which predicted operator the DNA "
            "complex used and whether that operator lies in the regulator's own promoter region, the "
            f"protein–DNA chain-pair ipTM (bold where the complex docks at ≥ {dk}; {docked} of "
            f"{len(rows)} do), the residues the model places at the regulatory metal site, and how many "
            "of those are functional SDPs of its cluster. The pages that follow draw each row in full. "
            f"Pooled over these regulators, {h1['in_class']} of {h1['n']} regulatory-metal ligands are "
            f"functional SDPs ({h1['expected']:.0f} expected; circular-shift P = {pval(h1['p_rotation'])}).")
    page = int(key.rsplit("_p", 1)[1])
    per = ALL_CANDIDATES_PER_PAGE
    on_page = rows[(page - 2) * per:(page - 1) * per]
    # named from the run Figure 4 was built on (numbers.json sdp.names), not the canonical survey,
    # which no longer holds the V. cholerae RFB16 regulators
    survey = N["sdp"].get("names") or {s["run_name"]: s for s in N["survey_table"]["rows"]}
    names = _and([f"*{v['genome']}* {survey[v['run_name']]['tf']} "
                  f"({survey[v['run_name']]['locus_tag']})" for _, v in on_page])
    return (
        f"**{names}.** " + common + " Each row carries the regulator's status (drawn in Figure 4, or "
        "folded and not drawn), how many AlphaFold 3 jobs it has and the best protein–DNA ipTM any of "
        "them reached on a predicted operator, then its genomic context with that operator, the DNA "
        "complex in the form that binds DNA (apo for ArsR/SmtB, metal-loaded for MerR and the Fur "
        "family), one protomer's DNA interface, the holo dimer in the complex's pose, and its metal "
        f"sites. Metal–ligand bonds are drawn to {S['cuts']['coordination_A']} Å and side chains within "
        f"{S['cuts']['dna_contact_A']} Å of the DNA are labelled whatever the ipTM.")


def sdp_summary(N: dict) -> dict:
    """The run document's Figure 4 section: what the SDP mapping shows, and why the Fur-family row is
    the regulator it is. {"heading", "intro", "bullets"}; every number from numbers.json "sdp"."""
    S = N["sdp"]
    R, cut = S["rows"], S["cuts"]
    dk = cut["docked_iptm"]
    rows, fur, others = _sdp_rows(N)
    h1, h1s, h2 = _pooled(N, "H1", "functional"), _pooled(N, "H1s", "functional"), _pooled(N, "H2", "structural")
    h3, h3d = _pooled(N, "H3", "structural"), _pooled(N, "H3", "structural", docked=True)
    h3xd = _pooled(N, "H3x", "structural", docked=True)
    docked = sorted(set(h3d["regulators"]) | set(h3xd["regulators"]))
    t1 = {k: R[k]["tests"]["H1 functional"] for k in R}
    in_fig = [SDP_NAME[k] for k in rows if R[k]["protein_dna_iptm"] >= dk]
    fur_family = [k for k in R if R[k]["family"] == "Fur"]
    survey = {k: next(x for x in N["survey_table"]["rows"] if x["run_name"] == R[k]["run_name"]) for k in R}

    def site(k: str) -> str:
        o = R[k]["operator"]
        return (f"rank-{o['rank']} operator "
                + ("in its own promoter region" if o["own_region"] else f"in the *{o['gene']}* promoter"))

    shared = len({site(k) for k in rows}) == 1

    def row(k: str) -> str:
        r = R[k]
        return (f"*{r['genome']}* {SDP_NAME[k]} ({r['family'].split('/')[0]} cluster {r['sdp_cluster']}; "
                f"protein–DNA ipTM {r['protein_dna_iptm']:.2f}" + ("" if shared else f"; {site(k)}") + ")")

    def strengths(k: str, o: str) -> list[str]:
        """What row k has over the other Fur-family candidate o, from the numbers above."""
        r, q = R[k], R[o]
        ag = {x: survey[x]["agreement"] for x in (k, o)}
        out = []
        if r["protein_dna_iptm"] > q["protein_dna_iptm"]:
            lead = ("the only Fur-family complex that docks" if r["protein_dna_iptm"] >= dk > q["protein_dna_iptm"]
                    else "the better-docked DNA complex" if q["protein_dna_iptm"] >= dk
                    else "the less weak DNA complex")
            out.append(f"{lead} ({r['protein_dna_iptm']:.2f} against {q['protein_dna_iptm']:.2f})")
        if (r["operator"]["rank"], -r["operator"]["score"]) < (q["operator"]["rank"], -q["operator"]["score"]):
            # q is quoted only when the regulon reports one for both sites
            qs = (f"; q = {pval(r['operator']['qvalue'])} against {pval(q['operator']['qvalue'])}"
                  if r["operator"]["qvalue"] is not None and q["operator"]["qvalue"] is not None else "")
            out.append(f"the stronger operator (rank {r['operator']['rank']}"
                       + (", in its own promoter region" if r["operator"]["own_region"] else "") + qs + ")")
        elif r["operator"]["own_region"] and not q["operator"]["own_region"]:
            out.append("an operator in its own promoter region")   # only when the clause above is absent
        if ag[k] > ag[o]:
            out.append(f"the more corroborated metal call (agreement {ag[k]:.2f} against {ag[o]:.2f})")
        if r["modeled"] == r["called"] and q["modeled"] != q["called"]:
            out.append(f"its native metal in the model ({_ox(r['called'])}, where {SDP_NAME[o]} needed a "
                       f"{_ox(q['modeled'])} stand-in for {_ox(q['called'])})")
        if "MN" in r["standin"]:
            out.append("a crystal structure of the same protein whose sites the model reproduces (2W57)")
        if "H2 structural" in r["tests"] and "H2 structural" not in q["tests"]:
            out.append("both Fur-family metal sites of the manuscript's Fig. 5C, with the structural one "
                       "confirmed by the SDPs")
        return out

    intro = [
        "Figure 4 shows one discovered regulator per family, coloured by the collaborator's "
        "specificity-determining positions (SDPs) for its SSN cluster. It replaces the earlier "
        "CueR/Fur/IdeR figure, whose DtxR-family regulator has no SDP table. Every DNA complex is modelled on "
        f"an operator the pipeline predicted; no curated site is used. The rows are {_and([row(k) for k in rows])}"
        + (f", each modelled on its {site(rows[0]).replace(' in its own', ', in its own')}." if shared else ".")
        + ((" The rows were chosen among the regulators whose complex docks (protein–DNA ipTM ≥ "
            f"{dk}) on an operator this run predicts, taking in each family one with no experimental "
            "characterisation of its own, so that the figure shows what the pipeline proposes for "
            "proteins nobody has studied: " + "; ".join(f"{SDP_NAME[k]}, {_char_note(R[k])}" for k in rows)
            + ".") if _uncharacterised(R, rows) else ""),
    ]
    bullets = [
        f"**Metal sites.** Across {w(len(h1['regulators']))} modelled regulators of the three families, "
        f"{h1['in_class']} of {h1['n']} residues coordinating the regulatory metal are functional SDPs, "
        f"against {h1['expected']:.0f} expected by chance (exact circular-shift P = {pval(h1['p_rotation'])}). "
        "In the figure rows: " + "; ".join(
            f"{SDP_NAME[k]} {_frac(t1[k])} ({_residues(t1[k])}; P = {pval(t1[k]['p_rotation'])})"
            for k in rows) + ". "
        f"In {SDP_NAME['RS12295']} Cys77 is functional, while Cys126, like the two cysteines of CueR, is "
        "conserved across the metal-sensing MerR clusters and so classed structural.",
        f"**Second shell.** Non-ligand residues within {cut['second_shell_A']:g} Å of the regulatory ion "
        f"are also enriched for functional SDPs ({h1s['in_class']}/{h1s['n']} against "
        f"{h1s['expected']:.0f} expected; P = {pval(h1s['p_rotation'])}). This test was added after the "
        "first folds were seen, so it is exploratory.",
        f"**Structural site.** Every ligand of the structural Zn site is a structural SDP, in the figure's "
        f"{SDP_NAME[fur]} ({_frac(R[fur]['tests']['H2 structural'])}: "
        f"{_residues(R[fur]['tests']['H2 structural'])}) and in each of the "
        f"{w(len(h2['regulators']))} Fur-family regulators that form one "
        f"({h2['in_class']}/{h2['n']} pooled, P = {pval(h2['p_rotation'])}), as the manuscript's Fig. 5C "
        "shows for scZur. Two of those regulators are in Fur clusters 4 and 9, which Fig. 5C does not "
        "draw, so the pattern holds outside the clusters it was read from.",
        f"**DNA interface.** Only {w(len(docked))} complexes dock at protein–DNA ipTM ≥ {dk}, "
        f"{w(len(in_fig))} of them in the figure ({_and(in_fig)})"
        + ("; no Fur-family complex does" if all(R[k]["protein_dna_iptm"] < dk for k in fur_family) else "")
        + f". In the {w(len(h3d['regulators']))} docked MerR and Fur-family complexes {_frac(h3d)} base "
        f"contacts are structural SDPs (P = {pval(h3d['p_rotation'])}); in the "
        f"{w(len(h3xd['regulators']))} docked ArsR complexes {_frac(h3xd)} are "
        f"(P = {pval(h3xd['p_rotation'])}). "
        + ("For MerR and the Fur family this agrees with a family-wide DNA-recognition role for the "
           f"structural SDPs, on {w(len(h3d['regulators']))} complexes only; the ArsR complexes do not "
           "show it. " if h3d["p_rotation"] < 0.05 <= h3xd["p_rotation"] else
           "The models neither support nor refute a family-wide DNA-recognition role for the structural "
           "SDPs. ")
        + f"The count over all MerR and Fur-family complexes ({_frac(h3)}, P = {pval(h3['p_rotation'])}) "
        "rests mostly on complexes below the cut, so it is not evidence.",
    ]
    arsr_alt = [k for k in R if not R[k]["figure_row"] and R[k]["family"].startswith("ArsR")]
    arsr = next((k for k in rows if R[k]["family"].startswith("ArsR")), None)
    if arsr and arsr_alt:
        a = arsr_alt[0]
        bullets.append(
            f"**ArsR row.** *{R[arsr]['genome']}* {SDP_NAME[arsr]} is drawn instead of "
            f"*{R[a]['genome']}* {SDP_NAME[a]}, of the same SDP cluster. Both dock on their rank-1 "
            f"own-promoter operator ({R[arsr]['protein_dna_iptm']:.2f} against "
            f"{R[a]['protein_dna_iptm']:.2f}), and {SDP_NAME[a]} is experimentally characterised while "
            f"{SDP_NAME[arsr]} is not: {_char_note(R[arsr])}. The row is therefore a prediction for an "
            f"unstudied ortholog, not a new kind of sensor. {SDP_NAME[a]} stays in the figure's "
            "configuration.")
    if others:
        def who_locus(k: str) -> str:
            return f"*{R[k]['genome']}* {SDP_NAME[k]} ({survey[k]['locus_tag']})"

        drawn = who_locus(fur)
        bullets.append(
            f"**Fur-family row.** Of the Fur-family regulators modelled on a predicted operator, "
            f"{drawn} is drawn. It has " + _and(strengths(fur, max(
                others, key=lambda o: R[o]["protein_dna_iptm"]))
                + (["no experimental characterisation of its own"] if R[fur].get("characterised") is False
                   else [])) + ". "
            + " ".join(
                (f"{who_locus(o)} has {_and(strengths(o, fur))}, and its complex docks too "
                 f"({R[o]['protein_dna_iptm']:.2f}), but it is experimentally characterised, so it was "
                 "not drawn." if R[o]["protein_dna_iptm"] >= dk and R[o].get("characterised") else
                 f"{who_locus(o)} has {_and(strengths(o, fur))}, but its complex reaches only "
                 f"{R[o]['protein_dna_iptm']:.2f}.") for o in others)
            + " Both alternatives stay in the figure's configuration, so the comparison can be redrawn.")
    bullets.append("**Reading.** AlphaFold 3 draws on the same evolutionary information as the SDPs, so this "
                   "co-localisation supports the SDP assignment rather than validating it independently.")
    return {"heading": "Figure 4: SDPs on the modelled regulators", "intro": intro, "bullets": bullets}


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


def site_only_sentence(comps: list[dict]) -> str:
    """The hatched segment of the inducer-class panel, with its count and families, from numbers.json."""
    n = sum(c.get("site_only", {}).get("organic", 0) for c in comps)
    fams = sorted({f for c in comps for f in c.get("site_only", {}).get("families", {})})
    if not n:
        return ""
    return (f"Hatched: organic calls in which MetalNet2 finds a conserved metal site ({n}, in "
            f"{', '.join(fams)}); the site is reported but not read as a sensing site (structural in GntR; "
            "overruled by a negative family coordination gate elsewhere), so the call stays organic. ")


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
            "which the inducer vocabulary classes as organic; 'non-metal:' followed by names is the same "
            "class-level call with the proposals it kept (the clade's label, then operon chemistry), none "
            "of which is resolved to one ligand; 'no call', class undetermined, is a regulator for which "
            "the pipeline made no inducer call at all. ‡ marks an organic call in which MetalNet2 finds "
            "a conserved metal site that is not read as a sensing site. "
            + _survey_table_columns(N, T["top_k"]) + " " + tail)


#: file stem -> legend function
LEGENDS = {
    "fig1_pipeline": legend_fig1,
    "fig2_benchmark": legend_fig2,
    "fig3_survey": legend_fig3,
    "fig_sdp_representatives": legend_fig4_sdp,
    **{k: (lambda N, p=k: legend_all_candidates(N, p)) for k, _, _, _ in all_candidates_figures()},
    "figS1_per_regulator": legend_figS1,
    "figS2_own_promoter": legend_figS2,
    "figS3_pwm_quality": legend_figS3,
    "figS4_ssn_clade": legend_figS4,
    "figS5_survey_metal_table": legend_figS5,
    "figS5bis_survey_table_1": lambda N: legend_figS5bis(N, 1),
    "figS5bis_survey_table_2": lambda N: legend_figS5bis(N, 2),
    "figS5bis_survey_table_3": lambda N: legend_figS5bis(N, 3),
}
#: the all-candidates supplement is appended last, after _and() and the legend function it uses exist
FIGURES += all_candidates_figures()


def figure_entries(N: dict | None = None) -> list[dict]:
    if N is None:
        N = json.loads((HERE / "numbers.json").read_text(encoding="utf-8"))
    return [{"key": k, "label": lab, "title": t, "main": m, "legend": LEGENDS[k](N)}
            for k, lab, t, m in FIGURES]


def sdp_provenance_line(N: dict) -> str:
    """Which run Figure 4's rows are read from: the canonical survey run, except regulators it does
    not hold (the V. cholerae RFB16 alternatives), which stay on the run they were folded for."""
    rows = N.get("sdp", {}).get("rows", {})
    if not rows:
        return ""
    return (" Figure 4 and Figure S6 are read from this run. Most of their AlphaFold 3 folds were made for "
            "operators predicted in an earlier run and are used because this run predicts the same operator "
            "at the same rank; the three *V. cholerae* regulators (now strain N16961) and EX350_RS11325 "
            "were refolded on this run's rank-1 operators.")


def main() -> int:
    N = json.loads((HERE / "numbers.json").read_text(encoding="utf-8"))
    L = ["# Figure legends", "", (f"Data-bearing figures use `numbers.json`. {provenance_line(N)} "
         "Figure 1 is a schematic of the v1.1 interfaces; no run data are plotted."
         + sdp_provenance_line(N)), ""]
    for e in figure_entries(N):
        L += [f"## {e['label']} | {e['title']}", "", f"`figures/{e['key']}.*`", "", e["legend"], ""]
    (HERE / "LEGENDS.md").write_text("\n".join(L), encoding="utf-8")
    print("  wrote LEGENDS.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
