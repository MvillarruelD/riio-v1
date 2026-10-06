"""
substrate_map.py -- transporter / efflux / resistance-enzyme  ->  cognate effector (ion or chemical).

Why this exists
---------------
The Ligify operon-chemistry route (`effector.ligify` + `effector.ligand`) infers the cognate ligand from
the *catalysed reactions* of neighbouring **enzymes**. That route is blind to the single most informative
gene class for metal-sensing regulators: the **transporters, efflux pumps and metallochaperones** the TF
co-regulates. A P1B-type ATPase, a CDF antiporter or an RND heavy-metal efflux complex carries no
small-molecule reaction to mine -- yet its *substrate specificity* names the inducer outright (a ZntA
exports Zn/Cd/Pb; a CopA exports Cu/Ag; an ArsB exports arsenite).

Gene Ontology cannot supply this: its molecular-function terms collapse to `metal ion binding`
(GO:0046872) / `cation transmembrane transport` -- true but ion-agnostic, and mostly IEA-transferred from
the same homology the SSN already exploits (circular). The *substrate resolution* lives in the curated
transporter classifications instead:

  * **TCDB** transporter families (e.g. 3.A.3.5 P-type heavy-metal ATPases; 2.A.4 CDF; 2.A.6.1 HME-RND),
  * **EC translocase** numbers (EC 7.2.2.- ion-translocating ATPases) and resistance-enzyme EC numbers,
  * **KEGG BRITE** ko02000 (transporters) / ko01504 (antimicrobial-resistance genes) and MetaCyc pathways.

This module distils those into an auditable rule table: each rule maps a transporter/enzyme class to the
ion/chemical it handles, matched by gene-name regex (most specific), EC number, or TCDB family (broadest).
`classify_gene()` returns the most specific `SubstrateHit`. It is consumed by `effector.regulon_enrichment`
(votes the inducer over a regulon / operon gene set) and by `report.ligand_regulon.ligand_genes` (flags
regulated-gene promoters at substrate resolution rather than the old hand keyword list).

Ion notation matches the rest of the pipeline (Zn2+, Cu+, Cd2+, Ni2+, Co2+, Hg2+, Pb2+, Mn2+, Fe, As3+,
Au+, CrO4). Run `python substrate_map.py` for a self-test.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

# canonical metal effectors (parity with effector.inducer / ssn_ligand_mine notation)
METAL_IONS = {"Zn2+", "Cu+", "Cd2+", "Ni2+", "Co2+", "Hg2+", "Pb2+", "Mn2+", "Fe", "As3+", "Au+", "CrO4"}


@dataclass(frozen=True)
class SubstrateRule:
    effector: str                 # canonical ion or chemical class the transporter/enzyme handles
    role: str                     # 'metal' | 'non-metal' | 'redox' | 'multidrug'
    substrate_class: str          # human-readable transporter/enzyme class (for provenance)
    specificity: int              # 3=gene-name (most specific) 2=EC 1=TCDB-family (broadest)
    ec: tuple = ()                # EC numbers (modern translocase 7.x AND legacy 3.6.3.x where relevant)
    tcdb: tuple = ()              # TCDB family prefixes (matched by startswith)
    name_patterns: tuple = ()     # regex over product + gene name (case-insensitive)
    source: str = ""              # provenance string (TCDB / EC / KEGG BRITE / MetaCyc)


@dataclass
class SubstrateHit:
    effector: str
    role: str
    substrate_class: str
    specificity: int
    matched_by: str               # 'name' | 'ec' | 'tcdb'
    detail: str                   # the token that matched (gene name fragment / EC / TC id)
    source: str
    gene: str = ""
    product: str = ""


# --------------------------------------------------------------------------- the rule table
# Ordered roughly specific->general within each metal; classify_gene scores by `specificity` so order is
# not load-bearing, but grouping by ion keeps it auditable. EC translocase numbers per IUBMB class 7.2.2
# (ion-translocating ATPases); legacy 3.6.3.x retained because many UniProt entries still carry them.
RULES: list[SubstrateRule] = [
    # ---- copper / silver (P1B-1 ATPase, multicopper oxidase, Cu chaperone, CusCBA RND) ----
    SubstrateRule("Cu+", "metal", "P1B-1-type Cu(+)/Ag(+)-translocating ATPase (CopA/GolT/CtpA)", 3,
                  ec=("7.2.2.8", "3.6.3.4"), tcdb=("3.A.3.5",),
                  name_patterns=(r"\bcopA\b", r"\bgolT\b", r"\bctpA\b",
                                 r"copper[- ]?(?:exporting|translocating|transporting).{0,12}ATPase",
                                 r"Cu\(\+?\).{0,12}ATPase"),
                  source="TCDB 3.A.3.5 / EC 7.2.2.8 / MetaCyc Cu efflux"),
    SubstrateRule("Cu+", "metal", "multicopper oxidase / cuprous oxidase (CueO/CopA-ox/CotA)", 3,
                  ec=("1.16.3.1", "1.10.3.-"),
                  name_patterns=(r"\bcueO\b", r"multicopper oxidase", r"cuprous oxidase",
                                 r"copper oxidase", r"\bcopA\b.*oxidase"),
                  source="EC 1.16.3.1 / MetaCyc cuprous-oxidation"),
    SubstrateRule("Cu+", "metal", "Cu(I) chaperone / Cu efflux RND (CopZ/Atx1/CusCBA/CusF)", 2,
                  tcdb=("2.A.6.1",),
                  name_patterns=(r"\bcopZ\b", r"\batx1\b", r"\bcusA\b", r"\bcusB\b", r"\bcusC\b",
                                 r"\bcusF\b", r"copper chaperone", r"copper.{0,8}efflux"),
                  source="TCDB 2.A.6.1 (HME-RND) / MetaCyc Cu chaperone"),
    SubstrateRule("Au+", "metal", "gold efflux determinant (GolS regulon: GolT/GolB)", 3,
                  name_patterns=(r"\bgolB\b", r"\bgolT\b", r"\bgesABC\b", r"gold[- ]resist"),
                  source="MetaCyc / Salmonella gol locus"),

    # ---- zinc / cadmium / lead (P1B-2 ATPase, CDF, ZnuABC, ZIP, CzcCBA) ----
    SubstrateRule("Zn2+", "metal", "P1B-2-type Zn(2+)/Cd(2+)/Pb(2+)-translocating ATPase (ZntA/CadA)", 3,
                  ec=("7.2.2.12", "7.2.2.21", "3.6.3.5", "3.6.3.3"), tcdb=("3.A.3.5",),
                  name_patterns=(r"\bzntA\b", r"\bcadA\b", r"\bpbrA\b",
                                 r"(?:zinc|cadmium|lead).{0,12}(?:translocating|exporting).{0,12}ATPase",
                                 r"Zn\(2\+\).{0,12}ATPase"),
                  source="TCDB 3.A.3.5 / EC 7.2.2.12 / MetaCyc Zn efflux"),
    SubstrateRule("Zn2+", "metal", "cation diffusion facilitator (CDF: CzcD/ZitB/YiiP/FieF)", 2,
                  tcdb=("2.A.4",),
                  name_patterns=(r"\bzitB\b", r"\bczcD\b", r"\byiiP\b", r"\bfieF\b",
                                 r"cation diffusion facilitator", r"\bCDF\b", r"zinc transporter ZnT"),
                  source="TCDB 2.A.4 (CDF) / KEGG ko02000"),
    SubstrateRule("Zn2+", "metal", "Zn(2+) ABC importer (ZnuABC) / ZIP uptake (ZupT)", 2,
                  tcdb=("3.A.1.15", "2.A.5"),
                  name_patterns=(r"\bznuA\b", r"\bznuB\b", r"\bznuC\b", r"\bzupT\b",
                                 r"zinc.{0,12}ABC", r"high-affinity zinc"),
                  source="TCDB 3.A.1.15 (ABC) / 2.A.5 (ZIP) / KEGG ko02000"),
    SubstrateRule("Cd2+", "metal", "cadmium-specific efflux / resistance determinant (CadD/CadX)", 3,
                  name_patterns=(r"\bcadD\b", r"\bcadC\b", r"cadmium resistance", r"cadmium efflux"),
                  source="MetaCyc Cd resistance / pI258 cad operon"),
    SubstrateRule("Pb2+", "metal", "lead resistance determinant (PbrABCD)", 3,
                  name_patterns=(r"\bpbrA\b", r"\bpbrB\b", r"\bpbrC\b", r"\bpbrD\b", r"\bpbrT\b",
                                 r"lead resistance", r"lead.{0,8}efflux"),
                  source="MetaCyc / Cupriavidus pbr operon"),

    # ---- cobalt / nickel (NiCoT, Ni-ABC, CnrCBA RND, urease accessory, RcnA) ----
    SubstrateRule("Ni2+", "metal", "nickel/cobalt transporter (NiCoT: HoxN/NikMNQO) + Ni ABC (NikABCDE)", 3,
                  tcdb=("2.A.52", "3.A.1.5"),
                  name_patterns=(r"\bnik[A-Z]\b", r"\bhoxN\b", r"\bnixA\b",
                                 r"nickel.{0,12}(?:transport|ABC|permease)", r"nickel/cobalt"),
                  source="TCDB 2.A.52 (NiCoT) / 3.A.1.5 / KEGG ko02000"),
    SubstrateRule("Ni2+", "metal", "urease (Ni metalloenzyme) + accessory (UreABC/UreEFGH)", 2,
                  ec=("3.5.1.5",),
                  name_patterns=(r"\bure[A-Z]\b", r"\burease\b"),
                  source="EC 3.5.1.5 / MetaCyc urea degradation (Ni cofactor)"),
    SubstrateRule("Co2+", "metal", "cobalt/nickel resistance RND + Co efflux (CnrCBA/RcnA/CorA)", 2,
                  tcdb=("2.A.6.1",),
                  name_patterns=(r"\bcnrA\b", r"\bcnrB\b", r"\bcnrC\b", r"\brcnA\b", r"\bcorA\b",
                                 r"cobalt.{0,12}(?:efflux|resistance|transport)"),
                  source="TCDB 2.A.6.1 / MetaCyc Co resistance"),

    # ---- mercury (mer operon: MerA reductase, MerB lyase, MerTPC transport) ----
    SubstrateRule("Hg2+", "metal", "mercuric reductase (MerA)", 3,
                  ec=("1.16.1.1",),
                  name_patterns=(r"\bmerA\b", r"mercuric reductase", r"mercury\(II\) reductase"),
                  source="EC 1.16.1.1 / MetaCyc mercury detoxification"),
    SubstrateRule("Hg2+", "metal", "organomercurial lyase (MerB) / Hg transport (MerTPCEF)", 3,
                  ec=("4.99.1.2",),
                  name_patterns=(r"\bmerB\b", r"\bmerT\b", r"\bmerP\b", r"\bmerC\b", r"\bmerE\b",
                                 r"\bmerF\b", r"alkylmercury lyase", r"mercur(?:y|ic).{0,12}transport"),
                  source="EC 4.99.1.2 / TCDB 9.A.2 (Mer) / MetaCyc"),

    # ---- arsenic (ars operon: ArsC reductase, ArsB/ACR3 efflux, ArsA ATPase) ----
    SubstrateRule("As3+", "metal", "arsenate reductase (ArsC)", 3,
                  ec=("1.20.4.1", "1.20.4.4"),
                  name_patterns=(r"\barsC\b", r"arsenate reductase"),
                  source="EC 1.20.4.1 / MetaCyc arsenate reduction"),
    SubstrateRule("As3+", "metal", "arsenite efflux permease (ArsB / ACR3) + ArsA ATPase", 3,
                  ec=("7.3.2.7",), tcdb=("2.A.45", "2.A.59"),
                  name_patterns=(r"\barsB\b", r"\barsA\b", r"\bacr3\b", r"\barsenical pump",
                                 r"arsenite.{0,12}(?:efflux|transport|permease)", r"\barsenite transporter"),
                  source="TCDB 2.A.45 (ArsB) / 2.A.59 (ACR3) / EC 7.3.2.7"),

    # ---- manganese / iron (NRAMP, Mn-ABC SitABCD, FeoB, ferric uptake) ----
    SubstrateRule("Mn2+", "metal", "Mn(2+)/Fe(2+) NRAMP transporter (MntH)", 2,
                  tcdb=("2.A.55",),
                  name_patterns=(r"\bmntH\b", r"\bnramp\b", r"natural resistance.{0,20}macrophage"),
                  source="TCDB 2.A.55 (NRAMP) / KEGG ko02000"),
    SubstrateRule("Mn2+", "metal", "Mn(2+) ABC importer (MntABC / SitABCD / PsaABC)", 2,
                  tcdb=("3.A.1.15",),
                  name_patterns=(r"\bmntA\b", r"\bmntB\b", r"\bmntC\b", r"\bsit[A-D]\b", r"\bpsa[A-C]\b",
                                 r"manganese.{0,12}ABC", r"manganese transport"),
                  source="TCDB 3.A.1.15 / KEGG ko02000"),
    SubstrateRule("Fe", "metal", "ferrous iron transport (FeoAB) / Fe ABC / siderophore (TonB-dependent)", 2,
                  tcdb=("2.A.108", "3.A.1.14"),
                  name_patterns=(r"\bfeoA\b", r"\bfeoB\b", r"\bfecA\b", r"\bfepA\b", r"\bfhuA\b",
                                 r"ferrous iron", r"ferric.{0,12}(?:transport|ABC|uptake)",
                                 r"siderophore", r"TonB-dependent.{0,12}receptor"),
                  source="TCDB 2.A.108 (FeoB) / 3.A.1.14 / MetaCyc Fe uptake"),

    # ---- chromate ----
    SubstrateRule("CrO4", "metal", "chromate ion transporter (ChrA / CHR family)", 2,
                  tcdb=("2.A.51",),
                  name_patterns=(r"\bchrA\b", r"chromate (?:transport|efflux|resistance)", r"\bCHR\b"),
                  source="TCDB 2.A.51 (CHR) / MetaCyc chromate resistance"),

    # ---- non-metal MerR-superfamily signals (multidrug, redox) -- keep the map from being metal-only ----
    SubstrateRule("multidrug", "multidrug", "multidrug efflux pump (RND AcrB / MFS Bmr-Blt / ABC)", 1,
                  tcdb=("2.A.6.2", "2.A.1.2", "2.A.1.3"),
                  name_patterns=(r"\bacrB\b", r"\bacrA\b", r"\bbmr\b", r"\bblt\b", r"\bemrB\b",
                                 r"\bmexB\b", r"multidrug.{0,12}(?:efflux|transport|resistance)",
                                 r"drug.{0,8}efflux"),
                  source="TCDB 2.A.6.2 (RND) / 2.A.1 (MFS) / KEGG ko01504"),
    SubstrateRule("superoxide/redox", "redox", "redox-stress determinant (SoxR regulon: SodA/Fpr/oxidative)", 2,
                  name_patterns=(r"\bsodA\b", r"\bsodB\b", r"\bfpr\b", r"superoxide dismutase",
                                 r"oxidative stress", r"ferredoxin.{0,8}reductase"),
                  source="MetaCyc oxidative-stress response / SoxRS regulon"),
]

# pre-split EC / TCDB indices for O(1)-ish lookup
_EC_INDEX: dict[str, list[SubstrateRule]] = {}
_TCDB_RULES = [r for r in RULES if r.tcdb]
for _r in RULES:
    for _ec in _r.ec:
        _EC_INDEX.setdefault(_ec, []).append(_r)
_COMPILED = [(r, [re.compile(p, re.I) for p in r.name_patterns]) for r in RULES]

# generic transporter/enzyme words that are not discriminative as ligand-gene keywords (dropped by
# keywords_for_effector so only gene names like zntA/cadA and element words like zinc survive).
_GENERIC_STEMS = {
    "translocating", "exporting", "transporting", "transport", "resistance", "efflux", "permease",
    "receptor", "dependent", "diffusion", "facilitator", "high", "affinity", "transporter", "oxidase",
    "reductase", "lyase", "macrophage", "natural", "oxidative", "stress", "dismutase", "ferredoxin",
    "atpase", "abc", "cation", "protein", "family", "membrane", "putative", "uncharacterized", "pump",
    # generic single words that, alone, match unrelated genes (e.g. 'chaperone' -> every DnaK/SurA) --
    # they are only meaningful as part of a phrase ('copper chaperone'), which is kept intact below
    "chaperone", "insertion", "assembly", "cofactor", "periplasmic", "factor", "inhibitor", "binding",
    "determinant", "regulator", "regulatory", "operon", "resist",
}


# --------------------------------------------------------------------------- classification
# a transcriptional regulator is never the substrate HANDLER -- guard against TF gene names whose stem
# collides with a transporter (zntR vs zntA, copR vs copA, arsR vs arsB) being voted as a transporter.
_REGULATOR_RE = re.compile(r"regulator|repressor|activator|transcription(?:al)? factor|"
                           r"DNA-binding|HTH[- ]|helix-turn-helix", re.I)


def is_regulator(product: str = "", name: str = "") -> bool:
    return bool(_REGULATOR_RE.search(f"{product or ''} {name or ''}"))


def _norm_ec(ec) -> list[str]:
    if not ec:
        return []
    if isinstance(ec, str):
        ec = [ec]
    return [e.strip() for e in ec if e and e.strip()]


def _ec_match(rule: SubstrateRule, ecs: list[str]) -> str | None:
    """Match an entry EC against a rule EC, honouring trailing '-' wildcards (e.g. 1.10.3.-)."""
    for want in rule.ec:
        for have in ecs:
            if want == have:
                return have
            if want.endswith(".-"):
                stem = want[:-1]                       # '1.10.3.'
                if have.startswith(stem):
                    return have
    return None


def classify_gene(product: str = "", name: str = "", *, ec=None, tcdb=None) -> SubstrateHit | None:
    """Best (most specific) substrate classification of one gene, or None.

    `product`/`name` are the GFF/UniProt free-text fields; `ec`/`tcdb` (str or list) are optional curated
    cross-references (from a UniProt entry) that raise confidence. Specificity ranks gene-name (3) > EC (2)
    > TCDB family (1); ties broken by rule order. Returns a `SubstrateHit` with full provenance.

    A gene whose *product* identifies it as a transcriptional regulator is never a substrate handler and
    returns None (so a TF like zntR/copR/arsR is not voted as its own cognate transporter)."""
    if is_regulator(product, ""):                  # guard on product only (a bare name stem is too weak)
        return None
    text = f"{product or ''} {name or ''}".strip()
    ecs = _norm_ec(ec)
    tcs = _norm_ec(tcdb)
    best: SubstrateHit | None = None

    def _consider(hit: SubstrateHit):
        nonlocal best
        if best is None or hit.specificity > best.specificity:
            best = hit

    # 1) gene-name / product regex (specificity 3)
    if text:
        for rule, pats in _COMPILED:
            for pat in pats:
                m = pat.search(text)
                if m:
                    _consider(SubstrateHit(rule.effector, rule.role, rule.substrate_class, 3,
                                           "name", m.group(0), rule.source, name, product))
                    break
    # 2) EC number (specificity 2)
    for e in ecs:
        for rule in _EC_INDEX.get(e, []):
            _consider(SubstrateHit(rule.effector, rule.role, rule.substrate_class, 2,
                                   "ec", e, rule.source, name, product))
        # wildcard EC rules (rule asks for 'x.y.z.-')
        for rule in RULES:
            hit_ec = _ec_match(rule, [e]) if any(w.endswith(".-") for w in rule.ec) else None
            if hit_ec:
                _consider(SubstrateHit(rule.effector, rule.role, rule.substrate_class, 2,
                                       "ec", hit_ec, rule.source, name, product))
    # 3) TCDB family prefix (specificity 1)
    for tc in tcs:
        for rule in _TCDB_RULES:
            if any(tc.startswith(fam) for fam in rule.tcdb):
                _consider(SubstrateHit(rule.effector, rule.role, rule.substrate_class, 1,
                                       "tcdb", tc, rule.source, name, product))
    return best


# bare element-name fallback (specificity 0): a generic product like "zinc transporter" with no curated
# gene name / EC / TCDB still names the ION. Kept SEPARATE from the curated RULES so a specific class
# always outranks it; consumed by regulon_enrichment as the tier below name/EC/TCDB and above ESM.
_ELEMENT_NAME = [
    (re.compile(r"\bzinc\b", re.I), "Zn2+"), (re.compile(r"\bcadmium\b", re.I), "Cd2+"),
    (re.compile(r"\bcopper\b|\bcupr(?:ic|ous)\b", re.I), "Cu+"), (re.compile(r"\bmercur", re.I), "Hg2+"),
    (re.compile(r"\barsen", re.I), "As3+"), (re.compile(r"\bnickel\b", re.I), "Ni2+"),
    (re.compile(r"\bcobalt\b", re.I), "Co2+"), (re.compile(r"\blead\b|\bplumb", re.I), "Pb2+"),
    (re.compile(r"\bmanganese\b", re.I), "Mn2+"), (re.compile(r"\bchromate\b|\bchromium\b", re.I), "CrO4"),
    (re.compile(r"\bferrous\b|\bferric\b|\biron\b|\bsiderophore\b", re.I), "Fe"),
    (re.compile(r"\bgold\b|\baurum\b", re.I), "Au+"),
]


def classify_generic(product: str = "", name: str = "") -> SubstrateHit | None:
    """Element-name fallback (specificity 0): name the ion from a bare metal word in the product when the
    curated `classify_gene` found nothing. Lower tier -- regulon_enrichment only uses it after classify_gene.
    Regulators are excluded (a 'zinc uptake regulator' is not a zinc handler)."""
    if is_regulator(product, ""):
        return None
    text = f"{product or ''} {name or ''}"
    for rx, ion in _ELEMENT_NAME:
        m = rx.search(text)
        if m:
            return SubstrateHit(ion, "metal", f"product names '{m.group(0)}'", 0, "element",
                                m.group(0), "product-name element token", name, product)
    return None


def keywords_for_effector(effector: str | None) -> list[str]:
    """All product-name keyword stems that flag a gene handling `effector` -- drop-in, substrate-resolved
    replacement for report.ligand_regulon._ELEMENT_KEYWORDS. Accepts an ion ('Zn2+') or an element token
    ('ZN'); returns lowercase literal stems extracted from the matching rules' name patterns."""
    if not effector:
        return []
    eff = effector.upper().rstrip("+-0123456789 ")          # 'ZN2+' -> 'ZN', 'Zn2+' -> 'ZN'
    out: set[str] = set()
    for rule in RULES:
        rule_el = rule.effector.upper().rstrip("+-0123456789 ")
        if rule_el != eff and rule.effector.upper() != effector.upper():
            continue
        for pat in rule.name_patterns:
            # strip regex word-boundary escapes (\b) so they don't leak a leading 'b' into the stems.
            clean = pat.replace(r"\b", " ").strip()
            if re.search(r"[.{}\[\]()+*?|]", clean):
                # regex-y pattern (e.g. 'copper.{0,8}efflux'): pull only discriminative literal stems.
                for tok in re.findall(r"[A-Za-z][A-Za-z0-9]{2,}", clean):
                    if tok.lower() not in _GENERIC_STEMS:
                        out.add(tok.lower())
            elif " " in clean:
                # MULTI-WORD descriptive phrase ('copper chaperone'): keep INTACT so it matches the phrase,
                # not a bare generic token ('chaperone' alone would flag every DnaK/SurA).
                out.add(clean.lower())
            elif clean.lower() not in _GENERIC_STEMS:
                # single discriminative gene-name stem ('copa', 'atx1', 'znta').
                out.add(clean.lower())
    return sorted(out)


# --------------------------------------------------------------------------- self-test
def _demo() -> None:
    cases = [
        # (product, name, ec, tcdb, expected effector, expected matched_by)
        ("Lead, cadmium, zinc and mercury-transporting ATPase", "zntA", None, None, "Zn2+", "name"),
        ("Copper-exporting P-type ATPase", "copA", "7.2.2.8", "3.A.3.5", "Cu+", "name"),
        ("Mercuric reductase", "merA", "1.16.1.1", None, "Hg2+", "name"),
        ("Arsenate reductase", "arsC", "1.20.4.1", None, "As3+", "name"),
        ("Arsenical pump membrane protein", "arsB", None, "2.A.45.1.1", "As3+", "name"),
        ("Nickel/cobalt transporter", "nikB", None, "2.A.52.1.1", "Ni2+", "name"),
        ("Cation diffusion facilitator family transporter", "czcD", None, "2.A.4.1.1", "Zn2+", "name"),
        ("Manganese transport protein", "mntH", None, "2.A.55.1.1", "Mn2+", "name"),
        ("hypothetical protein", "orf1", None, None, None, None),                # no signal
        ("MerR-family transcriptional regulator", "zntR", None, None, None, None),  # regulator guard (not a transporter)
        ("ArsR family transcriptional repressor", "arsR", None, None, None, None),  # regulator guard
        ("multidrug efflux RND transporter", "acrB", None, "2.A.6.2.1", "multidrug", "name"),
        # EC-only (no informative product name) still classifies via the EC index
        ("putative metal translocase", "yfgX", "7.2.2.12", None, "Zn2+", "ec"),
        # TCDB-only (broadest) -- product gives no hint, TC family carries the ion
        ("uncharacterized membrane protein", "ybxY", None, "3.A.3.5.7", "Cu+", "tcdb"),
    ]
    for product, name, ec, tcdb, exp_eff, exp_by in cases:
        hit = classify_gene(product, name, ec=ec, tcdb=tcdb)
        got_eff = hit.effector if hit else None
        got_by = hit.matched_by if hit else None
        ok = (got_eff == exp_eff)
        # for the TCDB-3.A.3.5 ambiguous case either Cu/Zn rule may win on a tie; accept any metal
        if name == "ybxY":
            ok = bool(hit and hit.matched_by == "tcdb" and hit.effector in METAL_IONS)
        assert ok, f"{name!r}: expected {exp_eff} via {exp_by}, got {got_eff} via {got_by} ({hit})"
        tag = (f"{hit.effector:9s} [{hit.matched_by}:{hit.detail}] {hit.substrate_class}"
               if hit else "(no substrate signal)")
        print(f"  {name:8s} {product[:46]:46s} -> {tag}")

    # keyword expansion is substrate-resolved (catches zntA/cadA/czcD/znuA for zinc, not a hand list)
    znkw = keywords_for_effector("Zn2+")
    assert "znta" in znkw and "cada" in znkw and "zitb" in znkw, f"zinc keywords thin: {znkw}"
    askw = keywords_for_effector("As3+")
    assert "arsb" in askw and "arsc" in askw and "acr3" in askw, f"arsenic keywords thin: {askw}"
    # element-name fallback (specificity 0) + regulator guard
    gen = classify_generic("cobalt-zinc-cadmium efflux system membrane protein", "czcA")
    assert gen and gen.specificity == 0 and gen.matched_by == "element", f"generic fallback failed: {gen}"
    assert classify_generic("zinc uptake transcriptional regulator", "zur") is None, "regulator not guarded"
    print(f"\n  generic fallback: czcA -> {gen.effector} (element '{gen.detail}'); Zur regulator guarded")

    print(f"\n  Zn2+ keywords: {znkw}")
    print(f"  As3+ keywords: {askw}")
    print("OK: transporter/efflux/enzyme classes mapped to ions at substrate resolution "
          "(name>EC>TCDB); keyword expansion substrate-resolved.")


if __name__ == "__main__":
    _demo()
