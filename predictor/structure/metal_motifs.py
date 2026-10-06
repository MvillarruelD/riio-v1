"""
metal_motifs.py -- the one catalogue of metal-coordination signatures the pipeline recognises.

Until 2026-09-02 the chemistry behind every metal call was spread over four places that nobody could
read side by side: the per-family branches of `structure/metal_site.coordination_gate`, the family lists
and precedence rules of `metal_site.inducer_class`, the per-clade `gate` strings in
`annotate/ssn_clusters`, and the regulon gene regexes in the analysis workspace. Nothing tied a rule to
its literature source, and nothing checked a stated prevalence against the data. That is not a
bookkeeping complaint: `metal_site` states the NikR H-x-H-x(4,7)-C site is present in "<=3.2 % of every
other family (Fur 3.2 %)" and it is in **60.6 % of Fur** -- harmless only because that gate happens to be
family-scoped.

So this module holds the signatures as DATA, each carrying its citation, the chemistry it implies, and
whether it may be used outside its family. `tools/build_metal_reference.py` measures every entry against
the vendored SSN member DBs -- per family AND per clade -- and writes both
`predictor/data/refs/metal_signatures.json` and `docs/METAL_SIGNATURES.md`. The prevalence numbers are
therefore GENERATED, never restated by hand, and the NikR/Fur error class cannot recur.

Two sources:

  * **Rondon, Antelo, Giedroc & Capdevila**, "Metallostasis" (Chem. Rev., in revision) section 3.2 --
    the per-family sensing chemistry: which residues, in which secondary-structure element, for which
    metal. Section numbers below are that manuscript's.
  * **Dudev & Lim**, Chem. Rev. 2014, 114, 538-556 -- the physics of selectivity: coordination number
    and geometry (their Table 1), ligand hardness (their Fig. 10 and section 6.3), and the
    Irving-Williams ordering. Used for `SELECTIVITY`, below.

WHAT THIS IS NOT. `SELECTIVITY` is a ligand-set heuristic for sanity-checking and tie-breaking. It is
not a thermodynamic calculation and cannot rank two ions that a site could plausibly bind: Dudev & Lim's
own conclusion is that the protein matrix and the cellular free-ion concentration decide the outcome,
and we model neither. A motif says "this looks like the site that senses X", never "X outcompetes Y here".

    python -m predictor.structure.metal_motifs        # self-test (offline)
"""
from __future__ import annotations

import re
from dataclasses import dataclass


@dataclass(frozen=True)
class MetalMotif:
    """One coordination signature, with everything needed to defend it in a Methods section."""

    name: str
    #: Canonical `ssn_clusters` family name this rule belongs to ("ArsR/SmtB"), or "" if family-agnostic.
    family: str
    #: What the site senses. Empty when the motif marks a class rather than an ion (e.g. an Fe-S cluster).
    implies: tuple = ()
    #: 'metal' | 'redox' | 'non-metal' | '' -- the class the motif argues for.
    implied_class: str = ""
    #: Python regex over the uppercase protein sequence. Empty when the rule is not expressible as one
    #: (see `implemented`), in which case the entry is DOCUMENTATION only.
    pattern: str = ""
    #: Coordination geometry and number, from the review / Dudev & Lim.
    geometry: str = ""
    coordination_number: str = ""
    ligand_set: str = ""
    #: 'family' -- valid ONLY when the query's family is `family`; measured to fire elsewhere.
    #: 'global' -- specific enough to be trusted without a family gate.
    scope: str = "family"
    #: May this motif, on its own, NAME an ion in the inducer call? Set from the measured per-clade
    #: concentration, not from how confident the literature sounds. A motif that marks a chemistry
    #: shared by several ions (a bare Cys pair) documents the site without naming the metal.
    names_ion: bool = False
    #: True when `pattern` is set and the code actually tests it. False entries are carried so the
    #: reference table is a complete account of the family's chemistry rather than only its automatable
    #: part -- a reader must be able to see what we do NOT detect.
    implemented: bool = True
    source: str = ""
    note: str = ""

    def search(self, seq: str):
        """First match in `seq`, or None. Returns None for documentation-only entries."""
        if not self.pattern:
            return None
        return re.search(self.pattern, (seq or "").upper())


# --------------------------------------------------------------------------- the catalogue
_R = "Rondon et al., Metallostasis (Chem. Rev.), "

MOTIFS: tuple = (
    # ---------------------------------------------------------------- ArsR/SmtB
    MetalMotif(
        name="ArsR_a5_DxHx10Hx2HE", family="ArsR/SmtB",
        implies=("Zn2+", "Ni2+"), implied_class="metal",
        pattern=r"D.H.{10}H.{2}[HE]",
        geometry="tetrahedral (Zn) / octahedral (Ni, completing with the N-terminal His-Gly)",
        coordination_number="4-6", ligand_set="Asp + 2-3 His (+/- Glu)",
        scope="global", names_ion=True,
        source=_R + "3.2.2.1",
        note="The interfacial alpha5 site shared by the CzrA/SmtB/ZiaR zinc sensors AND by NmtR (Ni), "
             "which completes an octahedron using its N-terminal His-Gly. It cannot separate Zn from Ni "
             "on sequence -- the review is explicit that they share the motif -- so it implies BOTH."),
    MetalMotif(
        name="ArsR_a3_CVC", family="ArsR/SmtB",
        implies=(), implied_class="metal",
        pattern=r"C[VA]C",
        geometry="trigonal planar (As) / varies", coordination_number="3-4",
        ligand_set="the alpha3 C-(V/A)-C thiol pair (+ a third Cys nearby or C-terminal)",
        scope="family", names_ion=False,
        source=_R + "3.2.2.3",
        note="The review names C-(V/A)-C as the canonical arsenite-sensor motif, and it does concentrate "
             "in the As(III) clade (measured: ArsR_c1 90 %). It DOES NOT identify arsenic on its own: it "
             "also fires on 71 % of ArsR_c3 and, decisively, on 41 % of ArsR_c2 -- the CzrA/SmtB/ZiaR "
             "ZINC clade. It is a sharper alpha3-thiol-site detector than a bare C-x-C (2.8 % vs 31.9 % "
             "out of family), and that is all it may be used for."),
    MetalMotif(
        name="ArsR_a3_thiol_CxC", family="ArsR/SmtB",
        implies=(), implied_class="metal",
        pattern=r"C.C",
        geometry="varies", coordination_number="3-4", ligand_set="Cys pair from alpha3 (+/- 2 N-terminal Cys)",
        scope="family", names_ion=False,
        source=_R + "3.2.2.1",
        note="The alpha3 thiol site. With two further N-terminal Cys this is the alpha3N site of the "
             "Cd/Pb sensors CadC and LmCadC. Deliberately does NOT name an ion: the review notes this "
             "Cys-rich site keeps only partial selectivity for Cd/Pb over Zn, and four Cys in the site "
             "actively degrades Pb specificity."),
    MetalMotif(
        name="ArsR_AfArsR_CCx5C", family="ArsR/SmtB",
        implies=("As3+",), implied_class="metal",
        pattern=r"CC.{5}C",
        geometry="trigonal planar", coordination_number="3", ligand_set="3 Cys",
        scope="family", names_ion=False,
        source=_R + "3.2.2.3",
        note="The second, distinct As(III) site of Acidithiobacillus ferrooxidans ArsR."),
    MetalMotif(
        name="ArsR_CmtR_a4_plus_tail", family="ArsR/SmtB",
        implies=("Cd2+",), implied_class="metal", pattern="", implemented=False,
        ligand_set="2 Cys from the alpha4 DNA-binding helix + 1 Cys from the C-terminal tail",
        scope="family", source=_R + "3.2.2.2",
        note="NOT DETECTED. The ligands come from three different structural elements, so no linear "
             "spacing identifies them; separating CmtR from any other multi-Cys ArsR needs an alignment."),

    # ---------------------------------------------------------------- MerR
    MetalMotif(
        name="MerR_ZntR_CCx3Hx4C", family="MerR",
        implies=("Zn2+",), implied_class="metal",
        pattern=r"CC.{3}H.{4}C",
        geometry="binuclear, both sites ~tetrahedral", coordination_number="4",
        ligand_set="Cys/Cys/His/Cys from the metal-binding loop + a Cys from the alpha5 N-terminus",
        scope="global", names_ion=True,
        source=_R + "3.2.3.2",
        note="The canonical MerR-family zinc efflux regulator (EcZntR)."),
    MetalMotif(
        name="MerR_SoxR_FeS", family="MerR",
        implies=(), implied_class="redox",
        pattern=r"C.{2}C.C.{4,6}C",
        geometry="[2Fe-2S] cluster", coordination_number="4", ligand_set="4 Cys (CIGCGCLSxxC)",
        scope="global", names_ion=False,
        source=_R + "3.2.3.6",
        note="The redox-sensing MerR branch: a sensory [2Fe-2S] cluster responding to superoxide, NO "
             "and redox-cycling agents. A cluster, not a sensed metal ion -- hence no `implies`."),
    MetalMotif(
        name="MerR_metal_CX7_8C", family="MerR",
        implies=(), implied_class="metal",
        pattern=r"C.{7,8}C",
        geometry="varies with the branch", coordination_number="2-4",
        ligand_set="the reactive Cys pair of the metal-binding loop",
        scope="family", names_ion=False,
        source=_R + "3.2.3",
        note="The existing MerR gate. It separates the metal-sensing clades from the non-metal ones "
             "(TnrA/GlnR/NmlR/YfmP) but cannot say WHICH metal -- the same pair serves Hg, Cd/Pb, Zn "
             "and Cu. `metal_site.coordination_gate` additionally requires it C-terminally."),
    MetalMotif(
        name="MerR_CueR_coinage", family="MerR",
        implies=("Cu+", "Ag+", "Au+"), implied_class="metal", pattern="", implemented=False,
        geometry="linear bis-thiolate", coordination_number="2", ligand_set="2 Cys from a SHORT MBL",
        scope="family", source=_R + "3.2.3.4",
        note="NOT DETECTED. Selectivity for the monovalent coinage metals comes from two NEGATIVE design "
             "elements -- a conserved Ser in alpha5 that excludes divalents, and an MBL shorter than the "
             "divalent sensors' -- neither of which is a presence/absence motif. Distinguishing CueR from "
             "a divalent sensor needs the alignment, not a regex."),
    MetalMotif(
        name="MerR_PbrR_CadR_trigonal", family="MerR",
        implies=("Pb2+", "Cd2+"), implied_class="metal", pattern="", implemented=False,
        geometry="trigonal pyramidal (Pb, with a stereochemically active lone pair)",
        coordination_number="3 (+1 Asn in CadR)",
        ligand_set="3 Cys: 2 from the MBL + 1 from the dimerisation helix",
        scope="family", source=_R + "3.2.3.1",
        note="NOT DETECTED. The third Cys comes from a different helix, so the ligand set has no fixed "
             "spacing. CadR adds an Asn as a fourth ligand that Pb's lone pair excludes in PbrR."),

    # ---------------------------------------------------------------- Fur
    MetalMotif(
        name="Fur_HHxHx2Cx2C", family="Fur",
        implies=(), implied_class="metal",
        pattern=r"HH.H.{2}C.{2}C",
        geometry="site 1 tetrathiolate (structural); sites 2/3 His/carboxylate (regulatory)",
        coordination_number="4 (site 1); 4-6 (site 2)",
        ligand_set="His from HHxH -> regulatory sites 2/3; Cx2C -> half of the structural S4 site 1",
        scope="global", names_ion=False,
        source=_R + "3.2.6",
        note="The conserved hinge motif of the family. It marks a Fur-fold metal architecture but not "
             "which metal: the SAME site 2 is tuned to Fe(II), Zn(II), Mn(II) or Ni(II) across "
             "Fur/Zur/Mur/Nur, and PerR uses it to hold the Fe that reacts with peroxide. The review "
             "presents this family precisely as the illustration that one site is retuned, not replaced."),
    MetalMotif(
        name="Fur_site2_His", family="Fur",
        implies=(), implied_class="metal",
        pattern=r"HH[HD]H",
        geometry="His/carboxylate", coordination_number="4-6", ligand_set="His-rich regulatory site 2",
        scope="family", names_ion=False,
        source=_R + "3.2.6",
        note="The regulatory site alone. Family-scoped and NOT safe outside it -- measured, it also "
             "fires on a substantial share of NikR."),
    MetalMotif(
        name="Fur_site1_Cx2C", family="Fur",
        implies=("Zn2+",), implied_class="metal",
        pattern=r"C.{2}C",
        geometry="tetrathiolate", coordination_number="4", ligand_set="2x Cx2C",
        scope="family", names_ion=False,
        source=_R + "3.2.6.1",
        note="The STRUCTURAL site: nearly always Zn, stabilising the dimer, and a Zn here behaves as "
             "apoprotein. Finding it says nothing about what the protein senses -- which is why it does "
             "not name an ion despite implying Zn occupancy."),

    # ---------------------------------------------------------------- NikR
    MetalMotif(
        name="NikR_His3Cys", family="NikR",
        implies=("Ni2+",), implied_class="metal",
        pattern=r"H.H.{4,7}C",
        geometry="square planar", coordination_number="4",
        ligand_set="His-x-His-x(4,7)-Cys (+ a His from the adjacent subunit)",
        scope="family", names_ion=True,
        source=_R + "3.2.8; Dudev & Lim 2014 section 6.4.2",
        note="FAMILY-SCOPED, and this is the entry that proves why the distinction matters: the code "
             "comment in metal_site.py claims this spacing occurs in <=3.2 % of other families, but "
             "measured on the current member DB it fires on a large fraction of Fur. Inside NikR it is "
             "the real square-planar Ni site (Dudev & Lim show the cavity disfavours Zn/Co/Cu, which "
             "prefer other geometries); outside it, it is noise."),

    # ---------------------------------------------------------------- CsoR / RcnR / FrmR
    MetalMotif(
        name="CsoR_WXYZ_metal_site", family="CsoR/FrmR",
        implies=(), implied_class="metal",
        pattern=r"C.{15,40}H.{3}C",
        geometry="trigonal S2N (Cu I) / varies", coordination_number="3-4",
        ligand_set="the W-X-Y-Z frame: X=Cys, Y=His', Z=Cys' (MtCsoR Cys36/His61'/Cys65')",
        scope="family", names_ion=False,
        source=_R + "3.2.4.1",
        note="The subunit-bridging metal site of this superfamily. Measured, it is an excellent METAL / "
             "NON-METAL discriminator within the family -- 95-99 % across CsoR_mtCsoR, CsoR_bsCsoR, "
             "CsoR_syInrS and CsoR_spNreA, against 0.0 % in both FrmR clades (formaldehyde) and "
             "0.3-1.4 % in the CstR persulfide clades. It does NOT name the metal: syInrS is a Ni/Co "
             "sensor and scores as highly as the Cu(I) clades, exactly as the shared W-X-Y-Z frame "
             "predicts. It also misses the Thermus x-C-H-H variant (10 % of CsoR_ttCsoR), which "
             "substitutes His for the second Cys and is still Cu-specific."),
    MetalMotif(
        name="CsoR_RcnR_Nterm_His", family="CsoR/FrmR",
        implies=("Ni2+", "Co2+"), implied_class="metal", pattern="", implemented=False,
        ligand_set="an N-terminal His donor added to the W-X-Y-Z frame (RcnR His3/Cys35/His60/His64)",
        scope="family", source=_R + "3.2.4.2",
        note="NOT DETECTED. RcnR is nested INSIDE the FrmR clade in our SSN, and the N-terminal His that "
             "distinguishes it is a single residue at the chain start. This is the known root cause of "
             "the RcnR->FrmR misassignment; a regex does not fix it."),
    MetalMotif(
        name="CsoR_FrmR_ligand_loss", family="CsoR/FrmR",
        implies=(), implied_class="non-metal", pattern="", implemented=False,
        ligand_set="Pro2 + Cys35 ONLY -- His64 replaced by Glu, N-terminal His absent",
        scope="family", source=_R + "3.2.4.2",
        note="NOT DETECTED as a motif, but its ABSENCE is what the live CsoR gate tests: "
             "`coordination_gate` requires >=2 His and >=1 Cys, which E. coli FrmR (1 His, 2 Cys) fails. "
             "Formaldehyde is sensed by a Pro2/Cys35 methylene bridge, not by a metal."),

    # ---------------------------------------------------------------- CopY
    MetalMotif(
        name="CopY_CxC", family="CopY",
        implies=("Cu+",), implied_class="metal",
        pattern=r"C.C",
        geometry="tetrahedral (Zn, resting) -> S4-Cu2 cluster (Cu-loaded)",
        coordination_number="4", ligand_set="CxC; a single Zn per dimer, displaced by two Cu(I)",
        scope="family", names_ion=True,
        source=_R + "3.2.5",
        note="The resting state holds a structural Zn that ACTIVATES DNA binding; rising Cu(I) displaces "
             "it to form the S4-Cu2 cluster that releases the operator. The BlaI/MecI branch of this "
             "family senses beta-lactam and is not a metal sensor at all."),

    # ---------------------------------------------------------------- Rrf2
    MetalMotif(
        name="Rrf2_FeS_Cys_triad", family="Rrf2",
        implies=(), implied_class="redox",
        pattern=r"C.{2,20}C.{2,20}C",
        geometry="[Fe-S] cluster", coordination_number="4",
        ligand_set="highly variable: IscR/RirA 3 Cys, NsrR 3 Cys + Asp, RsrR 1 His + 1 Glu + 2 Cys",
        scope="family", names_ion=False,
        source=_R + "3.2.9",
        note="Rrf2 senses cluster OCCUPANCY, NO, or thiol status -- never a free metal ion, which is why "
             "`coordination_gate` returns has_site=False for this family by design. The review is "
             "explicit that the ligand set is 'remarkably variable', so this pattern marks the Cys "
             "spacing and nothing more. SifR uses a SINGLE Cys to react with quinones and is missed."),

    # ---------------------------------------------------------------- MarR
    MetalMotif(
        name="MarR_AdcR_ZitR_linker", family="MarR/SlyA",
        implies=("Zn2+",), implied_class="metal", pattern="", implemented=False,
        geometry="two proximate pseudotetrahedral sites", coordination_number="4",
        ligand_set="ligands drawn from the unstructured region between alpha1 and alpha2",
        scope="family", source=_R + "3.2.11.1",
        note="NOT DETECTED. The discriminator is STRUCTURAL, not compositional: in the zinc-sensing "
             "MarRs (AdcR, ZitR) the ligands come from a disordered alpha1-alpha2 linker that is a "
             "continuous longer alpha2 helix in every other MarR. Detecting it needs the fold or an "
             "alignment. These are also the only characterised metal sensors in this large superfamily, "
             "and uniquely for MarR the metal ACTIVATES operator binding."),

    # ---------------------------------------------------------------- DtxR, TetR, GntR, LysR
    MetalMotif(
        name="DtxR_binuclear", family="DtxR/MntR",
        implies=("Fe2+", "Mn2+"), implied_class="metal", pattern="", implemented=False,
        geometry="binuclear; MntR's two sites are closer together and share ligands",
        coordination_number="6", scope="family", source=_R + "3.2.7",
        note="NOT DETECTED, and this is why the DtxR gate is silent: the site is assembled from "
             "dispersed His/Asp/Glu with no linear signature. Non-cognate metals (Fe, Zn, Co in MntR) "
             "bind MONONUCLEARLY and simply fail to trigger the allosteric response."),
    MetalMotif(
        name="TetR_SczA_dual_site", family="TetR/AcrR",
        implies=("Zn2+",), implied_class="metal", pattern="", implemented=False,
        geometry="two site pairs per dimer -- one tetrahedral, one octahedral",
        coordination_number="4 and 6", scope="family", source=_R + "3.2.12",
        note="NOT DETECTED. SczA is the lone metal-responsive TetR and activates a Zn efflux "
             "transporter; the rest of this very large family senses organic ligands in independently "
             "evolved C-terminal cavities."),
    MetalMotif(
        name="GntR_buried_structural_Zn", family="GntR",
        implies=("Zn2+",), implied_class="", pattern="", implemented=False,
        geometry="buried at the base of a solvent-accessible effector cavity",
        scope="family", source=_R + "3.2.10",
        note="NOT DETECTED, AND DELIBERATELY NOT A SENSOR CALL. ~70 % of GntR proteins carry a conserved "
             "Zn site, but the review states as an OPEN QUESTION whether it is regulatory or simply "
             "positions the carboxylate of the cognate organic acid. `metal_site.STRUCTURAL_SITE_FAMILIES` "
             "already encodes this caveat; finding a site here is a real observation, and 'therefore a "
             "metal sensor' is the step that does not follow."),
    MetalMotif(
        name="LysR_OxyR_peroxidatic_Cys", family="LysR-type (LTTR)",
        implies=(), implied_class="redox", pattern="", implemented=False,
        ligand_set="a peroxidatic Cys that forms a disulfide on H2O2 oxidation",
        scope="family", source=_R + "3.2.13",
        note="NOT DETECTED. Not a metal site at all -- OxyR senses H2O2 through Cys oxidation."),
)


# --------------------------------------------------------------------------- selectivity priors
#: Coordination number and geometry per ion (Dudev & Lim 2014, Table 1 -- CSD/PDB survey), with the
#: ligand preference their interaction energies establish (their Fig. 10 and section 6.3): among the
#: residues that coordinate transition metals, anionic Cys- is the MOST discriminative and strongly
#: prefers the soft Cu+/Hg2+/Cd2+; Asp-/Glu- carboxylates also discriminate; neutral His and the
#: backbone bind different metals with much more similar energies and so discriminate least.
#:
#: USE THIS TO SANITY-CHECK, NOT TO RANK. A geometry that contradicts an ion is evidence the site is
#: not that ion's; a geometry consistent with three ions does not order them. Dudev & Lim's own
#: conclusion is that the protein matrix and the cellular free-ion concentration decide the winner --
#: NikR binds Cu2+ more tightly than Ni2+ in vitro and is still a nickel sensor in vivo, because free
#: cytosolic Cu is ~10^-18 M. We model neither term.
SELECTIVITY: dict = {
    "Cu+":  {"cn": "2-4", "geometry": "linear / trigonal", "prefers": "Cys(thiolate)"},
    "Ag+":  {"cn": "2-3", "geometry": "linear / trigonal", "prefers": "Cys(thiolate)"},
    "Au+":  {"cn": "2",   "geometry": "linear",            "prefers": "Cys(thiolate)"},
    "Hg2+": {"cn": "2-4", "geometry": "linear / trigonal", "prefers": "Cys(thiolate)"},
    "Zn2+": {"cn": "4",   "geometry": "tetrahedral",       "prefers": "Cys / His / Asp-Glu"},
    "Cd2+": {"cn": "4-6", "geometry": "tetrahedral / octahedral", "prefers": "Cys(thiolate)"},
    "Pb2+": {"cn": "3-4", "geometry": "trigonal pyramidal (lone pair)", "prefers": "Cys(thiolate)"},
    "Ni2+": {"cn": "4 or 6", "geometry": "square planar or octahedral", "prefers": "His / Cys"},
    "Co2+": {"cn": "4-6", "geometry": "tetrahedral / octahedral", "prefers": "His / Cys"},
    "Fe2+": {"cn": "6",   "geometry": "octahedral",        "prefers": "His / Asp-Glu"},
    "Mn2+": {"cn": "6",   "geometry": "octahedral",        "prefers": "Asp-Glu / His"},
    "As3+": {"cn": "3",   "geometry": "trigonal pyramidal", "prefers": "Cys(thiolate)"},
}

#: Irving-Williams: the intrinsic affinity order of divalent first-row transition metals for ANY ligand
#: set (Dudev & Lim 2014 section 3). Every metalloregulator has to work against this ordering rather
#: than with it, which is why the cell holds free Zn at 10^-12-10^-15 M and free Cu at ~10^-18 M.
IRVING_WILLIAMS: tuple = ("Mg2+", "Mn2+", "Fe2+", "Co2+", "Ni2+", "Cu2+", "Zn2+")


# --------------------------------------------------------------------------- lookup
def motifs_for(family: str | None, *, implemented_only: bool = True) -> tuple:
    """Every catalogue entry usable for `family`: its own family rules plus the `global`-scope ones."""
    fam = family or ""
    out = []
    for m in MOTIFS:
        if implemented_only and not m.implemented:
            continue
        if m.family and m.family != fam and m.scope != "global":
            continue
        if m.scope == "global" and m.family and m.family != fam:
            # A global-scope motif is specific enough to trust anywhere, so it stays in the list even
            # when the query belongs to another family -- that is what `global` means.
            pass
        out.append(m)
    return tuple(out)


def match_motifs(seq: str, family: str | None) -> tuple:
    """The catalogue entries that FIRE on `seq`, given its family. Pure; no I/O."""
    seq = (seq or "").upper()
    return tuple(m for m in motifs_for(family) if m.search(seq))


def implied_ions(seq: str, family: str | None) -> dict:
    """{ion: [motif names]} from motifs that are allowed to name an ion (`names_ion`).

    An ion appears once per motif that argues for it. The caller decides what to do with a set of size
    >1 -- `effector.inducer` reports it as an unresolved candidate set rather than picking."""
    out: dict = {}
    for m in match_motifs(seq, family):
        if not m.names_ion:
            continue
        for ion in m.implies:
            out.setdefault(ion, []).append(m.name)
    return out


def by_name(name: str) -> MetalMotif | None:
    return next((m for m in MOTIFS if m.name == name), None)


# --------------------------------------------------------------------------- self-test
def _demo() -> None:
    names = [m.name for m in MOTIFS]
    assert len(names) == len(set(names)), "motif names must be unique"

    from predictor.annotate import ssn_clusters as _ssn
    fams = {m.family for m in MOTIFS if m.family}
    unknown = fams - set(_ssn.FAMILIES)
    assert not unknown, f"motif families must be canonical ssn_clusters names; got {unknown}"

    # every implemented entry must carry a compilable pattern; every documentation-only entry must not
    for m in MOTIFS:
        if m.implemented:
            assert m.pattern, f"{m.name}: implemented entries need a pattern"
            re.compile(m.pattern)
        else:
            assert not m.pattern, f"{m.name}: documentation-only entries must not carry a pattern"
        assert m.source, f"{m.name}: every entry needs a citation"
        assert m.note, f"{m.name}: every entry needs a note saying what it can and cannot do"
        if m.names_ion:
            assert m.implemented and m.implies, f"{m.name}: cannot name an ion without a pattern + ions"

    # a motif that names an ion must name one we have selectivity data for
    for m in MOTIFS:
        for ion in m.implies:
            assert ion in SELECTIVITY, f"{m.name}: {ion} missing from SELECTIVITY"

    # E. coli ZntR, the canonical MerR zinc sensor: its MBL motif must fire and must name Zn. Copied
    # verbatim from the shipped member DB record `MerR_9__P0ACS5`, so the test asserts against the same
    # sequence the pipeline would classify rather than against one written from memory.
    zntr = ("MYRIGELAKMAEVTPDTIRYYEKQQMMEHEVRTEGGFRLYTESDLQRLKFIRHARQLGFSLESIRELLSIRIDPEHHT"
            "CQESKGIVQERLQEVEARIAELQSMQRSLQRLNDACCGTAHSSVYCSILEALEQGASGVKSGC")
    hits = {m.name for m in match_motifs(zntr, "MerR")}
    assert "MerR_ZntR_CCx3Hx4C" in hits, hits
    ions = implied_ions(zntr, "MerR")
    assert "Zn2+" in ions, ions
    print(f"ZntR -> motifs {sorted(hits)}; ions {sorted(ions)}")

    # family scoping: a family-scoped motif must not be offered to another family
    nikr = by_name("NikR_His3Cys")
    assert nikr.scope == "family", "the NikR site is measured to fire on Fur; it must stay family-scoped"
    assert nikr not in motifs_for("Fur"), "family-scoped motifs must not leak into another family"
    assert nikr in motifs_for("NikR")

    # global-scope motifs are available regardless of family
    zn = by_name("MerR_ZntR_CCx3Hx4C")
    assert zn.scope == "global" and zn in motifs_for("Fur")

    n_impl = sum(1 for m in MOTIFS if m.implemented)
    n_name = sum(1 for m in MOTIFS if m.names_ion)
    print(f"OK: {len(MOTIFS)} signatures ({n_impl} detected, {len(MOTIFS)-n_impl} documented only); "
          f"{n_name} may name an ion; families all canonical.")


if __name__ == "__main__":
    _demo()
