"""
metal_site.py -- metal-coordination gate for inducer transfer (Phase 3 gate, Phase 2-adjacent).

SSN clusters map to inducers, but cluster -> *metal* is only safe when the coordination signature is
actually present (plan SS6: MerR clusters 3/5/10 are non-metal; "canonical" names mislead). This module
detects that signature so a metal label can be GATED on it -- the gate the unified inducer inference
(effector/inducer.py) applies before trusting a metal call.

Signatures:
  * MerR metal sensors -- a C-terminal reactive cysteine pair `C-X(7,8)-C` (CadR/ZntR/MerR/CueR clade);
    ABSENT in the non-metal TnrA/GlnR/NmlR subfamilies. This is the plan's explicit gate.
  * Fur -- a histidine/cysteine site-2 (`H..H..` + `C..C`), promiscuous -> advisory only.
  * ArsR/SmtB -- heterogeneous sites (alpha3 `ELCVCD..C`-like or alpha5 His); we report cysteine pairs
    but mark the call advisory, since ArsR metal-site location varies by subfamily.

Sequence detection is live here. Structural confirmation is `sensing_site_from_structure` below, which
reads a fold with its bound metals; the Folddisco constellation query this module once described was
measured over 139 candidates, answered none of them, and has been removed.

Run `python -m predictor.structure.metal_site --self-test`.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

_REPO = Path(__file__).resolve().parents[2]

# C-X(7,8)-C reactive cysteine pair (the MerR metal gate); also the generic divalent-Cys pair detector.
_CXC = re.compile(r"C.{7,8}C")
_CYS_PAIR = re.compile(r"C.{2,8}C")            # broader ArsR alpha3/alpha5-style Cys pairs
# NikR square-planar Ni(II) site: His87-x-His89-x(5)-Cys95 (E. coli). See coordination_gate for prevalence.
_NIKR_SITE = re.compile(r"H.H.{4,7}C")


def cys_pairs(seq: str, pattern=_CXC):
    """Start positions (0-based) of cysteine-pair matches; overlapping matches are found."""
    return [m.start() for m in re.finditer(rf"(?=({pattern.pattern}))", seq)]


def _catalogue_evidence(seq: str, family: str) -> dict:
    """The `metal_motifs` catalogue's read on this sequence: which signatures fire, and what they imply.

    Additive to the branches below, never a replacement for them. The branches decide `has_site` -- the
    metal/non-metal gate the inducer stage has always applied -- while this adds WHICH documented site
    was seen and, where a signature is concentrated enough in one clade to be diagnostic, which ion it
    argues for. Kept in one helper so every family branch reports it identically, including the families
    that have no branch at all."""
    try:
        from predictor.structure import metal_motifs as _mm
    except Exception:                                   # catalogue absent -> gate behaves exactly as before
        return {}
    hits = _mm.match_motifs(seq, family)
    if not hits:
        return {"motifs": [], "implied_ions": {}}
    return {"motifs": [m.name for m in hits],
            "implied_ions": _mm.implied_ions(seq, family)}


def coordination_gate(seq: str, family: str, *, fold_cif: str | None = None,
                      op_cif: str | None = None) -> dict:
    """Detect the metal-coordination signature. Returns dict(has_site, gate, positions, c_terminal,
    n_cys, advisory), plus `motifs` / `implied_ions` from the `metal_motifs` catalogue.

    `has_site` is what the inducer stage gates a *metal* label on, and it is decided exactly as before
    by the per-family branches in `_coordination_gate_impl`. `implied_ions` is the addition:
    {ion: [motif names]} from signatures measured to concentrate in a single SSN clade
    (`tools/build_metal_reference.py` computes that concentration; a motif that merely marks a chemistry
    shared by several ions contributes to `motifs` and names nothing). It is the only protein-level
    evidence about WHICH metal that does not come from the SSN clade -- which matters because on the
    last full run two thirds of candidates got no clade at all and therefore no ion.

    STRUCTURAL CONFIRMATION (production hook): when a `fold_cif` with bound metals is supplied, the
    sequence gate is upgraded with `structure_sites` -- every metal site typed structural-vs-regulatory
    (see `sensing_site_from_structure`). A two-site regulator (Fur) has a constitutive structural
    Zn(Cys4) cage and a separate regulatory site; only the REGULATORY-site metal is a candidate effector,
    so it is the one merged into `implied_ions`. A single site is passed through unchanged (the current
    single-site workflow is correct). Naming still respects the ceiling: `_finalize` promotes an
    implied ion only where the SSN clade is silent.

    Merged in a wrapper rather than inside each branch: there are seven return paths and one of them is
    the no-branch default, so adding the fields per branch would silently skip whichever was forgotten.
    """
    ev = _coordination_gate_impl(seq, family)
    ev.update(_catalogue_evidence(seq, family))
    if fold_cif:
        ev.update(_structure_confirmation(fold_cif, family, op_cif))
    return ev


# --------------------------------------------------------------------------- structural confirmation
# A structural site (constitutive Zn(Cys4) cage) is NOT the sensed effector; letting its metal name the
# inducer mislabels Fe/Mn/Ni sensors as Zn. Given a fold with bound metals we type every site and
# contribute ONLY the regulatory-site metal. AF3 folds often use a redox-inert surrogate (Co for the
# Fe/Mn/Ni response site); a surrogate names the metal CLASS and a Fe/Mn/Ni shortlist, never a single ion.
_METAL_ELEMENTS = {"ZN", "FE", "CO", "NI", "MN", "CU", "CD"}
_CANONICAL_ION = {"ZN": "Zn2+", "FE": "Fe2+", "MN": "Mn2+", "NI": "Ni2+", "CO": "Co2+",
                  "CU": "Cu+", "CD": "Cd2+"}
_SURROGATE_SHORTLIST = {"CO": ("Fe2+", "Mn2+", "Ni2+")}   # Co stands in for redox-labile divalents
_COORD_CUT = 2.8
_DNA_COMPS = {"DA", "DT", "DG", "DC"}


def _cif_atoms(text: str):
    """Minimal AF3/ModelCIF `_atom_site` reader -> list of (group, element, comp, chain, seq, x, y, z).

    Deliberately tolerant: any parse problem yields [] so the gate silently keeps its sequence-only
    behaviour rather than failing a run on a malformed fold.
    """
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        if lines[i].strip() == "loop_":
            j = i + 1
            cols = []
            while j < len(lines) and lines[j].lstrip().startswith("_"):
                cols.append(lines[j].strip())
                j += 1
            if cols and all(c.startswith("_atom_site.") for c in cols):
                idx = {c.split(".", 1)[1]: k for k, c in enumerate(cols)}
                need = ("group_PDB", "type_symbol", "label_comp_id", "label_asym_id",
                        "label_seq_id", "Cartn_x", "Cartn_y", "Cartn_z")
                if any(n not in idx for n in need):
                    return []
                out = []
                k = j
                while k < len(lines) and lines[k].strip() and not lines[k].lstrip().startswith(("_", "#", "loop_")):
                    p = lines[k].split()
                    if len(p) >= len(cols):
                        try:
                            s = p[idx["label_seq_id"]]
                            out.append((p[idx["group_PDB"]].upper(), p[idx["type_symbol"]].upper(),
                                        p[idx["label_comp_id"]].upper(), p[idx["label_asym_id"]],
                                        int(s) if s.lstrip("-").isdigit() else None,
                                        float(p[idx["Cartn_x"]]), float(p[idx["Cartn_y"]]), float(p[idx["Cartn_z"]])))
                        except (ValueError, IndexError):
                            pass
                    k += 1
                return out
            i = j
        else:
            i += 1
    return []


def sensing_site_from_structure(fold_cif: str, family: str, op_cif: str | None = None) -> dict:
    """Type every bound-metal site in a fold and return the sensed (regulatory) metal.

    Contract (the productionised prototype):
      * FIND every metal site (protein donor atoms within 2.8 A);
      * if MORE THAN ONE site, IDENTIFY the regulatory site -- the one coupled to the DNA-binding
        domain (ligands reaching the operator-contacting residues when `op_cif` is given, else the
        inter-domain-bridging / non-thiolate-cage site) -- and CALL its metal;
      * if ONE site, return it (single-site workflow; no invented second site).
    Returns {n_sites, mode, sensing_metal, sensing_element, structural_metals, chemistry, detail}.
    """
    import math
    atoms = _cif_atoms(fold_cif)
    metals = [a for a in atoms if a[0] == "HETATM" and a[1] in _METAL_ELEMENTS]
    donors = [a for a in atoms if a[0] == "ATOM" and a[1] in ("N", "O", "S")]
    if not metals:
        return {"n_metal_sites": 0, "mode": "no_metal", "sensing_metal": None}

    dna_seqs = None
    if op_cif:
        op_atoms = _cif_atoms(op_cif)
        dna = [a for a in op_atoms if a[2] in _DNA_COMPS]
        if dna:
            dna_seqs = set()
            for a in op_atoms:
                if a[0] == "ATOM" and a[4] is not None and any(
                        math.dist((a[5], a[6], a[7]), (d[5], d[6], d[7])) <= 10.0 for d in dna):
                    dna_seqs.add(a[4])

    from collections import Counter
    all_sites = []
    for m in metals:
        lig = [d for d in donors if math.dist((m[5], m[6], m[7]), (d[5], d[6], d[7])) <= _COORD_CUT]
        ncys = sum(1 for d in lig if d[2] == "CYS" and d[1] == "S")
        seqs = [d[4] for d in lig if d[4] is not None]
        near_dna = bool(dna_seqs) and any(s in dna_seqs for s in seqs)
        chains = Counter(d[3] for d in lig)
        prot_chain = chains.most_common(1)[0][0] if chains else m[3]
        all_sites.append({"element": m[1], "ncys": ncys, "nN": sum(d[1] == "N" for d in lig),
                          "nO": sum(d[1] == "O" for d in lig), "nS": sum(d[1] == "S" for d in lig),
                          "minseq": min(seqs) if seqs else None, "maxseq": max(seqs) if seqs else None,
                          "near_dna": near_dna, "prot_chain": prot_chain})
    # A homodimer repeats the same site set on each monomer -- classify ONE representative monomer
    # (the fullest), so `n_metal_sites` is per protomer and the single-site fallback is not tripped by
    # the partner chain's copy of the same site.
    by_mono: dict = {}
    for s in all_sites:
        by_mono.setdefault(s["prot_chain"], []).append(s)
    sites = max(by_mono.values(), key=len) if by_mono else all_sites
    if len(sites) == 1:
        s = sites[0]
        return {"n_metal_sites": 1, "mode": "single_site",
                "sensing_element": s["element"], "sensing_metal": _CANONICAL_ION.get(s["element"], s["element"]),
                "structural_metals": [], "chemistry": _chem(s), "detail": "one site; single-site workflow"}

    # >1 site: the regulatory site is the one coupled to the DNA-binding domain. Score each: a
    # DNA-contacting ligand (op fold) or an inter-domain span (ligands both N- and C-terminal of the
    # midpoint) marks the regulatory site; the compact thiolate cage is structural.
    span = max((s["maxseq"] or 0) for s in sites) or 150
    for s in sites:
        bridges = (s["minseq"] is not None and s["maxseq"] is not None
                   and s["minseq"] < 0.55 * span < s["maxseq"])
        s["reg_score"] = (2 if s["near_dna"] else 0) + (1 if bridges else 0) + (0.3 if s["ncys"] <= 1 else 0)
    order = sorted(sites, key=lambda s: -s["reg_score"])
    sensing, structural = order[0], order[1:]
    return {"n_metal_sites": len(sites), "mode": "multi_site",
            "sensing_element": sensing["element"],
            "sensing_metal": _CANONICAL_ION.get(sensing["element"], sensing["element"]),
            "structural_metals": [_CANONICAL_ION.get(s["element"], s["element"]) for s in structural],
            "chemistry": _chem(sensing),
            "detail": f"{len(sites)} sites; regulatory={_chem(sensing)}, "
                      f"structural={', '.join(_chem(s) for s in structural)}"}


def _chem(site: dict) -> str:
    if site["ncys"] >= 3:
        return f"{site['element']}·thiolate-cage(Cys{site['ncys']})"
    return f"{site['element']}·N/O(N{site['nN']}O{site['nO']}Cys{site['ncys']})"


def _structure_confirmation(fold_cif: str, family: str, op_cif: str | None) -> dict:
    """Merge structure-typed sensing-site evidence into the gate (see coordination_gate docstring)."""
    try:
        res = sensing_site_from_structure(fold_cif, family, op_cif)
    except Exception:
        return {}
    if not res.get("sensing_element"):
        return res if res else {}
    elem = res["sensing_element"]
    out = {"structure_sites": res, "has_site": True}
    shortlist = _SURROGATE_SHORTLIST.get(elem)
    implied = {}
    if shortlist:
        # a surrogate names the class + a Fe/Mn/Ni shortlist, never a single ion
        for ion in shortlist:
            implied[ion] = ["structure:regulatory_site(surrogate)"]
    else:
        implied[res["sensing_metal"]] = ["structure:regulatory_site"]
    # union with any sequence-motif implied ions already present is done by the caller wrapper via
    # dict.update on `ev`; here we only add the structure-implied set under its own key for provenance
    out["structure_implied_ions"] = implied
    return out


def _coordination_gate_impl(seq: str, family: str) -> dict:
    """The per-family gate branches. Use `coordination_gate`, which also attaches catalogue evidence."""
    seq = (seq or "").upper()
    n_cys = seq.count("C")
    fam = family or ""
    if "MerR" in fam:
        hits = cys_pairs(seq, _CXC)
        # the genuine MerR metal pair sits in the C-terminal dimerization helix (latter ~40% of the chain)
        c_term = [p for p in hits if p > 0.55 * len(seq)]
        return {"has_site": bool(c_term), "gate": "CX7-8C", "positions": hits,
                "c_terminal": bool(c_term), "n_cys": n_cys,
                "advisory": False if c_term else "no C-terminal CX(7-8)C -> treat as non-metal"}
    if "ArsR" in fam or "SmtB" in fam:
        # ArsR/SmtB has TWO metal-site chemistries: alpha3/alpha3N THIOL (Cys pairs; ArsR/CadC/SmtB) and
        # alpha5 HIS/carboxylate (CzrA His97/His100/Asp; NmtR His/Ni) -- the latter carries NO cysteines,
        # so a Cys-only gate would wrongly reject genuine Zn/Ni sensors. Check both.
        cys = cys_pairs(seq, _CYS_PAIR)
        n_his = seq.count("H")
        hxxh = bool(re.search(r"H..H", seq))            # alpha5 His pair signature
        thiol = bool(cys) and n_cys >= 2
        his_site = hxxh and n_his >= 2
        chem = "thiol(Cys)" if thiol else ("His/carboxylate" if his_site else "none")
        return {"has_site": bool(thiol or his_site), "gate": "ArsR-metal", "chemistry": chem,
                "positions": cys, "c_terminal": None, "n_cys": n_cys, "n_his": n_his, "hxxh": hxxh,
                "advisory": "ArsR sites are alpha3(Cys) OR alpha5(His/Asp); both chemistries checked"}
    if "NikR" in fam:
        # NikR's high-affinity Ni(II) site is square-planar His87/His89/Cys95 (+ His76' from the adjacent
        # subunit) in the C-terminal tetramerisation domain -- E. coli NikR numbering.
        #
        # THIS MOTIF IS FAMILY-SCOPED AND MUST STAY THAT WAY. An earlier version of this comment claimed
        # the His-x-His-x(4,7)-Cys spacing was "present in 85.5% of the vendored NikR member DB and in
        # <=3.2% of every other family's members (Fur 3.2%, MerR 1.1%, LysR 1.2%)". Re-measured against
        # the member DB this package actually ships (2026-09-02): NikR 85.4% -- that half was right --
        # but Fur 60.6%, LysR 8.7%, MerR 3.4%. The Fur figure was wrong by a factor of nineteen. It
        # changes no prediction, because this branch is only reached when the family is NikR, and that
        # is precisely the reason it survived unnoticed. Every prevalence in this module is now
        # generated and re-measured by `tools/build_metal_reference.py` into
        # `predictor/data/refs/metal_signatures.json`; see `predictor/structure/metal_motifs.py`.
        m = _NIKR_SITE.search(seq)
        return {"has_site": bool(m), "gate": "NikR-HxHxC", "chemistry": "His2/Cys (square-planar Ni)",
                "positions": [m.start()] if m else [], "c_terminal": None, "n_cys": n_cys,
                "n_his": seq.count("H"),
                "advisory": False if m else "no H-x-H-x(4-7)-C Ni site -> treat as non-metal"}
    if "CsoR" in fam or "FrmR" in fam:
        # The CsoR/RcnR/FrmR family shares one fold but splits on which metal-LIGATING residues survive:
        #   * CsoR (Cu(I))      -- Cys36, His61, Cys65 (trigonal S2N);
        #   * RcnR (Ni/Co)      -- His3 (N-term), Cys35, His60, His64;
        #   * FrmR (formaldehyde) -- Pro2 + Cys35 only; it has LOST His64 (Glu there) and the N-terminal His,
        #     and senses formaldehyde by a Pro2/Cys35 methylene bridge, not a metal (Osman et al. 2016,
        #     Nat Chem Biol 12:839; the gain-of-function FrmR(E64H) variant restores metal sensing).
        # So the discriminator is simply whether the family's His+Cys ligand set is still present. The rule
        # fires on 93.4% of the vendored CsoR member DB (E. coli RcnR: 5 His/1 Cys -> metal) and correctly
        # excludes E. coli FrmR (1 His/2 Cys). Scoped to this family: the same counts mean nothing elsewhere.
        n_his = seq.count("H")
        site = n_his >= 2 and n_cys >= 1
        return {"has_site": site, "gate": "CsoR-His/Cys", "chemistry": "His/Cys (Cu(I) or Ni/Co)" if site
                else "ligand set lost -> FrmR/aldehyde branch",
                "positions": cys_pairs(seq, _CYS_PAIR), "c_terminal": None, "n_cys": n_cys, "n_his": n_his,
                "advisory": False if site else "CsoR-family metal ligands absent -> non-metal (FrmR-type)"}
    if "Rrf2" in fam or "IscR" in fam or "NsrR" in fam:
        # Rrf2 is NOT a metal-ION sensor family: its members ligate an [Fe-S] CLUSTER (IscR Cys92/98/104,
        # NsrR Cys90/96/102) or sense thiol status (CymR), and respond to cluster occupancy / NO / cysteine.
        # `has_site` therefore reports False -- there is no metal-ion site to gate a Cu/Zn/Ni-style label on --
        # while `fes_triad` carries the positive redox evidence that `inducer_class` acts on.
        triad = fes_cluster_triad(seq)
        return {"has_site": False, "gate": "Rrf2-FeS", "chemistry": "[Fe-S] cluster (Cys triad)" if triad
                else "no Cys triad", "positions": [], "c_terminal": None, "n_cys": n_cys,
                "n_his": seq.count("H"), "fes_triad": triad,
                "advisory": "Rrf2 senses [Fe-S]/NO/thiol status -> redox, not a metal ion"}
    if "Fur" in fam:
        # Fur has TWO distinct sites (Rondon et al. Fig 4C): site 1 is STRUCTURAL and Cys-rich (`Cx2C`),
        # site 2 is the REGULATORY/allosteric metal-response site, a His/carboxylate `HH(H/D)H`. Gating a
        # *metal-response* claim on site 1's cysteines is therefore wrong on its own: the Mur (Mn(II))
        # clade carries site 2 in ~93% of members but `Cx2C` in only ~2%, so a `his and n_cys>=2` gate
        # rejects ~39/40 genuine Mn sensors. Accept EITHER chemistry (mirrors the ArsR dual-site gate).
        site2 = bool(re.search(r"HH[HD]H", seq))                  # regulatory His site (site 2)
        his_loose = bool(re.search(r"H.{1,3}H", seq))
        struct = his_loose and n_cys >= 2                          # site-1 Cys pair + His (Fur/Zur/PerR)
        chem = "His site-2" if site2 else ("His/Cys site-1" if struct else "none")
        return {"has_site": site2 or struct, "gate": "Fur-His/Cys", "chemistry": chem,
                "positions": cys_pairs(seq, _CYS_PAIR), "site2": site2,
                "c_terminal": None, "n_cys": n_cys,
                "advisory": "Fur is promiscuous -> advisory only; site1(Cys) OR site2(His) both accepted"}
    return {"has_site": None, "gate": "", "positions": [], "c_terminal": None, "n_cys": n_cys,
            "advisory": "no gate defined for this family"}


# --------------------------------------------------------------------------- Ligify gating (Phase 3.x)
# Families whose effector is an ORGANIC small molecule -> Ligify (operon-enzyme chemistry) is the right,
# and often only, inducer signal. TetR/AcrR is the archetype (tetracycline, bile acids, fatty acids, ...).
LIGAND_FAMILIES = ("TetR", "AcrR", "LysR", "GntR", "MarR", "SlyA", "AraC", "XylS", "IclR", "LacI", "GalR",
                   "LuxR", "FixJ", "SsuR", "VanR", "ArgR", "TrpR", "HxlR")
# Families that sense a METAL or a reactive species (RSS/ROS/RNS) -> the metabolite Ligify reads off a
# neighbouring enzyme is not the effector, so Ligify must be suppressed. 2026 drop added the new metal /
# metallocofactor sensors: NikR (Ni), CopY (Cu; the BlaI/MecI branch senses beta-lactam via signal
# transduction, also not a Ligify metabolite), Rrf2/IscR/NsrR ([Fe-S], NO, cysteine).
SENSOR_FAMILIES = ("ArsR", "SmtB", "MerR", "Fur", "Zur", "Nur", "DtxR", "MntR", "CsoR", "FrmR", "NiaR",
                   "NikR", "CopY")
# Families whose effector is a REACTIVE SPECIES or cofactor-occupancy signal rather than a metal ion. The
# 2026 drop initially filed Rrf2/IscR/NsrR under SENSOR_FAMILIES, which made every Rrf2 member classify as a
# metal sensor: the family prior fired before any redox evidence was consulted, and (worse) the resulting
# sensor_class="metal" tripped the metalloregulator guard in effector.inducer._finalize, actively demoting
# correct non-metal headlines. Rrf2 members ligate an [Fe-S] cluster (IscR, NsrR, RsrR) or read thiol status
# (CymR); none of them sense a free metal ion, so the family belongs in its own redox tier.
REDOX_FAMILIES = ("Rrf2", "IscR", "NsrR", "RsrR")

# Families where a detected metal site is, on the literature, more likely a STRUCTURAL cofactor than a
# regulatory one. Finding a site here is a real observation and must be reported as such -- it is the
# INTERPRETATION "therefore a metal sensor" that does not follow.
#
# GntR is the case that matters and it is the largest blind block we have (32 of 150 candidates).
# Rondon/Giedroc/Capdevila, Chem. Rev. metallostasis review, section 3.2.10:
#
#     "It has been proposed that ~70% of GntR proteins maintain a conserved, high binding site for
#      Zn(II). This metal binding site is buried at the bottom of a solvent accessible cavity, where it
#      facilitates the binding of small organic acids that harbor a carboxylate group ... The extent to
#      which metal binding is regulatory in this family of proteins or functions simply a structural
#      cofactor that enables molecular recognition of the cognate small molecule effector remains an
#      open question."
#
# Measured on run 6: MetalNet finds a site in 9 of 32 GntR candidates, and the residue sets it reports
# (E118/H164/H217/H239; H134/E175/H176/H198; E140/H144/H186) are His/Glu constellations in the
# C-terminal effector domain -- i.e. exactly the site the review describes. The site is real. Calling
# those regulators metal sensors on the strength of it is not supported.
STRUCTURAL_SITE_FAMILIES = {
    "GntR": ("~70% of GntR carry a conserved Zn site that is buried in the effector cavity and may "
             "serve to position a carboxylate-bearing organic acid rather than to sense metal "
             "(Chem. Rev. metallostasis review, 3.2.10 -- stated there as an open question)"),
}


def structural_site_caveat(family: str) -> str | None:
    """Why a metal site in this family may not be a regulatory one, or None if no caveat applies."""
    fam = family or ""
    for key, why in STRUCTURAL_SITE_FAMILIES.items():
        if key in fam:
            return why
    return None


def reactive_cys_cluster(seq: str, *, window: int = 16, min_cys: int = 3) -> bool:
    """>= `min_cys` cysteines within a `window`-residue span -> a compact reactive-Cys cluster. Used, like
    the metal gate, to tell a reactive-species sensor apart from an organic-ligand one."""
    pos = [i for i, a in enumerate((seq or "").upper()) if a == "C"]
    return any(pos[i + min_cys - 1] - pos[i] <= window for i in range(len(pos) - min_cys + 1))


# The [2Fe-2S] ferredoxin-like motif of the redox/reactive-species MerR sensors (SoxR C-X2-C-X-C-X(4-6)-C;
# also NsrR/IscR-adjacent). It is what tells a redox sensor apart from a genuine metal sensor: a benchmark
# against RegulonDB showed a plain Cys-density threshold cannot (Zur/CueR carry 4 Cys too), but this specific
# spacing matched SoxR and none of Cu/Zn/Fe/Mn sensors.
_2FE2S = re.compile(r"C.{2}C.C.{4,6}C")


def iron_sulfur_cluster(seq: str) -> bool:
    """SoxR-type [2Fe-2S] ferredoxin motif -> a reactive-species (RSS/ROS/RNS) sensor, not a metal-ion one."""
    return bool(_2FE2S.search((seq or "").upper()))


# The Rrf2-family [Fe-S] ligation triad: IscR Cys92/Cys98/Cys104, NsrR Cys90/Cys96/Cys102 -- i.e.
# C-x(4,6)-C-x(4,6)-C. Measured on the vendored SSN member DBs it is present in 74.0% of Rrf2 members and
# <=3.3% of every other family's (LysR 1.8%, CopY 3.3%, MerR 0.1%, Fur/DtxR/GntR/MarR/NikR/TetR 0.0%), and in
# Rrf2 it sits at a tightly conserved position (median 0.57 of the chain, IQR 0.55-0.58 -- exactly IscR's
# C92/162). It is used only to CORROBORATE a Rrf2-family redox call, never on its own: the same spacing does
# occur in Zur's structural Zn site, which is why the caller must scope it by family.
_FES_TRIAD = re.compile(r"C.{4,6}C.{4,6}C")


def fes_cluster_triad(seq: str) -> bool:
    """Rrf2-type [Fe-S] ligation triad C-x(4,6)-C-x(4,6)-C (IscR/NsrR). Family-scoped -- see `_FES_TRIAD`."""
    return bool(_FES_TRIAD.search((seq or "").upper()))


def metalnet_verdict(metalnet) -> bool | None:
    """Normalise whatever a caller passes for the MetalNet2 source to True / False / None.

    Accepts the `structure.metalnet.predict_site` result dict, a bare bool, or None. None means NO
    OPINION -- the engine is absent, the MSA was too shallow, or it could not score the sequence -- and
    must never be read as a negative: an abstention leaves the existing evidence untouched, whereas a
    false negative would actively demote a genuine metal sensor."""
    if metalnet is None or isinstance(metalnet, bool):
        return metalnet
    if isinstance(metalnet, dict):
        return metalnet.get("has_site") if metalnet.get("available") else None
    return None


def _metalnet_detail(metalnet) -> str:
    if not isinstance(metalnet, dict):
        return "MetalNet"
    n_res = metalnet.get("n_site_residues")
    res = ", ".join(metalnet.get("site_residues") or [])
    return f"MetalNet2 site: {n_res} co-evolving CHED residues ({res})" if n_res else "MetalNet2"


def inducer_class(seq: str, family: str, gate: dict | None = None, metalnet=None) -> tuple[str, str]:
    """Coarse effector class from offline signals: 'redox' (RSS/ROS/RNS via [2Fe-2S]/reactive-Cys), 'metal'
    (coordination site or metalloregulator family), or 'organic' (ligand family / no such site). Drives both
    the Ligify weight and the class-level inducer label.

    `metalnet` is the family-agnostic MetalNet2 verdict (see `structure.metalnet`). It is consulted only
    where our own chemistry is SILENT -- six of the twelve families, 82 of the 150 run-5 candidates -- and
    never to overturn a family-specific gate that actually ran. See rule 4b below."""
    g = gate if gate is not None else coordination_gate(seq, family)
    mn = metalnet_verdict(metalnet)
    fam = family or ""
    # Precedence is by STRENGTH OF EVIDENCE, not by convenience: a specific structural signature outranks a
    # family-scoped one, which outranks a bare family prior, which outranks the loose Cys-density heuristic.
    # 1. specific cofactor motif -- SoxR-type [2Fe-2S]
    if iron_sulfur_cluster(seq):
        return "redox", "[2Fe-2S] ferredoxin motif -> reactive-species (redox) sensor"
    # 2. redox families -- their effector is a reactive species / cofactor-occupancy signal, never a metal ion
    if any(f in fam for f in REDOX_FAMILIES):
        why = ("[Fe-S] Cys triad" if g.get("fes_triad") or fes_cluster_triad(seq) else "no Cys triad")
        return "redox", f"redox-sensor family ({fam}); {why} -> reactive-species / cluster-occupancy sensor"
    # 3. positive metal-coordination evidence
    if g.get("has_site") is True:
        return "metal", f"metal-coordination site ({g.get('gate') or 'gate'})"
    # 4. a sensor family whose gate ran and came back NEGATIVE. The gate exists precisely to catch the
    #    non-metal members inside a metalloregulator family -- the non-metal MerR clades (TnrA/GlnR/NmlR,
    #    plan SS6) and the CsoR-family FrmR/aldehyde branch. Falling through to the family prior below would
    #    discard that verdict and call them metal, which is the failure the gate was built to prevent.
    if g.get("has_site") is False and any(f in fam for f in SENSOR_FAMILIES):
        return "organic", (f"metalloregulator family ({fam}) but its coordination gate is NEGATIVE "
                           f"({g.get('advisory') or g.get('gate')}) -> non-metal member of a metal family")
    # 4b. MetalNet2, but ONLY where our own chemistry is silent. Deliberately placed after rule 4 and
    #     before the family prior: a family-specific gate that RAN carries more information about this
    #     protein than a family-agnostic classifier (the non-metal MerR clades and the CsoR-family FrmR
    #     branch are precisely the cases our gate exists to catch, and MetalNet cannot see them), while a
    #     SILENT gate carries none at all -- and silence is the majority case, six of twelve families.
    #     This is the rule that gives GntR/TetR/MarR/DtxR/CopY/LysR a structural opinion for the first time.
    if g.get("has_site") is None and mn is True:
        caveat = structural_site_caveat(fam)
        if caveat:
            # The site is real and is reported; what does not follow is "therefore a metal sensor".
            # Falling through leaves the family's own rule (organic-ligand, below) to set the class,
            # and `_finalize` states the site and this caveat side by side.
            pass
        else:
            return "metal", _metalnet_detail(metalnet) + " (family coordination gate silent)"
    # 5. family prior -- gate silent (no branch, or a dispersed site the sequence cannot see, e.g. MntR).
    #    A NEGATIVE MetalNet does not overturn this: a metalloregulator family whose site the sequence gate
    #    cannot see (MntR's dispersed His/Asp carboxylate site) is exactly the case MetalNet is also
    #    likeliest to miss, so a miss there is weak evidence. It is recorded as a dissent, not acted on.
    if any(f in fam for f in SENSOR_FAMILIES):
        mn_note = ("" if mn is None else
                   "; MetalNet corroborates" if mn else "; MetalNet found no site (dissent, not acted on)")
        return "metal", f"metalloregulator family ({fam}); coordination gate silent{mn_note}"
    # 6. organic-ligand family. Checked BEFORE the reactive-Cys heuristic: family assignment is real evidence,
    #    whereas ">=3 Cys in a 16-residue window" fires on 17.9% of MarR and 58.1% of ArsR members, so letting
    #    it outrank the family call mislabels ordinary ligand sensors as redox (E. coli MarR, whose Cys47/51/54
    #    are structural -- its effector is salicylate).
    if any(f in fam for f in LIGAND_FAMILIES):
        # With MetalNet wired in, an organic call for these families is no longer "our gate had no
        # chemistry here" but "neither method found a site" -- which is what makes the organic-family
        # guard in effector.inducer._finalize an assertion rather than an absence of evidence.
        if mn is False:
            return "organic", (f"organic-ligand family ({fam}); neither the coordination gate nor "
                               f"MetalNet found a metal site")
        if mn is True and structural_site_caveat(fam):
            return "organic", (f"organic-ligand family ({fam}); MetalNet DOES find a metal site, but "
                               f"in this family such a site is likely structural rather than "
                               f"regulatory -- reported, not treated as a sensor call")
        return "organic", f"organic-ligand family ({fam}), no metal/redox site"
    # 7. last-resort heuristic for an unknown family
    if reactive_cys_cluster(seq):
        return "redox", "reactive-Cys cluster (unknown family) -> reactive-species sensor"
    return "organic", ("no metal/redox signature" if mn is not False
                       else "no metal/redox signature; MetalNet found no site either")


def ligify_weight(seq: str, family: str, gate: dict | None = None, metalnet=None) -> tuple[float, str]:
    """(weight in [0,1], reason) for the Ligify inducer source, gated by structure then family.

    Ligify infers an ORGANIC effector from the chemistry of a neighbouring operon: right for ligand-sensing
    families, wrong for metal / reactive-species sensors. Suppressed (0.1) when the class is metal or redox;
    kept high (1.0) for a ligand family; family sets the prior when the structural signal is silent."""
    g = gate if gate is not None else coordination_gate(seq, family)
    mn = metalnet_verdict(metalnet)
    cls, why = inducer_class(seq, family, g, metalnet)
    # A structural-site family keeps its organic weight even when a site is found. Suppressing operon
    # chemistry x0.1 on a site the literature calls structural actively removes the best evidence this
    # regulator has: measured on run 6, all 9 GntR candidates with a MetalNet site had their Ligify
    # call down-weighted x0.1 by a Zn site that, per the review, may exist to position the organic acid.
    if structural_site_caveat(family) and g.get("has_site") is not True:
        return 1.00, (f"organic-ligand family ({family}); any metal site here is likely structural "
                      f"-> Ligify remains the signal")
    if cls in ("metal", "redox"):
        # a metalloregulator family with a SILENT gate and no MetalNet site is only advisory -> suppress
        # less hard. A site CONFIRMED by either method is direct evidence and suppresses fully.
        if (cls == "metal" and g.get("has_site") is not True and mn is not True
                and not iron_sulfur_cluster(seq)):
            return 0.30, f"{why} -> advisory"
        return 0.10, f"{why} -> not an organic ligand"
    fam = family or ""
    if any(f in fam for f in LIGAND_FAMILIES):
        return 1.00, f"organic-ligand family ({fam}) with no metal/redox site -> Ligify is the signal"
    if g.get("has_site") is False and any(f in fam for f in SENSOR_FAMILIES):
        # a NON-METAL member of a metalloregulator family (the MerR BmrR/MtaN/GyrI-like effector-binding
        # clades, ~61% of the MerR member DB; the CsoR-family FrmR branch). Its effector is an organic
        # molecule, so operon chemistry is the best remaining signal -- but the family still counsels
        # caution, hence not the full 1.00 an established ligand family gets.
        return 0.85, f"non-metal member of {fam} (coordination gate negative) -> lean on Ligify"
    return 0.70, "unknown family, no metal/redox signature -> lean on Ligify"


def _self_test() -> None:
    from predictor.annotate.af3_msa import parse_a3m, _op_msa

    def qseq(tf):
        a = _op_msa(tf)
        return "".join(c for c in next(parse_a3m(a))["seq"] if c.isalpha()).upper() if a else None

    cases = [("CueR_Ecoli", "MerR"), ("CzrA_Saureus", "ArsR/SmtB"), ("NmtR_Mtb", "ArsR/SmtB")]
    res = {}
    for tf, fam in cases:
        s = qseq(tf)
        if not s:
            continue
        g = coordination_gate(s, fam)
        res[tf] = g
        print(f"{tf:14s} ({fam:9s}) gate={g['gate']:10s} has_site={g['has_site']} "
              f"chem={g.get('chemistry','-'):16s} n_cys={g['n_cys']} n_his={g.get('n_his','-')}")
    # clear cases must pass: CueR (thiol CX7-8C) and CzrA (alpha5 His). NmtR's Ni site is formed by
    # sequence-DISPERSED residues -> the sequence gate cannot see it; only a structure can, which is a
    # known and documented limitation and the reason `sensing_site_from_structure` exists.
    if "CueR_Ecoli" in res:
        assert res["CueR_Ecoli"]["has_site"], "CueR thiol gate must pass"
    if "CzrA_Saureus" in res:
        assert res["CzrA_Saureus"]["has_site"], "CzrA alpha5-His gate must pass"
    if "NmtR_Mtb" in res and not res["NmtR_Mtb"]["has_site"]:
        print("  note: NmtR sequence-gate negative -> dispersed His/Asp Ni site; needs a fold")
    # a synthetic non-metal MerR (no C-terminal Cys pair) must fail the gate
    nonmetal = "M" + "A" * 60 + "L" * 60          # no cysteines at all
    g = coordination_gate(nonmetal, "MerR")
    assert g["has_site"] is False, "MerR without CX7-8C must fail the metal gate"

    # Fur regression: the Mur (Mn(II)) clade carries the REGULATORY His site-2 but no structural Cys
    # site-1, so a Cys-requiring gate rejected ~39/40 real Mn sensors. Site 2 alone must pass.
    mur_like = "M" + "A" * 40 + "HHDH" + "A" * 40          # site 2 present, zero cysteines
    g = coordination_gate(mur_like, "Fur")
    assert g["has_site"] and g["chemistry"] == "His site-2", f"Fur site-2-only (Mur) must pass: {g}"
    # real E. coli Fur (site 1 + site 2) must also pass
    g = coordination_gate("MTDNNTALKKAGLKVTLPRLKILEVLQEPDNHHVSAEDLYKRLIDMGEEIGLATVYRVLNQFDDAGIVTRHN"
                          "FEGGKSVFELTQQHHHDHLICLDCGKVIEFSDDSIEARQREIAAKHGIRLTNHSLYLYGHCAEGDCREDEHAHEGK", "Fur")
    assert g["has_site"], "E. coli Fur must pass the Fur gate"
    # a Fur-family sequence with neither chemistry must fail
    g = coordination_gate("M" + "A" * 100, "Fur")
    assert g["has_site"] is False, "Fur with no His/Cys site must fail the gate"
    # a synthetic MerR metal sensor with a C-terminal CX7C must pass
    metal = "M" + "A" * 70 + "C" + "X" * 7 + "C" + "A" * 10
    assert coordination_gate(metal, "MerR")["has_site"] is True, "C-terminal CX7C must pass the gate"
    print("OK: metal-coordination gate detects CX7-8C / Cys pairs and rejects non-metal MerR.")

    # --- Ligify gating: organic families trusted, metal/redox sensors suppressed ---
    assert ligify_weight("M" + "ARQLGID" * 8, "TetR/AcrR")[0] == 1.0, "TetR must trust Ligify"
    assert ligify_weight(metal, "MerR")[0] == 0.10, "MerR metal site must suppress Ligify"
    # A reactive-Cys cluster suppresses Ligify only for an UNKNOWN family. It deliberately no longer
    # overrides an organic-ligand family call: ">=3 Cys in 16 residues" fires on 17.9% of MarR and 4.7% of
    # LysR members (measured on the vendored SSN member DBs), so letting it win there mislabels ordinary
    # ligand sensors as redox and suppresses the operon-chemistry signal that is CORRECT for them -- the
    # E. coli MarR case, whose Cys47/51/54 are structural and whose effector is salicylate. Genuine redox
    # sensors are still caught by the two stronger routes (the [2Fe-2S] motif and REDOX_FAMILIES), both of
    # which are checked first, so nothing real is lost.
    assert ligify_weight("M" + "A" * 30 + "CAACAAC" + "A" * 40, "")[0] == 0.10, \
        "reactive-Cys cluster in an UNKNOWN family must suppress Ligify (RSS/ROS/RNS sensor)"
    assert ligify_weight("M" + "A" * 30 + "CAACAAC" + "A" * 40, "LysR-type (LTTR)")[0] == 1.00, \
        "organic-ligand family outranks the loose Cys-density heuristic -> Ligify stays trusted"
    # A metalloregulator family whose gate is SILENT (no branch / dispersed site, e.g. DtxR-MntR) keeps the
    # metal prior and only-advisory Ligify. One whose gate ran and said NO is a non-metal member of that
    # family (MerR's BmrR/MtaN clades, ~61% of the member DB) -> Ligify becomes the best remaining signal.
    assert ligify_weight("M" + "A" * 120, "DtxR/MntR")[0] == 0.30, "gate SILENT -> advisory Ligify"
    assert ligify_weight("M" + "A" * 120, "Fur")[0] == 0.85, "gate NEGATIVE -> non-metal member, use Ligify"
    assert not reactive_cys_cluster("M" + "C" + "A" * 50 + "C"), "2 distant Cys must NOT be a cluster"
    # SoxR-type [2Fe-2S] ferredoxin motif -> redox class, separable from metal-Cys sensors (RegulonDB iter.)
    soxr_like = "M" + "A" * 80 + "CAACACAAAAAC" + "A" * 20
    assert iron_sulfur_cluster(soxr_like), "[2Fe-2S] motif must be detected"
    assert inducer_class(soxr_like, "MerR")[0] == "redox", "SoxR-type Fe-S -> redox, not metal"
    assert inducer_class(metal, "MerR")[0] == "metal", "CX7-8C metal site -> metal"
    assert inducer_class("M" + "ARQLGID" * 8, "TetR/AcrR")[0] == "organic", "TetR -> organic"
    print("OK: Ligify gating trusts organic families and suppresses metal/redox sensors.")

    # --- 2026-drop families: the gate must produce EVIDENCE, not just inherit the family prior ---------
    # Rrf2 ([Fe-S]/NO/thiol status) is a REDOX family, not a metal one. Filing it under SENSOR_FAMILIES made
    # every member classify as metal, which additionally tripped the metalloregulator guard in
    # effector.inducer._finalize and demoted correct non-metal headlines. E. coli IscR C92/C98/C104.
    iscr = "M" + "A" * 88 + "CAAAAACAAAAAC" + "A" * 60
    assert fes_cluster_triad(iscr), "IscR-type C-x(4,6)-C-x(4,6)-C triad must be detected"
    assert inducer_class(iscr, "Rrf2")[0] == "redox", "Rrf2 family -> redox, never metal"
    assert coordination_gate(iscr, "Rrf2")["has_site"] is False, "Rrf2 has no metal-ION site to gate on"
    assert ligify_weight(iscr, "Rrf2")[0] == 0.10, "Rrf2 redox sensor must suppress Ligify"

    # NikR: square-planar His87/His89/Cys95 Ni site -> a CONFIRMED site, not the bare family prior.
    nikr = "M" + "A" * 84 + "HVHINHDDC" + "A" * 38
    g = coordination_gate(nikr, "NikR")
    assert g["has_site"] and g["gate"] == "NikR-HxHxC", f"NikR H-x-H-x(4-7)-C site must pass: {g}"
    assert inducer_class(nikr, "NikR")[1].startswith("metal-coordination site"), \
        "NikR must classify on site evidence, not the family prior"
    assert coordination_gate("M" + "A" * 130, "NikR")["has_site"] is False, "NikR without the site must fail"

    # CsoR/FrmR split on whether the family's His+Cys ligand set survives: RcnR (Ni/Co) keeps it, FrmR
    # (formaldehyde) has lost His64 + the N-terminal His and senses via Pro2/Cys35 (Osman 2016).
    rcnr = "MSHTIRDKQKLKARASKIQGQVVALKKMLDEPHECAAVLQQIAAIRGAVNGLMREVIKGHLTEHIVHQGDELKREEDLDVVLKVLDSYIK"
    frmr = "MPSTPEEKKKVLTRVRRIRGQIDALERSLEGDAECRAILQQIAAVRGAANGLMAEVLESHIRETFDRNDCYSREVSQSVDDTIELVRAYLK"
    assert coordination_gate(rcnr, "CsoR/FrmR")["has_site"] is True, "RcnR keeps the CsoR His/Cys ligands"
    assert coordination_gate(frmr, "CsoR/FrmR")["has_site"] is False, "FrmR has lost them -> non-metal"
    assert inducer_class(rcnr, "CsoR/FrmR")[0] == "metal", "RcnR -> metal (Ni/Co)"
    assert inducer_class(frmr, "CsoR/FrmR")[0] == "organic", "FrmR -> non-metal (formaldehyde)"

    # A metalloregulator family whose gate ran and said NO must not be rescued by the family prior -- the
    # whole point of the gate (plan SS6: the non-metal MerR clades TnrA/GlnR/NmlR).
    assert inducer_class(nonmetal, "MerR")[0] == "organic", \
        "a NEGATIVE coordination gate must outrank the metalloregulator family prior"
    print("OK: 2026-drop families (Rrf2 redox, NikR Ni site, CsoR/FrmR split) gate on evidence.")

    # --- MetalNet2 (2026-08): the family-AGNOSTIC gate, consulted only where our chemistry is silent ---
    # The coordination gate has branches for six of twelve families; for GntR/TetR/MarR/DtxR/CopY/LysR it
    # is silent on every member (82 of the 150 run-5 candidates). These assertions pin the precedence.
    mn_yes = {"available": True, "has_site": True, "n_site_residues": 4,
              "site_residues": ["C12", "C15", "H40", "C44"], "metal_type_top": "ZN"}
    mn_no = {"available": True, "has_site": False, "n_site_residues": 0}
    mn_abstain = {"available": False, "has_site": None, "reason": "engine unavailable"}
    assert metalnet_verdict(mn_yes) is True and metalnet_verdict(mn_no) is False
    assert metalnet_verdict(mn_abstain) is None, "an unavailable engine must abstain, not deny"
    assert metalnet_verdict(None) is None and metalnet_verdict(True) is True

    # A silent family where a MetalNet site DOES establish the metal class -- the 82-candidate gap
    # MetalNet was wired in to close. TetR stands for the general case.
    tetr = "M" + "ARQLGID" * 12                       # organic-ligand family, gate has no branch
    assert coordination_gate(tetr, "TetR/AcrR")["has_site"] is None, "TetR: our gate must be SILENT"
    assert inducer_class(tetr, "TetR/AcrR")[0] == "organic", "TetR alone stays organic"
    assert inducer_class(tetr, "TetR/AcrR", None, mn_yes)[0] == "metal", \
        "a MetalNet site where our chemistry is silent must establish the metal class (the 82-candidate gap)"
    assert "MetalNet" in inducer_class(tetr, "TetR/AcrR", None, mn_yes)[1], \
        "the reason must name the source"

    # GntR is the deliberate EXCEPTION, and this is what stops it regressing silently. ~70% of GntR
    # carry a conserved Zn site the literature says may be a structural cofactor for recognising the
    # organic acid rather than a sensed metal (STRUCTURAL_SITE_FAMILIES). GntR is also the largest
    # blind block we have -- 32 of 150 candidates, 9 with a MetalNet site -- so reading that site as a
    # sensor call would manufacture ~9 metalloregulators from chemistry that may not be regulatory.
    # The site is STATED in the report, never acted on.
    gntr = "M" + "ARQLGID" * 12
    assert coordination_gate(gntr, "GntR")["has_site"] is None, "GntR: our gate must be SILENT"
    assert inducer_class(gntr, "GntR")[0] == "organic", "GntR alone stays organic (unchanged behaviour)"
    assert inducer_class(gntr, "GntR", None, mn_yes)[0] == "organic", \
        "a GntR metal SITE must not become a metal CALL -- it is likely structural (see the review)"
    assert structural_site_caveat("GntR"), "and the caveat must be there to say so in the report"
    assert structural_site_caveat("TetR/AcrR") is None, "the exception is GntR, not every silent family"
    assert inducer_class(gntr, "GntR", None, mn_no)[0] == "organic", "MetalNet negative -> stays organic"
    assert "neither" in inducer_class(gntr, "GntR", None, mn_no)[1], \
        "with both methods run, the organic call must read as 'neither fired', not as absence of evidence"
    assert inducer_class(gntr, "GntR", None, mn_abstain)[0] == "organic", "abstention changes nothing"

    # MetalNet must NOT overturn a family-specific gate that actually ran -- in either direction. Those
    # verdicts encode chemistry MetalNet cannot see (the non-metal MerR clades; the FrmR branch).
    assert inducer_class(nonmetal, "MerR", None, mn_yes)[0] == "organic", \
        "a NEGATIVE family gate outranks a positive MetalNet"
    assert inducer_class(frmr, "CsoR/FrmR", None, mn_yes)[0] == "organic", "FrmR stays non-metal"
    assert inducer_class(metal, "MerR", None, mn_no)[0] == "metal", \
        "a POSITIVE family gate outranks a negative MetalNet"
    assert inducer_class(iscr, "Rrf2", None, mn_yes)[0] == "redox", \
        "a redox family stays redox: an [Fe-S] cluster is not a sensed metal ion"
    # a metalloregulator family whose gate is silent keeps its metal prior even if MetalNet misses --
    # MntR's dispersed carboxylate site is exactly what both methods are weakest at
    assert inducer_class("M" + "A" * 120, "DtxR/MntR", None, mn_no)[0] == "metal", \
        "a negative MetalNet must not demote a metalloregulator family (dissent is recorded, not acted on)"
    assert "dissent" in inducer_class("M" + "A" * 120, "DtxR/MntR", None, mn_no)[1]

    # Ligify weighting follows the same evidence: a MetalNet-confirmed site suppresses operon chemistry
    # as hard as our own gate does, instead of leaving it at the "advisory" 0.30.
    assert ligify_weight(tetr, "TetR/AcrR")[0] == 1.00, "TetR without MetalNet keeps Ligify (unchanged)"
    assert ligify_weight(tetr, "TetR/AcrR", metalnet=mn_yes)[0] == 0.10, \
        "a MetalNet-confirmed metal site must suppress the organic operon-chemistry source"
    # ...but NOT in GntR, and this is the second half of the same exception. If a GntR Zn site is a
    # structural cofactor for binding the organic acid, then suppressing the operon-chemistry source
    # on the strength of it would silence the evidence most likely to be RIGHT about this protein --
    # turning a caveat into an active harm rather than merely a missed opportunity.
    assert ligify_weight(gntr, "GntR")[0] == 1.00, "GntR without MetalNet keeps Ligify (unchanged)"
    assert ligify_weight(gntr, "GntR", metalnet=mn_yes)[0] == 1.00, \
        "a likely-STRUCTURAL GntR site must not suppress the organic source it exists to recognise"
    assert ligify_weight("M" + "A" * 120, "DtxR/MntR", metalnet=mn_yes)[0] == 0.10, \
        "MetalNet confirming a silent-gate metalloregulator upgrades advisory (0.30) to suppressed (0.10)"
    print("OK: MetalNet2 gates the 6 families our chemistry has no branch for, and never overrides one "
          "that ran.")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--seq"); ap.add_argument("--family", default="MerR")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args(argv)
    if a.self_test or not a.seq:
        _self_test()
        return
    s = a.seq
    if Path(a.seq).exists():
        s = "".join(l.strip() for l in Path(a.seq).read_text().splitlines() if not l.startswith(">"))
    print(coordination_gate(s, a.family))


if __name__ == "__main__":
    main()
