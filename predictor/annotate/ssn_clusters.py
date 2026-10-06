"""
ssn_clusters.py -- SSN (Sequence Similarity Network) cluster bookkeeping (Phase 0).

The Rondon/Antelo/Capdevila preprint clusters the ArsR, MerR and Fur metalloregulator families into
*isofunctional* groups (same inducer). The cluster FASTAs live under `predictor/data/ssn/clusters/`. This module
turns those raw FASTAs into the artifacts the rest of the feature consumes:

  * an **accession -> cluster** map (so a homolog hit can be tagged with its cluster),
  * a **cluster -> inducer** table with a *coordination-gate* flag (the metal label is only
    trustworthy when the coordination signature is present -- see the gate in `structure/metal_site`),
  * a **recognition-helix locator** (MerR `IRYYExxG`/`TLRYE`) used by the Test-1 pre-check and,
    later, as the Folddisco motif query.

Design notes
------------
* Header formats differ by family (verified against the data):
    MerR  `>A0A014LWC1 14|Erwinia mallotivora.|PF00376,PF09278`  -> acc = first token
    ArsR  `>tr|A0A011QX74|A0A011QX74_9LACT ... OS=.. OX=.. GN=..` -> acc = middle of tr|ACC|ENTRY
    Fur   uses BOTH styles (most clusters are MerR-style; mtZur is tr|ACC|ENTRY) -> parser handles both.
* Cluster -> inducer is CURATED here (the filenames carry a representative but not a machine-readable
  inducer). MerR labels come from the preprint (plan SS2.2); ArsR c1-c10 are **not inducer-labeled
  yet** (open task) -> inducer=None, role="unknown", and the gate is advisory only.
* Fur (added from the Rondon SSN drop) covers the 4 subgroups the preprint resolves: Fur/Mur (Fe/Mn),
  Zur/Nur (Zn/Ni), PerR (redox) and Irr (heme). Six labels are ANCHOR-VERIFIED -- the representative's
  own UniProt accession is a member of its cluster FASTA (ecFur P0A9A9, ecZur P0AC51, bsFur P54574,
  bsZur P54479, bsPerR P71086, cjFur P0C631).
* No external tool is needed for any function here. Assigning a *novel* sequence to a cluster is done
  by `assign_cluster_by_hits()` (pure: takes homolog hits the caller already computed), keeping the
  MMseqs/HMMER dependency at the call site (Phase 1).

Run `python -m predictor.annotate.ssn_clusters` for a self-test (parses the packaged FASTAs).
"""
from __future__ import annotations

import json
import re
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from predictor import resources

_REPO = Path(__file__).resolve().parents[2]
SSN_DIR = resources.SSN_DIR
_DATABASE = resources.SSN_DATABASE
_CACHE = resources.cache_path("ssn")
_ACC_MAP_PACKAGED = _DATABASE / "accession_to_cluster.json"
_ACC_MAP_CACHE = _CACHE / "accession_to_cluster.json"

# family canonical names (match family_db / FAMILY_PRIORS / operators_by_family spelling)
MERR = "MerR"
ARSR = "ArsR/SmtB"
FUR = "Fur"
# 2026 SSN drop: the remaining nine metalloregulator families now have cluster-level FASTAs. The canonical
# names MUST match packaged `pfam_to_family.json` so the evidence axes agree on one spelling. (The
# optional structure axis, `annotate/family_struct`, is not shipped in this release; if it is added,
# its reference-fold keys must use these same names.)
# family-classification axes (HMM / blastp / structure) and the SSN cluster layer all agree on spelling.
COPY = "CopY"
CSOR = "CsoR/FrmR"
DTXR = "DtxR/MntR"
GNTR = "GntR"
LYSR = "LysR-type (LTTR)"
MARR = "MarR/SlyA"
NIKR = "NikR"
RRF2 = "Rrf2"
TETR = "TetR/AcrR"

#: every family with an SSN cluster DB on disk. Callers that used to hardcode `[MERR, ARSR]` (or a
#: `"MerR" if fam == MERR else "ArsR"` binary) must use this + `family_tag()` instead, or a third
#: family silently lands in the second family's cache file.
FAMILIES = (MERR, ARSR, FUR, COPY, CSOR, DTXR, GNTR, LYSR, MARR, NIKR, RRF2, TETR)

#: filesystem-safe tag per family -- the single source of truth for cache/DB filenames
#: (`data/ssn/database/ssn_members_<tag>.fasta`, `curated_binders_<tag>.json`, ...). ArsR/SmtB -> "ArsR"
#: keeps the pre-existing on-disk names valid.
_TAGS = {MERR: "MerR", ARSR: "ArsR", FUR: "Fur",
         COPY: "CopY", CSOR: "CsoR", DTXR: "DtxR", GNTR: "GntR", LYSR: "LysR",
         MARR: "MarR", NIKR: "NikR", RRF2: "Rrf2", TETR: "TetR"}


def family_tag(family: str) -> str:
    """Filesystem-safe tag for `family` ('ArsR/SmtB' -> 'ArsR'). Unknown families degrade to a
    slugified name rather than silently colliding with another family's files."""
    if family in _TAGS:
        return _TAGS[family]
    return re.sub(r"[^A-Za-z0-9]+", "_", family or "unknown").strip("_") or "unknown"

# --------------------------------------------------------------------------- curated cluster table
# inducer ion notation matches the pipeline (CONFIGS effector="Zn2+", af3_export.map_effector_to_ion).
# role in {metal, non-metal, redox, multidrug, unknown}; gate names the coordination signature to
# verify before *trusting* a metal label (empty => no metal claim to gate).
@dataclass(frozen=True)
class ClusterInfo:
    cluster_id: str            # stable id, e.g. "MerR_9" / "ArsR_c1"
    family: str                # MERR / ARSR
    representative: str        # "ecZntR" / "" (from filename)
    inducer: str | None        # "Zn2+" / "nitrogen" / None(unknown)
    role: str                  # metal | non-metal | redox | multidrug | unknown
    gate: str                  # coordination signature to check, "" if none
    fasta: Path
    n_seqs: int = 0

    @property
    def is_metal(self) -> bool:
        return self.role == "metal"

    # ----------------------------------------------------------------- anchor set (vendored KB)
    # `representative` above is ONE hand-picked protein whose inducer became this cluster's label.
    # The anchor SET is every protein in the cluster with independently known biology, so a label can
    # finally carry support and disagreement instead of resting on a single accession. Loaded lazily
    # from packaged `data/refs/family_kb/`; absent files mean "no anchors known", never an error.
    @property
    def anchor_set(self):
        from predictor.annotate import family_kb
        return family_kb.anchor_set(self.cluster_id)

    @property
    def anchors(self) -> tuple:
        s = self.anchor_set
        return s.anchors if s else ()

    @property
    def n_independent_anchors(self) -> int:
        """Distinct SEQUENCES backing this cluster's label, not rows.

        The census holds 64 identical TrEMBL copies of E. coli CueR; a row count would report
        overwhelming support for one piece of evidence.
        """
        s = self.anchor_set
        return s.n_independent_anchors if s else 0

    @property
    def anchor_conflict(self) -> bool:
        """True when this cluster's anchors disagree about the inducer.

        Nine clusters do, and every one is real biology -- `ArsR_c2` carries Zn(II), Ni(II) (NmtR)
        AND Pb(II) (CadC) under a flat `Zn2+` label; `CsoR_ecFrmR` carries formaldehyde and Co2+
        because `ecRcnR` is nested inside it. A caller that needs one answer should treat these as
        low-confidence rather than trusting `inducer`.
        """
        s = self.anchor_set
        return bool(s and s.conflict)

    @property
    def label_confidence(self) -> str:
        """One of family_kb.TIERS (`high`, `tentative`, `mixed-anchor`, `unknown`, `structural-only`)."""
        s = self.anchor_set
        return s.tier if s else "unknown"


# MerR clusters: representative + inducer from the preprint (plan SS2.2 + Fig 3).
# Keyed by cluster integer parsed from the filename `MerR_cluster_<N>_<rep>.fasta`.
# The preprint's own curation sheet (`Table S1_And_Others.xlsx`, sheet "MerR": UniProt -> Cluster ->
# ligand) is the authority for these labels; each row below was cross-checked by looking the curated
# accession up in the cluster FASTA it is supposed to belong to (all 10 resolved consistently).
_MERR_TABLE = {
    1:  ("ppCadR", "Cd2+",          "metal",     "CX7-8C"),
    # Cluster 2 is deliberately POLYSPECIFIC. Two lines of evidence disagree about its representative and
    # neither is strong enough to overrule the other: the FASTA ships as `..._2_ecCadR.fasta` (Cd), while
    # the authors' curation sheet assigns the cluster exactly one characterized member, PbrR691 (Q1LKZ5,
    # Cupriavidus metallidurans, Pb(II)) -- verifiably in this FASTA, and their logo set carries
    # `cmPbrR_logo.svg`. CadR and PbrR are adjacent clades reading the same chemistry (soft thiophilic
    # divalents on the same C-terminal CX7-8C pair), so at this alignment score the SSN plausibly cannot
    # separate them and the members may genuinely be promiscuous. Labelling it `Cd2+/Pb2+` is the honest
    # call and costs nothing downstream: both ions map to the same AF3 surrogate (ZN), and
    # effector.inducer already prefers a query's nearest CURATED member ion over a lumped cluster label.
    # The legacy stem is kept so the loader's cluster-number parse is unaffected.
    2:  ("CadR/PbrR691", "Cd2+/Pb2+", "metal",   "CX7-8C"),
    3:  ("hiNmlR", "NO/RSS",        "non-metal", ""),
    4:  ("ecSoxR", "redox[2Fe-2S]", "redox",     ""),
    5:  ("bsTnrA", "nitrogen",      "non-metal", ""),
    6:  ("bcMerR1","Hg2+",          "metal",     "CX7-8C"),
    # Cluster 7 was previously "confirmed missing" -- it is the family's LARGEST cluster and simply was
    # not in the labelled drop. Recovered from the size-ranked partition (`SSN_MerR/MerR_FASTAs/
    # cluster_1.fasta`, 8942 seqs), identified by all 12 of the curation sheet's cluster-7 accessions
    # landing in it (CueR P0A9G4 Cu(I)/Ag(I), GolS A0A5Y3ZRF9 Au(I), hmrR, ...) and none elsewhere.
    # Polyspecific in the monovalent-coinage direction (Cu(I) dominant; Au(I)/Ag(I) minor).
    7:  ("ecCueR", "Cu+",           "metal",     "CX7-8C"),
    8:  ("bsYfmP", "multidrug",     "multidrug", ""),
    9:  ("ecZntR", "Zn2+",          "metal",     "CX7-8C"),
    10: ("bsGlnR", "nitrogen",      "non-metal", ""),
}

# ArsR clusters c1-c10: inducers assigned by VALIDATED-TF ANCHORS only.
# Keyword mining of UniProt descriptions is unreliable here (the family name "ArsR/arsenical" confounds it,
# and the residual signal is weak + Zn-diffuse), so only clusters that catch a characterized sensor are
# labelled; the rest stay unknown pending experimental anchors. gate = ArsR dual alpha3-thiol/alpha5-His.
_ARSR_GATE = "ArsR-metal"
_ARSR_TABLE = {
    "c2": ("CzrA/SmtB/ZiaR", "Zn2+", "metal"),       # HIGH: 3 independent Zn sensors all land here (support 1.0)
    "c7": ("NmtR(tentative)", "Ni2+", "metal"),      # TENTATIVE: single anchor (NmtR/Ni), and c7 also contains
                                                     # YgaV (thiol/RSS, NOT Ni) -> heterogeneous (retest finding).
}

# Fur clusters, keyed by the FASTA stem (the representative regulator). The preprint resolves the family
# into 4 subgroups and these 14 clusters populate all of them:
#   Fur/Mur  Fe(II)/Mn(II)   -- ecFur, bsFur, cjFur, baMur
#   Zur/Nur  Zn(II)/Ni(II)   -- ecZur, bsZur, mtZur, scZur, mtFurB, scNur
#   PerR     redox (H2O2)    -- bsPerR, liPerR, bbBosR
#   Irr      heme            -- bjIrr
# [ANCHOR] marks a label verified by finding the representative's own UniProt accession inside the
# cluster FASTA. PerR/Irr/BosR are "metalloforms that do not strictly respond to the metal" (the metal is
# redox-active or replaced by heme), so per the MerR convention they carry NO metal gate.
_FUR_GATE = "Fur-His/Cys"          # the gate name structure/metal_site.coordination_gate() returns
_FUR_TABLE = {
    "ecFur":  ("ecFur",  "Fe2+",        "metal",     _FUR_GATE),   # [ANCHOR] P0A9A9
    "bsFur":  ("bsFur",  "Fe2+",        "metal",     _FUR_GATE),   # [ANCHOR] P54574
    "cjFur":  ("cjFur",  "Fe2+",        "metal",     _FUR_GATE),   # [ANCHOR] P0C631
    "baMur":  ("baMur",  "Mn2+",        "metal",     _FUR_GATE),   # Mur clade, mostly Rhizobiales
    "ecZur":  ("ecZur",  "Zn2+",        "metal",     _FUR_GATE),   # [ANCHOR] P0AC51
    "bsZur":  ("bsZur",  "Zn2+",        "metal",     _FUR_GATE),   # [ANCHOR] P54479
    "mtZur":  ("mtZur",  "Zn2+",        "metal",     _FUR_GATE),
    "scZur":  ("scZur",  "Zn2+",        "metal",     _FUR_GATE),
    "mtFurB": ("mtFurB", "Zn2+",        "metal",     _FUR_GATE),   # FurB (Rv2359) IS a Zur (Maciag 2007)
    "scNur":  ("scNur",  "Ni2+",        "metal",     _FUR_GATE),   # Nur clade; Ni(II) is the family's novelty
    "bsPerR": ("bsPerR", "redox[H2O2]", "redox",     ""),          # [ANCHOR] P71086
    "liPerR": ("liPerR", "redox[H2O2]", "redox",     ""),
    "bbBosR": ("bbBosR", "redox[H2O2]", "redox",     ""),
    "bjIrr":  ("bjIrr",  "heme",        "non-metal", ""),
}

# The nine families added from the 2026 SSN drop. Keyed by the clade FASTA stem (the representative), like
# _FUR_TABLE. Labels are grounded in the SSN gene-name evidence + UniProt REST (COFACTOR/BINDING) + primary
# literature; 17/19 metal representatives are ANCHOR-VERIFIED (the reviewed rep accession is a member of its
# clade FASTA -- the Fur convention). Full provenance + refs live in analysis/ssn_final_state/SSN_MASTER.*.
# gate is a descriptive coordination-signature LABEL for is_metal clusters only (metal_site.coordination_gate
# stays sequence-driven and returns "" for these families -> their metal calls are SSN/anchor-based, reported
# "coordination gate silent"); redox/thiol chemistry is documented in SSN_MASTER.coord_residues.
#   tuple = (representative, inducer, role, gate)  role in {metal, non-metal, redox, multidrug, unknown}
_NEW_TABLES: dict[str, dict[str, tuple]] = {
    COPY: {   # PF03965 CopY/BlaI/MecI winged-HTH: one Cu(I) clade + a beta-lactam (BlaI/MecI) branch
        "c1": ("CopY",           "Cu+",          "metal",     "CopY-Cu(Cys)"),   # [ANCHOR] Q47839 (Cu; struct-Zn)
        "c2": ("BlaI/CopY-like", "beta-lactam",  "non-metal", ""),
        "c3": ("BlaI-like",      "beta-lactam",  "non-metal", ""),
        "c4": ("BlaI/MecI",      "beta-lactam",  "non-metal", ""),               # [ANCHOR] P0A042 / P68261
    },
    CSOR: {   # PF02583 CsoR-RcnR-CstR-FrmR-InrS superfamily; Cu(I)/Ni-Co/persulfide/formaldehyde
        "bsCsoR":       ("CsoR/RicR", "Cu+",           "metal",     "CsoR-CHC"),
        "mtCsoR":       ("CsoR",      "Cu+",           "metal",     "CsoR-CHC"),      # [ANCHOR] P9WP49
        "ttCsoR":       ("CsoR",      "Cu+",           "metal",     "CsoR-CHC"),
        "spNreA":       ("CsoR-like", "Cu+",           "metal",     "CsoR-CHC"),      # tentative (see MASTER)
        "ecRcnR":       ("RcnR",      "Ni2+/Co2+",     "metal",     "CsoR-His/Cys"),  # [ANCHOR] P64530
        "syInrS":       ("InrS",      "Ni2+",          "metal",     "CsoR-His/Cys"),
        "saCstR":       ("CstR",      "persulfide/RSS","redox",     ""),
        "scCstR":       ("CsoR/CstR", "Cu+/persulfide","redox",     ""),
        "ecFrmR":       ("FrmR/RcnR", "formaldehyde",  "non-metal", ""),             # [ANCHOR] P0AAP3
        "seFrmR_group": ("FrmR",      "formaldehyde",  "non-metal", ""),
    },
    DTXR: {   # PF01325 DtxR/MntR: Fe(II) (IdeR/DtxR) + Mn(II) (MntR/SloR)
        "IdeR":   ("IdeR/DtxR", "Fe2+", "metal", "DtxR-carboxylate"),   # [ANCHOR] P9WMH1
        "bsMntR": ("MntR",      "Mn2+", "metal", "DtxR-carboxylate"),   # [ANCHOR] P54512
        "ecMntR": ("MntR",      "Mn2+", "metal", "DtxR-carboxylate"),   # [ANCHOR] P0A9F1
        "mtMntR": ("MntR/SirR", "Mn2+", "metal", "DtxR-carboxylate"),   # Rv2788 (SirR renamed MntR)
        "smSloR": ("SloR",      "Mn2+", "metal", "DtxR-carboxylate"),
    },
    GNTR: {   # PF00392 GntR: organic effectors (not a metalloregulator family; here for completeness)
        "GntR_c1": ("GntR",      "gluconate",  "non-metal", ""),   # [ANCHOR] P10585
        "GntR_c2": ("CitO/GntR", "citrate",    "non-metal", ""),
        "GntR_c3": ("PdhR/LldR", "2-oxoacid",  "non-metal", ""),
        "GntR_c4": ("UxuR/ExuR", "hexuronate", "non-metal", ""),   # [ANCHOR] P39161
        "GntR_c5": ("LutR-like", "organic",    "non-metal", ""),
    },
    LYSR: {   # PF03466 LysR-type: OxyR (redox) + ModE (molybdate) + CysB/CmpR (sulfur/carbon)
        "LysR_c1": ("CysB/CmpR", "O-acetylserine", "non-metal", ""),
        "LysR_c2": ("OxyR",      "redox[H2O2]",    "redox",     ""),          # [ANCHOR] P0ACQ4
        "ModE":    ("ModE",      "molybdate",      "metal",     "ModE-Mop"),  # [ANCHOR] P0A9G8
    },
    MARR: {   # PF01047 MarR/SlyA: AdcR (Zn) + OhrR (redox) + MarR (multidrug)
        "c1":     ("AdcR",      "Zn2+",                 "metal",     "MarR-Zn(His/Glu)"),  # [ANCHOR] Q04I02
        "ecMarR": ("MarR",      "salicylate/multidrug", "multidrug", ""),                  # [ANCHOR] P27245
        "OhrR":   ("OhrR/SarZ", "redox[ROOH]",          "redox",     ""),                  # [ANCHOR] O34777
    },
    NIKR: {   # PF08753 NikR: Ni(II)
        "NikR_c1": ("NikR", "Ni2+", "metal", "NikR-His/Cys"),   # [ANCHOR] P0A6Z6
    },
    RRF2: {   # PF02082 Rrf2: IscR ([2Fe-2S]) / NsrR (NO) / CymR (cysteine) / RsrR (NAD redox)
        "ecIscR": ("IscR",      "redox[2Fe-2S]", "redox",     ""),   # [ANCHOR] P0AGK8
        "saCymR": ("CymR",      "cysteine/OAS",  "non-metal", ""),   # [ANCHOR] O34527
        "scNsrR": ("NsrR",      "NO",            "redox",     ""),
        "svRsrR": ("RsrR",      "redox[NAD(P)]", "redox",     ""),
        "tpIscR": ("IscR-like", "redox[2Fe-2S]", "redox",     ""),
    },
    TETR: {   # PF00440 TetR/AcrR: organic/multidrug (KstR2 steroid catabolism here)
        "c1_tetR": ("KstR2/TetR", "steroid-CoA/organic", "non-metal", ""),
    },
}

# recognition-helix signature. MerR N-terminal HTH `IRYYExxG` (also `TLRYE`); the data shows IRYYE /
# VRFYE variants -> anchor on R[YF][YF]E with a tolerant first residue. ArsR uses a different readout
# mode and is not motif-located here (returns None) -- consistent with plan SS7.
# Fur: the preprint reports `TVYRTLQ` as the conserved DBD motif. Measured column-wise on the preprint's
# OWN per-cluster MAFFT alignments (the same input its WebLogo figures are drawn from), that holds up: the
# `TVYR` core is essentially INVARIANT (per-column information content 4.2-4.3 bits of a 4.32 max) in all
# 10 canonical Fur/Zur/Mur/Nur clusters, and is carried by 96-100% of their members. The tail is what
# varies, so the published 7-mer is a representative consensus rather than a literal one:
#     mtFurB/mtZur/scZur  TVYRTLQ  (exact match -- Actinobacteria)
#     bsFur/cjFur/scNur   TVYRTL.  |  baMur TVYRTV  |  ecFur TVYR*V*LTQ  |  bsZur T*I*YR*N*L
# So anchor on the invariant core, not the full string (a `TVYRTLQ` literal would miss ecFur entirely).
# The core is absent from PerR/Irr/BosR -- and informatively so: bsPerR keeps `ATVY` (IC 4.3) but reads
# `ATVY*NN*` , i.e. it has lost the DNA-contacting arginine (the preprint's mgFur R57), consistent with its
# finding that PerR/Irr are the family's most divergent clade and do not read a canonical Fur box.
_RECOG = {
    MERR: re.compile(r"[ILVMTA]R[YF][YF]E"),
    FUR: re.compile(r"[ST][VIA]YR"),
}


# --------------------------------------------------------------------------- header parsing
def parse_accession(header: str) -> str:
    """Extract the UniProt accession from either family's FASTA header."""
    h = header[1:] if header.startswith(">") else header
    first = h.strip().split(None, 1)[0]
    if "|" in first:
        parts = first.split("|")
        if len(parts) >= 2 and parts[0] in ("tr", "sp"):
            return parts[1]
        return parts[0]
    return first


def parse_header(header: str) -> tuple[str, dict]:
    """Return (accession, meta). meta carries whatever the header exposes (organism, pfam, GN/OX)."""
    h = (header[1:] if header.startswith(">") else header).strip()
    acc = parse_accession(header)
    meta: dict = {}
    # ArsR-style UniProt description with OS=/OX=/GN= key-value tail
    if "OS=" in h or "OX=" in h or "GN=" in h:
        for key in ("OS", "OX", "GN", "PE", "SV"):
            m = re.search(rf"\b{key}=(.+?)(?=\s+\w\w?=|$)", h)
            if m:
                meta[key.lower()] = m.group(1).strip()
    # MerR-style `<num>|<organism>|<pfam,pfam>` tail
    elif "|" in h:
        tail = h.split(None, 1)[1] if " " in h else ""
        bits = tail.split("|")
        if len(bits) >= 3:
            meta["organism"] = bits[1].rstrip(".")
            meta["pfam"] = bits[2].split(",")
    return acc, meta


# --------------------------------------------------------------------------- enriched-label overrides
# The hardcoded tables above are the SEED. `tools/ssn_ligand_mine.py` enriches/corrects the cluster ->
# inducer labels (UniProt/PDB curated-binder mining + membership analysis + Folddisco) and writes them to
# packaged data/ssn/database/cluster_inducer_table.json. Loading that file as an OVERRIDE here gives the SSN its
# predictive power: assign_cluster -> info_for_cluster -> inducer inference now see the mined inducers, not
# just the 6 hand anchors. Absent file => pure seed defaults (so the module still works standalone).
_LABEL_TABLE_CACHE = _DATABASE / "cluster_inducer_table.json"


def _load_label_overrides() -> dict:
    if not _LABEL_TABLE_CACHE.exists():
        return {}
    try:
        rows = json.loads(_LABEL_TABLE_CACHE.read_text(encoding="utf-8"))
    except Exception:
        return {}
    return {r["cluster_id"]: r for r in rows if isinstance(r, dict) and "cluster_id" in r}


def _apply_override(c: ClusterInfo, ov: dict | None) -> ClusterInfo:
    if not ov:
        return c
    return ClusterInfo(
        cluster_id=c.cluster_id, family=c.family,
        representative=ov.get("representative") or c.representative,
        inducer=ov["inducer"] if ov.get("inducer") is not None else c.inducer,
        role=ov.get("role") or c.role,
        gate=ov["gate"] if ov.get("gate") is not None else c.gate,
        fasta=c.fasta, n_seqs=c.n_seqs)


# --------------------------------------------------------------------------- cluster loading
def _merr_clusters() -> list[ClusterInfo]:
    out = []
    for fa in sorted((SSN_DIR / "MerRs").glob("MerR_cluster_*.fasta")):
        m = re.search(r"cluster_(\d+)", fa.name)
        if not m:
            continue
        n = int(m.group(1))
        rep, inducer, role, gate = _MERR_TABLE.get(n, ("", None, "unknown", ""))
        out.append(ClusterInfo(f"MerR_{n}", MERR, rep, inducer, role, gate, fa))
    return out


def _arsr_clusters() -> list[ClusterInfo]:
    out = []
    for fa in sorted((SSN_DIR / "ArsRs").glob("c*.fasta"),
                     key=lambda p: int(re.search(r"c(\d+)", p.stem).group(1))):
        cid = fa.stem  # c1..c10
        rep, inducer, role = _ARSR_TABLE.get(cid, ("", None, "unknown"))
        out.append(ClusterInfo(f"ArsR_{cid}", ARSR, rep, inducer, role, _ARSR_GATE, fa))
    return out


def _fur_clusters() -> list[ClusterInfo]:
    out = []
    for fa in sorted((SSN_DIR / "Furs").glob("*.fasta")):
        stem = fa.stem                                     # ecFur, scNur, bjIrr, ...
        rep, inducer, role, gate = _FUR_TABLE.get(stem, (stem, None, "unknown", ""))
        out.append(ClusterInfo(f"Fur_{stem}", FUR, rep, inducer, role, gate, fa))
    return out


def _new_family_clusters(family: str) -> list[ClusterInfo]:
    """Stem-keyed loader for the nine 2026-drop families (same shape as `_fur_clusters`, parameterised).
    Clade FASTAs live under packaged `data/ssn/clusters/<tag>s/`; cluster_id is `<tag>_<stem>` with redundant
    leading `<tag>_` in the stem collapsed (so `GntR_c1.fasta` -> `GntR_c1`, not `GntR_GntR_c1`)."""
    tag = family_tag(family)
    table = _NEW_TABLES.get(family, {})
    out = []
    for fa in sorted((SSN_DIR / f"{tag}s").glob("*.fasta")):
        stem = fa.stem
        rep, inducer, role, gate = table.get(stem, (stem, None, "unknown", ""))
        cid_stem = stem[len(tag) + 1:] if stem.startswith(tag + "_") else stem
        out.append(ClusterInfo(f"{tag}_{cid_stem}", family, rep, inducer, role, gate, fa))
    return out


def iter_sequences(cluster, *, limit: int | None = None):
    """Yield (accession, sequence) for a ClusterInfo (or a raw FASTA Path). `limit` caps the count."""
    fa = cluster.fasta if isinstance(cluster, ClusterInfo) else Path(cluster)
    acc, buf, n = None, [], 0
    with fa.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if acc is not None:
                    yield acc, "".join(buf)
                    n += 1
                    if limit and n >= limit:
                        return
                acc, buf = parse_accession(line), []
            else:
                buf.append(line.strip())
    if acc is not None and (not limit or n < limit):
        yield acc, "".join(buf)


def _count_seqs(fa: Path) -> int:
    n = 0
    with fa.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def load_clusters(family: str | None = None, *, with_counts: bool = False,
                  apply_overrides: bool = True) -> list[ClusterInfo]:
    """All ClusterInfo for a family (or both). `with_counts` fills n_seqs (one pass per file).
    `apply_overrides` (default on) layers the enriched labels from cluster_inducer_table.json over the
    seed tables; pass False to get the raw hardcoded seed (used by the seed exporter / self-test)."""
    clusters: list[ClusterInfo] = []
    if family in (None, MERR):
        clusters += _merr_clusters()
    if family in (None, ARSR):
        clusters += _arsr_clusters()
    if family in (None, FUR):
        clusters += _fur_clusters()
    for _fam in (COPY, CSOR, DTXR, GNTR, LYSR, MARR, NIKR, RRF2, TETR):
        if family in (None, _fam):
            clusters += _new_family_clusters(_fam)
    if with_counts:
        clusters = [ClusterInfo(**{**c.__dict__, "n_seqs": _count_seqs(c.fasta)}) for c in clusters]
    if apply_overrides:
        ov = _load_label_overrides()
        clusters = [_apply_override(c, ov.get(c.cluster_id)) for c in clusters]
    return clusters


def cluster_table(family: str | None = None, *, apply_overrides: bool = False) -> list[dict]:
    """The cluster -> inducer artifact (with the gate flag), as plain dicts for JSON/CSV export.
    Defaults to the raw SEED (apply_overrides=False); the enriched/consumed table is maintained by
    tools/ssn_ligand_mine.py in the packaged cluster-inducer table."""
    rows = []
    for c in load_clusters(family, with_counts=True, apply_overrides=apply_overrides):
        rows.append(dict(cluster_id=c.cluster_id, family=c.family, representative=c.representative,
                         inducer=c.inducer, role=c.role, is_metal=c.is_metal, gate=c.gate,
                         n_seqs=c.n_seqs, fasta=str(c.fasta.relative_to(_REPO))))
    return rows


# --------------------------------------------------------------------------- accession -> cluster map
def build_accession_map(*, force: bool = False) -> dict[str, str]:
    """Map every member accession -> cluster_id. Read from packaged data, or rebuilt in the user cache.

    THE CLADES ARE NOT A PARTITION. 2,534 accessions belong to more than one, and three pairs are
    strict subsets: CsoR_ecRcnR (480) inside CsoR_ecFrmR (5,069), CsoR_seFrmR_group (1,973) inside
    the same, and Fur_mtZur (82) inside Fur_mtFurB (84). The first of those is the one that bites:
    RcnR and FrmR sense different things, so resolving a shared accession to the wrong side of that
    nesting is a wrong inducer, not a bookkeeping detail.

    The rule is MOST SPECIFIC CLADE WINS -- the smallest membership, ties broken by cluster_id so the
    result cannot depend on filesystem or iteration order. A protein belonging to both a broad clade
    and a narrow one inside it is better described by the narrow one; that is what drawing the narrow
    clade at all was meant to express.

    This replaces last-write-wins, which resolved overlaps by whichever clade load_clusters() read
    last. Measured over all 285,672 member accessions the two rules agree everywhere -- the old one
    reached the right answer, but by accident of ordering, and reordering load_clusters() would have
    changed predictions silently.
    """
    if not force:
        for cached in (_ACC_MAP_CACHE, _ACC_MAP_PACKAGED):
            if cached.exists():
                return json.loads(cached.read_text(encoding="utf-8"))

    members: dict[str, set[str]] = {}
    sizes: dict[str, int] = {}
    for c in load_clusters():
        n = 0
        with c.fasta.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith(">"):
                    members.setdefault(parse_accession(line), set()).add(c.cluster_id)
                    n += 1
        sizes[c.cluster_id] = n

    amap = {acc: min(cids, key=lambda cid: (sizes.get(cid, 0), cid))
            for acc, cids in members.items()}

    _CACHE.mkdir(parents=True, exist_ok=True)
    _ACC_MAP_CACHE.write_text(json.dumps(amap), encoding="utf-8")
    # Record WHICH accessions were ambiguous and how each was resolved, not just how many. A bare
    # count cannot be audited, and this is the file to read when a clade label looks wrong.
    multi = {a: sorted(c) for a, c in members.items() if len(c) > 1}
    if multi:
        (_CACHE / "accession_map_overlaps.json").write_text(json.dumps(
            {"rule": "most specific clade wins (smallest membership, ties by cluster_id)",
             "n_multi_clade": len(multi),
             "cluster_sizes": sizes,
             "resolved": {a: {"clades": c, "chosen": amap[a]} for a, c in sorted(multi.items())}},
            indent=1))
    return amap


def cluster_for_accession(acc: str, amap: dict[str, str] | None = None) -> str | None:
    amap = amap if amap is not None else build_accession_map()
    return amap.get(acc)


def info_for_cluster(cluster_id: str) -> ClusterInfo | None:
    for c in load_clusters():
        if c.cluster_id == cluster_id:
            return c
    return None


def assign_cluster_by_hits(hit_accessions, *, amap: dict[str, str] | None = None,
                           min_frac: float = 0.5):
    """Assign a query to a cluster from homolog hits the caller already computed (MMseqs/blastp).

    Pure majority vote over the clusters of the hit accessions. Returns (cluster_id, support_frac)
    or (None, 0.0) if no hit maps or the plurality is below `min_frac`. The MMseqs dependency stays
    at the call site (Phase 1 homolog_selection) -- this keeps ssn_clusters tool-free and testable.
    """
    amap = amap if amap is not None else build_accession_map()
    votes = Counter()
    mapped = 0
    for acc in hit_accessions:
        cid = amap.get(acc)
        if cid:
            votes[cid] += 1
            mapped += 1
    if not mapped:
        return None, 0.0
    cid, n = votes.most_common(1)[0]
    frac = n / mapped
    return (cid, frac) if frac >= min_frac else (None, frac)


# --------------------------------------------------------------------------- recognition helix
def recognition_region(seq: str, family: str, *, flank: int = 6) -> tuple[int, int, str] | None:
    """Locate the family recognition-helix motif. Returns (start, end, window) or None.

    Window is the motif +/- `flank` residues -- the slice the Test-1 pre-check compares within/across
    clusters and that later seeds the Folddisco per-monomer motif query.
    """
    pat = _RECOG.get(family)
    if pat is None:
        return None
    m = pat.search(seq)
    if not m:
        return None
    s = max(0, m.start() - flank)
    e = min(len(seq), m.end() + flank)
    return s, e, seq[s:e]


# --------------------------------------------------------------------------- self-test
def _demo() -> None:
    merr = load_clusters(MERR)
    arsr = load_clusters(ARSR)
    fur = load_clusters(FUR)
    print(f"MerR clusters: {len(merr)}  ArsR clusters: {len(arsr)}  Fur clusters: {len(fur)}")
    assert len(merr) == 10, f"expected 10 MerR clusters (7 recovered = CueR), got {len(merr)}"
    assert len(arsr) == 10, f"expected 10 ArsR clusters, got {len(arsr)}"
    assert len(fur) == 14, f"expected 14 Fur clusters, got {len(fur)}"

    # 2026 drop: the nine new families load with their expected clade counts (canonical names match
    # packaged pfam_to_family.json; the optional structure axis is not shipped here).
    new_counts = {COPY: 4, CSOR: 10, DTXR: 5, GNTR: 5, LYSR: 3, MARR: 3, NIKR: 1, RRF2: 5, TETR: 1}
    for fam, k in new_counts.items():
        cs = load_clusters(fam)
        assert len(cs) == k, f"{fam}: expected {k} clusters, got {len(cs)}"
    total = len(merr) + len(arsr) + len(fur) + sum(new_counts.values())
    assert len(load_clusters()) == total, f"load_clusters(None) must return all 12 families ({total})"
    assert total == 71, f"expected 71 clusters across 12 families, got {total}"

    # anchored metal/effector labels for the new families (seed layer; provenance + refs in
    # analysis/ssn_final_state/SSN_MASTER.*). 17/19 metal reps are anchor-verified in their clade FASTA.
    seed = {c.cluster_id: c for c in load_clusters(apply_overrides=False)}
    assert seed["NikR_c1"].inducer == "Ni2+" and seed["NikR_c1"].is_metal, seed["NikR_c1"]
    assert seed["CsoR_mtCsoR"].inducer == "Cu+" and seed["CsoR_ecRcnR"].inducer == "Ni2+/Co2+", "CsoR Cu/Ni"
    assert seed["DtxR_IdeR"].inducer == "Fe2+" and seed["DtxR_bsMntR"].inducer == "Mn2+", "DtxR Fe/Mn"
    assert seed["MarR_c1"].inducer == "Zn2+" and seed["MarR_c1"].is_metal, "AdcR -> Zn"
    assert seed["LysR_ModE"].inducer == "molybdate" and seed["Rrf2_ecIscR"].role == "redox", "ModE/IscR"
    assert seed["CopY_c1"].inducer == "Cu+" and seed["CopY_c4"].role == "non-metal", "CopY Cu + BlaI branch"
    # new metal clusters carry a descriptive gate; non-metal/redox carry none (the Fur/PerR convention)
    for c in load_clusters(apply_overrides=False):
        if c.family in new_counts:
            assert (c.gate != "") == c.is_metal, f"{c.cluster_id}: gate/is_metal mismatch ({c.gate!r})"

    # family tags must be distinct -- a collision silently clobbers another family's member DB
    tags = {family_tag(f) for f in FAMILIES}
    assert len(tags) == len(FAMILIES), f"family_tag collision: {tags}"
    assert family_tag(ARSR) == "ArsR", "ArsR/SmtB must keep its legacy on-disk tag"

    # cluster 7 (CueR/Cu+) was recovered from the size-ranked partition -> must load and be metal-gated
    cue = info_for_cluster("MerR_7")
    assert cue and cue.inducer == "Cu+" and cue.is_metal and cue.gate == "CX7-8C", cue

    # Fur: the 4 subgroups must all be represented, and only metal clusters carry the gate
    fur_ind = {c.cluster_id: (c.inducer, c.role, c.gate) for c in fur}
    assert fur_ind["Fur_ecFur"][0] == "Fe2+" and fur_ind["Fur_ecZur"][0] == "Zn2+", fur_ind
    assert fur_ind["Fur_scNur"][0] == "Ni2+", "Nur clade (Ni) must load -- recovered from scNur.fasta.aln"
    assert fur_ind["Fur_bjIrr"] == ("heme", "non-metal", ""), "heme regulator must not carry a metal gate"
    assert fur_ind["Fur_bsPerR"][1] == "redox" and fur_ind["Fur_bsPerR"][2] == "", fur_ind
    assert all(c.gate == _FUR_GATE for c in fur if c.is_metal), "every Fur metal cluster needs the gate"

    # curated labels wired correctly
    z = info_for_cluster("MerR_9")
    assert z and z.inducer == "Zn2+" and z.is_metal and z.gate == "CX7-8C", z
    nm = info_for_cluster("MerR_3")
    assert nm and nm.role == "non-metal" and nm.gate == "", "non-metal cluster must not carry a gate"
    a1_seed = next(c for c in load_clusters(ARSR, apply_overrides=False) if c.cluster_id == "ArsR_c1")
    assert a1_seed.inducer is None and a1_seed.role == "unknown", "seed ArsR c1 must be unlabeled"
    a1 = info_for_cluster("ArsR_c1")     # enriched (override layer): As3+ once ssn_ligand_mine has run
    print(f"  ArsR_c1 enriched inducer={a1.inducer} role={a1.role} (seed inducer={a1_seed.inducer})")

    # header parsing for all three styles
    assert parse_accession(">A0A014LWC1 14|Erwinia mallotivora.|PF00376,PF09278") == "A0A014LWC1"
    assert parse_accession(">A0A023CL38 140 bp") == "A0A023CL38", "DtxR minimal '>acc N bp' header"
    acc, meta = parse_header(">tr|A0A011QX74|A0A011QX74_9LACT X OS=Alkalibacterium sp. AK22 "
                             "OX=1229520 GN=ADIAL_0133 PE=4 SV=1")
    assert acc == "A0A011QX74" and meta.get("gn") == "ADIAL_0133" and meta.get("ox") == "1229520", meta

    # recognition motif present in a real ecZntR-cluster sequence; absent-locator returns None for ArsR
    znt = ("MYKIGQLAKLADVTPDTIRYYEKQQMMDHEIRSEGGFRLYSQNDLQRLKFIRYGRQLGFSLEAIRELLSIRVDPEHHTCQ")
    rr = recognition_region(znt, MERR)
    assert rr and "RYYE" in rr[2], rr
    assert recognition_region(znt, ARSR) is None, "ArsR recognition motif not defined yet"

    # Fur recognition helix: real E. coli Fur (P0A9A9) carries TVYR*V*LNQ -- the published TVYRTLQ
    # string does NOT match it, which is exactly why the locator anchors on the measured [ST][VIA]YR core.
    ecfur = ("MTDNNTALKKAGLKVTLPRLKILEVLQEPDNHHVSAEDLYKRLIDMGEEIGLATVYRVLNQFDDAGIVTRHNFEGGKSVFEL"
             "TQQHHHDHLICLDCGKVIEFSDDSIEARQREIAAKHGIRLTNHSLYLYGHCAEGDCREDEHAHEGK")
    fr = recognition_region(ecfur, FUR)
    assert fr and "TVYR" in fr[2], f"Fur recognition helix not located in E. coli Fur: {fr}"
    assert not re.search(r"TVYRTLQ", ecfur), "sanity: the published TVYRTLQ really is absent from ecFur"

    # accession map + majority vote (build is cached afterwards). Accessions below are the first
    # member of each cluster's FASTA: A0A009FWB9 -> ppCadR (MerR_1), A0A014LWC1 -> ecZntR (MerR_9).
    amap = build_accession_map(force=True)
    print(f"accession map: {len(amap)} sequences mapped")
    assert cluster_for_accession("A0A009FWB9", amap) == "MerR_1", "first ppCadR member -> MerR_1"
    assert cluster_for_accession("A0A014LWC1", amap) == "MerR_9", "first ecZntR member -> MerR_9"
    cid, frac = assign_cluster_by_hits(["A0A009FWB9", "A0A011QX74", "A0A009FWB9"], amap=amap)
    assert cid == "MerR_1" and frac >= 0.66, f"majority vote -> MerR_1, got {cid} ({frac:.2f})"
    print(f"majority vote demo -> {cid} (support {frac:.2f})")

    print("OK: SSN cluster parsing, curated table, header parsing, recognition locator, acc-map verified.")


def _write_table() -> Path:
    """Export the raw SEED cluster -> inducer table to a sidecar. The CONSUMED table
    (cluster_inducer_table.json) is owned by tools/ssn_ligand_mine.py (seed + mining + corrections); this
    exporter writes the seed separately so it can never clobber the enriched labels."""
    rows = cluster_table(apply_overrides=False)
    _CACHE.mkdir(parents=True, exist_ok=True)
    out = _CACHE / "cluster_inducer_table.seed.json"
    out.write_text(json.dumps(rows, indent=2))
    print(f"wrote {out.relative_to(_REPO)} (seed; enriched table is maintained by ssn_ligand_mine)")
    for r in rows:
        print(f"  {r['cluster_id']:9s} {str(r['representative']):8s} {str(r['inducer']):14s} "
              f"{r['role']:10s} gate={r['gate'] or '-':8s} n={r['n_seqs']}")
    return out


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "table":
        _write_table()
    else:
        _demo()
