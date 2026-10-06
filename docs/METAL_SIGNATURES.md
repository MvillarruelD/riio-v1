# Metal-coordination signatures

**Generated file — do not edit.** Regenerate with `python tools/build_metal_reference.py`; `--check` fails when it is stale. The signatures themselves live in `predictor/structure/metal_motifs.py`, and every prevalence below is measured at build time against the 288,212 sequences of the shipped SSN member databases (`predictor/data/ssn/database/ssn_members_<TAG>.fasta`).

This file exists because a prevalence stated by hand drifts. `structure/metal_site.py` records that the NikR H-x-H-x(4,7)-C site occurs in "≤3.2 % of every other family"; measured against the database this package ships, it occurs in **60.6 % of Fur**. Nothing re-checked the claim after it was written. Everything here is recomputed from the data instead.

## How to read it

- **in-family** — the share of that family's members carrying the motif. On its own this is often the wrong denominator: a motif specific to one clade looks weak against the whole family.
- **top clade** — the clade where it concentrates. This is usually the real signal.
- **max other family** — the highest rate in any family it does *not* belong to. This decides whether the motif can be trusted without a family gate.
- **scope** — `global` means the motif is specific enough to use anywhere (measured ≤ 2 % in every other family); `family` means it is only meaningful inside its own family and demonstrably fires elsewhere.
- **names ion** — whether this motif alone may name an ion in the inducer call. A motif marking a chemistry shared by several ions (a bare Cys pair) documents the site without naming the metal.

## Detected signatures

| motif | family | pattern | implies | class | in-family | top clade | max other family | scope | names ion |
|---|---|---|---|---|---:|---|---:|---|:-:|
| `ArsR_a5_DxHx10Hx2HE` | ArsR/SmtB | `D.H.{10}H.{2}[HE]` | Zn2+, Ni2+ | metal | 17.7 % | `ArsR_c5` 99.0 % | 0.3 % (NikR) | global | yes |
| `ArsR_a3_CVC` | ArsR/SmtB | `C[VA]C` | — | metal | 31.1 % | `ArsR_c1` 90.0 % | 2.8 % (MerR) | family | no |
| `ArsR_a3_thiol_CxC` | ArsR/SmtB | `C.C` | — | metal | 35.3 % | `ArsR_c1` 94.5 % | 31.9 % (CopY) | family | no |
| `ArsR_AfArsR_CCx5C` | ArsR/SmtB | `CC.{5}C` | As3+ | metal | 1.1 % | `ArsR_c4` 5.8 % | 0.8 % (Rrf2) | family | no |
| `MerR_ZntR_CCx3Hx4C` | MerR | `CC.{3}H.{4}C` | Zn2+ | metal | 3.2 % | `MerR_9` 60.3 % | 0.0 % (ArsR/SmtB) | global | yes |
| `MerR_SoxR_FeS` | MerR | `C.{2}C.C.{4,6}C` | — | redox | 19.2 % | `MerR_4` 99.4 % | 0.0 % (Fur) | global | no |
| `MerR_metal_CX7_8C` | MerR | `C.{7,8}C` | — | metal | 80.7 % | `MerR_4` 99.4 % | 35.2 % (LysR-type (LTTR)) | family | no |
| `Fur_HHxHx2Cx2C` | Fur | `HH.H.{2}C.{2}C` | — | metal | 37.3 % | `Fur_bsZur` 99.3 % | 0.0 % (MerR) | global | no |
| `Fur_site2_His` | Fur | `HH[HD]H` | — | metal | 52.3 % | `Fur_bsZur` 99.3 % | 12.5 % (NikR) | family | no |
| `Fur_site1_Cx2C` | Fur | `C.{2}C` | Zn2+ | metal | 74.0 % | `Fur_bbBosR` 100.0 % | 20.1 % (MerR) | family | no |
| `NikR_His3Cys` | NikR | `H.H.{4,7}C` | Ni2+ | metal | 85.4 % | `NikR_c1` 85.4 % | 60.6 % (Fur) | family | yes |
| `CsoR_WXYZ_metal_site` | CsoR/FrmR | `C.{15,40}H.{3}C` | — | metal | 45.1 % | `CsoR_mtCsoR` 99.3 % | 3.6 % (Fur) | family | no |
| `CopY_CxC` | CopY | `C.C` | Cu+ | metal | 31.9 % | `CopY_c1` 99.6 % | 35.3 % (ArsR/SmtB) | family | yes |
| `Rrf2_FeS_Cys_triad` | Rrf2 | `C.{2,20}C.{2,20}C` | — | redox | 86.2 % | `Rrf2_tpIscR` 100.0 % | 41.9 % (MerR) | family | no |

## Documented but not detected

These are real chemistry from the literature that no sequence pattern captures — the ligands come from different structural elements, or the discriminator is a *negative* design element, or the difference is a fold rather than a composition. They are listed so this table is a complete account of what each family senses, not only the automatable part.

| motif | family | implies | ligand set | why not detected |
|---|---|---|---|---|
| `ArsR_CmtR_a4_plus_tail` | ArsR/SmtB | Cd2+ | 2 Cys from the alpha4 DNA-binding helix + 1 Cys from the C-terminal tail | The ligands come from three different structural elements, so no linear spacing identifies them; separating CmtR from any other multi-Cys ArsR needs an alignment. |
| `MerR_CueR_coinage` | MerR | Cu+, Ag+, Au+ | 2 Cys from a SHORT MBL | Selectivity for the monovalent coinage metals comes from two NEGATIVE design elements -- a conserved Ser in alpha5 that excludes divalents, and an MBL shorter than the divalent sensors' -- neither of which is a presence/absence motif. Distinguishing CueR from a divalent sensor needs the alignment, not a regex. |
| `MerR_PbrR_CadR_trigonal` | MerR | Pb2+, Cd2+ | 3 Cys: 2 from the MBL + 1 from the dimerisation helix | The third Cys comes from a different helix, so the ligand set has no fixed spacing. CadR adds an Asn as a fourth ligand that Pb's lone pair excludes in PbrR. |
| `CsoR_RcnR_Nterm_His` | CsoR/FrmR | Ni2+, Co2+ | an N-terminal His donor added to the W-X-Y-Z frame (RcnR His3/Cys35/His60/His64) | RcnR is nested INSIDE the FrmR clade in our SSN, and the N-terminal His that distinguishes it is a single residue at the chain start. This is the known root cause of the RcnR->FrmR misassignment; a regex does not fix it. |
| `CsoR_FrmR_ligand_loss` | CsoR/FrmR | — | Pro2 + Cys35 ONLY -- His64 replaced by Glu, N-terminal His absent | coli FrmR (1 His, 2 Cys) fails. Formaldehyde is sensed by a Pro2/Cys35 methylene bridge, not by a metal. |
| `MarR_AdcR_ZitR_linker` | MarR/SlyA | Zn2+ | ligands drawn from the unstructured region between alpha1 and alpha2 | The discriminator is STRUCTURAL, not compositional: in the zinc-sensing MarRs (AdcR, ZitR) the ligands come from a disordered alpha1-alpha2 linker that is a continuous longer alpha2 helix in every other MarR. Detecting it needs the fold or an alignment. These are also the only characterised metal sensors in this large superfamily, and uniquely for MarR the metal ACTIVATES operator binding. |
| `DtxR_binuclear` | DtxR/MntR | Fe2+, Mn2+ | — | Non-cognate metals (Fe, Zn, Co in MntR) bind MONONUCLEARLY and simply fail to trigger the allosteric response. |
| `TetR_SczA_dual_site` | TetR/AcrR | Zn2+ | — | SczA is the lone metal-responsive TetR and activates a Zn efflux transporter; the rest of this very large family senses organic ligands in independently evolved C-terminal cavities. |
| `GntR_buried_structural_Zn` | GntR | Zn2+ | — | ~70 % of GntR proteins carry a conserved Zn site, but the review states as an OPEN QUESTION whether it is regulatory or simply positions the carboxylate of the cognate organic acid. `metal_site.STRUCTURAL_SITE_FAMILIES` already encodes this caveat; finding a site here is a real observation, and 'therefore a metal sensor' is the step that does not follow. |
| `LysR_OxyR_peroxidatic_Cys` | LysR-type (LTTR) | — | a peroxidatic Cys that forms a disulfide on H2O2 oxidation | Not a metal site at all -- OxyR senses H2O2 through Cys oxidation. |

## Per-clade detail

**`ArsR_a5_DxHx10Hx2HE`** — ArsR/SmtB

| clade | prevalence |
|---|---:|
| `ArsR_c5` | 99.0 % |
| `ArsR_c2` | 69.7 % |
| `ArsR_c7` | 0.3 % |
| `ArsR_c1` | 0.1 % |
| `ArsR_c8` | 0.0 % |
| `ArsR_c10` | 0.0 % |
| `ArsR_c3` | 0.0 % |
| `ArsR_c4` | 0.0 % |
| `ArsR_c6` | 0.0 % |
| `ArsR_c9` | 0.0 % |

**`ArsR_a3_CVC`** — ArsR/SmtB

| clade | prevalence |
|---|---:|
| `ArsR_c1` | 90.0 % |
| `ArsR_c3` | 71.2 % |
| `ArsR_c2` | 40.6 % |
| `ArsR_c4` | 1.7 % |
| `ArsR_c6` | 1.1 % |
| `ArsR_c8` | 0.2 % |
| `ArsR_c7` | 0.1 % |
| `ArsR_c5` | 0.1 % |
| `ArsR_c10` | 0.0 % |
| `ArsR_c9` | 0.0 % |

**`ArsR_a3_thiol_CxC`** — ArsR/SmtB

| clade | prevalence |
|---|---:|
| `ArsR_c1` | 94.5 % |
| `ArsR_c3` | 82.2 % |
| `ArsR_c2` | 42.4 % |
| `ArsR_c6` | 13.7 % |
| `ArsR_c4` | 10.7 % |
| `ArsR_c8` | 1.6 % |
| `ArsR_c7` | 0.4 % |
| `ArsR_c5` | 0.2 % |
| `ArsR_c9` | 0.1 % |
| `ArsR_c10` | 0.0 % |

**`ArsR_AfArsR_CCx5C`** — ArsR/SmtB

| clade | prevalence |
|---|---:|
| `ArsR_c4` | 5.8 % |
| `ArsR_c6` | 3.7 % |
| `ArsR_c1` | 0.3 % |
| `ArsR_c3` | 0.3 % |
| `ArsR_c2` | 0.1 % |
| `ArsR_c7` | 0.0 % |
| `ArsR_c10` | 0.0 % |
| `ArsR_c5` | 0.0 % |
| `ArsR_c8` | 0.0 % |
| `ArsR_c9` | 0.0 % |

**`MerR_ZntR_CCx3Hx4C`** — MerR

| clade | prevalence |
|---|---:|
| `MerR_9` | 60.3 % |
| `MerR_1` | 0.4 % |
| `MerR_2` | 0.1 % |
| `MerR_7` | 0.0 % |
| `MerR_10` | 0.0 % |
| `MerR_3` | 0.0 % |
| `MerR_4` | 0.0 % |
| `MerR_5` | 0.0 % |
| `MerR_6` | 0.0 % |
| `MerR_8` | 0.0 % |

**`MerR_SoxR_FeS`** — MerR

| clade | prevalence |
|---|---:|
| `MerR_4` | 99.4 % |
| `MerR_1` | 0.0 % |
| `MerR_10` | 0.0 % |
| `MerR_2` | 0.0 % |
| `MerR_3` | 0.0 % |
| `MerR_5` | 0.0 % |
| `MerR_6` | 0.0 % |
| `MerR_7` | 0.0 % |
| `MerR_8` | 0.0 % |
| `MerR_9` | 0.0 % |

**`MerR_metal_CX7_8C`** — MerR

| clade | prevalence |
|---|---:|
| `MerR_4` | 99.4 % |
| `MerR_7` | 99.0 % |
| `MerR_6` | 97.0 % |
| `MerR_2` | 90.9 % |
| `MerR_1` | 90.3 % |
| `MerR_9` | 83.1 % |
| `MerR_3` | 0.5 % |
| `MerR_10` | 0.0 % |
| `MerR_5` | 0.0 % |
| `MerR_8` | 0.0 % |

**`Fur_HHxHx2Cx2C`** — Fur

| clade | prevalence |
|---|---:|
| `Fur_bsZur` | 99.3 % |
| `Fur_scZur` | 97.7 % |
| `Fur_cjFur` | 96.8 % |
| `Fur_mtFurB` | 89.3 % |
| `Fur_mtZur` | 89.0 % |
| `Fur_ecFur` | 70.3 % |
| `Fur_scNur` | 45.2 % |
| `Fur_bsFur` | 27.4 % |
| `Fur_bsPerR` | 4.0 % |
| `Fur_baMur` | 3.1 % |
| `Fur_bbBosR` | 0.0 % |
| `Fur_bjIrr` | 0.0 % |
| `Fur_ecZur` | 0.0 % |
| `Fur_liPerR` | 0.0 % |

**`Fur_site2_His`** — Fur

| clade | prevalence |
|---|---:|
| `Fur_bsZur` | 99.3 % |
| `Fur_scZur` | 97.8 % |
| `Fur_cjFur` | 96.8 % |
| `Fur_ecFur` | 95.6 % |
| `Fur_baMur` | 91.6 % |
| `Fur_mtFurB` | 90.5 % |
| `Fur_mtZur` | 90.2 % |
| `Fur_liPerR` | 59.8 % |
| `Fur_bsFur` | 27.4 % |
| `Fur_ecZur` | 5.1 % |
| `Fur_bjIrr` | 1.8 % |
| `Fur_bsPerR` | 1.4 % |
| `Fur_scNur` | 0.8 % |
| `Fur_bbBosR` | 0.0 % |

**`Fur_site1_Cx2C`** — Fur

| clade | prevalence |
|---|---:|
| `Fur_bbBosR` | 100.0 % |
| `Fur_bsPerR` | 100.0 % |
| `Fur_bsZur` | 100.0 % |
| `Fur_mtFurB` | 100.0 % |
| `Fur_mtZur` | 100.0 % |
| `Fur_scNur` | 100.0 % |
| `Fur_scZur` | 100.0 % |
| `Fur_bsFur` | 99.9 % |
| `Fur_ecZur` | 99.9 % |
| `Fur_cjFur` | 99.7 % |
| `Fur_ecFur` | 72.1 % |
| `Fur_baMur` | 3.5 % |
| `Fur_bjIrr` | 0.4 % |
| `Fur_liPerR` | 0.0 % |

**`NikR_His3Cys`** — NikR

| clade | prevalence |
|---|---:|
| `NikR_c1` | 85.4 % |

**`CsoR_WXYZ_metal_site`** — CsoR/FrmR

| clade | prevalence |
|---|---:|
| `CsoR_mtCsoR` | 99.3 % |
| `CsoR_bsCsoR` | 98.2 % |
| `CsoR_syInrS` | 96.5 % |
| `CsoR_spNreA` | 95.0 % |
| `CsoR_ttCsoR` | 10.4 % |
| `CsoR_saCstR` | 1.4 % |
| `CsoR_scCstR` | 0.3 % |
| `CsoR_seFrmR_group` | 0.1 % |
| `CsoR_ecFrmR` | 0.0 % |
| `CsoR_ecRcnR` | 0.0 % |

**`CopY_CxC`** — CopY

| clade | prevalence |
|---|---:|
| `CopY_c1` | 99.6 % |
| `CopY_c2` | 0.2 % |
| `CopY_c3` | 0.0 % |
| `CopY_c4` | 0.0 % |

**`Rrf2_FeS_Cys_triad`** — Rrf2

| clade | prevalence |
|---|---:|
| `Rrf2_tpIscR` | 100.0 % |
| `Rrf2_scNsrR` | 99.7 % |
| `Rrf2_ecIscR` | 94.8 % |
| `Rrf2_saCymR` | 66.1 % |
| `Rrf2_svRsrR` | 22.4 % |

## Selectivity priors

From Dudev & Lim, *Chem. Rev.* 2014, 114, 538–556 (their Table 1 for coordination number and geometry; Fig. 10 and §6.3 for ligand preference).

**Use these to sanity-check, never to rank.** A geometry that contradicts an ion is evidence the site is not that ion's. A geometry consistent with three ions does not order them: Dudev & Lim's own conclusion is that the protein matrix and the cellular free-ion concentration decide the winner, and this pipeline models neither. NikR binds Cu(II) more tightly than Ni(II) *in vitro* and is still a nickel sensor *in vivo*, because free cytosolic Cu is ~10⁻¹⁸ M.

| ion | coordination number | geometry | preferred ligand |
|---|---|---|---|
| Cu+ | 2-4 | linear / trigonal | Cys(thiolate) |
| Ag+ | 2-3 | linear / trigonal | Cys(thiolate) |
| Au+ | 2 | linear | Cys(thiolate) |
| Hg2+ | 2-4 | linear / trigonal | Cys(thiolate) |
| Zn2+ | 4 | tetrahedral | Cys / His / Asp-Glu |
| Cd2+ | 4-6 | tetrahedral / octahedral | Cys(thiolate) |
| Pb2+ | 3-4 | trigonal pyramidal (lone pair) | Cys(thiolate) |
| Ni2+ | 4 or 6 | square planar or octahedral | His / Cys |
| Co2+ | 4-6 | tetrahedral / octahedral | His / Cys |
| Fe2+ | 6 | octahedral | His / Asp-Glu |
| Mn2+ | 6 | octahedral | Asp-Glu / His |
| As3+ | 3 | trigonal pyramidal | Cys(thiolate) |

Irving–Williams order of intrinsic affinity for any ligand set: Mg2+ < Mn2+ < Fe2+ < Co2+ < Ni2+ < Cu2+ < Zn2+. Every metalloregulator works *against* this ordering, which is why the cell holds free Zn at 10⁻¹²–10⁻¹⁵ M and free Cu at ~10⁻¹⁸ M.

