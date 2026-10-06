#!/usr/bin/env python
"""Panel E for the manuscript's four-genome figure — the PREDICTED counterpart of panel D.

Panel D (theirs)  one chip per FAMILY present in a genome that COULD sense an ion, two-tone by whether
                  the family carries a MetalNet CHDE site. Family-level, potential.
Panel E (this)    one chip per CANDIDATE REGULATOR actually predicted to sense that ion, two-tone by
                  whether our sequence coordination gate fired. Regulator-level, predicted.

Deliberately mirrors panel D so the two read as a matched pair: identical metal rows in the identical
groups (essential / toxic heavy metals / metalloids / oxyanions), identical genome column order, identical
chip geometry and two-tone blue encoding. Only the SEMANTICS of a chip differ -- family-that-could vs
regulator-that-is-predicted-to -- which is exactly the comparison the pair is meant to invite.

The inset bar reuses the manuscript's own genome colours (panel A/B) to show the inducer-class split,
i.e. how much of each genome's candidate repertoire is metal at all.

    py fig5_panelE.py
"""
import collections
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import fig_style as S                                              # noqa: E402
import matplotlib.pyplot as plt                                    # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle           # noqa: E402

# the report chain's output locations, so a figure is never drawn from an archived run
sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_paths as RP

TSV = (RP.COMPLEMENT_TSV if RP.COMPLEMENT_TSV.is_file()
       else RP.RESULTS / "per_regulator_predictions.tsv")
OUT = RP.FIGURES / "F23_panelE_inducer_per_regulator"

#: genome column order EXACTLY as the manuscript's panel D
GENOMES = [
    ("M. avium subsp. hominissuis", "M. avium"),
    ("M. tuberculosis H37Rv", "M. tuberculosis"),
    ("V. cholerae RFB16", "V. cholerae"),
    ("V. vulnificus NBRC15645", "V. vulnificus"),
]
#: the manuscript's genome palette, sampled from panel A/B of the published figure
GEN_COLOR = {"M. avium": "#5B9BD5", "M. tuberculosis": "#1F4E79",
             "V. cholerae": "#9E2A2B", "V. vulnificus": "#F4B183"}

#: metal rows + groups EXACTLY as panel D (ions we never call are kept, so the grids align row-for-row)
ION_GROUPS = [
    ("Essential metals", ["Fe", "Mn", "Zn", "Ni", "Co", "Cu"]),
    ("Toxic heavy metals", ["Cd", "Pb", "Hg", "Ag", "Au"]),
    ("Metalloids", ["As", "Sb", "Bi"]),
    ("Oxyanions", ["Mo", "W"]),
]

#: family -> the 3-letter chip code used in panel D
FAM_CODE = {"ArsR": "Ars", "MerR": "Mer", "Fur": "Fur", "CsoR": "Cso", "CopY": "Cop",
            "MarR": "Mar", "GntR": "Gnt", "TetR": "Tet", "DtxR": "Dtx", "Rrf2": "Rrf",
            "LysR": "Lys", "NikR": "Nik"}

#: panel D's two-tone blue
CHIP_DARK, CHIP_LIGHT = "#1F5C99", "#BDD7EE"

CLASS_COLOR = {"metal": "#1F5C99", "redox": "#E8A33D", "organic": "#D9D9D9"}
#: The survey's MetalNet CHDE count, over THEIR 150 candidates. Kept as the published comparator
#: and labelled as such -- it is not a measurement on this run's candidate set.
MS_CHDE_N, MS_CHDE_TOT = 41, 150
#: Row label for metal calls whose ion the evidence does not settle.
UNRESOLVED_ROW = "ion not resolved"


def element_of(inducer: str):
    """'Zn2+' -> 'Zn'; 'Cd2+/Pb2+' -> ['Cd','Pb']; non-ion -> []."""
    import re
    out = []
    for tok in (inducer or "").split("/"):
        m = re.match(r"\s*([A-Z][a-z]?)\d*[+-]", tok)
        if m:
            out.append(m.group(1))
    return out


def main():
    if not TSV.is_file():
        # One fail-closed idiom across the figure set: `require` names the missing input AND the
        # stage that writes it, so the message says what to run rather than only what is absent.
        S.require(False, f"the per-regulator complement table ({TSV.name})",
                  "fig5_complement.py, which runs immediately before this script")
    rows = list(csv.DictReader(TSV.open(encoding="utf-8"), delimiter="\t"))
    by_gen = collections.defaultdict(list)
    for r in rows:
        by_gen[r["genome"]].append(r)

    # ion x genome -> [(family_code, gate_positive)] for regulators predicted to sense that ion
    cell = collections.defaultdict(list)
    unresolved = 0
    for r in rows:
        if r["inducer_class"] != "metal":
            continue
        els = element_of(r["inducer"])
        if not els:
            # A metal call whose ION is not resolved. These had NO row at all, so the second most
            # common metal call in the run was invisible in the panel that reports the metal
            # repertoire. They get their own row, because "metal, ion undetermined" is a result.
            unresolved += 1
            cell[(UNRESOLVED_ROW, r["genome"])].append(
                (FAM_CODE.get(r["family"], r["family"][:3]), r["gate"] == "True"))
            continue
        for el in els:
            cell[(el, r["genome"])].append((FAM_CODE.get(r["family"], r["family"][:3]),
                                            r["gate"] == "True"))

    # Keep only ions this run actually called. The full table is 16 rows of which 8 were empty --
    # half the panel's height spent asserting that no regulator senses tungsten.
    seen = {k[0] for k in cell}
    ions, groups = [], []
    for gname, members in ION_GROUPS:
        keep = [m for m in members if m in seen]
        if not keep:
            continue
        groups.append((gname, len(ions), len(ions) + len(keep)))
        ions += keep
    if UNRESOLVED_ROW in seen:
        # no group band for this row: it is not a chemical group, and a rotated label here
        # collides with the Metalloids band directly above it
        groups.append(("", len(ions), len(ions) + 1))
        ions.append(UNRESOLVED_ROW)

    n_metal = sum(1 for r in rows if r["inducer_class"] == "metal")
    n_tot = len(rows)

    # ---------------------------------------------------------------- layout
    fig = plt.figure(figsize=(7.5, 5.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 2.35], wspace=0.05)
    axB = fig.add_subplot(gs[0, 0])     # inducer-class split, manuscript genome colours
    axD = fig.add_subplot(gs[0, 1])     # the panel-D mirror

    # ---- left: class split per genome, in the manuscript's genome colours -----------------------
    order = ["metal", "redox", "organic"]
    for i, (gkey, glab) in enumerate(GENOMES):
        rs = by_gen.get(gkey, [])
        c = collections.Counter(r["inducer_class"] for r in rs)
        y = len(GENOMES) - 1 - i
        left = 0.0
        for k in order:
            w = c[k]
            if w:
                axB.barh(y, w, left=left, height=0.52, color=CLASS_COLOR[k],
                         edgecolor="white", linewidth=0.6)
                if w >= 4:
                    axB.text(left + w / 2, y, str(w), ha="center", va="center", fontsize=6,
                             color="white" if k != "organic" else S.INK)
                left += w
        axB.text(left + 1.5, y, f"n={len(rs)}", va="center", fontsize=6.2, color=S.NEUTRAL_D)
        # genome identity tick in the manuscript's colour
        axB.add_patch(Rectangle((-6.4, y - 0.26), 1.5, 0.52, linewidth=0,
                                facecolor=GEN_COLOR[glab], clip_on=False))
    axB.set_ylim(-0.65, len(GENOMES) - 0.35)
    axB.set_xlim(0, 78)
    axB.set_yticks(range(len(GENOMES)))
    axB.set_yticklabels([lab for _, lab in GENOMES][::-1], fontstyle="italic", fontsize=6.8)
    axB.tick_params(axis="y", length=0, pad=12)
    axB.set_xlabel("candidate regulators", fontsize=6.5)
    axB.set_title("Predicted inducer class", fontsize=7.2, pad=6)
    for sp in ("top", "right", "left"):
        axB.spines[sp].set_visible(False)
    axB.set_xticks([0, 20, 40, 60])
    axB.tick_params(axis="x", labelsize=6)

    handles = [Rectangle((0, 0), 1, 1, facecolor=CLASS_COLOR[k], linewidth=0) for k in order]
    axB.legend(handles, ["metal ion", "reactive species", "organic ligand"],
               loc="upper center", bbox_to_anchor=(0.42, -0.10), ncol=1, fontsize=6.2,
               handlelength=1.0, handleheight=1.0, borderpad=0.0, labelspacing=0.35)

    # ---- right: the panel-D mirror --------------------------------------------------------------
    # Chips WRAP within a cell, 2 per line, as panel D does. Laying them in a single line overflowed
    # the column: 4 chips at 0.30 wide need 1.34 of a 1.0-wide column, so Cu/Cd bled into the next
    # genome and two labels collided.
    CW, CH, PER_LINE, MAX_CHIPS = 0.40, 0.26, 2, 6
    for ii, ion in enumerate(ions):
        for gi, (gkey, glab) in enumerate(GENOMES):
            chips = sorted(cell.get((ion, gkey), []), key=lambda t: (not t[1], t[0]))
            shown, extra = chips[:MAX_CHIPS], max(0, len(chips) - MAX_CHIPS)
            nlines = max(1, (len(shown) + PER_LINE - 1) // PER_LINE)
            for ci, (code, gate) in enumerate(shown):
                line, col = divmod(ci, PER_LINE)
                x = gi + 0.07 + col * (CW + 0.04)
                y = ii + (line - (nlines - 1) / 2) * (CH + 0.05)
                axD.add_patch(FancyBboxPatch((x, y - CH / 2), CW, CH,
                                             boxstyle="round,pad=0,rounding_size=0.05",
                                             linewidth=0,
                                             facecolor=CHIP_DARK if gate else CHIP_LIGHT))
                axD.text(x + CW / 2, y, code, ha="center", va="center", fontsize=5.2,
                         color="white" if gate else S.INK)
            if extra:
                axD.text(gi + 0.07 + 2 * (CW + 0.04) + 0.02, ii, f"+{extra}", ha="left",
                         va="center", fontsize=5.2, color=S.NEUTRAL_M)
    # row separators + group brackets
    for gname, lo, hi in groups:
        axD.axhline(hi - 0.5, color=S.NEUTRAL_L, linewidth=0.6)
        axD.text(-0.72, (lo + hi - 1) / 2, gname, rotation=90, va="center", ha="center",
                 fontsize=6.0, fontstyle="italic", color=S.INK)
    axD.set_xlim(-0.85, len(GENOMES))
    axD.set_ylim(len(ions) - 0.5, -0.9)
    axD.set_yticks(range(len(ions)))
    axD.set_yticklabels(ions, fontsize=6.6)
    axD.set_xticks([g + 0.5 for g in range(len(GENOMES))])
    axD.set_xticklabels([lab for _, lab in GENOMES], fontsize=6.8, fontstyle="italic")
    axD.xaxis.set_ticks_position("top")
    axD.tick_params(length=0, pad=2)
    for sp in ("top", "right", "bottom", "left"):
        axD.spines[sp].set_visible(False)

    lh = [FancyBboxPatch((0, 0), 1, 1, boxstyle="round,pad=0,rounding_size=0.1", linewidth=0,
                         facecolor=c) for c in (CHIP_DARK, CHIP_LIGHT)]
    axD.legend(lh, ["regulator with metal-coordination site",
                    "regulator predicted metal, site not detected"],
               loc="upper center", bbox_to_anchor=(0.5, -0.02), ncol=2, fontsize=6.2,
               handlelength=1.2, handleheight=1.2, borderpad=0.0, columnspacing=1.4)

    # ---------------------------------------------------------------- header
    fig.text(0.005, 0.985, "E", fontsize=11, fontweight="bold", va="top", color="black")
    fig.text(0.043, 0.983,
             "Predicted inducer per candidate regulator (predicted repertoire)", fontsize=8.2,
             va="top", color="black")
    fig.text(0.043, 0.951,
             f"metal call {n_metal}/{n_tot} ({100*n_metal/n_tot:.0f}%)  ·  "
             f"survey MetalNet CHDE {MS_CHDE_N}/{MS_CHDE_TOT} ({100*MS_CHDE_N/MS_CHDE_TOT:.0f}%)  ·  "
             f"each chip = one regulator, cf. panel D where each chip = one family",
             fontsize=6.4, va="top", color=S.NEUTRAL_D)

    for p in S.save(fig, OUT):
        print("wrote", p)
    print(f"\nmetal {n_metal}/{n_tot} ({100*n_metal/n_tot:.0f}%) vs MetalNet {MS_CHDE_N}/{MS_CHDE_TOT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
