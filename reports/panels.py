"""panels.py -- the composition panels Figure 2 (benchmark) and Figure 3 (survey) share.

One drawing routine per panel type, so the two benchmark genomes and the four survey genomes are
read in the same visual language: family x genome heat map, inducer-class composition, and the
per-genome metal panels (metal calls -> ion resolved -> operator supports the call). Every value
comes from a `composition` block of numbers.json.
"""
from __future__ import annotations

import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import MaxNLocator

from figure_style import (CLASS_COLOR, CLASS_LABEL, FS_BODY, FS_SMALL, C, ION_COLOR, OUTCOME_GLYPH,
                          ion_color, ion_label,
                          is_species)

CLASS_ORDER = ["metal", "redox", "organic", "undetermined"]
SHORT_FAM = {"LysR-type (LTTR)": "LysR"}
#: a shorter ion label for the narrow ion column, where the full one would overflow the bar
ION_SHORT = {"Cd2+/Pb2+": "Cd/Pb", "MoO42-": "MoO4", "As(III)": "As(III)"}
#: ions light enough to need dark text on top
_LIGHT_IONS = {"Zn2+", "Ni2+", "Co2+", "Mn2+"}

#: the four states of a census cell (discovered / known), in legend order
CENSUS_STATES = [("all", "all found", C["metal"], "white"),
                 ("some", "some found", C["metal_soft"], C["ink"]),
                 ("none", "none found", "white", C["ink"]),
                 ("absent", "family absent from the genome", C["empty"], C["ink"])]
_SURVEY_CMAP = LinearSegmentedColormap.from_list("survey", ["#E3ECF7", C["metal"]])
_OUTSIDE_CMAP = LinearSegmentedColormap.from_list("outside", [C["outside_soft"], C["outside"]])
_CELL = 0.9       # cell size in data units; the rest is the white gutter


def _cell(ax, x, y, fc, ec="none"):
    ax.add_patch(Rectangle((x - _CELL / 2, y - _CELL / 2), _CELL, _CELL, fc=fc, ec=ec, lw=0.5))


def family_heatmap(ax, families: list[str], rows: list[tuple[str, dict]], *,
                   totals: list[dict] | None = None, outside: list[dict] | None = None,
                   gap: float = 0.55) -> list[tuple[Patch, str]]:
    """Candidates per family (columns) and genome (rows): the twelve studied families, a gap, then
    every other family a candidate falls in, in a colour of its own.

    `rows` is [(label, composition)]. With `totals` (one {family: known} per row) each studied-family
    cell reads `found/known` and is coloured by its state (all / some / none found; family absent) --
    the benchmark census. Without it each cell is a count shaded by the count -- the survey. `outside`
    (one {family: n} per row) overrides the composition's `families_outside_study`. Empty cells are
    left blank. Returns the legend entries for the census states.
    """
    out_rows = outside if outside is not None else [c.get("families_outside_study", {}) for _, c in rows]
    tot_out = {}
    for o in out_rows:
        for f, n in o.items():
            tot_out[f] = tot_out.get(f, 0) + n
    extra = sorted(tot_out, key=lambda f: (-tot_out[f], f))
    xs = list(range(len(families))) + [len(families) + gap + k for k in range(len(extra))]
    found = np.array([[c["families"].get(f, 0) for f in families] for _, c in rows])
    vmax = max(int(found.max()), 1)
    omax = max([n for o in out_rows for n in o.values()] + [1])
    for i, (_, c) in enumerate(rows):
        for j, f in enumerate(families):
            d = int(found[i, j])
            if totals is not None:
                t = int(totals[i].get(f, 0))
                state = "absent" if not t else ("all" if d == t else ("some" if d else "none"))
                _, _, fc, tc = next(s for s in CENSUS_STATES if s[0] == state)
                _cell(ax, j, i, fc, ec=C["neutral_light"] if state == "none" else "none")
                if t:
                    ax.text(j, i, f"{d}/{t}", ha="center", va="center", fontsize=FS_SMALL, color=tc)
            elif d:
                v = d / vmax
                _cell(ax, j, i, _SURVEY_CMAP(0.08 + 0.92 * v))
                ax.text(j, i, str(d), ha="center", va="center", fontsize=FS_SMALL,
                        color="white" if v > 0.45 else C["ink"])
            else:
                _cell(ax, j, i, C["empty"])
        for k, f in enumerate(extra):
            n = int(out_rows[i].get(f, 0))
            x = xs[len(families) + k]
            if n:
                _cell(ax, x, i, _OUTSIDE_CMAP(0.15 + 0.85 * n / omax))
                ax.text(x, i, str(n), ha="center", va="center", fontsize=FS_SMALL, color=C["ink"])
            else:
                _cell(ax, x, i, C["empty"])
    # block headers above the columns
    head = dict(fontsize=FS_SMALL, color=C["neutral_dark"], va="bottom")
    ax.text(-0.5, -0.62, "the 12 studied families" + ("" if totals is not None else "  (blank: no candidate)"),
            ha="left", **head)
    if extra:
        ax.text(xs[len(families)] - 0.5, -0.62, "other families", ha="left",
                **head)
    ax.set_xticks(xs)
    ax.set_xticklabels([SHORT_FAM.get(f, f) for f in families + extra], rotation=35, ha="right",
                       fontsize=FS_BODY)
    for t in ax.get_xticklabels()[len(families):]:
        t.set_color(C["neutral_dark"])
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{lab} ({c['n']})" for lab, c in rows], fontsize=FS_BODY)
    for t, (lab, _) in zip(ax.get_yticklabels(), rows):
        t.set_style("italic" if is_species(lab) else "normal")
    ax.set_xlim(-0.55, xs[-1] + 0.55)
    ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.set_frame_on(False)
    ax.tick_params(length=0)
    return [(Patch(fc=fc, ec=C["neutral_light"] if key == "none" else "none", lw=0.5), lab)
            for key, lab, fc, _ in CENSUS_STATES]


def class_panel(ax, rows: list[tuple[str, dict]]) -> None:
    """Stacked inducer-class composition, one bar per genome; counts inside the segments."""
    bottoms = np.zeros(len(rows))
    for cls in CLASS_ORDER:
        counts = np.array([c["class_counts"].get(cls, 0) for _, c in rows])
        vals = 100 * counts / np.array([c["n"] for _, c in rows])
        ax.bar(range(len(rows)), vals, bottom=bottoms, width=0.66,
               color=CLASS_COLOR[cls], edgecolor="white", label=CLASS_LABEL[cls])
        for i, (v, n) in enumerate(zip(vals, counts)):
            if v >= 7:
                ax.text(i, bottoms[i] + v / 2, str(n), ha="center", va="center", fontsize=FS_SMALL,
                        color="white" if cls in ("metal", "redox") else C["ink"])
        bottoms += vals
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels([f"{lab}\n({c['n']})" for lab, c in rows], fontsize=FS_BODY)
    ax.set_ylim(0, 100)
    ax.set_ylabel("candidates (%)")


def class_legend(ax, **kw) -> None:
    ax.legend(handles=[Patch(fc=CLASS_COLOR[c], label=CLASS_LABEL[c]) for c in CLASS_ORDER],
              handlelength=1.0, columnspacing=1.0, **kw)


def _ion_order(ions) -> list[str]:
    return sorted(ions, key=lambda k: (list(ION_COLOR).index(k) if k in ION_COLOR else 99, k))


def metal_panels(fig, spec, rows: list[tuple[str, dict]], *, headroom: float = 1.28,
                 wspace: float = 0.35, label_ions: bool = True, ion_fs: float = 4.6) -> tuple[list, list[str]]:
    """One small panel per genome, each on its own y scale: metal calls | ion resolved (stacked by
    ion, each segment labelled) | operator supports the call (one colour).

    Every column carries its n/N: metal calls over the genome's candidates, ion-resolved over the
    metal calls, operator support over the ion-resolved calls (D3). Returns the axes and the ions
    drawn.
    """
    sub = spec.subgridspec(1, len(rows), wspace=wspace)
    axes, seen = [], {}
    for i, (lab, c) in enumerate(rows):
        ax = fig.add_subplot(sub[0, i])
        f, op = c["metal_funnel"], c["operator_supports_ion_call"]
        top = max(f["metal_calls"], 1)
        w = 0.74
        ax.bar(0, f["metal_calls"], width=w, color=C["metal"], edgecolor="white", linewidth=0.5)
        bottom = 0.0
        for ion in _ion_order(f["ions"]):
            n = f["ions"][ion]
            ax.bar(1, n, bottom=bottom, width=w, color=ion_color(ion), edgecolor="white", linewidth=0.5)
            if label_ions:
                ax.text(1, bottom + n / 2, ION_SHORT.get(ion, ion_label(ion)), ha="center", va="center",
                        fontsize=ion_fs, color=C["ink"] if ion in _LIGHT_IONS else "white")
            bottom += n
            seen[ion] = None
        ax.bar(2, op["n"], width=w, color=C["operator"], edgecolor="white", linewidth=0.5)
        for x, (a, b) in enumerate([(f["metal_calls"], c["n"]), (f["ion_resolved"], f["metal_calls"]),
                                    (op["n"], op["of_ion_resolved"])]):
            ax.text(x, a + top * 0.03, f"{a}/{b}", ha="center", va="bottom", fontsize=FS_SMALL)
        ax.set_xlim(-0.55, 2.55)
        ax.set_ylim(0, top * headroom)
        ax.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=4))
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(["metal", "ion", "oper."], fontsize=FS_SMALL)
        ax.tick_params(axis="x", length=0, pad=2)
        ax.tick_params(axis="y", labelsize=5.8)
        ax.set_title(lab, fontsize=FS_BODY, pad=3, style="italic" if is_species(lab) else "normal")
        if i == 0:
            ax.set_ylabel("regulators")
        axes.append(ax)
    return axes, _ion_order(seen)


def metal_key(ax, top_k: int, ions: list[str] | None = None, *, fontsize: float = FS_SMALL) -> None:
    """A blank axes holding the column key and what each column's n/N means. `ions` adds an ion
    legend (only needed when the segments are not labelled)."""
    ax.axis("off")
    handles = [Patch(fc=C["metal"], label="metal call"), Patch(fc=C["operator"], label="operator supports the call")]
    handles += [Patch(fc=ion_color(i), label=ion_label(i)) for i in ions or []]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(-0.05, 1.02), ncol=1,
              handlelength=0.9, handletextpad=0.4, fontsize=fontsize + 0.2)
    ax.text(-0.05, 0.02,
            "n/N above each column\n"
            "metal: metal calls / candidates\n"
            "ion: ion-resolved / metal calls,\n"
            "  each segment one ion\n"
            "operator: own-region operator in\n"
            f"  the top {top_k} / ion-resolved calls\n"
            "each genome on its own y scale",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=fontsize, color=C["neutral_dark"],
            linespacing=1.3)


def outcome_mark(ax, x: float, y: float, out: str, note: str = "", *, d: float = 6.4) -> None:
    """One outcome as a filled circle carrying its glyph; not applicable is a light dash.

    `d` is the circle diameter in points. Arial has no check/ballot glyphs, hence DejaVu Sans.
    """
    if out == "not_scored":
        ax.text(x, y, "–", ha="center", va="center", fontsize=d * 0.8, color=C["neutral_mid"])
        return
    ax.scatter([x], [y], s=d ** 2, color=C[out], lw=0, zorder=3)
    g = note or OUTCOME_GLYPH[out]
    ax.text(x, y, g, ha="center", va="center", fontsize=d * (0.62 if note else 0.72),
            fontfamily="DejaVu Sans", color="white" if out in ("correct", "incorrect") else C["ink"],
            zorder=4)


def outcome_legend_handles() -> tuple[list, list[str]]:
    from matplotlib.lines import Line2D
    mk = dict(marker="o", ls="none", markersize=5, markeredgewidth=0)
    return ([Line2D([], [], color=C["correct"], **mk), Line2D([], [], color=C["partial"], **mk),
             Line2D([], [], color=C["incorrect"], **mk),
             Line2D([], [], marker="$–$", ls="none", color=C["neutral_mid"], markersize=5)],
            ["yes / correct", "ion in the shortlist, not first", "no / incorrect", "not applicable"])
