"""build_survey_table.py -- Figure S5 (the survey's metal calls as a table) and Figure S5bis (every
survey candidate, one page per genome group).

A text figure: one block per survey genome, one row per regulator -- the family and SSN clade the
pipeline assigned, the inducer it called, its primary operator, whether an own-region operator
supports the call, and the gene that operator regulates. Every value is read from
numbers.json["survey_table"]. The survey genomes have no truth set, so nothing here is an accuracy
claim.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from figure_style import (CLASS_COLOR, CLASS_LABEL, FS_HEAD, C, MM, SURVEY_ORDER, W2, apply_style, ion_color,
                          ion_label, numbers, save, wrap)

MONO = "DejaVu Sans Mono"
#: table text on the shared type scale (figure_style.FS_SMALL); sequences a touch smaller, because a
#: monospace runs wider than Arial and a 38-bp operator must fit its column
FS, FS_SMALL, FS_MONO = 5.0, 5.0, 4.7
ROW, HEAD, GENOME_HEAD = 2.5, 6.5, 4.4      # mm
WIDTH = 183.0                               # mm
#: where a gene product starts inside the regulated-gene column (mm, and the label characters it allows)
PRODUCT_AT_MM, PRODUCT_AT_CHARS = 12.0, 13
#: Arial's average character width, in em, for fitting text to a column
CHAR_EM = 0.56
#: S5bis pages: which genomes share a page (kept in the survey order)
BIS_PAGES = [["M. tuberculosis"], ["M. avium"], ["V. cholerae", "V. vulnificus"]]


def _columns(with_class: bool) -> list[tuple[str, str, float]]:
    """(key, header, x in mm) in drawing order."""
    cols = [("tf", "TF", 0.0), ("family", "family", 15.5)]
    x = 32.5
    if with_class:
        cols.append(("cls", "class", x))
        x += 10.5
    cols += [("clade", "SSN clade", x), ("inducer", "inducer called", x + 19.5),
             ("operator", "primary operator (5'→3')", x + 49.0),
             ("support", "own-region\noperator", x + 91.0), ("target", "regulated gene", x + 101.5)]
    return cols


def _clip(text: str, n: int) -> str:
    return text if len(text) <= n else text[:n - 1] + "…"


def _inducer(ax, x, y, r: dict, n_chars: int) -> None:
    """The called inducer: a coloured ion, an unresolved metal with the channels that called it, or
    the ligand text."""
    if r["ion"]:
        ax.add_patch(Rectangle((x, y - 0.8), 1.6, 1.6, fc=ion_color(r["ion"]), ec="none"))
        ax.text(x + 2.4, y, ion_label(r["ion"]), fontsize=FS, va="center", ha="left")
    elif r["metal_call"]:
        why = ", ".join(r["metal_channels"])
        short = ", ".join(r["shortlist"])
        ax.text(x, y, "unresolved" + (f": {short}" if short else ""), fontsize=FS, va="center", ha="left")
        if why and not short:
            ax.text(x + 11.3, y, f"({why})", fontsize=FS_SMALL, va="center", ha="left",
                    color=C["neutral_dark"])
    else:
        # "non-metal effector (undetermined)" is a non-metal call that names no ligand
        vague = "undetermined" in r["inducer"] or r["cls"] == "undetermined"
        text = "non-metal, no ligand" if vague else r["inducer"]
        ax.text(x, y, _clip(text, n_chars), fontsize=FS, va="center", ha="left",
                color=C["neutral_dark"] if vague else C["ink"])


def _support(ax, x, y, r: dict) -> None:
    if not r["ion"]:
        ax.text(x, y, "·", fontsize=FS + 1, va="center", ha="center", color=C["neutral_mid"])
    elif r["own_region_top_k"]:
        ax.text(x, y, f"yes (#{r['own_region_rank']})", fontsize=FS, va="center", ha="center",
                color=C["ink"], fontweight="bold")
    else:
        ax.text(x, y, "no", fontsize=FS, va="center", ha="center", color=C["neutral_dark"])


def _target(ax, x, y, r: dict, col_mm: float) -> None:
    t = r["target"]
    if not t["genes"]:
        ax.text(x, y, "–", fontsize=FS, va="center", ha="left", color=C["neutral_mid"])
        return
    label = " / ".join(g["label"] for g in t["genes"])
    style = "italic" if t["rule"] == "nearest" else "normal"
    ax.text(x, y, label, fontsize=FS, va="center", ha="left", style=style,
            fontweight="bold" if any(g["self"] for g in t["genes"]) and len(t["genes"]) == 1 else "normal")
    # a product helps only when the gene has no name of its own
    unnamed = [g for g in t["genes"] if not g["self"] and not g["named"] and g["product"]]
    # the product sits in a fixed sub-column, so it never collides with a long label
    room = int((col_mm - PRODUCT_AT_MM - 1.5) / (FS_SMALL * CHAR_EM * 25.4 / 72))
    if unnamed and len(label) <= PRODUCT_AT_CHARS - 2 and room > 8:
        prod = unnamed[0]["product"].replace("hypothetical protein", "hypothetical")
        ax.text(x + PRODUCT_AT_MM, y, _clip(prod, room), fontsize=FS_SMALL, va="center", ha="left",
                color=C["neutral_dark"])


def _row(ax, cols, y, r: dict, stripe: bool) -> None:
    if stripe:
        ax.add_patch(Rectangle((-0.5, y - ROW / 2), WIDTH + 0.5, ROW, fc="#F4F6F7", ec="none", zorder=0))
    xs = {k: x for k, _, x in cols}
    keys = [k for k, _, _ in cols]
    nxt = {k: (xs[keys[i + 1]] if i + 1 < len(keys) else WIDTH) for i, k in enumerate(keys)}

    def room(k):  # characters that fit in a column at FS
        return int((nxt[k] - xs[k] - 1.0) / (FS * CHAR_EM * 25.4 / 72))

    ax.text(xs["tf"], y, r["tf"], fontsize=FS, va="center", ha="left",
            style="italic" if r["named"] else "normal")
    ax.text(xs["family"], y, r["family"].replace("LysR-type (LTTR)", "LysR"), fontsize=FS, va="center",
            ha="left")
    if "cls" in xs:
        ax.add_patch(Rectangle((xs["cls"], y - 0.8), 1.6, 1.6, fc=CLASS_COLOR[r["cls"]], ec="none"))
        ax.text(xs["cls"] + 2.4, y, CLASS_LABEL[r["cls"]].replace("undetermined", "undet."), fontsize=FS,
                va="center", ha="left")
    ax.text(xs["clade"], y, r["clade"] or "–", fontsize=FS, va="center", ha="left",
            color=C["ink"] if r["clade"] else C["neutral_mid"])
    _inducer(ax, xs["inducer"], y, r, room("inducer"))
    op = r["operator"]
    ax.text(xs["operator"], y, op["seq"] if op else "–", fontsize=FS_MONO, va="center", ha="left",
            fontfamily=MONO)
    _support(ax, xs["support"] + 4.0, y, r)
    _target(ax, xs["target"], y, r, WIDTH - xs["target"])


def table_figure(N: dict, genomes: list[str], rows: list[dict], *, with_class: bool, stem: str,
                 title: str, note: str) -> None:
    cols = _columns(with_class)
    by_genome = {g: [r for r in rows if r["genome"] == g] for g in genomes}
    # wrap() assumes 0.5 em per character; Arial runs wider, so wrap to a narrower measure
    foot = [wrap(t, WIDTH * MM * 0.86, FS_SMALL) for t in (N["survey_table"]["target_rule"], note)]
    n_foot = sum(f.count("\n") + 1 for f in foot)
    height = 1.5 + HEAD + sum(GENOME_HEAD + ROW * len(v) for v in by_genome.values()) + 2.4 * n_foot + 2.0
    fig = plt.figure(figsize=(W2 * 0.975, height * MM))  # the crop adds a margin: stay inside 183 mm
    ax = fig.add_axes([0, 0, 1, 1])
    # a little canvas beyond the content, so the cropped export stays inside 183 mm
    ax.set_xlim(-1.5, WIDTH + 1.5)
    ax.set_ylim(height, 0)
    ax.axis("off")
    # no in-figure title: the document prints it
    y = 1.5
    for k, lab, x in cols:
        center = k == "support"
        ax.text(x + (4.0 if center else 0), y + HEAD - 1.2, lab, fontsize=FS, fontweight="bold",
                va="bottom", ha="center" if center else "left", color=C["neutral_dark"], linespacing=1.0)
    y += HEAD
    ax.plot([-0.5, WIDTH], [y - 0.4] * 2, color=C["ink"], lw=0.6)
    for g, rs in by_genome.items():
        n_all = sum(r["genome"] == g for r in N["survey_table"]["rows"])
        n_metal = sum(r["genome"] == g and r["metal_call"] for r in N["survey_table"]["rows"])
        ax.text(0, y + GENOME_HEAD / 2 + 0.3, g, fontsize=FS_HEAD, fontweight="bold", style="italic", va="center")
        ax.text(len(g) * 6.4 * 0.52 * 25.4 / 72 + 2.0, y + GENOME_HEAD / 2 + 0.3,
                f"{n_metal} metal calls among {n_all} candidates", fontsize=FS, va="center",
                color=C["neutral_dark"])
        y += GENOME_HEAD
        for i, r in enumerate(rs):
            _row(ax, cols, y + ROW / 2, r, stripe=i % 2 == 1)
            y += ROW
        ax.plot([-0.5, WIDTH], [y + 0.1] * 2, color=C["neutral_light"], lw=0.5)
    y += 1.6
    for f in foot:
        ax.text(0, y, f, fontsize=FS_SMALL, va="top", ha="left", color=C["neutral_dark"], linespacing=1.3)
        y += 2.4 * (f.count("\n") + 1)
    save(fig, stem)


def main() -> int:
    apply_style()
    N = numbers()
    T = N["survey_table"]
    genomes = [g for g in SURVEY_ORDER if any(r["genome"] == g for r in T["rows"])]
    k = T["top_k"]
    note = (f"Operator: the primary operator (rank-1 rescan hit) on its own strand. Own-region operator: one "
            f"of the top {k} operators lies in the regulator's own flanking intergenic region (rank given), "
            "so the operator supports the call; scored only for calls that name one ion (·). Unresolved: a "
            "metal call naming no single ion, with the channels that called metal. No truth set exists "
            "for these genomes: nothing here is an accuracy claim.")
    metal = [r for r in T["rows"] if r["metal_call"]]
    table_figure(N, genomes, metal, with_class=False, stem="figS5_survey_metal_table",
                 title=f"The survey's {len(metal)} metal calls", note=note)
    for i, page in enumerate(BIS_PAGES, 1):
        gs = [g for g in page if g in genomes]
        rows = [r for r in T["rows"] if r["genome"] in gs]
        table_figure(N, gs, rows, with_class=True, stem=f"figS5bis_survey_table_{i}",
                     title=f"Every survey candidate ({i}/{len(BIS_PAGES)}): {', '.join(gs)}", note=note)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
