"""build_benchmark_figures.py -- Figure 2 (benchmark) and Figure S1 (every regulator, by axis).

Every axis is scored over ALL discovered TFs that carry truth for that axis, so denominators differ
by axis and are printed on every mark. Accuracy is claimed only here, for the two genomes with a
knowledge base. Panels a-c draw the benchmark genomes with the same routines as the survey figure
(panels.py), so the two figures read in one visual language. Every number is read from numbers.json.
"""
from __future__ import annotations

import cnsplots as cns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

from figure_style import (BENCH_GENOMES, CLASS_LABEL, FS_BODY, FS_HEAD, FS_SMALL, C, SHORT, W2, apply_style,
                          numbers, panel_head, save, type_scale, wrap)
from panels import (SHORT_FAM, class_legend, class_panel, family_heatmap, metal_panels, outcome_legend_handles,
                    outcome_mark)

MISS_REASON = {"profile_scope": "no clade profile reaches it", "threshold": "below the e-value cut",
               "curation": "dropped in curation"}


def _pct(a: int, b: int) -> float:
    return 100 * a / b if b else 0.0


def recovery_headline(N: dict) -> str:
    """The "did we miss any" answer as one sentence, every count read from numbers.json."""
    acct = N["benchmark"]["metal_sensor_accounting"]
    r, s = acct["ecoli_recovery"], acct["salmonella_recovery"]
    miss = "; ".join(f"missed: {tf}, {MISS_REASON.get(why, why.replace('_', ' '))}"
                     for tf, why in r["missed"].items())
    txt = (f"E. coli: {r['found']} of {r['in_study_families']} known metal sensors in the twelve "
           f"studied families recovered" + (f" ({miss})" if miss else ""))
    out = r["outside_study_families"]
    if out:
        txt += f"; {len(out)} more sit outside the twelve families"
        if r["outside_found"]:
            txt += f" ({len(r['outside_found'])} found anyway: {', '.join(r['outside_found'])})"
    return txt + (f".  SL1344: {s['metal_class']} of {s['sensors']} Osman sensors called metal, "
                  f"{s['exact_ion']} with the exact ion.")


def census_rows(N: dict) -> tuple[list[tuple[str, dict]], list[dict], list[dict]]:
    """Discovered/known per studied family plus the discovered TFs outside them, per genome, shaped
    like a composition block so the survey heat map draws it. Row n = every discovered TF."""
    rows, totals, outside = [], [], []
    for key, label in BENCH_GENOMES:
        d = N["benchmark"][key]["scores"]["discovery"]
        rows.append((SHORT[label], {"families": {f: v["discovered"] for f, v in d["by_family"].items()},
                                    "n": d["n_candidates"]}))
        totals.append({f: v["total"] for f, v in d["by_family"].items()})
        outside.append({f: len(v) for f, v in N["benchmark"][key]["discovered_outside_census"].items()})
    return rows, totals, outside


def panel_accuracy(ax, N: dict) -> None:
    """Accuracy per axis, each over the discovered TFs with truth for that axis."""
    e, s = N["benchmark"]["ecoli"]["scores"], N["benchmark"]["salmonella"]["scores"]

    def op(sc, k):
        b = sc["operator"]["all"]
        return b[k]["observed"], b["n"]

    axes_def = [("family", lambda sc: (sc["family"]["correct"], sc["family"]["n"])),
                ("inducer\nclass", lambda sc: (sc["class"]["correct"], sc["class"]["n"])),
                ("metal ion\n(exact)", lambda sc: (sc["ion"]["hit"], sc["ion"]["n"])),
                ("primary\noperator", lambda sc: op(sc, "1"))]
    x = np.arange(len(axes_def))
    w = 0.38
    for off, sc, key, lab in ((-w / 2, e, "ecoli", "E. coli"), (w / 2, s, "salm", "SL1344")):
        vals = [f(sc) for _, f in axes_def]
        ax.bar(x + off, [_pct(a, b) for a, b in vals], width=w, color=C[key], edgecolor="white", label=lab)
        for xi, (a, b) in zip(x + off, vals):
            ax.text(xi, _pct(a, b) + 2, f"{a}/{b}", ha="center", va="bottom", fontsize=FS_SMALL, rotation=90)
    ax.set_xticks(x)
    ax.set_xticklabels([lbl for lbl, _ in axes_def], fontsize=FS_BODY)
    ax.set_ylim(0, 132)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylabel("correct (% of TFs with truth)")
    ax.legend(loc="upper right", handlelength=1.0, ncol=2)


def panel_confusion(ax, N: dict) -> tuple[int, int]:
    """Class confusion over every discovered TF with a known effector, both genomes pooled."""
    rows = [(r["truth_class"], r["pred_class"]) for g in ("ecoli", "salmonella")
            for r in N["benchmark"][g]["per_tf"] if r.get("truth_class")]
    df = pd.DataFrame(rows, columns=["known", "predicted"]).replace(CLASS_LABEL)
    present = set(df.known) | set(df.predicted)
    order = [CLASS_LABEL[c] for c in ("metal", "redox", "organic", "undetermined") if CLASS_LABEL[c] in present]
    cns.confusionplot(data=df, x="predicted", y="known", x_order=order, y_order=order, cmap="Blues", ax=ax)
    ax.set_xlabel("predicted class")
    ax.set_ylabel("known class")
    return sum(1 for k, p in rows if k == p), len(rows)


# ---------------------------------------------------------------- panel f: metal-sensor accounting
#: column positions of the accounting table, in mm from the block's left edge
ACC_COLS = {"tf": 0.0, "ions": 9.0, "family": 29.0, "found": 46.0, "metal": 53.0, "ion": 59.5,
            "note": 63.5}


def _sensor_row(s: dict) -> tuple[list[str], str]:
    """(outcomes for found / metal class / exact ion, note) for one known metal sensor."""
    studied = s["family"] != "outside the 12"
    if s["discovered"] != "yes":
        found = "incorrect" if studied else "not_scored"
        note = MISS_REASON.get(s["miss_reason"], s["miss_reason"]) if studied else "family not studied"
        return [found, "not_scored", "not_scored"], note
    metal = "correct" if s["ion"] in ("hit", "shortlist", "class only") else "incorrect"
    ion = {"hit": "correct", "shortlist": "partial", "class only": "partial"}.get(s["ion"], "incorrect")
    note = "ion in the shortlist, not first" if s["ion"] == "shortlist" else \
        (f"called {s['called']}" if s["called"] else "")
    if not studied:
        note += " (found anyway)"
    return ["correct", metal, ion], note


def _accounting_block(ax, groups: list[tuple[str, list[dict]]], title: str, width_mm: float) -> None:
    """One genome's metal sensors as a compact table: TF, known ion(s), family, then found / metal
    class / exact ion as outcome marks and a short note. `groups` is [(heading, rows)]; each row
    is {tf, ions, family, outs, note, grey}."""
    head = dict(fontsize=FS_SMALL, color=C["neutral_dark"], va="bottom", fontweight="bold")
    y = -0.35
    for key, lab in (("tf", "TF"), ("ions", "known ion(s)"), ("family", "family"), ("note", "call / note")):
        ax.text(ACC_COLS[key], y, lab, ha="left", **head)
    for key, lab in (("found", "found"), ("metal", "metal"), ("ion", "ion")):
        ax.text(ACC_COLS[key], y, lab, ha="center", **head)
    ax.plot([0, width_mm], [y + 0.1] * 2, color=C["neutral_light"], lw=0.5)
    y = 0.35
    for heading, rows in groups:
        ax.text(0, y, heading, ha="left", va="center", fontsize=FS_SMALL, style="italic", color=C["neutral_dark"])
        y += 0.85
        for r in rows:
            ink = C["neutral_dark"] if r.get("grey") else C["ink"]
            ax.text(ACC_COLS["tf"], y, r["tf"], ha="left", va="center", fontsize=FS_SMALL, fontweight="bold",
                    color=ink)
            ax.text(ACC_COLS["ions"], y, r["ions"], ha="left", va="center", fontsize=FS_SMALL, color=ink)
            ax.text(ACC_COLS["family"], y, r["family"], ha="left", va="center", fontsize=FS_SMALL, color=ink,
                    style="italic" if r.get("grey") else "normal")
            for key, out in zip(("found", "metal", "ion"), r["outs"]):
                outcome_mark(ax, ACC_COLS[key], y, out, d=5.6)
            ax.text(ACC_COLS["note"], y, r["note"], ha="left", va="center", fontsize=FS_SMALL,
                    color=C["neutral_dark"])
            y += 1
        y += 0.25
    ax.set_xlim(-0.5, width_mm)
    ax.set_ylim(y - 0.3, -1.3)
    ax.axis("off")
    ax.set_title(title, loc="left", fontsize=FS_HEAD, fontweight="bold", pad=2)


def accounting_groups(N: dict) -> dict[str, list[tuple[str, list[dict]]]]:
    acct = N["benchmark"]["metal_sensor_accounting"]
    fam_of = {r["tf"]: r["truth_family"] for r in N["benchmark"]["ecoli"]["per_tf"]}

    def row(s):
        outs, note = _sensor_row(s)
        return {"tf": s["tf"], "ions": s["known_ions"].replace(";", ", "),
                "family": SHORT_FAM.get(s["family_name"], s["family_name"]),
                "outs": outs, "note": note, "grey": s["family"] == "outside the 12"}

    order = {"yes": 0, "no": 1}
    eco = sorted((s for s in acct["sensors"] if s["genome"] == "E. coli"),
                 key=lambda s: (order[s["discovered"]], s["tf"]))
    sal = sorted((s for s in acct["sensors"] if s["genome"] == "SL1344"), key=lambda s: s["tf"])
    fm = [{"tf": f["tf"], "ions": "none", "family": fam_of.get(f["tf"], ""),
           "outs": ["correct", "incorrect", "not_scored"],
           "note": f"false metal call; agreement {f['agreement']}"}
          for f in acct["false_metal"] if f["genome"] == "E. coli"]
    e_groups = [("in the 12 studied families", [row(s) for s in eco if s["family"] != "outside the 12"]),
                ("outside the 12 (family from Pfam)", [row(s) for s in eco if s["family"] == "outside the 12"])]
    if fm:
        e_groups.append(("not a metal sensor, called metal", fm))
    return {"ecoli": e_groups, "salmonella": [("in the 12 studied families", [row(s) for s in sal])]}


def fig_benchmark(N: dict, stem: str = "fig2_benchmark") -> None:
    """Figure 2. No in-figure title or caveat paragraph: the document prints the title, and the legend
    carries the recovery counts and the independence caveat (as Figure 1 does)."""
    with plt.rc_context(type_scale()):
        fig = plt.figure(figsize=(W2, 6.9))
        # panel a needs a deeper gap below it (rotated family names + the state legend) than the rows
        # below need between them, so the outer grid carries that gap and an inner grid the rest
        outer = fig.add_gridspec(2, 1, height_ratios=[0.6, 1.95], hspace=0.5,
                                 left=0.095, right=0.985, top=0.955, bottom=0.01)
        top_gs = outer[0].subgridspec(1, 3, wspace=0.5)
        low = outer[1].subgridspec(2, 3, height_ratios=[0.95, 1.0], hspace=0.42, wspace=0.5)

        ax = fig.add_subplot(top_gs[0, :2])
        rows, totals, outside = census_rows(N)
        handles = family_heatmap(ax, N["families"], rows, totals=totals, outside=outside)
        ax.legend([h for h, _ in handles], [lab for _, lab in handles], loc="upper left",
                  bbox_to_anchor=(-0.02, -0.6), ncol=4, handlelength=1.0, columnspacing=1.0)
        panel_head(ax, "a", "Discovery: found / known per studied family, and every other TF found",
                   pad_pt=12)

        bench = [(SHORT[lab], N["benchmark_genome_composition"][k]) for k, lab in BENCH_GENOMES]
        ax = fig.add_subplot(top_gs[0, 2])
        class_panel(ax, bench)
        class_legend(ax, loc="upper left", bbox_to_anchor=(-0.1, -0.36), ncol=2)
        panel_head(ax, "b", "Inducer class", pad_pt=12, letter_dx_pt=-28)

        axes, _ = metal_panels(fig, low[0, 0], bench, wspace=0.45)
        axes[0].legend(handles=[Patch(fc=C["metal"], label="metal call"),
                                Patch(fc=C["operator"], label="operator supports it")],
                       loc="upper left", bbox_to_anchor=(-0.45, -0.14), ncol=2, handlelength=0.9,
                       columnspacing=0.8)
        panel_head(axes[0], "c", "Metal calls, ions, operators", pad_pt=14, letter_dx_pt=-28)

        ax = fig.add_subplot(low[0, 1])
        panel_accuracy(ax, N)
        panel_head(ax, "d", "Accuracy, axis by axis", pad_pt=14, letter_dx_pt=-28)

        ax = fig.add_subplot(low[0, 2])
        diag, n = panel_confusion(ax, N)
        panel_head(ax, "e", f"Class calls: {diag}/{n} correct", pad_pt=14, letter_dx_pt=-40)

        groups = accounting_groups(N)
        n_sal = sum(len(r) for _, r in groups["salmonella"])
        sub = low[1, :].subgridspec(1, 2, wspace=0.08, width_ratios=[1.05, 1.0])
        ax = fig.add_subplot(sub[0, 0])
        n_known = sum(len(r) for h, r in groups["ecoli"] if not h.startswith("not a metal"))
        _accounting_block(ax, groups["ecoli"], f"E. coli K-12: all {n_known} known metal sensors", 95)
        panel_head(ax, "f", pad_pt=2)
        ax2 = fig.add_subplot(sub[0, 1])
        _accounting_block(ax2, groups["salmonella"], f"SL1344: the {n_sal} Osman 2019 sensors", 95)
        # both blocks share a row pitch: pad the shorter one's y range to the taller one's
        lo = max(ax.get_ylim()[0], ax2.get_ylim()[0])
        ax.set_ylim(lo, -1.3)
        ax2.set_ylim(lo, -1.3)
        h, lab = outcome_legend_handles()
        ax2.legend(h, lab, loc="lower left", bbox_to_anchor=(0.0, 0.02), ncol=2, handlelength=1.0,
                   handletextpad=0.4)
        save(fig, stem)


# ---------------------------------------------------------------- Figure S1
S1_COLS = [("family", "family"), ("class", "inducer\nclass"), ("ion", "metal\nion"), ("op1", "primary\noperator")]


def _outcome(row: dict, axis: str) -> tuple[str, str]:
    """(outcome, annotation) for one TF on one axis."""
    if axis == "family":
        v = row.get("family_ok")
        return ({"True": "correct", "False": "incorrect"}.get(v, "not_scored"), "")
    if axis == "class":
        v = row.get("class_ok")
        return ({"True": "correct", "False": "incorrect"}.get(v, "not_scored"), "")
    if axis == "ion":
        v = row.get("ion", "")
        return ({"hit": "correct", "shortlist": "partial", "class only": "partial",
                 "wrong ion": "incorrect", "miss": "incorrect"}.get(v, "not_scored"), "")
    if axis in ("op1", "op5"):
        v1, v5 = row.get("op_rec_1", ""), row.get("op_rec_5", "")
        if v1 == "":
            return "not_scored", ""
        rank = row.get("op_best_rank") or ""
        if axis == "op1":
            return ("correct" if v1 == "True" else "incorrect"), ""
        # rank shown only when it adds information beyond the top-1 column
        return (("correct" if v5 == "True" else "incorrect"),
                (f"#{rank}" if rank and v5 == "True" and str(rank) != "1" else ""))
    return "not_scored", ""


def has_truth_beyond_family(r: dict) -> bool:
    return any(_outcome(r, k)[0] != "not_scored" for k in ("class", "ion", "op1"))


def _matrix(ax, rows: list[dict], title: str, n_rows: int) -> None:
    """One mark per TF and axis, one row per TF with truth beyond its family."""
    dx = 1.0
    for i, r in enumerate(rows):
        for j, (key, _) in enumerate(S1_COLS):
            out, note = _outcome(r, key)
            outcome_mark(ax, j * dx, i, out, note, d=5.8)
        name = r["tf"] if len(r["tf"]) <= 16 else r["tf"][:15] + "…"
        star = "*" if r.get("anchored") == "yes" else ""
        weight = "bold" if r.get("truth_class") == "metal" else "normal"
        ax.text(-0.55, i, name + star, ha="right", va="center", fontsize=FS_SMALL, fontweight=weight)
    for j, (_, lbl) in enumerate(S1_COLS):
        ax.text(j * dx, -0.9, lbl, ha="center", va="bottom", fontsize=FS_SMALL, linespacing=1.05)
    ax.set_xlim(-2.9, (len(S1_COLS) - 1) * dx + 0.6)
    ax.set_ylim(n_rows - 0.4, -2.6)
    ax.axis("off")
    ax.set_title(title, loc="left", fontsize=FS_HEAD, fontweight="bold", pad=2)


def _order(rows: list[dict]) -> list[dict]:
    rank = {"metal": 0, "redox": 1, "organic": 2, "": 3}
    return sorted(rows, key=lambda r: (rank.get(r.get("truth_class", ""), 3), r["tf"].lower()))


def _family_only(ax, blocks: list[tuple[str, list[dict]]], width_in: float) -> None:
    """The TFs with no truth beyond their family, as name lists with the family tally."""
    ax.axis("off")
    y = 1.0
    for title, rows in blocks:
        ok = sum(r.get("family_ok") == "True" for r in rows)
        scored = sum(r.get("family_ok") in ("True", "False") for r in rows)
        unscored = [r["tf"] for r in rows if r.get("family_ok") not in ("True", "False")]
        ax.text(0, y, f"{title}: {len(rows)} TFs, family correct {ok}/{scored}", transform=ax.transAxes,
                fontsize=FS_BODY, fontweight="bold", va="top")
        y -= 0.075
        names = ", ".join(sorted((r["tf"] for r in rows if r.get("family_ok") in ("True", "False")),
                                 key=str.lower))
        txt = wrap(names, width_in, 5.3)
        ax.text(0, y, txt, transform=ax.transAxes, fontsize=FS_SMALL, va="top", linespacing=1.25)
        y -= 0.052 * (txt.count("\n") + 1) + 0.02
        if unscored:
            txt = wrap("family not scored (no census entry): " + ", ".join(unscored), width_in, 5.3)
            ax.text(0, y, txt, transform=ax.transAxes, fontsize=FS_SMALL, va="top", linespacing=1.25,
                    color=C["neutral_dark"])
            y -= 0.052 * (txt.count("\n") + 1) + 0.02
        y -= 0.05


def fig_per_regulator(N: dict) -> None:
    eco_all = _order(N["benchmark"]["ecoli"]["per_tf"])
    sal_all = sorted(N["benchmark"]["salmonella"]["per_tf"],
                     key=lambda r: (r.get("osman_sensor") != "yes", r["tf"].lower()))
    eco = [r for r in eco_all if has_truth_beyond_family(r)]
    sal = [r for r in sal_all if has_truth_beyond_family(r)]
    eco_f = [r for r in eco_all if not has_truth_beyond_family(r)]
    sal_f = [r for r in sal_all if not has_truth_beyond_family(r)]
    n = len(eco)
    height = 0.108 * (n + 3) + 0.85
    fig = plt.figure(figsize=(W2, height))
    top, bottom = 1 - 0.28 / height, 0.42 / height
    gs = fig.add_gridspec(1, 2, width_ratios=[0.8, 1.25], wspace=0.1, left=0.02, right=0.985,
                          top=top, bottom=bottom)
    _matrix(fig.add_subplot(gs[0, 0]), eco, f"E. coli K-12: the {len(eco)} of {len(eco_all)} TFs with truth "
            "beyond family", n + 3)
    right = gs[0, 1].subgridspec(2, 1, height_ratios=[len(sal) + 3, n - len(sal)], hspace=0.08)
    _matrix(fig.add_subplot(right[0, 0]), sal, f"SL1344: the {len(sal)} of {len(sal_all)} TFs with truth "
            "beyond family", len(sal) + 3)
    ax = fig.add_subplot(right[1, 0])
    ax.text(0, 1.02, "Family call only: no class, ion or operator truth", transform=ax.transAxes,
            fontsize=FS_HEAD, fontweight="bold", va="bottom")
    _family_only(ax, [("E. coli K-12", eco_f), ("SL1344", sal_f)], W2 * 0.5)
    h, _ = outcome_legend_handles()
    fig.legend(h, ["correct", "partial: right class, ion in the shortlist but not first", "incorrect",
                   "not scored: no truth for this TF on this axis"],
               loc="lower left", bbox_to_anchor=(0.02, 0.2 / height), ncol=4, handlelength=1.0,
               handletextpad=0.3, columnspacing=1.2)
    n_osman = N["benchmark"]["metal_sensor_accounting"]["salmonella_recovery"]["sensors"]
    fig.text(0.02, 0.06 / height, "Bold: metal sensors.  *: the TF's own sequence anchors a shipped SSN "
             "clade (curated support, not an independent prediction).  SL1344: class, ion and operator "
             f"truth exist only for the {n_osman} Osman 2019 sensors.", fontsize=FS_SMALL, color=C["neutral_dark"],
             va="bottom")
    save(fig, "figS1_per_regulator")


def main() -> int:
    apply_style()
    N = numbers()
    fig_benchmark(N)
    fig_per_regulator(N)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
