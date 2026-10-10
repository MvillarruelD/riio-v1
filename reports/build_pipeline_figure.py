"""Compact, shared Python renderer for Figure 1 and the run-scoped F0/F1 schematics.

Figure contract: RIIO connects four inputs to per-TF context, operator/inducer evidence,
and reviewable reports or AF3 submissions. This schematic-led composite merges biological
state drawings with the implemented v1.1 routing; icons illustrate states, not measured data.
Python only; NMI pastels; 177.8 mm canvas, >=7 pt text, editable SVG/PDF and 600 dpi rasters.
No prediction accuracy or neutrality claim is made. CLI/API syntax belongs in the legend.
"""
from __future__ import annotations

import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import (
    Circle,
    FancyArrowPatch,
    FancyBboxPatch,
    Polygon,
    Rectangle,
)

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42

INK, MUTED, RULE = "#272727", "#5E5E5E", "#BDBDC6"
ENTRY, ENTRY_SOFT = "#484878", "#E4E4F0"
OPER, OPER_SOFT = "#7884B4", "#E8ECF6"
IND, IND_SOFT = "#875F78", "#E4CCD8"
OUT, OUT_SOFT = "#484878", "#F0E9F0"
FS_TITLE, FS_LETTER, FS_HEAD, FS_BOX, FS_BODY = 9.0, 8.5, 8.0, 7.3, 7.0
W_MM = 177.8
MODES = [
    ("Single protein", "sequence, FASTA\nor protein accession"),
    ("Protein list", "multi-FASTA\nor accession table"),
    ("Genome accession", "RefSeq / GenBank\nassembly accession"),
    ("Local genome", "FASTA + GFF\nGenBank / proteome"),
]


class _Canvas:
    """Coordinates in mm, measured from the top left of the printed figure."""

    def __init__(self, height_mm):
        self.h = height_mm
        self.fig, self.ax = plt.subplots(figsize=(W_MM / 25.4, height_mm / 25.4))
        self.fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        self.ax.set(xlim=(0, W_MM), ylim=(0, height_mm))
        self.ax.axis("off")

    def y(self, top):
        return self.h - top

    def text(self, x, top, s, *, size=FS_BODY, weight="normal", color=INK, ha="left"):
        return self.ax.text(x, self.y(top), s, fontsize=size, fontweight=weight, color=color,
                            ha=ha, va="center", linespacing=1.1, zorder=5)

    def box(self, x, top, w, h, *, fill="white", edge=ENTRY, dashed=False, z=1):
        self.ax.add_patch(FancyBboxPatch(
            (x, self.y(top + h)), w, h, boxstyle="round,pad=0,rounding_size=1.2",
            facecolor=fill, edgecolor=edge, linewidth=0.6,
            linestyle=(0, (3, 2)) if dashed else "-", zorder=z))

    def heading(self, top, letter, title, *, note=""):
        self.text(2, top, letter, size=FS_LETTER, weight="bold")
        self.text(7, top, title, size=FS_HEAD, weight="bold")
        if note:
            self.text(W_MM - 2, top, note, color=MUTED, ha="right")

    def line(self, points, color=RULE, lw=0.6):
        xs, ys = zip(*[(x, self.y(top)) for x, top in points])
        self.ax.plot(xs, ys, color=color, lw=lw, zorder=2)

    def arrow(self, points, *, color=ENTRY, dashed=False):
        pts = [(x, self.y(top)) for x, top in points]
        if len(pts) > 2:
            xs, ys = zip(*pts[:-1])
            self.ax.plot(xs, ys, color=color, linewidth=0.75, zorder=2)
        self.ax.add_patch(FancyArrowPatch(
            pts[-2], pts[-1], arrowstyle="-|>", mutation_scale=7, color=color,
            linewidth=0.75, linestyle=(0, (3, 2)) if dashed else "-",
            shrinkA=0, shrinkB=0.6, zorder=2))

    def gene(self, x0, x1, top, *, reverse=False, color=RULE, h=0.85):
        head = min(1.0, (x1 - x0) * 0.35)
        yy = self.y(top)
        if reverse:
            pts = [(x1, yy - h), (x0 + head, yy - h), (x0, yy),
                   (x0 + head, yy + h), (x1, yy + h)]
        else:
            pts = [(x0, yy - h), (x1 - head, yy - h), (x1, yy),
                   (x1 - head, yy + h), (x0, yy + h)]
        self.ax.add_patch(Polygon(pts, facecolor=color, edgecolor=color, lw=0.45, zorder=4))

    def locus(self, x, top, w, *, color=ENTRY, promoter=False):
        self.line([(x, top), (x + w, top)])
        self.gene(x, x + w * 0.25, top)
        self.gene(x + w * 0.34, x + w * 0.56, top, reverse=True, color=color)
        self.gene(x + w * 0.72, x + w, top)
        if promoter:
            self.ax.add_patch(Rectangle((x + w * 0.58, self.y(top + 1.3)), w * 0.1, 2.6,
                                       facecolor=OPER, edgecolor="none", zorder=3))

    def proteins(self, x, top, *, n=4, color=ENTRY):
        for i, w in enumerate([9, 6.5, 8, 5.5][:n]):
            self.line([(x, top + i * 1.6), (x + w, top + i * 1.6)], color, lw=1.7)

    def ssn(self, x, top, *, color=ENTRY):
        for off, active in [(-4, False), (0, True), (4, False)]:
            pts = [(x + off + 1.5 * math.cos(k * math.pi / 3),
                    top + 1.8 * math.sin(k * math.pi / 3)) for k in range(6)]
            for k, p in enumerate(pts):
                self.line([p, pts[(k + 1) % 6]], color if active else RULE, lw=0.4)
                self.line([p, pts[(k + 2) % 6]], color if active else RULE, lw=0.4)
                self.ax.add_patch(Circle((p[0], self.y(p[1])), 0.3,
                                         facecolor=color if active else "white",
                                         edgecolor=color if active else RULE, lw=0.4, zorder=4))


def _operator_states(c, top, *, compact):
    """Restore the older locus, conserved-repeat, scan and ranked-site drawings."""
    centres = [52, 86, 120, 156]
    names = ["Homologue promoters", "Motif", "Genome rescan", "Ranked operators"]
    for x, title in zip(centres, names):
        c.text(x, top, title, size=FS_BOX, weight="bold", ha="center")
    for i in range(3):
        c.locus(43, top + 4.7 + i * 2, 18, color=RULE, promoter=True)
    # Aligned conserved half-sites across promoter sequences: no actual sequence is implied.
    for row in range(3):
        yy = top + 4.8 + row * 2
        c.line([(77, yy), (95, yy)])
        c.line([(79, yy), (83, yy)], OPER, lw=2)
        c.line([(89, yy), (93, yy)], OPER, lw=2)
    c.arrow([(79, top + 11), (83, top + 11)], color=OPER)
    c.arrow([(93, top + 11), (89, top + 11)], color=OPER)
    c.line([(111, top + 10.3), (129, top + 10.3)])
    c.line([(111, top + 6.8), (129, top + 6.8)], OPER, lw=0.4)
    for frac, h in [(0.08, 2), (0.24, 5), (0.42, 1.4), (0.59, 7), (0.76, 4.3), (0.92, 1)]:
        c.line([(111 + frac * 18, top + 10.3), (111 + frac * 18, top + 10.3 - h)],
               OPER if h >= 3.5 else RULE, lw=1.1)
    for rank, w in enumerate([14, 9, 6], 1):
        yy = top + 4.5 + (rank - 1) * 2.8
        c.text(144, yy, str(rank), color=OPER)
        c.line([(149, yy), (149 + w, yy)], OPER, lw=1.8)
    for x0, x1 in [(64, 73), (98, 107), (132, 140)]:
        c.arrow([(x0, top + 7), (x1, top + 7)], color=OPER)
    bodies = ["MMseqs2 orthologues", "inverted repeat → PWM",
              "default p ≤ 1e-4", "regulon + mode"]
    if compact:
        bodies = ["upstream DNA", "inverted repeat → PWM", "intergenic DNA", "regulon + mode"]
    for x, body in zip(centres, bodies):
        c.text(x, top + 15, body, color=MUTED, ha="center")


def draw_pipeline_figure(*, compact=False):
    """Both publication and graphical-abstract views use the same truthful stage graph."""
    c = _Canvas(130)
    c.text(2, 3.5, "RIIO: from protein or genome to operator, inducer and regulon",
           size=FS_TITLE, weight="bold")
    c.heading(10, "a", "Entry points", note="CLI · Python API · browser GUI")
    xs, width = [2, 46, 90, 134], 41.8
    centres = [x + width / 2 for x in xs]
    for i, (x, (title, body)) in enumerate(zip(xs, MODES)):
        c.box(x, 14, width, 16, fill=ENTRY_SOFT)
        if i < 2:
            c.proteins(x + 3, 20, n=1 if i == 0 else 3)
        else:
            c.locus(x + 3, 22.5, 10)
        c.text(x + 26, 18, title, size=FS_BOX, weight="bold", ha="center")
        c.text(x + 26, 24.5, body, ha="center", color=MUTED)

    # Supplied proteins bypass discovery; genome inputs retain the discovery/selection branch.
    c.text(90, 35, "b", size=FS_LETTER, weight="bold")
    c.text(95, 35, "Candidate discovery", size=FS_HEAD, weight="bold")
    c.box(90, 39, 85.8, 15, fill=ENTRY_SOFT, edge="none")
    c.proteins(94, 42, n=4)
    c.text(102, 50.5, "Proteome", size=FS_BOX, weight="bold", ha="center")
    c.arrow([(108, 47), (142, 47)])
    c.ssn(126, 43)
    c.text(126, 50.5, "SSN profiles", color=MUTED, ha="center")
    c.arrow([(132, 43), (142, 43)])
    c.text(144, 43, "BITACORA", size=FS_BOX, weight="bold")
    c.text(144, 48, "clade HMM + BLAST+")
    c.text(132.9, 58, "Candidate manifest + evidence → select TFs", ha="center", color=ENTRY)
    for cx in centres[2:]:
        c.line([(cx, 30), (cx, 32), (170, 32)], ENTRY)
    c.arrow([(170, 32), (170, 39)])
    c.arrow([(132.9, 54), (132.9, 55.8)])
    c.text(46, 43, "Supplied proteins skip discovery", size=FS_BOX, weight="bold", ha="center")
    c.text(46, 48, "all records submitted; no family filter", color=MUTED, ha="center")
    for cx in centres[:2]:
        c.line([(cx, 30), (cx, 32), (46, 32)], ENTRY)
    c.arrow([(46, 32), (46, 38)])

    # Only list/genome routes use the resumable runner. The single call goes directly to core.
    c.line([(centres[1], 51), (centres[1], 60), (132.9, 60), (132.9, 59.5)], ENTRY)
    c.arrow([(100, 60), (30, 60), (30, 67)])
    c.arrow([(centres[0], 51), (centres[0], 67)], color=ENTRY)
    c.text(5, 57, "single", color=MUTED)

    # Hero panel: actual biological states rather than a wall of command and file names.
    c.text(38, 62, "c", size=FS_LETTER, weight="bold")
    c.text(43, 62, "Per-TF inference", size=FS_HEAD, weight="bold")
    c.text(175.8, 62, "list/genome: resume / retry", color=MUTED, ha="right")
    c.box(2, 67, 31.5, 34, fill=ENTRY_SOFT)
    c.text(17.8, 70.5, "Resolve context", size=FS_BOX, weight="bold", ha="center")
    c.ssn(17.8, 75.5)
    c.text(17.8, 80, "family · SSN clade", ha="center")
    c.locus(6, 86, 23.5, promoter=True)
    c.text(17.8, 92, "locus + flanking DNA", ha="center")
    c.text(17.8, 97, "local or fetched", color=MUTED, ha="center")
    c.box(38, 67, 137.8, 19.5, fill=OPER_SOFT, edge="none")
    _operator_states(c, 69.5, compact=compact)
    c.arrow([(33.5, 76.5), (38, 76.5)], color=OPER)
    c.box(38, 89, 137.8, 12, fill=IND_SOFT, edge="none")
    channels = ["SSN clade\nlabel", "Cys / His\nmotif", "MetalNet2\nsite", "Nearby\ngenes", "Ligify\nchemistry"]
    for i, text in enumerate(channels):
        x = 39 + i * 20.4
        c.text(x + 9.3, 93, text, ha="center")
        c.line([(x + 9.3, 97), (x + 9.3, 99)], IND)
    c.line([(48.3, 99), (137, 99)], IND)
    c.arrow([(137, 99), (141, 99), (141, 94), (144, 94)], color=IND)
    c.text(159, 92, "Fused inducer call", size=FS_BOX, weight="bold", ha="center")
    c.text(159, 97, "class; ion if resolved", ha="center")
    c.arrow([(33.5, 94), (38, 94)], color=IND)
    c.text(2, 105, "Channels may abstain.", color=MUTED)
    c.text(175.8, 105, "Optional apo fold + QC (off by default).", color=MUTED, ha="right")

    c.heading(112, "d", "Outputs", note="reports rebuilt from saved data")
    groups = [("Per-TF bundle", "HTML + figures · dossier\nranked sites · run record"),
              ("Run report", "sortable HTML · S1 (XLSX / TSV)\ndiscovery evidence"),
              ("AlphaFold 3 hand-off", "apo / DNA / metal job files\nserver + local; submission only")]
    for x, (title, body) in zip([2, 61, 120], groups):
        c.box(x, 116, 55.8, 11, fill=OUT_SOFT, edge="none")
        c.text(x + 27.9, 118.5, title, size=FS_BOX, weight="bold", ha="center")
        c.text(x + 27.9, 123.5, body, color=MUTED, ha="center")
    c.arrow([(107, 101), (107, 110)], color=OUT)
    return c.fig


def main():
    try:
        from .figure_style import apply_style, save
    except ImportError:
        from figure_style import apply_style, save
    apply_style()
    save(draw_pipeline_figure(), "fig1_pipeline")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
