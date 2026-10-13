"""Run-scoped F0/F1 schematic; ``main()`` builds Figure 1 with ``build_fig1_artwork``.

Contract: RIIO combines TF sequence/context with complementary DNA-operator and
inducer evidence to produce interpretable regulatory hypotheses and auditable
exports. A schematic-led composite: context, a dominant biological evidence panel,
then synthesis. Drawings are conceptual, with no measured sequences, structures,
scores or prediction-accuracy claim. Python only; editable SVG/TrueType PDF,
177.8 x 120 mm, 6–7 pt ordinary text and 8 pt panel letters; muted semantic hues.
"""
from __future__ import annotations

import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import (
    Circle,
    Ellipse,
    FancyArrowPatch,
    PathPatch,
    Polygon,
    Rectangle,
)
from matplotlib.path import Path

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42

INK, MUTED, RULE = "#252B2C", "#657071", "#B8C1C2"
ENTRY, ENTRY_SOFT = "#45515B", "#DCE7EB"
# Reuse the paper's operator identity; pastel fills keep it subordinate to ink.
OPER, OPER_SOFT = "#2E8581", "#ABD3D0"
IND, IND_SOFT = "#795565", "#E7D7DF"
FS_LETTER, FS_HEAD, FS_BODY, FS_SMALL = 8.0, 7.0, 6.5, 6.0
W_MM, H_MM = 177.8, 120.0


class _Canvas:
    """Millimetres from the top-left of the printed figure."""

    def __init__(self):
        self.fig, self.ax = plt.subplots(figsize=(W_MM / 25.4, H_MM / 25.4))
        self.fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        self.ax.set(xlim=(0, W_MM), ylim=(0, H_MM))
        self.ax.axis("off")

    def y(self, top):
        return H_MM - top

    def text(self, x, top, s, *, size=FS_BODY, weight="normal", color=INK, ha="left"):
        return self.ax.text(x, self.y(top), s, fontsize=size, fontweight=weight,
                            color=color, ha=ha, va="center", linespacing=1.12, zorder=5)

    def heading(self, top, letter, title):
        self.text(2, top, letter, size=FS_LETTER, weight="bold")
        self.text(7, top, title, size=FS_HEAD, weight="bold")

    def line(self, points, color=RULE, lw=0.65, *, dashed=False, z=2):
        xs, ys = zip(*[(x, self.y(top)) for x, top in points])
        self.ax.plot(xs, ys, color=color, lw=lw, zorder=z,
                     linestyle=(0, (2, 2)) if dashed else "-")

    def arrow(self, points, *, color=ENTRY, dashed=False, lw=0.7):
        if len(points) > 2:
            self.line(points[:-1], color, lw, dashed=dashed)
        self.ax.add_patch(FancyArrowPatch(
            (points[-2][0], self.y(points[-2][1])),
            (points[-1][0], self.y(points[-1][1])), arrowstyle="-|>",
            mutation_scale=6, color=color, linewidth=lw, shrinkA=0,
            shrinkB=0.45, linestyle=(0, (2, 2)) if dashed else "-", zorder=2))

    def rectangle(self, x, top, w, h, *, fill=OPER_SOFT, edge="none", lw=0.4):
        self.ax.add_patch(Rectangle((x, self.y(top + h)), w, h,
                                    facecolor=fill, edgecolor=edge, lw=lw, zorder=3))

    def group(self, x, top, w, h, *, accent=RULE):
        """Open editorial grouping: a scientific header rule, never a UI card."""
        self.line([(x, top), (x + min(w, 26), top)], accent, lw=0.8)
        if w > 26:
            self.line([(x + 26, top), (x + w, top)], RULE, lw=0.35)
        self.line([(x, top + h), (x + w, top + h)], RULE, lw=0.35)

    def dna(self, x, top, w, *, color=OPER, highlight=None, amplitude=0.7):
        """Generic duplex texture; no nucleotide identities or distances implied."""
        for sign in [-1, 1]:
            points = [(x + w * i / 70, top + sign * amplitude * math.sin(i / 70 * 4 * math.pi))
                      for i in range(71)]
            self.line(points, color, lw=0.5 if sign == -1 else 0.65)
        for i in range(15):
            xx = x + w * i / 14
            yy = amplitude * math.sin(i / 14 * 4 * math.pi)
            self.line([(xx, top - yy), (xx, top + yy)], color, lw=0.32)
        if highlight:
            start, width = highlight
            self.rectangle(x + start, top - 1.1, width, 2.2,
                           fill=OPER_SOFT, edge="none")

    def protein(self, x, top, *, scale=1.0, color=ENTRY, fill=ENTRY_SOFT,
                detail=True, labelled=False):
        """A conceptual lobed protein object, not a fold or oligomeric model."""
        vertices = [(-5.7, -0.4), (-6.8, -2.8), (-4.4, -4.6), (-2.1, -3.6),
                    (-0.5, -4.8), (3.3, -4.3), (4.3, -2.1),
                    (6.6, -2.1), (7.1, 0.4), (4.9, 1.8),
                    (5.8, 4.0), (1.5, 4.9), (0.0, 3.4),
                    (-2.6, 5.0), (-6.0, 3.1), (-5.7, -0.4)]
        codes = [Path.MOVETO] + [Path.CURVE4] * 15
        pts = [(x + scale * dx, self.y(top + scale * dy)) for dx, dy in vertices]
        self.ax.add_patch(PathPatch(Path(pts, codes), facecolor=fill, edgecolor=color,
                                   lw=0.7, zorder=4))
        if detail:
            ribbon = ([(-2.8, -1.2), (-0.9, -3.4), (1.2, -1), (3.0, -2)]
                      if labelled else
                      [(-3.7, 0.8), (-2, -2.5), (-0.6, 2.6), (0.9, -0.8),
                       (1.8, -2), (2.1, 2.1), (3.3, 1.1)])
            pts = [(x + scale * dx, self.y(top + scale * dy)) for dx, dy in ribbon]
            self.ax.add_patch(PathPatch(Path(pts, [Path.MOVETO] + [Path.CURVE4] * (len(pts) - 1)),
                                       facecolor="none", edgecolor=color, lw=0.45, zorder=5))

    def file(self, x, top, *, tag="FASTA", w=9, h=10, color=ENTRY, fill="white"):
        fold = min(1.8, w * 0.2)
        vertices = [(x, self.y(top)), (x + w - fold, self.y(top)),
                    (x + w, self.y(top + fold)), (x + w, self.y(top + h)),
                    (x, self.y(top + h))]
        self.ax.add_patch(Polygon(vertices, facecolor=fill, edgecolor=RULE,
                                  lw=0.55, zorder=3))
        self.line([(x + w - fold, top), (x + w - fold, top + fold),
                   (x + w, top + fold)], RULE, lw=0.45, z=4)
        self.text(x + w / 2 - 0.2, top + 1.7, tag, size=FS_SMALL, ha="center")
        self.line([(x + 1, top + 3.1), (x + w - 1, top + 3.1)], color, lw=0.4, z=4)
        for row, span in enumerate([0.75, 0.55, 0.8]):
            self.line([(x + 1, top + 4.2 + row * 1.5),
                       (x + 1 + span * (w - 2), top + 4.2 + row * 1.5)],
                      color, lw=0.45, z=4)

    def circle(self, x, top, radius, *, fill="white", edge=ENTRY, lw=0.5):
        self.ax.add_patch(Circle((x, self.y(top)), radius, facecolor=fill,
                                edgecolor=edge, lw=lw, zorder=4))

    def gene(self, x0, x1, top, *, reverse=False, color=RULE, h=1):
        head = min(1.3, (x1 - x0) * 0.35)
        h *= 0.76
        yy = self.y(top)
        if reverse:
            pts = [(x1, yy - h), (x0 + head, yy - h), (x0, yy),
                   (x0 + head, yy + h), (x1, yy + h)]
        else:
            pts = [(x0, yy - h), (x1 - head, yy - h), (x1, yy),
                   (x1 - head, yy + h), (x0, yy + h)]
        fill = "#F3F5F5" if color == RULE else color
        edge = "#859193" if color == RULE else color
        self.ax.add_patch(Polygon(pts, facecolor=fill, edgecolor=edge, lw=0.4, zorder=4))

    def locus(self, x, top, w, *, highlight=True):
        self.line([(x, top), (x + w, top)], lw=0.6)
        self.gene(x, x + w * 0.22, top)
        self.gene(x + w * 0.40, x + w * 0.64, top, reverse=True,
                  color=ENTRY_SOFT if not highlight else ENTRY)
        self.gene(x + w * 0.73, x + w, top)

    def proteins(self, x, top, *, n=1, color=ENTRY):
        for i, width in enumerate([9, 7, 8][:n]):
            # A sequence glyph, not an inferred structure.
            self.line([(x, top + i * 1.6), (x + width, top + i * 1.6)], color, lw=1.7)
            for frac in [0.15, 0.45, 0.75]:
                self.line([(x + frac * width, top + i * 1.6 - 0.5),
                           (x + frac * width, top + i * 1.6 + 0.5)], "white", lw=0.4)

    def ssn(self, x, top, *, color=ENTRY, small=False):
        scale = 0.75 if small else 1.0
        offsets = [(-3.2, -1.5), (-1.5, -3), (0.1, -0.9),
                   (2.8, -2), (3.1, 1.3), (0.7, 2.5), (-2.6, 1.2)]
        pts = [(x + scale * xx, top + scale * yy) for xx, yy in offsets]
        for i, j in [(0, 1), (0, 2), (0, 6), (1, 2), (2, 3), (2, 4),
                     (2, 5), (2, 6), (3, 4), (4, 5), (5, 6)]:
            self.line([pts[i], pts[j]], color if 2 in (i, j) else RULE, lw=0.4)
        for i, p in enumerate(pts):
            self.circle(*p, 0.5 if small else 0.6,
                        fill=ENTRY if i == 2 else (IND_SOFT if color == IND else ENTRY_SOFT),
                        edge=ENTRY if i == 2 else color, lw=0.4)


def _inputs(c):
    c.heading(3.5, "a", "From proteins or genomes")
    centers = [20, 61, 105, 151]
    for x, title in zip(centers, ["Single protein", "Protein list", "Genome accession", "Local genome"]):
        c.text(x, 9.5, title, size=FS_HEAD, ha="center")
    c.file(10, 13, tag="FASTA", w=10, h=9)
    c.protein(27, 17.5, scale=0.6)
    c.text(20, 26, "sequence / FASTA / protein ID", size=FS_SMALL, ha="center")
    c.file(50, 13, tag="FASTA", w=10, h=9)
    c.file(59, 14.5, tag="table", w=10, h=9)
    for xx in [61, 63.5, 66]:
        c.line([(xx, 18.5), (xx, 22.5)], RULE, lw=0.4, z=4)
    c.text(61, 26, "multi-FASTA / accession table", size=FS_SMALL, ha="center")
    c.ax.add_patch(Ellipse((102, c.y(17.7)), 11, 8, fill=False,
                          edgecolor=ENTRY, lw=0.75, zorder=4))
    c.dna(97.8, 17.7, 8.4, color=ENTRY, amplitude=0.4)
    c.file(111, 13, tag="ID", w=8, h=9)
    c.text(105, 26, "assembly accession", size=FS_SMALL, ha="center")
    c.file(139, 13, tag="FASTA", w=10, h=9)
    c.file(148, 14.5, tag="GFF", w=10, h=9)
    c.locus(160, 18, 10, highlight=False)
    c.text(151, 26, "FASTA + GFF / GenBank / proteome", size=FS_SMALL, ha="center")
    # Supplied and discovery-selected sequences merge before per-TF inference.
    for x in centers[:2]:
        c.line([(x, 29), (x, 32), (40.5, 32)], ENTRY, lw=0.65)
    c.line([(40.5, 32), (40.5, 36), (16.5, 36), (16.5, 44)], ENTRY, lw=0.7)
    c.text(24, 39.5, "supplied TFs", size=FS_SMALL, color=MUTED)
    for x in centers[2:]:
        c.line([(x, 29), (x, 32), (128, 32)], ENTRY, lw=0.65)
    c.arrow([(128, 32), (128, 34)], color=ENTRY)
    c.text(128, 36, "Candidate discovery", weight="bold", ha="center")
    c.text(128, 40, "BITACORA · curated SSN profiles", size=FS_SMALL, ha="center")
    c.line([(102, 37.8), (62, 37.8), (62, 44), (16.5, 44)], ENTRY, lw=0.65)
    c.text(78, 41.5, "selected TFs", size=FS_SMALL, color=MUTED, ha="center")
    c.arrow([(16.5, 44), (16.5, 51)], color=ENTRY)


def _context_hub(c):
    c.group(3, 52, 27, 62, accent=ENTRY)
    c.text(16.5, 55.5, "TF context", weight="bold", ha="center")
    c.protein(16.5, 64, scale=0.8)
    c.line([(8, 70), (25, 70)], ENTRY, lw=0.5)
    for i in range(12):
        c.circle(8 + i * 1.55, 70, 0.5, fill=ENTRY_SOFT if i % 3 else ENTRY,
                 edge=ENTRY, lw=0.35)
    c.text(16.5, 74, "Sequence", size=FS_SMALL, ha="center")
    c.ssn(16.5, 82, small=False)
    c.text(16.5, 88, "Family + SSN clade", size=FS_SMALL, ha="center")
    c.locus(7, 97.5, 19)
    c.dna(7, 101, 19, color=RULE, amplitude=0.5)
    c.text(16.5, 107, "Locus + flanking DNA", size=FS_SMALL, ha="center")
    c.text(16.5, 111, "when available", size=FS_SMALL, color=MUTED, ha="center")
    c.arrow([(30, 65), (37, 65)], color=OPER)
    c.arrow([(30, 91), (38, 91)], color=IND)


def _operator_evidence(c):
    c.group(37, 52, 101, 29, accent=OPER)
    c.text(40, 55, "Operator evidence", weight="bold")
    centers = [49, 74, 99, 124]
    names = ["Homologue\npromoters", "Inverted-repeat\nPWM", "Intergenic\nscan", "Ranked\noperators"]
    for x, name in zip(centers, names):
        c.text(x, 60.5, name, size=FS_BODY, ha="center")
    for row in range(3):
        top = 66 + row * 3
        c.line([(41, top), (57, top)], lw=0.5)
        c.rectangle(42 + row * 0.3, top - 0.7, 6.8, 1.4)
        c.gene(50, 57, top, color=RULE, h=0.8)
        c.arrow([(43, top + 1), (45, top + 1)], color=OPER, lw=0.4)
        c.arrow([(48, top + 1), (46, top + 1)], color=OPER, lw=0.4)
    c.text(49, 77, "upstream DNA", size=FS_SMALL, color=MUTED, ha="center")
    for row in range(3):
        top = 65.5 + row * 1.5
        c.line([(66, top), (82, top)], lw=0.4)
        c.rectangle(67, top - 0.45, 3, 0.9, fill=OPER)
        c.rectangle(78, top - 0.45, 3, 0.9, fill=OPER)
    c.arrow([(67, 63.7), (70, 63.7)], color=OPER, lw=0.4)
    c.arrow([(81, 63.7), (78, 63.7)], color=OPER, lw=0.4)
    for col in range(10):
        for row in range(4):
            chosen = row == [0, 1, 2, 3, 0, 3, 2, 1, 0, 3][col]
            c.rectangle(66.5 + col * 1.55, 70 + row * 1.05, 1.35, 0.8,
                        fill=OPER if chosen else "#D8EAE8")
    c.text(74, 77, "half-sites → motif", size=FS_SMALL, color=MUTED, ha="center")
    for row, offsets in enumerate([[2, 11], [6], [1, 12]]):
        top = 66 + row * 3
        c.dna(91, top, 16, color=RULE, amplitude=0.45)
        for offset in offsets:
            c.rectangle(91 + offset, top - 0.8, 2.2, 1.6,
                        fill=OPER if row < 2 else OPER_SOFT)
    c.text(99, 77, "genome-wide", size=FS_SMALL, color=MUTED, ha="center")
    for rank in range(3):
        top = 66 + rank * 3
        c.text(115.3, top, str(rank + 1), size=FS_SMALL, color=MUTED)
        c.line([(119, top), (135, top)], lw=0.5)
        c.rectangle(120 + rank * 0.9, top - 0.75, 2.5, 1.5,
                    fill=OPER if rank == 0 else OPER_SOFT)
        c.gene(127, 135 - rank * 0.6, top, color=RULE, h=0.8)
    c.text(124, 77, "regulon + mode", size=FS_SMALL, color=MUTED, ha="center")
    for x0, x1 in [(58.5, 63), (84, 88.5), (109, 113.5)]:
        c.arrow([(x0, 68.5), (x1, 68.5)], color=OPER, lw=0.6)
    c.arrow([(138, 68.5), (147, 68.5)], color=OPER)


def _inducer_evidence(c):
    c.group(37, 84, 101, 30, accent=IND)
    c.text(40, 87, "Inducer evidence", weight="bold")
    centers = [47, 65, 83, 102, 123]
    c.ssn(47, 93, color=IND, small=True)
    c.line([(60, 93), (70, 93)], IND, lw=0.7)
    for x, letter in [(61, "C"), (65, "H"), (69, "C")]:
        c.circle(x, 93, 1.4, fill=IND_SOFT, edge="none")
        c.text(x, 93, letter, size=FS_SMALL, ha="center")
    for row in range(3):
        c.line([(78, 91.5 + row * 1.4), (88, 91.5 + row * 1.4)], IND_SOFT, lw=0.8)
        for xx in [80, 86]:
            c.circle(xx, 91.5 + row * 1.4, 0.4, fill=IND, edge="none")
    c.ax.add_patch(FancyArrowPatch((80, c.y(91)), (86, c.y(91)), arrowstyle="-",
                                   connectionstyle="arc3,rad=-0.5", lw=0.5,
                                   color=IND, zorder=2))
    c.line([(97, 93), (108, 93)])
    c.gene(97, 100, 93, color=RULE, h=0.7)
    c.gene(102, 104, 93, color=IND_SOFT, reverse=True, h=0.7)
    c.gene(106, 109, 93, color=RULE, h=0.7)
    pts = [(123 + 2 * math.cos(k * math.pi / 3), 93 + 2 * math.sin(k * math.pi / 3)) for k in range(6)]
    c.line(pts + pts[:1], IND, lw=0.6)
    c.line([(125, 93), (128, 91.5)], IND, lw=0.6)
    c.circle(128, 91.5, 0.5, fill=IND_SOFT, edge=IND)
    # All five channels receive the shared context, rather than an SSN-only feed.
    c.line([(38, 91), (38, 89.5), (123, 89.5)], IND, lw=0.4)
    for x in centers:
        c.line([(x, 89.5), (x, 90.5)], IND, lw=0.4)
    names = ["SSN clade", "Cys/His\nmotifs", "MetalNet2", "Nearby\ngenes", "Ligify"]
    for x, name in zip(centers, names):
        c.text(x, 97.5 if x != 83 else 96.5, name, size=FS_SMALL, ha="center")
        c.line([(x, 100.5), (x, 102)], IND, lw=0.4)
    c.line([(47, 102), (123, 102)], IND, lw=0.4)
    c.text(83, 100, "sequence + MSA", size=FS_SMALL, color=MUTED, ha="center")
    c.arrow([(90, 102), (90, 103.5)], color=IND, lw=0.5)
    c.text(66, 106.7, "Fused inducer call", weight="bold", ha="center")
    for x, label in [(101, "class"), (117, "ion?")]:
        c.text(x, 106.7, label, size=FS_SMALL, ha="center")
    c.line([(90, 103.5), (117, 103.5)], IND, lw=0.45)
    c.line([(101, 103.5), (101, 104.2)], IND, lw=0.45)
    c.line([(117, 103.5), (117, 104.2)], IND, lw=0.45)
    c.text(87.5, 112, "ion if resolved · channels may abstain", size=FS_SMALL,
           color=MUTED, ha="center")
    # Both class and conditional ion are fields of the fused call, not serial
    # decisions. Their shared output feeds the hypothetical inducer object.
    c.line([(101, 109.3), (101, 110.2), (128, 110.2)], IND, lw=0.45)
    c.line([(117, 109.3), (117, 110.2)], IND, lw=0.45)
    c.arrow([(128, 110.2), (144.5, 110.2), (144.5, 89.3),
             (165, 89.3), (168, 87.7)], color=IND)


def _synthesis(c):
    c.text(143, 46.5, "c", size=FS_LETTER, weight="bold")
    c.text(148, 47.5, "Regulatory hypothesis", size=FS_BODY, weight="bold")
    c.line([(147, 52), (175, 52)], RULE, lw=0.35)
    c.text(161, 57.5, "Operators + regulon", size=FS_BODY, ha="center")
    c.dna(148, 68.5, 28, color=RULE, amplitude=1.15)
    c.rectangle(151, 67.2, 5.2, 2.6, fill=OPER_SOFT, edge=OPER, lw=0.45)
    c.gene(158, 165, 68.5, color=RULE, h=1.2)
    c.gene(167, 174, 68.5, color=RULE, h=1.2)
    c.line([(158, 71), (158, 72), (174, 72), (174, 71)], RULE, lw=0.4)
    c.text(166, 75, "candidate targets", size=FS_SMALL, color=MUTED, ha="center")
    c.protein(154, 82, scale=1.25, labelled=True)
    c.text(154, 82, "TF", size=FS_BODY, weight="bold", ha="center")
    c.arrow([(153, 77), (153, 71)], color=OPER, dashed=True, lw=0.65)
    c.circle(168, 86, 1.7, fill=IND_SOFT, edge=IND)
    c.circle(171.8, 87.5, 0.7, fill="white", edge=IND)
    c.arrow([(165.7, 85), (161, 81.7)], color=IND, dashed=True, lw=0.65)
    c.text(168, 91, "Inducer", size=FS_SMALL, ha="center")
    c.text(161, 95, "TF family + clade", size=FS_SMALL, ha="center")
    c.line([(147, 98), (175, 98)], RULE, lw=0.45)
    c.text(161, 101.5, "Auditable outputs", size=FS_BODY, weight="bold", ha="center")
    c.file(147, 105, tag="HTML", w=9, h=9)
    c.text(159, 106.5, "Per-TF bundle", size=FS_SMALL)
    c.text(159, 110.5, "Run / discovery", size=FS_SMALL)
    c.text(161, 116, "S1 · AF3 input jobs", size=FS_SMALL, ha="center")


def draw_pipeline_figure(*, compact=False):
    """One connected graph for publication and run-scoped schematic exports."""
    c = _Canvas()
    _inputs(c)
    c.text(35, 47, "b", size=FS_LETTER, weight="bold")
    c.text(40, 47, "TF context and complementary evidence", weight="bold")
    _context_hub(c)
    _operator_evidence(c)
    _inducer_evidence(c)
    _synthesis(c)
    return c.fig


def main():
    """Figure 1 is the hand-built SVG; this matplotlib drawing stays for run-scoped F0/F1."""
    try:
        from . import build_fig1_artwork
    except ImportError:
        import build_fig1_artwork
    return build_fig1_artwork.main()


if __name__ == "__main__":
    raise SystemExit(main())
