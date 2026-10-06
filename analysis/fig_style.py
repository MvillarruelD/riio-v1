"""Shared publication-figure contract: restrained styling and deterministic exports.

Every primary figure is saved as editable SVG/PDF, a report-ready PNG and a 600-dpi LZW-compressed
TIFF.  Figure scripts should encode categories with text or geometry as well as colour.
"""
from __future__ import annotations
import os
import re
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- MANDATORY editable-SVG + font rules (first, always) ---
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams.update({
    "font.size": 7, "axes.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "figure.facecolor": "white", "savefig.facecolor": "white",
    # Nature-family conservative default for line/mixed artwork. Override only for quick iteration.
    "figure.dpi": 200, "savefig.dpi": float(os.environ.get("AN_DPI", 600)),
})

# Nature palette (references/api.md)
INK = "#272727"; NEUTRAL_D = "#4D4D4D"; NEUTRAL_M = "#767676"; NEUTRAL_L = "#CFCECE"; PANEL = "#F3F3F5"
# family identity — coherent blue/violet/teal (green & red reserved for directional verdicts)
FAM = {"ArsR": "#0F4D92", "ArsR/SmtB": "#0F4D92", "MerR": "#9A4D8E", "Fur": "#42949E"}
FAM_SOFT = {"ArsR": "#B4C0E4", "MerR": "#E7CBE2", "Fur": "#BFE0E2"}
# operator-hit verdict — directional
VERDICT = {"corroborated": "#2E9E44", "contradicted": "#E53935", "unsupported": "#E28E2C", "n/a": "#CFCECE"}
VERDICT_LABEL = {"corroborated": "corroborated", "contradicted": "contradicted",
                 "unsupported": "unsupported", "n/a": "non-metal"}
# DNA base colors for logos (classic, colour-blind-tolerant)
BASE = {"A": "#2E9E44", "C": "#0F4D92", "G": "#E8A33D", "T": "#E53935"}


#: Nature-family double-column width. The single number every primary figure is held to.
MAX_WIDTH_MM = 183.0
MAX_WIDTH_IN = MAX_WIDTH_MM / 25.4
PAD_INCHES = 0.04

# --------------------------------------------------------------------------- the width measurement
#: ONE definition of "how wide is this figure", used by `save()` here and by `figure_qa` as its gate.
#: They used to measure differently -- `save()` via `get_tightbbox` on the Agg renderer, `figure_qa`
#: via the SVG's own attributes -- and disagreed: on 2026-09-04 `save()` passed F4b silently while
#: the QA gate measured 184.8 mm and failed it. Agg and SVG lay text out with different metrics, so
#: a warning computed from the renderer is not a prediction of what lands on disk. Measuring the
#: written file makes the early warning and the gate the same check by construction.
_LENGTH = re.compile(r"^([0-9.]+)\s*([a-z%]*)$")
_UNIT_TO_MM = {"mm": 1.0, "cm": 10.0, "in": 25.4, "pt": 25.4 / 72, "pc": 25.4 / 6,
               "px": 25.4 / 96, "": 25.4 / 96}


def svg_size_mm(raw: str) -> tuple[float | None, float | None]:
    """(width, height) in millimetres from an SVG's own attributes, or (None, None)."""
    out: list[float | None] = []
    for name in ("width", "height"):
        match = re.search(rf'<svg[^>]*?\b{name}\s*=\s*"([^"]+)"', raw[:4000])
        if not match:
            out.append(None)
            continue
        parsed = _LENGTH.match(match.group(1).strip())
        if not parsed or parsed.group(2) not in _UNIT_TO_MM:
            out.append(None)
            continue
        out.append(float(parsed.group(1)) * _UNIT_TO_MM[parsed.group(2)])
    return out[0], out[1]


def final_width_in(fig) -> float:
    """The width the saved file will actually have, in inches.

    NOT `figsize`. Exports use `bbox_inches="tight"`, which expands the saved area to include
    anything drawn OUTSIDE the canvas -- a `suptitle` anchored past the axes, a legend placed with
    `bbox_to_anchor`, an over-long tick label. A figure declared at 7.2 in can therefore land at
    8.1 in on disk, and the overflow is invisible until someone measures the file. This measures it
    before the file is written.
    """
    fig.canvas.draw()
    bbox = fig.get_tightbbox(fig.canvas.get_renderer())
    return float(bbox.width) + 2 * PAD_INCHES


def save(fig, path_noext, formats=("svg", "pdf", "png", "tiff"), *, supplementary=False):
    """Save one figure in the project-wide publication formats and close it.

    SVG keeps text editable; PDF embeds TrueType text; PNG is used by the generated report; TIFF is
    the submission raster. The caller supplies a path without an extension.

    The width of the WRITTEN SVG is measured and reported, using the same `svg_size_mm` that
    `figure_qa.py` gates on -- so this warning and that gate are the same check by construction. The
    fix for an over-width figure is a layout change (bring the suptitle or legend inside the axes,
    shorten the labels, anchor an annotation `ha="right"`), never a silent rescale: rescaling a
    figure to fit also rescales its text, and 7 pt does not survive much of that.

    Pass `supplementary=True` for a figure the contract documents as full-width.
    """
    p = Path(path_noext)
    p.parent.mkdir(parents=True, exist_ok=True)
    canvas_mm = fig.get_size_inches()[0] * 25.4
    out = []
    for f in formats:
        fp = f"{p}.{f}"
        kwargs = {"bbox_inches": "tight", "pad_inches": PAD_INCHES}
        if f == "tiff":
            kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(fp, **kwargs)
        out.append(fp)
    plt.close(fig)

    # Measured AFTER writing, from the file itself. Measuring the figure object instead -- via
    # get_tightbbox on the Agg renderer -- is what let F4b through at 184.8 mm on 2026-09-04 while
    # the QA gate failed it: Agg and the SVG backend lay text out with different metrics, so the
    # renderer's answer is not a prediction of the file's. If the SVG was not requested there is
    # nothing to measure and the gate remains the only check.
    if not supplementary and "svg" in formats:
        try:
            width_mm, _ = svg_size_mm(Path(f"{p}.svg").read_text(encoding="utf-8", errors="replace"))
        except OSError:
            width_mm = None
        if width_mm and width_mm > MAX_WIDTH_MM + 0.5:
            print(f"    ! {p.name}: written width {width_mm:.1f} mm exceeds the "
                  f"{MAX_WIDTH_MM:.0f} mm publication limit "
                  f"(canvas is {canvas_mm:.1f} mm; the rest is content drawn outside it). "
                  f"Figure QA will fail until the layout is brought inside.")
    return out


# --------------------------------------------------------------------------- fail-closed inputs
def require(ok, what: str, produced_by: str, hint: str = "") -> None:
    """Stop with exit code 1 unless `ok`, naming the missing input and the stage that writes it.

    Why this exists rather than an early `return`. A figure script that prints "skipping" and exits
    0 tells the driver the stage succeeded. The report is then built from a figure set that is
    quietly one panel short -- and on a dirty output directory it is worse than that, because the
    missing panel is filled by an OLDER run's file and the report looks complete and is wrong.

    "Not enough data yet" is not a reason to succeed. The production driver refuses to run the
    figure stage until every bundle is complete, so by the time a figure script runs, missing input
    means a real defect upstream, not an interrupted batch.
    """
    if ok:
        return
    message = [f"cannot draw this figure: {what} is missing or empty.",
               f"  It is produced by: {produced_by}"]
    if hint:
        message.append(f"  {hint}")
    raise SystemExit("\n".join(message))


# --------------------------------------------------------------------------- family identity
#: The three families the focus figures (F3, F4a/F4b/F4c) are drawn for.
FOCUS_FAMILIES = ("ArsR", "MerR", "Fur")


def focus_family(predicted: str | None) -> str | None:
    """Map a PREDICTED family to its focus-family label, or None if it is not one of them.

    Take the family from the bundle's `dossier.json`, never from the run name. The run name encodes
    the family DISCOVERY assigned -- `MtubH37Rv__ArsR__NP_216235.1` -- and discovery matches a seed
    profile against part of a protein, so it is a hypothesis. The pipeline then classifies the full
    sequence and can disagree: measured on this run, 5 of 82 candidates discovered as ArsR were
    classified as MarR/SlyA, HxlR or IclR-type, and all of them declined to name a metal.

    The figures used to read the run-name string, so those five were drawn into the ArsR panels
    while `operator_predictions.csv` -- which takes `dossier["family"]` -- listed them as something
    else. Same run, two different answers, and the figure was the one that looked authoritative.

    Matching is on the head of the label before any "/", so "ArsR/SmtB" resolves to ArsR while
    "MarR/SlyA" correctly does not.
    """
    if not predicted:
        return None
    head = predicted.split("/")[0].strip().casefold()
    for fam in FOCUS_FAMILIES:
        if head == fam.casefold():
            return fam
    return None
