"""figure_style.py -- one style, one palette, one saver for every figure in reports/.

Built on cnsplots (Nature/Cell-style defaults: editable SVG/PDF text, open axes, 7-8 pt type),
with matplotlib GridSpec for composition. Sizes follow Nature: 183 mm double column, 89 mm single.

Outcome colours are Okabe-Ito (colour-blind safe) and every outcome also carries a glyph, so no
distinction rests on hue alone.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import cnsplots as cns  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
FIGDIR = HERE / "figures"
MM = 1 / 25.4
W2 = 183 * MM          # double column, inches
W1 = 89 * MM           # single column, inches
MAX_WIDTH_MM = 183.0   # figure_qa rejects anything wider

C = {
    # signal family: metal calls / metal sensors
    "metal": "#0F4D92",
    "metal_mid": "#3775BA",
    "metal_soft": "#B4C0E4",
    # the two genomes of the benchmark
    "ecoli": "#0F4D92",
    "salm": "#42949E",
    # inducer classes
    "redox": "#42949E",
    "organic": "#B9A17B",
    "undetermined": "#D8D8D8",
    # outcomes (Okabe-Ito): correct / partial / incorrect / not scored
    "correct": "#0072B2",
    "partial": "#E69F00",
    "incorrect": "#D55E00",
    "not_scored": "#E6E6E6",
    # neutrals
    "neutral_light": "#D8D8D8",
    "neutral_mid": "#A8A8A8",
    "neutral_dark": "#606060",
    "ink": "#272727",
    "chance": "#9E9E9E",
    # the operator path (Figure 1's lane colour), used wherever an operator result is drawn
    "operator": "#2E8581",
    # families outside the twelve studied ones
    "outside": "#B9A17B",
    "outside_soft": "#EFE7D8",
    "empty": "#F6F6F6",
}
CLASS_COLOR = {"metal": C["metal"], "redox": C["redox"], "organic": C["organic"],
               "undetermined": C["undetermined"]}
CLASS_LABEL = {"metal": "metal", "redox": "redox", "organic": "organic",
               "undetermined": "undetermined"}
#: One colour per ion, fixed across every figure so an ion reads the same everywhere. An ion not
#: listed here falls back to ION_FALLBACK rather than to a colour that means something else.
ION_COLOR = {
    "Zn2+": "#56B4E9",
    "Cu+": "#C96144",
    "Cd2+": "#8281B9",
    "Cd2+/Pb2+": "#4B4A8C",
    "As(III)": "#009E73",
    "Fe2+": "#9C5732",
    "Mn2+": "#CC79A7",
    "Ni2+": "#E69F00",
    "Co2+": "#F0C75E",
    "MoO42-": "#7A7A3A",
}
ION_FALLBACK = "#9AA3AF"
#: display labels (sub/superscripts as plain text keeps SVG text editable)
ION_LABEL = {"MoO42-": "MoO4 2-", "Cd2+/Pb2+": "Cd2+/Pb2+"}
#: The genomes, in the order every figure uses; benchmark genomes first.
BENCH_GENOMES = [("ecoli", "E. coli K-12"), ("salmonella", "S. Typhimurium SL1344")]
SURVEY_ORDER = ["M. tuberculosis", "M. avium", "V. cholerae", "V. vulnificus"]
SHORT = {"M. tuberculosis": "M. tb", "M. avium": "M. avium", "V. cholerae": "V. cho",
         "V. vulnificus": "V. vul", "E. coli K-12": "E. coli", "S. Typhimurium SL1344": "SL1344"}


def is_species(label: str) -> bool:
    """Italicise a genome label that is a binomial (not a strain code such as SL1344)."""
    return label.startswith(("M.", "V.", "E.", "S. "))


def ion_color(ion: str) -> str:
    return ION_COLOR.get(ion, ION_FALLBACK)


def ion_label(ion: str) -> str:
    return ION_LABEL.get(ion, ion)
OUTCOME_GLYPH = {"correct": "✓", "partial": "~", "incorrect": "✗", "not_scored": ""}


def apply_style() -> None:
    cns.settings.savefig_transparent = False
    cns.settings.savefig_dpi = 600
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "Liberation Sans"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        # the type scale (FS_* below): body text for axes, ticks and legends; headings for titles
        "font.size": 5.6,
        "axes.labelsize": 5.6,
        "axes.titlesize": 6.4,
        "legend.fontsize": 5.6,
        "legend.frameon": False,
        "xtick.labelsize": 5.6,
        "ytick.labelsize": 5.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.6,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "axes.grid": False,
    })


def numbers() -> dict:
    return json.loads((HERE / "numbers.json").read_text(encoding="utf-8"))


#: One type scale for every figure (Nature: 5-7 pt text, 8 pt bold panel letters). Figure 1 uses the
#: same values under its own names.
FS_LETTER = 8.0     # panel letters
FS_HEAD = 6.4       # panel titles and other headings
FS_BODY = 5.6       # axis labels, tick labels, legends, table text
FS_SMALL = 5.0      # annotations inside a panel: n/N labels, notes, cell text


def type_scale() -> dict:
    """rcParams that put axis text on the type scale; use as `with plt.rc_context(type_scale()):`."""
    return {"font.size": FS_BODY, "axes.labelsize": FS_BODY, "axes.titlesize": FS_HEAD,
            "xtick.labelsize": FS_BODY, "ytick.labelsize": FS_BODY, "legend.fontsize": FS_BODY}


def panel_head(ax, letter: str, title: str = "", *, pad_pt: float = 5.0, letter_dx_pt: float = -12.0) -> None:
    """A panel letter at a fixed offset from the axes' top-left corner, in points, and the title
    right after it, so every panel reads letter-then-title whatever each axes' size."""
    from matplotlib.transforms import ScaledTranslation
    fig = ax.figure

    def at(dx, dy):
        return ax.transAxes + ScaledTranslation(dx / 72, dy / 72, fig.dpi_scale_trans)

    ax.text(0, 1, letter, transform=at(letter_dx_pt, pad_pt), fontsize=FS_LETTER, fontweight="bold",
            ha="left", va="bottom", color=C["ink"])
    if title:
        ax.text(0, 1, title, transform=at(letter_dx_pt + 10.0, pad_pt), fontsize=FS_HEAD, fontweight="bold",
                ha="left",
                va="bottom", color=C["ink"])


def panel(ax, label: str, x: float = -0.08, y: float = 1.03) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="bottom", ha="left", color=C["ink"])


def _retry(write, attempts: int = 5, wait_s: float = 2.0) -> None:
    """Run one file write, retrying briefly on OSError.

    reports/ lives in Dropbox, which locks a file while it syncs it; overwriting a figure it is still
    uploading fails with "[Errno 22] Invalid argument". The lock clears within seconds.
    """
    import time
    for i in range(attempts):
        try:
            write()
            return
        except OSError:
            if i == attempts - 1:
                raise
            time.sleep(wait_s)


def save(fig, name: str) -> None:
    """SVG + PDF (editable text) + PNG (600 dpi) through cnsplots; TIFF at 600 dpi, LZW."""
    FIGDIR.mkdir(parents=True, exist_ok=True)
    plt.figure(fig.number)
    for ext in ("svg", "pdf", "png"):
        _retry(lambda: cns.savefig(str(FIGDIR / f"{name}.{ext}")))
    _retry(lambda: fig.savefig(FIGDIR / f"{name}.tiff", dpi=600, bbox_inches="tight",
                               pil_kwargs={"compression": "tiff_lzw"}))
    plt.close(fig)
    # The exports are cropped to their content, so a long line of text widens the page silently.
    # Measure the file, as figure_qa does, and say so here rather than at the QA gate.
    from PIL import Image
    with Image.open(FIGDIR / f"{name}.png") as im:
        w_mm = im.size[0] / (im.info.get("dpi", (600, 600))[0]) * 25.4
    flag = "" if w_mm <= MAX_WIDTH_MM else f"  !! wider than {MAX_WIDTH_MM:.0f} mm"
    print(f"  wrote figures/{name}.[svg|pdf|png|tiff]  ({w_mm:.0f} mm wide){flag}")


def wrap(text: str, width_in: float, fontsize: float) -> str:
    """Hard-wrap a line of figure text to fit `width_in` inches at `fontsize` pt (Arial average)."""
    import textwrap
    return textwrap.fill(text, width=max(20, int(width_in * 72 / (fontsize * 0.5))))


def frac(text: str) -> tuple[int, int]:
    """'7/17' -> (7, 17); tolerates trailing commentary."""
    a, b = text.split()[0].split("/")
    return int(a), int(b)
