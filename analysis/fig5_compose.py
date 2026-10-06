#!/usr/bin/env python
"""Compose the manuscript's four-genome figure (panels A-D) with our panel E into one 5-panel figure.

Panels A-D are taken from the copy embedded in `MainText__SSN_metalFamilies_v1.3.docx`
(word/media/image3.png, 796x977 px). That is the only copy available here, so **A-D in this composite
are limited to the resolution of the embedded raster** -- the output is a layout mock-up showing the
intended arrangement, not a print-ready master. To produce the final figure the authors should place
`F23_panelE_inducer_per_regulator.svg` (vector, editable text) beneath their own high-resolution A-D
source in Illustrator/Inkscape.

    py fig5_compose.py
"""
import sys
from pathlib import Path

from PIL import Image

# the report chain's output locations, so a figure is never drawn from an archived run
sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_paths as RP

HERE = Path(__file__).resolve().parent
PROJ = HERE.parent
#: The manuscript's own Figure-5 panels A-D. This used to point at
#: `%TEMP%/claude/ms_media/image3.png`, a scratch copy some earlier session unpacked -- so the figure
#: built once and then failed forever, because temp is swept. The .docx is a zip and the panel is
#: `word/media/image3.png` inside it, so extract it on demand from the manuscript we already have
#: and cache it in the run. Their artwork is never vendored into the repo.
MS_DIR = PROJ / "SSNs_metalloFamilies_2026"
MS_MEMBER = "word/media/image3.png"
AD = RP.RESULTS / "ms_figure5_panels_AD.png"
E = RP.FIGURES / "F23_panelE_inducer_per_regulator.png"
OUT = RP.FIGURES / "F24_fourgenome_ABCDE_composite.png"


def _extract_panels_ad() -> bool:
    """Pull panels A-D out of the newest manuscript .docx. False when it is not available."""
    import zipfile

    if AD.is_file() and AD.stat().st_size:
        return True
    docs = sorted(MS_DIR.glob("MainText__*.docx"))
    if not docs:
        print(f"no manuscript under {MS_DIR} -- cannot build the A-E composite")
        return False
    src = docs[-1]                      # newest version by name (v1.1 < v1.2 < v1.3)
    try:
        with zipfile.ZipFile(src) as z:
            if MS_MEMBER not in z.namelist():
                print(f"{src.name} has no {MS_MEMBER}")
                return False
            AD.parent.mkdir(parents=True, exist_ok=True)
            AD.write_bytes(z.read(MS_MEMBER))
    except Exception as exc:
        print(f"could not read {src.name}: {exc}")
        return False
    print(f"panels A-D extracted from {src.name} -> {AD.name}")
    return True

GAP_FRAC = 0.018          # vertical gap between the A-D block and panel E, as a fraction of width
PAD_FRAC = 0.012          # outer margin


def main():
    # The two missing-input cases here are NOT the same, and they must not be treated the same.
    #
    # Panels A-D come from the manuscript's own artwork, which this pipeline does not generate. A
    # run without the manuscript is a normal state, F24 is documented as non-primary (it embeds a
    # raster, so it is a layout mock-up rather than a submission master) and figure QA does not
    # check it. Failing here would stop the whole chain over a composite -- which is what happened
    # on the first production run. So: skip, say so plainly, and succeed.
    if not _extract_panels_ad():
        print("skipping F24 (non-primary): the manuscript artwork for panels A-D is not present.")
        print("  This is a normal state -- F24 is a layout mock-up, not a submission master.")
        return 0
    # Panel E is different. F23 is a PRIMARY figure produced by fig5_panelE.py, which runs
    # immediately before this script in the same stage. Its absence means that script failed, and
    # reporting success here would hide an upstream failure behind a skipped composite.
    if not E.is_file():
        print(f"cannot compose F24: missing {E}")
        print("  It is produced by: fig5_panelE.py, which runs immediately before this script.")
        print("  Panels A-D were found, so this is an upstream failure, not a missing manuscript.")
        return 1
    ad, e = Image.open(AD).convert("RGB"), Image.open(E).convert("RGB")

    # common width = panel E's native width (it is the higher-resolution of the two)
    W = e.width
    ad_h = round(ad.height * W / ad.width)
    ad_r = ad.resize((W, ad_h), Image.LANCZOS)

    gap, pad = round(W * GAP_FRAC), round(W * PAD_FRAC)
    H = pad + ad_h + gap + e.height + pad
    canvas = Image.new("RGB", (W + 2 * pad, H), "white")
    canvas.paste(ad_r, (pad, pad))
    canvas.paste(e, (pad, pad + ad_h + gap))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(OUT, dpi=(400, 400))
    print(f"wrote {OUT}  ({canvas.width}x{canvas.height})")
    print(f"  A-D  from {AD.name}: {ad.width}x{ad.height} -> upscaled to {W}x{ad_h}")
    print(f"       (upscale factor {W/ad.width:.1f}x -- mock-up only, see module docstring)")
    print(f"  E    {e.width}x{e.height} native")
    return 0


if __name__ == "__main__":
    sys.exit(main())
