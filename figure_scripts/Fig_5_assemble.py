#!/usr/bin/env python3
"""
============================================================================
Fig_5_assemble.py -- assemble Fig. 5 by stacking the two row images.
  row 1 (A, B)  ggplot output of figure_scripts/Fig_5.R
  row 2 (C)     the evolutionary timeline, exported from BioRender
Row 2 is scaled to row 1's width and given a panel letter drawn to match the
letters already in row 1: the glyph box of "A" is measured directly off row 1,
and "C" is rendered at the same cap height, the same left edge and the same
distance below the top of its block, so the three letters line up by
construction rather than by hand-tuned constants.
in : ./figures/parts/Fig_5_row1.png
     ./figures/parts/Fig_5_row2_timeline.png
out: ./figures/main/Fig_5.png
RUN: python3 figure_scripts/Fig_5_assemble.py   (figure_scripts/Fig_5.R calls it)
============================================================================
"""
import os
import sys

from PIL import Image, ImageDraw, ImageFont

# the repository root: one level above this file
SRC = os.environ.get("SRC", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ROW1 = os.path.join(SRC, "figures/parts/Fig_5_row1.png")
ROW2 = os.path.join(SRC, "figures/parts/Fig_5_row2_timeline.png")
OUT = os.path.join(SRC, "figures/main/Fig_5.png")

# row 1 is written by png(res = 200) in lib/fig_imperfect_treg.R; the composite
# inherits that convention so the figure keeps one physical size
DPI = 200
# The panel letter is set in Arial Bold, the face of the letters in row 1. The
# first of these that exists is used (macOS, Windows, common Linux locations);
# set PANEL_FONT to the path of a bold sans-serif .ttf to override.
_FONTS = [
    "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
    "C:/Windows/Fonts/arialbd.ttf",
    "/usr/share/fonts/truetype/msttcorefonts/Arial_Bold.ttf",
    "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
    "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
]
FONT = os.environ.get(
    "PANEL_FONT", next((f for f in _FONTS if os.path.exists(f)), _FONTS[0])
)
LETTER = os.environ.get("LETTER", "C")
# Optical correction. At equal cap height a round letter reads smaller than one
# with wide flat extremities, and Arial's built-in overshoot (1px here) is tuned
# for running text rather than for an isolated label. 1.00 makes "C" metrically
# identical to "A"; 1.04 makes it look identical. Set LETTER_SCALE=1.0 to revert.
LETTER_SCALE = float(os.environ.get("LETTER_SCALE", "1.04"))


def flatten(path):
    """Open a PNG and composite any alpha onto white."""
    im = Image.open(path)
    if im.mode in ("RGBA", "LA") or (im.mode == "P" and "transparency" in im.info):
        bg = Image.new("RGBA", im.size, (255, 255, 255, 255))
        im = Image.alpha_composite(bg, im.convert("RGBA"))
    return im.convert("RGB")


def letter_box(im, search_w=0.30, search_h=0.22, thresh=40):
    """Bounding box of the panel letter alone in the top-left of row 1.

    Two things have to be separated out. The axis text is grey30 while the
    panel title is pure black, so a threshold of 40 keeps only the title. The
    title is then "A   P(rescue) ...", so its full bounding box would span the
    body text too; the letter is recovered as the first contiguous run of
    inked columns, and its height is measured within those columns only.
    """
    w, h = im.size
    g = im.crop((0, 0, int(w * search_w), int(h * search_h))).convert("L")
    mask = g.point(lambda v: 255 if v < thresh else 0)
    px = mask.load()
    cols = [any(px[x, y] for y in range(mask.height)) for x in range(mask.width)]
    if not any(cols):
        raise SystemExit("[Fig. 5] could not find a panel letter in row 1")
    x0 = cols.index(True)
    x1 = x0
    while x1 + 1 < len(cols) and cols[x1 + 1]:
        x1 += 1
    rows = [y for y in range(mask.height) if any(px[x, y] for x in range(x0, x1 + 1))]
    return x0, rows[0], x1 + 1, rows[-1] + 1


def fit_font(path, target_cap, ref_letter="A"):
    """Font size at which REF_LETTER renders target_cap tall.

    Calibration has to use the same glyph that was measured off row 1. Matching
    the bounding box of "C" to the bounding box of "A" would set "C" too small,
    because round letters are drawn with overshoot above the cap line and below
    the baseline, so their box is taller than their cap height at equal size.
    """
    size, best = 8, None
    while size < 400:
        f = ImageFont.truetype(path, size)
        bb = f.getbbox(ref_letter)
        if (bb[3] - bb[1]) > target_cap:
            break
        best, size = f, size + 1
    if best is None:
        raise SystemExit("[Fig. 5] could not fit the panel letter")
    return best


def ink_bearing(font, letter, pad=200):
    """Distance from the text origin to the first inked column of LETTER.

    PIL's font.getbbox() reports 0 here for both "A" and "C", so the only
    reliable way to align the ink rather than the origin is to render once on
    a scratch canvas and measure where the glyph actually starts.
    """
    im = Image.new("L", (pad * 2, pad * 2), 255)
    ImageDraw.Draw(im).text((pad, pad), letter, font=font, fill=0, anchor="ls")
    bb = im.point(lambda v: 255 if v < 40 else 0).getbbox()
    return 0 if bb is None else bb[0] - pad


def main():
    for p in (ROW1, ROW2):
        if not os.path.exists(p):
            raise SystemExit("[Fig. 5] missing input: %s" % p)

    r1 = flatten(ROW1)
    r2 = flatten(ROW2)
    w = r1.width

    # scale row 2 to row 1's width, preserving aspect
    r2 = r2.resize((w, round(r2.height * w / r2.width)), Image.LANCZOS)

    # geometry of the existing panel letter, reused for the new one
    lx, ly, rx, by = letter_box(r1)
    cap = by - ly
    header = ly + cap + round(cap * 0.6)  # white strip above row 2 holding "C"

    out = Image.new("RGB", (w, r1.height + header + r2.height), "white")
    out.paste(r1, (0, 0))
    out.paste(r2, (0, r1.height + header))

    font = fit_font(FONT, round(cap * LETTER_SCALE))
    d = ImageDraw.Draw(out)
    # "A" has no descender, so the bottom of its box in row 1 is the baseline.
    # Setting the new letter on that baseline, with anchor "ls" (left, baseline),
    # is what a typesetter would do and keeps any overshoot where it belongs.
    # anchor "ls" aligns text ORIGINS, but the letters have different left side
    # bearings (1px for "A", 3px for "C" in Arial Bold at this size), so the
    # origin is shifted by the bearing to make the ink flush down the figure.
    # font.getbbox() reports 0 for both, so the bearing is measured by drawing.
    d.text((lx - ink_bearing(font, LETTER), r1.height + by), LETTER,
           font=font, fill=(0, 0, 0), anchor="ls")

    out.save(OUT, dpi=(DPI, DPI))
    print(
        "[Fig. 5] wrote %s (%d x %d) | row1 %dpx + header %dpx + row2 %dpx | "
        "letter '%s' cap %dpx at x=%d" % (
            os.path.relpath(OUT, SRC), out.width, out.height,
            r1.height, header, r2.height, LETTER,
            round(cap * LETTER_SCALE), lx,
        )
    )


if __name__ == "__main__":
    sys.exit(main())
