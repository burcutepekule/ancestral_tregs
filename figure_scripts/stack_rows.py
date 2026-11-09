#!/usr/bin/env python3
"""
============================================================================
stack_rows.py -- stack two separately rendered halves into one figure.

The same job Fig_5_assemble.py does, minus the panel-letter typesetting: here
both halves already carry their own letters (A, B from lib/fig_skeleton.R,
C, D from lib/fig_load_balance.R), so nothing is drawn on top. The bottom
half is resized to the top half's pixel width if the two differ, and the two
are butted together with no gap on a white canvas.

The callers render both halves at the SAME pixel width and scale each half's
device `res` in proportion to that width, so a resize here is a no-op and the
two halves keep the relative font sizes they have in the published figures.

in : argv[1] top half, argv[2] bottom half
out: argv[3]
env: DPI  the resolution written into the PNG header (default 250, the res of
          the phase-portrait half)
RUN: python3 figure_scripts/stack_rows.py top.png bottom.png out.png
============================================================================
"""
import os
import sys

from PIL import Image

DPI = float(os.environ.get("DPI", "250"))


def flatten(path):
    """Open a PNG and composite any alpha onto white."""
    im = Image.open(path)
    if im.mode in ("RGBA", "LA") or (im.mode == "P" and "transparency" in im.info):
        bg = Image.new("RGBA", im.size, (255, 255, 255, 255))
        im = Image.alpha_composite(bg, im.convert("RGBA"))
    return im.convert("RGB")


def main(argv):
    if len(argv) != 4:
        raise SystemExit("usage: stack_rows.py TOP.png BOTTOM.png OUT.png")
    top, bottom, out = argv[1], argv[2], argv[3]
    for p in (top, bottom):
        if not os.path.exists(p):
            raise SystemExit("[stack] missing input: %s" % p)

    t = flatten(top)
    b = flatten(bottom)
    if b.width != t.width:
        b = b.resize((t.width, round(b.height * t.width / b.width)), Image.LANCZOS)

    im = Image.new("RGB", (t.width, t.height + b.height), "white")
    im.paste(t, (0, 0))
    im.paste(b, (0, t.height))

    d = os.path.dirname(os.path.abspath(out))
    os.makedirs(d, exist_ok=True)
    im.save(out, dpi=(DPI, DPI))
    print("[stack] wrote %s (%d x %d) | top %dpx + bottom %dpx"
          % (out, im.width, im.height, t.height, b.height))


if __name__ == "__main__":
    sys.exit(main(sys.argv))
