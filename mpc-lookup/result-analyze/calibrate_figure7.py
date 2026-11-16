#!/usr/bin/env python3
"""Recalibrate the ``\\resizebox`` fractions and the panel gaps of the paper's
8-panel figure.

Why this step exists
--------------------
The paper includes each panel with
``\\resizebox{ <frac> \\textwidth}{!}{\\input{images/panel_x.pgf}}``, i.e. the
whole panel is scaled by ``scale = frac * \\textwidth / pgf_width``. draw.py
authors the matplotlib canvas at final print size (the frame is exactly
FRAME_W x FRAME_H inches), so ``frac`` has to equal
``pgf_width / \\textwidth`` for the scale to be 1 and for the fonts to be
printed at their nominal point size.

On top of that the four panels have **different decoration widths** (the
dual-axis panels are wider, the single-curve panels narrower), so ``\\hfill``
would share the free space out between boxes of unequal width and the frames
of the two rows would not line up. The paper therefore uses an explicit
``\\hspace*{...}``, whose value this script solves backwards from the
**measured offset of each frame** inside its pgf, so that the frame left edges
of the two rows coincide column by column.

When it must be re-run
----------------------
Re-run the script after changing any of the following, otherwise the columns
drift apart again and the panels come out at the wrong size:

* FRAME_W / FRAME_H in draw.py
* the data range, axis labels, tick labels or legend contents of any panel
  (all of them change the pgf width)
* the order or the number of the subfloats

Usage
-----
    python3 calibrate_figure7.py                # calibrate and write back
    python3 calibrate_figure7.py --dry-run      # print only, change nothing
    python3 calibrate_figure7.py --tex PATH     # point at another .tex file

The script rewrites two sets of numbers only; the caption, the sub-captions
and the labels are never touched.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Constants shared with draw.py and with the paper's layout
# ---------------------------------------------------------------------------
PAPER_TEXTWIDTH_IN = 7.0      # \textwidth of the IEEEtran conference class
FRAME_W_IN = 1.08             # keep in step with FRAME_W in draw.py; used only
                              # to print the "right pad" column, not in the
                              # calibration itself
ROW1 = ['a', 'b', 'c', 'd']   # top row, left to right
ROW2 = ['f', 'g', 'e', 'h']   # bottom row, left to right
MIN_GAP_IN = 0.006            # smallest gap a panel may end up with
MAX_ROW_IN = 6.985            # row width limit, a little under \textwidth so
                              # that rounding cannot cause an overfull box

PIC_RE = re.compile(
    r'\\pgfpathrectangle\{\\pgfpointorigin\}\{\\pgfqpoint\{([0-9.]+)in\}\{([0-9.]+)in\}\}')
AXES_RE = re.compile(
    r'\\pgfpathrectangle\{\\pgfqpoint\{([0-9.]+)in\}\{([0-9.]+)in\}\}'
    r'\{\\pgfqpoint\{([0-9.]+)in\}\{([0-9.]+)in\}\}')


def read_geometry(pgf: Path) -> tuple[float, float]:
    """Return (canvas width, left offset of the frame) in inches."""
    text = pgf.read_text(errors='replace')
    m = PIC_RE.search(text)
    if not m:
        raise SystemExit(f"Cannot parse the canvas size from {pgf}")
    width = float(m.group(1))
    m2 = AXES_RE.search(text)
    if not m2:
        raise SystemExit(f"Cannot parse the frame position from {pgf}")
    return width, float(m2.group(1))


def frame_positions(names, L, gaps, lead=0.0, W=None):
    """Accumulate in order; return the frame left edges and the row width."""
    xs, x = [], lead
    for i, n in enumerate(names):
        xs.append(x + L[n])
        x += W[n]
        if i < len(gaps):
            x += gaps[i]
    return xs, x


def solve_gaps(W, L):
    """Find the gaps of both rows so that the columns line up.

    The top gaps g1..g3 and the bottom gaps h1..h3 are six free variables, and
    alignment gives four equations (column 1 matches by construction), so a
    solution space exists. It is searched on a grid, taking the solution whose
    *smallest* gap is as large as possible (the most even-looking row) subject
    to the row width staying within MAX_ROW_IN.
    """
    best = None
    step = 0.0005
    for g1 in _frange(0.0, 0.10, step):
        for g2 in _frange(0.0, 0.10 - g1, step):
            for g3 in _frange(0.0, 0.10 - g1 - g2, step):
                xs1, end1 = frame_positions(ROW1, L, [g1, g2, g3], 0.0, W)
                if end1 > MAX_ROW_IN:
                    continue
                lead = xs1[0] - L[ROW2[0]]
                if lead < -1e-9:
                    continue
                h1 = xs1[1] - (lead + W[ROW2[0]] + L[ROW2[1]])
                h2 = xs1[2] - (lead + W[ROW2[0]] + h1 + W[ROW2[1]] + L[ROW2[2]])
                h3 = xs1[3] - (lead + W[ROW2[0]] + h1 + W[ROW2[1]] + h2
                               + W[ROW2[2]] + L[ROW2[3]])
                hs = [h1, h2, h3]
                if min(hs) < MIN_GAP_IN:
                    continue
                end2 = lead + sum(W[n] for n in ROW2) + sum(hs)
                if end2 > MAX_ROW_IN:
                    continue
                score = min([g1, g2, g3] + hs)
                if best is None or score > best['score']:
                    best = dict(score=score, g=[g1, g2, g3], h=hs,
                                lead=lead, end1=end1, end2=end2, xs1=xs1)
    return best


def _frange(start, stop, step):
    v = start
    while v <= stop + 1e-12:
        yield v
        v += step


def rewrite_tex(tex: Path, fracs: dict, gaps: list, dry_run: bool) -> None:
    """Replace the \\resizebox fractions and the \\hspace* values only."""
    src = tex.read_text()
    start = src.index('\\begin{figure*}[!htp]')
    end = src.index('\\end{figure*}', start) + len('\\end{figure*}')
    block = src[start:end]

    # 1) the fraction of every panel
    for n, frac in fracs.items():
        pat = re.compile(
            r'(\\resizebox\{\s*)[0-9.]+(\s*\\textwidth\}\{!\}\{\\input\{images/panel_'
            + n + r'\.pgf\}\})')
        block, k = pat.subn(lambda m: f"{m.group(1)}{frac:.4f}{m.group(2)}", block, count=1)
        if k != 1:
            raise SystemExit(f"No resizebox for panel_{n} in the figure block")

    # 2) the six \hspace* values, in row-major order
    hs_pat = re.compile(r'(\\hspace\*\{)[0-9.]+(in\})')
    found = hs_pat.findall(block)
    if len(found) != len(gaps):
        raise SystemExit(
            f"The figure block has {len(found)} \\hspace*, expected {len(gaps)}. "
            "If the gaps were changed back to \\hfill, restore the explicit "
            "\\hspace* first.")
    it = iter(gaps)
    block = hs_pat.sub(lambda m: f"{m.group(1)}{next(it):.4f}{m.group(2)}", block)

    if dry_run:
        print("\n[dry-run] figure block not written back")
        return
    tex.write_text(src[:start] + block + src[end:])
    print(f"\nWrote back {tex}")


def main() -> int:
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--tex', type=Path,
                    default=here / '../../../mplookup-paper/mplookup.tex',
                    help='path to the paper .tex file')
    ap.add_argument('--dry-run', action='store_true')
    args = ap.parse_args()

    W, L = {}, {}
    for n in 'abcdefgh':
        pgf = here / f'panel_{n}.pgf'
        if not pgf.exists():
            raise SystemExit(f"Missing {pgf}; run draw.py in result-analyze/ first")
        W[n], L[n] = read_geometry(pgf)

    print(f"{'panel':>7}{'canvas':>9}{'L offset':>10}{'right pad':>10}{'frac':>9}")
    for n in 'abcdefgh':
        print(f"{n:>7}{W[n]:>8.3f}in{L[n]:>9.3f}in"
              f"{W[n] - L[n] - FRAME_W_IN:>9.3f}in{W[n] / PAPER_TEXTWIDTH_IN:>9.4f}")

    fracs = {n: W[n] / PAPER_TEXTWIDTH_IN for n in 'abcdefgh'}
    best = solve_gaps(W, L)
    if best is None:
        print("\nNo feasible solution: the row is wider than the limit, or a "
              "gap came out negative.")
        print("Shrink FRAME_W in draw.py (by about 0.03in) and regenerate the "
              "pgf files first.")
        return 1

    print(f"\nTop row gaps  : {[round(x, 4) for x in best['g']]}"
          f"  row width {best['end1']:.4f}in")
    print(f"Bottom row lead: {best['lead']:.4f}in"
          f"  gaps {[round(x, 4) for x in best['h']]}"
          f"  row width {best['end2']:.4f}in")
    print(f"Smallest gap  : {best['score']:.4f}in")
    print(f"Frame left edges : {[round(x, 4) for x in best['xs1']]}in")

    rewrite_tex(args.tex, fracs, best['g'] + best['h'], args.dry_run)
    print("\nNext steps:")
    print("  1) cp panel_*.pgf ../../../mplookup-paper/images/")
    print("  2) cd ../../../mplookup-paper && pdflatex mplookup.tex twice")
    print("  3) Render page 12 with pdftoppm and check that the frame left "
          "edges of each column agree")
    return 0


if __name__ == '__main__':
    sys.exit(main())
