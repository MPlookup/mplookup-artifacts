"""
Experiment result analysis for Oblivious Lookup Argument.
Parses result files from mpc-lookup/results and generates figures.

Algorithm names:
  MPlookup  – the new O(N log²N) MPC algorithm
  Strawman  – the naive O(N²) baseline

Figure style: "slender" academic style (see style_guide/README.md).
Thin strokes (~0.9 pt), small markers (~3.6 pt), tab10 palette, compact
legends with light frames, faint grids, and high text-to-plot ratio.
"""

import re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import font_manager as _fm

from config import SAVE_FIG_FORMAT

# ---------------------------------------------------------------------------
# Slender academic style — single source of truth for all visual parameters
# ---------------------------------------------------------------------------

FONT_BASE    = 9.5      # axis labels / base text
FONT_TICK    = 8.5      # tick labels
FONT_LEGEND  = 8.5      # legend text (README rcParams)
FONT_STEPLG  = 8.5      # legend text for the dense 9-series panel

LW           = 0.9      # data line width
LW_REF       = 0.7      # reference / guide line width
MS           = 4.0      # marker size (README rcParams)
MSP          = 3.4      # marker size for dense multi-series panels
MARKER_EW    = 0.6      # marker edge width
SPINE_LW     = 0.6      # axes frame
TICK_LW      = 0.5      # tick marks
TICK_LEN     = 2.6      # tick length

GRID_COLOR   = '0.88'
GRID_LW      = 0.4

# ---------------------------------------------------------------------------
# Panel canvas — designed AT FINAL PRINT SIZE.
#
# The paper includes each panel with
#     \resizebox{<FRAC>\textwidth}{!}{\input{images/panel_x.pgf}}
# so the on-page width is FRAC * \textwidth and the whole picture is scaled by
#     scale = FRAC * \textwidth / pgf_width.
# If the matplotlib canvas is authored at exactly the width it will occupy on
# the page, then scale == 1 and every font is printed at its nominal point
# size. That is why the canvas below is sized from PAPER_TEXTWIDTH and the FRAC
# constants — not picked by eye.
#
#   panel_a/b/c/f/g/h : FRAC 0.28 -> 1.96 in on the page
#   panel_d/e         : FRAC 0.34 -> 2.38 in on the page
#
# bbox_inches='tight' + pad_inches=0.02 makes the exported pgf a little wider
# than figsize (ink such as tick labels sits outside the figure rectangle), so
# the canvas is authored _INK_OVERHANG narrower than the target print width;
# the constant is calibrated so the exported pgf width lands on the target width
# (measured: scale ends up ~1.02, i.e. printed sizes are 9.7/8.7 pt).
# ---------------------------------------------------------------------------
PAPER_TEXTWIDTH_IN = 7.0     # IEEEtran conference \textwidth
_INK_OVERHANG_IN   = 0.047   # pgf ends up this much wider than figsize

FRAC_STD_PANEL  = 0.28
FRAC_WIDE_PANEL = 0.34

# --- One frame geometry for every panel --------------------------------------
# The frame (the rectangle enclosed by the four spines) must come out the same
# size on paper in every panel. tight_layout cannot arrange that: the dual-axis
# panels carry a wider right-hand decoration, so their frame gets squeezed. The
# physical size of the frame is therefore set here directly, the canvas is given
# enough room on all four sides for ticks, axis labels and the legend, and
# bbox_inches='tight' trims the white left over. The exported pgf is then
# (frame width + this panel's own decoration width) wide, which differs per
# panel, and each panel's \resizebox fraction is calibrated from it so that the
# scale is always 1 and the frame lands on paper at exactly its nominal size.
FRAME_W = 1.06
FRAME_H = 0.985
_MARGIN_L = 0.70
_MARGIN_R = 0.62
_MARGIN_B = 0.95
_MARGIN_T = 0.16

PANEL_ASPECT = FRAME_W / FRAME_H

# tab10 palette (identical hues to the approved reference figures)
C_BLUE   = '#1f77b4'
C_ORANGE = '#ff7f0e'
C_GREEN  = '#2ca02c'
C_RED    = '#d62728'
C_PURPLE = '#9467bd'
C_BROWN  = '#8c564b'
C_PINK   = '#e377c2'
C_GRAY   = '#7f7f7f'
C_CYAN   = '#17becf'

# Per-step colours for the 9-step breakdown panel (a), (b) … (h)
STEP_PALETTE = [C_BLUE, C_ORANGE, C_GREEN, C_RED, C_PURPLE,
                C_BROWN, C_PINK, C_GRAY, C_CYAN]

# ---------------------------------------------------------------------------
# Quantity encoding — shared by the five "physical quantity" panels only,
# i.e. (a) time, (b) bytes+broadcasts, (c) speedup, (e) time, (f) bytes+
# broadcasts. The component panels (d), (g), (h) deliberately keep their own
# per-phase / per-step scheme (see STEP_PALETTE) and must not reuse these.
#
#   colour + marker  = which quantity   (time / bytes sent / broadcasts)
#                       — for *every* series, the strawman included
#   fill  + linestyle = which scheme     (MPlookup solid+filled,
#                                         Strawman dashed+hollow)
#
# "hollow + dashed" is the *exclusive* signature of the strawman baseline: no
# other series may use a dashed line or a hollow face, everything else is
# filled + solid. Without that rule the red hollow dashed triangle in (f)
# means "MPlookup's broadcasts" while the identical mark in (b) means "the
# strawman's broadcasts".
#
# The strawman keeps the colour of its quantity on purpose: in (b) the two
# baselines (bytes and broadcasts) land in the same region of the panel, and
# painting both grey made them impossible to tell apart.
#
# Time is orange rather than blue because time is the one quantity that never
# shares a panel with another: handing it the spare hue leaves blue and red for
# bytes vs broadcasts in (b)/(f), which is the comparison the eye has to make
# there. As a side effect no (colour, marker) pair here collides with the
# component panels' blue-circle / orange-square / green-triangle / red-diamond
# set.
Q_TIME  = dict(color=C_ORANGE, marker='o')  # time (incl. MPlookup total)
Q_BYTES = dict(color=C_BLUE,   marker='s')  # bytes sent
Q_BCAST = dict(color=C_RED,    marker='^')  # broadcasts

# Baseline / strawman series are drawn dashed with hollow markers so that the
# two curves can be told apart in greyscale print as well as in colour.
# Only the strawman may use this style — see the encoding rule above.
BASE_STYLE = dict(linestyle='--', markerfacecolor='white')

_BASE_RC = {
    'font.size': FONT_BASE,
    'axes.labelsize': FONT_BASE,
    'xtick.labelsize': FONT_TICK,
    'ytick.labelsize': FONT_TICK,
    'legend.fontsize': FONT_LEGEND,
    'axes.linewidth': SPINE_LW,
    'xtick.major.width': TICK_LW,
    'ytick.major.width': TICK_LW,
    'xtick.major.size': TICK_LEN,
    'ytick.major.size': TICK_LEN,
    'xtick.major.pad': 1.0,
    'ytick.major.pad': 1.2,
    'lines.linewidth': LW,
    'lines.markersize': MS,
    'lines.markeredgewidth': MARKER_EW,
    'patch.linewidth': 0.45,
    'figure.dpi': 160,
    'savefig.dpi': 300,
    'font.family': 'serif',
    'mathtext.fontset': 'cm',
    'axes.axisbelow': True,
}

if SAVE_FIG_FORMAT == "pgf":
    # https://matplotlib.org/stable/tutorials/text/pgf.html
    matplotlib.use("pgf")
    plt.rcParams.update(_BASE_RC)
    plt.rcParams.update(
        {
            "text.usetex": True,
            "pgf.texsystem": "pdflatex",
            # No font package here: with pgf.rcfonts = False the PGF output
            # inherits the fonts of the including document (IEEEtran/Times),
            # which keeps the figures serif and consistent with the paper.
            # Keep the preamble free of font packages (cmbright is not
            # installed and is sans-serif anyway).
            "pgf.preamble": "\n".join(
                [
                    r"\usepackage[utf8x]{inputenc}",
                    r"\usepackage[T1]{fontenc}",
                ]
            ),
            "pgf.rcfonts": False,
            "font.serif": [],  # use latex default
            "font.sans-serif": [],  # use latex default
            "font.monospace": [],  # use latex default
            "mathtext.fontset": "cm",
        }
    )
else:
    # Register Nimbus Roman (metrically identical to Times New Roman) so that
    # matplotlib can fall back to it when Times New Roman is not installed.
    import glob as _glob
    for _f in _glob.glob('/usr/share/fonts/opentype/urw-base35/NimbusRoman*.otf'):
        try:
            _fm.fontManager.addfont(_f)
        except Exception:
            pass

    plt.rcParams.update(_BASE_RC)
    plt.rcParams.update(
        {
            "text.usetex": False,
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Nimbus Roman", "Liberation Serif", "DejaVu Serif"],
        }
    )

RESULTS_DIR = Path(__file__).parent.parent / 'results'
FIGURES_DIR = Path(__file__).parent


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------
def _pow2_label(n):
    """Return '$2^{k}$' for n = 2^k."""
    k = int(round(np.log2(n)))
    return f'$2^{{{k}}}$'


def _install_pow2_formatter(ax, ns, all_labels):
    """Label power-of-two x ticks through a FuncFormatter.

    A plain `set_xticklabels([...])` installs a FixedFormatter, and matplotlib
    re-applies that formatter on every redraw — which silently undoes any later
    `set_text('')` thinning. Formatting on demand keeps full control: the
    `stride` in `state` can be raised after measuring the laid-out labels and
    survives every subsequent draw. All tick positions stay in place, so the
    grid and tick marks remain complete for every power of two.
    """
    ns_arr = np.asarray(ns, dtype=float)
    state = {'stride': 1}
    ax._pow2_tick_state = state

    # Index of each label inside the *visible* subsequence. Thinning must count
    # visible labels only, otherwise the labels that are hidden for n < 1024
    # still consume stride slots and too many visible labels get dropped.
    visible_pos = {}
    for j, i in enumerate(i for i in range(len(ns)) if all_labels or ns[i] >= 1024):
        visible_pos[i] = j

    def _fmt(x, _pos):
        i = int(np.argmin(np.abs(ns_arr - x)))
        j = visible_pos.get(i)
        if j is None:
            return ''
        return _pow2_label(ns[i]) if j % state['stride'] == 0 else ''

    ax.set_xticks(ns)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(_fmt))


def set_pow2_xticks(ax, ns):
    """Power-of-two x labels; hide labels for n < 1024."""
    _install_pow2_formatter(ax, ns, all_labels=False)


def set_pow2_xticks_all(ax, ns):
    """Power-of-two x labels for every n."""
    _install_pow2_formatter(ax, ns, all_labels=True)


# ---------------------------------------------------------------------------
# Final-geometry x-label thinning
# ---------------------------------------------------------------------------
_MIN_XLABEL_GAP = 4.0   # pt between neighbouring x tick labels


def refine_xtick_labels(fig, min_gap=_MIN_XLABEL_GAP):
    """Raise the label stride until neighbours are >= min_gap apart.

    Measured on the real laid-out geometry, so it adapts to whatever canvas
    size the panels are authored at: it keeps the largest number of labels that
    actually fit.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    def _min_gap(ax):
        spans = sorted((t.get_window_extent(renderer=renderer).x0,
                        t.get_window_extent(renderer=renderer).x1)
                       for t in ax.get_xticklabels() if t.get_text().strip())
        if len(spans) < 2:
            return float('inf')
        return min(spans[i + 1][0] - spans[i][1] for i in range(len(spans) - 1))

    for ax in fig.axes:
        state = getattr(ax, '_pow2_tick_state', None)
        if state is None:
            continue
        for stride in (1, 2, 3, 4, 6, 8):
            state['stride'] = stride
            fig.canvas.draw()
            if _min_gap(ax) >= min_gap:
                break


def merge_steps(steps_dict):
    """Re-map the 11-step result layout to the 9-step paper algorithm.

    Old steps 8 + 9 + 10 are summed into new step 8.
    Old step 11 becomes new step 9.
    Steps 1–7 are unchanged.
    """
    merged = {}
    for k in range(1, 8):
        if k in steps_dict:
            merged[k] = steps_dict[k]
    # New S8 = old S8 + S9 + S10 (combined time)
    t8 = sum(steps_dict.get(i, {}).get('time', 0.0) for i in [8, 9, 10])
    if any(i in steps_dict for i in [8, 9, 10]):
        merged[8] = {'name': 'Combine and sort write records', 'time': t8}
    # New S9 = old S11
    if 11 in steps_dict:
        merged[9] = steps_dict[11]
    return merged


# 9-step short labels used by all step breakdown figures
STEP_SHORT_LABELS = {
    1: 'S1: Sort query',
    2: 'S2: First-occurrence indicators',
    3: 'S3: Membership polynomial',
    4: 'S4: Poly eval at table',
    5: 'S5: Normalize indicators',
    6: 'S6: Sort unused elems',
    7: 'S7: Sort unfilled pos',
    8: 'S8: Combine & sort records',
    9: 'S9: Output permutation',
}

# Minimal step labels (just "S1" … "S9") for the compact step-breakdown legend.
STEP_TINY_LABELS = {k: f'S{k}' for k in range(1, 10)}


def _sci_label(x):
    """Format x as scientific notation."""
    float_digits_after_decimal = 0

    if x <= 0:
        return '0'
    exp = int(np.floor(np.log10(abs(x))))
    coef = x / (10.0 ** exp)
    coef = round(coef, float_digits_after_decimal)
    if coef == int(coef):
        coef = int(coef)

    if coef == 10:
        coef = 1
        exp += 1

    if coef == 1:
        return f'$10^{{{exp}}}$'
    return f'${coef:.{float_digits_after_decimal}f}\\times10^{{{exp}}}$'

_SCI_FORMATTER = ticker.FuncFormatter(lambda x, pos: _sci_label(x))


# ---------------------------------------------------------------------------
# Style helpers
# ---------------------------------------------------------------------------

def _grid(ax):
    """Faint major-grid only: structure without visual weight."""
    ax.grid(True, which='major', color=GRID_COLOR, linewidth=GRID_LW, alpha=1.0)


def _ticks(ax, labelsize=FONT_TICK):
    ax.tick_params(which='major', labelsize=labelsize,
                   width=TICK_LW, length=TICK_LEN, direction='out')
    ax.tick_params(which='minor', width=TICK_LW, length=TICK_LEN * 0.5)


def _spines(ax, lw=SPINE_LW):
    for side in ('left', 'right', 'top', 'bottom'):
        ax.spines[side].set_linewidth(lw)


def _axis_labels(ax, xlabel, ylabel, labelsize=FONT_BASE):
    # Empty strings are skipped so that dual-axis panels keep the coloured
    # y-labels they set on each twin axes.
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=labelsize, labelpad=0.6)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=labelsize, labelpad=0.9)


def _frame_legend(leg):
    """Apply the compact, lightly framed legend look."""
    leg.get_frame().set_linewidth(MARKER_EW)
    leg.get_frame().set_edgecolor('0.55')
    leg.get_frame().set_facecolor('white')
    leg.get_frame().set_alpha(0.92)
    return leg


def _legend(ax, handles=None, labels=None, fontsize=FONT_LEGEND, **kw):
    kw.setdefault('frameon', True)
    kw.setdefault('fancybox', False)
    kw.setdefault('handlelength', 1.2)
    kw.setdefault('handleheight', 0.7)
    kw.setdefault('handletextpad', 0.4)
    kw.setdefault('borderpad', 0.22)
    kw.setdefault('labelspacing', 0.16)
    kw.setdefault('borderaxespad', 0.22)
    if handles is None:
        leg = ax.legend(fontsize=fontsize, **kw)
    else:
        leg = ax.legend(handles, labels, fontsize=fontsize, **kw)
    return _frame_legend(leg)


# ---------------------------------------------------------------------------
# Legend placement that cannot cover the curves
# ---------------------------------------------------------------------------
LEGEND_LOCS = ('upper left', 'upper right', 'lower left', 'lower right',
               'upper center', 'lower center', 'center left', 'center right')


def _polyline_samples(ax, per_seg=10):
    """Densely sampled data polylines of one axes, in display coordinates."""
    out = []
    tr = ax.transData
    for ln in ax.get_lines():
        x = np.asarray(ln.get_xdata(), float)
        y = np.asarray(ln.get_ydata(), float)
        if x.size == 0:
            continue
        m = np.isfinite(x) & np.isfinite(y)
        x, y = x[m], y[m]
        if x.size == 0:
            continue
        xs, ys = [x[:1]], [y[:1]]
        for i in range(len(x) - 1):
            t = np.linspace(0.0, 1.0, per_seg)[1:]
            # interpolate in log space when possible: log axes look straight
            if x[i] > 0 and x[i + 1] > 0 and y[i] > 0 and y[i + 1] > 0:
                xi = np.exp(np.log(x[i]) * (1 - t) + np.log(x[i + 1]) * t)
                yi = np.exp(np.log(y[i]) * (1 - t) + np.log(y[i + 1]) * t)
            else:
                xi = x[i] * (1 - t) + x[i + 1] * t
                yi = y[i] * (1 - t) + y[i + 1] * t
            xs.append(xi); ys.append(yi)
        pts = np.column_stack([np.concatenate(xs), np.concatenate(ys)])
        out.append(tr.transform(pts))
    return out


def _legend_cover_count(fig, leg):
    """How many curve samples fall inside the legend box."""
    bb = leg.get_window_extent(fig.canvas.get_renderer())
    n = 0
    for a in fig.axes:
        for P in _polyline_samples(a):
            n += int(((P[:, 0] >= bb.x0) & (P[:, 0] <= bb.x1) &
                      (P[:, 1] >= bb.y0) & (P[:, 1] <= bb.y1)).sum())
    return n


def place_legend(ax, handles=None, labels=None, fontsize=FONT_LEGEND,
                 candidates=LEGEND_LOCS, **kw):
    """Build the legend and move it until it covers no curve sample.

    Fixed corner picks (loc='upper left') put the box straight on top of the
    rising curves once the panels get small. This measures, for every candidate
    corner, how many sampled points of every curve in the figure fall inside
    the legend box and keeps the first position that covers none (falling back
    to the least-bad one). Placement is therefore data-driven, not guessed.
    """
    leg = _legend(ax, handles, labels, fontsize=fontsize, loc=candidates[0], **kw)
    fig = ax.figure
    fig.canvas.draw()
    best_loc, best_n = candidates[0], _legend_cover_count(fig, leg)
    for loc in candidates[1:]:
        if best_n == 0:
            break
        leg.set_loc(loc)
        fig.canvas.draw()
        n = _legend_cover_count(fig, leg)
        if n < best_n:
            best_loc, best_n = loc, n
    leg.set_loc(best_loc)
    fig.canvas.draw()
    leg._cover_count = best_n
    return leg


def place_legend_below(ax, handles=None, labels=None, ncol=2,
                       fontsize=FONT_LEGEND, pad_in=0.03, **kw):
    """Legend as a compact strip just under the x-axis label.

    Inside a 1.2in frame the legend takes 11-30% of the plot area, so no corner
    is free and it lands on the curves (measured: up to 90 sampled points
    covered). Below the axes it cannot cover anything by construction, and the
    strip also reads better in a 4-up grid. The strip is positioned from the
    measured extent of the x-axis label, so it hugs the label for every panel.
    """
    if handles is None:
        handles, labels = ax.get_legend_handles_labels()
    kw.setdefault('frameon', False)
    kw.setdefault('fancybox', False)
    kw.setdefault('handlelength', 1.3)
    kw.setdefault('handleheight', 0.7)
    kw.setdefault('handletextpad', 0.4)
    kw.setdefault('borderpad', 0.05)
    kw.setdefault('labelspacing', 0.15)
    kw.setdefault('columnspacing', 0.9)
    kw.setdefault('borderaxespad', 0.0)
    leg = ax.legend(handles, labels, loc='upper center',
                    bbox_to_anchor=(0.5, 0.0), bbox_transform=ax.transAxes,
                    fontsize=fontsize, ncol=ncol, **kw)
    fig = ax.figure
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    xlab = ax.xaxis.label.get_window_extent(renderer=rend)
    y_px = xlab.y0 - pad_in * fig.dpi
    y_axes = ax.transAxes.inverted().transform((0.0, y_px))[1]
    leg.set_bbox_to_anchor((0.5, y_axes), transform=ax.transAxes)
    fig.canvas.draw()
    return leg


def _decade_ylim(vals, min_decades=2):
    """Expand a log axis to whole decades so it carries >= min_decades+1 ticks.

    The party-count panels hold either a constant (broadcasts are identical for
    every party count) or a narrow range, so a "tight around the data" limit
    showed a single labelled tick. Log axes label whole decades, hence the range
    is snapped outward to decade boundaries until at least `min_decades`
    decades are covered — 2 decades -> 3 labels.
    """
    v = [float(x) for x in vals if x is not None and np.isfinite(x) and x > 0]
    if not v:
        return None
    lo = int(np.floor(np.log10(min(v))))
    hi = int(np.ceil(np.log10(max(v))))
    while hi - lo < min_decades:
        lo -= 1
        if hi - lo < min_decades:
            hi += 1
    return 10.0 ** lo, 10.0 ** hi


def _pad_log_ylim(ax, min_pt=6.0):
    """Keep the extreme data point clear of the axis edge.

    Matplotlib's default log margin is 5% of the span; on a 7-decade axis that
    is only ~3.2 pt, i.e. less than the 4 pt marker diameter, so the marker of
    the lowest curve (Verify in the phase panel) is drawn right on the bottom
    spine and reads as if it were cut off. Expand the limits in log space until
    the lowest/highest point is at least `min_pt` from the edge.
    """
    if ax.get_yscale() != 'log':
        return
    vals = []
    for ln in ax.get_lines():
        v = np.asarray(ln.get_ydata(), float)
        v = v[np.isfinite(v) & (v > 0)]
        if v.size:
            vals.append(v)
    if not vals:
        return
    v = np.concatenate(vals)
    lo_d, hi_d = np.log10(v.min()), np.log10(v.max())
    cur_lo, cur_hi = np.log10(ax.get_ylim())
    span = max(cur_hi - cur_lo, 1e-9)
    need = min_pt / (FRAME_H * 72.0) * span
    ax.set_ylim(10 ** min(cur_lo, lo_d - need), 10 ** max(cur_hi, hi_d + need))


def _data_min(ax):
    """Smallest plotted y value (inf if the axes has no line data)."""
    vals = []
    for ln in ax.get_lines():
        v = np.asarray(ln.get_ydata(), float)
        v = v[np.isfinite(v)]
        if v.size:
            vals.append(v.min())
    return min(vals) if vals else float('inf')


def _pad_y_ticks(ax, min_pt=5.0):
    """Keep the outermost y tick labels clear of the frame.

    `_decade_ylim` snaps a log axis to whole decades so that it carries enough
    labelled ticks, but that also puts the first and last tick exactly on the
    bottom and top spine, so their labels read as if they were clipped. Expand
    the limits until the outermost labels sit at least `min_pt` inside the frame.
    Works for linear axes too (the speedup panel is linear and had the same
    problem).

    Two details cost time here. First, the tick locations are pinned with
    set_yticks before padding: expanding a log range otherwise lets LogLocator
    choose a finer step, which adds a new outermost tick that then needs padding
    again, and that feedback loop once blew a panel's range up from 1e-7..1e4 to
    1e-30..1e16. Second, the padding is driven by the *measured* label bounding
    boxes rather than by the tick positions: a label's box is not centred on its
    tick (it extends about 6.5pt below and 3.9pt above), so padding by a constant
    number of points leaves one side visibly closer to the frame.
    """
    is_log = ax.get_yscale() == 'log'
    to_u = np.log10 if is_log else (lambda v: v)
    fr_u = (lambda u: 10.0 ** u) if is_log else (lambda u: u)

    def limits():
        lo, hi = ax.get_ylim()
        if is_log and (lo <= 0 or hi <= 0):
            return None
        return to_u(np.float64(lo)), to_u(np.float64(hi))

    lh = limits()
    if lh is None:
        return
    ticks = [t for t in ax.get_yticks() if ax.get_ylim()[0] <= t <= ax.get_ylim()[1]]
    if is_log:
        ticks = [t for t in ticks if t > 0]
    if len(ticks) < 2:
        return
    ax.set_yticks(ticks)                      # pin, so the tick set is stable
    fig = ax.figure
    for _ in range(5):
        fig.canvas.draw()
        rend = fig.canvas.get_renderer()
        axbb = ax.get_window_extent(rend)
        lo, hi = ax.get_ylim()
        inrange = [(tick, lab) for tick, lab in zip(ax.get_yticks(), ax.get_yticklabels())
                   if lo <= tick <= hi and lab.get_text().strip()]
        if len(inrange) < 2:
            return
        boxes = [lab.get_window_extent(rend) for _, lab in inrange]
        clear_top = (axbb.y1 - max(b.y1 for b in boxes)) / fig.dpi * 72.0
        clear_bot = (min(b.y0 for b in boxes) - axbb.y0) / fig.dpi * 72.0
        if clear_top >= min_pt and clear_bot >= min_pt:
            return
        lo_u, hi_u = to_u(np.float64(lo)), to_u(np.float64(hi))
        span = max(float(hi_u - lo_u), 1e-12)
        per_pt = span / (FRAME_H * 72.0)
        grow_lo = max(0.0, min_pt - clear_bot) * per_pt
        grow_hi = max(0.0, min_pt - clear_top) * per_pt
        new_lo = float(lo_u - grow_lo)
        # On a linear axis whose data are all non-negative, do not push the
        # lower limit below 0 just to buy clearance: a quantity such as a
        # speedup cannot be negative, so stopping at 0 is the semantically right
        # answer (with the 0 tick sitting on the bottom spine). The test must
        # look at the data, not at the current limit: matplotlib's own autoscale
        # already produces a negative lower limit for this panel.
        if not is_log and new_lo < 0.0 and _data_min(ax) >= 0.0:
            new_lo = 0.0
        new_hi = float(hi_u + grow_hi)
        if abs(new_lo - lo_u) < 1e-9 and abs(new_hi - hi_u) < 1e-9:
            return          # the clamp already holds it still; retrying will not move it
        ax.set_ylim(fr_u(new_lo), fr_u(new_hi))


def pad_log_axes(fig, min_pt=6.0):
    """Apply the clearance rule to every log y-axis of a panel."""
    for a in fig.axes:
        _pad_log_ylim(a, min_pt)
        _pad_y_ticks(a, min_pt=5.0)


def _base_kw(color, marker):
    """Keyword arguments for a strawman / baseline series."""
    return dict(color=color, marker=marker, **BASE_STYLE)


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

def parse_duration_to_seconds(s):
    """Parse a duration string like '2.141s', '190.593ms', '65.082ms', '732ns', '1.362µs' into seconds."""
    s = s.strip()
    if s.endswith('ns'):
        return float(s[:-2]) * 1e-9
    if s.endswith('µs'):
        return float(s[:-2]) * 1e-6
    if s.endswith('ms'):
        return float(s[:-2]) * 1e-3
    if s.endswith('s'):
        return float(s[:-1])
    raise ValueError(f"Unknown duration format: {s!r}")


def parse_result_file(path):
    """Parse a result file and return a dict with all key metrics."""
    with open(path, 'r', errors='replace') as f:
        text = f.read()

    result = {}

    # n and m from summary line
    m = re.search(r'n=(\d+), m=(\d+)', text)
    if m:
        result['n'] = int(m.group(1))
        result['m'] = int(m.group(2))

    # MPC permutation time (new or naive)
    m = re.search(r'secure_oblivious_lookup_permutation: ([\d.]+)s', text)
    if m:
        result['mpc_perm_time'] = float(m.group(1))

    m = re.search(r'naive_secure_oblivious_lookup_permutation: ([\d.]+)s', text)
    if m:
        result['naive_mpc_perm_time'] = float(m.group(1))

    # Proof timing
    for key, label in [
        ('proof_setup', r'Proof setup:\s+([\d.]+\w+)'),
        ('proof_gen', r'Proof generation:\s+([\d.]+\w+)'),
        ('proof_ver', r'Proof verification:\s+([\d.]+\w+)'),
    ]:
        m = re.search(label, text)
        if m:
            result[key] = parse_duration_to_seconds(m.group(1))

    # Network stats
    m = re.search(r'bytes_sent:\s+(\d+)', text)
    if m:
        result['bytes_sent'] = int(m.group(1))
    m = re.search(r'bytes_recv:\s+(\d+)', text)
    if m:
        result['bytes_recv'] = int(m.group(1))
    m = re.search(r'broadcasts:\s+(\d+)', text)
    if m:
        result['broadcasts'] = int(m.group(1))

    # Step-by-step timing (MPlookup only)
    step_pattern = re.compile(
        r'End:\s+secure_oblivious_lookup_permutation step (\d+): ([^\n]+?) n=\d+ m=\d+ ([\d.]+\S+)'
    )
    steps = {}
    for sm in step_pattern.finditer(text):
        step_num = int(sm.group(1))
        step_name = sm.group(2).strip()
        step_time = parse_duration_to_seconds(sm.group(3))
        steps[step_num] = {'name': step_name, 'time': step_time}
    if steps:
        result['steps'] = steps

    return result


def load_all_results():
    """Load 2-party results (scalability experiments across different N values)."""
    new_results = {}
    naive_results = {}

    for fname in RESULTS_DIR.iterdir():
        if not fname.name.endswith('.txt'):
            continue
        # Only load 2-party files for the scalability figures
        if fname.name.startswith('result_naive_2_parties_'):
            data = parse_result_file(fname)
            n = data.get('n')
            if n is not None:
                naive_results[n] = data
        elif fname.name.startswith('result_2_parties_'):
            data = parse_result_file(fname)
            n = data.get('n')
            if n is not None:
                new_results[n] = data

    return new_results, naive_results


def load_multi_party_results(n_target):
    """Load MPlookup results for 2, 4, 8 … parties at a given N."""
    party_results = {}
    pattern = re.compile(r'^result_(\d+)_parties_n_(\d+)\.txt$')
    for fname in RESULTS_DIR.iterdir():
        m = pattern.match(fname.name)
        if m:
            parties = int(m.group(1))
            n = int(m.group(2))
            if n == n_target:
                data = parse_result_file(fname)
                party_results[parties] = data
    return party_results


# ---------------------------------------------------------------------------
# Panel draw helpers  (a) – (h)
# Each receives an axes object and draws into it.
# ---------------------------------------------------------------------------

def _panel_a_preproc_time(ax, new_ns, new_times, common_ns, naive_times):
    """(a) Preprocessing / MPC-permutation time vs n."""
    ax.plot(new_ns, new_times, **Q_TIME, label='MPlookup')
    ax.plot(common_ns, naive_times, **_base_kw(Q_TIME['color'], 'o'), label='Strawman')
    _axis_labels(ax, 'Input table size', 'Time (s)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.yaxis.set_minor_locator(ticker.NullLocator())
    set_pow2_xticks_all(ax, new_ns)
    _grid(ax)
    _spines(ax)
    _ticks(ax)
    place_legend_below(ax, ncol=1)


def _panel_b_comm_metrics(ax, new_ns, new_bytes, common_ns, naive_bytes, new_bc, naive_bc):
    """(b) Bytes sent (left y) + broadcasts (right y) vs n – dual y-axis."""
    # Left y-axis: bytes sent
    l1, = ax.plot(new_ns, new_bytes, **Q_BYTES, label='MPlookup')
    l2, = ax.plot(common_ns, naive_bytes,
                  **_base_kw(Q_BYTES['color'], 's'), label='Strawman')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.yaxis.set_minor_locator(ticker.NullLocator())
    ax.set_ylabel('Bytes sent', color=Q_BYTES['color'], labelpad=0.9)
    ax.tick_params(axis='y', labelcolor=Q_BYTES['color'])

    # Right y-axis: broadcasts
    ax2 = ax.twinx()
    l3, = ax2.plot(new_ns, new_bc, **Q_BCAST, label='MPlookup')
    l4, = ax2.plot(common_ns, naive_bc,
                   **_base_kw(Q_BCAST['color'], '^'), label='Strawman')
    ax2.set_yscale('log')
    ax2.yaxis.set_minor_locator(ticker.NullLocator())
    # The x tick labels belong to the primary axes; do not draw them twice.
    ax2.tick_params(axis='x', labelbottom=False)
    ax2.set_ylabel('Broadcasts', color=Q_BCAST['color'], labelpad=0.9)
    ax2.tick_params(axis='y', labelcolor=Q_BCAST['color'])

    set_pow2_xticks_all(ax, new_ns)
    _grid(ax)
    _spines(ax)
    _spines(ax2)
    _ticks(ax)
    _ticks(ax2)
    _axis_labels(ax, 'Input table size', '')  # y-labels set above (two axes)

    # Legend laid out as 2x2:
    #     row 1  [filled marker] Bytes sent, MPlookup  [hollow] Strawman
    #     row 2  [filled marker] Broadcasts, MPlookup  [hollow] Strawman
    # One row per quantity; the right column carries only the scheme name, and
    # the row position plus the marker shape tie it to the quantity. That is
    # more compact than repeating the quantity on all four entries. matplotlib
    # fills legends column-major, so the handle order has to be
    # (bytes/ours, broadcasts/ours, bytes/base, broadcasts/base).
    #
    # The negative columnspacing is deliberate. The pgf backend (text.usetex)
    # overestimates the width of a text run by about 1.4em, which pushes the
    # second column ~0.155in past the end of the first column's text; even with
    # handlelength, handletextpad and columnspacing all set to 0 the gap does not
    # close (measured). -1.1 cancels that excess, bringing the text gap back to
    # ~0.10in while the two markers still keep ~0.07in of clear space and do not
    # overlap. After touching this, re-render panel_b.pgf and check that the two
    # columns stay clear of each other.
    place_legend_below(ax,
                       [l1, l3, l2, l4],
                       ['Bytes sent, MPlookup', 'Broadcasts, MPlookup',
                        'Strawman', 'Strawman'],
                       ncol=2, fontsize=8.0,
                       handlelength=0.45, handletextpad=0.15, columnspacing=-1.1,
                       labelspacing=0.2)


def _panel_c_speedup(ax, common_ns, speedup_mpc, speedup_bytes, speedup_bc):
    """(c) Speedup of MPlookup over Strawman."""
    ax.plot(common_ns, speedup_mpc, **Q_TIME, label='Time')
    ax.plot(common_ns, speedup_bytes, **Q_BYTES, label='Bytes sent')
    ax.plot(common_ns, speedup_bc, **Q_BCAST, label='Broadcasts')
    ax.axhline(1.0, color='0.45', linewidth=LW_REF, linestyle=(0, (4, 2.5)), zorder=1)
    _axis_labels(ax, 'Input table size', 'Speedup $(\\times)$')
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.yaxis.set_minor_locator(ticker.FixedLocator([1]))
    set_pow2_xticks(ax, common_ns)
    _grid(ax)
    _spines(ax)
    _ticks(ax)
    # 3 items in 2 columns (2 rows), the same legend height as (d)/(h)
    place_legend_below(ax, ncol=2, fontsize=8.0,
                       handlelength=0.4, handletextpad=0.15, columnspacing=0.45)


def _panel_d_proof_vs_n(ax, proof_ns, t_setup, t_preproc, t_permvan, t_verify):
    """(d) Wall-clock time of each proof phase vs n – log-log."""
    ax.plot(proof_ns, t_setup,   color=C_BLUE,   marker='o', label='Setup')
    ax.plot(proof_ns, t_preproc, color=C_ORANGE, marker='s', label='Preprocess')
    ax.plot(proof_ns, t_permvan, color=C_GREEN,  marker='^', label='PermVanish')
    ax.plot(proof_ns, t_verify,  color=C_RED,    marker='D', label='Verify')
    ax.set_xscale('log')
    ax.set_yscale('log')
    _axis_labels(ax, 'Input table size', 'Time (s)')
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.yaxis.set_minor_locator(ticker.NullLocator())
    set_pow2_xticks_all(ax, proof_ns)
    _grid(ax)
    _spines(ax)
    _ticks(ax)
    place_legend_below(ax, ncol=2, markerscale=0.9)


def _panel_e_step_abs(ax, plot_ns, step_data, step_nums, colors_steps):
    """(e) Absolute time in each algorithmic step vs n – log-log, compact legend.

    Legend shows only the step label ("S1"…"S9") with a marker dot and no
    connecting line, laid out in rows of 5.
    """
    markers = ['o', 's', '^', 'D', 'x', 'v', 'p', '*', 'h']
    for j, snum in enumerate(step_nums):
        ax.plot(plot_ns, step_data[:, j],
                marker=markers[j % len(markers)], color=colors_steps[j],
                linewidth=LW, markersize=MSP,
                label=STEP_TINY_LABELS.get(snum, f'S{snum}'))
    ax.set_xscale('log')
    ax.set_yscale('log')
    _axis_labels(ax, 'Input table size', 'Time (s)')
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.yaxis.set_minor_locator(ticker.NullLocator())
    set_pow2_xticks_all(ax, plot_ns)
    _grid(ax)
    _spines(ax)
    _ticks(ax)

    # Compact legend: marker only (no line), 5 items per row, tight spacing.
    # Reorder handles so that row-major display (S1..S5 in row 1) is achieved
    # despite matplotlib's column-major legend filling.
    handles, labels = ax.get_legend_handles_labels()
    n_items = len(handles)
    # 9 curves: the column count follows the frame width; when the frame gets
    # too narrow use 3 columns (3x3) so the legend cannot overflow it. At 5
    # columns the legend is 2 rows tall, matching the other panels.
    ncol_leg = 5
    nrows_leg = -(-n_items // ncol_leg)  # ceil division
    new_order = [row * ncol_leg + col
                 for col in range(ncol_leg)
                 for row in range(nrows_leg)
                 if row * ncol_leg + col < n_items]
    handles = [handles[i] for i in new_order]
    labels  = [labels[i]  for i in new_order]
    leg = place_legend_below(ax, handles, labels, fontsize=8.0, ncol=ncol_leg,
                             handlelength=0.4, handletextpad=0.12,
                             columnspacing=0.4, labelspacing=0.08)
    for handle in leg.legend_handles:
        handle.set_linewidth(0)


def _panel_f_party_time(ax, parties, total_times):
    """MPlookup preprocessing time vs #parties (single series)."""
    ax.plot(parties, total_times, **Q_TIME, label='MPlookup')
    ax.set_xscale('log')
    ax.set_yscale('log')
    _axis_labels(ax, 'Number of parties', 'Time (s)')
    ax.set_xticks(parties)
    ax.xaxis.set_major_formatter(ticker.FixedFormatter([str(p) for p in parties]))
    ax.xaxis.set_minor_locator(ticker.NullLocator())

    # Snapped to decades for the same reason: the party counts span too little
    # for a range that hugs the data to carry more than two labelled ticks
    lim = _decade_ylim(total_times)
    if lim:
        ax.set_ylim(*lim)

    ax.yaxis.set_minor_locator(ticker.NullLocator())
    _grid(ax)
    _spines(ax)
    _ticks(ax)
    # Only one curve, but it still gets a legend entry so all eight panels look
    # uniform
    place_legend_below(ax, ncol=1)


def _panel_g_party_comm(ax, parties, bytes_vals, bc_vals):
    """(g) Bytes sent (left y) + broadcasts (right y) vs #parties – log-log, dual y-axis."""
    l1, = ax.plot(parties, bytes_vals, **Q_BYTES, label='Bytes sent')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xticks(parties)
    ax.xaxis.set_major_formatter(ticker.FixedFormatter([str(p) for p in parties]))
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    # Snapped to decades so that both y axes carry 3 ticks (a constant
    # broadcasts series would otherwise get a single tick)
    lim_b = _decade_ylim(bytes_vals)
    if lim_b:
        ax.set_ylim(*lim_b)

    ax.yaxis.set_minor_locator(ticker.NullLocator())
    ax.set_ylabel('Bytes sent', color=Q_BYTES['color'], labelpad=0.9)
    ax.tick_params(axis='y', labelcolor=Q_BYTES['color'])

    ax2 = ax.twinx()
    ax2.set_yscale('log')
    l2, = ax2.plot(parties, bc_vals, **Q_BCAST, label='Broadcasts')
    # y=0.45 lifts the rotated label clear of the topmost right-hand tick label;
    # x is left to matplotlib's automatic placement (a manual
    # set_label_coords here collides with the wide $10^{k}$ tick labels).
    # The x tick labels belong to the primary axes; do not draw them twice.
    ax2.tick_params(axis='x', labelbottom=False)
    ax2.set_ylabel('Broadcasts', color=Q_BCAST['color'], labelpad=0.9, y=0.45)

    # Broadcasts do not change with the party count; snap to decades for labels
    lim_bc = _decade_ylim(bc_vals)
    if lim_bc:
        ax2.set_ylim(*lim_bc)

    ax2.yaxis.set_minor_locator(ticker.NullLocator())
    ax2.tick_params(axis='y', labelcolor=Q_BCAST['color'])

    _axis_labels(ax, 'Number of parties', '')
    _grid(ax)
    _spines(ax)
    _spines(ax2)
    _ticks(ax)
    _ticks(ax2)

    place_legend_below(ax, [l1, l2], [l1.get_label(), l2.get_label()], ncol=1)


def _panel_h_party_proof(ax, parties, t_setup, t_preproc, t_permvan, t_verify):
    """(h) All four proof phases vs #parties – log-log."""
    ax.plot(parties, t_setup,   color=C_BLUE,   marker='o', label='Setup')
    ax.plot(parties, t_preproc, color=C_ORANGE, marker='s', label='Preprocess')
    ax.plot(parties, t_permvan, color=C_GREEN,  marker='^', label='PermVanish')
    ax.plot(parties, t_verify,  color=C_RED,    marker='D', label='Verify')
    ax.set_xscale('log')
    ax.set_yscale('log')
    _axis_labels(ax, 'Number of parties', 'Time (s)')
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.yaxis.set_minor_locator(ticker.NullLocator())
    ax.set_xticks(parties)
    ax.xaxis.set_major_formatter(ticker.FixedFormatter([str(p) for p in parties]))
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    _grid(ax)
    _spines(ax)
    _ticks(ax)
    place_legend_below(ax, ncol=2)


# ---------------------------------------------------------------------------
# Generate the 8 panel PDFs for LaTeX figure* inclusion
# ---------------------------------------------------------------------------

PANEL_PAD = 0.12          # tight_layout pad (tight: as little white as possible)


def fig_main_evaluation(new_results, naive_results, n_target=1024):
    r"""Save panel_a … panel_h (pdf or pgf) to FIGURES_DIR."""
    # ------------------------------------------------------------------ data
    common_ns   = sorted(set(new_results) & set(naive_results))
    new_ns      = sorted(new_results.keys())

    new_times   = [new_results[n]['mpc_perm_time']           for n in new_ns]
    naive_times = [naive_results[n]['naive_mpc_perm_time']   for n in common_ns]
    new_bytes   = [new_results[n]['bytes_sent']               for n in new_ns]
    naive_bytes = [naive_results[n]['bytes_sent']             for n in common_ns]
    new_bc      = [new_results[n]['broadcasts']               for n in new_ns]
    naive_bc    = [naive_results[n]['broadcasts']             for n in common_ns]

    speedup_mpc   = [naive_results[n]['naive_mpc_perm_time'] / new_results[n]['mpc_perm_time']
                     for n in common_ns]
    speedup_bytes = [naive_results[n]['bytes_sent'] / new_results[n]['bytes_sent']
                     for n in common_ns]
    speedup_bc    = [naive_results[n]['broadcasts'] / new_results[n]['broadcasts']
                     for n in common_ns]

    # Step-by-step absolute times
    plot_ns_steps = [n for n in sorted(new_results.keys()) if 'steps' in new_results[n]]
    sample_merged = merge_steps(new_results[plot_ns_steps[0]]['steps'])
    step_nums     = sorted(sample_merged.keys())
    colors_steps  = [STEP_PALETTE[i % len(STEP_PALETTE)] for i in range(len(step_nums))]
    step_data     = np.zeros((len(plot_ns_steps), len(step_nums)))
    for i, n in enumerate(plot_ns_steps):
        merged = merge_steps(new_results[n]['steps'])
        for j, snum in enumerate(step_nums):
            step_data[i, j] = merged.get(snum, {}).get('time', 0.0)

    # Proof phases vs n
    proof_ns  = sorted(n for n in new_results if new_results[n].get('proof_setup') is not None)
    t_setup   = [new_results[n]['proof_setup']   for n in proof_ns]
    t_preproc = [new_results[n]['mpc_perm_time'] for n in proof_ns]
    t_permvan = [new_results[n]['proof_gen']      for n in proof_ns]
    t_verify  = [new_results[n]['proof_ver']      for n in proof_ns]

    # Multi-party data
    party_results = load_multi_party_results(n_target)
    parties   = sorted(party_results.keys())
    p_times   = [party_results[p]['mpc_perm_time']              for p in parties]
    p_bytes   = [party_results[p].get('bytes_sent',   float('nan')) for p in parties]
    p_bc      = [party_results[p].get('broadcasts',   float('nan')) for p in parties]
    p_setup   = [party_results[p].get('proof_setup',  float('nan')) for p in parties]
    p_preproc = [party_results[p].get('mpc_perm_time',float('nan')) for p in parties]
    p_permvan = [party_results[p].get('proof_gen',    float('nan')) for p in parties]
    p_verify  = [party_results[p].get('proof_ver',    float('nan')) for p in parties]

    # ------------------------------------------------------------------ save
    # Two passes:
    #  1) draw once on a generous canvas and measure the box this panel actually
    #     needs (frame + decorations);
    #  2) rebuild the figure at "that width, plus one height shared by the whole
    #     group" and save it without any tight cropping.
    # Four panels side by side then have strictly equal canvas heights and
    # aligned frames, and the canvas width holds only the white space that is
    # needed, so a \resizebox fraction of W_i/\textwidth scales by exactly 1.
    def _make_panel_figure():
        W = _MARGIN_L + FRAME_W + _MARGIN_R
        H = _MARGIN_B + FRAME_H + _MARGIN_T
        fig = plt.figure(figsize=(W, H))
        ax = fig.add_axes([_MARGIN_L / W, _MARGIN_B / H, FRAME_W / W, FRAME_H / H])
        return fig, ax

    PAD = 0.02

    def _measure(draw_fn, args):
        fig, ax = _make_panel_figure()
        draw_fn(ax, *args)
        pad_log_axes(fig)
        refine_xtick_labels(fig)
        fig.canvas.draw()
        tb = fig.get_tightbbox(fig.canvas.get_renderer())     # in inches
        pos = ax.get_position()
        W, H = fig.get_size_inches()
        m = dict(w=tb.width, h=tb.height,
                 top_gap=tb.y1 - (pos.y0 + pos.height) * H,
                 left_gap=pos.x0 * W - tb.x0)
        plt.close(fig)
        return m

    jobs = [
        ('panel_a', _panel_a_preproc_time, (new_ns, new_times, common_ns, naive_times)),
        ('panel_b', _panel_b_comm_metrics, (new_ns, new_bytes, common_ns, naive_bytes, new_bc, naive_bc)),
        ('panel_c', _panel_c_speedup, (common_ns, speedup_mpc, speedup_bytes, speedup_bc)),
        ('panel_d', _panel_d_proof_vs_n, (proof_ns, t_setup, t_preproc, t_permvan, t_verify)),
        ('panel_e', _panel_e_step_abs, (plot_ns_steps, step_data, step_nums, colors_steps)),
        ('panel_f', _panel_f_party_time, (parties, p_times)),
        ('panel_g', _panel_g_party_comm, (parties, p_bytes, p_bc)),
        ('panel_h', _panel_h_party_proof, (parties, p_setup, p_preproc, p_permvan, p_verify)),
    ]

    metrics = {name: _measure(fn, args) for name, fn, args in jobs}
    # One shared height and one shared gap above the frame.
    # Matching the heights alone is not enough: the ink above each frame (top
    # tick labels and the like) differs in height, so the frame tops come out
    # ragged (measured: up to 0.055in apart). Taking the largest top gap of the
    # group makes every frame start at the same distance from the top of its
    # canvas, and LaTeX's bottom alignment of equal-height boxes then puts the
    # four panels of a row in a straight line.
    h_common = max(m['h'] for m in metrics.values()) + 2 * PAD
    top_gap_common = max(m['top_gap'] for m in metrics.values())

    extension = 'pgf' if SAVE_FIG_FORMAT == 'pgf' else SAVE_FIG_FORMAT
    for name, fn, args in jobs:
        m = metrics[name]
        w_i = m['w'] + 2 * PAD
        fig = plt.figure(figsize=(w_i, h_common))
        x0 = (PAD + m['left_gap']) / w_i
        y_top = 1.0 - (PAD + top_gap_common) / h_common
        ax = fig.add_axes([x0, y_top - FRAME_H / h_common,
                           FRAME_W / w_i, FRAME_H / h_common])
        fn(ax, *args)
        pad_log_axes(fig)
        refine_xtick_labels(fig)
        plt.savefig(FIGURES_DIR / f'{name}.{extension}')      # the canvas is the final size
        plt.close()
        print(f"Saved {name}.{extension}  ({w_i:.3f} x {h_common:.3f} in)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    print(f"Loading results from: {RESULTS_DIR.resolve()}")
    new_results, naive_results = load_all_results()
    print(f"MPlookup results: N = {sorted(new_results.keys())}")
    print(f"Strawman results: N = {sorted(naive_results.keys())}")

    fig_main_evaluation(new_results, naive_results, n_target=1024)

    print("\nAll figures saved.")
