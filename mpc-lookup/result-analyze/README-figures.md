# Evaluation figures (Figure 7 of the paper)

This directory produces the eight panels of the single evaluation figure in
`mplookup-paper/mplookup.tex`. Read this before touching `draw.py`.

## Files

| file | role |
|---|---|
| `draw.py` | parses `../results/*.txt` and writes `panel_a.pgf` … `panel_h.pgf` |
| `config.py` | output format switch, must be `'pgf'` for the paper |
| `calibrate_figure7.py` | writes the `\resizebox` fractions and `\hspace*` gaps back into the paper |
| `compute-evaluation-data.py` | turns raw experiment logs into the numbers quoted in the text |

## Panel → letter mapping in the figure

The letter is assigned by **position** in the LaTeX figure (subfig numbering),
not by the panel file name:

```
row 1 :  a  b  c  d      (vary input table size N)
row 2 :  f  g  e  h      (vary the number of parties)
```

So `panel_f.pgf` is printed as **(e)** and `panel_e.pgf` as **(g)**. Keep the
`\subfloat` order in the tex and the `ROW1`/`ROW2` constants in
`calibrate_figure7.py` consistent with this.

## Encoding contract (colour and marker)

The eight panels fall into two groups, and only one of them is colour-coded
centrally.

**Quantity panels — (a), (b), (c), (e), (f).** These five share a single
mapping, defined once at the top of `draw.py` as `Q_TIME` / `Q_BYTES` /
`Q_BCAST` and never spelled out per panel:

| channel | encodes |
|---|---|
| colour + marker | which quantity: time = orange `o`, bytes sent = blue `s`, broadcasts = red `^` |
| fill + linestyle | which scheme: MPlookup = filled + solid, Strawman = hollow + dashed |

Time is orange, not blue, so that the two quantities that share a panel in
(b)/(f) keep the blue/red contrast the whole figure is read with; an earlier
red/green pairing for bytes vs broadcasts looked harsh and is not worth it.

Rules that follow from it, and that a future edit must not break:

1. **Every** series takes its colour and marker from the quantity it plots,
   the strawman included. Do not give the baseline a private colour: in (b) the
   two baselines land in the same region of the panel, and painting them one
   colour made them impossible to tell apart.
2. **Hollow + dashed is the exclusive signature of the strawman.** No other
   series may use a dashed line or a hollow face. This is what stops a hollow
   dashed triangle meaning "MPlookup's broadcasts" in one panel and "the
   strawman's broadcasts" in the next.
3. Legends below the axes have `handlelength` < marker diameter, so the legend
   handles show the marker only. The dash is only visible on the curves. Keep
   that in mind when judging whether a panel is readable.

**Component panels — (d), (g), (h).** Setup / Preprocess / PermVanish / Verify
and the nine steps S1–S9 keep their own `STEP_PALETTE` scheme; they are a
decomposition, not a quantity, so they are deliberately exempt. Two
consequences worth knowing:

* (d)/(h)'s `Preprocess` series *is* the quantity plotted by (a) and (e) —
  `mpc_perm_time`, identical numbers. It is drawn orange `s` there and orange
  `o` here: the colour now agrees (same quantity) but the marker does not, so
  the read-across is half-visual. That is accepted, not a bug.
* With time on orange, no (colour, marker) pair in the quantity set
  (orange `o`, blue `s`, red `^`) equals one in the component set
  (blue `o`, orange `s`, green `^`, red `D`). Keep it that way if you
  recolour: shift either group, never both, and re-check the four pairs.

## Layout contract (why the code looks the way it does)

1. **The canvas is authored at final print size.** `FRAME_W` × `FRAME_H` (inches)
   is the size of the axes frame *on the paper*. Both are fixed, and the axes is
   placed with `add_axes` instead of `tight_layout`, so the frame is identical in
   every panel even though the decorations differ. `tight_layout` cannot do this:
   the dual-axis panels would come out narrower.
2. **Decorations are measured, not guessed.** The figure is drawn twice: once on
   a generous canvas to measure the ink extent of each panel, then rebuilt at
   exactly that size (plus a shared extra top gap so every frame top lines up).
   A shared canvas height keeps the rows bottom-aligned in LaTeX.
3. **Legends live below the axes** (`place_legend_below`). Inside a ~1.1in frame a
   legend takes 11–30% of the plot and lands on the curves; the highest measured
   overlap was 90 sampled points. Below the axis it cannot cover anything, and
   the strip is positioned from the measured extent of the x-axis label.
4. **Every panel has a legend**, even the single-series one, so the eight panels
   look uniform.
5. **x tick labels are thinned at final size** (`refine_xtick_labels`). The
   thinning has to go through a `FuncFormatter`: `set_xticklabels` installs a
   `FixedFormatter`, and matplotlib re-applies it on every redraw, which silently
   undoes a `set_text('')` pass.

## Regenerating and re-installing the figures

```bash
cd oblivious-lookup-argument-internal/mpc-lookup/result-analyze
python3 draw.py                       # ~1 min, runs pdflatex for the pgf backend
python3 calibrate_figure7.py          # writes fracs + gaps into ../../../mplookup-paper/mplookup.tex
cp panel_*.pgf ../../../mplookup-paper/images/
cd ../../../mplookup-paper && pdflatex mplookup.tex && pdflatex mplookup.tex
```

`calibrate_figure7.py` is idempotent and only rewrites two sets of numbers in the
figure block; it never touches captions or labels.

## When the calibration must be re-run

Any change that alters a panel's pgf width: `FRAME_W`/`FRAME_H`, axis limits,
axis labels, tick labels, legend contents, or the order/number of `\subfloat`s.
If the script reports "No feasible solution", the row no longer fits
`\textwidth`: decrease `FRAME_W` in `draw.py` by about 0.03in and start again.

A pure colour/marker change does *not* move pgf widths (a rotated y-label's
width is its font height, not the label's length, and the below-axes legend
stays narrower than the canvas), so re-running the calibration after one
normally reports the same numbers and leaves `mplookup.tex` untouched. Renaming
a legend entry can still widen its panel — check the reported fracs. A different
matplotlib version can instead change the exported canvas *height*, because
`bbox_inches='tight'` crops a little more or less white space from release to
release; that shortens the float and moves the caption up, so compare the page
count even when no width moved.

## When the experiment data changes

Everything downstream of `../results/*.txt` has to be refreshed together. Updating
only one of these leaves the paper in a state where the numbers and the figure no
longer agree.

1. `python3 draw.py` — regenerates the panels. Axis ranges, the tick set and the
   legend widths all derive from the data, so most panels change width and the
   fractions/gaps baked into `mplookup.tex` go stale. **Re-run the calibration**,
   this step is not optional.
2. `python3 compute-evaluation-data.py` — prints the whole "Implementation and
   Evaluation" section with every number recomputed from the data. Diff that
   output against `mplookup.tex`: the RQ paragraphs quote concrete figures
   ("224.40x", "947 s", "less than 17.6%", the doubling ratios) and they go stale
   silently, because nothing fails when they do.
3. Re-check how many pages the compiled paper takes, and where the body ends.
   This is the largest float in the paper, and caption lines and legend rows both
   add height, so a figure that grew can shift text onto following pages without
   any LaTeX warning. Whatever page limit the current venue imposes, compare the
   page count (and the page carrying the Acknowledgment heading) before and after
   the change rather than trusting the compile to be warning-free.

## Three traps that cost real time

* **Font metrics differ between the `pdf` and `pgf` backend.** `draw.py` has two
  rcParams branches. The `pgf` branch uses `text.usetex`, so text is measured with
  LaTeX fonts and comes out wider; a legend measured at 1.81in in `pdf` mode was
  2.05in in `pgf` mode. Always judge widths from the generated `.pgf` files, i.e.
  from `calibrate_figure7.py`, never from a `pdf`-mode dry run.
* **Explicit `\hspace*` is unforgiving.** Each "panel + gap" must be a single
  physical source line with no leading indentation, otherwise the whitespace
  after `%` becomes inter-word glue (~3.5pt) and the columns drift apart again.
* **The sub-caption inset is load-bearing.** The gaps between panels are only
  0.02-0.10in (they have to be, four panels plus their decorations must fit
  `\textwidth`), while each sub-caption is justified across its whole panel
  width. Without the `\leftskip=2mm\rightskip=2mm` on every caption the four
  captions in a row read as one continuous block. Removing it "to clean up"
  reintroduces exactly that.

## Verifying the result

Render page 12 of the paper at 600dpi and detect the vertical frame lines
(runs whose length is ≈`FRAME_H`):

1. all eight frames should be exactly `FRAME_W` wide;
2. for every column, the frame left edges of the two rows must match
   (raster tolerance ≤0.005in);
3. within a row the frame tops must match.

A helper for step 1–3 is described in the comment block above the figure in
`mplookup.tex`.
