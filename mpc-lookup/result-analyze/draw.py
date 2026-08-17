"""
Experiment result analysis for Oblivious Lookup Argument.
Parses result files from mpc-lookup/results and generates figures.

Algorithm names:
  MPlookup  – the new O(N log²N) MPC algorithm
  Strawman  – the naive O(N²) baseline
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

if SAVE_FIG_FORMAT == "pgf":
    # https://matplotlib.org/stable/tutorials/text/pgf.html
    matplotlib.use("pgf")
    plt.rcParams.update(
        {
            "text.usetex": True,
            "pgf.texsystem": "pdflatex",
            "pgf.preamble": "\n".join(
                [
                    "\\usepackage[utf8x]{inputenc}",
                    "\\usepackage[T1]{fontenc}",
                    "\\usepackage{cmbright}",
                ]
            ),
            "pgf.rcfonts": False,
            "font.serif": [],  # use latex default
            "font.sans-serif": [],  # use latex default
            "font.monospace": [],  # use latex default
            "font.family": "serif",
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


def set_pow2_xticks(ax, ns):
    """Set x-ticks as power-of-two labels with no rotation; hide labels for n < 1024."""
    ax.set_xticks(ns)
    ax.set_xticklabels([_pow2_label(n) if n >= 1024 else '' for n in ns])


def set_pow2_xticks_all(ax, ns):
    """Set x-ticks as power-of-two labels for all n (no hiding)."""
    ax.set_xticks(ns)
    ax.set_xticklabels([_pow2_label(n) for n in ns])


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

TRANSPARENCY = 0.7

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
    """Load MPlookup results for 2, 4, and 8 parties at a given N."""
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
    ax.plot(new_ns,    new_times,   'o-',  color='#2ca02c', linewidth=1.5, markersize=4,
        label='MPlookup', alpha=TRANSPARENCY)
    ax.plot(common_ns, naive_times, 'o-', color='#d62728',  linewidth=1.5, markersize=4,
        label='Strawman', alpha=TRANSPARENCY)
    ax.set_xlabel('Input table size', fontsize=9)
    ax.set_ylabel('Preprocess Time (s)', fontsize=9)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.yaxis.set_minor_locator(ticker.NullLocator())
    ax.legend(fontsize=7, bbox_to_anchor=(0.40, 1), bbox_transform=ax.transAxes)
    ax.grid(True, which='both', alpha=0.3)
    set_pow2_xticks_all(ax, new_ns)
    ax.tick_params(labelsize=8)


def _panel_b_comm_metrics(ax, new_ns, new_bytes, common_ns, naive_bytes, new_bc, naive_bc):
    """(b) Bytes sent (left y) + broadcasts (right y) vs n – dual y-axis.
    """
    # Left y-axis: bytes sent
    l1, = ax.plot(new_ns,    new_bytes,   f'v-',  color='steelblue',
                  linewidth=1.5, markersize=4, label='Bytes, MPlookup', alpha=TRANSPARENCY)
    l2, = ax.plot(common_ns, naive_bytes, 'v:', color='steelblue',
                  linewidth=1.5, markersize=4, label='Bytes, Strawman', alpha=TRANSPARENCY)
    ax.set_xlabel('Input table size', fontsize=9)
    ax.set_ylabel('Bytes Sent', fontsize=9, color='steelblue')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(axis='y', labelcolor='steelblue', labelsize=8)
    ax.tick_params(axis='x', labelsize=8)
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.yaxis.set_minor_locator(ticker.NullLocator())

    # Right y-axis: broadcasts
    ax2 = ax.twinx()
    l3, = ax2.plot(new_ns,    new_bc,   f'^-', color='coral',
                   linewidth=1.5, markersize=4, label='Broadcasts, MPlookup', alpha=TRANSPARENCY)
    l4, = ax2.plot(common_ns, naive_bc, '^:', color='coral',
                   linewidth=1.5, markersize=4, label='Broadcasts, Strawman', alpha=TRANSPARENCY)
    ax2.set_ylabel('Broadcasts', fontsize=9, color='coral')
    ax2.set_yscale('log')
    ax2.tick_params(axis='y', labelcolor='coral', labelsize=8)
    ax2.yaxis.set_minor_locator(ticker.NullLocator())

    # Split into two legends (one per y-axis) so each fits in white space.
    ax.legend([l1, l2], [l1.get_label(), l2.get_label()],
              fontsize=7, loc='upper left', ncol=1,
              bbox_to_anchor=(0.005, 0.9975), bbox_transform=ax.transAxes)
    ax2.legend([l3, l4], [l3.get_label(), l4.get_label()],
               fontsize=7, loc='upper right', ncol=1,
               bbox_to_anchor=(1.015, 0.185), bbox_transform=ax2.transAxes)
    ax.grid(True, which='both', alpha=0.3)
    set_pow2_xticks_all(ax, new_ns)


def _panel_c_speedup(ax, common_ns, speedup_mpc, speedup_bytes, speedup_bc):
    """(c) Speedup of MPlookup over Strawman."""
    ax.plot(common_ns, speedup_mpc,   'o-',  color='steelblue',  linewidth=1.5, markersize=4,
        label='Time', alpha=TRANSPARENCY)
    ax.plot(common_ns, speedup_bytes, 's--', color='darkorange', linewidth=1.5, markersize=4,
        label='Bytes Sent', alpha=TRANSPARENCY)
    ax.plot(common_ns, speedup_bc,    '^:',  color='green',      linewidth=1.5, markersize=4,
        label='Broadcasts', alpha=TRANSPARENCY)
    ax.axhline(1.0, color='gray', linewidth=1.0, linestyle='--', alpha=TRANSPARENCY)
    ax.set_xlabel('Input table size', fontsize=9)
    ax.set_ylabel('Speedup $({\\times})$', fontsize=9)
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.yaxis.set_minor_locator(ticker.FixedLocator([1]))
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)
    set_pow2_xticks(ax, common_ns)
    ax.tick_params(labelsize=8)


def _panel_d_proof_vs_n(ax, proof_ns, t_setup, t_preproc, t_permvan, t_verify):
    """(d) Wall-clock time of each proof phase vs n – log-log."""
    ax.plot(proof_ns, t_setup,   'o-',  color='#1f77b4', linewidth=1.5, markersize=4,
        label='Setup', alpha=TRANSPARENCY)
    ax.plot(proof_ns, t_preproc, 's-',  color='#ff7f0e', linewidth=1.5, markersize=4,
        label='Preprocess', alpha=TRANSPARENCY)
    ax.plot(proof_ns, t_permvan, '^-',  color='#2ca02c', linewidth=1.5, markersize=4,
        label='PermVanish', alpha=TRANSPARENCY)
    ax.plot(proof_ns, t_verify,  'D-',  color='#d62728', linewidth=1.5, markersize=4,
        label='Verify', alpha=TRANSPARENCY)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Input table size', fontsize=9)
    ax.set_ylabel('Time (s)', fontsize=9)
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.yaxis.set_minor_locator(ticker.NullLocator())
    ax.legend(fontsize=7)
    ax.grid(True, which='both', alpha=0.3)
    set_pow2_xticks_all(ax, proof_ns)
    ax.tick_params(labelsize=8)


def _panel_e_step_abs(ax, plot_ns, step_data, step_nums, colors_steps):
    """(e) Absolute time in each algorithmic step vs n – log-log, compact legend.

    Legend shows only the step label ("S1"…"S9") with a marker dot and no
    connecting line, laid out in rows of 5.
    """
    markers = ['o', 's', '^', 'D', 'x', 'v', 'p', '*', 'h']
    lines = []
    for j, snum in enumerate(step_nums):
        ln, = ax.plot(plot_ns, step_data[:, j],
                      marker=markers[j % len(markers)], color=colors_steps[j],
                      linewidth=1.5, markersize=4, alpha=TRANSPARENCY,
                      label=STEP_TINY_LABELS.get(snum, f'S{snum}'))
        lines.append(ln)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Input table size', fontsize=9)
    ax.set_ylabel('Time (s)', fontsize=9)
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.yaxis.set_minor_locator(ticker.NullLocator())
    ax.grid(True, which='both', alpha=0.3)
    set_pow2_xticks_all(ax, plot_ns)
    ax.tick_params(labelsize=8)

    # Compact legend: marker only (no line), 5 items per row, tight spacing.
    # Reorder handles so that row-major display (S1..S5 in row 1) is achieved
    # despite matplotlib's column-major legend filling.
    handles, labels = ax.get_legend_handles_labels()
    n_items = len(handles)
    ncol_leg = 5
    nrows_leg = -(-n_items // ncol_leg)  # ceil division
    new_order = [row * ncol_leg + col
                 for col in range(ncol_leg)
                 for row in range(nrows_leg)
                 if row * ncol_leg + col < n_items]
    handles = [handles[i] for i in new_order]
    labels  = [labels[i]  for i in new_order]
    leg = ax.legend(handles, labels, fontsize=7, loc='upper left', ncol=ncol_leg,
                    handlelength=0.8, handletextpad=0.3,
                    columnspacing=0.5, borderpad=0.4, labelspacing=0.2,
                    bbox_to_anchor=(0.01, 1.02), bbox_transform=ax.transAxes)
    for handle in leg.legend_handles:
        handle.set_linewidth(0)


def _panel_f_party_time(ax, parties, total_times):
    """(f) MPlookup preprocessing time vs #parties – no legend (single line)."""
    ax.plot(parties, total_times, 'o-', color='steelblue', linewidth=1.5, markersize=5,
            alpha=TRANSPARENCY)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Number of Parties', fontsize=9)
    ax.set_ylabel('Time (s)', fontsize=9)
    ax.set_xticks(parties)
    ax.xaxis.set_major_formatter(ticker.FixedFormatter([str(p) for p in parties]))
    ax.xaxis.set_minor_locator(ticker.NullLocator())

    # Set y-axis limits tight around the data
    y_min = min(t for t in total_times if t > 0) * 0.5
    y_max = max(total_times) * 2
    ax.set_ylim(y_min, y_max)

    ax.yaxis.set_minor_locator(ticker.NullLocator())
    ax.grid(True, which='both', alpha=0.3)
    ax.tick_params(labelsize=8)


def _panel_g_party_comm(ax, parties, bytes_vals, bc_vals):
    """(g) Bytes sent (left y) + broadcasts (right y) vs #parties – log-log, dual y-axis."""
    l1, = ax.plot(parties, bytes_vals, 's-', color='steelblue',
                  linewidth=1.5, markersize=5, label='Bytes Sent', alpha=TRANSPARENCY)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Number of Parties', fontsize=9)
    ax.set_ylabel('Bytes Sent', fontsize=9, color='steelblue')
    ax.set_xticks(parties)
    ax.xaxis.set_major_formatter(ticker.FixedFormatter([str(p) for p in parties]))
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    # Set y-axis limits tight around the data
    y_min_bytes = min(b for b in bytes_vals if b > 0) * 0.5
    y_max_bytes = max(bytes_vals) * 2
    ax.set_ylim(y_min_bytes, y_max_bytes)

    ax.yaxis.set_minor_locator(ticker.NullLocator())
    ax.tick_params(axis='y', labelcolor='steelblue', labelsize=8)
    ax.tick_params(axis='x', labelsize=8)

    # adjust the label location according to the drawing result!
    # ax.yaxis.set_label_coords(-0.25, 0.5)

    ax2 = ax.twinx()
    ax2.set_yscale('log')
    l2, = ax2.plot(parties, bc_vals, '^--', color='coral',
                   linewidth=1.5, markersize=5, label='Broadcasts', alpha=TRANSPARENCY)
    ax2.set_ylabel('Broadcasts', fontsize=9, color='coral')

    # Broadcasts are constant across party counts; use tight range
    bc_val = bc_vals[0]
    ax2.set_ylim(bc_val * 0.5, bc_val * 2)

    ax2.yaxis.set_minor_locator(ticker.NullLocator())
    ax2.tick_params(axis='y', labelcolor='coral', labelsize=8)

    # adjust the label location according to the drawing result!
    ax2.yaxis.set_label_coords(1.05, 0.5)

    ax.legend([l1, l2], [l1.get_label(), l2.get_label()],
              fontsize=7, loc='upper left', bbox_to_anchor=(0.01, 0.80),
              bbox_transform=ax.transAxes)
    ax.grid(True, which='both', alpha=0.3)


def _panel_h_party_proof(ax, parties, t_setup, t_preproc, t_permvan, t_verify):
    """(h) All four proof phases vs #parties – log-log."""
    ax.plot(parties, t_setup,   'o-',  color='#1f77b4', linewidth=1.5, markersize=5,
        label='Setup', alpha=TRANSPARENCY)
    ax.plot(parties, t_preproc, 's-',  color='#ff7f0e', linewidth=1.5, markersize=5,
        label='Preprocess', alpha=TRANSPARENCY)
    ax.plot(parties, t_permvan, '^-',  color='#2ca02c', linewidth=1.5, markersize=5,
        label='PermVanish', alpha=TRANSPARENCY)
    ax.plot(parties, t_verify,  'D-',  color='#d62728', linewidth=1.5, markersize=5,
        label='Verify', alpha=TRANSPARENCY)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Number of Parties', fontsize=9)
    ax.set_ylabel('Time (s)', fontsize=9)
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.yaxis.set_minor_locator(ticker.NullLocator())
    ax.set_xticks(parties)
    ax.xaxis.set_major_formatter(ticker.FixedFormatter([str(p) for p in parties]))
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.legend(fontsize=7, loc='upper left', bbox_to_anchor=(0.01, 0.80),
              bbox_transform=ax.transAxes)
    ax.grid(True, which='both', alpha=0.3)
    ax.tick_params(labelsize=8)


# ---------------------------------------------------------------------------
# Generate the 8 panel PDFs for LaTeX figure* inclusion
# ---------------------------------------------------------------------------

PANEL_PAD  = 0.4           # tight_layout pad


def fig_main_evaluation(new_results, naive_results, n_target=1024):
    r"""Save panel_a.pdf … panel_h.pdf to FIGURES_DIR."""
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
    colors_steps  = plt.cm.tab20(np.linspace(0, 1, len(step_nums)))
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
    def _save(name, figsize, draw_fn, *args, **kwargs):
        fig, ax = plt.subplots(figsize=figsize)
        draw_fn(ax, *args, **kwargs)
        plt.tight_layout(pad=PANEL_PAD)
        extension = 'pgf' if SAVE_FIG_FORMAT == 'pgf' else SAVE_FIG_FORMAT
        plt.savefig(FIGURES_DIR / f'{name}.{extension}', bbox_inches="tight")
        plt.close()
        print(f"Saved {name}.{extension}")

    _save('panel_a', (3.0, 2.4), _panel_a_preproc_time, new_ns, new_times, common_ns, naive_times)
    _save('panel_b', (3.0, 2.4), _panel_b_comm_metrics, new_ns, new_bytes, common_ns, naive_bytes, new_bc, naive_bc)
    _save('panel_c', (3.0, 2.4), _panel_c_speedup, common_ns, speedup_mpc, speedup_bytes, speedup_bc)
    _save('panel_d', (3.0, 2.4), _panel_d_proof_vs_n, proof_ns, t_setup, t_preproc, t_permvan, t_verify)
    _save('panel_e', (3.0, 2.4), _panel_e_step_abs, plot_ns_steps, step_data, step_nums, colors_steps)
    _save('panel_f', (3.0, 2.4), _panel_f_party_time, parties, p_times)
    _save('panel_g', (3.0, 2.4), _panel_g_party_comm, parties, p_bytes, p_bc)
    _save('panel_h', (3.0, 2.4), _panel_h_party_proof, parties, p_setup, p_preproc, p_permvan, p_verify)


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
