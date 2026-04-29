"""
compute-evaluation-data.py
==========================
Reproduces the full text of the "Implementation and Evaluation" section of the
paper, with every numerical value computed directly from the raw experiment
result files.

Run from any directory:
    python ref/experiments/result-analyze/compute-evaluation-data.py

The output is the evaluation section paragraphs in LaTeX source form, with all
numbers derived from the data.  The output should be identical to the
corresponding passages in mplookup.tex.
"""

import re
import sys
import math
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
RESULTS_DIR = Path(__file__).parent.parent / 'results'

# ---------------------------------------------------------------------------
# Parsing helpers (mirrors draw.py logic)
# ---------------------------------------------------------------------------

def parse_duration_to_seconds(s):
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
    with open(path, 'r', errors='replace') as f:
        text = f.read()

    result = {}

    m = re.search(r'n=(\d+), m=(\d+)', text)
    if m:
        result['n'] = int(m.group(1))
        result['m'] = int(m.group(2))

    m = re.search(r'secure_oblivious_lookup_permutation: ([\d.]+)s', text)
    if m:
        result['mpc_perm_time'] = float(m.group(1))

    m = re.search(r'naive_secure_oblivious_lookup_permutation: ([\d.]+)s', text)
    if m:
        result['naive_mpc_perm_time'] = float(m.group(1))

    for key, label in [
        ('proof_setup', r'Proof setup:\s+([\d.]+\w+)'),
        ('proof_gen',   r'Proof generation:\s+([\d.]+\w+)'),
        ('proof_ver',   r'Proof verification:\s+([\d.]+\w+)'),
    ]:
        m = re.search(label, text)
        if m:
            result[key] = parse_duration_to_seconds(m.group(1))

    m = re.search(r'bytes_sent:\s+(\d+)', text)
    if m:
        result['bytes_sent'] = int(m.group(1))
    m = re.search(r'broadcasts:\s+(\d+)', text)
    if m:
        result['broadcasts'] = int(m.group(1))

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


def merge_steps(steps_dict):
    """Re-map the 11-step result layout to the 9-step paper algorithm.

    Old S8+S9+S10 → new S8 (Combine and sort write records)
    Old S11       → new S9 (Construct output permutation)
    S1–S7 unchanged.
    """
    merged = {}
    for k in range(1, 8):
        if k in steps_dict:
            merged[k] = steps_dict[k]
    t8 = sum(steps_dict.get(i, {}).get('time', 0.0) for i in [8, 9, 10])
    if any(i in steps_dict for i in [8, 9, 10]):
        merged[8] = {'name': 'Combine and sort write records', 'time': t8}
    if 11 in steps_dict:
        merged[9] = steps_dict[11]
    return merged


def load_2party_results():
    new_results = {}
    naive_results = {}
    for fname in RESULTS_DIR.iterdir():
        if not fname.name.endswith('.txt'):
            continue
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


def load_multi_party_results(n_target=1024):
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
# Helpers
# ---------------------------------------------------------------------------

def pct(num, denom):
    return 100.0 * num / denom


def theoretical_doubling_ratio(n):
    """Expected T(2n)/T(n) for O(n log^2 n): 2*(log2(2n)/log2(n))^2."""
    k = math.log2(n)
    return 2.0 * ((k + 1) / k) ** 2


# ---------------------------------------------------------------------------
# Main computation
# ---------------------------------------------------------------------------

def main():
    new_results, naive_results = load_2party_results()
    party_results = load_multi_party_results(n_target=1024)

    common_ns      = sorted(set(new_results) & set(naive_results))
    new_ns         = sorted(new_results.keys())
    parties_sorted = sorted(party_results.keys())

    # ------------------------------------------------------------------
    # TeX formatting helpers
    # ------------------------------------------------------------------

    def tex_n(n):
        """Return integer/float as LaTeX with comma-thousands separators.

        e.g. 8192 → '8{,}192'
        """
        return f"{int(round(n)):,}".replace(',', '{,}')

    def oxford_join(items):
        """Join a list with Oxford comma: ['A','B','C'] → 'A, B, and C'."""
        if len(items) == 1:
            return items[0]
        if len(items) == 2:
            return f"{items[0]} and {items[1]}"
        return ', '.join(items[:-1]) + f", and {items[-1]}"

    # ------------------------------------------------------------------
    # RQ1 — crossover, speedup, doubling ratios, bytes, broadcasts
    # ------------------------------------------------------------------

    # First N where MPlookup is faster than Strawman
    crossover_n, crossover_n_below = None, None
    for i, n in enumerate(common_ns):
        if naive_results[n]['naive_mpc_perm_time'] / new_results[n]['mpc_perm_time'] > 1.0:
            crossover_n       = n
            crossover_n_below = common_ns[i - 1] if i > 0 else None
            break

    n_hi       = 8192
    t_new_hi   = new_results[n_hi]['mpc_perm_time']
    t_naive_hi = naive_results[n_hi]['naive_mpc_perm_time']
    speedup_hi = t_naive_hi / t_new_hi

    # Three consecutive doublings: 2048 → 4096 → 8192 → 16384
    trio_ns     = [2048, 4096, 8192, 16384]
    obs_ratios  = [new_results[trio_ns[i + 1]]['mpc_perm_time'] /
                   new_results[trio_ns[i]]['mpc_perm_time']
                   for i in range(len(trio_ns) - 1)]
    theo_ratios = [theoretical_doubling_ratio(trio_ns[i])
                   for i in range(len(trio_ns) - 1)]

    bs_speedup = (naive_results[n_hi]['bytes_sent'] /
                  new_results[n_hi]['bytes_sent'])

    # First N where MPlookup issues fewer broadcasts than Strawman
    bc_crossover_n = next(
        (n for n in common_ns
         if new_results[n]['broadcasts'] < naive_results[n]['broadcasts']),
        None)

    # ------------------------------------------------------------------
    # RQ2 — proof-phase breakdown
    # ------------------------------------------------------------------

    # Average verification time in ms, rounded to nearest integer
    ver_ms_vals   = [new_results[n]['proof_ver'] * 1000
                     for n in new_ns if 'proof_ver' in new_results[n]]
    ver_ms_approx = round(sum(ver_ms_vals) / len(ver_ms_vals))

    # Max fraction of total wall-clock time from PermVanish + Setup
    # (total = Preprocess + PermVanish + Setup + Verify)
    max_pct_nonpreproc = 0.0
    for n in new_ns:
        d = new_results[n]
        if 'proof_gen' not in d:
            continue
        total = (d['mpc_perm_time']
                 + d.get('proof_gen',   0.0)
                 + d.get('proof_setup', 0.0)
                 + d.get('proof_ver',   0.0))
        p = pct(d.get('proof_gen', 0.0) + d.get('proof_setup', 0.0), total)
        max_pct_nonpreproc = max(max_pct_nonpreproc, p)
    # Round up to one decimal place to get a safe upper-bound claim
    pct_upper = math.ceil(max_pct_nonpreproc * 10) / 10

    # ------------------------------------------------------------------
    # RQ3 — multi-party scaling at N=1024
    # ------------------------------------------------------------------

    base_time         = party_results[2]['mpc_perm_time']
    non_base          = [p for p in parties_sorted if p != 2]
    non_base_times    = [party_results[p]['mpc_perm_time'] for p in non_base]
    non_base_speedups = [party_results[p]['mpc_perm_time'] / base_time
                         for p in non_base]

    pg_2  = party_results[2].get('proof_gen', float('nan'))
    pg_16 = party_results.get(16, {}).get('proof_gen', float('nan'))

    # Pre-build joined TeX fragments
    obs_str      = ', '.join(f'${r:.2f}\\times$' for r in obs_ratios)
    theo_str     = ', '.join(f'${r:.2f}\\times$' for r in theo_ratios)
    times_tex    = [f'${tex_n(t)}$\\,s' for t in non_base_times]
    speedups_tex = [f'${sp:.1f}\\times$' for sp in non_base_speedups]

    # ------------------------------------------------------------------
    # Print evaluation section text
    # ------------------------------------------------------------------

    print("=" * 70)
    print("EVALUATION SECTION — all numbers computed from raw data")
    print("=" * 70)

    # --- RQ1 ---
    print()
    print(r"\subsection{RQ1: Scalability and Comparison}")
    print()
    print(r"Figures~\ref{fig:rq1_total_time} and~\ref{fig:rq1_comm} plot preprocessing time, bytes sent, and broadcast count.")
    print(r"We confirm \MPlookup scales as $O(N \log^2 N)$ and the strawman as $O(N^2)$, in agreement with Theorem~\ref{thm:complexity}.")
    print()
    print(r"At small $N$, the strawman is faster due to lower constant factors: each iteration uses one $\mathcal{F}_\textrm{Eq}$ call, whereas \MPlookup invokes sorting networks and polynomial arithmetic.")
    print(f"The crossover lies between $N = {tex_n(crossover_n_below)}$ and $N = {tex_n(crossover_n)}$.")
    print(f"At $N = {tex_n(n_hi)}$, \\MPlookup achieves a ${speedup_hi:.2f}\\times$ speedup, completing in ${tex_n(t_new_hi)}$\\,s vs.\\ ${tex_n(t_naive_hi)}$\\,s for the strawman.")
    print(f"From $N = {tex_n(trio_ns[0])}$ to ${tex_n(trio_ns[-1])}$, the observed doubling ratios {obs_str} closely match the theoretical {theo_str}.")
    print()
    print(f"Figure~\\ref{{fig:rq1_speedup}} shows speedups in preprocessing time, bytes sent, and broadcast count; the time speedup turns decisively in \\MPlookup's favour beyond $N = {tex_n(crossover_n)}$.")
    print(f"Communication savings are more pronounced: at $N = {tex_n(n_hi)}$, \\MPlookup sends approximately ${bs_speedup:.0f}\\times$ fewer bytes, as the strawman's $O(N^2)$ equality comparisons each require a full MPC round.")
    print(f"Broadcast count is higher for \\MPlookup at small $N$ due to sorting and subproduct tree rounds; at $N = {tex_n(bc_crossover_n)}$, \\MPlookup also issues fewer broadcasts.")

    # --- RQ2 ---
    print()
    print(r"\subsection{RQ2: Step-Level Performance Breakdown}")
    print()
    print(r"Figure~\ref{fig:rq2_steps} shows the time cost contributed by each of the nine steps of Algorithm~\ref{alg:mplookup-preprocessing} across all tested table sizes.")
    print(r"Figure~\ref{fig:rq2_nlogn} plots the time cost of each step divided by $N\log^2 N$ for Steps~1, 4, 6, 7, 8. For Steps~1, 6, 7, 8, flat curves confirm $O(N\log^2 N)$ scaling. For Step~4, the line is close to $O(N \log N)$ scaling. This is because both subproduct tree construction ($O(N \log N)$) and multi-point evaluation ($O(N \log^2 N)$) are included in Step~4, and the subproduct tree construction contributes more time cost.")
    print()
    print(r"Figure~\ref{fig:rq2_mplookup_phases} breaks down the four phases of \MPlookup, setup, proof generation consisting of $\Pi_\textrm{Preprocess}$ and $\Pi_\textrm{PermVanish}$, and finally verification.")
    print(f"$\\Pi_\\textrm{{Preprocess}}$ overwhelmingly dominates the total cost across all tested sizes, while verification remains essentially constant at approximately ${ver_ms_approx}$\\,ms regardless of $N$.")
    print(f"$\\Pi_\\textrm{{PermVanish}}$ and setup both grow with $N$ but contribute less than ${pct_upper:.1f}\\%$ of total time, confirming that preprocessing is the decisive bottleneck.")

    # --- RQ3 ---
    print()
    print(r"\subsection{RQ3: Scaling with the Number of Parties}")
    print()
    print(f"Figure~\\ref{{fig:rq3_parties}} shows preprocessing time at $N=1{{,}}024$, $M=512$: from ${tex_n(base_time)}$\\,s for 2 parties to {oxford_join(times_tex)}---speedups of {oxford_join(speedups_tex)} over the 2-party baseline.")
    print(r"Figure~\ref{fig:rq3_comm} shows bytes sent scales similarly to time, while broadcast count remains constant across party counts, consistent with the communication structure.")
    print()
    print(f"In Figure~\\ref{{fig:rq3_phases}}, the setup time and verification time remain stable because they are single-user protocols, while $\\Pi_\\textrm{{PermVanish}}$ grows from ${pg_2:.1f}$\\,s at 2 parties to ${pg_16:.1f}$\\,s at 16 parties as additional parties require more collaborative KZG commitment rounds.")
    print()
    print(r"We note that the evaluation results are measured with \textsf{CompatCircuit} as the ABB, whose $\mathcal{F}_\textrm{LT}$ operation requires $O(\log_2 p)$ sequential communication rounds per comparison over the BLS12-377 scalar field, where $\log_2 p \approx 253$.")
    print(r"This is the primary source of the large constant factor observed throughout.")

    print()
    print("=" * 70)


if __name__ == '__main__':
    main()
