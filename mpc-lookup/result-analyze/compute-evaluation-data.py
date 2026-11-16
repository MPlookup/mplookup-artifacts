"""
compute-evaluation-data.py
==========================
Reproduces the full text of the "Implementation and Evaluation" section of the
paper, with every numerical value computed directly from the raw experiment
result files.

Run from any directory:
    python3 mpc-lookup/result-analyze/compute-evaluation-data.py

Check the output against the paper (no printing, exit status 1 on mismatch):
    python3 mpc-lookup/result-analyze/compute-evaluation-data.py --check

Re-run this (and diff against mplookup.tex) whenever ../results/*.txt changes:
the evaluation paragraphs quote concrete numbers that go stale silently.

The output is the evaluation section paragraphs in LaTeX source form, with all
numbers derived from the data.  The output should be identical to the
corresponding passages in mplookup.tex, which --check verifies mechanically
(whitespace-insensitive, so ordinary re-wrapping in the .tex is fine); it also
re-checks the two table-size ranges quoted in "Experimental Setup", which are
not part of the printed output but are derived from the same data.
"""

import contextlib
import io
import os
import re
import sys
import math
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
RESULTS_DIR = Path(__file__).parent.parent / 'results'

# The scalability sweep in ../results is the two-party one (load_2party_results
# globs result_2_parties_*), so that is the default party count quoted in text.
SCALABILITY_PARTIES = 2
# Convenience default only: the authors' working copy has the paper next to
# this artifacts repository.  Never assume it -- --check tells the user when
# the path is missing instead of silently passing.
DEFAULT_PAPER_TEX = (Path(__file__).parent.parent.parent.parent
                     / 'mplookup-paper' / 'mplookup.tex')

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


def tex_n(n):
    """Return an integer as LaTeX with comma-thousands separators.

    e.g. 8192 -> '8{,}192'
    """
    return f"{int(round(n)):,}".replace(',', '{,}')


def tex_ratio(table_over_query):
    """Render a fixed query/table size ratio, e.g. 2.0 -> 'N/2'."""
    k = table_over_query
    if k < 1 or abs(k - round(k)) > 1e-9:
        raise ValueError(f'table/query ratio {k} is not an integer k >= 1')
    return f'N/{round(k)}'


# ---------------------------------------------------------------------------
# Main computation
# ---------------------------------------------------------------------------

def main():
    new_results, naive_results = load_2party_results()
    party_results = load_multi_party_results(n_target=1024)

    common_ns      = sorted(set(new_results) & set(naive_results))
    new_ns         = sorted(new_results.keys())
    parties_sorted = sorted(party_results.keys())

    # The two table-size sweeps and the query/table ratios come from the runs
    # themselves, so the ranges quoted in "Experimental Setup" cannot drift.
    def sweep(results):
        """Table sizes run, and the (constant) table-to-query ratio of a sweep."""
        pairs = [(d['n'], d.get('m')) for d in results.values() if 'n' in d]
        sizes = sorted({n for n, _ in pairs})
        ratios = {n / m for n, m in pairs if m}
        if len(ratios) != 1:
            raise ValueError(f'runs disagree on the query/table ratio: {sorted(ratios)}')
        return sizes, ratios.pop()

    new_sizes,   new_ratio   = sweep(new_results)
    naive_sizes, _           = sweep(naive_results)
    party_sizes, party_ratio = sweep(party_results)
    if len(party_sizes) != 1:
        raise ValueError(f'the multi-party runs disagree on N: {party_sizes}')
    party_n = party_sizes[0]
    party_m = int(round(party_n / party_ratio))

    def size_list(sizes):
        """'8, 16, 32, \\ldots, 32768' for the doubling sweep actually run."""
        head = ', '.join(str(n) for n in sizes[:3])
        return '\\{' + head + ', \\ldots, ' + tex_n(sizes[-1]) + '\\}'

    # ------------------------------------------------------------------
    # TeX formatting helpers
    # ------------------------------------------------------------------

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

    n_hi       = max(common_ns)
    t_new_hi   = new_results[n_hi]['mpc_perm_time']
    t_naive_hi = naive_results[n_hi]['naive_mpc_perm_time']
    speedup_hi = t_naive_hi / t_new_hi

    # Three consecutive doublings: 2048 → 4096 → 8192 → 16384 → 32768
    trio_ns     = [2048, 4096, 8192, 16384, 32768]
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
    # ... and the largest N below it, where MPlookup still issues more
    bc_crossover_below = None
    if bc_crossover_n is not None:
        i = common_ns.index(bc_crossover_n)
        bc_crossover_below = common_ns[i - 1] if i > 0 else None

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

    base_time       = party_results[2]['mpc_perm_time']
    non_base        = [p for p in parties_sorted if p != 2]
    non_base_times  = [party_results[p]['mpc_perm_time'] for p in non_base]
    # Ratio of the p-party time to the 2-party time; each ratio exceeds 1, i.e.
    # adding parties makes preprocessing slower -- these are slowdowns.
    non_base_ratios = [party_results[p]['mpc_perm_time'] / base_time
                       for p in non_base]

    pg_2  = party_results[2].get('proof_gen', float('nan'))
    pg_16 = party_results.get(16, {}).get('proof_gen', float('nan'))

    # Pre-build joined TeX fragments
    obs_str       = ', '.join(f'${r:.2f}\\times$' for r in obs_ratios)
    theo_str      = ', '.join(f'${r:.2f}\\times$' for r in theo_ratios)
    times_tex     = [f'${tex_n(t)}$\\,s' for t in non_base_times]
    slowdowns_tex = [f'${sp:.1f}\\times$' for sp in non_base_ratios]

    # ------------------------------------------------------------------
    # Print evaluation section text
    # ------------------------------------------------------------------

    # --- Experimental Setup ---
    # The prose is fixed; the two table-size ranges, the ratio M = N/k, and the
    # party counts all come from the runs above.
    print()
    print(r"\subsection{Experimental Setup}")
    print()
    print(r"Our experiments were conducted on a server instance on the ESXi 8.0 virtualization platform.")
    print(r"Denote the prover count as $n$.")
    print(r"The server instance simulating $n$ provers was allocated with up to $4n$ vCPU cores from Intel Xeon 4214R @ 2.4 GHz, $16n$ GB of RAM, running Ubuntu Linux 24.04 LTS. We set up $n$ provers in a local network.")
    print(r"Note that our implementation is \textit{not} parallelized, but more than one CPU core is needed to handle the communication.")
    print()
    print(r"We focus our evaluation on the preprocessing phase (Algorithm~\ref{alg:mplookup-preprocessing}) and the strawman baseline (Algorithm~\ref{alg:strawman-full}), as this phase dominates end-to-end proof generation time.")
    print(r"We also execute $\Pi_\textrm{PermVanish}$ and report its timing for completeness.")
    print(f"Unless stated otherwise, we use $n = {SCALABILITY_PARTIES}$ parties and set the query size to $M = {tex_ratio(new_ratio)}$, where $N$ is the table size.")
    print(f"For the scalability experiment we vary $N \\in {size_list(new_sizes)}$ for \\MPlookup and $N \\in {size_list(naive_sizes)}$ for the strawman.")
    print(f"For the multi-party experiment we fix $N = {tex_n(party_n)}$ and vary $n \\in \\{{{', '.join(str(p) for p in parties_sorted)}\\}}$.")

    # --- RQ1 ---
    print()
    print(r"\textbf{RQ1: Scalability and comparison.}")
    print(r"Figures~\ref{fig:rq1_total_time} and~\ref{fig:rq1_comm} plot preprocessing time, bytes sent, and broadcast count.")
    print(r"The strawman follows the predicted $O(N^2)$ curve, while \MPlookup's observed doubling ratios are close to linear over the tested range.")
    print()
    print(r"At small $N$, the strawman is faster due to lower constant factors. Each iteration uses one $\mathcal{F}_\textrm{Eq}$ call, whereas \MPlookup invokes sorting networks and polynomial arithmetic.")
    print(f"The crossover lies between $N = {tex_n(crossover_n_below)}$ and $N = {tex_n(crossover_n)}$.")
    print(f"At $N = {tex_n(n_hi)}$, \\MPlookup achieves a ${speedup_hi:.2f}\\times$ speedup, completing in ${tex_n(t_new_hi)}$\\,s vs.\\ ${tex_n(t_naive_hi)}$\\,s for the strawman.")
    print(f"From $N = {tex_n(trio_ns[0])}$ to ${tex_n(trio_ns[-1])}$, the observed doubling ratios are {obs_str}, compared to the theoretical $O(N\\log^2 N)$ ratios of {theo_str}.")
    print()
    print(f"Figure~\\ref{{fig:rq1_speedup}} shows speedups in preprocessing time, bytes sent, and broadcast count. The time speedup turns decisively in \\MPlookup's favour beyond $N = {tex_n(crossover_n)}$.")
    print(f"Communication savings are more pronounced. At $N = {tex_n(n_hi)}$, \\MPlookup sends approximately ${bs_speedup:.0f}\\times$ fewer bytes, as the strawman's $O(N^2)$ equality comparisons each require a full MPC round.")
    print(f"Broadcast count is higher for \\MPlookup only at $N = {tex_n(bc_crossover_below)}$, due to sorting and subproduct tree rounds; from $N = {tex_n(bc_crossover_n)}$ onward \\MPlookup also issues fewer broadcasts.")

    # --- RQ2 ---
    print()
    print(r"\textbf{RQ2: Step-level performance breakdown.}")
    # Same order as the paper: phases first (panel d), steps second (panel g).
    # The figure is read (d) before (g), so a paragraph that reversed them would
    # send the reader back up the page.
    print(r"Figure~\ref{fig:rq2_mplookup_phases} breaks down the four phases of \MPlookup, setup, proof generation consisting of $\Pi_\textrm{Preprocess}$ and $\Pi_\textrm{PermVanish}$, and finally verification.")
    print(f"$\\Pi_\\textrm{{Preprocess}}$ overwhelmingly dominates the total cost across all tested sizes, while verification remains essentially constant at approximately ${ver_ms_approx}$\\,ms regardless of $N$.")
    print(f"$\\Pi_\\textrm{{PermVanish}}$ and setup both grow with $N$ but contribute less than ${pct_upper:.1f}\\%$ of total time, confirming that preprocessing is the decisive bottleneck.")
    print()
    print(r"Figure~\ref{fig:rq2_steps} shows the time cost contributed by each of the nine steps of Algorithm~\ref{alg:mplookup-preprocessing} across all tested table sizes.")

    # --- RQ3 ---
    print()
    print(r"\textbf{RQ3: Scaling with the number of parties.}")
    print(f"Figure~\\ref{{fig:rq3_parties}} shows the preprocessing time at $N={tex_n(party_n)}$ and $M={tex_n(party_m)}$, which grows from ${tex_n(base_time)}$\\,s for 2 parties to {oxford_join(times_tex)}---slowdowns of {oxford_join(slowdowns_tex)} over the 2-party baseline.")
    print(r"Figure~\ref{fig:rq3_comm} shows bytes sent scales similarly to time, while broadcast count remains constant across party counts, consistent with the communication structure.")
    print()
    # Wording matches mplookup.tex, where this sentence was tightened by hand.
    # The numbers still come from the data.
    print(f"In Figure~\\ref{{fig:rq3_phases}}, setup and verification time remain stable as they are single-user protocols, while $\\Pi_\\textrm{{PermVanish}}$ grows from ${pg_2:.1f}$\\,s at 2 parties to ${pg_16:.1f}$\\,s at 16 parties since more parties mean more KZG commitment rounds.")

    print()


# ---------------------------------------------------------------------------
# --check: verify that the paper still quotes what the data says
# ---------------------------------------------------------------------------

def _normalize(s):
    """Collapse all whitespace so ordinary re-wrapping in the .tex is ignored."""
    return re.sub(r'\s+', ' ', s).strip()


def _gather_output():
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        main()
    return buf.getvalue()


def check_against_paper(path):
    """Return the generated passages that the paper does not contain verbatim."""
    tex = _normalize(path.read_text(encoding='utf-8'))
    problems = []

    for line in _gather_output().splitlines():
        line = _normalize(line)
        if not line or line.startswith(r'\subsection'):
            continue
        if line not in tex:
            problems.append(line)

    return problems


def paper_tex_path(argv):
    """Locate the paper .tex, or return None if it cannot be found.

    The paper is a separate checkout that is not part of the artifacts
    repository this script ships in, so its location is never derivable from
    __file__: it has to be supplied by the caller, either as `--check PATH`
    or through $MPLOOKUP_TEX.  The sibling path is only a convenience default
    for the authors' working copy.
    """
    explicit = [a for a in argv if a != '--check']
    if explicit:
        return Path(explicit[0])
    if os.environ.get('MPLOOKUP_TEX'):
        return Path(os.environ['MPLOOKUP_TEX'])
    return DEFAULT_PAPER_TEX


if __name__ == '__main__':
    if '--check' in sys.argv[1:]:
        path = paper_tex_path(sys.argv[1:])
        if path is None or not path.exists():
            print(f'error: no such file: {path}', file=sys.stderr)
            print('pass the path to the paper, e.g.\n'
                  '  compute-evaluation-data.py --check ../../mplookup-paper/mplookup.tex\n'
                  'or set $MPLOOKUP_TEX to it; this script cannot know where the '
                  'paper lives relative to the artifacts repository.',
                  file=sys.stderr)
            sys.exit(2)
        problems = check_against_paper(path)
        if problems:
            print(f'{len(problems)} passage(s) of {path.name} differ from the data:',
                  file=sys.stderr)
            for p in problems:
                print(f'  - {p}', file=sys.stderr)
            sys.exit(1)
        print(f'{path.name} is in sync with the evaluation data.')
    else:
        main()
