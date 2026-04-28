"""
compute_paper_stats.py
======================
Computes and prints every in-text numerical value used in the evaluation
section (§ "Implementation and Evaluation") of the paper.

Run from any directory:
    python ref/experiments/result-analyze/compute_paper_stats.py

Output is grouped by research question and labelled to match the exact
sentence in the TeX file where each number appears.
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


def fmt_int(x):
    """Format integer with comma thousands separator."""
    return f"{int(round(x)):,}"


def fmt_1f(x):
    return f"{x:.1f}"


def fmt_2f(x):
    return f"{x:.2f}"


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

    common_ns = sorted(set(new_results) & set(naive_results))
    new_ns    = sorted(new_results.keys())

    print("=" * 70)
    print("PAPER STATS — all values used in the evaluation section")
    print("=" * 70)

    # -----------------------------------------------------------------------
    # RQ1 — Crossover between Strawman and MPlookup
    # -----------------------------------------------------------------------
    print("\n--- RQ1: Scalability and Comparison against the Strawman ---\n")

    # Strawman-to-MPlookup time ratio at each common n (ratio < 1 → Strawman faster)
    print("Time ratio (Strawman / MPlookup) at each n:")
    crossover_n = None
    for n in common_ns:
        ratio = naive_results[n]['naive_mpc_perm_time'] / new_results[n]['mpc_perm_time']
        marker = " ← crossover" if crossover_n is None and ratio > 1.0 else ""
        if crossover_n is None and ratio > 1.0:
            crossover_n = n
        print(f"  n={n:6d}:  Strawman/MPlookup = {ratio:.4f}x{marker}")

    # Specific crossover sentences
    n_below = common_ns[common_ns.index(crossover_n) - 1] if crossover_n else None
    if n_below:
        r_below = naive_results[n_below]['naive_mpc_perm_time'] / new_results[n_below]['mpc_perm_time']
        r_cross = naive_results[crossover_n]['naive_mpc_perm_time'] / new_results[crossover_n]['mpc_perm_time']
        print(f"\nTeX sentence: crossover between n={n_below:,} and n={crossover_n:,}")
        print(f"  At n={n_below:,}: Strawman is {r_below:.2f}× the time of MPlookup  (TEX: ${{0:.2f}}\\times$)")
        print(f"  At n={crossover_n:,}: MPlookup is {r_cross:.2f}× faster             (TEX: ${{0:.2f}}\\times$)")

    # At n=8192
    n_highlight = 8192
    if n_highlight in new_results and n_highlight in naive_results:
        t_new   = new_results[n_highlight]['mpc_perm_time']
        t_naive = naive_results[n_highlight]['naive_mpc_perm_time']
        sp      = t_naive / t_new
        print(f"\nTeX sentence: At n={n_highlight:,}")
        print(f"  MPlookup time : {fmt_int(t_new)} s    (TEX: ${fmt_int(t_new)}$\\,s)")
        print(f"  Strawman time : {fmt_int(t_naive)} s  (TEX: ${fmt_int(t_naive)}$\\,s)")
        print(f"  Speedup       : {sp:.2f}×            (TEX: ${fmt_2f(sp)}\\times$)")

    # Doubling ratios for MPlookup (consecutive doublings across all available n)
    mplookup_sorted_ns = [n for n in new_ns if n >= 1024]
    print(f"\nMPlookup doubling ratios (observed vs theoretical) for n >= 1024:")
    doubling_data = []
    for i in range(len(mplookup_sorted_ns) - 1):
        n0 = mplookup_sorted_ns[i]
        n1 = mplookup_sorted_ns[i + 1]
        if n1 != 2 * n0:
            continue
        obs  = new_results[n1]['mpc_perm_time'] / new_results[n0]['mpc_perm_time']
        theo = theoretical_doubling_ratio(n0)
        doubling_data.append((n0, n1, obs, theo))
        print(f"  n={n0:6d}→{n1:6d}: observed={obs:.2f}×  theoretical={theo:.2f}×")

    # Pick the 3 doublings that the TeX sentence refers to.
    # The sentence says "from n=4096 to 16384" with 3 values — but that range
    # only has 2 doublings.  The data most naturally covering 3 doublings is
    # n=2048→4096→8192→16384.  We report that range (and note it in the output).
    trio_start = 2048
    trio = [(n0, n1, obs, theo) for (n0, n1, obs, theo) in doubling_data
            if n0 >= trio_start and n1 <= 16384]
    if len(trio) == 3:
        obs_vals  = [f"{obs:.2f}\\times" for _, _, obs, _ in trio]
        theo_vals = [f"{theo:.2f}\\times" for _, _, _, theo in trio]
        n_start = trio[0][0]
        n_end   = trio[-1][1]
        print(f"\nTeX sentence (3 doublings, n={n_start:,}→{n_end:,}):")
        print(f"  Observed  : {', '.join(obs_vals)}")
        print(f"  Theoretical: {', '.join(theo_vals)}")
        print(f"  (Update TeX range to 'from $n = {fmt_int(n_start)}$ to ${fmt_int(n_end)}$')")
    elif len(trio) != 3:
        print(f"\nWARNING: found {len(trio)} doublings in n={trio_start:,}–16384 (expected 3); adjust manually.")
        for n0, n1, obs, theo in trio:
            print(f"  n={n0}→{n1}: obs={obs:.2f} theo={theo:.2f}")

    # Bytes-sent speedup at n=8192
    if n_highlight in new_results and n_highlight in naive_results:
        bs_speedup = naive_results[n_highlight]['bytes_sent'] / new_results[n_highlight]['bytes_sent']
        print(f"\nTeX sentence: bytes speedup at n={n_highlight:,}")
        print(f"  Strawman bytes_sent : {naive_results[n_highlight]['bytes_sent']:,}")
        print(f"  MPlookup bytes_sent : {new_results[n_highlight]['bytes_sent']:,}")
        print(f"  Speedup             : {bs_speedup:.1f}×   (TEX: approximately ${bs_speedup:.0f}\\times$)")

    # Broadcast crossover (first n where MPlookup broadcasts < Strawman broadcasts)
    print(f"\nBroadcast crossover (first n where MPlookup < Strawman):")
    for n in common_ns:
        nb = new_results[n]['broadcasts']
        sb = naive_results[n]['broadcasts']
        marker = " ← crossover" if nb < sb else ""
        print(f"  n={n:6d}: MPlookup={nb:>14,}  Strawman={sb:>14,}{marker}")

    # PermVanish share of total time (max across all n)
    print(f"\nPermVanish (proof_gen) as % of total (mpc_perm_time + proof_gen), all n:")
    max_pct = 0.0
    for n in new_ns:
        d = new_results[n]
        if 'proof_gen' in d:
            total = d['mpc_perm_time'] + d['proof_gen']
            p = pct(d['proof_gen'], total)
            if p > max_pct:
                max_pct = p
            print(f"  n={n:6d}: {p:.4f}%")
    print(f"  Max across all n: {max_pct:.4f}% → TEX: 'less than {math.ceil(max_pct * 10) / 10:.1f}%'")

    # -----------------------------------------------------------------------
    # RQ2 — Step breakdown at n=1024 and n=16384
    # -----------------------------------------------------------------------
    print("\n--- RQ2: Step-Level Performance Breakdown ---\n")

    sort_steps = {1, 6, 7, 8}  # after merging

    for n_rq2 in [1024, 16384]:
        if n_rq2 not in new_results or 'steps' not in new_results[n_rq2]:
            print(f"  n={n_rq2}: step data not available")
            continue
        merged = merge_steps(new_results[n_rq2]['steps'])
        total  = new_results[n_rq2]['mpc_perm_time']
        sort_t = sum(merged[s]['time'] for s in sort_steps if s in merged)
        step4_t = merged.get(4, {}).get('time', 0.0)
        other_t = sum(merged[s]['time'] for s in [2, 3, 5, 9] if s in merged)

        print(f"n = {n_rq2:,} (total mpc_perm_time = {total:.3f} s):")
        for snum in sorted(merged.keys()):
            t = merged[snum]['time']
            print(f"  Step {snum}: {t:>12.3f} s  = {pct(t, total):.2f}%  [{merged[snum]['name']}]")
        print(f"  Sort steps (1,6,7,8) combined : {sort_t:.3f} s = {pct(sort_t, total):.2f}%")
        print(f"  Step 4 (poly eval)             : {step4_t:.3f} s = {pct(step4_t, total):.2f}%")
        print(f"  Steps 2,3,5,9 combined         : {other_t:.3f} s = {pct(other_t, total):.2f}%")
        print()

    # Step 8 at n=1024
    if 1024 in new_results and 'steps' in new_results[1024]:
        merged1024 = merge_steps(new_results[1024]['steps'])
        total1024  = new_results[1024]['mpc_perm_time']
        s8_t = merged1024.get(8, {}).get('time', 0.0)
        print(f"TeX sentence: Step 8 at n=1,024")
        print(f"  Step 8 time : {s8_t:.3f} s")
        print(f"  Step 8 share: {pct(s8_t, total1024):.1f}%   (TEX: approximately ${pct(s8_t, total1024):.0f}\\%$)")

    # Sort steps min share across all n with step data
    print(f"\nSort steps (1,6,7,8) share across all n:")
    for n in sorted(new_results.keys()):
        if 'steps' not in new_results[n]:
            continue
        merged = merge_steps(new_results[n]['steps'])
        total  = new_results[n]['mpc_perm_time']
        sort_t = sum(merged[s]['time'] for s in sort_steps if s in merged)
        print(f"  n={n:6d}: {pct(sort_t, total):.2f}%")

    # Step 4 range
    print(f"\nStep 4 (poly eval) share across all n:")
    for n in sorted(new_results.keys()):
        if 'steps' not in new_results[n]:
            continue
        merged = merge_steps(new_results[n]['steps'])
        total  = new_results[n]['mpc_perm_time']
        s4_t = merged.get(4, {}).get('time', 0.0)
        print(f"  n={n:6d}: {pct(s4_t, total):.2f}%")

    # Steps 2,3,5,9 max share
    print(f"\nSteps 2,3,5,9 combined share across all n:")
    max_minor = 0.0
    for n in sorted(new_results.keys()):
        if 'steps' not in new_results[n]:
            continue
        merged = merge_steps(new_results[n]['steps'])
        total  = new_results[n]['mpc_perm_time']
        other_t = sum(merged[s]['time'] for s in [2, 3, 5, 9] if s in merged)
        p = pct(other_t, total)
        if p > max_minor:
            max_minor = p
        print(f"  n={n:6d}: {p:.4f}%")
    print(f"  Max: {max_minor:.4f}% → TEX: 'less than {math.ceil(max_minor * 10) / 10:.1f}%'")

    # -----------------------------------------------------------------------
    # RQ3 — Multi-party scaling at n=1024
    # -----------------------------------------------------------------------
    print("\n--- RQ3: Scaling with the Number of Parties (n=1024) ---\n")

    parties_sorted = sorted(party_results.keys())
    base_time = party_results[2]['mpc_perm_time']

    print("Preprocessing time and speedup vs 2-party baseline:")
    for p in parties_sorted:
        t = party_results[p]['mpc_perm_time']
        sp = t / base_time
        print(f"  {p:2d} parties: {fmt_int(t)} s  (speedup over 2p: {sp:.2f}×)")

    print(f"\nTeX sentence: times and speedups")
    times_fmt = [f"${fmt_int(party_results[p]['mpc_perm_time'])}$\\,s" for p in parties_sorted]
    sp_fmt    = [f"${party_results[p]['mpc_perm_time'] / base_time:.1f}\\times$"
                 for p in parties_sorted if p != 2]
    print(f"  Times    : {', '.join(times_fmt)}")
    print(f"  Speedups : {', '.join(sp_fmt)}")

    # PermVanish vs parties
    print(f"\nPermVanish (proof_gen) across party counts:")
    for p in parties_sorted:
        pg = party_results[p].get('proof_gen', float('nan'))
        print(f"  {p:2d} parties: {pg:.4f} s  ≈ {fmt_1f(pg)} s")
    p2_pg  = party_results[2].get('proof_gen', float('nan'))
    p16_pg = party_results[16].get('proof_gen', float('nan')) if 16 in party_results else float('nan')
    print(f"\nTeX sentence: PermVanish from {fmt_1f(p2_pg)} s at 2 parties to {fmt_1f(p16_pg)} s at 16 parties")

    print("\n" + "=" * 70)
    print("END OF STATS")
    print("=" * 70)


if __name__ == '__main__':
    main()
