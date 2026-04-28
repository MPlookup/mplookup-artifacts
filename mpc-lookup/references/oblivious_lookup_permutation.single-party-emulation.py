# Secure Oblivious Permutation Algorithm Simulation in Python
# This script simulates the MPC algorithm in plaintext over a finite field F_p.
# All operations are performed modulo p, a large prime.
# The number of elements in t (N) and f (M) are configurable.
# Assumptions: M <= N, t has distinct elements, all elements are integers between 0 and p-1.
# Optimized to O(N log^2 N) time complexity using subproduct tree and iterative multipoint evaluation.

p: int = 998244353  # Prime modulus for the field, NTT-friendly

def add(a: int, b: int, mod: int) -> int:
    return (a + b) % mod

def sub(a: int, b: int, mod: int) -> int:
    return (a - b) % mod

def mul(a: int, b: int, mod: int) -> int:
    return (a * b) % mod

def inv(a: int, mod: int) -> int:
    return pow(a, mod - 2, mod)

def ntt(f: list[int], inverse: bool = False, mod: int = p) -> list[int]:
    n: int = len(f)
    if n == 1:
        return f
    # Bit reversal permutation
    logn: int = n.bit_length() - 1
    for i in range(n):
        j: int = 0
        temp: int = i
        for k in range(logn):
            j = (j << 1) | (temp & 1)
            temp >>= 1
        if i < j:
            f[i], f[j] = f[j], f[i]
    # Primitive root and computation
    g: int = 3
    ig: int = inv(g, mod)
    logn: int = n.bit_length() - 1
    len_: int = 2
    for stage in range(1, logn + 1):
        wlen: int = pow(ig if inverse else g, (mod - 1) // len_, mod)
        for i in range(0, n, len_):
            w: int = 1
            for jj in range(len_ // 2):
                u: int = f[i + jj]
                v: int = mul(f[i + jj + len_ // 2], w, mod)
                f[i + jj] = add(u, v, mod)
                f[i + jj + len_ // 2] = sub(u, v, mod)
                w = mul(w, wlen, mod)
        len_ *= 2
    if inverse:
        inv_n: int = inv(n, mod)
        for i in range(n):
            f[i] = mul(f[i], inv_n, mod)
    return f

def poly_mul(p1: list[int], p2: list[int], mod: int) -> list[int]:
    if not p1 or not p2:
        return []
    deg: int = len(p1) + len(p2) - 1
    sz: int = 1 << ((deg - 1).bit_length()) if deg > 0 else 0
    pp1: list[int] = p1 + [0] * (sz - len(p1))
    pp2: list[int] = p2 + [0] * (sz - len(p2))
    nt1: list[int] = ntt(pp1, mod=mod)
    nt2: list[int] = ntt(pp2, mod=mod)
    for i in range(sz):
        nt1[i] = mul(nt1[i], nt2[i], mod)
    res: list[int] = ntt(nt1, inverse=True, mod=mod)
    return res[:deg]

def poly_inv(f: list[int], deg: int, mod: int) -> list[int]:
    if deg <= 0:
        return []
    target_pow: int = 1 << ((deg - 1).bit_length()) if deg > 0 else 0
    logd: int = target_pow.bit_length() - 1
    g: list[int] = [inv(f[0], mod)]
    cur_deg: int = 1
    for level in range(logd):
        next_deg: int = cur_deg * 2
        ff: list[int] = f[:next_deg] + [0] * max(0, next_deg - len(f))
        tmp: list[int] = poly_mul(g, ff, mod)[:next_deg]
        tmp = poly_mul(tmp, g, mod)[:next_deg]
        new_g: list[int] = [0] * next_deg
        for i in range(cur_deg):
            new_g[i] = mul(2, g[i], mod)
        for i in range(next_deg):
            new_g[i] = sub(new_g[i], tmp[i], mod)
        g = new_g
        cur_deg = next_deg
    return g[:deg]

def poly_divmod(num: list[int], den: list[int], mod: int) -> tuple[list[int], list[int]]:
    m: int = len(num) - 1
    n: int = len(den) - 1
    if m < n:
        return [0], num
    deg_q: int = m - n + 1
    rev_num: list[int] = num[::-1]
    rev_den: list[int] = den[::-1]
    inv_rev_den: list[int] = poly_inv(rev_den, deg_q, mod)
    q_rev: list[int] = poly_mul(rev_num[:deg_q], inv_rev_den, mod)[:deg_q]
    q: list[int] = q_rev[::-1]
    # Compute remainder
    q_den: list[int] = poly_mul(q, den, mod)
    r_len: int = len(num)
    pad_len: int = max(0, r_len - len(q_den))
    q_den += [0] * pad_len
    q_den = q_den[:r_len]
    r: list[int] = [sub(num[i], q_den[i], mod) for i in range(r_len)]
    return q, r

def build_subproduct_tree(points: list[int], mod: int) -> list[list[list[int]]]:
    if not points:
        return [[[1]]]
    k: int = len(points)
    pow2: int = 1 << ((k - 1).bit_length()) if k > 0 else 0
    polys: list[list[int]] = [[sub(0, p, mod), 1] for p in points] + [[1]] * (pow2 - k)
    tree: list[list[list[int]]] = [polys[:]]
    while len(polys) > 1:
        new_polys: list[list[int]] = []
        for i in range(0, len(polys), 2):
            p1: list[int] = polys[i]
            p2: list[int] = polys[i + 1] if i + 1 < len(polys) else [1]
            product: list[int] = poly_mul(p1, p2, mod)
            new_polys.append(product)
        tree.append(new_polys)
        polys = new_polys
    tree = tree[::-1]  # tree[0]: root level, tree[-1]: leaves level
    return tree

def multipoint_eval(P: list[int], tree: list[list[list[int]]], mod: int) -> list[int]:
    if not tree or not tree[0]:
        return []
    num_levels: int = len(tree)
    reduced_P: list[list[list[int]]] = [[[] for _ in range(len(tree[lvl]))] for lvl in range(num_levels)]
    reduced_P[0][0] = P[:]
    for lvl in range(num_levels - 1):
        for nd in range(len(tree[lvl])):
            curr_P: list[int] = reduced_P[lvl][nd]
            left_idx: int = 2 * nd
            right_idx: int = 2 * nd + 1
            left_sub: list[int] = tree[lvl + 1][left_idx]
            right_sub: list[int] = tree[lvl + 1][right_idx]
            left_deg: int = len(left_sub) - 1
            right_deg: int = len(right_sub) - 1
            if left_deg > 0:
                rem: list[int] = poly_divmod(curr_P, left_sub, mod)[1]
                # Truncate to degree < deg(left_sub).  By the polynomial
                # remainder theorem, coefficients at indices >= left_deg are
                # exactly zero; removing them restores the O(n/2^d) size
                # invariant at depth d and keeps the overall cost O(n log^2 n).
                # This is safe even when t is secret-shared: deg(left_sub) is
                # determined solely by the public tree structure (level index
                # and table size n), so the truncation reveals nothing secret.
                reduced_P[lvl + 1][left_idx] = rem[:left_deg]
            else:
                reduced_P[lvl + 1][left_idx] = [0]  # padding node
            if right_deg > 0:
                rem = poly_divmod(curr_P, right_sub, mod)[1]
                reduced_P[lvl + 1][right_idx] = rem[:right_deg]  # same truncation
            else:
                reduced_P[lvl + 1][right_idx] = [0]  # padding node
    # Collect evaluations from real leaves
    leaf_level: int = num_levels - 1
    evals: list[int] = []
    for nd in range(len(tree[leaf_level])):
        sub: list[int] = tree[leaf_level][nd]
        if len(sub) > 1:  # Real leaf (degree 1)
            P_leaf: list[int] = reduced_P[leaf_level][nd]
            evals.append(P_leaf[0])
    return evals

def secure_oblivious_permutation(N: int, M: int, t: list[int], f: list[int], mod: int) -> tuple[list[int], list[int]]:
    # Step 1: Sort the query vector
    f_prime: list[int] = sorted(f)

    # Pad f_prime to size N with 0
    f_prime += [0] * (N - M)
    print("After Step 1 (sorted f_prime):", f_prime)

    # Step 2: Compute first-occurrence indicators
    f_double: list[int] = [0] * N
    f_double[0] = 1
    for i in range(1, N):
        is_before_M: int = 1 if i < M else 0
        eq: int = 1 if f_prime[i] == f_prime[i - 1] else 0
        f_double[i] = (1 - eq) * is_before_M
    print("After Step 2 (f_double):", f_double)

    # Step 3: Construct membership polynomial over f[0..M-1]
    f_tree: list[list[list[int]]] = build_subproduct_tree(f_prime[:M], mod)
    P_coeffs: list[int] = f_tree[0][0]
    print("After Step 3 (P_coeffs):", P_coeffs)

    # Step 4: Evaluate polynomial at table elements
    t_tree: list[list[list[int]]] = build_subproduct_tree(t, mod)
    p_values: list[int] = multipoint_eval(P_coeffs, t_tree, mod)
    print("After Step 4 (p_values):", p_values)

    # Step 5: Compute unused indicators (1 if not present in f, i.e., unused)
    unused: list[int] = [1 if p_values[j] != 0 else 0 for j in range(N)]
    print("After Step 5 (unused):", unused)

    # Step 6: Compact unused elements
    unusedValueRec: list[tuple[int, int, int]] = [(j, t[j], unused[j]) for j in range(N)]
    sorted_unused: list[tuple[int, int, int]] = sorted(unusedValueRec, key=lambda rec: (-rec[2], rec[0]))
    print("After Step 6 (unused_rec):", sorted_unused)

    # Step 7: Compact unfilled positions
    unfilledPositionRec: list[tuple[int, int]] = [(i, 1 - f_double[i]) for i in range(N)]
    sorted_unfilled: list[tuple[int, int]] = sorted(unfilledPositionRec, key=lambda rec: (-rec[1], rec[0]))
    print("After Step 7 (unfilled_rec):", sorted_unfilled)

    # Step 8: Create write records 1
    writeRec1: list[tuple[int, int, int]] = []
    for m in range(N):
        pos: int = sorted_unfilled[m][0]
        val: int = sorted_unused[m][1]
        valid: int = sorted_unfilled[m][1]
        writeRec1.append((pos, val, valid))
    print("After Step 8 (write_rec1):", writeRec1)

    # Step 9: Create write records 2
    writeRec2: list[tuple[int, int, int]] = []
    for k in range(N):
        pos: int = k
        val: int = f_prime[k]
        valid: int = f_double[k]
        writeRec2.append((pos, val, valid))
    print("After Step 9 (write_rec2):", writeRec2)

    # Step 10: Combine and sort write records
    allWrites: list[tuple[int, int, int]] = writeRec1 + writeRec2
    sorted_all: list[tuple[int, int, int]] = sorted(allWrites, key=lambda rec: (-rec[2], rec[0]))
    print("After Step 10 (all_writes):", sorted_all)

    # Step 11: Construct output permutation
    t_prime: list[int] = [sorted_all[r][1] for r in range(N)]
    print("After Step 11 (t_prime):", t_prime)

    # Return t' and f' (sorted query, size M)
    return t_prime, f_prime[:M]

# Example usage
if __name__ == "__main__":
    print("==== Example 1 ====")
    N: int = 8
    M: int = 4
    t: list[int] = [10, 20, 30, 0, 40, 50, 60, 70]
    f: list[int] = [30, 0, 30, 0]
    print("Input t:", t)
    print("Input f:", f)
    t_prime, f_prime = secure_oblivious_permutation(N, M, t, f, p)
    print("Permuted t':", t_prime)
    print("Sorted f':", f_prime)
    print()

    print("==== Example 2 ====")
    N: int = 8
    M: int = 4
    t: list[int] = [10, 20, 70, 0, 40, 50, 60, 30]
    f: list[int] = [30, 0, 30, 0]
    print("Input t:", t)
    print("Input f:", f)
    t_prime, f_prime = secure_oblivious_permutation(N, M, t, f, p)
    print("Permuted t':", t_prime)
    print("Sorted f':", f_prime)
    print()