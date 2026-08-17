//! Oblivious sorting operations for MPC fields.
//!
//! This module provides secure bitonic sorting algorithms that work
//! with MpcField values, maintaining obliviousness in MPC contexts.

use ark_ff::One;
use mpc_algebra::honest_but_curious::MpcField;
use mpc_algebra::{BooleanOps, DummyBooleanBeaverSource};
use mpc_net::{MpcMultiNet as Net, MpcNet};

type F = ark_bls12_377::Fr;
type MF = MpcField<F>;

/// Bitonic sorter for MpcField values using secure comparison.
/// 
/// This is an oblivious sorting algorithm that sorts in ascending order.
/// The input vector length must be a power of 2.
pub fn bitonic_sort(mut vec: Vec<MF>) -> Vec<MF> {
    let n = vec.len();

    if n < 2 || (n & (n - 1)) != 0 {
        panic!(
            "bitonic_sort: vector length must be a power of 2 and at least 2, got {}",
            n
        );
    }

    let mut k = 2;
    while k <= n {
        let mut j = k / 2;
        while j > 0 {
            let mut pairs = Vec::new();
            let mut indices = Vec::new();
            
            for i in 0..n {
                let i_xor_j = i ^ j;
                if i_xor_j > i {
                    let swap_if_greater = (i & k) == 0;
                    pairs.push((vec[i], vec[i_xor_j]));
                    indices.push((i, i_xor_j, swap_if_greater));
                }
            }
            
            let is_less_results = MF::batch_secure_cmp(&pairs);
            
            for (idx, &(i, i_xor_j, swap_if_greater)) in indices.iter().enumerate() {
                let x = vec[i];
                let y = vec[i_xor_j];
                let is_less = is_less_results[idx];
                
                let swap_condition = if swap_if_greater {
                    MF::Public(F::one()) - is_less
                } else {
                    is_less
                };
                
                let not_swap_condition = MF::Public(F::one()) - swap_condition;
                
                let new_x = swap_condition * y + not_swap_condition * x;
                let new_y = swap_condition * x + not_swap_condition * y;
                
                vec[i] = new_x;
                vec[i_xor_j] = new_y;
            }
            
            j /= 2;
        }
        k *= 2;
    }

    vec
}

/// Oblivious sort for (key, value1, value2) tuples.
///
/// Sorts by the key field in ascending order, moving associated values with it.
/// The input vector length must be a power of 2.
pub fn bitonic_sort_by_key(mut vec: Vec<(MF, MF, MF)>) -> Vec<(MF, MF, MF)> {
    let n = vec.len();

    if n < 2 || (n & (n - 1)) != 0 {
        panic!(
            "bitonic_sort_by_key: vector length must be a power of 2 and at least 2, got {}",
            n
        );
    }

    let mut k = 2;
    while k <= n {
        let mut j = k / 2;
        while j > 0 {
            let mut key_pairs = Vec::new();
            let mut indices = Vec::new();
            
            for i in 0..n {
                let i_xor_j = i ^ j;
                if i_xor_j > i {
                    let swap_if_greater = (i & k) == 0;
                    key_pairs.push((vec[i].0, vec[i_xor_j].0));
                    indices.push((i, i_xor_j, swap_if_greater));
                }
            }
            
            let is_less_results = MF::batch_secure_cmp(&key_pairs);
            
            for (idx, &(i, i_xor_j, swap_if_greater)) in indices.iter().enumerate() {
                let (key_i, val1_i, val2_i) = vec[i];
                let (key_j, val1_j, val2_j) = vec[i_xor_j];
                let is_less = is_less_results[idx];
                
                let swap_condition = if swap_if_greater {
                    MF::Public(F::one()) - is_less
                } else {
                    is_less
                };
                
                let not_swap_condition = MF::Public(F::one()) - swap_condition;
                
                let new_key_i = swap_condition * key_j + not_swap_condition * key_i;
                let new_val1_i = swap_condition * val1_j + not_swap_condition * val1_i;
                let new_val2_i = swap_condition * val2_j + not_swap_condition * val2_i;
                
                let new_key_j = swap_condition * key_i + not_swap_condition * key_j;
                let new_val1_j = swap_condition * val1_i + not_swap_condition * val1_j;
                let new_val2_j = swap_condition * val2_i + not_swap_condition * val2_j;
                
                vec[i] = (new_key_i, new_val1_i, new_val2_i);
                vec[i_xor_j] = (new_key_j, new_val1_j, new_val2_j);
            }
            
            j /= 2;
        }
        k *= 2;
    }

    vec
}

/// Radix sort for MpcField values using precomputed bits (bitonic network).
///
/// This is an optimized oblivious sorting algorithm that reduces the number of
/// sequential communication rounds from O(N log²N × B) to O(log²N × B) by
/// batch-decomposing all N elements once, then using precomputed bits for all
/// comparisons, batching all N/2 comparisons at each level simultaneously.
///
/// Performance:
/// - 1 batch broadcast for all N bit decompositions (vs. N broadcasts in bitonic)
/// - O(B) rounds per comparison level (vs. O(B × N/2) rounds in bitonic)
/// - Total: O(B × log²N) sequential rounds (vs. O(B × N/2 × log²N) in bitonic)
///
/// The input vector length must be a power of 2.
pub fn radix_sort(mut vec: Vec<MF>, bit_length: usize) -> Vec<MF> {
    let n = vec.len();

    if n < 2 || (n & (n - 1)) != 0 {
        panic!(
            "radix_sort: vector length must be a power of 2 and at least 2, got {}",
            n
        );
    }

    // Fast path: when all values are public, sort locally without MPC operations.
    // This preserves the Public variant through the output (the MPC path below
    // always returns Shared due to batch_b2a_bool).
    if vec.iter().all(|v| matches!(v, MF::Public(_))) {
        let mut plain: Vec<F> = vec.iter().map(|v| match v {
            MF::Public(x) => *x,
            _ => unreachable!(),
        }).collect();
        plain.sort();
        return plain.into_iter().map(MF::Public).collect();
    }

    // Step 1: Batch precompute all bits in one broadcast round, then truncate to
    // bit_length (values are in [0, 2^bit_length) so high bits are always 0).
    let full_bits = MF::batch_bit_decompose_bool(&vec);
    let mut bits: Vec<Vec<bool>> = full_bits
        .into_iter()
        .map(|mut b| {
            b.truncate(bit_length);
            b
        })
        .collect();

    let mut k = 2;
    while k <= n {
        let mut j = k / 2;
        while j > 0 {
            // Collect all pairs for this comparison level.
            let mut pair_indices: Vec<(usize, usize, bool)> = Vec::new();
            for i in 0..n {
                let ij = i ^ j;
                if ij > i {
                    let swap_if_greater = (i & k) == 0;
                    pair_indices.push((i, ij, swap_if_greater));
                }
            }
            let num_pairs = pair_indices.len();

            // Step 2: Batch compare all pairs using precomputed bits.
            let left_bits: Vec<Vec<bool>> =
                pair_indices.iter().map(|&(i, _, _)| bits[i].clone()).collect();
            let right_bits: Vec<Vec<bool>> =
                pair_indices.iter().map(|&(_, ij, _)| bits[ij].clone()).collect();
            let is_less_results = MF::batch_compare_bits(&left_bits, &right_bits);

            // Step 3: Compute boolean swap conditions (local, no communication).
            // swap_cond = is_less XOR swap_if_greater
            //   swap_if_greater=true  → swap when left >= right → swap_cond = NOT is_less
            //   swap_if_greater=false → swap when left < right  → swap_cond = is_less
            //
            // IMPORTANT: is_less is an XOR-shared boolean. XORing a public constant
            // (swig) must be applied to only ONE party's share, otherwise both parties
            // flip and the XOR of the two shares cancels out (true ^ true = false).
            // We apply the XOR to the king's share only.
            let swap_cond_bools: Vec<bool> = pair_indices
                .iter()
                .zip(is_less_results.iter())
                .map(|(&(_, _, swig), &is_less)| {
                    if swig { is_less ^ Net::am_king() } else { is_less }
                })
                .collect();

            // Step 4: Batch B2A conversion of swap conditions.
            let swap_cond_arith: Vec<MF> = MF::batch_b2a_bool(&swap_cond_bools);

            // Step 5: Compute arithmetic deltas (local, no communication).
            // delta[p] = vec[j] - vec[i]; new_vec[i] += swap_cond * delta; new_vec[j] -= ...
            let deltas: Vec<MF> = pair_indices
                .iter()
                .map(|&(i, ij, _)| vec[ij] - vec[i])
                .collect();

            // Step 6: Batch multiply swap conditions × deltas (2 broadcasts).
            let products = MF::batch_mul_vec(swap_cond_arith, deltas);

            // Step 7: Apply arithmetic updates (local).
            for (p, &(i, ij, _)) in pair_indices.iter().enumerate() {
                vec[i] = vec[i] + products[p];
                vec[ij] = vec[ij] - products[p];
            }

            // Step 8: Update precomputed bits with boolean MUX (2 broadcasts).
            // For each pair p and bit b:
            //   correction[p][b] = swap_cond_bool[p] AND (bits[i][b] XOR bits[ij][b])
            //   bits[i][b]  ^= correction[p][b]
            //   bits[ij][b] ^= correction[p][b]
            let mut and_lefts = Vec::with_capacity(num_pairs * bit_length);
            let mut and_rights = Vec::with_capacity(num_pairs * bit_length);
            for (p, &(i, ij, _)) in pair_indices.iter().enumerate() {
                for b in 0..bit_length {
                    and_lefts.push(swap_cond_bools[p]);
                    and_rights.push(bits[i][b] ^ bits[ij][b]);
                }
            }
            let mut beaver_source = DummyBooleanBeaverSource;
            let corrections =
                BooleanOps::batch_beaver_bitwise_and(&and_lefts, &and_rights, &mut beaver_source);
            for (p, &(i, ij, _)) in pair_indices.iter().enumerate() {
                for b in 0..bit_length {
                    let correction = corrections[p * bit_length + b];
                    bits[i][b] ^= correction;
                    bits[ij][b] ^= correction;
                }
            }

            j /= 2;
        }
        k *= 2;
    }

    vec
}

/// Radix sort for (key, value1, value2) tuples using precomputed bits (bitonic network).
///
/// Sorts by the key field in ascending order using the same batched-bits approach
/// as `radix_sort`. Only the key bits are precomputed; val1 and val2 are carried
/// alongside using arithmetic MUX.
///
/// The input vector length must be a power of 2.
pub fn radix_sort_by_key(
    mut vec: Vec<(MF, MF, MF)>,
    bit_length: usize,
) -> Vec<(MF, MF, MF)> {
    let n = vec.len();

    if n < 2 || (n & (n - 1)) != 0 {
        panic!(
            "radix_sort_by_key: vector length must be a power of 2 and at least 2, got {}",
            n
        );
    }

    // Fast path: when all values are public, sort locally without MPC operations.
    if vec.iter().all(|(k, v1, v2)| {
        matches!(k, MF::Public(_)) && matches!(v1, MF::Public(_)) && matches!(v2, MF::Public(_))
    }) {
        vec.sort_by_key(|(k, _, _)| match k {
            MF::Public(x) => *x,
            _ => unreachable!(),
        });
        return vec;
    }

    // Step 1: Batch precompute key bits in one broadcast round, then truncate.
    let keys: Vec<MF> = vec.iter().map(|(k, _, _)| *k).collect();
    let mut bits: Vec<Vec<bool>> = MF::batch_bit_decompose_bool(&keys)
        .into_iter()
        .map(|mut b| {
            b.truncate(bit_length);
            b
        })
        .collect();

    let mut k = 2;
    while k <= n {
        let mut j = k / 2;
        while j > 0 {
            // Collect all pairs for this comparison level.
            let mut pair_indices: Vec<(usize, usize, bool)> = Vec::new();
            for i in 0..n {
                let ij = i ^ j;
                if ij > i {
                    let swap_if_greater = (i & k) == 0;
                    pair_indices.push((i, ij, swap_if_greater));
                }
            }
            let num_pairs = pair_indices.len();

            // Step 2: Batch compare all pairs using precomputed key bits.
            let left_bits: Vec<Vec<bool>> =
                pair_indices.iter().map(|&(i, _, _)| bits[i].clone()).collect();
            let right_bits: Vec<Vec<bool>> =
                pair_indices.iter().map(|&(_, ij, _)| bits[ij].clone()).collect();
            let is_less_results = MF::batch_compare_bits(&left_bits, &right_bits);

            // Step 3: Boolean swap conditions (local).
            // IMPORTANT: XOR with public constant must be applied to king's share only.
            let swap_cond_bools: Vec<bool> = pair_indices
                .iter()
                .zip(is_less_results.iter())
                .map(|(&(_, _, swig), &is_less)| {
                    if swig { is_less ^ Net::am_king() } else { is_less }
                })
                .collect();

            // Step 4: Batch B2A conversion.
            let swap_cond_arith: Vec<MF> = MF::batch_b2a_bool(&swap_cond_bools);

            // Step 5: Compute arithmetic deltas for key, val1, val2 (local).
            // Layout: [delta_key_0, delta_val1_0, delta_val2_0, delta_key_1, ...]
            let mut all_conds = Vec::with_capacity(num_pairs * 3);
            let mut all_deltas = Vec::with_capacity(num_pairs * 3);
            for (p, &(i, ij, _)) in pair_indices.iter().enumerate() {
                let (ki, v1i, v2i) = vec[i];
                let (kj, v1j, v2j) = vec[ij];
                all_conds.push(swap_cond_arith[p]);
                all_conds.push(swap_cond_arith[p]);
                all_conds.push(swap_cond_arith[p]);
                all_deltas.push(kj - ki);
                all_deltas.push(v1j - v1i);
                all_deltas.push(v2j - v2i);
            }

            // Step 6: Batch multiply (2 broadcasts for all pairs × 3 fields).
            let products = MF::batch_mul_vec(all_conds, all_deltas);

            // Step 7: Apply updates (local).
            for (p, &(i, ij, _)) in pair_indices.iter().enumerate() {
                let dk = products[p * 3];
                let dv1 = products[p * 3 + 1];
                let dv2 = products[p * 3 + 2];
                vec[i].0 = vec[i].0 + dk;
                vec[i].1 = vec[i].1 + dv1;
                vec[i].2 = vec[i].2 + dv2;
                vec[ij].0 = vec[ij].0 - dk;
                vec[ij].1 = vec[ij].1 - dv1;
                vec[ij].2 = vec[ij].2 - dv2;
            }

            // Step 8: Update precomputed key bits with boolean MUX (2 broadcasts).
            let mut and_lefts = Vec::with_capacity(num_pairs * bit_length);
            let mut and_rights = Vec::with_capacity(num_pairs * bit_length);
            for (p, &(i, ij, _)) in pair_indices.iter().enumerate() {
                for b in 0..bit_length {
                    and_lefts.push(swap_cond_bools[p]);
                    and_rights.push(bits[i][b] ^ bits[ij][b]);
                }
            }
            let mut beaver_source = DummyBooleanBeaverSource;
            let corrections =
                BooleanOps::batch_beaver_bitwise_and(&and_lefts, &and_rights, &mut beaver_source);
            for (p, &(i, ij, _)) in pair_indices.iter().enumerate() {
                for b in 0..bit_length {
                    let correction = corrections[p * bit_length + b];
                    bits[i][b] ^= correction;
                    bits[ij][b] ^= correction;
                }
            }

            j /= 2;
        }
        k *= 2;
    }

    vec
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::UniformRand;
    use ark_std::test_rng;
    use mpc_algebra::Reveal;

    #[test]
    fn test_bitonic_sort_public_values() {
        let test_cases = vec![
            vec![F::from(3u64), F::from(1u64)],
            vec![F::from(4u64), F::from(2u64), F::from(3u64), F::from(1u64)],
            vec![
                F::from(8u64),
                F::from(3u64),
                F::from(6u64),
                F::from(1u64),
                F::from(7u64),
                F::from(5u64),
                F::from(2u64),
                F::from(4u64),
            ],
        ];

        for input in test_cases {
            let n = input.len();
            println!("Testing bitonic_sort with n={}", n);

            let mpc_input: Vec<MF> = input.iter().map(|&x| MF::Public(x)).collect();
            let sorted = bitonic_sort(mpc_input);
            let sorted_values: Vec<F> = sorted.iter().map(|x| x.reveal()).collect();

            let mut expected = input.clone();
            expected.sort();

            assert_eq!(
                sorted_values, expected,
                "Bitonic sort failed for input size {}",
                n
            );
            println!("  Passed for n={}", n);
        }
    }

    #[test]
    fn test_radix_sort_public_values() {
        let test_cases = vec![
            vec![F::from(3u64), F::from(1u64)],
            vec![F::from(4u64), F::from(2u64), F::from(3u64), F::from(1u64)],
            vec![
                F::from(8u64),
                F::from(3u64),
                F::from(6u64),
                F::from(1u64),
                F::from(7u64),
                F::from(5u64),
                F::from(2u64),
                F::from(4u64),
            ],
        ];

        for input in test_cases {
            let n = input.len();
            println!("Testing radix_sort with n={}", n);

            let mpc_input: Vec<MF> = input.iter().map(|&x| MF::Public(x)).collect();
            let sorted = radix_sort(mpc_input, 64);
            let sorted_values: Vec<F> = sorted.iter().map(|x| x.reveal()).collect();

            let mut expected = input.clone();
            expected.sort();

            assert_eq!(
                sorted_values, expected,
                "Radix sort failed for input size {}",
                n
            );
            println!("  Passed for n={}", n);
        }
    }

    #[test]
    fn test_radix_sort_by_key_public_values() {
        // Create (key, val1, val2) tuples with distinct keys
        let keys = vec![
            F::from(4u64),
            F::from(2u64),
            F::from(3u64),
            F::from(1u64),
        ];
        let vals: Vec<(F, F)> = keys
            .iter()
            .enumerate()
            .map(|(i, _)| (F::from(i as u64 * 10), F::from(i as u64 * 100)))
            .collect();

        let mpc_input: Vec<(MF, MF, MF)> = keys
            .iter()
            .zip(vals.iter())
            .map(|(&k, &(v1, v2))| (MF::Public(k), MF::Public(v1), MF::Public(v2)))
            .collect();

        let sorted = radix_sort_by_key(mpc_input, 64);

        // Keys should be sorted ascending; values should follow the key
        let sorted_keys: Vec<F> = sorted.iter().map(|(k, _, _)| k.reveal()).collect();
        let mut expected_keys = keys.clone();
        expected_keys.sort();

        assert_eq!(sorted_keys, expected_keys, "Keys not sorted correctly");

        // Verify that values followed their corresponding key
        for (i, &(_, v1, v2)) in sorted.iter().enumerate() {
            let key = sorted_keys[i];
            let orig_idx = keys.iter().position(|&k| k == key).unwrap();
            assert_eq!(
                v1.reveal(),
                vals[orig_idx].0,
                "val1 not correctly moved with key"
            );
            assert_eq!(
                v2.reveal(),
                vals[orig_idx].1,
                "val2 not correctly moved with key"
            );
        }

        println!("radix_sort_by_key test passed");
    }
}
