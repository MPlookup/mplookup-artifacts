//! Oblivious sorting operations for MPC fields.
//!
//! This module provides secure bitonic sorting algorithms that work
//! with MpcField values, maintaining obliviousness in MPC contexts.

use ark_ff::{One, Zero};
use mpc_algebra::honest_but_curious::MpcField;

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
}
