//! Secure oblivious lookup permutation algorithm.
//!
//! This module implements the main secure oblivious lookup permutation
//! algorithm that uses polynomial operations and sorting to create
//! a permutation satisfying the lookup relation.

use ark_ff::{BigInteger, FpParameters, One, PrimeField, Zero};
use ark_std::{end_timer, start_timer};
use mpc_algebra::honest_but_curious::MpcField;
use mpc_algebra::Reveal;

use super::polynomial::{build_subproduct_tree, multipoint_eval};
use super::sorting::{radix_sort, radix_sort_by_key};

type F = ark_bls12_377::Fr;
type MF = MpcField<F>;

// Sort key multiplier for oblivious sorting
const SORT_KEY_MULTIPLIER: u64 = 1u64 << 32;

/// Secure oblivious lookup permutation using MPC-friendly operations.
///
/// This function implements an oblivious lookup permutation that uses secure comparison
/// and sorting operations suitable for multi-party computation (MPC).
///
/// # Arguments
///
/// * `n` - Size of the table `t`
/// * `m` - Number of queries in `f`
/// * `t` - The table vector (MpcField values)
/// * `f` - The query vector (MpcField values)
/// * `debug` - Enable debug output
///
/// # Returns
///
/// A tuple `(t_prime, f_prime)` where:
/// * `t_prime` is a permuted version of `t`
/// * `f_prime` is a sorted version of `f`
pub fn secure_oblivious_lookup_permutation(
    n: usize,
    m: usize,
    t: &[MF],
    f: &[MF],
    debug: bool,
) -> (Vec<MF>, Vec<MF>) {
    // M and N must be power of two. M <= N. Otherwise, throw error.
    assert!(m <= n, "Number of queries m must be less than or equal to table size n");
    assert!(m.is_power_of_two(), "Number of queries m must be a power of two");
    assert!(n.is_power_of_two(), "Table size n must be a power of two");

    let timer_total =
        start_timer!(|| format!("secure_oblivious_lookup_permutation total n={} m={}", n, m));

    let unwrap_vec = |v: &[MF]| -> Vec<F> {
        let mut out = Vec::with_capacity(v.len());
        for x in v.iter() {
            let val = x.reveal();
            out.push(val);
        }
        out
    };

    let debug_print_vec = |label: &str, v: &[MF]| {
        if !debug {
            return;
        }
        let plain = unwrap_vec(v);
        println!("DEBUG {} (len={}):", label, plain.len());
        for i in 0..plain.len() {
            println!("{}[{}] = {}", label, i, plain[i]);
        }
    };

    let m_padded = m.next_power_of_two();
    let n_padded = n.next_power_of_two();

    let timer_step1 = start_timer!(|| format!(
        "secure_oblivious_lookup_permutation step 1: Sort the query vector n={} m={}",
        n, m
    ));
    let mut f_prime = f.to_vec();
    f_prime.resize(m_padded, MF::Public(F::zero()));
    f_prime = radix_sort(f_prime, <F as PrimeField>::Params::MODULUS.num_bits() as usize);
    f_prime.truncate(m);
    if debug {
        debug_print_vec("f_prime after sort", &f_prime);
        let plain = unwrap_vec(&f_prime);
        {
            let tmp = plain.clone();
            let mut tmp_sorted = tmp.clone();
            tmp_sorted.sort();
            if tmp != tmp_sorted {
                panic!("DEBUG: f_prime is not sorted (expected ascending)");
            } else {
                println!("DEBUG: f_prime is sorted (OK)");
            }
        }
    }

    let mut f_prime_full = f_prime.clone();
    f_prime_full.resize(n, MF::Public(F::zero()));
    end_timer!(timer_step1);

    let timer_step2 = start_timer!(|| format!(
        "secure_oblivious_lookup_permutation step 2: Compute first-occurrence indicators n={} m={}",
        n, m
    ));
    let mut f_double = vec![MF::Public(F::zero()); n];
    f_double[0] = MF::Public(F::one());
    for i in 1..n {
        let is_before_m = if i < m {
            MF::Public(F::one())
        } else {
            MF::Public(F::zero())
        };
        let eq = f_prime_full[i].secure_eq(&f_prime_full[i - 1]);
        f_double[i] = (MF::Public(F::one()) - eq) * is_before_m;
    }
    if debug {
        debug_print_vec("f_prime_full (after pad)", &f_prime_full);
        debug_print_vec("f_double after first-occurrence calc", &f_double);
        let (f_plain, f_double_plain) = (unwrap_vec(&f_prime_full), unwrap_vec(&f_double));
        {
            let mut expected = vec![F::zero(); n];
            expected[0] = F::one();
            for i in 1..n {
                let is_before_m_plain = if i < m { F::one() } else { F::zero() };
                let eq_plain = if f_plain[i] == f_plain[i - 1] {
                    F::one()
                } else {
                    F::zero()
                };
                expected[i] = (F::one() - eq_plain) * is_before_m_plain;
            }
            if expected != f_double_plain {
                panic!("DEBUG: f_double does not match expected first-occurrence indicators");
            } else {
                println!("DEBUG: f_double matches expected first-occurrence indicators (OK)");
            }
        }
    }
    end_timer!(timer_step2);

    let timer_step3 = start_timer!(|| format!(
        "secure_oblivious_lookup_permutation step 3: Construct membership polynomial n={} m={}",
        n, m
    ));
    let f_tree = build_subproduct_tree(&f_prime_full[0..m]);
    let p_coeffs = f_tree.last().unwrap()[0].clone();
    if debug {
        let coeffs_plain = unwrap_vec(&p_coeffs);
        println!("DEBUG p_coeffs (len={}):", coeffs_plain.len());
        for i in 0..coeffs_plain.len() {
            println!("p_coeffs[{}] = {}", i, coeffs_plain[i]);
        }
    }
    end_timer!(timer_step3);

    let timer_step4 = start_timer!(|| {
        format!("secure_oblivious_lookup_permutation step 4: Evaluate polynomial at table elements n={} m={}", n, m)
    });
    let t_tree = build_subproduct_tree(t);
    let p_values = multipoint_eval(&p_coeffs, &t_tree);
    if debug {
        debug_print_vec("p_values (eval at t)", &p_values);
    }
    end_timer!(timer_step4);

    let timer_step5 = start_timer!(|| format!(
        "secure_oblivious_lookup_permutation step 5: Compute unused indicators n={} m={}",
        n, m
    ));
    let unused: Vec<MF> = p_values
        .iter()
        .map(|&pv| {
            pv.secure_neq(&MF::Public(F::zero()))
        })
        .collect();
    if debug {
        debug_print_vec("unused flags", &unused);
        let (pv_plain, unused_plain) = (unwrap_vec(&p_values), unwrap_vec(&unused));
        {
            let mut expected_unused = Vec::with_capacity(pv_plain.len());
            for &v in pv_plain.iter() {
                expected_unused.push(if v != F::zero() { F::one() } else { F::zero() });
            }
            if expected_unused != unused_plain {
                panic!("DEBUG: unused flags do not match p_values != 0");
            } else {
                println!("DEBUG: unused flags validated against p_values (OK)");
            }
        }
    }
    end_timer!(timer_step5);

    let timer_step6 = start_timer!(|| format!(
        "secure_oblivious_lookup_permutation step 6: Compact unused elements n={} m={}",
        n, m
    ));
    let mut unused_value_rec: Vec<(MF, MF, MF)> = (0..n)
        .map(|j| {
            let large_val = MF::Public(F::from(SORT_KEY_MULTIPLIER));
            let index_val = MF::Public(F::from(j as u64));
            let sort_key = (MF::Public(F::one()) - unused[j]) * large_val + index_val;
            (sort_key, t[j], unused[j])
        })
        .collect();

    unused_value_rec.resize(
        n_padded,
        (
            MF::Public(F::zero()),
            MF::Public(F::zero()),
            MF::Public(F::zero()),
        ),
    );
    unused_value_rec = radix_sort_by_key(unused_value_rec, 34);
    unused_value_rec.truncate(n);
    if debug {
        println!("DEBUG: unused_value_rec after compaction (first 16 shown)");
        for i in 0..std::cmp::min(16, unused_value_rec.len()) {
            let (k, v, u) = unused_value_rec[i];
            let (kp, vp, up) = (unwrap_vec(&[k]), unwrap_vec(&[v]), unwrap_vec(&[u]));
            println!("  rec[{}] key={} val={} unused={}", i, kp[0], vp[0], up[0]);
        }
    }
    end_timer!(timer_step6);

    let timer_step7 = start_timer!(|| format!(
        "secure_oblivious_lookup_permutation step 7: Compact unfilled positions n={} m={}",
        n, m
    ));
    let mut unfilled_position_rec: Vec<(MF, MF, MF)> = (0..n)
        .map(|i| {
            let unfilled = MF::Public(F::one()) - f_double[i];
            let large_val = MF::Public(F::from(SORT_KEY_MULTIPLIER));
            let index_val = MF::Public(F::from(i as u64));
            let sort_key = (MF::Public(F::one()) - unfilled) * large_val + index_val;
            (sort_key, MF::Public(F::from(i as u64)), unfilled)
        })
        .collect();

    unfilled_position_rec.resize(
        n_padded,
        (
            MF::Public(F::zero()),
            MF::Public(F::zero()),
            MF::Public(F::zero()),
        ),
    );
    unfilled_position_rec = radix_sort_by_key(unfilled_position_rec, 34);
    unfilled_position_rec.truncate(n);
    if debug {
        println!("DEBUG: unfilled_position_rec after compaction (first 16 shown)");
        for i in 0..std::cmp::min(16, unfilled_position_rec.len()) {
            let (k, v, u) = unfilled_position_rec[i];
            let (kp, vp, up) = (unwrap_vec(&[k]), unwrap_vec(&[v]), unwrap_vec(&[u]));
            println!(
                "  rec[{}] key={} pos={} unfilled={}",
                i, kp[0], vp[0], up[0]
            );
        }
    }
    end_timer!(timer_step7);

    let timer_step8 = start_timer!(|| format!(
        "secure_oblivious_lookup_permutation step 8: Create write records 1 n={} m={}",
        n, m
    ));
    let write_rec1: Vec<(MF, MF, MF)> = (0..n)
        .map(|mm| {
            let pos = unfilled_position_rec[mm].1;
            let val = unused_value_rec[mm].1;
            let valid = unfilled_position_rec[mm].2;
            (pos, val, valid)
        })
        .collect();
    end_timer!(timer_step8);

    let timer_step9 = start_timer!(|| format!(
        "secure_oblivious_lookup_permutation step 9: Create write records 2 n={} m={}",
        n, m
    ));
    let write_rec2: Vec<(MF, MF, MF)> = (0..n)
        .map(|k| {
            let pos = MF::Public(F::from(k as u64));
            let val = f_prime_full[k];
            let valid = f_double[k];
            (pos, val, valid)
        })
        .collect();
    end_timer!(timer_step9);

    let timer_step10 = start_timer!(|| format!(
        "secure_oblivious_lookup_permutation step 10: Combine and sort write records n={} m={}",
        n, m
    ));
    let mut all_writes = Vec::with_capacity(2 * n);
    all_writes.extend_from_slice(&write_rec1);
    all_writes.extend_from_slice(&write_rec2);

    let mut all_writes_keyed: Vec<(MF, MF, MF)> = all_writes
        .iter()
        .map(|(pos, val, valid)| {
            let large_val = MF::Public(F::from(SORT_KEY_MULTIPLIER));
            let sort_key = (MF::Public(F::one()) - *valid) * large_val + *pos;
            (sort_key, *val, *valid)
        })
        .collect();

    let all_writes_len = all_writes_keyed.len();
    let all_writes_padded = all_writes_len.next_power_of_two();
    all_writes_keyed.resize(
        all_writes_padded,
        (
            MF::Public(F::zero()),
            MF::Public(F::zero()),
            MF::Public(F::zero()),
        ),
    );
    all_writes_keyed = radix_sort_by_key(all_writes_keyed, 34);
    if debug {
        println!("DEBUG: all_writes_keyed after sort (first 32 shown)");
        for i in 0..std::cmp::min(32, all_writes_keyed.len()) {
            let (k, v, valid) = all_writes_keyed[i];
            let (kp, vp, vld) = (unwrap_vec(&[k]), unwrap_vec(&[v]), unwrap_vec(&[valid]));
            println!(
                "  write[{}] key={} val={} valid={}",
                i, kp[0], vp[0], vld[0]
            );
        }
    }
    end_timer!(timer_step10);

    let timer_step11 = start_timer!(|| format!(
        "secure_oblivious_lookup_permutation step 11: Construct output permutation n={} m={}",
        n, m
    ));
    let t_prime: Vec<MF> = all_writes_keyed[0..n].iter().map(|r| r.1).collect();
    if debug {
        debug_print_vec("t_prime (final)", &t_prime);
        let (t_plain, f_plain, tprime_plain, fprime_plain) = (
            unwrap_vec(t),
            unwrap_vec(f),
            unwrap_vec(&t_prime),
            unwrap_vec(&f_prime),
        );
        println!("DEBUG: Running final verify_lookup_permutation_outputs check");
        verify_lookup_permutation_outputs_in_plaintext(n, m, &t_plain, &f_plain, &tprime_plain, &fprime_plain, debug);
        println!("DEBUG: final verify_lookup_permutation_outputs passed (OK)");
    }
    end_timer!(timer_step11);

    end_timer!(timer_total);

    (t_prime, f_prime)
}

/// Verify outputs of secure_oblivious_lookup_permutation.
///
/// Asserts the following properties:
/// 1) t' has length n and f' has length m
/// 2) t' is a permutation of the original table t
/// 3) f' is a permutation of the original queries f
/// 4) For each first occurrence in f', require f'[i] == t'[i]
pub fn verify_lookup_permutation_outputs_in_plaintext(
    n: usize,
    m: usize,
    t: &Vec<F>,
    f: &Vec<F>,
    t_prime: &Vec<F>,
    f_prime: &Vec<F>,
    debug: bool,
) {
    assert_eq!(t_prime.len(), n, "t' length mismatch");
    assert_eq!(f_prime.len(), m, "f' length mismatch");

    let mut t_sorted = t.clone();
    t_sorted.sort();
    let mut t_prime_sorted = t_prime.clone();
    t_prime_sorted.sort();
    assert_eq!(t_sorted, t_prime_sorted, "t' is not a permutation of t");

    let mut f_sorted = f.clone();
    f_sorted.sort();
    let mut f_prime_sorted = f_prime.clone();
    f_prime_sorted.sort();
    assert_eq!(f_sorted, f_prime_sorted, "f' is not a permutation of f");

    for i in 0..m {
        if i == 0 || f_prime[i] != f_prime[i - 1] {
            assert_eq!(
                f_prime[i], t_prime[i],
                "first-occurrence value mismatch at index {}",
                i
            );
        }
    }

    if debug {
        println!("All verification checks passed. This is a debug test. In real scenarios, values are secret-shared so a zero-knowledge proof (implemented in client.rs) is used to verify the result.")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::UniformRand;
    use ark_std::test_rng;
    use ark_std::rand::RngCore;

    #[test]
    fn test_secure_oblivious_lookup_permutation_local() {
        let n = 8;
        let m = 4;
        let t_plain = vec![
            F::from(10u64),
            F::from(20u64),
            F::from(30u64),
            F::from(0u64),
            F::from(40u64),
            F::from(50u64),
            F::from(60u64),
            F::from(70u64),
        ];
        let f_plain = vec![F::from(30u64), F::from(0u64), F::from(30u64), F::from(0u64)];

        let t: Vec<MF> = t_plain.iter().map(|&x| MF::Public(x)).collect();
        let f: Vec<MF> = f_plain.iter().map(|&x| MF::Public(x)).collect();

        let (t_prime, f_prime) = secure_oblivious_lookup_permutation(n, m, &t, &f, false);

        for (i, val) in t_prime.iter().enumerate() {
            println!("t'[{}] = {}", i, val);
        }
        for (i, val) in f_prime.iter().enumerate() {
            println!("f'[{}] = {}", i, val);
        }

        let t_prime_plain: Vec<F> = t_prime.iter().map(|x| x.reveal()).collect();
        let f_prime_plain: Vec<F> = f_prime.iter().map(|x| x.reveal()).collect();

        verify_lookup_permutation_outputs_in_plaintext(n, m, &t_plain, &f_plain, &t_prime_plain, &f_prime_plain, true);
    }

    #[test]
    fn test_secure_oblivious_lookup_permutation_large_local() {
        let mut rng = test_rng();
        let sizes = [8usize, 16, 32];

        for &n in sizes.iter() {
            let m = n / 2;

            println!("Testing n = {}, m = {}", n, m);

            let mut t_plain: Vec<F> = Vec::with_capacity(n);
            while t_plain.len() < n {
                let x = F::rand(&mut rng);
                if !t_plain.iter().any(|y| *y == x) {
                    t_plain.push(x);
                }
            }

            let mut f_plain: Vec<F> = (0..m)
                .map(|_| {
                    let idx = (rng.next_u64() as usize) % n;
                    t_plain[idx]
                })
                .collect();

            {
                let mut tmp = f_plain.clone();
                tmp.sort();
                let has_dup = tmp.windows(2).any(|w| w[0] == w[1]);
                if !has_dup && m >= 2 {
                    f_plain[m - 1] = f_plain[0];
                }
            }

            let t: Vec<MF> = t_plain.iter().map(|&x| MF::Public(x)).collect();
            let f: Vec<MF> = f_plain.iter().map(|&x| MF::Public(x)).collect();

            let (t_prime, f_prime) = secure_oblivious_lookup_permutation(n, m, &t, &f, false);

            let t_prime_plain: Vec<F> = t_prime.iter().map(|x| x.reveal()).collect();
            let f_prime_plain: Vec<F> = f_prime.iter().map(|x| x.reveal()).collect();

            verify_lookup_permutation_outputs_in_plaintext(
                n,
                m,
                &t_plain,
                &f_plain,
                &t_prime_plain,
                &f_prime_plain,
                true,
            );
        }
    }
}

/// Naive permutation method.
///
/// This function implements a straightforward oblivious permutation algorithm that
/// sorts the query vector and then constructs a permutation of the table vector.
/// It uses nested loops and is less efficient than the polynomial-based approach,
/// but serves as a baseline for comparison.
///
/// # Arguments
///
/// * `n` - Size of the table `t`
/// * `m` - Number of queries in `f`
/// * `t` - The table vector (MpcField values)
/// * `f` - The query vector (MpcField values)
/// * `debug` - Enable debug output
///
/// # Returns
///
/// A tuple `(t_prime, f_prime)` where:
/// * `t_prime` is a permuted version of `t`
/// * `f_prime` is a sorted version of `f`
pub fn naive_secure_oblivious_lookup_permutation(
    n: usize,
    m: usize,
    t: &[MF],
    f: &[MF],
    debug: bool,
) -> (Vec<MF>, Vec<MF>) {
    // M and N must be power of two. M <= N. Otherwise, throw error.
    assert!(m <= n, "Number of queries m must be less than or equal to table size n");
    assert!(m.is_power_of_two(), "Number of queries m must be a power of two");
    assert!(n.is_power_of_two(), "Table size n must be a power of two");

    let timer_total = start_timer!(|| format!("naive_secure_oblivious_lookup_permutation total n={} m={}", n, m));

    let unwrap_vec = |v: &[MF]| -> Vec<F> {
        let mut out = Vec::with_capacity(v.len());
        for x in v.iter() {
            let val = x.reveal();
            out.push(val);
        }
        out
    };

    let debug_print_vec = |label: &str, v: &[MF]| {
        if !debug {
            return;
        }
        let plain = unwrap_vec(v);
        println!("DEBUG {} (len={}):", label, plain.len());
        for i in 0..plain.len() {
            println!("{}[{}] = {}", label, i, plain[i]);
        }
    };

    // Step 1: Sort the query vector f to get f'
    let timer_step1 = start_timer!(|| "naive_secure_oblivious_lookup_permutation step 1: Sort query vector");
    let m_padded = m.next_power_of_two();
    let mut f_prime = f.to_vec();
    f_prime.resize(m_padded, MF::Public(F::zero()));
    f_prime = radix_sort(f_prime, <F as PrimeField>::Params::MODULUS.num_bits() as usize);
    f_prime.truncate(m);
    if debug {
        debug_print_vec("f_prime after sort", &f_prime);
    }
    end_timer!(timer_step1);

    // Step 2: Compute f'' (first-occurrence indicators)
    // f''[i] = 1 if f'[i] is not the same as f'[i-1], else 0
    let timer_step2 = start_timer!(|| "naive_secure_oblivious_lookup_permutation step 2: Compute first-occurrence indicators");
    let mut f_double = vec![MF::Public(F::zero()); n];
    f_double[0] = MF::Public(F::one());
    
    for i in 1..m {
        // Check if f'[i] == f'[i-1]
        let is_same = f_prime[i].secure_eq(&f_prime[i - 1]);
        // f''[i] = 1 - is_same
        f_double[i] = MF::Public(F::one()) - is_same;
    }
    
    if debug {
        debug_print_vec("f_double (first-occurrence indicators)", &f_double);
    }
    end_timer!(timer_step2);

    // Step 3: Initialize t' by filling with f' where f'' = 1
    // t'[i] = f''[i] * f'[i] for i < m
    // This fills positions where first occurrence happens
    let timer_step3 = start_timer!(|| "naive_secure_oblivious_lookup_permutation step 3: Initial filling of t'");
    let mut t_prime = vec![MF::Public(F::zero()); n];
    let mut fill_mask = vec![MF::Public(F::zero()); n];
    
    for i in 0..m {
        // t'[i] = f''[i] * f'[i]
        t_prime[i] = f_double[i] * f_prime[i];
        fill_mask[i] = f_double[i];
    }
    
    if debug {
        debug_print_vec("t_prime after initial fill", &t_prime);
        debug_print_vec("fill_mask after initial fill", &fill_mask);
    }
    end_timer!(timer_step3);

    // Step 4: Mark which elements in t have been used
    // We iterate through each query and mark corresponding table elements
    let timer_step4 = start_timer!(|| "naive_secure_oblivious_lookup_permutation step 4: Mark used table elements");
    let mut is_table_element_unused = vec![MF::Public(F::one()); n];
    
    for i in 0..m {
        let mut found = MF::Public(F::zero());
        
        for j in 0..n {
            // Condition 1: t[j] == f'[i]
            let condition1 = t[j].secure_eq(&f_prime[i]);
            
            // Condition 2: is_table_element_unused[j]
            let condition2 = is_table_element_unused[j];
            
            // Condition 3: f_double[i] (first occurrence)
            let condition3 = f_double[i];
            
            // All conditions must be true
            let condition = condition1 * condition2 * condition3;
            
            // notFound = 1 - found
            let not_found = MF::Public(F::one()) - found;
            
            // continue_wire = condition * not_found
            let continue_wire = condition * not_found;
            let not_continue = MF::Public(F::one()) - continue_wire;
            
            // Update is_table_element_unused[j] = not_continue * is_table_element_unused[j]
            is_table_element_unused[j] = not_continue * is_table_element_unused[j];
            
            // Update found = found + continue_wire
            found = found + continue_wire;
        }
    }
    
    if debug {
        debug_print_vec("is_table_element_unused after marking", &is_table_element_unused);
    }
    end_timer!(timer_step4);

    // Step 5: Fill remaining positions in t' with unused elements from t
    let timer_step5 = start_timer!(|| "naive_secure_oblivious_lookup_permutation step 5: Fill remaining positions");
    for i in 0..n {
        let mut found = MF::Public(F::zero());
        
        // Check if position i needs filling (fill_mask[i] == 0)
        let outer_continue = MF::Public(F::one()) - fill_mask[i];
        
        for j in 0..n {
            // Condition 1: outer_continue (position needs filling)
            let condition1 = outer_continue;
            
            // Condition 2: is_table_element_unused[j]
            let condition2 = is_table_element_unused[j];
            
            // Condition 3: not found yet
            let not_found = MF::Public(F::one()) - found;
            
            // inner_continue = all conditions
            let inner_continue = condition1 * condition2 * not_found;
            let not_inner_continue = MF::Public(F::one()) - inner_continue;
            
            // Update is_table_element_unused[j]
            is_table_element_unused[j] = is_table_element_unused[j] * not_inner_continue;
            
            // Update t'[i] = inner_continue * t[j] + not_inner_continue * t'[i]
            let left_item = inner_continue * t[j];
            let right_item = not_inner_continue * t_prime[i];
            t_prime[i] = left_item + right_item;
            
            // Update found
            found = found + inner_continue;
        }
        
        // Update fill_mask[i]
        fill_mask[i] = fill_mask[i] + found;
    }
    
    // Step 6: Return f_prime as the sorted queries
    // f_prime should just be the sorted queries (with duplicates preserved)
    // The client will pad it to length n using t_prime values for the proof system
    let timer_step6 = start_timer!(|| "naive_secure_oblivious_lookup_permutation step 6: Finalize f_prime");
    let f_prime_return = f_prime[0..m].to_vec();
    
    if debug {
        debug_print_vec("t_prime (final)", &t_prime);
        debug_print_vec("f_prime (sorted queries for return)", &f_prime_return);
        
        let (t_plain, f_plain, tprime_plain, fprime_return_plain) = (
            unwrap_vec(t),
            unwrap_vec(f),
            unwrap_vec(&t_prime),
            unwrap_vec(&f_prime_return),
        );
        println!("DEBUG: Running final verify_lookup_permutation_outputs check");
        verify_lookup_permutation_outputs_in_plaintext(n, m, &t_plain, &f_plain, &tprime_plain, &fprime_return_plain, debug);
        println!("DEBUG: final verify_lookup_permutation_outputs passed (OK)");
    }
    end_timer!(timer_step6);
    end_timer!(timer_step5);

    end_timer!(timer_total);

    (t_prime, f_prime_return)
}

#[cfg(test)]
mod naive_tests {
    use super::*;
    use ark_ff::UniformRand;
    use ark_std::test_rng;
    use ark_std::rand::RngCore;

    #[test]
    fn test_naive_secure_oblivious_lookup_permutation_local() {
        let n = 8;
        let m = 4;
        let t_plain = vec![
            F::from(10u64),
            F::from(20u64),
            F::from(30u64),
            F::from(0u64),
            F::from(40u64),
            F::from(50u64),
            F::from(60u64),
            F::from(70u64),
        ];
        let f_plain = vec![F::from(30u64), F::from(0u64), F::from(30u64), F::from(0u64)];

        let t: Vec<MF> = t_plain.iter().map(|&x| MF::Public(x)).collect();
        let f: Vec<MF> = f_plain.iter().map(|&x| MF::Public(x)).collect();

        let (t_prime, f_prime) = naive_secure_oblivious_lookup_permutation(n, m, &t, &f, false);

        for (i, val) in t_prime.iter().enumerate() {
            println!("t'[{}] = {}", i, val);
        }
        for (i, val) in f_prime.iter().enumerate() {
            println!("f'[{}] = {}", i, val);
        }

        let t_prime_plain: Vec<F> = t_prime.iter().map(|x| x.reveal()).collect();
        let f_prime_plain: Vec<F> = f_prime.iter().map(|x| x.reveal()).collect();

        // Verify with only the first m elements of f_prime
        let f_prime_plain_m = f_prime_plain[0..m].to_vec();
        verify_lookup_permutation_outputs_in_plaintext(n, m, &t_plain, &f_plain, &t_prime_plain, &f_prime_plain_m, true);
    }

    #[test]
    fn test_naive_secure_oblivious_lookup_permutation_large_local() {
        let mut rng = test_rng();
        let sizes = [8usize, 16, 32];

        for &n in sizes.iter() {
            let m = n / 2;

            println!("Testing n = {}, m = {}", n, m);

            let mut t_plain: Vec<F> = Vec::with_capacity(n);
            while t_plain.len() < n {
                let x = F::rand(&mut rng);
                if !t_plain.iter().any(|y| *y == x) {
                    t_plain.push(x);
                }
            }

            let mut f_plain: Vec<F> = (0..m)
                .map(|_| {
                    let idx = (rng.next_u64() as usize) % n;
                    t_plain[idx]
                })
                .collect();

            {
                let mut tmp = f_plain.clone();
                tmp.sort();
                let has_dup = tmp.windows(2).any(|w| w[0] == w[1]);
                if !has_dup && m >= 2 {
                    f_plain[m - 1] = f_plain[0];
                }
            }

            let t: Vec<MF> = t_plain.iter().map(|&x| MF::Public(x)).collect();
            let f: Vec<MF> = f_plain.iter().map(|&x| MF::Public(x)).collect();

            let (t_prime, f_prime) = naive_secure_oblivious_lookup_permutation(n, m, &t, &f, false);

            let t_prime_plain: Vec<F> = t_prime.iter().map(|x| x.reveal()).collect();
            let f_prime_plain: Vec<F> = f_prime.iter().map(|x| x.reveal()).collect();

            // Verify with only the first m elements of f_prime
            let f_prime_plain_m = f_prime_plain[0..m].to_vec();
            verify_lookup_permutation_outputs_in_plaintext(
                n,
                m,
                &t_plain,
                &f_plain,
                &t_prime_plain,
                &f_prime_plain_m,
                true,
            );
        }
    }
}

