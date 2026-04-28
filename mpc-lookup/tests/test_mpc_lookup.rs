//! Integration tests for MPC lookup functionality

use ark_ff::UniformRand;
use ark_std::test_rng;
use mpc_algebra::honest_but_curious::MpcField;
use mpc_algebra::Reveal;
use mpc_lookup::{bitonic_sort, secure_oblivious_lookup_permutation, verify_lookup_permutation_outputs_in_plaintext};

type F = ark_bls12_377::Fr;
type MF = MpcField<F>;

#[test]
fn test_batch_secure_cmp_public_values() {
    let pairs = vec![
        (MF::Public(F::from(1u64)), MF::Public(F::from(2u64))),
        (MF::Public(F::from(5u64)), MF::Public(F::from(3u64))),
        (MF::Public(F::from(7u64)), MF::Public(F::from(7u64))),
        (MF::Public(F::from(10u64)), MF::Public(F::from(20u64))),
    ];

    let results = MF::batch_secure_cmp(&pairs);

    assert_eq!(results.len(), 4);
    assert_eq!(results[0].reveal(), F::from(1u64));
    assert_eq!(results[1].reveal(), F::from(0u64));
    assert_eq!(results[2].reveal(), F::from(0u64));
    assert_eq!(results[3].reveal(), F::from(1u64));

    println!("test_batch_secure_cmp_public_values passed");
}
