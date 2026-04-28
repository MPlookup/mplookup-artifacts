//! Debug test to compare local vs MPC witness polynomial computation

use ark_bls12_377::{Bls12_377, Fr};
use ark_ff::{Field, One, UniformRand, Zero};
use ark_poly::univariate::DensePolynomial;
use ark_poly::{Polynomial, UVPolynomial};
use ark_poly_commit::marlin::marlin_pc::MarlinKZG10;
use ark_poly_commit::{LabeledCommitment, LabeledPolynomial, PolynomialCommitment};
use ark_std::test_rng;
use mpc_algebra::honest_but_curious::{MpcField, MpcPairingEngine};
use mpc_algebra::Reveal;
use std::iter::once;

type LocalPC = MarlinKZG10<Bls12_377, DensePolynomial<Fr>>;
type MpcPC = MarlinKZG10<MpcPairingEngine<Bls12_377>, DensePolynomial<MpcField<Fr>>>;

#[test]
fn test_compare_local_vs_mpc() {
    println!("\n=== Comparing Local vs MPC Polynomial Opening ===");
    
    let mut rng = test_rng();
    let max_degree = 48;
    let supported_degree = 15;
    
    // Create a simple polynomial with known coefficients for easier debugging
    let coeffs: Vec<Fr> = (0..=supported_degree).map(|i| Fr::from((i + 1) as u64)).collect();
    let poly_local = DensePolynomial::from_coefficients_vec(coeffs.clone());
    
    println!("Polynomial coefficients: [1, 2, 3, ..., {}]", supported_degree + 1);
    
    // Setup
    let pp = LocalPC::setup(max_degree, None, &mut rng).unwrap();
    let (ck_local, vk_local) = LocalPC::trim(&pp, max_degree, 0, Some(&[supported_degree])).unwrap();
    
    // Test point
    let point_local = Fr::from(7u64);
    let value_local = poly_local.evaluate(&point_local);
    
    println!("Evaluation point: 7");
    println!("Evaluation value: {}", value_local);
    
    // === LOCAL VERSION ===
    println!("\n--- LOCAL POLYNOMIAL COMMITMENT ---");
    
    let labeled_poly_local = LabeledPolynomial::new("test".to_string(), poly_local.clone(), Some(supported_degree), None);
    let (commitments_local, randomness_local) = LocalPC::commit(&ck_local, once(&labeled_poly_local), Some(&mut rng)).unwrap();
    
    println!("✓ Local commitment created");
    
    let proof_local = LocalPC::open(
        &ck_local,
        once(&labeled_poly_local),
        once(&commitments_local[0]),
        &point_local,
        Fr::one(),
        once(&randomness_local[0]),
        Some(&mut rng),
    ).unwrap();
    
    println!("✓ Local opening proof generated");
    
    let labeled_comm_local = LabeledCommitment::new(
        "test".to_string(),
        commitments_local[0].commitment().clone(),
        Some(supported_degree),
    );
    
    let result_local = LocalPC::check(
        &vk_local,
        once(&labeled_comm_local),
        &point_local,
        once(value_local),
        &proof_local,
        Fr::one(),
        None,
    ).unwrap();
    
    println!("Local pairing check result: {}", result_local);
    assert!(result_local, "Local version should pass");
    
    // === MPC VERSION ===
    println!("\n--- MPC POLYNOMIAL COMMITMENT ---");
    
    let ck_mpc = <ark_poly_commit::marlin::marlin_pc::CommitterKey<MpcPairingEngine<Bls12_377>> as Reveal>::from_public(ck_local);
    
    // Try with from_public
    println!("\n** Testing with from_public **");
    let coeffs_mpc_public: Vec<MpcField<Fr>> = coeffs.iter().map(|&c| MpcField::from_public(c)).collect();
    let poly_mpc_public = DensePolynomial::from_coefficients_vec(coeffs_mpc_public);
    
    let point_mpc = MpcField::from_public(point_local);
    let value_mpc = poly_mpc_public.evaluate(&point_mpc);
    let value_mpc_revealed = value_mpc.reveal();
    
    println!("MPC evaluation (from_public): {}", value_mpc_revealed);
    println!("Match with local: {}", value_mpc_revealed == value_local);
    
    let labeled_poly_mpc = LabeledPolynomial::new("test".to_string(), poly_mpc_public.clone(), Some(supported_degree), None);
    let (commitments_mpc, randomness_mpc) = MpcPC::commit(&ck_mpc, once(&labeled_poly_mpc), Some(&mut rng)).unwrap();
    
    println!("✓ MPC commitment created (from_public)");
    
    let proof_mpc = MpcPC::open(
        &ck_mpc,
        once(&labeled_poly_mpc),
        once(&commitments_mpc[0]),
        &point_mpc,
        MpcField::from_public(Fr::one()),
        once(&randomness_mpc[0]),
        Some(&mut rng),
    ).unwrap();
    
    println!("✓ MPC opening proof generated (from_public)");
    
    // Reveal and check
    let commitment_mpc_revealed = commitments_mpc[0].commitment().clone().reveal();
    let proof_mpc_revealed = proof_mpc.reveal();
    
    let labeled_comm_mpc = LabeledCommitment::new(
        "test".to_string(),
        commitment_mpc_revealed,
        Some(supported_degree),
    );
    
    let result_mpc = LocalPC::check(
        &vk_local,
        once(&labeled_comm_mpc),
        &point_local,
        once(value_mpc_revealed),
        &proof_mpc_revealed,
        Fr::one(),
        None,
    ).unwrap();
    
    println!("MPC pairing check result (from_public): {}", result_mpc);
    
    // Try with from_add_shared
    println!("\n** Testing with from_add_shared **");
    let coeffs_mpc_shared: Vec<MpcField<Fr>> = coeffs.iter().map(|&c| MpcField::from_add_shared(c)).collect();
    let poly_mpc_shared = DensePolynomial::from_coefficients_vec(coeffs_mpc_shared);
    
    let value_mpc_shared = poly_mpc_shared.evaluate(&point_mpc);
    let value_mpc_shared_revealed = value_mpc_shared.reveal();
    
    println!("MPC evaluation (from_add_shared): {}", value_mpc_shared_revealed);
    println!("Match with local: {}", value_mpc_shared_revealed == value_local);
    
    let labeled_poly_mpc_shared = LabeledPolynomial::new("test".to_string(), poly_mpc_shared.clone(), Some(supported_degree), None);
    let (commitments_mpc_shared, randomness_mpc_shared) = MpcPC::commit(&ck_mpc, once(&labeled_poly_mpc_shared), Some(&mut rng)).unwrap();
    
    println!("✓ MPC commitment created (from_add_shared)");
    
    let proof_mpc_shared = MpcPC::open(
        &ck_mpc,
        once(&labeled_poly_mpc_shared),
        once(&commitments_mpc_shared[0]),
        &point_mpc,
        MpcField::from_public(Fr::one()),
        once(&randomness_mpc_shared[0]),
        Some(&mut rng),
    ).unwrap();
    
    println!("✓ MPC opening proof generated (from_add_shared)");
    
    // Reveal and check
    let commitment_mpc_shared_revealed = commitments_mpc_shared[0].commitment().clone().reveal();
    let proof_mpc_shared_revealed = proof_mpc_shared.reveal();
    
    let labeled_comm_mpc_shared = LabeledCommitment::new(
        "test".to_string(),
        commitment_mpc_shared_revealed,
        Some(supported_degree),
    );
    
    let result_mpc_shared = LocalPC::check(
        &vk_local,
        once(&labeled_comm_mpc_shared),
        &point_local,
        once(value_mpc_shared_revealed),
        &proof_mpc_shared_revealed,
        Fr::one(),
        None,
    ).unwrap();
    
    println!("MPC pairing check result (from_add_shared): {}", result_mpc_shared);
    
    println!("\n=== SUMMARY ===");
    println!("Local: {}", if result_local { "✓ PASS" } else { "✗ FAIL" });
    println!("MPC (from_public): {}", if result_mpc { "✓ PASS" } else { "✗ FAIL" });
    println!("MPC (from_add_shared): {}", if result_mpc_shared { "✓ PASS" } else { "✗ FAIL" });
    
    assert!(result_mpc || result_mpc_shared, "At least one MPC version should pass");
}
