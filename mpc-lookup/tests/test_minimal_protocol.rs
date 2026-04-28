//! Minimal test to isolate the from_public vs from_add_shared difference

use ark_bls12_377::{Bls12_377, Fr};
use ark_ff::{Field, One, UniformRand};
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
fn test_witness_with_from_public() {
    println!("\n=== Minimal test: from_public ===");
    
    let mut rng = test_rng();
    let max_degree = 48;
    let supported_degree = 15;
    
    // Setup
    let pp = LocalPC::setup(max_degree, None, &mut rng).unwrap();
    let (ck_local, vk_local) = LocalPC::trim(&pp, max_degree, 0, Some(&[supported_degree])).unwrap();
    let ck_mpc = <ark_poly_commit::marlin::marlin_pc::CommitterKey<MpcPairingEngine<Bls12_377>> as Reveal>::from_public(ck_local);
    
    // Create polynomial with from_public
    let coeffs: Vec<Fr> = (0..=supported_degree).map(|_| Fr::rand(&mut rng)).collect();
    let coeffs_mpc: Vec<MpcField<Fr>> = coeffs.iter()
        .map(|&c| MpcField::from_public(c))
        .collect();
    let poly_mpc = DensePolynomial::from_coefficients_vec(coeffs_mpc);
    
    println!("Polynomial is_shared: {}", poly_mpc.is_shared());
    
    // Commit
    let labeled_poly_mpc = LabeledPolynomial::new("test".to_string(), poly_mpc.clone(), Some(supported_degree), None);
    let (commitments_mpc, randomness_mpc) = MpcPC::commit(&ck_mpc, once(&labeled_poly_mpc), Some(&mut rng)).unwrap();
    
    println!("✓ Commitment created");
    
    // Evaluate
    let point = Fr::rand(&mut rng);
    let point_mpc = MpcField::from_public(point);
    let value_mpc = poly_mpc.evaluate(&point_mpc);
    let value_local = value_mpc.reveal();
    
    println!("✓ Evaluated polynomial");
    
    // Compute witness polynomial manually to see what happens
    let divisor_coeffs: Vec<MpcField<Fr>> = vec![
        MpcField::from_public(-point),
        MpcField::from_public(Fr::one()),
    ];
    let divisor = DensePolynomial::from_coefficients_vec(divisor_coeffs);
    
    println!("Divisor is_shared: {}", divisor.is_shared());
    
    // Division should use plain algorithm since poly_mpc.is_shared() == false
    let witness = &poly_mpc / &divisor;
    
    println!("Witness polynomial degree: {}", witness.degree());
    println!("Witness is_shared: {}", witness.is_shared());
    
    // Try to open
    let opening_challenge = MpcField::from_public(Fr::one());
    
    println!("Attempting to open...");
    let result = MpcPC::open(
        &ck_mpc,
        once(&labeled_poly_mpc),
        once(&commitments_mpc[0]),
        &point_mpc,
        opening_challenge,
        once(&randomness_mpc[0]),
        Some(&mut rng),
    );
    
    match result {
        Ok(proof_mpc) => {
            println!("✓ Opening succeeded");
            
            let commitment_local = commitments_mpc[0].commitment().clone().reveal();
            let proof_local = proof_mpc.reveal();
            
            let labeled_comm_local = LabeledCommitment::new("test".to_string(), commitment_local, Some(supported_degree));
            
            let check_result = LocalPC::check(
                &vk_local,
                once(&labeled_comm_local),
                &point,
                once(value_local),
                &proof_local,
                Fr::one(),
                None,
            ).expect("Check returned error");
            
            println!("Pairing check result: {}", check_result);
            
            if check_result {
                println!("✓✓✓ TEST PASSED with from_public");
            } else {
                println!("✗✗✗ TEST FAILED with from_public - pairing check failed");
            }
            
            assert!(check_result, "Pairing check should pass");
        }
        Err(e) => {
            println!("✗ Opening failed: {:?}", e);
            panic!("Opening should succeed");
        }
    }
}

#[test]
fn test_witness_with_from_add_shared() {
    println!("\n=== Minimal test: from_add_shared ===");
    
    let mut rng = test_rng();
    let max_degree = 48;
    let supported_degree = 15;
    
    // Setup
    let pp = LocalPC::setup(max_degree, None, &mut rng).unwrap();
    let (ck_local, vk_local) = LocalPC::trim(&pp, max_degree, 0, Some(&[supported_degree])).unwrap();
    let ck_mpc = <ark_poly_commit::marlin::marlin_pc::CommitterKey<MpcPairingEngine<Bls12_377>> as Reveal>::from_public(ck_local);
    
    // Create polynomial with from_add_shared
    let coeffs: Vec<Fr> = (0..=supported_degree).map(|_| Fr::rand(&mut rng)).collect();
    let coeffs_mpc: Vec<MpcField<Fr>> = coeffs.iter()
        .map(|&c| MpcField::from_add_shared(c))
        .collect();
    let poly_mpc = DensePolynomial::from_coefficients_vec(coeffs_mpc);
    
    println!("Polynomial is_shared: {}", poly_mpc.is_shared());
    
    // Commit
    let labeled_poly_mpc = LabeledPolynomial::new("test".to_string(), poly_mpc.clone(), Some(supported_degree), None);
    let (commitments_mpc, randomness_mpc) = MpcPC::commit(&ck_mpc, once(&labeled_poly_mpc), Some(&mut rng)).unwrap();
    
    println!("✓ Commitment created");
    
    // Evaluate
    let point = Fr::rand(&mut rng);
    let point_mpc = MpcField::from_public(point);
    let value_mpc = poly_mpc.evaluate(&point_mpc);
    let value_local = value_mpc.reveal();
    
    println!("✓ Evaluated polynomial (revealed value may be 0 in single-party mode)");
    
    // Compute witness polynomial manually to see what happens
    let divisor_coeffs: Vec<MpcField<Fr>> = vec![
        MpcField::from_public(-point),
        MpcField::from_public(Fr::one()),
    ];
    let divisor = DensePolynomial::from_coefficients_vec(divisor_coeffs);
    
    println!("Divisor is_shared: {}", divisor.is_shared());
    
    // Division should use MPC algorithm since poly_mpc.is_shared() == true
    let witness = &poly_mpc / &divisor;
    
    println!("Witness polynomial degree: {}", witness.degree());
    println!("Witness is_shared: {}", witness.is_shared());
    
    // Try to open
    let opening_challenge = MpcField::from_public(Fr::one());
    
    println!("Attempting to open...");
    let result = MpcPC::open(
        &ck_mpc,
        once(&labeled_poly_mpc),
        once(&commitments_mpc[0]),
        &point_mpc,
        opening_challenge,
        once(&randomness_mpc[0]),
        Some(&mut rng),
    );
    
    match result {
        Ok(proof_mpc) => {
            println!("✓ Opening succeeded");
            
            let commitment_local = commitments_mpc[0].commitment().clone().reveal();
            let proof_local = proof_mpc.reveal();
            
            let labeled_comm_local = LabeledCommitment::new("test".to_string(), commitment_local, Some(supported_degree));
            
            let check_result = LocalPC::check(
                &vk_local,
                once(&labeled_comm_local),
                &point,
                once(value_local),
                &proof_local,
                Fr::one(),
                None,
            ).expect("Check returned error");
            
            println!("Pairing check result: {}", check_result);
            
            if check_result {
                println!("✓✓✓ TEST PASSED with from_add_shared");
            } else {
                println!("✗✗✗ TEST FAILED with from_add_shared - pairing check failed");
            }
            
            assert!(check_result, "Pairing check should pass");
        }
        Err(e) => {
            println!("✗ Opening failed: {:?}", e);
            panic!("Opening should succeed");
        }
    }
}
