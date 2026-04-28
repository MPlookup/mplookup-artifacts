//! Incremental tests to isolate MPC MarlinKZG10 verification failure
//!
//! Following the exact pattern from mpc-lookup/client.rs:
//! 1. Local setup (non-MPC)
//! 2. Convert to MPC types
//! 3. MPC prove
//! 4. Reveal and local verify
//!
//! Each test adds one more component to isolate where verification fails.

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

/// Test 1: Verify local MarlinKZG10 setup works
#[test]
fn test_1_local_setup_works() {
    println!("\n=== Test 1: Local MarlinKZG10 Setup ===");
    
    let mut rng = test_rng();
    let max_degree = 48; // 3 * domain_size for domain_size=16
    let supported_degree = 15; // domain_size - 1 for domain_size=16
    
    println!("Setting up with max_degree={}, supported_degrees=[{}]", max_degree, supported_degree);
    
    // This is how mpc-lookup does setup (client.rs line 161)
    let pp = LocalPC::setup(max_degree, None, &mut rng)
        .expect("Local setup failed");
    
    // This matches client.rs line 169
    let (ck, vk) = LocalPC::trim(&pp, max_degree, 0, Some(&[supported_degree]))
        .expect("Local trim failed");
    
    println!("✓ Local setup succeeded");
    println!("✓ Committer key created");
    println!("✓ Verifier key created");
}

/// Test 2: Commit to polynomial with MPC types using local setup
#[test]
fn test_2_mpc_commit_with_local_setup() {
    println!("\n=== Test 2: MPC Commit with Local Setup ===");
    
    let mut rng = test_rng();
    let max_degree = 48;
    let supported_degree = 15;
    
    // Step 1: Local setup (same as Test 1)
    let pp = LocalPC::setup(max_degree, None, &mut rng).unwrap();
    let (ck_local, vk_local) = LocalPC::trim(&pp, max_degree, 0, Some(&[supported_degree])).unwrap();
    
    println!("✓ Local setup complete");
    
    // Step 2: Convert to MPC types (matches client.rs line 176)
    let ck_mpc = <ark_poly_commit::marlin::marlin_pc::CommitterKey<MpcPairingEngine<Bls12_377>> as Reveal>::from_public(ck_local);
    
    println!("✓ Converted committer key to MPC type");
    
    // Step 3: Create MPC polynomial
    let coeffs: Vec<Fr> = (0..=supported_degree).map(|_| Fr::rand(&mut rng)).collect();
    let poly_local = DensePolynomial::from_coefficients_vec(coeffs.clone());
    
    // Share the polynomial coefficients (from_public is fine for test 2 since it doesn't do division)
    let coeffs_mpc: Vec<MpcField<Fr>> = coeffs.iter().map(|&c| MpcField::from_public(c)).collect();
    let poly_mpc = DensePolynomial::from_coefficients_vec(coeffs_mpc);
    
    println!("✓ Created MPC polynomial of degree {}", poly_mpc.degree());
    
    // Step 4: Commit with MPC
    let labeled_poly = LabeledPolynomial::new("test".to_string(), poly_mpc, Some(supported_degree), None);
    
    let (commitments, randomness) = MpcPC::commit(&ck_mpc, once(&labeled_poly), None)
        .expect("MPC commit failed");
    
    println!("✓ MPC commitment succeeded");
    println!("  Generated {} commitment(s)", commitments.len());
    
    // Step 5: Reveal commitment
    let commitment_local = commitments[0].commitment().clone().reveal();
    
    println!("✓ Revealed commitment to local type");
    println!("✓ Test 2 PASSED: MPC commit with local setup works!");
}

/// Test 3: Open polynomial with MPC, verify with local
#[test]
fn test_3_mpc_open_local_verify() {
    println!("\n=== Test 3: MPC Open + Local Verify ===");
    
    let mut rng = test_rng();
    let max_degree = 48;
    let supported_degree = 15;
    
    // Setup
    let pp = LocalPC::setup(max_degree, None, &mut rng).unwrap();
    let (ck_local, vk_local) = LocalPC::trim(&pp, max_degree, 0, Some(&[supported_degree])).unwrap();
    let ck_mpc = <ark_poly_commit::marlin::marlin_pc::CommitterKey<MpcPairingEngine<Bls12_377>> as Reveal>::from_public(ck_local);
    
    println!("✓ Setup complete");
    
    // Create and commit to polynomial
    let coeffs: Vec<Fr> = (0..=supported_degree).map(|_| Fr::rand(&mut rng)).collect();
    let poly_local = DensePolynomial::from_coefficients_vec(coeffs.clone());

    let coeffs_mpc: Vec<MpcField<Fr>> = coeffs.iter().map(|&c| MpcField::from_public(c)).collect();
    let poly_mpc = DensePolynomial::from_coefficients_vec(coeffs_mpc);
    
    let labeled_poly_mpc = LabeledPolynomial::new("test".to_string(), poly_mpc.clone(), Some(supported_degree), None);
    let (commitments_mpc, randomness_mpc) = MpcPC::commit(&ck_mpc, once(&labeled_poly_mpc), Some(&mut rng)).unwrap();
    
    println!("✓ MPC commitment created");
    
    // Choose random evaluation point
    let point = Fr::rand(&mut rng);
    let point_mpc = MpcField::from_public(point);
    
    // Evaluate polynomial
    let value_mpc = poly_mpc.evaluate(&point_mpc);
    let value_local = value_mpc.reveal();
    
    println!("✓ Evaluated polynomial at random point");
    println!("  Evaluation: {}", value_local);

    // Open with MPC (THIS IS THE CRITICAL STEP)
    // Following mpc-plonk pattern: PC::open with F::one() for single commitment
    let opening_challenge = MpcField::from_public(Fr::one());
    
    let proof_mpc = MpcPC::open(
        &ck_mpc,
        once(&labeled_poly_mpc),
        once(&commitments_mpc[0]),
        &point_mpc,
        opening_challenge,
        once(&randomness_mpc[0]),
        Some(&mut rng),
    ).expect("MPC open failed");
    
    println!("✓ MPC opening proof generated");
    
    // Reveal everything to local types
    let commitment_local = commitments_mpc[0].commitment().clone().reveal();
    let proof_local = proof_mpc.reveal();
    
    println!("✓ Revealed proof to local type");
    
    // Verify with LOCAL PC (THIS IS WHERE IT MIGHT FAIL)
    let labeled_comm_local = LabeledCommitment::new(
        "test".to_string(),
        commitment_local,
        Some(supported_degree),
    );
    
    let result = LocalPC::check(
        &vk_local,
        once(&labeled_comm_local),
        &point,
        once(value_local),
        &proof_local,
        Fr::one(), // Same opening challenge as used in open()
        None,
    ).expect("Check returned error");
    
    println!("\n=== VERIFICATION RESULT ===");
    println!("Pairing check result: {}", result);
    
    if result {
        println!("✓✓✓ Test 3 PASSED: MPC prove + local verify WORKS! ✓✓✓");
    } else {
        println!("✗✗✗ Test 3 FAILED: Pairing check returns false ✗✗✗");
        println!("This is the same failure as in mpc-lookup!");
    }
    
    assert!(result, "Pairing check should pass");
}

/// Test 4: Full flow with hiding commitments
#[test]
fn test_4_mpc_open_local_verify_hiding() {
    println!("\n=== Test 4: MPC Open + Local Verify (Hiding Commitments) ===");
    
    let mut rng = test_rng();
    let max_degree = 48;
    let supported_degree = 15;
    
    // Setup
    let pp = LocalPC::setup(max_degree, None, &mut rng).unwrap();
    let (ck_local, vk_local) = LocalPC::trim(&pp, max_degree, 0, Some(&[supported_degree])).unwrap();
    let ck_mpc = <ark_poly_commit::marlin::marlin_pc::CommitterKey<MpcPairingEngine<Bls12_377>> as Reveal>::from_public(ck_local);
    
    println!("✓ Setup complete");
    
    // Create polynomial
    let coeffs: Vec<Fr> = (0..=supported_degree).map(|_| Fr::rand(&mut rng)).collect();
    let coeffs_mpc: Vec<MpcField<Fr>> = coeffs.iter().map(|&c| MpcField::from_public(c)).collect();
    let poly_mpc = DensePolynomial::from_coefficients_vec(coeffs_mpc);
    
    let labeled_poly_mpc = LabeledPolynomial::new("test".to_string(), poly_mpc.clone(), Some(supported_degree), None);
    
    // Commit with RNG (hiding commitment)
    let (commitments_mpc, randomness_mpc) = MpcPC::commit(&ck_mpc, once(&labeled_poly_mpc), Some(&mut rng)).unwrap();
    
    println!("✓ MPC hiding commitment created");
    
    // Evaluate and open
    let point = Fr::rand(&mut rng);
    let point_mpc = MpcField::from_public(point);
    let value_mpc = poly_mpc.evaluate(&point_mpc);
    let value_local = value_mpc.reveal();
    
    let opening_challenge = MpcField::from_public(Fr::one());
    
    let proof_mpc = MpcPC::open(
        &ck_mpc,
        once(&labeled_poly_mpc),
        once(&commitments_mpc[0]),
        &point_mpc,
        opening_challenge,
        once(&randomness_mpc[0]),
        Some(&mut rng), // Provide RNG for opening too
    ).expect("MPC open failed");
    
    println!("✓ MPC opening proof generated");
    
    // Reveal and verify
    let commitment_local = commitments_mpc[0].commitment().clone().reveal();
    let proof_local = proof_mpc.reveal();
    
    let labeled_comm_local = LabeledCommitment::new("test".to_string(), commitment_local, Some(supported_degree));
    
    let result = LocalPC::check(
        &vk_local,
        once(&labeled_comm_local),
        &point,
        once(value_local),
        &proof_local,
        Fr::one(),
        None,
    ).expect("Check returned error");
    
    println!("\n=== VERIFICATION RESULT (Hiding) ===");
    println!("Pairing check result: {}", result);
    
    if result {
        println!("✓✓✓ Test 4 PASSED: Hiding commitments work! ✓✓✓");
    } else {
        println!("✗✗✗ Test 4 FAILED: Pairing check returns false even with hiding ✗✗✗");
    }
    
    assert!(result, "Pairing check should pass with hiding commitments");
}
