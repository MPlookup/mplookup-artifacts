// Test MarlinKZG10 polynomial commitment with plain Fr to verify basic functionality
// This isolates whether the pairing check works correctly before MPC integration

use ark_bls12_377::{Bls12_377, Fr};
use ark_ff::{One, UniformRand};
use ark_poly::{
    univariate::DensePolynomial, Polynomial, UVPolynomial,
};
use ark_poly_commit::{
    marlin_pc::MarlinKZG10, LabeledPolynomial, PolynomialCommitment,
};
use ark_std::test_rng;

type PC = MarlinKZG10<Bls12_377, DensePolynomial<Fr>>;

#[test]
fn test_marlin_kzg10_basic_commit_and_open() {
    println!("\n=== Testing Basic MarlinKZG10 with Fr ===");
    
    let rng = &mut test_rng();
    
    // Setup parameters
    let max_degree = 20;
    let supported_degree = 15;
    
    println!("Setting up universal params with max_degree={}", max_degree);
    let pp = PC::setup(max_degree, None, rng).expect("Setup failed");
    
    println!("Trimming to get committer and verifier keys with supported_degree={}", supported_degree);
    let (ck, vk) = PC::trim(
        &pp,
        max_degree,
        0, // no hiding bound
        Some(&[supported_degree]),
    ).expect("Trim failed");
    
    // Create a simple polynomial of degree 15
    println!("\nCreating polynomial of degree {}", supported_degree);
    let coeffs: Vec<Fr> = (0..=supported_degree)
        .map(|_| Fr::rand(rng))
        .collect();
    let poly = DensePolynomial::from_coefficients_vec(coeffs);
    
    println!("Polynomial degree: {}", poly.degree());
    
    // Create labeled polynomial with degree bound
    let labeled_poly = LabeledPolynomial::new(
        "test_poly".to_string(),
        poly.clone(),
        Some(supported_degree),
        None, // NO hiding bound - to avoid HidingBoundTooLarge error
    );
    
    // Commit to the polynomial
    println!("\nCommitting to polynomial...");
    let (commitments, randomness) = PC::commit(
        &ck,
        std::iter::once(&labeled_poly),
        Some(rng),
    ).expect("Commit failed");
    
    let commitment = commitments.first().expect("No commitment returned");
    println!("Commitment created successfully");
    
    // Choose a random evaluation point
    let point = Fr::rand(rng);
    let value = poly.evaluate(&point);
    
    println!("\nEvaluating polynomial at random point");
    println!("Point: {:?}", point);
    println!("Value: {:?}", value);
    
    // Create opening proof
    println!("\nGenerating opening proof...");
    let proof = PC::open(
        &ck,
        std::iter::once(&labeled_poly),
        std::iter::once(commitment),
        &point,
        Fr::one(), // opening challenge (must be 1 for single commitment)
        std::iter::once(&randomness[0]),
        Some(rng),
    ).expect("Open failed");
    
    println!("Opening proof generated successfully");
    
    // Verify the opening
    println!("\nVerifying opening proof...");
    let result = PC::check(
        &vk,
        std::iter::once(commitment),
        &point,
        std::iter::once(value),
        &proof,
        Fr::one(), // opening challenge (must be 1 for single commitment)
        None,
    ).expect("Check failed with error");
    
    println!("Verification result: {}", result);
    
    if result {
        println!("\n✓ SUCCESS: MarlinKZG10 pairing check PASSED with Fr");
        println!("This confirms MarlinKZG10 itself works correctly.");
        println!("The issue in mpc-lookup must be with MPC integration.");
    } else {
        println!("\n✗ FAILURE: MarlinKZG10 pairing check FAILED with Fr");
        println!("This would indicate a fundamental issue with MarlinKZG10 setup.");
    }
    
    assert!(result, "MarlinKZG10 verification should pass with Fr");
}

#[test]
fn test_marlin_kzg10_non_hiding() {
    println!("\n=== Testing MarlinKZG10 with Non-Hiding Commitments ===");
    
    let rng = &mut test_rng();
    
    let max_degree = 20;
    let supported_degree = 15;
    
    let pp = PC::setup(max_degree, None, rng).expect("Setup failed");
    let (ck, vk) = PC::trim(&pp, max_degree, 0, Some(&[supported_degree]))
        .expect("Trim failed");
    
    // Create polynomial
    let coeffs: Vec<Fr> = (0..=supported_degree)
        .map(|_| Fr::rand(rng))
        .collect();
    let poly = DensePolynomial::from_coefficients_vec(coeffs);
    
    let labeled_poly = LabeledPolynomial::new(
        "test_poly".to_string(),
        poly.clone(),
        Some(supported_degree),
        None, // NO hiding bound
    );
    
    // Commit WITHOUT rng (non-hiding)
    println!("Committing without RNG (non-hiding)...");
    let (commitments, randomness) = PC::commit(
        &ck,
        std::iter::once(&labeled_poly),
        None, // No RNG = non-hiding
    ).expect("Commit failed");
    
    let commitment = commitments.first().expect("No commitment returned");
    
    let point = Fr::rand(rng);
    let value = poly.evaluate(&point);
    
    // Open WITH rng
    println!("Opening with RNG...");
    let proof = PC::open(
        &ck,
        std::iter::once(&labeled_poly),
        std::iter::once(commitment),
        &point,
        Fr::one(),
        std::iter::once(&randomness[0]),
        Some(rng), // Provide RNG for opening
    ).expect("Open failed");
    
    // Verify
    let result = PC::check(
        &vk,
        std::iter::once(commitment),
        &point,
        std::iter::once(value),
        &proof,
        Fr::one(),
        None,
    ).expect("Check failed with error");
    
    println!("Non-hiding verification result: {}", result);
    
    if result {
        println!("✓ Non-hiding commitments work correctly");
    } else {
        println!("✗ Non-hiding commitments failed - this matches our mpc-lookup issue!");
    }
    
    assert!(result, "Non-hiding verification should pass");
}
