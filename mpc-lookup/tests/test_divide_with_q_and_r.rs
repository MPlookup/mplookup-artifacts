//! Test for divide_with_q_and_r function with MPC fields
//!
//! This test validates the correctness of polynomial division with both:
//! (a) Public coefficients (MpcField::from_public)
//! (b) Shared coefficients (MpcField::from_add_shared)
//!
//! The divisor is always public.

use ark_bls12_377::Fr;
use ark_ff::{Field, UniformRand, Zero};
use ark_poly::univariate::{DenseOrSparsePolynomial, DensePolynomial};
use ark_poly::{Polynomial, UVPolynomial};
use ark_std::test_rng;
use mpc_algebra::honest_but_curious::MpcField;
use mpc_algebra::Reveal;

/// Test divide_with_q_and_r with Public MpcField coefficients
#[test]
fn test_divide_with_q_and_r_public() {
    println!("\n=== Test: divide_with_q_and_r with Public coefficients ===");
    
    let _rng = test_rng();
    
    // Create a dividend polynomial: f(x) = 5x^3 + 4x^2 + 3x + 2
    let dividend_coeffs: Vec<Fr> = vec![
        Fr::from(2u64),
        Fr::from(3u64),
        Fr::from(4u64),
        Fr::from(5u64),
    ];
    
    // Create a divisor polynomial: g(x) = x + 1
    let divisor_coeffs: Vec<Fr> = vec![
        Fr::from(1u64),
        Fr::from(1u64),
    ];
    
    println!("Dividend polynomial degree: {}", dividend_coeffs.len() - 1);
    println!("Divisor polynomial degree: {}", divisor_coeffs.len() - 1);
    
    // Convert dividend to MPC field with PUBLIC sharing
    let dividend_coeffs_mpc: Vec<MpcField<Fr>> = dividend_coeffs
        .iter()
        .map(|&c| MpcField::from_public(c))
        .collect();
    let dividend_mpc = DensePolynomial::from_coefficients_vec(dividend_coeffs_mpc);
    
    // Divisor is also MPC field but with public values
    let divisor_coeffs_mpc: Vec<MpcField<Fr>> = divisor_coeffs
        .iter()
        .map(|&c| MpcField::from_public(c))
        .collect();
    let divisor_mpc = DensePolynomial::from_coefficients_vec(divisor_coeffs_mpc);
    
    // Convert to DenseOrSparsePolynomial for the function call
    let dividend_sparse: DenseOrSparsePolynomial<MpcField<Fr>> = (&dividend_mpc).into();
    let divisor_sparse: DenseOrSparsePolynomial<MpcField<Fr>> = (&divisor_mpc).into();
    
    // Perform division
    println!("Performing division with PUBLIC coefficients...");
    let result = dividend_sparse.divide_with_q_and_r(&divisor_sparse);
    
    assert!(result.is_some(), "Division should succeed");
    let (quotient_mpc, remainder_mpc) = result.unwrap();
    
    println!("✓ Division succeeded");
    println!("  Quotient degree: {}", quotient_mpc.degree());
    println!("  Remainder degree: {}", remainder_mpc.degree());
    
    // Reveal the MPC results
    let quotient = DensePolynomial::from_coefficients_vec(
        quotient_mpc.coeffs.iter().map(|c| c.reveal()).collect()
    );
    let remainder = DensePolynomial::from_coefficients_vec(
        remainder_mpc.coeffs.iter().map(|c| c.reveal()).collect()
    );
    let divisor = DensePolynomial::from_coefficients_vec(divisor_coeffs.clone());
    
    // Verify: dividend = quotient * divisor + remainder
    let product = &quotient * &divisor;
    let reconstructed = &product + &remainder;
    let original = DensePolynomial::from_coefficients_vec(dividend_coeffs.clone());
    
    println!("Verifying: dividend = quotient * divisor + remainder");
    
    // Compare coefficients, handling potential trailing zeros
    let max_len = reconstructed.coeffs.len().max(original.coeffs.len());
    for i in 0..max_len {
        let reconstructed_coeff = reconstructed.coeffs.get(i).copied().unwrap_or_else(Fr::zero);
        let original_coeff = original.coeffs.get(i).copied().unwrap_or_else(Fr::zero);
        assert_eq!(
            reconstructed_coeff, original_coeff,
            "Coefficient mismatch at position {}", i
        );
    }
    
    println!("✓✓✓ Test PASSED: divide_with_q_and_r with PUBLIC coefficients works correctly! ✓✓✓");
}

/// Test divide_with_q_and_r with Shared MpcField coefficients
#[test]
fn test_divide_with_q_and_r_shared() {
    println!("\n=== Test: divide_with_q_and_r with Shared coefficients ===");
    
    let _rng = test_rng();
    
    // Create a dividend polynomial: f(x) = 5x^3 + 4x^2 + 3x + 2
    let dividend_coeffs: Vec<Fr> = vec![
        Fr::from(2u64),
        Fr::from(3u64),
        Fr::from(4u64),
        Fr::from(5u64),
    ];
    
    // Create a divisor polynomial: g(x) = x + 1
    let divisor_coeffs: Vec<Fr> = vec![
        Fr::from(1u64),
        Fr::from(1u64),
    ];
    
    println!("Dividend polynomial degree: {}", dividend_coeffs.len() - 1);
    println!("Divisor polynomial degree: {}", divisor_coeffs.len() - 1);
    
    // Convert dividend to MPC field with SHARED (additive) sharing
    let dividend_coeffs_mpc: Vec<MpcField<Fr>> = dividend_coeffs
        .iter()
        .map(|&c| MpcField::from_add_shared(c))
        .collect();
    let dividend_mpc = DensePolynomial::from_coefficients_vec(dividend_coeffs_mpc);
    
    // Divisor is also MPC field but with public values
    let divisor_coeffs_mpc: Vec<MpcField<Fr>> = divisor_coeffs
        .iter()
        .map(|&c| MpcField::from_public(c))
        .collect();
    let divisor_mpc = DensePolynomial::from_coefficients_vec(divisor_coeffs_mpc);
    
    // Convert to DenseOrSparsePolynomial for the function call
    let dividend_sparse: DenseOrSparsePolynomial<MpcField<Fr>> = (&dividend_mpc).into();
    let divisor_sparse: DenseOrSparsePolynomial<MpcField<Fr>> = (&divisor_mpc).into();
    
    // Perform division
    println!("Performing division with SHARED coefficients...");
    let result = dividend_sparse.divide_with_q_and_r(&divisor_sparse);
    
    assert!(result.is_some(), "Division should succeed");
    let (quotient_mpc, remainder_mpc) = result.unwrap();
    
    println!("✓ Division succeeded");
    println!("  Quotient degree: {}", quotient_mpc.degree());
    println!("  Remainder degree: {}", remainder_mpc.degree());
    
    // Reveal the MPC results
    let quotient = DensePolynomial::from_coefficients_vec(
        quotient_mpc.coeffs.iter().map(|c| c.reveal()).collect()
    );
    let remainder = DensePolynomial::from_coefficients_vec(
        remainder_mpc.coeffs.iter().map(|c| c.reveal()).collect()
    );
    let divisor = DensePolynomial::from_coefficients_vec(divisor_coeffs.clone());
    
    // Verify: dividend = quotient * divisor + remainder
    let product = &quotient * &divisor;
    let reconstructed = &product + &remainder;
    let original = DensePolynomial::from_coefficients_vec(dividend_coeffs.clone());
    
    println!("Verifying: dividend = quotient * divisor + remainder");
    
    // Compare coefficients, handling potential trailing zeros
    let max_len = reconstructed.coeffs.len().max(original.coeffs.len());
    for i in 0..max_len {
        let reconstructed_coeff = reconstructed.coeffs.get(i).copied().unwrap_or_else(Fr::zero);
        let original_coeff = original.coeffs.get(i).copied().unwrap_or_else(Fr::zero);
        assert_eq!(
            reconstructed_coeff, original_coeff,
            "Coefficient mismatch at position {}", i
        );
    }
    
    println!("✓✓✓ Test PASSED: divide_with_q_and_r with SHARED coefficients works correctly! ✓✓✓");
}

/// Test with random polynomials to ensure robustness
#[test]
fn test_divide_with_q_and_r_random_public() {
    println!("\n=== Test: divide_with_q_and_r with random Public coefficients ===");
    
    let mut rng = test_rng();
    
    // Create random dividend polynomial of degree 10
    let dividend_degree = 10;
    let dividend_coeffs: Vec<Fr> = (0..=dividend_degree)
        .map(|_| Fr::rand(&mut rng))
        .collect();
    
    // Create random divisor polynomial of degree 3
    let divisor_degree = 3;
    let divisor_coeffs: Vec<Fr> = (0..=divisor_degree)
        .map(|_| Fr::rand(&mut rng))
        .collect();
    
    println!("Dividend polynomial degree: {}", dividend_degree);
    println!("Divisor polynomial degree: {}", divisor_degree);
    
    // Convert dividend to MPC field with PUBLIC sharing
    let dividend_coeffs_mpc: Vec<MpcField<Fr>> = dividend_coeffs
        .iter()
        .map(|&c| MpcField::from_public(c))
        .collect();
    let dividend_mpc = DensePolynomial::from_coefficients_vec(dividend_coeffs_mpc);
    
    // Divisor is also MPC field but with public values
    let divisor_coeffs_mpc: Vec<MpcField<Fr>> = divisor_coeffs
        .iter()
        .map(|&c| MpcField::from_public(c))
        .collect();
    let divisor_mpc = DensePolynomial::from_coefficients_vec(divisor_coeffs_mpc);
    
    // Convert to DenseOrSparsePolynomial for the function call
    let dividend_sparse: DenseOrSparsePolynomial<MpcField<Fr>> = (&dividend_mpc).into();
    let divisor_sparse: DenseOrSparsePolynomial<MpcField<Fr>> = (&divisor_mpc).into();
    
    // Perform division
    println!("Performing division with random PUBLIC coefficients...");
    let result = dividend_sparse.divide_with_q_and_r(&divisor_sparse);
    
    assert!(result.is_some(), "Division should succeed");
    let (quotient_mpc, remainder_mpc) = result.unwrap();
    
    println!("✓ Division succeeded");
    println!("  Quotient degree: {}", quotient_mpc.degree());
    println!("  Remainder degree: {}", remainder_mpc.degree());
    
    // Verify remainder degree is less than divisor degree
    assert!(
        remainder_mpc.degree() < divisor_mpc.degree(),
        "Remainder degree should be less than divisor degree"
    );
    
    // Reveal the MPC results
    let quotient = DensePolynomial::from_coefficients_vec(
        quotient_mpc.coeffs.iter().map(|c| c.reveal()).collect()
    );
    let remainder = DensePolynomial::from_coefficients_vec(
        remainder_mpc.coeffs.iter().map(|c| c.reveal()).collect()
    );
    let divisor = DensePolynomial::from_coefficients_vec(divisor_coeffs.clone());
    
    // Verify: dividend = quotient * divisor + remainder
    let product = &quotient * &divisor;
    let reconstructed = &product + &remainder;
    let original = DensePolynomial::from_coefficients_vec(dividend_coeffs.clone());
    
    println!("Verifying: dividend = quotient * divisor + remainder");
    
    // Compare coefficients, handling potential trailing zeros
    let max_len = reconstructed.coeffs.len().max(original.coeffs.len());
    for i in 0..max_len {
        let reconstructed_coeff = reconstructed.coeffs.get(i).copied().unwrap_or_else(Fr::zero);
        let original_coeff = original.coeffs.get(i).copied().unwrap_or_else(Fr::zero);
        assert_eq!(
            reconstructed_coeff, original_coeff,
            "Coefficient mismatch at position {}", i
        );
    }
    
    println!("✓✓✓ Test PASSED: divide_with_q_and_r with random PUBLIC coefficients works correctly! ✓✓✓");
}

/// Test with random polynomials and Shared coefficients
#[test]
fn test_divide_with_q_and_r_random_shared() {
    println!("\n=== Test: divide_with_q_and_r with random Shared coefficients ===");
    
    let mut rng = test_rng();
    
    // Create random dividend polynomial of degree 10
    let dividend_degree = 10;
    let dividend_coeffs: Vec<Fr> = (0..=dividend_degree)
        .map(|_| Fr::rand(&mut rng))
        .collect();
    
    // Create random divisor polynomial of degree 3
    let divisor_degree = 3;
    let divisor_coeffs: Vec<Fr> = (0..=divisor_degree)
        .map(|_| Fr::rand(&mut rng))
        .collect();
    
    println!("Dividend polynomial degree: {}", dividend_degree);
    println!("Divisor polynomial degree: {}", divisor_degree);
    
    // Convert dividend to MPC field with SHARED sharing
    let dividend_coeffs_mpc: Vec<MpcField<Fr>> = dividend_coeffs
        .iter()
        .map(|&c| MpcField::from_add_shared(c))
        .collect();
    let dividend_mpc = DensePolynomial::from_coefficients_vec(dividend_coeffs_mpc);
    
    // Divisor is also MPC field but with public values
    let divisor_coeffs_mpc: Vec<MpcField<Fr>> = divisor_coeffs
        .iter()
        .map(|&c| MpcField::from_public(c))
        .collect();
    let divisor_mpc = DensePolynomial::from_coefficients_vec(divisor_coeffs_mpc);
    
    // Convert to DenseOrSparsePolynomial for the function call
    let dividend_sparse: DenseOrSparsePolynomial<MpcField<Fr>> = (&dividend_mpc).into();
    let divisor_sparse: DenseOrSparsePolynomial<MpcField<Fr>> = (&divisor_mpc).into();
    
    // Perform division
    println!("Performing division with random SHARED coefficients...");
    let result = dividend_sparse.divide_with_q_and_r(&divisor_sparse);
    
    assert!(result.is_some(), "Division should succeed");
    let (quotient_mpc, remainder_mpc) = result.unwrap();
    
    println!("✓ Division succeeded");
    println!("  Quotient degree: {}", quotient_mpc.degree());
    println!("  Remainder degree: {}", remainder_mpc.degree());
    
    // Verify remainder degree is less than divisor degree
    assert!(
        remainder_mpc.degree() < divisor_mpc.degree(),
        "Remainder degree should be less than divisor degree"
    );
    
    // Reveal the MPC results
    let quotient = DensePolynomial::from_coefficients_vec(
        quotient_mpc.coeffs.iter().map(|c| c.reveal()).collect()
    );
    let remainder = DensePolynomial::from_coefficients_vec(
        remainder_mpc.coeffs.iter().map(|c| c.reveal()).collect()
    );
    let divisor = DensePolynomial::from_coefficients_vec(divisor_coeffs.clone());
    
    // Verify: dividend = quotient * divisor + remainder
    let product = &quotient * &divisor;
    let reconstructed = &product + &remainder;
    let original = DensePolynomial::from_coefficients_vec(dividend_coeffs.clone());
    
    println!("Verifying: dividend = quotient * divisor + remainder");
    
    // Compare coefficients, handling potential trailing zeros
    let max_len = reconstructed.coeffs.len().max(original.coeffs.len());
    for i in 0..max_len {
        let reconstructed_coeff = reconstructed.coeffs.get(i).copied().unwrap_or_else(Fr::zero);
        let original_coeff = original.coeffs.get(i).copied().unwrap_or_else(Fr::zero);
        assert_eq!(
            reconstructed_coeff, original_coeff,
            "Coefficient mismatch at position {}", i
        );
    }
    
    println!("✓✓✓ Test PASSED: divide_with_q_and_r with random SHARED coefficients works correctly! ✓✓✓");
}
