//! Test to check if plain division algorithm (TODO line 154) is correct

use ark_bls12_377::Fr;
use ark_ff::{Field, One, UniformRand, Zero};
use ark_poly::univariate::{DenseOrSparsePolynomial, DensePolynomial};
use ark_poly::{Polynomial, UVPolynomial};
use ark_std::test_rng;
use mpc_algebra::honest_but_curious::MpcField;
use mpc_algebra::Reveal;

#[test]
fn test_plain_division_algorithm_with_public() {
    println!("\n=== Test: Plain division algorithm (TODO line 154) ===");
    
    let mut rng = test_rng();
    
    // Create test polynomial: p(x) = x^3 + 2x^2 + 3x + 4
    let poly_coeffs: Vec<Fr> = vec![
        Fr::from(4u64),
        Fr::from(3u64),
        Fr::from(2u64),
        Fr::from(1u64),
    ];
    
    // Divisor: d(x) = x + 1
    let divisor_coeffs: Vec<Fr> = vec![
        Fr::from(1u64),
        Fr::from(1u64),
    ];
    
    // Expected quotient: q(x) = x^2 + x + 2
    // Expected remainder: r = 2
    // Verify: (x^2 + x + 2)(x + 1) + 2 = x^3 + x^2 + 2x + x^2 + x + 2 + 2 = x^3 + 2x^2 + 3x + 4 ✓
    
    println!("Testing with MpcField::from_public (uses plain division, line 154)");
    
    // Create polynomial with from_public
    let poly_mpc: Vec<MpcField<Fr>> = poly_coeffs.iter()
        .map(|&c| MpcField::from_public(c))
        .collect();
    let poly = DensePolynomial::from_coefficients_vec(poly_mpc);
    
    let divisor_mpc: Vec<MpcField<Fr>> = divisor_coeffs.iter()
        .map(|&c| MpcField::from_public(c))
        .collect();
    let divisor = DensePolynomial::from_coefficients_vec(divisor_mpc);
    
    println!("Polynomial is_shared: {}", poly.is_shared());
    println!("Divisor is_shared: {}", divisor.is_shared());
    
    // Perform division
    let poly_sparse: DenseOrSparsePolynomial<MpcField<Fr>> = (&poly).into();
    let divisor_sparse: DenseOrSparsePolynomial<MpcField<Fr>> = (&divisor).into();
    
    let result = poly_sparse.divide_with_q_and_r(&divisor_sparse);
    assert!(result.is_some(), "Division should succeed");
    
    let (quotient, remainder) = result.unwrap();
    
    println!("Quotient coefficients: {:?}", quotient.coeffs.iter().map(|c| c.reveal()).collect::<Vec<_>>());
    println!("Remainder coefficients: {:?}", remainder.coeffs.iter().map(|c| c.reveal()).collect::<Vec<_>>());
    
    // Check quotient
    let q_expected = vec![Fr::from(2u64), Fr::from(1u64), Fr::from(1u64)];
    for (i, (&expected, actual)) in q_expected.iter().zip(quotient.coeffs.iter()).enumerate() {
        let actual_revealed = actual.reveal();
        println!("Quotient[{}]: expected={:?}, actual={:?}", i, expected, actual_revealed);
        assert_eq!(expected, actual_revealed, "Quotient coefficient {} mismatch", i);
    }
    
    // Check remainder
    let r_expected = Fr::from(2u64);
    let r_actual = if remainder.coeffs.is_empty() {
        Fr::zero()
    } else {
        remainder.coeffs[0].reveal()
    };
    println!("Remainder: expected={:?}, actual={:?}", r_expected, r_actual);
    assert_eq!(r_expected, r_actual, "Remainder mismatch");
    
    // Verify: p(x) = q(x) * d(x) + r(x)
    let q_plain: Vec<Fr> = quotient.coeffs.iter().map(|c| c.reveal()).collect();
    let d_plain: Vec<Fr> = divisor_coeffs.clone();
    let r_plain: Vec<Fr> = if remainder.coeffs.is_empty() {
        vec![]
    } else {
        remainder.coeffs.iter().map(|c| c.reveal()).collect()
    };
    
    let q_poly = DensePolynomial::from_coefficients_vec(q_plain);
    let d_poly = DensePolynomial::from_coefficients_vec(d_plain);
    let r_poly = DensePolynomial::from_coefficients_vec(r_plain);
    
    let reconstructed = &(&q_poly * &d_poly) + &r_poly;
    let original = DensePolynomial::from_coefficients_vec(poly_coeffs);
    
    println!("\nVerification: q*d + r = p");
    for i in 0..original.coeffs.len() {
        let orig = original.coeffs.get(i).copied().unwrap_or(Fr::zero());
        let recon = reconstructed.coeffs.get(i).copied().unwrap_or(Fr::zero());
        println!("  Coeff[{}]: original={:?}, reconstructed={:?}", i, orig, recon);
        assert_eq!(orig, recon, "Verification failed at coefficient {}", i);
    }
    
    println!("\n✓✓✓ Plain division algorithm is CORRECT!");
}

#[test]
fn test_plain_division_algorithm_arithmetic() {
    println!("\n=== Test: Plain division with pure arithmetic (no MPC) ===");
    
    // Test the same division purely in Fr to verify expected results
    let poly_coeffs: Vec<Fr> = vec![
        Fr::from(4u64),
        Fr::from(3u64),
        Fr::from(2u64),
        Fr::from(1u64),
    ];
    
    let divisor_coeffs: Vec<Fr> = vec![
        Fr::from(1u64),
        Fr::from(1u64),
    ];
    
    let poly = DensePolynomial::from_coefficients_vec(poly_coeffs.clone());
    let divisor = DensePolynomial::from_coefficients_vec(divisor_coeffs.clone());
    
    let poly_sparse: DenseOrSparsePolynomial<Fr> = (&poly).into();
    let divisor_sparse: DenseOrSparsePolynomial<Fr> = (&divisor).into();
    
    let result = poly_sparse.divide_with_q_and_r(&divisor_sparse);
    assert!(result.is_some());
    
    let (quotient, remainder) = result.unwrap();
    
    println!("Quotient: {:?}", quotient.coeffs);
    println!("Remainder: {:?}", remainder.coeffs);
    
    // Verify
    let reconstructed = &(&quotient * &divisor) + &remainder;
    assert_eq!(reconstructed.coeffs, poly.coeffs);
    
    println!("✓ Pure arithmetic division works correctly");
}
