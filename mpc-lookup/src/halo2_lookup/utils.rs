//! # Utility Functions for Halo2 Lookup Argument
//!
//! Helper functions for polynomial operations and selector computation.
//!
//! ## Reusable Components from Existing Codebase
//!
//! **IMPORTANT**: Many utilities already exist in arkworks (algebra crate):
//!
//! ### Lagrange Basis (ALREADY EXISTS):
//! - **Location**: `algebra/poly/src/domain/radix2/mod.rs:119-180`
//! - **Function**: `Radix2EvaluationDomain::evaluate_all_lagrange_coefficients(x)`
//! - **Returns**: Vec<F> where `result[0]` is ℓ_0(x)
//! However, for API consistency with this module, we provide a simple wrapper below.
//!
//! ### Other Available Utilities:
//! - **Vanishing polynomial**: `domain.evaluate_vanishing_polynomial(x)` returns Z_H(x) = x^n - 1
//! - **FFT/IFFT**: `domain.fft(evals)` and `domain.ifft(coeffs)`
//! - **Polynomial shift**: `mpc_plonk::util::shift(f, a)` computes f(aX) from f(X)
//! - **Lagrange interpolation**: `mpc_plonk::util::interpolate(points)` for polynomial construction
//!
//! ### What This Module Implements:
//! - **compute_lagrange_basis_l0**: Wrapper for API consistency
//! - **compute_selector_polynomials**: Halo2-specific (q_blind, q_last) - trivial construction

use ark_ff::{FftField, Field, Zero, One};
use ark_std::vec::Vec;

/// Compute Lagrange basis polynomial ℓ_0(x).
///
/// ℓ_0(X) equals 1 at the first element of the domain (ω^0 = 1) and 0 at all other elements.
///
/// # Formula
///
/// ```text
/// ℓ_0(X) = (X^n - 1) / (n * (X - 1))
/// ```
///
/// where n is the domain size.
///
/// # Arguments
///
/// - `x`: Point to evaluate at
/// - `domain_size`: Size of the evaluation domain
/// - `omega`: Generator of the domain (unused but kept for API consistency)
///
/// # Returns
///
/// ℓ_0(x)
///
/// # Note on Existing Implementation
///
/// This function can also be computed using arkworks' built-in:
/// ```rust,ignore
/// use ark_poly::domain::radix2::Radix2EvaluationDomain;
/// use ark_poly::EvaluationDomain;
/// 
/// let domain = Radix2EvaluationDomain::new(domain_size).unwrap();
/// let l0_x = domain.evaluate_all_lagrange_coefficients(x)[0];
/// ```
///
/// However, for this module's API, implement the direct formula here.
///
/// # Implementation Algorithm
///
/// ```rust,ignore
/// // 1. Compute x^domain_size efficiently using pow()
/// let x_to_n = x.pow(&[domain_size as u64]);
/// 
/// // 2. Compute numerator = x^n - 1
/// let numerator = x_to_n - F::one();
/// 
/// // 3. Compute denominator = n * (x - 1)
/// let n = F::from(domain_size as u64);
/// let denominator = n * (x - F::one());
/// 
/// // 4. Return numerator / denominator
/// numerator * denominator.inverse().unwrap()
/// ```
pub fn compute_lagrange_basis_l0<F: FftField>(x: F, domain_size: usize, _omega: F) -> F {
    // Formula: ℓ_0(X) = (X^n - 1) / (n * (X - 1))
    
    // Special case: when x = 1, the formula is 0/0, but using L'Hôpital's rule we get 1
    if x == F::one() {
        return F::one();
    }
    
    // 1. Compute x^domain_size efficiently using pow()
    let x_to_n = x.pow(&[domain_size as u64]);
    
    // 2. Compute numerator = x^n - 1
    let numerator = x_to_n - F::one();
    
    // 3. Compute denominator = n * (x - 1)
    let n = F::from(domain_size as u64);
    let denominator = n * (x - F::one());
    
    // 4. Return numerator / denominator
    numerator * denominator.inverse().expect("denominator should be non-zero (x != 1)")
}

/// Compute selector polynomials q_blind and q_last.
///
/// These are Halo2-specific selectors that control constraint enforcement.
///
/// # q_blind Selector
///
/// - Equals 1 on the last `num_blinding_rows` rows
/// - Equals 0 on all other rows
/// - Used to disable constraints on blinding rows for zero-knowledge
///
/// # q_last Selector
///
/// - Equals 1 on the boundary row (at index num_usable_rows)
/// - Equals 0 on all other rows
/// - Used for idempotence check: q_last * (Z^2 - Z) = 0 ensures Z ∈ {0,1}
///
/// # Arguments
///
/// - `domain_size`: Size of the evaluation domain (must be power of 2)
/// - `num_blinding_rows`: Number of rows to blind for zero-knowledge
///
/// # Returns
///
/// (q_blind_values, q_last_values) as evaluation vectors of length `domain_size`
///
/// # Implementation Algorithm (SIMPLE - ~10 lines)
///
/// ```rust,ignore
/// // num_usable_rows is the count of usable rows AND the index of the boundary row
/// let num_usable_rows = domain_size - num_blinding_rows - 1;
/// 
/// // Create vectors filled with zeros
/// let mut q_blind = vec![F::zero(); domain_size];
/// let mut q_last = vec![F::zero(); domain_size];
/// 
/// // Set q_blind[i] = 1 for i >= num_usable_rows + 1 (last num_blinding_rows positions)
/// for i in (num_usable_rows + 1)..domain_size {
///     q_blind[i] = F::one();
/// }
/// 
/// // Set q_last[num_usable_rows] = 1 (the boundary row)
/// if num_usable_rows < domain_size {
///     q_last[num_usable_rows] = F::one();
/// }
/// 
/// (q_blind, q_last)
/// ```
///
/// # Example
///
/// For domain_size=8, num_blinding_rows=2:
/// - num_usable_rows = 8 - 2 - 1 = 5 (boundary row index)
/// - Usable data rows: 0-4 (5 rows, indices 0 to num_usable_rows-1)
/// - Boundary row: 5 (index num_usable_rows)
/// - q_blind = [0, 0, 0, 0, 0, 0, 1, 1] (1 on rows 6-7)
/// - q_last =  [0, 0, 0, 0, 0, 1, 0, 0] (1 on row 5)
pub fn compute_selector_polynomials<F: Field>(
    domain_size: usize,
    num_blinding_rows: usize,
) -> (Vec<F>, Vec<F>) {
    // According to Halo2 spec and LookupConfig::new():
    // num_usable_rows = domain_size - num_blinding_rows - 1
    // This value represents both:
    // - The COUNT of usable rows (rows 0 to num_usable_rows-1)
    // - The INDEX of the boundary row (row num_usable_rows)
    let num_usable_rows = domain_size.saturating_sub(num_blinding_rows).saturating_sub(1);
    
    // Create vectors filled with zeros
    let mut q_blind = vec![F::zero(); domain_size];
    let mut q_last = vec![F::zero(); domain_size];
    
    // Set q_blind[i] = 1 for i >= num_usable_rows + 1 (last num_blinding_rows positions)
    // Blinding rows start AFTER the boundary row
    for i in (num_usable_rows + 1)..domain_size {
        q_blind[i] = F::one();
    }
    
    // Set q_last[num_usable_rows] = 1 (the boundary row)
    // According to Halo2 spec, this is row u = domain_size - num_blinding_rows - 1
    // This is where we check Z² - Z = 0
    if num_usable_rows < domain_size {
        q_last[num_usable_rows] = F::one();
    }
    
    (q_blind, q_last)
}

/// Compute the Lagrange basis polynomial ℓ_0(X) as a polynomial.
///
/// ℓ_0(X) equals 1 at X=1 (first domain element) and 0 at all other domain elements.
///
/// # Formula
///
/// ℓ_0(X) = (X^n - 1) / (n * (X - 1))
///
/// We compute this as a polynomial by:
/// 1. Evaluating ℓ_0 at all domain points
/// 2. Using IFFT to convert to coefficient form
///
/// # Arguments
///
/// - `domain_size`: Size of the evaluation domain
/// - `omega`: Generator of the domain
///
/// # Returns
///
/// Coefficient representation of ℓ_0(X) as a DensePolynomial
pub fn compute_lagrange_l0_polynomial<F: FftField>(
    domain_size: usize,
    _omega: F,
) -> ark_poly::univariate::DensePolynomial<F> {
    use ark_poly::domain::radix2::Radix2EvaluationDomain;
    use ark_poly::EvaluationDomain;
    use ark_poly::univariate::DensePolynomial;
    
    // Create domain
    let domain = Radix2EvaluationDomain::<F>::new(domain_size)
        .expect("Failed to create domain for ℓ_0 computation");
    
    // Evaluate ℓ_0 at all domain points
    // At first point (1): ℓ_0 = 1
    // At all other points: ℓ_0 = 0
    let mut evaluations = vec![F::zero(); domain_size];
    evaluations[0] = F::one();
    
    // Convert from evaluations to coefficients using IFFT
    let mut coeffs = evaluations;
    domain.ifft_in_place(&mut coeffs);
    
    DensePolynomial { coeffs }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_377::Fr;
    use ark_poly::domain::radix2::Radix2EvaluationDomain;
    use ark_poly::EvaluationDomain;

    #[test]
    fn test_lagrange_basis_at_first_element() {
        let domain_size = 8;
        let domain = Radix2EvaluationDomain::<Fr>::new(domain_size).unwrap();
        let omega = domain.group_gen;

        // At first element (1), ℓ_0 should be 1
        let l0_at_first = compute_lagrange_basis_l0(Fr::one(), domain_size, omega);
        assert_eq!(l0_at_first, Fr::one());
    }

    #[test]
    fn test_lagrange_basis_at_other_elements() {
        let domain_size = 8;
        let domain = Radix2EvaluationDomain::<Fr>::new(domain_size).unwrap();
        let omega = domain.group_gen;

        // At other domain elements, ℓ_0 should be 0
        let l0_at_omega = compute_lagrange_basis_l0(omega, domain_size, omega);
        assert_eq!(l0_at_omega, Fr::zero());
    }

    #[test]
    fn test_selector_polynomials() {
        let domain_size = 8;
        let num_blinding_rows = 2;
        let (q_blind, q_last) = compute_selector_polynomials::<Fr>(domain_size, num_blinding_rows);

        // Check q_blind: should be 1 on last 2 rows
        assert_eq!(q_blind[0], Fr::zero());
        assert_eq!(q_blind[5], Fr::zero());
        assert_eq!(q_blind[6], Fr::one());
        assert_eq!(q_blind[7], Fr::one());

        // Check q_last: should be 1 only on row 5
        assert_eq!(q_last[5], Fr::one());
        for i in 0..domain_size {
            if i != 5 {
                assert_eq!(q_last[i], Fr::zero());
            }
        }
    }
}
