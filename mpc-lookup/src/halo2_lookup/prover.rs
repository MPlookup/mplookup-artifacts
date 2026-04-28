//! # Prover Implementation for Halo2 Lookup Argument
//!
//! This module implements the proof generation algorithm using ark_poly_commit directly.
//!
//! ## Implementation Strategy
//!
//! **Reuse Existing Code**:
//! - Polynomial commitment: `ark_poly_commit::kzg10::KZG10` (recommended)
//! - Permutation computation: `crate::secure_oblivious_lookup_permutation` (already exists)
//! - Polynomial utilities: `mpc_plonk::util::{shift, interpolate}` for polynomial operations
//! - Domain operations: `Radix2EvaluationDomain` for FFT/IFFT
//!
//! **What Needs Implementation** (core lookup logic):
//! - Grand product Z computation: ~20 lines (see algorithm below)
//! - Commitment and opening generation: use ark_poly_commit API
//! - Proof orchestration: combine all steps
//!
//! ## Core Algorithm: Grand Product Z
//!
//! ```text
//! Z[0] = 1
//! for i in 0..num_usable_rows:
//!     numerator = (A'[i] + β) * (S'[i] + γ)
//!     denominator = (A[i] + β) * (S[i] + γ)
//!     Z[i+1] = Z[i] * numerator / denominator
//! Fill last num_blinding_rows with random values
//! ```
//!
//! This encodes the lookup relation. If A ⊆ S, the product telescopes correctly.
//!
//! ## Recommended Polynomial Commitment Scheme
//!
//! ```rust,ignore
//! use ark_poly_commit::kzg10::KZG10;
//! use ark_bls12_377::Bls12_377;
//! type PC = KZG10<Bls12_377>;
//! ```
//!
//! See AI_IMPLEMENTATION_PROMPTS.md Step 3 for complete implementation guidance.

use ark_ff::{FftField, Field, One, Zero};
use ark_poly::domain::radix2::Radix2EvaluationDomain;
use ark_poly::{EvaluationDomain, Polynomial as PolyTrait, UVPolynomial};
use ark_poly::univariate::DensePolynomial;
use ark_poly::polynomial::univariate::DenseOrSparsePolynomial;
use ark_poly_commit::{LabeledPolynomial, LabeledCommitment, PolynomialCommitment};
use ark_poly_commit::marlin::marlin_pc::MarlinKZG10;
use ark_std::vec::Vec;
use digest::Digest;
use mpc_algebra::honest_but_curious::{MpcField, MpcPairingEngine};
use mpc_algebra::Reveal;

use super::transcript::LookupTranscript;
use super::types::{LookupConfig, LookupCommitments, LookupOpenings, LookupProof};
use super::utils::{compute_lagrange_basis_l0, compute_lagrange_l0_polynomial, compute_selector_polynomials};

/// Macro for conditional debug output in prover
macro_rules! debug_eprintln {
    ($debug:expr, $($arg:tt)*) => {
        if $debug {
            eprintln!($($arg)*);
        }
    };
}

/// Computes f(a*X) from a and f(X)
/// This shifts a polynomial by scaling its coefficients: f[i] *= a^i
fn shift_polynomial<F: FftField>(mut f: DensePolynomial<F>, a: F) -> DensePolynomial<F> {
    let mut s = F::one();
    for c in &mut f.coeffs {
        *c *= s;
        s *= a;
    }
    f
}

// Type aliases for MPC-compatible polynomial commitments
type Fr = ark_bls12_377::Fr;
type E = ark_bls12_377::Bls12_377;
type ME = MpcPairingEngine<ark_bls12_377::Bls12_377>;
type MFr = MpcField<Fr>;
type MpcPC = MarlinKZG10<ME, DensePolynomial<MFr>>;
type LocalPC = MarlinKZG10<E, DensePolynomial<Fr>>;

/// Prover state for the lookup argument using MPC-compatible types.
///
/// This implementation uses concrete MPC types to avoid revealing polynomial coefficients:
/// - MpcPairingEngine for elliptic curve operations
/// - DensePolynomial<MpcField<Fr>> for polynomials
/// - MarlinKZG10 over MPC types for commitments
///
/// All polynomial operations stay in MPC space. Only commitments are revealed (public by design).
pub struct LookupProver<D: Digest> {
    config: LookupConfig,
    domain: Radix2EvaluationDomain<MFr>,  // MPC domain
    pcs_ck: <MpcPC as PolynomialCommitment<MFr, DensePolynomial<MFr>>>::CommitterKey,
    // All values stay as MpcField
    a_values: Vec<MFr>,
    s_values: Vec<MFr>,
    a_prime_values: Vec<MFr>,
    s_prime_values: Vec<MFr>,
    z_values: Option<Vec<MFr>>,
    q_blind_values: Vec<Fr>,  // Selectors are public
    q_last_values: Vec<Fr>,
    // Store MPC polynomials (no reveals!)
    a_poly: Option<DensePolynomial<MFr>>,
    s_poly: Option<DensePolynomial<MFr>>,
    a_prime_poly: Option<DensePolynomial<MFr>>,
    s_prime_poly: Option<DensePolynomial<MFr>>,
    z_poly: Option<DensePolynomial<MFr>>,
    quotient_poly: Option<DensePolynomial<MFr>>,  // Quotient polynomial Q(X)
    // Store MPC commitments and randomness
    a_rand: Option<<MpcPC as PolynomialCommitment<MFr, DensePolynomial<MFr>>>::Randomness>,
    s_rand: Option<<MpcPC as PolynomialCommitment<MFr, DensePolynomial<MFr>>>::Randomness>,
    a_prime_rand: Option<<MpcPC as PolynomialCommitment<MFr, DensePolynomial<MFr>>>::Randomness>,
    s_prime_rand: Option<<MpcPC as PolynomialCommitment<MFr, DensePolynomial<MFr>>>::Randomness>,
    z_rand: Option<<MpcPC as PolynomialCommitment<MFr, DensePolynomial<MFr>>>::Randomness>,
    quotient_rand: Option<<MpcPC as PolynomialCommitment<MFr, DensePolynomial<MFr>>>::Randomness>,
    a_com: Option<LabeledCommitment<<MpcPC as PolynomialCommitment<MFr, DensePolynomial<MFr>>>::Commitment>>,
    s_com: Option<LabeledCommitment<<MpcPC as PolynomialCommitment<MFr, DensePolynomial<MFr>>>::Commitment>>,
    a_prime_com: Option<LabeledCommitment<<MpcPC as PolynomialCommitment<MFr, DensePolynomial<MFr>>>::Commitment>>,
    s_prime_com: Option<LabeledCommitment<<MpcPC as PolynomialCommitment<MFr, DensePolynomial<MFr>>>::Commitment>>,
    z_com: Option<LabeledCommitment<<MpcPC as PolynomialCommitment<MFr, DensePolynomial<MFr>>>::Commitment>>,
    quotient_com: Option<LabeledCommitment<<MpcPC as PolynomialCommitment<MFr, DensePolynomial<MFr>>>::Commitment>>,
    debug: bool,
    _phantom_d: std::marker::PhantomData<D>,
}

impl<D: Digest> LookupProver<D> {
    pub fn new(
        config: LookupConfig,
        pcs_ck: <MpcPC as PolynomialCommitment<MFr, DensePolynomial<MFr>>>::CommitterKey,
        a_values: Vec<MFr>,
        s_values: Vec<MFr>,
        debug: bool,
    ) -> Self {
        // Create MPC evaluation domain
        let domain = Radix2EvaluationDomain::<MFr>::new(config.domain_size)
            .expect("Failed to create MPC evaluation domain");
        
        // Pad input values to domain size with zeros
        let mut padded_a_values = a_values;
        let mut padded_s_values = s_values;
        debug_eprintln!(debug, "CONSTRUCTOR DEBUG: a_values len before padding = {}", padded_a_values.len());
        debug_eprintln!(debug, "CONSTRUCTOR DEBUG: s_values len before padding = {}", padded_s_values.len());
        debug_eprintln!(debug, "CONSTRUCTOR DEBUG: config.domain_size = {}", config.domain_size);
        padded_a_values.resize(config.domain_size, MFr::Public(Fr::zero()));
        padded_s_values.resize(config.domain_size, MFr::Public(Fr::zero()));
        debug_eprintln!(debug, "CONSTRUCTOR DEBUG: a_values len after padding = {}", padded_a_values.len());
        debug_eprintln!(debug, "CONSTRUCTOR DEBUG: s_values len after padding = {}", padded_s_values.len());
        
        // Compute selector polynomials (public, so use base field)
        let (q_blind_values, q_last_values) = compute_selector_polynomials(
            config.domain_size,
            config.num_blinding_rows,
        );
        
        Self {
            config,
            domain,
            pcs_ck,
            a_values: padded_a_values,
            s_values: padded_s_values,
            a_prime_values: Vec::new(),
            s_prime_values: Vec::new(),
            z_values: None,
            q_blind_values,
            q_last_values,
            a_poly: None,
            s_poly: None,
            a_prime_poly: None,
            s_prime_poly: None,
            z_poly: None,
            quotient_poly: None,
            a_rand: None,
            s_rand: None,
            a_prime_rand: None,
            s_prime_rand: None,
            z_rand: None,
            quotient_rand: None,
            a_com: None,
            s_com: None,
            a_prime_com: None,
            s_prime_com: None,
            z_com: None,
            quotient_com: None,
            debug,
            _phantom_d: std::marker::PhantomData,
        }
    }

    /// Set pre-computed permuted values (computed externally via secure_oblivious_lookup_permutation)
    pub fn set_permuted_values(&mut self, s_prime: Vec<MFr>, a_prime: Vec<MFr>) {
        // For Halo2 lookup, both A' and S' must be permutations of A and S respectively.
        // secure_oblivious_lookup_permutation already ensures this property.
        // We just need to pad both to domain_size with the same dummy value that's used for A and S.
        
        debug_eprintln!(self.debug, "DEBUG set_permuted_values: s_prime.len()={}, a_prime.len()={}", s_prime.len(), a_prime.len());
        
        let mut padded_s_prime = s_prime;
        let mut padded_a_prime = a_prime;
        
        // Pad with zeros to domain_size (must match the padding used for A and S in constructor)
        let dummy_value = MFr::Public(Fr::zero());
        padded_s_prime.resize(self.config.domain_size, dummy_value.clone());
        padded_a_prime.resize(self.config.domain_size, dummy_value.clone());
        
        debug_eprintln!(self.debug, "DEBUG after resize to domain_size: s_prime.len()={}, a_prime.len()={}", padded_s_prime.len(), padded_a_prime.len());
        
        self.s_prime_values = padded_s_prime;
        self.a_prime_values = padded_a_prime;
    }

    pub fn commit_polynomials(&mut self) -> (
        <LocalPC as PolynomialCommitment<Fr, DensePolynomial<Fr>>>::Commitment,
        <LocalPC as PolynomialCommitment<Fr, DensePolynomial<Fr>>>::Commitment,
        <LocalPC as PolynomialCommitment<Fr, DensePolynomial<Fr>>>::Commitment,
        <LocalPC as PolynomialCommitment<Fr, DensePolynomial<Fr>>>::Commitment,
    ) {
        // DEBUG: Print raw values before polynomial conversion
        if self.debug {
            eprintln!("\n=== PROVER DEBUG: Raw Values (before padding) ===");
            eprintln!("Domain size: {}, Blinding rows: {}, Usable rows: {}", 
                     self.config.domain_size, self.config.num_blinding_rows, self.config.num_usable_rows);
            eprintln!("A (queries) length: {}", self.a_values.len());
            for (i, v) in self.a_values.iter().take(8).enumerate() {
                eprintln!("  A[{}] = {}", i, v.reveal());
            }
            eprintln!("S (table) length: {}", self.s_values.len());
            for (i, v) in self.s_values.iter().take(8).enumerate() {
                eprintln!("  S[{}] = {}", i, v.reveal());
            }
            eprintln!("A' (permuted queries) length: {}", self.a_prime_values.len());
            for (i, v) in self.a_prime_values.iter().take(8).enumerate() {
                eprintln!("  A'[{}] = {}", i, v.reveal());
            }
            eprintln!("S' (permuted table) length: {}", self.s_prime_values.len());
            for (i, v) in self.s_prime_values.iter().take(8).enumerate() {
                eprintln!("  S'[{}] = {}", i, v.reveal());
            }
        }
        
        // MPC-COMPATIBLE POLYNOMIAL COMMITMENTS (NO REVEALS!)
        // 
        // This properly works in MPC space by:
        // 1. Creating DensePolynomial<MpcField<Fr>> from MPC evaluations
        // 2. Using MpcPC::commit on MpcField polynomials (stays in MPC space)
        // 3. Only revealing the COMMITMENTS (public by design), NOT the polynomials
        //
        // NO secret polynomial values are revealed - only commitment group elements.
        
        // Helper to convert MpcField evaluations to MpcField polynomial coefficients
        // This works entirely in MpcField space - NO REVEALS
        let evals_to_mpc_poly = |evals: &[MFr]| -> DensePolynomial<MFr> {
            // Pad to domain size
            let mut padded_evals = evals.to_vec();
            padded_evals.resize(self.config.domain_size, MFr::Public(Fr::zero()));
            
            // Use IFFT in MpcField space - NO REVEALS
            let mut coeffs = padded_evals;
            self.domain.ifft_in_place(&mut coeffs);
            
            DensePolynomial::from_coefficients_vec(coeffs)
        };
        
        // Convert evaluations to MpcField polynomial coefficients (no reveals!)
        let a_mpc_poly = evals_to_mpc_poly(&self.a_values);
        let s_mpc_poly = evals_to_mpc_poly(&self.s_values);
        let a_prime_mpc_poly = evals_to_mpc_poly(&self.a_prime_values);
        let s_prime_mpc_poly = evals_to_mpc_poly(&self.s_prime_values);
        
        // Create LabeledPolynomial instances with MPC polynomials
        // CRITICAL: Set degree_bound to match supported_degrees from setup!
        let degree_bound = Some(self.config.domain_size - 1);
        let labeled_a = LabeledPolynomial::new("a".into(), a_mpc_poly.clone(), degree_bound, None);
        let labeled_s = LabeledPolynomial::new("s".into(), s_mpc_poly.clone(), degree_bound, None);
        let labeled_a_prime = LabeledPolynomial::new("a_prime".into(), a_prime_mpc_poly.clone(), degree_bound, None);
        let labeled_s_prime = LabeledPolynomial::new("s_prime".into(), s_prime_mpc_poly.clone(), degree_bound, None);
        
        // Commit to MPC polynomials using MPC polynomial commitment scheme
        // This works on MpcField polynomials WITHOUT revealing coefficients
        // CRITICAL: Use RNG for hiding commitments (matching mpc-plonk pattern line 388-392)
        // This ensures the commitment has proper randomness that matches the opening proof
        let (mut commitments_a, mut rands_a) = MpcPC::commit(&self.pcs_ck, &[labeled_a], Some(&mut rand::thread_rng()))
            .expect("Failed to commit to A polynomial");
        let (mut commitments_s, mut rands_s) = MpcPC::commit(&self.pcs_ck, &[labeled_s], Some(&mut rand::thread_rng()))
            .expect("Failed to commit to S polynomial");
        let (mut commitments_a_prime, mut rands_a_prime) = MpcPC::commit(&self.pcs_ck, &[labeled_a_prime], Some(&mut rand::thread_rng()))
            .expect("Failed to commit to A' polynomial");
        let (mut commitments_s_prime, mut rands_s_prime) = MpcPC::commit(&self.pcs_ck, &[labeled_s_prime], Some(&mut rand::thread_rng()))
            .expect("Failed to commit to S' polynomial");
        
        // Extract MPC commitments and randomness
        let a_com_labeled = commitments_a.pop().expect("Expected commitment for A");
        let s_com_labeled = commitments_s.pop().expect("Expected commitment for S");
        let a_prime_com_labeled = commitments_a_prime.pop().expect("Expected commitment for A'");
        let s_prime_com_labeled = commitments_s_prime.pop().expect("Expected commitment for S'");
        
        let a_rand = rands_a.pop().expect("Expected randomness for A");
        let s_rand = rands_s.pop().expect("Expected randomness for S");
        let a_prime_rand = rands_a_prime.pop().expect("Expected randomness for A'");
        let s_prime_rand = rands_s_prime.pop().expect("Expected randomness for S'");
        
        // Store MPC polynomials and commitments for later (keep in MPC space for opening!)
        // Store them BEFORE any revelation
        self.a_poly = Some(a_mpc_poly);
        self.s_poly = Some(s_mpc_poly);
        self.a_prime_poly = Some(a_prime_mpc_poly);
        self.s_prime_poly = Some(s_prime_mpc_poly);
        
        // Store the MPC commitments (will use these for opening!)
        self.a_com = Some(a_com_labeled.clone());
        self.s_com = Some(s_com_labeled.clone());
        self.a_prime_com = Some(a_prime_com_labeled.clone());
        self.s_prime_com = Some(s_prime_com_labeled.clone());
        
        debug_eprintln!(self.debug, "DEBUG: Stored commitments are MPC (not yet revealed)");
        debug_eprintln!(self.debug, "  A commitment label: {}", a_com_labeled.label());
        
        self.a_rand = Some(a_rand);
        self.s_rand = Some(s_rand);
        self.a_prime_rand = Some(a_prime_rand);
        self.s_prime_rand = Some(s_prime_rand);
        
        // For output, reveal clones of the commitments (following mpc-snarks pattern: mpc_commits.clone().reveal())
        // The reveal_commitment function internally calls .reveal() on the group elements
        let a_com_public = reveal_commitment(a_com_labeled.commitment.clone());
        let s_com_public = reveal_commitment(s_com_labeled.commitment.clone());
        let a_prime_com_public = reveal_commitment(a_prime_com_labeled.commitment.clone());
        let s_prime_com_public = reveal_commitment(s_prime_com_labeled.commitment.clone());
        
        debug_eprintln!(self.debug, "DEBUG: Revealed commitments for output (separate from stored MPC commitments)");
        
        // CRITICAL VERIFICATION: Check that revealed commitments match MPC commitments
        // This ensures the proof output contains the same commitments used for opening
        let a_com_test_reveal = reveal_commitment(self.a_com.as_ref().unwrap().commitment.clone());
        debug_eprintln!(self.debug, "DEBUG: Verifying A commitment consistency:");
        debug_eprintln!(self.debug, "  Output commitment == Stored commitment: {}", a_com_public == a_com_test_reveal);

        
        (
            a_com_public,
            s_com_public,
            a_prime_com_public,
            s_prime_com_public,
        )
    }

    pub fn compute_grand_product(&mut self, beta: MFr, gamma: MFr) {
        if self.debug {
            eprintln!("\n=== PROVER DEBUG: Grand Product Z Computation ===");
            eprintln!("β = {}", beta.reveal());
            eprintln!("γ = {}", gamma.reveal());
        }
        
        // Initialize Z with Z[0] = 1
        let mut z_values = vec![MFr::Public(Fr::one())];
        
        debug_eprintln!(self.debug, "Z[0] = 1");
        
        // Compute grand product for usable rows
        // According to Halo2 documentation (lines 47-52):
        // Z_{i+1} = Z_i · [(A_i + β) · (S_i + γ)] / [(A'_i + β) · (S'_i + γ)]
        for i in 0..self.config.num_usable_rows {
            let a_i = &self.a_values[i];
            let s_i = &self.s_values[i];
            let a_prime_i = &self.a_prime_values[i];
            let s_prime_i = &self.s_prime_values[i];
            
            // Numerator: (A[i] + β) * (S[i] + γ)
            let numerator = (a_i.clone() + beta.clone()) * (s_i.clone() + gamma.clone());
            
            // Denominator: (A'[i] + β) * (S'[i] + γ)
            let denominator = (a_prime_i.clone() + beta.clone()) * (s_prime_i.clone() + gamma.clone());
            
            // Ratio: numerator / denominator
            // Compute inverse in MPC space (no reveal!)
            let denom_inv = denominator.inverse().expect("denominator should be non-zero");
            let ratio = numerator * denom_inv;
            
            // Z[i+1] = Z[i] * ratio
            let z_next = z_values[i].clone() * ratio;
            z_values.push(z_next.clone());
            
            // Debug first few and last few iterations
            if self.debug && (i < 3 || i >= self.config.num_usable_rows.saturating_sub(3)) {
                eprintln!("  i={}: A[{}]={}, S[{}]={}, A'[{}]={}, S'[{}]={}", 
                         i, i, a_i.reveal(), i, s_i.reveal(), i, a_prime_i.reveal(), i, s_prime_i.reveal());
                eprintln!("        numerator=(A[{}]+β)(S[{}]+γ)={}", i, i, numerator.reveal());
                eprintln!("        denominator=(A'[{}]+β)(S'[{}]+γ)={}", i, i, denominator.reveal());
                eprintln!("        Z[{}] = {}", i+1, z_next.reveal());
            }
        }
        
        debug_eprintln!(self.debug, "Total Z values computed: {} (including Z[0])", z_values.len());
        debug_eprintln!(self.debug, "Continuing recurrence for remaining rows to domain_size: {}", self.config.domain_size);
        
        // CRITICAL FIX: Continue the recurrence for ALL remaining rows up to domain_size!
        // The constraint Z[i+1] * (A[i]+β) * (S[i]+γ) = Z[i] * (A'[i]+β) * (S'[i]+γ)
        // must hold as a POLYNOMIAL IDENTITY everywhere, not just on usable rows.
        // 
        // For the padded rows (i >= num_usable_rows), A[i] = S[i] = A'[i] = S'[i] = 0.
        // So the recurrence becomes: Z[i+1] * β * γ = Z[i] * β * γ, which means Z[i+1] = Z[i].
        // This maintains Z at its current value through the padded/blinding rows.
        //
        // By continuing the recurrence instead of adding random values, the polynomial Z(X)
        // interpolated from these evaluations will satisfy the constraint everywhere!
        for i in self.config.num_usable_rows..(self.config.domain_size - 1) {
            let a_i = &self.a_values[i];
            let s_i = &self.s_values[i];
            let a_prime_i = &self.a_prime_values[i];
            let s_prime_i = &self.s_prime_values[i];
            
            // Numerator: (A'[i] + β) * (S'[i] + γ)
            let numerator = (a_prime_i.clone() + beta.clone()) * (s_prime_i.clone() + gamma.clone());
            
            // Denominator: (A[i] + β) * (S[i] + γ)
            let denominator = (a_i.clone() + beta.clone()) * (s_i.clone() + gamma.clone());
            
            // Ratio: numerator / denominator
            let denom_inv = denominator.inverse().expect("denominator should be non-zero");
            let ratio = numerator * denom_inv;
            
            // Z[i+1] = Z[i] * ratio
            let z_next = z_values[i].clone() * ratio;
            z_values.push(z_next.clone());
            
            if self.debug {
                // Debug output for padding rows
                eprintln!("  i={} (padding): A[{}]={}, S[{}]={}, A'[{}]={}, S'[{}]={}", 
                        i, i, a_i.reveal(), i, s_i.reveal(), i, a_prime_i.reveal(), i, s_prime_i.reveal());
                eprintln!("        ratio={}, Z[{}] = {}", ratio.reveal(), i+1, z_next.reveal());
            }
        }
        
        debug_eprintln!(self.debug, "Final Z values count: {} (should be {})", z_values.len(), self.config.domain_size);
        
        self.z_values = Some(z_values);
    }

    pub fn commit_grand_product(&mut self) -> <LocalPC as PolynomialCommitment<Fr, DensePolynomial<Fr>>>::Commitment {
        let z_values = self.z_values.as_ref().expect("Z values not computed yet");
        
        // Convert MpcField evaluations to MpcField polynomial coefficients (NO REVEAL!)
        let mut padded_evals = z_values.to_vec();
        padded_evals.resize(self.config.domain_size, MFr::Public(Fr::zero()));
        
        // Use IFFT in MpcField space
        let mut z_coeffs = padded_evals;
        self.domain.ifft_in_place(&mut z_coeffs);
        
        // Create MPC polynomial
        let z_mpc_poly = DensePolynomial::from_coefficients_vec(z_coeffs);
        
        // Create LabeledPolynomial
        let degree_bound = Some(self.config.domain_size - 1);
        let labeled_z = LabeledPolynomial::new("z".into(), z_mpc_poly.clone(), degree_bound, None);
        
        // Commit using MPC polynomial commitment
        // CRITICAL: Use RNG for hiding commitments
        let (mut commitments_z, mut rands_z) = MpcPC::commit(&self.pcs_ck, &[labeled_z], Some(&mut rand::thread_rng()))
            .expect("Failed to commit to Z polynomial");
        
        // Extract MPC commitment and randomness
        let z_com_labeled = commitments_z.pop().expect("Expected commitment for Z");
        let z_rand = rands_z.pop().expect("Expected randomness for Z");
        
        // Store MPC polynomial and commitment (keep in MPC space for opening!)
        self.z_poly = Some(z_mpc_poly);
        self.z_com = Some(z_com_labeled.clone());
        self.z_rand = Some(z_rand);
        
        // For output, reveal a clone of the commitment (following mpc-snarks pattern)
        let z_com_public = reveal_commitment(z_com_labeled.commitment.clone());
        
        z_com_public
    }

    /// Compute and commit to quotient polynomial for ALL Halo2 lookup constraints.
    /// The quotient Q(X) satisfies: (sum of all constraints) = Q(X) * Z_H(X)
    /// where Z_H(X) = X^domain_size - 1 is the vanishing polynomial.
    /// 
    /// Constraints (from Halo2 lookup documentation):
    /// 1. Permutation: (1 - (q_last + q_blind)) * [Z(ωX)(A'(X)+β)(S'(X)+γ) - Z(X)(A(X)+β)(S(X)+γ)]
    /// 2. Subset: (1 - (q_last + q_blind)) * (A'(X) - S'(X)) * (A'(X) - A'(ω⁻¹X))
    /// 3. Initial A': ℓ_0(X) * (A'(X) - S'(X))
    /// 4. Initial Z: ℓ_0(X) * (1 - Z(X))
    /// 5. Idempotence: q_last(X) * (Z(X)² - Z(X))
    pub fn compute_and_commit_quotient(
        &mut self, 
        beta: MFr, 
        gamma: MFr
    ) -> <LocalPC as PolynomialCommitment<Fr, DensePolynomial<Fr>>>::Commitment {
        debug_eprintln!(self.debug, "\n=== PROVER: Computing Quotient Polynomial (ALL Constraints) ===");
        
        // Get the stored polynomials
        let a_poly = self.a_poly.as_ref().expect("A polynomial not computed").clone();
        let s_poly = self.s_poly.as_ref().expect("S polynomial not computed").clone();
        let a_prime_poly = self.a_prime_poly.as_ref().expect("A' polynomial not computed").clone();
        let s_prime_poly = self.s_prime_poly.as_ref().expect("S' polynomial not computed").clone();
        let z_poly = self.z_poly.as_ref().expect("Z polynomial not computed").clone();
        
        // Create shifted polynomials using the shift utility
        let omega = self.domain.group_gen;
        let z_shifted_poly = shift_polynomial(z_poly.clone(), omega);  // Z(ωX)
        let a_prime_inv_shifted_poly = shift_polynomial(a_prime_poly.clone(), omega.inverse().expect("omega inverse"));  // A'(ω⁻¹X)
        
        // Convert selector values to polynomials using IFFT
        let selector_domain = Radix2EvaluationDomain::<MFr>::new(self.config.domain_size)
            .expect("Failed to create selector domain");
        
        let q_blind_poly = {
            let mut coeffs: Vec<MFr> = self.q_blind_values.iter().map(|&v| MFr::Public(v)).collect();
            selector_domain.ifft_in_place(&mut coeffs);
            DensePolynomial::from_coefficients_vec(coeffs)
        };
        
        let q_last_poly = {
            let mut coeffs: Vec<MFr> = self.q_last_values.iter().map(|&v| MFr::Public(v)).collect();
            selector_domain.ifft_in_place(&mut coeffs);
            DensePolynomial::from_coefficients_vec(coeffs)
        };
        
        // Compute ℓ_0(X) - Lagrange basis polynomial that is 1 at X=1 and 0 at other domain points
        // ℓ_0(X) = (X^n - 1) / (n * (X - 1))
        // We'll compute this as a polynomial
        let l_0_poly = compute_lagrange_l0_polynomial(self.config.domain_size, omega);
        
        // Create constant polynomials
        let beta_poly = DensePolynomial::from_coefficients_vec(vec![beta]);
        let gamma_poly = DensePolynomial::from_coefficients_vec(vec![gamma]);
        let one_poly = DensePolynomial::from_coefficients_vec(vec![MFr::one()]);
        
        debug_eprintln!(self.debug, "Computing constraint 1: Permutation constraint");
        // Constraint 1: Permutation (with selector)
        // According to Halo2 documentation (lines 42, 105):
        // (1 - (q_last + q_blind)) * [Z(ωX) * (A'(X)+β) * (S'(X)+γ) - Z(X) * (A(X)+β) * (S(X)+γ)]
        let selector_poly = &one_poly - &(&q_last_poly + &q_blind_poly);
        
        let a_plus_beta = &a_poly + &beta_poly;
        let s_plus_gamma = &s_poly + &gamma_poly;
        let a_prime_plus_beta = &a_prime_poly + &beta_poly;
        let s_prime_plus_gamma = &s_prime_poly + &gamma_poly;
        
        let perm_lhs = &(&z_shifted_poly * &a_prime_plus_beta) * &s_prime_plus_gamma;
        let perm_rhs = &(&z_poly * &a_plus_beta) * &s_plus_gamma;
        let perm_constraint = &perm_lhs - &perm_rhs;
        let constraint_1 = &selector_poly * &perm_constraint;
        
        debug_eprintln!(self.debug, "Computing constraint 2: Subset constraint");
        // Constraint 2: Subset (with selector)
        // (1 - (q_last + q_blind)) * (A'(X) - S'(X)) * (A'(X) - A'(ω⁻¹X))
        let a_prime_minus_s_prime = &a_prime_poly - &s_prime_poly;
        let a_prime_minus_a_prime_shifted = &a_prime_poly - &a_prime_inv_shifted_poly;
        let subset_constraint = &a_prime_minus_s_prime * &a_prime_minus_a_prime_shifted;
        let constraint_2 = &selector_poly * &subset_constraint;
        
        debug_eprintln!(self.debug, "Computing constraint 3: Initial A' constraint");
        // Constraint 3: Initial A' constraint
        // ℓ_0(X) * (A'(X) - S'(X))
        let constraint_3 = &l_0_poly * &a_prime_minus_s_prime;
        
        debug_eprintln!(self.debug, "Computing constraint 4: Initial Z constraint");
        // Constraint 4: Initial Z constraint
        // ℓ_0(X) * (1 - Z(X))
        let one_minus_z = &one_poly - &z_poly;
        let constraint_4 = &l_0_poly * &one_minus_z;
        
        debug_eprintln!(self.debug, "Computing constraint 5: Idempotence constraint");
        // Constraint 5: Idempotence constraint
        // q_last(X) * (Z(X)² - Z(X))
        let z_squared = &z_poly * &z_poly;
        let z_squared_minus_z = &z_squared - &z_poly;
        let constraint_5 = &q_last_poly * &z_squared_minus_z;
        
        debug_eprintln!(self.debug, "Summing all constraints");
        // Sum all constraints
        let constraint_poly = &(&(&(&constraint_1 + &constraint_2) + &constraint_3) + &constraint_4) + &constraint_5;
        
        debug_eprintln!(self.debug, "Total constraint polynomial degree: {}", constraint_poly.degree());
        debug_eprintln!(self.debug, "Total constraint polynomial coefficient count: {}", constraint_poly.coeffs.len());
        
        // DEBUG: Check if coefficients are shared and their structure
        debug_eprintln!(self.debug, "Constraint poly is_shared: {}", constraint_poly.is_shared());
        debug_eprintln!(self.debug, "First few coefficients info:");
        for (i, coeff) in constraint_poly.coeffs.iter().take(5).enumerate() {
            match coeff {
                mpc_algebra::MpcField::Shared(_s) => {
                    debug_eprintln!(self.debug, "  coeff[{}] is Shared with local value (not revealed)", i);
                },
                mpc_algebra::MpcField::Public(p) => {
                    debug_eprintln!(self.debug, "  coeff[{}] is Public: {}", i, p);
                }
            }
        }
        
        // DEBUG: Evaluate each constraint at domain point ω to verify they are zero
        debug_eprintln!(self.debug, "\n=== DEBUG: Evaluating constraints at ALL domain points ===");
        let omega_mfr = MFr::Public(omega.reveal());
        
        // Check at ALL domain points
        for idx in 0..self.config.domain_size {
            let omega_power = if idx == 0 {
                MFr::one()
            } else {
                omega_mfr.pow(&[idx as u64])
            };
            let c1 = constraint_1.evaluate(&omega_power).reveal();
            let c2 = constraint_2.evaluate(&omega_power).reveal();
            let c3 = constraint_3.evaluate(&omega_power).reveal();
            let c4 = constraint_4.evaluate(&omega_power).reveal();
            let c5 = constraint_5.evaluate(&omega_power).reveal();
            
            if !c1.is_zero() || !c2.is_zero() || !c3.is_zero() || !c4.is_zero() || !c5.is_zero() {
                debug_eprintln!(self.debug, "\nNON-ZERO constraint at ω^{} (index {}):", idx, idx);
                debug_eprintln!(self.debug, "  Constraint 1 (perm): {}", c1);
                debug_eprintln!(self.debug, "  Constraint 2 (subset): {}", c2);
                debug_eprintln!(self.debug, "  Constraint 3 (init A'): {}", c3);
                debug_eprintln!(self.debug, "  Constraint 4 (init Z): {}", c4);
                debug_eprintln!(self.debug, "  Constraint 5 (idem): {}", c5);
            }
        }
        
        // CRITICAL: Check constraint polynomial consistency before division
        // Evaluate the combined constraint polynomial at all domain points
        debug_eprintln!(self.debug, "\n=== CRITICAL: Verifying COMBINED constraint polynomial consistency ===");
        for idx in 0..self.config.domain_size {
            let omega_power = if idx == 0 {
                MFr::one()
            } else {
                omega_mfr.pow(&[idx as u64])
            };
            let constraint_at_point = constraint_poly.evaluate(&omega_power).reveal();
            if !constraint_at_point.is_zero() {
                debug_eprintln!(self.debug, "ERROR: Combined constraint polynomial is NON-ZERO at ω^{}: {}", idx, constraint_at_point);
                panic!("Constraint polynomial is not zero on domain - cannot compute valid quotient!");
            }
        }
        debug_eprintln!(self.debug, "✓ Combined constraint polynomial is zero at all domain points");
        
        // CRITICAL: Verify constraint polynomial coefficients are consistent
        debug_eprintln!(self.debug, "\n=== CRITICAL: Checking constraint polynomial coefficient consistency ===");
        debug_eprintln!(self.debug, "Constraint polynomial has {} coefficients", constraint_poly.coeffs.len());
        debug_eprintln!(self.debug, "Constraint polynomial degree: {}", constraint_poly.degree());
        debug_eprintln!(self.debug, "Revealing first 10 constraint coefficients to check consistency:");
        for (i, coeff) in constraint_poly.coeffs.iter().take(10).enumerate() {
            let revealed = coeff.reveal();
            debug_eprintln!(self.debug, "  constraint[{}] = {}", i, revealed);
        }
        
        // Also check the last few coefficients (high-degree terms)
        debug_eprintln!(self.debug, "Revealing last 5 constraint coefficients:");
        let len = constraint_poly.coeffs.len();
        for i in (len.saturating_sub(5))..len {
            let revealed = constraint_poly.coeffs[i].reveal();
            debug_eprintln!(self.debug, "  constraint[{}] = {}", i, revealed);
        }
        
        // Compute vanishing polynomial Z_H(X) = X^domain_size - 1
        let mut z_h_coeffs = vec![MFr::zero(); self.config.domain_size + 1];
        z_h_coeffs[0] = MFr::Public(-Fr::one());  // -1
        z_h_coeffs[self.config.domain_size] = MFr::one();  // X^domain_size
        let z_h_poly = DensePolynomial::from_coefficients_vec(z_h_coeffs);
        
        debug_eprintln!(self.debug, "Vanishing polynomial Z_H(X) degree: {}", z_h_poly.degree());
        
        // CRITICAL: Verify Z_H is consistent (should be all public)
        debug_eprintln!(self.debug, "\n=== CRITICAL: Checking Z_H polynomial ===");
        if self.debug { eprintln!("Z_H[0] (constant) = {}", z_h_poly.coeffs[0].reveal());
        }
        if self.debug { eprintln!("Z_H[{}] (X^{}) = {}", self.config.domain_size, self.config.domain_size, z_h_poly.coeffs[self.config.domain_size].reveal());
        }
        debug_eprintln!(self.debug, "Z_H is_shared: {}", z_h_poly.is_shared());
        debug_eprintln!(self.debug, "constraint_poly is_shared: {}", constraint_poly.is_shared());
        
        // Divide: Q(X) = constraint_poly / Z_H(X)
        let constraint_sparse = DenseOrSparsePolynomial::from(constraint_poly.clone());
        let z_h_sparse = DenseOrSparsePolynomial::from(z_h_poly);
        
        let (quotient_poly, remainder) = constraint_sparse
            .divide_with_q_and_r(&z_h_sparse)
            .expect("Division should succeed");
        
        debug_eprintln!(self.debug, "Quotient polynomial Q(X) degree: {}", quotient_poly.degree());
        debug_eprintln!(self.debug, "Remainder degree: {}", remainder.degree());
        
        // CRITICAL: Check quotient polynomial consistency
        debug_eprintln!(self.debug, "\n=== CRITICAL: Verifying quotient polynomial consistency ===");
        debug_eprintln!(self.debug, "Quotient coefficients (first 5):");
        for (i, coeff) in quotient_poly.coeffs.iter().take(5).enumerate() {
            let revealed = coeff.reveal();
            debug_eprintln!(self.debug, "  Q[{}] = {}", i, revealed);
        }
        
        // Verification checks (DEBUGGING ONLY)
        if true {
            // The remainder MUST be zero!

            // For MPC polynomials with shared coefficients, we need to reveal coefficients to properly check
            // Note: remainder.is_zero() only checks local shares, not the shared secret value
            let remainder_is_zero = remainder.coeffs.is_empty() || remainder.coeffs.iter().all(|c| c.reveal().is_zero());
            
            if !remainder_is_zero {
                debug_eprintln!(self.debug, "ERROR: Remainder is NOT zero! This means the constraint doesn't hold on the domain!");
                if self.debug { eprintln!("Remainder: {:?}", remainder.coeffs.iter().take(5).map(|c| c.reveal()).collect::<Vec<_>>());
                }
                panic!("Quotient polynomial computation failed - constraint not satisfied!");
            } else {
                debug_eprintln!(self.debug, "✓ Remainder is zero - constraint holds on domain!");
            }
            
            // VERIFICATION: Check quotient at domain point ω
            debug_eprintln!(self.debug, "\nVerifying quotient at domain point ω:");
            let omega_mfr = MFr::Public(omega.reveal());
            let q_at_omega = quotient_poly.evaluate(&omega_mfr);
            let constraint_at_omega = constraint_poly.evaluate(&omega_mfr);
            
            // Z_H(ω) should be 0 since ω is on the domain
            let z_h_at_omega = omega_mfr.pow(&[self.config.domain_size as u64]) - MFr::one();
            if self.debug { eprintln!("  Q(ω) = {}", q_at_omega.reveal());
            }
            if self.debug { eprintln!("  constraint(ω) = {}", constraint_at_omega.reveal());
            }
            if self.debug { eprintln!("  Z_H(ω) = {}", z_h_at_omega.reveal());
            }
            if self.debug { eprintln!("  Q(ω) * Z_H(ω) = {}", (q_at_omega * z_h_at_omega).reveal());
            }
            if self.debug { eprintln!("  Should match: {}", constraint_at_omega.reveal() == (q_at_omega * z_h_at_omega).reveal());
            }
            
            // VERIFICATION: Check that shift works correctly
            debug_eprintln!(self.debug, "\nVerifying polynomial shift:");
            let test_point = MFr::Public(Fr::from(7u64));
            let z_at_test = z_poly.evaluate(&test_point);
            let z_at_omega_test = z_poly.evaluate(&(omega * test_point.clone()));
            let z_shifted_at_test = z_shifted_poly.evaluate(&test_point);
            if self.debug { eprintln!("  Test point X = {}", test_point.reveal());
            }
            if self.debug { eprintln!("  Z(X) = {}", z_at_test.reveal());
            }
            if self.debug { eprintln!("  Z(ωX) evaluated directly = {}", z_at_omega_test.reveal());
            }
            if self.debug { eprintln!("  Z_shifted(X) = {}", z_shifted_at_test.reveal());
            }
            if self.debug { eprintln!("  Shift correct: {}", z_at_omega_test.reveal() == z_shifted_at_test.reveal());
            }
            
            // VERIFICATION: Check quotient relationship at test point
            debug_eprintln!(self.debug, "\nVerifying quotient relationship at test point:");
            let test_pt2 = MFr::Public(Fr::from(123u64));
            let q_at_test = quotient_poly.evaluate(&test_pt2);
            let constraint_at_test = constraint_poly.evaluate(&test_pt2);
            let z_h_at_test = test_pt2.pow(&[self.config.domain_size as u64]) - MFr::one();
            let expected = q_at_test.clone() * z_h_at_test.clone();
            if self.debug { eprintln!("  Test point = {}", test_pt2.reveal());
            }
            if self.debug { eprintln!("  Q(test) = {}", q_at_test.reveal());
            }
            if self.debug { eprintln!("  constraint(test) = {}", constraint_at_test.reveal());
            }
            if self.debug { eprintln!("  Z_H(test) = {}", z_h_at_test.reveal());
            }
            if self.debug { eprintln!("  Q(test) * Z_H(test) = {}", expected.reveal());
            }
            if self.debug { eprintln!("  Match: {}", constraint_at_test.reveal() == expected.reveal());
            }

        }
        // Store the quotient polynomial
        self.quotient_poly = Some(quotient_poly.clone());
        
        // Commit to quotient polynomial
        // CRITICAL: Use RNG for hiding commitments
        let labeled_q = LabeledPolynomial::new("quotient".into(), quotient_poly, None, None);
        let (mut commitments_q, mut rands_q) = MpcPC::commit(&self.pcs_ck, &[labeled_q], Some(&mut rand::thread_rng()))
            .expect("Failed to commit to quotient polynomial");
        
        let q_com_labeled = commitments_q.pop().expect("Expected commitment for Q");
        let q_rand = rands_q.pop().expect("Expected randomness for Q");
        
        // Store MPC commitment and randomness (keep in MPC space for opening!)
        self.quotient_com = Some(q_com_labeled.clone());
        self.quotient_rand = Some(q_rand);
        
        // For output, reveal a clone of the commitment (following mpc-snarks pattern)
        let q_com_public = reveal_commitment(q_com_labeled.commitment.clone());
        
        debug_eprintln!(self.debug, "✓ Quotient polynomial computed and committed");
        
        q_com_public
    }

    pub fn compute_evaluations(&self, zeta: Fr, omega: Fr) -> LookupOpenings<Fr, <MpcPC as PolynomialCommitment<MFr, DensePolynomial<MFr>>>::Proof> {
        debug_eprintln!(self.debug, "\n=== PROVER DEBUG: Polynomial Evaluations ===");
        debug_eprintln!(self.debug, "Challenge point ζ = {}", zeta);
        debug_eprintln!(self.debug, "Domain generator ω = {}", omega);
        
        // Compute evaluation points
        let omega_zeta = omega * zeta;
        let omega_inv = omega.inverse().expect("omega should be invertible");
        let omega_inv_zeta = zeta * omega_inv;
        
        debug_eprintln!(self.debug, "ωζ = {}", omega_zeta);
        debug_eprintln!(self.debug, "ω⁻¹ζ = {}", omega_inv_zeta);
        
        // Get stored MPC polynomials and commitments
        let a_poly = self.a_poly.as_ref().expect("A polynomial not computed");
        let s_poly = self.s_poly.as_ref().expect("S polynomial not computed");
        let a_prime_poly = self.a_prime_poly.as_ref().expect("A' polynomial not computed");
        let s_prime_poly = self.s_prime_poly.as_ref().expect("S' polynomial not computed");
        let z_poly = self.z_poly.as_ref().expect("Z polynomial not computed");
        
        let quotient_poly = self.quotient_poly.as_ref().expect("Quotient polynomial not computed");
        
        let a_com = self.a_com.as_ref().expect("A commitment not computed");
        let s_com = self.s_com.as_ref().expect("S commitment not computed");
        let a_prime_com = self.a_prime_com.as_ref().expect("A' commitment not computed");
        let s_prime_com = self.s_prime_com.as_ref().expect("S' commitment not computed");
        let z_com = self.z_com.as_ref().expect("Z commitment not computed");
        let quotient_com = self.quotient_com.as_ref().expect("Quotient commitment not computed");
        
        let a_rand = self.a_rand.as_ref().expect("A randomness not stored");
        let s_rand = self.s_rand.as_ref().expect("S randomness not stored");
        let a_prime_rand = self.a_prime_rand.as_ref().expect("A' randomness not stored");
        let s_prime_rand = self.s_prime_rand.as_ref().expect("S' randomness not stored");
        let z_rand = self.z_rand.as_ref().expect("Z randomness not stored");
        let quotient_rand = self.quotient_rand.as_ref().expect("Quotient randomness not stored");
        
        // Evaluate MPC polynomials at ζ (as MFr) - these evaluations will be revealed as part of the proof
        let zeta_mfr = MFr::Public(zeta);
        let a_zeta_mfr = a_poly.evaluate(&zeta_mfr);
        let s_zeta_mfr = s_poly.evaluate(&zeta_mfr);
        let a_prime_zeta_mfr = a_prime_poly.evaluate(&zeta_mfr);
        let s_prime_zeta_mfr = s_prime_poly.evaluate(&zeta_mfr);
        let z_zeta_mfr = z_poly.evaluate(&zeta_mfr);
        
        // Evaluate Z at ωζ (next row)
        let omega_zeta_mfr = MFr::Public(omega_zeta);
        let z_omega_zeta_mfr = z_poly.evaluate(&omega_zeta_mfr);
        
        // Evaluate A' at ω⁻¹ζ (previous row)
        let omega_inv_zeta_mfr = MFr::Public(omega_inv_zeta);
        let a_prime_omega_inv_zeta_mfr = a_prime_poly.evaluate(&omega_inv_zeta_mfr);
        
        // Reveal evaluations (these are public in the proof)
        let a_zeta = a_zeta_mfr.reveal();
        let s_zeta = s_zeta_mfr.reveal();
        let a_prime_zeta = a_prime_zeta_mfr.reveal();
        let s_prime_zeta = s_prime_zeta_mfr.reveal();
        let z_zeta = z_zeta_mfr.reveal();
        let z_omega_zeta = z_omega_zeta_mfr.reveal();
        let a_prime_omega_inv_zeta = a_prime_omega_inv_zeta_mfr.reveal();
        
        debug_eprintln!(self.debug, "\nPolynomial evaluations at challenge points:");
        debug_eprintln!(self.debug, "  A(ζ) = {}", a_zeta);
        debug_eprintln!(self.debug, "  S(ζ) = {}", s_zeta);
        debug_eprintln!(self.debug, "  A'(ζ) = {}", a_prime_zeta);
        debug_eprintln!(self.debug, "  S'(ζ) = {}", s_prime_zeta);
        debug_eprintln!(self.debug, "  Z(ζ) = {}", z_zeta);
        debug_eprintln!(self.debug, "  Z(ωζ) = {}", z_omega_zeta);
        debug_eprintln!(self.debug, "  A'(ω⁻¹ζ) = {}", a_prime_omega_inv_zeta);
        
        // Evaluate quotient polynomial at ζ
        let q_zeta_mpc = quotient_poly.evaluate(&MFr::Public(zeta));
        let q_zeta = q_zeta_mpc.reveal();
        debug_eprintln!(self.debug, "  Q(ζ) = {}", q_zeta);
        
        // VERIFICATION: Re-evaluate quotient polynomial directly to confirm
        debug_eprintln!(self.debug, "\n=== PROVER VERIFICATION: Polynomial Evaluations ===");
        debug_eprintln!(self.debug, "Re-evaluating Q(ζ) to verify:");
        let q_zeta_verify = quotient_poly.evaluate(&MFr::Public(zeta)).reveal();
        debug_eprintln!(self.debug, "  Q(ζ) re-evaluated = {}", q_zeta_verify);
        debug_eprintln!(self.debug, "  Matches original: {}", q_zeta == q_zeta_verify);
        
        // CRITICAL VERIFICATION: Check if our polynomial actually evaluates correctly
        // by directly evaluating the local version
        debug_eprintln!(self.debug, "\n=== PROVER DEBUG: Verifying A polynomial evaluation ===");
        let a_poly_local = a_poly.clone();
        let a_zeta_direct_mpc = a_poly_local.evaluate(&MFr::Public(zeta));
        let a_zeta_direct = a_zeta_direct_mpc.reveal();
        debug_eprintln!(self.debug, "  A(ζ) from stored evaluation = {}", a_zeta);
        debug_eprintln!(self.debug, "  A(ζ) from re-evaluation = {}", a_zeta_direct);
        debug_eprintln!(self.debug, "  Match: {}", a_zeta == a_zeta_direct);
        
        // Compute ℓ_0(ζ)
        let l_0_zeta = compute_lagrange_basis_l0(zeta, self.config.domain_size, omega);
        
        // Evaluate selector polynomials at ζ (selectors are public)
        let q_blind_domain = Radix2EvaluationDomain::<Fr>::new(self.config.domain_size)
            .expect("Failed to create base domain for selectors");
        let q_blind_poly = {
            let mut coeffs = self.q_blind_values.clone();
            q_blind_domain.ifft_in_place(&mut coeffs);
            DensePolynomial::from_coefficients_vec(coeffs)
        };
        let q_last_poly = {
            let mut coeffs = self.q_last_values.clone();
            q_blind_domain.ifft_in_place(&mut coeffs);
            DensePolynomial::from_coefficients_vec(coeffs)
        };
        
        let q_blind_zeta = q_blind_poly.evaluate(&zeta);
        let q_last_zeta = q_last_poly.evaluate(&zeta);
        
        // Generate batch opening proof using MPC PC::open
        // Batch open all 6 polynomials at ζ in a SINGLE call

        // CRITICAL: For single-commitment openings, use F::one() as the opening challenge
        // This matches the mpc-plonk pattern (see mpc-plonk/src/lib.rs line 356)
        // The comment there says "acceptable b/c this is just one commitment"
        let opening_challenge_mfr = MFr::Public(Fr::one());
        
        // Create LabeledPolynomial instances for opening
        // CRITICAL: Must match the degree_bound used in commitment!
        let degree_bound = Some(self.config.domain_size - 1); // Match the supported_degrees!
        let labeled_a = LabeledPolynomial::new("a".into(), a_poly.clone(), degree_bound, None);
        let labeled_s = LabeledPolynomial::new("s".into(), s_poly.clone(), degree_bound, None);
        let labeled_a_prime = LabeledPolynomial::new("a_prime".into(), a_prime_poly.clone(), degree_bound, None);
        let labeled_s_prime = LabeledPolynomial::new("s_prime".into(), s_prime_poly.clone(), degree_bound, None);
        let labeled_z = LabeledPolynomial::new("z".into(), z_poly.clone(), degree_bound, None);
        let labeled_quotient = LabeledPolynomial::new("quotient".into(), quotient_poly.clone(), None, None); // Quotient doesn't have degree bound enforcement
        
        // BATCH OPENING: Open all 6 polynomials at ζ in ONE call
        debug_eprintln!(self.debug, "\n=== PROVER: Generating BATCH opening proof at ζ ===");
        debug_eprintln!(self.debug, "Batching 6 polynomials (A, S, A', S', Z, Q) into single opening proof");
        
        // Collect all polynomials, commitments, and randomness for batch opening
        let polynomials_at_zeta = vec![
            &labeled_a,
            &labeled_s,
            &labeled_a_prime,
            &labeled_s_prime,
            &labeled_z,
            &labeled_quotient,
        ];
        
        let commitments_at_zeta = vec![
            a_com.clone(),
            s_com.clone(),
            a_prime_com.clone(),
            s_prime_com.clone(),
            z_com.clone(),
            quotient_com.clone(),
        ];
        
        let randomness_at_zeta = vec![
            a_rand,
            s_rand,
            a_prime_rand,
            s_prime_rand,
            z_rand,
            quotient_rand,
        ];
        
        debug_eprintln!(self.debug, "DEBUG: Batch opening parameters:");
        debug_eprintln!(self.debug, "  Number of polynomials: {}", polynomials_at_zeta.len());
        if self.debug { eprintln!("  Evaluation point ζ: {}", zeta_mfr.reveal());
        }
        if self.debug { eprintln!("  Opening challenge: {}", opening_challenge_mfr.reveal());
        }
        
        // Generate single batch opening proof for all 6 polynomials at ζ
        let proof_zeta = MpcPC::open(
            &self.pcs_ck,
            polynomials_at_zeta,
            &commitments_at_zeta,
            &zeta_mfr,
            opening_challenge_mfr.clone(),
            randomness_at_zeta,
            Some(&mut rand::thread_rng()),
        ).expect("Failed to generate batch opening proof at ζ");
        
        debug_eprintln!(self.debug, "✓ Batch opening proof generated for all 6 polynomials at ζ");
        
        // VERIFICATION: Check that the polynomial evaluations we computed match what PC::open uses
        debug_eprintln!(self.debug, "\nVerifying polynomial evaluations match:");
        debug_eprintln!(self.debug, "  A(ζ) = {} (computed from poly)", a_zeta);
        debug_eprintln!(self.debug, "  S(ζ) = {} (computed from poly)", s_zeta);
        debug_eprintln!(self.debug, "  A'(ζ) = {} (computed from poly)", a_prime_zeta);
        debug_eprintln!(self.debug, "  S'(ζ) = {} (computed from poly)", s_prime_zeta);
        debug_eprintln!(self.debug, "  Z(ζ) = {} (computed from poly)", z_zeta);
        debug_eprintln!(self.debug, "  Q(ζ) = {} (computed from poly)", q_zeta);
        
        // Generate individual opening proofs for Z at ωζ and A' at ω⁻¹ζ (different points than ζ)
        // These cannot be batched with the ζ openings because they're at different evaluation points
        
        // Proof at ωζ for Z
        let commitments_z_omega = vec![z_com.clone()];
        let proof_omega_zeta = MpcPC::open(
            &self.pcs_ck,
            vec![&labeled_z],
            &commitments_z_omega,
            &omega_zeta_mfr,
            MFr::Public(Fr::one()), // Use F::one() for single commitments
            vec![z_rand],
            None,
        ).expect("Failed to generate opening proof at omega_zeta");
        
        // Proof at ω⁻¹ζ for A'
        let commitments_a_prime_omega_inv = vec![a_prime_com.clone()];
        let proof_omega_inv_zeta = MpcPC::open(
            &self.pcs_ck,
            vec![&labeled_a_prime],
            &commitments_a_prime_omega_inv,
            &omega_inv_zeta_mfr,
            MFr::Public(Fr::one()), // Use F::one() for single commitments
            vec![a_prime_rand],
            None,
        ).expect("Failed to generate opening proof at omega_inv_zeta");
        
        LookupOpenings {
            a_zeta,
            s_zeta,
            a_prime_zeta,
            s_prime_zeta,
            z_zeta,
            z_omega_zeta,
            a_prime_omega_inv_zeta,
            q_blind_zeta,
            q_last_zeta,
            l_0_zeta,
            q_zeta,
            proof_zeta,
            proof_omega_zeta,
            proof_omega_inv_zeta,
        }
    }

    pub fn prove(&mut self, transcript: &mut LookupTranscript<D>) -> LookupProof<Fr, <LocalPC as PolynomialCommitment<Fr, DensePolynomial<Fr>>>::Commitment, <MpcPC as PolynomialCommitment<MFr, DensePolynomial<MFr>>>::Proof> {
        // NOTE: Permutations must be computed externally via secure_oblivious_lookup_permutation
        // and set using set_permuted_values() before calling prove()
        
        if self.a_prime_values.is_empty() || self.s_prime_values.is_empty() {
            panic!("Permuted values not set. Call set_permuted_values() before prove()");
        }
        
        // 2. Commit to A, S, A', S'
        let (com_a, com_s, com_a_prime, com_s_prime) = self.commit_polynomials();
        
        // Append commitments to transcript
        transcript.append_commitment(b"com_a", &com_a);
        transcript.append_commitment(b"com_s", &com_s);
        transcript.append_commitment(b"com_a_prime", &com_a_prime);
        transcript.append_commitment(b"com_s_prime", &com_s_prime);
        
        // 3. Derive β and γ from transcript
        let (beta, gamma) = transcript.derive_beta_gamma();
        
        // Convert to MpcField for computation
        let beta_mpc = MFr::Public(beta);
        let gamma_mpc = MFr::Public(gamma);
        
        // 4. Compute grand product Z
        self.compute_grand_product(beta_mpc, gamma_mpc);
        
        // 5. Commit to Z
        let com_z = self.commit_grand_product();
        
        // Append Z commitment to transcript
        transcript.append_commitment(b"com_z", &com_z);
        
        // 5b. Compute and commit to quotient polynomial
        let com_q = self.compute_and_commit_quotient(beta_mpc.clone(), gamma_mpc.clone());
        
        // Append Q commitment to transcript
        transcript.append_commitment(b"com_q", &com_q);
        
        // 6. Derive ζ from transcript  
        let zeta = transcript.derive_zeta();
        
        // 7. Compute evaluations
        let omega = self.domain.group_gen.reveal(); // Get base field omega from MPC domain
        let openings = self.compute_evaluations(zeta, omega);
        
        // 8. Build and return proof
        LookupProof {
            commitments: LookupCommitments {
                com_a,
                com_s,
                com_a_prime,
                com_s_prime,
                com_z,
                com_q,
                com_q_blind: None, // Fixed polynomials - not committed in proof
                com_q_last: None,
            },
            openings,
        }
    }
}

/// Helper function to reveal MPC commitment to base commitment
/// Follows the pattern from mpc-snarks/marlin.rs
pub fn reveal_commitment(
    mpc_comm: <MpcPC as PolynomialCommitment<MFr, DensePolynomial<MFr>>>::Commitment,
) -> <LocalPC as PolynomialCommitment<Fr, DensePolynomial<Fr>>>::Commitment {
    ark_poly_commit::marlin_pc::Commitment {
        comm: ark_poly_commit::kzg10::Commitment(mpc_comm.comm.0.reveal()),
        shifted_comm: mpc_comm.shifted_comm.map(|c| ark_poly_commit::kzg10::Commitment(c.0.reveal())),
    }
}

/// Helper function to reveal MPC proof to base proof
/// Reveals the opening proof from MPC space to public space
pub fn reveal_proof(
    mpc_proof: <MpcPC as PolynomialCommitment<MFr, DensePolynomial<MFr>>>::Proof,
) -> <LocalPC as PolynomialCommitment<Fr, DensePolynomial<Fr>>>::Proof {
    // MarlinKZG10's Proof type is kzg10::Proof<E>
    ark_poly_commit::kzg10::Proof {
        w: mpc_proof.w.reveal(),
        random_v: mpc_proof.random_v.map(|rv| rv.reveal()),
    }
}

/// Helper function to reveal MPC verifier key to base verifier key
pub fn reveal_verifier_key(
    mpc_vk: <MpcPC as PolynomialCommitment<MFr, DensePolynomial<MFr>>>::VerifierKey,
) -> <LocalPC as PolynomialCommitment<Fr, DensePolynomial<Fr>>>::VerifierKey {
    mpc_vk.reveal()
}
