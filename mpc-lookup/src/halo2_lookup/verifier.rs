//! # Verifier Implementation for Halo2 Lookup Argument
//!
//! This module implements the proof verification algorithm using ark_poly_commit directly.
//!
//! ## Verification Strategy
//!
//! The verifier checks 5 types of constraints:
//!
//! 1. **Permutation Constraint**: Grand product recurrence holds
//!    - `Z(ωζ) * (A'(ζ) + β) * (S'(ζ) + γ) = Z(ζ) * (A(ζ) + β) * (S(ζ) + γ)`
//!    - Active on usable rows only (disabled by q_blind and q_last)
//!
//! 2. **Subset Constraint**: A' is sorted/grouped
//!    - `(A'(ζ) - S'(ζ)) * (A'(ζ) - A'(ω⁻¹ζ)) = 0`
//!    - Means: A'(ζ) equals S'(ζ) OR equals previous A' value
//!
//! 3. **Initial Constraints**: Boundary conditions at first row
//!    - `ℓ_0(ζ) * (A'(ζ) - S'(ζ)) = 0` (A' matches S' at start)
//!    - `ℓ_0(ζ) * (1 - Z(ζ)) = 0` (Z starts at 1)
//!
//! 4. **Idempotence Constraint**: Z closure
//!    - `q_last(ζ) * (Z(ζ)² - Z(ζ)) = 0`
//!    - Ensures Z(last) ∈ {0, 1}
//!
//! 5. **Opening Proofs**: Polynomial evaluations are correct
//!    - Use ark_poly_commit::PolynomialCommitment::check()
//!
//! **All checks are simple field operations** (no exponentiations except in opening verification).
//!
//! See AI_IMPLEMENTATION_PROMPTS.md Step 4 for complete implementation with algorithms.

use ark_ff::FftField;
use ark_poly_commit::PolynomialCommitment;
use digest::Digest;

use super::transcript::LookupTranscript;
use super::types::{LookupChallenges, LookupOpenings, LookupProof, LookupVerificationKey};

/// Verifier for the lookup argument.
///
/// # Generic Parameters
///
/// - `F`: Field type
/// - `P`: Polynomial type
/// - `PC`: Polynomial commitment scheme from ark_poly_commit
/// - `D`: Digest for transcript
pub struct LookupVerifier<F, P, PC, D> 
where
    F: FftField,
    P: ark_poly::UVPolynomial<F>,
    PC: PolynomialCommitment<F, P>,
    D: Digest,
{
    vk: LookupVerificationKey<F, PC::Commitment>,
    pcs_params: PC::VerifierKey,
    debug: bool,
    _phantom_f: std::marker::PhantomData<F>,
    _phantom_p: std::marker::PhantomData<P>,
    _phantom_d: std::marker::PhantomData<D>,
}

impl<F, P, PC, D> LookupVerifier<F, P, PC, D> 
where
    F: FftField,
    P: ark_poly::UVPolynomial<F>,
    PC: PolynomialCommitment<F, P>,
    D: Digest,
{
    pub fn new(vk: LookupVerificationKey<F, PC::Commitment>, pcs_params: PC::VerifierKey, debug: bool) -> Self {
        Self { 
            vk, 
            pcs_params,
            debug,
            _phantom_f: std::marker::PhantomData,
            _phantom_p: std::marker::PhantomData,
            _phantom_d: std::marker::PhantomData,
        }
    }

    pub fn derive_challenges(
        &self,
        proof: &LookupProof<F, PC::Commitment, PC::Proof>,
        transcript: &mut LookupTranscript<D>,
    ) -> LookupChallenges<F> 
    where
        F: ark_ff::PrimeField + ark_ff::PubUniformRand,
    {
        // Append commitments in SAME order as prover
        transcript.append_commitment(b"com_a", &proof.commitments.com_a);
        transcript.append_commitment(b"com_s", &proof.commitments.com_s);
        transcript.append_commitment(b"com_a_prime", &proof.commitments.com_a_prime);
        transcript.append_commitment(b"com_s_prime", &proof.commitments.com_s_prime);
        
        let (beta, gamma) = transcript.derive_beta_gamma();
        
        transcript.append_commitment(b"com_z", &proof.commitments.com_z);
        
        // CRITICAL FIX: Append quotient commitment before deriving zeta (must match prover!)
        transcript.append_commitment(b"com_q", &proof.commitments.com_q);
        
        let zeta = transcript.derive_zeta();
        
        LookupChallenges::new(beta, gamma, zeta)
    }

    pub fn check_permutation_constraint(
        &self,
        openings: &LookupOpenings<F, PC::Proof>,
        challenges: &LookupChallenges<F>,
    ) -> bool {
        if self.debug {
            eprintln!("\n=== VERIFIER: Checking ALL Constraints via Quotient Polynomial ===");
        }
        
        // According to Halo2 lookup argument, the quotient polynomial Q(X) encodes ALL constraints:
        // 1. Permutation: (1 - (q_last + q_blind)) * [Z(ωX)(A'(X)+β)(S'(X)+γ) - Z(X)(A(X)+β)(S(X)+γ)]
        // 2. Subset: (1 - (q_last + q_blind)) * (A'(X) - S'(X)) * (A'(X) - A'(ω⁻¹X))
        // 3. Initial A': ℓ_0(X) * (A'(X) - S'(X))
        // 4. Initial Z: ℓ_0(X) * (1 - Z(X))
        // 5. Idempotence: q_last(X) * (Z(X)² - Z(X))
        //
        // All sum to: constraint_poly(X) = Q(X) * Z_H(X)
        //
        // At challenge point ζ, we verify: sum of all constraints = Q(ζ) * Z_H(ζ)
        
        if self.debug {
            eprintln!("β = {}", challenges.beta);
            eprintln!("γ = {}", challenges.gamma);
            eprintln!("ζ = {}", challenges.zeta);
            
            eprintln!("\nOpening values:");
            eprintln!("  A(ζ) = {}", openings.a_zeta);
            eprintln!("  S(ζ) = {}", openings.s_zeta);
            eprintln!("  A'(ζ) = {}", openings.a_prime_zeta);
            eprintln!("  S'(ζ) = {}", openings.s_prime_zeta);
            eprintln!("  Z(ζ) = {}", openings.z_zeta);
            eprintln!("  Z(ωζ) = {}", openings.z_omega_zeta);
            eprintln!("  A'(ω⁻¹ζ) = {}", openings.a_prime_omega_inv_zeta);
            eprintln!("  q_blind(ζ) = {}", openings.q_blind_zeta);
            eprintln!("  q_last(ζ) = {}", openings.q_last_zeta);
            eprintln!("  ℓ_0(ζ) = {}", openings.l_0_zeta);
            eprintln!("  Q(ζ) = {}", openings.q_zeta);
        }
        
        // Compute selector = 1 - (q_last + q_blind)
        let selector = F::one() - (openings.q_last_zeta + openings.q_blind_zeta);
        
        // Constraint 1: Permutation
        if self.debug {
            eprintln!("\nComputing constraint 1: Permutation");
        }
        let a_plus_beta = openings.a_zeta + challenges.beta;
        let s_plus_gamma = openings.s_zeta + challenges.gamma;
        let a_prime_plus_beta = openings.a_prime_zeta + challenges.beta;
        let s_prime_plus_gamma = openings.s_prime_zeta + challenges.gamma;
        
        let perm_lhs = openings.z_omega_zeta * a_prime_plus_beta * s_prime_plus_gamma;
        let perm_rhs = openings.z_zeta * a_plus_beta * s_plus_gamma;
        let perm_constraint = perm_lhs - perm_rhs;
        let constraint_1 = selector * perm_constraint;
        if self.debug {
            eprintln!("  Permutation constraint = {}", constraint_1);
        }
        
        // Constraint 2: Subset
        if self.debug {
            eprintln!("\nComputing constraint 2: Subset");
        }
        let a_prime_minus_s_prime = openings.a_prime_zeta - openings.s_prime_zeta;
        let a_prime_minus_a_prime_shifted = openings.a_prime_zeta - openings.a_prime_omega_inv_zeta;
        let subset_constraint = a_prime_minus_s_prime * a_prime_minus_a_prime_shifted;
        let constraint_2 = selector * subset_constraint;
        if self.debug {
            eprintln!("  Subset constraint = {}", constraint_2);
        }
        
        // Constraint 3: Initial A'
        if self.debug {
            eprintln!("\nComputing constraint 3: Initial A'");
        }
        let constraint_3 = openings.l_0_zeta * a_prime_minus_s_prime;
        if self.debug {
            eprintln!("  Initial A' constraint = {}", constraint_3);
        }
        
        // Constraint 4: Initial Z
        if self.debug {
            eprintln!("\nComputing constraint 4: Initial Z");
        }
        let one_minus_z = F::one() - openings.z_zeta;
        let constraint_4 = openings.l_0_zeta * one_minus_z;
        if self.debug {
            eprintln!("  Initial Z constraint = {}", constraint_4);
        }
        
        // Constraint 5: Idempotence
        if self.debug {
            eprintln!("\nComputing constraint 5: Idempotence");
        }
        let z_squared = openings.z_zeta * openings.z_zeta;
        let z_squared_minus_z = z_squared - openings.z_zeta;
        let constraint_5 = openings.q_last_zeta * z_squared_minus_z;
        if self.debug {
            eprintln!("  Idempotence constraint = {}", constraint_5);
        }
        
        // Sum all constraints
        let total_constraint = constraint_1 + constraint_2 + constraint_3 + constraint_4 + constraint_5;
        if self.debug {
            eprintln!("\nTotal constraint value = {}", total_constraint);
        }
        
        // Compute Z_H(ζ) = ζ^domain_size - 1
        let zeta_to_n = challenges.zeta.pow(&[self.vk.config.domain_size as u64]);
        let z_h_zeta = zeta_to_n - F::one();
        
        // The quotient check: total_constraint = Q(ζ) * Z_H(ζ)
        let expected_constraint = openings.q_zeta * z_h_zeta;
        
        if self.debug {
            eprintln!("\nQuotient polynomial verification:");
            eprintln!("  Sum of all constraints = {}", total_constraint);
            eprintln!("  Z_H(ζ) = ζ^{} - 1 = {}", self.vk.config.domain_size, z_h_zeta);
            eprintln!("  Q(ζ) * Z_H(ζ) = {}", expected_constraint);
            eprintln!("  Match: {}", total_constraint == expected_constraint);
        }
        
        if total_constraint == expected_constraint {
            if self.debug {
                eprintln!("✓ All constraints verified via quotient polynomial!");
            }
            true
        } else {
            if self.debug {
                eprintln!("✗ Quotient polynomial check failed!");
                eprintln!("  Difference: {}", total_constraint - expected_constraint);
            }
            false
        }
    }

    pub fn check_subset_constraint(
        &self,
        openings: &LookupOpenings<F, PC::Proof>,
    ) -> bool {
        if self.debug {
            eprintln!("\n=== VERIFIER: Subset Constraint Check ===");
        }
        
        // Active on usable rows
        let selector = F::one() - (openings.q_last_zeta + openings.q_blind_zeta);
        if self.debug {
            eprintln!("selector = {}", selector);
            eprintln!("q_last(ζ) = {}", openings.q_last_zeta);
            eprintln!("q_blind(ζ) = {}", openings.q_blind_zeta);
        }
        
        // Constraint: (A'(ζ) - S'(ζ)) * (A'(ζ) - A'(ω⁻¹ζ)) = 0
        // Meaning: A'(ζ) equals S'(ζ) OR equals previous A' value
        if self.debug {
            eprintln!("A'(ζ) = {}", openings.a_prime_zeta);
            eprintln!("S'(ζ) = {}", openings.s_prime_zeta);
            eprintln!("A'(ω⁻¹ζ) = {}", openings.a_prime_omega_inv_zeta);
        }
        
        let diff1 = openings.a_prime_zeta - openings.s_prime_zeta;
        let diff2 = openings.a_prime_zeta - openings.a_prime_omega_inv_zeta;
        if self.debug {
            eprintln!("A'(ζ) - S'(ζ) = {}", diff1);
            eprintln!("A'(ζ) - A'(ω⁻¹ζ) = {}", diff2);
        }
        
        let constraint_value = diff1 * diff2;
        if self.debug {
            eprintln!("constraint_value = {}", constraint_value);
        }
        
        let result = selector * constraint_value;
        if self.debug {
            eprintln!("selector * constraint_value = {}", result);
            eprintln!("is_zero? {}", result.is_zero());
        }
        
        result.is_zero()
    }

    pub fn check_initial_constraints(
        &self,
        openings: &LookupOpenings<F, PC::Proof>,
    ) -> bool {
        // At first row: ℓ_0(ζ) * (A'(ζ) - S'(ζ)) = 0
        let init_subset = openings.l_0_zeta * (openings.a_prime_zeta - openings.s_prime_zeta);

        // At first row: ℓ_0(ζ) * (1 - Z(ζ)) = 0 (ensures Z starts at 1)
        let init_z = openings.l_0_zeta * (F::one() - openings.z_zeta);
        
        if self.debug {
            eprintln!("DEBUG: l_0_zeta = {}", openings.l_0_zeta);
            eprintln!("DEBUG: a_prime_zeta = {}", openings.a_prime_zeta);
            eprintln!("DEBUG: s_prime_zeta = {}", openings.s_prime_zeta);
            eprintln!("DEBUG: z_zeta = {}", openings.z_zeta);
            eprintln!("DEBUG: init_subset = {}, is_zero = {}", init_subset, init_subset.is_zero());
            eprintln!("DEBUG: init_z = {}, is_zero = {}", init_z, init_z.is_zero());
        }
        
        init_subset.is_zero() && init_z.is_zero()
    }

    pub fn check_idempotence_constraint(
        &self,
        openings: &LookupOpenings<F, PC::Proof>,
    ) -> bool {
        // At last usable row: q_last(ζ) * (Z(ζ)² - Z(ζ)) = 0
        // Ensures Z(last) ∈ {0, 1}
        let z_squared = openings.z_zeta * openings.z_zeta;
        let constraint_value = openings.q_last_zeta * (z_squared - openings.z_zeta);
        
        constraint_value.is_zero()
    }

    pub fn verify_openings(
        &self,
        proof: &LookupProof<F, PC::Commitment, PC::Proof>,
        challenges: &LookupChallenges<F>,
    ) -> bool {
        use ark_poly_commit::LabeledCommitment;
        
        let zeta = challenges.zeta;
        let omega = self.vk.omega;
        let omega_zeta = omega * zeta;
        let omega_inv = omega.inverse().expect("omega should be invertible");
        let omega_inv_zeta = zeta * omega_inv;
        
        // CRITICAL: Use std::iter::once for single-polynomial openings at different points
        use std::iter::once;
        
        // Wrap commitments in LabeledCommitment for verification
        // CRITICAL: Must specify degree_bound matching setup's supported_degrees!
        let degree_bound = Some(self.vk.config.domain_size - 1);
        let labeled_a = LabeledCommitment::new("a".into(), proof.commitments.com_a.clone(), degree_bound);
        let labeled_s = LabeledCommitment::new("s".into(), proof.commitments.com_s.clone(), degree_bound);
        let labeled_a_prime = LabeledCommitment::new("a_prime".into(), proof.commitments.com_a_prime.clone(), degree_bound);
        let labeled_s_prime = LabeledCommitment::new("s_prime".into(), proof.commitments.com_s_prime.clone(), degree_bound);
        let labeled_z = LabeledCommitment::new("z".into(), proof.commitments.com_z.clone(), degree_bound);
        let labeled_quotient = LabeledCommitment::new("quotient".into(), proof.commitments.com_q.clone(), None); // Quotient has no enforced bound
        
        // BATCH VERIFICATION: Verify all 6 polynomial openings at ζ in ONE call
        if self.debug {
            eprintln!("\n=== VERIFIER: Checking BATCH Opening Proof at ζ ===");
            eprintln!("Batching verification of 6 polynomials (A, S, A', S', Z, Q)");
        }
        
        // Opening challenge for batch verification
        let opening_challenge = F::one();
        
        // Collect all commitments and evaluations for batch verification at ζ
        let commitments_at_zeta = vec![
            &labeled_a,
            &labeled_s,
            &labeled_a_prime,
            &labeled_s_prime,
            &labeled_z,
            &labeled_quotient,
        ];
        
        let evaluations_at_zeta = vec![
            proof.openings.a_zeta,
            proof.openings.s_zeta,
            proof.openings.a_prime_zeta,
            proof.openings.s_prime_zeta,
            proof.openings.z_zeta,
            proof.openings.q_zeta,
        ];
        
        if self.debug {
            eprintln!("Opening challenge: {}", opening_challenge);
            eprintln!("Number of polynomials in batch: {}", commitments_at_zeta.len());
            eprintln!("  ζ = {}", zeta);
            eprintln!("  A(ζ) = {}", proof.openings.a_zeta);
            eprintln!("  S(ζ) = {}", proof.openings.s_zeta);
            eprintln!("  A'(ζ) = {}", proof.openings.a_prime_zeta);
            eprintln!("  S'(ζ) = {}", proof.openings.s_prime_zeta);
            eprintln!("  Z(ζ) = {}", proof.openings.z_zeta);
            eprintln!("  Q(ζ) = {}", proof.openings.q_zeta);
        }
        
        // Batch verify all 6 polynomial openings at ζ
        let check_batch = PC::check(
            &self.pcs_params,
            commitments_at_zeta.into_iter(),
            &zeta,
            evaluations_at_zeta.into_iter(),
            &proof.openings.proof_zeta,
            opening_challenge,
            None
        );
        
        if self.debug {
            eprintln!("\nBatch opening proof result: {:?}", check_batch);
        }
        
        if !check_batch.unwrap_or(false) {
            if self.debug {
                eprintln!("✗ Batch opening proof verification FAILED");
            }
            return false;
        } else {
            if self.debug {
                eprintln!("✓ Batch opening proof verified successfully for all 6 polynomials at ζ");
            }
        }
        
        // Verify opening at ωζ for Z
        let opening_challenge_omega = F::one(); // Use F::one() for single commitments
        let check_omega = PC::check(
            &self.pcs_params,
            once(&labeled_z),  // Use once() like mpc-plonk
            &omega_zeta,
            once(proof.openings.z_omega_zeta),  // Use once() like mpc-plonk
            &proof.openings.proof_omega_zeta,
            opening_challenge_omega,
            None
        );
        
        if !check_omega.unwrap_or(false) {
            if self.debug {
                eprintln!("✗ Opening proof verification at ωζ FAILED");
            }
            return false;
        }
        
        // Verify opening at ω⁻¹ζ for A'
        let opening_challenge_omega_inv = F::one(); // Use F::one() for single commitments
        let check_omega_inv = PC::check(
            &self.pcs_params,
            once(&labeled_a_prime),  // Use once() like mpc-plonk
            &omega_inv_zeta,
            once(proof.openings.a_prime_omega_inv_zeta),  // Use once() like mpc-plonk
            &proof.openings.proof_omega_inv_zeta,
            opening_challenge_omega_inv,
            None
        );
        
        if !check_omega_inv.unwrap_or(false) {
            if self.debug {
                eprintln!("✗ Opening proof verification at ω⁻¹ζ FAILED");
            }
            return false;
        }
        
        true
    }

    pub fn verify(
        &self,
        proof: &LookupProof<F, PC::Commitment, PC::Proof>,
        transcript: &mut LookupTranscript<D>,
    ) -> bool 
    where
        F: ark_ff::PrimeField + ark_ff::PubUniformRand,
    {
        if self.debug {
            eprintln!("\n=== VERIFIER: Starting Verification ===");
        }
        
        // 1. Derive challenges (same as prover)
        let challenges = self.derive_challenges(proof, transcript);
        
        if self.debug {
            eprintln!("Challenges derived successfully");
            eprintln!("  β = {}", challenges.beta);
            eprintln!("  γ = {}", challenges.gamma);
            eprintln!("  ζ = {}", challenges.zeta);
        }
        
        // 2. First verify opening proofs to ensure all values are correct
        if self.debug {
            eprintln!("\nStep 1: Verifying opening proofs...");
        }
        if !self.verify_openings(proof, &challenges) {
            if self.debug {
                eprintln!("DEBUG: Opening proof verification failed");
            }
            return false;
        }
        if self.debug {
            eprintln!("✓ Opening proofs verified successfully");
        }
        
        // 3. Check all constraints via the quotient polynomial
        if self.debug {
            eprintln!("\nStep 2: Checking all constraints via quotient polynomial...");
        }
        if !self.check_permutation_constraint(&proof.openings, &challenges) {
            if self.debug {
                eprintln!("DEBUG: Quotient polynomial constraint check failed");
            }
            return false;
        }
        
        if self.debug {
            eprintln!("\n✓ All verification checks passed!");
        }
        
        true
    }
}
