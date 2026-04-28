//! # Core Data Structures for Halo2 Lookup Argument
//!
//! This module defines all the data structures needed for the lookup protocol.

use ark_ff::Field;

/// Configuration parameters for the Halo2 lookup argument.
///
/// These parameters define the size and security properties of the lookup protocol.
///
/// # Fields
///
/// - `domain_size`: Size of the evaluation domain (must be power of 2)
/// - `num_blinding_rows`: Number of rows filled with random values for zero-knowledge
/// - `num_usable_rows`: Index of the last usable row + 1 (= domain_size - num_blinding_rows - 1).
///   Usable rows are at indices 0 to num_usable_rows-1, boundary row is at index num_usable_rows.
///
/// # Example
///
/// ```rust,ignore
/// let config = LookupConfig::new(100, 50, 5);
/// assert_eq!(config.domain_size, 128); // Next power of 2 after max(100, 50)
/// assert_eq!(config.num_blinding_rows, 5);
/// assert_eq!(config.num_usable_rows, 122); // 128 - 5 - 1 (boundary row index)
/// ```
#[derive(Debug, Clone)]
pub struct LookupConfig {
    /// Size of the domain (must be power of 2)
    /// This should be ≥ max(|S|, |A|)
    pub domain_size: usize,

    /// Number of blinding rows for zero-knowledge (typically 3-6)
    /// The last `num_blinding_rows` rows are filled with random values
    pub num_blinding_rows: usize,

    /// Number of usable rows = domain_size - num_blinding_rows - 1
    pub num_usable_rows: usize,
}

impl LookupConfig {
    /// Create a new lookup configuration.
    ///
    /// # Arguments
    ///
    /// - `table_size`: Size of the lookup table S
    /// - `query_size`: Number of queries in A
    /// - `num_blinding_rows`: Number of blinding rows for ZK (typically 3-6)
    ///
    /// # Returns
    ///
    /// A LookupConfig with appropriate domain size (next power of 2)
    ///
    /// # Implementation Note
    ///
    /// The domain size is set to the next power of 2 that can accommodate the data.
    /// The usable rows (for constraint enforcement) are then reduced by blinding rows
    /// and one boundary row.
    pub fn new(table_size: usize, query_size: usize, num_blinding_rows: usize) -> Self {
        // Domain must accommodate:
        // - All data in usable rows: max(table_size, query_size)
        // - Plus one boundary row
        // - Plus blinding rows for zero-knowledge
        let min_data_size = table_size.max(query_size);
        let min_domain_size = min_data_size + num_blinding_rows + 1;
        let domain_size = min_domain_size.next_power_of_two();
        let num_usable_rows = domain_size - num_blinding_rows - 1;

        Self {
            domain_size,
            num_blinding_rows,
            num_usable_rows,
        }
    }
}

/// Challenges derived via Fiat-Shamir heuristic for the lookup argument.
///
/// These random field elements are computed by hashing the transcript of commitments,
/// ensuring the protocol is non-interactive and secure.
///
/// # Protocol Flow
///
/// 1. Prover commits to A, S, A', S'
/// 2. Derive β and γ from transcript
/// 3. Prover commits to Z
/// 4. Derive ζ from transcript
///
/// # Security
///
/// Challenges must be derived after the prover commits to prevent selective choice attacks.
/// The Fiat-Shamir transform ensures these challenges are as secure as random oracle outputs.
#[derive(Debug, Clone)]
pub struct LookupChallenges<F: Field> {
    /// Challenge β for grand product computation
    /// Used in: Z[i+1] = Z[i] * ((A'[i] + β) * (S'[i] + γ)) / ((A[i] + β) * (S[i] + γ))
    pub beta: F,

    /// Challenge γ for grand product computation
    /// Used in: Z[i+1] = Z[i] * ((A'[i] + β) * (S'[i] + γ)) / ((A[i] + β) * (S[i] + γ))
    pub gamma: F,

    /// Random evaluation point ζ for polynomial openings
    /// All polynomials are evaluated at this point for verification
    pub zeta: F,
}

impl<F: Field> LookupChallenges<F> {
    /// Create challenges with given values (primarily for testing).
    pub fn new(beta: F, gamma: F, zeta: F) -> Self {
        Self { beta, gamma, zeta }
    }
}

/// Polynomial commitments for the lookup proof.
///
/// This structure holds all the polynomial commitments created by the prover.
///
/// # Commitment Scheme
///
/// The generic parameter `G` represents the commitment type, which depends on the
/// polynomial commitment scheme used (e.g., KZG commitments are elliptic curve points).
///
/// # Fields
///
/// - `com_a`: Commitment to A(X) - query polynomial
/// - `com_s`: Commitment to S(X) - table polynomial
/// - `com_a_prime`: Commitment to A'(X) - permuted query polynomial
/// - `com_s_prime`: Commitment to S'(X) - permuted table polynomial
/// - `com_z`: Commitment to Z(X) - grand product polynomial
/// - `com_q_blind`: Optional commitment to q_blind(X) - blinding selector
/// - `com_q_last`: Optional commitment to q_last(X) - last row selector
///
/// Note: `com_q_blind` and `com_q_last` may be None if these are fixed/public polynomials.
#[derive(Debug, Clone)]
pub struct LookupCommitments<G> {
    /// Commitment to A(X) - query polynomial
    pub com_a: G,

    /// Commitment to S(X) - table polynomial
    pub com_s: G,

    /// Commitment to A'(X) - permuted query polynomial
    pub com_a_prime: G,

    /// Commitment to S'(X) - permuted table polynomial
    pub com_s_prime: G,

    /// Commitment to Z(X) - grand product polynomial
    pub com_z: G,

    /// Commitment to Q(X) - quotient polynomial for permutation constraint
    /// Q(X) satisfies: Z(ωX) * (A(X)+β) * (S(X)+γ) - Z(X) * (A'(X)+β) * (S'(X)+γ) = Q(X) * Z_H(X)
    pub com_q: G,

    /// Commitment to q_blind(X) - blinding selector (often fixed/public)
    pub com_q_blind: Option<G>,

    /// Commitment to q_last(X) - last row selector (often fixed/public)
    pub com_q_last: Option<G>,
}

/// Opening proofs for polynomial evaluations at challenge points.
///
/// These proofs convince the verifier that the claimed evaluations are correct
/// without revealing the entire polynomials.
///
/// # Evaluation Points
///
/// - ζ: Main evaluation point
/// - ωζ: Next row (for Z polynomial)
/// - ω^(-1)ζ: Previous row (for A' polynomial)
///
/// # Fields
///
/// All evaluations are at specific points derived from the challenge ζ and domain generator ω.
#[derive(Debug, Clone)]
pub struct LookupOpenings<F: Field, P> {
    /// Evaluation of A(X) at ζ
    pub a_zeta: F,

    /// Evaluation of S(X) at ζ
    pub s_zeta: F,

    /// Evaluation of A'(X) at ζ
    pub a_prime_zeta: F,

    /// Evaluation of S'(X) at ζ
    pub s_prime_zeta: F,

    /// Evaluation of Z(X) at ζ
    pub z_zeta: F,

    /// Evaluation of Z(X) at ωζ (next row)
    /// Used for: Z(ωζ) * (A'(ζ) + β) * (S'(ζ) + γ) = Z(ζ) * (A(ζ) + β) * (S(ζ) + γ)
    pub z_omega_zeta: F,

    /// Evaluation of A'(X) at ω^(-1)ζ (previous row)
    /// Used for: (A'(ζ) - S'(ζ)) * (A'(ζ) - A'(ω^(-1)ζ)) = 0
    pub a_prime_omega_inv_zeta: F,

    /// Evaluation of q_blind(X) at ζ
    pub q_blind_zeta: F,

    /// Evaluation of q_last(X) at ζ
    pub q_last_zeta: F,

    /// Evaluation of ℓ_0(X) at ζ (Lagrange basis for first row)
    /// Used for initial conditions: ℓ_0(ζ) * (A'(ζ) - S'(ζ)) = 0 and ℓ_0(ζ) * (1 - Z(ζ)) = 0
    pub l_0_zeta: F,

    /// Evaluation of Q(X) at ζ (quotient polynomial for permutation constraint)
    /// Used to verify: Z(ωζ) * (A(ζ)+β) * (S(ζ)+γ) - Z(ζ) * (A'(ζ)+β) * (S'(ζ)+γ) = Q(ζ) * Z_H(ζ)
    pub q_zeta: F,

    /// Opening proof for evaluations at ζ
    /// Proves correctness of: a_zeta, s_zeta, a_prime_zeta, s_prime_zeta, z_zeta, q_zeta, q_blind_zeta, q_last_zeta
    pub proof_zeta: P,

    /// Opening proof for evaluations at ωζ
    /// Proves correctness of: z_omega_zeta
    pub proof_omega_zeta: P,

    /// Opening proof for evaluations at ω^(-1)ζ
    /// Proves correctness of: a_prime_omega_inv_zeta
    pub proof_omega_inv_zeta: P,
}

/// Complete lookup argument proof.
///
/// This structure contains all information sent from prover to verifier.
///
/// # Protocol
///
/// The proof consists of:
/// 1. Commitments to all polynomials (A, S, A', S', Z, optionally q_blind, q_last)
/// 2. Evaluations of polynomials at challenge points
/// 3. Opening proofs that the evaluations are correct
///
/// # Generic Parameters
///
/// - `F`: Field type (e.g., BLS12-377 scalar field)
/// - `G`: Commitment type (depends on polynomial commitment scheme)
/// - `P`: Opening proof type (depends on polynomial commitment scheme)
#[derive(Debug, Clone)]
pub struct LookupProof<F: Field, G, P> {
    /// Polynomial commitments
    pub commitments: LookupCommitments<G>,

    /// Polynomial openings and proofs
    pub openings: LookupOpenings<F, P>,
}

/// Verification key for the lookup argument.
///
/// Contains public parameters and fixed polynomials needed for verification.
///
/// # Contents
///
/// - Domain parameters (size, generator ω)
/// - Fixed polynomial commitments (q_blind, q_last if applicable)
/// - Configuration parameters
///
/// # Usage
///
/// The verification key is generated during setup and distributed to all verifiers.
/// It must match the parameters used by the prover.
#[derive(Debug, Clone)]
pub struct LookupVerificationKey<F: Field, G> {
    /// Domain generator ω
    /// The multiplicative generator of the evaluation domain
    pub omega: F,

    /// Domain size (must be power of 2)
    pub domain_size: usize,

    /// Commitment to q_blind(X) if fixed
    /// If None, q_blind must be computed/committed by prover
    pub com_q_blind: Option<G>,

    /// Commitment to q_last(X) if fixed
    /// If None, q_last must be computed/committed by prover
    pub com_q_last: Option<G>,

    /// Configuration parameters
    pub config: LookupConfig,
}

#[cfg(test)]
mod tests {
    use ark_bls12_377::Fr;
    use super::*;
    use ark_ff::One;

    #[test]
    fn test_lookup_config_creation() {
        let config = LookupConfig::new(100, 50, 5);
        assert_eq!(config.domain_size, 128); // Next power of 2 after max(100, 50)
        assert_eq!(config.num_blinding_rows, 5);
        assert_eq!(config.num_usable_rows, 122); // 128 - 5 - 1
    }

    #[test]
    fn test_lookup_config_edge_cases() {
        // Test with power of 2: table_size=64 needs domain big enough for 64+3+1=68 -> 128
        let config1 = LookupConfig::new(64, 32, 3);
        assert_eq!(config1.domain_size, 128);  // Next power of 2 after 64+3+1=68
        assert_eq!(config1.num_usable_rows, 124); // 128 - 3 - 1

        // Test with small values: 1+2+1=4 -> domain=4
        let config2 = LookupConfig::new(1, 1, 2);
        assert_eq!(config2.domain_size, 4);  // Next power of 2 after 1+2+1=4
        assert_eq!(config2.num_usable_rows, 1); // 4 - 2 - 1

        // Test with large blinding: 16+10+1=27 -> 32
        let config3 = LookupConfig::new(16, 8, 10);
        assert_eq!(config3.domain_size, 32);  // Next power of 2 after 16+10+1=27
        assert_eq!(config3.num_usable_rows, 21); // 32 - 10 - 1
    }

    #[test]
    fn test_challenges_creation() {
        let beta = Fr::from(42u64);
        let gamma = Fr::from(123u64);
        let zeta = Fr::from(456u64);
        
        let challenges = LookupChallenges::new(beta, gamma, zeta);
        assert_eq!(challenges.beta, Fr::from(42u64));
        assert_eq!(challenges.gamma, Fr::from(123u64));
        assert_eq!(challenges.zeta, Fr::from(456u64));
    }
}
