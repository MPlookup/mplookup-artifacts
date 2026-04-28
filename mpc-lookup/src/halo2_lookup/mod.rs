//! # Halo2 Lookup Argument Module
//!
//! This module provides an implementation roadmap and placeholder structures for the Halo2 lookup
//! argument proof generation and verification in the MPC setting.
//!
//! ## Overview
//!
//! The Halo2 lookup argument is a zero-knowledge proof system that proves every value in a query
//! vector **A** (denoted as `f` in the existing code) exists in a table vector **S** (denoted as `t`).
//!
//! ### Current Integration
//!
//! The existing `mpc-lookup` crate already implements `secure_oblivious_lookup_permutation` which
//! computes the permuted vectors f' and t'. This module builds on top of that to add cryptographic
//! proofs.
//!
//! ## Module Structure
//!
//! - `types`: Core data structures (Proof, Commitments, Challenges, etc.) - **COMPLETE**
//! - `prover`: Proof generation implementation - **TODO**
//! - `verifier`: Proof verification implementation - **TODO**
//! - `transcript`: Fiat-Shamir transcript (wraps FiatShamirRng from mpc-plonk) - **TODO**
//! - `utils`: Helper functions (Lagrange basis, selectors, etc.) - **TODO**
//!
//! ## Reusable Components from Existing Codebase
//!
//! **IMPORTANT**: This implementation maximally reuses existing code:
//!
//! ### From ark_poly_commit (already in dependencies):
//! - **PolynomialCommitment trait**: Use directly for KZG10/IPA (no custom wrapper needed)
//! - **KZG10**: Recommended scheme compatible with mpc-plonk
//!
//! ### From mpc-plonk/src/util.rs:
//! - **FiatShamirRng**: Fiat-Shamir transcript with ChaChaRng - wrap in transcript.rs
//! - **shift()**: Polynomial shift f(aX) from f(X) - for Z(ωX), A'(ω^(-1)X)
//! - **interpolate()**: Lagrange interpolation - for polynomial construction
//!
//! ### From algebra/poly (arkworks):
//! - **Radix2EvaluationDomain**: FFT/IFFT and domain operations
//! - **evaluate_all_lagrange_coefficients()**: Computes ALL Lagrange basis polynomials at once
//! - **evaluate_vanishing_polynomial()**: Z_H(x) = x^n - 1
//!
//! ### From mpc-snarks/src/marlin.rs:
//! - **MPC reveal patterns**: How to reveal() commitments and proofs correctly
//!
//! **What still needs implementation** (core lookup logic only):
//! - Selector polynomials q_blind, q_last (trivial vector construction)
//! - Grand product Z computation (core algorithm, ~20 lines)
//! - Constraint checking functions (5 functions, field operations only)
//!
//! ## Usage Example
//!
//! ```rust,ignore
//! use mpc_lookup::halo2_lookup::*;
//! use ark_poly_commit::kzg10::KZG10;
//! use ark_bls12_377::{Bls12_377, Fr};
//! use blake2::Blake2s;
//!
//! // Setup configuration
//! let config = LookupConfig::new(table_size, query_size, 5);
//!
//! // Use KZG10 from ark_poly_commit directly
//! let pp = KZG10::<Bls12_377, DensePolynomial<Fr>>::setup(domain_size, false, &mut rng).unwrap();
//! let (ck, vk) = KZG10::trim(&pp, domain_size).unwrap();
//!
//! // Generate proof (integrates with existing permutation computation)
//! let mut transcript = LookupTranscript::<Blake2s>::new(b"lookup");
//! let mut prover = LookupProver::new(config, ck, a_values, s_values);
//! let proof = prover.prove(&mut transcript);
//!
//! // Verify proof
//! let mut transcript = LookupTranscript::<Blake2s>::new(b"lookup");
//! let verifier = LookupVerifier::new(verification_key, vk);
//! assert!(verifier.verify(&proof, &mut transcript));
//! ```
//!
//! ## Implementation Roadmap
//!
//! **For detailed step-by-step AI prompts, see AI_IMPLEMENTATION_PROMPTS.md**
//!
//! ### Phase 1: Utils & Transcript (0.5-1 day)
//! 1. Implement utils.rs (mostly wrappers around existing code)
//! 2. Wrap FiatShamirRng in transcript.rs
//!
//! ### Phase 2: Prover Core (7-10 days)
//! 1. Integrate with existing `secure_oblivious_lookup_permutation`
//! 2. Implement grand product Z computation (CORE ALGORITHM)
//! 3. Generate commitments using ark_poly_commit
//! 4. Compute evaluations and opening proofs
//!
//! ### Phase 3: Verifier Core (4-6 days)
//! 1. Re-derive challenges from proof (must match prover)
//! 2. Check 5 constraint types (permutation, subset, initial, idempotence)
//! 3. Verify opening proofs using ark_poly_commit
//!
//! ### Phase 4: Testing & Integration (5-7 days)
//! 1. End-to-end tests with valid and invalid queries
//! 2. MPC compatibility tests
//! 3. Performance benchmarks
//! 4. Examples and documentation
//!
//! **Total: 17-24 days (2.5-3.5 weeks)**
//!
//! ## Core Algorithm (Grand Product)
//!
//! The grand product Z encodes the lookup relation:
//!
//! ```text
//! Z[0] = 1
//! for i in 0..num_usable_rows:
//!     Z[i+1] = Z[i] * ((A'[i] + β)(S'[i] + γ)) / ((A[i] + β)(S[i] + γ))
//! ```
//!
//! If A ⊆ S, this product telescopes correctly. If not, the product won't satisfy constraints.
//!
//! ## References
//!
//! - [Halo2 Book - Lookup Argument](https://zcash.github.io/halo2/design/proving-system/lookup.html)
//! - [Halo2 Source - Prover](https://github.com/zcash/halo2/blob/main/halo2_proofs/src/plonk/lookup/prover.rs)
//! - [Halo2 Source - Verifier](https://github.com/zcash/halo2/blob/main/halo2_proofs/src/plonk/lookup/verifier.rs)
//! - [ark-poly-commit Documentation](https://docs.rs/ark-poly-commit/)
//! - **AI_IMPLEMENTATION_PROMPTS.md**: Step-by-step prompts for AI coding agents

pub mod types;
pub mod prover;
pub mod verifier;
pub mod transcript;
pub mod utils;
pub mod polynomial;
pub mod sorting;
pub mod permutation;

// Re-export main types
pub use types::{
    LookupConfig,
    LookupChallenges,
    LookupCommitments,
    LookupOpenings,
    LookupProof,
    LookupVerificationKey,
};

pub use prover::LookupProver;
pub use prover::{reveal_commitment, reveal_proof, reveal_verifier_key};
pub use verifier::LookupVerifier;
pub use transcript::LookupTranscript;

// Re-export utility functions
pub use utils::{
    compute_lagrange_basis_l0,
    compute_selector_polynomials,
};

// Re-export core lookup functions
pub use sorting::{bitonic_sort, bitonic_sort_by_key};
pub use permutation::{
    secure_oblivious_lookup_permutation,
    verify_lookup_permutation_outputs_in_plaintext,
};

// Re-export ark_poly_commit for convenience
// Developers should use ark_poly_commit::PolynomialCommitment directly
pub use ark_poly_commit::PolynomialCommitment;
