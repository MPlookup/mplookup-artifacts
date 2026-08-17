//! # MPC Lookup - Secure Oblivious Lookup Argument
//!
//! This crate provides a secure multi-party computation (MPC) implementation of
//! oblivious lookup arguments using the Halo2 protocol.
//!
//! ## Main Components
//!
//! - **Oblivious Sorting**: Secure bitonic sort for MPC field elements
//! - **Lookup Permutation**: Secure permutation algorithm for lookup arguments
//! - **Halo2 Proof System**: Zero-knowledge proof generation and verification
//!
//! ## Example
//!
//! ```rust,ignore
//! use mpc_lookup::{secure_oblivious_lookup_permutation, verify_lookup_permutation_outputs_in_plaintext};
//! use mpc_algebra::honest_but_curious::MpcField;
//!
//! // Run secure lookup permutation
//! let (t_prime, f_prime) = secure_oblivious_lookup_permutation(n, m, &t, &f, false);
//! ```

// Halo2 Lookup Argument Module - contains all core functionality
pub mod halo2_lookup;

// Re-export the main lookup functions at the crate root for backward compatibility
pub use halo2_lookup::sorting::{bitonic_sort, bitonic_sort_by_key, radix_sort, radix_sort_by_key};
pub use halo2_lookup::permutation::{
    secure_oblivious_lookup_permutation,
    naive_secure_oblivious_lookup_permutation,
    verify_lookup_permutation_outputs_in_plaintext,
};
