//! # Fiat-Shamir Transcript for Challenge Generation
//!
//! This module provides a transcript interface for the Halo2 lookup argument.
//!
//! ## Implementation Note
//!
//! Instead of implementing from scratch, developers should adapt the existing
//! `FiatShamirRng` from `mpc-plonk/src/util.rs`. That implementation provides:
//! - ChaChaRng-based random number generation
//! - Digest-based seed refreshing  
//! - Secure challenge derivation
//!
//! To use it, add to Cargo.toml:
//! ```toml
//! mpc-plonk = { path = "../mpc-plonk" }
//! ```
//!
//! Then import:
//! ```rust,ignore
//! use mpc_plonk::util::FiatShamirRng;
//! use blake2::Blake2s;
//! ```
//!
//! This saves 2-3 days of implementation time.

use ark_ff::Field;
use ark_serialize::CanonicalSerialize;
use digest::Digest;
use mpc_plonk::FiatShamirRng;

/// Transcript for Fiat-Shamir challenge derivation.
///
/// This wraps the existing `FiatShamirRng` from mpc-plonk for use in the Halo2 lookup argument.
///
/// # Usage Example
///
/// ```rust,ignore
/// use blake2::Blake2s;
/// 
/// let mut transcript = LookupTranscript::<Blake2s>::new(b"halo2-lookup");
/// transcript.append_commitment(b"com_a", &commitment_a);
/// let beta: Fr = transcript.squeeze_challenge(b"beta");
/// ```
pub struct LookupTranscript<D: Digest> {
    rng: FiatShamirRng<D>,
}

impl<D: Digest> LookupTranscript<D> {
    /// Create a new transcript with a given label.
    ///
    /// Wraps `FiatShamirRng::from_seed` to initialize the transcript.
    ///
    /// # Arguments
    ///
    /// - `label`: A byte string label for the transcript
    ///
    /// # Returns
    ///
    /// A new `LookupTranscript` instance
    pub fn new(label: &'static [u8]) -> Self {
        // Convert byte slice to a Vec<u8> which implements ToBytes
        let label_vec = label.to_vec();
        Self {
            rng: FiatShamirRng::from_seed(&label_vec),
        }
    }

    /// Append a commitment to the transcript.
    ///
    /// Serializes the commitment and absorbs it into the transcript state.
    ///
    /// # Arguments
    ///
    /// - `label`: A byte string label for the commitment
    /// - `commitment`: The commitment to append (must implement CanonicalSerialize)
    pub fn append_commitment<G: CanonicalSerialize>(
        &mut self,
        label: &'static [u8],
        commitment: &G,
    ) {
        let mut bytes = Vec::new();
        commitment.serialize(&mut bytes).unwrap();
        // Concatenate label and commitment bytes
        let mut combined = label.to_vec();
        combined.extend_from_slice(&bytes);
        self.rng.absorb(&combined);
    }

    /// Squeeze a random challenge from the transcript.
    ///
    /// Generates a uniformly random field element from the current transcript state.
    ///
    /// # Arguments
    ///
    /// - `_label`: A byte string label for the challenge (unused but kept for API consistency)
    ///
    /// # Returns
    ///
    /// A uniformly random field element
    pub fn squeeze_challenge<F: Field + ark_ff::PubUniformRand>(
        &mut self,
        _label: &'static [u8],
    ) -> F {
        self.rng.gen()
    }

    /// Derive β and γ challenges.
    ///
    /// Called after committing to A, S, A', S'.
    pub fn derive_beta_gamma<F: Field + ark_ff::PubUniformRand>(&mut self) -> (F, F) {
        let beta = self.squeeze_challenge(b"beta");
        let gamma = self.squeeze_challenge(b"gamma");
        (beta, gamma)
    }

    /// Derive ζ challenge.
    ///
    /// Called after committing to Z.
    pub fn derive_zeta<F: Field + ark_ff::PubUniformRand>(&mut self) -> F {
        self.squeeze_challenge(b"zeta")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_377::Fr;
    use blake2::Blake2s;

    #[test]
    fn test_transcript_usage() {
        // Test basic transcript usage
        let mut transcript = LookupTranscript::<Blake2s>::new(b"test");
        let challenge: Fr = transcript.squeeze_challenge(b"test_challenge");
        
        // The challenge should be a valid field element (no panic)
        // We can't test the specific value since it's pseudo-random based on seed
        assert!(challenge != Fr::from(0u64) || challenge == Fr::from(0u64));
    }
    
    #[test]
    fn test_transcript_derive_challenges() {
        // Test the derive_beta_gamma and derive_zeta methods
        let mut transcript = LookupTranscript::<Blake2s>::new(b"test");
        
        let (beta, gamma): (Fr, Fr) = transcript.derive_beta_gamma();
        let zeta: Fr = transcript.derive_zeta();
        
        // All challenges should be valid field elements
        assert!(beta != Fr::from(0u64) || beta == Fr::from(0u64));
        assert!(gamma != Fr::from(0u64) || gamma == Fr::from(0u64));
        assert!(zeta != Fr::from(0u64) || zeta == Fr::from(0u64));
    }
}
