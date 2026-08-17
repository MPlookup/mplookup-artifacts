//! Boolean secret sharing operations for MPC bit-level computations
//!
//! This module implements boolean secret sharing primitives including:
//! - EdaBits (Extended Doubly-Authenticated Bits)
//! - DaBits (Doubly-Authenticated Bits)
//! - Boolean operations (AND, OR, XOR, NOT)
//! - Bit adders with carry propagation
//! - B2A (Boolean-to-Arithmetic) conversion
//!
//! These primitives are used to implement secure comparison and bit decomposition
//! for secret-shared field elements.
//!
//! # Broadcast Optimization
//!
//! This module has been optimized to reduce the number of MPC broadcast operations:
//!
//! ## Single Operations (Original)
//! - `beaver_bitwise_and`: 2 broadcasts (one for each masked value)
//! - `beaver_bitwise_or`: 2 broadcasts (for the underlying AND operation)
//!
//! ## Batch Operations (Optimized)
//! - `batch_beaver_bitwise_and(n ops)`: 2 broadcasts total (all masked values batched)
//! - `batch_beaver_bitwise_or(n ops)`: 2 broadcasts total (all operations batched)
//!
//! ## Impact on Secure Comparison
//! For an n-bit comparison (e.g., 253 bits for BLS12-377):
//! - **Prefix-OR computation**: O(n) broadcasts (sequential dependency - cannot be batched)
//! - **Selection bits computation**: Batched to O(1) broadcasts
//! - **Final OR reduction**: Batched using tree reduction to O(log n) broadcasts
//! - **Overall**: O(n) broadcasts (dominated by prefix-OR), but with significant constant factor improvements
//! 
//! While we can't eliminate the O(n) sequential prefix-OR broadcasts, we achieve:
//! - Batched AND operations for selection bits: n operations in 2 broadcasts
//! - Batched OR operations for final reduction: tree-based reduction to log(n) broadcasts
//! - Batched AND operations in bit decomposition: 3 operations in 2 broadcasts per bit
//!
//! ## Impact on Bitonic Sort
//! For sorting n elements with O(n log² n) comparisons:
//! - Each comparison benefits from batched operations within selection and reduction
//! - Prefix-OR remains O(bit_length) but with optimized underlying operations
//! - Overall significant reduction in total broadcasts through batching within each comparison
//!
//! This optimization maintains all security guarantees while significantly improving
//! communication efficiency in MPC protocols.

use ark_ff::{PrimeField, Field, BigInteger};
use rand::Rng;
use mpc_net::{MpcNet, MpcMultiNet as Net};
use super::field::FieldShare;
use std::marker::PhantomData;
use derivative::Derivative;

/// EdaBit (Extended Doubly-Authenticated Bit) structure
/// 
/// Provides a random value `r` in both:
/// - Arithmetic domain: secret share `[r]_p ∈ F_p`
/// - Boolean domain: secret shares of individual bits `[r_i] ∈ {0,1}` for i ∈ [0, bit_length)
///
/// Used for bit decomposition of secret-shared field values.
#[derive(Clone, Debug)]
pub struct EdaBit<F: PrimeField, S: FieldShare<F>> {
    /// Arithmetic share of the random value in field F_p
    pub arith_share: S,
    /// Boolean shares of individual bits (LSB first)
    pub bool_shares: Vec<bool>,
    /// Phantom data for F
    _phantom: PhantomData<F>,
}

/// DaBit (Doubly-Authenticated Bit) structure using Prio+ protocol
///
/// Provides a random bit `r ∈ {0,1}` in both:
/// - Arithmetic domain: secret share `[r]_p ∈ F_p` (0 or 1)
/// - Boolean domain: secret share `[r]_2 ∈ {0,1}`
///
/// Used for B2A (Boolean-to-Arithmetic) conversion.
#[derive(Clone, Debug)]
pub struct DaBitPrioPlus<F: Field, S: FieldShare<F>> {
    /// Arithmetic share of the bit in field F_p (represents 0 or 1)
    pub arith_share: S,
    /// Boolean share of the bit
    pub bool_share: bool,
    /// Phantom data for F
    _phantom: PhantomData<F>,
}

/// Boolean secret sharing operations
pub struct BooleanOps;

impl BooleanOps {
    /// Bitwise NOT operation for boolean secret shares
    /// 
    /// In additive secret sharing over Z_2:
    /// - Party 0 (king): returns NOT of share
    /// - Other parties: return share unchanged
    #[inline]
    pub fn bitwise_not(share: bool) -> bool {
        if Net::am_king() {
            !share
        } else {
            share
        }
    }

    /// Bitwise XOR operation for boolean secret shares
    /// 
    /// XOR is locally computable: share1 XOR share2
    #[inline]
    pub fn bitwise_xor(left_share: bool, right_share: bool) -> bool {
        left_share ^ right_share
    }

    /// Secure bitwise AND using Beaver multiplication triple
    ///
    /// Uses the same pattern as field multiplication but in boolean domain.
    /// Requires communication for opening masked values.
    pub fn beaver_bitwise_and(
        left_share: bool,
        right_share: bool,
        beaver_source: &mut impl BooleanBeaverSource,
    ) -> bool {
        let (a_share, b_share, c_share) = beaver_source.bool_triple();
        
        // Compute and open d = left XOR a
        let d_share = left_share ^ a_share;
        let d = Self::open_bool(d_share);
        
        // Compute and open e = right XOR b
        let e_share = right_share ^ b_share;
        let e = Self::open_bool(e_share);
        
        // Compute result = c XOR (d AND b) XOR (e AND a) XOR (d AND e)
        let mut result = c_share;
        if d {
            result ^= b_share;
        }
        if e {
            result ^= a_share;
        }
        if Net::am_king() && d && e {
            result ^= true;
        }
        
        result
    }

    /// Secure bitwise OR derived from AND and XOR
    ///
    /// OR(a, b) = a XOR b XOR AND(a, b)
    pub fn beaver_bitwise_or(
        left_share: bool,
        right_share: bool,
        beaver_source: &mut impl BooleanBeaverSource,
    ) -> bool {
        let and_result = Self::beaver_bitwise_and(left_share, right_share, beaver_source);
        left_share ^ right_share ^ and_result
    }

    /// Batch secure bitwise AND using Beaver multiplication triples
    ///
    /// Performs multiple AND operations in parallel, batching the open operations
    /// to reduce communication rounds from O(n) to O(1).
    ///
    /// # Arguments
    /// * `left_shares` - Vector of left boolean shares
    /// * `right_shares` - Vector of right boolean shares
    /// * `beaver_source` - Source for boolean triples
    ///
    /// # Returns
    /// Vector of AND results corresponding to each pair
    pub fn batch_beaver_bitwise_and(
        left_shares: &[bool],
        right_shares: &[bool],
        beaver_source: &mut impl BooleanBeaverSource,
    ) -> Vec<bool> {
        assert_eq!(left_shares.len(), right_shares.len(), 
                   "Left and right share vectors must have the same length");
        
        let n = left_shares.len();
        let mut results = Vec::with_capacity(n);
        
        // Generate triples for all operations
        let mut triples = Vec::with_capacity(n);
        for _ in 0..n {
            triples.push(beaver_source.bool_triple());
        }
        
        // Compute all masked values
        let mut d_shares = Vec::with_capacity(n);
        let mut e_shares = Vec::with_capacity(n);
        
        for i in 0..n {
            let (a_share, b_share, _) = triples[i];
            d_shares.push(left_shares[i] ^ a_share);
            e_shares.push(right_shares[i] ^ b_share);
        }
        
        // Batch open all masked values at once (2 broadcasts instead of 2n)
        let d_opened = Self::batch_open_bool(&d_shares);
        let e_opened = Self::batch_open_bool(&e_shares);
        
        // Compute results locally
        for i in 0..n {
            let (a_share, b_share, c_share) = triples[i];
            let d = d_opened[i];
            let e = e_opened[i];
            
            let mut result = c_share;
            if d {
                result ^= b_share;
            }
            if e {
                result ^= a_share;
            }
            if Net::am_king() && d && e {
                result ^= true;
            }
            
            results.push(result);
        }
        
        results
    }

    /// Batch secure bitwise OR derived from batch AND and XOR
    ///
    /// OR(a, b) = a XOR b XOR AND(a, b)
    ///
    /// Performs multiple OR operations in parallel, reducing communication overhead.
    pub fn batch_beaver_bitwise_or(
        left_shares: &[bool],
        right_shares: &[bool],
        beaver_source: &mut impl BooleanBeaverSource,
    ) -> Vec<bool> {
        assert_eq!(left_shares.len(), right_shares.len(), 
                   "Left and right share vectors must have the same length");
        
        let and_results = Self::batch_beaver_bitwise_and(left_shares, right_shares, beaver_source);
        
        left_shares.iter()
            .zip(right_shares.iter())
            .zip(and_results.iter())
            .map(|((&l, &r), &and_res)| l ^ r ^ and_res)
            .collect()
    }

    /// Open a boolean secret share
    fn open_bool(share: bool) -> bool {
        use crate::channel::MpcSerNet;
        let all_shares = Net::broadcast(&share);
        // Fix for single-party mode: if broadcast returns empty (no other parties),
        // return the local share directly (king holds the full value in additive sharing).
        if all_shares.is_empty() {
            return share;
        }
        all_shares.into_iter().fold(false, |acc, s| acc ^ s)
    }

    /// Batch open multiple boolean secret shares
    /// 
    /// Opens multiple boolean shares in a single broadcast operation,
    /// reducing communication rounds from O(n) to O(1) where n is the number of shares.
    pub fn batch_open_bool(shares: &[bool]) -> Vec<bool> {
        use crate::channel::MpcSerNet;

        // Convert slice to Vec for broadcasting
        let shares_vec: Vec<bool> = shares.to_vec();

        // Broadcast all shares at once - returns Vec<Vec<bool>> where each inner Vec is from one party
        let all_party_shares: Vec<Vec<bool>> = Net::broadcast(&shares_vec);

        // Fix for single-party mode: if broadcast returns empty (no other parties),
        // return the local shares directly (king holds the full value in additive sharing).
        if all_party_shares.is_empty() {
            return shares_vec;
        }

        // XOR all shares for each position
        let n_shares = shares.len();
        let mut results = vec![false; n_shares];

        for party_shares in all_party_shares {
            for (i, &share) in party_shares.iter().enumerate() {
                results[i] ^= share;
            }
        }

        results
    }

    /// Bit adder: add two secret-shared bit sequences with carry propagation
    ///
    /// Adds `left_bits` + `right_bits` and returns result with optional overflow bit.
    ///
    /// # Arguments
    /// * `left_bits` - Secret-shared bits (LSB first)
    /// * `right_bits` - Secret-shared bits (LSB first)
    /// * `preserve_overflow` - If true, append carry bit to result
    /// * `beaver_source` - Source for boolean triples
    /// 
    /// # Optimization
    /// Batches AND operations within each bit position to reduce communication rounds.
    pub fn bits_add(
        left_bits: &[bool],
        right_bits: &[bool],
        preserve_overflow: bool,
        beaver_source: &mut impl BooleanBeaverSource,
    ) -> Vec<bool> {
        assert_eq!(left_bits.len(), right_bits.len(), "Bit sequences must have same length");
        
        let bit_count = left_bits.len();
        let mut result_bits = Vec::with_capacity(if preserve_overflow { bit_count + 1 } else { bit_count });
        let mut carry_share = false;
        
        for i in 0..bit_count {
            let a_bit_share = left_bits[i];
            let b_bit_share = right_bits[i];
            
            // Sum bit: a XOR b XOR carry
            result_bits.push(a_bit_share ^ b_bit_share ^ carry_share);
            
            // Update carry: (a AND b) OR (a AND carry) OR (b AND carry)
            // Batch the three AND operations
            let and_left = vec![a_bit_share, a_bit_share, b_bit_share];
            let and_right = vec![b_bit_share, carry_share, carry_share];
            let and_results = Self::batch_beaver_bitwise_and(&and_left, &and_right, beaver_source);
            
            let t1 = and_results[0]; // a AND b
            let t2 = and_results[1]; // a AND carry
            let t3 = and_results[2]; // b AND carry
            
            // Now compute (t1 OR t2) OR t3
            // First compute t1 OR t2, then OR with t3
            let t4 = Self::beaver_bitwise_or(t1, t2, beaver_source); // t1 OR t2
            carry_share = Self::beaver_bitwise_or(t4, t3, beaver_source); // (t1 OR t2) OR t3
        }
        
        if preserve_overflow {
            result_bits.push(carry_share);
        }
        
        result_bits
    }

    /// Bit adder with constant: add secret-shared bits to public constant
    ///
    /// More efficient than general bit adder since one operand is public.
    ///
    /// # Arguments
    /// * `left_bits` - Secret-shared bits (LSB first)
    /// * `right_const` - Public constant bits (LSB first)
    /// * `preserve_overflow` - If true, append carry bit to result
    /// * `beaver_source` - Source for boolean triples
    /// 
    /// # Optimization
    /// Batches operations where possible to reduce communication rounds.
    pub fn bits_add_const(
        left_bits: &[bool],
        right_const: &[bool],
        preserve_overflow: bool,
        beaver_source: &mut impl BooleanBeaverSource,
    ) -> Vec<bool> {
        assert_eq!(left_bits.len(), right_const.len(), "Bit sequences must have same length");
        
        let bit_count = left_bits.len();
        let mut result_bits = Vec::with_capacity(if preserve_overflow { bit_count + 1 } else { bit_count });
        let mut carry_share = false;
        
        for i in 0..bit_count {
            let a_bit_share = left_bits[i];
            // Convert public bit to secret share (only king has it)
            let b_bit = if Net::am_king() { right_const[i] } else { false };
            
            // Sum bit: a XOR b XOR carry
            result_bits.push(a_bit_share ^ b_bit ^ carry_share);
            
            // Update carry based on whether constant bit is 0 or 1
            if !right_const[i] {
                // If b=0: carry = a AND carry
                let t2 = Self::beaver_bitwise_and(a_bit_share, carry_share, beaver_source);
                carry_share = t2;
            } else {
                // If b=1: carry = a OR carry
                // We can batch the AND operation in the OR formula
                let t2 = Self::beaver_bitwise_and(a_bit_share, carry_share, beaver_source);
                // carry = a XOR carry XOR (a AND carry)
                carry_share = a_bit_share ^ carry_share ^ t2;
            }
        }
        
        if preserve_overflow {
            result_bits.push(carry_share);
        }
        
        result_bits
    }
}

/// Trait for providing boolean beaver triples (for AND operations)
pub trait BooleanBeaverSource: Clone {
    /// Generate a boolean multiplication triple (a, b, c) where c = a AND b
    fn bool_triple(&mut self) -> (bool, bool, bool);
}

/// Trait for providing EdaBits (Extended Doubly-Authenticated Bits)
pub trait EdaBitSource<F: PrimeField, S: FieldShare<F>>: Clone {
    /// Generate an EdaBit with specified bit length
    fn edabit(&mut self, bit_length: usize) -> EdaBit<F, S>;
}

/// Trait for providing DaBits (Doubly-Authenticated Bits)
pub trait DaBitSource<F: Field, S: FieldShare<F>>: Clone {
    /// Generate a DaBit
    fn dabit(&mut self) -> DaBitPrioPlus<F, S>;
}

/// Dummy boolean beaver source that generates trivial triples (1,1,1)
/// 
/// **WARNING: This is INSECURE and only for testing or trusted dealer scenarios.**
/// 
/// In a proper MPC protocol, beaver triples should be:
/// - Random values (a, b, c) where c = a AND b
/// - Generated securely in preprocessing phase
/// - Unknown to any single party
/// 
/// The current implementation generates predictable triples which violates
/// security requirements. For production use, implement proper preprocessing.
#[derive(Clone, Copy, Debug, Default)]
pub struct DummyBooleanBeaverSource;

impl BooleanBeaverSource for DummyBooleanBeaverSource {
    fn bool_triple(&mut self) -> (bool, bool, bool) {
        // TODO: Replace with secure preprocessing in production
        // For now, this generates insecure but functional triples
        // King generates (1,1,1), others generate (0,0,0)
        // This makes the AND computation work but is not cryptographically secure
        let a = if Net::am_king() { true } else { false };
        let b = if Net::am_king() { true } else { false };
        let c = if Net::am_king() { true } else { false };
        (a, b, c)
    }
}

/// Dummy EdaBit source following the multiplication triple pattern
/// 
/// **WARNING: This is INSECURE and only for testing or trusted dealer scenarios.**
/// 
/// In a proper MPC protocol, EdaBits should be:
/// - Random values with consistent arithmetic and boolean representations
/// - Generated securely in preprocessing phase
/// - Unknown to any single party
/// 
/// The current implementation generates predictable shares which violates
/// security requirements. For production use, implement proper preprocessing.
#[derive(Derivative)]
#[derivative(Default(bound = ""), Clone(bound = ""), Copy(bound = ""))]
pub struct DummyEdaBitSource<F, S> {
    _field: PhantomData<F>,
    _share: PhantomData<S>,
}

impl<F: PrimeField, S: FieldShare<F>> EdaBitSource<F, S> for DummyEdaBitSource<F, S> {
    fn edabit(&mut self, bit_length: usize) -> EdaBit<F, S> {
        // Follow multiplication triple pattern: no network transfer
        // Each party generates shares locally
        if Net::am_king() {
            // King generates random value and its bits
            let r = F::rand(&mut rand::thread_rng());
            let r_bigint = r.into_repr();
            let bits: Vec<bool> = r_bigint.to_bits_le()[..bit_length].to_vec();
            
            // Generate arithmetic share locally (king has value)
            let arith_share = S::from_add_shared(r);
            
            // Generate boolean shares locally (king has bit values)
            let bool_shares = bits;
            
            EdaBit {
                arith_share,
                bool_shares,
                _phantom: PhantomData,
            }
        } else {
            // Non-king parties generate zero shares locally
            let arith_share = S::from_add_shared(F::zero());
            let bool_shares = vec![false; bit_length];
            
            EdaBit {
                arith_share,
                bool_shares,
                _phantom: PhantomData,
            }
        }
    }
}

/// Dummy DaBit source following the multiplication triple pattern
/// 
/// **WARNING: This is INSECURE and only for testing or trusted dealer scenarios.**
/// 
/// In a proper MPC protocol, DaBits should be:
/// - Random bits with consistent arithmetic and boolean representations
/// - Generated securely in preprocessing phase
/// - Unknown to any single party
/// 
/// The current implementation generates predictable shares which violates
/// security requirements. For production use, implement proper preprocessing.
#[derive(Derivative)]
#[derivative(Default(bound = ""), Clone(bound = ""), Copy(bound = ""))]
pub struct DummyDaBitSource<F, S> {
    _field: PhantomData<F>,
    _share: PhantomData<S>,
}

impl<F: Field, S: FieldShare<F>> DaBitSource<F, S> for DummyDaBitSource<F, S> {
    fn dabit(&mut self) -> DaBitPrioPlus<F, S> {
        // Follow multiplication triple pattern: no network transfer
        // Each party generates shares locally
        if Net::am_king() {
            // King generates random bit
            let bit: bool = rand::random();
            let field_val = if bit { F::one() } else { F::zero() };
            
            // Generate arithmetic share locally (king has value)
            let arith_share = S::from_add_shared(field_val);
            
            // Generate boolean share locally (king has bit)
            let bool_share = bit;
            
            DaBitPrioPlus {
                arith_share,
                bool_share,
                _phantom: PhantomData,
            }
        } else {
            // Non-king parties generate zero shares locally
            let arith_share = S::from_add_shared(F::zero());
            let bool_share = false;
            
            DaBitPrioPlus {
                arith_share,
                bool_share,
                _phantom: PhantomData,
            }
        }
    }
}

/// B2A (Boolean-to-Arithmetic) conversion using daBit
///
/// Converts a boolean share [b]_2 to an arithmetic share [b]_p using a daBit.
///
/// Protocol:
/// 1. Consume a daBit ([r]_p, [r]_2)
/// 2. Compute and reveal delta = [b]_2 XOR [r]_2
/// 3. If delta = 0: return [r]_p
///    If delta = 1: return 1 - [r]_p
pub fn b2a<F: Field, S: FieldShare<F>>(
    bool_share: bool,
    dabit: DaBitPrioPlus<F, S>,
) -> S {
    // Compute delta = bool_share XOR dabit.bool_share
    let delta_share = bool_share ^ dabit.bool_share;
    
    // Reveal delta
    let delta = BooleanOps::open_bool(delta_share);
    
    // Return appropriate arithmetic share
    if !delta {
        dabit.arith_share
    } else {
        // Compute 1 - [r]_p
        let mut one_share = S::from_public(F::one());
        one_share.sub(&dabit.arith_share);
        one_share
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_377::Fr;
    use crate::share::add::AdditiveFieldShare;

    #[test]
    fn test_bitwise_not() {
        // Test that NOT works correctly for party 0 and others
        let share = true;
        let result = BooleanOps::bitwise_not(share);
        
        // Party 0 should flip, others should not
        if Net::am_king() {
            assert_eq!(result, false);
        } else {
            assert_eq!(result, true);
        }
    }

    #[test]
    fn test_bitwise_xor() {
        assert_eq!(BooleanOps::bitwise_xor(true, true), false);
        assert_eq!(BooleanOps::bitwise_xor(true, false), true);
        assert_eq!(BooleanOps::bitwise_xor(false, true), true);
        assert_eq!(BooleanOps::bitwise_xor(false, false), false);
    }
}
