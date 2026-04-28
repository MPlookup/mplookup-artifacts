use ark_ff::bytes::{FromBytes, ToBytes};
use ark_ff::prelude::*;
use ark_serialize::{
    CanonicalDeserialize, CanonicalDeserializeWithFlags, CanonicalSerialize,
    CanonicalSerializeWithFlags,
};
//use ark_poly::univariate::{DensePolynomial,DenseOrSparsePolynomial};
use core::ops::*;
use std::cmp::Ord;
use std::fmt::{Debug, Display};
use std::hash::Hash;
use rand::Rng;

use super::BeaverSource;
use crate::Reveal;

pub trait FieldShare<F: Field>:
    Clone
    + Copy
    + Display
    + Debug
    + Send
    + Sync
    + Eq
    + Hash
    + Ord
    + CanonicalSerialize
    + CanonicalDeserialize
    + CanonicalSerializeWithFlags
    + CanonicalDeserializeWithFlags
    + UniformRand
    + ToBytes
    + FromBytes
    + 'static
    + Reveal<Base = F>
{
    fn open(&self) -> F {
        <Self as Reveal>::reveal(*self)
    }

    fn map_homo<FF: Field, SS: FieldShare<FF>, Fun: Fn(F) -> FF>(self, f: Fun) -> SS {
        SS::from_add_shared(f(self.unwrap_as_public()))
    }

    fn batch_open(selfs: impl IntoIterator<Item = Self>) -> Vec<F> {
        selfs.into_iter().map(|s| s.open()).collect()
    }

    fn add(&mut self, other: &Self) -> &mut Self;

    fn sub(&mut self, other: &Self) -> &mut Self {
        let mut t = other.clone();
        t.neg();
        t.add(&self);
        *self = t;
        self
    }

    fn neg(&mut self) -> &mut Self {
        self.scale(&-<F as ark_ff::One>::one())
    }

    fn shift(&mut self, other: &F) -> &mut Self;

    fn scale(&mut self, other: &F) -> &mut Self;

    fn mul<S: BeaverSource<Self, Self, Self>>(self, other: Self, source: &mut S) -> Self {
        let (mut x, mut y, z) = source.triple();
        //println!("Triple:\n *{}\n *{}\n *{}", x, y, z);
        let s = self;
        let o = other;
        // output: z - open(s + x)y - open(o + y)x + open(s + x)open(o + y)
        //         xy - sy - xy - ox - yx + so + sy + xo + xy
        //         so
        let sx = {
            let mut t = s;
            t.add(&x).open()
        };
        let oy = {
            let mut t = o;
            t.add(&y).open()
        };
        let mut result = z;
        result.sub(y.scale(&sx)).sub(x.scale(&oy)).shift(&(sx * oy));
        #[cfg(debug_assertions)]
        {
            let a = s.reveal();
            let b = o.reveal();
            let r = result.reveal();
            if a * b != r {
                println!("Bad multiplication!.\n{}\n*\n{}\n=\n{}", a, b, r);
                panic!("Bad multiplication");
            }
        }
        result
    }

    fn batch_mul<S: BeaverSource<Self, Self, Self>>(
        xs: Vec<Self>,
        ys: Vec<Self>,
        source: &mut S,
    ) -> Vec<Self> {
        let ss = xs;
        let os = ys;
        let (xs, ys, zs) = source.triples(ss.len());
        // output: z - open(s + x)y - open(o + y)x + open(s + x)open(o + y)
        //         xy - sy - xy - ox - yx + so + sy + xo + xy
        //         so
        let sxs = Self::batch_open(ss.into_iter().zip(xs.iter()).map(|(mut s, x)| {
            s.add(x);
            s
        }));
        let oys = Self::batch_open(os.into_iter().zip(ys.iter()).map(|(mut o, y)| {
            o.add(y);
            o
        }));
        zs.into_iter()
            .zip(ys.into_iter())
            .zip(xs.into_iter())
            .enumerate()
            .map(|(i, ((mut z, mut y), mut x))| {
                z.sub(y.scale(&sxs[i]))
                    .sub(x.scale(&oys[i]))
                    .shift(&(sxs[i] * oys[i]));
                z
            })
            .collect()
    }

    fn inv<S: BeaverSource<Self, Self, Self>>(self, source: &mut S) -> Self {
        let (x, mut y) = source.inv_pair();
        let xa = x.mul(self, source).open().inverse().unwrap();
        *y.scale(&xa)
    }

    fn batch_inv<S: BeaverSource<Self, Self, Self>>(xs: Vec<Self>, source: &mut S) -> Vec<Self> {
        let (bs, cs) = source.inv_pairs(xs.len());
        cs.into_iter()
            .zip(
                Self::batch_open(Self::batch_mul(xs, bs, source))
                    .into_iter()
                    .map(|i| i.inverse().unwrap()),
            )
            .map(|(mut c, i)| {
                c.scale(&i);
                c
            })
            .collect()
    }

    fn div<S: BeaverSource<Self, Self, Self>>(self, other: Self, source: &mut S) -> Self {
        let o_inv = other.inv(source);
        self.mul(o_inv, source)
    }

    fn batch_div<S: BeaverSource<Self, Self, Self>>(
        xs: Vec<Self>,
        ys: Vec<Self>,
        source: &mut S,
    ) -> Vec<Self> {
        Self::batch_mul(xs, Self::batch_inv(ys, source), source)
    }

    fn partial_products<S: BeaverSource<Self, Self, Self>>(x: Vec<Self>, src: &mut S) -> Vec<Self> {
        let n = x.len();
        let (m, m_inv): (Vec<Self>, Vec<Self>) = (0..(n + 1)).map(|_| src.inv_pair()).unzip();
        let mx = Self::batch_mul(m[..n].iter().cloned().collect(), x, src);
        let mxm = Self::batch_mul(mx, m_inv[1..].iter().cloned().collect(), src);
        let mut mxm_pub = Self::batch_open(mxm);
        for i in 1..mxm_pub.len() {
            let last = mxm_pub[i - 1];
            mxm_pub[i] *= &last;
        }
        let m0 = vec![m[0]; n];
        let mms = Self::batch_mul(m0, m_inv[1..].iter().cloned().collect(), src);
        let mut mms_inv = Self::batch_inv(mms, src);
        //let mms_pub = Self::batch_open(mms);
        for i in 0..mxm_pub.len() {
            mms_inv[i].scale(&mxm_pub[i]);
        }
        debug_assert!(mxm_pub.len() == n);
        mms_inv
    }

    fn univariate_div_qr<'a>(
        _num: DenseOrSparsePolynomial<Self>,
        _den: DenseOrSparsePolynomial<F>,
    ) -> Option<(
        DensePolynomial<Self>,
        DensePolynomial<Self>,
    )> {
        todo!("Implement generic poly div")
    }

    /// Secure equality check for shared values using Fermat's little theorem.
    /// Computes: c = a - b, d = c^(p-2), e = c * d
    /// Returns the shared value e (0 if equal, 1 if not equal).
    fn secure_neq<S: BeaverSource<Self, Self, Self>>(&self, other: &Self, source: &mut S) -> Self
    where
        F: PrimeField,
    {
        // Step 1: c = a - b
        let mut c = self.clone();
        c.sub(other);
        
        // Step 2: d = c^(p-2) using Fermat's little theorem
        // For a prime field, a^(p-1) = 1 (for a != 0), so a^(p-2) = a^(-1)
        // If c = 0, then c^(p-2) = 0
        // Create a copy of MODULUS and subtract 2
        let two = <F as PrimeField>::BigInt::from(2u64);
        let mut modulus_minus_two = F::Params::MODULUS;
        modulus_minus_two.sub_noborrow(&two);
        let d = Self::pow_bigint(c, &modulus_minus_two, source);
        
        // Step 3: e = c * d
        // If c != 0: e = c * c^(-1) = 1
        // If c = 0: e = 0 * 0 = 0
        let e = c.mul(d, source);
        
        // Return e (0 means equal, 1 means not equal)
        e
    }
    
    /// Helper function to compute self^exp using binary exponentiation
    fn pow_bigint<S: BeaverSource<Self, Self, Self>>(
        base: Self,
        exp: &<F as PrimeField>::BigInt,
        source: &mut S,
    ) -> Self
    where
        F: PrimeField,
    {
        // Convert BigInt to bits for binary exponentiation
        let bits = exp.to_bits_be();
        
        // Start with result = 1
        let mut result = Self::from_public(F::one());
        let mut current = base;
        
        // Binary exponentiation from LSB to MSB
        for bit in bits.iter().rev() {
            if *bit {
                result = result.mul(current, source);
            }
            current = current.mul(current, source);
        }
        
        result
    }

    /// Secure comparison for shared values.
    /// 
    /// Returns a secret-shared value representing the less-than bit:
    /// - Returns Self (0) if self >= other
    /// - Returns Self (1) if self < other
    /// 
    /// This is similar to `secure_neq` and does NOT reveal the comparison result,
    /// maintaining privacy of the compared values.
    /// 
    /// # Protocol Overview
    /// 
    /// This protocol is implemented following Qi, Huayi, et al. "VDORAM: Towards a Random Access Machine 
    /// with Both Public Verifiability and Distributed Obliviousness." Cryptology ePrint Archive (2025).
    /// 
    /// The comparison works by:
    /// 1. Bit-decomposing both values using edaBits
    /// 2. Comparing bits from MSB to LSB to find first differing position
    /// 3. Using secure selection logic to extract the comparison bit
    /// 4. Returning the result as a secret-shared field element (0 or 1)
    /// 
    /// # Implementation
    /// 
    /// Uses trusted dealer model where king (party 0) generates:
    /// - edaBits for bit decomposition
    /// - daBits for B2A conversion
    /// - Boolean triples for secure AND operations
    /// 
    /// # Security Notes
    /// 
    /// Unlike the standard `Ord` trait which reveals ordering, this function:
    /// - Does NOT reveal whether self < other
    /// - Returns a secret share that can be used in further MPC computations
    /// - Leaks no information about the actual values or their ordering
    /// 
    /// # References
    /// - Qi, Huayi, et al. "VDORAM" (2025)
    /// - Rotaru & Wood "MArBled Circuits: Mixing Arithmetic and Boolean Circuits" (2019)
    /// - Escudero et al. "Improved Primitives for MPC over Mixed Arithmetic-Boolean Circuits" (2020)
    /// - Damgård et al. "Unconditionally Secure Constant-Rounds Multi-Party Computation" (2006)
    /// - Addanki et al. "Prio+: Privacy Preserving Aggregate Statistics via Boolean Shares" (2022)
    fn secure_cmp(&self, other: &Self) -> Self
    where
        F: PrimeField,
    {
        use crate::share::bit_ops::{DummyEdaBitSource, EdaBitSource, DummyDaBitSource, DaBitSource, b2a, BooleanOps, DummyBooleanBeaverSource};
        use rand::thread_rng;
        
        let rng = &mut thread_rng();
        let bit_length = F::Params::MODULUS.num_bits() as usize;
        
        // Bit decompose self and other
        let self_bits = Self::bit_decompose(self.clone(), bit_length, rng);
        let other_bits = Self::bit_decompose(other.clone(), bit_length, rng);
        
        // Compare bits from MSB to LSB and return secret-shared result
        Self::compare_bits(&self_bits, &other_bits, rng)
    }
    
    /// Bit decompose a secret-shared field element
    /// 
    /// Returns boolean shares of the bits (LSB first)
    fn bit_decompose<R: Rng>(value: Self, bit_length: usize, rng: &mut R) -> Vec<bool>
    where
        F: PrimeField,
    {
        use crate::share::bit_ops::*;
        
        // Generate edaBit using DummyEdaBitSource
        let edabit = DummyEdaBitSource::<F, Self>::default().edabit(bit_length);
        
        // Step 1: Compute [c] = [a] - [b]
        let mut c_share = value;
        c_share.sub(&edabit.arith_share);
        
        // Step 2: Reveal c
        let c = c_share.open();
        
        // Step 3: Compute e = c + 2^BitSize - p in extended ring
        // Break down for clarity: e = c + 2^l - p (mod 2^{l+1})
        let two_pow_l = num_bigint::BigUint::from(1u32) << bit_length;
        let two_pow_l_plus_1: num_bigint::BigUint = &two_pow_l << 1;
        let p = num_bigint::BigUint::from_bytes_le(&F::Params::MODULUS.to_bytes_le());
        let c_bigint = num_bigint::BigUint::from_bytes_le(&c.into_repr().to_bytes_le());
        
        // Compute e = (c + 2^l - p) mod 2^{l+1}
        let offset = if p < two_pow_l { 
            &two_pow_l - &p 
        } else { 
            two_pow_l_plus_1.clone() - (&p - &two_pow_l)
        };
        let e_bigint = (&c_bigint + offset) % &two_pow_l_plus_1;
        
        // Bit decompose e (public)
        let e_bytes = e_bigint.to_bytes_le();
        let mut e_bits = Vec::with_capacity(bit_length + 1);
        for i in 0..=bit_length {
            let byte_idx = i / 8;
            let bit_idx = i % 8;
            if byte_idx < e_bytes.len() {
                e_bits.push((e_bytes[byte_idx] >> bit_idx) & 1 == 1);
            } else {
                e_bits.push(false);
            }
        }
        
        // Step 4: Add [b_bits] + e_bits using bit adder
        let mut b_bits_with_overflow = edabit.bool_shares.clone();
        b_bits_with_overflow.push(false); // Add dummy bit
        
        let mut beaver_source = DummyBooleanBeaverSource;
        let d_prime_bits = BooleanOps::bits_add_const(
            &b_bits_with_overflow,
            &e_bits,
            true,
            &mut beaver_source,
        );
        
        // Step 5: Extract overflow bit q = d'[bit_length]
        let q_share = d_prime_bits[bit_length];
        
        // Step 6: Compute NOT q
        let not_q_share = BooleanOps::bitwise_not(q_share);
        
        // Step 7: Prepare p_bits and multiply with NOT q
        let p_bits: Vec<bool> = (0..bit_length)
            .map(|i| {
                let byte_idx = i / 8;
                let bit_idx = i % 8;
                let p_bytes = p.to_bytes_le();
                if byte_idx < p_bytes.len() {
                    (p_bytes[byte_idx] >> bit_idx) & 1 == 1
                } else {
                    false
                }
            })
            .collect();
        
        // Multiply each public p bit with [NOT q] using secure AND
        // Note: p_bit is public, so we share it and then AND with secret not_q_share
        let optional_p_bits: Vec<bool> = p_bits
            .iter()
            .map(|&p_bit| {
                use mpc_net::{MpcNet, MpcMultiNet as Net};
                // Convert public bit to secret share for secure multiplication
                let p_bit_share = if Net::am_king() { p_bit } else { false };
                // Perform secure AND between public bit (as share) and secret NOT q
                BooleanOps::beaver_bitwise_and(p_bit_share, not_q_share, &mut beaver_source)
            })
            .collect();
        
        // Step 8: Add d_prime_bits[0..bit_length] + optional_p_bits
        let d_prime_trimmed = d_prime_bits[..bit_length].to_vec();
        let h_bits = BooleanOps::bits_add(
            &d_prime_trimmed,
            &optional_p_bits,
            false,
            &mut beaver_source,
        );
        
        // Return the boolean bit shares
        h_bits
    }
    
    /// Compare two bit sequences (LSB first) and return secret-shared less-than bit
    /// 
    /// Uses comparison circuit to find first differing bit from MSB.
    /// Returns a secret-shared value:
    /// - 0 if left >= right
    /// - 1 if left < right
    /// 
    /// **SECURITY NOTE:** This function does NOT reveal the comparison result.
    /// The returned value is a secret share that maintains privacy.
    /// 
    /// **OPTIMIZATION:** Uses batched boolean operations to reduce communication rounds.
    fn compare_bits<R: Rng>(left_bits: &[bool], right_bits: &[bool], rng: &mut R) -> Self
    where
        F: PrimeField,
    {
        use crate::share::bit_ops::*;
        
        let bit_count = left_bits.len();
        assert_eq!(bit_count, right_bits.len());
        
        let mut beaver_source = DummyBooleanBeaverSource;
        
        // Step 1: XOR bits to find differing positions (MSB to LSB)
        // XOR is local, no communication needed
        let mut xor_bits = Vec::with_capacity(bit_count);
        for i in (0..bit_count).rev() {
            xor_bits.push(BooleanOps::bitwise_xor(left_bits[i], right_bits[i]));
        }
        
        // Step 2: Build prefix-OR to find first differing bit
        // prefix_or[i] indicates if any bit from MSB to position i differs
        // NOTE: This MUST be sequential due to data dependencies - cannot be batched
        let mut prefix_or = vec![false; bit_count];
        prefix_or[0] = xor_bits[0];
        for i in 1..bit_count {
            prefix_or[i] = BooleanOps::beaver_bitwise_or(
                prefix_or[i - 1],
                xor_bits[i],
                &mut beaver_source,
            );
        }
        
        // Step 3: Compute selection bits (indicates THE first differing bit)
        // selection[0] = xor[0] (MSB differs)
        // selection[i] = xor[i] AND NOT prefix_or[i-1] (this is first difference)
        let mut selection_bits = Vec::with_capacity(bit_count);
        selection_bits.push(xor_bits[0]);
        
        if bit_count > 1 {
            // Batch all AND operations for selection bits
            let mut and_left_shares = Vec::with_capacity(bit_count - 1);
            let mut and_right_shares = Vec::with_capacity(bit_count - 1);
            
            for i in 1..bit_count {
                let not_prev = BooleanOps::bitwise_not(prefix_or[i - 1]);
                and_left_shares.push(xor_bits[i]);
                and_right_shares.push(not_prev);
            }
            
            // Batch AND operations
            let and_results = BooleanOps::batch_beaver_bitwise_and(
                &and_left_shares,
                &and_right_shares,
                &mut beaver_source,
            );
            
            for i in 1..bit_count {
                selection_bits.push(and_results[i - 1]);
            }
        }
        
        // Step 4: Use selection to pick the right bit value
        // result = OR of (selection[i] AND right[MSB-i])
        // Batch all AND operations, then fold with OR
        let mut selected_and_left = Vec::with_capacity(bit_count);
        let mut selected_and_right = Vec::with_capacity(bit_count);
        
        for i in 0..bit_count {
            let bit_idx = bit_count - 1 - i; // Convert back to LSB indexing
            selected_and_left.push(selection_bits[i]);
            selected_and_right.push(right_bits[bit_idx]);
        }
        
        // Batch AND operations for selection
        let selected_results = BooleanOps::batch_beaver_bitwise_and(
            &selected_and_left,
            &selected_and_right,
            &mut beaver_source,
        );
        
        // Fold with OR using tree-based reduction to minimize communication rounds
        let mut current_level = selected_results;
        
        while current_level.len() > 1 {
            let mut next_level = Vec::new();
            let pairs_count = current_level.len() / 2;
            
            if pairs_count > 0 {
                let mut or_left = Vec::with_capacity(pairs_count);
                let mut or_right = Vec::with_capacity(pairs_count);
                
                for i in 0..pairs_count {
                    or_left.push(current_level[i * 2]);
                    or_right.push(current_level[i * 2 + 1]);
                }
                
                let or_results = BooleanOps::batch_beaver_bitwise_or(
                    &or_left,
                    &or_right,
                    &mut beaver_source,
                );
                
                next_level.extend(or_results);
            }
            
            // Handle odd element
            if current_level.len() % 2 == 1 {
                next_level.push(*current_level.last().unwrap());
            }
            
            current_level = next_level;
        }
        
        let result_share = current_level[0];
        
        // Step 5: Convert result to arithmetic domain (without revealing)
        let dabit_result = DummyDaBitSource::<F, Self>::default().dabit();
        let arith_result = b2a(result_share, dabit_result);
        
        // Return the secret-shared less-than bit (0 or 1)
        arith_result
    }

    /// Batch secure comparison for multiple pairs of secret-shared values.
    /// 
    /// Compares multiple pairs efficiently by batching bit decompositions
    /// and comparison operations together, reducing communication rounds.
    /// 
    /// For each pair (left[i], right[i]), returns a secret share of:
    /// - 1 if left[i] < right[i]
    /// - 0 if left[i] >= right[i]
    /// 
    /// # Performance
    /// This function batches the reveal operations in bit decomposition,
    /// reducing communication from O(n) sequential reveals to O(1) batched reveals
    /// for n pairs. Within each comparison, existing batched boolean operations are used.
    /// 
    /// # Security
    /// All operations maintain MPC security - no values are revealed during computation.
    fn batch_secure_cmp(left_values: &[Self], right_values: &[Self]) -> Vec<Self>
    where
        F: PrimeField,
    {
        use crate::share::bit_ops::*;
        use rand::thread_rng;
        
        assert_eq!(left_values.len(), right_values.len(), "batch_secure_cmp: input lengths must match");
        
        if left_values.is_empty() {
            return vec![];
        }
        
        let bit_length = F::Params::MODULUS.num_bits() as usize;
        let n_pairs = left_values.len();
        
        // Step 1: Generate edaBits for all values (2 * n_pairs values total)
        let mut all_edabits = Vec::with_capacity(2 * n_pairs);
        let mut all_c_shares = Vec::with_capacity(2 * n_pairs);
        
        for i in 0..n_pairs {
            // Generate edaBit for left value
            let left_edabit = DummyEdaBitSource::<F, Self>::default().edabit(bit_length);
            let mut left_c_share = left_values[i];
            left_c_share.sub(&left_edabit.arith_share);
            all_c_shares.push(left_c_share);
            all_edabits.push(left_edabit);
            
            // Generate edaBit for right value
            let right_edabit = DummyEdaBitSource::<F, Self>::default().edabit(bit_length);
            let mut right_c_share = right_values[i];
            right_c_share.sub(&right_edabit.arith_share);
            all_c_shares.push(right_c_share);
            all_edabits.push(right_edabit);
        }
        
        // Step 2: Batch open all c values (KEY OPTIMIZATION: single communication round)
        let all_c_values = Self::batch_open(all_c_shares);
        
        // Step 3: Bit decompose all values using the opened c values
        let mut all_left_bits = Vec::with_capacity(n_pairs);
        let mut all_right_bits = Vec::with_capacity(n_pairs);
        
        for i in 0..n_pairs {
            let left_idx = i * 2;
            let right_idx = i * 2 + 1;
            
            // Process left value bits
            let left_bits = Self::bit_decompose_from_opened(
                &all_c_values[left_idx],
                &all_edabits[left_idx].bool_shares,
                bit_length,
            );
            all_left_bits.push(left_bits);
            
            // Process right value bits
            let right_bits = Self::bit_decompose_from_opened(
                &all_c_values[right_idx],
                &all_edabits[right_idx].bool_shares,
                bit_length,
            );
            all_right_bits.push(right_bits);
        }
        
        // Step 4: Perform comparison for each pair using batched operations
        let mut rng = thread_rng();
        let mut results = Vec::with_capacity(n_pairs);
        for i in 0..n_pairs {
            let result = Self::compare_bits(&all_left_bits[i], &all_right_bits[i], &mut rng);
            results.push(result);
        }
        
        results
    }
    
    /// Helper function to complete bit decomposition after the value has been opened.
    /// 
    /// This is used by batch_secure_cmp to separate the open operation from the rest
    /// of the decomposition, enabling batched reveals for better communication efficiency.
    /// 
    /// # Arguments
    /// * `c` - The opened masked value (c = value - edaBit.arith_share)
    /// * `b_bits` - The boolean shares from the edaBit
    /// * `bit_length` - Number of bits to decompose
    /// 
    /// # Returns
    /// Boolean shares representing the bits of the original value (LSB first)
    fn bit_decompose_from_opened(c: &F, b_bits: &[bool], bit_length: usize) -> Vec<bool>
    where
        F: PrimeField,
    {
        use crate::share::bit_ops::*;
        
        // Step 3: Compute e = c + 2^BitSize - p in extended ring
        let two_pow_l = num_bigint::BigUint::from(1u32) << bit_length;
        let two_pow_l_plus_1: num_bigint::BigUint = &two_pow_l << 1;
        let p = num_bigint::BigUint::from_bytes_le(&F::Params::MODULUS.to_bytes_le());
        let c_bigint = num_bigint::BigUint::from_bytes_le(&c.into_repr().to_bytes_le());
        
        // Compute e = (c + 2^l - p) mod 2^{l+1}
        let offset = if p < two_pow_l { 
            &two_pow_l - &p 
        } else { 
            two_pow_l_plus_1.clone() - (&p - &two_pow_l)
        };
        let e_bigint = (&c_bigint + offset) % &two_pow_l_plus_1;
        
        // Bit decompose e (public)
        let e_bytes = e_bigint.to_bytes_le();
        let mut e_bits = Vec::with_capacity(bit_length + 1);
        for i in 0..=bit_length {
            let byte_idx = i / 8;
            let bit_idx = i % 8;
            if byte_idx < e_bytes.len() {
                e_bits.push((e_bytes[byte_idx] >> bit_idx) & 1 == 1);
            } else {
                e_bits.push(false);
            }
        }
        
        // Step 4: Add [b_bits] + e_bits using bit adder
        let mut b_bits_with_overflow = b_bits.to_vec();
        b_bits_with_overflow.push(false); // Add dummy bit
        
        let mut beaver_source = DummyBooleanBeaverSource;
        let d_prime_bits = BooleanOps::bits_add_const(
            &b_bits_with_overflow,
            &e_bits,
            true,
            &mut beaver_source,
        );
        
        // Step 5: Extract overflow bit q = d'[bit_length]
        let q_share = d_prime_bits[bit_length];
        
        // Step 6: Compute NOT q
        let not_q_share = BooleanOps::bitwise_not(q_share);
        
        // Step 7: Prepare p_bits and multiply with NOT q
        let p_bits: Vec<bool> = (0..bit_length)
            .map(|i| {
                let byte_idx = i / 8;
                let bit_idx = i % 8;
                let p_bytes = p.to_bytes_le();
                if byte_idx < p_bytes.len() {
                    (p_bytes[byte_idx] >> bit_idx) & 1 == 1
                } else {
                    false
                }
            })
            .collect();
        
        // Multiply each public p bit with [NOT q] using secure AND
        let optional_p_bits: Vec<bool> = p_bits
            .iter()
            .map(|&p_bit| {
                use mpc_net::{MpcNet, MpcMultiNet as Net};
                // Convert public bit to secret share for secure multiplication
                let p_bit_share = if Net::am_king() { p_bit } else { false };
                // Perform secure AND between public bit (as share) and secret NOT q
                BooleanOps::beaver_bitwise_and(p_bit_share, not_q_share, &mut beaver_source)
            })
            .collect();
        
        // Step 8: Add d_prime_bits[0..bit_length] + optional_p_bits
        let d_prime_trimmed = d_prime_bits[..bit_length].to_vec();
        let h_bits = BooleanOps::bits_add(
            &d_prime_trimmed,
            &optional_p_bits,
            false,
            &mut beaver_source,
        );
        
        // Return the boolean bit shares
        h_bits
    }

}

pub type DensePolynomial<T> = Vec<T>;
pub type SparsePolynomial<T> = Vec<(usize, T)>;
pub type DenseOrSparsePolynomial<T> = Result<DensePolynomial<T>, SparsePolynomial<T>>;

pub trait ExtFieldShare<F: Field>:
    Clone + Copy + Debug + 'static + Send + Sync + PartialEq + Eq
{
    type Base: FieldShare<F::BasePrimeField>;
    type Ext: FieldShare<F>;
}
