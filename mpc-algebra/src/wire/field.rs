use derivative::Derivative;
use log::debug;
use rand::Rng;
use zeroize::Zeroize;

use ark_ff::bytes::{FromBytes, ToBytes};
use ark_ff::prelude::*;
use ark_ff::{poly_stub, FftField};
use ark_serialize::{
    CanonicalDeserialize, CanonicalDeserializeWithFlags, CanonicalSerialize,
    CanonicalSerializeWithFlags, Flags, SerializationError,
};
use mpc_trait::MpcWire;

use std::fmt::{self, Debug, Display, Formatter};
use std::io::{self, Read, Write};
use std::iter::{Product, Sum};
use std::marker::PhantomData;
use std::ops::*;

use super::super::share::field::FieldShare;
use super::super::share::BeaverSource;
use crate::Reveal;
use mpc_net::{MpcNet, MpcMultiNet as Net};

#[derive(Clone, Copy, Hash, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum MpcField<F: Field, S: FieldShare<F>> {
    Public(F),
    Shared(S),
}

impl_basics_2!(FieldShare, Field, MpcField);

#[derive(Derivative)]
#[derivative(Default(bound = ""), Clone(bound = ""), Copy(bound = ""))]
pub struct DummyFieldTripleSource<T, S> {
    _scalar: PhantomData<T>,
    _share: PhantomData<S>,
}

impl<T: Field, S: FieldShare<T>> BeaverSource<S, S, S> for DummyFieldTripleSource<T, S> {
    #[inline]
    fn triple(&mut self) -> (S, S, S) {
        (
            S::from_add_shared(if Net::am_king() {
                T::one()
            } else {
                T::zero()
            }),
            S::from_add_shared(if Net::am_king() {
                T::one()
            } else {
                T::zero()
            }),
            S::from_add_shared(if Net::am_king() {
                T::one()
            } else {
                T::zero()
            }),
        )
    }
    #[inline]
    fn inv_pair(&mut self) -> (S, S) {
        (
            S::from_add_shared(if Net::am_king() {
                T::one()
            } else {
                T::zero()
            }),
            S::from_add_shared(if Net::am_king() {
                T::one()
            } else {
                T::zero()
            }),
        )
    }
}

impl<T: Field, S: FieldShare<T>> MpcField<T, S> {
    #[inline]
    pub fn inv(self) -> Option<Self> {
        match self {
            Self::Public(x) => x.inverse().map(MpcField::Public),
            Self::Shared(x) => Some(MpcField::Shared(
                x.inv(&mut DummyFieldTripleSource::default()),
            )),
        }
    }
    pub fn all_public_or_shared(v: impl IntoIterator<Item = Self>) -> Result<Vec<T>, Vec<S>> {
        let mut out_a = Vec::new();
        let mut out_b = Vec::new();
        for s in v {
            match s {
                Self::Public(x) => out_a.push(x),
                Self::Shared(x) => out_b.push(x),
            }
        }
        if out_a.len() > 0 && out_b.len() > 0 {
            panic!("Heterogeous")
        } else if out_a.len() > 0 {
            Ok(out_a)
        } else {
            Err(out_b)
        }
    }

    /// Secure comparison that works for both public and shared values.
    /// 
    /// Returns a field element representing the less-than bit:
    /// - Returns 0 if self >= other
    /// - Returns 1 if self < other
    /// 
    /// For **public values**, performs normal plaintext comparison and returns the result as a public field element.
    /// 
    /// For **shared values**, uses bit decomposition protocol to securely compare without revealing values.
    /// The result is returned as a **secret share** maintaining privacy.
    /// 
    /// For **heterogeneous values** (one Public, one Shared), automatically converts the Public value
    /// to Shared using `from_public()` before comparison, returning a Shared result.
    /// 
    /// # MPC-Safe Comparison Protocol
    /// 
    /// Implements secure comparison using:
    /// - edaBits for bit decomposition
    /// - Boolean secret sharing for bit-level operations
    /// - Bit adder circuits with carry propagation
    /// - Comparison circuit to extract less-than bit
    /// 
    /// ## Protocol Steps
    /// 
    /// ```text
    /// 1. Bit-decompose [a] and [b] using edaBits protocol
    /// 2. For each bit position i from MSB to LSB:
    ///    - Compute [xor_i] = [a_i] XOR [b_i]
    /// 3. Find first position where bits differ (using prefix-OR)
    /// 4. At that position, [b_i] gives the less-than bit
    /// 5. Return as secret-shared field element (no reveal)
    /// ```
    /// 
    /// See `FieldShare::secure_cmp()` documentation for detailed protocol description.
    /// 
    /// # Examples
    /// 
    /// ```no_run
    /// # use mpc_algebra::honest_but_curious::MpcField;
    /// # use ark_bls12_377::Fr;
    /// // Public values: returns public 0 or 1
    /// let a = MpcField::<Fr>::Public(Fr::from(10u64));
    /// let b = MpcField::<Fr>::Public(Fr::from(42u64));
    /// let is_less = a.secure_cmp(&b);  // Returns Public(Fr::from(1)) since 10 < 42
    /// 
    /// // Shared values: returns secret-shared 0 or 1
    /// # /*
    /// let x = MpcField::Shared(share1);
    /// let y = MpcField::Shared(share2);
    /// let is_less = x.secure_cmp(&y);  // Returns Shared([less_than_bit])
    /// // Can be used in further MPC computations without revealing the result
    /// # */
    /// ```
    #[inline]
    pub fn secure_cmp(&self, other: &Self) -> Self 
    where
        T: PrimeField,
        S: FieldShare<T>,
    {
        match (self, other) {
            (Self::Public(x), Self::Public(y)) => {
                // For public values, return 1 if x < y, else 0
                if x < y {
                    Self::Public(T::one())
                } else {
                    Self::Public(T::zero())
                }
            },
            (Self::Shared(x), Self::Shared(y)) => Self::Shared(x.secure_cmp(y)),
            (Self::Public(x), Self::Shared(y)) => {
                // Convert public value to shared and compare
                let x_shared = S::from_public(*x);
                Self::Shared(x_shared.secure_cmp(y))
            },
            (Self::Shared(x), Self::Public(y)) => {
                // Convert public value to shared and compare
                let y_shared = S::from_public(*y);
                Self::Shared(x.secure_cmp(&y_shared))
            }
        }
    }

    /// Batch secure comparison for multiple pairs.
    /// 
    /// Compares multiple pairs of values in a single batched operation.
    /// For each pair (left[i], right[i]), returns 1 if left[i] < right[i], else 0.
    /// 
    /// This function optimizes communication by batching the reveal operations
    /// during bit decomposition, reducing the number of communication rounds.
    /// 
    /// # Arguments
    /// * `pairs` - A slice of (left, right) tuples to compare
    /// 
    /// # Returns
    /// A vector of comparison results, where result[i] is 1 if left[i] < right[i], else 0
    /// 
    /// # Optimization
    /// This function achieves performance gains by:
    /// - For all-public values: uses fast path with no MPC operations
    /// - For shared values: batches all reveal operations in bit decomposition into a single round
    /// - Reduces n sequential reveals to 1 batched reveal for n comparisons
    /// - Comparison operations use batched boolean operations as in single secure_cmp
    #[inline]
    pub fn batch_secure_cmp(pairs: &[(Self, Self)]) -> Vec<Self>
    where
        T: PrimeField,
        S: FieldShare<T>,
    {
        if pairs.is_empty() {
            return vec![];
        }

        // Check if all values are public - fast path
        let all_public = pairs.iter().all(|(left, right)| {
            matches!(left, Self::Public(_)) && matches!(right, Self::Public(_))
        });

        if all_public {
            // Fast path: all public values
            return pairs.iter().map(|(left, right)| {
                if let (Self::Public(x), Self::Public(y)) = (left, right) {
                    if x < y {
                        Self::Public(T::one())
                    } else {
                        Self::Public(T::zero())
                    }
                } else {
                    unreachable!()
                }
            }).collect();
        }

        // Convert all values to shared and perform batched comparison
        let mut left_shares = Vec::with_capacity(pairs.len());
        let mut right_shares = Vec::with_capacity(pairs.len());

        for (left, right) in pairs {
            let left_share = match left {
                Self::Public(x) => S::from_public(*x),
                Self::Shared(x) => *x,
            };
            let right_share = match right {
                Self::Public(y) => S::from_public(*y),
                Self::Shared(y) => *y,
            };
            left_shares.push(left_share);
            right_shares.push(right_share);
        }

        // Batch compare all pairs
        let result_shares = S::batch_secure_cmp(&left_shares, &right_shares);

        // Convert back to MpcField
        result_shares.into_iter().map(Self::Shared).collect()
    }

    /// Batch bit decompose N values into boolean shares.
    ///
    /// Batches all reveals into a single communication round.
    /// Returns `bits[i][b]` = boolean share of bit `b` of `values[i]` (LSB first).
    /// All modulus bits are returned; callers that only need the low `n` bits may
    /// truncate each inner Vec.
    #[inline]
    pub fn batch_bit_decompose_bool(values: &[Self]) -> Vec<Vec<bool>>
    where
        T: PrimeField,
        S: FieldShare<T>,
    {
        if values.is_empty() {
            return vec![];
        }
        // Convert all to shares
        let shares: Vec<S> = values
            .iter()
            .map(|v| match v {
                Self::Public(x) => S::from_public(*x),
                Self::Shared(s) => *s,
            })
            .collect();
        S::batch_bit_decompose_elems(&shares)
    }

    /// Batch multiply two vectors of MpcField element-wise.
    ///
    /// All multiplications are batched into a single round of communication.
    #[inline]
    pub fn batch_mul_vec(xs: Vec<Self>, ys: Vec<Self>) -> Vec<Self>
    where
        T: PrimeField,
        S: FieldShare<T>,
    {
        assert_eq!(xs.len(), ys.len());
        if xs.is_empty() {
            return vec![];
        }
        // Check homogeneity
        let all_public = xs.iter().all(|x| matches!(x, Self::Public(_)))
            && ys.iter().all(|y| matches!(y, Self::Public(_)));
        if all_public {
            return xs
                .iter()
                .zip(ys.iter())
                .map(|(x, y)| match (x, y) {
                    (Self::Public(a), Self::Public(b)) => {
                        let mut r = *a;
                        r *= b;
                        Self::Public(r)
                    }
                    _ => unreachable!(),
                })
                .collect();
        }
        let x_shares: Vec<S> = xs
            .iter()
            .map(|v| match v {
                Self::Public(x) => S::from_public(*x),
                Self::Shared(s) => *s,
            })
            .collect();
        let y_shares: Vec<S> = ys
            .iter()
            .map(|v| match v {
                Self::Public(y) => S::from_public(*y),
                Self::Shared(s) => *s,
            })
            .collect();
        let result_shares =
            S::batch_mul(x_shares, y_shares, &mut DummyFieldTripleSource::default());
        result_shares.into_iter().map(Self::Shared).collect()
    }

    /// Batch compare N pairs of precomputed bit sequences.
    ///
    /// Returns `result[i]` = boolean share of (left[i] < right[i]).
    /// All comparisons are processed in O(bit_count) sequential rounds,
    /// regardless of N (N pairs are batched at each round).
    #[inline]
    pub fn batch_compare_bits(
        left_bits: &[Vec<bool>],
        right_bits: &[Vec<bool>],
    ) -> Vec<bool>
    where
        T: PrimeField,
        S: FieldShare<T>,
    {
        S::batch_compare_bits_shared(left_bits, right_bits)
    }

    /// Batch B2A conversion: convert N boolean shares to N arithmetic shares.
    ///
    /// Uses a single `batch_open_bool` call (1 broadcast) for all N conversions.
    #[inline]
    pub fn batch_b2a_bool(bool_shares: &[bool]) -> Vec<Self>
    where
        T: PrimeField,
        S: FieldShare<T>,
    {
        S::batch_b2a_elems(bool_shares)
            .into_iter()
            .map(Self::Shared)
            .collect()
    }
}

impl<T: PrimeField, S: FieldShare<T>> MpcField<T, S> {
    /// Secure inequality check that works for both public and shared values.
    /// 
    /// Returns a field element indicating whether values are not equal:
    /// - Returns 1 if values are NOT equal
    /// - Returns 0 if values ARE equal
    /// 
    /// For **public values**, performs normal inequality comparison and returns the result as a public field element.
    /// 
    /// For **shared values**, uses Fermat's little theorem to compute a shared result without revealing the values.
    /// 
    /// For **heterogeneous values** (one Public, one Shared), automatically converts the Public value
    /// to Shared using `from_public()` before comparison, returning a Shared result.
    #[inline]
    pub fn secure_neq(&self, other: &Self) -> Self {
        match (self, other) {
            (Self::Public(x), Self::Public(y)) => {
                // Return 0 if equal, 1 if not equal
                if x == y {
                    Self::Public(T::zero())
                } else {
                    Self::Public(T::one())
                }
            },
            (Self::Shared(x), Self::Shared(y)) => {
                Self::Shared(x.secure_neq(y, &mut DummyFieldTripleSource::default()))
            },
            (Self::Public(x), Self::Shared(y)) => {
                // Convert public value to shared and compare
                let x_shared = S::from_public(*x);
                Self::Shared(x_shared.secure_neq(y, &mut DummyFieldTripleSource::default()))
            },
            (Self::Shared(x), Self::Public(y)) => {
                // Convert public value to shared and compare
                let y_shared = S::from_public(*y);
                Self::Shared(x.secure_neq(&y_shared, &mut DummyFieldTripleSource::default()))
            }
        }
    }

    /// Secure equality check that works for both public and shared values.
    /// Returns 1 if equal, 0 if not equal (inverse of secure_neq).
    /// 
    /// For public values, performs normal equality comparison and returns Public(1) for equal, Public(0) for not equal.
    /// For shared values, computes 1 - secure_neq to get the equality indicator.
    /// 
    /// # Examples
    /// 
    /// ```no_run
    /// # use mpc_algebra::honest_but_curious::MpcField;
    /// # use ark_bls12_377::Fr;
    /// // Public values: returns 1 if equal, 0 if not equal
    /// let a = MpcField::<Fr>::Public(Fr::from(10u64));
    /// let b = MpcField::<Fr>::Public(Fr::from(10u64));
    /// let c = MpcField::<Fr>::Public(Fr::from(42u64));
    /// 
    /// let eq_result = a.secure_eq(&b);  // Returns Public(Fr::from(1))
    /// let neq_result = a.secure_eq(&c); // Returns Public(Fr::from(0))
    /// ```
    #[inline]
    pub fn secure_eq(&self, other: &Self) -> Self {
        let neq = self.secure_neq(other);
        match neq {
            Self::Public(x) => Self::Public(T::one() - x),
            Self::Shared(x) => {
                let mut one_share = S::from_public(T::one());
                one_share.sub(&x);
                Self::Shared(one_share)
            }
        }
    }
}
impl<'a, T: Field, S: FieldShare<T>> MulAssign<&'a MpcField<T, S>> for MpcField<T, S> {
    #[inline]
    fn mul_assign(&mut self, other: &Self) {
        match self {
            // for some reason, a two-stage match (rather than a tuple match) avoids moving
            // self
            MpcField::Public(x) => match other {
                MpcField::Public(y) => {
                    *x *= y;
                }
                MpcField::Shared(y) => {
                    let mut t = *y;
                    t.scale(x);
                    *self = MpcField::Shared(t);
                }
            },
            MpcField::Shared(x) => match other {
                MpcField::Public(y) => {
                    x.scale(y);
                }
                MpcField::Shared(y) => {
                    let t = x.mul(*y, &mut DummyFieldTripleSource::default());
                    *self = MpcField::Shared(t);
                }
            },
        }
    }
}
impl<T: Field, S: FieldShare<T>> One for MpcField<T, S> {
    #[inline]
    fn one() -> Self {
        MpcField::Public(T::one())
    }
}
impl<T: Field, S: FieldShare<T>> Product for MpcField<T, S> {
    #[inline]
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::one(), Mul::mul)
    }
}
impl<'a, T: Field, S: FieldShare<T> + 'a> Product<&'a MpcField<T, S>> for MpcField<T, S> {
    #[inline]
    fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::one(), |x, y| x.mul(y.clone()))
    }
}

impl<'a, T: Field, S: FieldShare<T>> DivAssign<&'a MpcField<T, S>> for MpcField<T, S> {
    #[inline]
    fn div_assign(&mut self, other: &Self) {
        match self {
            // for some reason, a two-stage match (rather than a tuple match) avoids moving
            // self
            MpcField::Public(x) => match other {
                MpcField::Public(y) => {
                    *x /= y;
                }
                MpcField::Shared(y) => {
                    let mut t = y.inv(&mut DummyFieldTripleSource::default());
                    t.scale(&x);
                    *self = MpcField::Shared(t);
                }
            },
            MpcField::Shared(x) => match other {
                MpcField::Public(y) => {
                    x.scale(&y.inverse().unwrap());
                }
                MpcField::Shared(y) => {
                    let src = &mut DummyFieldTripleSource::default();
                    *x = x.div(*y, src);
                }
            },
        }
    }
}

impl_ref_ops!(
    Mul,
    MulAssign,
    mul,
    mul_assign,
    Field,
    FieldShare,
    MpcField
);
impl_ref_ops!(
    Add,
    AddAssign,
    add,
    add_assign,
    Field,
    FieldShare,
    MpcField
);
impl_ref_ops!(
    Div,
    DivAssign,
    div,
    div_assign,
    Field,
    FieldShare,
    MpcField
);
impl_ref_ops!(
    Sub,
    SubAssign,
    sub,
    sub_assign,
    Field,
    FieldShare,
    MpcField
);

impl<T: Field, S: FieldShare<T>> MpcWire for MpcField<T, S> {
    #[inline]
    fn publicize(&mut self) {
        match self {
            MpcField::Shared(s) => {
                *self = MpcField::Public(s.open());
            }
            _ => {}
        }
        debug_assert!({
            let self_val = if let MpcField::Public(s) = self {
                s.clone()
            } else {
                unreachable!()
            };
            super::macros::check_eq(self_val.clone());
            true
        })
    }
    #[inline]
    fn is_shared(&self) -> bool {
        match self {
            MpcField::Shared(_) => true,
            MpcField::Public(_) => false,
        }
    }
}

impl<T: Field, S: FieldShare<T>> Reveal for MpcField<T, S> {
    type Base = T;
    #[inline]
    fn reveal(self) -> Self::Base {
        let result = match self {
            Self::Shared(s) => s.reveal(),
            Self::Public(s) => s,
        };
        super::macros::check_eq(result.clone());
        result
    }
    #[inline]
    fn from_public(b: Self::Base) -> Self {
        MpcField::Public(b)
    }
    #[inline]
    fn from_add_shared(b: Self::Base) -> Self {
        MpcField::Shared(S::from_add_shared(b))
    }
    #[inline]
    fn unwrap_as_public(self) -> Self::Base {
        match self {
            Self::Shared(s) => s.unwrap_as_public(),
            Self::Public(s) => s,
        }
    }
    #[inline]
    fn king_share<R: Rng>(f: Self::Base, rng: &mut R) -> Self {
        Self::Shared(S::king_share(f, rng))
    }
    #[inline]
    fn king_share_batch<R: Rng>(f: Vec<Self::Base>, rng: &mut R) -> Vec<Self> {
        S::king_share_batch(f, rng).into_iter().map(Self::Shared).collect()
    }
    fn init_protocol() {
        S::init_protocol()
    }
    fn deinit_protocol() {
        S::deinit_protocol()
    }
}

from_prim!(bool, Field, FieldShare, MpcField);
from_prim!(u8, Field, FieldShare, MpcField);
from_prim!(u16, Field, FieldShare, MpcField);
from_prim!(u32, Field, FieldShare, MpcField);
from_prim!(u64, Field, FieldShare, MpcField);
from_prim!(u128, Field, FieldShare, MpcField);

impl<T: PrimeField, S: FieldShare<T>> std::str::FromStr for MpcField<T, S> {
    type Err = T::Err;
    #[inline]
    fn from_str(s: &str) -> Result<Self, T::Err> {
        T::from_str(s).map(Self::Public)
    }
}

impl<F: PrimeField, S: FieldShare<F>> Field for MpcField<F, S> {
    type BasePrimeField = Self;
    #[inline]
    fn extension_degree() -> u64 {
        unimplemented!("extension_degree")
    }
    #[inline]
    fn from_base_prime_field_elems(_b: &[<Self as ark_ff::Field>::BasePrimeField]) -> Option<Self> {
        unimplemented!()
        // assert!(b.len() > 0);
        // let shared = b[0].is_shared();
        // assert!(b.iter().all(|e| e.is_shared() == shared));
        // let base_values = b.iter().map(|e| e.unwrap_as_public()).collect::<Vec<_>>();
        // F::from_base_prime_field_elems(&base_values).map(|val| Self::new(val, shared))
    }
    #[inline]
    fn double(&self) -> Self {
        Self::Public(F::from(2u8)) * self
    }
    #[inline]
    fn double_in_place(&mut self) -> &mut Self {
        *self *= Self::Public(F::from(2u8));
        self
    }
    #[inline]
    fn from_random_bytes_with_flags<Fl: Flags>(b: &[u8]) -> Option<(Self, Fl)> {
        F::from_random_bytes_with_flags(b).map(|(val, f)| (Self::Shared(S::from_public(val)), f))
    }
    #[inline]
    fn square(&self) -> Self {
        self.clone() * self
    }
    #[inline]
    fn square_in_place(&mut self) -> &mut Self {
        *self *= self.clone();
        self
    }
    #[inline]
    fn inverse(&self) -> Option<Self> {
        self.inv()
    }
    #[inline]
    fn inverse_in_place(&mut self) -> Option<&mut Self> {
        self.inv().map(|i| {
            *self = i;
            self
        })
    }
    #[inline]
    fn frobenius_map(&mut self, _: usize) {
        unimplemented!("frobenius_map")
    }

    fn batch_product_in_place(selfs: &mut [Self], others: &[Self]) {
        let selfs_shared = selfs[0].is_shared();
        let others_shared = others[0].is_shared();
        assert!(
            selfs.iter().all(|s| s.is_shared() == selfs_shared),
            "Selfs heterogenously shared!"
        );
        assert!(
            others.iter().all(|s| s.is_shared() == others_shared),
            "others heterogenously shared!"
        );
        if selfs_shared && others_shared {
            let sshares = selfs
                .iter()
                .map(|s| match s {
                    Self::Shared(s) => s.clone(),
                    Self::Public(_) => unreachable!(),
                })
                .collect();
            let oshares = others
                .iter()
                .map(|s| match s {
                    Self::Shared(s) => s.clone(),
                    Self::Public(_) => unreachable!(),
                })
                .collect();
            let nshares = S::batch_mul(sshares, oshares, &mut DummyFieldTripleSource::default());
            for (self_, new) in selfs.iter_mut().zip(nshares.into_iter()) {
                *self_ = Self::Shared(new);
            }
        } else {
            for (a, b) in ark_std::cfg_iter_mut!(selfs).zip(others.iter()) {
                *a *= b;
            }
        }
    }
    fn batch_division_in_place(selfs: &mut [Self], others: &[Self]) {
        let selfs_shared = selfs[0].is_shared();
        let others_shared = others[0].is_shared();
        assert!(
            selfs.iter().all(|s| s.is_shared() == selfs_shared),
            "Selfs heterogenously shared!"
        );
        assert!(
            others.iter().all(|s| s.is_shared() == others_shared),
            "others heterogenously shared!"
        );
        if selfs_shared && others_shared {
            let sshares = selfs
                .iter()
                .map(|s| match s {
                    Self::Shared(s) => s.clone(),
                    Self::Public(_) => unreachable!(),
                })
                .collect();
            let oshares = others
                .iter()
                .map(|s| match s {
                    Self::Shared(s) => s.clone(),
                    Self::Public(_) => unreachable!(),
                })
                .collect();
            let nshares = S::batch_div(sshares, oshares, &mut DummyFieldTripleSource::default());
            for (self_, new) in selfs.iter_mut().zip(nshares.into_iter()) {
                *self_ = Self::Shared(new);
            }
        } else {
            for (a, b) in ark_std::cfg_iter_mut!(selfs).zip(others.iter()) {
                *a *= b;
            }
        }
    }
    fn partial_products_in_place(selfs: &mut [Self]) {
        let selfs_shared = selfs[0].is_shared();
        assert!(
            selfs.iter().all(|s| s.is_shared() == selfs_shared),
            "Selfs heterogenously shared!"
        );
        if selfs_shared {
            let sshares = selfs
                .iter()
                .map(|s| match s {
                    Self::Shared(s) => s.clone(),
                    Self::Public(_) => unreachable!(),
                })
                .collect();
            for (self_, new) in selfs.iter_mut().zip(
                S::partial_products(sshares, &mut DummyFieldTripleSource::default()).into_iter(),
            ) {
                *self_ = Self::Shared(new);
            }
        } else {
            for i in 1..selfs.len() {
                let last = selfs[i - 1];
                selfs[i] *= &last;
            }
        }
    }
    fn has_univariate_div_qr() -> bool {
        true
    }
    fn univariate_div_qr<'a>(
        num: poly_stub::DenseOrSparsePolynomial<Self>,
        den: poly_stub::DenseOrSparsePolynomial<Self>,
    ) -> Option<(
        poly_stub::DensePolynomial<Self>,
        poly_stub::DensePolynomial<Self>,
    )> {
        use poly_stub::DenseOrSparsePolynomial::*;
        let shared_num = match num {
            DPolynomial(d) => Ok(d.into_owned().coeffs.into_iter().map(|c| match c {
                MpcField::Shared(s) => s,
                MpcField::Public(_) => panic!("public numerator"),
            }).collect()),
            SPolynomial(d) => Err(d.into_owned().coeffs.into_iter().map(|(i, c)| match c {
                MpcField::Shared(s) => (i, s),
                MpcField::Public(_) => panic!("public numerator"),
            }).collect()),
        };
        let pub_denom = match den {
            DPolynomial(d) => Ok(d.into_owned().coeffs.into_iter().map(|c| match c {
                MpcField::Public(s) => s,
                MpcField::Shared(_) => panic!("shared denominator"),
            }).collect()),
            SPolynomial(d) => Err(d.into_owned().coeffs.into_iter().map(|(i, c)| match c {
                MpcField::Public(s) => (i, s),
                MpcField::Shared(_) => panic!("shared denominator"),
            }).collect()),
        };
        S::univariate_div_qr(shared_num, pub_denom).map(|(q, r)| {
            (
                poly_stub::DensePolynomial {
                    coeffs: q.into_iter().map(|qc| MpcField::Shared(qc)).collect(),
                },
                poly_stub::DensePolynomial {
                    coeffs: r.into_iter().map(|rc| MpcField::Shared(rc)).collect(),
                },
            )
        })
    }
}

impl<F: PrimeField, S: FieldShare<F>> FftField for MpcField<F, S> {
    type FftParams = F::FftParams;
    #[inline]
    fn two_adic_root_of_unity() -> Self {
        Self::from_public(F::two_adic_root_of_unity())
    }
    #[inline]
    fn large_subgroup_root_of_unity() -> Option<Self> {
        F::large_subgroup_root_of_unity().map(Self::from_public)
    }
    #[inline]
    fn multiplicative_generator() -> Self {
        Self::from_public(F::multiplicative_generator())
    }
}

impl<F: PrimeField, S: FieldShare<F>> PrimeField for MpcField<F, S> {
    type Params = F::Params;
    type BigInt = F::BigInt;
    #[inline]
    fn from_repr(_r: <Self as PrimeField>::BigInt) -> Option<Self> {
        unimplemented!("No BigInt reprs for shared fields! (from_repr)")
        //F::from_repr(r).map(|v| Self::from_public(v))
    }
    // We're assuming that into_repr is linear
    #[inline]
    fn into_repr(&self) -> <Self as PrimeField>::BigInt {
        unimplemented!("No BigInt reprs for shared fields! (into_repr)")
        //self.unwrap_as_public().into_repr()
    }
}

impl<F: PrimeField, S: FieldShare<F>> SquareRootField for MpcField<F, S> {
    #[inline]
    fn legendre(&self) -> ark_ff::LegendreSymbol {
        todo!()
    }
    #[inline]
    fn sqrt(&self) -> Option<Self> {
        todo!()
    }
    #[inline]
    fn sqrt_in_place(&mut self) -> Option<&mut Self> {
        todo!()
    }
}

mod poly_impl {

    use crate::share::*;
    use crate::wire::*;
    use crate::Reveal;
    use ark_ff::PrimeField;
    use ark_poly::domain::{EvaluationDomain, GeneralEvaluationDomain};
    use ark_poly::evaluations::univariate::Evaluations;
    use ark_poly::univariate::DensePolynomial;

    impl<E: PrimeField, S: FieldShare<E>> Reveal for DensePolynomial<MpcField<E, S>> {
        type Base = DensePolynomial<E>;
        struct_reveal_simp_impl!(DensePolynomial; coeffs);
    }

    impl<F: PrimeField, S: FieldShare<F>> Reveal for Evaluations<MpcField<F, S>> {
        type Base = Evaluations<F>;

        fn reveal(self) -> Self::Base {
            Evaluations::from_vec_and_domain(
                self.evals.reveal(),
                GeneralEvaluationDomain::new(self.domain.size()).unwrap(),
            )
        }

        fn from_add_shared(b: Self::Base) -> Self {
            Evaluations::from_vec_and_domain(
                Reveal::from_add_shared(b.evals),
                GeneralEvaluationDomain::new(b.domain.size()).unwrap(),
            )
        }

        fn from_public(b: Self::Base) -> Self {
            Evaluations::from_vec_and_domain(
                Reveal::from_public(b.evals),
                GeneralEvaluationDomain::new(b.domain.size()).unwrap(),
            )
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_377::Fr;
    use crate::honest_but_curious::MpcField as HbcMpcField;

    #[test]
    fn test_secure_neq_public_values() {
        // Test equality with public values should work
        // Returns Public(0) for equal, Public(1) for not equal
        let a = HbcMpcField::<Fr>::Public(Fr::from(42u64));
        let b = HbcMpcField::<Fr>::Public(Fr::from(42u64));
        let c = HbcMpcField::<Fr>::Public(Fr::from(10u64));

        let eq_result = a.secure_neq(&b);
        let neq_result = a.secure_neq(&c);
        
        // Check that equal values return 0
        assert!(matches!(eq_result, HbcMpcField::Public(x) if x == Fr::zero()), 
                "Same public values should return 0");
        
        // Check that not-equal values return 1
        assert!(matches!(neq_result, HbcMpcField::Public(x) if x == Fr::one()), 
                "Different public values should return 1");
    }

    #[test]
    fn test_secure_eq_public_values() {
        // Test equality with public values should work
        // Returns Public(1) for equal, Public(0) for not equal (inverse of secure_neq)
        let a = HbcMpcField::<Fr>::Public(Fr::from(42u64));
        let b = HbcMpcField::<Fr>::Public(Fr::from(42u64));
        let c = HbcMpcField::<Fr>::Public(Fr::from(10u64));

        let eq_result = a.secure_eq(&b);
        let neq_result = a.secure_eq(&c);
        
        // Check that equal values return 1
        assert!(matches!(eq_result, HbcMpcField::Public(x) if x == Fr::one()), 
                "Same public values should return 1");
        
        // Check that not-equal values return 0
        assert!(matches!(neq_result, HbcMpcField::Public(x) if x == Fr::zero()), 
                "Different public values should return 0");
    }

    #[test]
    fn test_secure_cmp_public_values() {
        // Test comparison with public values should work
        // Returns 1 if a < b, else 0
        let a = HbcMpcField::<Fr>::Public(Fr::from(10u64));
        let b = HbcMpcField::<Fr>::Public(Fr::from(42u64));
        let c = HbcMpcField::<Fr>::Public(Fr::from(10u64));

        // 10 < 42 should return 1
        let result_ab = a.secure_cmp(&b);
        assert!(matches!(result_ab, HbcMpcField::Public(x) if x == Fr::one()), "10 < 42 should return 1");
        
        // 42 < 10 should return 0
        let result_ba = b.secure_cmp(&a);
        assert!(matches!(result_ba, HbcMpcField::Public(x) if x == Fr::zero()), "42 >= 10 should return 0");
        
        // 10 < 10 should return 0
        let result_ac = a.secure_cmp(&c);
        assert!(matches!(result_ac, HbcMpcField::Public(x) if x == Fr::zero()), "10 >= 10 should return 0");
    }

    #[test]
    #[ignore] // Requires MPC network setup
    fn test_secure_neq_shared_values() {
        use crate::share::add::AdditiveFieldShare;
        
        // Test that comparing shared values doesn't panic
        // Note: This test requires proper MPC network setup to run
        let a = HbcMpcField::<Fr>::Shared(AdditiveFieldShare::from_add_shared(Fr::from(42u64)));
        let b = HbcMpcField::<Fr>::Shared(AdditiveFieldShare::from_add_shared(Fr::from(42u64)));
        
        // This should return a Shared result without panicking
        let result = a.secure_neq(&b);
        
        // Verify it returns a Shared value
        assert!(matches!(result, HbcMpcField::Shared(_)), 
                "secure_neq on shared values should return Shared");
    }

    #[test]
    #[ignore] // Requires MPC network setup
    fn test_secure_cmp_shared_values() {
        use crate::share::add::AdditiveFieldShare;
        
        // Test that comparing shared values works
        // Note: This test requires proper MPC network setup to run
        let a = HbcMpcField::<Fr>::Shared(AdditiveFieldShare::from_add_shared(Fr::from(10u64)));
        let b = HbcMpcField::<Fr>::Shared(AdditiveFieldShare::from_add_shared(Fr::from(42u64)));
        
        // This should work now with the implementation
        let result = a.secure_cmp(&b);
        
        // In a real MPC setting, this would compute the comparison securely
        // The result should be a Shared value (0 or 1)
        assert!(matches!(result, HbcMpcField::Shared(_)),
                "secure_cmp on shared values should return Shared");
    }

    #[test]
    #[ignore] // Requires MPC network setup
    fn test_secure_neq_heterogeneous() {
        use crate::share::add::AdditiveFieldShare;
        
        // Test that comparing public with shared works by converting to shared
        let a_pub = HbcMpcField::<Fr>::Public(Fr::from(42u64));
        let b_shared = HbcMpcField::<Fr>::Shared(AdditiveFieldShare::from_add_shared(Fr::from(42u64)));
        
        // This should work now by converting public to shared
        let result = a_pub.secure_neq(&b_shared);
        
        // Verify it returns a Shared value
        assert!(matches!(result, HbcMpcField::Shared(_)), 
                "secure_neq on heterogeneous values (Public, Shared) should return Shared");
    }

    #[test]
    #[ignore] // Requires MPC network setup
    fn test_secure_cmp_heterogeneous() {
        use crate::share::add::AdditiveFieldShare;
        
        // Test that comparing public with shared works by converting to shared
        let a_pub = HbcMpcField::<Fr>::Public(Fr::from(10u64));
        let b_shared = HbcMpcField::<Fr>::Shared(AdditiveFieldShare::from_add_shared(Fr::from(42u64)));
        
        // This should work now by converting public to shared
        let result = a_pub.secure_cmp(&b_shared);
        
        // Verify it returns a Shared value
        assert!(matches!(result, HbcMpcField::Shared(_)),
                "secure_cmp on heterogeneous values (Public, Shared) should return Shared");
    }
}
