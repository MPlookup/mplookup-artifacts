//! Polynomial operations for oblivious lookup.
//!
//! This module provides efficient polynomial operations used in the
//! secure oblivious lookup permutation algorithm:
//! - Multiplication using FFT
//! - Inversion
//! - Division with remainder
//! - Subproduct tree construction
//! - Multipoint evaluation

use ark_ff::{FftField, One, Zero};
use ark_poly::domain::radix2::Radix2EvaluationDomain;
use ark_poly::EvaluationDomain;

pub(crate) fn poly_mul<F: FftField>(p1: &[F], p2: &[F]) -> Vec<F> {
    if p1.is_empty() || p2.is_empty() {
        return vec![];
    }
    let deg = p1.len() + p2.len() - 2;
    let len_product = deg + 1;
    
    let domain_size = len_product.next_power_of_two();
    let domain =
        Radix2EvaluationDomain::<F>::new(domain_size).expect("Failed to create domain");
    let mut coeffs1 = p1.to_vec();
    coeffs1.resize(domain_size, F::zero());
    domain.fft_in_place(&mut coeffs1);
    let mut coeffs2 = p2.to_vec();
    coeffs2.resize(domain_size, F::zero());
    domain.fft_in_place(&mut coeffs2);
    for i in 0..domain_size {
        coeffs1[i] *= coeffs2[i];
    }
    domain.ifft_in_place(&mut coeffs1);
    coeffs1.truncate(len_product);
    assert_eq!(coeffs1.len(), len_product);
    coeffs1
}

pub(crate) fn poly_inv<F: FftField>(f: &[F], deg: usize) -> Vec<F> {
    if deg == 0 {
        return vec![];
    }
    let mut g = vec![f[0].inverse().unwrap()];
    let mut cur_deg = 1;
    while cur_deg < deg {
        let next_deg = cur_deg * 2;
        let ff: Vec<F> = f
            .iter()
            .cloned()
            .chain(std::iter::repeat(F::zero()))
            .take(next_deg)
            .collect();
        let mut tmp = poly_mul(&g, &ff);
        tmp.truncate(next_deg);
        tmp = poly_mul(&g, &tmp);
        tmp.truncate(next_deg);
        let mut new_g = vec![F::zero(); next_deg];
        for i in 0..cur_deg {
            new_g[i] = g[i] + g[i];
        }
        for i in 0..next_deg {
            new_g[i] -= tmp[i];
        }
        g = new_g;
        cur_deg = next_deg;
    }
    g.into_iter().take(deg).collect()
}

pub(crate) fn poly_divmod<F: FftField>(num: &[F], den: &[F]) -> (Vec<F>, Vec<F>) {
    let m = num.len().checked_sub(1).unwrap_or(0);
    let n = den.len().checked_sub(1).unwrap_or(0);
    if m < n {
        return (vec![F::zero()], num.to_vec());
    }
    let deg_q = m - n + 1;
    let mut rev_num = num.to_vec();
    rev_num.reverse();
    let mut rev_den = den.to_vec();
    rev_den.reverse();
    let inv_rev_den = poly_inv(&rev_den, deg_q);
    let mut q_rev = poly_mul(&rev_num[0..deg_q], &inv_rev_den);
    q_rev.truncate(deg_q);
    let mut q = q_rev.clone();
    q.reverse();
    let mut q_den = poly_mul(&q, den);
    let r_len = num.len();
    q_den.resize(r_len, F::zero());
    let r: Vec<F> = num.iter().zip(q_den.iter()).map(|(&a, &b)| a - b).collect();
    (q, r)
}

pub(crate) fn build_subproduct_tree<F: FftField>(points: &[F]) -> Vec<Vec<Vec<F>>> {
    let k = points.len();
    if k == 0 {
        return vec![vec![vec![F::one()]]];
    }
    let pow2 = k.next_power_of_two();
    let mut polys: Vec<Vec<F>> = points.iter().map(|&p| vec![-p, F::one()]).collect();
    polys.resize_with(pow2, || vec![F::one()]);
    let mut tree = vec![polys.clone()];
    let mut polys = polys;
    while polys.len() > 1 {
        let mut new_polys = vec![];
        for i in (0..polys.len()).step_by(2) {
            let p1 = polys[i].clone();
            let p2 = if i + 1 < polys.len() {
                polys[i + 1].clone()
            } else {
                vec![F::one()]
            };
            let product = poly_mul(&p1, &p2);
            new_polys.push(product);
        }
        tree.push(new_polys.clone());
        polys = new_polys;
    }
    tree
}

pub(crate) fn multipoint_eval<F: FftField>(p: &[F], tree: &Vec<Vec<Vec<F>>>) -> Vec<F> {
    let num_levels = tree.len();
    let leaf_level = 0;
    let root_level = num_levels - 1;
    let mut reduced_p: Vec<Vec<Vec<F>>> = (0..num_levels).map(|_| vec![]).collect();
    for lvl in 0..num_levels {
        reduced_p[lvl] = vec![vec![]; tree[lvl].len()];
    }
    reduced_p[root_level][0] = p.to_vec();
    for lvl in (0..root_level).rev() {
        for parent_nd in 0..tree[lvl + 1].len() {
            let curr_p = reduced_p[lvl + 1][parent_nd].clone();
            let left_nd = parent_nd * 2;
            let right_nd = parent_nd * 2 + 1;
            if left_nd < tree[lvl].len() {
                let left_sub = tree[lvl][left_nd].clone();
                let left_deg = left_sub.len().saturating_sub(1);
                if left_deg > 0 {
                    let (_, mut rem) = poly_divmod(&curr_p, &left_sub);
                    // Truncate the remainder to degree < deg(left_sub).  By the
                    // polynomial remainder theorem the coefficients at indices
                    // >= left_deg are exactly zero; keeping them would leave
                    // curr_p at size O(n/2) at every tree depth, inflating each
                    // subsequent poly_divmod from O(n/2^d) to O(n/2) and
                    // making the overall cost O(n^2 log n) instead of
                    // O(n log^2 n).  Truncating to left_sub.len()-1 is safe
                    // even when t is secret-shared: the length of every
                    // subproduct-tree node is a fixed public value determined
                    // by the level index and the public table size n, so no
                    // secret information is revealed.
                    rem.truncate(left_sub.len().saturating_sub(1));
                    reduced_p[lvl][left_nd] = rem;
                } else {
                    reduced_p[lvl][left_nd] = curr_p.clone();
                }
            }
            if right_nd < tree[lvl].len() {
                let right_sub = tree[lvl][right_nd].clone();
                let right_deg = right_sub.len().saturating_sub(1);
                if right_deg > 0 {
                    let (_, mut rem) = poly_divmod(&curr_p, &right_sub);
                    rem.truncate(right_sub.len().saturating_sub(1)); // same truncation
                    reduced_p[lvl][right_nd] = rem;
                } else {
                    reduced_p[lvl][right_nd] = curr_p.clone();
                }
            }
        }
    }
    let mut evals = vec![];
    for nd in 0..tree[leaf_level].len() {
        if tree[leaf_level][nd].len() > 1 {
            let red = &reduced_p[leaf_level][nd];
            if red.is_empty() {
                evals.push(F::zero());
            } else {
                evals.push(red[0]);
            }
        }
    }
    evals
}
