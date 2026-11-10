//! BLS12-381 G1 subgroup membership check on OpenVM curve intrinsics.
//!
//! G1 has cofactor > 1, so a point that satisfies the curve equation may still lie
//! outside the prime-order subgroup. Ported from `openvm-curve-utils` in
//! axiom-crypto/openvm-reth-benchmark (#626); implements section 6 of
//! <https://eprint.iacr.org/2021/1130>.

use core::ops::{MulAssign, Neg};

use hex_literal::hex;
use openvm_ecc_guest::{weierstrass::WeierstrassPoint, Group};
use openvm_pairing::bls12_381::{Fp, G1Affine};

/// The BLS12-381 curve parameter `|u|`. `u` is negative; the sign is applied with `neg()`.
const X: [u64; 1] = [0xd201000000010000];

/// A non-trivial cube root of unity in Fq, used for the endomorphism `σ: (x, y) → (βx, y)`.
const BETA: Fp = Fp::from_const_bytes(hex!(
    "fefffeffffff012e02000a6213d817de8896f8e63ba9b3ddea770f6a07c669ba51ce76df2f67195f0000000000000000"
));

fn scalar_mul<const CHECK_SETUP: bool>(base: &G1Affine, scalar: &[u64]) -> G1Affine {
    let mut result = <G1Affine as WeierstrassPoint>::IDENTITY;
    let mut temp = base.clone();
    for limb in scalar {
        for bit_idx in 0..64u32 {
            if (limb >> bit_idx) & 1 == 1 {
                result.add_assign_impl::<CHECK_SETUP>(&temp);
            }
            temp.double_assign_impl::<CHECK_SETUP>();
        }
    }
    result
}

/// Returns `true` if `point` lies in the prime-order subgroup of G1.
///
/// The caller must already know `point` is on the curve; otherwise the result is meaningless.
pub fn g1_is_in_correct_subgroup(point: &G1Affine) -> bool {
    // `CHECK_SETUP = true` on the first call: G1Affine comes from `sw_declare`, which needs setup.
    let x_times_point = scalar_mul::<true>(point, &X);
    if point.eq(&x_times_point) && !point.is_identity() {
        return false;
    }
    let minus_x_squared_times_point = scalar_mul::<false>(&x_times_point, &X).neg();
    let mut endomorphism_point = point.clone();
    endomorphism_point.x_mut().mul_assign(&BETA);
    minus_x_squared_times_point.eq(&endomorphism_point)
}

#[cfg(test)]
mod tests {
    use super::*;
    use openvm_algebra_guest::IntMod;
    use openvm_ecc_guest::CyclicGroup;

    #[test]
    fn generator_in_subgroup() {
        assert!(g1_is_in_correct_subgroup(&G1Affine::GENERATOR));
    }

    #[test]
    fn identity_in_subgroup() {
        assert!(g1_is_in_correct_subgroup(
            &<G1Affine as WeierstrassPoint>::IDENTITY
        ));
    }

    #[test]
    fn rejects_on_curve_point_outside_subgroup() {
        let px = Fp::from_be_bytes(&[0u8; 48]).unwrap();
        let py = Fp::from_be_bytes(&hex!("1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaa9")).unwrap();
        // SAFETY: `from_xy` checks the curve equation and returns `None` if it fails.
        let point = unsafe { G1Affine::from_xy(px, py) }.expect("point should be on curve");
        assert!(!g1_is_in_correct_subgroup(&point));
    }
}
