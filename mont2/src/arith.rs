// Copyright 2021 Eric Schorn; Licensed under the 3-Clause BSD License.

use std::arch::x86_64::{_addcarryx_u64, _mulx_u64, _subborrow_u64};

#[rustfmt::skip]  // Save some vertical space
// BLS12-381 field prime modulus N, least significant limb first
const N: [u64; 6] = [
    0xb9fe_ffff_ffff_aaab, 0x1eab_fffe_b153_ffff, 0x6730_d2a0_f6b0_f624,
    0x6477_4b84_f385_12bf, 0x4b1b_a7b6_434b_acd7, 0x1a01_11ea_397f_e69a,
];

#[allow(clippy::needless_range_loop)]
// Assume properly reduced inputs and outputs
pub fn fe_add(result: &mut W6x64, a: &W6x64, b: &W6x64) {
    let mut sum = W6x64::default(); // Initializes to zero
    let mut carry = false;
    for i in 0..6 {
        let sum_car_a = a.v[i].overflowing_add(b.v[i]);
        let sum_car_b = sum_car_a.0.overflowing_add(u64::from(carry));
        sum.v[i] = sum_car_b.0;
        carry = sum_car_a.1 | sum_car_b.1
    }

    let mut trial = W6x64::default();
    let mut borrow = false;
    for i in 0..6 {
        // Note: a single overflowing_sub is sufficient because N[i]+borrow can never overflow
        let dif_bor = sum.v[i].overflowing_sub(N[i] + u64::from(borrow));
        trial.v[i] = dif_bor.0;
        borrow = dif_bor.1;
    }

    let select_sum = u64::from(borrow).wrapping_neg();
    for i in 0..6 {
        result.v[i] = (!select_sum & trial.v[i]) | (select_sum & sum.v[i]);
    }
}

#[rustfmt::skip]
// 2**384 - N, least significant limb first
const CORRECTION: [u64; 6] = [
    0x4601_0000_0000_5555, 0xe154_0001_4eac_0000, 0x98cf_2d5f_094f_09db,
    0x9b88_b47b_0c7a_ed40, 0xb4e4_5849_bcb4_5328, 0xe5fe_ee15_c680_1965,
];

#[allow(clippy::needless_range_loop)]
// Assume properly reduced inputs and outputs
pub fn fe_sub(result: &mut W6x64, a: &W6x64, b: &W6x64) {
    let mut diff = W6x64::default();
    let mut borrow_sub = false;
    for i in 0..6 {
        let dif_bor_a = a.v[i].overflowing_sub(b.v[i]);
        let dif_bor_b = dif_bor_a.0.overflowing_sub(u64::from(borrow_sub));
        diff.v[i] = dif_bor_b.0;
        borrow_sub = dif_bor_a.1 | dif_bor_b.1
    }

    let mask = u64::from(borrow_sub).wrapping_neg();
    let mut borrow_fix = false;
    for i in 0..6 {
        let dif_bor =
            // Note: a single overflowing_sub is sufficient because value+borrow can never overflow
            diff.v[i].overflowing_sub((mask & CORRECTION[i]) + u64::from(borrow_fix));
        result.v[i] = dif_bor.0;
        borrow_fix = dif_bor.1;
    }
}

#[rustfmt::skip]
// R^2 mod N, least significant limb first
const R_SQUARED: W6x64 = W6x64 {
    v: [0xf4df_1f34_1c34_1746, 0x0a76_e6a6_09d1_04f1, 0x8de5_476c_4c95_b6d5,
        0x67eb_88a9_939d_83c0, 0x9a79_3e85_b519_952d, 0x1198_8fe5_92ca_e3aa]
};

// Effectively a_mont = (a_norm * R) mod N
pub fn fe_to_mont(result: &mut W6x64, a: &[u64; 6]) {
    let a_w6x64 = W6x64 { v: *a };
    fe_mont_mul(&mut *result, &a_w6x64, &R_SQUARED);
}

#[rustfmt::skip]
// One, least significant limb first
const ONE: W6x64 = W6x64 {
    v: [0x0000_0000_0000_0001, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000,
        0x0000_0000_0000_0000, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000]
};

// Effectively a_norm = (a_mont * R^{-1}) mod N
pub fn fe_to_norm(result: &mut [u64; 6], a: &W6x64) {
    let mut result_w6x64 = W6x64::default();
    fe_mont_mul(&mut result_w6x64, &a, &ONE);
    *result = result_w6x64.v
}

#[derive(Default, Clone, Copy, Debug, PartialEq)] // Non constant-time Eq
#[repr(C)]
pub struct W6x64 {
    pub v: [u64; 6], // From least significant limb [0] to most significant [5]
}

const N_PRIME: u64 = 0x89f3_fffc_fffc_fffd;  // See constant.py

#[allow(clippy::cast_possible_truncation)]
// Effectively result_mont = (a_mont * b_mont * R^{-1}) mod N; Assume properly reduced input/output
pub fn fe_mont_mul(result: &mut W6x64, a: &W6x64, b: &W6x64) {
    let mut temp = [0_u64; 12];

    for i in 0..6 {
        let mut carry = 0_u64;
        for j in 0..6 {
            let hilo = u128::from(a.v[j]) * u128::from(b.v[i])
                + u128::from(temp[i + j])
                + u128::from(carry); // Note (2^64-1)*(2^64-1)+2*(2^64-1) = 2^128-1
            temp[i + j] = hilo as u64;
            carry = (hilo >> 64) as u64;
        }
        temp[i + 6] = temp[i + 6].wrapping_add(carry);

        let m: u64 = temp[i].wrapping_mul(N_PRIME);

        let mut carry = 0_u64;
        for j in 0..6 {
            let hilo =
                u128::from(m) * u128::from(N[j]) + u128::from(temp[i + j]) + u128::from(carry);
            temp[i + j] = hilo as u64;
            carry = (hilo >> 64) as u64;
        }
        temp[i + 6] = temp[i + 6].wrapping_add(carry);
    }

    let mut dec = [0_u64; 6];
    let mut borrow = false;
    for j in 0..6 {
        let (diff, borrow_tmp) = temp[j + 6].overflowing_sub(N[j] + u64::from(borrow));
        dec[j] = diff as u64;
        borrow = borrow_tmp;
    }

    let select_temp = u64::from(borrow).wrapping_neg();
    for j in 0..6 {
        result.v[j] = (select_temp & temp[j + 6]) | (!select_temp & dec[j]);
    }
}

macro_rules! full_add {
    ($carry_in:ident, $a:tt, $b:ident, $sum:ident, $carry_out:ident) => {
        let (sum0, carry0) = $a.overflowing_add($b);
        let (sum1, carry1) = sum0.overflowing_add(u64::from($carry_in));
        let $sum = sum1;
        let $carry_out = carry0 | carry1;
    };
}

// Same as above, but uses/expects previously declared (mut) sum variable
macro_rules! full_add2 {
    ($carry_in:ident, $a:ident, $b:ident, $sum:ident, $carry_out:ident) => {
        let (sum0, carry0) = $a.overflowing_add($b);
        let (sum1, carry1) = sum0.overflowing_add(u64::from($carry_in));
        $sum = sum1;
        let $carry_out = carry0 | carry1;
    };
}

macro_rules! full_sub {
    ($borrow_in:ident, $a:ident, $b:tt, $diff:tt, $borrow_out:ident) => {
        let (diff0, borrow0) = $a.overflowing_sub($b);
        let (diff1, borrow1) = diff0.overflowing_sub(u64::from($borrow_in));
        $diff = diff1;
        let $borrow_out = borrow0 | borrow1;
    };
}

macro_rules! mulx {
    ($a:tt, $b:tt, $hi:ident, $lo:ident) => {
        let hilo = u128::from($a) * u128::from($b);
        let $hi = (hilo >> 64) as u64;
        let $lo = hilo as u64;
    };
}

// Expected perf FOM = 12m + 6a + 13a/2 + 5*(12m + 13a/2) + 20 = 137/4.2GHz = 32.6ns; Act = 31.3nS
#[allow(clippy::similar_names, clippy::shadow_unrelated, unused_parens)]
pub fn fe_mont_mul_raw(result: &mut W6x64, a: &W6x64, b: &W6x64) {
    let (mut r10, mut r11, mut r12, mut r13, mut r14, mut r15);

    mulx!((a.v[0]), (b.v[0]), hi_a0b0, r10a);
    mulx!((a.v[1]), (b.v[0]), hi_a1b0, r11a);
    mulx!((a.v[2]), (b.v[0]), hi_a2b0, r12a);
    mulx!((a.v[3]), (b.v[0]), hi_a3b0, r13a);
    mulx!((a.v[4]), (b.v[0]), hi_a4b0, r14a);
    mulx!((a.v[5]), (b.v[0]), hi_a5b0, r15a);

    let (r11b, cf0_a) = r11a.overflowing_add(hi_a0b0);
    full_add!(cf0_a, r12a, hi_a1b0, r12b, cf0_b); // cf0_b -> carry flag, trace 0, step b
    full_add!(cf0_b, r13a, hi_a2b0, r13b, cf0_c);
    full_add!(cf0_c, r14a, hi_a3b0, r14b, cf0_d);
    full_add!(cf0_d, r15a, hi_a4b0, r15b, cf0_e);
    let r16a = hi_a5b0.wrapping_add(u64::from(cf0_e));

    let rdx = N_PRIME.wrapping_mul(r10a);

    mulx!((N[0]), rdx, hi_n0rdx, lo_n0rdx);
    mulx!((N[1]), rdx, hi_n1rdx, lo_n1rdx);
    mulx!((N[2]), rdx, hi_n2rdx, lo_n2rdx);
    mulx!((N[3]), rdx, hi_n3rdx, lo_n3rdx);
    mulx!((N[4]), rdx, hi_n4rdx, lo_n4rdx);
    mulx!((N[5]), rdx, hi_n5rdx, lo_n5rdx);

    let (r11c, cf1_a) = r11b.overflowing_add(hi_n0rdx);
    full_add!(cf1_a, r12b, hi_n1rdx, r12c, cf1_b);
    full_add!(cf1_b, r13b, hi_n2rdx, r13c, cf1_c);
    full_add!(cf1_c, r14b, hi_n3rdx, r14c, cf1_d);
    full_add!(cf1_d, r15b, hi_n4rdx, r15c, cf1_e);
    let r16b = hi_n5rdx.wrapping_add(u64::from(cf1_e)).wrapping_add(r16a);

    let (_, of1_a) = r10a.overflowing_add(lo_n0rdx);
    full_add2!(of1_a, r11c, lo_n1rdx, r10, of1_b); // of1_b -> overflow flag, trace 1, step b
    full_add2!(of1_b, r12c, lo_n2rdx, r11, of1_c);
    full_add2!(of1_c, r13c, lo_n3rdx, r12, of1_d);
    full_add2!(of1_d, r14c, lo_n4rdx, r13, of1_e);
    full_add2!(of1_e, r15c, lo_n5rdx, r14, of1_f);
    r15 = r16b.wrapping_add(u64::from(of1_f));

    for i in 1..6 {
        mulx!((a.v[0]), (b.v[i]), hi_a0bi, lo_a0bi);
        mulx!((a.v[1]), (b.v[i]), hi_a1bi, lo_a1bi);
        mulx!((a.v[2]), (b.v[i]), hi_a2bi, lo_a2bi);
        mulx!((a.v[3]), (b.v[i]), hi_a3bi, lo_a3bi);
        mulx!((a.v[4]), (b.v[i]), hi_a4bi, lo_a4bi);
        mulx!((a.v[5]), (b.v[i]), hi_a5bi, lo_a5bi);

        let res0 = r11.overflowing_add(hi_a0bi);
        let r11a = res0.0;
        let cf2_a = res0.1;
        full_add!(cf2_a, r12, hi_a1bi, r12a, cf2_b);
        full_add!(cf2_b, r13, hi_a2bi, r13a, cf2_c);
        full_add!(cf2_c, r14, hi_a3bi, r14a, cf2_d);
        full_add!(cf2_d, r15, hi_a4bi, r15a, cf2_e);
        let r16c = hi_a5bi.wrapping_add(u64::from(cf2_e));

        let res1 = r10.overflowing_add(lo_a0bi);
        let r10a = res1.0;
        let of2_a = res1.1;
        full_add!(of2_a, r11a, lo_a1bi, r11b, of2_b);
        full_add!(of2_b, r12a, lo_a2bi, r12b, of2_c);
        full_add!(of2_c, r13a, lo_a3bi, r13b, of2_d);
        full_add!(of2_d, r14a, lo_a4bi, r14b, of2_e);
        full_add!(of2_e, r15a, lo_a5bi, r15b, of2_f);
        let r16d = r16c.wrapping_add(u64::from(of2_f));

        let rdx = N_PRIME.wrapping_mul(r10a);

        mulx!((N[0]), rdx, hi_n0rdx, lo_n0rdx);
        mulx!((N[1]), rdx, hi_n1rdx, lo_n1rdx);
        mulx!((N[2]), rdx, hi_n2rdx, lo_n2rdx);
        mulx!((N[3]), rdx, hi_n3rdx, lo_n3rdx);
        mulx!((N[4]), rdx, hi_n4rdx, lo_n4rdx);
        mulx!((N[5]), rdx, hi_n5rdx, lo_n5rdx);

        let res2 = r11b.overflowing_add(hi_n0rdx);
        let r11c = res2.0;
        let cf3_a = res2.1;
        full_add!(cf3_a, r12b, hi_n1rdx, r12c, cf3_b);
        full_add!(cf3_b, r13b, hi_n2rdx, r13c, cf3_c);
        full_add!(cf3_c, r14b, hi_n3rdx, r14c, cf3_d);
        full_add!(cf3_d, r15b, hi_n4rdx, r15c, cf3_e);
        let r16e = hi_n5rdx.wrapping_add(u64::from(cf3_e)).wrapping_add(r16d);

        let res3 = r10a.overflowing_add(lo_n0rdx);
        let of3_a = res3.1;
        full_add2!(of3_a, r11c, lo_n1rdx, r10, of3_b);
        full_add2!(of3_b, r12c, lo_n2rdx, r11, of3_c);
        full_add2!(of3_c, r13c, lo_n3rdx, r12, of3_d);
        full_add2!(of3_d, r14c, lo_n4rdx, r13, of3_e);
        full_add2!(of3_e, r15c, lo_n5rdx, r14, of3_f);
        r15 = r16e.wrapping_add(u64::from(of3_f));
    }

    let res0 = r10.overflowing_sub(N[0]);
    let bor0 = res0.1;
    result.v[0] = res0.0;
    full_sub!(bor0, r11, (N[1]), (result.v[1]), bor1);
    full_sub!(bor1, r12, (N[2]), (result.v[2]), bor2);
    full_sub!(bor2, r13, (N[3]), (result.v[3]), bor3);
    full_sub!(bor3, r14, (N[4]), (result.v[4]), bor4);
    full_sub!(bor4, r15, (N[5]), (result.v[5]), bor5);

    // Mimic CMOV, but this is not constant-time!!
    if bor5 {
        result.v[0] = r10;
        result.v[1] = r11;
        result.v[2] = r12;
        result.v[3] = r13;
        result.v[4] = r14;
        result.v[5] = r15;
    }
}

#[allow(clippy::similar_names, clippy::too_many_lines)]
pub fn fe_mont_mul_intrinsics(result: &mut W6x64, a: &W6x64, b: &W6x64) {
    unsafe {
        let (mut r11, mut r12, mut r13, mut r14, mut r15, mut r16, mut rbx) = (0, 0, 0, 0, 0, 0, 0);

        let mut r10 = _mulx_u64(a.v[0], b.v[0], &mut r11);
        let mut rax = _mulx_u64(a.v[1], b.v[0], &mut r12);
        let tmp0 = r11.overflowing_add(rax);
        let cf0_a = u8::from(tmp0.1);
        r11 = tmp0.0;
        rax = _mulx_u64(a.v[2], b.v[0], &mut r13);
        let cf0_b = _addcarryx_u64(cf0_a, r12, rax, &mut r12);
        rax = _mulx_u64(a.v[3], b.v[0], &mut r14);
        let cf0_c = _addcarryx_u64(cf0_b, r13, rax, &mut r13);
        rax = _mulx_u64(a.v[4], b.v[0], &mut r15);
        let cf0_d = _addcarryx_u64(cf0_c, r14, rax, &mut r14);
        rax = _mulx_u64(a.v[5], b.v[0], &mut r16);
        let cf0_e = _addcarryx_u64(cf0_d, r15, rax, &mut r15);
        r16 = r16.wrapping_add(u64::from(cf0_e));

        let rdx = N_PRIME.wrapping_mul(r10);

        rax = _mulx_u64(N[0], rdx, &mut rbx);
        let tmp1 = r10.overflowing_add(rax);
        let cf1_a = u8::from(tmp1.1);
        let tmp2 = r11.overflowing_add(rbx);
        let of1_a = u8::from(tmp2.1);
        r11 = tmp2.0;
        r10 = _mulx_u64(N[1], rdx, &mut rbx);
        let cf1_b = _addcarryx_u64(cf1_a, r11, r10, &mut r10);
        let of1_b = _addcarryx_u64(of1_a, r12, rbx, &mut r12);
        r11 = _mulx_u64(N[2], rdx, &mut rbx);
        let cf1_c = _addcarryx_u64(cf1_b, r12, r11, &mut r11);
        let of1_c = _addcarryx_u64(of1_b, r13, rbx, &mut r13);
        r12 = _mulx_u64(N[3], rdx, &mut rbx);
        let cf1_d = _addcarryx_u64(cf1_c, r13, r12, &mut r12);
        let of1_d = _addcarryx_u64(of1_c, r14, rbx, &mut r14);
        r13 = _mulx_u64(N[4], rdx, &mut rbx);
        let cf1_e = _addcarryx_u64(cf1_d, r14, r13, &mut r13);
        let of1_e = _addcarryx_u64(of1_d, r15, rbx, &mut r15);
        r14 = _mulx_u64(N[5], rdx, &mut rbx);
        let cf1_f = _addcarryx_u64(cf1_e, r15, r14, &mut r14);
        r15 = r16.wrapping_add(u64::from(cf1_f));
        let _of1_f = _addcarryx_u64(of1_e, rbx, r15, &mut r15);

        for i in 1..6 {
            let mut r16 = 0;
            rax = _mulx_u64(a.v[0], b.v[i], &mut rbx);
            let tmp3 = r10.overflowing_add(rax);
            let cf2_a = tmp3.1 as u8;
            r10 = tmp3.0;
            let tmp4 = r11.overflowing_add(rbx);
            let of2_a = tmp4.1 as u8;
            r11 = tmp4.0;
            rax = _mulx_u64(a.v[1], b.v[i], &mut rbx);
            let cf2_b = _addcarryx_u64(cf2_a, r11, rax, &mut r11);
            let of2_b = _addcarryx_u64(of2_a, r12, rbx, &mut r12);
            rax = _mulx_u64(a.v[2], b.v[i], &mut rbx);
            let cf2_c = _addcarryx_u64(cf2_b, r12, rax, &mut r12);
            let of2_c = _addcarryx_u64(of2_b, r13, rbx, &mut r13);
            rax = _mulx_u64(a.v[3], b.v[i], &mut rbx);
            let cf2_d = _addcarryx_u64(cf2_c, r13, rax, &mut r13);
            let of2_d = _addcarryx_u64(of2_c, r14, rbx, &mut r14);
            rax = _mulx_u64(a.v[4], b.v[i], &mut rbx);
            let cf2_e = _addcarryx_u64(cf2_d, r14, rax, &mut r14);
            let of2_e = _addcarryx_u64(of2_d, r15, rbx, &mut r15);
            rax = _mulx_u64(a.v[5], b.v[i], &mut rbx);
            let cf2_f = _addcarryx_u64(cf2_e, r15, rax, &mut r15);
            let _of2_f = _addcarryx_u64(of2_e, rbx, r16, &mut r16);
            r16 = r16.wrapping_add(u64::from(cf2_f));

            let rdx = N_PRIME.wrapping_mul(r10);

            rax = _mulx_u64(N[0], rdx, &mut rbx);
            let tmp5 = r10.overflowing_add(rax);
            let cf3_a = u8::from(tmp5.1);
            let tmp6 = r11.overflowing_add(rbx);
            let of3_a = u8::from(tmp6.1);
            r11 = tmp6.0;
            r10 = _mulx_u64(N[1], rdx, &mut rbx);
            let cf3_b = _addcarryx_u64(cf3_a, r11, r10, &mut r10);
            let of3_b = _addcarryx_u64(of3_a, r12, rbx, &mut r12);
            r11 = _mulx_u64(N[2], rdx, &mut rbx);
            let cf3_c = _addcarryx_u64(cf3_b, r12, r11, &mut r11);
            let of3_c = _addcarryx_u64(of3_b, r13, rbx, &mut r13);
            r12 = _mulx_u64(N[3], rdx, &mut rbx);
            let cf3_d = _addcarryx_u64(cf3_c, r13, r12, &mut r12);
            let of3_d = _addcarryx_u64(of3_c, r14, rbx, &mut r14);
            r13 = _mulx_u64(N[4], rdx, &mut rbx);
            let cf3_e = _addcarryx_u64(cf3_d, r14, r13, &mut r13);
            let of3_e = _addcarryx_u64(of3_d, r15, rbx, &mut r15);
            r14 = _mulx_u64(N[5], rdx, &mut rbx);
            let cf3_f = _addcarryx_u64(cf3_e, r15, r14, &mut r14);
            r15 = r16.wrapping_add(u64::from(cf3_f));
            let _of3_f = _addcarryx_u64(of3_e, rbx, r15, &mut r15);
        }

        let res = r10.overflowing_sub(N[0]);
        result.v[0] = res.0;
        let bor1 = _subborrow_u64(u8::from(res.1), r11, N[1], &mut result.v[1]);
        let bor2 = _subborrow_u64(bor1, r12, N[2], &mut result.v[2]);
        let bor3 = _subborrow_u64(bor2, r13, N[3], &mut result.v[3]);
        let bor4 = _subborrow_u64(bor3, r14, N[4], &mut result.v[4]);
        let bor5 = _subborrow_u64(bor4, r15, N[5], &mut result.v[5]);

        // Mimic CMOV, but this is not constant-time!!
        if bor5 != 0 {
            result.v[0] = r10;
            result.v[1] = r11;
            result.v[2] = r12;
            result.v[3] = r13;
            result.v[4] = r14;
            result.v[5] = r15;
        }
    }
}
