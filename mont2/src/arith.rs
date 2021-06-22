use std::arch::x86_64::{_addcarryx_u64, _mulx_u64, _subborrow_u64};

#[derive(Default, Clone, Copy, Debug, PartialEq)] // Non constant-time Eq
#[repr(C)]
pub struct W6x64 {
    pub v: [u64; 6], // From least significant [0] to most significant [5]
}

#[rustfmt::skip]
// BLS12-381 field prime modulus N, least significant portion first
const N: [u64; 6] = [
    0xb9fe_ffff_ffff_aaab, 0x1eab_fffe_b153_ffff, 0x6730_d2a0_f6b0_f624,
    0x6477_4b84_f385_12bf, 0x4b1b_a7b6_434b_acd7, 0x1a01_11ea_397f_e69a,
];

#[allow(clippy::needless_range_loop)]
#[inline] // Assume properly reduced inputs and outputs
pub fn fe_add(result: &mut W6x64, a: &W6x64, b: &W6x64) {
    let mut sum = W6x64::default(); // Initializes to zero
    let mut carry = false;
    for i in 0..6 {
        let sum_car_a = a.v[i].overflowing_add(b.v[i]);
        let sum_car_b = sum_car_a.0.overflowing_add(if carry { 1 } else { 0 });
        sum.v[i] = sum_car_b.0;
        carry = sum_car_a.1 | sum_car_b.1
    }

    let mut trial = W6x64::default();
    let mut borrow = false;
    for i in 0..6 {
        // Note: a single overflowing_sub is sufficient because N[i]+borrow can never overflow
        let dif_bor = sum.v[i].overflowing_sub(N[i] + if borrow { 1 } else { 0 });
        trial.v[i] = dif_bor.0;
        borrow = dif_bor.1;
    }

    let select_sum = u64::from(borrow).wrapping_neg();
    for i in 0..6 {
        result.v[i] = (!select_sum & trial.v[i]) | (select_sum & sum.v[i]);
    }
}

#[rustfmt::skip]
// 2**384 - N, least significant portion first
const CORRECTION: [u64; 6] = [
    0x4601_0000_0000_5555, 0xe154_0001_4eac_0000, 0x98cf_2d5f_094f_09db,
    0x9b88_b47b_0c7a_ed40, 0xb4e4_5849_bcb4_5328, 0xe5fe_ee15_c680_1965,
];

#[allow(clippy::needless_range_loop)]
#[inline] // Assume properly reduced inputs and outputs
pub fn fe_sub(result: &mut W6x64, a: &W6x64, b: &W6x64) {
    let mut diff = W6x64::default();
    let mut borrow_sub = false;
    for i in 0..6 {
        let dif_bor_a = a.v[i].overflowing_sub(b.v[i]);
        let dif_bor_b = dif_bor_a.0.overflowing_sub(if borrow_sub { 1 } else { 0 });
        diff.v[i] = dif_bor_b.0;
        borrow_sub = dif_bor_a.1 | dif_bor_b.1
    }

    let mask = u64::from(borrow_sub).wrapping_neg();
    let mut borrow_fix = false;
    for i in 0..6 {
        let dif_bor =
            // Note: a single overflowing_sub is sufficient because value+borrow can never overflow
            diff.v[i].overflowing_sub((mask & CORRECTION[i]) + if borrow_fix { 1 } else { 0 });
        result.v[i] = dif_bor.0;
        borrow_fix = dif_bor.1;
    }
}

#[rustfmt::skip]
// R^2 mod N, least significant portion first
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
// One, least significant portion first
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

const N_PRIME: u64 = 0x89f3_fffc_fffc_fffd;

#[allow(clippy::cast_possible_truncation)]
// Effectively result_mont = (a_mont * b_mont * R^{-1}) mod N; Assume properly reduced inputs and outputs
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
        temp[i + 6] += carry;

        let m: u64 = temp[i].wrapping_mul(N_PRIME);

        let mut carry = 0_u64;
        for j in 0..6 {
            let hilo =
                u128::from(m) * u128::from(N[j]) + u128::from(temp[i + j]) + u128::from(carry);
            temp[i + j] = hilo as u64;
            carry = (hilo >> 64) as u64;
        }
        temp[i + 6] += carry;
    }

    let mut dec = [0_u64; 6];
    let mut borrow = 0_u64;
    for j in 0..6 {
        let (diff, borrow_t0) = temp[j + 6].overflowing_sub(N[j] + borrow);
        dec[j] = diff as u64;
        borrow = u64::from(borrow_t0);
    }

    let select_temp = borrow.wrapping_neg();
    for j in 0..6 {
        result.v[j] = (select_temp & temp[j + 6]) | (!select_temp & dec[j]);
    }
}

// Expected perf FOM = 6*(6+6) single-cycle mul and 6*(2*13) half-cycle adds = 150 cycles * (1/4.5GHz) = 36nS; Actual = 35nS
#[rustfmt::skip]
pub fn fe_mont_mul_raw(result: &mut W6x64, a: &W6x64, b: &W6x64) {
    let hilo = u128::from(a.v[0]) * u128::from(b.v[0]); let hi_a0b0 = (hilo >> 64) as u64; let r10a = hilo as u64;
    let hilo = u128::from(a.v[1]) * u128::from(b.v[0]); let hi_a1b0 = (hilo >> 64) as u64; let r11a = hilo as u64;
    let hilo = u128::from(a.v[2]) * u128::from(b.v[0]); let hi_a2b0 = (hilo >> 64) as u64; let r12a = hilo as u64;
    let hilo = u128::from(a.v[3]) * u128::from(b.v[0]); let hi_a3b0 = (hilo >> 64) as u64; let r13a = hilo as u64;
    let hilo = u128::from(a.v[4]) * u128::from(b.v[0]); let hi_a4b0 = (hilo >> 64) as u64; let r14a = hilo as u64;
    let hilo = u128::from(a.v[5]) * u128::from(b.v[0]); let hi_a5b0 = (hilo >> 64) as u64; let r15a = hilo as u64;

    let res = r11a.overflowing_add(hi_a0b0); let r11b = res.0; let cf0_a = res.1;
    let res1 = r12a.overflowing_add(hi_a1b0); let res2 = res1.0.overflowing_add(u64::from(cf0_a)); let r12b = res2.0; let cf0_b = res1.1 | res2.1;
    let res1 = r13a.overflowing_add(hi_a2b0); let res2 = res1.0.overflowing_add(u64::from(cf0_b)); let r13b = res2.0; let cf0_c = res1.1 | res2.1;
    let res1 = r14a.overflowing_add(hi_a3b0); let res2 = res1.0.overflowing_add(u64::from(cf0_c)); let r14b = res2.0; let cf0_d = res1.1 | res2.1;
    let res1 = r15a.overflowing_add(hi_a4b0); let res2 = res1.0.overflowing_add(u64::from(cf0_d)); let r15b = res2.0; let cf0_e = res1.1 | res2.1;
    let r16a = hi_a5b0.wrapping_add(u64::from(cf0_e));

    let rdx = N_PRIME.wrapping_mul(r10a);

    let hilo = u128::from(N[0]) * u128::from(rdx); let hi_n0rdx = (hilo >> 64) as u64; let lo_n0rdx = hilo as u64;
    let hilo = u128::from(N[1]) * u128::from(rdx); let hi_n1rdx = (hilo >> 64) as u64; let lo_n1rdx = hilo as u64;
    let hilo = u128::from(N[2]) * u128::from(rdx); let hi_n2rdx = (hilo >> 64) as u64; let lo_n2rdx = hilo as u64;
    let hilo = u128::from(N[3]) * u128::from(rdx); let hi_n3rdx = (hilo >> 64) as u64; let lo_n3rdx = hilo as u64;
    let hilo = u128::from(N[4]) * u128::from(rdx); let hi_n4rdx = (hilo >> 64) as u64; let lo_n4rdx = hilo as u64;
    let hilo = u128::from(N[5]) * u128::from(rdx); let hi_n5rdx = (hilo >> 64) as u64; let lo_n5rdx = hilo as u64;

    let res = r11b.overflowing_add(hi_n0rdx); let r11c = res.0; let cf1_a = res.1;
    let res1 = r12b.overflowing_add(hi_n1rdx); let res2 = res1.0.overflowing_add(u64::from(cf1_a)); let r12c = res2.0; let cf1_b = res1.1 | res2.1;
    let res1 = r13b.overflowing_add(hi_n2rdx); let res2 = res1.0.overflowing_add(u64::from(cf1_b)); let r13c = res2.0; let cf1_c = res1.1 | res2.1;
    let res1 = r14b.overflowing_add(hi_n3rdx); let res2 = res1.0.overflowing_add(u64::from(cf1_c)); let r14c = res2.0; let cf1_d = res1.1 | res2.1;
    let res1 = r15b.overflowing_add(hi_n4rdx); let res2 = res1.0.overflowing_add(u64::from(cf1_d)); let r15c = res2.0; let cf1_e = res1.1 | res2.1;
    let r16b = r16a.wrapping_add(hi_n5rdx);

    let res = r10a.overflowing_add(lo_n0rdx); let of1_a = res.1 as u64;
    let res1 = r11c.overflowing_add(lo_n1rdx); let res2 = res1.0.overflowing_add(u64::from(of1_a)); let mut r10 = res2.0; let of1_b = res1.1 | res2.1;
    let res1 = r12c.overflowing_add(lo_n2rdx); let res2 = res1.0.overflowing_add(u64::from(of1_b)); let mut r11 = res2.0; let of1_c = res1.1 | res2.1;
    let res1 = r13c.overflowing_add(lo_n3rdx); let res2 = res1.0.overflowing_add(u64::from(of1_c)); let mut r12 = res2.0; let of1_d = res1.1 | res2.1;
    let res1 = r14c.overflowing_add(lo_n4rdx); let res2 = res1.0.overflowing_add(u64::from(of1_d)); let mut r13 = res2.0; let of1_e = res1.1 | res2.1;
    let res1 = r15c.overflowing_add(lo_n5rdx); let res2 = res1.0.overflowing_add(u64::from(of1_e)); let mut r14 = res2.0; let of1_f = res1.1 | res2.1;

    let mut r15 = r16b.wrapping_add(u64::from(cf1_e) + u64::from(of1_f));

    for i in 1..6 {
        let hilo = u128::from(a.v[0]) * u128::from(b.v[i]); let hi_a0vi = (hilo >> 64) as u64; let lo_a0vi = hilo as u64;
        let hilo = u128::from(a.v[1]) * u128::from(b.v[i]); let hi_a1vi = (hilo >> 64) as u64; let lo_a1vi = hilo as u64;
        let hilo = u128::from(a.v[2]) * u128::from(b.v[i]); let hi_a2vi = (hilo >> 64) as u64; let lo_a2vi = hilo as u64;
        let hilo = u128::from(a.v[3]) * u128::from(b.v[i]); let hi_a3vi = (hilo >> 64) as u64; let lo_a3vi = hilo as u64;
        let hilo = u128::from(a.v[4]) * u128::from(b.v[i]); let hi_a4vi = (hilo >> 64) as u64; let lo_a4vi = hilo as u64;
        let hilo = u128::from(a.v[5]) * u128::from(b.v[i]); let hi_a5vi = (hilo >> 64) as u64; let lo_a5vi = hilo as u64;

        let res0 = r11.overflowing_add(hi_a0vi); let r11a = res0.0; let cf2_a = res0.1;
        let res0 = r12.overflowing_add(hi_a1vi); let res1 = res0.0.overflowing_add(u64::from(cf2_a)); let r12a = res1.0; let cf2_b = res0.1 | res1.1;
        let res0 = r13.overflowing_add(hi_a2vi); let res1 = res0.0.overflowing_add(u64::from(cf2_b)); let r13a = res1.0; let cf2_c = res0.1 | res1.1;
        let res0 = r14.overflowing_add(hi_a3vi); let res1 = res0.0.overflowing_add(u64::from(cf2_c)); let r14a = res1.0; let cf2_d = res0.1 | res1.1;
        let res0 = r15.overflowing_add(hi_a4vi); let res1 = res0.0.overflowing_add(u64::from(cf2_d)); let r15a = res1.0; let cf2_e = res0.1 | res1.1;
        let r16c = hi_a5vi.wrapping_add(u64::from(cf2_e));

        let res8 = r10.overflowing_add(lo_a0vi); let r10a = res8.0; let of2_a = res8.1;
        let res8 = r11a.overflowing_add(lo_a1vi); let res9 = res8.0.overflowing_add(u64::from(of2_a)); let r11a = res9.0; let of2_b = res8.1 | res9.1;
        let res8 = r12a.overflowing_add(lo_a2vi); let res9 = res8.0.overflowing_add(u64::from(of2_b)); let r12b = res9.0; let of2_c = res8.1 | res9.1;
        let res8 = r13a.overflowing_add(lo_a3vi); let res9 = res8.0.overflowing_add(u64::from(of2_c)); let r13b = res9.0; let of2_d = res8.1 | res9.1;
        let res8 = r14a.overflowing_add(lo_a4vi); let res9 = res8.0.overflowing_add(u64::from(of2_d)); let r14b = res9.0; let of2_e = res8.1 | res9.1;
        let res8 = r15a.overflowing_add(lo_a5vi); let res9 = res8.0.overflowing_add(u64::from(of2_e)); let r15b = res9.0; let of2_f = res8.1 | res9.1;
        let r16d = r16c.wrapping_add(of2_f as u64);

        let rdx = N_PRIME.wrapping_mul(r10a);

        let hilo = u128::from(N[0]) * u128::from(rdx); let hi_n0rdx = (hilo >> 64) as u64; let lo_n0rdx = hilo as u64;
        let hilo = u128::from(N[1]) * u128::from(rdx); let hi_n1rdx = (hilo >> 64) as u64; let lo_n1rdx = hilo as u64;
        let hilo = u128::from(N[2]) * u128::from(rdx); let hi_n2rdx = (hilo >> 64) as u64; let lo_n2rdx = hilo as u64;
        let hilo = u128::from(N[3]) * u128::from(rdx); let hi_n3rdx = (hilo >> 64) as u64; let lo_n3rdx = hilo as u64;
        let hilo = u128::from(N[4]) * u128::from(rdx); let hi_n4rdx = (hilo >> 64) as u64; let lo_n4rdx = hilo as u64;
        let hilo = u128::from(N[5]) * u128::from(rdx); let hi_n5rdx = (hilo >> 64) as u64; let lo_n5rdx = hilo as u64;

        let res0 = r11a.overflowing_add(hi_n0rdx); let r11b = res0.0; let cf3_a = res0.1;
        let res0 = r12b.overflowing_add(hi_n1rdx); let res1 = res0.0.overflowing_add(u64::from(cf3_a)); let r12c = res1.0; let cf3_b = res0.1 | res1.1;
        let res0 = r13b.overflowing_add(hi_n2rdx); let res1 = res0.0.overflowing_add(u64::from(cf3_b)); let r13c = res1.0; let cf3_c = res0.1 | res1.1;
        let res0 = r14b.overflowing_add(hi_n3rdx); let res1 = res0.0.overflowing_add(u64::from(cf3_c)); let r14c = res1.0; let cf3_d = res0.1 | res1.1;
        let res0 = r15b.overflowing_add(hi_n4rdx); let res1 = res0.0.overflowing_add(u64::from(cf3_d)); let r15c = res1.0; let cf3_e = res0.1 | res1.1;
        let r16e = r16d.wrapping_add(hi_n5rdx);

        let res8 = r10a.overflowing_add(lo_n0rdx); let of3_a = res8.1;
        let res8 = r11b.overflowing_add(lo_n1rdx); let res9 = res8.0.overflowing_add(u64::from(of3_a)); r10 = res9.0; let of3_b = res8.1 | res9.1;
        let res8 = r12c.overflowing_add(lo_n2rdx); let res9 = res8.0.overflowing_add(u64::from(of3_b)); r11 = res9.0; let of3_c = res8.1 | res9.1;
        let res8 = r13c.overflowing_add(lo_n3rdx); let res9 = res8.0.overflowing_add(u64::from(of3_c)); r12 = res9.0; let of3_d = res8.1 | res9.1;
        let res8 = r14c.overflowing_add(lo_n4rdx); let res9 = res8.0.overflowing_add(u64::from(of3_d)); r13 = res9.0; let of3_e = res8.1 | res9.1;
        let res8 = r15c.overflowing_add(lo_n5rdx); let res9 = res8.0.overflowing_add(u64::from(of3_e)); r14 = res9.0; let of3_f = res8.1 | res9.1;

        r15 = r16e.wrapping_add(u64::from(cf3_e) + u64::from(of3_f));
    }

    let res0 = r10.overflowing_sub(N[0]); let bor0 = res0.1;
    let res1a = r11.overflowing_sub(N[1]); let res1b = res1a.0.overflowing_sub(u64::from(bor0)); let bor1 = res1a.1 | res1b.1;
    let res2a = r12.overflowing_sub(N[2]); let res2b = res2a.0.overflowing_sub(u64::from(bor1)); let bor2 = res2a.1 | res2b.1;
    let res3a = r13.overflowing_sub(N[3]); let res3b = res3a.0.overflowing_sub(u64::from(bor2)); let bor3 = res3a.1 | res3b.1;
    let res4a = r14.overflowing_sub(N[4]); let res4b = res4a.0.overflowing_sub(u64::from(bor3)); let bor4 = res4a.1 | res4b.1;
    let res5a = r15.overflowing_sub(N[5]); let res5b = res5a.0.overflowing_sub(u64::from(bor4)); let bor5 = res5a.1 | res5b.1;

    if bor5 {
        result.v[0] = r10;
        result.v[1] = r11;
        result.v[2] = r12;
        result.v[3] = r13;
        result.v[4] = r14;
        result.v[5] = r15;
    } else {
        result.v[0] = res0.0;
        result.v[1] = res1b.0;
        result.v[2] = res2b.0;
        result.v[3] = res3b.0;
        result.v[4] = res4b.0;
        result.v[5] = res5b.0;
    }
}

pub fn fe_mont_mul_intrinsics(result: &mut W6x64, a: &W6x64, b: &W6x64) {
    unsafe {
        let (mut r10, mut r11, mut r12, mut r13, mut r14, mut r15, mut rbx) = (0, 0, 0, 0, 0, 0, 0);
        let mut rax;

        for i in 0..6 {
            let mut rbp = 0;
            rax = _mulx_u64(a.v[0], b.v[i], &mut rbx);
            //eprintln!("b.v[0] = {:x?}", b.v[0]);
            let car_0a = _addcarryx_u64(0_u8, r10, rax, &mut r10);
            //eprintln!("r10 {:x?}", r10);
            let ovf_0a = _addcarryx_u64(0_u8, r11, rbx, &mut r11);
            rax = _mulx_u64(a.v[1], b.v[i], &mut rbx);
            let car_1a = _addcarryx_u64(car_0a, r11, rax, &mut r11);
            let ovf_1a = _addcarryx_u64(ovf_0a, r12, rbx, &mut r12);
            rax = _mulx_u64(a.v[2], b.v[i], &mut rbx);
            let car_2a = _addcarryx_u64(car_1a, r12, rax, &mut r12);
            let ovf_2a = _addcarryx_u64(ovf_1a, r13, rbx, &mut r13);
            rax = _mulx_u64(a.v[3], b.v[i], &mut rbx);
            let car_3a = _addcarryx_u64(car_2a, r13, rax, &mut r13);
            let ovf_3a = _addcarryx_u64(ovf_2a, r14, rbx, &mut r14);
            rax = _mulx_u64(a.v[4], b.v[i], &mut rbx);
            let car_4a = _addcarryx_u64(car_3a, r14, rax, &mut r14);
            let ovf_4a = _addcarryx_u64(ovf_3a, r15, rbx, &mut r15);
            rax = _mulx_u64(a.v[5], b.v[i], &mut rbx);
            let car_5a = _addcarryx_u64(car_4a, r15, rax, &mut r15);
            let _ovf_5a = _addcarryx_u64(ovf_4a, rbx, rbp, &mut rbp);
            let _car_5a = _addcarryx_u64(car_5a, 0, rbp, &mut rbp);

            let rdx = N_PRIME.wrapping_mul(r10);

            rax = _mulx_u64(N[0], rdx, &mut rbx);
            let car_0b = _addcarryx_u64(0_u8, r10, rax, &mut rax);
            let ovf_0b = _addcarryx_u64(0_u8, r11, rbx, &mut r11);
            r10 = _mulx_u64(N[1], rdx, &mut rbx);
            let car_1b = _addcarryx_u64(car_0b, r11, r10, &mut r10);
            let ovf_1b = _addcarryx_u64(ovf_0b, r12, rbx, &mut r12);
            r11 = _mulx_u64(N[2], rdx, &mut rbx);
            let car_2b = _addcarryx_u64(car_1b, r12, r11, &mut r11);
            let ovf_2b = _addcarryx_u64(ovf_1b, r13, rbx, &mut r13);
            r12 = _mulx_u64(N[3], rdx, &mut rbx);
            let car_3b = _addcarryx_u64(car_2b, r13, r12, &mut r12);
            let ovf_3b = _addcarryx_u64(ovf_2b, r14, rbx, &mut r14);
            r13 = _mulx_u64(N[4], rdx, &mut rbx);
            let car_4b = _addcarryx_u64(car_3b, r14, r13, &mut r13);
            let ovf_4b = _addcarryx_u64(ovf_3b, r15, rbx, &mut r15);
            r14 = _mulx_u64(N[5], rdx, &mut rbx);
            let car_5b = _addcarryx_u64(car_4b, r15, r14, &mut r14);
            let _ovf_5b = _addcarryx_u64(car_5b, rbp, 0_u64, &mut r15);
            let _ovf_5b = _addcarryx_u64(ovf_4b, rbx, r15, &mut r15);
        }

        let res = r10.overflowing_sub(N[0]);
        result.v[0] = res.0;
        let bor1 = _subborrow_u64(u8::from(res.1), r11, N[1], &mut result.v[1]);
        let bor2 = _subborrow_u64(bor1, r12, N[2], &mut result.v[2]);
        let bor3 = _subborrow_u64(bor2, r13, N[3], &mut result.v[3]);
        let bor4 = _subborrow_u64(bor3, r14, N[4], &mut result.v[4]);
        let bor5 = _subborrow_u64(bor4, r15, N[5], &mut result.v[5]);

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
