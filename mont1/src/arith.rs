#[repr(C)]
#[derive(Default, Clone, Copy, Debug, PartialEq)]
pub struct W6x64 {
    pub v: [u64; 6], // From least significant [0] to most significant [5]
}

// BLS12-381 field prime, least significant portion first
#[rustfmt::skip]
const N: [u64; 6] = [
    0xb9fe_ffff_ffff_aaab, 0x1eab_fffe_b153_ffff, 0x6730_d2a0_f6b0_f624,
    0x6477_4b84_f385_12bf, 0x4b1b_a7b6_434b_acd7, 0x1a01_11ea_397f_e69a,
];

#[allow(clippy::needless_range_loop)]
#[inline]
pub fn fe_add(result: &mut W6x64, a: &W6x64, b: &W6x64) {
    let mut sum = W6x64::default();
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
        let dif_bor = sum.v[i].overflowing_sub(N[i] + if borrow { 1 } else { 0 });
        trial.v[i] = dif_bor.0;
        borrow = dif_bor.1;
    }

    let select_sum = u64::from(borrow).wrapping_neg();
    for i in 0..6 {
        result.v[i] = (!select_sum & trial.v[i]) | (select_sum & sum.v[i]);
    }
}

// 2**384 - N, least significant portion first
#[rustfmt::skip]
const CORRECTION: [u64; 6] = [
    0x4601_0000_0000_5555, 0xe154_0001_4eac_0000, 0x98cf_2d5f_094f_09db,
    0x9b88_b47b_0c7a_ed40, 0xb4e4_5849_bcb4_5328, 0xe5fe_ee15_c680_1965,
];

#[allow(clippy::needless_range_loop)]
#[inline]
pub fn fe_sub(result: &mut W6x64, a: &W6x64, b: &W6x64) {
    let mut diff = W6x64::default();
    let mut borrow_diff = false;
    for i in 0..6 {
        let dif_bor_a = a.v[i].overflowing_sub(b.v[i]);
        let dif_bor_b = dif_bor_a.0.overflowing_sub(if borrow_diff { 1 } else { 0 });
        diff.v[i] = dif_bor_b.0;
        borrow_diff = dif_bor_a.1 | dif_bor_b.1
    }

    let mask = u64::from(borrow_diff).wrapping_neg();
    let mut borrow_fix = false;
    for i in 0..6 {
        let dif_bor =
            diff.v[i].overflowing_sub((mask & CORRECTION[i]) + if borrow_fix { 1 } else { 0 });
        result.v[i] = dif_bor.0;
        borrow_fix = dif_bor.1;
    }
}

// R^2 mod N, least significant portion first
#[rustfmt::skip]
const R_SQUARED: W6x64 = W6x64 {
    v: [0xf4df_1f34_1c34_1746, 0x0a76_e6a6_09d1_04f1, 0x8de5_476c_4c95_b6d5,
        0x67eb_88a9_939d_83c0, 0x9a79_3e85_b519_952d, 0x1198_8fe5_92ca_e3aa]
};

pub fn fe_to_mont(result: &mut W6x64, a: &[u64; 6]) {
    let w6x64_a = W6x64 { v: *a };
    fe_mont_mul(&mut *result, &w6x64_a, &R_SQUARED);
}

// One, least significant portion first
#[rustfmt::skip]
const ONE: W6x64 = W6x64 {
    v: [0x0000_0000_0000_0001, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000,
        0x0000_0000_0000_0000, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000]
};

pub fn fe_to_norm(result: &mut [u64; 6], a: &W6x64) {
    let mut w6x64_result = W6x64::default();
    fe_mont_mul(&mut w6x64_result, &a, &ONE);
    *result = w6x64_result.v
}

const N_PRIME: u64 = 0x89f3_fffc_fffc_fffd;

#[allow(clippy::cast_possible_truncation)]
#[inline]
pub fn fe_mont_mul(result: &mut W6x64, a: &W6x64, b: &W6x64) {
    let mut temp = [0_u64; 12];

    for i in 0..6 {
        let mut carry = 0_u64;
        for j in 0..6 {
            let hilo = u128::from(a.v[j]) * u128::from(b.v[i])
                + u128::from(temp[i + j])
                + u128::from(carry);
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
