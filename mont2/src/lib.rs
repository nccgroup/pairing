// Copyright 2021 Eric Schorn; Licensed under the 3-Clause BSD License.

#![deny(clippy::all)]
#![deny(clippy::pedantic)]
#![deny(clippy::cargo)]

pub mod arith;

// Benchmark with `cargo bench` but prefer the following:
//   RUSTFLAGS="--emit asm -C target-cpu=native" cargo bench
//   RUSTFLAGS="--emit asm -C target-feature=+bmi2" cargo bench

extern "C" {
    pub fn fe_mont_mul_asm(result: &mut u64, a: &u64, b: &u64);
}

#[cfg(test)]
#[macro_use]
extern crate lazy_static;

#[cfg(test)]
mod tests {
    use crate::arith::{
        fe_add, fe_mont_mul, fe_mont_mul_intrinsics, fe_mont_mul_raw, fe_sub, fe_to_mont,
        fe_to_norm, W6x64,
    };
    use crate::fe_mont_mul_asm;
    use num_bigint::BigUint;
    use num_traits::Num;
    use rand::Rng;
    use std::convert::TryInto;

    lazy_static! {
        static ref R: BigUint = BigUint::from(1_u64) << 384;
        static ref MODULUS: BigUint = BigUint::from_str_radix(
            "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab",
            16).unwrap();
    }

    fn rnd_big_mod_n() -> BigUint {
        let mut rnd_bytes = [0_u8; 64];
        rand::thread_rng().fill(&mut rnd_bytes[..]);
        BigUint::from_bytes_le(&rnd_bytes) % &(*MODULUS)
    }

    fn big_to_6u64(x: &BigUint) -> [u64; 6] {
        let mut bytes = [0_u8; 48];
        bytes[0..(((7 + x.bits()) / 8) as usize)].clone_from_slice(&x.to_bytes_le());
        let mut result = [0_u64; 6];
        for i in 0..6 {
            result[i] = u64::from_le_bytes(bytes[i * 8..(i + 1) * 8].try_into().unwrap());
        }
        result
    }

    lazy_static! { static ref R_INV: BigUint = BigUint::from_str_radix(
        "14fec701e8fb0ce9ed5e64273c4f538b1797ab1458a88de9343ea97914956dc87fe11274d898fafbf4d38259380b4820",
        16).unwrap();
    }

    // Reference for Montgomery multiplication with BigUint
    fn mont_mul_biguint(result: &mut BigUint, a: &BigUint, b: &BigUint) {
        *result = (a * b * &(*R_INV)) % &(*MODULUS);
    }

    #[test]
    fn test_fe_add() {
        let mut actual_mont = W6x64::default();
        let mut actual_norm = [0_u64; 6];
        let mut a_mont = W6x64::default();
        let mut b_mont = W6x64::default();

        for _i in 0..1_000_000 {
            let a_big = rnd_big_mod_n();
            let b_big = rnd_big_mod_n();
            fe_to_mont(&mut a_mont, &big_to_6u64(&a_big));
            fe_to_mont(&mut b_mont, &big_to_6u64(&b_big));

            let expected = (&a_big + &b_big) % &(*MODULUS);
            fe_add(&mut actual_mont, &a_mont, &b_mont);

            fe_to_norm(&mut actual_norm, &actual_mont);
            assert_eq!(big_to_6u64(&expected), actual_norm);
        }
    }

    #[test]
    fn test_fe_sub() {
        let mut actual_mont = W6x64::default();
        let mut actual_norm = [0_u64; 6];
        let mut a_mont = W6x64::default();
        let mut b_mont = W6x64::default();

        for _i in 0..1_000_000 {
            let a_big = rnd_big_mod_n();
            let b_big = rnd_big_mod_n();
            fe_to_mont(&mut a_mont, &big_to_6u64(&a_big));
            fe_to_mont(&mut b_mont, &big_to_6u64(&b_big));

            let expected = (&a_big + &(*MODULUS) - &b_big) % &(*MODULUS);
            fe_sub(&mut actual_mont, &a_mont, &b_mont);

            fe_to_norm(&mut actual_norm, &actual_mont);
            assert_eq!(big_to_6u64(&expected), actual_norm);
        }
    }

    #[test]
    fn test_fe_to_mont_to_norm() {
        let mut a_mont = W6x64::default();
        let mut actual_norm = [0_u64; 6];

        for _i in 0..1_000_000 {
            let a_big = rnd_big_mod_n();
            fe_to_mont(&mut a_mont, &big_to_6u64(&a_big));
            fe_to_norm(&mut actual_norm, &a_mont);
            assert_eq!(big_to_6u64(&a_big), actual_norm);
        }
    }

    #[test]
    fn test_fe_mont_mul() {
        let mut actual_mont = W6x64::default();
        let mut exp_mont = BigUint::default();
        let mut a_mont = W6x64::default();
        let mut b_mont = W6x64::default();

        for _i in 0..1_000_000 {
            let a_big = rnd_big_mod_n();
            let b_big = rnd_big_mod_n();
            fe_to_mont(&mut a_mont, &big_to_6u64(&a_big));
            fe_to_mont(&mut b_mont, &big_to_6u64(&b_big));

            fe_mont_mul(&mut actual_mont, &a_mont, &b_mont);

            mont_mul_biguint(&mut exp_mont, &(a_big * &(*R)), &(b_big * &(*R)));
            assert_eq!(big_to_6u64(&exp_mont), actual_mont.v);
        }
    }

    #[test]
    fn test_fe_mont_mul_intrinsics() {
        let mut actual_mont = W6x64::default();
        let mut exp_mont = BigUint::default();
        let mut a_mont = W6x64::default();
        let mut b_mont = W6x64::default();

        for _i in 0..1_000_000 {
            let a_big = rnd_big_mod_n();
            let b_big = rnd_big_mod_n();
            fe_to_mont(&mut a_mont, &big_to_6u64(&a_big));
            fe_to_mont(&mut b_mont, &big_to_6u64(&b_big));

            fe_mont_mul_intrinsics(&mut actual_mont, &a_mont, &b_mont);

            mont_mul_biguint(&mut exp_mont, &(a_big * &(*R)), &(b_big * &(*R)));
            assert_eq!(big_to_6u64(&exp_mont), actual_mont.v);
        }
    }

    #[test]
    fn test_fe_mont_mul_raw() {
        let mut actual_mont = W6x64::default();
        let mut exp_mont = BigUint::default();
        let mut a_mont = W6x64::default();
        let mut b_mont = W6x64::default();

        for _i in 0..1_000_000 {
            let a_big = rnd_big_mod_n();
            let b_big = rnd_big_mod_n();
            fe_to_mont(&mut a_mont, &big_to_6u64(&a_big));
            fe_to_mont(&mut b_mont, &big_to_6u64(&b_big));

            fe_mont_mul_raw(&mut actual_mont, &a_mont, &b_mont);

            mont_mul_biguint(&mut exp_mont, &(a_big * &(*R)), &(b_big * &(*R)));
            assert_eq!(big_to_6u64(&exp_mont), actual_mont.v);
        }
    }

    #[test]
    fn test_fe_mont_mul_asm() {
        let mut actual_mont = W6x64::default();
        let mut exp_mont = BigUint::default();
        let mut a_mont = W6x64::default();
        let mut b_mont = W6x64::default();

        for _i in 0..1_000_000 {
            let a_big = rnd_big_mod_n();
            let b_big = rnd_big_mod_n();
            fe_to_mont(&mut a_mont, &big_to_6u64(&a_big));
            fe_to_mont(&mut b_mont, &big_to_6u64(&b_big));

            unsafe {
                fe_mont_mul_asm(&mut actual_mont.v[0], &a_mont.v[0], &b_mont.v[0]);
            }

            mont_mul_biguint(&mut exp_mont, &(a_big * &(*R)), &(b_big * &(*R)));
            assert_eq!(big_to_6u64(&exp_mont), actual_mont.v);
        }
    }
}
