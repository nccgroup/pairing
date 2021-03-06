// Copyright 2021 Eric Schorn; Licensed under the 3-Clause BSD License.

use criterion::{criterion_group, criterion_main, Criterion};
use mont2::arith::{fe_add, fe_mont_mul, fe_mont_mul_intrinsics, fe_mont_mul_raw, fe_sub, W6x64};
use mont2::fe_mont_mul_asm;
use num_bigint::BigUint;
use num_traits::Num;
use std::time::Duration;

#[macro_use]
extern crate lazy_static;

// RUSTFLAGS="--emit asm -C target-cpu=native" cargo bench
// RUSTFLAGS="--emit asm -C target-feature=+bmi2" cargo bench
// see, e.g.: mont2/target/release/deps/mont2-eb826ed791da90e8.s

// Arbitrary input and expected values for 1000 iterations to ensure functionality
#[rustfmt::skip]
const X: W6x64 = W6x64 {
    v: [0xc34110121829fa85, 0xc42f61586f13abac, 0x5a98f20b2164430a,
        0xcdd6beb839ca6556, 0xdacae65ae941e8e8, 0xf594a44cbdf0ae1]
};

#[rustfmt::skip]
const Y: W6x64 = W6x64 {
    v: [0xf4921aadbbf08d96, 0x9f5973902a56b682, 0x4b86761f89b618b2,
        0xca440e25b9c201dd, 0xd3caeb49dc668726, 0x416ce3c635e5e23],
};

#[rustfmt::skip]
const EXP_SUM: W6x64 = W6x64 {
    v: [0xbd3d31dc0303fa06, 0x704875edd38742a3, 0x60549f5927a1c745,
        0xe234ee37eb7d3cee, 0xa13832d81ab0d5c5, 0x10fdd2f03da8f7ca],
};

#[rustfmt::skip]
const EXP_DIFF: W6x64 = W6x64 {
    v: [0xeb500a9ba3c63dbc, 0xf9d612366c970ad5, 0x581e56b55f02cbcb,
        0x60e49af2737caf46, 0x441baca536704b15, 0xebe95e1d0ff39dc],
};

#[rustfmt::skip]
const EXP_PROD: W6x64 = W6x64 {
    v: [0xb54cf29498954919, 0x8f2491ddb5cef751, 0xb155fe8acce5c7d3,
        0x448683648418e8dd, 0xf3599187e803fc7e, 0x1118bd439ac24052],
};

lazy_static! { static ref EXPECTED: BigUint = BigUint::from_str_radix(
    "169d18ab74c03e6199a9ec1869d2a2a0d53be1749c6acd5028310a17f06383087d69cb203aa01ae0a73a546f5db98555",
    16).unwrap();
}

// Montgomery addition x1000 written in Rust
fn add_rust(x: &W6x64, y: &W6x64, expected: &W6x64) {
    let mut xx = x.clone();
    let mut yy = y.clone();
    let mut result = W6x64::default();
    for _i in 0..1_000 {
        fe_add(&mut result, &xx, &yy);
        yy = xx;
        xx = result;
    }
    assert_eq!(&result, expected);
}

// Montgomery subtraction x1000 written in Rust
fn sub_rust(x: &W6x64, y: &W6x64, expected: &W6x64) {
    let mut xx = x.clone();
    let mut yy = y.clone();
    let mut result = W6x64::default();
    for _i in 0..1_000 {
        fe_sub(&mut result, &xx, &yy);
        yy = xx;
        xx = result;
    }
    assert_eq!(&result, expected);
}

lazy_static! { static ref MODULUS: BigUint = BigUint::from_str_radix(
    "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab",
    16).unwrap();
}

// Reference multiplication with BigUint
fn mul_biguint(result: &mut BigUint, a: &BigUint, b: &BigUint) {
    *result = (a * b) % &(*MODULUS);
}

// Multiplication x1000 with BigUint
fn mul_big(x: &BigUint, y: &BigUint, expected: &BigUint) {
    let mut xx = x.clone();
    let mut yy = y.clone();
    for _i in 0..1_000 {
        let mut result = BigUint::default();
        //let result = (&xx * &yy) % modulus;
        mul_biguint(&mut result, &xx, &yy);
        yy = xx;
        xx = result;
    }
    assert_eq!(&xx, expected)
}

// Montgomery multiplication x1000 written in Rust (mont1 blog)
fn mul_rust(x: &W6x64, y: &W6x64, expected: &W6x64) {
    let mut xx = x.clone();
    let mut yy = y.clone();
    let mut result = W6x64::default();
    for _i in 0..1_000 {
        fe_mont_mul(&mut result, &xx, &yy);
        yy = xx;
        xx = result;
    }
    assert_eq!(&result, expected);
}

// Montgomery multiplication x1000 written in Assembly
fn mul_asm(x: &W6x64, y: &W6x64, expected: &W6x64) {
    let mut xx = x.clone();
    let mut yy = y.clone();
    let mut result = W6x64::default();
    for _i in 0..1_000 {
        unsafe {
            fe_mont_mul_asm(&mut result.v[0], &xx.v[0], &yy.v[0]);
        }
        yy = xx;
        xx = result;
    }
    assert_eq!(&result, expected);
}

// Montgomery multiplication x1000 written in Rust with intrinsics
fn mul_rust_intrinsics(x: &W6x64, y: &W6x64, expected: &W6x64) {
    let mut xx = x.clone();
    let mut yy = y.clone();
    let mut result = W6x64::default();
    for _i in 0..1_000 {
        fe_mont_mul_intrinsics(&mut result, &xx, &yy);
        yy = xx;
        xx = result;
    }
    assert_eq!(&result, expected);
}

// Montgomery multiplication x1000 written in Rust (backported asm)
fn mul_rust_raw(x: &W6x64, y: &W6x64, expected: &W6x64) {
    let mut xx = x.clone();
    let mut yy = y.clone();
    let mut result = W6x64::default();
    for _i in 0..1_000 {
        fe_mont_mul_raw(&mut result, &xx, &yy);
        yy = xx;
        xx = result;
    }
    assert_eq!(&result, expected);
}

// Harness for addition with inputs and expected result
pub fn bench_add(c: &mut Criterion) {
    c.bench_function("1. Addition X 1000 iterations", |b| b.iter(|| add_rust(&X, &Y, &EXP_SUM)));
}

// Harness for subtraction with inputs and expected result
pub fn bench_sub(c: &mut Criterion) {
    c.bench_function("2. Subtraction X 1000 iterations", |b| {
        b.iter(|| sub_rust(&X, &Y, &EXP_DIFF))
    });
}

// Harness for multiplication by BigUint with inputs and expected result
pub fn bench_mul_big(c: &mut Criterion) {
    let x = BigUint::from(u128::MAX);
    let y = BigUint::from(u64::MAX);
    c.bench_function("3. Multiplication by BigUint X 1000 iterations", |b| {
        b.iter(|| mul_big(&x, &y, &(*EXPECTED)))
    });
}

pub fn bench_mul_rust(c: &mut Criterion) {
    c.bench_function("4. Multiplication in Rust (mont1 blog) X 1000 iterations", |b| {
        b.iter(|| mul_rust(&X, &Y, &EXP_PROD))
    });
}

pub fn bench_mul_rust_raw(c: &mut Criterion) {
    c.bench_function("5. Multiplication in flat Rust X 1000 iterations", |b| {
        b.iter(|| mul_rust_raw(&X, &Y, &EXP_PROD))
    });
}

pub fn bench_mul_rust_intrinsics(c: &mut Criterion) {
    c.bench_function("6. Multiplication in Rust with intrinsics X 1000 iterations", |b| {
        b.iter(|| mul_rust_intrinsics(&X, &Y, &EXP_PROD))
    });
}

pub fn bench_mul_asm(c: &mut Criterion) {
    c.bench_function("7. Multiplication in Rust with assembly X 1000 iterations", |b| {
        b.iter(|| mul_asm(&X, &Y, &EXP_PROD))
    });
}

// Run all seven harnesses
criterion_group! {
    name = benches;
    config = Criterion::default().measurement_time(Duration::new(60, 0));
    targets = bench_add, bench_sub, bench_mul_big, bench_mul_rust, bench_mul_rust_raw,
    bench_mul_rust_intrinsics, bench_mul_asm
}
criterion_main!(benches);
