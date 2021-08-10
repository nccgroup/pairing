// Copyright 2021 Eric Schorn; Licensed under the 3-Clause BSD License.

extern crate cc;

fn main() {
    cc::Build::new().file("src/mont_mul_asm.S").compile("mont_mul_asm");
}
