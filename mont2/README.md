# Optimizing Pairing-Based Cryptography: Montgomery Arithmetic in Rust


This subdirectory contains the working code corresponding to the blog
post covering Montgomery Multiplication in Assembly.

> **Please note:** This code is for educational purposes, has not undergone
> a security audit and is not suitable for production. Use at your own risk.

* The blog post: <https://research.nccgroup.com/2021/09/10/optimizing-pairing-based-cryptography-montgomery-multiplication-in-assembly/>

> This is the second blog post in a new code-centric series about selected 
> optimizations found in pairing-based cryptography. The first post in this 
> series covered modular Montgomery arithmetic in Rust from start to finish. 
> This second post takes the Montgomery multiplication algorithm developed 
> in Rust even further to seek the maximum performance a modern x86-64 machine 
> can deliver from an implementation hand-written in assembly language. 
> Several specialized instructions and advanced micro-architectural features 
> enabling increased parallelism result in the Montgomery multiplication 
> routine running more than 15X faster than a generic Big Integer implementation.

This code runs on Unbuntu, Mac OS and Windows. After installing Rust, git and clang
run:

~~~
$ git clone https://github.com/nccgroup/pairing.git
$ cd pairing/mont2
$ cargo test
$ RUSTFLAGS="--emit asm -C target-cpu=native" cargo bench
$ RUSTFLAGS="--emit asm -C target-feature=+bmi2" cargo bench
~~~

The arithmetic routines can be found in <https://github.com/nccgroup/pairing/blob/main/mont2/src/arith.rs>

The Montgomery multiplication assembly code can be found in <https://github.com/nccgroup/pairing/blob/main/mont2/src/mont_mul_asm.S>

The benchmarking code can be found in <https://github.com/nccgroup/pairing/blob/main/mont2/src/bench.rs>

Example results

~~~
Ubuntu 20.04.2 LTS on Intel Core i7-7700K CPU @ 4.20GHz with Rust version 1.53

Addition X 1000 iterations                                 [3.4908 us 3.4921 us 3.4934 us]  // native
Subtraction X 1000 iterations                              [3.1235 us 3.1238 us 3.1240 us]  // native
Multiplication by BigUint X 1000 iterations:               [493.13 us 493.15 us 493.17 us]  // bmi2
Multiplication in Rust (mont1 blog) X 1000 iterations:     [46.052 us 46.061 us 46.074 us]  // native
Multiplication in Rust with intrinsics X 1000 iterations:  [31.461 us 31.467 us 31.473 us]  // native
Multiplication in RAW Rust X 1000 iterations:              [32.039 us 32.060 us 32.084 us]  // bmi2
Multiplication in Rust with assembly X 1000 iterations:    [28.924 us 28.935 us 28.945 us]  // 'bmi2'
~~~

---
Copyright 2021 Eric Schorn; Licensed under the 3-Clause BSD License.
