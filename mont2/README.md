# Optimizing Pairing-Based Cryptography: Montgomery Arithmetic in Rust


This subdirectory contains the working code corresponding to the blog
post covering xxxxx

> **Please note:** This code is for educational purposes, has not undergone
> a security audit and is not suitable for production. Use at your own risk.

* The blog post: <https://blog.link.here>

> Blog post description here

This code runs on Unbuntu, Mac OS and Windows. After installing Rust, run:

~~~
$ git clone https://github.com/nccgroup/pairing.git
$ cd pairing/mont2
$ cargo test
$ cargo bench
~~~

The arithmetic routines can be found in <https://github.com/nccgroup/pairing/blob/main/mont1/src/arith.rs>

The benchmarking code can be found in <https://github.com/nccgroup/pairing/blob/main/mont1/src/bench.rs>

---
Copyright 2021 Eric Schorn; Licensed under the 3-Clause BSD License.
