# Optimizations for Pairing-Based Cryptography

This repository contains the working Rust code that corresponds to an
ongoing series of blog posts covering optimizations in pairing-based
cryptography using BLS12-381 as the reference curve.

> **Please note:** This code is for educational purposes, has not undergone
> a security audit and is not suitable for production. Use at your own risk.


## Montgomery arithmetic (1 of 2)

This post will cover modular Montgomery arithmetic from start to
finish, including context, alternatives, theory and practical working
code in Rust running **9X faster** than a generic Big Integer
implementation. The next blog post will further optimize the
(relatively) heavyweight multiplication routine in bare-metal x86-64
assembly language.

* The blog post: <https://research.nccgroup.com/2021/06/09/optimizing-pairing-based-cryptography-montgomery-arithmetic-in-rust/>
* The full code and README is in the `mont1` subdirectory: <https://github.com/nccgroup/pairing/tree/main/mont1>


## Montgomery arithmetic (2 of 2)

This second post takes the Montgomery multiplication algorithm developed
in Rust even further to seek the maximum performance a modern x86-64 machine
can deliver from an implementation hand-written in assembly language.
Several specialized instructions and advanced micro-architectural features
enabling increased parallelism result in the Montgomery multiplication
routine running more than 15X faster than a generic Big Integer implementation.

* The blog post: <https://coming.soon/>
* The full code and README is in the `mont2` subdirectory: <https://github.com/nccgroup/pairing/tree/main/mont2>

---

Copyright 2021 Eric Schorn; Licensed under the 3-Clause BSD License.
