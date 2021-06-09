# Optimizing Pairing-Based Cryptography: Montgomery Arithmetic in Rust


This subdirectory contains the working code corresponding to the blog
post covering Montgomery Arithmetic in Rust.

> **Please note:** This code is for educational purposes, has not undergone
> a security audit and is not suitable for production. Use at your own risk.

* The blog post: <https://research.nccgroup.com/2021/06/09/optimizing-pairing-based-cryptography-montgomery-arithmetic-in-rust/>

> This is the first code-centric blog post in a new series about selected 
> optimizations found in pairing-based cryptography. Pairing operations are 
> foundational to the BLS Signatures central to Ethereum 2.0, zero-knowledge 
> arguments central to Zcash and Filecoin, and a wide variety of other 
> emerging applications. A prior blog series implemented the entire pairing 
> operation from start to finish in roughly 200 lines of straightforward 
> Haskell. This series will dive deeply into individual optimizations that 
> drive raw performance in more operational systems. This post will cover 
> modular Montgomery arithmetic from start to finish, including context, 
> alternatives, theory and practical working code in Rust running **9X faster** 
> than a generic Big Integer implementation. The next blog post will further 
> optimize the (relatively) heavyweight multiplication routine in bare-metal 
> x86-64 assembly language.

This code runs on Unbuntu, Mac OS and Windows. After installing Rust, run:

~~~
$ git clone https://github.com/nccgroup/pairing.git
$ cd pairing/mont1
$ cargo test
$ cargo bench
~~~

The arithmetic routines can be found in <https://github.com/nccgroup/pairing/blob/main/mont1/src/arith.rs>

The benchmarking code can be found in <https://github.com/nccgroup/pairing/blob/main/mont1/src/bench.rs>

---
Copyright 2021 Eric Schorn; Licensed under the 3-Clause BSD License.
