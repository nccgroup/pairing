# Optimizing Pairing-Based Cryptography: Montgomery Arithmetic in Rust

[![cargo test mont1](https://github.com/integritychain/pairing/actions/workflows/mont1.yml/badge.svg)](https://github.com/integritychain/pairing/actions/workflows/mont1.yml)

This subdirectory contains the working code corresponding to the blog
post covering [Montgomery Arithmetic in
Rust](https://research.nccgroup.com/2020/08/13/pairing-over-bls12-381-part-3-pairing/).

> This is the first code-centric blog post in a new series about selected optimizations found in pairing-based cryptography. Pairing operations are foundational to the BLS Signatures[^1] central to Ethereum 2.0, zero-knowledge arguments central to Zcash and Filecoin[^2], and a wide variety of other emerging applications. A prior blog series implemented the entire pairing operation from start to finish in roughly 200 lines of straightforward Haskell. This series will dive deeply into individual optimizations that drive raw performance in more operational systems. This post will cover modular Montgomery arithmetic[^3] from start to finish, including context, alternatives, theory and practical working code in Rust running **9X faster** than an generic Big Integer implementation. The next blog post will further optimize the (relatively) heavyweight multiplication routine in bare-metal x86-64 assembly language.

This code runs on Unbuntu, Mac OS and Windows. To install and run:

~~~
$ git clone https://github.com/eschorn1/ff_12381.git
$ cd pairing/mont1
$ cargo test
$ cargo bench
~~~

Copyright 2021 Eric Schorn; Licensed under the 3-Clause BSD License.
