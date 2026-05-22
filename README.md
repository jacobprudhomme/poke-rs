# poke-rs

Rust implementations of POKÉ and INKE, compact and efficient public-key encryption schemes from higher-dimensional isogenies. These are meant to be cryptographic-quality reference implementations of the algorithms, both to give a real-world-ish estimate of their efficiency, as well as to give people a clean, correct and documented implementation they can refer to. However&mdash;and this should go without saying&mdash;**DO NOT USE THEM IN A PRODUCTION ENVIRONMENT**!

These were implemented as one part of my master's thesis within [LASEC](https://lasec.epfl.ch) at EPFL: "Gotta Study 'Em All: An efficiency and security analysis of the POKÉ family of PKEs from higher-dimensional isogenies". This project was co-supervised by Dr. Andrea Basso at IBM Research Zürich and Prof. Péter Kutas at Eötvös Loránd University (ELTE).

This implementation relies heavily on the great libraries [fp2](https://github.com/GiacomoPope/fp2) and [isogeny_rs](https://github.com/GiacomoPope/isogeny_rs) by Dr. Giacomo Pope, for the former's finite field arithmetic, and the latter's 1D/2D-isogeny computations and other useful functionality (big number arithmetic, pairings, discrete logarithms, etc.).

Source literature:

- [POKÉ: A Compact and Efficient PKE from Higher-dimensional Isogenies](https://eprint.iacr.org/2024/624)
- [INKE: Fast Isogeny-Based PKE Using Intermediate Curves](https://eprint.iacr.org/2025/1458)

In the future, I plan to also include an implementation of PIKE:

- [PIKE: Faster Isogeny-Based Public Key Encryption with Pairing-Assisted Decryption](https://eprint.iacr.org/2026/473)

## Using the library

## Examples

Examples can be run using the following command:

```bash
cargo run --example {poke,inke}_{i,iii,v}
```

The trailing roman numeral indicates the NIST security level to run the protocol at.

## Tests

For the moment, the only functional tests of public interest are the `encryption` ones. All the other ones are meant to test internal choices for how to implement things, comparing the different methods to make sure they act functionally the same.

To run them, execute

```bash
cargo test --test encryption
```

## Benchmarks

Similar to the tests, many of the exposed benchmarks are not really of public interest and exist to compare different methods of implementation. The ones of interest are `poke` and `inke` (as you could have probably guessed).

To run them, execute

```bash
cargo bench --bench {poke,inke}
```

### Results

All benchmarks were collected on an 2021 M1 MacBook Pro, 14", at 1000 iterations. For Rust, I used Criterion. For C code coming from other implementations, I used their comparable standard benchmarking code.

| Key Generation | POKÉ (ms) | INKE (ms) |
|----------------|-----------|-----------|
| Level I        | 71.181    | 63.450    |
| Level III      | 288.57    | 216.88    |
| Level V        | 705.21    | 570.59    |

| Encryption | POKÉ (ms) | INKE (ms) |
|------------|-----------|-----------|
| Level I    | 8.8790    | 10.152    |
| Level III  | 32.415    | 30.793    |
| Level V    | 71.928    | 73.465    |

| Decryption | POKÉ (ms) | INKE (ms) |
|------------|-----------|-----------|
| Level I    | 11.889    | 5.0418    |
| Level III  | 43.272    | 15.204    |
| Level V    | 94.475    | 35.880    |

## TODO

- [ ] Replace modular reduction from [rug](https://docs.rs/rug/latest/rug) with my own constant-time implementation (this is the only non-constant-time part remaining that can currently be made constant-time)
- [ ] Add complete test suite
- [ ] Integrate reusable parts of library into [isogeny_rs](https://github.com/GiacomoPope/isogeny_rs)
- [ ] Split INKE out into its own package
- [ ] Publish crate
- [ ] Experiment with other implementation options. Such as...
  - [ ] Evaluating 2D-isogenies using composition with dual (i.e. no discrete logarithms)
  - [ ] Optimizing implementation of discrete logarithm-solving
  - [ ] Other things when I remember them (I should have written a list...)
