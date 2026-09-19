# SageMath implementation of Δ-KLPT⁺ and Delfar

This repository is a proof-of-concept implementation for the paper
[Δ-KLPT⁺: Improved Norm Bounds for the Quaternion Isogeny Path Problem
and an Application to Ring Signatures](https://eprint.iacr.org/2026/XXXX).
This contains a SageMath implementation of a new KLPT-type algorithm Δ-KLPT⁺ and of Delfar,
a fully anonymous ring signature scheme based on it. It contains the three
components implemented for the paper:

- Δ-SigningKLPT⁺, the variant of Δ-KLPT⁺ used for signing, together with its
  subroutine `IdealNormReduce`;
- SQIsign based on Δ-SigningKLPT⁺;
- the ring signature scheme Delfar.

## Requirements

- SageMath 10.x (`sage` on your `PATH`). The timings in the paper were measured
  with SageMath 10.9.

No additional Python packages are required.

## Layout

| Path | Contents |
|------|----------|
| `sqisign.py` | `SQIsign(f, c, lam)`: SQIsign based on Δ-SigningKLPT⁺, with `Keygen`, `Sign`, `Verify` |
| `ring_sqisign.py` | `RingSQIsign(f, c, lam, n_parties)`: the ring signature scheme Delfar, with `Keygen`, `Sign`, `Verify` for a ring of `n_parties` public keys |
| `parameters.py` | Parameter sets `{f, c, lam}`; `p = c * 2^f - 1` |
| `quaternion.py`, `lattice.py` | Quaternion ideals, `IdealNormReduce`, Δ-SigningKLPT⁺ (`DeltaKLPT_plus`), lattice enumeration |
| `special_curve.py`, `montgomery.py`, `util.py` | Curve `E0`, ideal-to-isogeny, x-only arithmetic |
| `theta_structures/`, `theta_isogenies/`, `utilities/` | Dimension-2 theta isogeny machinery, third-party code (see below) |

## Usage

All entry points are `.sage` scripts. Run each with `sage <script>.sage` from
this directory, so that the `.py` modules are importable.

### `test_sqisign.sage`

Runs Keygen / Sign / Verify of SQIsign based on Δ-SigningKLPT⁺ 10 times for
every parameter set in `parameters.py` and asserts that all signatures verify.

```sh
sage test_sqisign.sage
```

### `test_ring_sqisign.sage`

Same as above for Delfar with a ring of 3 parties. For each trial a random
party is chosen as the signer.

```sh
sage test_ring_sqisign.sage
```

### `benchmark.sage`

Measures the average Keygen / Sign / Verify time of SQIsign based on
Δ-SigningKLPT⁺ (`base`) and of Delfar (`ring`) over all parameter sets.

```sh
sage benchmark.sage
sage benchmark.sage --num-trials 5 --tests base
sage benchmark.sage --num-trials 5 --num-parties 4 --tests ring
```

| Option | Default | Meaning |
|--------|---------|---------|
| `--num-trials N` | 10 | trials per parameter set |
| `--num-parties N` | 3 | ring size for the `ring` benchmark |
| `--tests {base,ring} ...` | both | which benchmarks to run: `base` is SQIsign based on Δ-SigningKLPT⁺, `ring` is Delfar |

Some Sage launchers pass a literal `--` to the script. If arguments are not
recognised, insert `--` before them:

```sh
sage benchmark.sage -- --num-trials 5
```

### `count_iter_INR.sage`

Counts the average number of `IdealNormReduce` calls in `DeltaKLPT_plus` for
each parameter set. Each option of `DeltaKLPT_plus` switches one of the four
optimizations of the paper, in the order in which they are listed there:

| Option | Optimization |
|--------|--------------|
| `initial_reduce` | reducing the common norm before the randomization phase |
| `shortest_alpha` | reducing the norm of the randomization quaternion γ |
| `search_prime` | enumerating candidates in the final `IdealNormReduce` call |
| `prime1mod4` | mitigating the condition on the final common norm N |

The script runs four settings: no optimization, the first two, the first three,
and all four. The first three are the default of `DeltaKLPT_plus` and the
setting used in Delfar. Takes `--num-trials` like `benchmark.sage`.

```sh
sage count_iter_INR.sage --num-trials 20
```

## Reproducing the tables of the paper

The tables of the paper report averages over 100 executions for the parameter
sets with f = 324, 500 and 664. The scripts loop over all parameter sets in
`parameters.py`, so they also print results for f = 248.

Running times of Delfar for a ring of n = 5 parties:

```sh
sage benchmark.sage --num-trials 100 --num-parties 5
```

This prints the times of both `base` and `ring`. In the paper, the KeyGen and
Verify times of `ring` are divided by n. The Sign time is reported as the `base`
signing time, which equals that of a ring with n = 1, plus n - 1 times the
simulator time. The simulator time is the difference between the `ring` and
`base` signing times divided by n - 1.

Numbers of `IdealNormReduce` calls:

```sh
sage count_iter_INR.sage --num-trials 100
```

Both runs take several hours: one Delfar signature takes from several seconds to
about a minute, depending on the parameter set.

## Third-party code

The directories `theta_structures/`, `theta_isogenies/` and `utilities/` are
unmodified copies of the corresponding directories of `Theta-SageMath` in
[ThetaIsogenies/two-isogenies](https://github.com/ThetaIsogenies/two-isogenies),
the code accompanying the paper
[An Algorithmic Approach to (2, 2)-isogenies in the Theta Model and Applications to Isogeny-based Cryptography](https://eprint.iacr.org/2023/1747)
by Pierrick Dartois, Luciano Maino, Giacomo Pope and Damien Robert.
That code is distributed under the MIT License
(Copyright (c) 2023 Pierrick Dartois, Luciano Maino, Giacomo Pope and Damien Robert);
a copy of the license is included as `LICENSE` in each of the three directories.

## Use of AI tools
A part of the code in this repository was generated with the assistance of Claude, an AI assistant developed by Anthropic.
In particular, the following files were almost entirely generated by Claude:

- `montgomery.py`: x-only arithmetic on Montgomery curves,
- `lattice.py`: lattice reduction and enumeration,
- `util.py`: utility functions, including a deterministic pseudorandom generator,
- `benchmark.sage` and `count_iter_INR.sage`: the measurement scripts used to
  produce the tables of the paper,
- `README.md`: the documentation of the repository.

The remaining files, which in particular contain the main algorithms
proposed in the paper (`IdealNormReduce`, Δ-KLPT⁺ and Δ-SigningKLPT⁺),
were written by the authors.
For these files, Claude was used for partial assistance, namely
debugging, performance optimization, and the implementation of some
individual functions.
