# newKLPT — SQIsign and Ring SQIsign in SageMath

A SageMath implementation of the SQIsign signature scheme and a ring-signature
variant built on it. The quaternion side uses a variant of KLPT
(`DeltaKLPT_plus`) that fixes the norm of the response ideal; the isogeny side
uses x-only Montgomery arithmetic and 2-dimensional theta isogenies.

## Requirements

- SageMath 10.x (`sage` on your `PATH`)

No additional Python packages are required.

## Layout

| Path | Contents |
|------|----------|
| `sqisign.py` | `SQIsign(f, c, lam)` with `Keygen`, `Sign`, `Verify` |
| `ring_sqisign.py` | `RingSQIsign(f, c, lam, n_parties)` with the same API |
| `parameters.py` | Parameter sets `{f, c, lam}`; `p = c * 2^f - 1` |
| `quaternion.py`, `lattice.py` | Quaternion ideals, KLPT variants, lattice enumeration |
| `special_curve.py`, `montgomery.py`, `util.py` | Curve `E0`, ideal-to-isogeny, x-only arithmetic |
| `theta_structures/`, `theta_isogenies/`, `utilities/` | Dimension-2 theta isogeny machinery |

## Usage

All entry points are `.sage` scripts. Run each with `sage <script>.sage` from
the repository root, so that the `.py` modules are importable.

### `test_sqisign.sage`

Runs Keygen / Sign / Verify 10 times for every parameter set in
`parameters.py` and asserts that all signatures verify.

```sh
sage test_sqisign.sage
```

### `test_ring_sqisign.sage`

Same as above for the ring signature with 3 parties. For each trial a random
party is chosen as the signer.

```sh
sage test_ring_sqisign.sage
```

### `benchmark.sage`

Measures average Keygen / Sign / Verify time for SQIsign and Ring SQIsign over
all parameter sets.

```sh
sage benchmark.sage
sage benchmark.sage --num-trials 5 --tests base
sage benchmark.sage --num-trials 5 --num-parties 4 --tests ring
```

| Option | Default | Meaning |
|--------|---------|---------|
| `--num-trials N` | 10 | trials per parameter set |
| `--num-parties N` | 3 | ring size for the `ring` benchmark |
| `--tests {base,ring} ...` | both | which benchmarks to run |

Some Sage launchers pass a literal `--` to the script. If arguments are not
recognised, insert `--` before them:

```sh
sage benchmark.sage -- --num-trials 5
```

### `count_iter_INR.sage`

Counts the average number of iterations inside `DeltaKLPT_plus` for each
parameter set, under four combinations of the options `initial_reduce`,
`shortest_alpha`, `search_prime` and `prime1mod4`. Takes `--num-trials` like
`benchmark.sage`.

```sh
sage count_iter_INR.sage --num-trials 20
```
