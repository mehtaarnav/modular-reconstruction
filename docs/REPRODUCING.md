# Reproducing the certified reconstruction

This guide concerns **one fixed instance**: PR #2, Trial 0. It does not change the input distribution, search for new parameter regimes, or reproduce an unproved worst-case decoding claim.

[Back to project overview](../README.md) · [Established results and certificates](../crt_codes/results/ESTABLISHED_RESULTS.md)

## Environment

The successful GitHub Actions jobs used Ubuntu, **Python 3.11**, and:

```bash
python -m pip install "sympy==1.14.0" "pytest==9.2.0" "fpylll==0.6.4" "cysignals==1.12.5"
```

`fpylll` requires native lattice libraries and may have different installation requirements outside Linux. If installation fails, use the archived GitHub Actions certificate links below; do not treat an unrun script as a successful reproduction.

## First, check the exact arithmetic

```bash
python -m pytest crt_codes/test_exhaustive.py crt_codes/test_integer_affine.py -q
```

These tests cover modular Vandermonde inversion, CRT-subset reconstruction, unimodular integer-kernel saturation, exact rational polytope bounds, and noncompletion behavior when a declared enumeration budget is exceeded.

## Reconstruct Trial 0

The fixed parameters are:

| Parameter | Value |
| :--- | :--- |
| Cyclotomic field | `Q(zeta_11)`, degree 10 |
| Full blocks | `n = 20` |
| Rational prime generator | `bits = 24`, code seed `2` |
| Coefficient bound | `|c_j| < 2^71` |
| Galois interpolation order | `s = 7` |
| Received-word seed | `20261121` |
| Selected instance | `trial = 0` |
| Required agreement | `>= 6` **whole blocks** |

### Scalar interpolation and same-word comparison

```bash
python -m crt_codes.scalar_affine_trial0 --output scalar-affine-trial0.json
```

The script deterministically constructs the same received word, builds the 11-dimensional weighted scalar congruence lattice, verifies each exact reduced relation and its height bound, and computes the complete integer solution list. It also reruns the Galois certificate on **that same received word** for an exact affine-space comparison. No planted message or reference candidate is used to choose interpolation relations.

Expected mathematical checks:

| Check | Expected result |
| :--- | :--- |
| Reduced scalar basis vectors examined | 11 |
| Scalar relations passing exact height test | 8 |
| Scalar rational equation rank / nullity | 8 / 2 |
| Scalar–Galois rational affine spaces | Equal |
| Certified integer-coordinate box | `[-1,0] × [-1,0]` |
| Box combinations / admissible integer messages | 4 / 3 |
| Messages satisfying six-block threshold | 3 |
| Independent complete reference count | 3 |

**Received-word fingerprint** (SHA-256 of canonical JSON blocks):

```text
cbcd857846e26de3c42e70498d04957c08ddbe38f7206adc92b6ff79eebf5fcb
```

The full output contains the eleven reduced relations, their exact integer coefficients and height values, the unimodular integer-kernel certificate, exact rational LP bounds, and stage timings.

### Galois Stage B certificate by itself

```bash
python -m crt_codes.stage_b_certificate --output stage-b-certificate.json
```

This independently rebuilds the six selected Galois interpolation relations and verifies the sufficient Stage A inequality:

```text
U_max < 2^117 < 2^138 <= P_min,6
```

It then certifies the saturated integer kernel, rational LP box, and complete three-message list. Stage-level timings are machine-dependent and are **not** a general speed comparison.

## Archived certificates

The original GitHub-hosted runs and their downloadable JSON artifacts are preserved:

| Run | Evidence |
| :--- | :--- |
| [Three-codeword adversarial audit](https://github.com/mehtaarnav/quantum_shenanigans/actions/runs/44875984991) | Empirical three-trial comparison and independent reference |
| [Instrumented Trial 0](https://github.com/mehtaarnav/quantum_shenanigans/actions/runs/44877406000) | Lattice construction, LLL, extraction, verification timings |
| [Exact Galois Stage B](https://github.com/mehtaarnav/quantum_shenanigans/actions/runs/44879104769) | Full integer matrix, unimodular certificate, LP vertices, exact list |
| [Exact scalar/Galois comparison](https://github.com/mehtaarnav/quantum_shenanigans/actions/runs/44881119774) | All scalar relations, same rational plane, same certified list |

GitHub may eventually expire workflow artifacts according to repository retention settings; the mathematical summaries remain in [`crt_codes/results/`](../crt_codes/results/).

## Boundaries of the result

The exhaustive CRT-subset decoder is a small-instance independent oracle **for correctness checks**, not an efficient general replacement for list decoding. The original Galois candidate extractor contains a limited-radius enumeration; the **separate Stage B certificate** establishes completeness for this input without relying on that cutoff.

The [established-results ledger](../crt_codes/results/ESTABLISHED_RESULTS.md) distinguishes the fixed-instance certificates from general theorems, inherited literature, historical experiments, and unresolved heuristic predictions.
