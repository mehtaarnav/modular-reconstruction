# Contributing

Corrections and small, reproducible improvements are welcome. This is a research repository: accuracy and a clear record of assumptions matter more than adding features quickly.

## Before opening a pull request

- **Code changes:** include a minimal reproduction, relevant tests, and the precise environment when behavior depends on lattice libraries.
- **Mathematical claims:** distinguish a general proof from a certificate for one input, an empirical observation, or a heuristic. State assumptions and link to the exact code or source used.
- **Experimental results:** record parameters, seeds, comparison criteria, and whether the decoder uses any information about the transmitted message.
- **Documentation:** keep the [established results ledger](crt_codes/results/ESTABLISHED_RESULTS.md) as the source of truth; exploratory notes should not silently become proved claims.

The current exact-arithmetic tests can be run with:

```bash
python -m pip install "sympy==1.14.0" "pytest==9.2.0"
python -m pytest crt_codes/test_exhaustive.py crt_codes/test_integer_affine.py -q
```

The [reproduction guide](docs/REPRODUCING.md) covers the pinned optional `fpylll` dependencies and the fixed certified Trial-0 instance.

When reporting an issue, include the observed and expected behavior, the full command, Python version, and relevant output. Please do not publish a decoder-radius or runtime claim based on a single successful reconstruction.
