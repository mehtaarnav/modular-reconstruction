# Three-codeword adversarial audit (reproducible protocol)

## Question

Does the Galois-folded lattice decoder recover **all three** independently
chosen messages when the received word consists of six full blocks from
each message and two unrelated blocks? At M=11, phi=10, n=20, s=7, the
project's calculated Johnson agreement threshold is about 7.9 blocks;
each planted message has at least six matching blocks.

This is a **stress test**, not a proof of a worst-case beyond-Johnson list
decoding radius. Even a positive outcome on several seeds would not prove
one. Nor is a full-block invertible Vandermonde transform a new code family
relative to interleaved CRT codes under the block-error metric. The possible
contribution is in the algorithm and a justified decoding guarantee.

## Run

From the repo root:

```sh
python -m pip install sympy pytest
python -m pytest crt_codes/test_exhaustive.py -q
python -m crt_codes.triad_audit --trials 3 --skip-lattice --output crt_codes/results/triad-reference.json

# Also evaluate the lattice candidate generator:
python -m pip install fpylll cysignals
python -m crt_codes.triad_audit --trials 3 --output crt_codes/results/triad-full.json
```

The first two commands require no fpylll; `folded.py` lazily imports it only
when its lattice decoder is called. `triad_audit.py` always prints JSON.
Results files are intentionally **not pre-populated**: do not confuse a
proposed stress test with empirical output. CI artifacts carry actual logs.

## Why the reference decoder is complete

For a rational prime p congruent to 1 mod M, let x_0,...,x_(phi-1)
be the *distinct* M-th roots in its full orbit. The full block satisfies
y_j = sum_(t=0)^(phi-1) c_t x_j^t (mod p). The Vandermonde determinant
prod_(i<j)(x_j-x_i) is nonzero modulo p, so Gaussian elimination recovers
all coefficient residues c_t mod p.

Let H=2^coeff_bits be the strict coefficient bound |c_t|<H. Choose the
smallest k such that the product of the **smallest k** block primes is
strictly greater than 2(H-1). This ensures that every k-subset has
a large enough CRT modulus to recover each signed coefficient uniquely.
Enumerating every k-subset, reconstructing every coefficient by CRT,
keeping only candidates in (-H,H), and re-encoding every candidate is a
**complete reference list** for all codewords with at least k agreeing blocks.

The audit refuses to claim completeness when the requested agreement
threshold is below k. Its worst-case cost is combinatorial in n and k.
It must NOT be used to claim an efficient, polynomial-time decoder.

## Output interpretation

- `reference.planted_found`: must equal 3. This checks the reference and
  code construction, not the proposed decoder.
- `reference.list_size`: number of *all* messages with >=6 whole-block
  agreements, not just the three planted ones.
- `lattice.planted_found`: of the three planted messages, how many survive
  the current `decode` implementation and a separate exact re-encoding
  filter at >=6 agreeing blocks.
- `lattice.missed_from_reference`: complete-list candidates that the
  lattice decoder failed to return.
- `lattice.spurious_vs_reference`: should be zero after verification.

This protocol gives the lattice decoder no hidden original message or
oracle-computed agreement mass. The old CRT baseline in
`experiment.py` *does* compute `log_b` from the transmitted message;
its measurements should be labeled **oracle-aided**. A fixed-input
interleaved CRT baseline remains necessary before making novelty claims.

## Caveats

1. The saved 7/7 experiment used only three trials. It is not a theorem.
2. The heuristic radius is (s/(s+1))*(1-R) as in RESEARCH.md; the earlier
   `folded.py` docstring used a different expression and has been aligned.
3. The code's `_integer_points` routine enumerates a radius-2 box in an
   LLL-reduced kernel basis. Failure to recover a valid codeword can result
   from this *enumeration cap*, not failure of interpolation or decoding in
   principle. Record and analyze that separately.
4. Lattice `decode` discards inconsistent short-vector equations. That
   does not certify the surviving system contains all valid codewords.
5. A promising further experiment is adversarially selecting received
   blocks to maximize the nullity of the sigma-linear interpolation system.
   Do not scale to larger M until the small case passes a careful audit.
