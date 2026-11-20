# CRT codes: integer "multiplicity" and the Johnson-to-capacity gap

Alpöge's RH note (alpo.ge/remark.pdf, 2026-11-20) names the obstacle his proof
gets around: Stepanov's method needs *multiplicity*, and "multiplicity is hard to
find over the integers, corresponding to something like differentiating whole
numbers". His substitute is a chain of about √X theta lifts.

Coding theory has the same obstacle. For polynomials:

- Reed–Solomon codes list-decode up to the **Johnson bound** (Guruswami–Sudan).
- Multiplicity codes and folded RS codes, which carry derivative or twisted
  information per coordinate, reach **list-decoding capacity** (Guruswami–Wang,
  Kopparty, Guruswami–Rudra).

The integer analogue of Reed–Solomon is the **Chinese Remainder code**: the
message is m < K and the codeword is (m mod q_1, …, m mod q_n). With prime-power
moduli q_i = p_i^s, knowing m mod p^s is the integer counterpart of knowing
f, f′, …, f^(s−1) at a point. The best known CRT decoders (Boneh;
Guruswami–Sahai–Sudan) stop at the Johnson bound.

This directory is the measurement harness for that gap.

- `code.py`: CRT codes with prime or prime-power moduli, and the error channel.
- `decode.py`: the lattice list-decoder (Coppersmith/Howgrave-Graham "small root
  modulo an unknown divisor" with multiplicity h). This is the baseline.
- `experiment.py`: sweeps the error count against the unique, Johnson and
  capacity radii.
- `results/baseline.md`: tables from `python3 -m crt_codes.experiment --trials 12`.

## Baseline result (n = 24, k = 6, rate 1/4)

| Setting | Decodes reliably up to | Unique | Johnson | Capacity |
|---|---|---|---|---|
| Prime moduli, full errors (h = 4, 6) | 11 errors | 9 | 12 | 18 |
| Prime-power moduli p², full errors | 11 errors | 9 | 12 | 18 |
| p², errors corrupting only the high p-adic digit | 22 errors | – | – | – |

The lattice decoder reaches the Johnson bound. Prime-power ("integer
multiplicity") moduli give **no gain against full errors**. Derivatives close
exactly this gap for polynomial codes, and there is no derivative of an integer.
A decoder that beats Johnson here is the target. It needs an integer substitute
for differentiation, which is the role theta lifts play in Alpöge's argument.

## Run

```bash
pip install fpylll cysignals sympy
python3 -m pytest crt_codes -q
python3 -m crt_codes.experiment --trials 12
```

## Stage 2: Galois-folded cyclotomic codes (`folded.py`)

See `RESEARCH.md` for the literature, construction, heuristic radius
(s/(s+1))·(1 − R) and current results. On ℚ(ζ_11) with n = 20 and rate 0.156,
the decoder corrects 13 and 14 random errors against a Johnson radius of 12.09.
Worst-case two-codeword tests are in `beyond_johnson.py`.
