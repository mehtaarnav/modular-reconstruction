# Integer codes past the Johnson bound via Galois folding

## The question

For polynomials, list decoding up to capacity (1 − R) is now essentially solved:

- folded Reed–Solomon codes and multiplicity codes (Guruswami–Rudra; Guruswami–Wang; Kopparty);
- folded codes from cyclotomic function fields, which use Artin–Frobenius automorphisms
  (Guruswami, arXiv:0811.4139, ANT 2010), and from function-field towers (Guruswami–Xing);
- plain Reed–Solomon codes over prime fields on *every* evaluation set
  (Brakensiek–Chen–Putterman–Zhang–Zheng, arXiv:2609.08005, Sept 2026). They add Hasse
  derivatives of the unknown message as phantom interpolation variables, tied down by a
  local Taylor identity.

For **integers**, the analogue of Reed–Solomon is the Chinese Remainder code: m ↦ (m mod p_i).
The best known worst-case decoders (Boneh; Guruswami–Sahai–Sudan 2000) stop at the Johnson
bound. Every capacity-achieving polynomial construction relies on derivatives or on a global
twist X ↦ γX, and integers have neither. Alpöge's RH note (2026-11-20) identifies this as the
obstacle in his own setting: multiplicity over ℤ "corresponds to something like
differentiating whole numbers".

Related integer work found:

- Interleaved CRT codes decoded beyond half the minimum distance with LLL, for random errors
  (Li–Sidorenko–Nielsen, ISIT 2013), with a probabilistic failure analysis (Abbondati et al.,
  2023; 2025).
- Number-field codes (Lenstra; Guruswami 2003), with a list decoder by Biasse–Quintin
  (arXiv:1107.2321).
- Guruswami's cyclotomic-function-field paper names the number-field analogue ℚ(ω) and
  leaves it open.

None of these give a worst-case decoder for an integer code family beyond Johnson.

## Construction: Galois-folded cyclotomic codes

Take K = ℚ(ζ_M) with M prime, φ = M − 1, and σ: ζ ↦ ζ^g for a generator g of (ℤ/M)^×.
Then Gal(K/ℚ) = ⟨σ⟩, and the fixed field of σ is ℚ.

- **Message:** α ∈ ℤ[ζ] with integer coefficients |c_t| < H.
- **Coordinates:** primes p_1, …, p_n ≡ 1 (mod M), each splitting into φ degree-1 primes, which
  form one σ-orbit.
- **Symbol at p:** y_j = α(ω^{g^j}) mod p for j = 0..φ−1, where ω is a primitive M-th root
  of unity mod p. This is α reduced at every prime above p.

Since σ^l(α) mod 𝔓_j = y_{j+l}, each block is a folded-RS block whose *entire* orbit is used,
with cyclic windows. Polynomials cannot do this: F_q^× is a single orbit of X ↦ γX, so FRS
must cut that orbit into blocks and pays the factor m/(m − s + 1). Number fields provide one
full orbit per rational prime.

(Up to a change of basis in each block, this is an interleaved CRT code on the φ coefficients
of α. The difference from Li–Sidorenko–Nielsen is the decoder, which uses the σ-structure,
and the error model, which is worst-case.)

## Decoder (Guruswami–Wang linear-algebraic template, with heights in place of degrees)

1. **Interpolation, one LLL call.** Find short A_0, …, A_s ∈ ℤ[ζ] with, for every block i and
   every window j,

   A_0(x) + Σ_{l=1..s} A_l(x) · y_{i, j+l−1} ≡ 0 (mod p_i),  x = ω_i^{g^j}.

   This is a lattice of index ∏ p_i^φ in ℤ^{(s+1)φ}. A_0 gets the extra height allowance H·φ·s.
2. **Vanishing.** E = A_0 + Σ A_l σ^{l−1}(α) lies in every prime above every agreeing block, so
   N(E) is divisible by ∏_{agree} p_i^φ. If that product exceeds the archimedean size of E,
   then E = 0.
3. **σ-linear algebra.** E = 0 is a σ-difference equation of order s − 1 over K, whose constants
   are ℚ. Its solutions form an affine ℚ-space of dimension ≤ s − 1 (Casoratian argument).
   Every sufficiently short lattice vector gives such an equation, so intersect them.
4. **List recovery.** Enumerate small integer points of the remaining affine space with one
   embedding LLL plus Euclid on the inhomogeneous column, then filter by agreement.

**Heuristic radius.** Minkowski plus the weighting gives E = 0 whenever the agreement fraction
a satisfies

  a > 1/(s+1) + R·s/(s+1),  i.e.  errors < (s/(s+1))·(1 − R),

which tends to capacity 1 − R as s → φ grows. Johnson is 1 − √R. LLL loses a factor that
grows with the dimension (s+1)φ; in practice it is small here.

## Results so far (`folded.py`, `beyond_johnson.py`)

| Field | n | rate R | s | Johnson radius (errors) | predicted (s/(s+1))(1−R)·n | observed (random errors, all trials) |
|---|---|---|---|---|---|---|
| ℚ(ζ_5), φ=4 | 16 | 0.138 | 3 | < 10.06 | 10.3 | 10 (at Johnson, **not** past it) |
| ℚ(ζ_11), φ=10 | 20 | 0.156 | 7 | < 12.1 | 14.8 | **14** (13 and 14 ok; see `results/beyond_johnson.txt` for 15 and 16) |

**Worst-case two-codeword mixtures** at ℚ(ζ_5): received words that agree with α on a blocks
and with β on b blocks. For (a, b) ∈ {(8,8), (7,9), (7,7), (6,10), (6,6)}, both α and β
satisfied the interpolated σ-linear system in every trial, and both were recovered in 30 of
30 trials. All these agreements are at or above the Johnson threshold (5.95), so this checks
list recovery, not radius. The ℚ(ζ_11) mixtures with agreements *below* the Johnson threshold are in
`results/beyond_johnson.txt`.

## What is and isn't established

- **Heuristic, not proved:** the radius uses Minkowski's bound for an LLL-found vector. A proof
  needs either an explicit LLL-loss accounting or a Siegel/Bost-slope argument for existence.
  The latter is the Arakelov tool that also appears in Alpöge's proof (note, step 13).
- **List size is the open part.** The output is an affine ℚ-space of dimension ≤ s − 1. In the
  experiments, intersecting several short vectors cut it to a point, or to a line when two
  codewords both agreed heavily. A worst-case polynomial list bound would need pruning, in the
  role subspace-evasive sets play for folded RS. That is the next theoretical target.
- **Novelty is checked only by search.** Searches covered Guruswami (function fields),
  Li–Sidorenko–Nielsen (interleaved CRT, random errors), Biasse–Quintin (number-field codes),
  BCPZZ (RS capacity), and general queries for folded or capacity results for CRT codes and
  number-field codes. They found no worst-case decoder past Johnson for an integer code family.
  This needs a check by people in the area before any claim is made publicly.

## Next steps

1. Scale up: ℚ(ζ_23) (φ = 22) with s up to about 20, to show the radius approaching capacity.
2. Adversarial stress tests beyond two-codeword mixtures, e.g. many codewords, and errors built
   to enlarge the kernel.
3. Prove the interpolation step with an explicit Siegel lemma (or Bost slope inequality) for
   the weighted lattice.
4. A list-size bound: pruning with the Frobenius at an auxiliary prime ℓ (σ ≡ ℓ-th power mod ℓ
   for σ = Frob_ℓ), the number-field counterpart of Guruswami–Rudra's root-finding over an
   extension field.
5. Plain CRT codes: can the phantom-derivative interpolation of arXiv:2609.08005 be given an
   integer counterpart? This is the version of the problem stated in Alpöge's note.
