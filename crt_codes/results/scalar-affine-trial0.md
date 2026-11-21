# Fixed Trial-0: 11-dimensional scalar affine relation vs Galois interpolation

**Finding:** the scalar comparator produces **eight** height-certified independent affine relations and the **same rational affine solution space** as the six Galois interpolation relations. Its certified integer extraction finds **exactly three** admissible messages, all with at least six full-block agreements; the list equals the independently complete reference. This falsifies **Galois-linearity being necessary to recover this particular received word**. No general speedup or radius theorem follows.

**Provenance:** same PR #2 Trial 0, `M=11,n=20,bits=24,coeff_bits=71,H=2^71,code_seed=2,received_seed=20261121,agreement>=6`. SHA256 of canonical JSON received word: `cbcd857846e26de3c42e70498d04957c08ddbe38f7206adc92b6ff79eebf5fcb`. One CI run: [GitHub Actions 44881119774](https://github.com/mehtaarnav/quantum_shenanigans/actions/runs/44881119774). **Full exact machine-readable JSON certificate**: [workflow artifact 14697369959](https://github.com/mehtaarnav/quantum_shenanigans/actions/runs/44881119774/artifacts/14697369959). The artifact includes all reduced vectors, all integer coefficients, exact U values, CRT residues, unimodular transformation, integer solution basis, rational LP extrema and timing diagnostics.

**Source of formulas:** [`crt_codes/scalar_affine_trial0.py`](https://github.com/mehtaarnav/quantum_shenanigans/blob/audit/pr2-scalar-affine-comparator/crt_codes/scalar_affine_trial0.py). Source of original code: [`crt_codes/folded.py:39–76,118–190`](https://github.com/mehtaarnav/quantum_shenanigans/blob/5ae11eb6bd0bf5ba35b887d672ced724c9ce5c18/crt_codes/folded.py#L118-L190). Certified integer kernel/extraction functions: [`crt_codes/integer_affine.py`](https://github.com/mehtaarnav/quantum_shenanigans/blob/audit/pr2-scalar-affine-comparator/crt_codes/integer_affine.py).

## 1. Exact 11-dimensional congruence lattice

For every rational prime `p_i`, invert the ten-symbol Vandermonde block over `F_(p_i)` to get residue `r_(i,j)` of message coefficient `c_j`. CRT the `j`th stream across all 20 primes to get `R_j (mod N)`, where `N=prod_i p_i` (461 bits).

Search for a scalar *affine functional* on message coefficient vectors:

```
       a0 + sum_{j=0}^9 aj*cj = 0  over Z.
```

The modular interpolation condition is `a0 + sum_j aj*Rj == 0 (mod N)`. The 11-dimensional weighted congruence lattice is exactly

```
L = {(a0, H*a1, ..., H*a10) : a0 + sum_j aj*Rj == 0 (mod N)}
```

with explicit integral row basis `(N,0,...,0)` and `(-Rj,0,...,H at position j,...,0)`. Its determinant is `N*H^10`. Apply the default fpylll LLL to this **11×11** basis, unweight each returned basis row *exactly*, verify `N` congruence and each of the 20 individual prime congruences, and compute

```
       U = |a0| + (H-1)*sum_j |aj|.
```

The exact comparison threshold is

```
P_min6 = 412140041510669591866160348537616348856951.
```

If `U<P_min6`, every admissible message agreeing on six full blocks satisfies the integer equality, because the left side is divisible by the six agreement primes and bounded strictly below their product. No planted messages or reference candidates were consulted when choosing relations.

## 2. All eleven reduced basis vectors, exact height classification

| Reduced row | Retained? | Exact U |
|---:|:---:|---:|
| 0 | Yes | 75737313737571461262952966364 |
| 1 | Yes | 105344766295861970673791493629 |
| 2 | Yes | 111719813862786044752287328928 |
| 3 | Yes | 109193419286162267741582901833 |
| 4 | Yes | 117115156764528948107621828783 |
| 5 | Yes | 115447327651099055906245071606 |
| 6 | Yes | 112527321774211299815502376078 |
| 7 | Yes | 114185005432957118482610275013 |
| 8 | No | 602736970313061231885698060724854890261617 |
| 9 | No | 551320676588714646841471372290303390266200 |
| 10 | No | 703972933385729269880360039709346253673941 |

The largest **retained** U was `117115156764528948107621828783`, substantially smaller than `P_min6`. Each of the eleven reduced rows passed all congruence checks. The last three were rejected **solely on the explicit height inequality**.

The eight selected **exact integer affine relation coefficients** `(a0 ; [a1,...,a10])` are:

```text
0:  1131722296477526105658182005 ; [4537117,-2027036,3262056,-501425,-3341868,131354,-2403168,-4196283,3512174,-7684216]
1: -3341657349427358362468914170 ; [-6934922,-4220950,-1780707,-3166966,-5608224,1201109,-5193641,7173133,6176224,-1744121]
2: 10914519007631869340482591487 ; [5003847,-2423991,3400405,-4217625,5378612,-2994935,3777455,7262989,-7937388,295456]
3: -1106896445730701160591933552 ; [7052354,-1363494,6534264,6397726,-4191860,-3189041,-10273843,-2190823,-693150,3889868]
4: -11215349335822585091057935672 ; [-11353450,2607916,6999218,1115957,8528173,910213,5380835,-2987053,1704517,-3262981]
5: 11162944987669499689559015168 ; [-850658,-4356948,7225760,2795065,812575,-3332885,5196820,7796751,11325776,-472916]
6: -3219080223363010190921846751 ; [-259734,4851935,6180219,-10610891,554364,-13000650,1924004,-3910920,-4446454,-554670]
7: -1850175591405258892347192616 ; [3833279,11471234,-778575,303383,4541705,-3891528,-5058151,3017486,14243647,-436663]
```

## 3. Rational rank and affine solution space

Eight certified scalar relations have coefficient matrix rank **8** over Q, leaving rational affine dimension **2** in ten message coordinates. Rebuild the six Galois sigma-linear relations independently on the *identical received word* and test equality by exact rational matrix ranks:

| Check | Result |
|---|---:|
| Scalar equation rank | 8 |
| Galois equation rank | 8 |
| Rank of all scalar + Galois equations stacked | 8 |
| Rank of stacked augmented matrix | 8 |
| Both individual affine systems consistent | Yes |
| **Same affine solution space over Q** | **Yes** |
| Common rational dimension | 2 |

The *saturated integer kernel* of the scalar system was computed with an exact unimodular column-Hermite transform (`det(V)=±1`, `T*V=[C|0]`), then reduced/centered by exact rank-two Gauss operations. **It yields the same centered particular integer solution and the same two integer kernel basis vectors** as the Galois Stage-B certificate [`stage-b-trial0.md`](https://github.com/mehtaarnav/quantum_shenanigans/blob/207c2fb05c43df07b8ec824431702382eb7adfbc/crt_codes/results/stage-b-trial0.md). All ten message coefficient bounds were translated into 20 exact inequalities, followed by exact-rational vertex enumeration.

Certified coordinate bounds: `t1 ∈ [-1,0]`, `t2 ∈ [-1,0]`. This outer box has **four** integer pairs. Exactly **three** satisfy the coefficient box: `(-1,0)`, `(0,-1)` and `(0,0)`. The excluded point `(-1,-1)` violates message coefficient bounds.

All **three** exact scalar candidates have full-block agreement **6/20**; the independent exhaustive CRT reference returns the **identical three messages**. No candidate omission, timeout or heuristic truncation occurred. The reusable extraction code reports `complete=false` rather than issuing a false claim for oversized/unsupported polytopes; this run is fully tractable.

## 4. Stage-level wall-clock comparison (same GitHub runner, one fixed word)

| Step | Scalar affine | Galois |
|---|---:|---:|
| Invert received blocks, CRT streams, construct interpolation lattice | 0.002058 s | 0.011008 s |
| **Main interpolation LLL** | **0.000655 s** (11×11) | **13.752036 s** (90×90 padded) |
| Exact relation verification and height selection | 0.000250 s | 0.000108 s (short-vector collection only) |
| Equation system construction/rank | included in extraction | 0.575012 s |
| Exact saturated integer extraction, rational LP, enumeration | 0.022300 s | 0.052435 s original extraction stage; separate exact Stage B certification was also run |
| Exact block-agreement filtering | 0.000344 s | separately certified |
| **Core decoding** | **about 0.0256 s** (sum of scalar stages including exact extraction) | **14.3906 s** (original decoder, excluding separate full Stage-B certification) |

Same-run Galois full certificate (including its separate setup/reference and exact checks): 14.4814 s. Complete independent reference validation time: 0.0170 s. The difference is substantial **on this one input**, but the operations are not asymptotically equivalent, and timings must not be extrapolated to a general speed advantage.

## 5. Claim ledger

| Status | Statement |
|---|---|
| **Proved algebraically** | Exact scalar congruence lattice, determinant `N*H^10`, and U-based interpolation vanishing condition for all messages with at least six fully agreeing blocks. |
| **Certified for this input** | Eight exact height-valid scalar relations, rank eight, saturated integer kernel, exact affine plane, complete four-pair integer search and complete qualifying list of exactly three messages. |
| **Reproduced** | Original Galois decoder's same-word affine plane and qualified list are identical to the scalar comparator's; measured interpolation LLL times 0.000655s vs 13.752s. |
| **Falsified for this input** | Galois automorphism relations are *necessary* to produce a certified complete list on PR #2 Trial-0. Scalar affine relations suffice. |
| **Unresolved** | Whether scalar affine interpolation finds useful certificates on other received words; comparative worst-case decoding guarantees, list-size bounds and asymptotic runtime. No new inputs were tested. |

**Conclusion:** the original Trial-0 successful Galois decoding does not isolate any benefit due to Galois automorphisms. The 11D scalar affine construction already supplies exactly the same certified affine plane and complete three-message list, with a much smaller interpolation lattice on the fixed input. This does *not* imply general equivalence of the two relation families or decoders.
