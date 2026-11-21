# Exact Stage B certificate — PR #2 Trial 0 only

**Scope and provenance.** This report certifies the *existing* Trial-0 received word, not a new test regime: `M=11`, `n=20`, `bits=24`, `coeff_bits=71`, `s=7`, code seed `2`, received seed `20261121`, and full-block agreement threshold `6`. The original decoded six auxiliary vectors are preserved in the artifact. **One** CI reproduction: [run 44879104769](https://github.com/mehtaarnav/quantum_shenanigans/actions/runs/44879104769), including the downloadable `stage-b-certificate.json` artifact containing the entire integer matrix `T`, right-hand side `rhs`, exact unimodular transformation `V`, the selected auxiliary vectors, rational LP vertices, candidates, timing details, and code-instance metadata.

Source definitions:
- Cyclotomic field, primes, message bounds and encoding: [`crt_codes/folded.py:39–76`](https://github.com/mehtaarnav/quantum_shenanigans/blob/5ae11eb6bd0bf5ba35b887d672ced724c9ce5c18/crt_codes/folded.py#L39-L76).
- Extracted equations: [`crt_codes/folded.py:164–190`](https://github.com/mehtaarnav/quantum_shenanigans/blob/5ae11eb6bd0bf5ba35b887d672ced724c9ce5c18/crt_codes/folded.py#L164-L190), independently rebuilt in [`stage_b_certificate.py:41–59`](https://github.com/mehtaarnav/quantum_shenanigans/blob/714e159bba2af10254b27b02aff65df4f6943bda/crt_codes/stage_b_certificate.py#L41-L59).
- Saturated integer kernel by explicit column-unimodular Euclidean/Hermite elimination: [`integer_affine.py:19–112`](https://github.com/mehtaarnav/quantum_shenanigans/blob/714e159bba2af10254b27b02aff65df4f6943bda/crt_codes/integer_affine.py#L19-L112).
- Exact rank-two Gauss reduction, centering, rational LP, integer enumeration: [`integer_affine.py:115–260`](https://github.com/mehtaarnav/quantum_shenanigans/blob/714e159bba2af10254b27b02aff65df4f6943bda/crt_codes/integer_affine.py#L115-L260).
- Exact whole-block filtering and complete reference comparison: [`stage_b_certificate.py:181–221`](https://github.com/mehtaarnav/quantum_shenanigans/blob/714e159bba2af10254b27b02aff65df4f6943bda/crt_codes/stage_b_certificate.py#L181-L221).

## Stage A: exact uniform certificate

The six selected auxiliaries were checked at **all 200 interpolation congruences each**. With `H=2**71`, let `B_l=max_t |a_(l,t)|` and
```
U_A = B_0 + 4*10*(H - 1)*sum(B_l for l=1,...,7).
```
For each input with `|c_t|<=H-1`, the coefficient bound `||A_0+Σ_l A_l sigma^(l-1)(alpha_c)||_infinity <= U_A` follows from `sigma`'s at-most-twofold coefficient expansion and at most `2*phi` signed terms in each product coefficient after cyclotomic reduction. The maximum bound over all six selected auxiliaries, and the product of the six smallest primes, are **exactly**

```
U_max   = 166123626202756937457210496757940727
2**117  = 166153499473114484112975882535043072
2**138  = 348449143727040986586495598010130648530944
P_min,6 = 412140041510669591866160348537616348856951
```

Thus `U_max < 2**117 < 2**138 <= P_min,6`. Every message with at least six entire agreeing blocks satisfies **all six selected sigma-linear equations**, because every coefficient of its residual is divisible by the agreement-prime product and strictly smaller in absolute value. This certificate applies to *this received word* and these selected interpolation vectors only. [Source: `stage_b_certificate.py:62–104`](https://github.com/mehtaarnav/quantum_shenanigans/blob/714e159bba2af10254b27b02aff65df4f6943bda/crt_codes/stage_b_certificate.py#L62-L104).

## Stage B: exact integer kernel saturation

The extracted integer system `T*x=rhs` has shape **60×10**, exact rank **8**, integer nullity **2**. Column-Hermite-style elimination uses elementary unimodular column operations, constructing a **10×10** integer matrix `V` and checking:
- `det(V) = +1`;
- `T*V = [C | 0_(60×2)]`, with `rank(C)=8`;
- `T*x0=rhs`, `T*b1=T*b2=0`.

Since `V` is unimodular, its final two columns form a **saturated** basis of `ker_Z(T)`, and rank-eight uniqueness of the first eight transformed coordinates proves the affine solution set exactly equals `x=x0+t1*b1+t2*b2` for **all** `t1,t2 ∈ Z`. An exact `GL_2(Z)` Gauss change and integral centering preserve saturation. The full `V` and `T,rhs` are in the JSON artifact; the final centered vectors are:

```text
x0 = [
 -828649802069662333576, -1271129642427256802731,
 -623314204177818072439, -702641666241335307047,
 1710947054425240092502, -320152011377461650559,
 -1434835713324810360781, -1538151554043433249678,
 558252923409689554722, 568859911045947385363
]
b1 = [
 -236285373907600196198, -3308282403115859050702,
 -2461733346264268678763, 414645705511447484873,
 3038748813913881422090, -2620288381565678428504,
 -1284991020155376603465, -1400299186487361550861,
 733821663516656472588, -203292704940382506093
]
b2 = [
 711159013783786911102, -2642617947098160344645,
 149933098514155870039, -2860985177897567283102,
 2372072875615061825485, 1263119231138483830899,
 -925042588950427438261, -2024113908455670529250,
 1773609431996124201204, 2562620157101410629403
]
```

## Exact rational polytope and coordinate bounds

For all ten coefficients, enforce the twenty inequalities
`-H+1 <= x0[t] + t1*b1[t] + t2*b2[t] <= H-1`.

Every pairwise intersection of constraint boundaries was solved over exact `Fraction` arithmetic; feasible vertices were checked against **all twenty** inequalities. Extremizing each coordinate over this bounded rational polygon gives:

| LP variable | Exact continuous minimum | Exact continuous maximum | Certified integer interval |
|---|---|---|---|
| `t1` | `-5836894459193502025139808959595595183320746 / 5551588114555466846952768623088864533501089` | `1132560130339254853208509870628040558229167 / 2010763412331023353095196322330721100156670` | **[-1, 0]** |
| `t2` | `-11642120966660417022778719283162872422567182 / 10560697102685547108810544268560323109892689` | `1061901777391546391353350753636817573117037 / 1935477029399080589096976123688882537982317` | **[-1, 0]** |

Taking `ceil(min)` and `floor(max)` therefore leaves precisely **four integer coordinate pairs**. The complete inequality sweep covers all four (without any radius cutoff):
- `(t1,t2)=(-1,-1)`: **rejected by the original coefficient box**, with violations at coefficient indices 1 and 4 (zero-based);
- `(-1,0)`, `(0,-1)`, `(0,0)`: **all three admissible**.

Thus **four outer-box combinations were covered**, **three actual admissible integer points were constructed and exactly checked**, and no admissible common solution of the selected equations was omitted. The explicit search budget in the reusable implementation is not binding here: the certified outer box size is **4**, far below the declared budget.

## Exact agreement filter and independent reference

All three admissible integer points have whole-block agreement **6/20**. Independent `exhaustive_decode(code,received,6)` enumerates all 4-subsets of matching blocks and reconstructs **exactly the same three codewords**:
- Stage B admissible common solutions: **3**.
- Stage B qualifying messages (agreement ≥6): **3**.
- Reference qualifying messages: **3**.
- Missing or extra messages: **0**.

The planted messages are used *only* to validate the result, not to construct `T`, solve the integer system, derive the coordinate bounds, or enumerate candidates. The original heuristic radius-two decoder also returned the same three in this run, but the **new** candidate enumeration has a separate exact completeness certificate.

## Stage-level wall-clock timings

One GitHub Ubuntu hosted runner, one reproduction; all values seconds. Runtime variability does not imply comparative speed guarantees.

| Step | Seconds |
|---|---:|
| Code setup | 0.000562 |
| Existing Trial-0 received word generation | 0.000799 |
| Original lattice decoder total | **21.982017** |
| — Constraint CRT / interpolation lattice construction | 0.021386 |
| — Main interpolation LLL | **20.814396** |
| — Short vector collection | 0.000157 |
| — Sigma-linear equation/rank operations | 1.033623 |
| — Integer extraction in original decoder | 0.098914 |
| Rebuild all selected equations exactly | 0.007195 |
| Stage A exact uniform-bound certificate | 0.010383 |
| Unimodular column-Hermite elimination | 0.005075 |
| Exact rank-two Gauss reduction and centering | 0.001025 |
| Check `det(V)`, `T*V`, and rank | 0.101877 |
| Exact rational LP | 0.006664 |
| Complete polytope integer enumeration | **0.000097** |
| Exact whole-block filtering | 0.000724 |
| Independent reference decoder | 0.032740 |

The main interpolation LLL dominates; do not claim speedups from these timings. Full raw timing data is in the artifact.

## Updated claim ledger

| Class | Statement |
|---|---|
| **Proved** | Unimodular transformation, rank 8, saturated integer kernel, exact affine solution lattice, exact rational LP bounds, and exhaustive enumeration of all admissible solutions for the **fixed selected equation system**. |
| **Proved for this input** | All qualifying messages necessarily satisfy all six selected sigma equations by the numeric Stage A bound; Stage B exhaustively enumerates their admissible integer solutions; after exact block verification the complete list for this particular received word is exactly three. |
| **Reproduced** | Same original Trial-0 seed/code; three outputs from old heuristic decoder, three from new exact enumeration, three from independent complete reference; all agree exactly. |
| **Conditional** | General completeness for other received words if Stage A height and interpolation certificates hold and a tractable, fully enumerated saturated integer solution polytope can be certified. |
| **Heuristic** | The original general decoding-radius prediction and the original `radius=2` shortcut outside this input. |
| **Unresolved** | General polynomial-time Stage B extraction/list size, worst-case Stage A height guarantees, scaling and performance relative to known methods. No baseline, sweep, or new error regime was run. |

**Conclusion:** Stage B is **fully certified for the existing PR #2 Trial-0 instance**. This is not a worst-case theorem and does not establish that Stage B is tractable for other inputs.
