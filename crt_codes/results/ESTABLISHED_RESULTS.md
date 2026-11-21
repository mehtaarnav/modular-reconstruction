# Established results and certificate ledger

**Status (2026-11-21):** This is the canonical account of what has been established in the CRT/cyclotomic research. It supersedes any *current-status or novelty* implications in the exploratory [RESEARCH.md](../RESEARCH.md), [crt_codes/README.md](../README.md), and saved output snapshots. Those files remain useful historical records, **not proofs of a new decoding radius**.

**Scope:** Only the already executed PR #2 Trial-0 instance has full Stage A + Stage B + independent-reference certificates. The original three-trial workflow has successful outputs, but trials 1 and 2 do **not** have the same full certificate. No new received word, parameter sweep, obstruction counterexample, or random-noise rank theorem is presented here.

## 1. Instance and definitions

The fixed received word is PR #2 **trial 0**: field \(K=\mathbb Q(\zeta_{11})\), \(\phi(11)=10\); \(n=20\) distinct primes \(p_i\equiv1\pmod {11}\); \(\sigma(\zeta)=\zeta^2\); code seed **2**; received seed **20261121**; prime-size parameter **bits=24**; coefficient-bound parameter **coeff_bits=71**; Galois interpolation order **s=7**. The bound is \(H=2^{71}=2361183241434822606848\); the messages are all ten-tuples \(c\in\mathbb Z^{10}\) satisfying \(|c_j|\le H-1\), with cardinality \((2H-1)^{10}\).

Each block is the complete vector of ten evaluations
\[
C_i(c)=\left(\sum_{j=0}^{9}c_j x_{i\ell}^{\,j}\bmod p_i\right)_{\ell=0}^{9},
\quad x_{i\ell}=\omega_i^{2^\ell}\in\mathbb F_{p_i}.
\]
A **matching block** means exact equality of all ten symbols, and \(a(c,r)=\#\{i:C_i(c)=r_i\}\). The certified list threshold is **\(a(c,r)\ge6\) matching blocks**.

- Source definition: [\`folded.py\` lines 39–76, pinned at PR #2](https://github.com/mehtaarnav/quantum_shenanigans/blob/5ae11eb6bd0bf5ba35b887d672ced724c9ce5c18/crt_codes/folded.py#L39-L76).
- Exact experiment/seed: [\`triad_audit.py\` lines 27–66, pinned](https://github.com/mehtaarnav/quantum_shenanigans/blob/5ae11eb6bd0bf5ba35b887d672ced724c9ce5c18/crt_codes/triad_audit.py#L27-L66).
- Full 20-prime list and received-word fingerprint: [scalar comparator JSON artifact](https://github.com/mehtaarnav/quantum_shenanigans/actions/runs/44881119774/artifacts/14697369959); SHA-256 of canonical received JSON: \`cbcd857846e26de3c42e70498d04957c08ddbe38f7206adc92b6ff79eebf5fcb\`.

## 2. General mathematical facts (not decoding-radius theorems)

**Blockwise equivalence to interleaved CRT.** Each full cyclotomic block is an invertible Vandermonde transform of the ten coefficient residues \(c_j\bmod p_i\). Thus in the full-block Hamming metric this construction is blockwise equivalent to ten interleaved CRT streams. It is **not** a newly proved code family distinct from interleaved CRT.

**Four-block reconstruction and distance.** Let \(q\) be the smallest index such that \(\prod_{i=1}^{q}p_i>2(H-1)\). For the fixed primes, \(q=4\): the product of the three smallest primes is below \(2(H-1)\), and that of the four smallest is above it. Every four complete matching blocks uniquely determine every signed coefficient. The exact minimum whole-block distance of this particular code is **17**, since distinct messages can coincide on at most three blocks, and a valid difference in a constant coefficient equal to \(p_1p_2p_3\) realizes three coincidences.

**Exact reference-list completeness.** For every received word with these parameters, any qualifying codeword with at least six matching blocks is reconstructed from at least one four-subset of those matching blocks. The reference decoder enumerates all \(\binom{20}{4}=4845\) four-subsets, reconstructs signed coefficients by modular Vandermonde inversion + CRT, checks the original box, and **re-encodes and tests all 20 blocks**. This proves completeness **for the specified agreement threshold**, not efficiency as \(n\) grows. See [\`exhaustive.py\` lines 16–103](https://github.com/mehtaarnav/quantum_shenanigans/blob/5ae11eb6bd0bf5ba35b887d672ced724c9ce5c18/crt_codes/exhaustive.py#L16-L103).

**Scalar height-certificate soundness.** For arbitrary \(n\) distinct primes, \(m\) interleaved coefficient streams with \(|c_j|\le H-1\), and target \(t\) full matching blocks, define \(N=\prod_i p_i\), \(P_t=\prod_{i=1}^{t}p_i\) (sorted smallest), and CRT residues \(R_j\bmod N\). Every exact modular relation
\[
a_0+\sum_{j=1}^m a_jR_j\equiv0\pmod N
\quad\text{with}\quad
U(a)=|a_0|+(H-1)\sum_{j=1}^m|a_j|<P_t
\]
is an **integer equality** for *every* qualifying message: its value is divisible by the product of any \(t\) agreeing primes yet has absolute value \(<P_t\). This is a general deterministic *soundness* theorem for retained relations, not a guarantee that enough independent relations exist.

**Relation lattice and prior-art identification.** The weighted hyperplane lattice has dimension \(m+1\), basis \((N,0,\dots,0)\), \((-R_j,He_j)\), and determinant \(NH^m\). This is the degree-one/multiplicity-one specialization of the Cohn–Heninger approximate-common-divisor interpolation basis. The interleaved CRT locator lattice \(\mathcal E=\{(\lambda,u):u_j-\lambda R_j\equiv0\bmod N\}\) satisfies \(\mathcal E=N\mathcal H^*\); with complementary \(H\)-weights, \(\mathcal E D_E=NH(\mathcal H D_H)^*\), and its determinant is \(HN^m\). Thus there is **published interpolation prior art** and an exact duality, not a novel primitive. [Literature and proof note](../notes/scalar-affine-duality-literature.md).

**One-relation existence is not a rank theorem.** Minkowski yields a nonzero certified scalar relation if
\[
P_t^{m+1}>(m+1)!\,N(H-1)^m.
\]
For standard \(3/4\)-LLL on the implemented \(NH^m\)-determinant lattice, a sufficient condition to obtain **at least one** certified row is
\[
P_t>
\sqrt{1+m(1-1/H)^2}\;2^{m/4}(NH^m)^{1/(m+1)}.
\]
Both statements concern one relation; **neither forces independent rank**, a small output list, or polynomial-time extraction. Details and limitations: [theory note](../notes/scalar-affine-duality-literature.md).

**Conditional extraction correctness.** If the certified equation system \(Tc=b\) is consistent, exact unimodular column-Hermite/Smith operations recover a particular integer solution \(c_0\) and a **saturated** kernel \(B\) such that all integer solutions are \(c=c_0+Bz\), \(z\in\mathbb Z^k\), \(k=m-\operatorname{rank}_{\mathbb Q}T\). Translating the coefficient bounds into exact polytope inequalities, bounding them with rational LP, and **enumerating every feasible integer point** before exact agreement filtering returns the complete qualifying list. If the polytope is too large or the implementation cannot handle its dimension, it must report noncompletion; there is **no general polynomial extraction guarantee**.

## 3. Certified for *exactly* PR #2 Trial 0

### Independent reference

The complete four-block-subset reference found **exactly three** codewords with \(\ge6\) matching blocks. Each has exactly six matching blocks, and there were zero additional qualifying messages. The original PR #2 three-trial CI succeeded, but only **trial 0** was subjected to the detailed exact stage-level proof and comparator.

[PR #2 reference run](https://github.com/mehtaarnav/quantum_shenanigans/actions/runs/44875984991) · [PR #2 fixed-trial protocol](../results/triad-protocol.md)

### Galois Stage A — exact numerical certificate

The six selected interpolation vectors passed all **\(20\cdot10=200\)** congruences each. For \(A_0,\ldots,A_7\), define \(B_\ell=\max_j|a_{\ell j}|\) and
\[
U_A=B_0+4\cdot10\,(H-1)\sum_{\ell=1}^{7}B_\ell.
\]
The recorded maximum over the six vectors was
\[
U_{\max}=166123626202756937457210496757940727.
\]
With the product of the six smallest primes
\[
P_{\min,6}=412140041510669591866160348537616348856951,
\]
the **integer certificate** is
\[
\boxed{U_{\max}<2^{117}<2^{138}\le P_{\min,6}.}
\]
Consequently **every** message with \(\ge6\) agreeing blocks satisfies **all six** selected Galois-linear equations, not merely the three planted messages. This is a fixed-received-word certificate, not a universal LLL height guarantee.

[Certificate implementation, lines 62–104](https://github.com/mehtaarnav/quantum_shenanigans/blob/714e159bba2af10254b27b02aff65df4f6943bda/crt_codes/stage_b_certificate.py#L62-L104) · [Single-run JSON artifact](https://github.com/mehtaarnav/quantum_shenanigans/actions/runs/44879104769/artifacts/14696999034)

### Galois Stage B — certified complete integer extraction

The Galois equation system \(Tc=b\) had shape **60 × 10**, rational rank **8**, and integer nullity **2**. The exact unimodular transform \(V\) satisfied \(\det V=+1\), \(TV=[C\mid0]\), and \(Tc_0=b\). Its two kernel columns, after exact unimodular reduction, span **all** integer solutions. Exact rational LP on the original coefficient box gave \(z_1,z_2\in[-1,0]\). Of **four** outer-box pairs, one violated coefficient limits; **three** were admissible, passed the six-block agreement test, and matched the complete reference list.

[Detailed Galois exact Stage-B report](../results/stage-b-trial0.md) · [Integer-kernel and LP code](../integer_affine.py) · [JSON artifact](https://github.com/mehtaarnav/quantum_shenanigans/actions/runs/44879104769/artifacts/14696999034)

### Scalar affine comparator — same solution plane, same complete list

The **11-dimensional scalar affine** lattice uses weights \((1,H,\ldots,H)\). All eleven LLL-reduced basis rows satisfied the exact modular congruences; **eight** passed \(U(a)<P_{\min,6}\). Their exact rational coefficient rank was **8**, leaving nullity **2**. The scalar and Galois systems' stacked coefficient and augmented matrices both had rank **8**, proving their **entire rational affine solution spaces coincide**, not merely the final three-message lists.

An independent exact saturated-kernel/LP extraction for the scalar system again gave \(z_1,z_2\in[-1,0]\), **four** bounded pairs, **three** admissible messages, and **exact equality** with the reference. It follows that **Galois structure was not necessary to recover or certify the list on this instance**.

The two procedures had markedly different **single-run** timing on one GitHub runner: scalar core **about 0.0256 s**, Galois original decoder **about 14.3906 s**; their interpolation LLL stages were **0.000655 s** (11 × 11) versus **13.752036 s** (90 × 90). These are profiling observations for one fixed instance, **not** a speedup theorem.

[Full scalar relation table and exact affine equality](../results/scalar-affine-trial0.md) · [Scalar implementation](../scalar_affine_trial0.py) · [Exact JSON artifact](https://github.com/mehtaarnav/quantum_shenanigans/actions/runs/44881119774/artifacts/14697369959)

## 4. Claim ledger: keep the levels separate

| Class | Established statement / limitation |
|---|---|
| **Proved generally** | Modular Vandermonde/CRT inversion, conditional \(q\)-block uniqueness, scalar \(U<P_t\) soundness, scalar weighted lattice determinant and duality, one-relation Minkowski/LLL sufficient bounds, and correctness of **fully completed** saturated-polytope extraction. |
| **Certified only for Trial 0** | Galois and scalar height bounds, both systems' rank **8** and identical affine rational plane, unimodular saturated integer kernels, exact LP box of size **4**, exactly **3** admissible/qualifying messages, complete independent agreement-reference equality. |
| **Reproduced empirically** | Original three-trial PR #2 recovery; one-run timing comparison for the two implementations. Trials 1–2 do **not** have full Stage A/B certification. |
| **Prior art / interpretation** | Cyclotomic full blocks are interleaved CRT under an invertible transform; scalar linear interpolation specializes Cohn–Heninger and is dual to the interleaved CRT locator lattice. Claims of a wholly new code or scalar interpolation primitive are **retired**. |
| **Heuristic, not a theorem** | Galois decoding-radius expression \(\frac{s}{s+1}(1-R)\), any random-noise covolume/rank estimate, and extrapolation of single-run runtimes. **No overwhelming-probability rank statement is established.** |
| **Open / not established** | A general guarantee of enough independent certified relations, polynomial list-size or complete-extraction complexity, or any novel worst-case beyond-Johnson decoding guarantee. |

### Obstruction construction — narrow statement only, *not* a certified decoding counterexample

The separate proposed obstruction argument's **limited relation-level conclusion** is that **all height-certified scalar relations leave the final message-coordinate axis unconstrained** (their coefficients in that coordinate vanish). This does **not** establish that the qualifying-message set itself is confined to the zero message, nor does it establish a rank-gap counterexample for an instance with a unique qualifying message.

**Still to be established:** selecting received residues for which **zero is the unique qualifying message** while retaining that unconstrained relation axis. No such completed residue construction or certified counterexample is claimed in this ledger. This proposal is **excluded from the proved-results table** above.

Likewise, a random-noise **covolume calculation is heuristic only**. It provides no rigorous overwhelming-probability assertion about the rank of the certified relation system; no random-noise rank claim appears among the proven statements.

### Historical numbers needing care

- The saved **7.91** Johnson agreement threshold is \(20\sqrt{R_{\mathrm{proxy}}}\), with \(R_{\mathrm{proxy}}=0.156246475000973\), **not** a distance-derived Johnson theorem. The exact whole-block minimum distance here is **17**, giving the corresponding large-alphabet Johnson agreement proxy \(\sqrt{20(20-17)}=\sqrt{60}\approx7.746\).
- The saved [\`beyond_johnson.txt\`](../results/beyond_johnson.txt) shows **14.38** predicted errors, whereas the formula in the research note/\`beyond_johnson.py\` computes about **14.77**. They must not be silently conflated.
- The historical scalar CRT benchmark [\`experiment.py:49–59\`](https://github.com/mehtaarnav/quantum_shenanigans/blob/5ae11eb6bd0bf5ba35b887d672ced724c9ce5c18/crt_codes/experiment.py#L49-L59) derives its \`log_b\` decoder input from the **hidden transmitted message**. Its recorded results are **oracle-assisted**, not a fair stand-alone decoder performance guarantee.

## 5. Certificate index and reproduction (only the already fixed instance)

| Artifact | Purpose |
|---|---|
| [PR #2, original three-trial workflow](https://github.com/mehtaarnav/quantum_shenanigans/actions/runs/44875984991) | Empirical recovery outcomes; **not** full theorem certificates for all trials |
| [PR #3, instrumented Trial 0](https://github.com/mehtaarnav/quantum_shenanigans/actions/runs/44877406000) | Stage-level Galois timings and instance-specific height verification |
| [PR #4, exact Galois Stage B](https://github.com/mehtaarnav/quantum_shenanigans/actions/runs/44879104769) | Full input, 60×10 integer system, \`det V=1\`, exact kernel, rational LP, list and reference equality |
| [PR #5, scalar/Galois same-word comparison](https://github.com/mehtaarnav/quantum_shenanigans/actions/runs/44881119774) | Eleven scalar basis relations, eight certified height rows, exact affine equality and timings |
| [Literature/theorem note](../notes/scalar-affine-duality-literature.md) | Published interpolation identity/duality and four-stage conditional correctness theorem |

Run an **existing** fixed instance only, from repository root:
\`\`\`sh
python -m pip install "sympy==1.14.0" "fpylll==0.6.4" "cysignals==1.12.5"
python -m crt_codes.stage_b_certificate --output stage-b-certificate.json
python -m crt_codes.scalar_affine_trial0 --output scalar-affine-trial0.json
\`\`\`
The scalar comparator invokes the fixed Galois certificate internally for same-word comparison. These two commands are independent reproducibility entrypoints, not new input regimes or parameter sweeps.

**No further investigation is initiated by this consolidation.**
