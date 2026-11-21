# Scalar affine-list decoder: bounded literature and theorem audit

**Scope.** Theory and literature only. No new received words, code parameters, seeds, experiments, or implementations. This note considers the fixed prior Trial-0 comparator as motivation, but proves statements for arbitrary n, s, H, t and pairwise distinct primes. Original exact comparator: [scalar_affine_trial0.py](https://github.com/mehtaarnav/quantum_shenanigans/blob/3dc49606cb45406a90e037b31357a49cd5e8f00a/crt_codes/scalar_affine_trial0.py); original complete certificate: [Trial-0 report](https://github.com/mehtaarnav/quantum_shenanigans/blob/3dc49606cb45406a90e037b31357a49cd5e8f00a/crt_codes/results/scalar-affine-trial0.md).

## Result / novelty assessment

**Outcome (a), exact equivalence of the interpolation *lattice* to an existing published specialization:** Our 11-dimensional construction is the degree-one, multiplicity-one case of Cohn–Heninger's *Approximate common divisors via lattices* (ANTS X, 2012), with equal coefficient scales X_j=H. It is also the **scaled dual** of Li–Sidorenko–Nielsen's (ISIT 2013) weighted interleaved CRT locator lattice. These are precise identities of Z-lattices, not merely analogies.

A distinct implementation decision remains: retain **every independently height-certified** linear relation produced by reduction; check the exact coefficient matrix rank rather than assuming independent roots; compute a **saturated** integer affine kernel using unimodular transforms; impose exact coefficient-box inequalities; and enumerate all feasible integer points, returning an explicit noncompletion result if infeasible at the declared resource budget. The inspected papers do not specify this precise complete finite-box list pipeline, but Cohn–Heninger already use multiple short equations, an l1 height test, dependency-aware root solving, and additional relations. Consequently this pipeline must **not** be claimed novel without a broader prior-art investigation. There is no new decoding-radius or polynomial-time list-size theorem here.

Bibliographic corrections:
* Bleichenbacher–Kiayias–Yung (2003), *Decoding of Interleaved Reed Solomon Codes over Noisy Data*, appeared at **ICALP 2003**, not PKC 2003: https://www.research.ed.ac.uk/en/publications/decoding-of-interleaved-reed-solomon-codes-over-noisy-data/
* The four-author Abbondati–Afflatet–Guerrini–Lebreton work is **ITW 2023**: https://esploro.umontpellier.fr/esploro/outputs/conferenceProceeding/Probabilistic-Analysis-of-LLL-based-Decoder-of/99102343009311
* The **Journal of Symbolic Computation 2026**, vol. 132, article 102481 work on simultaneous rational number codes has three authors Abbondati–Guerrini–Lebreton and is not the ITW 2023 paper: https://www.sciencedirect.com/science/article/pii/S074771712500063X

## Common mathematical setup

Let s>=1, H>=2, n>=1, and 1<=t<=n. Let p_1<...<p_n be distinct rational primes, N=prod_i p_i, and R_j in Z/NZ be the CRT lift of the received residue in stream j, for j=1,...,s. Messages are integer vectors c=(c_1,...,c_s) with |c_j|<=H-1. Their block at prime p_i is (c_j mod p_i)_{j=1}^s. Write a(c,r)=number of whole blocks in agreement, and P_t=prod_{i=1}^t p_i (smallest t-prime product).

If a(c,r)>=t and S is any t agreeing blocks, then for **every** stream j,
R_j-c_j == 0 mod P_S, with P_S=prod_{i in S} p_i >= P_t. In particular,
gcd(N,R_1-c_1,...,R_s-c_s)>=P_S>=P_t.
The last inequality connects whole-block agreement to the large common-divisor objective in Cohn–Heninger; the converse generally fails for unequal primes because the gcd may have few large prime factors.

Define
Hcal={a=(a_0,...,a_s) in Z^(s+1): a_0+sum_j a_j R_j == 0 mod N}.
A basis is N*e_0 and (-R_j,e_j), j=1,...,s; det(Hcal)=N.
The actual interpolation LLL lattice is Hcal * D_H, D_H=diag(1,H,...,H), with det=N*H^s.
For a relation a use **exact message-box** certificate
U(a)=|a_0|+(H-1)*sum_j |a_j|.
An integer relation is certified if its congruence holds and U(a)<P_t.

## Exact comparison of relation constructions

| Work | Actual relation mechanism | Dimension/index or determinant | Output / completeness assumption |
|---|---|---|---|
| **Our scalar affine-list** | a_0+sum_j a_j R_j==0 mod N; each certified a gives a_0+sum_j a_j c_j=0 for all qualifying c | d=s+1, det weighted N H^s | Retain several certified affine equations, exact rank, saturated integer kernel, exact bounded polytope and agreement verification. Worst-case list enumeration potentially exponential. |
| **Bleichenbacher–Kiayias–Yung, ICALP 2003** | Over a finite field use a common error-locator polynomial E(X), with E(x_i)y_(ij)=Q_j(x_i), degrees bounded by error count and message degrees | Simultaneous polynomial interpolation/key equations over F_q; **not** an integer lattice of determinant N H^s or H N^s | Probabilistic simultaneous **polynomial** recovery for interleaved RS; not a deterministic CRT integer-list certificate. |
| **Cohn–Heninger, ANTS X 2012** | Construct polynomials from prod_j (X_j x_j - R_j)^(e_j) N^ell, with degree and multiplicity budgets. The *degree=1,multiplicity=1* generators are precisely N and (H x_j-R_j) | **Exactly our Hcal*D_H** at degree/multiplicity one; det N H^s. General lattice dimension binomial(degree+s,s), not necessarily s+1 | Uses weighted l1 for vanishing; combines **multiple relations**, sometimes adds more if dependent; general root finding assumes suitable algebraic independence. Their general algorithm includes our relation lattice as a special case. |
| **Li–Sidorenko–Nielsen, ISIT 2013** | Common error-locator lambda and u_j=lambda c_j, with u_j==lambda R_j mod N; basis (1,R) and (0,N e_j) | d=s+1, unweighted det N^s, equal-height normalization weight diag(H,1,...,1) has det H N^s | Recovers locator or reports failure; short-vector analysis and channel-probability model, not a worst-case all-message-list certificate. |
| **Abbondati–Afflatet–Guerrini–Lebreton, ITW 2023** | Analyzes the LSN LLL-based ICR locator decoder | Same **locator-side** module; no new hyperplane lattice identity required | Provides rigorous *probabilistic failure bounds* for random-error model; no claim of arbitrary-adversarial complete lists. |
| **Abbondati–Guerrini–Lebreton, ISSAC 2024** | Simultaneous rational reconstruction: write f_j/g as received residues and seek phi=lambda*g, psi_j=lambda*f_j, with phi R_j==psi_j mod N | Locator lattice enlarged to accommodate shared unknown denominator, same key-equation direction; polynomial-field analogue included | Analysis of decoding probability/uniqueness, not exact finite-box enumeration of every admissible vector in a positive-dimensional rational affine space. |
| **Abbondati–Guerrini–Lebreton, JSC 2026** | Extends simultaneous rational number code recovery to prime powers, denominator bad primes and hybrid errors | Extensions of simultaneous numerator/denominator/key-equation framework | Probabilistic reconstruction in richer modular model, not the same pure integer coefficient-box list guarantee. |

Sources:
* BKY: https://www.research.ed.ac.uk/en/publications/decoding-of-interleaved-reed-solomon-codes-over-noisy-data/
* Cohn–Heninger, especially §2, equations for generating products and l1 coefficient norms: https://mathweb.ucsd.edu/~kedlaya/ants10/cohn/paper.pdf
* Li–Sidorenko–Nielsen, equation (5), Algorithm 1 and Corollary 1: https://backend.orbit.dtu.dk/ws/files/119987681/ISIT.ICRT.pdf
* Abbondati et al., ITW 2023: https://esploro.umontpellier.fr/esploro/outputs/conferenceProceeding/Probabilistic-Analysis-of-LLL-based-Decoder-of/99102343009311
* Abbondati et al., ISSAC 2024: https://doi.org/10.1145/3666000.3669686 ; accessible project slides https://www.cirm-math.fr/RepOrga/3047/Slides/Slides_Matteo_Abbondati.pdf
* Abbondati et al., JSC 2026: https://www.sciencedirect.com/science/article/pii/S074771712500063X
* Abbott, *Fault-Tolerant Modular Reconstruction of Rational Numbers* (JSC 2017): https://arxiv.org/abs/1303.2965 (scalar rational reconstruction with faulty residues, distinct from s-stream affine-list recovery).

### Two exact equivalences

**(I) Cohn–Heninger degree-one basis identity.**
Cohn–Heninger generate the polynomial lattice by products
prod_j (X_j x_j-R_j)^(e_j) N^ell,
with sum e_j <= degree_budget and ell=max(multiplicity_budget-sum e_j,0).
Taking both budgets equal to one and X_j=H leaves exactly
N, H*x_1-R_1, ..., H*x_s-R_s.
Their coefficient lattice, in monomial order (1,x_1,...,x_s), is therefore exactly Hcal*D_H.
Their generic weighted l1 polynomial-height test is |a_0|+H sum|a_j|. Ours uses the *sharper exact admissible bound* H-1 in U(a), without changing the LLL lattice. Their paper explicitly states that the l1 norm bounds evaluations and uses several short polynomials to reconstruct common roots. **Do not claim this interpolation method or "multiple short vanishing relations" is new.**

**(II) Interleaved-CRT weighted-dual identity.**
Define the locator lattice
Ecal={(lambda,u_1,...,u_s) in Z^(s+1): u_j-lambda R_j ==0 mod N}.
Its basis is (1,R_1,...,R_s) and (0,N e_j), so det(Ecal)=N^s.
For any v=(lambda,u), v lies in N*(Hcal)^* iff v·(N e_0) divisible by N (automatic) and v·(-R_j,e_j)==u_j-lambda R_j==0 mod N for all j.
Therefore **Ecal=N*(Hcal)^*** as exact integral lattices.

Let D_E=diag(H,1,...,1). Since D_H*D_E=H I,
**Ecal D_E = N H (Hcal D_H)^***.
The locator-side weighted determinant is H N^s; the hyperplane-side weighted determinant is N H^s. Duality does not imply LLL outputs or end-to-end decoding guarantees coincide.

## Proposition: strongest unconditional general guarantee supported by the scalar method

Let all assumptions in the common setup hold. Write d=s+1, and let A be **any finite set of exact modular relation vectors from Hcal** satisfying U(a)<P_t. Build integer equations
T c=b by taking each retained vector a=(a_0,a_1,...,a_s) as row (a_1,...,a_s), with right side -a_0.

**1. Soundness for arbitrary received words.**
Every admissible c with a(c,r)>=t satisfies T c=b. Proof: choose t agreeing primes S. For every retained a,
P_S divides a_0+sum_j a_j c_j because its reduction modulo every agreeing p_i equals the interpolation congruence. Its absolute value is at most U(a)<P_t<=P_S; hence it is zero. This proof needs neither randomness nor LLL nor a rank assumption. If the equations are inconsistent, no qualifying c exists.

**2. Existence by geometry of numbers.**
Rescale instead by D_{H-1}=diag(1,H-1,...,H-1). The resulting full-rank lattice has determinant
N(H-1)^s, and U(a)=||a D_{H-1}||_1.
The l1 cross-polytope of radius P_t has volume 2^d P_t^d/d!.
Minkowski's convex-body theorem therefore yields a **nonzero certified relation** whenever

  P_t^d > d! N(H-1)^s.

Strict inequality permits a slightly smaller body, giving U(a)<P_t. Because pure constant relations have |a_0|>=N>=P_t, the resulting relation has at least one nonzero message coefficient. The determinant alone guarantees **at least one** relation, not enough independent relations to decode or bound the list polynomially.

**3. LLL-computable certificate.**
Use the actual weighted lattice L_H=Hcal*D_H of determinant D=N H^s.
For a standard 3/4-LLL-reduced first basis row v,
||v||_2 <= 2^(s/4) D^(1/d).
Since
U(a) <= sqrt(1+s((H-1)/H)^2) * ||(a_0,H a_1,...,H a_s)||_2,
the **sufficient deterministic condition**

  P_t > sqrt(1+s((H-1)/H)^2) * 2^(s/4) * (N H^s)^(1/(s+1))

guarantees that LLL's first vector is height-certified. The 3/4 guarantee is conservative relative to stronger reductions; all actual output vectors still must be checked exactly for congruence and U<P_t. LLL has polynomial bit complexity in s and the input bit lengths for fixed reduction constants. Crucially, this does **not** certify independent second/third/... relations. Selecting many reduced basis rows is a valid *conditional* strategy; the observed rank must be computed.

For the fixed prior Trial-0 values, s=10, H=2^71, N<2^461 and P_6>=2^138. The displayed LLL upper bound is below 2^111, strictly less than P_6. Thus even a standard worst-case LLL bound guarantees **at least one** certified scalar relation on *every* received word with those modulus and box parameters; Trial 0's observed **eight independent** certified relations do **not** follow from that bound.

**4. Rank and complete integer extraction.**
Let r=rank_Q(T), k=s-r, and check rank(T)=rank([T|b]) before assuming consistency. When inconsistent, the qualifying list is empty. Otherwise rational solutions form an affine Q-space of dimension k. An SNF or column-HNF computation with explicit unimodular transformations either proves that no integer solution exists, or gives an exact integer particular solution c0 and a **saturated** integer kernel basis b_1,...,b_k with

  {c in Z^s: T c=b} = c0 + sum_{j=1}^k Z b_j.

Now impose all 2s exact inequalities

  -(H-1) <= c0_i + sum_j z_j b_(j,i) <= H-1.

They define a bounded rational polytope P in R^k. Solve exact rational LPs to find each coordinate's continuous minimum and maximum, and round outward using ceil(min), floor(max) to obtain a certified outer integer box. An exhaustive search of that box, with all inequalities and exact whole-block agreement checked, returns **exactly all qualifying messages**. This is an **output-correctness theorem if and only if the full enumeration completes**; any declared resource budget must report noncompletion, not silently truncate.

Standard exact HNF/SNF, rank, LLL and rational LP computations can be performed in polynomial bit complexity. Naive complete enumeration has up to

  B = prod_{j=1}^k max(0, floor(z_j,max)-ceil(z_j,min)+1)

integer coordinate tuples; time can be O(s B) arithmetic operations plus polynomial preprocessing and CRT verification. Although fixed-dimensional integer counting/feasibility has faster methods, **enumerating the full list** is inherently output-sensitive. It can be exponential in k*log H and no uniform polynomial-time list-size or complete-decoding guarantee follows.

Since k columns of a rank-k integer kernel basis have a nonsingular k×k coordinate minor, projection onto those k original message coordinates is injective on the affine solution set. Hence the number of *admissible integer common solutions* has the unconditional coarse bound (2H-1)^k. For k=0 there is at most one candidate. For k>0 that bound may still be exponential in the bit length of H.

**5. Independent complete CRT fallback and list bound.**
If q=min{j:prod_(i=1)^j p_i>2(H-1)} exists and t>=q, each qualifying message is uniquely identified by any q matching whole blocks. Exhaustive q-subset CRT reconstruction is therefore a complete, independent decoder with

  list size <= binomial(n,q)

and O(binomial(n,q)*poly(n,s,log N,log H)) bit complexity. This does not prove uniform polynomial time when q grows with input size, but is a useful deterministic completeness fallback independent of short relations.

### What the theorem does NOT establish

* Geometry of numbers plus LLL gives a **single** relation under an explicit sufficient threshold, not rank r=s, r=s-O(1), or a polynomial list size.
* The stacked equations need not be zero-dimensional. Algebraic independence is neither assumed nor an automatic consequence of multiple LLL vectors.
* Exact rank and saturation certify all integer points of the resulting affine lattice **if** full bounded enumeration is completed; they do not make exponential enumeration polynomial.
* Cohn–Heninger already present multiple l1-certified relations, and the scalar interpolation lattice is their degree-one specialization. The precise code-level pipeline might be an explicit certifying implementation detail, but no originality claim for it is justified without deeper prior-art work.

## Final classification

**(a) Explicit equivalence to published work** for the interpolation lattice and l1 short-polynomial relation mechanism (Cohn–Heninger, degree/multiplicity one), and exact weighted-lattice **duality** with the interleaved CRT locator lattice (Li–Sidorenko–Nielsen). There is only a **narrow unresolved comparison** regarding whether this specific positive-dimensional saturated-kernel/LP exhaustive finite-box output-certification pipeline has appeared in the modular-reconstruction literature.

The correct research statement is a deterministic **certificate/extraction proposition**, not a new interpolation primitive, a beyond-Johnson radius theorem, or an observed speedup. No new experiments were run.
