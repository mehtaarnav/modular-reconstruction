"""Galois-folded cyclotomic number-field codes and a linear-algebraic lattice decoder.

Field: K = Q(zeta_M), M prime, phi = M - 1, sigma: zeta -> zeta^g with g a
generator of (Z/M)^*, so Gal(K/Q) = <sigma> and the fixed field of sigma is Q.

Message: alpha in Z[zeta] with integer coefficients |c_t| <= H (power basis).
Block at a prime p = 1 mod M: fix a primitive M-th root omega mod p. The phi
primes above p are P_j = (p, zeta - omega^{g^j}), and the block symbol is

        y_j = alpha(omega^{g^j}) mod p,      j = 0..phi-1,

i.e. alpha reduced at the full Galois orbit of primes above p. Since
sigma^l(alpha) mod P_j = y_{j+l}, this is folded Reed-Solomon structure, with
each rational prime supplying one full orbit. (Polynomials can't do that,
because F_q^* has only one orbit under x -> gamma x.)

Decoder, the Guruswami-Wang linear-algebraic template with heights in place of
degrees: find short A_0..A_s in Z[zeta] such that, for every block i and
window j (cyclic),

        A_0(x) + sum_{l=1..s} A_l(x) * y_{i, j+l-1} = 0 mod p_i,  x = omega_i^{g^j}.

Then E = A_0 + sum A_l sigma^{l-1}(alpha) vanishes mod every prime above every
agreeing block. If it is small enough it is 0 in K, so alpha solves a sigma-linear
equation whose solution set is an affine Q-space of dimension <= s - 1.
Intersecting over several short vectors shrinks that space. Heuristic radius:
heuristic error fraction < (s/(s+1))*(1-R).
"""

from __future__ import annotations

import math
import random
import time
from fractions import Fraction

from sympy import Matrix, isprime, primitive_root


class CycloCode:
    def __init__(self, M: int, n: int, bits: int, coeff_bits: int, seed: int = 0):
        assert isprime(M)
        self.M, self.phi = M, M - 1
        self.g = int(primitive_root(M))
        rng = random.Random(seed)
        self.H = 1 << coeff_bits  # |c_t| < H
        ps, p = [], (1 << (bits - 1)) + rng.randrange(1 << (bits - 3))
        p -= p % M
        while len(ps) < n:
            p += M
            if isprime(p + 1):
                ps.append(p + 1)
        self.primes = ps
        self.omegas = []
        for q in ps:
            while True:
                w = pow(rng.randrange(2, q - 1), (q - 1) // M, q)
                if w != 1:
                    break
            self.omegas.append(w)
        # evaluation points per block: omega^{g^j}
        self.points = [[pow(w, pow(self.g, j, M), q) for j in range(self.phi)] for w, q in zip(self.omegas, ps)]

    @property
    def n(self):
        return len(self.primes)

    def rate(self) -> float:
        msg_bits = self.phi * math.log2(2 * self.H)
        code_bits = sum(self.phi * math.log2(q) for q in self.primes)
        return msg_bits / code_bits

    def random_message(self, rng) -> list[int]:
        return [rng.randrange(-self.H + 1, self.H) for _ in range(self.phi)]

    def encode(self, c: list[int]) -> list[list[int]]:
        return [[_ev(c, x, q) for x in pts] for pts, q in zip(self.points, self.primes)]

    def corrupt(self, word, e, rng):
        word = [list(b) for b in word]
        for i in rng.sample(range(self.n), e):
            q = self.primes[i]
            word[i] = [rng.randrange(q) for _ in range(self.phi)]
        return word

    # --- arithmetic in Z[zeta]/(Phi_M), power basis of length phi ---
    def reduce(self, poly: list) -> list:
        """Reduce a polynomial (any length) mod Phi_M = 1 + x + ... + x^{M-1}."""
        M, phi = self.M, self.phi
        r = [0] * M
        for i, c in enumerate(poly):
            r[i % M] += c  # x^M = 1
        top = r[M - 1]  # x^{M-1} = -(1 + ... + x^{M-2})
        return [r[i] - top for i in range(phi)]

    def mul(self, a, b):
        out = [0] * (2 * self.phi)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    out[i + j] += x * y
        return self.reduce(out)

    def sigma(self, a, times=1):
        e = pow(self.g, times, self.M)
        out = [0] * self.M
        for i, c in enumerate(a):
            out[(i * e) % self.M] += c
        return self.reduce(out)


def _ev(c, x, q):
    v = 0
    for coef in reversed(c):
        v = (v * x + coef) % q
    return v


def decode(code: CycloCode, received, s: int, n_vectors: int = 6, profile: dict | None = None):
    """Return (candidates, info). Candidates are integer coefficient vectors."""
    from fpylll import LLL, IntegerMatrix
    t_start = time.perf_counter()
    phi, n = code.phi, code.n
    D = (s + 1) * phi  # unknowns: coefficients of A_0..A_s
    # Rows of the constraint map Z^D -> prod_i F_{p_i}^{phi}, combined by CRT.
    N = math.prod(code.primes)
    rows = [[0] * D for _ in range(phi)]
    for j in range(phi):
        for d in range(D):
            l, t = divmod(d, phi)  # coefficient t of A_l
            # value of this unknown's contribution at block i, window j
            residues = []
            for i in range(n):
                q, x = code.primes[i], code.points[i][j]
                base = pow(x, t, q)
                if l == 0:
                    residues.append(base)
                else:
                    residues.append(base * received[i][(j + l - 1) % phi] % q)
            rows[j][d] = _crt(residues, code.primes)
    # Height weighting: A_0 may be ~H*phi*s larger than A_l (l >= 1).
    W0, W1 = 1, code.H * phi * s
    big = N * (1 << 20)
    dim = D + phi
    B = IntegerMatrix(dim, dim)
    for d in range(D):
        B[d, d] = W0 if d < phi else W1
        for j in range(phi):
            B[d, D + j] = rows[j][d] * big
    for j in range(phi):
        B[D + j, D + j] = N * big
    t_construct = time.perf_counter()
    LLL.reduction(B)
    t_reduce = time.perf_counter()
    vecs = []
    for r in range(dim):
        if all(B[r, D + j] == 0 for j in range(phi)):
            v = [B[r, d] // (W0 if d < phi else W1) for d in range(D)]
            if any(v[phi:]):
                vecs.append(v)
        if len(vecs) >= n_vectors:
            break
    t_vectors = time.perf_counter()
    # Each short vector gives a sigma-linear system for alpha. Vectors that are
    # not short enough give systems the true alpha does not satisfy, so add them
    # one at a time and keep only those consistent with what is already stacked.
    eqs, rhs = [], []
    selected_vectors = []
    skipped = 0
    used = 0
    for v in vecs:
        A = [v[l * phi:(l + 1) * phi] for l in range(s + 1)]
        # Matrix of alpha -> sum_l A_l * sigma^{l-1}(alpha) in the power basis.
        cols = []
        for t in range(phi):
            e_t = [0] * phi
            e_t[t] = 1
            acc = [0] * phi
            for l in range(1, s + 1):
                term = code.mul(A[l], code.sigma(e_t, l - 1))
                acc = [a + b for a, b in zip(acc, term)]
            cols.append(acc)
        new_eqs = eqs + [[cols[t][row] for t in range(phi)] for row in range(phi)]
        new_rhs = rhs + [-A[0][row] for row in range(phi)]
        Mx = Matrix(new_eqs)
        if Mx.rank() != Mx.row_join(Matrix(new_rhs)).rank():
            skipped += 1
            continue  # inconsistent: this vector was not short enough
        eqs, rhs = new_eqs, new_rhs
        selected_vectors.append(v)
        used += 1
        if Mx.rank() == phi:
            break
    t_equations = time.perf_counter()
    if not eqs:
        return [], {"vectors": len(vecs), "used": 0, "nullity": None}
    Mx = Matrix(eqs)
    nullity = phi - Mx.rank()
    point_profile = {} if profile is not None else None
    cands = _integer_points(eqs, rhs, phi, code.H, radius=2, profile=point_profile) if nullity <= 3 else []
    t_end = time.perf_counter()
    if profile is not None:
        profile.update({
            "timings_seconds": {
                "constraint_crt_and_lattice_construction": t_construct - t_start,
                "interpolation_lll": t_reduce - t_construct,
                "short_vector_collection": t_vectors - t_reduce,
                "sigma_linear_equation_build_and_rank": t_equations - t_vectors,
                "integer_point_extraction_total": t_end - t_equations,
                "decode_total": t_end - t_start,
            },
            "integer_point_stages": point_profile,
            "lattice": {"N_bits": N.bit_length(), "N": str(N), "D": D,
                         "embedded_rows": dim, "congruences": n * phi,
                         "crt_aggregated_rows": phi, "W0": W0,
                         "W1": W1, "big": str(big), "n_vectors_cap": n_vectors,
                         "requested_s": s, "coefficient_degree_bound": phi - 1},
            "vectors_count": len(vecs), "accepted_count": used,
            "rejected_inconsistent_count": skipped,
            "rank": phi - nullity,
            "short_vector_heights": [{"A0_max_abs": max(map(abs, v[:phi])),
                                      "other_A_max_abs": max(map(abs, v[phi:]))}
                                     for v in vecs],
            "selected_auxiliary_vectors": selected_vectors,  # internal audit only
        })
    return cands, {"vectors": len(vecs), "used": used, "nullity": nullity}


def _integer_points(eqs, rhs, phi, H, radius=2, profile: dict | None = None):
    """Small integer solutions x (|x_t| < H) of eqs * x = rhs.

    One LLL on the embedding [eqs^T*W | I | 0 ; -rhs^T*W | 0 | 1] gives both a
    reduced basis of the integer kernel and a short particular solution. We then
    enumerate x0 + sum t_i b_i with |t_i| <= radius. This is list recovery when
    the solution space is a line or plane, e.g. when two codewords both agree
    heavily with the received word."""
    from fpylll import LLL, IntegerMatrix
    t_start = time.perf_counter()
    m = len(eqs)
    big = 1 << 200
    dim = phi + 1
    B = IntegerMatrix(dim, m + phi + 1)
    for t in range(phi):
        for r in range(m):
            B[t, r] = eqs[r][t] * big
        B[t, m + t] = 1
    for r in range(m):
        B[phi, r] = -rhs[r] * big
    B[phi, m + phi] = 1
    t_embedding = time.perf_counter()
    LLL.reduction(B)
    t_lll = time.perf_counter()
    rows = [[B[i, m + t] for t in range(phi)] + [B[i, m + phi]]
            for i in range(dim) if not any(B[i, r] for r in range(m))]
    # Euclid on the last column: leave one row with c = gcd (= +-1 if an integer
    # solution exists) and make every other row's c zero (those span the kernel).
    while sum(1 for r in rows if r[-1] != 0) > 1:
        nz = sorted((r for r in rows if r[-1] != 0), key=lambda r: abs(r[-1]))
        piv = nz[0]
        for r in nz[1:]:
            q = r[-1] // piv[-1]
            for t in range(phi + 1):
                r[t] -= q * piv[t]
    kernel = [r[:phi] for r in rows if r[-1] == 0]
    piv = [r for r in rows if r[-1] != 0]
    if not piv or abs(piv[0][-1]) != 1:
        return []
    x0 = [piv[0][-1] * v for v in piv[0][:phi]]
    if kernel:
        K = IntegerMatrix.from_matrix(kernel)
        LLL.reduction(K)
        kernel = [[K[i, t] for t in range(phi)] for i in range(K.nrows)]
        # Babai-style size reduction of x0 against the reduced kernel basis.
        for b in reversed(kernel):
            bb = sum(v * v for v in b)
            if bb:
                c = round(sum(x * v for x, v in zip(x0, b)) / bb)
                x0 = [x - c * v for x, v in zip(x0, b)]
    t_kernel = time.perf_counter()
    out = set()
    import itertools
    for ts in itertools.product(range(-radius, radius + 1), repeat=len(kernel)):
        x = [x0[t] + sum(tt * kv[t] for tt, kv in zip(ts, kernel)) for t in range(phi)]
        if all(abs(v) < H for v in x):
            out.add(tuple(x))
    t_end = time.perf_counter()
    if profile is not None:
        profile.update({
            "embedding_construction": t_embedding - t_start,
            "embedding_lll": t_lll - t_embedding,
            "euclidean_kernel_and_babai": t_kernel - t_lll,
            "bounded_grid_enumeration": t_end - t_kernel,
            "total": t_end - t_start,
            "nullity_integer_kernel": len(kernel),
            "radius": radius,
            "grid_points_examined": (2 * radius + 1) ** len(kernel),
            "bounded_points_found": len(out),
            "embedding_big_bits": big.bit_length(),
            "embedding_rows": dim, "embedding_columns": m + phi + 1,
        })
    return [list(x) for x in out]


def agreement(code: CycloCode, c, received) -> int:
    """Number of blocks where the codeword of c equals the received block."""
    w = code.encode(c)
    return sum(a == b for a, b in zip(w, received))


def _crt(res, mods):
    R, N = 0, 1
    for r, q in zip(res, mods):
        t = ((r - R) * pow(N, -1, q)) % q
        R, N = R + N * t, N * q
    return R % N
