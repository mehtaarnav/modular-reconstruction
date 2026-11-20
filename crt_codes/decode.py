"""Lattice list-decoder for CRT codes (Boneh 2000; Guruswami-Sahai-Sudan 2000),
via Coppersmith / Howgrave-Graham "small roots modulo an unknown divisor".

If the received word agrees with the codeword of m on a set S, then
P_S = prod_{i in S} q_i divides both N and (m - R), where R = CRT(received).
So m is a root of f(x) = x - R modulo an unknown divisor b = P_S of N, with
|m| < K. Howgrave-Graham's lattice, with multiplicity h, recovers every such m
as long as (asymptotically, h -> infinity)

        log K < (log b)**2 / log N,     i.e.   log b > sqrt(log K * log N).

That is the Johnson bound for CRT codes. List-decoding *capacity* would be
log b > log K. The gap between the two is the open territory.
"""

from __future__ import annotations

import math

from fpylll import LLL, IntegerMatrix
from sympy import Poly, symbols

_x = symbols("x")


def _poly_mul(a: list[int], b: list[int]) -> list[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                out[i + j] += ai * bj
    return out


def list_decode(moduli, K: int, received: list[int], log_b: float, h: int) -> set[int]:
    """All m in [0, K) whose agreement product with `received` is (heuristically)
    at least exp(log_b), using multiplicity h. Returns the candidates found; the
    caller filters by actual agreement."""
    N = math.prod(moduli)
    R = 0
    Nacc = 1
    for r, q in zip(received, moduli):
        t = ((r - R) * pow(Nacc, -1, q)) % q
        R, Nacc = R + Nacc * t, Nacc * q
    R %= N

    beta = log_b / math.log(N)
    beta = min(max(beta, 1e-6), 1.0)
    # Number of extra shifts x^j f^h (May's choice t ~ h(1/beta - 1)).
    t = max(1, math.ceil(h * (1 / beta - 1)))
    d = h + t
    X = K  # bound on the root

    f = [-R, 1]  # x - R, low degree first
    fpow = [[1]]
    for _ in range(h):
        fpow.append(_poly_mul(fpow[-1], f))

    rows = []
    for i in range(h):  # N^(h-i) f^i
        rows.append([c * N ** (h - i) for c in fpow[i]])
    for j in range(t):  # x^j f^h
        rows.append([0] * j + fpow[h])
    M = IntegerMatrix(d, d)
    for r_idx, row in enumerate(rows):
        for c_idx, c in enumerate(row):
            M[r_idx, c_idx] = c * X**c_idx
    LLL.reduction(M)

    found: set[int] = set()
    for r_idx in range(min(3, d)):
        coeffs = [M[r_idx, c] for c in range(d)]
        if not any(coeffs):
            continue
        poly = []
        for c_idx, c in enumerate(coeffs):
            q, rem = divmod(c, X**c_idx)
            if rem:
                break
            poly.append(q)
        else:
            P = Poly(list(reversed(poly)), _x)
            if P.degree() < 1:
                continue
            for root in P.ground_roots():
                if root.is_integer and 0 <= int(root) < K:
                    found.add(int(root))
    return found
