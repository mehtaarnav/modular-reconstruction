"""Independent, exact CRT-subset reference decoder for cyclotomic block codes.

This is a *combinatorial* reference, not an efficient asymptotic decoder.
Each full cyclotomic block is a Vandermonde transform of the coefficient
residues modulo its rational prime. Invert this transform, enumerate enough
blocks for unique coefficient reconstruction, and verify candidates against
ALL received blocks. Unlike a lattice heuristic, it does not need the
transmitted message, a guessed message, or an oracle-supplied agreement set.
"""
from __future__ import annotations

import itertools
import math


def interpolate_block(points: list[int], values: list[int], prime: int) -> list[int]:
    """Recover power-basis coefficients over F_p from all orbit evaluations."""
    d = len(points)
    if len(values) != d or len(set(points)) != d:
        raise ValueError("Expected a full block at distinct evaluation points")
    a = [[pow(x, t, prime) for t in range(d)] + [y % prime]
         for x, y in zip(points, values)]
    for c in range(d):
        pivot = next((i for i in range(c, d) if a[i][c] % prime), None)
        if pivot is None:
            raise ValueError("Block Vandermonde matrix is singular")
        a[c], a[pivot] = a[pivot], a[c]
        inv = pow(a[c][c], -1, prime)
        a[c] = [(z * inv) % prime for z in a[c]]
        for i in range(d):
            if i == c:
                continue
            t = a[i][c]
            if t:
                a[i] = [(u - t * v) % prime for u, v in zip(a[i], a[c])]
    return [a[i][-1] for i in range(d)]


def min_reconstruction_blocks(code) -> int:
    """Minimum k for *every* k-subset to reconstruct |c_t| < H uniquely."""
    product = 1
    for k, p in enumerate(sorted(code.primes), 1):
        product *= p
        if product > 2 * (code.H - 1):
            return k
    raise ValueError("Combined moduli do not identify every message coefficient")


def block_coefficients(code, received) -> list[list[int]]:
    if len(received) != code.n:
        raise ValueError("Incorrect block count")
    return [interpolate_block(code.points[i], received[i], code.primes[i])
            for i in range(code.n)]


def exhaustive_decode(code, received, min_agreement: int):
    """Return all valid messages with agreement >= min_agreement and audit stats.

    Completeness follows if min_agreement >= k, where k is the number of
    full blocks needed to reconstruct a message uniquely, because at least
    one enumerated k-subset then consists entirely of agreeing blocks.

    Complexity is O(binomial(n,k)); only use as a small-instance reference.
    """
    if not 0 <= min_agreement <= code.n:
        raise ValueError("Invalid agreement threshold")
    k = min_reconstruction_blocks(code)
    if min_agreement < k:
        raise ValueError("Completeness cannot be guaranteed at this threshold")

    residue = block_coefficients(code, received)
    candidates = set()
    checked = set()
    subsets = 0
    for ids in itertools.combinations(range(code.n), k):
        subsets += 1
        moduli = [code.primes[i] for i in ids]
        N = math.prod(moduli)
        weights = [(N // p) * pow(N // p, -1, p) for p in moduli]
        coeff = []
        for t in range(code.phi):
            x = sum(residue[i][t] * w for i, w in zip(ids, weights)) % N
            if x > N // 2:
                x -= N
            if not -code.H < x < code.H:
                break
            coeff.append(x)
        if len(coeff) != code.phi:
            continue
        key = tuple(coeff)
        if key in checked:
            continue
        checked.add(key)
        encoded = code.encode(coeff)
        agreement = sum(a == b for a, b in zip(encoded, received))
        if agreement >= min_agreement:
            candidates.add(key)

    return [list(c) for c in sorted(candidates)], {
        "complete": True, "required_blocks": k,
        "subsets": subsets, "checked": len(checked),
        "list_size": len(candidates),
    }
