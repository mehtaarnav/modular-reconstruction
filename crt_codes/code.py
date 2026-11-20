"""Chinese Remainder (CRT) codes and an error channel.

A CRT code is the integer analogue of a Reed-Solomon code: the message is an
integer 0 <= m < K and the codeword is (m mod q_1, ..., m mod q_n) for pairwise
coprime moduli q_i. With q_i = p_i**s (prime powers) one gets the integer
analogue of a *multiplicity code*: knowing m mod p**s is the integer counterpart
of knowing f, f', ..., f^(s-1) at a point (f mod (x - a)**s).
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass

from sympy import nextprime


def primes_near(n: int, bits: int, seed: int = 0) -> list[int]:
    """n distinct primes of roughly `bits` bits, in increasing order."""
    rng = random.Random(seed)
    start = (1 << (bits - 1)) + rng.randrange(1 << (bits - 2))
    out, p = [], start
    while len(out) < n:
        p = nextprime(p)
        out.append(p)
    return out


def crt(residues: list[int], moduli: list[int]) -> tuple[int, int]:
    """Return (R, N) with R = residues[i] mod moduli[i] for all i and N = prod(moduli)."""
    R, N = 0, 1
    for r, q in zip(residues, moduli):
        # Solve R + N*t = r (mod q).
        t = ((r - R) * pow(N, -1, q)) % q
        R, N = R + N * t, N * q
    return R % N, N


@dataclass(frozen=True)
class CRTCode:
    moduli: tuple[int, ...]
    K: int  # messages are integers in [0, K)

    @staticmethod
    def make(n: int, k: int, bits: int, s: int = 1, seed: int = 0) -> "CRTCode":
        """n moduli p_i**s with p_i ~ 2**bits; K = product of the k smallest moduli
        (so "rate" is about k/n, as for Reed-Solomon)."""
        ps = primes_near(n, bits, seed)
        moduli = tuple(p**s for p in ps)
        K = math.prod(moduli[:k])
        return CRTCode(moduli, K)

    @property
    def n(self) -> int:
        return len(self.moduli)

    @property
    def N(self) -> int:
        return math.prod(self.moduli)

    def encode(self, m: int) -> list[int]:
        assert 0 <= m < self.K
        return [m % q for q in self.moduli]

    def agreement_weight(self, m: int, received: list[int]) -> float:
        """log of prod over i of gcd(m - r_i, q_i): full agreement at i contributes
        log q_i, agreement only mod p**j (j < s) contributes j*log p."""
        return sum(math.log(math.gcd(m - r, q)) for r, q in zip(received, self.moduli))


def corrupt(code: CRTCode, word: list[int], e: int, rng: random.Random) -> tuple[list[int], list[int]]:
    """Replace e random coordinates by uniformly random *different* symbols."""
    word = list(word)
    pos = rng.sample(range(code.n), e)
    for i in pos:
        q = code.moduli[i]
        v = rng.randrange(q - 1)
        word[i] = v if v < word[i] else v + 1
    return word, pos
