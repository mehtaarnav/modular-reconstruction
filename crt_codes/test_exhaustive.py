"""Independent validation of full-block interpolation and subset completeness.

These tests need only the Python standard library (not sympy/fpylll).
"""
import random

import pytest

from crt_codes.exhaustive import (
    block_coefficients, exhaustive_decode,
    interpolate_block, min_reconstruction_blocks,
)


class ToyCode:
    phi = 3
    H = 20
    primes = [11, 13, 17, 19, 23, 29, 31]
    points = [[1, 2, 3] for _ in primes]
    n = len(primes)

    def encode(self, c):
        return [[sum(v * pow(x, t, p) for t, v in enumerate(c)) % p
                 for x in pts]
                for p, pts in zip(self.primes, self.points)]


def test_vandermonde_round_trip():
    code = ToyCode()
    c = [8, -14, 19]
    assert block_coefficients(code, code.encode(c)) == [
        [x % p for x in c] for p in code.primes
    ]


def test_reconstruction_requires_smallest_product_bound():
    code = ToyCode()
    assert min_reconstruction_blocks(code) == 2
    assert 11 * 13 > 2 * (code.H - 1)


def test_exact_two_message_list_and_no_spurious():
    code = ToyCode()
    a, b = [2, -3, 5], [-7, 11, -2]
    wa, wb = code.encode(a), code.encode(b)
    received = wa[:3] + wb[3:6] + [[3, 4, 5]]
    found, stats = exhaustive_decode(code, received, 3)
    assert {tuple(x) for x in found} == {tuple(a), tuple(b)}
    assert stats["complete"]
    assert stats["subsets"] == 21


def test_threshold_below_reconstruction_not_claimed_complete():
    code = ToyCode()
    with pytest.raises(ValueError, match="Completeness"):
        exhaustive_decode(code, code.encode([0, 0, 0]), 1)


def test_random_single_message():
    code = ToyCode()
    rng = random.Random(42)
    for _ in range(10):
        c = [rng.randrange(-19, 20) for _ in range(3)]
        received = code.encode(c)
        found, _ = exhaustive_decode(code, received, code.n)
        assert found == [c]
