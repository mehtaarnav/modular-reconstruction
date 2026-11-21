import math
import random

from crt_codes.code import CRTCode, corrupt, crt
from crt_codes.decode import list_decode


def test_crt_roundtrip():
    code = CRTCode.make(n=8, k=3, bits=12, seed=3)
    m = 123456789 % code.K
    R, N = crt(code.encode(m), list(code.moduli))
    assert R == m and N == code.N


def test_decodes_beyond_unique_radius():
    code = CRTCode.make(n=20, k=5, bits=16, seed=1)
    rng = random.Random(0)
    m = rng.randrange(code.K)
    r, _ = corrupt(code, code.encode(m), 9, rng)  # unique radius is 7.5
    assert m in list_decode(code.moduli, code.K, r, code.agreement_weight(m, r) * 0.999, h=4)


def test_prime_power_partial_agreement_counts():
    code = CRTCode.make(n=6, k=2, bits=8, s=2, seed=5)
    m = 4242 % code.K
    w = code.encode(m)
    p = round(math.sqrt(code.moduli[0]))
    w[0] = (w[0] + p) % code.moduli[0]  # still agrees mod p, not mod p^2
    full = sum(math.log(q) for q in code.moduli)
    assert abs(code.agreement_weight(m, w) - (full - math.log(p))) < 1e-9


def test_folded_decodes_past_johnson():
    """Q(zeta_11), n=20: Johnson radius is 12.09 errors; decode 13. (~30 s)"""
    import math
    import random
    from crt_codes.folded import CycloCode, decode
    code = CycloCode(M=11, n=20, bits=24, coeff_bits=71, seed=2)
    rng = random.Random(1)
    e = 13
    assert e > code.n * (1 - math.sqrt(code.rate()))
    c = code.random_message(rng)
    r = code.corrupt(code.encode(c), e, rng)
    cands, _ = decode(code, r, 7)
    assert c in cands


def test_folded_two_codeword_mixture():
    import random
    from crt_codes.folded import CycloCode, decode
    code = CycloCode(M=5, n=16, bits=24, coeff_bits=50, seed=2)
    rng = random.Random(7)
    al, be = code.random_message(rng), code.random_message(rng)
    wa, wb = code.encode(al), code.encode(be)
    r = wa[:8] + wb[8:]
    cands, _ = decode(code, r, 3)
    assert al in cands and be in cands
