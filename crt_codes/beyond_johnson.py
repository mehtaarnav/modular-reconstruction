"""Decisive tests for the Galois-folded cyclotomic decoder (crt_codes/folded.py).

1. Random errors past the Johnson radius.
2. Worst-case: a received word that mixes two codewords alpha and beta block by
   block, each with agreement *below* the Johnson threshold n*sqrt(R) but above
   the predicted threshold n*(1/(s+1) + R). Both must come back in the list.

Usage: python3 -m crt_codes.beyond_johnson [--trials 3]
"""
import argparse, math, random, time
from .folded import CycloCode, decode, agreement


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--trials", type=int, default=3); a = ap.parse_args()
    M, n, bits, cb, s = 11, 20, 24, 71, 7
    code = CycloCode(M=M, n=n, bits=bits, coeff_bits=cb, seed=2)
    rng = random.Random(3); R = code.rate()
    print(f"M={M} phi={code.phi} n={n} rate={R:.3f} s={s}")
    print(f"  Johnson: errors < {n*(1-math.sqrt(R)):.2f}  (agreement > {n*math.sqrt(R):.2f})")
    print(f"  predicted: errors < {n*(s/(s+1))*(1-R):.2f}  (agreement > {n*(1/(s+1)+R*s/(s+1)):.2f})")
    print(f"  capacity: errors < {n*(1-R):.2f}", flush=True)
    for e in (15, 16):
        ok = 0
        for _ in range(a.trials):
            c = code.random_message(rng); r = code.corrupt(code.encode(c), e, rng)
            cands, _ = decode(code, r, s); ok += c in cands
        print(f"random errors e={e}: {ok}/{a.trials}", flush=True)
    for aa, ab in ((7, 7), (6, 6), (7, 6)):
        res = []
        for _ in range(a.trials):
            al, be = code.random_message(rng), code.random_message(rng)
            wa, wb = code.encode(al), code.encode(be)
            idx = list(range(n)); rng.shuffle(idx); r = [None] * n
            for i in idx[:aa]: r[i] = wa[i]
            for i in idx[aa:aa + ab]: r[i] = wb[i]
            for i in idx[aa + ab:]: r[i] = [rng.randrange(code.primes[i]) for _ in range(code.phi)]
            t0 = time.time(); cands, info = decode(code, r, s)
            thr = min(aa, ab)
            good = [x for x in cands if agreement(code, x, r) >= thr]
            res.append(("A" if al in good else "-") + ("B" if be in good else "-") + f"/list{len(good)}/null{info['nullity']}")
        print(f"two-codeword mixture, agreements alpha={aa} beta={ab} (Johnson needs >{n*math.sqrt(R):.2f}): {' '.join(res)}", flush=True)


if __name__ == "__main__":
    main()
