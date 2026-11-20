"""Measure where the lattice CRT list-decoder stops, against the three radii
that matter: unique decoding, the Johnson bound, and list-decoding capacity.

Usage:  python3 -m crt_codes.experiment [--trials 20] [--out crt_codes/results/baseline.md]
"""

from __future__ import annotations

import argparse
import math
import random
import time
from pathlib import Path

from .code import CRTCode, corrupt
from .decode import list_decode


def radii(code: CRTCode, k: int) -> dict[str, float]:
    """Error counts (equal-size moduli approximation) at the three radii."""
    n = code.n
    return {
        "unique": (n - k) / 2,
        "johnson": n - math.sqrt(k * n),
        "capacity": n - k,
    }


def corrupt_high_digit(code: CRTCode, word, e, rng, p_of):
    """Errors that keep m mod p but corrupt m mod p**s (only the 'derivative' digits)."""
    word = list(word)
    pos = rng.sample(range(code.n), e)
    for i in pos:
        q, p = code.moduli[i], p_of[i]
        new = (word[i] % p) + p * rng.randrange(1, q // p)
        word[i] = new % q
    return word, pos


def run(n, k, bits, s, h, errors, trials, mode, seed=0):
    code = CRTCode.make(n=n, k=k, bits=bits, s=s, seed=seed)
    p_of = [round(q ** (1 / s)) for q in code.moduli]
    rng = random.Random(seed + 1)
    rows = []
    for e in errors:
        ok = 0
        lists = 0
        t0 = time.time()
        for _ in range(trials):
            m = rng.randrange(code.K)
            w = code.encode(m)
            if mode == "full":
                r, _ = corrupt(code, w, e, rng)
            else:
                r, _ = corrupt_high_digit(code, w, e, rng, p_of)
            log_b = code.agreement_weight(m, r) * 0.999
            S = list_decode(code.moduli, code.K, r, log_b=log_b, h=h)
            ok += m in S
            lists += len(S)
        rows.append((e, ok / trials, lists / trials, (time.time() - t0) / trials))
    return code, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=20)
    ap.add_argument("--out", default="crt_codes/results/baseline.md")
    a = ap.parse_args()

    n, k = 24, 6
    configs = [
        ("prime moduli, s=1, 16-bit primes, full errors", dict(bits=16, s=1, mode="full")),
        ("prime-power moduli p^2, 8-bit primes, full errors", dict(bits=8, s=2, mode="full")),
        ("prime-power moduli p^2, 8-bit primes, high-digit-only errors", dict(bits=8, s=2, mode="high")),
    ]
    errors = list(range(6, 23))
    lines = [f"# CRT list-decoding baseline (n={n}, k={k}, rate {k/n:.2f})", ""]
    for h in (2, 4, 6):
        for title, cfg in configs:
            code, rows = run(n, k, cfg["bits"], cfg["s"], h, errors, a.trials, cfg["mode"])
            rad = radii(code, k)
            lines += [
                f"## {title}, multiplicity h={h}",
                "",
                f"Radii in errors: unique {rad['unique']:.1f}, Johnson {rad['johnson']:.1f}, capacity {rad['capacity']:.1f}",
                "",
                "| errors | success | avg list size | sec/trial |",
                "|---:|---:|---:|---:|",
            ]
            lines += [f"| {e} | {sr:.2f} | {ls:.2f} | {t:.3f} |" for e, sr, ls, t in rows]
            lines.append("")
            last = max((e for e, sr, _, _ in rows if sr >= 0.95), default=None)
            print(f"h={h} | {title}: decodes reliably up to e={last}  "
                  f"(unique {rad['unique']:.1f}, Johnson {rad['johnson']:.1f}, capacity {rad['capacity']:.1f})", flush=True)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text("\n".join(lines))
    print("wrote", a.out)


if __name__ == "__main__":
    main()
