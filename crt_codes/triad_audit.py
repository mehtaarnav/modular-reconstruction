"""Adversarial three-codeword block-mixture benchmark for Galois-folded CRT.

Usage:
    python -m crt_codes.triad_audit --trials 3 --seed 20261121
    python -m crt_codes.triad_audit --trials 3 --skip-lattice

For each trial, 6 blocks come from each of 3 independent codewords.
The remaining 2 blocks contain independent random symbols. We use a
complete exhaustive subset decoder as ground truth for >=6 agreements.
No original message is passed into either decoder.

Exhaustive reference is exponential in k; this is a fixed-size falsification
harness, NOT an efficient alternative to folded list decoding.
"""
from __future__ import annotations

import argparse
import json
import math
import random
import time
from pathlib import Path

from .exhaustive import exhaustive_decode, min_reconstruction_blocks


def make_triad(code, rng, agreements=6):
    """Construct three planted codewords with disjoint supporting blocks."""
    if 3 * agreements > code.n:
        raise ValueError("Not enough blocks for disjoint planted supports")
    planted = []
    while len(planted) < 3:
        x = code.random_message(rng)
        if x not in planted:
            planted.append(x)
    words = [code.encode(x) for x in planted]
    positions = list(range(code.n))
    rng.shuffle(positions)
    received = [None] * code.n
    support = []
    for j in range(3):
        subset = positions[j * agreements:(j + 1) * agreements]
        support.append(sorted(subset))
        for i in subset:
            received[i] = list(words[j][i])
    for i in positions[3 * agreements:]:
        received[i] = [rng.randrange(code.primes[i]) for _ in range(code.phi)]
    realized = [sum(a == b for a, b in zip(word, received)) for word in words]
    return planted, received, support, realized


def run(trials=3, seed=20261121, skip_lattice=False, n=20,
        M=11, bits=24, coeff_bits=71, s=7):
    # CycloCode currently imports fpylll at module import time. Keep the
    # pure-Python exhaustive module independent so its tests still work.
    from .folded import CycloCode, decode, agreement
    code = CycloCode(M=M, n=n, bits=bits, coeff_bits=coeff_bits, seed=2)
    rate = code.rate()
    johnson = n * math.sqrt(rate)
    predicted_a = n * (1 / (s + 1) + rate * s / (s + 1))
    if min_reconstruction_blocks(code) > 6:
        raise ValueError("Need more than six full blocks to reconstruct")
    rows = []
    rng = random.Random(seed)
    for trial in range(trials):
        planted, received, support, realized = make_triad(code, rng)
        start = time.perf_counter()
        reference, audit = exhaustive_decode(code, received, 6)
        ref_secs = time.perf_counter() - start
        ref_set = {tuple(v) for v in reference}
        planted_set = {tuple(v) for v in planted}
        assert planted_set <= ref_set, "The complete reference missed a planted message"
        row = {
            "trial": trial, "seed": seed, "planted_agreement": realized,
            "supports": support, "reference": {
                **audit, "planted_found": sum(tuple(v) in ref_set for v in planted),
                "seconds": round(ref_secs, 3),
            }, "lattice": None,
        }
        if not skip_lattice:
            start = time.perf_counter()
            raw, info = decode(code, received, s)
            secs = time.perf_counter() - start
            verified = {tuple(v) for v in raw
                        if agreement(code, v, received) >= 6}
            row["lattice"] = {
                **info, "raw_candidates": len(raw),
                "verified_candidates": len(verified),
                "planted_found": sum(tuple(v) in verified for v in planted),
                "missed_from_reference": len(ref_set - verified),
                "spurious_vs_reference": len(verified - ref_set),
                "seconds": round(secs, 3),
            }
        rows.append(row)
    return {
        "parameters": {
            "M": M, "n": n, "phi": code.phi, "bits": bits,
            "coeff_bits": coeff_bits, "s": s, "rate": rate,
            "seed": seed, "trials": trials,
            "agreement_threshold": 6, "johnson_agreement": johnson,
            "heuristic_agreement": predicted_a,
            "skipped_lattice": skip_lattice,
        }, "trials": rows,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--trials", type=int, default=3)
    p.add_argument("--seed", type=int, default=20261121)
    p.add_argument("--skip-lattice", action="store_true")
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    result = run(args.trials, args.seed, args.skip_lattice)
    out = json.dumps(result, sort_keys=True, indent=2)
    print(out)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(out + "\n")


if __name__ == "__main__":
    main()
