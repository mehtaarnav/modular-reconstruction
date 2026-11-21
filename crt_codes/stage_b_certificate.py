"""Certify Stage B on ONE existing PR #2 Trial-0 instance.

Usage:
  python -m crt_codes.stage_b_certificate --output stage-b-certificate.json

This is NOT a sweep. It uses the exact original PR #2 seed and code:
M=11, n=20, bits=24, coeff_bits=71, code seed=2,
received seed=20261121, trial=0, s=7, agreement>=6.

The original lattice decoder runs once to produce the six interpolation
vectors. A separate exact column-Hermite computation certifies integer
kernel saturation, then exact rational LP and exhaustive integer search.
There is no use of the planted messages to generate candidates.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import random
import time
from fractions import Fraction
from pathlib import Path

from sympy import Matrix

from .exhaustive import exhaustive_decode
from .folded import CycloCode, decode, agreement, _ev
from .integer_affine import (
    column_hermite_solution, gauss_reduce_two, center_particular,
    build_box_inequalities, exact_rank_two_lp, complete_integer_points,
    exact_dot, matmul,
)
from .triad_audit import make_triad

M, N_BLOCKS, BITS, COEFF_BITS, CODE_SEED = 11, 20, 24, 71, 2
RECEIVED_SEED, S, MIN_AGREEMENT = 20261121, 7, 6


def exact_equations(code, vectors, s):
    """Reimplement decode's sigma-linear equations using exact arithmetic."""
    phi = code.phi
    T, rhs = [], []
    for v in vectors:
        A = [v[j * phi:(j + 1) * phi] for j in range(s + 1)]
        columns = []
        for t in range(phi):
            unit = [0] * phi
            unit[t] = 1
            acc = [0] * phi
            for j in range(1, s + 1):
                term = code.mul(A[j], code.sigma(unit, j - 1))
                acc = [u + w for u,w in zip(acc, term)]
            columns.append(acc)
        for i in range(phi):
            T.append([columns[t][i] for t in range(phi)])
            rhs.append(-A[0][i])
    return T, rhs


def certify_selected_auxiliaries(code, vectors, received, s):
    """Exact modular checks plus uniform coefficient-height certificate."""
    phi, H = code.phi, code.H
    pmin = math.prod(sorted(code.primes)[:MIN_AGREEMENT])
    aux_certificates = []
    for v in vectors:
        A = [v[j * phi:(j + 1) * phi] for j in range(s + 1)]
        checks = 0
        for i, p in enumerate(code.primes):
            for j, x in enumerate(code.points[i]):
                value = _ev(A[0], x, p)
                for ell in range(1, s + 1):
                    y = received[i][(j + ell - 1) % phi]
                    value += _ev(A[ell], x, p) * y
                if value % p:
                    raise ArithmeticError("Selected auxiliary violates interpolation")
                checks += 1
        B0 = max(abs(z) for z in A[0])
        Bpos = [max(abs(z) for z in Al) for Al in A[1:]]
        # sigma maps |c_t|<=H-1 into |coeff|<=2(H-1), and
        # each coefficient of mul(A_l,sigma(c)) after Phi_11 reduction
        # contains at most 2*phi such signed products, hence 4*phi.
        U = B0 + 4 * phi * (H - 1) * sum(Bpos)
        aux_certificates.append({
            "interpolation_congruences_checked": checks,
            "B0": B0, "B1_through_B7": Bpos,
            "uniform_coefficient_bound_U": U, "U_bits": U.bit_length(),
            "U_below_Pmin6": U < pmin,
        })
    Umax = max(cert["uniform_coefficient_bound_U"] for cert in aux_certificates)
    if not (Umax < 2**117 < 2**138 <= pmin):
        raise ArithmeticError("Requested exact Stage A inequality not satisfied")
    return {
        "auxiliaries": aux_certificates,
        "maximum_uniform_U": Umax,
        "U_below_2pow117": Umax < 2**117,
        "two_pow_117": str(2**117),
        "two_pow_138": str(2**138),
        "P_min6": pmin,
        "P_min6_bits": pmin.bit_length(),
        "certificate": "U_max < 2^117 < 2^138 <= P_min6",
        "holds_for_every_message_with_at_least_six_full_matching_blocks": True,
    }


def frac_repr(value: Fraction):
    return {"numerator": value.numerator, "denominator": value.denominator}


def run_once():
    stages = {}
    start = time.perf_counter()
    code = CycloCode(M=M, n=N_BLOCKS, bits=BITS,
                     coeff_bits=COEFF_BITS, seed=CODE_SEED)
    stages["code_construction"] = time.perf_counter() - start
    rng = random.Random(RECEIVED_SEED)
    start = time.perf_counter()
    planted, received, supports, realized = make_triad(code, rng, agreements=6)
    stages["received_word_generation"] = time.perf_counter() - start

    # Exactly one LLL/decoder run on the original fixed Trial-0 word.
    profile = {}
    start = time.perf_counter()
    heuristic_candidates, decoder_info = decode(code, received, S, profile=profile)
    stages["original_decoder_total"] = time.perf_counter() - start
    vectors = profile.pop("selected_auxiliary_vectors")
    if len(vectors) != decoder_info["used"]:
        raise ArithmeticError("Selected vector count mismatch")
    if len(vectors) != 6:
        raise ArithmeticError("Expected exactly six accepted auxiliary vectors")

    start = time.perf_counter()
    T, rhs = exact_equations(code, vectors, S)
    stages["exact_sigma_equation_reconstruction"] = time.perf_counter() - start

    start = time.perf_counter()
    stage_a = certify_selected_auxiliaries(code, vectors, received, S)
    stages["stage_A_exact_certificate"] = time.perf_counter() - start

    start = time.perf_counter()
    x0_uncentered, basis_uncentered, V, form = column_hermite_solution(T, rhs)
    stages["unimodular_column_hermite_elimination"] = time.perf_counter() - start
    if form["rank"] != 8 or form["nullity"] != 2:
        raise ArithmeticError("Unexpected extracted rational rank/nullity")

    start = time.perf_counter()
    basis, U, gauss_iterations = gauss_reduce_two(basis_uncentered)
    x0, center_shift = center_particular(x0_uncentered, basis)
    # Update right unimodular transformation with the GL_2(Z) kernel basis change.
    Vred = [row[:] for row in V]
    for i in range(code.phi):
        for j in range(2):
            Vred[i][form["rank"] + j] = sum(
                V[i][form["rank"] + h] * U[j][h] for h in range(2)
            )
    assert basis == [
        [Vred[i][form["rank"] + j] for i in range(code.phi)]
        for j in range(2)
    ]
    assert [exact_dot(row, x0) for row in T] == rhs
    assert all(all(exact_dot(row, v) == 0 for row in T) for v in basis)
    stages["exact_kernel_gauss_reduction_and_centering"] = time.perf_counter() - start

    start = time.perf_counter()
    det_V = int(Matrix(Vred).det(method="domain-ge"))
    transformed = matmul(T, Vred)
    if abs(det_V) != 1:
        raise ArithmeticError("V is not unimodular")
    if any(row[col] for row in transformed for col in range(form["rank"], code.phi)):
        raise ArithmeticError("T V has a nonzero kernel column")
    if int(Matrix(T).rank()) != form["rank"]:
        raise ArithmeticError("Exact rank certificate failed")
    stages["unimodular_determinant_and_rank_verification"] = time.perf_counter() - start

    start = time.perf_counter()
    inequalities = build_box_inequalities(x0, basis, code.H)
    polytope = exact_rank_two_lp(inequalities)
    stages["exact_rational_LP_vertex_enumeration"] = time.perf_counter() - start

    start = time.perf_counter()
    all_integer, count = complete_integer_points(
        x0, basis, code.H, inequalities, polytope["bounds"],
        max_box_cardinality=2_000_000,
    )
    stages["certified_integer_polytope_enumeration"] = time.perf_counter() - start

    if not count["complete"]:
        qualifying = None
        comparison = {"complete": False, "reason": count["reason"]}
    else:
        start = time.perf_counter()
        qualifying = sorted(
            {tuple(c) for _, _, c in all_integer
             if agreement(code, c, received) >= MIN_AGREEMENT}
        )
        stages["exact_agreement_filtering"] = time.perf_counter() - start
        start = time.perf_counter()
        reference, reference_info = exhaustive_decode(
            code, received, MIN_AGREEMENT
        )
        stages["independent_complete_reference"] = time.perf_counter() - start
        reference_set = {tuple(c) for c in reference}
        if set(qualifying) != reference_set:
            raise ArithmeticError("Complete Stage B list disagrees with reference")
        comparison = {
            "complete": True,
            "reference_list": reference,
            "reference_stats": reference_info,
            "stageB_qualifying_list": [list(c) for c in qualifying],
            "matched_reference_exactly": True,
            "qualifying_count": len(qualifying),
            "original_radius_two_decoder_returned": len(heuristic_candidates),
            "original_decoder_qualifying": sum(
                agreement(code, c, received) >= MIN_AGREEMENT
                for c in heuristic_candidates
            ),
            "planted_recovered_for_validation_only": sum(
                tuple(c) in set(qualifying) for c in planted
            ),
        }
    kernel_system = {
        "T": T, "rhs": rhs,
        "x0": x0, "saturated_integer_kernel_basis": basis,
        "V_unimodular": Vred,
        "determinant_V": det_V,
        "rank": form["rank"], "integer_nullity": form["nullity"],
        "column_Hermite_pivot_rows": form["pivot_rows"],
        "column_operations": form["elementary_unimodular_column_operations"],
        "integer_basis_Gauss_transform": U,
        "integer_basis_Gauss_iterations": gauss_iterations,
        "integer_center_shift": list(center_shift),
        "certificates": {
            "T_times_x0_equals_rhs": True,
            "T_times_each_kernel_vector_zero": True,
            "abs_det_V_equals_one": True,
            "last_two_columns_T_times_V_zero": True,
            "rank_T_equals_eight": True,
        },
    }
    exact_lp_certificate = {
        "inequalities": inequalities,
        "minima_continuous": [frac_repr(q) for q in polytope["continuous_min"]],
        "maxima_continuous": [frac_repr(q) for q in polytope["continuous_max"]],
        "certified_coordinate_bounds": polytope["bounds"],
        "feasible_rational_vertices": [
            [frac_repr(q) for q in vertex]
            for vertex in polytope["vertices"]
        ],
        "vertex_count": len(polytope["vertices"]),
        "outer_box_cardinality": polytope["outer_box_cardinality"],
        "continuous_polytope_is_bounded": True,
        "integer_search": count,
        "all_integer_solutions": [
            {"coordinates": [u,v], "coefficient_vector": c}
            for u,v,c in all_integer
        ] if all_integer is not None else None,
    }
    raw_inputs = json.dumps(
        {"T":T,"rhs":rhs},separators=(",",":"),sort_keys=True
    ).encode()
    return {
        "instance": {
            "code": {"M":M,"phi":code.phi,"n":code.n,"bits":BITS,
                     "coeff_bits":COEFF_BITS,"H":code.H,"s":S,
                     "code_seed":CODE_SEED,"primes":code.primes,
                     "g":code.g},
            "received_seed":RECEIVED_SEED,"trial":0,
            "minimum_full_block_agreement":MIN_AGREEMENT,
            "received_blocks":received,
            "support_indices":supports,"planted_agreements":realized,
            "no_planted_messages_used_to_build_any_certificate":True,
            "T_rhs_canonical_SHA256": hashlib.sha256(raw_inputs).hexdigest(),
        },
        "stage_A":stage_a,
        "stage_B_integer_system":kernel_system,
        "stage_B_polytope":exact_lp_certificate,
        "comparison":comparison,
        "timings_seconds":stages,
        "original_decoder_info":decoder_info,
        "original_decoder_profile_timings":profile.get("timings_seconds"),
        "selected_interpolation_vectors":vectors,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path, default=Path("stage-b-certificate.json"))
    args = p.parse_args()
    result = run_once()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    b = result["stage_B_polytope"]
    print(json.dumps({
        "run": result["instance"]["trial"],
        "stage_A": result["stage_A"]["certificate"],
        "stage_A_U_max": str(result["stage_A"]["maximum_uniform_U"]),
        "P_min6": str(result["stage_A"]["P_min6"]),
        "T_rank": result["stage_B_integer_system"]["rank"],
        "integer_nullity": result["stage_B_integer_system"]["integer_nullity"],
        "det_V": result["stage_B_integer_system"]["determinant_V"],
        "integer_particular_x0": result["stage_B_integer_system"]["x0"],
        "saturated_kernel_basis": result["stage_B_integer_system"]["saturated_integer_kernel_basis"],
        "continuous_lp_min": b["minima_continuous"],
        "continuous_lp_max": b["maxima_continuous"],
        "exact_coordinate_bounds": b["certified_coordinate_bounds"],
        "outer_box_cardinality": b["outer_box_cardinality"],
        "enumeration": b["integer_search"],
        "comparison": {
            "complete": result["comparison"]["complete"],
            "qualifying_count": result["comparison"].get("qualifying_count"),
            "matched_reference": result["comparison"].get("matched_reference_exactly"),
            "reference_count": result["comparison"].get("reference_stats",{}).get("list_size"),
        },
        "timings_seconds": result["timings_seconds"],
        "original_decode_profile": result["original_decoder_profile_timings"],
        "certificate_artifact": str(args.output),
    },indent=2,sort_keys=True))


if __name__ == "__main__":
    main()
