"""11-dimensional scalar affine-relation comparator, PR #2 Trial 0 only.

No new received word, parameter sweep, decoding baseline, or error regime.
No planted messages or reference candidates enter relation selection.

Received block i is inverted by the 10x10 Vandermonde transform to obtain
r_(i,j) = c_j (mod p_i) for each coefficient stream. For each j, CRT gives
R_j (mod N), N=prod_i p_i. The 11-dimensional congruence lattice is

  L = {(a0, H*a1, ..., H*a10) in Z^11:
             a0 + sum_j a_j*R_j == 0 (mod N)}

with explicit row generators (N,0,...,0) and
(-R_j,0,...,H in j-th coordinate,...,0), determinant N*H^10.
For each LLL-reduced basis row, undo the weights exactly, verify the
congruence, and retain iff U=|a0|+(H-1)*sum_j|a_j| < P_min,6.

If a message c agrees in at least 6 whole blocks, then for each selected
relation a0 + sum_j a_j*c_j is divisible by those 6 primes. Its absolute
value is at most U < P_min,6, hence is zero over Z. Thus the selected
relations yield certified affine linear equations on ALL qualifying messages.

The integer affine solution set is computed by a unimodular column-Hermite
reduction, then (when nullity=2) the exact-rational rank-two LP and
nontruncated integer-point enumeration from the previous Stage B audit.
If enumeration is too large, the report is explicit about noncompletion.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import random
import time
from pathlib import Path

from fpylll import IntegerMatrix, LLL
from sympy import Matrix

from .exhaustive import block_coefficients, exhaustive_decode
from .folded import CycloCode, agreement, _crt
from .integer_affine import (
    column_hermite_solution, gauss_reduce_two, center_particular,
    build_box_inequalities, exact_rank_two_lp, complete_integer_points,
    exact_dot, matmul,
)
from .stage_b_certificate import (
    M, N_BLOCKS, BITS, COEFF_BITS, CODE_SEED, RECEIVED_SEED,
    S, MIN_AGREEMENT, run_once as galois_stage_b_run_once,
)
from .triad_audit import make_triad


def scalar_lattice(code, received):
    """Return exact 11D weighted congruence lattice and inverse-CRT data."""
    if code.phi != 10 or code.n != 20:
        raise ValueError("This comparator is restricted to existing Trial-0 parameters")
    modular_coeffs = block_coefficients(code, received)
    N = math.prod(code.primes)
    residues = [
        _crt([modular_coeffs[i][j] for i in range(code.n)], code.primes)
        for j in range(code.phi)
    ]
    H = code.H
    B = IntegerMatrix(code.phi + 1, code.phi + 1)
    B[0, 0] = N
    for j, Rj in enumerate(residues):
        B[j+1, 0] = -Rj
        B[j+1, j+1] = H
    return B, residues, modular_coeffs, N


def select_scalar_relations(B, residues, modular_coeffs, code, N):
    """Check EVERY reduced basis row before choosing by the exact U bound."""
    pmin6 = math.prod(sorted(code.primes)[:MIN_AGREEMENT])
    H = code.H
    all_relations = []
    accepted = []
    for i in range(B.nrows):
        weighted = [int(B[i, j]) for j in range(B.ncols)]
        if any(weighted[j+1] % H for j in range(code.phi)):
            raise ArithmeticError("Reduced vector has nonintegral unweighted coefficients")
        a0 = weighted[0]
        a = [weighted[j+1] // H for j in range(code.phi)]
        global_mod = (a0 + exact_dot(a, residues)) % N
        if global_mod != 0:
            raise ArithmeticError("Reduced vector fails the exact CRT congruence")
        modular_checks = sum(
            (a0 + exact_dot(a, modular_coeffs[k])) % code.primes[k] == 0
            for k in range(code.n)
        )
        if modular_checks != code.n:
            raise ArithmeticError("Reduced vector fails a blockwise congruence")
        U = abs(a0) + (H - 1) * sum(abs(z) for z in a)
        qualifying = U < pmin6
        entry = {
            "reduced_row": i,
            "a0": a0, "a": a,
            "U_exact": U,
            "passes_exact_CRT_congruence": True,
            "passes_all_blockwise_congruences": True,
            "U_below_P_min6": qualifying,
            "weighted_L2_norm_sq": sum(z*z for z in weighted),
        }
        all_relations.append(entry)
        if qualifying:
            accepted.append(entry)
    return all_relations, accepted, pmin6


def exact_affine_extraction(T, rhs, H, max_outer_box=2_000_000):
    """Fully saturate the integer kernel; enumerate if certified tractable.

    No silent radius cutoff. A large or unsupported polytope is explicitly
    reported as unfinished with its exact known algebraic structure.
    """
    if not T:
        return {
            "finished": False,
            "reason": "No certified scalar relations: full Z^10 is unconstrained",
            "rank": 0, "integer_nullity": 10,
        }, None
    start = time.perf_counter()
    x0raw, basis_raw, V, form = column_hermite_solution(T, rhs)
    times = {"unimodular_column_elimination": time.perf_counter() - start}
    rank = form["rank"]
    nullity = form["nullity"]
    start = time.perf_counter()
    Vdet = int(Matrix(V).det(method="domain-ge"))
    if abs(Vdet) != 1:
        raise ArithmeticError("Scalar integer kernel is not saturated")
    if int(Matrix(T).rank()) != rank:
        raise ArithmeticError("Incorrect scalar relation rank")
    if matmul(T, V) != [
        [sum(T[i][k]*V[k][j] for k in range(len(V))) for j in range(len(V))]
        for i in range(len(T))
    ]:
        raise ArithmeticError("Matrix multiplication error")
    if any(z for row in matmul(T, V) for z in row[rank:]):
        raise ArithmeticError("Scalar T*V kernel columns nonzero")
    times["unimodular_det_and_rank_check"] = time.perf_counter() - start
    result = {
        "finished": False,
        "T": T, "rhs": rhs, "rank": rank,
        "integer_nullity": nullity,
        "x0_uncentered": x0raw, "saturated_integer_kernel_basis": basis_raw,
        "V_unimodular": V, "determinant_V": Vdet,
        "unimodular_column_operations": form["elementary_unimodular_column_operations"],
        "timings_seconds": times,
        "certified_saturation": True,
    }
    if nullity != 2:
        result["reason"] = (
            "Exact integer parameterization complete; this audit's exact "
            "rational polytope enumerator handles only rank-two kernels, "
            f"but integer nullity is {nullity}. No truncation performed."
        )
        return result, None

    start = time.perf_counter()
    basis, Uchange, gauss_iterations = gauss_reduce_two(basis_raw)
    x0, translate = center_particular(x0raw, basis)
    if [exact_dot(row, x0) for row in T] != rhs:
        raise ArithmeticError("Centered particular solution is not exact")
    if any(any(exact_dot(row,b) for row in T) for b in basis):
        raise ArithmeticError("Reduced kernel is not exact")
    if abs(int(Matrix(Uchange).det())) != 1:
        raise ArithmeticError("Kernel change of basis not unimodular")
    times["exact_gauss_reduction_and_centering"] = time.perf_counter() - start
    result.update({
        "x0": x0, "saturated_integer_kernel_basis": basis,
        "kernel_unimodular_change": Uchange,
        "kernel_center_shift": list(translate),
        "gauss_iterations": gauss_iterations,
    })
    start = time.perf_counter()
    inequalities = build_box_inequalities(x0, basis, H)
    poly = exact_rank_two_lp(inequalities)
    times["exact_rational_LP"] = time.perf_counter() - start
    bounds = poly["bounds"]
    result.update({
        "polytope_inequalities": inequalities,
        "LP_bounds_exact_fraction": {
            "min": [
                {"num": x.numerator, "den": x.denominator}
                for x in poly["continuous_min"]
            ],
            "max": [
                {"num": x.numerator, "den": x.denominator}
                for x in poly["continuous_max"]
            ],
        },
        "certified_integer_coordinate_bounds": bounds,
        "LP_feasible_vertex_count": len(poly["vertices"]),
        "outer_box_size": poly["outer_box_cardinality"],
    })
    start = time.perf_counter()
    points, count = complete_integer_points(
        x0, basis, H, inequalities, bounds,
        max_box_cardinality=max_outer_box,
    )
    times["integer_polytope_enumeration"] = time.perf_counter() - start
    result["enumeration"] = count
    if points is None:
        result["reason"] = (
            "Certified coordinate box has too many integer pairs for "
            "complete explicit enumeration on this run; no partial list "
            "is represented as complete."
        )
        return result, None
    result["finished"] = True
    result["all_integer_polytope_points"] = [
        {"coordinates": [t1, t2], "message": c}
        for t1, t2, c in points
    ]
    return result, [c for _, _, c in points]


def compare_affine_rational_spaces(T_scalar, rhs_scalar, T_galois, rhs_galois):
    """Exact algebraic test for equality or containment of affine Q-spaces."""
    if not T_scalar:
        return {"scalar_rank": 0, "galois_rank": int(Matrix(T_galois).rank()),
                "equal": False, "relationship": "scalar constraints absent"}
    Msc, Mg = Matrix(T_scalar), Matrix(T_galois)
    bsc, bg = Matrix(rhs_scalar), Matrix(rhs_galois)
    rank_s = int(Msc.rank())
    rank_g = int(Mg.rank())
    consistent_s = int(Msc.row_join(bsc).rank()) == rank_s
    consistent_g = int(Mg.row_join(bg).rank()) == rank_g
    stacked = Msc.col_join(Mg)
    stacked_augmented = stacked.row_join(bsc.col_join(bg))
    rank_union = int(stacked.rank())
    consistent_union = int(stacked_augmented.rank()) == rank_union
    same = consistent_s and consistent_g and consistent_union and rank_s == rank_g == rank_union
    if not (consistent_s and consistent_g):
        relation = "one affine system is inconsistent"
    elif not consistent_union:
        relation = "disjoint affine spaces"
    elif same:
        relation = "equal over Q (hence equal admissible integer slices)"
    elif rank_s == rank_union:
        relation = "scalar affine space is a proper subset of Galois affine space"
    elif rank_g == rank_union:
        relation = "Galois affine space is a proper subset of scalar affine space"
    else:
        relation = "affine spaces intersect, but neither contains the other"
    return {
        "scalar_rank": rank_s,
        "galois_rank": rank_g,
        "stacked_rank": rank_union,
        "stacked_augmented_rank": int(stacked_augmented.rank()),
        "both_individually_consistent": consistent_s and consistent_g,
        "jointly_consistent": consistent_union,
        "equal": same,
        "relationship": relation,
        "intersection_rational_dimension": 10 - rank_union if consistent_union else None,
    }


def run():
    times = {}
    start = time.perf_counter()
    code = CycloCode(M=M, n=N_BLOCKS, bits=BITS, coeff_bits=COEFF_BITS,
                     seed=CODE_SEED)
    times["code_setup"] = time.perf_counter() - start
    rng = random.Random(RECEIVED_SEED)
    start = time.perf_counter()
    _unused_planted, received, _unused_support, _unused_realized = make_triad(
        code, rng, agreements=6
    )
    times["fixed_received_word_generation"] = time.perf_counter() - start

    start = time.perf_counter()
    B, residues, modular_coeffs, N = scalar_lattice(code, received)
    times["inverse_vandermonde_CRT_lattice_construction"] = time.perf_counter() - start
    expected_det_abs = N * (code.H ** code.phi)

    start = time.perf_counter()
    LLL.reduction(B)
    times["scalar_LLL_11x11"] = time.perf_counter() - start
    # Every LLL row is still in the explicit 11-D weighted congruence lattice.
    start = time.perf_counter()
    all_relations, accepted, pmin6 = select_scalar_relations(
        B, residues, modular_coeffs, code, N
    )
    times["scalar_exact_relation_and_height_verification"] = time.perf_counter() - start
    if expected_det_abs <= 0:
        raise ArithmeticError("Invalid lattice determinant")
    scalar_T = [v["a"] for v in accepted]
    scalar_rhs = [-v["a0"] for v in accepted]
    start = time.perf_counter()
    extraction, candidate_points = exact_affine_extraction(
        scalar_T, scalar_rhs, code.H
    )
    times["scalar_integer_parameterization_LP_enumeration"] = time.perf_counter() - start

    scalar_qualifying = None
    if candidate_points is not None:
        start = time.perf_counter()
        scalar_qualifying = sorted({
            tuple(c) for c in candidate_points
            if agreement(code, c, received) >= MIN_AGREEMENT
        })
        times["scalar_exact_whole_block_agreement_filter"] = time.perf_counter() - start

    # The independent reference is used ONLY AFTER the relations are selected
    # and the scalar integer solution set / post-filter candidates are fixed.
    start = time.perf_counter()
    reference, reference_info = exhaustive_decode(code, received, MIN_AGREEMENT)
    times["reference_decode_validation_only"] = time.perf_counter() - start
    ref_set = {tuple(c) for c in reference}
    if scalar_qualifying is not None and set(scalar_qualifying) != ref_set:
        raise ArithmeticError(
            "Scalar complete-parameterization list differs from exact reference"
        )

    # The earlier exact Galois Stage-B certificate is independently rerun
    # ON THE SAME fixed received word, for same-run stage-level comparison.
    start = time.perf_counter()
    galois_cert = galois_stage_b_run_once()
    times["galois_complete_certificate_same_runner"] = time.perf_counter() - start
    if galois_cert["instance"]["received_blocks"] != received:
        raise ArithmeticError("Galois comparator did not use identical received word")
    gsys = galois_cert["stage_B_integer_system"]
    affine_comparison = compare_affine_rational_spaces(
        scalar_T, scalar_rhs, gsys["T"], gsys["rhs"]
    )
    galois_list = {
        tuple(c) for c in galois_cert["comparison"].get("stageB_qualifying_list", [])
    }
    if galois_list != ref_set:
        raise ArithmeticError("Previously certified Galois list disagrees with reference")
    if scalar_qualifying is not None and set(scalar_qualifying) != galois_list:
        raise ArithmeticError("Scalar list differs from Galois list")
    input_fingerprint = hashlib.sha256(
        json.dumps(received, sort_keys=True, separators=(",",":")).encode()
    ).hexdigest()
    relation_fingerprint = hashlib.sha256(
        json.dumps({"T": scalar_T, "rhs":scalar_rhs},
                   sort_keys=True,separators=(",",":")).encode()
    ).hexdigest()
    return {
        "scope": {
            "instance": "PR #2 Trial 0 only",
            "M": M, "phi": code.phi, "n":code.n,
            "code_seed": CODE_SEED, "received_seed":RECEIVED_SEED,
            "bits": BITS, "coeff_bits":COEFF_BITS,
            "H": code.H, "s_galois": S, "minimum_agreement":MIN_AGREEMENT,
            "primes":code.primes,
            "received_blocks":received,
            "received_SHA256":input_fingerprint,
            "planted_messages_not_used_to_select_relations":True,
            "reference_list_not_used_to_select_relations":True,
        },
        "scalar": {
            "lattice_dimension":11,
            "weight_vector":[1]+[code.H]*code.phi,
            "N":N, "N_bits":N.bit_length(),
            "lattice_determinant_abs":expected_det_abs,
            "residue_coeffs_after_block_vandermonde":modular_coeffs,
            "crt_received_coefficients":residues,
            "all_eleven_reduced_basis_relations":all_relations,
            "retained_height_certified_relations":accepted,
            "retained_relation_count":len(accepted),
            "P_min6":pmin6,
            "relations_integer_matrix_SHA256":relation_fingerprint,
            "rank_Q":extraction["rank"],
            "nullity_Q":extraction["integer_nullity"],
            "exact_integer_extraction":extraction,
            "qualifying_list": [list(c) for c in scalar_qualifying]
                if scalar_qualifying is not None else None,
            "qualifying_count":len(scalar_qualifying)
                if scalar_qualifying is not None else None,
            "qualifying_list_equal_complete_reference":
                set(scalar_qualifying) == ref_set
                if scalar_qualifying is not None else None,
        },
        "galois": {
            "same_received_SHA256": input_fingerprint,
            "original_decoder_info":galois_cert["original_decoder_info"],
            "original_decoder_profile":galois_cert["original_decoder_profile_timings"],
            "stageB_exact_polytope_bounds": galois_cert["stage_B_polytope"][
                "certified_coordinate_bounds"
            ],
            "stageB_integer_nullity":gsys["integer_nullity"],
            "stageB_complete":galois_cert["comparison"]["complete"],
            "qualifying_count":len(galois_list),
            "qualifying_list_matches_reference":galois_list == ref_set,
            "all_timings_seconds":galois_cert["timings_seconds"],
        },
        "affine_space_comparison":affine_comparison,
        "reference": {
            "list_size":len(ref_set), "stats":reference_info,
            "used_for_comparison_only":True
        },
        "timings_seconds":times,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,default=Path("scalar-affine-trial0.json"))
    args=p.parse_args()
    result=run()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    extraction=result["scalar"]["exact_integer_extraction"]
    print(json.dumps({
        "received_SHA256": result["scope"]["received_SHA256"],
        "dimension": result["scalar"]["lattice_dimension"],
        "N_bits": result["scalar"]["N_bits"],
        "L2_basis_U": [v["U_exact"] for v in result["scalar"][
            "all_eleven_reduced_basis_relations"]],
        "P_min6": result["scalar"]["P_min6"],
        "certified_relation_count":result["scalar"]["retained_relation_count"],
        "scalar_rank":result["scalar"]["rank_Q"],
        "scalar_nullity":result["scalar"]["nullity_Q"],
        "integer_kernel_saturated":extraction.get("certified_saturation"),
        "certified_coordinate_bounds": extraction.get(
            "certified_integer_coordinate_bounds"),
        "outer_box_size":extraction.get("outer_box_size"),
        "integer_enumeration": extraction.get("enumeration"),
        "complete":extraction["finished"],
        "incomplete_reason":extraction.get("reason"),
        "qualifying_count":result["scalar"]["qualifying_count"],
        "reference_count":result["reference"]["list_size"],
        "affine_space_comparison":result["affine_space_comparison"],
        "timings_seconds":result["timings_seconds"],
        "scalar_extraction_timings":extraction.get("timings_seconds"),
        "galois_decode_timings":result["galois"]["original_decoder_profile"],
        "certificate_artifact":str(args.output)
    },sort_keys=True,indent=2))


if __name__=="__main__":
    main()
