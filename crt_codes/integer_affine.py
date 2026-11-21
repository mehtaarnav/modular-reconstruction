"""Exact integer affine-lattice certification (column-Hermite Euclidean operations).

No floating point, candidate cutoffs, or unsaturated rational nullspaces.
Given integer T (m x n), b (m), produce a unimodular V (n x n)
with T*V = [C | 0], where C has full column rank r. Solving C*y=b
over Z yields x0=V*(y,0); the final n-r columns of V form a saturated
Z-basis for ker_Z(T).

Rank-2 kernel basis reduction uses exact Gauss operations, also
unimodular. Polytope bounds come from ALL rational line intersections.
"""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from math import prod


def extended_gcd(a: int, b: int) -> tuple[int, int, int]:
    """Return g>0 and u,v with u*a+v*b=g."""
    if a == b == 0:
        raise ValueError("gcd(0,0) is not a pivot")
    ar, br, u0, u1, v0, v1 = a, b, 1, 0, 0, 1
    while br:
        q = ar // br
        ar, br = br, ar - q * br
        u0, u1 = u1, u0 - q * u1
        v0, v1 = v1, v0 - q * v1
    if ar < 0:
        return -ar, -u0, -v0
    return ar, u0, v0


def exact_dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def matmul(A, B):
    return [[exact_dot(row, column) for column in zip(*B)] for row in A]


def column_hermite_solution(T: list[list[int]], rhs: list[int]):
    """Column Hermite-style elimination, with explicit unimodular V.

    The form T*V=[C|0] need not be the *canonical* HNF, but every column
    operation is a GL_n(Z) operation, so the saturation certificate is
    independent of canonical pivot choices.
    """
    if not T or not T[0] or len(T) != len(rhs):
        raise ValueError("Invalid or mismatched integer linear system")
    n = len(T[0])
    if any(len(row) != n or any(type(v) is not int for v in row) for row in T):
        raise ValueError("T must be rectangular with integer entries")
    if any(type(v) is not int for v in rhs):
        raise ValueError("rhs must be integral")
    A = [row[:] for row in T]
    V = [[int(i == j) for j in range(n)] for i in range(n)]
    rank = 0
    pivots = []
    elementary_ops = 0

    def col_op(i, j, s, t, u, v):
        nonlocal elementary_ops
        if s * v - t * u not in (-1, 1):
            raise ArithmeticError("Column operation is not unimodular")
        for matrix in (A, V):
            for row in matrix:
                ai, aj = row[i], row[j]
                row[i] = s * ai + t * aj
                row[j] = u * ai + v * aj
        elementary_ops += 1

    for row_id in range(len(T)):
        if rank == n:
            break
        pivot_col = next((j for j in range(rank, n) if A[row_id][j]), None)
        if pivot_col is None:
            continue
        if pivot_col != rank:
            col_op(rank, pivot_col, 0, 1, 1, 0)
        for j in range(rank + 1, n):
            a, b = A[row_id][rank], A[row_id][j]
            if not b:
                continue
            g, s, t = extended_gcd(a, b)
            col_op(rank, j, s, t, -b // g, a // g)
            assert A[row_id][rank] == g and A[row_id][j] == 0
        pivots.append(row_id)
        rank += 1

    if matmul(T, V) != A:
        raise ArithmeticError("Unimodular reduction certificate failed")
    if any(A[i][j] for i in range(len(T)) for j in range(rank, n)):
        raise ArithmeticError("The alleged kernel columns are not zero")
    y = [0] * n
    for col, row_id in enumerate(pivots):
        numerator = rhs[row_id] - sum(A[row_id][j] * y[j] for j in range(col))
        d = A[row_id][col]
        if not d or numerator % d:
            raise ValueError("The integer system has no solution")
        y[col] = numerator // d

    x0 = [sum(V[i][j] * y[j] for j in range(rank)) for i in range(n)]
    if [exact_dot(row, x0) for row in T] != rhs:
        raise ValueError("No integer solution (inconsistent remaining rows)")
    basis = [[V[i][j] for i in range(n)] for j in range(rank, n)]
    if any(any(exact_dot(row, b) for row in T) for b in basis):
        raise ArithmeticError("Kernel identity failed")
    return x0, basis, V, {
        "rank": rank, "nullity": n - rank, "pivot_rows": pivots,
        "elementary_unimodular_column_operations": elementary_ops,
    }


def gauss_reduce_two(basis):
    """Exact Gauss reduction in rank 2, returning basis and GL_2(Z) change."""
    if len(basis) != 2:
        raise ValueError("This exact reduction handles rank two")
    B = [list(v) for v in basis]
    U = [[1, 0], [0, 1]]
    iterations = 0
    while True:
        iterations += 1
        if iterations > 100000:
            raise ArithmeticError("Exact Gauss reduction did not terminate")
        if exact_dot(B[1], B[1]) < exact_dot(B[0], B[0]):
            B.reverse()
            U.reverse()
        q = round(Fraction(exact_dot(B[0], B[1]), exact_dot(B[0], B[0])))
        if q == 0:
            break
        B[1] = [u - q * v for u, v in zip(B[1], B[0])]
        U[1] = [u - q * v for u, v in zip(U[1], U[0])]
    if abs(U[0][0] * U[1][1] - U[0][1] * U[1][0]) != 1:
        raise ArithmeticError("Kernel reduction lost unimodularity")
    if B != matmul(U, basis):
        raise ArithmeticError("Kernel reduction change of basis mismatch")
    return B, U, iterations


def center_particular(x0, basis):
    """Translate integer x0 by kernel vectors, using exact rational projection."""
    if len(basis) != 2:
        raise ValueError("This centering handles rank two")
    a = exact_dot(basis[0], basis[0])
    b = exact_dot(basis[0], basis[1])
    c = exact_dot(basis[1], basis[1])
    det = a * c - b * b
    if det <= 0:
        raise ArithmeticError("Kernel basis not independent")
    h0, h1 = exact_dot(basis[0], x0), exact_dot(basis[1], x0)
    shift = (
        round(Fraction(-c * h0 + b * h1, det)),
        round(Fraction(b * h0 - a * h1, det)),
    )
    return [
        x0[i] + sum(shift[j] * basis[j][i] for j in range(2))
        for i in range(len(x0))
    ], shift


def ceil_exact(q: Fraction) -> int:
    return -((-q.numerator) // q.denominator)


def floor_exact(q: Fraction) -> int:
    return q.numerator // q.denominator


def build_box_inequalities(x0, basis, H: int):
    """(a,b,c) means a*t1 + b*t2 <= c, all entries integers."""
    if H < 1 or len(basis) != 2:
        raise ValueError("Expected positive H and rational nullity two")
    inequalities = []
    for i in range(len(x0)):
        a, b = basis[0][i], basis[1][i]
        inequalities.extend((
            (a, b, H - 1 - x0[i]),
            (-a, -b, H - 1 + x0[i]),
        ))
    return inequalities


def exact_rank_two_lp(inequalities):
    """Exact continuous 2D LP by rational feasible vertex enumeration."""
    vertices = set()
    for (a,b,c),(d,e,f) in combinations(inequalities, 2):
        det = a * e - b * d
        if not det:
            continue
        u = Fraction(c * e - b * f, det)
        v = Fraction(a * f - c * d, det)
        if all(p * u + q * v <= rhs for p,q,rhs in inequalities):
            vertices.add((u,v))
    if not vertices:
        raise ValueError("No feasible vertices: empty or lower-dimensional polytope")
    mins = [min(z[j] for z in vertices) for j in range(2)]
    maxs = [max(z[j] for z in vertices) for j in range(2)]
    bounds = [(ceil_exact(mins[j]), floor_exact(maxs[j])) for j in range(2)]
    return {
        "bounds": bounds, "continuous_min": mins,
        "continuous_max": maxs, "vertices": sorted(vertices),
        "outer_box_cardinality": prod(max(0, hi - lo + 1) for lo,hi in bounds),
    }


def complete_integer_points(x0, basis, H, inequalities, bounds,
                            max_box_cardinality=2_000_000):
    """Exact rank-two polytope enumeration; no silent truncation.

    If certified outer box is too large, return a clear refusal with the
    exact box size. Never represent a truncated list as exhaustive.
    """
    (lo1, hi1), (lo2, hi2) = bounds
    box_size = max(0, hi1 - lo1 + 1) * max(0, hi2 - lo2 + 1)
    if box_size > max_box_cardinality:
        return None, {
            "complete": False, "reason": "outer box exceeds declared enumeration budget",
            "outer_box_cardinality": box_size,
            "budget": max_box_cardinality,
            "feasible_integer_points_checked": 0,
        }
    found = []
    axis_one_rows = 0
    pair_checks = 0
    for t1 in range(lo1, hi1 + 1):
        axis_one_rows += 1
        low, high = lo2, hi2
        infeasible = False
        for a,b,c in inequalities:
            rhs = c - a * t1
            if b == 0:
                if rhs < 0:
                    infeasible = True
                    break
            elif b > 0:
                high = min(high, floor_exact(Fraction(rhs, b)))
            else:
                low = max(low, ceil_exact(Fraction(rhs, b)))
            if low > high:
                infeasible = True
                break
        if infeasible:
            continue
        for t2 in range(low, high + 1):
            pair_checks += 1
            x = [
                x0[i] + t1 * basis[0][i] + t2 * basis[1][i]
                for i in range(len(x0))
            ]
            if not all(-H + 1 <= value <= H - 1 for value in x):
                raise ArithmeticError("Inequality-sweep error")
            found.append((t1, t2, x))
    return found, {
        "complete": True, "outer_box_cardinality": box_size,
        "axis_one_rows_visited": axis_one_rows,
        "feasible_integer_points_checked": pair_checks,
        "all_polytope_integer_points": len(found),
        "budget": max_box_cardinality,
    }
