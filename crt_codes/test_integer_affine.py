"""Exact arithmetic unit tests. No additional received-word experiments."""
from itertools import product

from sympy import Matrix

from crt_codes.integer_affine import (
    column_hermite_solution, gauss_reduce_two, center_particular,
    build_box_inequalities, exact_rank_two_lp, complete_integer_points,
    exact_dot, matmul, extended_gcd,
)


def test_signed_extended_gcd_and_unimodular_diagonalization():
    for a, b in ((45, -66), (-45, 66), (-45, -66), (0, -13), (17, 0)):
        g, u, v = extended_gcd(a,b)
        assert g > 0
        assert a*u+b*v == g


def test_column_hermite_certifies_saturated_kernel_and_integer_solution():
    T = [[2, 4, 0, 0], [0, 6, 9, 0], [4, 8, 0, 0]]
    hidden_solution = [-1, 2, 0, 1]
    rhs = [exact_dot(row, hidden_solution) for row in T]
    x0, basis, V, data = column_hermite_solution(T, rhs)
    assert data["rank"] == 2 and data["nullity"] == 2
    assert abs(Matrix(V).det()) == 1
    assert [exact_dot(row, x0) for row in T] == rhs
    assert all([exact_dot(row, vec) for row in T] == [0]*len(T) for vec in basis)
    assert all(all(row[j] == 0 for j in range(2,4)) for row in matmul(T,V))
    inverse_V = Matrix(V).inv()
    assert all(value.q == 1 for value in inverse_V)
    for candidate in product(range(-2,3), repeat=4):
        if [exact_dot(row, candidate) for row in T] != rhs:
            continue
        coords = inverse_V * (Matrix(candidate)-Matrix(x0))
        assert all(v == 0 for v in coords[:2,0])
        assert all(v.q == 1 for v in coords[2:,0])
        assert list(candidate) == [
            x0[i]+sum(int(coords[2+j])*basis[j][i] for j in range(2))
            for i in range(4)
        ]


def test_exact_gauss_unimodular_change_and_centering():
    old = [[1000, 0, 0], [1001, 1, 0]]
    reduced, U, _ = gauss_reduce_two(old)
    assert abs(Matrix(U).det()) == 1
    assert matmul(U, old) == reduced
    x0 = [1001*2000+1000*3000, 2000, 0]
    centered, shift = center_particular(x0, reduced)
    assert all(isinstance(v, int) for v in centered)
    assert centered == [x0[i]+sum(shift[j]*reduced[j][i] for j in range(2)) for i in range(3)]


def test_exact_rational_lp_outer_box_and_complete_integer_count():
    x0, basis, H = [0, 0, 0], [[1,1,0], [1,-1,0]], 3
    inequalities = build_box_inequalities(x0, basis, H)
    poly = exact_rank_two_lp(inequalities)
    assert poly["bounds"] == [(-2,2),(-2,2)]
    assert poly["outer_box_cardinality"] == 25
    found, stats = complete_integer_points(x0,basis,H,inequalities,poly["bounds"])
    assert stats["complete"] is True
    assert stats["feasible_integer_points_checked"] == 13
    assert len(found) == 13
    assert {tuple((u,v)) for u,v,_ in found} == {
        (u,v) for u,v in product(range(-2,3), repeat=2)
        if abs(u+v)<=2 and abs(u-v)<=2
    }


def test_declared_budget_does_not_pretend_to_be_complete():
    x0, basis, H = [0,0,0], [[1,0,0],[0,1,0]], 10000
    inequalities = build_box_inequalities(x0,basis,H)
    poly = exact_rank_two_lp(inequalities)
    found, stats = complete_integer_points(x0,basis,H,inequalities,poly["bounds"],
                                          max_box_cardinality=1)
    assert found is None
    assert stats["complete"] is False
    assert stats["outer_box_cardinality"] == (2*H-1)**2
