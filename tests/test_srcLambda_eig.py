#!/usr/bin/env python
"""
Unit test: compare srcLambda eigenvalue computation against numpy (LAPACK dgeev).
Focus on disc > 0 case with dominant complex conjugate pair.
"""
import numpy as np

def cubic_eig_max(A):
    """
    Replicate the Fortran computeSrcLambda logic for a 3x3 matrix A.
    Returns max(0, max real part of eigenvalues).
    """
    A11, A12, A13 = A[0, 0], A[0, 1], A[0, 2]
    A21, A22, A23 = A[1, 0], A[1, 1], A[1, 2]
    A31, A32, A33 = A[2, 0], A[2, 1], A[2, 2]

    # Characteristic polynomial coefficients
    a = A11 + A22 + A33  # trace
    b = A11*A22 + A22*A33 + A33*A11 - A12*A21 - A23*A32 - A31*A13
    c = A11*(A22*A33 - A23*A32) - A12*(A21*A33 - A23*A31) + A13*(A21*A32 - A22*A31)

    # Depressed cubic: t³ + pt + qc = 0, λ = t + a/3
    oneThird = 1.0 / 3.0
    p = b - a*a*oneThird
    qc = -2*a*a*a/27.0 + a*b*oneThird - c
    disc = qc*qc/4.0 + p*p*p/27.0

    if disc <= 0:
        # Three real roots
        sqrtP = np.sqrt(-p*oneThird)
        phi = 0.0
        if sqrtP > 1e-30:
            phi = np.arccos(np.clip(-qc/(2*sqrtP**3), -1, 1))
            t = 2 * sqrtP * np.cos(phi * oneThird)
        else:
            t = 0.0
        lambdaMax = t + a * oneThird
        t = 2 * sqrtP * np.cos((phi + 2*np.pi) * oneThird)
        lambdaMax = max(lambdaMax, t + a * oneThird)
        t = 2 * sqrtP * np.cos((phi + 4*np.pi) * oneThird)
        lambdaMax = max(lambdaMax, t + a * oneThird)
    else:
        # One real root, two complex conjugates
        r = np.sqrt(disc)
        t = np.sign(-qc/2 + r) * np.abs(-qc/2 + r)**oneThird + \
            np.sign(-qc/2 - r) * np.abs(-qc/2 - r)**oneThird
        lambdaMax = t + a * oneThird
        # Complex pair real part: (a - λ_real) / 2
        lambdaMax = max(lambdaMax, (a - lambdaMax) / 2)

    return max(0.0, lambdaMax)


def reference_eig_max(A):
    """Reference: numpy/LAPACK eigenvalues."""
    eigs = np.linalg.eigvals(A)
    return max(0.0, np.max(eigs.real))


def test_case(name, A, tol=1e-10):
    """Test one matrix."""
    our_val = cubic_eig_max(A)
    ref_val = reference_eig_max(A)
    eigs = np.linalg.eigvals(A)
    passed = np.abs(our_val - ref_val) < tol
    status = "PASS" if passed else "FAIL"
    print(f"{status}: {name}")
    print(f"  Eigenvalues: {eigs}")
    print(f"  Reference max(0, Re): {ref_val:.10f}")
    print(f"  Our result:           {our_val:.10f}")
    print(f"  Error:                {abs(our_val - ref_val):.2e}")
    return passed


if __name__ == "__main__":
    print("="*60)
    print("srcLambda eigenvalue unit test vs LAPACK dgeev")
    print("="*60)

    all_passed = True

    # Case 1: disc > 0, complex pair has LARGER real part than real root
    # This is the bug case - complex pair dominates
    print("\n--- Case 1: disc>0, complex pair dominates ---")
    A1 = np.array([
        [5.0, 1.0, 0.0],
        [0.0, 2.0, 3.0],
        [0.0, -3.0, 2.0]
    ])
    all_passed &= test_case("Complex pair dominates", A1)

    # Case 2: disc > 0, real root dominates
    print("\n--- Case 2: disc>0, real root dominates ---")
    A2 = np.array([
        [10.0, 1.0, 0.0],
        [0.0, -1.0, 2.0],
        [0.0, -2.0, -1.0]
    ])
    all_passed &= test_case("Real root dominates", A2)

    # Case 3: Three real roots (disc <= 0)
    print("\n--- Case 3: disc<=0, three real roots ---")
    A3 = np.array([
        [3.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
        [0.0, 0.0, 1.0]
    ])
    all_passed &= test_case("Diagonal (3 real)", A3)

    # Case 4: All negative eigenvalues -> srcLambda = 0
    print("\n--- Case 4: All negative eigenvalues ---")
    A4 = np.array([
        [-5.0, 1.0, 0.0],
        [0.0, -2.0, 1.0],
        [0.0, 0.0, -3.0]
    ])
    all_passed &= test_case("All negative", A4)

    # Case 5: Typical source Jacobian pattern
    print("\n--- Case 5: Typical SA-γ-Reθt source Jacobian ---")
    A5 = np.array([
        [-2.5, 0.1, 0.0],
        [0.5, 1.2, -0.3],
        [0.0, 0.0, -0.8]
    ])
    all_passed &= test_case("Typical source Jacobian", A5)

    # Case 6: Near-singular (sqrtP ~ 0) - detailed output
    print("\n--- Case 6: Near degenerate (detailed) ---")
    np.random.seed(42)
    A6 = np.array([
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0]
    ]) + 1e-12 * np.random.randn(3, 3)
    our_val = cubic_eig_max(A6)
    ref_val = reference_eig_max(A6)
    err = our_val - ref_val
    tol = 1e-5
    print(f"  Our value:     {our_val:.15f}")
    print(f"  dgeev value:   {ref_val:.15f}")
    print(f"  Error:         {err:+.2e} ({'over' if err > 0 else 'under'}estimate)")
    print(f"  Tolerance:     {tol:.0e}")
    print(f"  Status:        {'PASS' if abs(err) < tol else 'FAIL'}")
    all_passed &= abs(err) < tol

    # Case 7: R̃e_θt relaxation cell - qq33 dominates, A_source = -qq has negative diagonal
    # Typical: qq = [[small, 0, 0], [0, small, 0], [0, 0, large_positive]]
    # A_source = -qq → diagonal mostly negative → λ_max ≈ 0 or small
    print("\n--- Case 7: Reθt relaxation cell (expect srcLambda ≈ 0) ---")
    A7 = np.array([
        [-0.1, 0.01, 0.0],
        [0.0, -0.2, 0.0],
        [0.0, 0.0, -5.0]  # Strong relaxation term
    ])
    all_passed &= test_case("Reθt relaxation", A7)
    print(f"  Expected: srcLambda = 0 (all eigenvalues negative)")

    # Case 8: All-stable cell (destruction-dominated SA + stable γ + stable Reθt)
    print("\n--- Case 8: All-stable cell (expect srcLambda = 0) ---")
    A8 = np.array([
        [-2.0, 0.1, 0.0],
        [0.05, -1.5, 0.02],
        [0.0, 0.01, -3.0]
    ])
    all_passed &= test_case("All stable", A8)
    eigs8 = np.linalg.eigvals(A8)
    print(f"  Eigenvalues: {eigs8.real}")
    print(f"  All Re(λ) < 0: {all(eigs8.real < 0)}")

    print("\n" + "="*60)
    if all_passed:
        print("ALL TESTS PASSED")
    else:
        print("SOME TESTS FAILED")
    print("="*60)
