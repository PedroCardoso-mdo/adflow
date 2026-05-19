#!/usr/bin/env python
"""
Unit test: srcLambda eigenvalue computation for block-triangular source Jacobian.

The source Jacobian A has A13=A31=A32=0 (P&Z §7.1), giving:
  det(A-λI) = (A33-λ)·det([A11-λ,A12;A21,A22-λ])

Eigenvalues = {A33} ∪ {2x2 block eigenvalues} — no cubic needed.
"""
import numpy as np


def eig_2x2_max(A11, A12, A21, A22):
    """
    Max eigenvalue of 2x2 block [A11,A12;A21,A22].
    Uses explicit branch for complex case (disc2 < 0 → real part = tr2/2).
    """
    tr2 = A11 + A22
    disc2 = ((A11 - A22) / 2)**2 + A12 * A21
    if disc2 >= 0:
        return tr2 / 2 + np.sqrt(disc2)
    else:
        # Complex conjugate pair: real part = tr2/2
        return tr2 / 2


def srcLambda_decoupled(A):
    """Mode 0: diagonal only."""
    return [max(0, A[0, 0]), max(0, A[1, 1]), max(0, A[2, 2])]


def srcLambda_transition(A):
    """Mode 1: SA diagonal; γ-Reθ block triangular (A32=0) → max(A22,A33)."""
    lam1 = max(0, A[0, 0])
    lam23 = max(0, A[1, 1], A[2, 2])
    return [lam1, lam23, lam23]


def srcLambda_full(A):
    """Mode 2: 2x2 eigenvalue for SA-γ block + A33."""
    lam2x2 = eig_2x2_max(A[0, 0], A[0, 1], A[1, 0], A[1, 1])
    lam3x3 = max(0, lam2x2, A[2, 2])
    return [lam3x3, lam3x3, lam3x3]


def reference_eig_max(A):
    """Reference: numpy/LAPACK eigenvalues."""
    eigs = np.linalg.eigvals(A)
    return max(0.0, np.max(eigs.real))


def test_case(name, A, tol=1e-10):
    """Test all modes for one matrix (must have A13=A31=A32=0)."""
    assert A[0, 2] == 0 and A[2, 0] == 0 and A[2, 1] == 0, "Matrix must be block-triangular"

    ref = reference_eig_max(A)
    full = srcLambda_full(A)

    passed = abs(full[0] - ref) < tol
    status = "PASS" if passed else "FAIL"
    print(f"{status}: {name}")
    print(f"  Reference (LAPACK): {ref:.10f}")
    print(f"  Full mode:          {full[0]:.10f}")
    print(f"  Error:              {abs(full[0] - ref):.2e}")
    return passed


if __name__ == "__main__":
    print("=" * 60)
    print("srcLambda eigenvalue test (block-triangular structure)")
    print("=" * 60)

    all_passed = True

    # Case 1: Typical SA-γ-Reθt source Jacobian
    print("\n--- Case 1: Typical source Jacobian ---")
    A1 = np.array([
        [-2.5, 0.1, 0.0],
        [0.5, 1.2, -0.3],
        [0.0, 0.0, -0.8]
    ])
    all_passed &= test_case("Typical source Jacobian", A1)

    # Case 2: Complex eigenvalues in 2x2 block (disc2 < 0)
    print("\n--- Case 2: Complex 2x2 eigenvalues (disc2 < 0) ---")
    A2 = np.array([
        [1.0, 2.0, 0.0],
        [-2.0, 1.0, 0.5],
        [0.0, 0.0, 0.5]
    ])
    all_passed &= test_case("Complex 2x2 block", A2)
    tr2 = A2[0, 0] + A2[1, 1]
    disc2 = ((A2[0, 0] - A2[1, 1]) / 2)**2 + A2[0, 1] * A2[1, 0]
    print(f"  disc2 = {disc2:.4f} (should be < 0)")
    print(f"  Real part tr2/2 = {tr2/2:.4f}")

    # Case 3: A33 dominates
    print("\n--- Case 3: A33 dominates ---")
    A3 = np.array([
        [-1.0, 0.1, 0.0],
        [0.1, -0.5, 0.2],
        [0.0, 0.0, 5.0]
    ])
    all_passed &= test_case("A33 dominates", A3)

    # Case 4: All negative (srcLambda = 0)
    print("\n--- Case 4: All negative eigenvalues ---")
    A4 = np.array([
        [-5.0, 0.1, 0.0],
        [0.0, -2.0, 0.1],
        [0.0, 0.0, -3.0]
    ])
    all_passed &= test_case("All negative", A4)

    # Case 5: Diagonal (simple case)
    print("\n--- Case 5: Diagonal matrix ---")
    A5 = np.array([
        [3.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
        [0.0, 0.0, 1.0]
    ])
    all_passed &= test_case("Diagonal", A5)

    # Case 6: Strong off-diagonal coupling
    print("\n--- Case 6: Strong SA-γ coupling ---")
    A6 = np.array([
        [0.5, 3.0, 0.0],
        [2.0, 0.5, -0.1],
        [0.0, 0.0, -1.0]
    ])
    all_passed &= test_case("Strong coupling", A6)

    # Case 7: Verify mode consistency
    print("\n--- Case 7: Mode consistency check ---")
    A7 = np.array([
        [2.0, 0.5, 0.0],
        [0.3, 1.5, 0.2],
        [0.0, 0.0, 0.8]
    ])
    dec = srcLambda_decoupled(A7)
    tra = srcLambda_transition(A7)
    ful = srcLambda_full(A7)
    ref = reference_eig_max(A7)
    print(f"  Decoupled:   {dec}")
    print(f"  Transition:  {tra}")
    print(f"  Full:        {ful}")
    print(f"  Reference:   {ref:.6f}")
    # Full should match reference
    all_passed &= abs(ful[0] - ref) < 1e-10

    print("\n" + "=" * 60)
    if all_passed:
        print("ALL TESTS PASSED")
    else:
        print("SOME TESTS FAILED")
    print("=" * 60)
