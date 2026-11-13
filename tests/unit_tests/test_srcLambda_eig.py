import unittest
import numpy as np

# Closed-form spectral bound used by computeSrcLambda (src/turbulence/saGammaRetheta.F90) for the
# source-term time-step restriction of the SA-noft2-Gamma-Retheta model. The 3x3 source Jacobian A of
# (nuTilde, gamma, reThetat) has A13 = A31 = A32 = 0, so its eigenvalues are A33 and the eigenvalues
# of the 2x2 block [A11, A12; A21, A22], and no cubic has to be solved.
#
# The routine is not exposed to python: these tests check the python transcription of the formula
# below against LAPACK, i.e. the algebra, not the compiled code.


def srcLambdaFull(A):
    """Largest non-negative real part of the eigenvalues of A (mode "full" of computeSrcLambda)."""
    tr2 = A[0, 0] + A[1, 1]
    disc2 = ((A[0, 0] - A[1, 1]) / 2) ** 2 + A[0, 1] * A[1, 0]
    if disc2 >= 0:
        lam2x2 = tr2 / 2 + np.sqrt(disc2)
    else:
        # complex conjugate pair: the real part is tr2 / 2
        lam2x2 = tr2 / 2
    return max(0.0, lam2x2, A[2, 2])


def srcLambdaTransition(A):
    """Mode "transition": SA on its own; the gamma-reThetat block is triangular because A32 = 0."""
    return max(0.0, A[0, 0]), max(0.0, A[1, 1], A[2, 2])


def srcLambdaReference(A):
    return max(0.0, np.max(np.linalg.eigvals(A).real))


class TestSrcLambdaEig(unittest.TestCase):
    N_PROCS = 1

    cases = {
        "typical": [[-2.5, 0.1, 0.0], [0.5, 1.2, -0.3], [0.0, 0.0, -0.8]],
        "complex 2x2 pair": [[1.0, 2.0, 0.0], [-2.0, 1.0, 0.5], [0.0, 0.0, 0.5]],
        "A33 dominates": [[-1.0, 0.1, 0.0], [0.1, -0.5, 0.2], [0.0, 0.0, 5.0]],
        "all negative": [[-5.0, 0.1, 0.0], [0.0, -2.0, 0.1], [0.0, 0.0, -3.0]],
        "diagonal": [[3.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 1.0]],
        "strong SA-gamma coupling": [[0.5, 3.0, 0.0], [2.0, 0.5, -0.1], [0.0, 0.0, -1.0]],
        "all positive": [[2.0, 0.5, 0.0], [0.3, 1.5, 0.2], [0.0, 0.0, 0.8]],
    }

    def test_full_mode(self):
        for name, A in self.cases.items():
            A = np.array(A)
            np.testing.assert_allclose(srcLambdaFull(A), srcLambdaReference(A), rtol=1e-12, atol=1e-12, err_msg=name)

    def test_transition_mode(self):
        # the gamma-reThetat block [A22, A23; 0, A33] on its own
        for name, A in self.cases.items():
            A = np.array(A)
            lamSA, lamTrans = srcLambdaTransition(A)
            self.assertEqual(lamSA, max(0.0, A[0, 0]), name)
            np.testing.assert_allclose(lamTrans, srcLambdaReference(A[1:, 1:]), rtol=1e-12, atol=1e-12, err_msg=name)


if __name__ == "__main__":
    unittest.main()
