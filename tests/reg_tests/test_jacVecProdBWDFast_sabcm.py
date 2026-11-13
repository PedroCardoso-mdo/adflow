# built-ins
import unittest
import os
import copy
from parameterized import parameterized_class
import numpy as np
from mpi4py import MPI

# MACH classes
from adflow import ADFLOW

from reg_default_options import adflowDefOpts

from reg_transition import sabcm_test_params, setSABCMTestState


baseDir = os.path.dirname(os.path.abspath(__file__))

# four cases: smooth and original SA-BCM, each at a turbulent and at a transitional state
test_params = copy.deepcopy(sabcm_test_params)


@parameterized_class(test_params)
class TestJacVecBWDFastSABCM(unittest.TestCase):
    """
    Tests that given a flow state the reverse and reverse-fast jacobian vector products of the SA-BCM
    model are consistent with each other and with the forward mode (dot-product test).

    No reference file is needed, so this inherits unittest.TestCase and not RegTest.
    """

    N_PROCS = 2

    def setUp(self):
        if not hasattr(self, "name"):
            # return immediately when the setup method is being called on the based class and NOT the
            # classes created using parametrized
            # this will happen when testing, but will hopefully be fixed down the line
            return

        options = copy.copy(adflowDefOpts)
        options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
        options.update(self.options)

        # Create the solver
        self.CFDSolver = ADFLOW(options=copy.deepcopy(options), debug=True)

        self.ap = copy.deepcopy(self.aero_prob)

        # propagates the values from the restart file throughout the code and builds the test state
        setSABCMTestState(self.CFDSolver, self.ap, self.nuTildeScale)

    # ------------------- Derivative routine checks ----------------------------
    def test_BWD(self):
        dwBar = self.CFDSolver.getStatePerturbation(314)

        wBar = self.CFDSolver.computeJacobianVectorProductBwd(resBar=dwBar, wDeriv=True)
        wBarfast = self.CFDSolver.computeJacobianVectorProductBwdFast(resBar=dwBar)

        np.testing.assert_allclose(wBar, wBarfast, atol=1e-16, err_msg="w wrt res")

    def test_dot_product(self):
        wDot = self.CFDSolver.getStatePerturbation(314)
        resBar = self.CFDSolver.getStatePerturbation(271)

        resDot = self.CFDSolver.computeJacobianVectorProductFwd(wDot=wDot, residualDeriv=True)
        wBar = self.CFDSolver.computeJacobianVectorProductBwd(resBar=resBar, wDeriv=True)
        wBarFast = self.CFDSolver.computeJacobianVectorProductBwdFast(resBar=resBar)

        fwdProd = self.CFDSolver.comm.allreduce(np.dot(resBar, resDot), op=MPI.SUM)
        bwdProd = self.CFDSolver.comm.allreduce(np.dot(wBar, wDot), op=MPI.SUM)
        bwdFastProd = self.CFDSolver.comm.allreduce(np.dot(wBarFast, wDot), op=MPI.SUM)

        np.testing.assert_allclose(bwdProd, fwdProd, rtol=1e-10, err_msg="BWD dot product")
        np.testing.assert_allclose(bwdFastProd, fwdProd, rtol=1e-10, err_msg="BWDFast dot product")


if __name__ == "__main__":
    unittest.main()
