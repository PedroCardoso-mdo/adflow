# built-ins
import unittest
import os
import copy
from parameterized import parameterized_class

# MACH classes
from adflow import ADFLOW, ADFLOW_C

from reg_default_options import adflowDefOpts, defaultAeroDVs
import reg_test_utils as utils
from baseclasses.testing import getTol

import reg_test_classes
from reg_transition import sabcm_test_params, setSABCMTestState


baseDir = os.path.dirname(os.path.abspath(__file__))

# four cases: smooth and original SA-BCM, each at a turbulent and at a transitional state
test_params = copy.deepcopy(sabcm_test_params)


@parameterized_class(test_params)
class TestJacVecFwdSABCM(reg_test_classes.RegTest):
    """
    Tests that given a flow state the residuals and the FWD jacobian vector products of the SA-BCM model
    agree with the values recorded in the ref file.
    """

    N_PROCS = 2

    def setUp(self):
        if not hasattr(self, "name"):
            # return immediately when the setup method is being called on the based class and NOT the
            # classes created using parametrized
            # this will happen when testing, but will hopefully be fixed down the line
            return

        super().setUp()

        options = copy.copy(adflowDefOpts)
        options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
        options.update(self.options)

        # Create the solver
        self.CFDSolver = ADFLOW(options=copy.deepcopy(options), debug=True)

        self.ap = copy.deepcopy(self.aero_prob)
        # add the default dvs to the problem
        for dv in defaultAeroDVs:
            self.ap.addDV(dv)

        # propagates the values from the restart file throughout the code and builds the test state
        setSABCMTestState(self.CFDSolver, self.ap, self.nuTildeScale)

    def test_residuals(self):
        utils.assert_residuals_allclose(self.handler, self.CFDSolver, self.ap, tol=1e-10)

        # The model only changes the turbulence equation, so check its residual on its own as well
        rtol, atol = getTol(tol=1e-10)
        res = self.CFDSolver.getResidual(self.ap) / self.CFDSolver.getFreeStreamResidual(self.ap)
        nw = self.CFDSolver.adflow.flowvarrefstate.nw
        self.handler.root_print("Norm of turbulence residual")
        self.handler.par_add_norm("Norm of turbulence residual", res.reshape(-1, nw)[:, -1], rtol=rtol, atol=atol)

    # ------------------- Derivative routine checks ----------------------------
    def test_wDot(self):
        utils.assert_fwd_mode_wdot_allclose(self.handler, self.CFDSolver, self.ap, seed=314, tol=5e-9)

    def test_xVDot(self):
        utils.assert_fwd_mode_xVDot_allclose(self.handler, self.CFDSolver, self.ap, seed=314, tol=1e-10)


@parameterized_class(test_params)
class TestJacVecFwdSABCMCS(reg_test_classes.CmplxRegTest):
    """
    Tests the FWD jacobian vector products of the SA-BCM model against complex step: the CS values are
    compared with the AD values the real build trained into the ref file.
    """

    N_PROCS = 2

    h = 1e-40

    def setUp(self):
        if not hasattr(self, "name"):
            # return immediately when the setup method is being called on the based class and NOT the
            # classes created using parametrized
            # this will happen when training, but will hopefully be fixed down the line
            return

        super().setUp()

        options = copy.copy(adflowDefOpts)
        options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
        options.update(self.options)

        self.CFDSolver = ADFLOW_C(options=copy.deepcopy(options), debug=True)

        self.ap = copy.deepcopy(self.aero_prob)
        # add the default dvs to the problem
        for dv in defaultAeroDVs:
            self.ap.addDV(dv)

        # propagates the values from the restart file throughout the code and builds the test state
        setSABCMTestState(self.CFDSolver, self.ap, self.nuTildeScale)

    # ------------------- Derivative routine checks ----------------------------
    def cmplx_test_wDot(self):
        if not hasattr(self, "name"):
            return

        wDot = self.CFDSolver.getStatePerturbation(314)

        resDot, funcsDot, fDot = self.CFDSolver.computeJacobianVectorProductFwd(
            wDot=wDot, residualDeriv=True, funcDeriv=True, fDeriv=True, mode="CS", h=self.h
        )

        rtol, atol = getTol(tol=5e-9)

        self.handler.root_print("||dR/dw * wDot||")
        self.handler.par_add_norm("||dR/dw * wDot||", resDot, rtol=rtol, atol=atol)

        self.handler.root_print("dFuncs/dw * wDot")
        self.handler.root_add_dict("dFuncs/dw * wDot", funcsDot, rtol=rtol, atol=atol)

        self.handler.root_print("||dF/dw * wDot||")
        self.handler.par_add_norm("||dF/dw * wDot||", fDot, rtol=rtol, atol=atol)

    def cmplx_test_xVDot(self):
        if not hasattr(self, "name"):
            return

        xVDot = self.CFDSolver.getSpatialPerturbation(314)

        resDot, funcsDot, fDot = self.CFDSolver.computeJacobianVectorProductFwd(
            xVDot=xVDot, residualDeriv=True, funcDeriv=True, fDeriv=True, mode="CS", h=self.h
        )

        rtol, atol = getTol(tol=1e-10)

        self.handler.root_print("||dR/dXv * xVDot||")
        self.handler.par_add_norm("||dR/dXv * xVDot||", resDot, rtol=rtol, atol=atol)

        # These can be finiky sometimes so a bigger tolerance.
        self.handler.root_print("dFuncs/dXv * xVDot")
        self.handler.root_add_dict("dFuncs/dXv * xVDot", funcsDot, rtol=rtol * 10, atol=atol * 10)

        self.handler.root_print("||dF/dXv * xVDot||")
        self.handler.par_add_norm("||dF/dXv * xVDot||", fDot, rtol=rtol, atol=atol)


if __name__ == "__main__":
    unittest.main()
