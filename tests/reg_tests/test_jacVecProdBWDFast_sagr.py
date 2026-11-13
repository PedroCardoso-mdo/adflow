# built-ins
import unittest
import os
import copy
from parameterized import parameterized_class
import numpy as np

# MACH classes
from adflow import ADFLOW

from reg_default_options import adflowDefOpts
import reg_test_utils as utils

import reg_test_classes
import reg_transition
from reg_transition import ap_sagr_tut_wing, sagrOptions, sagrCrossflowOptions, transitionAeroDVs


baseDir = os.path.dirname(os.path.abspath(__file__))

test_params = [
    {
        "name": "sagr_tut_wing",
        "options": copy.deepcopy(sagrOptions),
        "ref_file": "jacvecbwd_sagr_tut_wing.json",
        "aero_prob": copy.deepcopy(ap_sagr_tut_wing),
        "dot_product_pairs": reg_transition.couplingDotProductPairs["sagr"],
        "N_PROCS": 1,
    },
    # same state, crossflow transition term switched on
    {
        "name": "sagr_crossflow_tut_wing",
        "options": copy.deepcopy(sagrCrossflowOptions),
        "ref_file": "jacvecbwd_sagr_crossflow_tut_wing.json",
        "aero_prob": copy.deepcopy(ap_sagr_tut_wing),
        "dot_product_pairs": reg_transition.couplingDotProductPairs["sagr"],
        "N_PROCS": 1,
    },
]


@parameterized_class(test_params)
class TestJacVecBWDFastSAGR(reg_test_classes.RegTest):
    """
    Tests that given a flow state the reverse-fast (state only) jacobian vector products of the
    SA-noft2-Gamma-Retheta model agree with the full reverse mode. Same checks as TestJacVecBWDFast
    (test_jacVecProdBWDFast.py) plus the products seeded one residual row block at a time.
    """

    # nothing is recorded in the ref file, so there is nothing to train
    no_train = True

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
        # add the aero dvs to the problem
        for dv in transitionAeroDVs:
            self.ap.addDV(dv)

        # propagates the values from the restart file throughout the code
        self.CFDSolver.getResidual(self.ap)

    # ------------------- Derivative routine checks ----------------------------
    def test_BWD(self):
        dwBar = self.CFDSolver.getStatePerturbation(314)

        wBar = self.CFDSolver.computeJacobianVectorProductBwd(
            resBar=dwBar,
            wDeriv=True,
        )

        wBarfast = self.CFDSolver.computeJacobianVectorProductBwdFast(resBar=dwBar)

        np.testing.assert_allclose(wBar, wBarfast, atol=1e-16, err_msg="w wrt res")

    def test_BWD_block_rows(self):
        # seed one residual row block at a time so each transposed coupling column is checked separately
        reg_transition.assert_bwdfast_blocks_allclose(self.CFDSolver, seed=314, atol=1e-16)

    def test_repeated_calls(self):
        dwBar = self.CFDSolver.getStatePerturbation(314)

        wBarfast1 = self.CFDSolver.computeJacobianVectorProductBwdFast(resBar=dwBar)
        wBarfast2 = self.CFDSolver.computeJacobianVectorProductBwdFast(resBar=dwBar)

        np.testing.assert_allclose(wBarfast1, wBarfast2, atol=1e-16, err_msg="w wrt res double call")


@parameterized_class(test_params)
class TestDotProductsSAGR(reg_test_classes.RegTest):
    """
    Dot-product tests between the forward and reverse modes of the SA-noft2-Gamma-Retheta model: the
    set of reg_test_utils.assert_dot_products_allclose plus blockwise products that isolate the coupling
    between the SA and the transition equations.

    Dot-product tests show that the forward and reverse modes are consistent, not that they are correct;
    the complex-step classes in test_jacVecProdFWD_sagr.py and test_adjoint_sagr.py check correctness.
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
        # add the aero dvs to the problem
        for dv in transitionAeroDVs:
            self.ap.addDV(dv)

        # propagates the values from the restart file throughout the code
        self.CFDSolver.getResidual(self.ap)

    def test_dot_products(self):
        utils.assert_dot_products_allclose(self.handler, self.CFDSolver, tol=2e-10)

    def test_coupling_dot_products(self):
        reg_transition.assert_coupling_dot_products_allclose(
            self.handler, self.CFDSolver, self.dot_product_pairs, seed=314, tol=2e-10
        )


if __name__ == "__main__":
    unittest.main()
