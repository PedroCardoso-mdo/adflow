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
import reg_sagr
from reg_sagr import ap_sagr_flatplate, sagrBaseOptions, sagrAeroDVs


baseDir = os.path.dirname(os.path.abspath(__file__))

test_params = [
    {
        "name": "sagr_flatplate",
        "options": copy.deepcopy(sagrBaseOptions),
        "ref_file": "jacvecbwd_sagr_flatplate.json",
        "aero_prob": copy.deepcopy(ap_sagr_flatplate),
        "N_PROCS": 1,
    },
]


@parameterized_class(test_params)
class TestJacVecBWDFastSAGR(reg_test_classes.RegTest):
    """
    Tests that the reverse-fast state-only mode (Tapenade _fast_b routines,
    computeJacobianVectorProductBwdFast) is consistent with the full reverse
    mode (_b) for the SA-noft2-Gamma-Retheta model. Mirrors TestJacVecBWDFast
    (test_jacVecProdBWDFast.py), extended with row-block-seeded checks that
    isolate the transition equations.

    Note (docs/audits/sst_dev_lessons.md, watch item 1): the
    autoEditReverseFast.py push/pop stripping that broke SST's _fast_b is
    still active on this branch. A failure here that is isolated to the
    gamma/reThetat row seeds is the signature of that stripping bug.
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
        # add the SA-GR aero dvs to the problem
        for dv in sagrAeroDVs:
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
        # seed one residual row block at a time (meanflow / nuTilde / gamma /
        # reThetat) so each transposed coupling column is checked separately
        reg_sagr.assert_bwdfast_blocks_allclose(self.CFDSolver, seed=314, atol=1e-16)

    def test_repeated_calls(self):
        dwBar = self.CFDSolver.getStatePerturbation(314)

        wBarfast1 = self.CFDSolver.computeJacobianVectorProductBwdFast(resBar=dwBar)
        wBarfast2 = self.CFDSolver.computeJacobianVectorProductBwdFast(resBar=dwBar)

        np.testing.assert_allclose(wBarfast1, wBarfast2, atol=1e-16, err_msg="w wrt res double call")


@parameterized_class(test_params)
class TestDotProductsSAGR(reg_test_classes.RegTest):
    """
    Transpose (dot-product) tests between the forward (_d) and reverse (_b)
    modes for the SA-noft2-Gamma-Retheta model: the full set from
    reg_test_utils.assert_dot_products_allclose plus blockwise products that
    isolate the SA<->transition coupling entries from audit 06.

    Transpose tests prove forward/reverse *consistency*, not correctness —
    correctness comes from the CS classes in test_jacVecProdFWD_sagr.py and
    test_adjoint_sagr.py (sst_dev post-mortem).
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
        # add the SA-GR aero dvs to the problem
        for dv in sagrAeroDVs:
            self.ap.addDV(dv)

        # propagates the values from the restart file throughout the code
        self.CFDSolver.getResidual(self.ap)

    def test_dot_products(self):
        utils.assert_dot_products_allclose(self.handler, self.CFDSolver, tol=2e-10)

    def test_coupling_dot_products(self):
        reg_sagr.assert_coupling_dot_products_allclose(self.handler, self.CFDSolver, seed=314, tol=2e-10)


if __name__ == "__main__":
    unittest.main()
