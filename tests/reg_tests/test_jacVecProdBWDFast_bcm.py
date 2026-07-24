# built-ins
import unittest
import os
import copy
from parameterized import parameterized_class
import numpy as np

# MACH classes
from adflow import ADFLOW

from reg_default_options import adflowDefOpts, defaultAeroDVs
import reg_test_utils as utils
import reg_test_classes
import reg_bcm
from reg_bcm import ap_bcm_tut_wing, bcmBaseOptionsSmooth, bcmBaseOptionsHard

baseDir = os.path.dirname(os.path.abspath(__file__))

test_params = [
    {
        "name": "bcm_smooth_tut_wing",
        "options": copy.deepcopy(bcmBaseOptionsSmooth),
        "ref_file": "jacvecbwd_bcm_smooth_tut_wing.json",
        "aero_prob": copy.deepcopy(ap_bcm_tut_wing),
        "N_PROCS": 2,
    },
    {
        "name": "bcm_hard_tut_wing",
        "options": copy.deepcopy(bcmBaseOptionsHard),
        "ref_file": "jacvecbwd_bcm_hard_tut_wing.json",
        "aero_prob": copy.deepcopy(ap_bcm_tut_wing),
        "N_PROCS": 2,
    },
]


@parameterized_class(test_params)
class TestJacVecBWDFastBCM(reg_test_classes.RegTest):
    """
    _b vs _fast_b consistency for the SA-BCM case, both SABCM_Exp variants. Mirrors
    test_jacVecProdBWDFast.py.

    NOTE: if a future guard commit (e.g. the divide-by-zero guards on sibling branch
    sa-bcm-timing, see docs/adjoint-trace.md) went through autoEditReverseFast.py, that script
    is known to strip push/pop in ways that broke another model's _fast_b upstream. Treat any
    failure here as an autoEditReverseFast.py stripping suspect *before* suspecting the SA-BCM
    model itself -- check docs/adjoint-trace.md for the current guard/rerun status first.
    """

    N_PROCS = 2

    def setUp(self):
        if not hasattr(self, "name"):
            return

        super().setUp()

        options = copy.copy(adflowDefOpts)
        options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
        options.update(self.options)

        self.CFDSolver = ADFLOW(options=copy.deepcopy(options), debug=True)

        self.ap = copy.deepcopy(self.aero_prob)
        for dv in defaultAeroDVs:
            self.ap.addDV(dv)

        self.CFDSolver.getResidual(self.ap)

    def test_BWD(self):
        dwBar = self.CFDSolver.getStatePerturbation(314)

        wBar = self.CFDSolver.computeJacobianVectorProductBwd(resBar=dwBar, wDeriv=True)
        wBarFast = self.CFDSolver.computeJacobianVectorProductBwdFast(resBar=dwBar)

        np.testing.assert_allclose(wBar, wBarFast, atol=1e-16, err_msg="w wrt res")

    def test_BWD_block_rows(self):
        # seed one residual row block at a time (meanflow / nuTilde) so each transposed
        # coupling column is checked separately
        reg_bcm.assert_bwdfast_blocks_allclose(self.CFDSolver, seed=314, atol=1e-16)

    def test_BWD_nutilde_seed(self):
        """Seed resBar on the nuTilde row only -- isolates the SA-BCM-modified production term's
        column in the transposed Jacobian. A failure isolated to this seed (vs. the full-state
        seed above passing) points at the use_SABCM-gated reverse code specifically."""
        nw, blocks = reg_bcm.getStateBlocks(self.CFDSolver)
        dwBarFull = self.CFDSolver.getStatePerturbation(314)
        dwBar = reg_bcm.maskStateVector(dwBarFull, nw, blocks["nuTilde"])

        wBar = self.CFDSolver.computeJacobianVectorProductBwd(resBar=dwBar, wDeriv=True)
        wBarFast = self.CFDSolver.computeJacobianVectorProductBwdFast(resBar=dwBar)

        np.testing.assert_allclose(wBar, wBarFast, atol=1e-16, err_msg="w wrt res, nuTilde-row seed")

    def test_repeated_calls(self):
        dwBar = self.CFDSolver.getStatePerturbation(314)

        wBarFast1 = self.CFDSolver.computeJacobianVectorProductBwdFast(resBar=dwBar)
        wBarFast2 = self.CFDSolver.computeJacobianVectorProductBwdFast(resBar=dwBar)

        np.testing.assert_allclose(wBarFast1, wBarFast2, atol=1e-16, err_msg="w wrt res double call")


@parameterized_class(test_params)
class TestDotProductsBCM(reg_test_classes.RegTest):
    """
    Fwd/rev transpose consistency for the SA-BCM case, both variants. Proves consistency, not
    correctness -- fwd and rev can both be wrong identically and still pass this test.
    Correctness comes from TestJacVecFwdBCMCS (complex build) only -- see
    tests/reg_tests/README_BCM.md's verification ladder.
    """

    N_PROCS = 2

    def setUp(self):
        if not hasattr(self, "name"):
            return

        super().setUp()

        options = copy.copy(adflowDefOpts)
        options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
        options.update(self.options)

        self.CFDSolver = ADFLOW(options=copy.deepcopy(options), debug=True)

        self.ap = copy.deepcopy(self.aero_prob)
        for dv in defaultAeroDVs:
            self.ap.addDV(dv)

        self.CFDSolver.getResidual(self.ap)

    def test_dot_products(self):
        utils.assert_dot_products_allclose(self.handler, self.CFDSolver, seed=314, tol=1e-8)


if __name__ == "__main__":
    unittest.main()
