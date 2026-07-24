# built-ins
import unittest
import os
import copy
from parameterized import parameterized_class
import numpy as np

# MACH classes
from adflow import ADFLOW

from reg_default_options import adflowDefOpts
from reg_bcm import ap_bcm_tut_wing, bcmBaseOptionsSmooth, bcmBaseOptionsHard, bcmAeroDVs


baseDir = os.path.dirname(os.path.abspath(__file__))

# Variable layout (NKSolvers.F90:getStates): variable-fastest, w(i,j,k,1:nw).
# SA-BCM is plain-SA-shaped: nw = 6, l (0-based) = 5 = nuTilde. No gamma/reThetat --
# SA-BCM adds no new transported state, see reg_bcm.py module docstring.
VAR_NAMES = ["rho", "rhou", "rhov", "rhow", "rhoE", "nuTilde"]

test_params = [
    {
        "name": "bcm_smooth_tut_wing",
        "options": copy.deepcopy(bcmBaseOptionsSmooth),
        "aero_prob": copy.deepcopy(ap_bcm_tut_wing),
        "N_PROCS": 2,
    },
    {
        "name": "bcm_hard_tut_wing",
        "options": copy.deepcopy(bcmBaseOptionsHard),
        "aero_prob": copy.deepcopy(ap_bcm_tut_wing),
        "N_PROCS": 2,
    },
]


@parameterized_class(test_params)
class TestBlocketteResidualBCM(unittest.TestCase):
    """
    Block vs blockette residual-operator consistency for SA-BCM, both SABCM_Exp variants.

    ADflow builds the RANS residual through two interchangeable paths inside
    blocketteRes (src/NKSolver/blockette.F90), selected by `useBlockettes`:

      useBlockettes = False -> blockResCore     -> saSource
                               (reference, src/turbulence/sa.F90, Tapenade-differentiated)
      useBlockettes = True  -> blocketteResCore -> saSource
                               (hand-maintained duplicate, blockette.F90:1004-1250)

    Unlike SA-GR (force-disabled after test_blockette_sagr.py caught drift), nothing
    force-disables useBlockettes for SA-BCM, and reg_bcm.py deliberately keeps it True
    by default (matching the previously-used AR5 run script). This test is the check
    that answers whether that default is actually safe -- see dev/diag_blockette_bcm.py
    for the raw-output version of the same check, run first.

    Self-contained (no stored reference) -> inherits unittest.TestCase, not RegTest.
    """

    N_PROCS = 2

    def setUp(self):
        if not hasattr(self, "name"):
            # base (non-parameterized) class -- skip, mirrors the other suites
            return

        options = copy.copy(adflowDefOpts)
        options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
        options.update(self.options)
        # We only want the residual operator, not a solve.
        options["ncycles"] = 1
        options["usenksolver"] = False
        options["useanksolver"] = False

        self.CFDSolver = ADFLOW(options=copy.deepcopy(options), debug=True)

        self.ap = copy.deepcopy(self.aero_prob)
        for dv in bcmAeroDVs:
            self.ap.addDV(dv)

        # propagate the restart state throughout the code
        self.CFDSolver.getResidual(self.ap)

    def _residual(self, useBlockettes, state):
        # force the desired path, restore the identical owned state, evaluate
        self.CFDSolver.setOption("useBlockettes", useBlockettes)
        self.assertEqual(self.CFDSolver.getOption("useBlockettes"), useBlockettes)
        self.CFDSolver.setStates(state)
        return self.CFDSolver.getResidual(self.ap).copy()

    def test_block_vs_blockette_residual(self):
        nw = int(self.CFDSolver.adflow.flowvarrefstate.nw)
        self.assertEqual(nw, 6, "expected nw = 6 for SA-BCM (plain-SA state layout)")

        # A ~1% perturbation of the converged state makes the residual O(0.1) so
        # every source/diffusion term is active -- a converged state has residual
        # ~1e-9 and would hide term-level mismatches. Same rationale as SA-GR's
        # equivalent test.
        state0 = self.CFDSolver.getStates()
        rng = np.random.default_rng(12345)
        state = state0 * (1.0 + 1.0e-2 * rng.standard_normal(state0.shape))

        res_block = self._residual(False, state)
        res_blockette = self._residual(True, state)

        rb = res_block.reshape(-1, nw)
        rk = res_blockette.reshape(-1, nw)

        # Per-variable check: a localized mismatch (e.g. only nuTilde, or only near
        # the transition front) could be diluted by the correct meanflow rows in a
        # single whole-vector assert.
        for l in range(nw):
            np.testing.assert_allclose(
                rk[:, l],
                rb[:, l],
                rtol=1e-7,
                atol=1e-10,
                err_msg="block vs blockette residual mismatch on %s (SABCM_Exp=%s)"
                % (VAR_NAMES[l], self.options["sabcm_exp"]),
            )


if __name__ == "__main__":
    unittest.main()
