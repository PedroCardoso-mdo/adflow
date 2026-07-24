# built-ins
import unittest
import os
import copy
from parameterized import parameterized_class
import numpy as np

# MACH classes
from adflow import ADFLOW

from reg_default_options import adflowDefOpts
from reg_sagr import ap_sagr_tut_wing, sagrBaseOptions, sagrAeroDVs


baseDir = os.path.dirname(os.path.abspath(__file__))

# Variable layout (NKSolvers.F90:getStates): variable-fastest, w(i,j,k,1:nw).
# For SA-GR nw = 8; l (0-based) = 5/6/7 = nuTilde/gamma/reThetat.
VAR_NAMES = ["rho", "rhou", "rhov", "rhow", "rhoE", "nuTilde", "gamma", "reThetat"]

test_params = [
    {
        "name": "sagr_tut_wing",
        "options": copy.deepcopy(sagrBaseOptions),
        "aero_prob": copy.deepcopy(ap_sagr_tut_wing),
        "N_PROCS": 2,
    },
]


@parameterized_class(test_params)
class TestBlocketteResidualSAGR(unittest.TestCase):
    """
    Block vs blockette residual-operator consistency for the
    SA-noft2-Gamma-Retheta model.

    ADflow builds the RANS residual through two interchangeable paths inside
    blocketteRes (src/NKSolver/blockette.F90), selected by `useBlockettes`:

      useBlockettes = False -> blockResCore     -> saGammaReTheta_block
                               (reference, src/turbulence/saGammaRetheta.F90)
      useBlockettes = True  -> blocketteResCore -> inlined saGammaRethetaSource
                               /Advection/Viscous/ResScale (cache-blocked copies)

    pyADflow force-disables blockettes for SA-GR (_updateTurbResScale) because
    the inlined kernels had drifted from the model. This test locks in that the
    two paths return the SAME residual for the SAME state w -- for every
    variable, including the two transition variables gamma (itu2) and
    reThetat (itu3).

    Regression guard: a sign-flip in the inlined first-order upwind advection
    (blockette saGammaRethetaAdvection) made the transition residuals differ
    from the block path by a factor ~2 while the meanflow stayed correct; the
    per-variable assert below is what pins that down (see
    docs/task-log/2026-07-24-blockette-sagr-residual-sync.md).

    Self-contained (no stored reference) -> inherits unittest.TestCase, not
    RegTest.
    """

    N_PROCS = 2

    def setUp(self):
        if not hasattr(self, "name"):
            # base (non-parameterized) class -- skip, mirrors the other suites
            return

        options = copy.copy(adflowDefOpts)
        options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
        options.update(self.options)
        # Exercise the crossflow D_scf source branch in BOTH paths (crossflow is
        # a standing part of the model; inert on this ~2D wing but the code runs).
        options["transitioncrossflow"] = True
        # We only want the residual operator, not a solve.
        options["ncycles"] = 1
        options["usenksolver"] = False
        options["useanksolver"] = False

        self.CFDSolver = ADFLOW(options=copy.deepcopy(options), debug=True)

        self.ap = copy.deepcopy(self.aero_prob)
        for dv in sagrAeroDVs:
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
        self.assertEqual(nw, 8, "expected nw = 8 for SA-GR")

        # A ~1% perturbation of the converged state makes the residual O(0.1)
        # so every source/advection/diffusion term is active -- a converged
        # state has residual ~1e-9 and would hide term-level mismatches.
        state0 = self.CFDSolver.getStates()
        rng = np.random.default_rng(12345)
        state = state0 * (1.0 + 1.0e-2 * rng.standard_normal(state0.shape))

        res_block = self._residual(False, state)
        res_blockette = self._residual(True, state)

        rb = res_block.reshape(-1, nw)
        rk = res_blockette.reshape(-1, nw)

        # Per-variable check: the transition-variable mismatch this guards
        # against left the meanflow correct, so a single whole-vector assert
        # could be dominated by the (correct) meanflow. rtol=1e-7 sits ~3
        # orders above the measured tiled-summation roundoff (~1e-10 rel);
        # atol covers the near-zero nuTilde residual.
        for l in range(nw):
            np.testing.assert_allclose(
                rk[:, l],
                rb[:, l],
                rtol=1e-7,
                atol=1e-11,
                err_msg=f"blockette != block residual for variable '{VAR_NAMES[l]}'",
            )


if __name__ == "__main__":
    unittest.main()
