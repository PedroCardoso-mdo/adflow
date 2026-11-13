# built-ins
import unittest
import os
import copy
from parameterized import parameterized_class
import numpy as np

# MACH classes
from adflow import ADFLOW

from reg_default_options import adflowDefOpts
from reg_transition import ap_sagr_tut_wing, sagrOptions, transitionAeroDVs


baseDir = os.path.dirname(os.path.abspath(__file__))

# State layout: variable-fastest, w(i,j,k,1:nw), with nw = 8 for this model
VAR_NAMES = ["rho", "rhou", "rhov", "rhow", "rhoE", "nuTilde", "gamma", "reThetat"]

test_params = [
    {
        "name": "sagr_tut_wing",
        "options": copy.deepcopy(sagrOptions),
        "aero_prob": copy.deepcopy(ap_sagr_tut_wing),
        "N_PROCS": 2,
    },
]


@parameterized_class(test_params)
class TestBlocketteResidualSAGR(unittest.TestCase):
    """
    Tests that the two residual paths of blocketteRes (src/NKSolver/blockette.F90) return the same
    residual for the same state with the SA-noft2-Gamma-Retheta model:

      useBlockettes = False -> blockResCore, which calls the model routines in
                               src/turbulence/saGammaRetheta.F90 (reference)
      useBlockettes = True  -> blocketteResCore, which has its own cache-blocked copies of the
                               source, advection and viscous kernels

    pyADflow switches blockettes off for the transition models, so this test sets the option after
    the solver is created. The comparison is done per variable: an error in the transition kernels
    leaves the mean flow residual unchanged and would be hidden in a whole-vector norm.

    No reference file is needed, so this inherits unittest.TestCase and not RegTest.
    """

    N_PROCS = 2

    def setUp(self):
        if not hasattr(self, "name"):
            # return immediately when the setup method is being called on the based class and NOT the
            # classes created using parametrized
            return

        options = copy.copy(adflowDefOpts)
        options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
        options.update(self.options)
        # run the crossflow source branch in both paths (the term is active on this swept wing: see the
        # sagr_crossflow_tut_wing derivative tests)
        options["transitioncrossflow"] = True
        # only the residual is evaluated, there is no solve
        options["ncycles"] = 1
        options["usenksolver"] = False
        options["useanksolver"] = False

        self.CFDSolver = ADFLOW(options=copy.deepcopy(options), debug=True)

        self.ap = copy.deepcopy(self.aero_prob)
        for dv in transitionAeroDVs:
            self.ap.addDV(dv)

        # propagates the values from the restart file throughout the code
        self.CFDSolver.getResidual(self.ap)

    def _residual(self, useBlockettes, state):
        # select the residual path, set the same state and evaluate
        self.CFDSolver.setOption("useBlockettes", useBlockettes)
        self.assertEqual(self.CFDSolver.getOption("useBlockettes"), useBlockettes)
        self.CFDSolver.setStates(state)
        return self.CFDSolver.getResidual(self.ap).copy()

    def test_block_vs_blockette_residual(self):
        nw = int(self.CFDSolver.adflow.flowvarrefstate.nw)
        self.assertEqual(nw, 8, "expected nw = 8 for SA-noft2-Gamma-Retheta")

        # Perturb the converged state by about 1 % so the residual is O(0.1) and every source,
        # advection and diffusion term contributes; at the converged state the terms cancel.
        state0 = self.CFDSolver.getStates()
        rng = np.random.default_rng(12345)
        state = state0 * (1.0 + 1.0e-2 * rng.standard_normal(state0.shape))

        res_block = self._residual(False, state)
        res_blockette = self._residual(True, state)

        rb = res_block.reshape(-1, nw)
        rk = res_blockette.reshape(-1, nw)

        # The two paths sum the fluxes in a different order: the round-off difference is about 1e-10
        # relative. atol covers the near-zero nuTilde residual.
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
