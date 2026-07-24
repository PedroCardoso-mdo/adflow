# built-ins
import unittest
import os
import copy
from parameterized import parameterized_class
import numpy as np

# MACH classes
from adflow import ADFLOW, ADFLOW_C

from reg_default_options import adflowDefOpts, defaultAeroDVs
import reg_test_utils as utils
from baseclasses.testing import getTol

import reg_test_classes
from reg_bcm import ap_bcm_tut_wing, bcmBaseOptionsSmooth, bcmBaseOptionsHard, bcmAeroDVs
import reg_bcm

baseDir = os.path.dirname(os.path.abspath(__file__))

test_params = [
    {
        "name": "bcm_smooth_tut_wing",
        "options": copy.deepcopy(bcmBaseOptionsSmooth),
        "ref_file": "jacvecfwd_bcm_smooth_tut_wing.json",
        "aero_prob": copy.deepcopy(ap_bcm_tut_wing),
        "N_PROCS": 2,
    },
    {
        "name": "bcm_hard_tut_wing",
        "options": copy.deepcopy(bcmBaseOptionsHard),
        "ref_file": "jacvecfwd_bcm_hard_tut_wing.json",
        "aero_prob": copy.deepcopy(ap_bcm_tut_wing),
        "N_PROCS": 2,
    },
]


@parameterized_class(test_params)
class TestJacVecFwdBCM(reg_test_classes.RegTest):
    """
    Tests that given a converged use_SABCM=True flow state, the FWD jacobian vector products
    agree with the previous values recorded in the ref file. Mirrors test_jacVecProdFWD.py's
    TestJacVecFwd, plus SA-BCM-specific checks. Runs both SABCM_Exp variants (smooth/hard).

    Unlike SA-GR, SA-BCM adds no new state variable (nw=6, unchanged from plain SA), so the
    "coupling blocks" here are just meanflow<->nuTilde (2x2), not a 4x4 gamma/reThetat split --
    see reg_bcm.py's getStateBlocks. Even without a new state, tTgamma's dependence on
    rho/rlv/d2wall/chi(nuTilde) makes this cross-coupling nontrivial and worth locking down.
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

        # propagates the values from the restart file throughout the code
        self.CFDSolver.getResidual(self.ap)

    # ------------------- Derivative routine checks ----------------------------
    def test_wDot(self):
        utils.assert_fwd_mode_wdot_allclose(self.handler, self.CFDSolver, self.ap, seed=314, tol=5e-9)

    def test_xVDot(self):
        utils.assert_fwd_mode_xVDot_allclose(self.handler, self.CFDSolver, self.ap, seed=314, tol=1e-10)

    def test_xDvDot(self):
        utils.assert_fwd_mode_xDvDot_allclose(self.handler, self.CFDSolver, self.ap, seed=1.0, tol=1e-10)

    def test_xDvDot_bcm_rows(self):
        """SA-BCM-relevant freestream DVs only (alpha/mach/P/T -- the ones feeding Re_theta,
        Re_theta_c, term2 through vorticity and viscosity ratios)."""
        reg_bcm.assert_bcm_xdvdot_allclose(self.handler, self.CFDSolver, self.ap, seed=1.0, tol=1e-10)

    def test_coupling_blocks(self):
        # per-block dR[row]/dw[col] norms: meanflow<->nuTilde (SA-BCM has no gamma/reThetat)
        reg_bcm.assert_coupling_blocks_allclose(self.handler, self.CFDSolver, seed=314, tol=5e-9)


@parameterized_class(test_params)
class TestJacVecFwdBCMFD(reg_test_classes.RegTest):
    """
    Self-contained AD-vs-FD check for the SA-BCM case (both variants) -- no ref files needed.

    NOTE: expect FD noise near the tanh blend kink (manuscript Eq. 8, around term1+term2 == S0,
    smooth variant only) and near the KS max() aggregation (Eq. 10, both variants). Do not
    loosen tolerances to compensate -- that's exactly the tolerance-inflation failure mode
    CLAUDE.md/README_BCM.md warn against. Use TestJacVecFwdBCMCS (complex build) for the
    decisive check instead.
    """

    N_PROCS = 2
    no_train = True

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

    def test_wDot(self):
        wDot = self.CFDSolver.getStatePerturbation(321)

        resDot, funcsDot, fDot = self.CFDSolver.computeJacobianVectorProductFwd(
            wDot=wDot, residualDeriv=True, funcDeriv=True, fDeriv=True
        )
        resDot_FD, funcsDot_FD, fDot_FD = self.CFDSolver.computeJacobianVectorProductFwd(
            wDot=wDot, residualDeriv=True, funcDeriv=True, fDeriv=True, mode="FD", h=1e-8
        )

        np.testing.assert_allclose(resDot_FD, resDot, rtol=8e-4, err_msg="residual")

        for func in funcsDot:
            np.testing.assert_allclose(funcsDot_FD[func], funcsDot[func], rtol=1e-5, err_msg=func)

        np.testing.assert_allclose(fDot_FD, fDot, rtol=5e-4, err_msg="forces")

    def test_wDot_bcm_columns(self):
        """Mask the wDot seed to just nuTilde (the only state variable SA-BCM's modified
        production term acts through) so a coupling error in the modified term isn't diluted by
        mean-flow entries."""
        nw, blocks = reg_bcm.getStateBlocks(self.CFDSolver)
        wDotFull = self.CFDSolver.getStatePerturbation(321)
        wDot = reg_bcm.maskStateVector(wDotFull, nw, blocks["nuTilde"])

        resDot = self.CFDSolver.computeJacobianVectorProductFwd(wDot=wDot, residualDeriv=True)
        resDot_FD = self.CFDSolver.computeJacobianVectorProductFwd(wDot=wDot, residualDeriv=True, mode="FD", h=1e-8)

        np.testing.assert_allclose(resDot_FD, resDot, rtol=8e-4, err_msg="residual, nuTilde-only seed")

    def test_xDvDot(self):
        step_size = {"alpha": 1e-4, "mach": 1e-5, "P": 1e-1, "T": 1e-4}

        for aeroDV in self.ap.DVs.values():
            key = aeroDV.key
            if key not in step_size:
                continue
            xDvDot = {key: 1.0}

            resDot, funcsDot, fDot = self.CFDSolver.computeJacobianVectorProductFwd(
                xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True
            )
            resDot_FD, funcsDot_FD, fDot_FD = self.CFDSolver.computeJacobianVectorProductFwd(
                xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, mode="FD", h=step_size[key]
            )

            np.testing.assert_allclose(resDot_FD, resDot, atol=5e-5, err_msg=f"residual wrt {key}")

            for func in funcsDot:
                if np.abs(funcsDot[func]) <= 1e-16:
                    np.testing.assert_allclose(
                        funcsDot_FD[func], funcsDot[func], atol=5e-5, err_msg=f"{func} wrt {key}"
                    )
                else:
                    np.testing.assert_allclose(
                        funcsDot_FD[func], funcsDot[func], rtol=1e-3, err_msg=f"{func} wrt {key}"
                    )

            np.testing.assert_allclose(fDot_FD, fDot, atol=5e-7, err_msg=f"forces wrt {key}")


@parameterized_class(test_params)
class TestJacVecFwdBCMCS(reg_test_classes.CmplxRegTest):
    """
    Decisive AD-vs-complex-step check for the SA-BCM case (both variants). Records results under
    the same handler keys the real-build TestJacVecFwdBCM trained, so running it against those
    refs compares AD vs CS key-by-key -- this is the check that actually answers the manuscript's
    "adjoint results need revalidation" caveat, not the FD class above.
    """

    N_PROCS = 2
    h = 1e-40

    def setUp(self):
        if not hasattr(self, "name"):
            return
        super().setUp()

        options = copy.copy(adflowDefOpts)
        options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
        options.update(self.options)

        self.ap = copy.deepcopy(self.aero_prob)
        for dv in defaultAeroDVs:
            self.ap.addDV(dv)

        self.CFDSolver = ADFLOW_C(options=options, debug=True)
        self.CFDSolver.getResidual(self.ap)

    def cmplx_test_wDot(self):
        if not hasattr(self, "name"):
            return

        wDot = self.CFDSolver.getStatePerturbation(314)

        resDot_CS, funcsDot_CS, fDot_CS = self.CFDSolver.computeJacobianVectorProductFwd(
            wDot=wDot, residualDeriv=True, funcDeriv=True, fDeriv=True, mode="CS", h=self.h
        )

        rtol, atol = getTol(tol=1e-9)

        self.handler.root_print("||dR/dw * wDot||")
        self.handler.par_add_norm("||dR/dw * wDot||", resDot_CS, rtol=rtol, atol=atol)

        self.handler.root_print("dFuncs/dw * wDot")
        self.handler.root_add_dict("dFuncs/dw * wDot", funcsDot_CS, rtol=rtol, atol=atol)

        self.handler.root_print("||dF/dw * wDot||")
        self.handler.par_add_norm("||dF/dw * wDot||", fDot_CS, rtol=rtol, atol=atol)

    def cmplx_test_xVDot(self):
        if not hasattr(self, "name"):
            return

        xVDot = self.CFDSolver.getSpatialPerturbation(314)
        rtol, atol = getTol(tol=1e-10)

        resDot_cs, funcsDot_cs, fDot_cs = self.CFDSolver.computeJacobianVectorProductFwd(
            xVDot=xVDot, residualDeriv=True, funcDeriv=True, fDeriv=True, mode="CS", h=self.h
        )

        self.handler.root_print("||dR/dXv * xVDot||")
        self.handler.par_add_norm("||dR/dXv * xVDot||", resDot_cs, rtol=rtol, atol=atol)

        self.handler.root_print("dFuncs/dXv * xVDot")
        self.handler.root_add_dict("dFuncs/dXv * xVDot", funcsDot_cs, rtol=rtol * 10, atol=atol * 10)

        self.handler.root_print("||dF/dXv * xVDot||")
        self.handler.par_add_norm("||dF/dXv * xVDot||", fDot_cs, rtol=rtol, atol=atol)

    def cmplx_test_xDvDot(self):
        if not hasattr(self, "name"):
            return
        rtol, atol = getTol(tol=1e-10)

        for aeroDV in self.ap.DVs.values():
            key = aeroDV.key
            self.handler.root_print("  -> %s" % key)
            xDvDot = {key: 1.0}

            resDot_cs, funcsDot_cs, fDot_cs = self.CFDSolver.computeJacobianVectorProductFwd(
                xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, mode="CS", h=self.h
            )

            self.handler.root_print("||dR/d%s||" % key)
            self.handler.par_add_norm("||dR/d%s||" % key, resDot_cs, rtol=rtol, atol=atol)

            self.handler.root_print("dFuncs/d%s" % key)
            self.handler.root_add_dict("dFuncs/d%s" % key, funcsDot_cs, rtol=rtol, atol=atol)

            self.handler.root_print("||dF/d%s||" % key)
            self.handler.par_add_norm("||dF/d%s||" % key, fDot_cs, rtol=rtol, atol=atol)

    def cmplx_test_xDvDot_bcm_rows(self):
        if not hasattr(self, "name"):
            return
        reg_bcm.assert_bcm_xdvdot_allclose(self.handler, self.CFDSolver, self.ap, seed=1.0, mode="CS", h=self.h, tol=1e-10)

    def cmplx_test_coupling_blocks(self):
        if not hasattr(self, "name"):
            return
        # same keys as TestJacVecFwdBCM.test_coupling_blocks -> the ref file trained by the AD
        # run is compared block-by-block against CS
        reg_bcm.assert_coupling_blocks_allclose(self.handler, self.CFDSolver, seed=314, mode="CS", h=self.h, tol=5e-9)


if __name__ == "__main__":
    unittest.main()
