# built-ins
import unittest
import os
import copy
from parameterized import parameterized_class
import numpy as np

# MACH classes
from adflow import ADFLOW, ADFLOW_C

from reg_default_options import adflowDefOpts
import reg_test_utils as utils
from baseclasses.testing import getTol

import reg_test_classes
import reg_sagr
from reg_sagr import ap_sagr_flatplate, sagrBaseOptions, sagrAeroDVs


baseDir = os.path.dirname(os.path.abspath(__file__))

test_params = [
    {
        "name": "sagr_flatplate",
        "options": copy.deepcopy(sagrBaseOptions),
        "ref_file": "jacvecfwd_sagr_flatplate.json",
        "aero_prob": copy.deepcopy(ap_sagr_flatplate),
        "N_PROCS": 1,
    },
]


@parameterized_class(test_params)
class TestJacVecFwdSAGR(reg_test_classes.RegTest):
    """
    Tests that given a flow state the FWD (tangent, Tapenade _d) jacobian
    vector products for the SA-noft2-Gamma-Retheta model agree with the
    values recorded in the ref file. Mirrors TestJacVecFwd
    (test_jacVecProdFWD.py) with additional per-variable (gamma, reThetat)
    and coupling-block checks from audit 06.
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
        # add the SA-GR aero dvs to the problem (mach/P/T exercise the
        # uInf/muInf vorticity-limiter path, audit-06 F1)
        for dv in sagrAeroDVs:
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

    # ------------------- SA-GR specific checks --------------------------------
    def test_coupling_blocks(self):
        # per-block dR[row]/dw[col] norms: gamma/reThetat in isolation plus
        # the SA<->transition cross blocks (audit 06 section 4)
        reg_sagr.assert_coupling_blocks_allclose(self.handler, self.CFDSolver, seed=314, tol=5e-9)

    def test_xDvDot_transition_rows(self):
        # transition-row norms of dR/d(mach, P, T, alpha): farfield/wInf BC
        # chain and the uInf/muInf limiter path (audit-06 F1)
        reg_sagr.assert_transition_xdvdot_allclose(self.handler, self.CFDSolver, self.ap, seed=1.0, tol=1e-10)


@parameterized_class(test_params)
class TestJacVecFwdSAGRFD(reg_test_classes.RegTest):
    """
    Tests that given a flow state the FWD jacobian vector products for the
    SA-noft2-Gamma-Retheta model agree with FD. Mirrors TestJacVecFwdFD.

    Note (audit sst_dev_lessons, watch item 2): the SA-GR residual has
    one-sided kinks (vorticity cap, smoothMinMax blend points); expect FD
    noise there rather than AD error. Do not loosen tolerances to pass —
    use the CS class below for a decisive check.
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
        for dv in sagrAeroDVs:
            self.ap.addDV(dv)

        # propagates the values from the restart file throughout the code
        self.CFDSolver.getResidual(self.ap)

    # ------------------- Derivative routine checks ----------------------------
    def test_wDot(self):
        # perturb each input and check that the outputs match the FD to with in reason
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

    def test_wDot_transition_columns(self):
        # same FD comparison but seeding gamma and reThetat columns in
        # isolation, so a coupling-block error is not diluted by the mean
        # flow entries
        nw, blocks = reg_sagr.getStateBlocks(self.CFDSolver)
        wDotFull = self.CFDSolver.getStatePerturbation(321)

        for colName in ["nuTilde", "gamma", "reThetat"]:
            wDot = reg_sagr.maskStateVector(wDotFull, nw, blocks[colName])

            resDot = self.CFDSolver.computeJacobianVectorProductFwd(wDot=wDot, residualDeriv=True)
            resDot_FD = self.CFDSolver.computeJacobianVectorProductFwd(
                wDot=wDot, residualDeriv=True, mode="FD", h=1e-8
            )

            np.testing.assert_allclose(resDot_FD, resDot, rtol=8e-4, err_msg="residual wrt %s" % colName)

    def test_xDvDot(self):
        # perturb each input and check that the outputs match the FD to with in reason
        step_size = {
            "alpha": 1e-4,
            "mach": 1e-5,
            "P": 1e-1,
            "T": 1e-4,
        }

        for aeroDV in self.ap.DVs.values():
            key = aeroDV.key
            xDvDot = {key: 1.0}

            resDot, funcsDot, fDot = self.CFDSolver.computeJacobianVectorProductFwd(
                xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True
            )

            resDot_FD, funcsDot_FD, fDot_FD = self.CFDSolver.computeJacobianVectorProductFwd(
                xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, mode="FD", h=step_size[key]
            )

            # the tolerances here are loose becuase different ouputs have different optimal steps
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
class TestJacVecFwdSAGRCS(reg_test_classes.CmplxRegTest):
    """
    Tests the SA-GR jacobian vector products against complex step. Mirrors
    TestJacVecFwdCS: the CS values are recorded under the same handler keys
    the real-build (AD) classes train into the ref file, so running this on
    the complex build compares AD against CS.
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

        self.ap = copy.deepcopy(self.aero_prob)
        # add the SA-GR aero dvs to the problem
        for dv in sagrAeroDVs:
            self.ap.addDV(dv)

        self.CFDSolver = ADFLOW_C(options=options, debug=True)

        # propagates the values from the restart file throughout the code
        self.CFDSolver.getResidual(self.ap)

    # ------------------- Derivative routine checks ----------------------------
    def cmplx_test_wDot(self):
        if not hasattr(self, "name"):
            return

        # perturb each input and check that the outputs match the FD to with in reason
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

        # perturb each input and check that the outputs match the FD to with in reason
        xVDot = self.CFDSolver.getSpatialPerturbation(314)

        rtol, atol = getTol(tol=1e-10)

        resDot_cs, funcsDot_cs, fDot_cs = self.CFDSolver.computeJacobianVectorProductFwd(
            xVDot=xVDot, residualDeriv=True, funcDeriv=True, fDeriv=True, mode="CS", h=self.h
        )

        self.handler.root_print("||dR/dXv * xVDot||")
        self.handler.par_add_norm("||dR/dXv * xVDot||", resDot_cs, rtol=rtol, atol=atol)

        # These can be finiky sometimes so a bigger tolerance.
        self.handler.root_print("dFuncs/dXv * xVDot")
        self.handler.root_add_dict("dFuncs/dXv * xVDot", funcsDot_cs, rtol=rtol * 10, atol=atol * 10)

        self.handler.root_print("||dF/dXv * xVDot||")
        self.handler.par_add_norm("||dF/dXv * xVDot||", fDot_cs, rtol=rtol, atol=atol)

    def cmplx_test_xDvDot(self):
        if not hasattr(self, "name"):
            return
        rtol, atol = getTol(tol=1e-10)

        # perturb each input and check that the outputs match the FD to with in reason
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

    # ------------------- SA-GR specific checks --------------------------------
    def cmplx_test_coupling_blocks(self):
        if not hasattr(self, "name"):
            return
        # same keys as TestJacVecFwdSAGR.test_coupling_blocks -> the ref file
        # trained by the AD run is compared block-by-block against CS
        reg_sagr.assert_coupling_blocks_allclose(
            self.handler, self.CFDSolver, seed=314, mode="CS", h=self.h, tol=5e-9
        )

    def cmplx_test_xDvDot_transition_rows(self):
        if not hasattr(self, "name"):
            return
        # decisive check for audit-06 F1: after the Tapenade rerun with
        # uInf/muInf active, the gamma/reThetat row norms of dR/d(mach,P,T)
        # from AD must match these CS values
        reg_sagr.assert_transition_xdvdot_allclose(
            self.handler, self.CFDSolver, self.ap, seed=1.0, mode="CS", h=self.h, tol=1e-10
        )


if __name__ == "__main__":
    unittest.main()
