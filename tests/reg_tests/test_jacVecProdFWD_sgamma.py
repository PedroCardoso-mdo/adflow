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
import reg_transition
from reg_transition import ap_sgamma_tut_wing, sgammaOptions, sgammaArtViscOptions, transitionAeroDVs


baseDir = os.path.dirname(os.path.abspath(__file__))

test_params = [
    {
        "name": "sgamma_tut_wing",
        "options": copy.deepcopy(sgammaOptions),
        "ref_file": "jacvecfwd_sgamma_tut_wing.json",
        "aero_prob": copy.deepcopy(ap_sgamma_tut_wing),
        "N_PROCS": 1,
    },
    # same state, artificial viscosity on gamma switched on
    {
        "name": "sgamma_artvisc_tut_wing",
        "options": copy.deepcopy(sgammaArtViscOptions),
        "ref_file": "jacvecfwd_sgamma_artvisc_tut_wing.json",
        "aero_prob": copy.deepcopy(ap_sgamma_tut_wing),
        "N_PROCS": 1,
    },
]


@parameterized_class(test_params)
class TestJacVecFwdSGAMMA(reg_test_classes.RegTest):
    """
    Tests that given a flow state the FWD jacobian vector products of the SA-noft2-Gamma model
    agree with the values recorded in the ref file. Same checks as TestJacVecFwd (test_jacVecProdFWD.py)
    plus the products split by variable block (mean flow, nuTilde, gamma).
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

    # ------------------- Derivative routine checks ----------------------------
    def test_wDot(self):
        utils.assert_fwd_mode_wdot_allclose(self.handler, self.CFDSolver, self.ap, seed=314, tol=5e-9)

    def test_xVDot(self):
        utils.assert_fwd_mode_xVDot_allclose(self.handler, self.CFDSolver, self.ap, seed=314, tol=1e-10)

    def test_xDvDot(self):
        utils.assert_fwd_mode_xDvDot_allclose(self.handler, self.CFDSolver, self.ap, seed=1.0, tol=1e-10)

    # ------------------- Transition model checks ------------------------------
    def test_coupling_blocks(self):
        # norm of every dR[row block]/dw[column block] product
        reg_transition.assert_coupling_blocks_allclose(self.handler, self.CFDSolver, seed=314, tol=5e-9)

    def test_xDvDot_transition_rows(self):
        # norm of dR/d(alpha, mach, P, T) split by residual row block
        reg_transition.assert_transition_xdvdot_allclose(self.handler, self.CFDSolver, self.ap, seed=1.0, tol=1e-10)


@parameterized_class(test_params)
class TestJacVecFwdSGAMMAFD(reg_test_classes.RegTest):
    """
    Tests that given a flow state the FWD jacobian vector products of the SA-noft2-Gamma model
    agree with FD. Same checks as TestJacVecFwdFD (test_jacVecProdFWD.py).

    The residual has one-sided switches (limiters, blending functions) and spans about 13 orders of
    magnitude, so the residual comparisons of the state products are limited by FD noise and are marked expectedFailure.
    The complex-step class below is the step-free check.
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
    # The residual comparisons sweep the FD step (assert_fd_allclose_hsweep) and report the best step
    # reached. No step meets the tolerance on this case, hence expectedFailure; an "unexpected success"
    # means the decorator can be dropped.
    @unittest.expectedFailure
    def test_wDot(self):
        # perturb each input and check that the outputs match the FD to with in reason
        wDot = self.CFDSolver.getStatePerturbation(321)

        resDot, funcsDot, fDot = self.CFDSolver.computeJacobianVectorProductFwd(
            wDot=wDot, residualDeriv=True, funcDeriv=True, fDeriv=True
        )
        _, funcsDot_FD, fDot_FD = self.CFDSolver.computeJacobianVectorProductFwd(
            wDot=wDot, residualDeriv=True, funcDeriv=True, fDeriv=True, mode="FD", h=1e-8
        )

        reg_transition.assert_fd_allclose_hsweep(
            lambda h: self.CFDSolver.computeJacobianVectorProductFwd(
                wDot=wDot, residualDeriv=True, mode="FD", h=h
            ),
            resDot,
            rtol=8e-4,
            hPrimary=1e-8,
            comm=self.CFDSolver.comm,
            root_print=self.handler.root_print,
            err_msg="residual (dR/dw * wDot)",
        )

        for func in funcsDot:
            np.testing.assert_allclose(funcsDot_FD[func], funcsDot[func], rtol=1e-5, err_msg=func)

        np.testing.assert_allclose(fDot_FD, fDot, rtol=5e-4, err_msg="forces")

    @unittest.expectedFailure
    def test_wDot_transition_columns(self):
        # same comparison seeding one turbulence variable at a time, so that an error in one block is
        # not diluted by the mean flow entries
        nw, blocks = reg_transition.getStateBlocks(self.CFDSolver)
        wDotFull = self.CFDSolver.getStatePerturbation(321)

        for colName in ["nuTilde", "gamma"]:
            wDot = reg_transition.maskStateVector(wDotFull, nw, blocks[colName])

            resDot = self.CFDSolver.computeJacobianVectorProductFwd(wDot=wDot, residualDeriv=True)

            reg_transition.assert_fd_allclose_hsweep(
                lambda h, wDot=wDot: self.CFDSolver.computeJacobianVectorProductFwd(
                    wDot=wDot, residualDeriv=True, mode="FD", h=h
                ),
                resDot,
                rtol=8e-4,
                hPrimary=1e-8,
                comm=self.CFDSolver.comm,
                root_print=self.handler.root_print,
                err_msg="residual wrt %s" % colName,
            )

    # the FD of the aero DV products meets the tolerance for this model: no expectedFailure
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

            _, funcsDot_FD, fDot_FD = self.CFDSolver.computeJacobianVectorProductFwd(
                xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, mode="FD", h=step_size[key]
            )

            # the tolerances here are loose becuase different ouputs have different optimal steps
            reg_transition.assert_fd_allclose_hsweep(
                lambda h, xDvDot=xDvDot: self.CFDSolver.computeJacobianVectorProductFwd(
                    xDvDot=xDvDot, residualDeriv=True, mode="FD", h=h
                ),
                resDot,
                atol=5e-5,
                hPrimary=step_size[key],
                comm=self.CFDSolver.comm,
                root_print=self.handler.root_print,
                err_msg=f"residual wrt {key}",
            )

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
class TestJacVecFwdSGAMMACS(reg_test_classes.CmplxRegTest):
    """
    Tests the FWD jacobian vector products of the SA-noft2-Gamma model against complex step.
    Same checks as TestJacVecFwdCS (test_jacVecProdFWD.py): the CS values are compared with the AD
    values the real build trained into the ref file.
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
        # add the aero dvs to the problem
        for dv in transitionAeroDVs:
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

    # ------------------- Transition model checks ------------------------------
    def cmplx_test_coupling_blocks(self):
        if not hasattr(self, "name"):
            return
        # same keys as TestJacVecFwdSGAMMA.test_coupling_blocks
        reg_transition.assert_coupling_blocks_allclose(
            self.handler, self.CFDSolver, seed=314, mode="CS", h=self.h, tol=5e-9
        )

    def cmplx_test_xDvDot_transition_rows(self):
        if not hasattr(self, "name"):
            return
        # same keys as TestJacVecFwdSGAMMA.test_xDvDot_transition_rows
        reg_transition.assert_transition_xdvdot_allclose(
            self.handler, self.CFDSolver, self.ap, seed=1.0, mode="CS", h=self.h, tol=1e-10
        )


if __name__ == "__main__":
    unittest.main()
