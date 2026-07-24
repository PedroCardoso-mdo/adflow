# built-ins
import unittest
import numpy
import os
import copy
from collections import defaultdict
from parameterized import parameterized_class
from mpi4py import MPI

# MACH classes
from pygeo import DVGeometry
from pyspline import Curve
from idwarp import USMesh

from adflow import ADFLOW
from adflow import ADFLOW_C
from idwarp import USMesh_C

import reg_test_utils as utils

from reg_default_options import adflowDefOpts, defaultAeroDVs, IDWarpDefOpts
from reg_bcm import ap_bcm_tut_wing, bcmBaseOptionsSmooth, bcmBaseOptionsHard, bcmFFDFile
import reg_test_classes

baseDir = os.path.dirname(os.path.abspath(__file__))


def getDVGeo(ffdFile, isComplex=False):
    # Identical to the sibling SA-GR repo's test_adjoint_sagr.py DVGeo: SA-BCM now runs on the
    # SAME tutorial-wing mesh + FFD (see reg_bcm.py case-inputs note), so this exercises the
    # exact same twist/span (global, ref-axis) and shape (local) design variables, driving the
    # full dR/dXv path through the SA-BCM-modified residual (wall distance, metrics, vorticity).
    DVGeo = DVGeometry(ffdFile, isComplex=isComplex)

    nTwist = 6
    DVGeo.addRefAxis(
        "wing",
        Curve(
            x=numpy.linspace(5.0 / 4.0, 1.5 / 4.0 + 7.5, nTwist),
            y=numpy.zeros(nTwist),
            z=numpy.linspace(0, 14, nTwist),
            k=2,
        ),
    )

    def twist(val, geo):
        for i in range(nTwist):
            geo.rot_z["wing"].coef[i] = val[i]

    def span(val, geo):
        C = geo.extractCoef("wing")
        s = geo.extractS("wing")
        for i in range(len(C)):
            C[i, 2] += s[i] * val[0]
        geo.restoreCoef(C, "wing")

    DVGeo.addGlobalDV("twist", [0] * nTwist, twist, lower=-10, upper=10, scale=1.0)
    DVGeo.addGlobalDV("span", [0], span, lower=-10, upper=10, scale=1.0)
    DVGeo.addLocalDV("shape", lower=-0.5, upper=0.5, axis="y", scale=10.0)

    return DVGeo


test_params = [
    {
        "name": "bcm_smooth_tut_wing",
        "options": copy.deepcopy(bcmBaseOptionsSmooth),
        "ref_file": "adjoint_bcm_smooth_tut_wing.json",
        "aero_prob": copy.deepcopy(ap_bcm_tut_wing),
        # cd is the transition-sensitive functional; never verify cl alone (CLAUDE.md /
        # README_BCM.md, mirrors the sst_dev post-mortem SA-GR's suite guards against too)
        "evalFuncs": ["cl", "cd", "cmz", "drag"],
        "N_PROCS": 2,
    },
    {
        "name": "bcm_hard_tut_wing",
        "options": copy.deepcopy(bcmBaseOptionsHard),
        "ref_file": "adjoint_bcm_hard_tut_wing.json",
        "aero_prob": copy.deepcopy(ap_bcm_tut_wing),
        "evalFuncs": ["cl", "cd", "cmz", "drag"],
        "N_PROCS": 2,
    },
]


@parameterized_class(test_params)
class TestAdjointBCM(reg_test_classes.RegTest):
    """
    Tests that total sensitivities calculated by solving the use_SABCM=True adjoint are correct,
    for both SABCM_Exp variants. Mirrors test_adjoint.py's TestAdjoint / the sibling SA-GR repo's
    TestAdjointSAGR.
    """

    N_PROCS = 2

    options = None
    ap = None
    ref_file = None

    def setUp(self):
        if not hasattr(self, "name"):
            return

        super().setUp()

        options = copy.copy(adflowDefOpts)
        options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
        options.update(self.options)

        self.ffdFile = bcmFFDFile

        mesh_options = copy.copy(IDWarpDefOpts)
        mesh_options.update({"gridFile": options["gridfile"]})

        self.ap = copy.deepcopy(self.aero_prob)
        self.ap.evalFuncs = self.evalFuncs

        for dv in defaultAeroDVs:
            self.ap.addDV(dv)

        self.CFDSolver = ADFLOW(options=options, debug=True)

        self.CFDSolver.setMesh(USMesh(options=mesh_options))
        self.CFDSolver.setDVGeo(getDVGeo(self.ffdFile, isComplex=False), pointSetKwargs={"embTol": 1e-12, "eps": 1e-14})

        # propagates the values from the restart file throughout the code
        self.CFDSolver.getResidual(self.ap)

    def test_residuals(self):
        utils.assert_residuals_allclose(self.handler, self.CFDSolver, self.ap, tol=1e-10)

    def test_adjoint(self):
        utils.assert_adjoint_sens_allclose(self.handler, self.CFDSolver, self.ap, tol=1e-10)
        self.assert_adjoint_failure()

    def test_adjoint2(self):
        utils.assert_adjoint2_sens_allclose(self.handler, self.CFDSolver, self.ap, tol=1e-10)
        self.assert_adjoint_failure()

    def test_adjoint_states(self):
        utils.assert_adjoint_states_allclose(self.handler, self.CFDSolver, self.ap, tol=1e-10)
        self.assert_adjoint_failure()


@parameterized_class(test_params)
class TestCmplxStepBCM(reg_test_classes.CmplxRegTest):
    """
    Decisive check for the manuscript's caveat: adjoint totals (trained by TestAdjointBCM)
    compared against complex-step re-convergence, for both variants. This is the check that
    either confirms the fix or reproduces the reported bug -- see docs/current-task.md.
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

        self.ffdFile = bcmFFDFile

        mesh_options = copy.copy(IDWarpDefOpts)
        mesh_options.update({"gridFile": options["gridfile"]})

        self.ap = copy.deepcopy(self.aero_prob)
        self.ap.evalFuncs = self.evalFuncs

        for dv in defaultAeroDVs:
            self.ap.addDV(dv)

        # Complex-build solver overrides -- same rationale as the sibling SA-GR repo's
        # TestCmplxStepSAGR: the complexify build excludes the Tapenade AD routines, so the
        # AD-based preconditioner (ANK/NKADPC=True) is unavailable; re-converge with FD-colored
        # PC on the decoupled ANK->NK path instead.
        options["ankadpc"] = False
        options["nkadpc"] = False
        options["ankcoupledswitchtol"] = 1e-16

        self.CFDSolver = ADFLOW_C(options=options, debug=True)

        self.CFDSolver.setMesh(USMesh_C(options=mesh_options))
        self.CFDSolver.setDVGeo(getDVGeo(self.ffdFile, isComplex=True), pointSetKwargs={"embTol": 1e-12, "eps": 1e-14})

        self.CFDSolver.getResidual(self.ap)

    def cmplx_test_aero_dvs(self):
        if not hasattr(self, "name"):
            return

        rtol = 1e-8
        atol = 5e-10

        funcsSens = defaultdict(lambda: {})

        for dv in ["alpha", "mach"]:
            setattr(self.ap, dv, getattr(self.ap, dv) + self.h * 1j)

            self.CFDSolver.resetFlow(self.ap)
            self.CFDSolver(self.ap, writeSolution=False)
            self.assert_solution_failure()

            funcs = {}
            self.CFDSolver.evalFunctions(self.ap, funcs)
            setattr(self.ap, dv, getattr(self.ap, dv) - self.h * 1j)

            for f in self.ap.evalFuncs:
                key = self.ap.name + "_" + f
                dv_key = dv + "_" + self.ap.name
                funcsSens[key][dv_key] = numpy.imag(funcs[key]) / self.h

        if MPI.COMM_WORLD.rank == 0:
            print("====================================")
            print(self.name, funcsSens)
            print("====================================")

        self.handler.root_add_dict("Eval Functions Sens:", funcsSens, rtol=rtol, atol=atol)

    def cmplx_test_geom_dvs(self):
        if not hasattr(self, "name"):
            return

        funcsSens = defaultdict(lambda: {})

        xRef = {"twist": [0.0] * 6, "span": [0.0], "shape": numpy.zeros(72, dtype="D")}

        rtol = 5e-9
        atol = 5e-9

        for dv in ["span", "twist", "shape"]:
            xRef[dv][0] += self.h * 1j

            self.CFDSolver.resetFlow(self.ap)
            self.CFDSolver.DVGeo.setDesignVars(xRef)
            self.CFDSolver(self.ap, writeSolution=False)
            self.assert_solution_failure()

            funcs = {}
            self.CFDSolver.evalFunctions(self.ap, funcs)

            xRef[dv][0] -= self.h * 1j

            for f in self.ap.evalFuncs:
                key = self.ap.name + "_" + f
                dv_key = dv
                funcsSens[key][dv_key] = numpy.imag(funcs[key]) / self.h

                err_msg = "Failed value for: {}".format(key + " " + dv_key)

                ref_val = self.handler.db["Eval Functions Sens:"][key][dv_key]
                if not isinstance(ref_val, float):
                    ref_val = ref_val.flatten()[0]

                numpy.testing.assert_allclose(funcsSens[key][dv_key], ref_val, atol=atol, rtol=rtol, err_msg=err_msg)

        if MPI.COMM_WORLD.rank == 0:
            print("====================================")
            print(self.name, funcsSens)
            print("====================================")


if __name__ == "__main__":
    unittest.main()
