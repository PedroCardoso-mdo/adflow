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

# MACH testing class
from adflow import ADFLOW

from adflow import ADFLOW_C
from idwarp import USMesh_C


import reg_test_utils as utils

from reg_default_options import adflowDefOpts, IDWarpDefOpts

import reg_test_classes
from reg_sagr import ap_sagr_tut_wing, sagrBaseOptions, sagrAeroDVs, sagrFFDFile

baseDir = os.path.dirname(os.path.abspath(__file__))


def getDVGeo(ffdFile, isComplex=False):
    # Identical to test_adjoint.py's tutorial-wing DVGeo: the SA-GR case now
    # runs on the SAME tutorial-wing mesh + FFD as the SA adjoint test (see
    # reg_sagr.py case-inputs note), so we exercise the exact same twist/span
    # (global, ref-axis) and shape (local) design variables. These drive the
    # full dR/dXv path through the transition residuals (wall distance,
    # metrics, vorticity limiter) just as they do for plain SA.
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
        "name": "sagr_tut_wing",
        "options": copy.deepcopy(sagrBaseOptions),
        "ref_file": "adjoint_sagr_tut_wing.json",
        "aero_prob": copy.deepcopy(ap_sagr_tut_wing),
        # cd is the transition-sensitive functional; do NOT reduce this to
        # cl only (sst_dev post-mortem: adjoint regression covered only cl
        # and the verification hole was never closed)
        "evalFuncs": ["cl", "cd", "cmz", "drag"],
        "N_PROCS": 2,
    },
]


@parameterized_class(test_params)
class TestAdjointSAGR(reg_test_classes.RegTest):
    """
    Tests that total sensitivities calculated by solving the SA-GR adjoint
    (reverse mode, Tapenade _b routines, 8-state adjoint system) are correct.
    Mirrors TestAdjoint (test_adjoint.py), tutorial-wing RANS case.
    """

    N_PROCS = 2

    options = None
    ap = None
    ref_file = None

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

        self.ffdFile = sagrFFDFile

        mesh_options = copy.copy(IDWarpDefOpts)
        mesh_options.update({"gridFile": options["gridfile"]})

        self.ap = copy.deepcopy(self.aero_prob)

        # Setup aeroproblem
        self.ap.evalFuncs = self.evalFuncs

        # add the SA-GR aero dvs to the problem
        for dv in sagrAeroDVs:
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
class TestCmplxStepSAGR(reg_test_classes.CmplxRegTest):
    """
    Complex-step verification of the SA-GR total sensitivities: re-converge
    the complexified solver with a 1e-40j perturbation on each DV and compare
    imag(f)/h against the adjoint totals stored in the ref file by
    TestAdjointSAGR. Mirrors TestCmplxStep (test_adjoint.py).
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

        self.ffdFile = sagrFFDFile

        mesh_options = copy.copy(IDWarpDefOpts)
        mesh_options.update({"gridFile": options["gridfile"]})

        self.ap = copy.deepcopy(self.aero_prob)

        # Setup aeroproblem
        self.ap.evalFuncs = self.evalFuncs

        # add the SA-GR aero dvs to the problem
        for dv in sagrAeroDVs:
            self.ap.addDV(dv)

        # Complex-build solver overrides. The complexify build deliberately
        # excludes the Tapenade AD routines (audit 08: complexify runs over all
        # files MINUS adjoint/output{Forward,Reverse,ReverseFast}), so the
        # AD-based preconditioner used by the real solve (ANK/NKADPC=True in
        # sagrBaseOptions) is unavailable in the complex build -- attempting it
        # aborts with "Forward AD routines are not complexified". The coupled
        # ANK path also warns it may diverge on the stiff transition sources.
        # So re-converge the complex primal with FD-colored PC on the decoupled
        # (DADI-turb) ANK->NK path instead. Verified (2026-07-23) to reach the
        # same steady state the real AD-PC ladder converges to. This only
        # affects HOW the complex solver iterates to R(w)=0; the converged
        # state (hence the CS-vs-adjoint comparison) is unchanged.
        options["ankadpc"] = False
        options["nkadpc"] = False
        options["ankcoupledswitchtol"] = 1e-16  # never couple -> stay decoupled

        self.CFDSolver = ADFLOW_C(options=options, debug=True)

        self.CFDSolver.setMesh(USMesh_C(options=mesh_options))
        self.CFDSolver.setDVGeo(getDVGeo(self.ffdFile, isComplex=True), pointSetKwargs={"embTol": 1e-12, "eps": 1e-14})

        # propagates the values from the restart file throughout the code
        self.CFDSolver.getResidual(self.ap)

    def cmplx_test_aero_dvs(self):
        if not hasattr(self, "name"):
            # return immediately when the setup method is being called on the based class and NOT the
            # classes created using parametrized
            # this will happen when training, but will hopefully be fixed down the line
            return

        rtol = 1e-8
        atol = 5e-10

        funcsSens = defaultdict(lambda: {})

        # mach is the priority DV: it drives uInf/muInf through the P&Z
        # Eq. 52-53 vorticity limiter (audit-06 F1) and the farfield
        # wInf(itu2/itu3) ghost state
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
            # return immediately when the setup method is being called on the based class and NOT the
            # classes created using parametrized
            # this will happen when training, but will hopefully be fixed down the line
            return

        # redo the setup for a cmplx test
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
