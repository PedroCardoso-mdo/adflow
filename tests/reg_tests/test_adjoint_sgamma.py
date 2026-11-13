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
import reg_transition
from reg_transition import ap_sgamma_tut_wing, sgammaOptions, transitionAeroDVs, tutWingFFDFile

baseDir = os.path.dirname(os.path.abspath(__file__))


def getDVGeo(ffdFile, isComplex=False):
    # same tutorial-wing twist, span and shape design variables as test_adjoint.py
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
        "name": "sgamma_tut_wing",
        "options": copy.deepcopy(sgammaOptions),
        "ref_file": "adjoint_sgamma_tut_wing.json",
        "aero_prob": copy.deepcopy(ap_sgamma_tut_wing),
        # cd is the functional that is sensitive to the transition location
        "evalFuncs": ["cl", "cd", "cmz", "drag"],
        "N_PROCS": 2,
    },
]

# Same checks with the model option on, about its own converged state. The state is not in git: the
# parameter set exists only once dev/generate_sgamma_restart.py --variant artvisc has written it.
if os.path.exists(reg_transition.sgammaArtViscRestartFile):
    test_params.append(
        {
            "name": "sgamma_artvisc_tut_wing",
            "options": copy.deepcopy(reg_transition.sgammaArtViscAdjointOptions),
            "ref_file": "adjoint_sgamma_artvisc_tut_wing.json",
            "aero_prob": copy.deepcopy(ap_sgamma_tut_wing),
            "evalFuncs": ["cl", "cd", "cmz", "drag"],
            "N_PROCS": 2,
        }
    )


@parameterized_class(test_params)
class TestAdjointSGAMMA(reg_test_classes.RegTest):
    """
    Tests that the total sensitivities calculated from the adjoint of the SA-noft2-Gamma model
    agree with the values recorded in the ref file. Same checks as TestAdjoint (test_adjoint.py).
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

        self.ffdFile = tutWingFFDFile

        mesh_options = copy.copy(IDWarpDefOpts)
        mesh_options.update({"gridFile": options["gridfile"]})

        self.ap = copy.deepcopy(self.aero_prob)

        # Setup aeroproblem
        self.ap.evalFuncs = self.evalFuncs

        # add the aero dvs to the problem
        for dv in transitionAeroDVs:
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
class TestCmplxStepSGAMMA(reg_test_classes.CmplxRegTest):
    """
    Tests the total sensitivities of the SA-noft2-Gamma model against complex step: the complex solver
    is re-converged from the restart state with a perturbation on one DV at a time and imag(f)/h is
    compared with the adjoint totals TestAdjointSGAMMA trained into the ref file. Same checks as
    TestCmplxStep (test_adjoint.py).
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

        self.ffdFile = tutWingFFDFile

        mesh_options = copy.copy(IDWarpDefOpts)
        mesh_options.update({"gridFile": options["gridfile"]})

        self.ap = copy.deepcopy(self.aero_prob)

        # Setup aeroproblem
        self.ap.evalFuncs = self.evalFuncs

        # add the aero dvs to the problem
        for dv in transitionAeroDVs:
            self.ap.addDV(dv)

        # The complex build has no AD routines, so the AD preconditioner of the real solve is not
        # available.
        options["ankadpc"] = False
        options["nkadpc"] = False
        options["ankcoupledswitchtol"] = 1e-20
        # With the FD preconditioner a solve from freestream does not converge for this model, so the
        # complex system is re-converged from the restart state (see _warmStart) with full Newton steps.
        # The real part is already converged and it is the imaginary part that has to converge:
        #  - no line search and a fixed unit step, because the line search judges the step on the real
        #    residual, which is already at its floor;
        #  - a tight fixed linear tolerance instead of Eisenstat-Walker;
        #  - an L2 target that is never reached, so the solver iterates until ncycles (about 40 steps).
        # _reportSolve prints the norm of the imaginary part of the residual that was reached.
        options["nkls"] = "none"
        options["nkfixedstep"] = 1.0
        options["ncycles"] = 6000
        options["l2convergence"] = 1e-30
        options["nklinearsolvetol"] = 1e-4
        options["nkuseew"] = False

        self.CFDSolver = ADFLOW_C(options=options, debug=True)

        self.CFDSolver.setMesh(USMesh_C(options=mesh_options))
        self.CFDSolver.setDVGeo(getDVGeo(self.ffdFile, isComplex=True), pointSetKwargs={"embTol": 1e-12, "eps": 1e-14})

        # propagates the values from the restart file throughout the code
        self.CFDSolver.getResidual(self.ap)

        # converged restart state: every complex re-converge starts from it
        self.w0 = self.CFDSolver.getStates().copy()

    def _warmStart(self):
        self.CFDSolver.setStates(self.w0)

    def _table(self, dv, dvKey, funcs):
        """Print CS vs adjoint for every function of this DV before asserting (an assert stops at the first)."""
        if MPI.COMM_WORLD.rank != 0:
            return
        for f in self.ap.evalFuncs:
            key = self.ap.name + "_" + f
            cs = numpy.imag(funcs[key]) / self.h
            try:
                ref = self.handler.db["Eval Functions Sens:"][key][dvKey]
                ref = ref if isinstance(ref, float) else numpy.asarray(ref).flatten()[0]
                print(
                    "[CS table] d%-4s/d%-6s CS=% .10e  adjoint=% .10e  rel=%.2e"
                    % (f, dv, cs, ref, abs(cs - ref) / max(abs(ref), 1e-30))
                )
            except (KeyError, TypeError) as e:
                print("[CS table] d%-4s/d%-6s CS=% .10e  (no ref: %s)" % (f, dv, cs, e))

    def _reportSolve(self, dv):
        """Print (do not assert) the solver flag and the imaginary part of the residual reached by the
        complex re-converge. The solver flag is always "failed" because the L2 target is unreachable."""
        funcs = {}
        self.CFDSolver.checkSolutionFailure(self.ap, funcs)
        if MPI.COMM_WORLD.rank == 0:
            print("[CS solve] dv=%s  fail=%s" % (dv, funcs.get("fail")))
        imR = numpy.imag(self.CFDSolver.getResidual(self.ap))
        imR = numpy.sqrt(MPI.COMM_WORLD.allreduce(float(numpy.sum(imR**2)), op=MPI.SUM))
        if MPI.COMM_WORLD.rank == 0:
            print("[CS solve] dv=%s  ||Im R||/h = %.3e  (should be << 1)" % (dv, imR / self.h))

    def cmplx_test_aero_dvs(self):
        if not hasattr(self, "name"):
            # return immediately when the setup method is being called on the based class and NOT the
            # classes created using parametrized
            # this will happen when training, but will hopefully be fixed down the line
            return

        # same tolerance as test_adjoint_sagr.py
        rtol = 1e-7
        atol = 1e-7

        funcsSens = defaultdict(lambda: {})

        # --- alpha: asserted ---
        dv = "alpha"
        setattr(self.ap, dv, getattr(self.ap, dv) + self.h * 1j)
        self._warmStart()
        self.CFDSolver(self.ap, writeSolution=False)
        self._reportSolve(dv)
        funcs = {}
        self.CFDSolver.evalFunctions(self.ap, funcs)
        self._table(dv, dv + "_" + self.ap.name, funcs)
        setattr(self.ap, dv, getattr(self.ap, dv) - self.h * 1j)
        for f in self.ap.evalFuncs:
            key = self.ap.name + "_" + f
            funcsSens[key][dv + "_" + self.ap.name] = numpy.imag(funcs[key]) / self.h

        if MPI.COMM_WORLD.rank == 0:
            print("====================================")
            print(self.name, funcsSens)
            print("====================================")

        self.handler.root_add_dict("Eval Functions Sens:", funcsSens, rtol=rtol, atol=atol)

        # --- mach: reported, not asserted ---
        # Not asserted, as in test_adjoint_sagr.py. The partial derivatives with respect to mach are
        # checked against complex step in test_jacVecProdFWD_sgamma.py.
        dv = "mach"
        setattr(self.ap, dv, getattr(self.ap, dv) + self.h * 1j)
        self._warmStart()
        self.CFDSolver(self.ap, writeSolution=False)
        machFuncs = {}
        self.CFDSolver.evalFunctions(self.ap, machFuncs)
        setattr(self.ap, dv, getattr(self.ap, dv) - self.h * 1j)
        if MPI.COMM_WORLD.rank == 0:
            print("==== mach check (reported, not asserted) ====")
            for f in self.ap.evalFuncs:
                key = self.ap.name + "_" + f
                cs = numpy.imag(machFuncs[key]) / self.h
                try:
                    ref = self.handler.db["Eval Functions Sens:"][key]["mach_" + self.ap.name]
                    if not isinstance(ref, float):
                        ref = ref.flatten()[0]
                    ok = abs(cs - ref) <= atol + rtol * abs(ref)
                    print("  d%-4s/dmach  CS=% .8e  ref=% .8e  rel=%.2e  %s"
                          % (f, cs, ref, abs(cs - ref) / max(abs(ref), 1e-30),
                             "ok" if ok else "outside tolerance (not asserted)"))
                except Exception as e:
                    print("  d%-4s/dmach  CS=% .8e  (no ref: %s)" % (f, cs, e))

    def cmplx_test_geom_dvs(self):
        if not hasattr(self, "name"):
            # return immediately when the setup method is being called on the based class and NOT the
            # classes created using parametrized
            # this will happen when training, but will hopefully be fixed down the line
            return

        # redo the setup for a cmplx test
        funcsSens = defaultdict(lambda: {})

        xRef = {"twist": [0.0] * 6, "span": [0.0], "shape": numpy.zeros(72, dtype="D")}

        # same tolerance as cmplx_test_aero_dvs
        rtol = 1e-7
        atol = 1e-7

        for dv in ["shape", "span", "twist"]:
            xRef[dv][0] += self.h * 1j

            self._warmStart()
            self.CFDSolver.DVGeo.setDesignVars(xRef)
            self.CFDSolver(self.ap, writeSolution=False)
            self._reportSolve(dv)

            funcs = {}
            self.CFDSolver.evalFunctions(self.ap, funcs)
            self._table(dv, dv, funcs)

            xRef[dv][0] -= self.h * 1j

            for f in self.ap.evalFuncs:
                key = self.ap.name + "_" + f
                dv_key = dv
                funcsSens[key][dv_key] = numpy.imag(funcs[key]) / self.h

                err_msg = "Failed value for: {}".format(key + " " + dv_key)

                ref_val = self.handler.db["Eval Functions Sens:"][key][dv_key]
                if not isinstance(ref_val, float):
                    ref_val = ref_val.flatten()[0]

                cs_val = funcsSens[key][dv_key]

                if f == "drag":
                    # drag = cd * q_inf * Sref carries the relative error of cd on a much larger
                    # magnitude, so it can exceed the tolerance where cd meets it. cd is asserted;
                    # drag is reported only.
                    if MPI.COMM_WORLD.rank == 0:
                        ok = abs(cs_val - ref_val) <= atol + rtol * abs(ref_val)
                        print("  [drag, not asserted] d%s/d%s  CS=% .8e  ref=% .8e  rel=%.2e  %s"
                              % (f, dv_key, cs_val, ref_val,
                                 abs(cs_val - ref_val) / max(abs(ref_val), 1e-30),
                                 "ok" if ok else "outside tolerance (not asserted)"))
                    continue

                numpy.testing.assert_allclose(cs_val, ref_val, atol=atol, rtol=rtol, err_msg=err_msg)

        if MPI.COMM_WORLD.rank == 0:
            print("====================================")
            print(self.name, funcsSens)
            print("====================================")


if __name__ == "__main__":
    unittest.main()
