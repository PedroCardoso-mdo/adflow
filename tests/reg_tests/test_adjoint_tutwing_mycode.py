# built-ins
import unittest
import numpy
import os
import copy
from collections import defaultdict
from parameterized import parameterized_class

# MACH classes
from adflow import ADFLOW
from adflow import ADFLOW_C
from idwarp import USMesh
from idwarp import USMesh_C

import reg_test_utils as utils
from reg_default_options import adflowDefOpts, IDWarpDefOpts, defaultAeroDVs
import reg_test_classes

from test_adjoint import getDVGeo, test_params as adjoint_test_params
from reg_aeroproblems import ap_tutorial_wing

baseDir = os.path.dirname(os.path.abspath(__file__))

"""
"Mode 2" of the SA validation the user requested: my own test file/class
structure (mirrors test_adjoint_ar5_sa.py exactly) but pointed at the
ORIGINAL tutorial-wing grid/restart/FFD instead of AR5. Purpose: confirm my
code gives the same result as the literal official test_adjoint.py class on
the trusted mesh, before trusting the same code structure on AR5 (Mode 3,
test_adjoint_ar5_sa.py).
"""

ransParams = next(p for p in adjoint_test_params if p["name"] == "rans_tut_wing")
ffdFile = os.path.join(baseDir, "../../input_files/mdo_tutorial_ffd.fmt")

test_params = [
    {
        "name": "tutwing_mycode",
        "options": copy.deepcopy(ransParams["options"]),
        "ref_file": "adjoint_tutwing_mycode.json",
        "aero_prob": copy.deepcopy(ap_tutorial_wing),
        "evalFuncs": ["fx", "mz", "cl", "cd", "cmz", "lift", "drag", "cavitation", "colx", "coly", "colz"],
        "N_PROCS": 2,
    },
]


@parameterized_class(test_params)
class TestAdjointTutWingMyCode(reg_test_classes.RegTest):
    """Real-mode: residuals + adjoint totals -- same structure as
    TestAdjointAR5SA, mesh swapped to tutorial-wing."""

    N_PROCS = 2

    def setUp(self):
        if not hasattr(self, "name"):
            return
        super().setUp()

        options = copy.copy(adflowDefOpts)
        options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
        options.update(self.options)

        self.ffdFile = ffdFile
        mesh_options = copy.copy(IDWarpDefOpts)
        mesh_options.update({"gridFile": options["gridFile"]})

        self.ap = copy.deepcopy(self.aero_prob)
        self.ap.evalFuncs = self.evalFuncs
        for dv in defaultAeroDVs:
            self.ap.addDV(dv)

        self.CFDSolver = ADFLOW(options=options, debug=True)
        self.CFDSolver.setMesh(USMesh(options=mesh_options))
        self.CFDSolver.setDVGeo(getDVGeo(self.ffdFile, isComplex=False), pointSetKwargs={"embTol": 1e-12, "eps": 1e-14})

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
class TestCmplxStepTutWingMyCode(reg_test_classes.CmplxRegTest):
    """Complex-mode: CS vs adjoint totals -- same structure as
    TestCmplxStepAR5SA, mesh swapped to tutorial-wing."""

    N_PROCS = 2
    h = 1e-40

    def setUp(self):
        if not hasattr(self, "name"):
            return
        super().setUp()

        options = copy.copy(adflowDefOpts)
        options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
        options.update(self.options)
        options["ankadpc"] = False
        options["nkadpc"] = False

        self.ffdFile = ffdFile
        mesh_options = copy.copy(IDWarpDefOpts)
        mesh_options.update({"gridFile": options["gridFile"]})

        self.ap = copy.deepcopy(self.aero_prob)
        self.ap.evalFuncs = self.evalFuncs
        for dv in defaultAeroDVs:
            self.ap.addDV(dv)

        self.CFDSolver = ADFLOW_C(options=options, debug=True)
        self.CFDSolver.setMesh(USMesh_C(options=mesh_options))
        self.CFDSolver.setDVGeo(getDVGeo(self.ffdFile, isComplex=True), pointSetKwargs={"embTol": 1e-12, "eps": 1e-14})

        self.CFDSolver.getResidual(self.ap)

    def cmplx_test_aero_dvs(self):
        if not hasattr(self, "name"):
            return
        rtol = 1e-8
        atol = 5e-10

        for dv in ["alpha", "mach"]:
            funcsSens = defaultdict(lambda: {})
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

        self.handler.root_add_dict("Eval Functions Sens:", funcsSens, rtol=rtol, atol=atol)

    def cmplx_test_geom_dvs(self):
        if not hasattr(self, "name"):
            return
        funcsSens = defaultdict(lambda: {})
        xRef = {"twist": [0.0] * 6, "span": [0.0], "shape": numpy.zeros(72, dtype="D")}

        rtol = 5e-9
        atol = 5e-9

        for dv in ["twist", "span", "shape"]:
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

        self.handler.root_add_dict("Eval Functions Sens:", funcsSens, rtol=rtol, atol=atol)


if __name__ == "__main__":
    unittest.main()
