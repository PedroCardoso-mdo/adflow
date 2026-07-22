#!/usr/bin/env python
"""Stage 1 of the plain-SA partial-derivative isolation ladder (see
docs/current-task.md): reverse-fast (_fast_b) vs full-reverse (_b)
consistency, using the raw ADflow API directly (no reg_sagr.py
block-splitting helpers -- that's "my script", tested separately at stages
2/3) on ADflow's own SA tutorial-wing mesh+restart (same fixture as
test_functionals.py's rans_tut_wing case). No upstream test file ships a
RANS/SA case for this check (test_jacVecProdBWDFast.py only has an Euler
param) so this is a minimal standalone equivalent.

Usage: mpirun -n 2 /home/mdo/packages_v2/mach/bin/python sanity_check_bwdfast_stage1.py
"""
import os
import copy

from mpi4py import MPI
from adflow import ADFLOW
import numpy

# dev/ is one level below reg_tests/; make the shared reg_* fixtures
# and sibling campaign modules importable when run directly (mpirun python).
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from reg_default_options import adflowDefOpts, defaultAeroDVs
from reg_aeroproblems import ap_tutorial_wing

baseDir = os.path.dirname(os.path.abspath(__file__))
rank = MPI.COMM_WORLD.rank

options = copy.copy(adflowDefOpts)
options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
options.update(
    {
        "gridfile": os.path.join(baseDir, "../../input_files/mdo_tutorial_rans_scalar_jst.cgns"),
        "restartfile": os.path.join(baseDir, "../../input_files/mdo_tutorial_rans_scalar_jst.cgns"),
        "mgcycle": "sg",
        "equationtype": "RANS",
        "smoother": "DADI",
        "cfl": 1.5,
        "cflcoarse": 1.25,
        "resaveraging": "never",
        "nsubiter": 3,
        "nsubiterturb": 3,
        "ncyclescoarse": 100,
        "ncycles": 1000,
        "usenksolver": True,
        "l2convergence": 1e-14,
        "l2convergencecoarse": 1e-4,
        "nkswitchtol": 1e-3,
        "adjointl2convergence": 1e-14,
        "frozenturbulence": False,
    }
)

CFDSolver = ADFLOW(options=options, debug=True)
ap = copy.deepcopy(ap_tutorial_wing)
for dv in defaultAeroDVs:
    ap.addDV(dv)

CFDSolver.getResidual(ap)

dwBar = CFDSolver.getStatePerturbation(314)
wBar = CFDSolver.computeJacobianVectorProductBwd(resBar=dwBar, wDeriv=True)
wBarFast = CFDSolver.computeJacobianVectorProductBwdFast(resBar=dwBar)

numpy.testing.assert_allclose(wBar, wBarFast, atol=1e-16, err_msg="w wrt res")

wBarFast2 = CFDSolver.computeJacobianVectorProductBwdFast(resBar=dwBar)
numpy.testing.assert_allclose(wBarFast, wBarFast2, atol=1e-16, err_msg="w wrt res, double call")

if rank == 0:
    print("\n=== STAGE 1 BWD-vs-BWDFast PASSED ===\n", flush=True)
