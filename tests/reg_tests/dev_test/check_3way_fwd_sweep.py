#!/usr/bin/env python
"""FD step-size sweep for the 3-way forward-mode check (real build only --
AD and CS were already confirmed to agree to full precision at both stages
in check_3way_fwd.py, so sweeping FD against AD here is equivalent to
sweeping against CS, without needing the complex build again).

Compares two FD formulas across an array of h:
  - one-sided (the built-in mode="FD" in computeJacobianVectorProductFwd,
    see src/adjoint/adjointDebug.F90:computeMatrixFreeProductFwdFD --
    (R(w+h*wDot) - R(w))/h, first-order accurate)
  - centered (manual, via getStates/setStates --
    (R(w+h*wDot) - R(w-h*wDot))/(2h), second-order accurate)
for the same wDot/probe seeds used in check_3way_fwd.py, so the AD value
printed there is the reference both are compared against.

Usage: mpirun -n 2 /home/mdo/packages_v2/mach/bin/python check_3way_fwd_sweep.py --mesh {tutorial,ar5}
"""
import argparse
import os
import copy

from mpi4py import MPI
import numpy
from adflow import ADFLOW

from reg_default_options import adflowDefOpts, defaultAeroDVs

parser = argparse.ArgumentParser()
parser.add_argument("--mesh", choices=["tutorial", "ar5"], default="tutorial")
parser.add_argument(
    "--turbmodel",
    choices=["sa", "sagr"],
    default="sa",
    help="sa = plain SA (nw=6); sagr = SA-noft2-Gamma-Retheta (nw=8, AR5 SA-GR "
    "restart). sagr forces the AR5 SA-GR case regardless of --mesh.",
)
args = parser.parse_args()

baseDir = os.path.dirname(os.path.abspath(__file__))
comm = MPI.COMM_WORLD
rank = comm.rank

options = copy.copy(adflowDefOpts)
options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])

if args.turbmodel == "sagr":
    import reg_sagr

    options.update(copy.deepcopy(reg_sagr.sagrBaseOptions))
    ap = copy.deepcopy(reg_sagr.ap_sagr_ar5_wing)
    for dv in reg_sagr.sagrAeroDVs:
        ap.addDV(dv)
elif args.mesh == "ar5":
    from test_adjoint_ar5_sa import ar5SAOptions, ap_ar5_sa

    options.update(copy.deepcopy(ar5SAOptions))
    ap = copy.deepcopy(ap_ar5_sa)
else:
    from reg_aeroproblems import ap_tutorial_wing

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
    ap = copy.deepcopy(ap_tutorial_wing)

if args.turbmodel != "sagr":  # sagr already added its own DVs above
    for dv in defaultAeroDVs:
        ap.addDV(dv)

CFDSolver = ADFLOW(options=options, debug=True)
CFDSolver.getResidual(ap)

# same seeds as check_3way_fwd.py -- AD/CS reference from that script applies here
wDot = CFDSolver.getStatePerturbation(321)
probe = CFDSolver.getStatePerturbation(999)


def reducedDot(vec):
    local = numpy.sum(vec * probe)
    total = comm.reduce(local)
    return total


resDot_AD = CFDSolver.computeJacobianVectorProductFwd(wDot=wDot, residualDeriv=True, mode="AD")
scalarAD = reducedDot(resDot_AD)

w0 = CFDSolver.getStates()

hVals = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12]

if rank == 0:
    print("\n=== FD step-size sweep, mesh=%s, turbmodel=%s ===" % (args.mesh, args.turbmodel), flush=True)
    print("AD reference (== CS, confirmed in check_3way_fwd.py) = %.10e\n" % scalarAD, flush=True)
    print("%-10s %-18s %-12s %-18s %-12s" % ("h", "one-sided FD", "rel_err", "centered FD", "rel_err"), flush=True)

for h in hVals:
    # one-sided (built-in)
    resDot_FD1 = CFDSolver.computeJacobianVectorProductFwd(wDot=wDot, residualDeriv=True, mode="FD", h=h)
    scalarFD1 = reducedDot(resDot_FD1)

    # centered (manual, via getStates/setStates)
    CFDSolver.setStates(w0 + h * wDot)
    resPlus = CFDSolver.getResidual(ap)
    CFDSolver.setStates(w0 - h * wDot)
    resMinus = CFDSolver.getResidual(ap)
    CFDSolver.setStates(w0)

    scalarPlus = reducedDot(resPlus)
    scalarMinus = reducedDot(resMinus)

    if rank == 0:
        scalarFD2 = (scalarPlus - scalarMinus) / (2 * h)
        relErr1 = abs(scalarFD1 - scalarAD) / abs(scalarAD)
        relErr2 = abs(scalarFD2 - scalarAD) / abs(scalarAD)
        print(
            "%-10.1e %-18.10e %-12.3e %-18.10e %-12.3e" % (h, scalarFD1, relErr1, scalarFD2, relErr2), flush=True
        )
