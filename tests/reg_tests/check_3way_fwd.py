#!/usr/bin/env python
"""The 3-way forward-mode check (see docs/current-task.md /
docs/ADFLOW_BASE/debugging_derivatives.md): perturb one input direction
(wDot, a fixed state perturbation) and compute the same directional
derivative three independent ways at a single fixed restart state (no
reconvergence):
  - analytic FWD routine  (computeJacobianVectorProductFwd mode="AD")
  - finite difference     (mode="FD", h=1e-8)
  - complex step          (mode="CS", h=1e-40, complex build only)

Run once with --build real (prints AD and FD) and once with --build complex
(prints CS) -- ADFLOW and ADFLOW_C can't be instantiated in the same
process. Both builds reduce dR/dw*wDot against the same fixed probe vector
(getStatePerturbation with a second, independent seed) to a single
MPI-reduced scalar, so the three numbers are directly comparable across the
two separate runs/processes without needing to save/compare full arrays.

Usage:
  mpirun -n 2 <real  python> check_3way_fwd.py --mesh {tutorial,ar5} --build real
  mpirun -n 2 <cmplx python> check_3way_fwd.py --mesh {tutorial,ar5} --build complex
"""
import argparse
import os
import copy

from mpi4py import MPI
import numpy

from reg_default_options import adflowDefOpts, defaultAeroDVs

parser = argparse.ArgumentParser()
parser.add_argument("--mesh", choices=["tutorial", "ar5"], default="tutorial")
parser.add_argument("--build", choices=["real", "complex"], required=True)
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

if args.build == "real":
    from adflow import ADFLOW as ADFLOWClass
else:
    from adflow import ADFLOW_C as ADFLOWClass

CFDSolver = ADFLOWClass(options=options, debug=True)
CFDSolver.getResidual(ap)

# same seeds on both builds -> directly comparable reduced scalar
wDot = CFDSolver.getStatePerturbation(321)
probe = CFDSolver.getStatePerturbation(999)


def reducedDot(vec):
    local = numpy.sum(vec * probe)
    total = comm.reduce(local)
    return total


if args.build == "real":
    resDot_AD = CFDSolver.computeJacobianVectorProductFwd(wDot=wDot, residualDeriv=True, mode="AD")
    resDot_FD = CFDSolver.computeJacobianVectorProductFwd(wDot=wDot, residualDeriv=True, mode="FD", h=1e-8)

    scalarAD = reducedDot(resDot_AD)
    scalarFD = reducedDot(resDot_FD)

    if rank == 0:
        relErr = abs(scalarFD - scalarAD) / max(abs(scalarAD), 1e-300)
        print("\n=== 3-way fwd check, mesh=%s, turbmodel=%s, build=real ===" % (args.mesh, args.turbmodel), flush=True)
        print("  AD = %.10e" % scalarAD, flush=True)
        print("  FD = %.10e   rel_err(FD,AD) = %.3e" % (scalarFD, relErr), flush=True)
        print("  (CS: run --build complex and compare against AD above)", flush=True)
else:
    resDot_CS = CFDSolver.computeJacobianVectorProductFwd(wDot=wDot, residualDeriv=True, mode="CS", h=1e-40)
    scalarCS = reducedDot(resDot_CS)

    if rank == 0:
        print("\n=== 3-way fwd check, mesh=%s, turbmodel=%s, build=complex ===" % (args.mesh, args.turbmodel), flush=True)
        print("  CS = %.10e" % numpy.real(scalarCS), flush=True)
        print("  (compare against AD/FD from the --build real run)", flush=True)
