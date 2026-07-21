#!/usr/bin/env python
"""Real-mode FD step-size sweep for dcd/dmach on the AR5 plain-wing mesh,
plain SA (test_adjoint_ar5_sa.py's case). Reconstructed 2026-07-21 -- the
original ad-hoc version this was based on was never committed; only its
log survived (see docs/current-task.md's h-sweep table, and
/tmp/sweep_h_fd.log for the run that produced it).

Purpose: separate "the adjoint is wrong" from "this mesh can't be converged
tightly enough for FD/CS to resolve the true derivative" for the AR5 CS
mismatch (dcd/dmach: CS=-231.1 vs adjoint=-0.537, see docs/current-task.md).
Each +h/-h solve does a full resetFlow + cold re-solve (ncycles budget
below), matching TestCmplxStepAR5SA's own resetFlow-based methodology so
the comparison is apples-to-apples with the CS check.

Usage: mpirun -np 2 /home/mdo/packages_v2/mach/bin/python sweep_h_fd.py [--ncycles N]
"""
import argparse
import copy
import os

from mpi4py import MPI

from adflow import ADFLOW
from idwarp import USMesh

from reg_default_options import adflowDefOpts, IDWarpDefOpts
from test_adjoint_ar5_sa import ar5SAOptions, ap_ar5_sa, ar5FFDFile, getDVGeo

baseDir = os.path.dirname(os.path.abspath(__file__))
rank = MPI.COMM_WORLD.rank

parser = argparse.ArgumentParser()
parser.add_argument("--ncycles", type=int, default=5000)
parser.add_argument("--hlist", type=float, nargs="+", default=[1e-3, 1e-4, 1e-5, 1e-6])
args = parser.parse_args()

options = copy.copy(adflowDefOpts)
options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
options.update(copy.deepcopy(ar5SAOptions))
options["ncycles"] = args.ncycles

mesh_options = copy.copy(IDWarpDefOpts)
mesh_options.update({"gridFile": options["gridFile"]})

ap = copy.deepcopy(ap_ar5_sa)
ap.evalFuncs = ["cl", "cd", "cmz", "drag"]
for dv in ["alpha", "mach", "P", "T"]:
    ap.addDV(dv)

CFDSolver = ADFLOW(options=options, debug=True)
CFDSolver.setMesh(USMesh(options=mesh_options))
CFDSolver.setDVGeo(getDVGeo(ar5FFDFile, isComplex=False), pointSetKwargs={"embTol": 1e-12, "eps": 1e-14})

# baseline solve (establishes restart-consistent state) + adjoint dcd/dmach
CFDSolver(ap)
funcsSens = {}
CFDSolver.evalFunctionsSens(ap, funcsSens, evalFuncs=ap.evalFuncs)
adj = funcsSens["ar5_plain_wing_sa_cd"]["mach_ar5_plain_wing_sa"]
if rank == 0:
    print(f"\n=== adjoint dcd/dmach = {adj: .6f} ===\n")

x0 = ap.mach
for h in args.hlist:
    setattr(ap, "mach", x0 + h)
    CFDSolver.resetFlow(ap)
    CFDSolver(ap, writeSolution=False)
    failP = ap.solveFailed or ap.fatalFail
    funcsP = {}
    CFDSolver.evalFunctions(ap, funcsP)
    cdP = funcsP["ar5_plain_wing_sa_cd"]

    setattr(ap, "mach", x0 - h)
    CFDSolver.resetFlow(ap)
    CFDSolver(ap, writeSolution=False)
    failM = ap.solveFailed or ap.fatalFail
    funcsM = {}
    CFDSolver.evalFunctions(ap, funcsM)
    cdM = funcsM["ar5_plain_wing_sa_cd"]

    setattr(ap, "mach", x0)

    fd = (cdP - cdM) / (2 * h)
    if rank == 0:
        print(f"h={h:.0e}  failP={bool(failP)} failM={bool(failM)}  cdP={cdP:.8f}  cdM={cdM:.8f}  FD={fd:.6f}  adj={adj:.6f}")
