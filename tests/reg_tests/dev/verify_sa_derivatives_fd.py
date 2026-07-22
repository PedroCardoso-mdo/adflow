#!/usr/bin/env python
"""Standalone FD-vs-adjoint total-derivative verification for the PLAIN SA
model (mdo_tutorial_rans_scalar_jst.cgns, the well-established tutorial-wing
case) using the existing tutorial FFD (input_files/mdo_tutorial_ffd.fmt).

Purpose: validate the derivative-verification METHOD itself -- central FD on
re-converged functionals vs. the adjoint's evalFunctionsSens totals -- on a
model that's known-good, before trusting the same method on SA-Gamma-Retheta
(which doesn't have a working FFD yet, see reg_sagr.py's sagrFFDFile note).

Mirrors the mechanics of test_adjoint.py's TestCmplxStep (same rans_tut_wing
options block, same getDVGeo/FFD setup, same DVGeo.setDesignVars pattern for
geometric DVs) but with real central finite differences instead of complex
step, since this is the FD pass requested first.

Usage: mpirun -np 4 /home/mdo/packages_v2/mach/bin/python verify_sa_derivatives_fd.py
"""
import copy
import os

import numpy as np
from mpi4py import MPI

from adflow import ADFLOW
from idwarp import USMesh

# dev/ is one level below reg_tests/; make the shared reg_* fixtures
# and sibling campaign modules importable when run directly (mpirun python).
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from reg_default_options import adflowDefOpts, defaultAeroDVs, IDWarpDefOpts
from reg_aeroproblems import ap_tutorial_wing
from test_adjoint import getDVGeo, test_params as adjoint_test_params

baseDir = os.path.dirname(os.path.abspath(__file__))
rank = MPI.COMM_WORLD.rank

# Reuse the already-validated "rans_tut_wing" options block from
# test_adjoint.py verbatim -- this is the established SA RANS case the
# upstream adjoint regression tests already trust.
ransParams = next(p for p in adjoint_test_params if p["name"] == "rans_tut_wing")

options = copy.copy(adflowDefOpts)
options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
options.update(copy.deepcopy(ransParams["options"]))

ffdFile = os.path.join(baseDir, "../../input_files/mdo_tutorial_ffd.fmt")
mesh_options = copy.copy(IDWarpDefOpts)
mesh_options.update({"gridFile": options["gridFile"]})

ap = copy.deepcopy(ap_tutorial_wing)
ap.evalFuncs = ["cl", "cd"]
for dv in defaultAeroDVs:
    ap.addDV(dv)

CFDSolver = ADFLOW(options=options, debug=True)
CFDSolver.setMesh(USMesh(options=mesh_options))
CFDSolver.setDVGeo(getDVGeo(ffdFile, isComplex=False), pointSetKwargs={"embTol": 1e-12, "eps": 1e-14})


def report(name, fd, adj, rtol):
    fd = float(np.atleast_1d(fd).flatten()[0])
    adj = float(np.atleast_1d(adj).flatten()[0])
    denom = max(abs(adj), 1e-12)
    err = abs(fd - adj) / denom
    status = "PASS" if err < rtol else "FAIL"
    if rank == 0:
        print(f"[{status}] {name:14s} FD={fd: .8e}  adjoint={adj: .8e}  rel_err={err:.3e}  (rtol={rtol:.0e})")
    return status == "PASS"


results = []

# ---- baseline solve + adjoint sensitivities ----
if rank == 0:
    print("=== baseline solve ===")
CFDSolver(ap)
funcsSens = {}
CFDSolver.evalFunctionsSens(ap, funcsSens, evalFuncs=ap.evalFuncs)

if rank == 0:
    print("funcsSens keys:", {k: list(v.keys()) for k, v in funcsSens.items()})

# ---- aero DVs: alpha, mach (central FD, re-converged each side) ----
if rank == 0:
    print("\n=== aero DVs (alpha, mach) ===")
h_aero = {"alpha": 1e-4, "mach": 1e-5}
for dv in ["alpha", "mach"]:
    h = h_aero[dv]
    x0 = getattr(ap, dv)

    setattr(ap, dv, x0 + h)
    CFDSolver(ap, writeSolution=False)
    funcsP = {}
    CFDSolver.evalFunctions(ap, funcsP)

    setattr(ap, dv, x0 - h)
    CFDSolver(ap, writeSolution=False)
    funcsM = {}
    CFDSolver.evalFunctions(ap, funcsM)

    setattr(ap, dv, x0)
    CFDSolver(ap, writeSolution=False)  # restore baseline state for the next DV

    for f in ap.evalFuncs:
        key = ap.name + "_" + f
        dvkey = dv + "_" + ap.name
        fd = (funcsP[key] - funcsM[key]) / (2 * h)
        adj = funcsSens[key][dvkey]
        results.append(report(f"d{f}/d{dv}", fd, adj, rtol=1e-3))

# ---- geometric DVs via FFD: twist, span, shape (central FD) ----
if rank == 0:
    print("\n=== geometric DVs via FFD (twist, span, shape) ===")
h_geo = {"twist": 1e-3, "span": 1e-4, "shape": 1e-4}
xRef = {"twist": [0.0] * 6, "span": [0.0], "shape": np.zeros(72)}

for dv in ["twist", "span", "shape"]:
    h = h_geo[dv]

    xRef[dv][0] = h
    CFDSolver.DVGeo.setDesignVars(xRef)
    CFDSolver(ap, writeSolution=False)
    funcsP = {}
    CFDSolver.evalFunctions(ap, funcsP)

    xRef[dv][0] = -h
    CFDSolver.DVGeo.setDesignVars(xRef)
    CFDSolver(ap, writeSolution=False)
    funcsM = {}
    CFDSolver.evalFunctions(ap, funcsM)

    xRef[dv][0] = 0.0
    CFDSolver.DVGeo.setDesignVars(xRef)
    CFDSolver(ap, writeSolution=False)  # restore

    for f in ap.evalFuncs:
        key = ap.name + "_" + f
        fd = (funcsP[key] - funcsM[key]) / (2 * h)
        adj_arr = np.atleast_1d(np.asarray(funcsSens[key][dv]))
        if rank == 0:
            print(f"  (debug) funcsSens[{key}][{dv}] shape={adj_arr.shape}")
        adj = adj_arr.flatten()[0]
        results.append(report(f"d{f}/d{dv}[0]", fd, adj, rtol=1e-2))

if rank == 0:
    n_pass = sum(results)
    print(f"\n{n_pass}/{len(results)} checks passed")

# ---- step-size sweep on "shape" (local DV, index 0) to separate a real
# derivative bug from FD truncation/cancellation noise ----
if rank == 0:
    print("\n=== shape[0] step-size sweep ===")
adj_arr = np.atleast_1d(np.asarray(funcsSens["mdo_tutorial_cl"]["shape"]))
adj_cl_shape0 = adj_arr.flatten()[0]
for h in [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]:
    xRef["shape"][0] = h
    CFDSolver.DVGeo.setDesignVars(xRef)
    CFDSolver(ap, writeSolution=False)
    funcsP = {}
    CFDSolver.evalFunctions(ap, funcsP)

    xRef["shape"][0] = -h
    CFDSolver.DVGeo.setDesignVars(xRef)
    CFDSolver(ap, writeSolution=False)
    funcsM = {}
    CFDSolver.evalFunctions(ap, funcsM)

    fd = (funcsP["mdo_tutorial_cl"] - funcsM["mdo_tutorial_cl"]) / (2 * h)
    if rank == 0:
        err = abs(fd - adj_cl_shape0) / max(abs(adj_cl_shape0), 1e-12)
        print(f"  h={h:.0e}  FD_dcl/dshape[0]={fd: .8e}  adjoint={adj_cl_shape0: .8e}  rel_err={err:.3e}")

xRef["shape"][0] = 0.0
CFDSolver.DVGeo.setDesignVars(xRef)
CFDSolver(ap, writeSolution=False)
