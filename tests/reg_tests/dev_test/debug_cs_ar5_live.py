#!/usr/bin/env python
"""Standalone, non-testflo replica of TestCmplxStepAR5SA.cmplx_test_aero_dvs
(mach DV only), run directly under mpirun so ADflow's own solver log streams
to the terminal in real time -- testflo always buffers subprocess output
until the test finishes, even with -s, so it can't be used to watch a
stalling/diverging solve live.

Usage:
  mpirun -n 2 /home/mdo/packages_v2/mach/bin/python debug_cs_ar5_live.py [--no-restart] [--cold-start] [--ank-first] [--dv {alpha,mach,shape}] [--index N]

--dv shape --index N perturbs FFD local-shape control point N (of 112, see
getDVGeo in test_adjoint_ar5_sa.py) by h*1j via DVGeo.setDesignVars, one
point at a time, instead of alpha/mach. Compared against the trained
adjoint reference (refs/adjoint_ar5_sa.json, "shape" key, index N) at the
end automatically.

--ank-first (only with --cold-start): ADflow has no built-in "detect a
stall, then switch to NK" logic for plain SA (that only exists for SA-GR's
transitionNK/transitionNKStallRtolCap, which is transition-specific -- see
CLAUDE.md). So this is done as two explicit phases instead: phase 1 runs
with useNKSolver=False (ANK/CANK/CSANK ladder only) up to --ank-cycles
iterations; phase 2 re-enables NK (nkswitchtol forced to 1.0 so it engages
immediately) and continues the solve from wherever phase 1 left off (no
resetFlow -- the state carries over). If phase 1 already converged, phase 2
is a fast no-op; if it stalled, NK gets an actual shot at finishing it.

--no-restart drops restartfile so ADflow cold-starts from freestream instead
of the trained restart -- isolates whether the CS re-solve's non-convergence
is specific to the restart state or happens from a cold start too. Uses the
same (small) ncycles budget as the restart case, so it is NOT expected to
converge -- it's just for watching where a bare cold start with the CS
test's restart-tuned options gets stuck.

--cold-start does the real thing: sets mach = x0 + h*1j from the very start
(iteration 0) and runs the full validated production ladder (options
transplanted from run_strategy.py / STRATEGY.md, adapted from
SA-Gamma-Retheta to plain SA -- the SA-GR-only transition options are
dropped, turbResScale is left unset since _updateTurbResScale() already
auto-defaults it to 10000.0 for plain SA, pyADflow.py:6673) directly on the
perturbed problem. Because h is infinitesimal, the real part of the
converged solution is indistinguishable from the unperturbed case, and the
imaginary part accumulates the exact derivative through the whole nonlinear
solve -- one solve, no separate baseline+resolve needed. Read the derivative
straight off the converged complex functionals: Im(f)/h.
"""
import argparse
import copy
import os

import numpy
from mpi4py import MPI

from adflow import ADFLOW_C
from idwarp import USMesh_C

from reg_default_options import adflowDefOpts, IDWarpDefOpts
from test_adjoint_ar5_sa import ar5SAOptions, ap_ar5_sa, ar5FFDFile, getDVGeo

parser = argparse.ArgumentParser()
parser.add_argument("--no-restart", action="store_true")
parser.add_argument("--cold-start", action="store_true")
parser.add_argument("--ank-first", action="store_true")
parser.add_argument("--ank-cycles", type=int, default=100000)
parser.add_argument("--l2convergence", type=float, default=1e-20)
parser.add_argument("--nkswitchtol", type=float, default=1e-5)
parser.add_argument("--nk-cycles", type=int, default=200000)
parser.add_argument("--dv", choices=["alpha", "mach", "shape"], default="alpha")
parser.add_argument("--index", type=int, default=0, help="FFD shape point index (only for --dv shape)")
parser.add_argument("--no-nk", action="store_true", help="disable NK entirely -- ANK/CANK/CSANK ladder only")
parser.add_argument("--cycles", type=int, default=400000, help="ncycles for the single continuous cold-start solve")
args = parser.parse_args()

baseDir = os.path.dirname(os.path.abspath(__file__))
rank = MPI.COMM_WORLD.rank

options = copy.copy(adflowDefOpts)
options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
options.update(copy.deepcopy(ar5SAOptions))
options["ankadpc"] = False
options["nkadpc"] = False
# ar5SAOptions' l2convergence=1e-14 was tuned for perturbation-off-restart,
# not a cold start -- raised to production's own 1e-12 (run_strategy.py) so
# the run can actually satisfy it instead of chasing an unreachable target
# until it hits the ncycles ceiling.
options["l2convergence"] = args.l2convergence

if args.cold_start:
    options.pop("restartfile", None)
    # full validated ladder (run_strategy.py), adapted to plain SA:
    # useBlockettes/acousticScaleFactor differ from ar5SAOptions's defaults
    # and are treated as load-bearing there; turbResScale is intentionally
    # left alone (auto-defaults correctly for "SA").
    options["useBlockettes"] = False
    options["acousticScaleFactor"] = 0.2
    options["ncycles"] = args.cycles
    # skip CANK (coupled ANK) and CSANK/SANK (second-order switch) for this
    # SA run -- stay segregated first-order ANK until NK takes over. Set
    # the coupling/second-order tolerances low enough that ANK will never
    # reach them before nkswitchtol (below) preempts.
    options["ankcoupledswitchtol"] = 1e-16
    options["anksecondordswitchtol"] = 1e-16
    options["nkswitchtol"] = args.nkswitchtol
    if args.no_nk:
        options["usenksolver"] = False
    # this run uses turbulencemodel="SA" (plain SA), so every branch/dev
    # addition made to NKSolvers.F90 on this feature branch is already a
    # structural no-op here -- they're all guarded by
    # "turbModel == spalartallmarasnoft2gammaretheta .and. transitionNK
    # .and. transitionNKActive" (or a subset thereof), grep-confirmed with
    # no exceptions. transitionNK=False here is belt-and-suspenders: makes
    # that explicit rather than relying on the turbulence-model mismatch
    # alone, so this is provably stock/vanilla ADflow NK.
    options["transitionnk"] = False
    if args.ank_first:
        options["usenksolver"] = False
        options["ncycles"] = args.ank_cycles
    if rank == 0:
        print("\n=== --cold-start: freestream + full production ladder ===\n", flush=True)
elif args.no_restart:
    options.pop("restartfile", None)
    if rank == 0:
        print("\n=== --no-restart: cold-starting from freestream, no restartfile ===\n", flush=True)

mesh_options = copy.copy(IDWarpDefOpts)
mesh_options.update({"gridFile": options["gridFile"]})

ap = copy.deepcopy(ap_ar5_sa)
ap.evalFuncs = ["cl", "cd", "cmz", "drag"]
for dv in ["alpha", "mach", "P", "T"]:
    ap.addDV(dv)

CFDSolver = ADFLOW_C(options=options, debug=True)
CFDSolver.setMesh(USMesh_C(options=mesh_options))
CFDSolver.setDVGeo(getDVGeo(ar5FFDFile, isComplex=True), pointSetKwargs={"embTol": 1e-12, "eps": 1e-14})

h = 1e-10
dv = args.dv
xRef = None
if dv == "shape":
    xRef = {"shape": numpy.zeros(112, dtype="D")}
    x0 = 0.0
else:
    x0 = getattr(ap, dv)

if args.cold_start:
    # perturb by h*1j BEFORE the first solve -- one cold-start solve, done
    # directly on the perturbed problem; no separate baseline needed. Same
    # residual-triggered single-solve approach as alpha/mach -- no
    # iteration-count-gated phases (--ank-first is still there but unused
    # by default, per "i dont want it by number of iter").
    if dv == "shape":
        xRef["shape"][args.index] = h * 1j
        CFDSolver.DVGeo.setDesignVars(xRef)
        if rank == 0:
            print(f"\n=== cold-start CS solve: shape[{args.index}] = {h}j ===\n", flush=True)
    else:
        setattr(ap, dv, x0 + h * 1j)
        if rank == 0:
            print(f"\n=== cold-start CS solve: {dv} = {x0} + {h}j ===\n", flush=True)

    if args.ank_first:
        if rank == 0:
            print(f"\n=== phase 1: ANK only, up to {args.ank_cycles} cycles, NK disabled ===\n", flush=True)
        CFDSolver(ap, writeSolution=False)

        phase1Funcs = {}
        CFDSolver.checkSolutionFailure(ap, phase1Funcs)
        if rank == 0:
            print(f"\n=== phase 1 done: fail={phase1Funcs.get('fail')} -- now enabling NK ===\n", flush=True)

        # phase 2: hand off to NK from wherever phase 1 left off (no
        # resetFlow -- reuse the current state, converged or stalled). Give
        # NK a big budget here -- this phase is the one that needs to
        # actually drive the residual down, not just engage promptly.
        if rank == 0:
            print(f"\n=== phase 2: NK, up to {args.nk_cycles} cycles ===\n", flush=True)
        CFDSolver.setOption("useNKSolver", True)
        CFDSolver.setOption("nkswitchtol", 1.0)
        CFDSolver.setOption("ncycles", args.nk_cycles)
        CFDSolver(ap, writeSolution=False)
    else:
        CFDSolver(ap, writeSolution=False)
else:
    CFDSolver.getResidual(ap)
    if dv == "shape":
        xRef["shape"][args.index] = h * 1j
        CFDSolver.DVGeo.setDesignVars(xRef)
        if rank == 0:
            print(f"\n=== CS solve: shape[{args.index}] = {h}j ===\n", flush=True)
    else:
        setattr(ap, dv, x0 + h * 1j)
        if rank == 0:
            print(f"\n=== CS solve: {dv} + {h}j ===\n", flush=True)

    CFDSolver.resetFlow(ap)
    CFDSolver(ap, writeSolution=False)

funcs = {}
CFDSolver.checkSolutionFailure(ap, funcs)
if rank == 0:
    print(f"\n=== solve done: fail={funcs.get('fail')} ===\n", flush=True)

CFDSolver.evalFunctions(ap, funcs)
if dv != "shape":
    setattr(ap, dv, x0)

if rank == 0:
    import json

    refFile = os.path.join(baseDir, "refs", "adjoint_ar5_sa.json")
    adjSens = json.load(open(refFile))["Eval Functions Sens:"]
    for f in ap.evalFuncs:
        key = ap.name + "_" + f
        val = numpy.imag(funcs[key]) / h
        dvKey = dv + "_" + ap.name if dv != "shape" else "shape"
        adjRaw = adjSens.get(key, {}).get(dvKey)
        if dv == "shape" and isinstance(adjRaw, dict) and "__ndarray__" in adjRaw:
            adjVal = numpy.array(adjRaw["__ndarray__"]).flatten()[args.index]
        else:
            adjVal = adjRaw
        label = f"d{f}/d{dv}" + (f"[{args.index}]" if dv == "shape" else "")
        if adjVal is not None:
            adjVal = float(numpy.atleast_1d(adjVal).flatten()[0])
            err = abs(val - adjVal) / max(abs(adjVal), 1e-12)
            print(f"{label} (CS) = {val:.6e}   adjoint = {adjVal:.6e}   rel_err = {err:.3e}")
        else:
            print(f"{label} (CS) = {val:.6e}   adjoint = N/A")
