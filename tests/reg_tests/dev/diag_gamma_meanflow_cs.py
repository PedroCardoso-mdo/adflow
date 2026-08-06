#!/usr/bin/env python
"""Localize the 2.15e-8 mean-flow leak. Seed gamma-only and reThetat-only via
CS; extract the mean-flow rows of the CS product. Report per-component norm,
nonzero-cell count, and -- decisively -- whether the gamma-seed and
reThetat-seed mean-flow vectors are ELEMENT-WISE IDENTICAL (=> shared-scratch
artifact, not a physical coupling which would differ by seed)."""
import os, sys, copy

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mpi4py import MPI
import numpy as np
from adflow import ADFLOW_C
import reg_sagr
from reg_default_options import adflowDefOpts

comm = MPI.COMM_WORLD
rank = comm.rank
baseDir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
options = copy.copy(adflowDefOpts)
options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
options.update(copy.deepcopy(reg_sagr.sagrBaseOptions))
ap = copy.deepcopy(reg_sagr.ap_sagr_ar5_wing)
for dv in reg_sagr.sagrAeroDVs:
    ap.addDV(dv)
CFD = ADFLOW_C(options=options, debug=True)
CFD.getResidual(ap)

nw, blocks = reg_sagr.getStateBlocks(CFD)
wDotFull = CFD.getStatePerturbation(314)


def mask(vec, offs):
    out = np.zeros_like(vec)
    for o in offs:
        out[o::nw] = vec[o::nw]
    return out


mf_vecs = {}
compNames = ["rho", "rhou", "rhov", "rhow", "rhoE"]
for col in ["gamma", "reThetat"]:
    wDot = mask(wDotFull, blocks[col])
    resDot = np.real(CFD.computeJacobianVectorProductFwd(wDot=wDot, residualDeriv=True, mode="CS", h=1e-40))
    mf = mask(resDot, blocks["meanflow"])
    mf_vecs[col] = mf
    percomp = []
    for l in blocks["meanflow"]:
        c = resDot[l::nw]
        percomp.append(np.sqrt(comm.reduce(np.sum(c ** 2)) or 0.0))
    locmax = np.max(np.abs(mf)) if mf.size else 0.0
    gmax = comm.reduce(locmax, op=MPI.MAX)
    gnnz = comm.reduce(int(np.sum(np.abs(mf) > 1e-14)))
    if rank == 0:
        print("\n=== seed %s -> meanflow rows (CS) ===" % col)
        print("  max=%.6e  nonzero_cells(>1e-14)=%d" % (gmax, gnnz))
        print("  per-comp norm %s = %s" % (compNames, ["%.4e" % v for v in percomp]))

# decisive: are the two mean-flow leak vectors element-wise identical?
diff = mf_vecs["gamma"] - mf_vecs["reThetat"]
loc_absmax_g = np.max(np.abs(mf_vecs["gamma"])) if mf_vecs["gamma"].size else 0.0
loc_absmax_diff = np.max(np.abs(diff)) if diff.size else 0.0
g = comm.reduce(loc_absmax_g, op=MPI.MAX)
d = comm.reduce(loc_absmax_diff, op=MPI.MAX)
if rank == 0:
    print("\n=== gamma-seed vs reThetat-seed meanflow leak: element-wise identical? ===")
    print("  max|gamma_leak|              = %.6e" % g)
    print("  max|gamma_leak - ret_leak|   = %.6e" % d)
    print("  => %s" % ("IDENTICAL (shared-scratch artifact, seed-independent)"
                        if d < 1e-16 * max(g, 1e-300) else "DIFFER (seed-dependent -> possible real coupling)"))
