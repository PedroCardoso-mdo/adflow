#!/usr/bin/env python
"""Option (b): find a state reset that clears the cross-call CS residue while
preserving the converged (real) linearization point. After a meanflow-column
CS warm-up (which injects ~2e-8 into subsequent zero mean-flow rows), apply
each reset and measure ||dR[meanflow]/dw[gamma]|| (should be exact 0)."""
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
w0 = np.real(CFD.getStates()).copy()  # the true real converged state


def mask(vec, offs):
    out = np.zeros_like(vec)
    for o in offs:
        out[o::nw] = vec[o::nw]
    return out


def cs_col(col):
    wDot = mask(wDotFull, blocks[col])
    return np.real(CFD.computeJacobianVectorProductFwd(wDot=wDot, residualDeriv=True, mode="CS", h=1e-40))


def mf_leak(resDot):
    return np.sqrt(comm.reduce(np.sum(mask(resDot, blocks["meanflow"]) ** 2)) or 0.0)


def trial(tag, reset):
    _ = cs_col("meanflow")   # inject residue
    if reset is not None:
        reset()
    leak = mf_leak(cs_col("gamma"))
    if rank == 0:
        print("  %-46s -> %.6e %s" % (tag, leak, "(CLEARED)" if leak < 1e-12 else ""))


def reset_setstates_real():
    CFD.setStates(np.real(CFD.getStates()))


def reset_setstates_w0():
    CFD.setStates(w0)


def reset_setstates_w0_getres():
    CFD.setStates(w0)
    CFD.getResidual(ap)


if rank == 0:
    print("\nreset between meanflow-CS warm-up and gamma-CS:")
trial("(none)", None)
trial("setStates(real(getStates()))", reset_setstates_real)
trial("setStates(w0_real)", reset_setstates_w0)
trial("setStates(w0_real) + getResidual", reset_setstates_w0_getres)
