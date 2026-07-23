#!/usr/bin/env python
"""Compute CS for ALL 16 dR[row]/dw[col] blocks and compare each to the AD
ref (jacvecfwd_sagr_tut_wing.json). Reports every block's pass/fail under the
official tol (rtol=atol=5e-9), so we know exactly which per-block checks pass."""
import os, sys, copy, json

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mpi4py import MPI
import numpy as np
from adflow import ADFLOW_C
import reg_sagr
from reg_default_options import adflowDefOpts

comm = MPI.COMM_WORLD
rank = comm.rank
baseDir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ref = json.load(open(os.path.join(baseDir, "refs", "jacvecfwd_sagr_tut_wing.json")))
options = copy.copy(adflowDefOpts)
options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
options.update(copy.deepcopy(reg_sagr.sagrBaseOptions))
ap = copy.deepcopy(reg_sagr.ap_sagr_ar5_wing)
for dv in reg_sagr.sagrAeroDVs:
    ap.addDV(dv)
CFD = ADFLOW_C(options=options, debug=True)
CFD.getResidual(ap)

nw, blocks = reg_sagr.getStateBlocks(CFD)
wDotFull = CFD.getStatePerturbation(314)  # same seed as the official test


def mask(vec, offs):
    out = np.zeros_like(vec)
    for o in offs:
        out[o::nw] = vec[o::nw]
    return out


rtol = atol = 5e-9
rows = []
for col in blocks:
    wDot = mask(wDotFull, blocks[col])
    resDot = CFD.computeJacobianVectorProductFwd(wDot=wDot, residualDeriv=True, mode="CS", h=1e-40)
    for row in blocks:
        rd = mask(np.real(resDot), blocks[row])
        cs = np.sqrt(comm.reduce(np.sum(rd ** 2)) or 0.0)
        if rank == 0:
            key = "||dR[%s]/dw[%s] * eDot||" % (row, col)
            ad = ref[key]
            ok = abs(cs - ad) <= atol + rtol * abs(ad)
            rows.append((key, ad, cs, ok))

if rank == 0:
    print("\n%-42s %-16s %-16s %s" % ("block", "AD(ref)", "CS", "pass(5e-9)"))
    for key, ad, cs, ok in rows:
        print("%-42s %-16.6e %-16.6e %s" % (key, ad, cs, "PASS" if ok else "*** FAIL ***"))
    nf = sum(1 for *_, ok in rows if not ok)
    print("\n%d/16 blocks fail" % nf)
