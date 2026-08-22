#!/usr/bin/env python
"""Dev diagnostic: SA-BCM block vs blockette residual-operator consistency, RAW output.

ADflow builds the RANS residual through two interchangeable paths inside blocketteRes
(src/NKSolver/blockette.F90), selected by `useBlockettes`:

    useBlockettes = False -> blockResCore     -> saSource (reference, src/turbulence/sa.F90,
                                                  Tapenade-differentiated)
    useBlockettes = True  -> blocketteResCore -> saSource (hand-maintained duplicate,
                                                  src/NKSolver/blockette.F90:1004-1250)

Unlike SA-GR (whose blockette copy WAS found to have drifted, and is now force-disabled by
pyADflow), NOTHING force-disables useBlockettes for SA-BCM -- and reg_bcm.py's
bcmBaseOptions deliberately keeps useBlockettes=True (matching the previously-used AR5
script). That means production SA-BCM runs, right now, use the hand-synced blockette.F90
copy by default, and it has never been checked against sa.F90. This script is that check,
with the per-variable diffs printed to stdout (not swallowed by testflo) so a mismatch is
immediately legible. test_blockette_bcm.py is the registered/asserted version of the same
check, added after this script's output has been inspected once.

Run for BOTH SA-BCM variants (SABCM_Exp False/True) -- the KS-smoothmax blend differs
between them and could plausibly drift differently in the hand-maintained copy.

examples:
  mpirun -np 2 --bind-to core /home/mdo/packages_v2/mach/bin/python dev/diag_blockette_bcm.py --variant smooth
  mpirun -np 2 --bind-to core /home/mdo/packages_v2/mach/bin/python dev/diag_blockette_bcm.py --variant hard
"""

import argparse
import copy
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from adflow import ADFLOW

from reg_default_options import adflowDefOpts
from reg_bcm import ap_bcm_tut_wing, bcmBaseOptionsSmooth, bcmBaseOptionsHard, bcmAeroDVs

baseDir = os.path.dirname(os.path.abspath(__file__))

VAR_NAMES = ["rho", "rhou", "rhov", "rhow", "rhoE", "nuTilde"]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--variant", choices=["smooth", "hard"], required=True)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--perturbation", type=float, default=1.0e-2, help="relative state perturbation")
    args = parser.parse_args()

    variantOptions = bcmBaseOptionsSmooth if args.variant == "smooth" else bcmBaseOptionsHard

    options = copy.copy(adflowDefOpts)
    options["outputdirectory"] = os.path.join(baseDir, "../output_files")
    options.update(copy.deepcopy(variantOptions))
    # exercise a state that isn't the (near-zero-residual) converged restart, so every
    # source/diffusion term is active -- mirrors test_blockette_sagr.py's rationale
    options["ncycles"] = 1
    options["usenksolver"] = False
    options["useanksolver"] = False

    CFDSolver = ADFLOW(options=copy.deepcopy(options), debug=True)

    ap = copy.deepcopy(ap_bcm_tut_wing)
    for dv in bcmAeroDVs:
        ap.addDV(dv)

    CFDSolver.getResidual(ap)  # propagate the restart state

    nw = int(CFDSolver.adflow.flowvarrefstate.nw)
    if nw != 6:
        raise RuntimeError("expected nw=6 for SA-BCM (plain-SA state layout); got nw=%d" % nw)

    state0 = CFDSolver.getStates()
    rng = np.random.default_rng(args.seed)
    state = state0 * (1.0 + args.perturbation * rng.standard_normal(state0.shape))

    def residual(useBlockettes):
        CFDSolver.setOption("useBlockettes", useBlockettes)
        assert CFDSolver.getOption("useBlockettes") == useBlockettes
        CFDSolver.setStates(state)
        return CFDSolver.getResidual(ap).copy()

    print("==== block vs blockette, SA-BCM variant=%s (SABCM_Exp=%s) ====" % (args.variant, variantOptions["sabcm_exp"]))
    res_block = residual(False)
    res_blockette = residual(True)

    rb = res_block.reshape(-1, nw)
    rk = res_blockette.reshape(-1, nw)

    print("%-10s %14s %14s %14s" % ("var", "max|block|", "max|blockette|", "max|diff|"))
    worst = (-1.0, None)
    for l in range(nw):
        mb = np.max(np.abs(rb[:, l])) if rb.shape[0] else 0.0
        mk = np.max(np.abs(rk[:, l])) if rk.shape[0] else 0.0
        diff = np.max(np.abs(rk[:, l] - rb[:, l])) if rb.shape[0] else 0.0
        relDiff = diff / max(mb, 1e-300)
        print("%-10s %14.6e %14.6e %14.6e  (rel=%.3e)" % (VAR_NAMES[l], mb, mk, diff, relDiff))
        if diff > worst[0]:
            worst = (diff, VAR_NAMES[l])

    print("\nworst absolute diff: %.6e on %s" % worst)
    if worst[0] > 1e-7 * max(np.max(np.abs(rb)), 1e-300):
        print(
            "!! LARGE MISMATCH -- blockette.F90's saSource has drifted from sa.F90 for this "
            "variant, same failure mode SA-GR hit before being force-disabled. Do NOT trust "
            "useBlockettes=True runs of this variant until this is fixed."
        )
    else:
        print("block and blockette agree to floating-point roundoff for this variant.")


if __name__ == "__main__":
    main()
