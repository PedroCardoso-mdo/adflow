#!/usr/bin/env python
"""Generate the SA-noft2-Gamma-Retheta restart CGNS for the derivative
regression tests (test_*_sagr.py).

2026-07-19: switched from the SA RANS tutorial-wing grid to the AR5
plain-wing case (`ar5_plain_wing_vol_L3.cgns`, copied into input_files/,
gitignored like the rest of that directory) and the validated
ANK->CANK->NK production ladder (STRATEGY.md /
`.../3D_Plain_Wing/best_strategie/run_strategy.py`), so derivatives are
validated about a real NK-converged state on a real 3D wing rather than the
small ANK-only tutorial-wing case.

The written file contains all 8 states: the SA-GR restart variable set
includes Intermittency and ReThetat (outputMod.F90 solNames(itu2/itu3)),
and variableReading.F90 reads them back on restart.

Run on the REAL (not complex) build, e.g. from tests/reg_tests:

    python generate_sagr_restart.py

then report the Mach/Tu used so reg_sagr.py can be pointed at the output
(the test AeroProblem must match the conditions converged here exactly) --
already done for the --mach/--tu/--alpha defaults below, which match
ap_sagr_ar5_wing in reg_sagr.py.
"""

# built-ins
import argparse
import copy
import os

# MACH classes
from adflow import ADFLOW
from baseclasses import AeroProblem

# dev/ is one level below reg_tests/; make the shared reg_* fixtures
# and sibling campaign modules importable when run directly (mpirun python).
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from reg_default_options import adflowDefOpts
from reg_sagr import sagrBaseOptions

baseDir = os.path.dirname(os.path.abspath(__file__))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gridFile",
        default=os.path.join(baseDir, "../../input_files/ar5_plain_wing_vol_L3.cgns"),
        help="grid to converge (default: the AR5 plain-wing case)",
    )
    parser.add_argument(
        "--output",
        default=os.path.join(baseDir, "../../input_files/ar5_plain_wing_sagr_dp.cgns"),
        help="grid+solution CGNS to write (becomes the tests' restartfile)",
    )
    parser.add_argument("--mach", type=float, default=0.2, help="Mach number (matches ap_sagr_ar5_wing)")
    parser.add_argument("--alpha", type=float, default=0.0, help="angle of attack [deg] (AR5 test value: 0.0)")
    parser.add_argument("--tu", type=float, default=0.0025, help="freestream turbulence intensity (fraction)")
    parser.add_argument(
        "--l2",
        type=float,
        default=1e-14,
        help="L2 convergence target -- kept tight/effectively disabled. L2Convergence is checked "
        "every iteration regardless of solver phase, and is LOOSER than ANKCoupledSwitchTol/"
        "nkswitchtol (1e-5) at any value above 1e-5, so a loose target (e.g. 1e-3) stops the run "
        "during plain ANK, before CANK/NK ever engage. Convergence depth is controlled by "
        "--ncycles instead (see its help).",
    )
    parser.add_argument(
        "--ncycles",
        type=int,
        default=3594,
        help="max iterations -- caps the run to exactly reproduce "
        "full_ladder_1e-5_production/run_full.log at outer iter 282 / Iter Tot 3594 "
        "(totalRes=9.9223314015988036E-04), the specific point requested for this restart",
    )
    parser.add_argument(
        "--restartfile",
        default=None,
        help="restart from this CGNS instead of a cold start on --gridFile "
        "(e.g. the existing ar5_plain_wing_sagr_dp.cgns, to continue "
        "converging it further). --gridFile is still used as the grid.",
    )
    parser.add_argument(
        "--nk",
        action="store_true",
        help="switch on the NK solver path (usenksolver=True, useanksolver=False) "
        "instead of the ANK-only path in sagrBaseOptions. Useful for driving a "
        "restart the rest of the way to convergence once ANK has it in the NK "
        "basin of attraction.",
    )
    parser.add_argument(
        "--crossflow",
        action="store_true",
        help="enable transitionCrossflow (D_scf). Off by default here: the "
        "helicity term divides by max(velMag, 1e-10) TWICE (paper Eqs. 23-26, "
        "faithfully implemented), which blows up in near-wall/near-stagnation "
        "3D cells where velMag is genuinely tiny -- unvalidated per CLAUDE.md "
        "rule 3 and e5cd58cd's own 'validation pending' note. sagrBaseOptions "
        "keeps it True only for adjoint-block test coverage on the ~2D "
        "tutorial mesh, where D_scf is identically zero anyway.",
    )
    args = parser.parse_args()

    options = copy.copy(adflowDefOpts)
    options.update(copy.deepcopy(sagrBaseOptions))
    options["transitioncrossflow"] = args.crossflow

    options["gridfile"] = args.gridFile
    if args.restartfile is not None:
        # continue converging an existing solution (e.g. more iterations on
        # the current ar5_plain_wing_sagr_dp.cgns)
        options["restartfile"] = args.restartfile
    else:
        # cold start on the AR5 grid (freestream); matches the validated
        # production run (full_ladder_1e-5_production) -- no restart
        options.pop("restartfile", None)

    options["turbintensityinf"] = args.tu
    options["l2convergence"] = args.l2
    options["ncycles"] = args.ncycles
    options["outputdirectory"] = os.path.dirname(os.path.abspath(args.output))
    # solutionprecision/writevolumesolution/writesurfacesolution/
    # volumevariables/surfacevariables all already set in sagrBaseOptions
    # (double precision + full transitionDebug field set for validation)

    if args.nk:
        options["useanksolver"] = False
        options["usenksolver"] = True

    # AR5 plain-wing AeroProblem, matching ap_sagr_ar5_wing in reg_sagr.py —
    # these numbers must stay mirrored there
    ap = AeroProblem(
        name="ar5_plain_wing",
        alpha=args.alpha,
        mach=args.mach,
        # explicit P/T (standard atmosphere sea level, numerically identical
        # to altitude=0.0) -- see reg_sagr.py's ap_sagr_ar5_wing comment
        P=101325.0,
        T=288.15,
        areaRef=0.1,
        chordRef=1.0,
        beta=0.0,
        xRef=0.0,
        yRef=0.0,
        zRef=0.0,
        evalFuncs=["cl", "cd", "cmz", "drag"],
    )

    CFDSolver = ADFLOW(options=options)
    CFDSolver(ap)

    funcs = {}
    CFDSolver.checkSolutionFailure(ap, funcs)
    if funcs["fail"] and CFDSolver.comm.rank == 0:
        print(
            "WARNING: SA-GR solve on %s did not converge cleanly (l2convergence "
            "not reached within ncycles). Writing the output anyway for "
            "inspection -- do NOT use it as a test restart until you get a "
            "clean convergence." % args.gridFile
        )

    CFDSolver.evalFunctions(ap, funcs)
    if CFDSolver.comm.rank == 0:
        print("Converged functions:", funcs)

    CFDSolver.writeVolumeSolutionFile(args.output, writeGrid=True)

    surfaceOutput = os.path.splitext(args.output)[0] + "_surf.cgns"
    CFDSolver.writeSurfaceSolutionFile(surfaceOutput)

    if CFDSolver.comm.rank == 0:
        print("Wrote restart:", args.output)
        print("Wrote surface (with cf):", surfaceOutput)
        if not funcs["fail"]:
            print(
                "Now confirm in tests/reg_tests/reg_sagr.py: sagrRestartFile points at "
                "this file, and ap_sagr_ar5_wing matches mach=%g, alpha=%g, altitude=0, "
                "chordRef=1.0, areaRef=0.1 with turbintensityinf=%g."
                % (args.mach, args.alpha, args.tu)
            )


if __name__ == "__main__":
    main()
