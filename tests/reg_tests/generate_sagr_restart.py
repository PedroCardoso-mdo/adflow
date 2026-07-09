#!/usr/bin/env python
"""Generate the SA-noft2-Gamma-Retheta restart CGNS for the derivative
regression tests (test_*_sagr.py).

The upstream SA/Euler restart files (mdo_tutorial_*.cgns) are NOT generated
by any script in the repo — they are downloaded pre-made by
input_files/get-input-files.sh. This script is the SA-GR equivalent: it
takes the SA RANS tutorial-wing grid, converges it with the
SA-noft2-Gamma-Retheta model at a chosen Mach (P and T are kept at the SA
test values, so Re scales with Mach), and writes a grid+solution CGNS that
the tests use as ``restartfile``.

The written file contains all 8 states: the SA-GR restart variable set
includes Intermittency and ReThetat (outputMod.F90 solNames(itu2/itu3)),
and variableReading.F90 reads them back on restart.

Run on the REAL (not complex) build, e.g. from tests/reg_tests:

    python generate_sagr_restart.py --mach 0.25 --tu 0.003

then report the Mach/Tu used so reg_sagr.py can be pointed at the output
(the test AeroProblem must match the conditions converged here exactly).
"""

# built-ins
import argparse
import copy
import os

# MACH classes
from adflow import ADFLOW
from baseclasses import AeroProblem

from reg_default_options import adflowDefOpts
from reg_sagr import sagrBaseOptions

baseDir = os.path.dirname(os.path.abspath(__file__))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gridFile",
        default=os.path.join(baseDir, "../../input_files/mdo_tutorial_rans_scalar_jst.cgns"),
        #default=os.path.join(baseDir, "../../input_files/volumeMesh_NLF_L2.cgns"),
        #default=os.path.join(baseDir, "../../input_files/vmesh_L2.cgns"),
        help="grid to converge (default: the SA RANS tutorial wing)",
    )
    parser.add_argument(
        "--output",
        default=os.path.join(baseDir, "../../input_files/mdo_tutorial_sagr_scalar_jst.cgns"),
        help="grid+solution CGNS to write (becomes the tests' restartfile)",
    )
    parser.add_argument("--mach", type=float, default=0.12, help="Mach number (Re scales with it, P/T fixed)")
    parser.add_argument("--alpha", type=float, default=1.8, help="angle of attack [deg] (SA test value: 1.8)")
    parser.add_argument("--tu", type=float, default=0.015, help="freestream turbulence intensity (fraction)")
    parser.add_argument("--l2", type=float, default=1e-14, help="L2 convergence target")
    parser.add_argument("--ncycles", type=int, default=40000, help="max iterations")
    parser.add_argument(
        "--restartfile",
        default=None,
        help="restart from this CGNS instead of a cold start on --gridFile "
        "(e.g. the existing mdo_tutorial_sagr_scalar_jst.cgns, to continue "
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
    # options["transitioncrossflow"] = args.crossflow  # not available in this ADflow install

    options["gridfile"] = args.gridFile
    if args.restartfile is not None:
        # continue converging an existing solution (e.g. more iterations on
        # the current mdo_tutorial_sagr_scalar_jst.cgns)
        options["restartfile"] = args.restartfile
    else:
        # cold start on the tutorial grid; the smoke-converged SA-GR state is
        # what gets written out
        options.pop("restartfile", None)

    options["turbintensityinf"] = args.tu
    options["l2convergence"] = args.l2
    options["ncycles"] = args.ncycles
    # tied to whatever Mach this run actually uses (polar-script convention),
    # not a fixed constant
    options["acousticscalefactor"] = args.mach
    options["outputdirectory"] = os.path.dirname(os.path.abspath(args.output))
    # full precision so the tests linearize about exactly this state
    options["solutionprecision"] = "double"
    options["writevolumesolution"] = False
    options["writesurfacesolution"] = True
    options["surfacevariables"] = ["cp", "vx", "vy", "vz", "mach", "cf", "cfx", "cfy", "cfz"]

    if args.nk:
        options["useanksolver"] = False
        options["usenksolver"] = True

    # tutorial-wing AeroProblem (reg_aeroproblems.ap_tutorial_wing) with only
    # Mach/alpha overridden — these numbers must be mirrored in reg_sagr.py
    ap = AeroProblem(
        name="mdo_tutorial",
        alpha=args.alpha,
        mach=args.mach,
        P=20000.0,
        T=220.0,
        areaRef=45.5,
        chordRef=3.25,
        beta=0.0,
        R=287.87,
        xRef=0.0,
        yRef=0.0,
        zRef=0.0,
        evalFuncs=["cl", "cd"],
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
                "Now set in tests/reg_tests/reg_sagr.py: sagrGridFile/sagrRestartFile "
                "to this file, and the AeroProblem to mach=%g, alpha=%g, P=20000, "
                "T=220, chordRef=3.25, areaRef=45.5 with turbintensityinf=%g."
                % (args.mach, args.alpha, args.tu)
            )


if __name__ == "__main__":
    main()
