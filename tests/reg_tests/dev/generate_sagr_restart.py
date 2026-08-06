#!/usr/bin/env python
"""Generate the SA-noft2-Gamma-Retheta restart CGNS for the derivative
regression tests (test_*_sagr.py).

2026-07-19: switched from the SA RANS tutorial-wing grid to the AR5
plain-wing case, so derivatives were validated about a real NK-converged
state on a real 3D wing rather than the small ANK-only tutorial-wing case.

2026-07-23: switched BACK to the tutorial-wing grid
(`mdo_tutorial_rans_scalar_jst.cgns`, the standard ADflow test asset already
in `input_files/` -- no separate copy needed, unlike the AR5 files). Reason:
the AR5 case never got past a chronic quasi-stall (`Step` pinned ~0.01,
`docs/current-task.md`'s "Step 3 (AR5) CS check fails hard") no matter the
iteration budget. The SAME tutorial-wing mesh, run with the SA-GR model at
Mach=0.15 (down from the mesh's stock 0.8 -- see run_sagr.py's Reynolds-
number note in
`.../3D_Plain_Wing/pedro_test/run_sagr.py`) instead of AR5's Mach=0.2,
converges cleanly with no stall. `--mach`/`--alpha`/the AeroProblem geometry
below now match that tutorial-wing case, not AR5's.

The written file contains all 8 states: the SA-GR restart variable set
includes Intermittency and ReThetat (outputMod.F90 solNames(itu2/itu3)),
and variableReading.F90 reads them back on restart.

Run on the REAL (not complex) build, e.g. from tests/reg_tests:

    python generate_sagr_restart.py

then report the Mach/Tu used so reg_sagr.py can be pointed at the output
(the test AeroProblem must match the conditions converged here exactly) --
already done for the --mach/--tu/--alpha defaults below, which match
ap_sagr_tut_wing in reg_sagr.py.
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
        default=os.path.join(baseDir, "../../../input_files/mdo_tutorial_rans_scalar_jst.cgns"),
        help="grid to converge (default: the standard ADflow tutorial-wing case)",
    )
    parser.add_argument(
        "--output",
        default=os.path.join(baseDir, "../../../input_files/mdo_tutorial_sagr_dp.cgns"),
        help="grid+solution CGNS to write (becomes the tests' restartfile)",
    )
    parser.add_argument("--mach", type=float, default=0.15, help="Mach number (matches ap_sagr_tut_wing)")
    parser.add_argument("--alpha", type=float, default=1.8, help="angle of attack [deg] (tutorial-wing test value)")
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
        default=10000,
        help="max iterations -- this tutorial-wing/M=0.15 case converges cleanly "
        "(no AR5-style stall), so this is a generous cap, not a stall workaround. "
        "Confirmed 2026-07-23: naturally hit L2Convergence (didn't need the cap) "
        "at outer iter 83 / Iter Tot 3676, scaledTotalRes=1.09e-9, fail=False.",
    )
    parser.add_argument(
        "--restartfile",
        default=None,
        help="restart from this CGNS instead of a cold start on --gridFile "
        "(e.g. the existing mdo_tutorial_sagr_dp.cgns, to continue "
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
        help="enable transitionCrossflow (D_scf). Off by default here (matches "
        "sagrBaseOptions as of 2026-07-23): D_scf is identically zero on this "
        "~2D tutorial wing regardless, so this flag is mostly moot on this "
        "mesh -- kept for parity with the AR5-mesh workflow, where it mattered "
        "(the helicity term divides by max(velMag, 1e-10) TWICE, Eqs. 23-26, "
        "and is unvalidated on a real 3D wing per CLAUDE.md rule 3).",
    )
    args = parser.parse_args()

    options = copy.copy(adflowDefOpts)
    options.update(copy.deepcopy(sagrBaseOptions))
    options["transitioncrossflow"] = args.crossflow

    options["gridfile"] = args.gridFile
    if args.restartfile is not None:
        # continue converging an existing solution (e.g. more iterations on
        # the current mdo_tutorial_sagr_dp.cgns)
        options["restartfile"] = args.restartfile
    else:
        # cold start on the tutorial-wing grid (freestream); M=0.15 is far
        # less transonic than the mesh's stock M=0.8 restart, so this
        # survives a cold start through the ANK->CANK->NK ladder (unlike a
        # cold M=0.8 start, which NaNs -- see pedro_test/run_fresh.log)
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

    # tutorial-wing AeroProblem, matching ap_sagr_tut_wing in reg_sagr.py —
    # these numbers must stay mirrored there
    ap = AeroProblem(
        name="mdo_tutorial_sagr",
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
                "this file, and ap_sagr_tut_wing matches mach=%g, alpha=%g, P=20000, T=220, "
                "chordRef=3.25, areaRef=45.5 with turbintensityinf=%g."
                % (args.mach, args.alpha, args.tu)
            )


if __name__ == "__main__":
    main()
