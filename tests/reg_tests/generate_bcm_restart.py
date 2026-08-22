#!/usr/bin/env python
"""Generate the SA-BCM restart CGNS(s) for the derivative regression tests (test_*_bcm.py).

Cold-starts from the SAME tutorial-wing grid SA-GR uses (mdo_tutorial_rans_scalar_jst.cgns,
the standard ADflow test asset in input_files/) at the SAME flow conditions
(mach=0.15, alpha=1.8) as ap_sagr_tut_wing/ap_bcm_tut_wing, for apples-to-apples
comparability -- NOT the manuscript's NLF-0416 case (see reg_bcm.py module docstring).

Thin wrapper around dev/run_bcm_case.py's option set (kept in sync by hand -- see that
file's buildOptions()) plus a final CGNS restart write, since the case is expected NOT
to converge cleanly on the first attempt; use dev/run_bcm_case.py directly first to
iterate on options with raw solver output visible, then come back here once a variant
is known to run to actually produce the restart the registered suite depends on.

Writes ONE restart per SABCM_Exp variant (the two blends converge to different states):
  --variant smooth  (SABCM_Exp=False) -> input_files/mdo_tutorial_bcm_smooth_dp.cgns
  --variant hard    (SABCM_Exp=True)  -> input_files/mdo_tutorial_bcm_hard_dp.cgns

Run on the REAL (not complex) build, e.g. from tests/reg_tests:
    mpirun -np 2 --bind-to core python generate_bcm_restart.py --variant smooth
    mpirun -np 2 --bind-to core python generate_bcm_restart.py --variant hard
"""

import argparse
import copy
import os
import sys

from adflow import ADFLOW
from baseclasses import AeroProblem

from reg_default_options import adflowDefOpts

baseDir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(baseDir, "dev"))
from run_bcm_case import buildOptions  # noqa: E402  (canonical SA-BCM option dict, kept in sync there)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--variant", choices=["smooth", "hard"], required=True)
    parser.add_argument(
        "--gridFile",
        default=os.path.join(baseDir, "../../input_files/mdo_tutorial_rans_scalar_jst.cgns"),
        help="grid to converge (default: the standard ADflow tutorial-wing case, same as SA-GR's)",
    )
    parser.add_argument("--output", default=None, help="default: input_files/mdo_tutorial_bcm_<variant>_dp.cgns")
    parser.add_argument("--mach", type=float, default=0.15, help="matches ap_bcm_tut_wing")
    parser.add_argument("--alpha", type=float, default=1.8, help="matches ap_bcm_tut_wing")
    parser.add_argument(
        "--tu",
        type=float,
        default=0.5,
        help="freestream turbulence intensity Tu_inf (%%) -- feeds SABCM_TU. Not a fraction, and "
        "not the SA-GR-repo-only 'turbintensityinf' option (doesn't exist in this repo).",
    )
    parser.add_argument("--l2", type=float, default=1e-14)
    parser.add_argument("--ncycles", type=int, default=10000)
    parser.add_argument(
        "--restartfile",
        default=None,
        help="restart from this CGNS instead of a cold start on --gridFile (continue converging)",
    )
    parser.add_argument(
        "--nk",
        action="store_true",
        help="switch on the NK solver path only (usenksolver=True, useanksolver=False) "
        "instead of the ANK+NK ladder in buildOptions() -- useful once ANK has the case "
        "in the NK basin of attraction",
    )
    args = parser.parse_args()

    output = args.output or os.path.join(
        baseDir, "../../input_files/mdo_tutorial_bcm_%s_dp.cgns" % args.variant
    )

    outputDir = os.path.dirname(os.path.abspath(output))
    options = buildOptions(args.variant, outputDir)
    options["gridfile"] = args.gridFile
    if args.restartfile is not None:
        options["restartfile"] = args.restartfile
    else:
        # cold start on the tutorial-wing grid, same rationale as SA-GR's
        # generate_sagr_restart.py: M=0.15 survives a cold start through
        # ANK->NK, unlike the mesh's stock M=0.8 restart.
        options.pop("restartfile", None)

    options["sabcm_tu"] = args.tu
    options["l2convergence"] = args.l2
    options["ncycles"] = args.ncycles
    if args.nk:
        options["useanksolver"] = False
        options["usenksolver"] = True

    ap = AeroProblem(
        name="mdo_tutorial_bcm",
        alpha=args.alpha,
        mach=args.mach,
        P=20000.0,
        T=220.0,
        R=287.87,
        areaRef=45.5,
        chordRef=3.25,
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
            "WARNING: SA-BCM (%s) solve on %s did NOT converge cleanly (l2convergence not "
            "reached within ncycles). Writing the output anyway for inspection -- do NOT use "
            "it as a test restart until you get a clean convergence. Go back to "
            "dev/run_bcm_case.py --variant %s to iterate on options with raw solver output "
            "visible." % (args.variant, args.gridFile, args.variant)
        )

    CFDSolver.evalFunctions(ap, funcs)
    if CFDSolver.comm.rank == 0:
        print("Converged functions:", funcs)

    CFDSolver.writeVolumeSolutionFile(output, writeGrid=True)
    surfaceOutput = os.path.splitext(output)[0] + "_surf.cgns"
    CFDSolver.writeSurfaceSolutionFile(surfaceOutput)

    if CFDSolver.comm.rank == 0:
        print("Wrote restart:", output)
        print("Wrote surface (with cf):", surfaceOutput)
        if not funcs["fail"]:
            print(
                "Now confirm in tests/reg_tests/reg_bcm.py: bcmRestartFile%s points at this "
                "file, and ap_bcm_tut_wing matches mach=%g, alpha=%g, P=20000, T=220, "
                "chordRef=3.25, areaRef=45.5 with SABCM_TU=%g."
                % (args.variant.capitalize(), args.mach, args.alpha, args.tu)
            )


if __name__ == "__main__":
    main()
