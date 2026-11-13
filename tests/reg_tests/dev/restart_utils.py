"""Common part of the generate_<model>_restart.py scripts."""

# built-ins
import argparse
import copy
import os
import sys

# MACH classes
from adflow import ADFLOW

# the fixtures live one level up, in tests/reg_tests
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from reg_default_options import adflowDefOpts  # noqa: E402
import reg_transition  # noqa: E402

cases = {
    "sagr": (reg_transition.sagrOptions, reg_transition.ap_sagr_tut_wing, reg_transition.sagrRestartFile),
    "sgamma": (reg_transition.sgammaOptions, reg_transition.ap_sgamma_tut_wing, reg_transition.sgammaRestartFile),
    # states for the total-derivative tests of a model option (--variant)
    "sagr_crossflow": (
        reg_transition.sagrCrossflowAdjointOptions,
        reg_transition.ap_sagr_tut_wing,
        reg_transition.sagrCrossflowRestartFile,
    ),
    "sgamma_artvisc": (
        reg_transition.sgammaArtViscAdjointOptions,
        reg_transition.ap_sgamma_tut_wing,
        reg_transition.sgammaArtViscRestartFile,
    ),
}


def generateRestart(model, description=None):
    # --variant selects the state of a model option: sagr crossflow, sgamma artvisc
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--variant", default=None)
    variant = pre.parse_known_args()[0].variant
    if variant is not None:
        model = model + "_" + variant
    modelOptions, aeroProblem, restartFile = cases[model]

    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--variant", default=None, help="model option: crossflow (sagr) or artvisc (sgamma)")
    parser.add_argument("--gridFile", default=reg_transition.tutWingGridFile, help="grid to converge")
    parser.add_argument("--output", default=restartFile, help="volume solution (with grid) to write")
    parser.add_argument(
        "--restartFile",
        default=None,
        help="continue from this solution instead of starting from freestream",
    )
    parser.add_argument("--nCycles", type=int, default=modelOptions["ncycles"], help="iteration limit")
    args = parser.parse_args()

    options = copy.copy(adflowDefOpts)
    options.update(copy.deepcopy(modelOptions))
    options["gridfile"] = args.gridFile
    if args.restartFile is not None:
        options["restartfile"] = args.restartFile
    else:
        options.pop("restartfile", None)
    options["ncycles"] = args.nCycles
    options["outputdirectory"] = os.path.dirname(os.path.abspath(args.output))

    ap = copy.deepcopy(aeroProblem)

    CFDSolver = ADFLOW(options=options)
    CFDSolver(ap)

    funcs = {}
    CFDSolver.checkSolutionFailure(ap, funcs)
    failed = funcs["fail"]
    CFDSolver.evalFunctions(ap, funcs)
    if CFDSolver.comm.rank == 0:
        print("Functions:", funcs)

    CFDSolver.writeVolumeSolutionFile(args.output, writeGrid=True)
    surfaceOutput = os.path.splitext(args.output)[0] + "_surf.cgns"
    CFDSolver.writeSurfaceSolutionFile(surfaceOutput)

    if CFDSolver.comm.rank == 0:
        print("Wrote restart:", args.output)
        print("Wrote surface solution:", surfaceOutput)
        if failed:
            print(
                "WARNING: the solve did not reach l2convergence within nCycles. "
                "Do not use this file as the restart state of the tests."
            )
