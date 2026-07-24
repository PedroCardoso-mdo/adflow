#!/usr/bin/env python
"""Dev diagnostic: interactive SA-BCM run driver with RAW ADflow output visible.

This is the entry point for standing up SA-BCM's tutorial-wing/M=0.15 case for the
first time. Unlike the registered testflo suite, this does NOT capture stdout and does
NOT assume convergence -- it runs `solver(ap)` directly so every ANK/NK iteration line
streams to the terminal live, so you can watch it stall/diverge/converge and adjust
flags between attempts. Expect this NOT to converge cleanly on the first try.

Option dict below is adapted from a previously-used AR5-wing aeroOptions dict (see the
"AR5 REFERENCE (unused)" block at the bottom of this file for the original, kept for
provenance), retargeted at SA-GR's tutorial-wing mesh at mach=0.15 -- NOT the AR5 wing,
and NOT the AR5 script's liftIndex/areaRef/chordRef, which don't apply here. Differences
from the AR5 script, decided explicitly (see plan/session notes, not re-litigated here):
  - SABCM_Const2  = 0.02  (repo default), not the AR5 script's 2.0
  - useft2SA      = False for both SA-BCM and the plain-SA path
  - useBlockettes = True  (kept, matching AR5 script -- see reg_bcm.py module docstring
                            for why this is left on rather than defaulted off)

Two SA-BCM variants, both runnable from this one script:
  --variant smooth   SABCM_Exp=False (tanh blend, manuscript default)
  --variant hard      SABCM_Exp=True  (exp-sqrt blend, Mura & Cakmakcioglu original)

Cold-starts from the raw tutorial-wing grid (mdo_tutorial_rans_scalar_jst.cgns) -- the
same grid SA-GR cold-starts from -- rather than requiring a pre-converged restart. Once
a variant is known to run, generate_bcm_restart.py writes the restart the registered
suite (reg_bcm.py) actually uses; this script is the exploratory step before that.

examples:
  mpirun -np 2 --bind-to core /home/mdo/packages_v2/mach/bin/python dev/run_bcm_case.py --variant smooth
  mpirun -np 2 --bind-to core /home/mdo/packages_v2/mach/bin/python dev/run_bcm_case.py --variant hard
  # once a variant runs, continue from a partial restart instead of cold-starting again:
  mpirun -np 2 --bind-to core /home/mdo/packages_v2/mach/bin/python dev/run_bcm_case.py --variant smooth --restartfile /tmp/bcm_smooth_partial.cgns
"""

import argparse
import copy
import os
import sys

# dev/ is one level below reg_tests/; make the shared reg_* fixtures importable
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from adflow import ADFLOW
from baseclasses import AeroProblem

from reg_default_options import adflowDefOpts

baseDir = os.path.dirname(os.path.abspath(__file__))


def buildOptions(variant, outputDir):
    """The canonical SA-BCM option dict, adapted from the AR5 aeroOptions snippet.
    reg_bcm.py's bcmBaseOptionsSmooth/Hard mirror this exactly -- keep them in sync
    by hand if either changes."""
    gridFile = os.path.join(baseDir, "../../../input_files/mdo_tutorial_rans_scalar_jst.cgns")

    options = copy.copy(adflowDefOpts)
    options.update(
        {
            # ---- I/O ----
            "gridfile": gridFile,
            "outputdirectory": outputDir,
            "monitorvariables": ["resrho", "cl", "cd", "cdp", "cdv", "yplus", "resturb"],
            "surfacevariables": ["rho", "p", "vx", "vy", "vz", "cp", "mach", "cf", "cfx", "cfy", "cfz", "yplus"],
            "volumevariables": ["eddyratio"],
            "printiterations": True,
            "printtiming": True,
            # ---- physics ----
            "equationtype": "RANS",
            "turbulencemodel": "SA",
            "useft2sa": False,  # decided: off for both SA-BCM and plain-SA baseline
            "useblockettes": True,  # decided: kept True -- see reg_bcm.py module docstring
            "eddyvisinfratio": 1e-10,
            # ---- SA-BCM ----
            "use_sabcm": variant != "sa",  # variant="sa" -> plain SA baseline, SABCM off
            "sabcm_exp": variant == "hard",  # False=smooth (tanh), True=hard (exp-sqrt)
            "sabcm_const1": 0.002,
            "sabcm_const2": 0.02,  # decided: repo default, NOT the AR5 script's 2.0
            "sabcm_tu": 0.5,
            "sabcm_s0_tanh": 0.5,
            "sabcm_fsmooth": 0.08,
            "sabcm_maxsmooth": 50.0,
            "useapproxwalldistance": True,  # hard-required by inputParamRoutines.F90 when use_SABCM=True
            "frozenturbulence": False,
            # ---- startup smoother ----
            "smoother": "DADI",
            "mgcycle": "sg",
            "acousticscalefactor": 0.15,  # matches SA-GR's tutorial-wing/M=0.15 case
            "infchangecorrection": True,
            "resaveraging": "alternate",
            # ---- ANK ----
            "useanksolver": True,
            "ankcfl0": 1.0,
            "ankchartimesteptype": "VLR",
            "anksecondordswitchtol": 1e-3,
            "ankadpc": True,
            "ankcfllimit": 3.0e5,
            # ---- NK ----
            "usenksolver": True,
            "nkadpc": True,
            # nkswitchtol left at the ADflow default (the AR5 snippet's vars.nkswitchtol
            # was an unresolved variable -- not hardcoded here on purpose)
            # ---- termination ----
            "l2convergence": 1e-14,
            "ncycles": 10000,
            "solutionprecision": "double",
            "writevolumesolution": True,
            "writesurfacesolution": True,
        }
    )
    return options


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--variant",
        choices=["smooth", "hard", "sa"],
        required=True,
        help="smooth=SABCM_Exp=False, hard=SABCM_Exp=True, sa=use_SABCM off (plain-SA baseline, "
        "everything else identical) for apples-to-apples comparison",
    )
    parser.add_argument(
        "--restartfile",
        default=None,
        help="restart from this CGNS instead of a cold start on the tutorial-wing grid",
    )
    parser.add_argument("--mach", type=float, default=0.15, help="matches SA-GR's tutorial-wing case")
    parser.add_argument("--alpha", type=float, default=1.8, help="matches SA-GR's tutorial-wing case")
    parser.add_argument(
        "--tu",
        type=float,
        default=0.5,
        help="freestream turbulence intensity Tu_inf (%%) -- feeds SABCM_TU (Re_theta_c "
        "correlation). Unlike SA-GR's 'turbintensityinf' (a fraction, not a valid option in "
        "this repo), SA-BCM's only freestream-Tu knob is SABCM_TU, in percent.",
    )
    parser.add_argument("--ncycles", type=int, default=10000)
    parser.add_argument("--output", default=None, help="write a grid+solution CGNS here when done (optional)")
    args = parser.parse_args()

    outputDir = os.path.join(baseDir, "../output_files")
    options = buildOptions(args.variant, outputDir)
    options["ncycles"] = args.ncycles
    options["sabcm_tu"] = args.tu
    if args.restartfile is not None:
        options["restartfile"] = args.restartfile
    else:
        options.pop("restartfile", None)  # cold start, like SA-GR's generate_sagr_restart.py default

    ap = AeroProblem(
        name="mdo_tutorial_bcm_%s" % args.variant,
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

    print("==== SA-BCM dev run: variant=%s (SABCM_Exp=%s) ====" % (args.variant, options["sabcm_exp"]))
    CFDSolver = ADFLOW(options=options, debug=True)
    CFDSolver(ap)  # raw ADflow stdout streams here, uncaptured

    funcs = {}
    CFDSolver.checkSolutionFailure(ap, funcs)
    if funcs["fail"]:
        print("==== SA-BCM (%s) did NOT converge cleanly within ncycles=%d ====" % (args.variant, args.ncycles))
    else:
        print("==== SA-BCM (%s) converged ====" % args.variant)

    CFDSolver.evalFunctions(ap, funcs)
    print("Functions:", funcs)

    if args.output is not None:
        CFDSolver.writeVolumeSolutionFile(args.output, writeGrid=True)
        print("Wrote:", args.output)


if __name__ == "__main__":
    main()


# --------------------------------------------------------------------------
# AR5 REFERENCE (unused) -- the original AR5-wing aeroOptions dict this script's
# option set was adapted from, kept verbatim for provenance/comparison. NOT wired
# into buildOptions() above -- liftIndex/areaRef/chordRef/mesh are AR5-wing-specific
# and don't apply to the tutorial-wing case this script actually runs.
# --------------------------------------------------------------------------
_AR5_REFERENCE_AERO_OPTIONS = {
    # I/O Parameters
    "gridFile": None,  # vars.gridFile -- AR5 wing mesh
    "outputDirectory": None,  # vars.output
    "monitorvariables": ["resrho", "cl", "cd", "cdp", "cdv", "yplus", "resturb"],
    "surfaceVariables": ["rho", "p", "vx", "vy", "vz", "cp", "mach", "cf", "cfx", "cfy", "cfz", "yplus", "blank"],
    "volumeVariables": ["eddyratio"],
    # Physics Parameters
    "equationType": None,  # vars.equationType
    "useft2SA": False,
    "useBlockettes": True,
    "ANKCFL0": 1.0,
    "ANKCharTimeStepType": "VLR",
    "eddyVisInfRatio": 1e-10,
    # Transition Model Parameters
    "use_SABCM": True,
    "SABCM_Const2": 2.0,  # differs from repo default (0.02) -- NOT reused, see module docstring
    "SABCM_TU": 0.5,
    # Solver Parameters
    "smoother": "DADI",
    "MGCycle": "sg",
    "infchangecorrection": True,
    "acousticScaleFactor": None,  # vars.mach
    # ANK Solver Parameters
    "useANKSolver": True,
    "ANKSecondOrdSwitchTol": 1e-3,
    "ANKADPC": True,
    "ANKCFLLimit": 3.0e5,
    # NK Solver Parameters
    "useNKSolver": True,
    "nkswitchtol": None,  # vars.nkswitchtol -- unresolved in the original snippet
    "NKADPC": True,
    # Termination Criteria
    "L2Convergence": None,  # vars.L2Convergence
    "nCycles": None,  # vars.nCycles
    # Adjoint Parameters
    "setMonitor": True,
    "printIterations": True,
    "printTiming": True,
    "adjointMonitorStep": 10,
    "verifyState": True,
    "verifySpatial": True,
    "ADPC": True,
    "skipAfterFailedAdjoint": True,
    "adjointMaxIter": None,  # vars.adjointMaxIter
    "adjointL2Convergence": None,  # vars.adjointL2Convergence
    "liftIndex": None,  # AR5-wing-specific, not reused (tutorial wing uses pyADflow default)
}
