"""Dev diagnostic: solve each functional's adjoint IN ISOLATION (fresh AP each time).

Context: diag_full_derivatives_bcm.py solves cl, cd, cmz, drag on the SAME AeroProblem, in
that order. skipAfterFailedAdjoint=True means once cmz's KSP diverges, adjointFailed=True is
latched on that AP and drag's solveAdjoint is skipped entirely (never attempted, just zeroed
by pyADflow) -- see pyADflow.py:4090-4096. So the "drag adjoint failed" conclusion in the
Verification_tuturial_mesh bundle is actually "drag was never tried", not "drag diverged".

This script builds a fresh solver/AP per functional and calls evalFunctionsSens with
evalFuncs=[<one func>] so each one gets its own solveAdjoint call with adjointFailed=False,
independent of the others' outcome.

usage:
  mpirun -np 12 --bind-to core python dev/diag_adjoint_isolated_bcm.py --variant hard
"""
import os
import sys
import argparse

import numpy
from mpi4py import MPI

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from diag_full_derivatives_bcm import buildSolver, EVAL_FUNCS  # noqa: E402

COMM = MPI.COMM_WORLD
RANK = COMM.rank


def rprint(*a, **kw):
    if RANK == 0:
        print(*a, **kw)
        sys.stdout.flush()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--variant", choices=["smooth", "hard", "sa"], required=True)
    p.add_argument("--funcs", default="cl,cd,cmz,drag", help="comma list, each solved in ISOLATION")
    p.add_argument("--restartfile", default=None)
    args = p.parse_args()

    funcs = args.funcs.split(",")

    for f in funcs:
        rprint("\n" + "=" * 100)
        rprint("==== isolated adjoint solve: variant=%s  func=%s  (fresh AP, no prior history) ===="
               % (args.variant, f))
        rprint("=" * 100)

        solver, ap, nShape = buildSolver(args.variant, complex_build=False, restartfile=args.restartfile)

        funcsSens = {}
        solver.evalFunctionsSens(ap, funcsSens, evalFuncs=[f])

        fail = {}
        solver.checkAdjointFailure(ap, fail)
        rprint("\n---- %s: adjoint fail flag = %s ----" % (f, fail.get("fail")))

        if RANK == 0:
            key = ap.name + "_" + f
            d = funcsSens[key]
            for dv in [k for k in d if k.endswith(ap.name)]:
                rprint("  d%s/d%-20s % .10e" % (f, dv, numpy.real(numpy.atleast_1d(d[dv]).flatten()[0])))
            rprint("  d%s/dspan[0]              % .10e" % (f, numpy.real(numpy.atleast_1d(d["span"]).flatten()[0])))
            sh = numpy.real(numpy.atleast_1d(d["shape"]).flatten())
            for j in (0, 71):
                if j < len(sh):
                    rprint("  d%s/dshape[%d]             % .10e" % (f, j, sh[j]))

        del solver


if __name__ == "__main__":
    main()
