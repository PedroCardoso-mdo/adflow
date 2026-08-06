"""Find the MINIMUM complex-solve depth (iters + L2) needed for the OFFICIAL
adjoint-vs-CS test points to pass their inherited tolerances -- WITHOUT touching
the tolerances.

Motivation: TestCmplxStepSAGR (test_adjoint_sagr.py) fails cmplx_test_aero_dvs /
cmplx_test_geom_dvs not because the adjoint is wrong but because the complex
build's FD-PC re-converge was too shallow (~1e-8) at the default budget. The
complex build CAN go deeper -- it just needs more iterations (10k runs reach
~1e-10). This tool measures how deep is *enough*.

Method (staged / warm-continue): for each official DV point, resetFlow once,
apply the 1e-40j perturbation once, then call the solver in CHUNK-iter bursts.
Because l2convergence is pinned to 1e-20 the solve never early-exits, so each
burst adds exactly CHUNK iters and continues from the previous state. After every
burst we evaluate the CS derivative for all funcs and test it against the trained
adjoint ref at the OFFICIAL tolerance. We report the first (cumulative iters,
scaled L2) at which every func passes.

Official points + tolerances (from test_adjoint_sagr.py):
  aero: alpha, mach          rtol=1e-8, atol=5e-10
  geom: span[0], twist[0], shape[0]   rtol=5e-9, atol=5e-9

Run (complex build), e.g.:
  mpirun -np 2 --bind-to core /home/mdo/packages_v2/mach/bin/python \
      dev/diag_min_iters.py --chunk 500 --maxit 12000
"""
import os
import sys
import argparse

import numpy
from mpi4py import MPI

_HERE = os.path.dirname(os.path.abspath(__file__))
RT = _HERE if os.path.isfile(os.path.join(_HERE, "reg_sagr.py")) else os.path.dirname(_HERE)
sys.path.insert(0, RT)

# reuse the exact build/ref helpers the full-derivative dev tool uses
from diag_full_derivatives import buildSolver, loadRef, refLookup, EVAL_FUNCS, NTWIST, rprint

COMM = MPI.COMM_WORLD


def _default_ref():
    for c in (os.path.join(RT, "refs", "adjoint_sagr_tut_wing.json"),
              os.path.join(RT, "inputs", "adjoint_sagr_tut_wing.json")):
        if os.path.isfile(c):
            return c
    return os.path.join(RT, "refs", "adjoint_sagr_tut_wing.json")


# (name, kind, idx, rtol, atol) -- the official test points
POINTS = [
    ("alpha", "aero", 0, 1e-8, 5e-10),
    ("mach", "aero", 0, 1e-8, 5e-10),
    ("span", "geom", 0, 5e-9, 5e-9),
    ("twist", "geom", 0, 5e-9, 5e-9),
    ("shape", "geom", 0, 5e-9, 5e-9),
]


def passes(cs, ref, rtol, atol):
    """numpy.allclose element semantics for a single scalar."""
    if ref is None:
        return False
    return abs(cs - ref) <= atol + rtol * abs(ref)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--chunk", type=int, default=500, help="iters per burst")
    p.add_argument("--maxit", type=int, default=12000, help="give-up cap per point")
    p.add_argument("--ref", default=_default_ref())
    args = p.parse_args()
    h = 1e-40

    solver, ap, nShape = buildSolver(complex_build=True)
    ref = loadRef(args.ref, ap.name)
    # never early-exit; burst size controls granularity
    solver.setOption("l2convergence", 1e-20)
    solver.setOption("ncycles", args.chunk)

    rprint("\n==== MIN-ITERS TO PASS OFFICIAL TOLERANCES (chunk=%d, maxit=%d) ====" % (args.chunk, args.maxit))
    rprint("     ref=%s  nShape=%d\n" % (args.ref, nShape))

    summary = []
    for name, kind, idx, rtol, atol in POINTS:
        rprint("\n#### point: %s%s  (tol rtol=%.0e atol=%.0e) ####" % (
            name, "" if kind == "aero" else "[%d]" % idx, rtol, atol))

        # reference values for the 4 funcs at this point
        refs = {}
        for f in EVAL_FUNCS:
            if kind == "aero":
                refs[f] = refLookup(ref, ap.name, f, name + "_" + ap.name)
            else:
                refs[f] = refLookup(ref, ap.name, f, name, idx=idx)

        # reset once, perturb once
        solver.resetFlow(ap)
        if kind == "aero":
            setattr(ap, name, getattr(ap, name) + h * 1j)
        else:
            xRef = {"twist": numpy.zeros(NTWIST, dtype="D"), "span": numpy.zeros(1, dtype="D"),
                    "shape": numpy.zeros(nShape, dtype="D")}
            xRef[name][idx] += h * 1j
            solver.DVGeo.setDesignVars(xRef)

        r0_init = None
        cum = 0
        hit = None
        while cum < args.maxit:
            solver(ap, writeSolution=False)  # continues from current state
            cum += args.chunk
            r0, _, rf = solver.getResNorms()
            if r0_init is None:
                r0_init = numpy.real(r0)
            absres = numpy.real(rf)
            scaled = absres / r0_init if r0_init > 0 else absres

            funcs = {}
            solver.evalFunctions(ap, funcs)
            cs = {f: numpy.imag(funcs[ap.name + "_" + f]) / h for f in EVAL_FUNCS}
            ok = all(passes(cs[f], refs[f], rtol, atol) for f in EVAL_FUNCS)
            worst = max((abs(cs[f] - refs[f]) / max(abs(refs[f]), 1e-30)) for f in EVAL_FUNCS if refs[f] is not None)
            rprint("  iters=%5d  L2=%.2e  absRes=%.2e  worstRel=%.2e  %s" % (
                cum, scaled, absres, worst, "PASS" if ok else "fail"))
            if ok:
                hit = (cum, scaled, absres, worst)
                break

        # unperturb aero (geom xRef is local; reset next loop)
        if kind == "aero":
            setattr(ap, name, getattr(ap, name) - h * 1j)

        if hit:
            rprint("  -> PASSES official tol at iters=%d, L2=%.2e (worstRel=%.2e)" % (hit[0], hit[1], hit[3]))
            summary.append((name if kind == "aero" else "%s[%d]" % (name, idx), hit[0], hit[1], hit[3], "PASS"))
        else:
            rprint("  -> did NOT pass within %d iters (best worstRel above)" % args.maxit)
            summary.append((name if kind == "aero" else "%s[%d]" % (name, idx), args.maxit, scaled, worst, "FAIL"))

    rprint("\n================= MIN-ITERS SUMMARY =================")
    rprint("%-12s %8s %10s %11s %6s" % ("point", "min-iters", "L2", "worstRel", "res"))
    for nm, it, l2, wr, st in summary:
        rprint("%-12s %8d %10.2e %11.2e %6s" % (nm, it, l2, wr, st))
    if summary:
        worst_it = max(summary, key=lambda x: x[1] if x[4] == "PASS" else 1e18)
        rprint("\n=> to pass ALL official points: ncycles >= %d, and L2 ~ %.1e"
               % (max(s[1] for s in summary if s[4] == "PASS") if any(s[4] == "PASS" for s in summary) else args.maxit,
                  max(s[2] for s in summary if s[4] == "PASS") if any(s[4] == "PASS" for s in summary) else 0.0))


if __name__ == "__main__":
    main()
