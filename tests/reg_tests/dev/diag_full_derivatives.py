"""Dev diagnostic: SA-GR full total-derivative check with ADflow output VISIBLE.

Unlike testflo (which swallows per-test stdout), this prints the ADflow flow-
solver and adjoint-KSP convergence live, and runs the full df/dx check over
ALL design variables the adjoint/CS tests use -- aero (alpha, mach, P, T) and
geometric (span, twist, and MULTIPLE individual local `shape` components) --
not just the reduced set the regression test exercises.

It mirrors test_adjoint_sagr.py's setup exactly (same mesh, FFD, options,
tutorial-wing twist/span/shape DVGeo), including the complex-build solver
overrides (ANK/NKADPC off + no coupled ANK: the AD preconditioner is not in
the complexify build -- see test_adjoint_sagr.py).

Two modes, matched to the installed build:

  --mode adjoint   REAL build (ADFLOW). Solve the adjoint for every func and
                   print the full df/dx incl. EVERY twist/span/shape component,
                   with the adjoint KSP trace visible. Writes the totals to
                   --out (default refs/adjoint_sagr_tut_wing.json is the trained
                   ref, so by default `cs` mode compares against the same
                   values the regression test does).

  --mode cs        COMPLEX build (ADFLOW_C). For each selected DV -- aero
                   (alpha, mach, P, T), span, twist, and shape (ALL FFD points
                   by default) -- reset the flow to a clean restart, perturb by
                   1e-40j, re-converge (flow-solver output visible), print
                   imag(f)/h, and at the end print a TABLE tabulating
                   adjoint-vs-CS absolute / relative / % error against --ref
                   (default: the trained adjoint ref json).

Each complex re-converge is a full nonlinear solve, and the primal is forced to
grind all `--ncycles` iters (default 5000, `--l2 1e-20` disables early exit) to
push the FD-PC complex solve past its ~1e-8 stall. `resetFlow` runs before every
DV so there is no cross-DV contamination. Sweeping ALL shape points is the
default and is expensive -- narrow it with `--shape a,b,c` when iterating.

examples:
  # real build: dump the adjoint totals for all DVs (writes json for `cs` to read)
  mpirun -np 2 --bind-to core /home/mdo/packages_v2/mach/bin/python \
      dev/diag_full_derivatives.py --mode adjoint --out /tmp/adj_all.json

  # complex build: CS-vs-adjoint over aero + span + ALL twist + ALL shape,
  #                deep 5000-iter re-converge, final abs/rel/% error table
  mpirun -np 2 --bind-to core /home/mdo/packages_v2/mach/bin/python \
      dev/diag_full_derivatives.py --mode cs --ref /tmp/adj_all.json

  # cheaper: just a few shape points, shallower budget
  mpirun -np 2 --bind-to core /home/mdo/packages_v2/mach/bin/python \
      dev/diag_full_derivatives.py --mode cs --shape 0,24,48,71 --ncycles 2000
"""
import os
import sys
import copy
import json
import argparse

import numpy
from mpi4py import MPI

# Resolve the config/ref root so this same file runs in TWO layouts:
#   - repo:       tests/reg_tests/dev/diag_full_derivatives.py  (RT = tests/reg_tests)
#   - standalone: <folder>/diag_full_derivatives.py beside reg_sagr.py + inputs/
# (a self-contained copy, e.g. the Verification_tuturial_mesh archive that must
#  run on another machine). Detection: if reg_sagr.py sits next to this file, we
#  are standalone; otherwise assume the repo dev/ layout.
_HERE = os.path.dirname(os.path.abspath(__file__))
if os.path.isfile(os.path.join(_HERE, "reg_sagr.py")):
    RT = _HERE  # standalone folder
else:
    RT = os.path.dirname(_HERE)  # repo: dev/ -> tests/reg_tests
sys.path.insert(0, RT)


def _default_ref():
    """Trained adjoint ref, wherever this layout keeps it (repo refs/ or the
    standalone inputs/)."""
    for c in (os.path.join(RT, "refs", "adjoint_sagr_tut_wing.json"),
              os.path.join(RT, "inputs", "adjoint_sagr_tut_wing.json")):
        if os.path.isfile(c):
            return c
    return os.path.join(RT, "refs", "adjoint_sagr_tut_wing.json")

from pygeo import DVGeometry
from pyspline import Curve
from reg_default_options import adflowDefOpts, IDWarpDefOpts
from reg_sagr import ap_sagr_tut_wing, sagrBaseOptions, sagrAeroDVs, sagrFFDFile

COMM = MPI.COMM_WORLD
RANK = COMM.rank
EVAL_FUNCS = ["cl", "cd", "cmz", "drag"]
NTWIST = 6


def rprint(*a):
    if RANK == 0:
        print(*a, flush=True)


def getDVGeo(ffdFile, isComplex):
    """Identical twist/span/shape ref-axis DVGeo to test_adjoint_sagr.py."""
    DVGeo = DVGeometry(ffdFile, isComplex=isComplex)
    DVGeo.addRefAxis(
        "wing",
        Curve(
            x=numpy.linspace(5.0 / 4.0, 1.5 / 4.0 + 7.5, NTWIST),
            y=numpy.zeros(NTWIST),
            z=numpy.linspace(0, 14, NTWIST),
            k=2,
        ),
    )

    def twist(val, geo):
        for i in range(NTWIST):
            geo.rot_z["wing"].coef[i] = val[i]

    def span(val, geo):
        C = geo.extractCoef("wing")
        s = geo.extractS("wing")
        for i in range(len(C)):
            C[i, 2] += s[i] * val[0]
        geo.restoreCoef(C, "wing")

    DVGeo.addGlobalDV("twist", [0] * NTWIST, twist, lower=-10, upper=10, scale=1.0)
    DVGeo.addGlobalDV("span", [0], span, lower=-10, upper=10, scale=1.0)
    DVGeo.addLocalDV("shape", lower=-0.5, upper=0.5, axis="y", scale=10.0)
    return DVGeo


def buildSolver(complex_build, ncycles=None, l2=None):
    """Mirror test_adjoint_sagr.py setUp. Returns (solver, ap, nShape).

    ncycles / l2 (when given) override the primal iteration budget and the
    L2Convergence target. For the complex re-converge we set l2 absurdly tight
    (1e-20) so the run NEVER early-exits on L2Convergence and instead grinds all
    `ncycles` iterations -- the point being to push the FD-PC complex solve as
    deep as it can go past its ~1e-8 stall (see complex-build-no-ad-pc caveat).
    """
    if complex_build:
        from adflow import ADFLOW_C as ADF
        from idwarp import USMesh_C as MESH
    else:
        from adflow import ADFLOW as ADF
        from idwarp import USMesh as MESH

    options = copy.copy(adflowDefOpts)
    options["outputdirectory"] = os.path.join(RT, options["outputdirectory"])
    options.update(copy.deepcopy(sagrBaseOptions))
    # make ADflow talk
    options["printiterations"] = True

    # deeper primal: run all iters (no early L2 exit) to push the solve past
    # the FD-PC stall -- applied to both modes so adjoint linearizes about the
    # same deep state the CS solve reaches.
    if ncycles is not None:
        options["ncycles"] = int(ncycles)
    if l2 is not None:
        options["l2convergence"] = float(l2)

    if complex_build:
        # complex build has no AD routines -> AD preconditioner unavailable;
        # coupled ANK warns it may diverge on the stiff sources. FD-PC decoupled.
        options["ankadpc"] = False
        options["nkadpc"] = False
        options["ankcoupledswitchtol"] = 1e-16

    mesh_options = copy.copy(IDWarpDefOpts)
    mesh_options.update({"gridFile": options["gridfile"]})

    ap = copy.deepcopy(ap_sagr_tut_wing)
    ap.evalFuncs = EVAL_FUNCS
    for dv in sagrAeroDVs:
        ap.addDV(dv)

    solver = ADF(options=options, debug=True)
    solver.setMesh(MESH(options=mesh_options))
    solver.setDVGeo(getDVGeo(sagrFFDFile, isComplex=complex_build), pointSetKwargs={"embTol": 1e-12, "eps": 1e-14})
    solver.getResidual(ap)
    nShape = solver.DVGeo.getNDV() - NTWIST - 1  # total - twist(6) - span(1)
    return solver, ap, nShape


# ---------------------------------------------------------------- ref helpers
def _asarray(node):
    """Decode a BaseRegTest-serialized value (float / list / __ndarray__ dict)."""
    if isinstance(node, dict):
        if "__ndarray__" in node:
            return numpy.array(node["__ndarray__"]).flatten()
        if "value" in node:
            return _asarray(node["value"])
    return numpy.atleast_1d(numpy.array(node, dtype=float)).flatten()


def loadRef(reffile, apname):
    """{func: {dv_or_'aeroname': ndarray}} from the adjoint ref json, or None."""
    if not reffile or not os.path.isfile(reffile):
        return None
    efs = json.load(open(reffile)).get("Eval Functions Sens:", {})
    return efs  # indexed as efs[apname+'_'+f][dvkey]


def refLookup(ref, apname, f, dvkey, idx=0):
    if ref is None:
        return None
    try:
        return _asarray(ref[apname + "_" + f][dvkey])[idx]
    except Exception:
        return None


def cmpRow(label, cs, ref):
    if ref is None:
        return "%-34s  CS=% .10e   (no ref)" % (label, cs)
    ad = abs(cs - ref)
    rd = ad / max(abs(ref), 1e-30)
    return "%-34s  adj=% .10e  CS=% .10e  |Δ|=%.2e  rel=%.2e" % (label, ref, cs, ad, rd)


def resolve_indices(spec, n):
    """Parse a --shape/--twist spec into a list of 0-based indices.
    "all" -> every index [0, n); otherwise comma-separated ints (clamped to n)."""
    spec = str(spec).strip().lower()
    if spec in ("", "all"):
        return list(range(n))
    return [int(x) for x in spec.split(",") if x != "" and 0 <= int(x) < n]


def print_table(rows):
    """rows: list of (label, cs, ref, l2). Prints an aligned adjoint-vs-CS table
    with absolute, relative, and % error plus the L2 norm the complex primal
    re-converged to for that DV (totalRfinal/totalR0), then worst-case lines."""
    hdr = "%-24s %18s %18s %11s %11s %11s %11s" % (
        "d(func)/d(dv)", "adjoint", "complex-step", "|abs err|", "rel err", "% err", "L2 reached",
    )
    rprint("\n" + "=" * len(hdr))
    rprint("  SUMMARY TABLE: adjoint vs complex-step  (L2 reached = CS re-converge totalRfinal/totalR0)")
    rprint("=" * len(hdr))
    rprint(hdr)
    rprint("-" * len(hdr))
    worst_rel = (-1.0, "")
    n_noref = 0
    for label, cs, ref, l2 in rows:
        l2s = "%.2e" % l2 if l2 is not None else "-"
        if ref is None:
            rprint("%-24s %18.10e %18s %11s %11s %11s %11s" % (label, cs, "(no ref)", "-", "-", "-", l2s))
            n_noref += 1
            continue
        ad = abs(cs - ref)
        rd = ad / max(abs(ref), 1e-30)
        rprint("%-24s %18.10e %18.10e %11.2e %11.2e %10.4f%% %11s" % (label, ref, cs, ad, rd, rd * 100.0, l2s))
        if rd > worst_rel[0]:
            worst_rel = (rd, label)
    rprint("-" * len(hdr))
    if worst_rel[0] >= 0:
        rprint("worst relative error: %.3e (%.4f%%) at %s" % (worst_rel[0], worst_rel[0] * 100.0, worst_rel[1]))
    if n_noref:
        rprint("(%d row(s) had no reference value)" % n_noref)


# --------------------------------------------------------------------- modes
def run_adjoint(args):
    solver, ap, nShape = buildSolver(complex_build=False, ncycles=args.ncycles, l2=args.l2)
    shape_idx = resolve_indices(args.shape, nShape)
    rprint("\n==== ADJOINT MODE (real build) — nShape(local)=%d, printing %d shape comp ====" % (nShape, len(shape_idx)))
    funcsSens = {}
    solver.evalFunctionsSens(ap, funcsSens)  # adjoint KSP trace prints here
    fail = {}
    solver.checkAdjointFailure(ap, fail)
    rprint("\n==== adjoint fail flag: %s ====" % fail.get("fail"))

    if RANK == 0:
        for f in EVAL_FUNCS:
            key = ap.name + "_" + f
            rprint("\n--- df/d(x) for %s ---" % key)
            d = funcsSens[key]
            for dv in [k for k in d if k.endswith(ap.name)]:  # aero
                rprint("  %-28s % .10e" % (dv, numpy.real(numpy.atleast_1d(d[dv]).flatten()[0])))
            rprint("  span                         % .10e" % numpy.real(numpy.atleast_1d(d["span"]).flatten()[0]))
            tw = numpy.real(numpy.atleast_1d(d["twist"]).flatten())
            for i in range(len(tw)):
                rprint("  twist[%d]                     % .10e" % (i, tw[i]))
            sh = numpy.real(numpy.atleast_1d(d["shape"]).flatten())
            for j in shape_idx:
                if j < len(sh):
                    rprint("  shape[%d]%s% .10e" % (j, " " * (21 - len(str(j))), sh[j]))
        # dump for cs mode to read
        out = {"Eval Functions Sens:": {}}
        for f in EVAL_FUNCS:
            key = ap.name + "_" + f
            out["Eval Functions Sens:"][key] = {}
            for dv, v in funcsSens[key].items():
                arr = numpy.real(numpy.atleast_1d(numpy.array(v)).flatten())
                out["Eval Functions Sens:"][key][dv] = {
                    "__ndarray__": [arr.tolist()], "dtype": "float64", "shape": [1, arr.size]
                }
        json.dump(out, open(args.out, "w"))
        rprint("\nwrote adjoint totals -> %s" % args.out)


def run_cs(args):
    h = 1e-40
    solver, ap, nShape = buildSolver(complex_build=True, ncycles=args.ncycles, l2=args.l2)
    ref = loadRef(args.ref, ap.name)
    shape_idx = resolve_indices(args.shape, nShape)
    twist_idx = resolve_indices(args.twist, NTWIST)
    rprint("\n==== CS MODE (complex build) — ref=%s  nShape(local)=%d ====" % (args.ref, nShape))
    rprint("     primal budget: ncycles=%s  l2convergence=%s (force all iters)" % (args.ncycles, args.l2))
    rprint("     sweeping %d shape + %d twist components" % (len(shape_idx), len(twist_idx)))

    rows = []  # (label, cs, ref, l2) for the final table

    # --- warm start ---------------------------------------------------------
    # Default: resetFlow (freestream) before each DV -> the perturbed complex
    # solve must re-converge from scratch, and on stiff directions (mach) the
    # FD-PC stalls the REAL residual at ~4e-5, capping the derivative at ~1e-3.
    # --warm instead deep-converges the UNPERTURBED base once, snapshots its
    # (essentially real) state, and re-seats that state before applying each
    # 1e-40j perturbation. NK then only has to propagate the infinitesimal
    # imaginary part from an already-deep real state, rather than drag the real
    # residual down from freestream -- the test of whether mach's ~1e-3 is an
    # FD-PC-from-freestream artifact or a true floor.
    base_w = None
    if args.warm:
        rprint("\n#### --warm: deep-converging the UNPERTURBED base once ####")
        solver.resetFlow(ap)
        solver(ap, writeSolution=False)
        base_w = solver.getStates().copy()
        r0, _, rf = solver.getResNorms()
        rprint("  [base L2 reached: %.3e]" % (numpy.real(rf) / max(numpy.real(r0), 1e-300)))

    def seat_base():
        """Start point for the next perturbed solve: warm (deep base state) or
        cold (freestream)."""
        if args.warm and base_w is not None:
            solver.setStates(base_w.copy())
        else:
            solver.resetFlow(ap)

    def l2reached():
        """L2 norm the primal actually re-converged to: totalRfinal/totalR0."""
        r0, _, rf = solver.getResNorms()
        r0 = numpy.real(r0)
        return numpy.real(rf) / r0 if r0 > 0 else None

    def do_aero(dv):
        rprint("\n#### re-converging complex solver: aero DV %s ####" % dv)
        setattr(ap, dv, getattr(ap, dv) + h * 1j)
        seat_base()  # warm (deep base) or cold (freestream) start
        solver(ap, writeSolution=False)
        l2 = l2reached()
        funcs = {}
        solver.evalFunctions(ap, funcs)
        setattr(ap, dv, getattr(ap, dv) - h * 1j)
        rprint("  [L2 reached: %.3e]" % l2 if l2 is not None else "  [L2 reached: n/a]")
        for f in EVAL_FUNCS:
            cs = numpy.imag(funcs[ap.name + "_" + f]) / h
            r = refLookup(ref, ap.name, f, dv + "_" + ap.name)
            lbl = "d%s/d%s" % (f, dv)
            rprint("  " + cmpRow(lbl, cs, r))
            rows.append((lbl, cs, r, l2))

    def do_geom(dvname, idx, arrlen):
        rprint("\n#### re-converging complex solver: geom DV %s[%d] ####" % (dvname, idx))
        xRef = {"twist": numpy.zeros(NTWIST, dtype="D"), "span": numpy.zeros(1, dtype="D"),
                "shape": numpy.zeros(nShape, dtype="D")}
        xRef[dvname][idx] += h * 1j
        seat_base()  # warm (deep base) or cold (freestream) start
        solver.DVGeo.setDesignVars(xRef)
        solver(ap, writeSolution=False)
        l2 = l2reached()
        funcs = {}
        solver.evalFunctions(ap, funcs)
        rprint("  [L2 reached: %.3e]" % l2 if l2 is not None else "  [L2 reached: n/a]")
        for f in EVAL_FUNCS:
            cs = numpy.imag(funcs[ap.name + "_" + f]) / h
            r = refLookup(ref, ap.name, f, dvname, idx=idx)
            lbl = "d%s/d%s[%d]" % (f, dvname, idx)
            rprint("  " + cmpRow(lbl, cs, r))
            rows.append((lbl, cs, r, l2))

    if not args.skip_aero:
        aero_dvs = [d.strip() for d in args.aero.split(",") if d.strip()] if args.aero else ["alpha", "mach", "P", "T"]
        for dv in aero_dvs:
            do_aero(dv)
    if not args.skip_geom:
        if not args.skip_span:
            do_geom("span", 0, 1)
        if not args.skip_twist:
            for i in twist_idx:
                do_geom("twist", i, NTWIST)
        for j in shape_idx:
            do_geom("shape", j, nShape)

    print_table(rows)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--mode", choices=["adjoint", "cs"], required=True)
    p.add_argument("--ref", default=_default_ref(),
                   help="adjoint totals json to compare CS against (cs mode)")
    p.add_argument("--out", default=_default_ref(),
                   help="where adjoint mode writes the totals")
    p.add_argument("--shape", default="all",
                   help="'all' (every FFD shape point, default) or comma-separated indices")
    p.add_argument("--twist", default="all",
                   help="'all' (default) or comma-separated twist component indices")
    p.add_argument("--ncycles", type=int, default=5000,
                   help="primal iteration budget for the (complex) re-converge (default 5000)")
    p.add_argument("--l2", type=float, default=1e-20,
                   help="L2Convergence target; 1e-20 forces all ncycles iters (default 1e-20)")
    p.add_argument("--aero", default="", help="comma list of aero DVs to sweep (default: alpha,mach,P,T)")
    p.add_argument("--warm", action="store_true",
                   help="warm-start: deep-converge the unperturbed base once, then re-seat that "
                        "state before each perturbation (NK only propagates the imaginary part) "
                        "instead of resetFlow-to-freestream. Tests the mach FD-PC-floor hypothesis.")
    p.add_argument("--skip-aero", action="store_true")
    p.add_argument("--skip-geom", action="store_true")
    p.add_argument("--skip-span", action="store_true", help="skip the span DV (geom mode)")
    p.add_argument("--skip-twist", action="store_true", help="skip all twist DVs (geom mode)")
    args = p.parse_args()

    if args.mode == "adjoint":
        run_adjoint(args)
    else:
        run_cs(args)


if __name__ == "__main__":
    main()
