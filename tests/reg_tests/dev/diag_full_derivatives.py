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

  --mode cs        COMPLEX build (ADFLOW_C). For each selected DV -- aero, span,
                   selected twist components, and selected shape components --
                   perturb by 1e-40j, re-converge (flow-solver output visible),
                   print imag(f)/h, and tabulate adjoint-vs-CS abs/rel diff
                   against --ref (default: the trained adjoint ref json, which
                   already stores every shape/twist component).

Each complex re-converge is a full nonlinear solve (~5-10 min), so the DV set
is configurable; mind the cost when sweeping many shape components.

examples:
  # real build: dump the adjoint totals for all DVs (writes json for `cs` to read)
  mpirun -np 2 --bind-to core /home/mdo/packages_v2/mach/bin/python \
      dev/diag_full_derivatives.py --mode adjoint --out /tmp/adj_all.json

  # complex build: CS-vs-adjoint for aero + span + twist 0,3 + shape 0,24,48,71
  mpirun -np 2 --bind-to core /home/mdo/packages_v2/mach/bin/python \
      dev/diag_full_derivatives.py --mode cs --ref /tmp/adj_all.json \
      --twist 0,3 --shape 0,24,48,71
"""
import os
import sys
import copy
import json
import argparse

import numpy
from mpi4py import MPI

# tests/reg_tests on the path (this file lives in tests/reg_tests/dev/)
RT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, RT)

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


def buildSolver(complex_build):
    """Mirror test_adjoint_sagr.py setUp. Returns (solver, ap, nShape)."""
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


# --------------------------------------------------------------------- modes
def run_adjoint(args):
    solver, ap, nShape = buildSolver(complex_build=False)
    rprint("\n==== ADJOINT MODE (real build) — nShape(local)=%d ====" % nShape)
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
            for j in args.shape:
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
    solver, ap, nShape = buildSolver(complex_build=True)
    ref = loadRef(args.ref, ap.name)
    rprint("\n==== CS MODE (complex build) — ref=%s  nShape(local)=%d ====" % (args.ref, nShape))

    results = []  # (label, func, cs, ref)

    def do_aero(dv):
        rprint("\n#### re-converging complex solver: aero DV %s ####" % dv)
        setattr(ap, dv, getattr(ap, dv) + h * 1j)
        solver.resetFlow(ap)
        solver(ap, writeSolution=False)
        funcs = {}
        solver.evalFunctions(ap, funcs)
        setattr(ap, dv, getattr(ap, dv) - h * 1j)
        for f in EVAL_FUNCS:
            cs = numpy.imag(funcs[ap.name + "_" + f]) / h
            r = refLookup(ref, ap.name, f, dv + "_" + ap.name)
            results.append(("d%s/d%s" % (f, dv), cmpRow("d%s/d%s" % (f, dv), cs, r)))

    def do_geom(dvname, idx, arrlen):
        rprint("\n#### re-converging complex solver: geom DV %s[%d] ####" % (dvname, idx))
        xRef = {"twist": numpy.zeros(NTWIST, dtype="D"), "span": numpy.zeros(1, dtype="D"),
                "shape": numpy.zeros(nShape, dtype="D")}
        xRef[dvname][idx] += h * 1j
        solver.resetFlow(ap)
        solver.DVGeo.setDesignVars(xRef)
        solver(ap, writeSolution=False)
        funcs = {}
        solver.evalFunctions(ap, funcs)
        for f in EVAL_FUNCS:
            cs = numpy.imag(funcs[ap.name + "_" + f]) / h
            r = refLookup(ref, ap.name, f, dvname, idx=idx)
            lbl = "d%s/d%s[%d]" % (f, dvname, idx)
            results.append((lbl, cmpRow(lbl, cs, r)))

    if not args.skip_aero:
        for dv in ["alpha", "mach", "P", "T"]:
            do_aero(dv)
    if not args.skip_geom:
        do_geom("span", 0, 1)
        for i in args.twist:
            do_geom("twist", i, NTWIST)
        for j in args.shape:
            do_geom("shape", j, nShape)

    rprint("\n================= SUMMARY: adjoint vs complex-step =================")
    for _, row in results:
        rprint(row)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--mode", choices=["adjoint", "cs"], required=True)
    p.add_argument("--ref", default=os.path.join(RT, "refs", "adjoint_sagr_tut_wing.json"),
                   help="adjoint totals json to compare CS against (cs mode)")
    p.add_argument("--out", default=os.path.join(RT, "refs", "adjoint_sagr_tut_wing.json"),
                   help="where adjoint mode writes the totals")
    p.add_argument("--shape", default="0,24,48,71",
                   help="comma-separated local shape component indices to sweep")
    p.add_argument("--twist", default="0,3", help="comma-separated twist component indices")
    p.add_argument("--skip-aero", action="store_true")
    p.add_argument("--skip-geom", action="store_true")
    args = p.parse_args()
    args.shape = [int(x) for x in args.shape.split(",") if x != ""]
    args.twist = [int(x) for x in args.twist.split(",") if x != ""]

    if args.mode == "adjoint":
        run_adjoint(args)
    else:
        run_cs(args)


if __name__ == "__main__":
    main()
