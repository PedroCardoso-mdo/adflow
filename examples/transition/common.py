"""Plumbing shared by the transition examples: paths, command line, solve and report.

Nothing here decides physics or solver settings: those come from recipes.py and from the
case script that calls these helpers.
"""
import argparse
import ast
import hashlib
import json
import math
import os
import time

import numpy as np

import recipes

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
# Meshes are not in git. Default location: <repo>/showcase/meshes ; override with TRANSITION_MESH_ROOT.
MESH_ROOT = os.environ.get("TRANSITION_MESH_ROOT", os.path.join(REPO, "showcase", "meshes"))
FFD_DIR = os.path.join(HERE, "ffd")


def meshPath(*parts):
    return os.path.join(MESH_ROOT, *parts)


def baseParser(doc, models=recipes.MODELS, defaultModel="gr"):
    """Arguments every example understands."""
    p = argparse.ArgumentParser(description=doc, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--model", choices=models, default=defaultModel,
                   help="gr = SA-noft2-Gamma-Retheta, sgamma = SA-noft2-Gamma, bcm / bcmo = SA-BCM smooth / original, sa = plain SA")
    p.add_argument("--mesh", default=None, help="CGNS volume mesh (default: the validated mesh under showcase/meshes)")
    p.add_argument("--outputDir", default=None, help="output directory (default: ./output_<case>)")
    p.add_argument("--nCycles", type=int, default=None,
                   help="cap on 'Iter Tot' (default: the recipe value). 4 = two solver iterations, for a smoke test")
    p.add_argument("--nuinf", type=float, default=None, help="eddyVisInfRatio (default: recipes.NUINF of the model)")
    p.add_argument("--recipe", default=None, choices=sorted(recipes.RECIPES), help="solver recipe (default: the validated one)")
    p.add_argument("--opt", action="append", default=[], metavar="KEY=VAL", help="override one ADflow option (repeatable)")
    p.add_argument("--writeVolume", action="store_true", help="also write the volume solution")
    return p


def applyOverrides(opts, args):
    """--nCycles, --writeVolume and --opt KEY=VAL on top of the recipe (keys are case-insensitive)."""
    over = {}
    for kv in args.opt:
        k, v = kv.split("=", 1)
        try:
            v = ast.literal_eval(v)
        except (ValueError, SyntaxError):
            pass
        over[k] = v
    if args.nCycles is not None:
        over["nCycles"] = args.nCycles
    if args.writeVolume:
        over["writeVolumeSolution"] = True
    lower = {k.lower(): k for k in opts}
    for k, v in over.items():
        opts[lower.get(k.lower(), k)] = v
    return over


def md5(path, n=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(n), b""):
            h.update(chunk)
    return h.hexdigest()


def provenance(comm, gridFile):
    """Which ADflow ran (git stamp written at install time), on how many ranks, on which mesh."""
    import adflow

    info = {"adflowFile": adflow.__file__, "nProcs": comm.size, "gridFile": os.path.abspath(gridFile)}
    try:
        from adflow import gitInfo

        info["adflowGit"] = {k: getattr(gitInfo, k) for k in dir(gitInfo) if not k.startswith("_")}
    except ImportError:
        info["adflowGit"] = None
    if comm.rank == 0:
        info["gridMd5"] = md5(gridFile)
    return info


def makeSolver(opts, comm):
    from adflow import ADFLOW

    if comm.rank == 0:
        os.makedirs(opts["outputDirectory"], exist_ok=True)
    comm.barrier()
    return ADFLOW(options=opts, comm=comm)


def solve(CFD, ap, comm, label=None, coldStart=False):
    """Run one solve and return the numbers: cl, cd, cdp, cdv, residual drop, wall time."""
    if coldStart:
        CFD.resetFlow(ap)          # free-stream start, as every validated run
    t0 = time.time()
    CFD(ap)
    wall = time.time() - t0
    funcs = {}
    CFD.evalFunctions(ap, funcs, evalFuncs=["cl", "cd", "cdp", "cdv"])
    CFD.checkSolutionFailure(ap, funcs)
    r0, _, rFinal = (float(v) for v in CFD.getResNorms())
    rel = rFinal / r0 if r0 > 0 else float("nan")
    target = CFD.getOption("L2Convergence")
    rec = dict(
        label=label or ap.name,
        alpha=float(np.real(ap.alpha)),
        cl=float(np.real(funcs[ap.name + "_cl"])),
        cd=float(np.real(funcs[ap.name + "_cd"])),
        cdp=float(np.real(funcs[ap.name + "_cdp"])),
        cdv=float(np.real(funcs[ap.name + "_cdv"])),
        totalR0=r0,                      # residual of the free-stream state (iteration 0 of the log)
        totalRFinal=rFinal,
        resRel=rel,
        orders=-math.log10(rel) if rel > 0 else float("nan"),
        converged=bool(rel <= 1.2 * target and not funcs.get("fail", False)),
        fail=bool(funcs.get("fail", False)),
        wallTimeSeconds=wall,
    )
    if comm.rank == 0:
        print(
            f"RESULT {rec['label']}: alpha={rec['alpha']:.4f} cl={rec['cl']:.6f} cd={rec['cd'] * 1e4:.3f} counts "
            f"(cdp {rec['cdp'] * 1e4:.3f}, cdv {rec['cdv'] * 1e4:.3f})  res/r0={rel:.2e} ({rec['orders']:.2f} orders, "
            f"target {target:g}) converged={rec['converged']}  totalR0={r0:.16e}  wall={wall:.0f}s",
            flush=True,
        )
    return rec


def writeResult(comm, opts, case, records, prov, overrides=None):
    """<outputDir>/result.json: the numbers and everything needed to reproduce them."""
    if comm.rank != 0:
        return

    def conv(o):
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        if isinstance(o, np.ndarray):
            return o.tolist()
        return str(o)

    path = os.path.join(opts["outputDirectory"], "result.json")
    with open(path, "w") as f:
        json.dump(dict(case=case, results=records, overrides=overrides or {}, options=opts, provenance=prov),
                  f, indent=1, default=conv)
    print(f"wrote {path}", flush=True)


def banner(comm, text):
    if comm.rank == 0:
        print("=" * 100 + "\n" + text + "\n" + "=" * 100, flush=True)
