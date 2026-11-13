#!/usr/bin/env python
"""2D shape optimisation of the NACA 0012 with a transition model: single-point and multipoint.

    minimise   sum_i w_i cd_i
    w.r.t.     20 FFD shape variables (10 x 2 x 2 box, vertical motion, +-0.05) + one alpha per point
    subject to cl_i = cl*_i,  volume >= initial,  0.5 <= t/t0 <= 3 at 100 stations,
               leading / trailing edge fixed, both span planes move together (2D)
    optimiser  SLSQP

Trimmed from the campaign-24 driver (L/24_minmodel_2d/common/{run_opt.py,problem.py}); every
design point is a cold start from free stream.

Mesh    showcase/meshes/naca0012/mesh_L0/volumeMesh.cgns (193 224 cells); FFD ffd/naca0012_ffd_10x2x2.xyz
Flow    M 0.1252, sea level (Re 2.92e6); a point is "Tu[%]:cl", e.g. 0.3:0.5 (article point A)

Reference numbers (cd in counts, own model, L0; showcase/_manifest_2d.md 24-E/F/G)
    Single point A = Tu 0.3 %, cl 0.5, from the NACA 0012:
        scalar   GR 65.1 -> 31.6 (-51.5 %) | SA-gamma 68.4 -> 32.3 | SA-BCM (original) 65.5 -> 31.3
        matrix   GR 66.0 -> 30.7 (-53.6 %) | SA-gamma 72.0 -> 46.5 | SA-BCM (original) 68.3 -> 48.8
        Status: GR final; SA-gamma and SA-BCM provisional (matrix SA-BCM still running, job 2218879).
        Stopped by the stall rule (< 0.1 count over 10 major iterations) or the 48 h limit,
        never by the SLSQP criterion.
    Old point Tu 0.5 %, cl 0.3: GR 62.19 -> 31.24 (-49.8 %), final.
    Multipoint 0.5:0.3,1.0:0.3,0.5:0.4 (job 1980228, 3 x 42 ranks, GR, scalar): 91 major iterations,
        per point -41.5 / -40.7 / -43.4 %, cl constraints met to 2e-5. Final.
    Multipoint 0.3:0.5,0.6:0.5,0.3:0.7 in matrix (mp7): running; the GR design points stalled at
        1-2e-5 until ANKSecondOrdSwitchTol was raised to 1e-4 (recipe gr_2d_matrix_opt, jobs 2219449 /
        2219697: cd 47.34 counts at Tu 0.3 %, cl 0.5 on that shape). Status: open.
    Before trusting an optimum of SA-gamma or SA-BCM, re-analyse it with GR: the gain often does
    not survive (--model gr --startShapeTxt xstar.shape.txt --startAlpha <a> --baselineOnly evaluates
    a shape at the start alpha, without trim).

Cluster cost: 64 ranks per single-point optimisation, GR 0.42 h per major iteration, 90 % of the
    gain after ~24 h; multipoint 32-42 ranks per point, 48 h.

Commands
    mpirun -np 64 python opt_airfoil.py --model gr --points 0.3:0.5 --outputDir c2_gr
    mpirun -np 64 python opt_airfoil.py --model gr --points 0.3:0.5 --dissipation matrix --outputDir c2_gr_mtx
    mpirun -np 126 python opt_airfoil.py --model gr --points 0.3:0.5,0.6:0.5,0.3:0.7 --dissipation matrix
    mpirun -np 6 --bind-to core python opt_airfoil.py --model gr --points 0.3:0.5 --baselineOnly --nCycles 4   # smoke
"""
import os

import numpy as np
from mpi4py import MPI

import common
import recipes

# geometry and optimisation set-up (L/24_minmodel_2d/config/setup.py:118-132)
MACH, ALTITUDE, AREA_REF, CHORD_REF = 0.125191453, 0.0, 0.1, 1.0
Z_MIN, Z_MAX = 0.0, 0.1
SHAPE_BOUND = 0.05
ALPHA_BOUNDS = (0.0, 10.0)
THICKNESS = (0.5, 3.0)
N_CHORD, N_SPAN = 100, 2
OBJ_SCALE = 1e4
DIVERGED_CD = 0.05          # 500 counts: stands in for the cd of a diverged solve (see --divergedTol)

p = common.baseParser(__doc__)
p.add_argument("--points", default="0.3:0.5", help="'Tu:cl,Tu:cl,...' with Tu in percent")
p.add_argument("--weights", default=None, help="comma list, default equal")
p.add_argument("--dissipation", choices=["scalar", "matrix"], default="scalar")
p.add_argument("--level", default="L0")
p.add_argument("--ffd", default=os.path.join(common.FFD_DIR, "naca0012_ffd_10x2x2.xyz"))
p.add_argument("--startShapeTxt", default=None, help="start from this shape vector (np.savetxt) instead of the NACA 0012")
p.add_argument("--startAlpha", type=float, default=2.7, help="alpha at cl 0.3; other cl are shifted by the lift slope")
p.add_argument("--shapeScale", type=float, default=200.0, help="optimiser scaling of the shape variables")
p.add_argument("--acc", type=float, default=None, help="SLSQP ACC (default 1e-8 scalar, 1e-7 matrix)")
p.add_argument("--maxIt", type=int, default=200)
p.add_argument("--hotStart", default=None, help="opt.hst of the previous round (in another directory)")
p.add_argument("--divergedTol", type=float, default=1e-5,
               help="a failed solve with res/r0 above this gets cd = 500 counts so that SLSQP backs off "
                    "(pySLSQP ignores the fail flag); GR in matrix used 2e-5")
p.add_argument("--baselineOnly", action="store_true", help="build the problem, evaluate the start design, stop")
args = p.parse_args()

points = [tuple(float(v) for v in tok.split(":")) for tok in args.points.split(",")]
nPts = len(points)
weights = [1.0 / nPts] * nPts if args.weights is None else [float(w) for w in args.weights.split(",")]
assert len(weights) == nPts and abs(sum(weights) - 1.0) < 1e-6
world = MPI.COMM_WORLD
assert world.size % nPts == 0, f"{world.size} ranks cannot be split equally over {nPts} points"

from multipoint import multiPointSparse  # noqa: E402

MP = multiPointSparse(world)
MP.addProcessorSet("cruise", nMembers=nPts, memberSizes=[world.size // nPts] * nPts)
comm, _, _, _, ptID = MP.createCommunicators()
tuPt, clPt = points[ptID]

grid = args.mesh or common.meshPath("naca0012", f"mesh_{args.level}", "volumeMesh.cgns")
outDir = args.outputDir or f"output_opt_airfoil_{args.model}_{args.dissipation}"
recipe = args.recipe or recipes.recipe2D(args.model, args.dissipation, purpose="opt")
opts = recipes.options(args.model, recipe, args.dissipation, tuPt, grid, outDir, nuinf=args.nuinf, liftIndex=2)
overrides = common.applyOverrides(opts, args)
if world.rank == 0:
    os.makedirs(outDir, exist_ok=True)
common.banner(world, f"NACA 0012 optimisation | {recipes.MODEL_LABEL[args.model]} | {args.dissipation} dissipation | "
                     f"points (Tu %, cl) {points} weights {weights} | {world.size // nPts} ranks per point\n"
                     f"recipe {recipe}: {recipes.DESCRIPTION[recipe]}\noverrides {overrides}")

from baseclasses import AeroProblem  # noqa: E402
from idwarp import USMesh  # noqa: E402
from pygeo import DVConstraints, DVGeometry  # noqa: E402
from pyoptsparse import OPT, Optimization  # noqa: E402

CFD = common.makeSolver(opts, comm)
CFD.addSlices("z", [0.5 * (Z_MIN + Z_MAX)], sliceType="absolute")

# ---- flow conditions: one AeroProblem per point, alpha is a design variable of each
aps = []
for i, (tu, cl) in enumerate(points):
    a0 = args.startAlpha + (cl - 0.3) / 0.1095
    ap = AeroProblem(name=f"cruise{i}", alpha=a0, mach=MACH, altitude=ALTITUDE, areaRef=AREA_REF,
                     chordRef=CHORD_REF, evalFuncs=["cl", "cd"])       # two adjoints per gradient
    ap.addDV("alpha", value=a0, lower=ALPHA_BOUNDS[0], upper=ALPHA_BOUNDS[1], scale=1.0)
    aps.append(ap)
ap = aps[ptID]

# ---- geometry: FFD shape variables
DVGeo = DVGeometry(args.ffd)
DVGeo.addLocalDV("shape", lower=-SHAPE_BOUND, upper=SHAPE_BOUND, axis="y", scale=args.shapeScale)
if args.startShapeTxt:
    DVGeo.setDesignVars({"shape": np.loadtxt(args.startShapeTxt).ravel()})
CFD.setDVGeo(DVGeo)

# ---- geometric constraints
DVCon = DVConstraints()
DVCon.setDVGeo(DVGeo)
DVCon.setSurface(CFD.getTriangulatedMeshSurface())
lIndex = DVGeo.getLocalIndex(0)
# leading and trailing edge fixed (upper and lower control point move in opposite directions)
DVCon.addLeTeConstraints(0, indSetA=[lIndex[0, 0, 0], lIndex[-1, 0, 0]], indSetB=[lIndex[0, 1, 0], lIndex[-1, 1, 0]])
# 2D: the two span planes of the FFD move together
indA = [lIndex[i, j, 0] for j in range(2) for i in range(lIndex.shape[0])]
indB = [lIndex[i, j, 1] for j in range(2) for i in range(lIndex.shape[0])]
DVCon.addLinearConstraintsShape(indA, indB, factorA=1.0, factorB=-1.0, lower=0, upper=0)
eps = 1e-4
leList = [[eps, 0.0, Z_MIN + eps], [eps, 0.0, Z_MAX - eps]]
teList = [[1.0 - eps, 0.0, Z_MIN + eps], [1.0 - eps, 0.0, Z_MAX - eps]]
DVCon.addVolumeConstraint(leList, teList, N_SPAN, N_CHORD, lower=1.0, scaled=True)
DVCon.addThicknessConstraints2D(leList, teList, N_SPAN, N_CHORD, lower=THICKNESS[0], upper=THICKNESS[1])
CFD.setMesh(USMesh(options={"gridFile": grid}, comm=comm))


def cruiseFuncs(x):
    """One design point: deform, cold start, solve, return cl / cd and the geometric constraints."""
    DVGeo.setDesignVars(x)
    ap.setDesignVars(x)
    CFD.resetFlow(ap)
    CFD(ap)
    funcs = {}
    DVCon.evalFunctions(funcs)
    CFD.evalFunctions(ap, funcs)
    CFD.checkSolutionFailure(ap, funcs)
    split = {}
    CFD.evalFunctions(ap, split, evalFuncs=["cdp", "cdv"])     # reported only, no adjoint
    r0, _, rf = CFD.getResNorms()
    rel = float(rf / r0)
    if comm.rank == 0:
        print(f"[point {ptID}: Tu {tuPt} %, cl* {clPt}] alpha={float(np.real(ap.alpha)):.5f} "
              f"cl={funcs[ap.name + '_cl']:.8f} cd={funcs[ap.name + '_cd'] * 1e4:.4f} counts "
              f"(cdp {split[ap.name + '_cdp'] * 1e4:.4f}, cdv {split[ap.name + '_cdv'] * 1e4:.4f}) "
              f"res/r0={rel:.2e} totalR0={float(r0):.16e} fail={funcs.get('fail', False)}", flush=True)
    if funcs.get("fail", False) and not rel <= args.divergedTol and not args.baselineOnly:
        funcs[ap.name + "_cd"] = DIVERGED_CD
    return funcs


def cruiseFuncsSens(x, funcs):
    funcsSens = {}
    DVCon.evalFunctionsSens(funcsSens)
    CFD.evalFunctionsSens(ap, funcsSens)
    CFD.checkAdjointFailure(ap, funcsSens)
    return funcsSens


def objCon(funcs, printOK):
    funcs["obj"] = sum(w * funcs[a.name + "_cd"] for w, a in zip(weights, aps))
    for a, (_, cl) in zip(aps, points):
        funcs["cl_con_" + a.name] = funcs[a.name + "_cl"] - cl
    if printOK:
        print("objective and constraints:", {k: funcs[k] for k in sorted(funcs) if k == "obj" or "cruise" in k}, flush=True)
    return funcs


optProb = Optimization(f"naca0012_{args.model}_{nPts}pt", MP.obj, comm=world)
optProb.addObj("obj", scale=OBJ_SCALE)
for a in aps:
    a.addVariablesPyOpt(optProb)
DVGeo.addVariablesPyOpt(optProb)
DVCon.addConstraintsPyOpt(optProb)
for a in aps:
    optProb.addCon("cl_con_" + a.name, lower=0.0, upper=0.0, scale=1.0, wrt=["alpha_" + a.name, "shape"])
MP.setProcSetObjFunc("cruise", cruiseFuncs)
MP.setProcSetSensFunc("cruise", cruiseFuncsSens)
MP.setObjCon(objCon)
MP.setOptProb(optProb)

if args.baselineOnly:
    # the start design only: no optimiser, no adjoint
    funcs, fail = MP.obj(optProb.getDVs())
    if world.rank == 0:
        print(f"BASELINE obj={funcs['obj'] * 1e4:.4f} counts fail={fail}", flush=True)
        common.writeResult(world, opts, dict(script="opt_airfoil.py", mode="baselineOnly", model=args.model,
                                             dissipation=args.dissipation, points=points, weights=weights, recipe=recipe),
                           [{k: float(np.real(v)) for k, v in funcs.items() if np.size(v) == 1}],
                           common.provenance(world, grid), overrides)
else:
    acc = args.acc if args.acc is not None else (1e-7 if args.dissipation == "matrix" else 1e-8)
    opt = OPT("SLSQP", options={"ACC": acc, "MAXIT": args.maxIt, "IPRINT": 1, "IOUT": 60,
                                "IFILE": os.path.join(outDir, "SLSQP.out")})
    sol = opt(optProb, MP.sens, storeHistory=os.path.join(outDir, "opt.hst"), hotStart=args.hotStart)
    if world.rank == 0:
        print(sol)
        np.savetxt(os.path.join(outDir, "xstar.shape.txt"), np.asarray(sol.xStar["shape"]).ravel())
        with open(os.path.join(outDir, "xstar.alpha.txt"), "w") as f:
            for a in aps:
                f.write(f"{float(np.ravel(sol.xStar['alpha_' + a.name])[0]):.6f}\n")
