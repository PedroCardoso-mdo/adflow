#!/usr/bin/env python
"""3D wing shape optimisation set-ups on the AR5 wing: the GR attempt and the SA-BCM article case.

    minimise cd  w.r.t. FFD shape variables + alpha  subject to the lift and thickness constraints, SLSQP

--case gr_c376   GR on the 397-wrap mesh "c376" (13 684 496 cells) -- an ATTEMPT, not finished
    Mesh   showcase/meshes/ar5_wing/chord_ladder_c376/corrected_foil_ar5_vol_c376.cgns (321 MB)
           (the coarsest mesh of the chordwise ladder that resolves the transition front)
    FFD    ffd/ar5_corrected_foil_ffd_7x13x2.xyz : 182 shape variables along z, bounds +-0.02 m
    Flow   M 0.08164, 304.8 m, Tu 0.5 % (also 0.4 and 0.25 %), start alpha 1.1325 deg, cl = 0.8932
    Model  GR, crossflow off, scalar vis4 0.01, eddyVisInfRatio 5e-7; primal to 1e-7, adjoint to 1e-3
    Constraints  volume >= initial, thickness 0.9-1.5 x initial at 8 x 25 stations on the wing
           (y 0.05 .. 2.90 m; the winglet is held by the bounds only), leading / trailing edge fixed
    Results (H/08_optimization/3d_c376; 4 nodes x 32 ranks, 48 h per job, ~5 h per primal and
           ~2 h per adjoint pair):
        Tu 0.4 % (jobs 1830863, 1843539, 1844851, 1849924): cd 287.0 -> 269.0 counts (-6.3 %) in
            6 major iterations, but cl 0.88728 (-0.66 % of the target): part of the gain is less lift.
        Tu 0.5 % (jobs 1931137, 1939946): 287.1 -> 269.5 (-6.1 %), cl 0.88834, 12 major iterations;
            the last evaluation diverged.
        Fully turbulent SA control (job 1927255): 362.8 -> 345.2 (-4.9 %), cl 0.89320, converged.
    Status: attempt. The lift constraint is not closed, no optimum was re-solved cold, there is no
           adjoint-versus-complex-step check on this mesh; the campaign stopped when the quota ran out.

--case bcm_article   SA-BCM article optimisation (943 104 cells)
    Mesh   showcase/meshes/ar5_wing/ar5_asa_L2_943104cells.cgns
    FFD    ffd/ar5_article_ffd_12x29x2.xyz, section-normal shape variables, bounds +-0.05
    Flow   M 0.133, 304.8 m, Tu 0.25 %, start alpha 1.057684 deg, cl >= 0.8932
    Constraints  thickness >= 0.9 x initial at 19 x 6 stations, leading / trailing edge fixed (no volume constraint)
    Results (showcase/bcm_article/cluster_package/06_opt_3d/TABLE_3D_TU025.md; 64 ranks, 12.3 h):
        baseline SA-BCM  cd 0.02445674 (cl 0.894395, alpha fixed)
        C2 optimum       cd 0.02405300 (-1.65 %; -1.4 % against a trimmed baseline), 11 + 19 SLSQP iterations
        C3 from the SA optimum  cd 0.02397400 (-1.97 %), SLSQP Inform 0
        SA: 0.03244747 -> 0.03134485 (-3.40 %)
    Status: reproducible. The article's central claims about the smooth form (convergence and
           derivative accuracy against the original model) were later withdrawn; see README.md.

Commands
    mpirun -np 128 python opt_ar5_wing.py --case gr_c376 --tu 0.5           # 4 nodes x 32 ranks
    mpirun -np 64 python opt_ar5_wing.py --case bcm_article
    mpirun -np 6 --bind-to core python opt_ar5_wing.py --case bcm_article --baselineOnly --nCycles 5   # smoke
    mpirun -np 6 --bind-to core python opt_ar5_wing.py --case gr_c376 --baselineOnly --nCycles 4 \\
        --mesh ../../showcase/meshes/ar5_wing/corrected_foil_ar5_vol_L4.cgns    # smoke on a coarse mesh of the family
"""
import os

import numpy as np
from mpi4py import MPI

import common
import recipes

AREA_REF, CHORD_REF, ALTITUDE, CL_STAR = 2.1691, 0.6091, 304.8, 0.8932

p = common.baseParser(__doc__, models=("gr", "bcm", "bcmo", "sa"), defaultModel=None)
p.add_argument("--case", choices=["gr_c376", "bcm_article"], required=True)
p.add_argument("--tu", type=float, default=None, help="turbulence intensity [%%] (default 0.5 gr_c376, 0.25 bcm_article)")
p.add_argument("--alpha0", type=float, default=None, help="start alpha [deg]")
p.add_argument("--hotStart", default=None, help="opt.hst of the previous round (in another directory)")
p.add_argument("--baselineOnly", action="store_true", help="build the problem, evaluate the start design, stop")
args = p.parse_args()

world = MPI.COMM_WORLD
from multipoint import multiPointSparse  # noqa: E402

MP = multiPointSparse(world)
MP.addProcessorSet("cruise", nMembers=1, memberSizes=world.size)
comm, _, _, _, _ = MP.createCommunicators()

if args.case == "gr_c376":
    # H/08_optimization/3d_c376/{vars_master.py,solver_setup.py,run_opt.py}
    model = args.model or "gr"
    grid = args.mesh or common.meshPath("ar5_wing", "chord_ladder_c376", "corrected_foil_ar5_vol_c376.cgns")
    ffd = os.path.join(common.FFD_DIR, "ar5_corrected_foil_ffd_7x13x2.xyz")
    mach, alpha0, tu = 0.08164, 1.1325, 0.5
    recipe = args.recipe or "gr_3d_wing_ar5_opt"
    dissipation = "scalar"                                     # vis4 0.01
    nuinf = (0.009 if model == "sa" else 5e-7) if args.nuinf is None else args.nuinf
    extra = dict(volumeVariables=["intermittency", "fonset", "eddyratio"] if model == "gr" else ["eddyratio"])
    alphaBounds, optAcc, optMaxIt, objScale, clScale, clUpper = (-3.0, 6.0), 1e-6, 60, 1e4, 1.0, 0.0
else:
    # showcase/bcm_article/cluster_package/06_opt_3d/02_c2_opt_tu025/{run.py,vars_master.py}
    model = args.model or "bcm"
    grid = args.mesh or common.meshPath("ar5_wing", "ar5_asa_L2_943104cells.cgns")
    ffd = os.path.join(common.FFD_DIR, "ar5_article_ffd_12x29x2.xyz")
    mach, alpha0, tu = 0.133, 1.057684, 0.25
    recipe = args.recipe or ("sa_3d_wing" if model == "sa" else "bcm_3d_wing")
    dissipation = "scalar_default"                             # ADflow default vis4 0.0156
    nuinf = (1e-10 if model == "sa" else 1e-7) if args.nuinf is None else args.nuinf
    extra = dict(acousticScaleFactor=mach)
    alphaBounds, optAcc, optMaxIt, objScale, clScale, clUpper = (0.0, 15.0), 1e-8, 500, 1e2, 50.0, None
tu = tu if args.tu is None else args.tu
alpha0 = alpha0 if args.alpha0 is None else args.alpha0
outDir = args.outputDir or f"output_opt_ar5_{args.case}"

opts = recipes.options(model, recipe, dissipation, tu, grid, outDir, nuinf=nuinf, liftIndex=3,
                       surfaceVariables=["cp", "mach", "cf", "cfx", "cfy", "cfz", "yplus", "blank"], **extra)
if model == "sa" and args.case == "gr_c376":
    # fully turbulent control on the GR case: the recipe keys that only exist for 3 equations go
    opts.pop("turbResScale", None)
overrides = common.applyOverrides(opts, args)
if world.rank == 0:
    os.makedirs(outDir, exist_ok=True)
common.banner(world, f"AR5 wing optimisation, case {args.case} | {recipes.MODEL_LABEL[model]} | M {mach} Tu {tu} % "
                     f"alpha0 {alpha0} cl* {CL_STAR}\nmesh {grid}\nrecipe {recipe}: {recipes.DESCRIPTION[recipe]}"
                     f"\noverrides {overrides}")

from baseclasses import AeroProblem  # noqa: E402
from idwarp import USMesh  # noqa: E402
from pygeo import DVConstraints, DVGeometry  # noqa: E402
from pyoptsparse import OPT, Optimization  # noqa: E402

CFD = common.makeSolver(opts, comm)
ap = AeroProblem(name="cruise", alpha=alpha0, mach=mach, altitude=ALTITUDE, areaRef=AREA_REF, chordRef=CHORD_REF,
                 evalFuncs=["cl", "cd"])
ap.addDV("alpha", value=alpha0, lower=alphaBounds[0], upper=alphaBounds[1], scale=1.0)

DVGeo = DVGeometry(ffd)
DVCon = DVConstraints()
if args.case == "gr_c376":
    # ---- shape: every FFD control point along the lift direction (z)
    DVGeo.addLocalDV("shape", lower=-0.02, upper=0.02, axis="z", scale=25.0)
    CFD.setDVGeo(DVGeo)
    DVCon.setDVGeo(DVGeo)
    surf = CFD.getTriangulatedMeshSurface()
    DVCon.setSurface(surf)
    lIndex = DVGeo.getLocalIndex(0)
    indA, indB = [], []
    for j in range(lIndex.shape[1]):                 # leading and trailing edge of every span station
        indA += [lIndex[0, j, 0], lIndex[-1, j, 0]]
        indB += [lIndex[0, j, 1], lIndex[-1, j, 1]]
    DVCon.addLeTeConstraints(0, indSetA=indA, indSetB=indB)
    # ---- constraint box on the main wing, with the LE / TE curves measured on the wall mesh
    # (getTriangulatedMeshSurface returns p0 and two edge VECTORS; the FFD box is larger than the wing)
    p0, v1, v2 = (np.asarray(a) for a in surf[:3])
    pts = np.vstack([p0, p0 + v1, p0 + v2])
    leList, teList = [], []
    for yS in np.linspace(0.05, 2.90, 8):
        band = np.abs(pts[:, 1] - yS) < 0.01 * pts[:, 1].max()
        xl, xt = pts[band, 0].min(), pts[band, 0].max()
        zc = 0.5 * (pts[band, 2].min() + pts[band, 2].max())
        leList.append([xl + 0.03 * (xt - xl), yS, zc])
        teList.append([xt - 0.03 * (xt - xl), yS, zc])
    DVCon.addVolumeConstraint(leList, teList, 8, 25, lower=1.0, scaled=True)
    DVCon.addThicknessConstraints2D(leList, teList, 8, 25, lower=0.9, upper=1.5)
    meshOptions = {"gridFile": grid}
else:
    # ---- shape: control points along the local section normal
    DVGeo.addLocalSectionDV("shape", secIndex="j", lower=-0.05, upper=0.05, axis=1, scale=10)
    CFD.setDVGeo(DVGeo)
    DVCon.setDVGeo(DVGeo)
    DVCon.setSurface(CFD.getTriangulatedMeshSurface())
    m = 0.02
    leList = [[1.5546 + m, 0.0002, 0.30705], [1.5747 + m, 1.1900, 0.30608], [1.5933 + m, 3.2900, 0.30129]]
    teList = [[2.2538 - m, 0.0002, 0.27605], [2.2537 - m, 1.1900, 0.27608], [2.0711 - m, 3.2900, 0.29653]]
    DVCon.addThicknessConstraints2D(leList, teList, nSpan=19, nChord=6, lower=0.9, scaled=True)
    DVCon.addLeTeConstraints(0, "iHigh")
    meshOptions = {"gridFile": grid, "LdefFact": 100.0}
CFD.setMesh(USMesh(options=meshOptions, comm=comm))


def cruiseFuncs(x):
    DVGeo.setDesignVars(x)
    ap.setDesignVars(x)
    CFD.resetFlow(ap)                               # cold start at every design point
    CFD(ap)
    funcs = {}
    DVCon.evalFunctions(funcs)
    CFD.evalFunctions(ap, funcs)
    CFD.checkSolutionFailure(ap, funcs)
    split = {}
    CFD.evalFunctions(ap, split, evalFuncs=["cdp", "cdv"])
    r0, _, rf = CFD.getResNorms()
    if comm.rank == 0:
        print(f"[design point] alpha={float(np.real(ap.alpha)):.5f} cl={funcs['cruise_cl']:.8f} cd={funcs['cruise_cd']:.8f} "
              f"(cdp {split['cruise_cdp']:.8f}, cdv {split['cruise_cdv']:.8f}) res/r0={float(rf / r0):.2e} "
              f"totalR0={float(r0):.16e} fail={funcs.get('fail', False)}", flush=True)
    if args.case == "bcm_article" and funcs.get("fail", False) and not args.baselineOnly:
        funcs["cruise_cd"] = max(funcs["cruise_cd"] * 10.0, 1.0)    # article script: make the line search back off
    return funcs


def cruiseFuncsSens(x, funcs):
    funcsSens = {}
    DVCon.evalFunctionsSens(funcsSens)
    CFD.evalFunctionsSens(ap, funcsSens)
    CFD.checkAdjointFailure(ap, funcsSens)
    return funcsSens


def objCon(funcs, printOK):
    funcs["obj"] = funcs["cruise_cd"]
    funcs["cl_con"] = funcs["cruise_cl"] - CL_STAR
    if printOK:
        print("objective and lift constraint:", funcs["obj"], funcs["cl_con"], flush=True)
    return funcs


optProb = Optimization(f"ar5_{args.case}", MP.obj, comm=world)
optProb.addObj("obj", scale=objScale)
ap.addVariablesPyOpt(optProb)
DVGeo.addVariablesPyOpt(optProb)
DVCon.addConstraintsPyOpt(optProb)
optProb.addCon("cl_con", lower=0.0, upper=clUpper, scale=clScale)     # equality (GR) or cl >= cl* (article)
MP.setProcSetObjFunc("cruise", cruiseFuncs)
MP.setProcSetSensFunc("cruise", cruiseFuncsSens)
MP.setObjCon(objCon)
MP.setOptProb(optProb)

if args.baselineOnly:
    funcs, fail = MP.obj(optProb.getDVs())
    if world.rank == 0:
        print(f"BASELINE cd={funcs['obj']:.8f} cl-cl*={funcs['cl_con']:.3e} nDV={DVGeo.getNDV() + 1} fail={fail}", flush=True)
        common.writeResult(world, opts, dict(script="opt_ar5_wing.py", mode="baselineOnly", case=args.case, model=model,
                                             recipe=recipe, tuPercent=tu, nShapeDV=DVGeo.getNDV()),
                           [{k: float(np.real(v)) for k, v in funcs.items() if np.size(v) == 1}],
                           common.provenance(world, grid), overrides)
else:
    opt = OPT("SLSQP", options={"ACC": optAcc, "MAXIT": optMaxIt, "IPRINT": 1, "IOUT": 60,
                                "IFILE": os.path.join(outDir, "SLSQP.out")})
    sol = opt(optProb, MP.sens, storeHistory=os.path.join(outDir, "opt.hst"), hotStart=args.hotStart)
    if world.rank == 0:
        print(sol)
