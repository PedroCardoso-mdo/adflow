#!/usr/bin/env python
"""AR5 wing + winglet analysis: GR (campaign 11 refinement study) and SA-BCM (article wing).

What it demonstrates
    The transition models on a real 3D wing. Two validated set-ups, selected by --model:

    --model gr    "corrected foil" AR5 wing + winglet, 5-level mesh family, crossflow OFF.
                  M 0.08164, altitude 304.8 m (Re 1.127e6 on the mean chord 0.6091 m), alpha 1.8 deg,
                  Tu 0.25 %, ADflow default scalar dissipation (vis4 0.0156) and default
                  eddyVisInfRatio 0.009 -- that is what the refinement study ran.
    --model bcm   SA-BCM article wing (943 104 cells). M 0.133, altitude 304.8 m, alpha 1.057684 deg,
                  Tu 0.25 %, eddyVisInfRatio 1e-7, acousticScaleFactor = Mach, default vis4.
                  (bcmo = original SA-BCM, sa = fully turbulent, on the same mesh and recipe.)

Meshes (showcase/meshes/ar5_wing/)
    corrected_foil_ar5_vol_L4.cgns   459 452 cells        corrected_foil_ar5_vol_L2.cgns  1 783 442
    corrected_foil_ar5_vol_L3.cgns   904 134              corrected_foil_ar5_vol_L1.cgns  3 667 320
    (L0, 7 415 408 cells, is on the cluster only and needs 2 nodes)
    ar5_asa_L2_943104cells.cgns      943 104 cells (SA-BCM article)

Reference numbers
    GR refinement, alpha 1.8, Tu 0.25 % (jobs 1810888-92 + no-NK restarts 1810916-25, 64 ranks):
        level  orders  cl        cd        cdp       cdv       transition onset x/c (eta 0.05)
        L4     8.00    0.90955   0.045928  0.042343  0.003584  0.23
        L3     8.07    0.93011   0.041680  0.038444  0.003236  0.28
        L2     8.04    0.94849   0.038856  0.035616  0.003240  0.28
        L1     8.05    0.95822   0.035907  0.032626  0.003281  0.33
        L0     7.61    0.96054   0.033905  0.030429  0.003476  0.33   (not converged)
        Status: forces converged on 4 of 5 levels, but NOT mesh-converged (cd still moves 5.6 %
        from L1 to L0), and the transition front needs >= 397 wrap points (this family has 181):
        on the finer chordwise meshes it sits at x/c 0.62, not 0.33. Use as a solver demonstration.
        L1: 1 h 34 min on 64 ranks for both legs. The second leg is a restart without NK:
            --restartFile <outputDir>/ar5_gr_000_vol.cgns --opt useNKSolver=False
    SA-BCM article baseline (job 1904793, 64 ranks, 163 iterations, 312 s, converged to 1e-6):
        cl 0.894395  cd 0.02445674  (cdp 0.02122303, cdv 0.00323371)
        Fully turbulent SA on the same wing, trimmed: cl 0.893208, cd 0.03244747.
        Status: reproducible numbers; the article's claims about the smooth form were later
        withdrawn (see README.md).

Three models on one mesh (--dissipation matrix)
    The 397-wrap mesh c376 (13 684 496 cells; the coarsest mesh of the chordwise ladder on which the
    transition front has settled, x/c 0.62), matrix dissipation, M 0.08164, Tu 0.4 %, alpha 1.13249 deg,
    for --model gr | sgamma | bcm | bcmo | sa. Needs 2 cluster nodes (about 4 GB per rank at 128 ranks).
    Free-stream nu_tilde: the validated value of each model (recipes.NUINF).
    Runs of August 2026 at this condition (128 ranks on 2 nodes), with other nu_tilde values
    (GR 0.009, SA-BCM 5e-7; pass --nuinf to repeat them):
        GR, crossflow off  cl 0.89970  cd 0.028352  7.7 orders  (restart from a GR solution, 3.8 h)
        GR, crossflow on   cl 0.89983  cd 0.028349  8.7 orders  (restart from the run above, 3.9 h)
        SA-BCM smooth      cl 0.89735  cd 0.028290  6.5 orders  (cold start, 9.2 h, stopped)
    SA-gamma has no validated 3D recipe: its run is a first attempt.

Commands
    mpirun -np 64 python run_ar5_wing.py --model gr --level L1 --writeVolume
    mpirun -np 64 python run_ar5_wing.py --model bcm
    mpirun -np 6 --bind-to core python run_ar5_wing.py --model gr --level L4 --nCycles 4     # smoke
    mpirun -np 6 --bind-to core python run_ar5_wing.py --model bcm --nCycles 5               # smoke
    srun -n 128 python run_ar5_wing.py --model gr --dissipation matrix --writeVolume          # 2 nodes
"""
from mpi4py import MPI

import common
import recipes

AREA_REF, CHORD_REF, ALTITUDE = 2.1691, 0.6091, 304.8      # half-model area, mean aerodynamic chord

p = common.baseParser(__doc__, models=("gr", "sgamma", "bcm", "bcmo", "sa"))
p.add_argument("--level", default="L1", choices=["L1", "L2", "L3", "L4"], help="GR mesh level")
p.add_argument("--dissipation", default="scalar_default", choices=["scalar_default", "matrix"],
               help="matrix = the three-model benchmark on the 13.7 M mesh (see the docstring)")
p.add_argument("--alpha", type=float, default=None)
p.add_argument("--tu", type=float, default=None, help="turbulence intensity [%%] (default 0.25; 0.4 with --dissipation matrix)")
p.add_argument("--crossflow", action="store_true", help="GR only: switch the crossflow term on")
p.add_argument("--restartFile", default=None, help="volume solution to restart from")
args = p.parse_args()

comm = MPI.COMM_WORLD
if args.dissipation == "matrix":
    # One mesh, one condition, every model: the 397-wrap mesh c376 (13 684 496 cells), M 0.08164,
    # Tu 0.4 %, alpha 1.13249 deg (the angle that trims GR to cl 0.8932 in scalar dissipation).
    # H/08_optimization/3d_c376/run_c376_mtx_cross.py (GR), L/15_sabcm_polars/jobs/job_c376_bcm_mtx.sh (SA-BCM).
    grid = args.mesh or common.meshPath("ar5_wing", "chord_ladder_c376", "corrected_foil_ar5_vol_c376.cgns")
    mach, alpha = 0.08164, 1.1324921152893268 if args.alpha is None else args.alpha
    tu = 0.4 if args.tu is None else args.tu
    recipe = args.recipe or {"gr": "gr_3d_wing_matrix", "sgamma": "sgamma_3d_wing_matrix"}.get(args.model, "bcm_3d_wing_matrix")
    nuinf = recipes.NUINF[args.model] if args.nuinf is None else args.nuinf   # validated value of each model
    extra = {}
    if args.model in ("gr", "sgamma"):
        extra.update(volumeVariables=["intermittency", "fonset"])
    if args.model == "gr":
        extra.update(transitionCrossflow=args.crossflow)
    if args.model in ("bcm", "bcmo"):
        extra.update(acousticScaleFactor=0.15)
elif args.model == "gr":
    # L/11_ar5_corrected_foil/common/run_ar5cf.py:55-61, batch_jobs/ar5cf.sh
    grid = args.mesh or common.meshPath("ar5_wing", f"corrected_foil_ar5_vol_{args.level}.cgns")
    mach, alpha = 0.08164, 1.8 if args.alpha is None else args.alpha
    tu = 0.25 if args.tu is None else args.tu
    recipe = args.recipe or "gr_3d_wing_ar5"
    nuinf = 0.009 if args.nuinf is None else args.nuinf      # ADflow default: what campaign 11 ran
    extra = dict(volumeVariables=["intermittency", "fonset"], transitionCrossflow=args.crossflow)
elif args.model == "sgamma":
    raise SystemExit("--model sgamma on this wing exists only as --dissipation matrix")
else:
    # showcase/bcm_article/cluster_package/06_opt_3d/02_c2_opt_tu025/{run.py,vars_master.py}
    grid = args.mesh or common.meshPath("ar5_wing", "ar5_asa_L2_943104cells.cgns")
    mach, alpha = 0.133, 1.057684 if args.alpha is None else args.alpha
    tu = 0.25 if args.tu is None else args.tu
    recipe = args.recipe or ("sa_3d_wing" if args.model == "sa" else "bcm_3d_wing")
    nuinf = (1e-10 if args.model == "sa" else 1e-7) if args.nuinf is None else args.nuinf
    extra = dict(acousticScaleFactor=mach)
outDir = args.outputDir or f"output_ar5_{args.model}"

opts = recipes.options(args.model, recipe, args.dissipation, tu, grid, outDir, nuinf=nuinf,
                       liftIndex=3, surfaceVariables=["cp", "mach", "cf", "cfx", "cfy", "cfz", "yplus", "blank"], **extra)
if args.restartFile:
    opts["restartFile"] = args.restartFile
overrides = common.applyOverrides(opts, args)

common.banner(comm, f"AR5 wing | {recipes.MODEL_LABEL[args.model]} | M {mach} alpha {alpha} Tu {tu} % | {args.dissipation} | "
                    f"nu_inf ratio {nuinf:g}\nmesh {grid}\nrecipe {recipe}: {recipes.DESCRIPTION[recipe]}\noverrides {overrides}")

from baseclasses import AeroProblem  # noqa: E402

CFD = common.makeSolver(opts, comm)
ap = AeroProblem(name=f"ar5_{args.model}", alpha=alpha, mach=mach, altitude=ALTITUDE,
                 areaRef=AREA_REF, chordRef=CHORD_REF, evalFuncs=["cl", "cd", "cdp", "cdv"])
records = [common.solve(CFD, ap, comm)]
common.writeResult(comm, opts, dict(script="run_ar5_wing.py", model=args.model, recipe=recipe, mach=mach,
                                    tuPercent=tu, dissipation=args.dissipation, reynolds=float(ap.re * CHORD_REF)),
                   records, common.provenance(comm, grid), overrides)
