#!/usr/bin/env python
"""2D airfoil analysis with any of the transition models (fixed alpha, one or several).

What it demonstrates
    The three transition models (and plain SA) on the same meshes, conditions and dissipation:
    NACA 0012 at the article condition, and the NLF(1)-0416 / S809 validation polars of
    Piotrowski & Zingg (2020).

Meshes (showcase/meshes/<airfoil>/mesh_<level>/volumeMesh.cgns; pyHyp O-grids, span 0.1, far field 500 c)
    L0 = 193 224 cells (used below), L1 = 98 736, L-1 = 381 240, L-2 = 753 989.

Conditions
    naca0012  M 0.1252, sea level (Re 2.92e6), Tu 0.5 %   [article point A: --tu 0.3 --alpha 4.52]
    nlf0416   M 0.10, Re 4.0e6, T 288.15 K, Tu 0.15 %
    s809      M 0.10, Re 2.0e6, T 288.15 K, Tu 0.07 %

Reference numbers (cd in counts; every one converged to res/r0 <= 1e-8 unless stated)
    NACA 0012 L0, alpha 2.7, scalar (jobs 1962375-79, 32 ranks):
        GR     cl 0.30113  cd 62.19 (cdp 14.24, cdv 47.95)   299 s
        SA-g   cl 0.30239  cd 67.13 (cdp 14.71, cdv 52.42)   506 s
        BCM    cl 0.30417  cd 66.01 (cdp 14.39, cdv 51.62)   444 s
    NACA 0012 L0, alpha 4.52, Tu 0.3 %, scalar (job 2197279, 8 ranks, ~30 min):
        GR 65.06 (cl 0.5000) | SA-g 68.83 | BCM original 65.76 | SA 102.02
    NLF0416 L0, alpha 4:
        GR scalar  cl 0.9631 cd 65.58 (job 1961835, 32 ranks, 572 s)   paper 65.40, experiment 64.36
        GR matrix  cl 0.9611 cd 65.22 (job 2201491, 10 ranks, 2798 s)
        SA-g scalar 69.28 (nuinf 2.79e-7, job 1961971) | BCM 64.87 (job 1961959) | BCM original 64.73 | SA 112.76
    S809 L0, alpha 4:
        GR scalar  cl 0.6549 cd 66.42 (job 1961835, 893 s) | GR matrix cl 0.6557 cd 63.86 (job 2201434, 4219 s)
    Polars (alpha -6..6 for NLF0416, 0..9 for S809): GR 23/23 and SA-BCM 23/23 points at 1e-8 in both
    schemes; SA-gamma 16-17/23 in scalar (the rest at 3-6 orders). Full tables: showcase/_manifest_2d.md.
    Status: validated.

Cluster cost: one point = 8-32 ranks, 5-70 min (128-core node: 4 x 32 or 12 x 10 points at once).

Commands
    mpirun -np 32 python run_airfoil.py --airfoil naca0012 --model gr --alpha 2.7
    mpirun -np 32 python run_airfoil.py --airfoil nlf0416 --model gr --dissipation matrix --alpha 4
    mpirun -np 32 python run_airfoil.py --airfoil s809 --model bcmo --alpha 0 2 4 6 8
    mpirun -np 6 --bind-to core python run_airfoil.py --airfoil nlf0416 --model gr --alpha 4 --nCycles 4   # smoke
"""
from mpi4py import MPI

import common
import recipes

AREA_REF, CHORD_REF, LIFT_INDEX, SLICE_Z = 0.1, 1.0, 2, 0.05   # chord 1, span 0.1 (z), lift in y
AIRFOILS = {
    # L/24_minmodel_2d/config/setup.py:63-65
    "naca0012": dict(mach=0.125191453, altitude=0.0, tu=0.5, alpha=[2.7]),
    "nlf0416": dict(mach=0.1, reynolds=4.0e6, T=288.15, tu=0.15, alpha=[4.0]),
    "s809": dict(mach=0.1, reynolds=2.0e6, T=288.15, tu=0.07, alpha=[4.0]),
}

p = common.baseParser(__doc__)
p.add_argument("--airfoil", choices=sorted(AIRFOILS), default="naca0012")
p.add_argument("--level", default="L0", help="mesh level of the family (L0 = 193 224 cells)")
p.add_argument("--dissipation", choices=["scalar", "matrix"], default="scalar")
p.add_argument("--alpha", type=float, nargs="+", default=None, help="angle(s) of attack [deg]; each one is a cold start")
p.add_argument("--tu", type=float, default=None, help="free-stream turbulence intensity [%%] (default: the case value)")
args = p.parse_args()

comm = MPI.COMM_WORLD
case = AIRFOILS[args.airfoil]
tu = case["tu"] if args.tu is None else args.tu
alphas = case["alpha"] if args.alpha is None else args.alpha
grid = args.mesh or common.meshPath(args.airfoil, f"mesh_{args.level}", "volumeMesh.cgns")
outDir = args.outputDir or f"output_{args.airfoil}_{args.model}_{args.dissipation}"
recipe = args.recipe or recipes.recipe2D(args.model, args.dissipation,
                                         purpose="analysis" if args.airfoil == "naca0012" else "polar")

opts = recipes.options(args.model, recipe, args.dissipation, tu, grid, outDir, nuinf=args.nuinf, liftIndex=LIFT_INDEX)
overrides = common.applyOverrides(opts, args)

common.banner(comm, f"{args.airfoil} {args.level} | {recipes.MODEL_LABEL[args.model]} | {args.dissipation} dissipation | "
                    f"Tu {tu} % | nu_inf ratio {opts['eddyVisInfRatio']:g}\nrecipe {recipe}: {recipes.DESCRIPTION[recipe]}"
                    f"\noverrides {overrides}")

from baseclasses import AeroProblem  # noqa: E402

CFD = common.makeSolver(opts, comm)
CFD.addSlices("z", [SLICE_Z], sliceType="absolute")     # cp(x), cf(x): transition location

records = []
for alpha in alphas:
    kw = dict(name=f"{args.airfoil}_a{alpha:g}", alpha=alpha, mach=case["mach"], areaRef=AREA_REF,
              chordRef=CHORD_REF, evalFuncs=["cl", "cd", "cdp", "cdv"])
    if "altitude" in case:
        kw["altitude"] = case["altitude"]
    else:
        kw.update(reynolds=case["reynolds"], reynoldsLength=CHORD_REF, T=case["T"])
    ap = AeroProblem(**kw)
    records.append(common.solve(CFD, ap, comm, coldStart=True))

common.writeResult(comm, opts, dict(script="run_airfoil.py", airfoil=args.airfoil, level=args.level, model=args.model,
                                    dissipation=args.dissipation, tuPercent=tu, recipe=recipe,
                                    mach=case["mach"], reynolds=float(ap.re * CHORD_REF)),
                   records, common.provenance(comm, grid), overrides)
