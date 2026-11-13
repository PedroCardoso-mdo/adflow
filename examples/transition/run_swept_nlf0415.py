#!/usr/bin/env python
"""NLF(2)-0415 infinite swept wing: crossflow transition with the GR model (paper section V.C).

What it demonstrates
    transitionCrossflow = True reproduces the transition location against Reynolds number of
    Piotrowski & Zingg (2020) -- with MATRIX dissipation. With scalar dissipation the same case
    converges but shows no crossflow front (x_tr ~0.83 at every Re): --dissipation scalar is the
    counter-example.

Mesh
    showcase/meshes/nlf0415_swept/mesh_L0/volumeMesh_periodic.cgns : 64 800 cells (541 x 121 nodes,
    one span cell, translational-periodic span faces). L1 129 200, L2 259 200, L3 516 800 cells.
    The geometry is unswept; the 45 deg sweep is the flow angle beta.

Conditions (job 1826778)
    M 0.10, T 288.15 K, alpha -4 deg, beta 45 deg, Tu 0.20 %, roughness 3.3e-6 m,
    Re 1.8e6 .. 3.5e6 on the NORMAL chord (reynoldsLength 1.0), eddyVisInfRatio 2.79e-7.

Reference numbers, matrix dissipation, L0 (job 1826778, 9 ranks per case, all at res/r0 <= 2e-9)
    Re/1e6   x_tr    paper   cl         cd (cdp, cdv) [counts]    iterations   wall
    1.796    0.680   0.710   -0.09793   71.08 (17.63, 53.45)      591          4221 s
    2.000    0.604   0.623   -0.10154   72.85 (17.76, 55.09)      615          4201 s
    2.204    0.531   0.541   -0.10399   73.98 (17.74, 56.24)      634          4290 s
    2.370    0.484   0.486   -0.10552   74.69 (17.66, 57.03)      627          3648 s
    2.498    0.447   0.444   -0.10658   75.24 (17.57, 57.67)      633          3276 s
    3.000    0.339   0.349   -0.10888   76.06 (16.96, 59.11)      495          2173 s
    3.498    0.266   0.266   -0.10984   75.76 (15.91, 59.84)      582          4304 s
    Mean |x_tr error| 1.9 % of chord against the paper. Status: validated.
    Scalar counter-example: converged ~8 orders, x_tr ~0.83-0.84, no Re trend. Those runs used the
    Reynolds number on the streamwise chord (--reynoldsLength 1.41421356). Status: negative result.

Cluster cost: 7 Reynolds numbers x 9 ranks on one node, 36-72 min each.

Commands
    mpirun -np 9 python run_swept_nlf0415.py --reynolds 2.0e6
    mpirun -np 9 python run_swept_nlf0415.py --reynolds 2.0e6 --dissipation scalar
    mpirun -np 6 --bind-to core python run_swept_nlf0415.py --reynolds 2.0e6 --nCycles 4      # smoke
"""
from mpi4py import MPI

import common
import recipes

p = common.baseParser(__doc__, models=("gr",))
p.add_argument("--level", default="L0")
p.add_argument("--dissipation", choices=["matrix", "scalar"], default="matrix")
p.add_argument("--reynolds", type=float, default=2.0e6)
p.add_argument("--reynoldsLength", type=float, default=1.0, help="1.0 = normal chord (the paper's definition)")
p.add_argument("--alpha", type=float, default=-4.0)
p.add_argument("--beta", type=float, default=45.0, help="sweep, realised by the flow angle")
p.add_argument("--tu", type=float, default=0.20, help="turbulence intensity [%%]")
p.add_argument("--roughness", type=float, default=3.3e-6, help="transitionRoughnessHeight [m]")
p.add_argument("--crossflow", type=int, default=1, choices=[0, 1])
args = p.parse_args()

comm = MPI.COMM_WORLD
grid = args.mesh or common.meshPath("nlf0415_swept", f"mesh_{args.level}", "volumeMesh_periodic.cgns")
outDir = args.outputDir or f"output_swept_nlf0415_{args.dissipation}"
recipe = args.recipe or ("gr_3d_crossflow_swept" if args.dissipation == "matrix" else "gr_3d_crossflow_swept_scalar")

opts = recipes.options(
    "gr", recipe, args.dissipation, args.tu, grid, outDir, nuinf=args.nuinf,
    liftIndex=2,
    transitionCrossflow=bool(args.crossflow),
    transitionRoughnessHeight=args.roughness,
    # crossflow diagnostics on the surface and in the volume (helicity criterion hcf against 0.1066)
    surfaceVariables=["cp", "mach", "cf", "yplus", "hcf", "rescf", "dscf"],
    volumeVariables=["intermittency", "fonset", "dscf", "rescf", "hcf"],
)
overrides = common.applyOverrides(opts, args)

common.banner(comm, f"NLF(2)-0415 swept {args.level} | GR, crossflow {'ON' if args.crossflow else 'OFF'} | "
                    f"{args.dissipation} dissipation | Re {args.reynolds:.4g} on L = {args.reynoldsLength}\n"
                    f"recipe {recipe}: {recipes.DESCRIPTION[recipe]}\noverrides {overrides}")

from baseclasses import AeroProblem  # noqa: E402

CFD = common.makeSolver(opts, comm)
CFD.addSlices("z", [0.05], sliceType="absolute")

ap = AeroProblem(
    name=f"nlf0415_re{args.reynolds / 1e6:.3f}", alpha=args.alpha, beta=args.beta, mach=0.10,
    reynolds=args.reynolds, reynoldsLength=args.reynoldsLength, T=288.15,
    areaRef=0.1, chordRef=1.0, evalFuncs=["cl", "cd", "cdp", "cdv"],
)
records = [common.solve(CFD, ap, comm)]
common.writeResult(comm, opts, dict(script="run_swept_nlf0415.py", level=args.level, dissipation=args.dissipation,
                                    recipe=recipe, reynolds=args.reynolds, reynoldsLength=args.reynoldsLength,
                                    crossflow=bool(args.crossflow)),
                   records, common.provenance(comm, grid), overrides)
