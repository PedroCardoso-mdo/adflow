#!/usr/bin/env python
"""TU Braunschweig sickle wing on the paper authors' medium grid, GR with crossflow (paper section V.D).

What it demonstrates
    A geometrically swept 3D wing where crossflow moves the transition front at the kinks.
    On the authors' own grid, GR + crossflow + matrix dissipation gives the paper's forces.

Mesh
    showcase/meshes/sickle_wing/sickle_paper_L0.cgns : 4 404 960 cells, 18 blocks (merged from the
    authors' 592-block PLOT3D grid), O-grid, far field ~50 c. The grid was sent privately by the
    authors: check with them before redistributing it.

Conditions (run v14)
    M 0.16, altitude 3270.3675 m (Re 2.75e6 on chord 1.0 m), alpha -2.6 deg, Tu 0.17 %,
    roughness 1.2225e-5 m (the paper's 9.78e-6 m on its 0.8 m chord), eddyVisInfRatio 5e-7,
    areaRef 1.3072 (measured half-wing planform), liftIndex 3.

Reference numbers (job script job_sickle_paper_v14_cfon.sh, binary 36b01134, 64 ranks on one node)
    CL -0.066371      (paper Fig. 13b: -0.0640 .. -0.0660)
    CD  64.72 counts  (cdp 16.89, cdv 47.83; paper 60 .. 63)
    Convergence: 6.25 orders on totalR after 6426 iterations and 17.7 h; the 1e-8 target was NOT
    reached (res/r0 3.4e-8 on density). The residual sat flat for ~15 h while the crossflow front
    settled, then fell 20x in 80 iterations.
    Status: worth showing, not fully converged; the per-station x_tr comparison with the paper
    is qualitative only.

Cluster cost: 64 ranks, one 128-core node (memory-bound: ~100 GB), 17.7 h.

Commands
    mpirun -np 64 python run_sickle_wing.py
    mpirun -np 6 --bind-to core python run_sickle_wing.py --nCycles 4       # smoke (needs ~60 GB)
"""
from mpi4py import MPI

import common
import recipes

p = common.baseParser(__doc__, models=("gr",))
p.add_argument("--alpha", type=float, default=-2.6)
p.add_argument("--tu", type=float, default=0.17, help="turbulence intensity [%%]")
p.add_argument("--roughness", type=float, default=1.2225e-5, help="transitionRoughnessHeight [m]")
p.add_argument("--crossflow", type=int, default=1, choices=[0, 1])
args = p.parse_args()

comm = MPI.COMM_WORLD
grid = args.mesh or common.meshPath("sickle_wing", "sickle_paper_L0.cgns")
outDir = args.outputDir or "output_sickle_wing"
recipe = args.recipe or "gr_3d_crossflow_sickle"

opts = recipes.options(
    "gr", recipe, "matrix", args.tu, grid, outDir,
    nuinf=5e-7 if args.nuinf is None else args.nuinf,      # v14 value (the swept wing used 2.79e-7)
    liftIndex=3,
    transitionCrossflow=bool(args.crossflow),
    transitionRoughnessHeight=args.roughness,
    surfaceVariables=["cp", "mach", "cf", "yplus", "hcf", "rescf", "dscf"],
    volumeVariables=["intermittency", "fonset", "dscf", "rescf", "hcf"],
)
overrides = common.applyOverrides(opts, args)

common.banner(comm, f"Sickle wing, authors' grid | GR, crossflow {'ON' if args.crossflow else 'OFF'} | matrix dissipation\n"
                    f"recipe {recipe}: {recipes.DESCRIPTION[recipe]}\noverrides {overrides}")

from baseclasses import AeroProblem  # noqa: E402

CFD = common.makeSolver(opts, comm)
ap = AeroProblem(name="sickle", alpha=args.alpha, mach=0.16, altitude=3270.3675,
                 areaRef=1.3072, chordRef=1.0, evalFuncs=["cl", "cd", "cdp", "cdv"])
records = [common.solve(CFD, ap, comm)]
common.writeResult(comm, opts, dict(script="run_sickle_wing.py", recipe=recipe, crossflow=bool(args.crossflow),
                                    reynolds=float(ap.re * 1.0)),
                   records, common.provenance(comm, grid), overrides)
