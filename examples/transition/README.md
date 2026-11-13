# Transition-model examples

Runnable examples for the three laminar–turbulent transition models of this ADflow fork, each with
the solver recipe that produced the validated result.

| Key      | Model                                                | Selected by                                                                           |
| -------- | ---------------------------------------------------- | ------------------------------------------------------------------------------------- |
| `gr`     | γ–R̅eθt of Piotrowski & Zingg (2020), the main model | `turbulenceModel="SA-noft2-Gamma-Retheta"`; crossflow with `transitionCrossflow=True` |
| `sgamma` | one-equation γ                                       | `turbulenceModel="SA-noft2-Gamma"`                                                    |
| `bcm`    | SA-BCM, smooth form                                  | `turbulenceModel="SA"`, `useSABCM=True`                                               |
| `bcmo`   | SA-BCM, original form                                | as `bcm` with `SABCMSmooth=False`                                                     |
| `sa`     | plain SA, fully turbulent reference                  | `turbulenceModel="SA"`                                                                |

## Files

| File                                       | What it is                                                                       |
| ------------------------------------------ | -------------------------------------------------------------------------------- |
| `recipes.py`                               | every validated option set as an explicit dict, each with the job that proved it |
| `common.py`                                | command line, solve, `RESULT` line, `result.json`                                |
| `run_airfoil.py`                           | 2D analysis, any model (NACA 0012, NLF(1)-0416, S809; one alpha or a polar)      |
| `run_swept_nlf0415.py`                     | NLF(2)-0415 infinite swept wing with crossflow, matrix and scalar dissipation    |
| `run_sickle_wing.py`                       | sickle wing on the paper authors' grid, with crossflow                           |
| `run_ar5_wing.py`                          | AR5 wing analysis: GR on the refinement family, SA-BCM on the article mesh       |
| `opt_airfoil.py`                           | 2D shape optimisation, single-point and multipoint (FFD, lift by alpha, SLSQP)   |
| `opt_ar5_wing.py`                          | 3D wing optimisation set-ups: the GR attempt (c376) and the SA-BCM article case  |
| `ffd/`                                     | the three FFD boxes (small text files)                                           |
| `job_template.sh`                          | Slurm template: one 128-core node, one `srun` per case                           |
| `smoke_all.sh`                             | two-iteration smoke test of every example on a workstation                       |
| `smoke_compare.py`, `smoke_reference.json` | comparison of the smoke runs with the first iterations of the cluster runs       |

Each script starts with a header: what the case shows, the reference numbers with the convergence
level actually reached, the cluster cost and the exact command. `python <script> --help` prints it.

Meshes are not in git. The scripts look for them in `<repo>/showcase/meshes` (or
`$TRANSITION_MESH_ROOT`); `--mesh` takes any other file.

## The examples

cd in drag counts (1e-4) unless written as a coefficient. "1e-8" = total residual reduced by 8 orders
from free stream. Full tables, figures and dead ends: `showcase/_manifest_2d.md`, `showcase/_manifest_3d.md`.

| Example                                     | Case                                                        | Model(s)                          | Mesh (under `showcase/meshes/`), cells                                         | Cluster cost                      | Expected numbers                                                         | Status                                                            |
| ------------------------------------------- | ----------------------------------------------------------- | --------------------------------- | ------------------------------------------------------------------------------ | --------------------------------- | ------------------------------------------------------------------------ | ----------------------------------------------------------------- |
| `run_airfoil.py`                            | NACA 0012, α 2.7°, Tu 0.5 %, scalar                         | GR / SA-γ / SA-BCM                | `naca0012/mesh_L0/volumeMesh.cgns`, 193 224                                    | 32 ranks, 5–8 min                 | cd 62.19 / 67.13 / 66.01 counts; cl 0.3011 / 0.3024 / 0.3042             | validated, 1e-8 (jobs 1962375–79)                                 |
| `run_airfoil.py`                            | NACA 0012, α 4.52°, Tu 0.3 % (article point A), scalar      | GR / SA-γ / SA-BCM orig. / SA     | same                                                                           | 8 ranks, ~30 min                  | cd 65.06 / 68.83 / 65.76 / 102.02; cl 0.500                              | validated, 1e-8 (job 2197279)                                     |
| `run_airfoil.py`                            | NLF(1)-0416, α 4°, Re 4e6, Tu 0.15 %                        | GR scalar / GR matrix             | `nlf0416/mesh_L0/volumeMesh.cgns`, 193 224                                     | 32 ranks 10 min / 10 ranks 47 min | cd 65.58 / 65.22 (paper 65.40, experiment 64.36); cl 0.963 / 0.961       | validated, 1e-8 (jobs 1961835, 2201491)                           |
| `run_airfoil.py`                            | NLF(1)-0416, α 4°                                           | SA-γ / SA-BCM / SA-BCM orig. / SA | same                                                                           | 8–32 ranks, 7–37 min              | cd 69.28 / 64.87 / 64.73 / 112.76                                        | validated, 1e-8 (jobs 1961971, 1961959, 2168686, 1962367)         |
| `run_airfoil.py`                            | S809, α 4°, Re 2e6, Tu 0.07 %                               | GR scalar / GR matrix             | `s809/mesh_L0/volumeMesh.cgns`, 193 224                                        | 32 ranks 15 min / 14 ranks 70 min | cd 66.42 / 63.86; cl 0.655 / 0.656                                       | validated, 1e-8 (jobs 1961835, 2201434)                           |
| `run_airfoil.py --alpha ...`                | NLF0416 (α −6…6) and S809 (α 0…9) polars                    | all                               | same L0 meshes                                                                 | one point per 10–32 ranks         | tables in `showcase/_manifest_2d.md` (24-B, 24-C, 12-A)                  | GR and SA-BCM 23/23 points at 1e-8; SA-γ 16–17/23                 |
| `run_swept_nlf0415.py`                      | NLF(2)-0415 infinite swept wing, crossflow, 7 Re            | GR, matrix                        | `nlf0415_swept/mesh_L0/volumeMesh_periodic.cgns`, 64 800                       | 9 ranks per Re, 36–72 min         | x_tr 0.680 … 0.266 (paper 0.710 … 0.266); Re 2e6: cl −0.10154, cd 72.85  | validated, 1e-8 (job 1826778)                                     |
| `run_swept_nlf0415.py --dissipation scalar` | same wing, scalar dissipation                               | GR, scalar                        | same                                                                           | 9 ranks, 7–14 min                 | x_tr ≈ 0.83 at every Re: no crossflow front                              | negative result (counter-example)                                 |
| `run_sickle_wing.py`                        | Sickle wing, authors' grid, crossflow                       | GR, matrix                        | `sickle_wing/sickle_paper_L0.cgns`, 4 404 960                                  | 64 ranks, 17.7 h, ~100 GB         | CL −0.0664 (paper −0.064…−0.066); CD 64.7 counts (paper 60…63)           | worth showing; 6.25 orders, 1e-8 not reached                      |
| `run_ar5_wing.py --model gr`                | AR5 wing + winglet, α 1.8°, Tu 0.25 %                       | GR, scalar (ADflow default)       | `ar5_wing/corrected_foil_ar5_vol_L{4,3,2,1}.cgns`, 0.46 / 0.90 / 1.78 / 3.67 M | 64 ranks; L1 1 h 34 min           | L1: cl 0.95822, cd 0.035907 (cdp 0.032626, cdv 0.003281)                 | converged 1e-8 on L4–L1; not mesh-converged, front under-resolved |
| `run_ar5_wing.py --model bcm`               | AR5 wing of the SA-BCM article, α 1.0577°                   | SA-BCM (also original, SA)        | `ar5_wing/ar5_asa_L2_943104cells.cgns`, 943 104                                | 64 ranks, 5 min                   | cl 0.894395, cd 0.02445674 (cdp 0.02122303, cdv 0.00323371)              | reproducible; converged to 1e-6 (job 1904793)                     |
| `opt_airfoil.py`                            | NACA 0012 single-point optimisation at A (Tu 0.3 %, cl 0.5) | GR scalar / GR matrix             | `naca0012/mesh_L0/volumeMesh.cgns`, 193 224                                    | 64 ranks, 24–48 h                 | cd 65.1 → 31.6 (−51.5 %) / 66.0 → 30.7 (−53.6 %)                         | final (stall rule)                                                |
| `opt_airfoil.py`                            | same, other models                                          | SA-γ, SA-BCM orig.                | same                                                                           | 64 ranks, 48 h                    | scalar 68.4 → 32.3, 65.5 → 31.3; matrix 72.0 → 46.5, 68.3 → 48.8         | provisional                                                       |
| `opt_airfoil.py --points a,b,c`             | multipoint 0.5:0.3, 1.0:0.3, 0.5:0.4                        | GR scalar                         | same                                                                           | 3 × 42 ranks, 48 h                | −41.5 / −40.7 / −43.4 % per point                                        | final (job 1980228)                                               |
| `opt_airfoil.py --points a,b,c`             | multipoint 0.3:0.5, 0.6:0.5, 0.3:0.7, matrix                | GR matrix                         | same                                                                           | 3 × 42 ranks, 48 h                | cd 47.34 at Tu 0.3 %, cl 0.5 on the multipoint shape with the new recipe | open (running)                                                    |
| `opt_ar5_wing.py --case gr_c376`            | AR5 wing, 182 shape variables, cl 0.8932                    | GR scalar (SA control)            | `ar5_wing/chord_ladder_c376/corrected_foil_ar5_vol_c376.cgns`, 13 684 496      | 4 nodes × 32 ranks, 3 × 48 h      | cd 287.1 → 269.5 counts (−6.1 %), cl 0.888; SA 362.8 → 345.2 (−4.9 %)    | attempt: lift constraint not closed, unfinished                   |
| `opt_ar5_wing.py --case bcm_article`        | AR5 wing of the SA-BCM article, cl ≥ 0.8932                 | SA-BCM (SA)                       | `ar5_wing/ar5_asa_L2_943104cells.cgns`, 943 104                                | 64 ranks, 12.3 h                  | cd 0.024457 → 0.024053 (−1.65 %); SA −3.40 %                             | reproducible; article claims withdrawn (see below)                |

Not runnable on a workstation: the c376 mesh (13.7 M cells) and the AR5 L0 mesh (7.4 M cells, cluster only).

### Things to know before quoting a number

- **Free-stream ν̃ is a boundary condition of each model** (`recipes.NUINF`): GR 2.79e-7, SA-γ 1e-7,
  SA-BCM 1e-9, SA 9e-3. With the ADflow default 9e-3 a transition model loses its laminar region
  (NLF0416 α 4°, GR: +17 % cd). The AR5 refinement study still ran at 9e-3 and the sickle wing at 5e-7.
- **Crossflow needs matrix dissipation.** With scalar dissipation the swept wing converges and shows no
  crossflow front.
- **AR5 wing (GR):** the solver converges, but the transition front of this 181-wrap mesh family is a
  mesh artefact (x/c 0.33 against 0.62 from 397 wrap points up). The 3D GR optimisation is an attempt.
- **SA-BCM article:** the numbers reproduce, but the article's central claims were withdrawn after
  later tests (the original model converges with the same recipe; its derivative errors came from a
  stalled primal). Do not present the smooth form as more reliable than the original.
- **SA-γ**: the final form is the first onset form (smooth min/max, P&Z F_turb) with artificial viscosity on γ (`sgammaArtVisc`).
  The stored runs do not all use it: the scalar optima ran with the tanh onset and no artificial viscosity,
  the matrix optimisations with coefficient 5, the mesh study and derivative checks with 10.
- Optima of SA-γ and SA-BCM must be re-analysed with GR before they are trusted (the gain often does
  not survive).

## Recipes

`recipes.py` is the single place for solver settings. `recipes.options(model, recipe, dissipation,
tuPercent, gridFile, outputDir)` returns the complete ADflow option dict.

| Recipe                           | Ladder                                                                    | Proved by                                                |
| -------------------------------- | ------------------------------------------------------------------------- | -------------------------------------------------------- |
| `gr_2d_scalar`                   | ANK → CANK 1e-5 → CSANK 1e-6, no NK                                       | campaign 24, jobs 1962375–79, 1988112                    |
| `gr_2d_scalar_polar`             | ANK → CANK 1e-2 → CSANK 1e-5, no NK                                       | 24-B, jobs 1961834–39                                    |
| `gr_2d_matrix`                   | ANK → CANK 1e-2 → CSANK 1e-4, no NK                                       | 12-A jobs 1827893 / 1828015; 24-C jobs 2201434 / 2201491 |
| `gr_2d_matrix_opt`               | second order 1e-4, coupling 1e-5, `ANKMaxIter` 80, no NK                  | jobs 2219449, 2219697 (2026-11-12)                       |
| `gr_3d_crossflow_swept`          | ANK → CANK 1e-2 → CSANK 1e-5 → NK 1e-6, `ANKCFLLimit` 1e8                 | job 1826778                                              |
| `gr_3d_crossflow_swept_scalar`   | ANK → CANK 1e-3 → CSANK 1e-5 → NK 1e-6                                    | `job_swept_L0_resweep_ladder.sh`                         |
| `gr_3d_crossflow_sickle`         | ANK → CANK 1e-3 → CSANK 1e-6, no NK, per-node damping, `ANKCFLMinCap` 1e3 | `job_sickle_paper_v14_cfon.sh`                           |
| `gr_3d_wing_ar5`                 | ANK → CANK 1e-4 → CSANK 3e-5 → NK 1e-6, then restart without NK           | jobs 1810888–92, 1810916–25                              |
| `gr_3d_wing_ar5_opt`             | ANK → CANK 1e-2 → CSANK 1e-3, no NK, 1e-7; adjoint 1e-3                   | jobs 1844851, 1931137, 1939946                           |
| `sgamma_2d`                      | ANK → second order 1e-2, never coupled, no NK                             | jobs 1962375–79, 1961969–74                              |
| `sgamma_2d_tuned`                | coupling 1e-5, second order 1e-4                                          | jobs 2197279, 2199671, 2203553                           |
| `bcm_2d` (= `sa_2d`)             | ANK → second order 1e-3, no NK, 3 turbulence sub-iterations               | jobs 1962377, 1961957–62, 2168686, 1962365–70            |
| `bcm_2d_tuned` (= `sa_2d_tuned`) | `bcm_2d` + VLR, line-search tolerances, CFL 1e8, second order 1e-4        | jobs 2198926, 2203553                                    |
| `bcm_3d_wing`                    | ANK from CFL 1 → second order 1e-3, target 1e-6                           | jobs 1904652, 1904793                                    |
| `sa_3d_wing`                     | as `bcm_3d_wing`, CFL limit 3e5                                           | article `08_base_sa`, `07_c1_sa`                         |

ANK = approximate Newton–Krylov with segregated turbulence; CANK = flow and turbulence coupled
(`ANKCoupledSwitchTol`); CSANK = coupled with the exact matrix-free Jacobian (`ANKSecondOrdSwitchTol`);
NK = Newton–Krylov (`NKSwitchTol`). All switch levels are relative to the free-stream residual.

Setting found on 2026-11-12 for GR optimisation solves with matrix dissipation:
`ANKSecondOrdSwitchTol = 1e-4` (recipe `gr_2d_matrix_opt`). It switches the exact Jacobian on above the
1–2e-5 level where the laminar-bubble stall sits; the multipoint shape that never converged from free
stream reaches 1e-8 (cd 47.34 counts at Tu 0.3 %, cl 0.5; jobs 2219449, 2219697).

For the transition models pyADflow sets by itself, unless the user passes them:
`globalPreconditioner="field split"`, `ILUFill=3`, `ASMOverlap=3`, `adjointSolver="LGMRES"`,
`adjointSubspaceSize=400`, `NKLSRelax=True`. They are written out in `recipes.py`
(`TRANSITION_PYADFLOW_DEFAULTS`). The per-node damping of the sickle recipe is selected by
`ANKTransitionGlobalLambda=False` alone.

### Where the sources disagreed

| Item                         | Source A                                                       | Source B                                                                       | Used here                                      |
| ---------------------------- | -------------------------------------------------------------- | ------------------------------------------------------------------------------ | ---------------------------------------------- |
| Sickle-wing ladder           | earlier project notes: CANK 1e-2 → CSANK 1e-4                   | job script v14: CANK 1e-3 → CSANK 1e-6 + per-node damping + `ANKCFLMinCap` 1e3 | the job script (the only finished run)         |
| Campaign-12 `*_final` series | `14_nlf_2d_solver_ab/PURPOSE.md`: scalar vis4 0.01, ν̃∞ 1.5e-3 | job scripts: matrix vis4 0.04, ν̃∞ 5e-7                                        | the job scripts                                |
| Matrix coefficients          | `12_paper_conditions/PURPOSE.md`: Vl = Vn = 0                  | job headers: Vl 0 (`epsShear`), Vn 0.25 (`epsAcoustic`)                        | Vl 0, Vn 0.25                                  |
| GR polar recipe, scalar L0   | `recipes.py` "polar" = "opt" (CANK 1e-5)                       | `result.json` of jobs 1961834–39: "polar_old" (CANK 1e-2 → CSANK 1e-5), NK off | what `result.json` records                     |
| c376 optimisation ladder     | `vars_master.py` defaults 1e-4 / 3e-5 (job 1830863)            | option dump of jobs 1844851, 1931137, 1939946: 1e-2 / 1e-3                     | 1e-2 / 1e-3 (the later jobs)                   |
| SA-γ ν̃∞                     | `setup.py`: 1e-7 (NACA runs)                                   | L0 polars of the figures: 2.79e-7                                              | default 1e-7; `--nuinf 2.79e-7` for the polars |
| Swept-wing Reynolds length   | `13_nlf0415_swept/PURPOSE.md`: √2                              | final job 1826778: 1.0 (normal chord)                                          | 1.0                                            |
| SA-BCM article adjoint PC    | manifest: fill not listed                                      | `run.py`: `ILUFill` 4                                                          | `ILUFill` 4                                    |
| SA on the article wing       | SA value ν̃∞ 9e-3                                              | `08_base_sa/run.py`: 1e-10, CFL limit 3e5                                      | what the script ran                            |

## How to run on a cluster

`job_template.sh` is a Slurm script for one 128-core node with several cases side by side:

```bash
cd examples/transition
# edit the three "environment" lines and the account in job_template.sh
mkdir -p logs && sbatch job_template.sh
```

- One `srun` per case, each in its own subshell with its own `timeout -s USR2`.
- ADflow traps SIGUSR2: it writes the solution and leaves the solve, and the script still writes
  `result.json`. Stop a job with `scancel --signal=USR2 <jobid>`, not with a plain `scancel`.
- `OMP_NUM_THREADS=1`. Fill the node (it is billed whole), but for 3D cases prefer more ranks per
  case to more cases: ADflow is memory-bound.
- Never set `writeSolutionEachIter=True`.
- Every log must start with the ADflow banner (`ADflow … | git <hash> (<branch>) | installed …`).

Memory guide (GR, 64 ranks): AR5 L1 (3.67 M cells) ≈ 105 GB; sickle wing (4.4 M) ≈ 100 GB without NK;
AR5 L0 (7.4 M) needs 2 nodes; the c376 optimisation 4 nodes × 32 ranks.

## Two-iteration smoke test on a workstation

```bash
source <python environment>            # mpi4py, baseclasses, pygeo, idwarp, pyoptsparse, multipoint
export OMP_NUM_THREADS=1 PYTHONPATH=<repo>      # only if this build is not the installed adflow
cd examples/transition
./smoke_all.sh                         # all cases, 6 ranks, one at a time
SKIP_LARGE=1 ./smoke_all.sh            # without the sickle wing (4.4 M cells, ~60 GB) and AR5 L1
./smoke_all.sh swept                   # only the cases whose name contains "swept"
```

A single case by hand:

```bash
mpirun -np 6 --bind-to core python run_airfoil.py --airfoil nlf0416 --model gr --alpha 4 --nCycles 4
mpirun -np 6 --bind-to core python opt_airfoil.py --model gr --points 0.3:0.5 --baselineOnly --nCycles 4
```

`--nCycles 4` stops after the first two solver iterations (`nCycles` counts "Iter Tot"). The
optimisation scripts take `--baselineOnly`: the problem is built and the start design evaluated, with
no optimiser iteration and no adjoint. `smoke_all.sh` ends with `smoke_compare.py`, which writes
`showcase/smoke_runs/SMOKE_REPORT.md`: the residual of the free-stream state (iteration 0) must agree
with the cluster run of the same case to round-off; iterations 1–2 agree to a few digits only, because
the preconditioner depends on the rank count.

Last run (2026-11-12, 6 ranks, build stamp 2f3d507d): 30 cases ran, all 29 with a cluster twin agree
at iteration 0 to 1e-12 or better; 2D cases take 7–14 s, AR5 L1 2 min, the sickle wing 2 min 15 s
(about 48 GB). Details and caveats: `showcase/smoke_runs/SMOKE_REPORT.md`.

What a smoke run does not test: optimiser iterations, adjoints, restarts, the switch levels of the
ladders and the NK phase.
