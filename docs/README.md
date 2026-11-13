# ADflow with laminar–turbulent transition models

This fork of [ADflow](https://github.com/mdolab/adflow) adds three transition models to the
Spalart–Allmaras (SA) RANS solver, the solver and adjoint extensions they need, and the
Tapenade-differentiated code that makes them usable in gradient-based optimisation.
With plain SA the code behaves as upstream (checked, see [verification.md](verification.md)).

## The three models

| Short name     | How to select                                                      | Transport equations | Reference                                                                              | Status                                                                                                                                         |
| -------------- | ------------------------------------------------------------------ | ------------------- | -------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| GR             | `turbulenceModel="SA-noft2-Gamma-Retheta"`                         | 3 (ν̃, γ, R̅eθt)      | Piotrowski and Zingg, AIAA Journal, 2020 (doi:10.2514/1.J059784)                       | Validated in analysis (2D airfoils against the paper). Used in optimisation: 2D finished, 3D an unfinished attempt.                            |
| GR + crossflow | as GR, with `transitionCrossflow=True`                             | 3                   | same paper, Eqs. 15–26                                                                 | Validated on the NLF(2)-0415 infinite swept wing with matrix dissipation. Sickle wing: residual reduced 6.25 orders. Not used in optimisation. |
| SA-γ           | `turbulenceModel="SA-noft2-Gamma"`                                 | 2 (ν̃, γ)            | Menter et al. 2015 one-equation γ model, coupled to SA after Nichols 2019 and Lee 2021 | Final form: first onset form + `sgammaArtVisc`. Validated in 2D analysis; optimisation results mix settings (see validation.md) |
| SA-BCM         | `turbulenceModel="SA"`, `useSABCM=True` (`SABCMSmooth` True/False) | 1 (ν̃), algebraic γ  | Mura and Cakmakcioglu, AIAA 2020-2714                                                  | Validated in 2D analysis. Used in optimisation; its optima must be re-analysed with GR. Earlier claims about the smooth form were withdrawn.   |

GR is the main model and the reference against which the other two are judged.

## What else the fork adds

| Area    | Addition                                                                                                                                                                                                    |
| ------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Solver  | Coupled ANK for the transition unknowns, column scaling of the turbulence Krylov systems, block-coupled turbulence DADI, source-term time-step restriction, per-node update damping, relaxed NK line search |
| Scheme  | Run-time coefficients of the matrix-dissipation eigenvalue limiters (`epsAcoustic`, `epsShear`)                                                                                                             |
| Adjoint | Tapenade forward, reverse and fast-reverse code for both transported models and for SA-BCM; field-split preconditioner; LGMRES                                                                              |
| Output  | Transition diagnostics as volume and surface variables; a build stamp (git hash, branch, install time) printed at the top of every run                                                                      |

## Quick start

```bash
# 1. Build and install as upstream ADflow (two notes in developer.md: no OpenMP, Python version).
# 2. Meshes are not in git: put them under showcase/meshes or set TRANSITION_MESH_ROOT.
cd examples/transition
python run_airfoil.py --help         # case description, reference numbers, exact command
mpirun -np 8 python run_airfoil.py --airfoil nlf0416 --model gr --alpha 4
#    reference: cd 65.58 counts, cl 0.963, residual reduced 8 orders (32 ranks, 10 min)
mpirun -np 6 python run_airfoil.py --airfoil nlf0416 --model gr --alpha 4 --nCycles 4   # 2-iteration smoke run
./smoke_all.sh                       # every example for 2 iterations, compared with the reference logs
sbatch job_template.sh               # Slurm template: several cases on one node
```

In a script of your own, take the option dictionary from the recipes instead of writing it by hand:

```python
import recipes   # examples/transition/recipes.py
options = recipes.options("gr", "gr_2d_scalar", "scalar", tuPercent=0.15,
                          gridFile="mesh.cgns", outputDir="out")
```

Two settings matter more than any other: `eddyVisInfRatio` must be set per model (the ADflow default 0.009
removes the laminar region) and `ANKCFLLimit` must be raised to 1e8. Both are explained in [recipes.md](recipes.md).

## Pages

| Page                                           | Content                                                                                      |
| ---------------------------------------------- | -------------------------------------------------------------------------------------------- |
| [models.md](models.md)                         | Equations implemented, constants, boundary and free-stream conditions, model options, limits |
| [options.md](options.md)                       | Every option added over upstream: default, meaning, which model reads it                     |
| [recipes.md](recipes.md)                       | Validated solver settings per model and case type, and what fails without them               |
| [validation.md](validation.md)                 | What was computed, against what, with the numbers and a status label per case                |
| [verification.md](verification.md)             | Derivative checks, comparison with upstream, smoke tests, and the state of each check |
| [developer.md](developer.md)                   | Code map, non-dimensionalisation, AD wiring and Tapenade, differences from upstream, build   |
| [tried_and_rejected.md](tried_and_rejected.md) | Ideas and options that were tested and dropped, with the evidence                            |

Other material in the repository: `examples/transition/README.md` (runnable examples),
`showcase/_manifest_2d.md` and `showcase/_manifest_3d.md` (inventories of every campaign with file-level
provenance), `tests/reg_tests/README_transition.md` (derivative test suite), `doc/` (upstream ADflow documentation).

## Conventions used in these pages

- cd in drag counts (1 count = 1e-4), with the pressure / viscous split (cdp / cdv) where it was recorded.
- "Residual reduced N orders" refers to the total residual relative to the free-stream state.
- Campaign names (`12_paper_conditions`, `24_minmodel_2d`, …) are the folder names used in the two manifests.
  Cluster job numbers appear only as provenance.
- ANK = approximate Newton–Krylov with segregated turbulence. CANK = flow and turbulence in one Krylov system
  (below `ANKCoupledSwitchTol`). CSANK = coupled, with the second-order matrix-free Jacobian (below
  `ANKSecondOrdSwitchTol`). NK = Newton–Krylov (below `NKSwitchTol`). All switch levels are relative to the
  free-stream residual.
- Status labels: **validated** (compared with a reference and converged), **converged N orders** (no
  reference comparison, or target missed), **provisional** (will be redone), **attempt** (unfinished).
