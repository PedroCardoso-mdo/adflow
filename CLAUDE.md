# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

# ADflow transition models (branch `transition-models`)

## Models on this branch

| Model                    | Selected by                                        | State (nwt)    | Code                                |
|--------------------------|----------------------------------------------------|----------------|-------------------------------------|
| SA-GR (P&Z 2020 γ-Re̅θt)  | `turbulenceModel = "SA-noft2-Gamma-Retheta"`       | ν̃, γ, Re̅θt (3) | `src/turbulence/saGammaRetheta.F90` |
| SA-sγ (one-eq. γ)        | `turbulenceModel = "SA-noft2-Gamma"` (id 9)        | ν̃, γ (2)       | `src/turbulence/saGamma.F90`        |
| SA-BCM (algebraic γ)     | `use_SABCM = True` on plain SA (`SABCM_*` options) | ν̃ (1)          | `src/turbulence/sa.F90` (protected) |

**Status: finished.** SA-GR converges, is validated against the paper, and
passes the full derivative ladder (dot-product, `_b` vs `_fast_b`, AD vs
complex step, total adjoint vs CS) on the standard ADflow tutorial-wing mesh.
SA-GR is the main model and the only one used for optimisation. SA-sγ design
docs live outside the repo in `/home/mdo/Desktop/Run/MDO_PhD/Transition/sa_sgamma/`.

Solver paths, all runtime-selectable: DADI (`TurbDADICoupled` full / decoupled
/ transition-only), turb-ANK, and coupled ANK (flow + turbulence in one
Newton-Krylov system).

## Hard rules

1. **Do not modify the SA model for GR / SA-sγ work.** Transition is a
   modifier — γ multiplies SA production only (Eq. 41); the coupling lives in
   `saGammaRetheta.F90` / `saGamma.F90`. SA-BCM itself lives in `sa.F90`, so
   the `guard_protected_files.sh` hook asks for approval on every `sa.F90` edit.
2. **Crossflow is OFF by default** (`transitionCrossflow = False`): with it on,
   the tutorial wing plateaus at ~3e-2. Do not flip the default without explicit
   user instruction. Multigrid is not supported for transition.
3. **All transition diagnostics go to the volume CGNS** via the
   `transitionDebug` array. No ASCII debug files, no per-cell printouts.
4. **First-order upwind** for γ and Re̅θt convection (paper §IV.A).
5. **Never hand-edit Tapenade output** (`src/adjoint/{outputForward,
   outputReverse,outputReverseFast,temp_*}/`; `settings.json` prompts on it).
   Regenerate with Tapenade. Hand-written adjoint code (`adjoint*.F90`,
   `masterRoutines.F90`, `fortranPC.F90`, `autoEdit/*.py`, `Makefile_tapenade`)
   is fine to edit.
6. **Every run log must identify its ADflow revision.** The root proc prints
   `| ADflow … | git <hash> (<branch>) | installed <timestamp>` (+ `*** DIRTY
   WORKING TREE ***`), stamped into the gitignored `adflow/gitInfo.py` at
   `pip install` time. A log without the banner is not a valid result —
   reinstall and rerun. Never silence it.
7. **Paper wins** when paper and code disagree on physics. Source of truth:
   `docs/SA_GAMMA_RETHETHA_BASE/SAGR_01_paper_piotrowski_zingg_2020.md`
   (SA-BCM: `docs/papers/`).
8. **Units are p-ρ non-dimensional.** Velocity normalizes to M·√γ (not 1),
   viscosities are ratios to μ_∞ (`rlv`, `rev`), `1/Re` is NOT absorbed into the
   viscosity. Read `docs/ADFLOW_BASE/ADFLOW_08_nondimensionalization.md` before
   touching any equation with velocity, viscosity, time scale or Re.
9. **Token discipline:** load only the files `docs/README.md` routes you to.

Be concise in replies — in wording, not in effort.

## Docs

`docs/README.md` is the index and routing table. Everything else in `docs/`
is one of: `CORE_BASE/CORE_01_architecture.md` (code map, all options,
design rationale, known limitations, SA-sγ and SA-BCM parts),
`CORE_BASE/CORE_02_convergence_strategy.md` (recipes per model/case),
`ADFLOW_BASE/ADFLOW_08/09` (units; AD wiring and maintenance),
`VERIFICATION/VERIF_00/02/06` (derivative ladder; rotating frame; solver
fidelity F-tags cited from code), and reference texts (paper, `papers/`,
`MICHAEL_PIOTROWSKI/`). Upstream ADflow docs are in `doc/*.rst`.

## File locations (most-edited)

| What            | Where                                                                  |
|-----------------|------------------------------------------------------------------------|
| Residuals       | `src/turbulence/saGammaRetheta.F90` (GR), `src/turbulence/saGamma.F90` (SA-sγ) |
| Helpers         | `src/turbulence/turbUtils.F90`: `reThetaTCorrelation`, `flengthCorrelation`, `rethetacCorrelation`, `smoothMinMax` |
| Constants       | `src/modules/paramTurb.F90`                                            |
| Input params    | `src/modules/inputParam.F90` (+ `src/f2py/adflow.pyf` — a variable missing there is a silent no-op) |
| Solvers         | `src/NKSolver/NKSolvers.F90` (ANK / NK / turbKSP), `src/solver/blockette.F90` |
| BCs             | `src/turbulence/turbBCRoutines.F90`                                    |
| Init            | `src/initFlow/initializeFlow.F90`                                      |
| Output dispatch | `grep -rn 'case ("eddy")' src/`                                        |
| Python wrapper  | `adflow/pyADflow.py`                                                   |
| AD generated    | `src/adjoint/output{Forward,Reverse,ReverseFast}/`                     |

## AD-relevant code (triggers `TAPENADE NEEDED`)

Any change to differentiated residual math needs a Tapenade regen before the
adjoint is trusted: `saGammaRetheta.F90` (Source, Viscous, correlations),
`saGamma.F90` (`sgSource`, `sgViscous`), `sa.F90`, the helpers in
`turbUtils.F90`, and anything else in `src/adjoint/Makefile_tapenade`. Code
behind `#ifndef USE_TAPENADE` (LHS/qq blocks, DADI-only) does not trigger it.
Say `TAPENADE NEEDED` in the wrap-up; the user runs Tapenade, then `/build`.

## Build, test, run

- **Build:** `/build real|complex|both` (user skill). It removes `build/`
  before `pip install` — tests and runs import adflow from site-packages, so
  an un-reinstalled build is silently stale.
- **Derivative verification:** `adflow-derivative-verification` skill →
  `tests/reg_tests/run_{sagr,sgamma,bcm}_tests.sh [blockette|real|cs|adjoint|train|genw]`.
  Single test: `testflo <file>.py:<Class>.<test>` (`testflo --dryrun` lists them).
- **Every ADflow solve goes to Deucalion (sbatch), never the workstation** —
  smoke tests included. Skills: `deucalion-hpc`, `hpc-pack`. Run folders live in
  `/home/mdo/Desktop/Run/MDO_PhD/Transition/gama_rethetha/` (index: its
  `README.md`; HPC recipe: `HPC_HOWTO.md`), each with a `PURPOSE.md`.
- When a task is done, invoke the `wrapup` skill (proposes doc updates, asks
  before commit/push).
