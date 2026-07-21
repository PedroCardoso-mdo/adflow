# 00 — Inventory: SA-γ-Re̅θt (SA-sLM2015) Implementation Map

Date: 2026-07-07. Branch: `sa_gamma_rethetha`. This is a file map only — no analysis or conclusions.

## 1. Transition-model source files (primal)

| File | Role |
|------|------|
| `src/turbulence/saGammaRetheta.F90` | Main model module `saGammaReTheta`: `saGammaReTheta_block` driver, `Source`, `Viscous`, `ResScale`, `saGammaReThetaSolve` (DD-ADI), `evalSrcJacBlock` (3×3 source Jacobian), `computeSrcLambda`. SA-production coupling lives here (γ clamped to [0,1] multiplies SA production only; `gammaForSA` debug slot 30). |
| `src/turbulence/turbUtils.F90` | Correlation functions (`reThetaTCorrelation`, etc.) and `smoothMinMax` helpers (~lines 2279–2412), plus generic turb utilities. |
| `src/turbulence/turbAPI.F90` | Dispatch layer: `turbSolveDDADI` and residual entry points select SA / SST / SA-γ-Re̅θt (`saGammaRetheta_block`) by `turbModel`. |
| `src/turbulence/turbBCRoutines.F90` | Turbulence BCs incl. transition-variable wall / farfield ghost-cell treatment. |
| `src/turbulence/sa.F90` | Baseline SA model — **unmodified by design** (γ coupling is applied inside `saGammaRetheta.F90`, not here). |
| `src/turbulence/SST.F90` | Menter SST model (reference/template for 2-eq machinery, no transition coupling). |
| `src/bcdata/BCData.F90` | Farfield/inflow BC data arrays; extended for transition variables. |
| `src/initFlow/initializeFlow.F90` | `referenceState` sets `wInf(itu2/itu3)` (γ∞ = 1, Re̅θt∞ via `reThetaTCorrelation`); flow-field init for transition vars. |
| `src/initFlow/variableReading.F90` | Restart-file reading of the transition variables. |
| `src/modules/paramTurb.F90` | Model constants (LM2015 / SA-sLM2015 correlation constants). |
| `src/modules/inputParam.F90` | Runtime options for the transition model (solver-path selection, etc.). |
| `src/inputParam/inputParamRoutines.F90` | Option parsing / validation for the transition options. |
| `src/modules/constants.F90` | State-vector indices: `itu1` (SA ν̃), `itu2`, `itu3` (=8), `itu4` (=9) — γ and Re̅θt slots; `spalartAllmarasNoft2GammaRetheta` model enum. |
| `src/modules/block.F90` | Block-level `transitionDebug` array (volume-CGNS diagnostics) and state allocation. |
| `src/modules/blockPointers.F90` | Pointer aliases to the block arrays incl. `transitionDebug`. |
| `src/modules/cgnsNames.f90` / `src/modules/extraOutput.f90` / `src/output/outputMod.F90` | CGNS variable names and output dispatch for γ, Re̅θt, and transition-debug fields. |
| `src/NKSolver/blockette.F90` | Blockette (cache-blocked) residual path: `saGammaRethetaSource/Advection/Viscous/ResScale` (lines ~6903–7559) and dispatch at ~630/822. |
| `src/adjoint/masterRoutines.F90` | Hand-written AD driver; routes transition-model residual through the AD master routines. |
| `python/pyADflow.py` | Python wrapper: option plumbing and output for the transition model. |

## 2. Tapenade-differentiated files

Suffix conventions: `_d` = forward (tangent) mode; `_b` = reverse (adjoint) mode; `_fast_b` = reverse "fast" state-only mode. Finalized files live in `src/adjoint/output{Forward,Reverse,ReverseFast}/`; raw Tapenade intermediates in `src/adjoint/temp_{forward,reverse,reverse_fast}/`. All are generated — never hand-edited.

| Primal file | Forward | Reverse | Reverse-fast |
|-------------|---------|---------|--------------|
| `saGammaRetheta.F90` | `outputForward/saGammaRetheta_d.f90` | `outputReverse/saGammaRetheta_b.f90` | `outputReverseFast/saGammaRetheta_fast_b.f90` |
| `turbBCRoutines.F90` | `outputForward/turbBCRoutines_d.f90` | `outputReverse/turbBCRoutines_b.f90` | `outputReverseFast/turbBCRoutines_fast_b.f90` |
| `turbUtils.F90` | `outputForward/turbUtils_d.f90` | `outputReverse/turbUtils_b.f90` | `outputReverseFast/turbUtils_fast_b.f90` |
| `initializeFlow.F90` | `outputForward/initializeFlow_d.f90` | `outputReverse/initializeFlow_b.f90` | `outputReverseFast/initializeFlow_fast_b.f90` |
| `BCData.F90` | `outputForward/BCData_d.f90` | `outputReverse/BCData_b.f90` | (temp only: `temp_reverse_fast/BCData_fast_b.f90`) |
| `sa.F90` | `outputForward/sa_d.f90` | `outputReverse/sa_b.f90` | `outputReverseFast/sa_fast_b.f90` |

Matching intermediates exist under `src/adjoint/temp_forward/`, `temp_reverse/`, `temp_reverse_fast/` (plus differentiated module shells `block_d.f90`, `blockPointers_d.f90`, etc.). Branch-specific AD wiring is documented in `docs/VERIFICATION/adjoint-trace.md`.

## 3. Nondimensionalization

| File | Role |
|------|------|
| `src/initFlow/initializeFlow.F90` — `subroutine referenceState` (lines 10–190) | Computes the p-ρ non-dimensional reference state: `pInf`, `rhoInf`, `uInf = M√(γ pInf/rhoInf)`, `muInf`, `rGas`, and free-stream state `wInf`. |
| `src/modules/flowVarRefState.F90` | Module holding the reference quantities (`pRef`, `rhoRef`, `muInf`, `gammaInf`, `wInf`, `Re`, …) used everywhere. |
| `docs/nondimensionalization.md` | KB doc explaining the p-ρ scaling (velocity → M√γ, viscosities as ratios to μ∞, 1/Re not absorbed). |

Transition variables through it: `wInf(itu2)` (intermittency γ∞) and `wInf(itu3)` (Re̅θt∞, set via `reThetaTCorrelation` at free-stream turbulence intensity) are set in `referenceState` under `case (spalartallmarasnoft2gammaretheta)` (initializeFlow.F90:140–144); γ is already dimensionless and Re̅θt is a Reynolds-number-like dimensionless variable, so neither is rescaled by pRef/rhoRef — but every velocity/viscosity entering the source terms (`uInf`, `muInf`, `rlv`, `rev`, time scales) is in p-ρ units.

## 4. Solver hooks (DADI / ANK)

| File | Role |
|------|------|
| `src/turbulence/turbAPI.F90` — `turbSolveDDADI` | Decoupled DADI entry: selects `saSolve` / SST solve / `saGammaReThetaSolve` by `turbModel`. |
| `src/turbulence/saGammaRetheta.F90` — `saGammaReThetaSolve` | The model's DD-ADI (additive-form) implicit solve, incl. the 3×3 / decoupled block variants; `evalSrcJacBlock` supplies the source Jacobian. |
| `src/turbulence/sa.F90` / `src/turbulence/SST.F90` | SA and SST DD-ADI solves (same pattern, dispatched from `turbAPI`). |
| `src/NKSolver/NKSolvers.F90` | ANK and NK solvers; coupled-ANK path takes flow+turbulence in one Newton–Krylov system (note at line ~2414 contrasts its Jacobian form with the model's DD-ADI). |
| `src/NKSolver/blockette.F90` | Fast residual evaluation used by ANK/NK: per-model residual kernels incl. the `saGammaRetheta*` set; also calls `saGammaRetheta_block(.true.)` on the non-blockette path. |
| `src/solver/solvers.F90`, `src/solver/smoothers.F90`, `src/solver/residuals.F90` | RK/DADI flow-side smoothers and master residual assembly that call into the turbulence layer. |

## 5. Reference papers

| Reference | Path / status |
|-----------|---------------|
| Piotrowski & Zingg 2020, SA-sLM2015 (γ-Re̅θt) — physics source of truth | `docs/SA_GAMMA_RETHETHA_BASE/Piotrowski_Zingg_2020_SA-sLM2015_clean (1).md` (full text, per CLAUDE.md rule 9) and PDF `doc/PiotrowskiandZingg2020 (2) (1).pdf`. |
| Original Langtry–Menter 2009 correlation paper | **Not present** in the repo — supply if needed. |
| ADflow paper (Mader et al. 2020 / Kenway et al.) | Not present as text/PDF; only cited in `doc/citation.rst` / `doc/citations.bib`. Supply the PDF if needed. |
| Local theory KB (derived, not primary) | `docs/ADFLOW_BASE/ADFLOW_00…05` and `docs/architecture.md`. |

## 6. Existing SA/SST test suites (templates)

| File | Role |
|------|------|
| `tests/reg_tests/test_solve.py` | Forward-solve regression tests (SA and SST cases). |
| `tests/reg_tests/test_functionals.py` | Functional evaluation + partial-derivative checks. |
| `tests/reg_tests/test_adjoint.py` | Reverse-mode (adjoint) total-derivative tests; includes complex-step DVGeo comparison (`isComplex=True`, ~line 334). |
| `tests/reg_tests/test_jacVecProdFWD.py` | Forward-mode (tangent) Jacobian-vector product tests. |
| `tests/reg_tests/test_jacVecProdBWDFast.py` | Reverse-fast (`_fast_b`) transpose-product tests. |
| `tests/reg_tests/reg_default_options.py`, `reg_aeroproblems.py`, `reg_test_classes.py`, `reg_test_utils.py` | Shared fixtures: default options, aero problems, base classes, ref-file utilities (`tests/reg_tests/refs/`). |
| `tests/unit_tests/` | Unit tests (`test_basics.py`, `test_files.py`, …). |
| `tests/test_srcLambda_eig.py` | Branch-local test of the source-Jacobian eigenvalue machinery (`computeSrcLambda`). |
| `src_cs/` (repo root) | Complexified build tree used for complex-step verification. |
