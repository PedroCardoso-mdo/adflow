# 08 — Derivative-test preparation for SA-γ-R̃e_θt (pre-complex-build)

**Date:** 2026-07-08. Branch `sa_gamma_rethetha`. Preparation only — nothing
was built or executed (the complex-step ADflow build is not installed yet).
Follows audits `06_adjoint_wiring.md` and `adjoint_audit_2026-07-07.md`
(the dated visual-wiring pass) plus `sst_dev_lessons.md`.

Note: prior audits live in `docs/audits/` (the old top-level `audits/` was
deleted when docs were consolidated), so this file is
`docs/audits/08_test_prep.md` rather than `audits/08_test_prep.md`.

---

## 1. The existing derivative test suite (what was located, exactly)

All derivative verification lives in `tests/reg_tests/`. **There are no SST
derivative tests on this branch** — the fork predates upstream PR #331
(`sst_dev`); `grep -ri sst tests/` is empty. The SA coverage (via the
`rans_tut_wing` cases on `mdo_tutorial_rans_scalar_jst.cgns`, SA being the
RANS default) is therefore the template. The three modes and where each is
verified:

| Mode | Generated code | Python entry point | Test file / class | Verification method |
|---|---|---|---|---|
| Forward (tangent) | `outputForward/*_d.f90` | `computeJacobianVectorProductFwd(wDot/xVDot/xDvDot, ...)` | `test_jacVecProdFWD.py`: `TestJacVecFwd` (ref file), `TestJacVecFwdFD` (vs FD, h=1e-8), `TestJacVecFwdCS` (vs complex step, h=1e-40, `ADFLOW_C`) | norms vs trained ref; FD rtol 8e-4 (res) / 1e-5 (funcs); CS records into the same ref keys |
| Reverse (adjoint) | `outputReverse/*_b.f90` | `computeJacobianVectorProductBwd(resBar/fBar/funcsBar, ...)`; totals via `evalFunctionsSens` | `test_functionals.py` (`test_jac_vec_prod_bwd`, `test_dot_products`), `test_adjoint.py`: `TestAdjoint` (adjoint totals vs ref) + `TestCmplxStep` (totals vs re-converged complex solve, aero DVs + DVGeo `isComplex=True` geometric DVs) | dot-product (transpose) tests tol 1e-10–2e-10; adjoint totals tol 1e-10; CS rtol 1e-8 (aero) / 5e-9 (geom) |
| **"Fast" mode** | `outputReverseFast/*_fast_b.f90` | **`computeJacobianVectorProductBwdFast(resBar=...)`** — this is ADflow's precise name for the lower-cost path: a **state-only reverse mode** (dR/dw^T products only, no mesh/DV/function seeds), used to assemble the adjoint/PC matrices cheaply | `test_jacVecProdBWDFast.py`: `TestJacVecBWDFast` | consistency vs the full reverse mode, `atol=1e-16`, plus a repeated-call determinism check |

Shared fixtures: `reg_default_options.py`, `reg_aeroproblems.py`,
`reg_test_classes.py` (`RegTest` trains/compares a `BaseRegTest` JSON ref in
`refs/`; `CmplxRegTest` never trains — its `cmplx_test_*` methods compare
against the ref the real build trained), `reg_test_utils.py` (the
`assert_*_allclose` helpers).

Ladder implied by the suite (and by the `sst_dev` post-mortem): FD sanity →
CS (decisive) → fwd/rev dot products → `_fast_b` vs `_b` consistency →
adjoint totals vs CS.

## 2. New SA-GR test scaffolding (written, not executed)

Four new files in `tests/reg_tests/`, following the upstream structure and
naming (`parameterized_class` cases, `RegTest`/`CmplxRegTest` bases, one JSON
ref per file in `refs/`):

| File | Mirrors | Classes |
|---|---|---|
| `reg_sagr.py` | `reg_default_options.py` + `reg_aeroproblems.py` + `reg_test_utils.py` | shared SA-GR options/AeroProblem + block-masking helpers + 4 new assert helpers |
| `test_jacVecProdFWD_sagr.py` | `test_jacVecProdFWD.py` | `TestJacVecFwdSAGR`, `TestJacVecFwdSAGRFD`, `TestJacVecFwdSAGRCS` |
| `test_jacVecProdBWDFast_sagr.py` | `test_jacVecProdBWDFast.py` (+ dot products from `test_functionals.py`) | `TestJacVecBWDFastSAGR`, `TestDotProductsSAGR` |
| `test_adjoint_sagr.py` | `test_adjoint.py` | `TestAdjointSAGR`, `TestCmplxStepSAGR` |

Ref files (to be trained on the real build): `jacvecfwd_sagr_flatplate.json`,
`jacvecbwd_sagr_flatplate.json`, `adjoint_sagr_flatplate.json`.

### 2.1 How the extensions map to the brief

**New state variables in isolation.** The flattened state/residual vectors
are variable-fastest (`w(i,j,k,1:nw)`, verified in
`NKSolvers.F90:getStates`), and `flowvarrefstate.nw/nt1/nt2` are f2py-exposed,
so `reg_sagr.getStateBlocks`/`maskStateVector` split any state-shaped vector
into `meanflow` / `nuTilde` / `gamma` / `reThetat` blocks identically on the
real and complex builds (a guard raises unless nw=8 with 3 turb variables).
`test_wDot_transition_columns` (FD), the column-seeded coupling tests, and
the row-seeded `_fast_b` tests all use these masks.

**Coupling partials flagged in audit 06 (§4).**
`assert_coupling_blocks_allclose` records the norm of every 4×4 block of
dR/dw (column-masked forward product, row-masked result), under identical
handler keys in the AD and CS classes — so training on the real build and
running `cmplx_test_coupling_blocks` compares each block AD-vs-CS, including:
d(R_ν̃)/dγ (the Eq. 41 production coupling), d(R_γ)/dν̃ (r_T chains),
d(R_γ)/dR̃ (correlation chains), d(R_θ)/d(meanflow), d(R_θ)/dν̃ (nonzero
only via crossflow D_scf — finding F5's residual-level counterpart), and the
structurally-zero d(R_ν̃)/dR̃ block as a regression guard.
`assert_coupling_dot_products_allclose` runs blockwise transpose tests on the
same pairs, isolating each off-diagonal block in forward-vs-reverse
consistency. `assert_bwdfast_blocks_allclose` seeds `resBar` one row block at
a time so a `_fast_b` failure localizes to an equation (targets the
`autoEditReverseFast.py` stripping watch item from `sst_dev_lessons.md`).

**Transition-model BC derivatives.** Two mechanisms:
(i) wall/farfield ghost-cell chains (ν̃ antisymmetric; γ, R̃ Neumann at
walls; `bvt = wInf` at farfield) sit inside every dR/dw product because
`applyAllTurbBC` runs inside the residual — the column/row-masked products
cover them; (ii) the flow-condition side (farfield ghost `wInf(itu2/itu3)`
and the uInf/muInf vorticity-limiter path — **audit-06 F1**) is targeted by
`assert_transition_xdvdot_allclose`: transition-row norms of dR/d(mach, P, T,
alpha), recorded by both the AD and CS classes.
`cmplx_test_xDvDot_transition_rows` is the decisive post-Tapenade-rerun check
for F1. The adjoint CS class prioritizes `mach` for the same reason, and
`evalFuncs` includes `cd`/`cmz` (not `cl` only — sst_dev post-mortem).

### 2.2 Syntactic status

All four files pass `python -m py_compile`; `reg_sagr.py` imports cleanly in
this environment (baseclasses present). The three test modules import up to
`ModuleNotFoundError: parameterized` — the identical import the existing
upstream tests use, i.e. an environment gap, not a code defect (install
`parameterized`, or run under the MACH test environment/`testflo`, which
already runs the upstream suite). `BaseRegTest.par_add_norm/par_add_sum/
root_add_dict` signatures were checked against the installed baseclasses.

**Placeholders that must be filled before execution** (all in `reg_sagr.py`,
single location): grid/restart CGNS (`sagr_flatplate_t3a.cgns`), FFD file
(`sagr_flatplate_ffd.fmt`), and the AeroProblem flow conditions
(mach/Re/T/Tu placeholder values for a T3-series flat plate). The restart
must be a converged SA-GR solution (the suite's pattern: derivatives are
linearized about a physical state loaded via `restartfile`).

## 3. Expected complex-build flags / configuration

From `src_cs/build/Makefile{,1}` and `config/`:

1. `config/config.mk` must define `COMPLEXIFY_DIR` with
   `COMPLEXIFY_INCLUDE_FLAGS`/`COMPLEXIFY_LINKER_FLAGS`
   (`-L$(COMPLEXIFY_DIR)/lib -lcomplexify`) — the `complexify` package
   (pip `complexify`, also in `setup.py` extras `[complex]`) provides both
   the source translator and the runtime library.
2. Source generation: `make complexify` via `src_cs/build/Makefile1` runs the
   `complexify` translator over every file in `src/build/fileList` **minus
   `adjoint/output{Forward,Reverse,ReverseFast}`** (the Tapenade files are
   excluded from the complex build — so the complex solver contains no AD
   code, only the complexified primal; this is why the CS classes record into
   ref files rather than comparing AD-vs-CS in-process).
3. Compilation: `src_cs/build/Makefile` with `-DUSE_COMPLEX` and
   `FF90_PRECISION_FLAGS`, linking `libadflow_cs.so` and moving it into
   `adflow/` next to `libadflow.so`. `ADFLOW_C` (`adflow/pyADflow_C.py`)
   loads `libadflow_cs` and sets `dtype="D"`.
4. Because SA-GR modifies primal files (`saGammaRetheta.F90`, `turbUtils.F90`,
   `turbBCRoutines.F90`, `initializeFlow.F90`, `blockette.F90`, module files),
   all of these get complexified — no SA-GR-specific build change is
   expected. Watch: `complexify` may need updating on new intrinsics
   (`tanh`, `log`, `sinh` are already used elsewhere; `smoothMinMax`
   and correlations use only exp/log/tanh/min/max — should pass).
5. Python side: complex idwarp (`USMesh_C`) and complex-capable pygeo
   (`DVGeometry(isComplex=True)`) are required for `TestCmplxStepSAGR`
   (same as upstream `TestCmplxStep`).
6. Run pattern: real build → `testflo ... --train` (or `train_*` entries) to
   write refs → complex build → run `cmplx_test_*` methods.

**Pre-existing blocker to sequence first (audit-06 F1):** the current AD
files still have `vortlimd = 0.0_8`; `Makefile_tapenade` was fixed in
21f29567 but Tapenade has not been rerun. Train the FWD refs only *after*
the Tapenade rerun + `make`, otherwise `cmplx_test_xDvDot_transition_rows`
and the mach-DV adjoint CS check will (correctly) fail against refs trained
from the stale AD.

## 4. Checklist — needed from you to proceed

1. **Complex ADflow install location / status**: where `libadflow_cs.so`
   will live (expected: built in `src_cs/build`, moved to `adflow/`), and
   whether `complexify` (translator + runtime lib) is installed;
   `COMPLEXIFY_DIR` value for `config/config.mk`.
2. **Complex dependency stack**: confirm complex idwarp (`USMesh_C`) and
   pygeo with `isComplex=True` support are available in the test
   environment; also `parameterized` and `testflo` (the standard runners for
   this suite) in whatever env the tests will run.
3. **Test case inputs**: a converged (or smoke-converged) SA-GR CGNS
   grid+restart for a transition case (flat plate T3-series or the P&Z
   validation case) + its exact flow conditions (mach, Re, T, Tu∞,
   chordRef/areaRef) to replace the placeholders in
   `tests/reg_tests/reg_sagr.py`; an FFD `.fmt` file embedding that surface
   for the geometric-DV adjoint/CS tests (or say so and I drop
   `TestCmplxStepSAGR.cmplx_test_geom_dvs` / DVGeo wiring).
4. **Tolerance policy confirmation**: I seeded the SA values used elsewhere —
   ref-comparison 5e-9 (fwd wDot) / 1e-10 (xVDot, xDvDot, bwd, adjoint),
   dot products 2e-10 (the RANS value), `_fast_b` vs `_b` atol 1e-16, FD
   rtol 8e-4/1e-5, CS-vs-adjoint rtol 1e-8 (aero) / 5e-9 (geom), CS h=1e-40.
   Confirm, or supply the tolerances you want for the transition model
   (per sst_dev lessons I will *not* inflate a failing tolerance; masked/
   per-block metrics exist precisely to localize instead).
5. **Sequencing**: confirm the Tapenade rerun (F1 fix) + `make` + commit of
   the six AD files happens before training the forward-mode refs.
6. **Solver settings for the case**: whether the restart converges tightly
   enough with the options in `reg_sagr.sagrBaseOptions` (DADI +
   NK, `nkswitchtol` 1e-3), or the smoke-script options you already use for
   this case should replace them.

**Stop point reached** — nothing built, nothing run, per the brief.
