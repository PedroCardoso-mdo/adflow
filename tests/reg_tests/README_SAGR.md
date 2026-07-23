# SA-noft2-Gamma-Retheta (SA-GR) derivative tests

Test scaffolding for verifying the SA-GR transition model's AD derivatives
(Tapenade forward `_d`, reverse `_b`, reverse-fast `_fast_b`) against finite
differences and the complex-step (CS) build. Written per
`docs/audits/08_test_prep.md`; background in `docs/audits/06_adjoint_wiring.md`.

**2026-07-23: case is the SA tutorial-wing grid** (`mdo_tutorial_sagr_dp.cgns`,
an NK-converged state on the standard ADflow tutorial-wing mesh at Mach=0.15).
It was briefly on the AR5 plain-wing case (2026-07-19) but that stalled
chronically — see `reg_sagr.py`'s header comments for the full rationale.
The old `_flatplate` labels (class-name suffixes, ref JSON filenames) were a
legacy misnomer from before either switch and have now been **renamed to
`_tut_wing`** (refs, test `name`s, `ap_sagr_tut_wing`) to stop misleading —
the case has never been a flat plate. `ap_sagr_ar5_wing` remains as a
backward-compat alias for the dev/ diagnostic scripts.

## File tree

```
tests/reg_tests/
├── reg_sagr.py                      # shared config + helpers (everything SA-GR-specific lives here)
├── generate_sagr_restart.py         # one-off: produce the restart CGNS the tests linearize about
├── test_jacVecProdFWD_sagr.py       # forward (tangent) mode: AD vs ref / FD / CS
├── test_jacVecProdBWDFast_sagr.py   # reverse-fast vs reverse consistency + dot products
├── test_adjoint_sagr.py             # adjoint totals vs ref + complex-step totals
└── refs/                            # trained JSON reference files land here
    ├── jacvecfwd_sagr_tut_wing.json    (created by training)
    ├── jacvecbwd_sagr_tut_wing.json    (created by training)
    └── adjoint_sagr_tut_wing.json      (created by training)
```

They mirror the upstream SA suite one-to-one:

| SA-GR file | Mirrors | Verifies |
|---|---|---|
| `test_jacVecProdFWD_sagr.py` | `test_jacVecProdFWD.py` | `computeJacobianVectorProductFwd` (the `outputForward/*_d.f90` code) |
| `test_jacVecProdBWDFast_sagr.py` | `test_jacVecProdBWDFast.py` + dot products from `test_functionals.py` | `computeJacobianVectorProductBwdFast` (`outputReverseFast/*_fast_b.f90`) against the full reverse mode, and fwd↔rev transpose consistency |
| `test_adjoint_sagr.py` | `test_adjoint.py` | total sensitivities from the adjoint solve (`outputReverse/*_b.f90`) vs complex step |
| `reg_sagr.py` | `reg_default_options.py` + `reg_aeroproblems.py` + `reg_test_utils.py` | — (fixtures + assert helpers) |

## What each file contains

### `reg_sagr.py`
- **Case config**: grid/restart/FFD paths, `sagrBaseOptions` (solver options:
  DADI+NK, no wall functions, no ft2, Tu∞ > 0, `eddyVisInfRatio=1e-10`),
  the AeroProblem, and the aero DVs (`alpha, mach, P, T` — mach/P/T drive
  uInf/muInf and the farfield `wInf(itu2/itu3)`, i.e. audit-06 finding F1).
  **The grid/AeroProblem values are placeholders until
  `generate_sagr_restart.py` has been run** — the test conditions must match
  the restart's exactly.
- **Block helpers**: `getStateBlocks` / `maskStateVector` split any flattened
  state/residual vector into `meanflow` / `nuTilde` / `gamma` / `reThetat`
  blocks (the vectors are variable-fastest, `w(i,j,k,1:nw)`, nw=8).
- **Assert helpers** (the SA-GR-specific extensions):
  - `assert_coupling_blocks_allclose` — norms of all 16 `dR[row]/dw[col]`
    blocks (column-masked forward product, row-masked result). Covers every
    SA↔transition cross term from audit 06 §4, incl. the structurally-zero
    `dR[nuTilde]/dw[reThetat]` as a regression guard.
  - `assert_transition_xdvdot_allclose` — transition-row norms of
    `dR/d(mach,P,T,alpha)`: the targeted farfield-BC / vorticity-limiter (F1)
    check.
  - `assert_bwdfast_blocks_allclose` — `_b` vs `_fast_b` with `resBar` seeded
    one equation-row block at a time (localizes the `autoEditReverseFast.py`
    stripping hazard, `docs/audits/sst_dev_lessons.md` watch item 1).
  - `assert_coupling_dot_products_allclose` — blockwise transpose tests
    isolating each off-diagonal Jacobian block in fwd-vs-rev consistency.

### `test_jacVecProdFWD_sagr.py`
- `TestJacVecFwdSAGR` (real build, ref-based): wDot/xVDot/xDvDot norms +
  the coupling-block and transition-row-xDvDot norms.
- `TestJacVecFwdSAGRFD` (real build, self-contained): AD vs FD, incl.
  gamma/reThetat columns in isolation. Kinks (vorticity cap, smoothMinMax)
  cause FD noise — don't loosen tolerances; CS is the decisive check.
- `TestJacVecFwdSAGRCS` (complex build, `cmplx_test_*`): records CS values
  under the same handler keys the real-build class trained, so running it
  compares AD vs CS key-by-key. `cmplx_test_xDvDot_transition_rows` is the
  decisive post-Tapenade-rerun check for audit-06 F1.

### `test_jacVecProdBWDFast_sagr.py`
- `TestJacVecBWDFastSAGR`: full + row-block-seeded `_b` vs `_fast_b`
  (atol 1e-16) + repeated-call determinism.
- `TestDotProductsSAGR`: the standard w/Xv→R/F dot products + the blockwise
  coupling dot products. Transpose tests prove *consistency*, not
  correctness — correctness comes from the CS classes.

### `test_adjoint_sagr.py`
- `TestAdjointSAGR` (real build): residuals + `evalFunctionsSens` totals vs
  ref, for `cl, cd, cmz, drag` (cd/cmz deliberately included — sst_dev
  post-mortem: never verify `cl` only). Uses `mdo_tutorial_ffd.fmt` + the
  twist/span/shape ref-axis DVGeo the SA `test_adjoint.py` uses. Needs
  `adjointMaxIter=3000` (set in `reg_sagr.py`): the SA-GR 8-state adjoint KSP
  takes ~583–597 iters to reach `adjointL2Convergence=1e-14`; the default 500
  cap clips it just short and mis-flags a converging adjoint as failed.
- `TestCmplxStepSAGR` (complex build): re-converges the complex solver with a
  1e-40j perturbation per DV (`alpha`, `mach`, and `twist`/`span`/`shape` via
  DVGeo/idwarp) and compares imag(f)/h against the adjoint totals in the ref.
  **Complex-build caveat:** the complexify build excludes the Tapenade AD
  routines, so the AD preconditioner the real solve relies on is unavailable —
  `setUp` overrides `ankadpc/nkadpc=False` + `ankcoupledswitchtol=1e-16` so the
  complex primal re-converges with FD-colored PC on the decoupled ANK→NK path
  (otherwise it aborts with "Forward AD routines are not complexified"). That
  FD-PC path stalls ~1e-8 (STRATEGY.md: FD-PC weak on SA-GR), so the complex
  functions — and hence CS-vs-adjoint agreement — cap at ~1e-8 absolute. The
  adjoint totals are correct to that level; whether to loosen the inherited
  `5e-9` tol or push the complex solve deeper is an **open item** (see
  `docs/task-log/2026-07-23-sagr-full-adjoint-test.md`).

### `dev/diag_full_derivatives.py`
Dev diagnostic (NOT a testflo test) for the full df/dx check with ADflow's
flow-solver / adjoint-KSP output VISIBLE, over ALL aero + geom DVs incl.
multiple individual `shape` components. `--mode adjoint` (real build: dump all
adjoint totals) writes a json that `--mode cs` (complex build: per-DV complex
re-converge, tabulated vs the adjoint) reads. Each CS re-converge is a full
nonlinear solve (~8 min), so the DV/shape set is configurable (`--shape`,
`--twist`, `--skip-aero/geom`).

### `generate_sagr_restart.py`
The upstream SA restarts are downloaded pre-made (`input_files/
get-input-files.sh`); this script is the SA-GR equivalent. It converges the
**AR5 plain-wing grid** (`ar5_plain_wing_vol_L3.cgns`, copied into
`input_files/`, gitignored like the rest of that directory) with SA-GR
through the validated ANK->CANK->NK production ladder (STRATEGY.md /
`.../3D_Plain_Wing/best_strategie/run_strategy.py`) and writes a
grid+solution CGNS containing all 8 states (the SA-GR restart set includes
Intermittency/ReThetat — `outputMod.F90`). Still writes the file with a
WARNING if `L2Convergence` isn't reached within `ncycles` (expected here:
this case hits a reproducible NK stall, so `ncycles` — not
`L2Convergence` — controls how deep the restart is; see the script's own
`--l2`/`--ncycles` help text).

## How to run (in order)

Prereqs: `testflo` and `parameterized` in the test env
(`pip install .[testing]`), and `input_files/get-input-files.sh` run once.

```bash
# 0. (once, real build) generate the restart the tests linearize about --
#    defaults already match ap_sagr_ar5_wing (mach=0.2, alpha=0.0, tu=0.0025),
#    so no flags are needed unless deviating from that case
cd tests/reg_tests
python generate_sagr_restart.py
#    -> confirm reg_sagr.py's sagrRestartFile/ap_sagr_ar5_wing still match
#       the script's printed summary (see script output)

# 1. (real build) self-contained FD sanity check — no ref files needed
testflo test_jacVecProdFWD_sagr.py:TestJacVecFwdSAGRFD -v

# 2. (real build) reverse-fast consistency — no ref files needed
testflo test_jacVecProdBWDFast_sagr.py -v

# 3. (real build) TRAIN the reference files (writes refs/*.json)
#    !! only after the Tapenade rerun for audit-06 F1 (vortlimd) + make !!
testflo test_jacVecProdFWD_sagr.py test_adjoint_sagr.py test_jacVecProdBWDFast_sagr.py -m "train*"

# 4. (real build) regression run against the trained refs
testflo test_jacVecProdFWD_sagr.py test_adjoint_sagr.py test_jacVecProdBWDFast_sagr.py

# 5. (complex build: libadflow_cs.so present, complex idwarp/pygeo) the
#    decisive AD-vs-CS comparison
testflo test_jacVecProdFWD_sagr.py test_adjoint_sagr.py -m "cmplx_test*"
```

**Building the complex lib:** use the `/build` skill, or
`PETSC_ARCH=complex-debug make -f Makefile_CS` (the complex build MUST point at
the complex PETSc arch; the default `real-debug` arch gives COMPLEX→REAL PETSc
type-mismatch errors in fortranPC/adjointAPI). After any Fortran rebuild,
`pip install .` into the mach env (site-packages, not `./adflow`, is what the
tests import — see [[tests-load-site-packages-not-repo]]). The whole suite is
also driven by `./run_sagr_tests.sh {all,real,cs,adjoint,train,genw}`.

Single test example:
`testflo test_jacVecProdFWD_sagr.py:TestJacVecFwdSAGR_0_sagr_tut_wing.test_coupling_blocks`
(use `testflo --dryrun` to list exact signatures).

## Verification ladder / how to read failures

1. **FD class fails, CS passes** → FD step/kink noise, not an AD bug.
2. **CS coupling-block key fails** → that specific `dR[row]/dw[col]` block is
   wrong in the forward code; audit 06 §4 lists where each block lives.
3. **`_fast_b` block-row test fails on gamma/reThetat seeds only** → suspect
   `autoEditReverseFast.py` push/pop stripping before suspecting the model.
4. **Dot products fail** → fwd and rev disagree with each other (bookkeeping);
   they can also both be wrong identically and still pass — hence CS.
5. **`cmplx_test_xDvDot_transition_rows` / mach-DV adjoint CS fails** →
   audit-06 F1 territory (uInf/muInf limiter path) — check the Tapenade
   regeneration happened before training.

Never fix a failing comparison by inflating its tolerance
(`docs/audits/sst_dev_lessons.md` §3).
