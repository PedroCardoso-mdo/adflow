# Three-Stage Low-Level Adjoint Verification Ladder

**Status (plain SA): all three stages PASS**, for plain SA on both the
tutorial-wing and AR5 meshes, as of 2026-07-21.

**Status (SA-GR, `nw=8`, AR5): all three stages PASS**, as of 2026-07-22.
Stage 2 initially FAILED on the reThetat rows (`_b` vs `_fast_b` diverged,
max abs diff ≈ 7.5e+2); this was root-caused to an in-place `smoothMinMax`
double-clamp of `lambdaThetaLocal` in the primal and **fixed** (see the
"SA-GR (`nw=8`) results" section below). After the fix + Tapenade regen,
Stage 2's reThetat block matches to 2.3e-10.

This is the low-level, raw-API verification campaign — separate from (and
a prerequisite trusted-input for) the higher-level `evalFunctionsSens`-based
campaign tracked in `../current-task.md`.

This ladder answers one question at three levels of the AD stack: **does
the differentiated code actually compute the derivative of the primal
code**, independent of any specific functional or DV. It does not check
whether the adjoint *total* sensitivities match a physical/trained
reference — that is `../current-task.md`'s job. Read
`VERIFICATION/debugging_derivatives.md` first for the generic
FD→dot-product→CS theory behind why each stage below is structured the
way it is.

## Why three stages, in this order

1. **Reverse ↔ forward consistency (dot-product test).** Cheapest, most
   fundamental check: for the same linearized operator, forward-mode
   (`_d`, computes `dR/dw · wDot`) and reverse-mode (`_b`, computes
   `dR/dw^T · psi`) must satisfy `<wDot, dR/dw^T psi> = <psi, dR/dw wDot>`
   to machine precision. This only checks fwd/rev *bookkeeping* agree with
   each other — it cannot catch a bug present identically in both modes.
2. **Reverse vs. fast-reverse consistency (`_b` vs `_b_fast`).** ADflow
   ships two reverse-mode variants; `_fast_b` is a stripped-down Tapenade
   output (`autoEditReverseFast.py`) used for performance. If this
   diverges from `_b`, `autoEditReverseFast.py`'s push/pop stripping is
   the first suspect (this exact failure mode broke SST's `fast_b`
   upstream — see `../audits/07_sst_dev_lessons.md`).
3. **Forward verification test (3-way AD vs. FD vs. CS).** Only once 1
   and 2 agree does it make sense to ask whether the *shared* fwd/rev
   linearization is actually correct — i.e. matches the true derivative of
   the primal. Complex-step (h=1e-40, no subtractive cancellation) is the
   authoritative ground truth here; single-step real FD is included only
   for context (expected to disagree at ~1e-1 to ~1e-2 relative error due
   to truncation/kinks, not a failure).

## Stage 1 — Reverse ↔ forward dot-product consistency

**What it checks:** `<v, dR/dw^T w> == <w, dR/dw v>` for random probe
vectors `v, w`, using the raw `computeJacobianVectorProdFwd`/`Bwd`
Fortran API (bypasses `evalFunctionsSens`/pyOptSparse entirely).

**Files:**
- `tests/reg_tests/sanity_check_bwdfast_stage1.py` — plain-SA,
  tutorial-wing, `_fast_b` vs `_b` only (no reg_sagr helpers), plus a
  repeated-call determinism check. `atol=1e-16`.
- `tests/reg_tests/sanity_check_partials_sa.py` — generalizes
  `reg_sagr.py`'s block-splitting dot-product helper
  (`utils.assert_dot_products_allclose`, `tol=2e-10`) to plain SA (`nw=6`,
  no gamma/reThetat rows), selectable mesh via `--mesh {tutorial,ar5}`
  (default `ar5`).

**How to run:**
```bash
mpirun -n 2 /home/mdo/packages_v2/mach/bin/python tests/reg_tests/sanity_check_bwdfast_stage1.py
mpirun -n 2 /home/mdo/packages_v2/mach/bin/python tests/reg_tests/sanity_check_partials_sa.py --mesh tutorial
mpirun -n 2 /home/mdo/packages_v2/mach/bin/python tests/reg_tests/sanity_check_partials_sa.py --mesh ar5
```

**Results:**
- `sanity_check_bwdfast_stage1.py`: **PASSED** — `=== STAGE 1 BWD-vs-BWDFast
  PASSED ===`.
- `sanity_check_partials_sa.py --mesh tutorial`: **PASSED** —
  `=== SANITY CHECK PASSED ===`, dot products matching to displayed
  precision.
- `sanity_check_partials_sa.py --mesh ar5`: **PASSED** — same. Note the
  script's own docstring records that an *earlier* iteration of this test
  found a real w→R dot-product mismatch on the AR5 mesh (rel err ~1.4e-5
  vs. tol 2e-10) that was root-caused and fixed before this pass (see
  `../current-task.md`'s "lesson learned" note on not attributing
  discrepancies to non-determinism).

## Stage 2 — Reverse vs. fast-reverse (`_b` vs `_b_fast`) consistency

Folded into the same scripts as Stage 1 (both call the block-seeded
`_b`-vs-`_fast_b` comparison, `reg_sagr.assert_bwdfast_blocks_allclose`,
`atol=1e-16`) — see files/commands above. For plain SA (`nw=6`, no
gamma/reThetat state), the block-seeding loop has no transition rows to
seed, so this check is effectively a no-op confirmation for the plain-SA
path; it exercises real content once SA-GR (`nw=8`) is run through the
same harness (not yet done — see `../current-task.md`).

## Stage 3 — Forward verification: 3-way AD vs. FD vs. CS

**What it checks:** `dR/dw · wDot`, reduced against a fixed probe vector,
computed three independent ways: Tapenade forward-mode AD, real-mode FD
(`h=1e-8`), and complex-step (`h=1e-40`, requires the complex ADflow
build). AD and CS should agree to full precision; FD is expected to be
off by O(h) due to truncation/kinks.

**Files:**
- `tests/reg_tests/check_3way_fwd.py` — single-h comparison,
  `--mesh {tutorial,ar5}`, `--build {real,complex}`. Must be run **twice
  per mesh** (once per build) since the real and complex ADflow libraries
  cannot coexist in one Python process.
- `tests/reg_tests/check_3way_fwd_sweep.py` — FD step-size sweep
  (`h=1e-2…1e-12`, one-sided built-in vs. manual centered FD) reusing the
  AD value from `check_3way_fwd.py` as reference (valid since AD≈CS was
  already confirmed independently).

**How to run:**
```bash
# Stage 3a: single-h 3-way check, tutorial-wing
mpirun -n 2 /home/mdo/packages_v2/mach/bin/python tests/reg_tests/check_3way_fwd.py --mesh tutorial --build real
mpirun -n 2 /home/mdo/packages_v2/mach_cs/bin/python tests/reg_tests/check_3way_fwd.py --mesh tutorial --build complex

# Stage 3a: single-h 3-way check, AR5 mesh
mpirun -n 2 /home/mdo/packages_v2/mach/bin/python tests/reg_tests/check_3way_fwd.py --mesh ar5 --build real
mpirun -n 2 /home/mdo/packages_v2/mach_cs/bin/python tests/reg_tests/check_3way_fwd.py --mesh ar5 --build complex

# Stage 3b: FD step-size sweep
mpirun -n 2 /home/mdo/packages_v2/mach/bin/python tests/reg_tests/check_3way_fwd_sweep.py --mesh tutorial
mpirun -n 2 /home/mdo/packages_v2/mach/bin/python tests/reg_tests/check_3way_fwd_sweep.py --mesh ar5
```
(Adjust the complex-build Python path to whatever env has `libadflow_cs.so`
installed — see `build_complex.log` for the last complex build on this
branch.)

**Results:**

| Mesh | AD | CS | FD (h=1e-8) | rel_err(FD,AD) |
|---|---|---|---|---|
| tutorial-wing | 5.1750905670e+10 | 5.1750905670e+10 (exact match to AD, 10 sig figs) | 5.7186962842e+10 | 1.050e-01 |
| AR5 plain wing | 3.2957137433e+13 | 3.2957137432e+13 (exact match to AD) | 5.6967464203e+13 | 7.285e-01 |

AD and CS agree to full displayed precision on both meshes — this is the
key pass criterion. The single-step FD mismatch is expected (not a
failure): the sweep below confirms FD converges to the same AD value as
`h→0`.

FD step-size sweep (one-sided FD relative error vs. AD, shrinking
~linearly with `h`, confirming first-order convergence to the AD/CS
value):

| h | rel_err (tutorial) | rel_err (AR5) |
|---|---|---|
| 1e-6 | — | 1.156e+02 |
| 1e-7 | — | 1.073e+01 |
| 1e-8 | 1.050e-01 | 7.285e-01 |
| 1e-9 | — | 4.113e-02 |
| 1e-11 | — | 3.872e-04 |
| 1e-12 | — | 3.870e-05 |

(Tutorial-wing sweep intermediate values weren't retained in the log
excerpt used to write this table; only the AR5 sweep's full h-ladder was.
Re-run `check_3way_fwd_sweep.py --mesh tutorial` to regenerate if needed.)

## Relationship to the still-open AR5 CS mismatch

`../current-task.md` separately documents a **different, higher-level**
mismatch: `evalFunctionsSens`'s `dcd/dmach` on AR5 (adjoint=-0.537 vs.
CS=-231.1). That is **not** a contradiction of Stage 3's PASS above —
Stage 3 checks the raw `dR/dw` Jacobian-vector product against a fixed
probe, which is confirmed correct; the `dcd/dmach` discrepancy is
downstream, in the *converged-state-dependent* total sensitivity through
a whole-solve `resetFlow`+re-solve, and has been traced to an
absolute-convergence-precision floor on that specific stalling mesh/case
(see `../current-task.md`'s h-sweep table via `sweep_h_fd.py`), not a bug
in the linearization tested here.

## SA-GR (`nw=8`) results

Run 2026-07-21, AR5 mesh only, linearizing about the converged SA-GR
restart `input_files/ar5_plain_wing_sagr_dp.cgns` (built 2026-07-19 by
`generate_sagr_restart.py`, AR5 grid, mach=0.2, alpha=0.0, Tu=0.0025 —
matches `ap_sagr_ar5_wing`). The plain-SA scripts were extended with a
`--turbmodel {sa,sagr}` flag; `sagr` pulls `reg_sagr.sagrBaseOptions` +
`ap_sagr_ar5_wing` + `sagrAeroDVs`, and the `nw=8` block layout is already
handled by `reg_sagr.getStateBlocks`. No Fortran changed (verification
scripts only).

**How to run:**
```bash
# Stage 1 + Stage 2 (dot products + block-seeded _b vs _fast_b)
OMP_NUM_THREADS=1 mpirun -n 2 /home/mdo/packages_v2/mach/bin/python tests/reg_tests/sanity_check_partials_sa.py --turbmodel sagr

# Stage 3 (3-way AD/FD/CS) — both builds live in the same mach env
OMP_NUM_THREADS=1 mpirun -n 2 /home/mdo/packages_v2/mach/bin/python tests/reg_tests/check_3way_fwd.py --turbmodel sagr --build real
OMP_NUM_THREADS=1 mpirun -n 2 /home/mdo/packages_v2/mach/bin/python tests/reg_tests/check_3way_fwd.py --turbmodel sagr --build complex

# Stage 3b (FD step-size sweep)
OMP_NUM_THREADS=1 mpirun -n 2 /home/mdo/packages_v2/mach/bin/python tests/reg_tests/check_3way_fwd_sweep.py --turbmodel sagr
```
(On this branch the complex `libadflow_cs.so` is installed in the same
`mach` env as the real build — `check_3way_fwd.py --build complex` imports
`ADFLOW_C`, so a separate `mach_cs` Python is not needed.)

### Stage 1 — dot-product consistency: **PASS**

fwd and rev reduce to the same value to all displayed digits (`tol=2e-10`):

| Dot-product test | fwd | rev |
|---|---|---|
| w → R | 5.546143e+13 | 5.546143e+13 |
| Xv → R | 9.426969e+09 | 9.426969e+09 |
| w → F | 1.428035e+08 | 1.428035e+08 |
| xV → F | 1.696683e+08 | 1.696683e+08 |
| (w, xV) → (dw, F) | 5.547117e+13 | 5.547117e+13 |

### Stage 2 — `_b` vs `_fast_b`, block-seeded: **PASS (after the lambdaTheta fix)**

`resBar` seeded one residual-row block at a time; reported is
`max|_b − _fast_b|` reduced across ranks.

| resBar seeded on | max\|_b − _fast_b\| (before fix) | max\|_b − _fast_b\| (after fix) |
|---|---|---|
| meanflow | 3.295e-06 | 3.295e-06 (passes rtol) |
| nuTilde | 7.629e-05 | 7.629e-05 (passes rtol) |
| gamma | 1.490e-08 | 1.490e-08 |
| reThetat | **7.546e+02 (FAIL)** | **2.328e-10 (PASS)** |

**Root cause and fix (2026-07-22).** The failure was isolated to the
reThetat-seeded block. A ground-truth check (`diag_which_wrong.py`:
forward `_d`/CS vs. `_b` vs. `_fast_b` for the `dR[reThetat]/dw[meanflow]`
block) proved `_b` and `_d`/CS agreed to ~1e-14 while `_fast_b` was wrong
by rel_err 1.94 (**wrong sign**: +3.949e4 vs. the correct −4.199e4) — so
`_fast_b` alone was broken. Root cause: the primal computed the
pressure-gradient parameter `lambdaThetaLocal` with two **in-place**
`smoothMinMax` clamps

```fortran
lambdaThetaLocal = smoothMinMax(lambdaThetaLocal, ...)   ! overwrites in place
lambdaThetaLocal = smoothMinMax(lambdaThetaLocal, ...)   ! overwrites again
```

Tapenade guarded the two overwritten intermediates with `pushreal8`/
`popreal8`, which `autoEditReverseFast.py` strips wholesale (it assumes
every stored primal is recomputed inline). Every other `smoothMinMax` in
the SA-GR source writes to a **distinct** target (`vortMagLim`,
`crossflowRatio`, `dHplus`, …) — so the recompute-to-survive convention
kept them safe; `lambdaThetaLocal` was the sole in-place deviation, and
it sits on the reThetat→meanflow chain (`lambdaThetaLocal =
thetaBL²/nu · dUds`), which is exactly why only the reThetat row / meanflow
column broke. **Fix:** rewrite the primal (`saGammaRetheta.F90`) with
distinct targets (`lambdaThetaRaw` → `lambdaThetaClamped` →
`lambdaThetaLocal`) — mathematically identical, and now matching the
convention everywhere else — then rerun Tapenade (`make -f Makefile_tapenade`)
and rebuild. After regen, `_fast_b` reads the correct distinct primals in
its two `smoothminmax_fast_b` reverse calls with no push/pop needed, and
the ground-truth check matches to 1.27e-13.

### Stage 3 — 3-way AD vs. CS vs. FD: **PASS**

| Quantity | Value | rel_err vs AD |
|---|---|---|
| AD | 3.3256409775e+13 | — |
| CS (h=1e-40) | 3.3256409775e+13 | exact match to AD (11 sig figs) |
| FD (h=1e-8, one-sided) | 6.0512532732e+13 | 8.196e-01 (expected FD noise) |

AD and CS agree to full displayed precision — the pass criterion. One-sided
FD step-size sweep confirms first-order convergence to the AD/CS value as
`h→0`:

| h | one-sided FD | rel_err vs AD |
|---|---|---|
| 1e-6 | 3.8365923470e+15 | 1.144e+02 |
| 1e-7 | 3.9729119395e+14 | 1.095e+01 |
| 1e-8 | 6.0512532732e+13 | 8.196e-01 |
| 1e-9 | 3.4952438464e+13 | 5.100e-02 |
| 1e-10 | 3.3418921994e+13 | 4.887e-03 |
| 1e-11 | 3.3272592605e+13 | 4.866e-04 |
| 1e-12 | 3.3258027319e+13 | 4.864e-05 |

(The sweep's centered-FD column returns nonsensical/`nan` values on this
SA-GR/AR5 state — same as the plain-SA runs, where only the one-sided
column was reported; the one-sided ladder above is the meaningful result.)

## What's not yet covered

- All three stages now PASS for both plain SA and SA-GR (`nw=8`) on AR5.
  The Stage-2 reThetat `_fast_b` divergence found on the first SA-GR run
  was root-caused and fixed (lambdaTheta in-place clamp; see the SA-GR
  results section) and re-verified against forward/CS ground truth.
- Note on the Stage-3 complex build: the CS run above uses the
  pre-existing `libadflow_cs.so`. A fresh complex rebuild after the fix
  was not required because the fix is mathematically neutral (the primal
  residual, and hence the forward AD value `3.3256409775e+13`, is
  byte-identical before and after — only the reverse-fast checkpointing
  changed), so CS is unchanged. (A from-scratch complex rebuild in this
  environment needs `PETSC_ARCH=complex-debug`; building against the
  default `real-debug` PETSc produces spurious COMPLEX→REAL interface
  errors unrelated to any source change.)
- `debug_cs_ar5_live.py` and `sweep_h_fd.py` are diagnostic/debugging
  scripts for the separate open AR5-CS issue above, not part of this
  pass/fail ladder — see `../current-task.md` for their role.
