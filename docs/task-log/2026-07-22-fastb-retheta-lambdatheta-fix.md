# Fix reverse-fast `dR[reThetat]/dw[meanflow]` (lambdaTheta in-place clamp) — 2026-07-22

**Problem:** The first-ever SA-GR run of the 3-stage adjoint ladder
(2026-07-21) found Stage 2 (`_b` vs `_fast_b`) FAILING on the reThetat
rows: `max|_b − _fast_b| ≈ 7.5e+2`. The reverse-fast product
`computeJacobianVectorProductBwdFast` is the operator the adjoint KSP
applies every iteration (`adjointAPI.F90:dRdwTMatMult`), so this was a
real defect in the production adjoint operator for the transition→mean-flow
coupling, not a test-only artifact.

**Diagnosis (ground-truth check):** `diag_which_wrong.py` compared the
`dR[reThetat]/dw[meanflow]` block computed three ways via the transpose
identity: forward `_d` (= CS, proven by Stage 3), full reverse `_b`, and
reverse-fast `_fast_b`. Result: `_d` = `_b` = −4.199e4 (agree to ~1e-14);
`_fast_b` = +3.949e4 — **wrong sign**, rel_err 1.94. So `_fast_b` alone was
broken. Root cause: the primal computed the pressure-gradient parameter
with two **in-place** `smoothMinMax` clamps —
`lambdaThetaLocal = smoothMinMax(lambdaThetaLocal, …)` twice. Tapenade
guards the overwritten intermediates with `pushreal8`/`popreal8`, which
`autoEditReverseFast.py` strips wholesale (it assumes each stored primal is
recomputed inline). Every *other* `smoothMinMax` in the SA-GR source writes
to a distinct target (`vortMagLim`, `crossflowRatio`, `dHplus`, …) — the
recompute-to-survive convention that makes the strip safe — so
`lambdaThetaLocal` was the sole in-place deviation, and it sits on the
reThetat→meanflow chain (`lambdaThetaLocal = thetaBL²/nu · dUds`), which is
exactly why only the reThetat row / meanflow column broke.

**Solution:** Rewrote the primal to use distinct targets, matching the
convention everywhere else in the file (mathematically identical):
- `src/turbulence/saGammaRetheta.F90` — added locals `lambdaThetaRaw`,
  `lambdaThetaClamped`; replaced the two in-place clamps with
  `lambdaThetaRaw = thetaBL²/nu·dUds` →
  `lambdaThetaClamped = smoothMinMax(lambdaThetaRaw, …)` →
  `lambdaThetaLocal = smoothMinMax(lambdaThetaClamped, …)`.
- Reran Tapenade (`make -f Makefile_tapenade default`,
  `TAPENADE_HOME=/home/mdo/packages_v2/tapenade_3.16`, headless), which
  regenerated `outputForward/{saGammaRetheta,turbUtils}_d.f90`,
  `outputReverse/{…}_b.f90`, `outputReverseFast/{…}_fast_b.f90`. The
  regenerated `_fast_b` now reads `lambdathetaclamped` (V1) and
  `lambdathetaraw` (V0) in its two `smoothminmax_fast_b` reverse calls,
  with zero push/pop needed.

**Files created/touched and why:**
- `src/turbulence/saGammaRetheta.F90` — the actual fix (distinct-target
  clamps + two new locals). Hand-written transition model, editable (not
  SA, not a Tapenade file).
- `src/adjoint/output{Forward,Reverse,ReverseFast}/{saGammaRetheta,turbUtils}_{d,b,fast_b}.f90`
  — Tapenade-regenerated from the fixed primal (the legitimate rule-6 path:
  regenerate, never hand-edit).
- `docs/VERIFICATION/three-stage-verification.md` — Stage 2 SA-GR updated
  FAIL → PASS with the root-cause/fix writeup and before/after table.
- `docs/current-task.md` — Stage-2 open item marked DONE.
- `tests/reg_tests/sanity_check_partials_sa.py`, `check_3way_fwd.py`,
  `check_3way_fwd_sweep.py` — the `--turbmodel sagr` harness (from the
  2026-07-21 task) used to reproduce the failure and confirm the fix.

**Verification (all on the real build unless noted):**
- Ground truth: `_fast_b` `dR[reThetat]/dw[meanflow]` vs forward `_d`
  1.94 (wrong sign) → **1.27e-13**.
- Stage 2 SA-GR reThetat block `max|_b − _fast_b|`: 7.546e+02 → **2.328e-10**;
  `sanity_check_partials_sa.py --turbmodel sagr` prints `SANITY CHECK PASSED`.
- Stage 3 SA-GR: AD = CS = **3.3256409775e+13** (the forward AD value is
  byte-identical before/after, confirming the primal change is
  mathematically neutral). CS re-confirmed against a **freshly rebuilt**
  `libadflow_cs.so` (`git clean -fdx src_cs/` +
  `PETSC_ARCH=complex-debug make -f Makefile_CS`), so nothing is stale.
- Plain-SA Stage 2 (regression): still `SANITY CHECK PASSED`.
- Build: real `make` + `pip install` clean; Tapenade regen exit 0.

**Follow-ups:**
- Complex-build recipe (verified this session): `git clean -fdx src_cs/`
  then `PETSC_ARCH=complex-debug make -f Makefile_CS`, then
  `pip install . --no-build-isolation`. The `PETSC_ARCH=complex-debug`
  export is mandatory — the default `real-debug` (or stale `.mod` from a
  prior real-arch attempt) yields spurious COMPLEX→REAL interface errors
  in `fortranPC.F90`/`adjointAPI.F90`/`surfaceUtils.F90`/etc. that are
  toolchain artifacts, not source problems.
- SA-GR's higher-level `evalFunctionsSens`-based adjoint/CS validation on
  AR5 still not started — tracked in `../current-task.md`.
