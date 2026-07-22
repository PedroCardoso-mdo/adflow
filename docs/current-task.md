# Current Task

> Update this file's content whenever the active task changes (rule: **one
> task per session**, CLAUDE.md #1). When a task finishes, move its summary
> into a new file under `task-log/` (see `task-log/README.md` for the
> template and index) and replace this file's body with the next task, or
> with "No task in progress" if idle.

## Active Task: SA / SA-GR derivative verification on the AR5 mesh

**Started:** 2026-07-20
**Status:** in progress — validated through Step 2 (aero DVs); Step 3 (AR5)
blocked on a convergence-robustness issue, not a derivative bug.

**PARTIALS DONE (2026-07-22).** The raw-API 3-stage verification ladder
(dot-product, fast-reverse, 3-way AD/FD/CS) is now a **registered testflo
suite** driven by `tests/reg_tests/run_sagr_tests.sh`, passing for plain SA
and SA-GR (`nw=8`) on AR5 with **crossflow always ON** (real 13/13, CS 10/10).
So the *partials* — `dR/dw` (all SA↔transition coupling blocks incl. the
D_scf crossflow block), `dR/dXv`, `dR/d{aeroDV}`, output partials — are
validated against complex-step in all three AD modes. See
`VERIFICATION/three-stage-verification.md` §"Canonical way to run" and
`task-log/2026-07-22-sagr-test-suite-standardization.md`.

**What remains (this task): the TOTAL sensitivity `dF/dX`** through the
assembled adjoint solve — the "complete-mode" `test_adjoint`. That is a
*different* thing from the validated partials and is still blocked below,
because this AR5 state does not converge deeply enough (absolute
precision floor). Validated partials ≠ validated gradient.

### Objective

Validate that SA (then SA-GR) derivatives are correct on the AR5 plain-wing
mesh + AR5 FFD (`input_files/ar5_plain_wing_vol_L3.cgns`,
`ar5_plain_wing_ffd_L3.xyz`), using the same rigorous methodology already
trusted on the tutorial-wing case. 3-step process, each via the adjoint
(real, vs. trained reference) + complex-step (vs. adjoint) route:

1. Baseline: official `test_adjoint.py` class, tutorial-wing mesh.
2. Same methodology, my own test file, tutorial-wing mesh (proves my code
   matches the baseline before trusting it on a new mesh).
3. Same methodology, my own test file, AR5 mesh (the actual goal).

### Tests that have passed

- `test_adjoint.py:TestAdjoint_3_rans_tut_wing` — 4/4 (real adjoint totals
  vs. trained reference, tutorial-wing, plain SA).
- `test_adjoint.py:TestCmplxStep_3_rans_tut_wing` — 2/2 (CS vs. adjoint
  totals, tutorial-wing, plain SA, incl. FFD shape DV).
- `test_adjoint_tutwing_mycode.py:TestAdjointTutWingMyCode` — real adjoint,
  my own test file, tutorial-wing. Matches the baseline.
- `test_adjoint_tutwing_mycode.py:TestCmplxStepTutWingMyCode.cmplx_test_aero_dvs`
  — CS vs. adjoint for alpha/mach, tutorial-wing, my own test file.
- `test_adjoint_ar5_sa.py:TestAdjointAR5SA` — real adjoint, AR5 mesh, plain
  SA.

### Still open (pick up here)

- `cmplx_test_geom_dvs` (FFD shape DV) on tutorial-wing (Step 2) fails on
  the `cavitation` functional specifically (98.6% mismatch) — cl/cd/cmz not
  yet confirmed to pass in isolation. Cavitation uses indicator/threshold
  functions, plausibly CS-unfriendly by nature — isolate before concluding
  anything.
- **Step 3 (AR5) CS check fails hard**: `dcd/dmach` CS=-231.1 vs.
  adjoint=-0.537 (~429x off). Root cause **narrowed down 2026-07-20** via a
  real-build FD step-size sweep (`sweep_h_fd.py`, ADPC on, both `+h`/`-h`
  solves confirmed converged, `fail=False`):

  | h | FD dcd/dmach | notes |
  |---|---|---|
  | 1e-3 | -1.73 | `+h` solve failed to converge |
  | **1e-4** | **-0.666** | both converged; closest to adjoint's -0.537 |
  | 1e-5 | +19.4 | both "converged" but nonsensical |
  | 1e-6 | -1.96 | both "converged" but still off |

  Non-monotonic, erratic vs. h — the classic signature of an **absolute
  precision floor**, not a derivative bug: this mesh's residual creeps very
  slowly (`Step` pinned ~0.01, the same quasi-stall documented elsewhere on
  this branch) and never reaches machine-zero within a practical `ncycles`
  budget, so `cd` itself is only stable to ~6 significant digits. For
  `h=1e-5`/`1e-6` the expected FD signal (~5e-6/5e-7 in `cd`) is at or below
  that noise floor; for `h=1e-4` the signal (~5e-5) rises above it, giving
  the much saner -0.666. **Working hypothesis: the adjoint (-0.537) is
  correct; the mesh just can't be converged tightly enough, in absolute
  terms, for small-h FD or CS (h=1e-40) to resolve the true derivative** —
  CS needs even more absolute precision than FD, which is consistent with
  CS failing even harder (-231 vs. FD's -1.7 to +19).

  **Tried and ruled out (2026-07-20): "just give it more cycles" does not
  fix this.** Raising `ncycles` for `h=1e-4` from 5000 (both sides
  converged, FD=-0.666, closest match) to 15000 made the `+h` side
  **fail** (`fail=True`) and land on a *worse* FD value (-1.406), despite
  identical `h`, options, and `l2convergence=1e-14` (already tight — the
  target isn't the bottleneck, reaching it is). Same case, same options,
  different outcome with more cycles: this is genuine, chaotic
  floating-point sensitivity right at this stall's knife-edge (same
  phenomenon documented in `nk_switch_crossing_test/PURPOSE.md`), not a
  simple iteration-budget problem. Pushing `ncycles` further is not
  expected to help reliably.

  **Where this leaves the task**: the adjoint (-0.537) is trusted (matches
  the already-validated methodology from Steps 1-2). FD/CS validation of
  this *specific* AR5+SA state via `resetFlow`+re-solve is not currently
  practical because the perturbed re-solve lands in/near the same
  chaotic quasi-stall this branch has documented repeatedly. Real fix
  would be a genuine convergence-robustness improvement for this stall
  (see `nk_switch_crossing_test/run_stall_fix_organic.log` — unvalidated
  candidate fix exists, `transitionNKStallRtolCap`, but that's gated to
  SA-GR/`transitionNK` only, not plain SA) — a separate, bigger task, not
  a quick follow-up.
- **DONE (2026-07-21):** the raw, low-level rung-1 check flagged above
  (large mismatch on the nuTilde residual row) was root-caused and fixed,
  then re-verified via a proper 3-stage ladder (reverse↔forward
  dot-product consistency, reverse vs. fast-reverse consistency, 3-way
  AD/FD/CS forward check) — **all three stages PASS**, both meshes, plain
  SA. Full writeup: `docs/VERIFICATION/three-stage-verification.md`;
  case file: `task-log/2026-07-21-three-stage-adjoint-verification.md`.
  This does not touch the AR5 `dcd/dmach` CS mismatch below, which is a
  separate, converged-state-dependent issue (see that doc's "Relationship
  to the still-open AR5 CS mismatch" section).
- **DONE (2026-07-21):** ran the 3-stage ladder for SA-GR (`nw=8`,
  gamma/reThetat rows), AR5 mesh, via a new `--turbmodel sagr` flag on
  `sanity_check_partials_sa.py`/`check_3way_fwd.py`/`check_3way_fwd_sweep.py`.
  **Stage 1 (dot products) PASS; Stage 3 (AD/FD/CS) PASS** (AD = CS =
  3.3256409775e+13 exact); **Stage 2 (`_b` vs `_fast_b`) FAILS on the
  reThetat rows** (max\|_b − _fast_b\| ≈ 7.5e+2, ~44% elements mismatched;
  meanflow/nuTilde/gamma pass). Full writeup:
  `VERIFICATION/three-stage-verification.md` ("SA-GR (`nw=8`) results");
  case file: `task-log/2026-07-21-three-stage-adjoint-verification-sagr.md`.
- **DONE (2026-07-22): Stage-2 `_b`-vs-`_fast_b` reThetat divergence
  root-caused and FIXED.** Ground-truth check proved `_fast_b` (not `_b`)
  was wrong for `dR[reThetat]/dw[meanflow]` (rel_err 1.94, wrong sign).
  Cause: the primal computed `lambdaThetaLocal` with two **in-place**
  `smoothMinMax` clamps; Tapenade's `push/popreal8` for the overwritten
  intermediates is stripped by `autoEditReverseFast.py`, leaving `_fast_b`
  reading a stale primal — the sole in-place `smoothMinMax` in the source
  (all others use distinct targets). Fix: rewrote the primal with distinct
  targets (`lambdaThetaRaw`→`lambdaThetaClamped`→`lambdaThetaLocal`,
  `saGammaRetheta.F90`), reran Tapenade (regenerated `_d`/`_b`/`_fast_b`)
  and rebuilt. Now all three stages PASS for SA-GR (reThetat block
  `_b`-vs-`_fast_b` 7.5e+2 → 2.3e-10; `_fast_b` matches forward/CS ground
  truth to 1.27e-13; Stage 3 AD=CS=3.3256409775e+13 unchanged). Case file:
  `task-log/2026-07-22-fastb-retheta-lambdatheta-fix.md`.
- SA-GR's own adjoint/CS validation on the AR5 mesh (the higher-level
  `evalFunctionsSens`-based route, distinct from the raw-API ladder above)
  hasn't started — `test_adjoint_sagr.py`/`reg_sagr.py` were repointed at
  AR5 earlier but never run through this same 3-step process.

### Key files

- `tests/reg_tests/test_adjoint_ar5_sa.py` — Step 3 (real: `TestAdjointAR5SA`;
  CS: `TestCmplxStepAR5SA`).
- `tests/reg_tests/test_adjoint_tutwing_mycode.py` — Step 2.
- `tests/reg_tests/generate_sagr_restart.py`, `reg_sagr.py` — SA-GR side of
  the same AR5 case switch (grid/FFD/options done; not yet validated).

### Lesson learned this session (see memory: `nondeterminism-scope-limited.md`)

Do not attribute derivative-check discrepancies to floating-point
non-determinism. Every real discrepancy found so far had a concrete, fixable
cause: a stale `resDot` reused after a state-perturbing FD probe, a missing
`defaultAeroDVs` entry, a mismatched `N_PROCS`. Keep debugging until a
concrete cause is found or ruled out.
