# Current Task

> Update this file's content whenever the active task changes (rule: **one
> task per session**, CLAUDE.md #1). When a task finishes, move its summary
> into a new file under `task-log/` (see `task-log/README.md` for the
> template and index) and replace this file's body with the next task, or
> with "No task in progress" if idle.

## Active Task: gammaForSA clamp AD-vs-CS mismatch (tutorial-wing/M=0.15 case)

**Started:** 2026-07-23. **Status: RESOLVED 2026-07-23.** The full SA-GR
partial suite now passes: `./run_sagr_tests.sh all` → real (Stage 1
dot-products, Stage 2 `_b`↔`_fast_b`, Stage 3 AD/FD) **PASS** and complex
(Stage 3 CS ground truth) **PASS**.

**Root cause was a build/install gotcha, NOT a physics/AD bug.** The
`rsaGRgammaForSAMargin=1e-3` clamp fix (both `Source` ~L533 and
`evalSrcJacBlock` ~L2533 in `saGammaRetheta.F90`) plus the regenerated
Tapenade files were correct and *were* compiled into the freshly-built libs
in the repo tree (`./adflow/libadflow{,_cs}.so`). But those libs were never
`pip install`ed into the mach env, and **the tests import `adflow` from
site-packages** (their CWD is `tests/reg_tests/`, not the repo root, so
`./adflow` is not on `sys.path`). Site-packages still held the pre-fix
15:30 build with the *bare* clamp `min(max(gamma,0),1)`. So every test run
— across sessions — executed the stale bare-clamp binary and produced the
*identical* CS=0-at-`gamma==1` failure, which read as "the fix didn't take."
Confirmed decisively: a runtime probe bumping `gamma==1` cells to 1.0005
(strictly between bare bound 1.0 and margin bound 1.001) gave `dR/dgamma==0`
under the site-packages lib (bare) and `dR/dgamma!=0` under the repo lib via
`PYTHONPATH` (margin).

**Fix applied:** `pip install . --no-deps` into `/home/mdo/packages_v2/mach`
(refreshes site-packages `.so`s to the margin build), then
`./run_sagr_tests.sh train` to regenerate the JSON refs against the fixed
libs (the old refs were trained on the stale build; the dot-product test
passing after retrain confirms the forward `_d` = reverse `_b` identity
holds — no `_d`/`_b` inconsistency). Refs updated:
`refs/jacvec{fwd,bwd}_sagr_flatplate.json`.

**Standing lesson (see memory):** after any Fortran rebuild, `pip install`
into the mach env before running `tests/reg_tests/` — building `./adflow`
alone does not reach the tests. See [[tests-load-site-packages-not-repo]].

<details><summary>Original investigation notes (pre-resolution)</summary>

Condensed summary:

- **Case switch (2026-07-23):** SA-GR test suite (`reg_sagr.py`,
  `dev/generate_sagr_restart.py`) moved from the AR5 mesh (chronic
  quasi-stall, see the superseded task below) to the standard ADflow
  tutorial-wing mesh at Mach=0.15 (mid-transition, not full-laminar/
  full-turbulent), which converges cleanly to `L2Convergence` with no
  stall. New restart: `input_files/mdo_tutorial_sagr_dp.cgns`.
- **Bug found (confirmed via per-cell AD-vs-CS diffing):**
  `dR[nuTilde]/dw[gamma]` mismatches CS by ~2e-8 relative, concentrated
  100% at cells where `gamma==1.0` bit-exact.
  `saGammaRetheta.F90`'s `gammaForSA = min(max(gamma,0),1)` (2 occurrences:
  `Source` ~line 533, PC-only `evalSrcJacBlock` ~line 2529) ties exactly at
  gamma's own natural saturation values (0/1), which real converged states
  hit routinely — unlike the file's other clamps (`rsaGRgammaLo/Hi`,
  `rsaGRreThetaLo`), which are arbitrary safety floors never hit exactly
  (verified). Tapenade's tangent picks the wrong branch at the tie.
- **Fix (per user's explicit instruction, NOT smoothed):** padded the clamp
  bounds via a new `rsaGRgammaForSAMargin=1e-3` constant
  (`paramTurb.F90`) so gamma=0/1 sit strictly inside the unclamped
  pass-through region — a pure divergence safeguard, not a physics change.
  Applied to both occurrences.
- **Tapenade rerun + rebuild:** done, with the user's explicit one-time
  go-ahead. Regenerated exactly `outputForward/Reverse/ReverseFast
  saGammaRetheta_{d,b,fast_b}.f90`. **First incremental `make -j` rebuild
  was a red herring** — looked successful but per-cell diffing showed zero
  behavior change; a **full clean rebuild** (`make clean && make -j` real,
  `Makefile_CS clean` + rebuild complex) was required to actually pick up
  the change. No `git clean` needed for either.
- **Still unresolved:** after the clean rebuild + retrained JSON refs, the
  CS test still fails with numbers that print identically to before. Some
  evidence suggests AD's tangent output changed at gamma==1 cells, but this
  was not rigorously confirmed cell-for-cell. The specific worst cell
  (rank 1, cell index 2420 in the `wDot=getStatePerturbation(314)` masked-
  to-gamma-column setup) still shows AD nonzero, CS exactly 0.0 even after
  the clean rebuild — unexplained. Leading hypothesis (unverified): something
  in how the complex build's `min`/`max` (via the `complexify` library) or
  its primal evaluation of the padded clamp handles this construct
  differently from the real build. Next step: isolate with `gamma=1.0005`
  vs `gamma=1.05` synthetic states to directly test clamp branch selection
  in both builds independent of the full solve.

</details>

---

## Previous task (superseded 2026-07-23 by the mesh switch above): SA / SA-GR derivative verification on the AR5 mesh

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
