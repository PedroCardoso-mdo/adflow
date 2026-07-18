# Current Task

> Update this file's content whenever the active task changes (rule: **one
> task per session**, CLAUDE.md #1). When a task finishes, move its summary
> into a new file under `task-log/` (see `task-log/README.md` for the template and index) and replace this file's body with
> the next task, or with "No task in progress" if idle.

## Active Task: NK paper-faithful convergence push (Eq. 59, Algorithm 2, Eq. 58 S_r/S_a)

**Started:** 2026-07-16
**Status:** implementation + smoke verification done; production run next

### Objective

Get NK as close as possible to the paper's (Piotrowski & Zingg 2020) §IV.B
Newton–Krylov–Schur algorithm, then run a real production convergence
attempt via the validated `best_strategie/run_strategy.py nk` stage. This
directory (`sa_gamma_rethetha_paper_solver`) is dedicated to closing the
open items in `adflow-vs-paper-solver.md` §6. This pauses (does not
replace) the adjoint-verification task below.

### Context (read only what's listed — CLAUDE.md rule 8)

- `docs/adflow-vs-paper-solver.md` §2, §4, §5, §6 — mechanism detail for
  everything below, the S_r/S_a smoke-test results, and the `.pyf` bug.
- `docs/architecture.md` Part 2 — every new option's exact semantics.

### What was done

1. **Eq. 59 reactivation-on-backtrack inside NK** — `NKStep` updates
   `noBacktrackCount` after the line search; `FormJacobianNK` +
   new `applyNKSrcDtDiagonal` inject the source-dt diagonal onto the
   already-assembled, already column-scaled `dRdwPre` (never touches
   `setupStateResidualMatrix`, shared with the frozen adjoint).
2. **Faithful Algorithm 2 inside NK** — new `applyNKAlgorithm2Damping`,
   called from `NKStep` after the line search accepts `work`: per-node
   bounds-triggered exponential back-off on γ/Re̅θt only, reusing DD-ADI's
   existing constants/options.
3. **Eq. 58 S_r geometric row scaling** — new factor in `setRVec`
   (`volRef**(5/3)` flow rows, `volRef**(2/3)` turb rows), gated by new
   option `transitionRowVolScale` (default `False`).
4. **Eq. 58 S_a autoscale proxy** — new `computeNKResidualAutoscale`
   (called once per Jacobian reform) + `nkAutoScaleFac` (moved to
   `inputIteration` so `utils.F90` can share it without a circular
   `use`), gated by new option `transitionResidualAutoscale` (default
   `False`). Not verified against Osusky & Zingg's actual method (no
   published formula available) — a same-intent proxy only.
5. New monitor variable `"scaledtotalr"` (`sumAllResidualsScaled`,
   `utils.F90`) exposing the S_r/S_a-scaled residual for visibility;
   deliberately does not feed `totalR`/switch tolerances.
6. New master switch `transitionNK` (default `True`) gating the whole
   column-scaling + Eq. 59 bundle across NK/ANK/turbKSP, replacing bare
   `turbModel==SA-GR` checks — `.and. transitionNK` added at all 14 gate
   sites (verified via grep only the gamma/Re̅θt physicality-bounds
   checks remain unguarded, correctly).
7. **Found and partially fixed a real infra bug**: `src/f2py/adflow.pyf`
   is hand-maintained, not auto-regenerated, and was stale — new
   `inputIteration` variables were invisible to real Python↔Fortran
   communication despite `setOption`/read-back appearing to work (f2py's
   `fortran` objects silently accept phantom attribute names). Fixed for
   `transitionNK`/`transitionRowVolScale`/`transitionResidualAutoscale`
   (added to `.pyf`, verified via Fortran-side debug print). **Not**
   fixed for pre-existing `transitionSrcDtRestrict` and siblings — out of
   this session's scope, flagged in `adflow-vs-paper-solver.md` §5/§6 and
   `architecture.md`.

### Smoke test results (`nk_eq59_reactivation_test/`, 10-iteration windows)

- Items 1+2 (Eq. 59 + Algorithm 2): clean, matches/improves on baseline,
  no NaN. `run_10iter_step1_eq59_transitionNK.log`,
  `run_10iter_step2_algorithm2.log`.
- Item 3 (S_r), once genuinely wired: **severe linear-solve stall** (lin
  res pinned ~1.0 for 9 iterations). `run_10iter_step5_Sr_Sa_genuinely_wired_STALL.log`.
  Stays off by default.
- Item 4 (S_a alone), once genuinely wired: no stall, real progress,
  noisier steps than baseline, not clearly better.
  `run_10iter_step4_Sa_alone_genuinely_wired.log`. Stays off by default.
- Earlier "positive" results for items 3/4 (`run_10iter_step3_Sr_geometric_ON.log`,
  `run_10iter_step4_Sa_autoscale_ON.log`) predate the `.pyf` fix and
  reflect **floating-point recompilation noise, not real effects** —
  kept in the folder for the record but not evidence of anything.

### Known simplification vs. ANK

ANK's coupled-path Eq. 59 reactivation also resets on
`totalR > ANK_secondOrdSwitchTol * totalR0`. `ANK_secondOrdSwitchTol` is
module-local inside `module ANKSolver`, which `use`s `module NKSolver` —
not the reverse — so `NKStep` cannot see it without a circular `use` or
relocating the variable (deferred; a residual rise already shows up as a
backtracked step in practice).

### Definition of done for this task

Per CLAUDE.md: code compiles, smoke test runs N iterations without NaN —
satisfied for items 1-4. Remaining: production run via
`best_strategie/run_strategy.py nk` with the validated-safe options only
(`transitionNK` default True; S_r/S_a left off), logged to
`best_strategie/logs/`. Report the log/final residual as facts (vs. the
documented 3.3e-9 record / ~5e-9 saturation point) — no pass/fail
interpretation from me, per CLAUDE.md rule 7.

When done, log to `task-log/` per its template and update the index; then
resume the adjoint task below, or start the next open item (§6: fix the
`.pyf` bug for the remaining options, or investigate S_r's stall) in a
fresh session.

---

## Paused Task: Verify SA-GR adjoint derivatives (partials campaign)

**Started:** 2026-07-09
**Status:** paused (see Active Task above)

### Objective

Confirm the SA-γ-Re̅θt adjoint (`db_forward`/`_d`, `_b`, `_fast_b` Tapenade
output) produces correct derivatives before treating the adjoint as done.
This is the verification step the A4 audit (`audits/design-decisions.md`,
`audits/adjoint_audit_2026-07-07.md`) deferred — code review only, no
partials were actually run yet.

### Context (read only what's listed — CLAUDE.md rule 8)

- `docs/adjoint-trace.md` — adjoint/AD touchpoints on this branch.
- `docs/audits/adjoint_audit_2026-07-07.md` — pre-partials visual audit,
  the `vortlimd=0` finding, and the watch items below.
- `docs/audits/sst_dev_lessons.md` — SST post-mortem: what a verification
  campaign missed last time (rtol inflation masking a real bug, cd/cm/DVs
  not checked, only cl).
- `docs/TODO.md` §"Adjoint / partials" — the checklist this task works off.

### Working files (this session, `tests/reg_tests/`)

- `generate_sagr_restart.py` — produces a converged SA-GR restart solution
  to use as the base state for partials (dot-product/FD/CS tests need a
  physically converged point, not a fresh initialization).
- `reg_sagr.py` — registration/harness glue for the SA-GR regression case.
- `test_adjoint_sagr.py` — dR/dw and dR/dx adjoint correctness tests.
- `test_jacVecProdFWD_sagr.py` — forward-mode (`_d`) Jacobian-vector product
  tests.
- `test_jacVecProdBWDFast_sagr.py` — reverse-mode-fast (`_b_fast`) tests;
  per `TODO.md` item (b), if this one fails, suspect
  `autoEditReverseFast.py` push/pop stripping first (this broke SST's
  fast_b upstream, and it's still active on this branch).
- `README_SAGR.md` — how to run the above and what each checks.

### Checklist (from `TODO.md`, don't duplicate detail here — read it there)

- [ ] Rerun Tapenade to pick up `uInf`/`muInf` as active in the
      `saGammaRetheta%Source` head (currently `vortlimd` hard-zero).
- [ ] Decide frozen vs. differentiated vorticity limiter cap (default =
      differentiate).
- [ ] Run dot-product test; if BWDFast fails, check `autoEditReverseFast.py`
      stripping before anything else.
- [ ] Expect FD noise (not AD error) near `smoothMinMax`/vortLim blend
      points — use complex-step or a mask there, never inflate global rtol
      (this is exactly what went wrong in the SST post-mortem).
- [ ] Validate cd/cm and all DVs with CS, not just cl.

### Definition of done for this task

Per CLAUDE.md: adjoint dot-product/FD/CS tests pass (or documented,
understood failures with a fix plan) — physics correctness beyond
derivative correctness is not in scope. When done, log it as a new case file in `task-log/` and add it to `task-log/README.md`'s index.
