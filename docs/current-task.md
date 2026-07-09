# Current Task

> Update this file's content whenever the active task changes (rule: **one
> task per session**, CLAUDE.md #1). When a task finishes, move its summary
> into a new file under `task-log/` (see `task-log/README.md` for the template and index) and replace this file's body with
> the next task, or with "No task in progress" if idle.

## Task: Verify SA-GR adjoint derivatives (partials campaign)

**Started:** 2026-07-09
**Status:** in progress

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
