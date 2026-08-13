# Current Task

> Update this file's content whenever the active task changes (rule: **one
> task per session**, CLAUDE.md #1). When a task finishes, move its summary
> into a new file under `task-log/` (see `task-log/README.md` for the
> template and index) and replace this file's body with the next task, or
> with "No task in progress" if idle.

## No task in progress (as of 2026-08-12)

The previous "Active Task" (gammaForSA clamp AD-vs-CS mismatch) was RESOLVED
on 2026-07-23 and its full record lives in
[`task-log/2026-07-23-sagr-full-adjoint-test.md`](../task-log/2026-07-23-sagr-full-adjoint-test.md)
and [`../task-log/2026-07-24-sagr-cs-tolerance-nonblocking.md`](../task-log/2026-07-24-sagr-cs-tolerance-nonblocking.md)
(CS tolerance settled at 5e-8, mach/drag non-blocking). This file had
retained the resolved task plus ~140 lines of a superseded July partials
campaign; that content was already duplicated by the task-log entries and
`VERIFICATION/three-stage-verification.md`, and was removed 2026-08-12
(recoverable from git history).

## Where the post-July work is recorded (no task-log entries were written)

| Work (date) | Record |
|---|---|
| Crossflow default flip → False + psi-history adjoint diagnostic (2026-07-27, `a34441c9`) | `architecture.md` Part 2, `VERIFICATION/adjoint-trace.md` |
| λθ-clamp switchable, unclamped tested & REJECTED, default reverted (2026-08-03, `be9d6d1d`/`efed31cf`) | `VERIFICATION/paper-validation-campaign.md` |
| Matrix-dissipation eigenvalue limiters Vn/Vl + Tapenade regen (2026-08-04, `2c4ce2c1`) | `architecture.md` §"Matrix-dissipation eigenvalue limiters" |
| S809/NLF0416 paper-validation campaign (2026-08-03/06, paused) | `VERIFICATION/paper-validation-campaign.md` |
| AR5 corrected-foil refinement: premature-NK rule, EW mitigation falsified, ILU study (2026-08-07/12, `ef9fc10d`/`54c43475`) | `convergence-strategy.md` |

## Known open items (pointers, not a plan)

- Deep-NK preconditioner wall — code items in
  `SA_GAMMA_RETHETHA_BASE/adflow-vs-paper-solver.md` §8; option-level levers
  exhausted (`convergence-strategy.md`).
- `.pyf` no-ops: `ANKPhysicalLSTolReTheta` / `omegaMinGamma`
  (`architecture.md`, Turb-ANK options).
- Mach CS derivative doesn't settle (non-blocking) —
  `../task-log/2026-07-24-sagr-cs-tolerance-nonblocking.md`.
- Rotating-frame path not exercised by any test —
  `VERIFICATION/rotating-frame-audit.md` §5.
- Physics gap vs paper (+20/30 counts above bucket) — campaign paused,
  next step defined in `VERIFICATION/paper-validation-campaign.md`.
