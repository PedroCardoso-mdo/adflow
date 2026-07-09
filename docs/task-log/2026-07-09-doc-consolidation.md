# Doc consolidation (findings/ redistribution) — 2026-07-09

**Problem:** `docs/` had grown too many small files; `docs/findings/`
(`A1` folded into `A_confirmacao.md`, `A2_convergencia.md`,
`A3_coerencia.md`, `A_confirmacao.md`, `D1_transitionRefLength.md`)
duplicated content already merged into `architecture.md`/
`nondimensionalization.md`, and scattered the rest across 4 overlapping
files with no single "decisions" reference. `task-log.md` was also a single
ever-growing file rather than one case per finished task.

**Solution:**
- Folded all still-unique content from `findings/*` into
  `docs/audits/design-decisions.md`, organized by topic instead of by
  audit-task ID, with the already-redundant `D1_transitionRefLength`
  content dropped (fully covered by `architecture.md` Part 2 +
  `nondimensionalization.md` §5 already). Placed under `audits/` since
  that's where the underlying investigation happened, and marked
  explicitly as **not a defining/normative file** — a memory of past
  discussion, not a spec; code/paper win if it ever disagrees.
  `docs/audits/{00_inventory,06_adjoint_wiring,07_sst_dev_lessons,
  08_test_prep,adjoint_audit_2026-07-07}.md` were left untouched (out of
  scope).
- Updated cross-references in `TODO.md`, `architecture.md`,
  `nondimensionalization.md`, `README.md` to point at
  `audits/design-decisions.md` instead of `findings/`.
- Deleted `docs/findings/`.
- Converted the finished-task log from one flat file into a directory
  (`docs/task-log/`) with an index (`README.md`) and one file per task —
  this file is the first case.

**Files created/touched and why:**
- `docs/audits/design-decisions.md` — single canonical log of resolved
  audit questions (A1-A3), condensed ~3x from the source material;
  organized by theme (nondim, convergence strategy, code coherence) so a
  future reader can find "why is X the way it is" without knowing which
  audit task found it. Explicitly a memory/log, not a source of truth.
- `docs/current-task.md` — always holds exactly the one task in progress,
  per CLAUDE.md's one-task-per-session rule; meant to be overwritten each
  session rather than accumulated.
- `docs/task-log/README.md` — index + template for the finished-task log.
- `docs/task-log/<date>-<slug>.md` (this file) — one case per finished
  task, so the log grows by adding files, not by growing one file
  indefinitely (the same sprawl problem `findings/` had, at a smaller
  scale).

**Verification:** doc-only change, no code/build impact.

**Follow-ups:** none.
