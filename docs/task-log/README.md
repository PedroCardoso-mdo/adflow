# Finished Task Log — SA-γ-Re̅θt

One file per finished task (a "case"), added here over time as work closes —
not appended into one ever-growing file. When a task in `../current-task.md`
finishes: write `<YYYY-MM-DD>-<slug>.md` in this directory using the
template below, add a row to the index, then reset `../current-task.md` to
the next task.

## Index

| Date | Task | File |
|---|---|---|
| 2026-07-09 | Doc consolidation (`findings/` redistribution) | [2026-07-09-doc-consolidation.md](2026-07-09-doc-consolidation.md) |
| 2026-07-21 | Three-stage low-level adjoint verification (plain SA) | [2026-07-21-three-stage-adjoint-verification.md](2026-07-21-three-stage-adjoint-verification.md) |
| 2026-07-21 | Three-stage low-level adjoint verification (SA-GR, nw=8) — Stages 1&3 PASS, Stage 2 FAILS on reThetat | [2026-07-21-three-stage-adjoint-verification-sagr.md](2026-07-21-three-stage-adjoint-verification-sagr.md) |
| 2026-07-22 | Fix reverse-fast `dR[reThetat]/dw[meanflow]` (lambdaTheta in-place clamp) — Stage 2 SA-GR now PASSES | [2026-07-22-fastb-retheta-lambdatheta-fix.md](2026-07-22-fastb-retheta-lambdatheta-fix.md) |

## Template

```
# <task name> — <YYYY-MM-DD>

**Problem:** what was broken/missing/in question, 1-3 sentences.

**Solution:** what was actually changed (files:lines), 1-5 bullets.

**Files created/touched and why:**
- `path/to/file` — purpose in this task, and why it's structured the way
  it is (strategy), not just what it contains.

**Verification:** what was run to confirm done-ness (compiles, smoke test
N iters no NaN, TAPENADE NEEDED?, etc.) per CLAUDE.md's definition of done.

**Follow-ups:** anything punted to `../TODO.md`, with the item name.
```
