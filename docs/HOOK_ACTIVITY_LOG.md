# Hook Activity Log

Auto-appended by the project hooks in `.claude/hooks/` — newest entries at the
bottom. Two hooks write here:

- **auto-pip** (`auto_pip_after_make.sh`, PostToolUse/Bash): logs each time it
  reinstalls adflow into the mach env after a `make`, so site-packages never
  drifts from `./adflow` (the stale-install trap that broke the gammaForSA
  clamp task).
- **wrap-up** (`propose_wrapup.sh`, Stop): logs each time it prompts for a
  doc-update + commit/push at task end, with a snapshot of the uncommitted
  task-relevant change-set.

---
- **2026-07-23 18:45:44** — wrap-up prompt fired (change-set 0ceeb08fc831); uncommitted task-relevant paths:
    ?? docs/HOOK_ACTIVITY_LOG.md
