# Hook Activity Log

Auto-appended by the project hooks in `.claude/hooks/` — newest entries at the
bottom. Three hooks write here:

- **auto-pip** (`auto_pip_after_make.sh`, PostToolUse/Bash): logs each time it
  reinstalls adflow into the mach env after a `make`, so site-packages never
  drifts from `./adflow` (the stale-install trap that broke the gammaForSA
  clamp task).
- **wrap-up** (`propose_wrapup.sh`, Stop): logs each time it prompts for a
  doc-update + commit/push at task end, with a snapshot of the uncommitted
  task-relevant change-set.
- **guard** (`guard_protected_files.sh`, PreToolUse/Edit|Write|NotebookEdit):
  logs each time it hard-blocks an edit to the SA model `src/turbulence/sa.F90`
  (CLAUDE.md rule 2).

---
- **2026-07-23 18:45:44** — wrap-up prompt fired (change-set 0ceeb08fc831); uncommitted task-relevant paths:
    ?? docs/HOOK_ACTIVITY_LOG.md
- **2026-07-23 19:15:00** — wrap-up prompt fired (change-set 7a8095be7783); uncommitted task-relevant paths:
     M CLAUDE.md
     M docs/HOOK_ACTIVITY_LOG.md
     M tests/reg_tests/reg_sagr.py
     M tests/reg_tests/run_sagr_tests.sh
     M tests/reg_tests/test_adjoint_sagr.py
