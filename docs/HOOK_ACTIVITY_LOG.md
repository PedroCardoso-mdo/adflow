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
- **2026-07-23 19:27:38** — wrap-up prompt fired (change-set 275b966324ca); uncommitted task-relevant paths:
    M  CLAUDE.md
    M  docs/HOOK_ACTIVITY_LOG.md
     M src/turbulence/saGammaRetheta.F90
     M tests/reg_tests/README_SAGR.md
     M tests/reg_tests/dev/diag_all16_blocks_cs.py
    R  tests/reg_tests/refs/jacvecbwd_sagr_flatplate.json -> tests/reg_tests/refs/jacvecbwd_sagr_tut_wing.json
    R  tests/reg_tests/refs/jacvecfwd_sagr_flatplate.json -> tests/reg_tests/refs/jacvecfwd_sagr_tut_wing.json
     M tests/reg_tests/reg_sagr.py
     M tests/reg_tests/run_sagr_tests.sh
     M tests/reg_tests/test_adjoint_sagr.py
     M tests/reg_tests/test_jacVecProdBWDFast_sagr.py
     M tests/reg_tests/test_jacVecProdFWD_sagr.py
- **2026-07-23 19:31:51** — wrap-up prompt fired (change-set d4a9587c7f1e); uncommitted task-relevant paths:
    M  CLAUDE.md
    MM docs/HOOK_ACTIVITY_LOG.md
     M docs/current-task.md
     M docs/task-log/README.md
     M src/turbulence/saGammaRetheta.F90
     M tests/reg_tests/README_SAGR.md
     M tests/reg_tests/dev/diag_all16_blocks_cs.py
    R  tests/reg_tests/refs/jacvecbwd_sagr_flatplate.json -> tests/reg_tests/refs/jacvecbwd_sagr_tut_wing.json
    R  tests/reg_tests/refs/jacvecfwd_sagr_flatplate.json -> tests/reg_tests/refs/jacvecfwd_sagr_tut_wing.json
     M tests/reg_tests/reg_sagr.py
     M tests/reg_tests/run_sagr_tests.sh
     M tests/reg_tests/test_adjoint_sagr.py
     M tests/reg_tests/test_jacVecProdBWDFast_sagr.py
     M tests/reg_tests/test_jacVecProdFWD_sagr.py
    ?? docs/task-log/2026-07-23-sagr-full-adjoint-test.md
- **2026-07-23 19:38:23** — wrap-up prompt fired (change-set 0191b1a7c039); uncommitted task-relevant paths:
     M docs/HOOK_ACTIVITY_LOG.md
     M docs/current-task.md
     M docs/task-log/README.md
     M src/turbulence/saGammaRetheta.F90
     M tests/reg_tests/README_SAGR.md
     M tests/reg_tests/dev/diag_all16_blocks_cs.py
     M tests/reg_tests/reg_sagr.py
     M tests/reg_tests/run_sagr_tests.sh
     M tests/reg_tests/test_adjoint_sagr.py
     M tests/reg_tests/test_jacVecProdBWDFast_sagr.py
     M tests/reg_tests/test_jacVecProdFWD_sagr.py
    ?? docs/task-log/2026-07-23-sagr-full-adjoint-test.md
- **2026-07-23 19:40:07** — wrap-up prompt fired (change-set 0933e396cac1); uncommitted task-relevant paths:
     M docs/HOOK_ACTIVITY_LOG.md
     M docs/current-task.md
     M docs/task-log/README.md
     M src/turbulence/saGammaRetheta.F90
     M tests/reg_tests/README_SAGR.md
     M tests/reg_tests/dev/diag_all16_blocks_cs.py
     M tests/reg_tests/reg_sagr.py
     M tests/reg_tests/run_sagr_tests.sh
     M tests/reg_tests/test_adjoint_sagr.py
     M tests/reg_tests/test_jacVecProdBWDFast_sagr.py
     M tests/reg_tests/test_jacVecProdFWD_sagr.py
    ?? docs/task-log/2026-07-23-sagr-full-adjoint-test.md
    ?? src/adjoint/preprocess_forward/
- **2026-07-23 19:45:47** — auto-pip: reinstalled adflow (complex-step (CS) build) into mach env after `cd /home/mdo/MDOLab_3_v2/adflow_sa_gamma_rethetha_paper_solver echo "=== auto-pip: syntax check ==="` (site-packages now matches ./adflow)
- **2026-07-23 19:46:09** — wrap-up prompt fired (change-set 84836aa8b066); uncommitted task-relevant paths:
     M docs/current-task.md
     M docs/task-log/README.md
     M src/adjoint/outputForward/saGammaRetheta_d.f90
     M src/adjoint/outputReverse/saGammaRetheta_b.f90
     M src/adjoint/outputReverseFast/saGammaRetheta_fast_b.f90
     M src/turbulence/saGammaRetheta.F90
     M tests/reg_tests/README_SAGR.md
     M tests/reg_tests/dev/diag_all16_blocks_cs.py
     M tests/reg_tests/reg_sagr.py
     M tests/reg_tests/run_sagr_tests.sh
     M tests/reg_tests/test_adjoint_sagr.py
     M tests/reg_tests/test_jacVecProdBWDFast_sagr.py
     M tests/reg_tests/test_jacVecProdFWD_sagr.py
    ?? docs/task-log/2026-07-23-sagr-full-adjoint-test.md
- **2026-07-23 19:49:48** — auto-pip: reinstalled adflow (real build) into mach env after `cd /home/mdo/MDOLab_3_v2/adflow_sa_gamma_rethetha_paper_solver chmod +x .claude/hooks/announce_long_` (site-packages now matches ./adflow)
- **2026-07-23 20:30:14** — wrap-up prompt fired (change-set 735b17185e43); uncommitted task-relevant paths:
     M docs/current-task.md
     M docs/task-log/README.md
     M src/adjoint/outputForward/saGammaRetheta_d.f90
     M src/adjoint/outputReverse/saGammaRetheta_b.f90
     M src/adjoint/outputReverseFast/saGammaRetheta_fast_b.f90
     M src/turbulence/saGammaRetheta.F90
     M tests/reg_tests/README_SAGR.md
     M tests/reg_tests/dev/diag_all16_blocks_cs.py
     M tests/reg_tests/reg_sagr.py
     M tests/reg_tests/run_sagr_tests.sh
     M tests/reg_tests/test_adjoint_sagr.py
     M tests/reg_tests/test_jacVecProdBWDFast_sagr.py
     M tests/reg_tests/test_jacVecProdFWD_sagr.py
    ?? docs/task-log/2026-07-23-sagr-full-adjoint-test.md
    ?? tests/reg_tests/refs/adjoint_sagr_tut_wing.json
- **2026-07-23 20:40:23** — auto-pip: reinstalled adflow (complex-step (CS) build) into mach env after `cd /home/mdo/MDOLab_3_v2/adflow_sa_gamma_rethetha_paper_solver git add .claude/skills/build/SKILL.md` (site-packages now matches ./adflow)
- **2026-07-23 21:07:16** — auto-pip: reinstalled adflow (real build) into mach env after `cd /home/mdo/MDOLab_3_v2/adflow_sa_gamma_rethetha_paper_solver/tests/reg_tests; grep -n "complex bui` (site-packages now matches ./adflow)
