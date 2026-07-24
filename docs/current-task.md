# Current Task

> Overwrite this file's body whenever the active task changes (CLAUDE.md: one task per
> session). When a task finishes, write its summary as a new file in `task-log/` (see
> `task-log/README.md`), add it to that index, and replace this file's body with the next task.

## Task: Stand up the SA-BCM derivative-test harness

**Started:** 2026-07-09
**Status:** in progress

### Objective

Get `tests/reg_tests/*_bcm.py` running end-to-end so the manuscript's caveat — "adjoint-based
sensitivity results must be revalidated" — can actually be checked, per the verification ladder
in `tests/reg_tests/README_BCM.md`.

### Context (read only what's listed — CLAUDE.md rule 8)
- `docs/architecture.md` — runtime options, code path, paper-symbol lookup.
- `tests/reg_tests/README_BCM.md` — how to run the ladder and read failures.
- `tests/reg_tests/reg_bcm.py` — fixtures and assert helpers.

### Working files (this session)
- `tests/reg_tests/reg_bcm.py` — shared config/fixtures.
- `tests/reg_tests/generate_bcm_restart.py` — one-time restart CGNS generator.
- `tests/reg_tests/test_jacVecProdFWD_bcm.py`, `test_jacVecProdBWDFast_bcm.py`,
  `test_adjoint_bcm.py` — the actual test classes.

### Checklist
- [x] Discover SA-BCM's actual shape (no new state vars — modifier inside SA's own equation).
- [x] Confirm AD is currently in sync with `sa.F90` (no Tapenade rerun needed yet).
- [x] Grid/FFD: settled on SA-GR's tutorial wing wholesale (mach=0.15, alpha=1.8), not a new
      NACA0012/NLF-0416 mesh — see `reg_bcm.py` and `docs/task-log/` (prior session).
- [x] Harness code written (`reg_bcm.py`, `dev/run_bcm_case.py`, `dev/diag_blockette_bcm.py`,
      `generate_bcm_restart.py`, `test_*_bcm.py`, `run_bcm_tests.sh`) — see
      `tests/reg_tests/README_BCM.md`. Not yet run.
- [x] Paper-vs-code verification (both variants) + Tapenade diff-variable audit — see
      `docs/task-log/2026-07-24-paper-verification-and-nk-relax.md`. **Result: faithful, no bugs
      found; the only deviations are the two already-documented, deliberate smoothings.**
- [x] `NKLSRelax` option added (off by default) for the one generic NK convergence idea worth
      trying cheaply — see same task-log entry.
- [x] Real build + mach install done (this repo's `adflow` is now what's active in the shared
      mach venv — re-`pip install` the sibling SA-GR repo before switching back to that work).
- [ ] Complex build (`make -f Makefile_CS PETSC_ARCH=complex-debug`) — **never done in this
      repo**, needed for the CS ground-truth stage.
- [ ] Populate `input_files/` (`./get-input-files.sh` or copy the two tutorial-wing stock files
      from the sibling SA-GR repo — see `README_BCM.md` step 0).
- [ ] Run `dev/run_bcm_case.py` (raw log, expect it not to converge first try) for both variants.
- [ ] Run `dev/diag_blockette_bcm.py` — safety-critical, `useBlockettes=True` is the default here.
- [ ] Run `generate_bcm_restart.py`, then `run_bcm_tests.sh blockette` → `train` → full ladder.

### Definition of done
Per CLAUDE.md: compiles, smoke test N iters no NaN, plus the harness itself runs end-to-end
(even if it reveals the manuscript's suspected bug — finding it, not fixing physics, is this
task's job; see checklist above for what's still missing before that's possible).
