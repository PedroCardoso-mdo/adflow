# Current Task

> Overwrite this file's body whenever the active task changes (CLAUDE.md: one task per
> session). When a task finishes, write its summary as a new file in `task-log/` (see
> `task-log/README.md`), add it to that index, and replace this file's body with the next task.

## Task: Finish the 3D original-vs-smooth derivative table (adjoint half)

**Started:** 2026-08-25
**Status:** in progress — waiting on Tapenade rerun; CS sweeps running

See `task-log/2026-08-25-3d-original-vs-smooth-table.md` ("Remaining") for the exact command
sequence: Tapenade → make + pip → partials gate → adjoint refs (sa/smooth/hard) → train/all →
`dev/make_side_by_side_table.py` → results in `tests/reg_tests/dev/results/` and copy to
`/home/mdo/Desktop/Run/MDO_PhD/Transition/SA_BCM/Verification_tuturial_mesh/`.

### Context
- `docs/task-log/2026-08-25-3d-original-vs-smooth-table.md`
- `tests/reg_tests/dev/diag_full_derivatives_bcm.py`, `tests/reg_tests/dev/make_side_by_side_table.py`
- `tests/reg_tests/reg_bcm.py` (P=15 kPa), `tests/reg_tests/generate_bcm_restart.py`
