# Current Task: OpenMP Validation

## Objective
Verify that existing OpenMP parallelization in blockette.F90 is working and achieving speedup.

## Plan

| Phase | Description | Status |
|-------|-------------|--------|
| 1 | Create validation test script | ✅ Done |
| 2 | Run 1-thread baseline | ✅ Done (422.9s) |
| 3 | Run 4-thread test | ✅ Done (328.7s) |
| 4 | Analyze speedup | ✅ Done - 1.29x BUT WRONG RESULTS |
| 5 | Document results / known problems | ✅ Done |
| 6 | Fix race condition in blockette.F90 | ✅ Done (commit dd79abae) |

## STATUS: COMPLETE

The OpenMP race condition was fixed (commit `dd79abae`). 1-thread and N-thread
runs now produce matching results (CL=0.8850, CD=0.0328) with measurable
speedup: 1.52x at 4 threads, 1.76x at 14 threads, up to 2.58x hybrid (4 MPI x
3 OMP). Full details in `docs/OPENMP_WORK_REPORT.md`; root cause / fix archived
in `docs/fixed_problems/OPENMP_RACE_CONDITION.md`.

### Remaining open follow-up (not part of original plan)

- The `4x2` hybrid config (4 MPI x 2 OMP) still gives wrong results — possible
  additional race in the multi-MPI + multi-OMP path. Needs investigation.
- Consider OMP for the non-blockette path (`masterRoutines.F90`) and the
  adjoint/Tapenade interaction. See "Next Steps" in `docs/OPENMP_WORK_REPORT.md`.

## Test Setup

- **Script**: `Test_OpenMp/validate_openmp.py`
- **Live log**: `tail -f Test_OpenMp/test_logs/live.log`
- **Logs saved to**: `Test_OpenMp/test_logs/<timestamp>/`

## Binding Strategy (Hybrid MPI+OpenMP)

```bash
OMP_NUM_THREADS=N OMP_PROC_BIND=close mpirun --map-by slot:PE=N --bind-to core -np 1 python analysis.py
```

- `--map-by slot:PE=N` allocates N cores per MPI rank
- `--bind-to core` binds threads to those cores
- Prevents thread migration, preserves memory locality

## Expected Outcome

- **Speedup > 1.0x** → OpenMP working
- **Speedup ~1.0x** → OpenMP not working, investigate
- **CL/CD must match** between 1-thread and 4-thread runs

## If OpenMP Not Working

Check:
1. Loop indices `ii, jj, kk` privatization in blockette.F90:449
2. THREADPRIVATE declarations complete (blockette.F90:64-67)
3. Compiler flags: `-fopenmp` in FF90_FLAGS and LINKER_FLAGS (config.mk)
4. Runtime: `OMP_NUM_THREADS` actually seen by process

## Key Files

| File | Purpose |
|------|---------|
| `src/NKSolver/blockette.F90` | OpenMP parallel region (lines 449-852) |
| `config.mk` | Build flags with -fopenmp |
| `docs/OPENMP.md` | OpenMP analysis |
| `Test_OpenMp/validate_openmp.py` | Validation test |

## Last Updated
2026-06-08 - Plan concluded; Phase 6 fix verified, docs reconciled.
