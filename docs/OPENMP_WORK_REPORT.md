# OpenMP Work Report - 2026-06-01/02

## Summary

Successfully fixed OpenMP race condition in blockette.F90 and validated hybrid MPI+OpenMP scaling.

## Problem Found

OpenMP was coded but producing **wrong results** with multiple threads.

**Symptom:** 1-thread CL=0.0328, 4-thread CL=0.0094 (completely different)

**Root cause:** Missing THREADPRIVATE declarations for module variables written inside parallel region.

## Fix Applied

File: `src/NKSolver/blockette.F90`

```fortran
! Added line 68-69:
!$OMP THREADPRIVATE(singleHaloStart, doubleHaloStart, nodeStart)
!$OMP THREADPRIVATE(dtl, sFaceI, sFaceJ, sFaceK)

! Modified line 451 - added ii,jj,kk to private clause:
!$OMP parallel do private(i,j,k,l,tCopyStart,ii,jj,kk) collapse(2) reduction(+:copyTime)
```

**Note:** Loop indices `ii,jj,kk` cannot be THREADPRIVATE (OpenMP rule), must be in `private()` clause.

**Commit:** `dd79abae`

## Validation Results

### After Fix - Results Match

| Config | Iter | CL | CD | Time |
|--------|------|----|----|------|
| 1x1 | 1120 | 0.8850 | 0.0328 | 834s |
| 1x2 | 1120 | 0.8850 | 0.0328 | 642s |
| 1x4 | 1120 | 0.8850 | 0.0328 | 548s |
| 1x8 | 1120 | 0.8850 | 0.0328 | 503s |
| 1x14 | 1120 | 0.8850 | 0.0328 | 474s |
| 2x2 | 1417 | 0.8850 | 0.0328 | 495s |
| 2x4 | 1417 | 0.8850 | 0.0328 | 455s |
| 2x7 | 1417 | 0.8850 | 0.0328 | 415s |
| 4x1 | 1225 | 0.8850 | 0.0328 | 373s |
| 4x3 | 1225 | 0.8850 | 0.0328 | 323s |

**Speedup (1x1 baseline = 834s):**
- 1x4: 1.52x
- 1x14: 1.76x
- 4x3: 2.58x

### Known Issue

4x2 configuration failed with wrong results - potential additional race condition in multi-MPI + multi-OMP. Needs investigation.

## System Configuration

- **CPU:** AMD Ryzen 9 7950X (16 cores, SMT off)
- **1 socket, 1 NUMA node**
- **Binding command:**
```bash
OMP_NUM_THREADS=N OMP_PROC_BIND=close mpirun --map-by slot:PE=N --bind-to core -np M python -u analysis.py
```

## Code Architecture

### When useBlockettes=true (default)
- `blockette.F90` handles residuals with OpenMP parallelization
- SA routines (`saSource`, `saViscous`) are inside the OMP parallel region
- **OpenMP works here** (after fix)

### When useBlockettes=false
- `masterRoutines.F90` loops over blocks
- Calls original `sa.F90` routines (no OpenMP)
- For derivatives: calls Tapenade-generated `sa_d.f90`, `sa_b.f90`

### Tapenade and OpenMP

Tapenade **does not preserve** OpenMP directives. Options:
1. Post-process generated files after `./AD_I.sh` (fragile)
2. Add OMP to block loop in `masterRoutines.F90` (safer)
3. Use blockettes mode (already has OMP)

### Reverse Mode Challenge

Tapenade reverse mode generates accumulations like:
```fortran
wb(i,j,k,:) = wb(i,j,k,:) + ...
```

These cause race conditions with naive OpenMP. Solution: thread-local arrays with reduction at end.

## Files Created/Modified

| File | Change |
|------|--------|
| `src/NKSolver/blockette.F90` | Fixed THREADPRIVATE + private clause |
| `Test_OpenMp/validate_openmp.py` | MPI x OMP scaling test script |
| `docs/CURRENT_TASK.md` | Task tracking |
| `docs/known_problems/OPENMP_RACE_CONDITION.md` | Problem documentation |
| `docs/OPENMP.md` | OpenMP analysis |

## Test Infrastructure

- **Test script:** `Test_OpenMp/validate_openmp.py`
- **Live monitoring:** `tail -f Test_OpenMp/test_logs/live.log`
- **Logs saved to:** `Test_OpenMp/test_logs/<timestamp>/`
- **Uses unbuffered output:** `python -u` + `tee`

## Next Steps

1. Investigate 4x2 failure
2. Consider adding OMP to block loop in `masterRoutines.F90` for non-blockette mode
3. For adjoint: evaluate post-processing Tapenade output vs block-level parallelization
