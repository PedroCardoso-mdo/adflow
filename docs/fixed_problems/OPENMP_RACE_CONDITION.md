# OpenMP Race Condition in Blockette

## Problem

Running with multiple OpenMP threads produced **different results** than
single-threaded execution.

Test run 2026-06-01 (before fix):

| Metric | 1 Thread | 4 Threads |
|--------|----------|-----------|
| resrho | 8.88E-01 | 2.44E-03 |
| CL     | 0.0328   | 0.0094    |
| CD     | 1498     | 0.46      |

Results should be identical within floating-point tolerance; the large
difference indicated corrupted shared data.

## Root Cause

In `src/NKSolver/blockette.F90` the OpenMP parallel region (lines 449-852)
wrote several module-level variables that were **not** declared THREADPRIVATE,
so all threads shared a single copy and clobbered each other:

- `singleHaloStart`, `doubleHaloStart`, `nodeStart`
- `dtl`, `sFaceI`, `sFaceJ`, `sFaceK`

In addition, the inner loop indices `ii, jj, kk` were neither THREADPRIVATE
(not allowed for loop variables) nor listed in the `private()` clause, so they
were shared as well.

## Solution

- Added the missing module variables to THREADPRIVATE declarations
  (`blockette.F90:68-69`).
- Added `ii, jj, kk` to the `private()` clause of the parallel-do
  (`blockette.F90:451`) — loop indices cannot be THREADPRIVATE, so they must be
  privatized on the directive.

After the fix, 1-thread and N-thread runs produce matching CL/CD (CL=0.8850,
CD=0.0328) and the parallel region scales (1.52x at 4 threads, 1.76x at 14
threads; up to 2.58x with hybrid 4 MPI x 3 OMP). See `docs/OPENMP_WORK_REPORT.md`
for the full scaling table.

## Remaining Open Item

The `4x2` hybrid configuration (4 MPI x 2 OMP) still produced wrong results,
suggesting an additional race in the multi-MPI + multi-OMP path. Tracked as a
next step in `docs/OPENMP_WORK_REPORT.md`, not yet investigated.

## Related Commits

| Commit | Description |
|--------|-------------|
| `dd79abae` | Fix OpenMP race condition in blockette.F90 |
| `b6dc4a87` | Add OpenMP work report documentation |

## Files Changed

- `src/NKSolver/blockette.F90` - Added THREADPRIVATE declarations for
  singleHaloStart/doubleHaloStart/nodeStart and dtl/sFaceI/sFaceJ/sFaceK; added
  ii/jj/kk to the parallel-do private clause.
