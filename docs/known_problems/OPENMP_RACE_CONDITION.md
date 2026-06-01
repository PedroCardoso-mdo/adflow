# OpenMP Race Condition in Blockette

## Problem

Running with multiple OpenMP threads produces **different results** than single-threaded execution.

## Evidence

Test run 2026-06-01:

| Metric | 1 Thread | 4 Threads |
|--------|----------|-----------|
| resrho | 8.88E-01 | 2.44E-03 |
| CL | 0.0328 | 0.0094 |
| CD | 1498 | 0.46 |

Results should be identical (within floating-point tolerance). This large difference indicates corrupted data.

## Likely Cause

Race condition in `src/NKSolver/blockette.F90` parallel region (lines 449-852).

Possible issues:
1. **Missing private variables** - Loop indices `ii, jj, kk` not in private clause
2. **Shared arrays being written** - Arrays not in THREADPRIVATE may be shared
3. **Reduction missing** - Accumulations without proper reduction clause

## Investigation Steps

1. Add `ii, jj, kk` to private clause explicitly
2. Review all variables used inside parallel region
3. Check if any module variables are written but not THREADPRIVATE
4. Test with `OMP_NUM_THREADS=2` to see if problem scales

## Related Files

- `src/NKSolver/blockette.F90:449` - Parallel do directive
- `src/NKSolver/blockette.F90:64-67` - THREADPRIVATE declarations
- `docs/OPENMP.md` - OpenMP analysis
