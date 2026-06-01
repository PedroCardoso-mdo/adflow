# NaN Issue with >64 Ranks

## Problem

SA-BCM branch crashes with NaN in lift/drag at iteration 0 when using more than 64 MPI ranks, even with SA-BCM and ANK profiling disabled.

## Root Cause

Uninitialized Fortran module variables in `inputPhysics`:
- `use_SABCM`
- `SABCM_Exp`
- `SABCM_Const1`, `SABCM_Const2`, `SABCM_TU`, `SABCM_S0_tanh`, `SABCM_fsmooth`, `SABCM_maxsmooth`

These variables contained garbage values until Python's `setOption` was called. With more MPI ranks, different ranks had different garbage values, causing inconsistent behavior and NaN propagation during MPI reductions for force computation.

Additionally, `Tgamma` pointer in `setPointers` used level `mm` instead of level `1`, potentially accessing unallocated memory on coarse grid levels.

## Solution

1. Added Fortran default initialization for all SABCM variables in `inputParamRoutines.F90`
2. Fixed `Tgamma` pointer to use level 1 (fine grid) in `setPointers`

## Related Commits

| Commit | Description |
|--------|-------------|
| `09dc94f7` | Fix uninitialized SABCM variables causing NaN with >64 ranks |
| `e4b0b6ae` | Add guards to prevent division by zero in SA-BCM |

## Files Changed

- `src/inputParam/inputParamRoutines.F90` - Added SABCM variable defaults
- `src/utils/utils.F90` - Fixed Tgamma pointer level
- `src/turbulence/sa.F90` - Added division guards (preventive)
- `src/NKSolver/blockette.F90` - Added division guards (preventive)
