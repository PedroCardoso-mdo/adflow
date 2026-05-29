# NaN Issue with >64 Ranks

## Problem

SA-BCM + timing crashes with NaNs when using more than 64 MPI ranks.

## Analysis Summary

**Timing code (ankProfiling.F90):** Appears safe - no >64 rank specific issues found.

**SA-BCM code (sa.F90):** Contains potential NaN sources that may be triggered by domain decomposition changes:

### Suspect Lines in sa.F90

1. **Line 323:** `Re_vorty = sqrtVort * rho / rlv`
   - Division by laminar viscosity `rlv`
   - With more ranks, smaller blocks may have cells with `rlv=0`

2. **Line 327:** `tterm1 = (Re_theta - Re_theta_c) / (Re_theta_c * SABCM_Const1)`
   - Division by product of two terms

3. **Line 335:** `exp(stransition - k_max)`
   - Large exponential; comment says `stransition` "can be huge"
   - Log-sum-exp smoothing may have edge cases

4. **Line 345:** Division by `SABCM_fsmooth`

## Hypotheses

1. **Uninitialized `rlv`:** With more ranks, some processors may have blocks with uninitialized or zero viscosity
2. **Edge cells:** Smaller per-rank block sizes expose boundary cells with invalid values
3. **Halo exchange timing:** Race condition or incomplete exchange before SA-BCM evaluation
4. **Memory layout:** Different array indexing with finer decomposition

## Debugging Steps

1. Add NaN checks after line 323, 327, 335, 345 in sa.F90
2. Add guards: `if (rlv < 1e-30) rlv = 1e-30`
3. Run with 65 ranks and print first NaN location
4. Check halo exchange completion before turbulence evaluation

## Quick Fix to Test

Add to `saSource` before line 323:
```fortran
if (rlv < 1.0e-30_realType) then
    rlv = 1.0e-30_realType
endif
```
