# Blockette SA-Gamma-ReTheta Analysis and Optimization Plan

## Context

The SA-gamma-rethetha transition model is already implemented in `src/NKSolver/blockette.F90` (lines 6903-7559). However, there are critical issues:

1. **OpenMP correctness bug** - Missing THREADPRIVATE declarations cause wrong results with multiple threads
2. **Performance** - Helper functions are called rather than inlined, adding function call overhead

This plan addresses both issues while keeping the original blockette.F90 intact for SA-only cases.

---

## CRITICAL FINDING: OpenMP Issues Affect BOTH SA and SA-Gamma-ReTheta

The turbulence routines (saSource, saGammaRethetaSource, etc.) run **serially** within `blocketteResCore()`. The parallel region is at a **higher level** (block-level parallelism). Both SA-only and SA-gamma-rethetha have **identical OpenMP safety issues**.

### Issue A: sa.F90 Module Variables (HIGH RISK)

**File:** `src/turbulence/sa.F90`, line 8

Variables written **every call** without THREADPRIVATE:
- `cv13`
- `kar2Inv`
- `cw36`
- `cb3Inv`

**Written at:**
- `saSource`: Lines 1017-1020
- `saViscous`: Lines 1196-1199
- `saGammaRethetaSource`: Lines 6948-6951

**Impact:** When blocketteRes is called from parallel region, each thread overwrites these shared scalars → race condition.

### Issue B: blockette.F90 Face Velocity Arrays

**File:** `src/NKSolver/blockette.F90`, lines 53-55

Variables NOT in THREADPRIVATE:
- `sFaceI`, `sFaceJ`, `sFaceK`

**Written at:** Lines 581, 589, 597, 602-604

### Issue C: blockette.F90 Halo Start Variables

**File:** `src/NKSolver/blockette.F90`, line 18

Variables NOT in THREADPRIVATE:
- `singleHaloStart`, `doubleHaloStart`, `nodeStart`

**Written at:** Lines 376-378, 434-436

### Helper Functions: SAFE

All turbUtils.F90 helper functions are **thread-safe** (pure functions):
- `smoothMinMax` - no module state
- `flengthCorrelation` - no module state
- `rethetacCorrelation` - no module state
- `reThetaTCorrelation` - no module state

---

## Part 1: Current Implementation Analysis

### File Structure

| Routine | Lines | Status |
|---------|-------|--------|
| `saGammaRethetaSource` | 6903-7135 (233 lines) | Exists, needs OpenMP fix |
| `saGammaRethetaAdvection` | 7137-7288 (152 lines) | Exists |
| `saGammaRethetaViscous` | 7290-7537 (248 lines) | Exists |
| `saGammaRethetaResScale` | 7539-7559 (21 lines) | Exists |

### Helper Functions (Called, Not Inlined)

| Function | Called At | Source |
|----------|-----------|--------|
| `smoothMinMax` | Lines 7082, 7115, 7116 | turbUtils.F90:2371-2412 |
| `flengthCorrelation` | Line 7091 | turbUtils.F90:2329-2352 |
| `rethetacCorrelation` | Line 7086 | turbUtils.F90:2354-2369 |
| `reThetaTCorrelation` | Line 7118 | turbUtils.F90:2279-2327 |

**Note:** fv1, fv2, ft2 are already computed inline (good).

---

## Part 2: OpenMP Bug Analysis

### Symptom
- 1-thread: CL = 0.0328
- 4-threads: CL = 0.0094 (completely wrong)

### Root Cause
Module variables written inside parallel region without THREADPRIVATE declaration.

### Current THREADPRIVATE Declarations (Lines 64-67)
```fortran
!$OMP THREADPRIVATE(nx, ny, nz, il, jl, kl, ie, je, ke, ib, jb, kb)
!$OMP THREADPRIVATE(w, p, gamma, ss, x, rlv, rev, vol, aa, radI, radJ, radK)
!$OMP THREADPRIVATE(dss, volRef, d2wall, iblank, porI, porJ, porK, fw, dw)
!$OMP THREADPRIVATE(sI, sJ, sK, ux, uy, uz, vx, vy, vz, wx, wy, wz, qx, qy, qz)
```

### Missing Variables

| Variable | Declared At | Written At | Issue |
|----------|-------------|------------|-------|
| `singleHaloStart` | Line 18 | Lines 376-378, 434-436 | **MISSING THREADPRIVATE** |
| `doubleHaloStart` | Line 18 | Lines 376-378, 434-436 | **MISSING THREADPRIVATE** |
| `nodeStart` | Line 18 | Lines 376-378, 434-436 | **MISSING THREADPRIVATE** |
| `ii, jj, kk` | Line 21 | Loop counters 355-357 | **MISSING from private()** |

### Current Parallel Directive (Line 354)
```fortran
!$OMP parallel do private(i,j,k,l) collapse(2)
```

**Problem:** `ii, jj, kk` are loop indices but declared as module variables (line 21). They are NOT in the private() clause.

---

## Part 3: Implementation Plan

### Task 1: Fix sa.F90 Module Variables (CRITICAL - affects both SA and SA-GR)

**File:** `src/turbulence/sa.F90`

**Option A (Preferred):** Make constants PARAMETER instead of computed:
```fortran
! Replace computed values with compile-time constants
real(kind=realType), parameter :: cv13 = 7.1_realType**3
real(kind=realType), parameter :: kar2Inv = 1.0_realType / (0.41_realType**2)
real(kind=realType), parameter :: cw36 = 0.3_realType**6
real(kind=realType), parameter :: cb3Inv = 1.0_realType / 0.622_realType
```

**Option B:** Add THREADPRIVATE to sa.F90:
```fortran
!$OMP THREADPRIVATE(cv13, kar2Inv, cw36, cb3Inv)
```

### Task 2: Fix blockette.F90 THREADPRIVATE

**File:** `src/NKSolver/blockette.F90`

**Action 1:** Add missing THREADPRIVATE after line 67:
```fortran
!$OMP THREADPRIVATE(singleHaloStart, doubleHaloStart, nodeStart)
!$OMP THREADPRIVATE(sFaceI, sFaceJ, sFaceK)
```

**Action 2:** Modify line 354 to add loop indices to private clause:
```fortran
!$OMP parallel do private(i,j,k,l,ii,jj,kk) collapse(2)
```

**Note:** Loop indices (`ii,jj,kk`) must be in `private()` clause, NOT THREADPRIVATE (OpenMP rule: loop counters of collapsed loops must be private, not threadprivate).

### Task 2: Inline Helper Functions (Performance)

**Goal:** Eliminate function call overhead for hot-path correlations.

**Functions to inline in saGammaRethetaSource:**

#### 2a. Inline `smoothMinMax` (3 calls)
Replace function call with:
```fortran
! Instead of: result = smoothMinMax(g1, g2, p)
! Inline:
diff = g1 - g2
phi = half * (g1 + g2 + diff * tanh(p * diff))
```

#### 2b. Inline `rethetacCorrelation` (1 call)
Replace function call with:
```fortran
! Instead of: reThetaC = rethetacCorrelation(reThetaTilde)
! Inline:
if (reThetaTilde <= 1860.0_realType) then
    reThetaC = reThetaTilde - 396.035e-2_realType + ...  ! polynomial
else
    reThetaC = reThetaTilde - 593.11_realType + ...      ! different polynomial
end if
```

#### 2c. Inline `flengthCorrelation` (1 call)
Replace function call with:
```fortran
! Instead of: fLength = flengthCorrelation(reThetaTilde)
! Inline (Eqs. 49-50 from paper):
if (reThetaTilde < 400.0_realType) then
    fLength = 398.189e-1_realType - ...
else if (reThetaTilde < 596.0_realType) then
    fLength = 263.404_realType - ...
else if (reThetaTilde < 1200.0_realType) then
    fLength = 0.5_realType - ...
else
    fLength = 0.3188_realType
end if
```

#### 2d. Keep `reThetaTCorrelation` as function call
- Called only once per cell
- Complex logic with lambdaTheta dependency
- Inlining would add ~50 lines for marginal benefit

### Task 3: Verify Other Parallel Regions

Search for any other `!$OMP parallel` directives and verify all written variables are properly scoped.

---

## Part 4: Verification

### Test 1: OpenMP Correctness
```bash
# Run with 1 thread
OMP_NUM_THREADS=1 python run_case.py
# Record CL value

# Run with 4 threads
OMP_NUM_THREADS=4 python run_case.py
# CL should match 1-thread result (within tolerance)
```

### Test 2: Performance
```bash
# Before changes - measure time
time OMP_NUM_THREADS=4 python run_case.py

# After inlining - measure time
time OMP_NUM_THREADS=4 python run_case.py
# Should see improvement (estimate 5-10% in turbulence kernel)
```

### Test 3: Compilation
```bash
cd build && make -j
# Must compile without errors
```

---

## Part 5: Summary of Changes

| File | Change | Priority |
|------|--------|----------|
| `src/turbulence/sa.F90` line 8 | Make cv13, kar2Inv, cw36, cb3Inv PARAMETER or THREADPRIVATE | **Critical** |
| `src/NKSolver/blockette.F90` line 68 | Add `!$OMP THREADPRIVATE(singleHaloStart, doubleHaloStart, nodeStart)` | **Critical** |
| `src/NKSolver/blockette.F90` line 68 | Add `!$OMP THREADPRIVATE(sFaceI, sFaceJ, sFaceK)` | **Critical** |
| `src/NKSolver/blockette.F90` line 354 | Change `private(i,j,k,l)` to `private(i,j,k,l,ii,jj,kk)` | **Critical** |
| `src/NKSolver/blockette.F90` lines 7082, 7115, 7116 | Inline `smoothMinMax` | Performance |
| `src/NKSolver/blockette.F90` line 7086 | Inline `rethetacCorrelation` | Performance |
| `src/NKSolver/blockette.F90` line 7091 | Inline `flengthCorrelation` | Performance |

---

## Appendix A: SA vs SA-Gamma-ReTheta Comparison

| Aspect | SA-only | SA-Gamma-ReTheta |
|--------|---------|------------------|
| Equations | 1 (nu_tilde) | 3 (nu_tilde, gamma, Re_theta) |
| sa.F90 bug | **YES** - writes cv13 etc. | **YES** - same bug |
| Helper functions | None | 4 (all thread-safe) |
| OpenMP risk | **Same** | **Same** |

**Conclusion:** Both paths have identical OpenMP bugs. Fixing sa.F90 fixes both.

---

## Appendix B: Why Not a Separate File?

Considered creating `blockette_sa_gamma_rethetha.F90` but rejected because:
1. SA-gamma-rethetha routines already exist in blockette.F90
2. Would duplicate ~3000 lines of shared infrastructure
3. OpenMP threadprivate arrays are shared between SA and SA-gamma-rethetha
4. Dispatch logic already handles model selection
5. Maintenance burden of two files

**Recommendation:** Fix in place, keep single file.
