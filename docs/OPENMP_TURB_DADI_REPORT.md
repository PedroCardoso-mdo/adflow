# Report — OpenMP for SA-γ-Reθ DADI + Constant Bug Fix

_Date: 2026-06-03_

## 1. Objective

Parallelize the SA-noft2-Gamma-Retheta turbulence DADI residual/Jacobian
assembly with OpenMP so idle cores are used during the turbulence solve, and
measure the resulting speedup via a 1-vs-12-thread comparison.

---

## 2. Work done

### 2.1 OpenMP on the DADI turbulence assembly
**File:** `src/turbulence/saGammaRetheta.F90` (commit `abefec77`)

- Added `!$OMP parallel do collapse(3)` to the primal loops of `Source()`
  (1 loop) and `Viscous()` (3 directional loops).
- Pragmas placed only in the `#else` (primal) branch; the
  `#ifdef TAPENADE_REVERSE` paths are untouched → differentiated math
  unchanged → **no Tapenade regeneration required**.
- Per-cell scalar temporaries are `private`; loop-invariant constants and the
  global per-block arrays (`scratch`, `qq`, `w`, metrics) stay shared with
  disjoint per-cell writes (no race).
- DADI tridiagonal solve (`saGammaReThetaSolve`) left serial (inherently
  sequential).
- Always-on via `OMP_NUM_THREADS`; no new Python option. Purely additive
  (58 insertions, 0 deletions).

### 2.2 OpenMP thread-safety in blockette
**File:** `src/NKSolver/blockette.F90` (commit `d2b837b6`)

- Added missing `THREADPRIVATE` for `singleHaloStart, doubleHaloStart,
  nodeStart, dtl, sFaceI/J/K`.
- Added `ii,jj,kk` to the `private()` clause of the block parallel loop.

---

## 3. Critical bug found and fixed

### 3.1 Symptom
The `useBlockettes=True` SA-γ-Reθ run (config 1) would not converge —
stuck past 14,000 iterations — while the known-good serial reference
(`out_serial.txt`, blockettes off, 1 thread) converges in ~1200 s.

### 3.2 Root cause
In commit `d2b837b6`, converting the `sa.F90` module constants to `parameter`s,
**wrong values were hardcoded** by assuming textbook SA constants instead of
reading this codebase's `paramTurb`:

| constant | definition | hardcoded (WRONG) | correct |
|----------|------------|-------------------|---------|
| `cv13`   | `rsaCv1**3` (7.1³) | 357.911 | 357.911 ✓ |
| `kar2Inv`| `1/rsaK²` (1/0.41²) | 5.9488399 | 5.948840 ✓~ |
| `cw36`   | `rsaCw3**6` | **0.000729** (0.3⁶) | **64.0** (2.0⁶) |
| `cb3Inv` | `1/rsaCb3`  | **1.6077** (1/0.622) | **1.5** (1/0.6667) |

`cw36` was wrong by ~88000×, `cb3Inv` by 7%. These drive the SA destruction
(`fw`) and diffusion terms → broken SA physics.

### 3.3 Why only the blockette path was hit
- The `sa` module constants are consumed by the **blockette** SA-GR routines
  (`useBlockettes=True`).
- The **DADI path** (`saGammaRetheta.F90`, used when `useBlockettes=False`)
  computes these constants *locally at runtime* → was always correct. The
  1200 s serial reference (blk=False) was never affected.
- This bug also affects the **plain SA model** (any code using the `sa`
  module), so the fix matters beyond transition.

### 3.4 Fix
**File:** `src/turbulence/sa.F90` — derive the constants from `paramTurb`,
giving values that are both exact and OpenMP-safe:
```fortran
use paramTurb, only: rsaCv1, rsaK, rsaCw3, rsaCb3
real(kind=realType), parameter :: cv13    = rsaCv1**3
real(kind=realType), parameter :: kar2Inv = one/(rsaK**2)
real(kind=realType), parameter :: cw36    = rsaCw3**6
real(kind=realType), parameter :: cb3Inv  = one/rsaCb3
```
Rebuilt and reinstalled into the venv (site-packages `.so` md5 verified to
match the repo build).

---

## 4. Test harness

**Location:** `TestOpenMp/` (run via `validate_openmp.py`)

- Case: `volumeMesh_L2.cgns`, SA-noft2-Gamma-Retheta, DADI overlay
  (`ANKUseTurbDADI=True`), alpha=0, converge to `L2=1e-4`, max 25000 cycles.
- 4 configs (1 rank): (1 thread, blk=False), (12 threads, blk=False),
  (1 thread, blk=True), (12 threads, blk=True).
- `ADFLOW_ALLOW_SAGR_BLOCKETTES=1` bypasses the pyADflow guard that normally
  force-disables blockettes for SA-GR (gated edit in `pyADflow.py`).
- Logs: per-config under `TestOpenMp/test_logs/<timestamp>/`, plus a
  `summary.txt` table (ranks, omp, requested/effective blockettes, solve
  wall time, speedup, CL/CD).

### What each comparison isolates
| Comparison | Meaning |
|---|---|
| **blk=False, 1 vs 12 threads** | **Clean turbulence-OpenMP speedup** — the only parallel region active. Equivalent setup to the serial reference. |
| blk=True, 1 vs 12 threads | Flow blockette OpenMP + turbulence OpenMP combined. |
| 1 vs 12 threads (same blk) | **Correctness gate** — CL/CD and iteration count must match. |

Caveat: total wall time includes the serial flow solve, so the blk=False
ratio is a *lower bound* (Amdahl) on the turbulence-kernel speedup — a ratio
>1 still proves the OpenMP works. An undiluted kernel measurement would need a
dedicated timer around `saGammaRetheta_block`.

---

## 5. Status

- **OpenMP code:** implemented, builds clean, importable.
- **Constant bug:** identified and fixed; rebuilt + reinstalled.
- **Validation (blk=False, omp=1):** converging **healthily** post-fix (CFL
  ramping normally, residuals dropping) — final time-to-1e-4 pending.
- **Full 4-config matrix:** to be (re)run after validation; summary table will
  report time / ranks / threads / effective blockettes / speedup / CL-CD.

## 6. Pending actions
1. Confirm blk=False/omp=1 converges to 1e-4 in a time comparable to the
   1200 s reference.
2. Commit the `sa.F90` constant fix (correctness-critical, affects SA + SA-GR).
3. Run the full 4-config matrix; report the speedup table.
4. (Optional) Add a turbulence-kernel timer for the undiluted speedup figure.

---

## 7. Files touched
| File | Change | Tapenade? |
|------|--------|-----------|
| `src/turbulence/saGammaRetheta.F90` | OpenMP pragmas (primal only) | No |
| `src/turbulence/sa.F90` | constants → derived parameters (bug fix) | No (math identical to original runtime values) |
| `src/NKSolver/blockette.F90` | THREADPRIVATE + private fixes | No |
| `adflow/pyADflow.py` | env-gated SA-GR blockette guard | n/a (Python) |
| `TestOpenMp/run.py`, `validate_openmp.py` | env-driven timing harness | n/a |
