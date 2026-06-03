# MAYBE / TODO — Coarser OpenMP for SA-γ-Reθ DADI assembly

_Status: idea, not yet implemented. Captured 2026-06-03._

## Current state (what's implemented)

`src/turbulence/saGammaRetheta.F90` parallelizes the DADI residual/Jacobian
assembly with **4 separate** `!$OMP parallel do collapse(3)` regions:

| Lines | Routine | Region |
|-------|---------|--------|
| 290–765 | `Source()` | source residual + Jacobian (1 loop) |
| 840–1033 | `Viscous()` | k-direction diffusion |
| 1045–1236 | `Viscous()` | j-direction diffusion |
| 1248–1438 | `Viscous()` | i-direction diffusion |

**Not parallelized in the DADI path:**
- Advection — `turbAdvection()` in `turbUtils.F90` (shared with kw/SST/etc.; left untouched on purpose).
- Unsteady term — `unsteadyTurbTerm()` in `turbUtils.F90`.
- `saGammaReThetaSolve()` — DD-ADI tridiagonal solve; inherently sequential, stays serial.

Call order in `saGammaReTheta_block`:
`Source → turbAdvection → unsteadyTurbTerm → Viscous → ResScale → saGammaReThetaSolve`.

## How blockette does it (the better pattern)

`src/NKSolver/blockette.F90` uses **ONE** coarse `!$OMP parallel do collapse(2)`
(line 356) over the 8³ sub-blocks. Inside that single region each thread runs the
whole chain serially on its own block:
```
saSource / saAdvection / saViscous / saResScale            (SA)
saGammaRethetaSource / ...Advection / ...Viscous / ...ResScale  (SA-GR)
```
So advection + unsteady are parallelized "for free" (they're inside the block
region), with **one fork/join** instead of four.

## The idea

Refactor the DADI assembly toward blockette's single-region style:

1. Wrap the **whole assembly phase** (Source + Advection + Viscous, and unsteady)
   in **one** `!$OMP parallel` region instead of 4 separate `parallel do`s —
   lower fork/join overhead and more work in parallel.
2. Give SA-γ-Reθ its **own parallel advection** rather than threading the shared
   `turbAdvection`/`unsteadyTurbTerm` in `turbUtils.F90`. Per project rule, do
   NOT modify the shared turbUtils routines — write an SA-GR-specific path or an
   inlined advection inside the parallel region.
3. Keep the DD-ADI solve serial (must run after the full block is assembled).

## Why it matters / caveats

- **Speedup:** advection + unsteady are currently serial → an un-parallelized
  chunk of per-iteration turbulence cost (on top of the serial flow solve for
  blk=False). Parallelizing them raises the achievable speedup.
- **Overhead:** 4 fork/joins per turbulence solve → 1.
- **Constraint:** the DADI solve needs the fully-assembled `scratch`/`qq` for the
  whole block before the tridiagonal sweeps, so only the *assembly* can be
  single-region; the solve stays separate and serial.
- **Hard rule:** never touch SA-only or shared turbUtils code — see
  CLAUDE.md rule #2. Any advection parallelization must be SA-GR-local.

## Prerequisite

First confirm the **current** OpenMP scope is race-free: compare blk=False
1×1 vs 1×12 iter-by-iter on the same build. Fix any missing `private` before
expanding the parallel region.
