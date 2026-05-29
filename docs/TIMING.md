# Timing Instrumentation

## Overview

Timing profiling tracks 55 sections using MPI timing (`mpi_wtime()`). Statistics gathered via `MPI_Allreduce` for max/sum across ranks.

## Key File

`src/NKSolver/ankProfiling.F90`

## Implementation

- Uses `mpi_wtime()` for wall-clock timing
- Collects max and sum via `MPI_Allreduce`
- Average computed as: `avg = sum / nProc`
- Report output on rank 0 only

## Data Types

| Type | Size | Notes |
|------|------|-------|
| `intType` | 4 bytes | Handles up to 2B ranks |
| `alwaysRealType` | 8 bytes (double) | Timing values |
| MPI datatype | `MPI_DOUBLE` | For timing floats |

## >64 Rank NaN Issue Analysis

The timing code itself appears safe:
- Division by `nProc` casts to double first (lines 340-343)
- All MPI datatypes correct
- Integer overflow not expected

**Likely cause is NOT in timing code** - suspect SA-BCM numerical issues that manifest with more ranks due to:
- Different domain decomposition exposing edge cases
- Different cell counts per rank affecting local values
- Potential uninitialized memory in smaller per-rank blocks

## Sections Tracked

Major timing categories:
- ANK solver iterations
- Linear solver (PETSc)
- Preconditioner setup/apply
- Residual evaluation
- Matrix-free products
- Communication/halo exchange

## Report Output (ankProfReport)

- Runs after MPI collectives complete
- "Unaccounted" time computed via subtraction
- May show small negative values (cosmetic, not NaN)
