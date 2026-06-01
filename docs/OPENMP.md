# OpenMP Parallelization Analysis

## Overview

OpenMP parallelization exists in `src/NKSolver/blockette.F90` for cache-blocked residual computations.

## Build Configuration

OpenMP is enabled in `config.mk`:
- `FF90_FLAGS`: includes `-fopenmp` (line 35)
- `LINKER_FLAGS`: includes `-fopenmp` (line 45)

Note: `FF77_FLAGS` does NOT include `-fopenmp`.

## Code Structure

### Location
`src/NKSolver/blockette.F90`

### THREADPRIVATE Variables (lines 64-67)
Module-level arrays declared thread-private for safe parallel access:
```fortran
!$OMP THREADPRIVATE(nx, ny, nz, il, jl, kl, ie, je, ke, ib, jb, kb)
!$OMP THREADPRIVATE(w, p, gamma, ss, x, rlv, rev, vol, aa, radI, radJ, radK)
!$OMP THREADPRIVATE(dss, volRef, d2wall, iblank, porI, porJ, porK, fw, dw)
!$OMP THREADPRIVATE(sI, sJ, sK, ux, uy, uz, vx, vy, vz, wx, wy, wz, qx, qy, qz)
```

### Parallel Region (lines 449-852)
```fortran
!$OMP parallel do private(i,j,k,l,tCopyStart) collapse(2) reduction(+:copyTime)
do kk = 2, bkl, BS
    do jj = 2, bjl, BS
        do ii = 2, bil, BS
            ! ... blockette operations ...
        end do
    end do
end do
!$OMP END PARALLEL DO
```

## Potential Issues

1. **Loop indices `ii, jj, kk`**: Not explicitly in private clause. Compiler should auto-privatize loop variables, but worth verifying.

2. **`collapse(2)` on triple loop**: Only outer two loops (kk, jj) are collapsed. The innermost (ii) runs sequentially within each thread.

3. **No runtime verification**: No code to confirm OpenMP is actually running with multiple threads.

## Testing OpenMP

```bash
export OMP_NUM_THREADS=4
# Run solver and compare timing with OMP_NUM_THREADS=1
```

## Next Steps

1. Add debug output to verify thread count at runtime
2. Confirm `ii, jj, kk` are properly privatized
3. Benchmark with different thread counts to demonstrate speedup
