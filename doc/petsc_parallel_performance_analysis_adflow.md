# PETSc Parallel Performance Analysis in ADflow

## Scope and intent
This document is a code-driven analysis of how PETSc is used in ADflow, how PETSc executes in parallel, where PETSc interacts with ADflow routines, and what performance costs are visible vs hidden inside PETSc internals.

It focuses on the ANK/NK linear-solver path and matrix assembly path in the current implementation.

## 1) Where PETSc enters ADflow

### 1.1 Primary call sites in ADflow
- Solver setup and matrix-free operator setup are in `src/NKSolver/NKSolvers.F90`.
- Jacobian/PC assembly and KSP/PC configuration helpers are in `src/adjoint/adjointUtils.F90`.
- AMG shell preconditioner and multilevel KSPs are in `src/solver/amg.F90`.

### 1.2 Main PETSc objects used
- `Vec`: distributed state/residual/search direction vectors.
- `Mat`: preconditioner matrix, timestep matrix, matrix-free `MatMFFD` operators.
- `KSP`: Krylov linear solvers for NK/ANK and multilevel smoothers.
- `PC`: ASM, shell PC, and nested KSP-based preconditioning.

## 2) ADflow-PETSc interaction graph

### 2.1 Matrix-free matvec path (most frequent in ANK/NK)
In ANK, PETSc matvec callbacks are ADflow residual evaluations plus a timestep-matrix term:
1. PETSc requests `MatMult` of the matrix-free operator.
2. ADflow callback computes residual-like action (`FormFunction_mf`).
3. ADflow applies `MatMultAdd(timeStepMat, inVec, rVec, rVec)`.
4. PETSc continues GMRES orthogonalization and convergence checks.

Key code locations:
- `MatMFFDSetFunction(...)` and `MatMFFDSetBase(...)` in `src/NKSolver/NKSolvers.F90`.
- `KSPSolve(ANK_KSP, ...)` in `src/NKSolver/NKSolvers.F90`.
- Callback `FormFunction_mf` and `FormFunction_mf_turb` in `src/NKSolver/NKSolvers.F90`.

### 2.2 Matrix-based PC/Jacobian path
ADflow assembles block entries through repeated residual perturbations (FD/AD + coloring), inserts into PETSc matrices, and finalizes with assembly:
1. `setupStateResidualMatrix(...)` computes blocks and calls `MatSetValuesBlocked(...)` repeatedly.
2. `MatAssemblyBegin/End(...)` finalizes distributed matrix structure/values.
3. `MatAXPY(...)` adds timestep contributions to preconditioner matrix where used.
4. KSP/PC are configured (`setupStandardKSP` or multigrid path).

Key code locations:
- `setupStateResidualMatrix` in `src/adjoint/adjointUtils.F90`.
- `FormJacobianANK` in `src/NKSolver/NKSolvers.F90`.
- `setupStandardKSP` and `setupStandardMultigrid` in `src/adjoint/adjointUtils.F90`.

### 2.3 AMG shell preconditioner path
When multigrid shell PC is used:
1. PETSc invokes `PCShell` apply.
2. ADflow `applyShellPC` executes multilevel operations (`MGPreCon`) and/or level KSP solves.
3. Level smoothers/solves are PETSc `KSP` objects on level matrices.

Key code locations:
- `applyShellPC`, `setupShellPC`, `MGPreCon` in `src/solver/amg.F90`.

## 3) PETSc parallel model as used by ADflow

### 3.1 Communicator and partitioning
- ADflow uses `ADFLOW_COMM_WORLD` when creating PETSc vectors/matrices/KSP.
- Vector and matrix ownership are distributed by local cell counts and block size (`nState`).
- Blocked insertion (`MatSetValuesBlocked`) maps local cell/block contributions to global rows via ADflow global indices.

### 3.2 PETSc data movement categories in this code
- Point-to-point neighbor exchange handled by ADflow (halo routines) before/around residual evaluations.
- PETSc collectives/reductions inside KSP and norms/dots.
- PETSc assembly communication during `MatAssemblyBegin/End`.
- PETSc preconditioner-local solves (ASM sub-KSPs) and their setup/overlap handling.

---

## 3b) Deep-dive: how ADflow provisions PETSc with parallel layout, connectivity, and vectors

This section traces the exact code mechanisms by which ADflow tells PETSc what it owns, where off-processor entries go, and how data moves between ADflow's block arrays and PETSc's distributed Vecs and Mats.

### 3b.1 Global index assignment (`setGlobalCellsAndNodes` in `src/preprocessing/preprocessingAPI.F90`)

Before any PETSc object is created, ADflow must assign each cell a unique global integer index that identifies its row in the distributed matrix and its position in distributed Vecs. This is done once during preprocessing by `setGlobalCellsAndNodes`.

**Step 1 — local cell count.**
Each process counts its owned interior cells across all local blocks:
```fortran
nCellsLocal(level) = 0
do nn = 1, nDom
    call setPointers(nn, level, 1)
    nCellsLocal(level) = nCellsLocal(level) + nx * ny * nz
end do
```
`nx*ny*nz = (il-1)*(jl-1)*(kl-1)` which is the number of interior (non-halo) cells in a block.

**Step 2 — MPI prefix sum for the per-process offset.**
```fortran
call mpi_allreduce(nCellsLocal, nCellsGlobal, 1, adflow_integer, mpi_sum, ...)
call mpi_gather(nCellsLocal, 1, adflow_integer, nCells, 1, adflow_integer, 0, ...)
! root computes prefix sum:
nCellOffset(1) = 0
do nn = 2, nProc
    nCellOffset(nn) = nCellOffset(nn-1) + nCells(nn-1)
end do
call mpi_scatter(nCellOffset, 1, ..., nCellOffsetLocal(level), 1, ...)
```
After this, each process has `nCellOffsetLocal` = the global index of its first cell.

**Step 3 — per-block sub-offset.**
Within a process with multiple blocks, a further offset is accumulated:
```fortran
nCellBlockOffset(1) = nCellOffsetLocal(level)
do nn = 2, nDom
    nCellBlockOffset(nn) = nCellBlockOffset(nn-1) + nx_prev * ny_prev * nz_prev
end do
```

**Step 4 — per-cell global index.**
Each interior cell at `(i,j,k)` in block `nn`, spectral instance `sps`, gets:
```fortran
globalCell(i,j,k) = nCellBlockOffset(nn) * nTimeIntervalsSpectral &
                    + nx * ny * nz * (sps - 1)                     &
                    + (i-2) + (j-2)*nx + (k-2)*nx*ny
```
This is a 0-based integer. The spectral instances of a given block are laid out contiguously in global index space. Within a spectral instance the ordering is row-major `i -> j -> k`.

**Step 5 — propagate indices into halo cells.**
PETSc must be told about off-processor neighbor indices when inserting off-diagonal blocks. To make these available locally, ADflow runs a standard cell-halo exchange on the integer `globalCell` array:
```fortran
! Assign globalCell pointers as intCommVars to every block, then:
call wHalo1to1IntGeneric(1, level, sps, commPatternCell_2nd, internalCell_2nd)
```
After this, halo cells in every block hold the `globalCell` value of the owning process. This lets any stencil-walking code on this process find the exact global PETSc row index of any neighbor, even if that neighbor belongs to a different process.

**Performance note:** `setGlobalCellsAndNodes` runs once at startup. Its MPI cost is O(nProc) integer reductions — negligible compared to solver cost. Its only lasting memory cost is the per-block `globalCell(0:ib, 0:jb, 0:kb)` integer array.

---

### 3b.2 How PETSc learns the vector size and block structure

When `setupANKsolver` (in `src/NKSolver/NKSolvers.F90`) creates the state and residual vectors, the call sequence is:
```fortran
nDimW = nState * nCellsLocal(1) * nTimeIntervalsSpectral

call VecCreate(ADFLOW_COMM_WORLD, wVec, ierr)
call VecSetSizes(wVec, nDimW, PETSC_DECIDE, ierr)   ! local=nDimW, global=auto-sum
call VecSetBlockSize(wVec, nw, ierr)                ! block size = nState per cell
call VecSetType(wVec, VECMPI, ierr)
```

- `nDimW` is the **local** size in scalar DOFs: `nState` flow variables times the number of local cells times spectral instances.
- `PETSC_DECIDE` tells PETSc to compute the global size by summing `nDimW` across all processes via an internal `MPI_Allreduce`. No explicit global size is ever passed from ADflow code.
- `VecSetBlockSize` registers `nw = nState` as the block size. This enables blocked scatter, blocked norm computation, and consistent interpretation during `MatMult`.
- The same `nDimW` local size is reused for the residual Vec `rVec` and for the matrix-free operator `MatCreateMFFD`.

The turbulence solver creates an independent parallel Vec with a separate scalar-DOF count:
```fortran
nDimWTurb = nStateTurb * nCellsLocal(1) * nTimeIntervalsSpectral
call VecSetSizes(wVecTurb, nDimWTurb, PETSC_DECIDE, ierr)
call VecSetBlockSize(wVecTurb, nStateTurb, ierr)
```

**Ownership range:** Once created, PETSc internally records that this process owns rows `[iStart, iEnd)` of the global Vec, where `iStart` is determined by the prefix-sum of all `nDimW` values. This range is exactly consistent with `nCellOffsetLocal * nState`, so cell `globalCell(i,j,k)` always maps to a local row from this process's perspective.

---

### 3b.3 How PETSc gets the sparsity pattern (`statePreAllocation` in `src/adjoint/adjointUtils.F90`)

PETSc's `MatCreateBAIJ` requires the number of nonzero blocks per row split into diagonal and off-diagonal counts. ADflow computes this by walking the finite-difference stencil on every cell using `statePreAllocation`.

**Input:** a stencil array of `(di, dj, dk)` offsets defining which cells influence cell `(i,j,k)`.

**Key range check:** The on-processor global row range is:
```fortran
iRowStart = flowDoms(1, 1, 1)%globalCell(2, 2, 2)          ! first owned cell of block 1
iRowEnd   = flowDoms(nDom, 1, nTIS)%globalCell(il, jl, kl) ! last owned cell of last block
```

For each owned cell, the routine loops over stencil neighbors, reads their `globalCell` value (which is already populated in halo cells by step 3b.1), and classifies:
```fortran
if (gc >= iRowStart .and. gc <= iRowEnd) then
    onProc(ii) = onProc(ii) + 1   ! diagonal block
else
    offProc(ii) = offProc(ii) + 1 ! off-diagonal block
end if
```

**Overset handling:** for fringe (`iblank == -1`) cells, the stencil entries are replaced by the 8 donor cells accessed via `gInd`, which holds the donor `globalCell` indices. A `unique()` call deduplicates donors that appear via multiple fringes.

**Transposed matrix** (for the adjoint): off-processor nonzero counts cannot be determined locally, so ADflow uses a PETSc Vec as a temporary distributed accumulator. For each off-proc neighbor `gc`, it calls `VecSetValue(offProcVec, gc, 1.0, ADD_VALUES)` and then after `VecAssemblyBegin/End` reads back the local portion to recover the transposed row nonzero counts.

**Matrix creation:**
```fortran
call MatCreateBAIJ(ADFLOW_COMM_WORLD, nState,   &
                   nDimW, nDimW,                &   ! local rows/cols [scalar DOF]
                   PETSC_DETERMINE, PETSC_DETERMINE, &  ! global sizes auto-computed
                   0, nnzDiagonal, 0, nnzOffDiag, &
                   matrix, ierr)
```
- `nState` = block size (state DOFs per cell, e.g. 5 for Euler or 6 for RANS SA).
- `nnzDiagonal(i)` and `nnzOffDiag(i)` are in **blocks** per block-row, not scalars.
- PETSc pre-allocates exactly the right memory for local diagonal blocks and off-processor off-diagonal blocks based on this array, avoiding costly dynamic memory during assembly.

**`MAT_ROW_ORIENTED = .FALSE.`:** Because Fortran arrays are column-major but PETSc expects row-major block layout by default, this option is set to allow passing dense block submatrices directly in Fortran's column-major order without transposing.

---

### 3b.4 How ADflow copies data into and out of PETSc Vecs

PETSc Vecs are distributed objects, but each process's local portion is a contiguous array in memory. ADflow accesses it via `VecGetArrayF90`:

```fortran
call VecGetArrayF90(wVec, wvec_pointer, ierr)
ii = 1
do nn = 1, nDom
    do sps = 1, nTimeIntervalsSpectral
        call setPointers(nn, 1, sps)
        do k = 2, kl
            do j = 2, jl
                do i = 2, il
                    do l = 1, nw
                        wvec_pointer(ii) = w(i, j, k, l)
                        ii = ii + 1
                    end do
                end do
            end do
        end do
    end do
end do
call VecRestoreArrayF90(wVec, wvec_pointer, ierr)
```

Critical observation: the loop ordering `nn → sps → k → j → i → l` is **identical** to the ordering used in `setGlobalCellsAndNodes` to assign `globalCell`. So the sequential position `ii - 1` of the first variable `l=1` of cell `(i,j,k)` in block `nn` equals exactly `(globalCell(i,j,k) - nCellOffsetLocal) * nState`. There is no explicit index remapping — the loop order is the remapping.

`VecGetArrayF90` gives a zero-copy Fortran pointer to the Vec's internal buffer (local portion). The copy from `w(i,j,k,l)` loops are the only data movement cost. `VecRestoreArrayF90` is mandatory and signals to PETSc that the pointer is no longer in use.

For the residual Vec `rVec`, `setRVec` applies a volume scaling (`dw / volRef`) during the copy, and also computes the L2 norm in the same pass via a local accumulator followed by `MPI_Allreduce`.

---

### 3b.5 How matrix values are inserted with global row/column indices

Once the state/Jacobian blocks `blk(nState, nState)` are computed for a given cell `(i,j,k)`, they are inserted into the PETSc matrix using the `globalCell` index as the block row/column key:

```fortran
irow = globalCell(i, j, k)   ! 0-based global block row (from preprocessing)
! cols(1:N) = globalCell of each stencil neighbor (read from halo-populated array)
call MatSetValuesBlocked(matrix, 1, cols, 1, irow, blk, ADD_VALUES, ierr)
```

`MatSetValuesBlocked` accepts 0-based block indices. PETSc internally:
1. Checks if `irow` is owned by this process (within its local row range).
2. Checks if each column index `cols(j)` is in the diagonal range (on-process) or off-diagonal (cross-process).
3. Diagonal entries are stored directly into the local BAIJ diagonal block storage.
4. Off-diagonal entries are queued in a staging buffer.

During `MatAssemblyBegin`, PETSc initiates non-blocking MPI sends of the off-diagonal staging buffer to the target process that owns those rows. `MatAssemblyEnd` completes the communication and finalizes local storage. This is the phase where matrix topology information (what column indices were actually used) is communicated and structure is finalized.

For the timestep matrix `timeStepMat`, only diagonal blocks are inserted (`irow = icol = globalCell(i,j,k)`), so `MatAssemblyBegin/End` for that matrix involves no off-processor communication.

---

### 3b.6 How PETSc tracks connectivity for matrix-free operations

For the matrix-free operator (`MatMFFD`), PETSc does **not** store any sparsity structure. Its only record of "connections" is the callback function pointer and the base state/residual vectors set via `MatMFFDSetBase`. The operator's parallel topology is implicitly encoded in the residual evaluation callback `FormFunction_mf`:

```fortran
! Inside KSPSolve, for each matvec PETSc calls:
call FormFunction_mf(snes, vecIn, vecOut, ctx, ierr)
    ! which does:
    call setW(vecIn)           ! Vec -> ADflow block arrays (VecGetArrayF90 loop)
    call computeResidualNK()   ! ADflow residual: all halo exchanges + kernel loops
    call setRVecANK(vecOut)    ! ADflow block arrays -> Vec (VecGetArrayF90 loop)
```

PETSc treats the entire callback as an opaque function `F(x)`. The parallel communication (halo exchange) happens inside `computeResidualNK` using ADflow's own MPI routines, completely outside PETSc's visibility. PETSc only adds the finite-difference perturbation to the vector before calling the function and subtracts the base residual afterward, which are both local Vec operations on the locally owned portion.

---

### 3b.7 AMG coarse-level parallel layout

Each AMG coarse level creates its own independent distributed matrix and vectors following the same pattern:

```fortran
! ncoarse = number of local coarse cells at this level
call MatCreateBAIJ(adflow_comm_world, bs, ncoarse*bs, ncoarse*bs, &
                   PETSC_DETERMINE, PETSC_DETERMINE, 0, nnzOn, 0, nnzOff, A(lvl), ierr)
call VecCreateMPI(adflow_comm_world, ncoarse*bs, PETSC_DETERMINE, res(lvl), ierr)
```

The coarse global indices are computed independently by AMG using accumulated coarse cell offsets, following the same prefix-sum pattern as the fine-level `setGlobalCellsAndNodes`. Coarse-to-fine restriction and fine-to-coarse prolongation are done in ADflow's `restrictVec`/`prolongVec` routines by direct local-to-local array indexing through precomputed `coarseIndices(nn,lvl)` maps, with no cross-process communication required at the interpolation stage.

---

### 3b.8 Summary: what PETSc can and cannot infer about the mesh

| Information | How PETSc gets it | Code site |
|---|---|---|
| Local DOF count | `VecSetSizes(wVec, nDimW, PETSC_DECIDE)` | `NKSolvers.F90` setup |
| Global DOF count | Internal `MPI_Allreduce` of local sizes | inside PETSc |
| Block size (state vars per cell) | `VecSetBlockSize(wVec, nState)` | `NKSolvers.F90` setup |
| Matrix diagonal/off-diag nnz | `nnzDiagonal`, `nnzOffDiag` arrays at `MatCreateBAIJ` | `statePreAllocation` in `adjointUtils.F90` |
| Which global row each local cell maps to | 0-based `globalCell(i,j,k)` computed by prefix sum | `setGlobalCellsAndNodes` in `preprocessingAPI.F90` |
| Neighbor connectivity (sparsity) | Stencil walk + `globalCell` read from halo-populated array | `statePreAllocation` |
| Off-proc neighbor global indices (halo cells) | Integer halo exchange `wHalo1to1IntGeneric` on `globalCell` | `setGlobalCellsAndNodes` |
| Physical mesh geometry, cell neighbors for residual | Not passed to PETSc; hidden inside ADflow callback | `FormFunction_mf`, halo routines |
| Overset donor connectivity | Resolved by ADflow before passing `globalCell` of donors | `statePreAllocation` via `gInd` |

PETSc has no knowledge of the mesh geometry, block topology, or overset interpolation stencils. Its only inputs are: local size, block size, per-row nonzero counts, and at assembly time the actual row/column index pairs with values. All mesh knowledge is encoded by ADflow before any PETSc call is made.

## 4) What ADflow can measure directly vs what is internal to PETSc

### 4.1 Directly observable from ADflow-side instrumentation
- Time spent in ADflow callbacks and routines around PETSc calls.
- Wall time across explicit PETSc API calls (`KSPSolve`, `MatAssembly*`, `MatMultAdd`, `Vec*`).
- Custom staged timers in ADflow halo path and residual kernels.

### 4.2 Not directly observable without PETSc internal instrumentation
- Exact breakdown of `KSPSolve` internal time into:
  - basis orthogonalization kernels,
  - reduction latency vs local flops,
  - preconditioner apply internals,
  - MPI progress behavior.
- Exact split inside `MatAssembly*` between local compaction/sorting and remote communication.
- Internal memory traffic, temporary buffer behavior, and per-kernel cache behavior in PETSc internals.
- Time spent in PETSc code paths not surfaced as separate API calls.

Conclusion: ADflow timers are accurate for API-boundary wall time attribution, but not for exact PETSc internal micro-breakdowns unless PETSc event logging is enabled.

## 5) What is actually done inside key PETSc routines in this workflow

### 5.1 `KSPSolve`
In this workflow, `KSPSolve` typically performs:
- iterative Krylov steps (GMRES or configured type),
- repeated operator applications (`MatMult` callback for matrix-free),
- residual norms and dot products (global reductions),
- preconditioner applications (`PCApply`), possibly nested KSP/ASM/shell PC.

Performance implications:
- Global reductions are synchronization points that grow costly with strong scaling.
- Matrix-free mode multiplies callback cost by Krylov iteration count.

### 5.2 `MatMFFDSetBase` and matrix-free `MatMult`
- PETSc stores base state/residual data needed for finite-difference matrix-free action.
- During `MatMult`, PETSc invokes ADflow callback to evaluate function action around the base.

Performance implications:
- Expensive because each matvec triggers residual path + communication.
- Solver convergence quality directly controls callback count.

### 5.3 `MatAssemblyBegin/End`
- Finalizes distributed sparse matrix entries inserted by many `MatSetValuesBlocked` calls.
- Performs ownership reconciliation and communication for off-process entries.

Performance implications:
- Collectively synchronized phase.
- Often communication-heavy at scale.

### 5.4 `PCASM*` and sub-KSP setup
- Builds additive Schwarz decomposition, overlap, and local sub-KSPs.
- Local preconditioners (e.g., ILU) are configured per subdomain KSP.

Performance implications:
- Setup costs can be nontrivial.
- Apply cost depends on overlap, fill level, and local iteration counts.

### 5.5 `PCShell` (`applyShellPC`)
- PETSc calls user-provided apply routine.
- In ADflow this dispatches to multilevel logic and possibly level solves.

Performance implications:
- Control is ADflow-side, but each nested PETSc solve has its own reductions/synchronization.

## 6) End-to-end parallel critical paths in ANK/NK

### 6.1 Matrix-free ANK/NK step critical path
1. Optional Jacobian/PC refresh (`setupStateResidualMatrix` + assembly).
2. Base residual setup (`FormFunction_mf` + `MatMFFDSetBase`).
3. `KSPSolve` loop:
   - repeated `MatMult` callbacks into ADflow residual path,
   - `MatMultAdd(timeStepMat, ...)`,
   - PETSc reductions and PC applications.
4. Nonlinear update checks and possible line-search/backtracking residual recomputes.

### 6.2 Why scaling degrades
- ADflow side: halo exchanges and residual recomputation cost.
- PETSc side: reduction/synchronization frequency in Krylov.
- Periodic assembly phases add collective communication spikes.

## 7) Performance model for PETSc-coupled ADflow cost

A practical step model:

$$
T_{step} \approx N_{mv} \cdot T_{mv\_callback} + T_{ksp\_internal} + T_{pc\_setup} + T_{assembly} + T_{other}
$$

with

- $N_{mv}$: number of matrix-vector applications in Krylov,
- $T_{mv\_callback}$: ADflow matrix-free callback cost (residual path + halo + local kernels + copy/writeback),
- $T_{ksp\_internal}$: PETSc-internal Krylov overhead (orthogonalization + reductions + control),
- $T_{pc\_setup}$: preconditioner setup and factor/setup costs,
- $T_{assembly}$: matrix assembly communication/closure.

Interpretation:
- Reducing Krylov iterations has multiplicative benefit by reducing callback invocations and PETSc global sync events.
- Kernel micro-optimizations help, but iteration economics and synchronization frequency usually dominate at scale.

## 8) Recommended PETSc-focused profiling strategy (complete, non-hand-wavy)

### 8.1 Keep ADflow timers (already useful)
- Preserve current ANK and whalo staged timers.
- Treat them as code-region wall-time attribution.

### 8.2 Add PETSc event-level profiling for internals
Enable PETSc logging and correlate with ADflow timers:
- `-log_view` (or equivalent in run setup) to capture PETSc event costs.
- Inspect events around KSP, MatMult, MatAssembly, PCApply, Vec operations.

This is the required step to separate:
- ADflow callback cost,
- PETSc internal solver/reduction cost,
- assembly/internal communication cost.

### 8.3 Build a closure table per run
For each configuration, report:
1. ADflow total step time
2. ADflow residual/halo/copy/prep timers
3. PETSc KSP/PC/Mat events
4. closure residual (unattributed)

This gives a defensible, complete performance accounting.

## 9) High-impact optimization levers specific to PETSc coupling

1. Reduce Krylov iterations (`N_mv`) via stronger/updated preconditioning and robust tolerances.
2. Reduce matrix-free callback cost per matvec (residual path cost and halo critical path).
3. Minimize assembly frequency and assembly overhead when refresh is not needed every step.
4. Tune ASM overlap/fill/local iterations for best setup-vs-apply tradeoff.
5. For shell/multigrid path, tune level smoother settings and outer AMG iterations to lower total KSP work.

## 10) What this analysis does and does not claim

This analysis is implementation-accurate at the ADflow/PETSc interface and identifies where costs enter and synchronize.

It does not claim exact internal PETSc micro-breakdowns without PETSc event logs. Those internals are only partially observable from ADflow-side timers alone.

That distinction is essential for correct interpretation in performance papers.


\subsection{Linear Solver}
In both the startup and Newton phases, a Krylov solver, such as the generalized minimal residual (GMRES) method, is used to approximate the solution of sparse linear systems of the form 
\begin{equation}
    A x = b
\end{equation}
\cite{nk_parallel}. In parallel implementations, the unknowns (\(x\)) and their corresponding equations are distributed across processes through domain decomposition. In ADFlow, each process handles one or more blocks of the domain. For a given process \(i\), the unknowns in the linear system can be categorized as follows:


\begin{itemize}
    \item \textbf{Internal unknowns:} Variables that only appear in equations within process \( i \).
    \item \textbf{Internal-interface unknowns:} Variables assigned to process \( i \) but coupled with variables on other processes.
    \item \textbf{External-interface unknowns:} Variables belonging to other processes but appearing in equations on process \( i \).
\end{itemize}

Considering summation-by-parts techniques, internal- and external-interface unknowns represent nodes shared between adjacent blocks \cite{nk_parallel}.

If the global linear system is grouped by subdomain, the equations for process \( i \) can be expressed as:
\[
A_i x_i + E_i y_{i,\text{ext}} = b_i,
\]
where \( x_i \) and \( b_i \) are the unknowns and right-hand sides local to process \( i \), and \( y_{i,\text{ext}} \) are the external-interface unknowns coupled with \( x_i \). Grouping the unknowns this way results in a sparse block structure for the global system, with internal-interface unknowns listed last in each subdomain. This ordering improves interprocess communication efficiency and reduces indirect addressing overhead during matrix-vector multiplication \cite{nk_parallel}.

The local equations for process \( i \) can be further partitioned as:
\[
\begin{bmatrix}
B_i & F_i \\
E_i & C_i
\end{bmatrix}
\begin{bmatrix}
u_i \\
y_i
\end{bmatrix}
=
\begin{bmatrix}
f_i \\
g_i
\end{bmatrix}
+ \sum_{j \in N_i} E_{ij} y_j,
\]
where \( u_i \) are internal variables, \( y_i \) are internal-interface variables, and \( f_i, g_i \) are the corresponding partitions of \( b_i \). The set \( N_i \) contains the neighboring subdomains of \( i \) \cite{nk_parallel}.

For Krylov-subspace solvers, efficient parallelization of inner products, matrix-vector products, and preconditioners is critical. Inner products are straightforward, computed by summing local products using an MPI reduction. Matrix-vector products involve both local computations and contributions from neighboring subdomains. Communication time for interface variables can be overlapped with local computations using nonblocking MPI communications \cite{nk_parallel}.

Preconditioning is a notably challenging aspect of Newton–Krylov solvers. Although efficient serial preconditioners, such as ILU, are well-established, their parallel counterparts often encounter difficulties, including inefficiencies arising from idle processes and significant communication overhead. An effective parallel preconditioning strategy necessitates a careful balance between scalability and performance. This balance involves optimizing computational efficiency within individual processes while simultaneously minimizing interprocess communication to ensure the overall performance of the solver remains robust \cite{nk_parallel}.

ADFlow uses a parallel preconditioner based on the additive Schwarz method \cite{ADFLOW_adjoint}. The underlying methods require inversion of local submatrices. The preconditioner uses an incomplete lower/upper factorization of the local submatrix of the modified Jacobian. The blocks are composed of the flow unknowns at each node. Notice that the factorization itself does not require interprocessor communication because it is applied to the local submatrices only.

In ADFlow the additive Schwarz preconditioning is essentially a block Jacobi iteration where the equation is solved using a direct method \cite{nk_parallel}.
