# ADflow Parallelization and Runtime Execution Report

Date: 2026-03-14
Codebase analyzed: `adflow- 2.12.2 dev 128`

## 1. Scope and Intent

This report documents the parallel execution path of ADflow from Python entry to preprocessing, flow solve, and adjoint solve, with explicit routine names and where MPI/OpenMP are used.

Main goals covered:
- Domain decomposition and mesh partitioning path.
- Halo exchange timing and implementation details.
- Boundary condition ordering and treatment of boundary/halo cells.
- Residual assembly path and global residual reduction path.
- NK/ANK/SANK behavior in serial and distributed settings.
- Adjoint path and matrix-free parallel behavior.
- MPI semantics (nonblocking point-to-point, collectives, synchronization points).
- OpenMP placement and practical impact.
- Optimization-focused interpretation of likely runtime bottlenecks.

Notes:
- The solver is primarily MPI-parallel. OpenMP usage is limited.
- "SNK" in user terminology corresponds to SANK labeling in this code (second-order ANK phase), not a separate standalone solver module.

## 2. High-Level Runtime Map (Python -> Fortran)

### 2.1 Initialization and preprocessing path

From `adflow/pyADflow.py`:
1. `self.adflow.partitioning.partitionandreadgrid(False)`
2. `self.adflow.preprocessingapi.preprocessing()`
3. Overset callbacks/cutting
4. `self.adflow.preprocessingapi.preprocessingoverset(flag, famList)`
5. `self.adflow.initializeflow.updatebcdataalllevels()`
6. `self.adflow.initializeflow.initflow()`

Key references:
- `adflow/pyADflow.py:230`
- `adflow/pyADflow.py:236`
- `adflow/pyADflow.py:329`
- `adflow/pyADflow.py:333`
- `adflow/pyADflow.py:334`

### 2.2 Solve path

From `adflow/pyADflow.py`:
- Steady/time-spectral: `self.adflow.solvers.solver()`
- Unsteady: `self.adflow.solvers.solverunsteadyinit()` then stepping

Key references:
- `adflow/pyADflow.py:1284`
- `adflow/pyADflow.py:1287`

Fortran entry:
- `src/solver/solvers.F90:4` `subroutine solver`
- `src/solver/solvers.F90:892` `subroutine solveState`

### 2.3 Adjoint path

From `adflow/pyADflow.py`:
- Setup: `createpetscvars`, `setupallresidualmatricesfwd`, `setuppetscksp`
- Solve: `self.adflow.adjointapi.solveadjoint(RHS, psi, True)`

Key references:
- `adflow/pyADflow.py:4036`
- `adflow/pyADflow.py:4040`
- `adflow/pyADflow.py:4041`
- `adflow/pyADflow.py:4101`

Fortran adjoint entry points:
- `src/adjoint/adjointAPI.F90:619` `setupAllResidualMatricesfwd`
- `src/adjoint/adjointAPI.F90:661` `solveAdjoint`
- `src/adjoint/adjointAPI.F90:1007` `dRdwTMatMult` (matrix-free transpose product)

## 3. Mesh Division and Partitioning (Domain Decomposition)

Top-level partition routine:
- `src/partitioning/partitioning.F90:5` `partitionAndReadGrid(partitionOnly)`

Execution order:
1. `determineGridFileNames`
2. `readBlockSizes`
3. `determineNeighborIDs`
4. `determineSections`
5. `loadBalanceGrid`
6. `initFineGridIblank`
7. `allocCoorFineGrid`
8. `readGrid`
9. spectral transforms and coordinates (`timePeriodSpectral`, `timeRotMatricesSpectral`, `fineGridSpectralCoor`)

Load balancing and splitting routines:
- `src/partitioning/loadBalance.F90:5` `loadBalanceGrid`
- `src/partitioning/loadBalance.F90:409` `blockDistribution`
- `src/partitioning/loadBalance.F90:880` `splitBlocksLoadBalance`
- `src/partitioning/loadBalance.F90:2433` `graphPartitioning`
- `src/partitioning/loadBalance.F90:2790` `splitBlock`

Data structures that carry partition metadata:
- `src/partitioning/partitionMod.F90`:
  - `distributionBlockType`
  - `splitCGNSType`
  - `subblocksOfCGNSType`
  - global arrays `blocks`, `part`, etc.

Interpretation:
- The partitioning logic supports splitting original CGNS blocks.
- Local subdomains are mapped into `flowDoms` with local dimensions and neighbor metadata.
- Boundary subfaces and internal interfaces are explicitly remapped post-partition.

## 4. Communication Model and Halo Exchange

## 4.1 Main halo APIs

- `src/utils/haloExchange.F90:5` `whalo1(...)` first halo exchange
- `src/utils/haloExchange.F90:109` `whalo2(...)` second halo exchange
- `src/utils/haloExchange.F90:553` `whalo1to1RealGeneric(...)`
- `src/utils/haloExchange.F90:1330` `wOverset(...)`
- `src/utils/haloExchange.F90:2225` `resHalo1(...)`

Each `whalo*` call does three logical parts:
1. exchange for 1-to-1 interfaces,
2. exchange for overset interfaces,
3. orphan averaging and consistency corrections.

## 4.2 MPI primitives used in halos

Inside generic halo routines (`whalo1to1RealGeneric`, similarly for integer and overset variants):
- `mpi_isend(...)`
- `mpi_irecv(...)`
- local-copy section (same-rank donor -> halo)
- `mpi_waitany(...)` loops for receives
- `mpi_waitany(...)` loops for sends

Key reference example:
- `src/utils/haloExchange.F90:620`
- `src/utils/haloExchange.F90:644`
- `src/utils/haloExchange.F90:688`
- `src/utils/haloExchange.F90:716`

Semantics:
- Point-to-point exchange is nonblocking at posting time (`Isend/Irecv`).
- Completion is synchronized before use through `Waitany` loops.
- So communication is not globally synchronous by default, but each phase is locally synchronized before consuming new halo data.

## 4.3 Where halos happen in runtime

### During explicit unsteady RK stage:
- Early pressure-only exchange:
  - `whalo1(currentLevel, 1, 0, .true., .false., .false.)`
- Apply BCs
- Full state exchange:
  - `whalo2(currentLevel, 1, nw, .true., .true., .true.)`

Reference:
- `src/solver/solvers.F90:729`
- `src/solver/solvers.F90:740`

### During RK/DADI smoothers and blockette residual path:
- After local state update and BC application:
  - `whalo1` or `whalo2` depending on second-halo need.

Reference:
- `src/solver/smoothers.F90` `executeDADIStep`

### During NK/ANK blockette-based residual evaluation:
- `blocketteRes` explicitly calls `whalo2(...)` before residual core loops.

Reference:
- `src/NKSolver/blockette.F90` in `blocketteRes`

## 5. Boundary Conditions and Boundary Cell Treatment

Global BC driver:
- `src/solver/BCRoutines.F90:15` `applyAllBC(secondHalo)`
- block-level dispatcher:
  - `src/solver/BCRoutines.F90:57` `applyAllBC_block(secondHalo)`

BC ordering in `applyAllBC_block` is explicit and important:
1. Symmetry (`bcSymm1stHalo`, optional `bcSymm2ndHalo`)
2. Polar symmetry (`bcSymmPolar*`)
3. NS adiabatic wall (`bcNSWallAdiabatic`)
4. NS isothermal wall (`bcNSWallIsoThermal`)
5. Farfield (`bcFarfield`)
6. Subsonic outflow (`bcSubsonicOutflow`)
7. Subsonic inflow (`bcSubsonicInflow`)
8. Extrapolation and supersonic outflow (`bcExtrap`)
9. Euler wall (`bcEulerWall`)
10. Supersonic inflow (`bcSupersonicInflow`)

Boundary cells are not threaded via OpenMP in this module. They are handled by block loops and BC subface pointer ranges. BC values are applied after state update and again after halo exchange in overset cases where interpolated halo states require BC reapplication.

## 6. Residual Computation and Global Residual Assembly

Primary residual routines:
- `src/solver/residuals.F90:4` `residual_block`
- `src/solver/residuals.F90:1028` `residual`

Execution pattern:
- `residual` loops over spectral instances and local domains.
- For each domain, `setPointers` then `residual_block`.
- `residual_block` computes inviscid central flux + dissipation (scalar/matrix/upwind) + optional viscous terms + source additions.

Global residual monitor assembly:
- done in `convergenceInfo` (`src/solver/solvers.F90`)
- uses local accumulation from `sumResiduals` / `sumAllResiduals`
- then global reductions:
  - summation monitors: `mpi_allreduce(..., mpi_sum, ...)`
  - max monitors: `mpi_allreduce(..., mpi_max, ...)`

References:
- `src/solver/solvers.F90:1349`
- `src/solver/solvers.F90:1538`
- `src/solver/solvers.F90:1544`

Interpretation:
- The "general solution residual" is globally joined by `MPI_Allreduce` over monitor arrays.
- Convergence checks are globally consistent by construction.

## 7. Multigrid and Smoother Runtime Sequence

Multigrid cycle driver:
- `src/solver/multiGrid.F90:825` `executeMGCycle`

Core pattern inside one nonlinear iteration:
1. set `currentLevel = groundLevel`
2. follow cycle instructions (`cycling` array): smooth / restrict / prolong
3. before smoothing when needed:
   - `timeStep(.false.)`
   - `initRes`
   - `sourceTerms`
   - `residual`
4. smoother call:
   - `RungeKuttaSmoother` or `DADISmoother`
5. at cycle end, recompute residual for new solution

Smoothers:
- RK path: `RungeKuttaSmoother` -> repeated `executeRkStage`
- DADI path: `DADISmoother` -> `executeDADIStep`

Both apply this logic:
- update local states,
- compute pressure/viscosities,
- optional early pressure halo,
- apply BCs,
- full halo exchange.

## 8. NK / ANK / SANK (SNK) Detailed Behavior

## 8.1 Solver mode switching

Top-level switch is in `solveState`:
- coarse levels: MG smoother only
- fine level: RK startup, then ANK and/or NK based on options and residual thresholds

Reference:
- `src/solver/solvers.F90:892`

## 8.2 NK (Newton-Krylov)

Key routines:
- `NKStep` (`src/NKSolver/NKSolvers.F90:512`)
- `FormJacobianNK`
- `computeResidualNK` -> `blocketteRes`
- `setRVec` (packs global residual vector and computes global norms)
- matrix-free operator setup via `MatMFFDSetBase`

NK flow (per step):
1. possibly setup solver/preconditioner
2. compute residual vector (`blocketteRes` path)
3. optionally update Jacobian/preconditioner based on lag policy
4. run KSP linear solve (`KSPSolve`)
5. line search (none/cubic/nonmonotone)
6. update state vector and residual vector

Parallel details:
- residual vector is assembled by local loops over local cells;
- global residual split norms (flow/turb) use `mpi_allreduce` in `setRVec`;
- matrix-vector products invoke distributed matrix-free callbacks.

## 8.3 ANK (Approximate Newton-Krylov)

Key routines:
- `ANKStep` (`src/NKSolver/NKSolvers.F90:3629`)
- `FormJacobianANK`
- `computeTimeStepMat`
- `computeUnsteadyResANK`
- optional turbulence KSP: `ANKTurbSolveKSP`

ANK flow highlights:
- adaptive CFL policy controls pseudo-time/preconditioner behavior.
- can run uncoupled or coupled flow+turb variants.
- uses approximate flux mode above `ANK_secondOrdSwitchTol`.
- transitions to second-order mode below threshold.

SANK/CSANK in logs:
- not separate modules; labels represent ANK phase variants:
  - `SANK`: second-order ANK
  - `CSANK`: coupled second-order ANK
  - `CANK`: coupled ANK (not yet second-order threshold)

## 8.4 Single-rank vs multi-rank behavior

Single rank:
- same routines, no inter-rank halo traffic, only local copies.
- MPI collectives become local operations with negligible cost.

Multi-rank:
- halos and global reductions dominate communication phases.
- KSP and matrix-free products depend on distributed vectors and communication patterns.
- preconditioner and residual recomputation cadence strongly affects total communication volume.

## 9. Adjoint Parallel Execution

Main adjoint solve:
- `solveAdjoint` in `src/adjoint/adjointAPI.F90:661`

Pattern:
1. allocate/zero derivative storage
2. place arrays into PETSc vectors
3. compute initial residual norms
4. set KSP tolerances
5. solve linear system with PETSc KSP
6. update adjoint state and check convergence

Matrix-free transpose product:
- `dRdwTMatMult` -> `computeMatrixFreeProductBwdFast`
- this is the central distributed operator for adjoint KSP when matrix-free is active.

Adjoint residual/Jacobian assembly setup:
- `setupAllResidualMatricesfwd` can explicitly assemble forward-mode matrices (transpose usage for adjoint).

Adjoint communication touchpoints:
- reductions for timing/statistics in `adjointAPI`.
- forward/reverse master routines use halo exchanges and `whalo2_d` / `whalo2_b` for derivative/adjoint data movement.

References:
- `src/adjoint/masterRoutines.F90:157` `whalo2(...)`
- `src/adjoint/masterRoutines.F90:493` `whalo2_d(...)`
- `src/adjoint/masterRoutines.F90:836` `whalo2_b(...)`

## 10. MPI Semantics: Buffered vs Synchronous

Observed behavior in this code:
- Point-to-point halos are generally nonblocking (`MPI_Isend`, `MPI_Irecv`), then explicitly completed with `MPI_Waitany`.
- Collectives (`MPI_Allreduce`, `MPI_Reduce`, `MPI_Alltoall`, `MPI_Bcast`) are used for global synchronization of norms, metadata, and overset structures.
- `MPI_Barrier` is used at specific I/O or staging points, not as the primary halo mechanism.

Practical meaning:
- Data transfer is posted asynchronously, but each phase waits before data consumption, so phases are effectively synchronized at completion points.
- Whether the MPI library uses internal eager buffering is implementation-dependent; code correctness does not rely on buffering semantics.

## 11. OpenMP Usage and Threading Reality

OpenMP in this branch is concentrated in blockette residual kernels:
- `src/NKSolver/blockette.F90:64-67` `THREADPRIVATE` declarations
- `src/NKSolver/blockette.F90:354` `!$OMP parallel do ... collapse(2)`

Interpretation:
- OpenMP is used to parallelize blockette tile loops inside NK/ANK residual evaluation.
- Core solver modules (`residuals.F90`, `smoothers.F90`, `BCRoutines.F90`) are largely MPI-only in this snapshot.
- Therefore, total scaling is predominantly MPI scaling with targeted thread acceleration in blockette paths.

## 12. Where Runtime Is Likely Spent

High-cost categories likely dominating wall clock:
1. Residual/flux kernels (`residual_block`, blockette cores).
2. Halo exchange phases (`whalo2` especially for wide stencils and second halos).
3. Global reductions each nonlinear iteration (residual monitors, KSP-related norms).
4. Preconditioner/Jacobian refreshes in NK/ANK (`FormJacobian*`, AMG/KSP setup).
5. Overset communication and connectivity updates for overset cases.

Adjoint runs often cost similar or more than primal due to:
- large Krylov iterations,
- transpose/matrix-free products,
- stricter tolerances for stable sensitivities.

## 13. Optimization-Focused Recommendations (Based on Current Structure)

## 13.1 Parallel efficiency

1. Reduce unnecessary global reductions in monitoring if possible.
2. Increase work per communication event (fewer halo phases per equivalent progress).
3. Tune NK/ANK jacobian lag and CFL adaptation to reduce expensive preconditioner rebuilds.
4. For overset-heavy cases, profile `updateOversetConnectivity` cadence and communication volume.

## 13.2 OpenMP and hybrid parallelism

1. Expand OpenMP beyond blockette loops into major residual kernels if thread-safe.
2. Evaluate MPI rank/thread balance (fewer ranks, more threads) where network contention dominates.

## 13.3 Numerical and data-layout considerations

1. Consider reduced precision only after sensitivity/robustness validation (do not assume safe globally).
2. Focus first on memory locality and cache reuse in flux loops and block traversal.
3. Minimize branch-heavy logic in hottest loops only when profiling confirms branch misprediction cost.

## 13.4 Immediate profiling plan

1. Instrument per-phase timers around:
   - `whalo1/whalo2`
   - `residual` / `blocketteRes`
   - `FormJacobianNK/ANK`
   - `KSPSolve`
2. Collect MPI time breakdown (point-to-point vs collectives).
3. Compare scaling curves for:
   - pure MPI
   - MPI+OpenMP (if enabled)
4. Separate primal and adjoint scaling studies.

## 14. Compact End-to-End Flow Map (Steady / Time-Spectral)

1. Python constructor -> partition -> preprocessing -> overset preprocessing -> init flow.
2. `solvers.solver` starts MG level loop.
3. `solveState` initializes residual baseline via `computeResidualNK` and `timeStep`.
4. Nonlinear loop chooses RK/DADI, ANK, NK according to thresholds/options.
5. Each update path performs BC application + halo exchanges + residual recomputation.
6. `convergenceInfo` aggregates residual/function monitors with `MPI_Allreduce`.
7. Postprocessing/write paths execute after convergence/exit.

## 15. Compact End-to-End Flow Map (Adjoint)

1. Python `_setupAdjoint` allocates PETSc vectors/matrices and KSP.
2. Optionally assemble residual matrices (`setupAllResidualMatricesfwd`).
3. `solveAdjoint` builds residual, sets tolerances, runs KSP solve.
4. Matrix-free products call transpose product callbacks (`dRdwTMatMult`).
5. Halo/reduction logic inside master derivative routines ensures distributed consistency.
6. Adjoint convergence/failure status returned to Python for sensitivity workflow.

## 16. Final Clarifications for Your Questions

- Where are halos exchanged: in `whalo1/whalo2` calls from smoothers, unsteady RK stages, blockette residual paths, and derivative/adjoint master paths.
- Are MPI communications buffered or synchronous: nonblocking post + explicit waits in point-to-point phases; collectives are synchronization points by definition.
- How are residuals joined globally: local sums -> `MPI_Allreduce` in monitor/residual-vector routines.
- How are neighboring cells accounted for in matrix-free NK/ANK: through halo-updated state fields prior to residual/Jv evaluations; matrix-free products use those consistent distributed states.
- Where OpenMP works: primarily `blockette.F90` residual tiles in this branch.

---

Companion document for optimization planning:
`doc/parallelization_hotspot_roadmap.md`
