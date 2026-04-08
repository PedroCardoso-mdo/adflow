# ADflow Hotspot Ranking and Optimization Roadmap

Date: 2026-03-15
Companion to: `doc/parallelization_report.md`
Deep-dive companion: `doc/solver_pass_nk_ank_sa_pc_optimization.md`
Codebase: `adflow- 2.12.2 dev 128`

## 1. Purpose

This document ranks likely runtime hotspots by routine and provides:
- estimated computational complexity,
- expected communication/computation ratio,
- scaling risk indicators,
- optimization actions prioritized for practical speedup.

Important: Rankings are model-based from code structure and communication patterns, not from measured runtime traces. Use this as a high-quality starting hypothesis before final tuning.

## 2. Rating Definitions

Communication/computation ratio:
- `Very Low`: mostly local arithmetic, minimal MPI pressure.
- `Low`: local dominated, occasional collective/halo overhead.
- `Medium`: communication visible and can affect scaling.
- `High`: communication can dominate at moderate/high core counts.
- `Very High`: strong-scaling limiter likely dominated by communication and synchronization.

Complexity notation:
- `N_cells_local`: owned cells on one rank.
- `N_faces_local`: owned and relevant halo/boundary faces on one rank.
- `N_halo`: halo points exchanged per iteration.
- `N_iter_nl`: nonlinear iterations.
- `N_iter_ksp`: Krylov iterations.
- `N_ranks`: MPI ranks.

## 3. Ranked Hotspot Table

## Rank 1
Routine(s):
- `src/NKSolver/blockette.F90` `blocketteRes`, `blocketteResCore`
- `src/solver/residuals.F90` `residual_block`

Phase:
- Per-iteration residual evaluation in RK/DADI/NK/ANK loops.

Estimated complexity:
- `O(N_cells_local * stencil_ops + N_faces_local * flux_ops)` per residual call.

Communication/computation ratio:
- `Medium` in isolation; `High` when coupled tightly with frequent halo calls.

Why high impact:
- Called repeatedly in every nonlinear iteration and inside line searches/KSP matrix-free contexts.
- Includes expensive inviscid + viscous + turbulence related operations.

Scaling risks:
- Memory bandwidth pressure.
- Cache misses on large block strides.
- Insufficient threading coverage outside blockette paths.

Top optimization actions:
1. Improve memory locality and data reuse in inner loops.
2. Fuse/streamline repeated derived quantity calculations where numerically safe.
3. Expand OpenMP to additional residual kernels beyond blockette where race-free.

## Rank 2
Routine(s):
- `src/utils/haloExchange.F90` `whalo2`, `whalo1`, `whalo1to1RealGeneric`, overset exchange variants

Phase:
- State synchronization across partitions.

Estimated complexity:
- Packing/unpacking: `O(N_halo)`
- Network: `O(messages + bytes)` dependent on partition surface area and rank graph.

Communication/computation ratio:
- `High` to `Very High` at scale.

Why high impact:
- Appears in nearly all update paths (smoothers, unsteady RK stages, blockette/NK/ANK).
- Uses nonblocking point-to-point then explicit waits; completion points synchronize progress.

Scaling risks:
- Surface/volume effect under strong scaling.
- Increased MPI latency cost for many small neighbor messages.
- Overset/periodic correction overhead after receives.

Top optimization actions:
1. Reduce halo exchange frequency where mathematically acceptable.
2. Increase work per exchange (algorithmic batching).
3. Revisit partition quality to reduce edge cuts and message count.

## Rank 3
Routine(s):
- `src/solver/solvers.F90` `convergenceInfo` (`MPI_Allreduce` monitors)
- `src/NKSolver/NKSolvers.F90` `setRVec` global norm reduction

Phase:
- Global residual/monitor assembly each iteration.

Estimated complexity:
- Local accumulation: `O(N_cells_local)`
- Global sync: `O(log N_ranks)` latency-sensitive collective cost.

Communication/computation ratio:
- `High` at large rank counts.

Why high impact:
- Every nonlinear iteration uses global checks.
- Collectives enforce global synchronization and can expose load imbalance.

Scaling risks:
- Collective latency accumulation.
- Waiting for slowest rank.

Top optimization actions:
1. Reduce frequency/width of monitor reductions during early iterations.
2. Separate strict convergence reductions from optional telemetry.
3. Profile rank imbalance and rebalance domain workload.

## Rank 4
Routine(s):
- `src/NKSolver/NKSolvers.F90` `NKStep`, `ANKStep`, `FormJacobianNK`, `FormJacobianANK`, `KSPSolve` calls

Phase:
- Nonlinear Newton-like updates and linear solves.

Estimated complexity:
- Per nonlinear step: `O(N_iter_ksp * cost(MatVec + PCApply)) + O(PC_build_if_triggered)`.

Communication/computation ratio:
- `Medium` to `High` (depends on KSP iteration count and PC freshness).

Why high impact:
- Jacobian/preconditioner refreshes are expensive.
- KSP iterations can balloon when CFL/step/lag settings are suboptimal.

Scaling risks:
- Frequent PC rebuilds.
- Line-search retries increase extra residual evaluations.

Top optimization actions:
1. Tune `jacobianLag`, step factors, and CFL evolution.
2. Monitor linear residual trajectories to reduce failed/repeated solves.
3. Use separate tuning for coupled and uncoupled ANK turbulence modes.

## Rank 5
Routine(s):
- `src/partitioning/loadBalance.F90` partitioning/splitting and graph partitioning

Phase:
- Startup preprocessing.

Estimated complexity:
- One-time partition effort, roughly superlinear in block/interface graph size depending on splitter and graph partition stage.

Communication/computation ratio:
- `Low` during solve (startup only), but high strategic importance.

Why high impact:
- Partition quality directly controls halo traffic and load imbalance for entire run.

Scaling risks:
- Poor partitions cause chronic communication overhead and collective waiting.

Top optimization actions:
1. Validate partition quality metrics before long production runs.
2. Compare alternative partition settings for edge-cut minimization.
3. Favor decomposition that preserves block locality for dominant communication patterns.

## Rank 6
Routine(s):
- Overset communication/connectivity routines (`src/overset/*`, `wOverset*` in halo exchange)

Phase:
- Overset updates and fringe/orphan handling.

Estimated complexity:
- Problem- and overlap-dependent; often significant communication and metadata movement.

Communication/computation ratio:
- `High` for overset-heavy cases.

Why high impact:
- Additional communication layers beyond standard block interfaces.
- Connectivity update costs can be nontrivial in moving-grid setups.

Scaling risks:
- Communication graph complexity across overlap regions.

Top optimization actions:
1. Minimize overset update frequency when physically/algorithmically acceptable.
2. Review overlap region size and fringe population.
3. Profile overset-specific phases separately from baseline halo exchange.

## Rank 7
Routine(s):
- Turbulence update paths (`turbSolveDDADI`, turbulence residuals and BC paths)

Phase:
- RANS/turbulence equation updates.

Estimated complexity:
- `O(N_cells_local * turb_ops)` per turbulence update.

Communication/computation ratio:
- `Medium` (can rise with frequent coupled updates).

Why high impact:
- Adds substantial arithmetic and sometimes extra solve loops.
- Strongly coupled to ANK mode behavior and convergence smoothness.

Scaling risks:
- Increased nonlinear iteration count if turbulence solve is stiff.

Top optimization actions:
1. Tune coupled vs uncoupled turbulence strategy in ANK.
2. Inspect turbulence-specific linear residuals and iteration caps.
3. Verify turbulence BC/update cadence to avoid unnecessary recomputation.

## Rank 8
Routine(s):
- Adjoint matrix-free products and KSP solve path
  - `src/adjoint/adjointAPI.F90` `solveAdjoint`, `dRdwTMatMult`
  - `src/adjoint/masterRoutines.F90` halo derivative exchanges

Phase:
- Sensitivity/adjoint solve stage.

Estimated complexity:
- Similar structure to primal Krylov solve: `O(N_iter_adj_ksp * cost(MatVec_T + PCApply_T))`.

Communication/computation ratio:
- `Medium` to `High`.

Why high impact:
- Adjoint solve can match or exceed primal solve cost.
- Tight tolerances and transpose products increase iteration pressure.

Scaling risks:
- KSP stagnation or high iteration counts.
- Halo/collective synchronization within derivative paths.

Top optimization actions:
1. Track adjoint iteration histories and convergence reason codes.
2. Separate primal and adjoint tuning (they need different tolerances/PC behavior).
3. Validate matrix-free vs assembled options for target case sizes.

## 4. Cross-Cutting Bottleneck Map

Dominant bottleneck classes expected:
1. Arithmetic kernel throughput: residual and flux loops.
2. Communication latency/bandwidth: halo and collectives.
3. Synchronization: frequent global reductions and wait points.
4. Solver configuration sensitivity: NK/ANK/adjoint KSP behavior.
5. Partition quality dependency: edge cuts and load imbalance.

## 5. Optimization Roadmap (Priority Order)

## Stage A: Measurement Baseline (must do first)

1. Add per-phase timers around:
   - `blocketteRes` / `residual`
   - `whalo1` / `whalo2`
   - `FormJacobianNK` / `FormJacobianANK`
   - `KSPSolve` (primal and adjoint)
   - `convergenceInfo` reduction section
2. Capture per-rank max/min/avg timing for imbalance visibility.
3. Record iteration counters: nonlinear, KSP, line-search retries, PC rebuild count.

Deliverable:
- one CSV or markdown timing table per run.

## Stage B: Communication and Decomposition

1. Evaluate partition quality and communication volume for production meshes.
2. Compare decomposition alternatives and rank counts.
3. Reduce unnecessary monitor reductions and output frequency in early convergence phases.

Expected gain:
- Better strong-scaling efficiency and reduced idle synchronization.

## Stage C: Kernel and Memory Throughput

1. Profile cache/memory behavior in residual/blockette kernels.
2. Improve data locality and limit redundant derived quantity recomputation.
3. Expand thread-level parallel coverage where safe.

Expected gain:
- Lower time per residual call and improved node-level utilization.

## Stage D: NK/ANK/Adjoint Solver Tuning

1. Tune Jacobian lag and CFL adaptation windows.
2. Reduce failed line-search backtracks.
3. Tune KSP tolerance strategy to avoid over-solving early and under-solving late.

Expected gain:
- Fewer expensive PC rebuilds and lower total linear iterations.

## Stage E: Overset/Turbulence-Specific Refinement

1. Profile overset update phases independently.
2. Tune turbulence coupling mode and update cadence.
3. Validate improvements on both primal and adjoint workloads.

Expected gain:
- Case-dependent but potentially large for overset RANS workloads.

## 6. Initial KPI Set for Decision Making

Track these KPIs per run and per rank configuration:
1. `% time in residual kernels`
2. `% time in halo exchange`
3. `% time in collectives`
4. `avg and max KSP iterations (NK/ANK/adjoint)`
5. `PC rebuild count per 100 nonlinear iterations`
6. `parallel efficiency` vs baseline rank/thread setting
7. `time-to-solution` at fixed convergence criteria

## 7. Practical Next Experiments

Recommended first experiment matrix:
1. `MPI ranks`: low/medium/high (e.g., 32/64/128)
2. `threads per rank`: 1 vs 2-4 (if OpenMP enabled and safe)
3. `NK/ANK settings`: default vs moderate jacobian lag increase
4. `monitor frequency`: default vs reduced early-iteration reductions

For each experiment, collect:
- total runtime,
- residual and halo timing share,
- KSP iteration counts,
- convergence robustness (no deterioration in final residual target).

## 8. Usage Guidance

Use this ranking as a living artifact:
- update the ranking after real profiling data,
- promote/demote hotspots based on measured percentages,
- tie every optimization to one KPI and one acceptance criterion.

If the next step is useful, a third document can be generated with a concrete case-sheet template for your runs (machine, mesh size, rank/thread map, timers, KPIs, and conclusions) so every optimization trial is directly comparable.
