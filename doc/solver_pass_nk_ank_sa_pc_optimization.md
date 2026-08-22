# ADflow Solver-Pass Deep Dive: NK/ANK, SA, PC, and Performance Optimization

Date: 2026-03-15
Related docs:
- `doc/parallelization_report.md`
- `doc/parallelization_hotspot_roadmap.md`

## 1. What This Document Is

This is a practical, detailed walkthrough of the code path that matters most for runtime speed in ADflow, focused on:
- solver pass behavior for NK and ANK,
- turbulence pass behavior (especially SA),
- preconditioner (PC) setup/apply cost,
- coding style and architecture effects on performance,
- what is already done well,
- what can be optimized next,
- how to test optimizations correctly.

This is written for users who are not performance experts yet.

## 2. Big Picture of the Expensive Path

For steady/time-spectral solve, one expensive nonlinear iteration typically looks like this:
1. choose update mode in `solveState`:
   - RK/DADI multigrid cycle, or
   - ANK step, or
   - NK step
2. update states
3. apply BCs
4. halo exchange (`whalo1`/`whalo2`)
5. compute residual again
6. reduce global residual metrics (`MPI_Allreduce`)

Main dispatch:
- `src/solver/solvers.F90:892` `solveState`

NK path highlights:
- `src/NKSolver/NKSolvers.F90:512` `NKStep`
- matrix-free residual via `computeResidualNK` -> `blocketteRes`
- Krylov solve with PETSc `KSPSolve`

ANK path highlights:
- `src/NKSolver/NKSolvers.F90:3629` `ANKStep`
- PC matrix form in `FormJacobianANK`
- timestep diagonal terms in `computeTimeStepMat`
- optional turbulence KSP update in `ANKTurbSolveKSP`

## 3. NK Pass (What Happens and Why It Costs Time)

## 3.1 Sequence inside NKStep

Key routine:
- `src/NKSolver/NKSolvers.F90:512` `NKStep`

Costly operations per NK step:
1. possibly build/update Jacobian PC (`FormJacobianNK`), controlled by lag.
2. compute matrix-free base residual (`formFunction_mf`, `MatMFFDSetBase`).
3. solve linear system (`KSPSolve`).
4. line search can trigger extra residual evaluations.
5. update vectors and monitor linear convergence.

Why expensive:
- Each KSP iteration requires matrix-vector products that depend on fresh residual evaluations.
- If line search backtracks, you pay extra residual calls.
- If PC is stale, KSP iterations rise sharply.

## 3.2 What is already good

- Uses matrix-free operators to avoid explicit huge Jacobian assembly in main operator path.
- Uses lagged PC strategy (`NK_jacobianLag`) to avoid rebuilding every nonlinear step.
- Uses Eisenstat-Walker-like tolerance logic for adaptive linear solve tolerance.

## 3.3 Risks and improvement levers

- Too-frequent PC rebuild: strong cost spikes.
- Too-infrequent PC rebuild: KSP iteration count explodes.
- Line-search settings can force too many extra residual evaluations.

Practical tuning knobs:
- `NK_jacobianLag`
- `NK_rtolInit`, EW settings
- line-search mode (`noLineSearch`, `cubicLineSearch`, `nonMonotoneLineSearch`)

## 4. ANK Pass (What Happens and Why It Costs Time)

## 4.1 Sequence inside ANKStep

Key routine:
- `src/NKSolver/NKSolvers.F90:3629` `ANKStep`

Typical ANK step cost stack:
1. decide coupled/uncoupled mode and second-order switch.
2. maybe rebuild ANK PC (`FormJacobianANK`).
3. update timestep matrix (`computeTimeStepMat`).
4. KSP solve for update direction.
5. physicality checks and potential backtracking.
6. residual recompute for next iteration.
7. turbulence update path (DADI or turbulence KSP) if uncoupled RANS.

PC form routines:
- `src/NKSolver/NKSolvers.F90:1935` `FormJacobianANK`
- `src/NKSolver/NKSolvers.F90:2041` `computeTimeStepMat`

## 4.2 SANK (SNK) meaning in this code

There is no separate standalone SNK module here.
- `SANK` labels in iteration type are ANK running in second-order mode.
- `CSANK` is coupled + second-order ANK.

This matters because optimization is in ANK settings, not in a separate solver.

## 4.3 What is already good

- Adaptive CFL and step logic to avoid instability.
- Ability to switch from approximate to second-order behavior.
- Optional coupled/uncoupled turbulence treatment.

## 4.4 Risks and improvement levers

- Frequent PC refresh in difficult flow regimes can dominate runtime.
- Backtracking loops can multiply cost per nonlinear iteration.
- Turbulence coupling mode can help or hurt depending on case stiffness.

Primary knobs to tune:
- `ANK_jacobianLag`
- `ANK_CFL*` parameters
- `ANK_secondOrdSwitchTol`
- coupled vs uncoupled turbulence mode

## 5. SA and Turbulence Path (RANS Cost Driver)

Main turbulence orchestrator:
- `src/turbulence/turbAPI.F90:4` `turbSolveDDADI`

SA kernels:
- `src/turbulence/sa.F90:88` `saSource`
- `src/turbulence/sa.F90:396` `saViscous`
- plus `turbAdvection` (`src/turbulence/turbUtils.F90:825`)

Why SA can be expensive:
- heavy gradient/stress computations in cell loops,
- additional halo exchanges of turbulence variables,
- extra subiterations (`nSubIterTurb`) in decoupled turbulence solve.

What is already good:
- model-specific routines are separated (clean enough for targeted profiling).
- turbulence solve can be decoupled to control stability.

Optimization opportunities:
1. reduce unnecessary turbulence subiterations when convergence allows.
2. ensure turbulence halo exchanges happen only when needed.
3. verify order settings (`approxSA`, turbulence order) for runtime/accuracy trade.

## 6. Preconditioner (PC) and AMG Commentary

PC-related routines:
- `src/NKSolver/NKSolvers.F90:1105` `applyPC`
- `src/NKSolver/NKSolvers.F90:1169` `applyAdjointPC`
- `src/solver/amg.F90:75` `setupAMG`
- `src/solver/amg.F90:433` `destroyAMG`

What is done well:
- PC paths are modular and PETSc-driven.
- AMG infrastructure exists and coarse-level support is integrated for ANK MG preconditioning.

Performance concerns:
- PC setup cost can be very high (especially if rebuilt frequently).
- applying PC too weakly increases KSP iterations.
- AMG level count and smoother choices can be expensive if overtuned.

Practical guidance:
1. Treat PC rebuild frequency as a first-class optimization target.
2. Track and report `PC setup time` separately from `KSP solve time`.
3. Optimize for total time, not for lowest KSP iteration count alone.

## 6A. How the Linear Solver Itself Works in Parallel

This section answers the specific question: what happens inside the linear solver on many MPI ranks?

The short version is:
1. PETSc creates distributed vectors and matrices.
2. Each MPI rank owns the unknowns for its local CFD cells.
3. GMRES (or another PETSc Krylov method) runs globally across all ranks.
4. Each Krylov iteration needs:
    - a matrix-vector product,
    - one or more global dot products / norms,
    - a preconditioner application.
5. In ADflow, the matrix-vector product is often matrix-free, so it calls back into CFD residual code.
6. The preconditioner is usually approximate and more local than the full operator.

## 6A.1 What each MPI rank owns

For NK, PETSc vectors are created in `setupNKsolver`:
- `src/NKSolver/NKSolvers.F90:84` `setupNKsolver`

Local vector size:
- `nDimW = nw * nCellsLocal(1) * nTimeIntervalsSpectral`

Important point:
- each rank stores only its local state unknowns,
- but those unknowns correspond to cells whose residual depends on neighbor-cell halo data,
- therefore the linear operator is distributed, not independent per rank.

PETSc objects created in parallel:
- `VecSetType(..., VECMPI)` for distributed vectors,
- `KSPCreate(ADFLOW_COMM_WORLD, ...)` for a communicator-wide Krylov solver,
- matrix-free operator `dRdw`,
- stored approximate PC matrix `dRdwPre`.

The adjoint does the same idea in `createPETScVars`:
- `src/adjoint/adjointAPI.F90:1130` `createPETScVars`

## 6A.2 What linear system is being solved

For NK, the nonlinear update comes from a linear system that looks conceptually like:

$$
J(w) \Delta w = -R(w)
$$

but ADflow often solves a pseudo-transient modified version. In `NKMatMult` the shell matrix adds:

$$
y \leftarrow Jx + \frac{1}{\mathrm{NK\_CFL}}x
$$

Reference:
- `src/NKSolver/NKSolvers.F90:244` `NKMatMult`

So the linear operator used by the Krylov method is not just the raw Jacobian. It is a stabilized operator.

For ANK, the PC matrix also gets a diagonal timestep contribution through `computeTimeStepMat` and `MatAXPY` in `FormJacobianANK`.

## 6A.3 What happens in one Krylov iteration

Take one GMRES iteration for NK as the model.

Step 1: PETSc asks for `y = A x`
- for NK this goes through the shell/matrix-free path:
   - `NKMatMult` or `MatMFFD` path,
   - which ultimately uses `FormFunction_mf`, `computeResidualNK`, and `setRVec`.

Important routines:
- `src/NKSolver/NKSolvers.F90:437` `FormFunction_mf`
- `src/NKSolver/NKSolvers.F90:1084` `computeResidualNK`
- `src/NKSolver/NKSolvers.F90:1262` `setRVec`

What this means physically:
- the linear algebra layer asks for a Jacobian-vector action,
- ADflow computes it by perturbation/matrix-free residual machinery,
- each rank evaluates its local residual contribution,
- but local residual evaluation requires halo-consistent neighbor state.

Step 2: PETSc performs orthogonalization and residual updates
- GMRES needs dot products and norms.
- Those are global operations across ranks.
- PETSc implements them with MPI reductions internally.

This is one of the key synchronization costs in parallel linear solves.

Step 3: PETSc applies the preconditioner
- This is not the full exact operator.
- It is usually an ASM or AMG-based approximate solve.

Step 4: PETSc decides whether to continue
- based on linear residual norms,
- again requiring globally consistent information.

So even if the CFD residual computation is local-heavy, the Krylov method itself necessarily includes repeated global synchronizations.

## 6A.4 Where the matrix-vector product gets its parallelism from

The matrix-vector product is parallel because the underlying residual evaluation is parallel.

In NK:
- `FormFunction_mf` does:
   1. `setW(wVec)`
   2. `computeResidualNK(...)`
   3. `setRVec(rVec)`

Inside `computeResidualNK`:
- `blocketteRes` is called.

Inside `blocketteRes`:
1. local state/BC preparation,
2. `whalo2(...)` exchange of required state ranges,
3. blockette residual loops,
4. source terms / wall / force accumulation if needed.

That means a single Krylov matvec is not a tiny algebraic operation here. It is a distributed CFD residual pass.

That is the most important fact for performance understanding.

## 6A.5 Where MPI enters inside the linear solve

MPI enters in two places:

1. Inside operator/preconditioner application
- halo exchange for state consistency,
- neighbor communication via `MPI_Isend`, `MPI_Irecv`, `MPI_Waitany`.

2. Inside PETSc Krylov operations
- global norms,
- dot products,
- convergence tests,
- orthogonalization steps.

The first is PDE/operator communication.
The second is linear algebra synchronization.

Both matter.

At low rank counts, operator cost often dominates.
At high rank counts, GMRES global reductions can become a major limiter.

## 6A.6 How the preconditioner works in parallel

### ASM path

Common setup routine:
- `src/adjoint/adjointUtils.F90:1374` `setupStandardKSP`

Hierarchy created by the code:
1. outer KSP object,
2. global PC of type ASM,
3. subdomain KSP/PC on each local block,
4. local factorization, usually ILU.

Important detail:
- ASM overlap is configurable.
- overlap 0 is closer to block Jacobi.
- larger overlap improves robustness but costs more communication/work.

What parallel ASM means:
- each rank mostly solves local preconditioning problems,
- with limited neighbor overlap,
- so the PC is cheaper and more local than the full global solve.

This is why the preconditioner is useful: it avoids doing a full global exact solve each Krylov iteration.

### MG / Shell-PC path

Multigrid setup:
- `src/adjoint/adjointUtils.F90:1564` `setupStandardMultigrid`
- `src/solver/amg.F90:75` `setupAMG`

This path uses a shell preconditioner:
- PETSc calls an apply function,
- ADflow AMG data structures do the coarse/fine smoothing hierarchy.

Parallel meaning:
- coarse indices and coarse matrices are distributed too,
- the coarser the level, the more communication/synchronization sensitivity grows.

## 6A.7 Why linear solver scaling can get bad

The linear solver can scale poorly for several reasons:

1. Each matvec is expensive
- because it is really a distributed residual evaluation.

2. Each GMRES iteration requires global reductions
- norms and orthogonalization are synchronization-heavy.

3. Weak preconditioner quality
- raises iteration count,
- multiplying both matvec and reduction costs.

4. Strong scaling effect
- as cells per rank shrink, halo and global sync overhead becomes a larger fraction.

5. Load imbalance
- every reduction waits for the slowest rank.

## 6A.8 What is local and what is global

Mostly local:
- packing local vectors,
- local residual loops,
- local ILU subsolves,
- much of blockette arithmetic.

Neighbor-coupled:
- halo exchange before/inside residual-consistent operator applications,
- overlapped ASM regions.

Fully global:
- GMRES norms and dot products,
- convergence checks,
- some monitor reductions.

That local/neighbor/global split is the right mental model for the linear solver.

## 6A.9 What to optimize first in the linear solver

If your goal is faster linear solve time, optimize in this order:

1. Reduce Krylov iteration count sensibly
- better PC settings,
- better lag/CFL settings,
- avoid over-stale PCs.

2. Reduce cost per matvec
- optimize residual kernels and halo behavior.

3. Reduce synchronization cost
- especially high-rank GMRES reduction overhead.

4. Improve partition quality
- fewer interface cuts, better balance.

Do not start with GPU or precision changes before you know which of these dominates.

## 6A.10 A simple mental picture

One GMRES iteration in ADflow is approximately:

1. Every rank takes its local piece of `x`.
2. Neighbor data is exchanged so the operator is valid near partition boundaries.
3. Every rank computes its local piece of `y = A x`.
4. PETSc performs global reductions for GMRES algebra.
5. Every rank applies a mostly local approximate preconditioner.
6. PETSc repeats until the linear tolerance is reached.

That is how the linear solver works in parallel in practice.

## 7. Code Style and Architecture: Performance Impact

## 7.1 Style observations

Good:
- clear physics/solver modular decomposition,
- named routines for specific algorithmic phases,
- many loops already in cache-friendlier i-j-k order (i inner).

Costs from style/architecture:
- heavy use of global module state and pointer switching (`setPointers`) increases hidden side effects and can limit compiler optimization.
- repeated derived quantity recalculations across routine boundaries can cost bandwidth.
- many runtime condition branches inside critical paths may reduce vectorization opportunities.

Recommendation:
- keep correctness architecture, but isolate hottest kernels and simplify data/branch behavior there.

## 7.2 Branching (`if`) comments

General rule for hot loops:
- Avoid deep conditional trees inside innermost loops if possible.
- Move case selection (e.g., model/discretization choices) outside loops and run specialized loops.

Current code already partially does this at routine-level with model-specific routines, but some inner-loop branch effects remain.

## 7.3 RAM access/data layout comments

Fortran arrays are column-major and first index contiguous, so i-inner loops are generally correct for locality.

Still important:
- reduce indirect addressing and pointer chasing in hottest loops,
- improve temporal locality (reuse loaded fields before eviction),
- be careful with large scratch arrays that exceed cache.

## 8. Float Type / Precision Comments

Precision is centrally controlled by:
- `src/modules/precision.F90`

Current defaults and options:
- default double precision (`real(kind=8)`),
- optional single precision via `USE_SINGLE_PRECISION`,
- optional quadruple (rarely practical for performance).

Practical recommendation:
1. Keep double precision for baseline production and optimization tuning.
2. If testing single precision, do staged validation:
   - primal convergence behavior,
   - force/drag/lift deltas,
   - adjoint consistency/sensitivity deltas.
3. Mixed precision is usually safer than full single for CFD + adjoint, but requires focused implementation work.

Do not switch to single globally without strict regression checks.

## 9. OpenMP Status and What to Do

Observed OpenMP usage is concentrated in blockette:
- `src/NKSolver/blockette.F90:64-67` threadprivate data
- `src/NKSolver/blockette.F90:354` parallel-do collapse over tiles

Interpretation:
- thread parallelism exists mainly in blockette residual kernels.
- many other expensive paths are still MPI-dominant.

Recommendations:
1. Keep OpenMP in blockette and tune thread count/affinity first.
2. Add OpenMP incrementally to other hotspot kernels only after race-analysis.
3. Validate reproducibility and convergence when changing thread model.

## 10. GPU: Should You Use It?

Short answer:
- Potentially yes long-term, but not the first optimization step for this codebase as-is.

Why:
- current design has many tightly-coupled Fortran kernels, pointer/global state, and frequent MPI synchronization.
- direct GPU offload of only small parts may underperform due to transfer and orchestration overhead.

A practical GPU path would require:
1. identify top 2-3 arithmetic kernels with low branching and high arithmetic intensity,
2. refactor data layout for contiguous batch operations,
3. keep data resident on device across multiple kernel stages,
4. minimize host-device synchronization points.

Recommendation today:
- first optimize CPU path and communication efficiency.
- then evaluate GPU for selected kernels (residual/flux pieces) with a prototype branch.

## 11. What Is Already Optimized vs What Still Needs Work

Already done well:
1. modular solver architecture with explicit NK/ANK/turb/adjoint paths.
2. matrix-free nonlinear/adjoint operator usage.
3. nonblocking MPI halo pattern with explicit completion.
4. some OpenMP acceleration in blockette.

Still open for significant speedup:
1. communication frequency and collective synchronization load.
2. PC rebuild policy and KSP iteration stability.
3. kernel memory-bandwidth efficiency in residual and turbulence loops.
4. broader thread-level parallel coverage.
5. partition quality tuning for communication surface minimization.

## 12. Detailed Test and Validation Methodology (Step-by-Step)

## 12.1 Baseline profiling setup

Collect per-run metrics for:
1. total wall time,
2. time in:
   - residual/blockette,
   - halo exchange,
   - collectives/reductions,
   - PC setup,
   - KSP solve,
3. iteration counters:
   - nonlinear iterations,
   - KSP iterations per nonlinear step,
   - line-search backtracks,
   - PC rebuild count.

## 12.2 Experiment matrix (minimum)

Run at least:
1. Rank scaling: 1, 2, 4, 8, ... target production count.
2. Thread scaling: 1, 2, 4 threads/rank (if OpenMP build).
3. NK/ANK parameter sweep:
   - jacobian lag low/medium/high,
   - CFL growth/cutback variants.
4. Turbulence mode sweep:
   - coupled vs uncoupled ANK turbulence update.

## 12.3 Accuracy regression checks

For each performance variant, compare against trusted baseline:
1. final residual norms,
2. integrated functionals (CL, CD, CM) and relative error,
3. adjoint residual behavior and gradient differences (if adjoint used).

Accept only variants that preserve required accuracy and robustness.

## 12.4 Decision protocol

Use this order:
1. reject unstable/non-convergent variants,
2. among stable variants, rank by time-to-solution,
3. if two are close, choose the one with lower communication sensitivity (better scaling potential).

## 13. Fast-Start Optimization Plan (Recommended)

Phase 1 (1-2 weeks):
1. instrument timers and gather baseline on 2-3 representative meshes.
2. tune ANK/NK lag and CFL settings for best total time.
3. tune rank/thread mapping for current hardware.

Phase 2 (2-4 weeks):
1. optimize residual + turbulence kernel memory behavior.
2. reduce monitor/collective overhead where safe.
3. improve partition settings for communication reduction.

Phase 3 (longer term):
1. broaden OpenMP coverage carefully.
2. prototype GPU for one high-intensity kernel chain.
3. evaluate mixed precision strategy with strict regression testing.

## 14. Final Practical Advice

If your immediate goal is faster runs with low development risk, do this first:
1. Profile and tune NK/ANK PC rebuild cadence and KSP behavior.
2. Tune MPI decomposition and thread placement.
3. Reduce communication/synchronization overhead before touching precision or GPU.

If your long-term goal is major acceleration, prepare for kernel refactoring (data layout and side-effect reduction) before serious GPU effort.
