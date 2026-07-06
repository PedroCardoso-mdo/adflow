# ADflow — Solvers, Developer Guide, Adjoint, Differentiation & Options (compiled reference)

> **Provenance.** Compiled verbatim from the official ADflow documentation (MDO Lab), version slug `latest`, retrieved 6 July 2026, from:
> - Solvers — https://mdolab-adflow.readthedocs-hosted.com/en/latest/solvers.html
> - Developers Guide — https://mdolab-adflow.readthedocs-hosted.com/en/latest/devguide.html
> - Adjoint overview — https://mdolab-adflow.readthedocs-hosted.com/en/latest/adjoint.html
> - Differentiation overview — https://mdolab-adflow.readthedocs-hosted.com/en/latest/autodiff.html
> - Options — https://mdolab-adflow.readthedocs-hosted.com/en/latest/options.html
>
> **Fidelity note.** Prose, option names, types, defaults, and enumerated values are reproduced as published. Two things to be aware of, flagged rather than silently "fixed":
> 1. The **Options** page truncated on fetch inside `sepSensorKsRho`. The final ~10 separation/cavitation-sensor options are listed by name (from the page's own contents) but their descriptions are marked **[description not captured — see source URL]**. Nothing was invented.
> 2. A few documented **internal inconsistencies** are preserved and flagged inline with `⚠` (e.g. `CFL` default field says `1.7` while its prose says `1.5`).

---

# PART 1 — SOLVERS

This section contains some practical information about the available solver algorithms in ADflow, along with best practices. It assumes users started experimenting with ADflow following the tutorial. The additional options mentioned here can be added to baseline runscripts to improve performance.

ADflow can switch between solver algorithms during a solution. This is controlled by the *relative convergence* metric — the ratio of the current residual norm to the initial residual norm. "2 orders of magnitude relative convergence" means the L2 norm of the current residual vector is two orders lower than the initial norm. All three solver algorithms can technically be used in a single simulation; typically users use either multigrid or approximate Newton–Krylov (ANK) for startup, combined with Newton–Krylov (NK) for the final stages.

If all three solvers are used, ADflow uses multigrid until relative convergence reaches `ANKSwitchTol`; then ANK is used to reach `NKSwitchTol`; after that the solver switches to NK. Besides the switch tolerances, users must set `useANKSolver` and `useNKSolver` to `True`.

For RANS, the startup stage is the initial 4–6 orders of convergence; the terminal stage is the rest until `L2Convergence`. Reaching steady state for RANS is equivalent to reducing the residual norm by 6–12 orders. Tighter tolerances (e.g. `'L2Convergence': 1e-12`) yield more accurate solutions and gradients at higher computational cost.

## Multigrid

Multigrid is the baseline solver. It uses either Runge–Kutta (RK) or the Diagonalized Diagonally-Dominant Alternating Direction Implicit (D3ADI) algorithm as the smoother. The smoothers only update flow variables (density, momentum components, energy); a Diagonalized ADI (DADI) method updates the turbulence model for RANS.

In full multigrid, the solution initializes on the coarsest level; smoothers run to reach `L2ConvergenceCoarse`, then move to the next finer level, repeating until the finest level, where the prescribed `MGCycle` is performed.

Performance depends on a few options. More coarse grid levels generally improves performance (practically ~3–5). The number of levels is the leading number in `MGCycle` (e.g. `'3w'` = a *w* cycle with 3 levels). The second critical parameter is the CFL number (time-step size). Larger CFL → larger steps → faster convergence, but both smoothers are *conditionally stable* (a CFL limit exists before instability). Reducing CFL is more stable but slower.

Two CFL numbers: `CFL` (finest level) and `CFLCoarse` (coarser levels). RK is fully explicit; defaults `'CFL': 1.7`, `'CFLCoarse': 1.0` are reasonable. With `'smoother': 'DADI'` (D3ADI) a higher `CFL` helps, but performance doesn't generally improve with higher CFL because the factorization breaks down — practical D3ADI CFL is ~3–5.

Multigrid startup works well with multiblock meshes (coarse levels usually available). With overset meshes coarse levels may be unavailable → use ANK for startup. If the flowfield has separated regions, multigrid will likely stall → use ANK.

## Approximate Newton–Krylov (ANK)

ANK is one of ADflow's fully implicit methods, developed for efficient overset startup; it also converges difficult, heavily separated cases. Enable with `'useANKSolver': True`. Default `'ANKSwitchTol'` is `1.0`, so ADflow starts with ANK if enabled. However, if coarse levels are available and `MGCycle` is not `'sg'` (single grid), a full multigrid startup may run first, with ANK activated at the finest level. `MGStartLevel` controls the initializing grid level (default `-1` = coarsest available). ANK can also follow multigrid: set `ANKSwitchTol` to e.g. `1e-2` to switch after 2 orders.

ANK has many tunable parameters and modes. Defaults suit typical transonic aeronautical applications. For a quick solution just set `'useANKSolver': True`. For RANS, set a relatively high turbulence sub-iteration count `nSubiterTurb` between `3` and `7`.

### ANK Solver Algorithm

ANK uses backward Euler; being fully implicit it is unconditionally stable (no CFL stability limit). CFL has important outcomes: small CFL (~1–10) → small steps, favorable stability, but convergence stagnates after a few orders; high CFL (~100–1e5) → approaches Newton's method with favorable nonlinear convergence, but may stagnate if far from the solution. ADflow starts small (`'ANKCFL0': 5.0`) and adaptively ramps CFL — **pseudo-transient continuation**.

Each nonlinear iteration determines state updates by *inexactly* solving a large linear system — reducing the linear residual by a factor of `ANKLinearSolveTol` (default `0.05`). Exact solves are avoided because many cheaper nonlinear iterations are usually better (RANS startup often needs tens–hundreds of nonlinear iterations). A higher linear tolerance may destabilize when the solution vector is too far from the actual solution. The linear systems contain a Jacobian and a time-stepping term (from backward Euler); the time-stepping term uses a global CFL applied per cell; the Jacobian holds partials of residuals w.r.t. states per cell.

The linear solver is GMRES. A matrix-based preconditioner improves GMRES; the actual linear system is never formed — a matrix-free approach is used (GMRES needs only matrix-vector products). The preconditioner is based on a *first order* Jacobian with memory-reducing approximations: a 7-cell stencil vs the 33-cell stencil of the second order scheme. The matrix-free operations also contain approximations; the resulting Jacobian is somewhere between full-second and first order — hence *approximate* Newton–Krylov.

Matrix-free ops let the approximate Jacobian be up to date every nonlinear iteration, but keeping an up-to-date preconditioner every iteration is impractical (forming/factorizing is costly), so the preconditioner is **lagged** automatically. Note: the ANK CFL is only updated on iterations where the preconditioner is updated — these are marked with a `*` leading the iteration type.

After inexactly solving, an update vector is obtained but usually not taken fully; it is relaxed by a factor between 0 and 1. To determine it, first a **physicality check**: limit density and energy change per cell to `ANKPhysicalLSTol` fraction (default `0.2` = 20%). For turbulence, only updates that reduce the turbulence variable are checked, limited to 99% via `'ANKPhysicalLSTolTurb': 0.99`. Then a **backtracking line search** finds a step giving a reduction in the *unsteady* residual norm (different from the printed steady residual — the steady/total norm can increase while unsteady decreases). It starts from the physicality step size and traces back. The solver then updates the state and repeats until convergence or `NKSwitchTol`.

### ANK Matrix-Free Operations

Matrix-free ops allow modifying the Jacobian formulation on the fly at no extra memory cost. Defaults contain approximations vs exact residual routines. Users can switch to an exact Jacobian at runtime via `ANKSecondOrdSwitchTol` — a relative convergence limit above which the approximate Jacobian is used, and below which the exact formulation is used. E.g. `'ANKSecondOrdSwitchTol': 1e-2` uses approximate for the first 2 orders, then exact. This only changes the implicit system handling, not the baseline residuals — effect is only on nonlinear convergence rate and per-iteration cost.

The approximate Jacobian is better conditioned (easier to solve) but reduces update accuracy (slower nonlinear convergence). The exact Jacobian is harder to solve but yields better nonlinear convergence — tradeoff between per-iteration cost and convergence rate.

The second order switch defaults to `1e-16` (disabled). For many cases the approximate Jacobian is faster for the first 3–4 orders; switching to second order after helps. Watch the linear residual per iteration: a second order Jacobian creates a harder linear system that might fail and destabilize. If the linear tolerance can't be met (linear residual above the 0.05 default), stay with the approximate formulation and disable the second order switch. Optimal switch point is case dependent. The solver prints `S` before `ANK` when using the second order Jacobian.

### ANK Turbulence Coupling

Turbulence models can be hard to converge, so ADflow solves turbulence separately from the flow variables (density, momentum, energy). A decoupled algorithm updates flow variables, then performs turbulence sub-iterations before moving on. Best diagnostic: print the turbulence residual norm via `'resturb'` in `monitorVariables`.

Decoupled ANK has two turbulence solvers:
- **`turbDADI`** — DADI-based; uses a large number of sub-iterations. `nSubiterTurb` sets the count after each flow update (recommend 3–7; difficult cases up to 10).
- **`turbKSP`** — an isolated ANK just for turbulence, with its own matrices/linear system. Enable by `'ANKUseTurbDADI': False`. Debug info via `'ANKTurbKSPDebug': True` (prints linear convergence, step size, GMRES iterations, etc. — not pretty but useful). Recommend 1 sub-iteration (more powerful but expensive); more via `ANKNSubiterTurb` (>1).

For smaller multiblock cases (<1M cells) use `turbDADI`. For larger overset cases (>1M cells) `turbKSP` performs better.

ANK can also couple turbulence to the flow (**coupled ANK**, iterations denoted `C` before `ANK`). Coupled considers turbulence–flow coupling and is expected to converge better in late startup; but running decoupled for the first 4–5 orders is almost always better. Enable via a target relative convergence in `ANKCoupledSwitchTol` (default `1e-16`); e.g. `'ankcoupledswitchtol': 1e-4` switches after 4 orders. In coupled mode, turbulence and flow are updated together with no turbulence sub-iterations.

Two important aspects of converging SA:
1. The turbulence residual drops a few orders in the first 1–2 orders of relative convergence, then *increases* until about 4 orders (a *hill*), then usually decreases monotonically. If coupled switch happens before the turbulence goes over the hill, the solver may stall or perform badly. Use decoupled before the hill, coupled after.
2. Turbulence residual scaling: flow variables are normalized w.r.t. free-stream reference values (density and velocity of 1 = free-stream), keeping density/momentum/energy residuals similar in magnitude. The turbulence variable can span orders even after normalization, so it isn't normalized the same way — its residual norm is usually 4–5 orders lower than the flow residuals. To avoid numerical issues in coupled ANK (and NK), the turbulence residual norm is scaled by `1e4`. Coupled solvers perform well only when the printed turbulence residual norm is ~3–5 orders below the flow residual norms. If the difference is much larger/smaller than `1e4`, the scaling is off and the coupled linear systems become difficult. Users can manually set `turbResScale`, but the recommendation is to **not** modify it and instead continue with decoupled ANK until the scaling is naturally achieved.

### Interpreting the Output (ANK)

- **`Iter Tot`** — cumulative linear iterations + residual evaluations. With ANK: total GMRES iterations (flow linear system in decoupled, or coupled system in coupled) plus residual evals for the line search per nonlinear iteration.
- **`Iter Type`** — solver type; last three chars always `ANK`. Leading chars: `*` = preconditioner updated that iteration; `C` = coupled; `S` = second order.
- **`CFL`** — CFL for this iteration; only updated on preconditioner-update iterations.
- **`Step`** — relaxation factor. `1.0` = full update; less = relaxed.
- **`Lin Res`** — relative convergence achieved by the linear solver (default target 0.05). GMRES iterations are capped for cost; if it runs out, this exceeds the default (linear solution "failed", but a partial update is still usable).

With `turbKSP` under decoupled ANK, `'ANKTurbKSPDebug': True` prints turbulence info between iteration outputs (turbulence output first, then the default output — same nonlinear iteration). It prints `LIN RES, ITER, INITRES, REASON, STEP` followed by 5 numbers:
- **`LIN RES`** — relative convergence of the turbulence linear solver only.
- **`ITER`** — iterations the linear solver took.
- **`INITRES`** — initial linear residual norm (developer use).
- **`REASON`** — termination reason: `2` = desired relative convergence reached; `-3` = ran out of iterations. See PETSc `KSPConvergedReason`.
- **`STEP`** — relaxation factor for the turbulence update only (unrelated to the flow update).

### Expected Performance (ANK)

Rules of thumb — concerns: nonlinear convergence, number of nonlinear iterations, cumulative linear iterations, step size, linear residual.

- ANK should reduce total residual by 4–5 orders in ~100 iterations for simple cases (wing-only, wing-body). Complex geometries (wing-body-tail, nacelles) may reach ~200. If it takes more than a few hundred nonlinear iterations for 4–5 orders, diagnose (see Troubleshooting).
- Cumulative linear iterations at 4–5 orders should be a few thousand. Simple cases: 2–3 thousand; complex: up to 5 thousand. Maxing out linear iterations each nonlinear iteration adds load.
- Step size: repeatedly very small steps (<0.1) badly hurt nonlinear convergence. Starting small or short transients are okay; tens of iterations with small steps usually signals a problem.
- Linear residual: target 0.05, but not always reached — okay if it sits around 0.1. Repeatedly above 0.5 indicates the linear system is too stiff, the preconditioner/linear solver is too weak, or both.

### Troubleshooting (ANK)

The solver is tuned for transonic aeronautical cases. Common failure modes and fixes below (kept updated as new issues arise).

**Very Small Step Sizes.** Usually with coupled ANK — reduce the coupled switch tolerance so the solver converges tighter before switching. If small even in default (decoupled) ANK and CFL hit its limit, reduce `ANKCFLLimit` from `1e5` (slower but more stable; don't go below a few hundred). If that fails, switch to second order right before small steps start via `ANKSecondOrdSwitchTol` (exact implicit → more accurate updates). If the problem occurs before max CFL and second order doesn't help, relax the step-size algorithms: increase `ANKUnsteadyLSTol` from `1.0` (e.g. `1.5`) to allow unsteady residual increase (may help over *hills*, but risks divergence), or increase `ANKPhysicalLSTol` from `0.2` toward (but not above) 1 for more aggressive density/energy updates (less robust; >1 risks negative density/energy). If it continues, report to developers.

**High Linear Residuals.** Not reaching 0.05 is usually okay; above 0.5 is a problem; above 0.9 the solver practically stalls. Automatic mitigation: if the linear residual exceeds `ANKLinResMax`, the solver reduces CFL until linear convergence drops below it.

**Large Number of Nonlinear Iterations.** If linear systems converge to 0.05 and steps are full, yet nonlinear convergence is slow/stalled, activate second order, coupled, or both. Record the relative convergence where it stalls and use a slightly higher value for the second order/coupled switches.

**Turbulence Residuals not Converging.** Flow variables converge but the turbulence residual stalls high (want it ~4 orders below the mean-flow residuals). Because multiple turbulence solvers exist, try switching turbulence methods, or increase turbulence sub-iterations (if not coupled). With `turbKSP`, set `'ANKTurbKSPDebug': True` and monitor. Fixes here apply to the standalone `turbKSP` solver (same default algorithm as ANK).

**Special Cases.**
- Complex full-aircraft (tail, nacelle, pylon) are hard: heavy overset use introduces inter-block couplings in the global Jacobians and large–small cell volume couplings that worsen conditioning; configs like strut-braced wings or nacelles often have early-optimization separation. Expect more nonlinear iterations and more expensive linear solves.
- Actuator regions (powered fan momentum) struggle early — use ADflow's feature to gradually ramp momentum source terms as the solver converges.
- Massive separation — use a higher turbulence sub-iteration number and start with `turbDADI` or `turbKSP`. (Wing-body configs at 90° AoA have been converged.)
- Supersonic — increase `ANKPhysicalLSTol` from 0.2 to 0.4–0.6 to help when steps are tiny (moving shocks settle faster with larger physical changes).
- Very low Mach — ideally use an incompressible code; experience exists with wind turbines/automobiles. Typical issue: many nonlinear iterations from sub-optimal convergence rate. Experiment with the second order switch and a lower linear solver tolerance.

## Newton–Krylov (NK)

NK is recommended for the final stages. Best nonlinear convergence when the initial guess is near the *basin of attraction* (well-behaved cases can drop 2–3 orders in one nonlinear iteration). Used far from the solution, NK stalls or performs badly. Enable with `'useNKSolver': True` and set `NKSwitchTol` (case dependent; RANS best-case ~`1e-4`; very difficult cases down to `1e-8`; typically `1e-5`–`1e-6`). Fewer tunables than ANK, mostly linear-solver related.

### NK Solver Algorithm

NK solves the nonlinear governing equations by Newton's method — a large linear system per iteration for the state update, solved with GMRES (Krylov subspace → Newton–*Krylov*). All state variables are coupled, using the default scaling from the Turbulence Coupling section. Equivalent to Euler's method with an infinite time step → no time step in the linear system; the implicit component is only the Jacobian. A matrix-based preconditioner (approximate Jacobian) is still used, but the main driver is the exact matrix-free residual operations → always solving for the exact Jacobian. After solving, a **cubic line search** (default) guarantees residual-norm reduction. Line search set via `NKLS` — default `'cubic'`; `'non-monotone'` relaxes the decrease criterion; `'none'` takes the full step each iteration.

### Selecting the Linear Solver Tolerance (Eisenstat–Walker)

ADflow uses the Eisenstat–Walker (EW) algorithm to avoid *over-solving* the linear system per nonlinear iteration. The Newton update is a linearization about the current state: far from the solution, nonlinear gains are limited; near it, the Newton update yields several orders in one iteration. Linear solves are expensive; tighter tolerances cost more. Tradeoff: over-solving → better convergence but expensive; under-solving → cheaper but less accurate update, slower nonlinear convergence. EW monitors linear and nonlinear convergence rates and picks the optimal linear tolerance for the next iteration.

Practical behavior: the solver starts with `NKLinearSolveTol` = `0.3`. To use a constant linear tolerance, set this option **and** disable EW via `'NKUseEW': False`. With default `'NKUseEW': True`, only the first NK iteration targets 0.3; afterward the solver monitors nonlinear convergence and sets the next tolerance. Read `Lin Res` in the output. Lower chosen tolerance → satisfactory nonlinear convergence, improve with a tighter linear solve (desired; NK gradually lowers the tolerance; expect 2–3 orders nonlinear convergence per iteration). If nonlinear convergence is unsatisfactory, EW picks a *larger* tolerance to avoid over-solving (state not close; take more low-cost iterations). Usually fine — it eventually lowers again. If it doesn't after a handful of iterations, lower `NKSwitchTol` and retry (ANK handles transients better; switching to NK later helps). Hard-coded upper limit for the linear tolerance is `0.8`: consistently solving to 0.8 means the state is far away → retry with lower `NKSwitchTol`.

### NK Linear Solver Performance

The above assumes the linear solver isn't limiting (target reached every iteration) — often untrue for difficult cases, especially large overset meshes (weak/outdated preconditioner + running out of GMRES iterations). The hard-coded upper limit is 0.8; `Lin Res` above it means the linear solver is failing. As EW picks lower tolerances, the linear solver can also fail — usually okay; monitor the change in `Iter Tot` (cumulative linear iterations). If the change between nonlinear iterations exceeds `NKSubspaceSize`, the linear solver is failing to reach tolerance.

To strengthen the linear solver (each costs memory, CPU, or both; the adjoint solver is usually the memory bottleneck — see Performance — so NK has room):
- `NKPCILUFill` — higher ILU fill (more memory, more per-iteration compute).
- `NKASMOverlap` — more overlap between parallel subdomains (more communication/memory); useful with very many processors when aggressive decomposition breaks the linear solver.
- `NKOuterPreconIts` and `NKInnerPreconIts` — more global/local preconditioner iterations (no memory cost, more compute). Handle all these with care — small changes cause large linear-solver changes; avoid going past 3–4 for these parameters.
- `NKSubspaceSize` — larger GMRES subspace (default 60; more memory and increasing compute due to orthogonalization).

By default the NK preconditioner is lagged 20 nonlinear iterations (`NKJacobianLag`). If the first iteration is good but performance degrades after a few, reduce this — but forming/factorizing preconditioners is expensive, so some lag is recommended.

The preconditioner basis matrix is an approximate Jacobian fully formed by finite-differences with an efficient coloring algorithm. If FD is inaccurate, the preconditioner may not help even with strong tuning; set `NKADPC: True` to use forward-mode AD for the basis matrix (analytical partials in the approximate Jacobian; considerably more expensive per preconditioner).

### NK Troubleshooting

Three main NK failure modes below. In most cases, simply reduce the NK switch tolerance and converge further with ANK (more tunables → more likely to fix). Failure modes tend to be coupled, complicating NK diagnosis.

**Failed Linear Solutions.** The linear solver may miss the prescribed tolerance. Because EW varies the target, diagnosis isn't straightforward. Check linear iterations within a nonlinear iteration (difference in total iteration number between current and previous). If it reached `NKSubspaceSize`, the linear solver likely failed. Default subspace 60; GMRES is used without restarts. Note: due to the post-solve line search, the reported total-iteration change can exceed 60 (ADflow counts each line-search iteration as a linear iteration — costs are similar, ~one residual evaluation).

**Very Small Step Sizes.** Common — the solver can't take a meaningful step, so state changes are tiny. Ideal fix: reduce NK switch tolerance and retry. Occurs if NK starts before transients settle, or flow/turbulence residuals aren't scaled properly — ANK handles both better (recommended). Alternatively relax line search via `'NKLS': 'non-monotone'`, or disable it with `'none'` — not advised (usually diverges or produces NaNs). Even when it works, it's usually slower than converging a few more orders with ANK and retrying NK.

**EW Algorithm Stalling.** EW may consistently pick very large linear tolerances, capping NK's potential — due to unsatisfactory nonlinear convergence (multiple possible causes). Easier to return to ANK and switch to NK later. To force a constant linear tolerance, set `NKUseEW: False` and `NKLinearSolveTol`. This may add unnecessary cost — the lack of nonlinear convergence might be from small step sizes, but the solver will keep solving linear systems to tight tolerances until the max iteration limit.

---

# PART 2 — DEVELOPERS GUIDE

This guide is intended for developers of ADflow. Higher-level descriptions of the code are found in its subsections.

## Adjoint (subsection index)
- Adjoint overview → Part 3 below
  - Derivative seed manipulation
  - Python interface for derivatives

## Differentiation (subsection index)
- Differentiation overview → Part 4 below
  - Automatic Differentiation (AD)
  - Manual Differentiation

## Extra Stuff
- Triangles of zipper mesh live only in the root processor.
- `surfaceIntegration.F90`: Takes care of forces and flow-through integrations. User-defined surfaces can only be used for flow-through integrations.
- We only assemble the full Jacobian for the preconditioner. The adjoint is matrix-free. In the future we need a matrix-free preconditioner to avoid memory limitations.
- Overset interpolation is in the `wOversetGeneric` subroutine, located in `src/utils/haloExchange.F90`.

---

# PART 3 — ADJOINT OVERVIEW

*This section is based on a session given by Gaetan Kenway and Charles Mader on 03-22-2017.*

## Derivative seed manipulation

The workhorse routines for derivative seed manipulation are the master routines in `src/adjoint/masterRoutines.F90`: `master`, `master_d`, and `master_b`.

### MASTER subroutine

Inputs/outputs of the master routine (inputs are implicit to the routine):

```
            +--------------------+
      w --->|                    |---> R (dw)
            |                    |
      x --->|                    |---> forces (or tractions)
            |       MASTER       |
   x_bc --->|                    |---> funcs
            |                    |
x_extra --->|                    |
            +--------------------+
                ^
                |
             fam_list
```

Inputs (all variables that matter for the flow solution):
- `w` — the flow states.
- `x` — the volume node coordinates.
- `x_bc` — coordinates of nodes in boundary conditions such as pressure and temperature faces.
- `x_extra` — aeroProblem variables such as alpha, Mach, reference point coordinates, and free-stream pressure, temperature, and density.
- `fam_list` — list of families to subdivide the surface integration of output quantities.

All inputs except `fam_list` are active variables (should be differentiated). Every function gets an associated familyList; the default is used unless defined specifically.

Outputs (expected results from a flow solver):
- `R` (`dw` in Fortran) — flow solution residuals in every cell.
- `forces` — nodal forces or tractions (forces/area), depending on the `forcesAsTractions` flag. Includes the zipper mesh nodes.
- `funcs` — ADflow functions requested by the user; integrated quantities like `lift`, `cd`, `mass flow`, definable on arbitrary integration planes. User-defined functions are not included here (handled on the Python side).

The master function performs all operations to turn solver inputs into outputs — a collection of functions that "reproduce" the solver residual and function calculation. The call order can differ from the analysis code; there are also preprocessing calls. This code is **not** used during the actual flow solution; it is a reference to develop the differentiated code, which is done by hand.

To do this, master must touch practically all ADflow modules (hence its extensive `use` list). Anyone implementing new modules that change the residual computation **should** also add them to master and its differentiated versions. Preprocessing routines (e.g. face areas, cell volumes) are included in master because they affect the linearization.

### MASTER_B subroutine

Reverse AD version of master. These differentiated versions are done manually for two reasons:
- Tapenade cannot handle pointers such as those in `setPointers` and communication calls (MPI and PETSc).
- More efficient code by avoiding local calls to the non-linear code to push/pop values during the reverse AD run, since most necessary information is already stored during the forward pass.

**Important:** The linearization of the low speed preconditioner is not implemented because it needs the residuals before the preconditioner update, which are not stored during the forward pass. All other residual operations are incremental (they just add to the residual), so the reverse code needs no intermediate residual values to accumulate seeds.

**Note:** `master_b` is highly optimized to call only the functions needed to get the requested derivatives.

### MASTER_STATE_B subroutine

Reverse-mode AD code (e.g. via Tapenade) can be 4×–6× slower than the original due to extra data, push/pop calls, and local forward-pass recomputation of intermediate states. The adjoint solution calls the reverse code several times (to compute dRdw vector products), so it must be as fast as possible.

This subroutine is an optimized `master_b` for dRdw products assuming **frozen spatial variables** (nodal volume coordinates `x`). It is the "fast" reverse code used during adjoint solves, calling specialized "fast" versions of the differentiated modules that also assume frozen spatial variables (fewer active variables). Python never sees the fast mode — only the Fortran adjoint solver does.

## Python interface for derivatives

The main Python methods for derivative seed manipulation are `computeJacobianVectorProductFwd` and `computeJacobianVectorProductBwd` in `pyADflow.py`. These eventually call the Fortran master routines.

### ComputeJacobianVectorProductBwd

Can receive and return multiple subsets of seeds. The X̄_v seeds are typically only for derivative verification. X̄_s seeds are available only if the mesh object is included (`ADFLOW.setMesh`). X̄_Dv seeds need both the mesh object and the DVGeo object (`ADFLOW.setDVGeo`); X̄_Dv includes X̄_DvAero (design variables from the AeroProblem).

All native cost functions (`cl`, `cd`, `lift`, `drag`, …) are applied at all design families first, but family subsets can be specified (e.g. wing/fuselage drag breakdown). One of the first procedures assembles `funcsBar`, a matrix of derivative seeds for all functions (columns) and every family subset (rows), ensuring the smallest number of reverse passes.

Users can supply a custom function of the intrinsic ADflow functions; ADflow uses complex-step to get its Jacobian, potentially reducing the number of adjoints. E.g. instead of separate adjoints for `cl` and `cd` to get L/D sensitivity, solve one adjoint for cl/cd:

```
def userFunc(funcs):
    funcs['LD'] = funcs['cl']/funcs['cd']
```

BC definitions and sensitivities are stored under `AeroProblem` to ease different BCs for multipoint cases. The reverse master (`master_b`) is wrapped by `computeMatrixFreeProductBwd` in `adjointAPI.F90`, which allocates memory etc. before derivatives are computed.

---

# PART 4 — DIFFERENTIATION OVERVIEW

## Automatic Differentiation (AD)

When running Tapenade in reverse mode, ADflow explicitly asks to use "d" to flag derivative seeds instead of the usual "b" identifier, so the same variable definition works for both forward and reverse modes. However, subroutine names still use "_b" for their reverse counterparts.

Tapenade should be called in a single command to generate all AD code — a single head ensures correct dependencies. Almost all variables are treated as input and output simultaneously so seeds are not zeroed after subroutine calls (zeroing vectors is costly). Therefore, all seeds must be set to zero at the start of the reverse AD pass — done in `src/adjoint/adjointUtils.F90` at `zeroADSeeds` (invoked in `adjointAPI.F90`).

If developers make important changes to ADflow modules, the differentiated code must be updated. Steps:
- Open a terminal at `src/adjoint`.
- Run:

```
$ make -f Makefile_tapenade ad_forward
$ make -f Makefile_tapenade ad_reverse
$ make -f Makefile_tapenade ad_reverse_fast
```

## Manual Differentiation

Reverse-code generation is not fully automatic. Some subroutines are manually differentiated for two reasons:
- Tapenade cannot handle pointers and communication calls (MPI and PETSc).
- More efficient code by avoiding local calls to the non-linear code to push/pop values during the reverse AD run (info already stored during the forward pass).

Examples of manually differentiated subroutines: `master`, `haloExchange`, `surfaceIntegration`, and `getSolution`.

If developers change existing base modules, no manual differentiation is needed — just run the AD command above. If **new modules** are added, or the residual computation drastically changes, manual differentiation will be necessary (especially in `master`, `master_d`, and `master_b`).

**Note:** If a routine is differenced by hand, it should be labeled clearly at the top as manually differentiated. These routines mostly call code that has been differentiated by Tapenade.

### Differentiation of communication calls

Tapenade does not differentiate functions with MPI and PETSc calls, so they must be manually differentiated — for every MPI call in the non-linear (original) code, add a corresponding MPI or PETSc call to transfer seeds. Reverse counterparts of some communication procedures:
- MPI Sends become MPI Receives, and vice-versa.
- MPI Reduces become MPI Broadcasts, and vice-versa.
- All PETSc forward scatters with insert flag become PETSc reverse scatters with accumulate flag.

---

# PART 5 — OPTIONS

Format: **`optionName`** — `type` = default — then description and (where applicable) enumerated values. Grouping headers are added for navigation only; option order and content follow the source page.

## I/O, grids, output

**`gridFile`** — `str` = `default.cgns` — The grid file to use. Must be a multiblock structured or overset CGNS file containing all block-to-block and boundary condition information.

**`restartFile`** — `str or list or None` = `None` — A single string or list of strings pointing to CGNS volume solution file(s) written by ADflow. Steady-state restart needs one file (not required as a single-item list). Unsteady restart typically needs a list of 2 items for second order restart.

**`meshSurfaceFamily`** — `str or None` = `None` — Custom family identifying the surface mesh to ADflow (sets `self.meshFamilyGroup`). If None, all wall surfaces are the surface mesh. Used for mesh warping and surface output.

**`designSurfaceFamily`** — `str or None` = `None` — Custom family defining the "design" surface (may be modified geometrically; sets `self.designFamilyGroup`). If None, all wall surfaces are the design surface. Distinguishes design nodes from regular nodes (regular nodes don't change with design but are still in `meshFamilyGroup` for warping).

**`closedSurfaceFamilies`** — `list or None` = `None` — Set of surfaces forming a closed surface; used in overset initialization to find wall cells. If None, all wall surfaces are used.

**`storeRindLayer`** — `bool` = `True` — Write halo/"rind cell" info into CGNS files. Required for some postprocessors (e.g. Tecplot) to compute contour lines correctly.

**`outputDirectory`** — `str` = `./` — Directory where output files are written.

**`outputSurfaceFamily`** — `str` = `allSurfaces` — The family included in surface output files.

**`writeSurfaceSolution`** — `bool` = `True` — Write surface solution automatically after each solution call.

**`writeVolumeSolution`** — `bool` = `True` — Write volume solution automatically after each solution call.

**`writeSolutionEachIter`** — `bool` = `False` — Write solution files after every nonlinear iteration (files depend on `writeTecplotSurfaceSolution`, `writeVolumeSolution`, `writeSurfaceSolution`). Useful for debugging/developing the solver.

**`writeTecplotSurfaceSolution`** — `bool` = `False` — Write a surface output in Tecplot format (less tweaking in Tecplot). Includes the zipper mesh and only active overset cells if a zipper mesh was created.

**`nSaveVolume`** — `int` = `1` — In unsteady runs, write volume output every this many unsteady iterations.

**`nSaveSurface`** — `int` = `1` — In unsteady runs, write surface output every this many unsteady iterations.

**`solutionPrecision`** — `str` = `single` — Precision when writing volume solution files. `single`: best if not using restart files (half the size; enough for visualization). `double`: best for restart files.

**`gridPrecision`** — `str` = `double` — Precision when writing the volume mesh in volume output files. `double`: preferred, especially RANS grids (lower precision harms boundary-layer cells). `single`: not typically used.

**`solutionPrecisionSurface`** — `str` = `single` — Precision for surface solution files. `single`: smaller file, enough for visualization. `double`: preferred for numerical analyses.

**`gridPrecisionSurface`** — `str` = `single` — Precision for the surface mesh in surface output files. `single`: smaller. `double`: preferred for numerical analyses.

**`isosurface`** — `dict` = `{}` — Type and values for isosurfaces; any `volumeVariables` may be used. Example `{"vx":-0.001, "shock":1.0}` places an isosurface at ~0 x-velocity and at shock sensor 1.

**`isoVariables`** — `list` = `[]` — State variables included on the isosurfaces requested via `isosurface`.

**`viscousSurfaceVelocities`** — `bool` = `True` — RANS/laminar NS only. Write surface velocities as the value on the first off-wall cell (surface velocities won't be zero). Enables oil-flow surface streamline postprocessing.

## Discretization & dissipation

**`discretization`** — `str` = `central plus scalar dissipation` — Discretization method (default recommended for robustness/speed at the cost of accuracy). `central plus scalar dissipation`: central FV with JST scalar dissipation. `central plus matrix dissipation`: central FV with JST matrix dissipation (may help poor meshes; minimal improvement on well-posed pyHyp meshes; for convergence issues lower `vis4` to 0.1 and lower CFL). `upwind`: upwind scheme.

**`coarseDiscretization`** — `str` = `central plus scalar dissipation` — Discretization for the coarse grid; generally same as `discretization`. Same options/notes as above.

**`limiter`** — `str` = `van Albada` — Flux limiter for the `upwind` scheme. `van Albada`; `minmod`; `no limiter`.

**`smoother`** — `str` = `DADI` — Smoother for the multigrid solver. `DADI`: D3ADI scheme (typically faster than RK, may be less robust). `Runge-Kutta`: five-stage fourth-order low-memory explicit RK time stepping.

**`dissipationScalingExponent`** — `float` = `0.67` — Exponent in the JST dissipation scheme. Typically unchanged; 2/3 is the theoretical best for an orthogonal 3-D grid.

**`acousticScaleFactor`** — `float` = `1.0` — Multiplies the acoustic contribution in the spectral-radius computation. Only for `central plus scalar dissipation`. Affects dissipation. 1.0 = standard JST weighting. <1.0 → more accurate for low Mach but too small slows convergence (setting it to the freestream Mach tends to work well for low Mach). >1.0 should not be used (compromises accuracy). See Seraj2023b.

**`vis4`** — `float` = `0.0156` — Coefficient of fourth-order dissipation (scalar and matrix JST). Default recommended if convergence is attainable; may be raised to 0.02–0.025 to help convergence at the cost of accuracy.

**`vis2`** — `float` = `0.25` — Coefficient of second-order dissipation (scalar and matrix JST). Only on at shocks; may be 0.0 for entirely subsonic simulations.

**`vis2Coarse`** — `float` = `0.5` — `vis2` used for the coarse grid (typically larger than `vis2`). Default 0.5 usually sufficient.

**`useDissContinuation`** — `bool` = `False` — Dissipation-based continuation (only `central plus scalar dissipation`). Adds second-order dissipation near shocks early, tapering off as a sigmoid of relative convergence. Helps converge supersonic and some transonic flows; no effect on purely subsonic. Solution must be well converged so the added dissipation reduces to a non-affecting value (for default `dissContMagnitude`/`dissContMidpoint`, L2 of 1e-6 suffices). See Seraj2023b.

**`dissContMagnitude`** — `float` = `1.0` — Magnitude of the additional dissipation. For large enough `dissContMidpoint`/`dissContSharpness` (defaults), initial dissipation = this value plus `vis2`. Values between a fifth and half of the freestream Mach typically work well.

**`dissContMidpoint`** — `float` = `3.0` — Orders of relative convergence at which the additional dissipation is half of `dissContMagnitude`. Default works well; should be at most 4.0 to avoid affecting a converged result for L2 of 1e-8.

**`dissContSharpness`** — `float` = `3.0` — Sharpness of the continuation sigmoid. Larger = steeper drop from initial to final dissipation. Should be at least 2.0 to avoid affecting a converged result for L2 of 1e-8.

## Equations, flow type, turbulence & SA constants

**`equationType`** — `str` = `RANS` — `RANS`; `Euler`; `laminar NS`.

**`equationMode`** — `str` = `steady` — `steady` (extensively tested); `unsteady` (not well tested); `time spectral` (extensively tested).

**`flowType`** — `str` = `external` — `external` (extensively tested); `internal` (not well tested).

**`turbulenceModel`** — `str` = `SA` — `SA`: Spalart-Allmaras — recommended for external aero and the **only** turbulence model that has been differentiated and the only one tested. Others: `SA-Edwards`; `k-omega Wilcox`; `k-omega modified`; `k-tau`; `Menter SST`; `v2f`.

**`turbulenceOrder`** — `str` = `first order` — `first order`: recommended (adjoint systems much easier to solve). `second order`: not typically used.

**`turbResScale`** — `float or list or None` = `None` — Affects how the total residual is scaled; set automatically by turbulence model. Defaults usually sufficient. Values can be a float scalar up to a 4-element list. Defaults/types: `SA` — float scalar, Default `10e4`; `SA-Edwards` — NOT IMPLEMENTED; `k-omega Wilcox` — NOT IMPLEMENTED; `k-omega modified` — NOT IMPLEMENTED; `k-tau` — NOT IMPLEMENTED; `Menter SST` — float list of 2, Default `[1e3, 1e-6]`; `v2f` — NOT IMPLEMENTED.

**`meshMaxSkewness`** — `float` = `1.0` — ADflow errors and fails if mesh skewness exceeds this. Only used when `useSkewnessCheck` is active.

**`useSkewnessCheck`** — `bool` = `False` — Compute per-cell skewness and error if above `meshMaxSkewness`. See `printBadlySkewedCells`.

**`turbulenceProduction`** — `str` = `strain` — Production term in the k equation (only for models with k). `strain`; `vorticity`; `Kato-Launder`.

**`useQCR`** — `bool` = `False` — Use the QCR2000 version of Spalart-Allmaras.

**`useRotationSA`** — `bool` = `False` — Use the rotation correction with Spalart-Allmaras.

**`useft2SA`** — `bool` = `True` — Include the ft2 term in Spalart-Allmaras.

**`eddyVisInfRatio`** — `float` = `0.009` — Free-stream eddy viscosity value. See turbmodels.larc.nasa.gov/spalart.html.

**`SAKappa`** — `float` = `0.41` — Von Karman constant in SA.

**`SAcb1`** — `float` = `0.1355` — SA constant cb1 (coefficient on the production term).

**`SAcb2`** — `float` = `0.622` — SA constant cb2 (coefficient on the cross-diffusion production term).

**`SAsigma`** — `float` = `0.66666666667` — SA constant sigma (cb3 in some formulations), the diffusion coefficient.

**`SAcv1`** — `float` = `7.1` — SA constant cv1 (used in the viscous damping function fv1).

**`SAcw2`** — `float` = `0.3` — SA constant cw2 (used in the destruction term function fw).

**`SAcw3`** — `float` = `2.0` — SA constant cw3 (used in the destruction term function fw).

**`SAct1`** — `float` = `1.0` — SA constant ct1 (used in the trip term ft1).

**`SAct2`** — `float` = `2.0` — SA constant ct2 (used in the trip term ft2).

**`SAct3`** — `float` = `1.2` — SA constant ct3 (used in the laminar suppression function ft2).

**`SAct4`** — `float` = `0.5` — SA constant ct4 (used in the laminar suppression function ft2).

**`SAcrot`** — `float` = `2.0` — SA constant crot (scaling factor for the rotation correction term).

## Wall treatment & misc physics

**`useWallFunctions`** — `bool` = `False` — Use wall functions. Generally not recommended (potentially very poor drag). Required routines are **not differentiated**, so wall-function simulations cannot be used for optimization.

**`useApproxWallDistance`** — `bool` = `True` — Cheap wall-distance calculation. Exact wall distances computed at initialization; parametric location of the closest wall point stored per cell. After deformation, the spatial search is skipped and the distance between new parametric location and new cell center is used. Faster; enables efficient wall-distance updates for aerostructural analysis.

**`updateWallAssociations`** — `bool` = `False` — Update wall associations even with approximate wall distance. Off by default because the update itself cannot be differentiated. For large mesh changes users may want it for more accurate wall distance.

**`eulerWallTreatment`** — `str` = `linear pressure extrapolation` — Inviscid wall BC implementation (default usually fine). `linear pressure extrapolation` (works with adjoint); `constant pressure extrapolation` (works with adjoint); `quadratic pressure extrapolation` (not tested); `normal momentum` (not tested).

**`viscWallTreatment`** — `str` = `constant pressure extrapolation` — Viscous wall BC implementation (default usually fine). `constant pressure extrapolation` (standard); `linear pressure extrapolation` (not typically used).

**`liftIndex`** — `int` = `2` — Coordinate index for the lift direction. `2`: y-axis; `3`: z-axis.

**`lowSpeedPreconditioner`** — `bool` = `False` — Low-speed preconditioner (previously for very low Mach; phased out in favor of ANK).

**`wallDistCutoff`** — `float` = `1e+20` — Distance cutoff beyond which wall distance is not computed (set to a large constant).

**`nearWallDist`** — `float` = `0.1` — Distance used to flag a cell as "nearWall".

**`infChangeCorrection`** — `bool` = `True` — Adjust the flowfield everywhere if the AoA changes between optimization iterations. When AoA changes after a converged run, residuals can stay low (initially only farfield boundary cells), confusing L2-based solver logic. This does a simple state-vector adjustment for free-stream changes (method via `infChangeCorrectionType`), run only if the free-stream state-change norm exceeds `infChangeCorrectionTol`. Supersedes the `RKReset` approach.

**`infChangeCorrectionTol`** — `float` = `1e-12` — Tolerance deciding whether to run the correction. Compute the change in the free-stream state vector and its norm; if above this, correct. Default suits the `"offset"` update; with `"rotate"`, use a higher tolerance (the rotation may be ill-conditioned for small state changes).

**`infChangeCorrectionType`** — `str` = `offset` — Correction method when `infChangeCorrection` is True. `offset` (default): add the delta free-stream state to every cell. `rotate`: rotate cell velocities by the free-stream velocity rotation and add the delta free-stream velocity to cell density and energy.

**`cavitationNumber`** — `float` = `1.4` — The −Cp value that triggers the cavitation sensor.

**`cpMinRho`** — `float` = `100.0` — rho parameter for the KS aggregation computing minimum Cp. Higher rho → closer to the actual min Cp, at the cost of a more nonlinear function.

## Iteration counts, CFL, multigrid, RK smoothing

**`nCycles`** — `int` = `2000` — Max number of "iterations". Counted differently per solver: multigrid → fine-grid cycles; ANK → steps plus KSP iterations; NK → function evaluations (mat-vec products or line search), each ≈ one residual evaluation. Terminates once the combined count exceeds `nCycles`.

**`timeLimit`** — `float` = `-1.0` — If positive, the max seconds the solver runs. Triggered internally by setting total iterations to max, so the solver reports analysis failed. Users must manage the fail flags. Useful when some cases take extremely long in optimizations.

**`nCyclesCoarse`** — `int` = `500` — Max iterations on the coarse grid during full-multigrid startup.

**`nSubiterTurb`** — `int` = `3` — Iterations of the turbulent DADI solver (RANS only). Increasing to 5–7 may lower total time for some cases.

**`nSubiter`** — `int` = `1` — RK or D3ADI iterations per nonlinear iteration (between turbulence updates). No effect for NK or ANK.

**`CFL`** — `float` = `1.7` — CFL number for Runge-Kutta simulations; the main speed/robustness parameter. Lower = more robust but slower. ⚠ *The description text states "The default value of 1.5 is a good place to start," which is inconsistent with the field default of 1.7 shown here — reproduced as published.* Some experimentation is needed for the max CFL.

**`CFLCoarse`** — `float` = `1.0` — CFL on the coarse multigrid grids; often desirable somewhat lower than the fine-grid CFL.

**`MGCycle`** — `str` = `sg` — Multigrid cycle type. Grid dimensions must permit the requested level. `sg` = single grid (no multigrid); `3w` = 3 levels with a 'w' cycle; `3v` = 'v' cycle; etc.

**`MGStartLevel`** — `int` = `-1` — Starting grid level for "full multigrid startup". Can significantly cut simulation time (good start from coarse approximate solutions). `-1` = coarsest level. For RANS, starting on the coarsest grid is often impossible if it has very few cells. Should be `1` for unsteady.

**`resAveraging`** — `str` = `alternate` — Frequency of residual averaging (smoothing) for the RK solver. `alternate`: every second RK stage (saves computation, minimal convergence impact). `never`; `always` (every RK stage).

**`smoothParameter`** — `float` = `1.5` — Parameter in residual smoothing; typically unchanged.

**`CFLLimit`** — `float` = `1.5` — Max CFL runnable without residual smoothing. If actual CFL is lower, no smoothing is applied regardless of `resAveraging`.

**`useBlockettes`** — `bool` = `True` — Speeds residual calculations via smaller data subsets ("blockettes") to reduce cache misses. Speedup ~2–3× (problem dependent). Not exhaustively tested and has caused errors before — verify it doesn't change your results.

**`useLinResMonitor`** — `bool` = `False` — Enable the linear residual monitor for NK and ANK (output resembles the adjoint solver's linear residual output).

**`restrictionRelaxation`** — `float` = `0.8` — Relaxation for the multigrid restriction operation (0.0–1.0). 1.0 = no relaxation. Default 0.80 works well across many cases.

## Overset

**`backgroundVolScale`** — `float` = `1.0` — Modifies the perceived quality of background cells; higher discourages background-mesh use during implicit hole cutting.

**`oversetProjTol`** — `float` = `1e-12` — Tolerance when projecting search points onto a block's surface definition.

**`overlapFactor`** — `float` = `0.9` — How much to derate cell volumes when seeking a donor (should be <1.0). With 0.9, a volume-1.0 cell won't take a volume-0.99 donor (0.99 > 1.0×0.9). Creates smoother cuts in similar-cell-size overlaps. Changing it generally doesn't change overlap much (valid-donor checks exist).

**`oversetLoadBalance`** — `bool` = `True` — Use load balancing for the implicit hole-cutting connectivity algorithm (Kenway2017a).

**`debugZipper`** — `bool` = `False` — Verbose printing/writing to debug zipper-mesh issues.

**`zipperSurfaceFamily`** — `str or None` = `None` — Surface families used when creating zipper meshes.

**`cutCallback`** — `function or None` = `None` — User Python function to explicitly blank cells by an arbitrary pattern (e.g. blank cells on the wrong side of a symmetry plane). Flagged cells are set as flood seeds. Signature: `def cutCallback(xCen, CGNSZoneNameIDs, cellIDs, flags):` filling `flags` (NumPy array). `xCen` = cell-center coordinates; `CGNSZoneNameIDs`, `cellIDs` help define the pattern. `flags` initialized to zero; 1 tags a cell for flooding. Example: `flags[xCen[:, 1] < 0] = 1` (flag cells with negative y-center for an x-z symmetry plane keeping the +y half).

**`explicitSurfaceCallback`** — `function or None` = `None` — User Python callback passing info for the explicit surface blanking method (requires understanding of overset flooding). For configs where automatic flooding is extremely hard, this surface-based blanking tags cells *inside* a user-provided closed geometry to reduce flood seeds/flooded cells. Cells with at least one vertex inside the surface are explicitly blanked (-4), adding two fringe-cell layers that prevent flooding in most cases. The callback receives `CGNSZoneNameIDs` (needed to map volume blocks to surfaces). It returns a dict keyed by bookkeeping names; each entry has `surfFile` (a closed plot3d surface; normals point outward; may be open at a symmetry plane), `blockIDs` (blocks to search/blank with this surface), optional `coordXfer` (a node coordinate transform; signature per DVGeometry's `addPointSet`), and optional `kMin` (only blank cells with k-index above this, to avoid tagging near-wall cells). See the source page for the full worked wing/fuselage example.

**`oversetUpdateMode`** — `str` = `frozen` — Overset connectivity update after warping. `frozen`: not updated after initialization (cheapest; whole mesh warped together, e.g. USMesh in IDWarp). `fast`: only weights updated, donors unchanged (whole mesh warped together; small design changes only — large changes raise an error from nonphysical weights). `full`: recomputed from scratch (component meshes warped independently, e.g. MultiUSMesh; hole cutting likely fails for large changes; even if it works, aero derivatives are inaccurate because hole-cutting routines are not differentiated).

**`nRefine`** — `int` = `10` — Max refinement-loop iterations in the implicit hole-cutting method.

**`nFloodIter`** — `int` = `-1` — Flooding iterations per hole-cutting iteration (debugging meshes that flood completely on the first iteration). `-1` = run until flooded cells stop changing; a positive value stops flooding at that count.

**`useZipperMesh`** — `bool` = `True` — Use a zipper mesh. Required to accurately compute integrated quantities (lift, drag) when surface meshes overlap.

**`useOversetWallScaling`** — `bool` = `False` — Modify perceived cell quality during hole cutting to favor cells aligned with viscous walls.

**`selfZipCutoff`** — `float` = `120.0` — Cutoff angle (degrees) in the self-zip step of zipper-mesh gap triangulation. Triangulation triggers if the angle between adjacent substrings is below the cutoff. Only for the first self-zip iteration at each node; 90° for subsequent iterations.

**`oversetPriority`** — `dict` = `{}` — Modify perceived cell quality per block. Keys = block IDs, values = quality multipliers. Higher discourages, lower encourages that block. May be required to get flooding working.

**`recomputeOverlapMatrix`** — `bool` = `True` — Recompute the domain overlap matrix on a full overset update. The matrix defines domain connections for the donor search (a domain must overlap the current one to be searched). It changes on full updates; if outdated, donor search can fail where overlapping domains aren't represented. Set False for efficiency if the user knows it won't change.

**`oversetDebugPrint`** — `bool` = `False` — Debug printout when hole cutting fails. Not useful anymore; instead obtain the volume solution after failure and plot cells with incomplete connectivities (iblank = -5).

## Unsteady & time-spectral

**`timeIntegrationScheme`** — `str` = `BDF` — Time integration for unsteady analysis. `BDF`: backward differentiation formula (second order BDF is currently the only option known to work). `explicit RK`; `implicit RK`.

**`timeAccuracy`** — `int` = `2` — Order of the time integration scheme. `2` (second order); `1`; `3`.

**`nTimeStepsCoarse`** — `int` = `48` — Coarse-mesh time steps (only for periodic problems where full multigrid is possible).

**`nTimeStepsFine`** — `int` = `400` — Time steps in an unsteady simulation.

**`deltaT`** — `float` = `0.01` — Time step for unsteady simulation.

**`useALE`** — `bool` = `True` — Use the Arbitrary Lagrangian-Eulerian formulation for unsteady mesh deformations.

**`useGridMotion`** — `bool` = `False` — Whether a rigid body grid motion has been specified.

**`coupledSolution`** — `bool` = `False` — ADflow coupled to other solvers, but only in unsteady mode. Previously for aerothermoelastic analysis.

**`timeIntervals`** — `int` = `1` — Number of "spectral instances" for a time spectral simulation. Only meaningful when `equationMode` is `time spectral`.

**`alphaMode`** — `bool` = `False` — Specified alpha motion for time spectral.

**`betaMode`** — `bool` = `False` — Specified beta motion for time spectral. Untested.

**`machMode`** — `bool` = `False` — Specified Mach number motion for time spectral. Untested.

**`pMode`** — `bool` = `False` — Specified p-motion (roll) for time spectral. Untested.

**`qMode`** — `bool` = `False` — Specified q-motion (pitch) for time spectral.

**`rMode`** — `bool` = `False` — Specified r-motion (yaw) for time spectral. Untested.

**`altitudeMode`** — `bool` = `False` — Specified h-variation motion for time spectral. Untested.

**`windAxis`** — `bool` = `False` — Wind axis for time spectral stability derivative computations.

**`alphaFollowing`** — `bool` = `True` — Whether alpha follows the body in p, q, r mode.

**`TSStability`** — `bool` = `False` — Compute time spectral stability info from a time-spectral CFD solution.

**`useTSInterpolatedGridVelocity`** — `bool` = `False` — Use spectral differentiation to compute grid velocity for the time spectral equation.

**`useExternalDynamicMesh`** — `bool` = `False` — Use an externally provided deformed mesh for the time spectral equation (each time instance gets its own deformed grid; aerodynamic or aeroelastic).

## Convergence

**`L2Convergence`** — `float` = `1e-08` — Desired convergence factor. For multigrid, relative to the initial fine-grid residual (actual convergence vs freestream residual may be 1–2 orders lower). For NK, reference is the freestream residual.

**`L2ConvergenceRel`** — `float` = `1e-16` — Typically only with an aerostructural solver. Relative tolerance vs the current starting point.

**`L2ConvergenceCoarse`** — `float` = `0.01` — Convergence factor on coarse grids during multigrid startup. Most benefit after 2–3 orders, so typically 1e-2 to 1e-3.

**`maxL2DeviationFactor`** — `float` = `1.0` — If the solver runs out of iterations, the max factor the residual may exceed the `L2Convergence` target and still be "converged".

## Newton–Krylov (NK) options

**`useNKSolver`** — `bool` = `False` — Enable the NK solver. If False, all "NK" options have no effect.

**`NKSwitchTol`** — `float` = `1e-05` — Relative tolerance before switching to NK. Must be low enough that difficult transients have passed. If NK stalls, lower it (run startup longer before switching).

**`NKSubspaceSize`** — `int` = `60` — GMRES subspace size for NK. Increase for difficult problems at the expense of memory.

**`NKLinearSolveTol`** — `float` = `0.3` — Initial tolerance for the Newton linear system. Only used for the first iteration if `NKUseEW` is True.

**`NKUseEW`** — `bool` = `True` — Use the Eisenstat-Walker algorithm for per-iteration linear convergence. If False, always target `NKLinearSolveTol`. If True, only the initial iteration targets it; later iterations pick a tolerance from the previous nonlinear convergence (avoids wasted computation far from the solution; tightens near it). If working, NK converges in a few to tens of iterations. It may repeatedly pick very low tolerances and loop — then disable this (or better, switch to NK later).

**`NKADPC`** — `bool` = `False` — Use the AD version of the NK preconditioner (vs finite-difference). AD improves preconditioner accuracy and linear performance at higher preconditioner cost.

**`NKViscPC`** — `bool` = `False` — Use the full viscous stencil for the NK preconditioner (memory/cost penalty, more accurate). A smaller-stencil preconditioner is usually good enough.

**`NKGlobalPreconditioner`** — `str` = `additive Schwarz` — Like `globalPreconditioner` but for NK. `additive Schwarz` (restricted ASM); `multigrid` (AMG).

**`NKASMOverlap`** — `int` = `1` — Overlap levels in the additive Schwarz preconditioner for NK. More = stronger but costlier/more memory. Typical 1 (easy) up to 2–3 (difficult).

**`NKASMOverlapCoarse`** — `int` = `0` — Same as `NKASMOverlap` but coarse levels when using `multigrid` for `NKGlobalPreconditioner`.

**`NKPCILUFill`** — `int` = `2` — Levels of fill for the local ILU factorization in NK. Typical 1 (easy) up to 3 (difficult). More = stronger preconditioner (fewer linear iterations) but costlier/more memory.

**`NKPCILUFillCoarse`** — `int` = `0` — Same as `NKPCILUFill` but coarse levels (multigrid `NKGlobalPreconditioner`).

**`NKJacobianLag`** — `int` = `20` — Frequency of NK preconditioner reformation (Jacobian "lagged" by this many iterations). Simple problems may never reform; difficult ones may need lower (more assemblies, which are costly).

**`applyPCSubspaceSize`** — `int` = `10` — Only for aerostructural analysis. Subspace size and total iterations when ADflow only preconditions residuals via `globalNKPrecon`.

**`NKInnerPreconIts`** — `int` = `1` — Local preconditioning iterations for NK. Leave at 1 unless very difficult.

**`NKInnerPreconItsCoarse`** — `int` = `1` — Same as `NKInnerPreconIts` but coarse levels (multigrid `NKGlobalPreconditioner`).

**`NKOuterPreconIts`** — `int` = `1` — Global preconditioning iterations for NK. Typical 1–3.

**`NKAMGLevels`** — `int` = `2` — Like `adjointAMGLevels` but for NK.

**`NKAMGNSmooth`** — `int` = `1` — Like `adjointAMGNSmooth` but for NK.

**`NKLS`** — `str` = `cubic` — NK line search. `cubic` (cubic interpolation); `none` (no line search); `non-monotone`.

**`NKFixedStep`** — `float` = `0.25` — Step size for NK with no line search (`"NKLS": "none"`).

**`RKReset`** — `bool` = `False` — Run `nRKReset` RK startup iterations when using NK and restarting from a converged solution (lets residual increase; prevents NK stalling).

**`nRKReset`** — `int` = `5` — Number of RK startup iterations when `RKReset` is True.

## Approximate Newton–Krylov (ANK) options

**`useANKSolver`** — `bool` = `True` — Enable the ANK solver.

**`ANKUseTurbDADI`** — `bool` = `True` — Solve turbulence with the DADI solver (decoupled ANK only). If False, an internal ANK-like solver (`turbKSP`) is used in decoupled mode.

**`ANKUseApproxSA`** — `bool` = `False` — ANK switches from approximate to exact SA implicit formulation when it switches first→second order (printed SANK or CSANK). Set False to force the approximate SA treatment always. No effect on the final solution, but can help convergence on challenging cases.

**`ANKSwitchTol`** — `float` = `1000.0` — Relative convergence before switching to ANK. Default ensures ANK is used from the first iteration, even after design changes.

**`ANKSubspaceSize`** — `int` = `-1` — Positive value sets the ANK subspace size; otherwise the max iteration value is used.

**`ANKMaxIter`** — `int` = `40` — Max linear iterations per ANK step.

**`ANKLinearSolveTol`** — `float` = `0.05` — Linear solver tolerance for ANK.

**`ANKLinearSolveBuffer`** — `float` = `0.01` — *Do not modify unless you understand PETSc linear tolerances.* PETSc's absolute linear convergence must be lower than the final residual convergence norm (the linear residual is a linearization of the nonlinear residual; solving to an absolute linear tolerance doesn't guarantee the same nonlinear reduction). This option offsets the ANK absolute linear tolerance and matters only very near the final L2 target. E.g. `0.01` forces tighter linear convergence when absolute linear convergence is reached before the relative target. Lowering to 1e-4/1e-6 may give "extra" nonlinear convergence in the final iteration at more linear cost; raising to 1e-1 reduces this buffer (ends closer to the target) but too large may cause issues in the final iterations needing several small Newton steps.

**`ANKLinResMax`** — `float` = `0.1` — Adaptively adjusts CFL so the ANK linear residual doesn't exceed this. Target linear convergence may not be reached per iteration (fine), but CFL is managed so the post-iteration linear residual stays below this — the max allowed linear residual per iteration.

**`ANKGlobalPreconditioner`** — `str` = `additive Schwarz` — Like `globalPreconditioner` but for ANK. `additive Schwarz` (restricted ASM); `multigrid` (AMG).

**`ANKASMOverlap`** — `int` = `1` — Like `NKASMOverlap` but for ANK.

**`ANKASMOverlapCoarse`** — `int` = `0` — Coarse-level version (multigrid `ANKGlobalPreconditioner`).

**`ANKPCILUFill`** — `int` = `2` — Like `NKPCILUFill` but for ANK.

**`ANKPCILUFillCoarse`** — `int` = `0` — Coarse-level version (multigrid `ANKGlobalPreconditioner`).

**`ANKJacobianLag`** — `int` = `10` — Max nonlinear iterations between preconditioner updates (an adaptive algorithm may update more often; this is the upper limit). See Yildirim2019b.

**`ANKInnerPreconIts`** — `int` = `1` — Like `NKInnerPreconIts` but for ANK.

**`ANKInnerPreconItsCoarse`** — `int` = `1` — Coarse-level version (multigrid `ANKGlobalPreconditioner`).

**`ANKOuterPreconIts`** — `int` = `1` — Like `NKOuterPreconIts` but for ANK.

**`ANKAMGLevels`** — `int` = `2` — Like `adjointAMGLevels` but for ANK.

**`ANKAMGNSmooth`** — `int` = `1` — Like `adjointAMGNSmooth` but for ANK.

**`ANKCFL0`** — `float` = `5.0` — Initial CFL for ANK.

**`ANKCFLMin`** — `float` = `1.0` — Coefficient in the ANK minimum-CFL algorithm (the minimum CFL rises as the sim converges; this is the initial value). See Yildirim2019b.

**`ANKCFLLimit`** — `float` = `100000.0` — Upper CFL limit for ANK. Larger = better nonlinear convergence at the cost of robustness and expensive linear solves.

**`ANKCFLFactor`** — `float` = `10.0` — Multiplication factor in the ANK CFL ramping algorithm. See Yildirim2019b.

**`ANKCFLExponent`** — `float` = `0.5` — Exponent in the ANK CFL ramping algorithm. See Yildirim2019b.

**`ANKCFLCutback`** — `float` = `0.5` — Cutback factor for the ANK CFL when adaptively reduced.

**`ANKCFLReset`** — `bool` = `True` — Reset the ANK CFL to `ANKCFL0` at the first ANK iterations. Each aeroproblem saves its last ANK CFL; if False, subsequent solutions with the same aeroproblem start at that last CFL; if True, always reset to `ANKCFL0`.

**`ANKStepFactor`** — `float` = `1.0` — Initial step size taken by ANK (initializes the line searches; actual step can be lower but not higher).

**`ANKStepMin`** — `float` = `0.01` — Minimum ANK step size. If a smaller step is required, the CFL is cut back.

**`ANKConstCFLStep`** — `float` = `0.4` — If the ANK step size is below this, don't ramp CFL (keep constant). If the last step size is larger (progress is being made), ramp the CFL.

**`ANKPhysicalLSTol`** — `float` = `0.2` — ANK physical line search parameter. 0.2 = physical parameters (density, energy) may change 20% per iteration per cell; the line search adjusts step size to obey this globally.

**`ANKPhysicalLSTolTurb`** — `float` = `0.99` — Physicality check for the turbulence model in ANK. The turbulence variable may increase freely but only decrease by this fraction per iteration. 0.99 = up to 99% decrease per iteration; since the initial value is positive, this effectively prevents negatives.

**`ANKUnsteadyLSTol`** — `float` = `1.0` — Reduction factor in the ANK unsteady line search. 1.0 = accept any step that doesn't increase the unsteady residual norm. See Yildirim2019b.

**`ANKSecondOrdSwitchTol`** — `float` = `1e-16` — Relative convergence target where ANK switches to a second-order-accurate implicit formulation. For transonic, set just after shocks settle (~1e-4 to 1e-6).

**`ANKCoupledSwitchTol`** — `float` = `1e-16` — Relative convergence target where ANK switches to a coupled turbulence formulation (flow and turbulence solved together — better nonlinear convergence, less robust). Use after 5–6 orders of nonlinear convergence.

**`ANKTurbCFLScale`** — `float` = `1.0` — Scaling of the CFL used for the separate turbulence solver in ANK.

**`ANKUseFullVisc`** — `bool` = `True` — Use full viscous terms in the ANK implicit formulation (True = the R1 residual approximation level). See Yildirim2019b.

**`ANKPCUpdateTol`** — `float` = `0.5` — If ANK converges by this amount relative to the last preconditioner-update iteration, update the preconditioner again.

**`ANKPCUpdateCutoff`** — `float` = `1e-16` — Cutoff below which the PC-update tolerance is adjusted based on nonlinear convergence. When L2 is above this option, ANK PC updates trigger if relative convergence since the last update is below `ANKPCUpdateTol`; when L2 is below this option, `ANKPCUpdateTolAfterCutoff` is used. Useful for CSANK, which can converge L2 by multiple orders per iteration near the end (no need to update PC every iteration).

**`ANKPCUpdateTolAfterCutoff`** — `float` = `0.0001` — ANK PC update tolerance used after reaching `ANKPCUpdateCutoff` convergence (before the cutoff, `ANKPCUpdateTol` is used).

**`ANKADPC`** — `bool` = `False` — Use the AD-based preconditioner for ANK. The finite-difference preconditioner is usually good enough and cheaper.

**`ANKNSubiterTurb`** — `int` = `1` — Turbulent subiterations per ANK iteration if ANK turbKSP is used (`"ANKUseTurbDADI": False`). Keep at 1 for most cases (turbKSP is more effective and expensive).

**`ANKTurbKSPDebug`** — `bool` = `False` — Debug printout from the ANK turbKSP solver. If turbDADI is disabled, turbKSP is used; without this option True you get no turbulence-solver diagnostics.

**`ANKUseMatrixFree`** — `bool` = `True` — Use the matrix-free Jacobian for ANK. If False, use the approximate Jacobian (formulated for the preconditioner matrix) as the implicit formulation.

**`ANKCharTimeStepType`** — `str` = `None` — Characteristic-based preconditioning of the ANK time-step matrix; accelerates low-Mach convergence (speedup depends on how low local Mach is). Mainly useful for freestream Mach ≤ 0.4; usually combine with `acousticScaleFactor`. See Seraj2023b. `None`: no characteristic time stepping. `Turkel`: Turkel preconditioner (all flows, applied only to locally subsonic flow). `VLR`: van Leer-Lee-Roe (locally subsonic and supersonic; mainly useful subsonic; usually slightly better than Turkel).

## Partitioning & load balancing

**`blockSplitting`** — `bool` = `True` — Whether blocks may be split for better load balancing.

**`loadImbalance`** — `float` = `0.1` — Allowable load imbalance between processors when mapping blocks. 0.1 = 10%.

**`loadBalanceIter`** — `int` = `10` — METIS graph partitioning iterations. More = better load balancing but tends to split blocks more (communication penalty).

**`partitionOnly`** — `bool` = `False` — Only run partitioning, not the flow solution (to check load balancing without a CFD solve).

**`partitionLikeNProc`** — `int` = `-1` — If positive, partition as if using that many processors. Recreate a partition-related issue seen on many procs using fewer. Purely a debug option; partitions will exceed procs, giving a load imbalance ~1.0 (bad for performance).

## Printing, monitoring & solution variables

**`numberSolutions`** — `bool` = `True` — Attach the AeroProblem numbering to the grid solution file.

**`writeSolutionDigits`** — `int` = `3` — Digits in solution output filenames (4 → 0023, 3 → 023).

**`printIterations`** — `bool` = `True` — Print monitoring values at each iteration.

**`printTiming`** — `bool` = `True` — Print time for initialization, flow solution, function evaluations, and adjoint solution.

**`printIntro`** — `bool` = `True` — Print the ADflow intro message during initialization.

**`printAllOptions`** — `bool` = `True` — Print all options during initialization.

**`setMonitor`** — `bool` = `True` — Monitor the adjoint iterations.

**`printWarnings`** — `bool` = `True` — Print warnings (e.g. bad-quality volumes).

**`printNegativeVolumes`** — `bool` = `False` — Print block indices, cell-center coordinates, and volume for each negative-volume cell.

**`printBadlySkewedCells`** — `bool` = `False` — Print block indices, cell-center coordinates, and skewness for each cell above `meshMaxSkewness` (only when `useSkewnessCheck` is active).

**`monitorVariables`** — `list` = `['cpu', 'resrho', 'resturb', 'cl', 'cd']` — Variables whose convergence is monitored. Options: `resrho` (density residual); `resmom` (momentum residuals); `resrhoe` (total energy residual); `resturb` (turbulence residuals); `cl`; `clp` (pressure part of cl); `clv` (viscous part of cl); `cd`; `cdp`; `cdv`; `cfx`, `cfy`, `cfz` (force coefficients); `cmx`, `cmy`, `cmz` (moment coefficients); `hdiff` (max relative difference between H and Hinf); `mach` (max Mach); `yplus` (max y+); `eddyv` (max eddy/laminar viscosity ratio).

**`surfaceVariables`** — `list` = `['cp', 'vx', 'vy', 'vz', 'mach']` — Variables written to the CGNS surface solution file: `rho`, `p`, `temp`, `vx`, `vy`, `vz`, `cp`, `ptloss` (relative total pressure loss), `mach`, `blank` (cell iblank values). For viscous flows also: `cf` (skin friction magnitude), `cfx`, `cfy`, `cfz`, `ch` (Stanton number), `yplus`, `forceInDragDir`, `forceInLiftDir` (force coefficient per unit area projected in drag/lift direction).

**`volumeVariables`** — `list` = `['resrho']` — Variables (in addition to restart variables) written to the CGNS volume solution file: `mx`, `my`, `mz` (momentum), `rhoe` (total energy), `temp`, `vort` (vorticity magnitude), `vortx`, `vorty`, `vortz`, `cp`, `mach`, `macht` (turbulent Mach), `ptloss`, `eddy` (eddy viscosity), `eddyratio`, `dist` (wall distance to nearest viscous wall), `resrho`, `resmom`, `resrhoe`, `resturb`, `blank`.

**`storeConvHist`** — `bool` = `True` — Save convergence history into Fortran arrays (accessible via `getConvergenceHistory`).

**`forcesAsTractions`** — `bool` = `True` — Return tractions instead of forces.

## Adjoint solver & preconditioner

**`adjointL2Convergence`** — `float` = `1e-06` — Adjoint tolerance relative to the residual for a zero initial guess. With zipper meshes the adjoint can stall early; running until the absolute residual reaches a minimum should still give accurate derivatives.

**`adjointL2ConvergenceRel`** — `float` = `1e-16` — Adjoint tolerance relative to the residual at the start of the adjoint call (including a possible non-zero restart).

**`adjointL2ConvergenceAbs`** — `float` = `1e-16` — Adjoint absolute tolerance.

**`adjointDivTol`** — `float` = `100000.0` — Relative amount the adjoint residual may increase before the method is declared diverging.

**`adjointMaxL2DeviationFactor`** — `float` = `1.0` — Like `maxL2DeviationFactor` but for the adjoint; applied only to `adjointL2Convergence`.

**`approxPC`** — `bool` = `True` — Use the approximate Jacobian for the adjoint preconditioner.

**`ADPC`** — `bool` = `False` — Use AD for the adjoint preconditioner.

**`viscPC`** — `bool` = `False` — Keep cross-derivative terms in the adjoint preconditioner.

**`useDiagTSPC`** — `bool` = `True` — Include off-time-instance terms in the time spectral adjoint preconditioner.

**`restartAdjoint`** — `bool` = `True` — Restart the adjoint from the previous solution.

**`adjointSolver`** — `str` = `GMRES` — Linear solver for the adjoint (see PETSc). `GMRES` (best performance); `TFQMR` (transpose-free QMR); `Richardson` (preconditioned Richardson); `BCGS` (BiCGStab); `IBCGS` (improved BiCGStab).

**`adjointMaxIter`** — `int` = `500` — Max iterations for the adjoint solution.

**`adjointSubspaceSize`** — `int` = `100` — Krylov subspace size for the adjoint solution.

**`GMRESOrthogonalizationType`** — `str` = `modified Gram-Schmidt` — Orthogonalization for GMRES (mostly affects the adjoint; also used in ANK/NK linear solvers). `modified Gram-Schmidt` (best accuracy vs complex step, and performance); `CGS never refine` (fast, not most accurate); `CGS refine if needed` (slow, inaccurate); `CGS always refine` (slow, inaccurate).

**`adjointMonitorStep`** — `int` = `10` — Write the adjoint residual norm every this many iterations.

**`dissipationLumpingParameter`** — `float` = `6.0` — Scaling for dissipation lumping in the approximate adjoint preconditioner.

**`preconditionerSide`** — `str` = `right` — Side to apply the adjoint preconditioner. `right`; `left`.

**`matrixOrdering`** — `str` = `RCM` — Matrix ordering for the adjoint preconditioner (see PETSc). `RCM` (reverse Cuthill-McKee); `natural`; `nested dissection`; `one way dissection`; `quotient minimum degree`.

**`globalPreconditioner`** — `str` = `additive Schwarz` — Global preconditioner for the adjoint system. `additive Schwarz` (restricted ASM); `multigrid` (AMG — much faster than ASM for large meshes; smoothing and coarse solves use standard ASM+ILU).

**`localPreconditioner`** — `str` = `ILU` — Local preconditioner for the adjoint. `ILU` (incomplete LU — currently the only supported option).

**`ILUFill`** — `int` = `2` — Levels of fill for the local ILU factorization in the adjoint. Typical 1 (easy) up to 3 (difficult). More = stronger (fewer iterations) but costlier/more memory.

**`ILUFillCoarse`** — `int` = `0` — Same as `ILUFill` but coarse levels (multigrid `globalPreconditioner`).

**`ASMOverlap`** — `int` = `1` — Overlap levels in the additive Schwarz preconditioner for the adjoint. Typical 1 (easy) up to 2–3 (difficult).

**`ASMOverlapCoarse`** — `int` = `0` — Same as `ASMOverlap` but coarse levels (multigrid `globalPreconditioner`).

**`innerPreconIts`** — `int` = `1` — Local preconditioning iterations for the adjoint. More may help difficult problems (more time each).

**`innerPreconItsCoarse`** — `int` = `1` — Same as `innerPreconIts` but coarse levels (multigrid `globalPreconditioner`).

**`outerPreconIts`** — `int` = `3` — Global preconditioning iterations for the adjoint. More may help difficult problems (more time each). Default sufficient for most.

**`adjointAMGLevels`** — `int` = `2` — Levels for the adjoint solver's algebraic multigrid preconditioner.

**`adjointAMGNSmooth`** — `int` = `1` — Smoothing iterations per level for the adjoint solver's AMG preconditioner.

**`applyAdjointPCSubspaceSize`** — `int` = `20` — Krylov subspace size for the adjoint preconditioner.

**`frozenTurbulence`** — `bool` = `False` — Use the frozen-turbulence assumption in the adjoint (neglect the turbulence-model linearization). Only Spalart-Allmaras is ADed. May help convergence for high transonic flows, but the resulting sensitivity is **less accurate**.

**`useMatrixFreedrdw`** — `bool` = `True` — Use matrix-free routines for mat-vec products with the full Jacobian. If False, form the full matrix exactly (matrix-based mat-vec). **Should always be True for overset meshes**, otherwise gradients are inaccurate (known bug — github.com/mdolab/adflow/issues/204). Matrix-free has lower memory; runtime is a tradeoff between Jacobian assembly and adjoint solution times. See Kenway2019a.

**`skipAfterFailedAdjoint`** — `bool` = `True` — If an adjoint fails in the current sensitivity evaluation, skip the rest for efficiency (use `checkAdjointFailure` for the correct fail flag). Next design evaluation retries. If False, all adjoints are solved to whatever tolerance possible and partially converged solutions are used for total derivatives (user decides if acceptable).

## Adjoint debug (Tapenade verification)

**`firstRun`** — `bool` = `True` — Adjoint debugging only. Setting False turns on the Tapenade debugger.

**`verifyState`** — `bool` = `True` — Adjoint debugging only. Verifies dRdw.

**`verifySpatial`** — `bool` = `True` — Adjoint debugging only. Verifies dRdx.

**`verifyExtra`** — `bool` = `True` — Adjoint debugging only. Verifies dIda.

## Separation & cavitation sensors

**`computeSepSensorKs`** — `bool` = `False` — Compute KS-based separation cost functions (`sepSensorKs`, `sepSensorKsArea`). If not True, the code returns zero for these. When True there is a small overhead (first finds the conventional max sensor value over the surface to compute the aggregated max).

**`sepSensorKsRho`** — `float` = `1000.0` — rho parameter for the KS aggregation computing the maximum separation-sensor value in `sepSensorKs`/`sepSensorKsArea`. A higher rho approaches the actual max separation sensor more accurately… **[description truncated on source fetch — see source URL for the remainder]**

> **⚠ Gap — not captured (page truncated on fetch).** The following options appear in the source page's table of contents after `sepSensorKsRho` but their descriptions were not retrievable through the fetcher. They are listed here for completeness; consult the source Options page for details. Nothing below was invented.
> - **`sepSensorOffset`** — [description not captured]
> - **`sepSensorKsOffset`** — [description not captured]
> - **`sepSensorSharpness`** — [description not captured]
> - **`sepSensorKsSharpness`** — [description not captured]
> - **`sepSensorKsPhi`** — [description not captured]
> - **`cavSensorOffset`** — [description not captured]
> - **`cavSensorSharpness`** — [description not captured]
> - **`cavExponent`** — [description not captured]
> - **`computeCavitation`** — [description not captured]

---

*End of compiled reference. Source: ADflow documentation (MDO Lab), retrieved 6 July 2026. Referenced internal papers cited by the docs include Yildirim2019b (ANK startup), Kenway2017a (parallel overset), Kenway2019a (effective adjoint approaches), and Seraj2023b (dissipation/time-step scaling for low/high Mach).*
