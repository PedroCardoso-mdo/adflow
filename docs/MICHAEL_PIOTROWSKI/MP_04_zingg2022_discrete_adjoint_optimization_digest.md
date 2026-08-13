# MP_04 digest — Piotrowski & Zingg (2022), Discrete-adjoint shape optimization

Companion to [`MP_04_zingg2022_discrete_adjoint_optimization_full.md`](MP_04_zingg2022_discrete_adjoint_optimization_full.md).

The optimization paper: SA-sLM2015 fully coupled to a discrete-adjoint
gradient-based ASO framework (Jetstream), applied to lift-constrained drag
minimizations. Relevant to us for **adjoint verification practice** and for
**what free-transition optimization demands of the mesh** — not for the primal
solver, which is covered better by MP_03/MP_06.

## 1. Adjoint verification — their ladder vs ours

Their partial derivatives are **hand-linearized analytically**, with some
entries by complex step. Verification, in their words:

- Flow Jacobian `∂R/∂Q` — verified against **complex step**.
- Metric linearization `∂R/∂C` — verified against **complex step**.
- Total gradient — verified against **second-order finite difference**.
- `∂R/∂X` (design variables incl. **α and sideslip**, the latter formed by
  complex step) is verified *concurrently* with `∂R/∂C`, because design
  variables are built from the B-spline control points.

Two structural points worth checking against our own adjoint:

- **The residual depends on the control points through the off-wall spacing**,
  not only through the grid metrics — "the residuals depend on the B-spline
  control points through the grid metrics **and the off-wall spacing**, the
  latter of which is used in the turbulence and transition model source terms."
  Wall distance `d` appears throughout the transition source terms (`Re_v`,
  `F_θt`, `δ`), so any geometric-derivative path that treats `d` as constant is
  wrong. Worth confirming in `VERIF_00_three_stage_verification.md` territory.
- **Their verification case is chosen to contain every relevant flow feature**:
  a blunt-TE RAE2822 extruded at 25° sweep with periodic BCs, run with
  **SA-sLM2015cc including crossflow** — i.e. the verification case exercises
  the crossflow correlations and compressibility corrections, not just 2D TS
  transition. Our dot-product/fast_b tests should be checked for the same
  coverage (crossflow terms live or dead in the tested configuration).

The adjoint systems themselves: flow adjoint by **GCROT** (simplified/flexible
variant), mesh adjoint by preconditioned CG, design updates by **SNOPT**.

## 2. The finding that costs the most compute: streamwise resolution

> "increasing the streamwise resolution, and therefore better resolving the
> large gradients in the boundary-layer transition region, **improves the
> capability of the optimization algorithm to delay boundary-layer
> transition**."

And the scaling law behind it:

> "**As the transition length decreases with increasing Reynolds number, finer
> streamwise mesh spacings are required**… Despite increasing the streamwise
> grid resolution, drag minimizations at the Dash8-Q400 and sweep-corrected
> 737-800 design conditions **fail to delay transition as far aft** as the
> stability-analysis-based framework of Rashad & Zingg."

The failure mode when streamwise resolution is short is *not* a crash — the
optimizer simply **cannot find the laminar-extent gain**, because the
transition-region gradients are under-resolved and so is the gradient of drag
w.r.t. shape. This is a silent accuracy failure that looks like a converged
optimization.

**For us:** any AR5/sickle optimization campaign is subject to this. It also
re-frames the mesh-family work (`04_mesh_families`, `05_mesh_independence`) —
the relevant refinement direction for transition is **streamwise near the
front**, and the requirement **tightens with Reynolds number**. Cf. the
sickle-wing chord-resolution misreading already recorded in memory.

## 3. Multimodality

The Cessna 172R conditions produced **at least two distinct local minima** — a
thin reflexed trailing edge with steep recovery, versus an aft-loaded design;
the aft-loaded one wins. Free-transition airfoil design spaces are multi-modal,
so a single gradient-based run is not the answer. Relevant if we ever report
optimization results.

## 4. Optimality convergence at high Re — an open problem, not a solved one

From the thesis recommendations (MP_06) on the same results: the optimizer
"struggles to reduce optimality at higher Reynolds numbers", and increasing
transition length does not fix it. Their suspected cause is **noise in the
gradient**, specifically from "different formulations of the **local pressure
gradient parameter**" (`λ_θ`). Worth remembering before we attribute a noisy
SA-GR gradient to our AD wiring — the model itself has a known noisy input.

## 5. Solver settings restated here

Same architecture as MP_03: SBP + SAT, matrix dissipation (Swanson–Turkel),
first-order upwind on turbulence/transition, fully-coupled Newton-Krylov-Schur
with PTC, approximate-Schur PC. Nothing new; MP_03 and MP_06 are the better
sources.
