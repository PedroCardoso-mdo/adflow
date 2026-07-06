# ADFLOW 03 — Concordance (paper math ↔ code option ↔ docs)

**Use for:** you see an option and want the algorithm behind it, or you read the math and want the flag. This is the bridge between `01`/`02` (theory) and `05` (options).

**How to read:** "Paper concept" = Yildirim2019 / Kenway2019 notation. "ADflow option" = name in `05`. "Default" = code default. `⇄` = "same thing, different name." **Verify option names/defaults against `05`/the running code** before trusting them in a script — defaults drift.

---

## A. Flow solver / ANK crosswalk (Yildirim2019 ⇄ Docs)

| Paper concept (Yildirim2019) | ADflow option (`05`) | Default | Notes / gotcha |
|---|---|---|---|
| Matrix-free operator level `R_m`; `R1`=approximate, `R0`=exact | `ANKSecondOrdSwitchTol` | `1e-16` (≈ off) | Rel-conv threshold: **above** → approximate (`R1`); **below** → exact/2nd-order (`R0`). Output prints `S` when 2nd-order active. **The code exposes a 2-way switch (approx↔exact); the paper's `R2` is the *preconditioner* Jacobian, not a run-time matrix-free choice.** |
| Preconditioner `M = I/Δt + ∂R2/∂Q` (first-order, 7-block) | (automatic; built from R2) + local `ILUFill`, `ASMOverlap` | ILU 2, overlap 1 | You don't pick R2; you tune its ILU/ASM strength. |
| Adaptive **PC lagging** (`η_pc=0.5`, `n_lag=10`) | automatic; iterations marked `*` | — | Not "Jacobian lagging" — the matrix-free op stays current. Refresh triggers: GMRES cap hit after a prior success, or `η_pc` reached, or `n_lag` cap. |
| Linear tolerance `η_lin` `(12)` | `ANKLinearSolveTol` | `0.05` | Fixed by design (EW deliberately **not** used in ANK). |
| GMRES iteration cap `n_lin` `(§3.5)` | ANK linear iteration limit (see `05`) | 40 | Hitting the cap = "linear failed" but partial step is still taken. |
| Pseudo-time init `CFL(0)=5` | `ANKCFL0` | `5.0` | — |
| CFL ceiling `CFL_max=1e5` | `ANKCFLLimit` | `1e5` | Lower it (not below a few hundred) to stabilize tiny-step cases. |
| mRDM/SER CFL ramp `(18)–(20)` (`α=10, β=0.5, θ_red=0.5`) | automatic PTC | — | **CFL updates only on PC-refresh iterations** (the `*` ones); other iterations hold CFL. |
| Physicality `θ_phys,ρ=θ_phys,E=0.2` `(15)` | `ANKPhysicalLSTol` | `0.2` | Raise toward (not past) 1.0 for supersonic/moving shocks (0.4–0.6). |
| Physicality `θ_phys,ν̃=0.99` (neg. updates only) | `ANKPhysicalLSTolTurb` | `0.99` | Lets SA grow, forbids negative `ν̃`. |
| Backtracking acceptance `R_u ≤ ‖R0‖` `(17)`, `θ_bt=0.7` | `ANKUnsteadyLSTol` (accept factor); `θ_bt` internal | `1.0` | **`ANKUnsteadyLSTol=1.0` ⇄ the paper's exact criterion (17).** Raise to ~1.5 to *allow* an unsteady-residual increase (helps over "hills," risks divergence). |
| Decoupled turbulence (block Gauss–Seidel, Alg 3) | default; solver select `ANKUseTurbDADI` | `True` (`turbDADI`) | `False` → `turbKSP` (isolated ANK for turbulence; better for >1M / overset). |
| Turbulence sub-iterations | `nSubiterTurb` (turbDADI) / `ANKNSubiterTurb` (turbKSP) | 3–7 rec. | Massive-separation cases up to 10. |
| Coupled turbulence (standalone Newton) | `ANKCoupledSwitchTol` | `1e-16` (≈ off) | Switch **after** the SA "hill" (~4 orders); before it → stall. |
| Turbulence residual scale `×1e4` (§2.4) | `turbResScale` | `10e4` (SA) | **Numeric mismatch to flag:** paper says ×`1e4`; code default `10e4`. Recommendation: **don't touch it**; keep decoupled ANK until printed turb residual sits ~3–5 orders below flow residuals (that's the sign the scaling is right). |
| Convergence tolerance `η_abs=1e-12` `(7)` | `L2Convergence` | 1e-12-class | Final solution/gradient accuracy vs cost. |
| Switch to NK at `η_ANK=1e-5` | `NKSwitchTol` | ~1e-5 | Lower (1e-6…1e-8) for hard cases; premature switch fails (see `04`). |
| NK: exact `R0`, no time term, EW tolerance | `useNKSolver`, `NKUseEW`, `NKLinearSolveTol`, `NKLS` | EW on, 0.3, cubic | EW upper limit hard-coded 0.8; consistently hitting 0.8 ⇒ state too far ⇒ lower `NKSwitchTol`. |
| Startup selection (ANK vs MG) | `useANKSolver`, `ANKSwitchTol`, `MGStartLevel`, `MGCycle` | ANK 1.0 | Overset / separated flow ⇒ ANK (MG has no coarse levels / stalls). |

## B. Adjoint / differentiation crosswalk (Kenway2019 ⇄ Docs)

| Paper concept (Kenway2019) | ADflow option (`05`) | Default | Notes / gotcha |
|---|---|---|---|
| Jacobian-free GMRES adjoint `(26)` (reverse-AD mat-vec, no stored `[∂R/∂w]ᵀ`) | `useMatrixFreedrdw` | `True` | **Must be `True` for overset** (matrix-based drdw gives wrong overset gradients — known bug #204). Matrix-free = lower memory; runtime = assembly vs solve trade-off. |
| Right-preconditioned adjoint `(27)` | `preconditionerSide` | `right` | — |
| Nested PC: **ASM** global | `globalPreconditioner` | `additive Schwarz` | `multigrid` (AMG) much faster for large meshes. |
| Nested PC: **ILU** local | `localPreconditioner`, `ILUFill`, `ILUFillCoarse` | ILU, 2, 0 | Fill 1 (easy) → 3 (hard); more = stronger but costlier. |
| ASM overlap | `ASMOverlap`, `ASMOverlapCoarse` | 1, 0 | 1 (easy) → 2–3 (hard). |
| **Outer** Richardson (ASM) = 3 iters | `outerPreconIts`, `outerPreconItsCoarse` | 3 | Outer > inner in effectiveness (cross-block info). |
| **Inner** Richardson (ILU) = 2 iters | `innerPreconIts`, `innerPreconItsCoarse` | 1 | (Paper reports 2; code default 1 — tune per problem.) |
| Matrix ordering (RCM) | `matrixOrdering` | `RCM` | PETSc orderings. |
| Krylov restart size | `applyAdjointPCSubspaceSize` | 20 | — |
| **No** frozen turbulence (differentiate `ν̃`) for accurate gradients | `frozenTurbulence` | `False` | `True` speeds hard transonic adjoints but gives **less accurate** sensitivities. Only SA is AD'd. |
| Source-transformation AD (Tapenade), master routine (Algs 1–4) | (build system) | — | AD-tool verification via `firstRun/verifyState(dRdw)/verifySpatial(dRdx)/verifyExtra(dIda)`. |
| FD-Jacobian variant ≈ 0.1% gradient error | — (implementation choice) | — | Prefer Jacobian-free/AD for production; complex-step for reference. |
| Skip remaining adjoints after a failure | `skipAfterFailedAdjoint`, `checkAdjointFailure` | `True` | `False` uses partially-converged adjoints for totals (user judges acceptability). |
| Dissipation lumping in approximate adjoint PC | `dissipationLumpingParameter` | 6.0 | — |
| Adjoint AMG controls | `adjointAMGLevels`, `adjointAMGNSmooth` | 2, 1 | Used when `globalPreconditioner=multigrid`. |

## C. Ambiguities & gotchas (read before trusting a mapping)

1. **`R1`/`R2`/`R0` are not three run-time knobs.** In the code, the ANK matrix-free operator toggles between **approximate (`R1`)** and **exact (`R0`)** via `ANKSecondOrdSwitchTol`. **`R2` is the preconditioner's first-order Jacobian** (7-block stencil), tuned indirectly through ILU/ASM — you never select `R2` as the mat-vec operator at run time. Confusing these is the #1 mismatch between reading the paper and reading a runscript.
2. **"Second order" vs "exact."** Docs/output say "second order" (`S` flag) for what the paper calls the exact/`R0` Jacobian. Same object.
3. **Steady vs unsteady residual.** The line search reduces the **unsteady** residual (`ANKUnsteadyLSTol` governs its acceptance); the printed L2 is the **steady/total** residual and *can rise* while unsteady falls. Do not read a transient total-residual bump as divergence.
4. **`turbResScale`**: paper `1e4`, code default `10e4` for SA. Treat as "auto-set, don't touch" unless the printed turb/flow residual gap is far from ~4 orders.
5. **CFL cadence:** ANK CFL changes **only** on preconditioner-refresh iterations (`*`). If CFL "looks stuck," check whether the PC is being refreshed, not the CFL logic.
6. **Matrix-free is mandatory on overset** for *both* primal (`useMatrixFreedrdw` context) and correct adjoints. On overset, `R2`/matrix-based preconditioning ignores overset stencils by design — accepted, because the matrix-free operator carries the exact coupling.
7. **Default-value drift:** several defaults above (e.g. inner Richardson 1 vs paper's 2, `ANKSecondOrdSwitchTol` off-by-default) reflect code choices that differ from the paper's *experiments*. When a number matters, the code/`05` wins over the paper.
