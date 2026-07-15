# ADflow solver algorithms vs. Piotrowski & Zingg (2020) — where they differ and why it matters

Purpose: explain why the paper converges SA-sLM2015 cases to machine zero in
few iterations while ADflow's solver hierarchy (ANK/SANK/CANK/CSANK/NK)
stalls or oscillates on the same physics. Grounded in the paper's §IV
(Solution Methodology) and in the empirical campaign of 2026-07-14/15
(`3D_Plain_Wing/long_overnight/DECISIONS.md`, `ab_matrix/RESULTS.md`).

## 1. The two architectures side by side

| Aspect | Paper (Newton–Krylov–Schur, Diablo) | ADflow on this branch |
|---|---|---|
| Phase 1 (globalization) | **Fully coupled approximate-Newton**: implicit Euler + local Δt, first-order analytic Jacobian, from iteration 1 | Segregated ANK (turb via DADI or turbKSP) and/or coupled CANK (first-order Jacobian); coupling entered at `ANKCoupledSwitchTol` |
| Phase 2 (endgame) | **Inexact Newton**: true second-order Jacobian (analytic or FD), switch at rel. residual ~1e-5, converges to machine zero | NK (matrix-free second-order), switch at `nkswitchtol`; CSANK exists as an in-ANK second-order mode |
| Linear solver | GMRES + **approximate-Schur** parallel PC | GMRES + ILU-based global PC (FD-coloring or AD-assembled `ANKADPC`/`NKADPC`) |
| Step control | Per-node, per-variable **bounds-triggered damping** (Alg. 2) + unsteady backtracking line search | **One global λ** for the whole field (physical LS + unsteady LS pick the worst-cell step) |
| Source stiffness | **Eq. 59 source-term Δt restriction** (λ_src·Δt ≤ 0.9, QR eigenvalues of the 3×3 source Jacobian); active in phase 1, deactivated after 5 clean Newton iters, **reactivated on backtracking** | Ported to this branch (`transitionSrcDtRestrict`, coupled path, incl. the 5-iter deactivation) — but **no reactivation inside NK** (ADflow NK has no pseudo-transient machinery at all) |
| Linear-system scaling | Eq. 58 row+column+auto scaling; fixed normalizations ν̃_max=1e3, γ_max=10, Re̅θt_max=1e4 | `turbResScale` (residual/row side); NK **column scaling** added 2026-07-13 on this branch; no residual autoscaling |
| Monitored residual | Partially scaled residual S_r·R | Unscaled per-equation L2 norms + `totalRes` |
| Convection of SA/γ/Re̅θt | First-order upwind (deliberately dissipative; helps keep variables in bounds) | Same on this branch (paper §IV.A followed) |

## 2. The three differences that actually explain the iteration gap

### a) Step control: global scalar vs. local damping (the big one)

ADflow's ANK/CANK picks a single λ ∈ [0.01, 1] for **all** cells and
equations, set by the most violating cell through the physical/unsteady line
searches. At a moving transition front, a handful of cells always want a
tiny update, so the entire 175k-cell update gets throttled to λ = 0.01–0.05
— observed for hours in `full_dp.log` iters 330–480 (retheta residual
oscillating 50 → 620 while CFL was cut 9.5e3 → 1.18e3).

The paper never damps globally. Algorithm 2 takes the full Newton update
everywhere and then, **only at nodes where a transition variable leaves its
bounds** (γ outside [1e-10, 2], Re̅θt < 20), multiplies that node's
transition update by 0.99^m until it re-enters. With Eq. 59 active the
bounds are "rarely reached", i.e. the damping is a rarely-firing safety
valve — 99.9% of the field takes λ = 1.

Caution from this branch: a naive "per-cell λ" port (damping front cells
every iteration, not bounds-triggered) was tested 2026-07-14 (`paper_mimic`
run 2) and was WORSE than global λ — it desynchronizes the update. The
paper's version is different in both trigger (bounds violation only) and
scope (transition variables only, per node). A faithful re-implementation is
still an open code item.

### b) Newton with a parachute: source-Δt reactivation

The paper switches to Newton at rel. 1e-5 and survives because when a
Newton step needs backtracking, the source-term time stepping (Eq. 59)
switches back ON — Newton temporarily becomes pseudo-transient again, takes
a safe step, then hands back to pure Newton. ADflow's NK has no equivalent:
a rejected step just backtracks the same global direction, which is exactly
the 200-evals/iteration, λ = 0.01–0.03, lin-res 0.8 thrash measured at rel.
~1e-6 in `nk_now.log` (iter 14) and `mimic_run5`. This is why on this branch
NK only survives when engaged deep (~1e-6 rel., front nearly settled), while
the paper engages at 1e-5 routinely.

### c) Linear solve quality at the front

The paper: approximate-Schur PC + *stricter* linear tolerances for
transition cases (§IV.B: "often requiring stricter linear solver tolerances
relative to fully turbulent simulations"), full second-order Jacobian
available analytically. ADflow: the FD-colored PC is unusable for SA-GR
Newton (lin res 0.99, `mimic_run4`); the AD-assembled PC (`NKADPC`) fixes
engagement (lin res 0.3) but degrades to 0.8 as the front moves and the
lagged PC (NKJacobianLag 20) goes stale. Note: naively bundling stronger-PC
options (`NKViscPC`, `NKJacobianLag 3`, ILU fill 3, outer its 2, ASM overlap
2) made the first step catastrophically wrong on 2026-07-15 (`nk_pc1.log`)
— untangle one knob at a time.

## 3. Differences that turned out to matter less (on this case)

- **Coupled-from-iteration-1**: the paper does it, but with Alg. 2 + Eq. 59
  + their scaling. In ADflow, early coupling (CANK at 1e-2 or even 1e-4)
  reproducibly ends in the global-λ stagnation of §2a. ADflow's own docs
  prescribe the opposite for small coupled steps: couple later. Empirically
  (2026-07-15), segregated ANK+DADI with the wall-BC fix converges the
  transition residuals fast from freestream (retheta res 3.4e4 → 35 in ~60
  iters), so the coupled mid-phase may be unnecessary here.
- **Iteration counting**: one paper iteration = one Newton step with a
  well-solved linear system (many inner GMRES its). Raw outer-iteration
  comparisons overstate the gap; wall-time comparisons are the honest ones.
  The structural gaps of §2 are real either way.

## 4. What this branch already ported from the paper

- Eq. 59 source-term Δt restriction in the coupled ANK path
  (`transitionSrcDtRestrict`, srcLambda frozen per iteration, §IV.B.3
  five-clean-iterations deactivation counter). Missing: reactivation inside
  NK (needs NK-side machinery that does not exist in ADflow).
- First-order upwind convection for γ/Re̅θt (and SA via
  `transitionUseApproxSA`).
- NK column scaling (analogue of the S_c part of Eq. 58) — 2026-07-13.
- Bounds/clipping philosophy: ADflow clips transition variables; the paper
  explicitly warns hard clipping "can potentially lead to stalling" and
  prefers Alg. 2 damping. Worth revisiting where the branch clips.

## 5. Open code items, in impact order (2026-07-15 view)

1. **Faithful Algorithm 2** (bounds-triggered per-node damping of transition
   updates, γ ∈ [1e-10, 2], Re̅θt ≥ 20, θ_fac = 0.99) as an alternative to
   the global physical LS for the transition rows.
2. **Eq. 59 reactivation in NK** (on backtracking / residual rise above the
   phase-switch level) — would let NK engage shallower than 1e-6.
3. **Residual autoscaling / monitored scaled residual** (Eq. 58 S_a, S_r)
   — affects switch tolerances and totalRes composition; today's `totalRes`
   mixes equations with arbitrary relative weights (`turbResScale`).
4. Stronger NK PC for the front — one knob at a time from a
   double-precision checkpoint (`NKJacobianLag` first; `NKViscPC` suspect).
