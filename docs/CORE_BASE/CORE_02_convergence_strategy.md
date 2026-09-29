# Converging transition cases (SA-GR, SA-sγ)

Validated recipes per model/case, the options that must be set, known limits,
and what has been tried and falsified. Option defaults quoted here are the
`adflow/pyADflow.py` defaults; code wins if they ever disagree.

Run-side deliverable for the 3D plain wing (recipe, per-phase restarts,
phase-entry runner `run_strategy.py`):
`~/Desktop/Run/MDO_PhD/Transition/gama_rethetha/03_convergence_strategy/3d_plain_wing/best_strategy/`.

## Recipes

### SA-GR, 3D wing (plain wing, AR5 family, tutorial wing)

| Phase | Enter at (rel totalRes) | Key options | Behaviour |
|---|---|---|---|
| ANK segregated | start | `ANKUseTurbDADI: True` | flow converges to rel ~1e-5; Re̅θt residual parks (expected). Faster variant: SANK via `ANKSecondOrdSwitchTol 1e-4`, never couple (needs its own leg) |
| CANK | `ANKCoupledSwitchTol: 1e-5` | `ANKADPC: True`, LS below | full steps to rel ~1e-7; kills the Re̅θt residual in ~20 iterations |
| CSANK | `ANKSecondOrdSwitchTol: 1e-6` | same | to rel ~3.5e-8 (one order past CANK) |
| NK | `nkswitchtol` ≈ 4e-8 (just above CSANK's floor) | `NKADPC: True`, `NKSubspaceSize` 200–300 | one full step at entry, then the deep-NK wall below rel ~5e-9 |

On large meshes NK at the correct point merely re-attains CSANK's depth;
`useNKSolver: False` and letting CSANK finish is equally good.

### SA-GR, sickle wing (crossflow)

CANK at rel 1e-2 → CSANK at 1e-4, **no NK**, `eddyVisInfRatio 5e-7`
(default 0.009). This is the one case where early coupling is validated; the
"do not couple early" rule below is from the plain wing.

### SA-GR, 2D optimisation (NACA 0012, L0/L1)

`L2Convergence 1e-8` + `nkswitchtol 1e-8` + `ANKCFLLimit 1e8` (default 1e5).
CSANK carries every solve to 1e-8 in 230–350 iterations and NK never engages.
With `nkswitchtol 1e-6` NK pinned at `Step 0.00` on almost every solve and
the optimiser carried on regardless, because **pySLSQP ignores the `fail`
flag**. Always check `fail` in the history. The only remaining `fail`s are far
line-search trials, which the optimiser rejects anyway.

### SA-sγ (`SA-noft2-Gamma`), tutorial wing

This recipe is used by `tests/reg_tests/reg_sgamma.py`:
- ANK → SANK early: `ANKSecondOrdSwitchTol 1e-2`.
- **Never couple** (`ANKCoupledSwitchTol 1e-20`): CANK pins the step.
- `ANKCFLLimit 1e8`, `ANKCFL0 5`.
- LS `ANKUnsteadyLSTol 2.0`, `ANKPhysicalLSTol 0.8`, `ANKPhysicalLSTolTurb 0.99`.
- `ANKADPC`/`NKADPC: True`, `NKSubspaceSize 300`, `nkswitchtol 1e-8`.

## Non-negotiable options (SA-GR)

- `ANKUnsteadyLSTol: 1.5`, `ANKPhysicalLSTol: 0.5` (defaults 1.0 / 0.2).
  These are the fix for coupled-phase step collapse (steps 0.01 → 1.00). The
  defaults stagnate, and 2.0 / 0.7 gives bit-identical results.
- `solutionPrecision: "double"`. Single-precision restarts truncate the
  transition front and poison every restart.
- `ANKADPC` / `NKADPC: True` (default False). The FD-coloured PC is unusable
  for SA-GR Newton phases (lin res ~0.99).
- `turbResScale` is left at its auto value: SA-GR `[1e4, 0.1, 1e-4]`, SA-sγ
  `[1e4, 0.1]`.
- First-order transition advection (`transitionFirstOrderUpwind`, default on).
- `NKLSRelax` (default True): NK cubic line search with Armijo alpha 1e-3
  and turb-blowup pre-limit factor 3.0. Off gives upstream 1e-2 / 2.0, where
  NK sits at minlambda on SA-GR. Do not relax further: 1e-4, or a factor of
  5.0, lets through steps that crash, because NK has no ρ/E physicality check.

## Rules and diagnostics

- **NK engaged too early is the most expensive mistake.** With
  `nkswitchtol 1e-6`, every AR5 level where NK engaged stalled at rel ~1e-6
  (`Step = 0.00`, CFL `----`, totalRes creeping up in the 12th digit).
  Restarting with `useNKSolver: False` reached 8 orders in 2–27 iterations.
  If NK pins at `Step = 0.00`, suspect early engagement first. Never restart
  NK at the same threshold.
- **Do not couple early (plain wing).** CANK at rel 1e-2 or 1e-4 stagnates,
  even with the LS relaxations, and 1e-5 is optimal there. The sickle recipe
  above is the exception.
- **Do not hold CANK past ~1e-7 or CSANK past ~3.5e-8.** A front adjustment
  kicks the residual and the phase oscillates permanently. Hand over before
  the kick.
- **A crossflow plateau can be the front settling, not a stall.** The sickle
  medium grid sat flat at ~4.8 orders for 15 h with clean health indicators
  (Step 1.00, lin res ~0.05), then dropped 20× in 80 iterations. Before
  killing a flat crossflow run, check whether γ and Re̅θt are still
  evolving.
- **A genuine stall without crossflow** shows as a few γ cells pinned at their
  bound, whose updates Algorithm-2 damping deletes (`damp=`/`wf=` on
  STALLDIAG). `srcDtDeactivateIters` (default 5) never reactivates Eq. 59 on
  a flat floor. Raising it to 100000 moves the pinned residual slightly but
  does not break the floor.
- After a deep restart the ANK CFL re-ramps from
  `ANKCFLMin·(totalR0/totalR)^0.5`. Raise `ANKCFLMin` if the run is stranded.
- `ANKNSubiterTurb` only acts in the turbKSP path. With `ANKUseTurbDADI: True`
  only `nSubiterTurb` (default 3) reaches DADI.
- `scancel --signal=USR2` ends the current *solve*, not the job. A runner with
  several staged `CFDSolver(ap)` calls continues to the next one, so follow
  with a plain `scancel`.

## Known limits

- **Deep-NK wall** below rel ~5e-9: the NK linear solve saturates (lin res
  0.8–0.99, GMRES exhausted, 60–200 evals/iteration). It depends on how
  settled the field is: a fully settled field reached rel 6.4e-11 before
  hitting the same wall. No option moves it (table below). It needs PC code
  work.
- NK has no ρ/E physicality check (ANK has one); this is what caps `NKLSRelax`.
- **Memory is driven by ILU fill, not the NK subspace.** At 7.42M cells on one
  node, ILU(1) fits at 3 GB/rank, while ILU(2) and ILU(3) OOM during CANK,
  before NK engages. More ranks means more total memory (per-rank PC slices
  and halos), so scale ranks with nodes.

## Falsified levers

| Lever | Result |
|---|---|
| `NKUseEW False` + `NKLinearSolveTol 0.3` | GMRES still returns 0.998. The PC, not EW, starves deep NK |
| `ANKPCILUFill` / `NKPCILUFill` 3 | Worse. NK stalls (lin res 0.996) where ILU(2) gives 0.76. Keep 2 (ILU(1): slower, no stall) |
| `NKJacobianLag 5` | No effect. PC quality, not staleness, is the limit |
| LS 2.0 / 0.7 vs 1.5 / 0.5 (SA-GR) | Bit-identical |
| Early NK (`nkswitchtol` 1e-5 / 1e-6) | Stalls production runs (see rules) |
| Bundle of stronger-PC NK knobs at once | First step blew the transition residual up 125×. Change one knob at a time |
| `transitionRowVolScale` (Eq. 58 S_r) | Stalls the NK linear solve. Keep off |
| `transitionResidualAutoscale` (Eq. 58 S_a) | Marginal: progresses but noisier. Keep off |
| `ANKCFLLimit` as a lever on the plain wing | Not one: collapses coincide with front adjustments, not CFL growth |

Untried, no memory cost: `NKOuterPreconIts` / `NKInnerPreconIts`.

## ADflow vs the paper's solver (P&Z §IV)

| Aspect | Paper (Newton–Krylov–Schur) | ADflow on this branch |
|---|---|---|
| Globalisation | Fully coupled approximate Newton from iteration 1 | Segregated ANK (DADI or turbKSP) → coupled CANK at `ANKCoupledSwitchTol` |
| Endgame | Inexact Newton at rel ~1e-5, to machine zero | Matrix-free NK at `nkswitchtol`; CSANK as an in-ANK second-order mode |
| Linear PC | Approximate Schur | ILU/ASM, FD-coloured or AD-assembled (`ANKADPC`/`NKADPC`) |
| Step control | Per-node bounds-triggered damping (Alg. 2) + backtracking | Global λ line search. Alg. 2 per-node damping on γ/Re̅θt in DADI and NK (`applyNKAlgorithm2Damping`); no ρ/E equivalent |
| Source stiffness | Eq. 59 Δt restriction, off after 5 clean Newton steps, back on after a backtrack | Same (`transitionSrcDtRestrict`, `srcDtDeactivateIters`). Additive in DADI; MAX form in turbKSP/CANK; additive diagonal in NK (`applyNKSrcDtDiagonal`, preconditioner only, not the true residual) |
| Scaling | Eq. 58 row + column + auto | `turbResScale` rows + NK/turbKSP column scaling. Eq. 58 S_r/S_a exist but stay off |
| Convection of ν̃/γ/Re̅θt | First-order upwind | Same |

The iteration gap is structural: global λ instead of per-node damping, and a
weaker PC. One paper iteration is one well-solved Newton step, so compare wall
time rather than outer iteration counts.
