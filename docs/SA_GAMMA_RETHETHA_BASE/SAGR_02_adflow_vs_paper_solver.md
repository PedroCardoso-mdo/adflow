# ADflow solver algorithms vs. Piotrowski & Zingg (2020) — where they differ and why it matters

Purpose: explain why the paper converges SA-sLM2015 cases to machine zero in
few iterations while ADflow's solver hierarchy (ANK/SANK/CANK/CSANK/NK)
stalls or oscillates on the same physics. Grounded in the paper's §IV
(Solution Methodology) and in the empirical campaign of 2026-07-14/15
(`03_convergence_strategy/3d_plain_wing/_old/campaign_2026-07-14_to_16/` — `TESTS_AND_CONCLUSIONS.md`,
`long_overnight/DECISIONS.md`, `archive/ab_matrix/RESULTS.md`; deliverable
ladder + restart states in `03_convergence_strategy/3d_plain_wing/best_strategy/`).

## 1. The two architectures side by side

| Aspect | Paper (Newton–Krylov–Schur, Diablo) | ADflow on this branch |
|---|---|---|
| Phase 1 (globalization) | **Fully coupled approximate-Newton**: implicit Euler + local Δt, first-order analytic Jacobian, from iteration 1 | Segregated ANK (turb via DADI or turbKSP) and/or coupled CANK (first-order Jacobian); coupling entered at `ANKCoupledSwitchTol` |
| Phase 2 (endgame) | **Inexact Newton**: true second-order Jacobian (analytic or FD), switch at rel. residual ~1e-5, converges to machine zero | NK (matrix-free second-order), switch at `nkswitchtol`; CSANK exists as an in-ANK second-order mode |
| Linear solver | GMRES + **approximate-Schur** parallel PC | GMRES + ILU-based global PC (FD-coloring or AD-assembled `ANKADPC`/`NKADPC`) |
| Step control | Per-node, per-variable **bounds-triggered damping** (Alg. 2) + unsteady backtracking line search | **Global λ** line search (`nkls`: none/cubic/non-monotone, `NKSolvers.F90:682-691`) **plus** Algorithm 2's per-node γ/Re̅θt damping, ported into NK 2026-07-16 (`applyNKAlgorithm2Damping`, called right after the LS accepts a step — §4). Remaining gap: no equivalent per-node physicality check for ρ/E (§7) — only γ/Re̅θt are covered |
| Source stiffness | **Eq. 59 source-term Δt restriction** (λ_src·Δt ≤ 0.9, QR eigenvalues of the 3×3 source Jacobian); active in phase 1, deactivated after 5 clean Newton iters, **reactivated on backtracking** | Ported to this branch (`transitionSrcDtRestrict`, coupled path, incl. the 5-iter deactivation) **and** reactivation-on-backtrack inside NK, ported 2026-07-16 (`noBacktrackCount`/`srcDtDeactivateIters`, §4) |
| Linear-system scaling | Eq. 58 row+column+auto scaling; fixed normalizations ν̃_max=1e3, γ_max=10, Re̅θt_max=1e4 | `turbResScale` (residual/row side); NK **column scaling** added 2026-07-13 on this branch; no residual autoscaling |
| Monitored residual | Partially scaled residual S_r·R | Unscaled per-equation L2 norms + `totalRes` |
| Convection of SA/γ/Re̅θt | First-order upwind (deliberately dissipative; helps keep variables in bounds) | Same on this branch (paper §IV.A followed) |

## 2. The three differences that actually explain the iteration gap

### a) Step control: global scalar vs. local damping (RESOLVED for NK 2026-07-16 — historical account below)

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
scope (transition variables only, per node). **This is the version that was
faithfully ported into NK on 2026-07-16** (`applyNKAlgorithm2Damping`, §4) —
bounds-triggered, γ/Re̅θt only, applied after the global LS accepts a step,
confirmed to fire silently (no "exhausted" warnings) in smoke testing,
matching the paper's description of it as a rarely-firing safety valve.

**So this is no longer the open gap it once was.** The global-λ line search
(`nkls`) and the per-node Algorithm 2 damping now coexist in NK exactly as
in the paper: LS picks one global step, then Algorithm 2 individually backs
off any cell whose γ/Re̅θt left bounds. If NK is still converging poorly,
the remaining suspects are the ones actually still open (§7, §8): NK has
**no** per-node physicality check for density/energy (unlike ANK's
`physicalityCheckANK`), so an over-permissive line-search `alpha` can accept
a globally-stepped update that damages ρ/E with nothing to catch it before
Algorithm 2 even runs on γ/Re̅θt; and PC staleness at the moving front
(§8.4) has since been partially tested — JacLag 5 had no effect; the PC's quality is the limit (convergence-strategy.md, 2026-08-07/12).

### b) Newton with a parachute: source-Δt reactivation (RESOLVED 2026-07-16 — historical account below)

The paper switches to Newton at rel. 1e-5 and survives because when a
Newton step needs backtracking, the source-term time stepping (Eq. 59)
switches back ON — Newton temporarily becomes pseudo-transient again, takes
a safe step, then hands back to pure Newton. At the time this was written,
ADflow's NK had no equivalent: a rejected step just backtracked the same
global direction, which is exactly the 200-evals/iteration, λ = 0.01–0.03,
lin-res 0.8 thrash measured at rel. ~1e-6 in `nk_now.log` (iter 14) and
`mimic_run5`. This was why NK only survived when engaged deep (~1e-6 rel.,
front nearly settled).

**Reactivation-on-backtrack was ported into NK 2026-07-16** (same
`noBacktrackCount`/`srcDtDeactivateIters` globals ANK uses, §4), and §7
records that NK *survived* (no crash) at the paper's ~1e-5 engagement point
in ONE 35-iteration run on a 175k-cell restart (`run_early_engage.log`).

> **⚠️ Correction (2026-08-12, per commits `ef9fc10d`/`54c43475`):** "no
> crash" is not "converges". On the AR5 corrected-foil family (0.46M–7.42M
> cells) even `nkswitchtol = 1e-6` — LATER than 1e-5 — engaged NK ~1.5
> orders early and stalled every engaged level at `Step = 0.00`. Early NK
> engagement is NOT a solved problem: the validated handover is at CSANK's
> floor (~4.2e-8), or NK disabled entirely. See
> `docs/CORE_BASE/CORE_02_convergence_strategy.md` ("Engaging NK too early is the most
> expensive mistake").

If convergence is still poor from a shallow restart, look at the other two
open items (§7's missing physicality check, §8's PC quality) rather than
this one.

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
  five-clean-iterations deactivation counter).
- Eq. 59 reactivation-on-backtrack **inside NK** — 2026-07-16.
  `NKStep` (`src/NKSolver/NKSolvers.F90`)
  now updates the same `noBacktrackCount`/`srcDtDeactivateIters` globals
  ANK uses: a backtracked step (`stepMonitor < 1`) or outright line-search
  failure resets the counter, reactivating the restriction for the next
  `srcDtDeactivateIters` clean Newton steps. Unlike ANK there is no
  totalR-vs-`ANKSecondOrdSwitchTol` leg (NK-only module, can't see
  `ANKSolver`'s module-local tol; a backtrack already catches a residual
  rise in practice). The restriction itself is injected as an additive
  diagonal (`applyNKSrcDtDiagonal`, new subroutine) on the **already
  assembled, already column-scaled** `dRdwPre`, mirroring exactly where
  ANKStep's `timeStepMat` gets `MatAXPY`'d in — it never touches
  `setupStateResidualMatrix` (shared with the frozen adjoint) or the true
  residual (`FormFunction_mf`/`computeResidualNK`), so line-search
  acceptance and the adjoint are both untouched.
- **Faithful Algorithm 2 inside NK** — 2026-07-16. New
  `applyNKAlgorithm2Damping` (`NKSolvers.F90`), called from `NKStep` after
  the line search accepts `work` and before it becomes the new state:
  per-node exponential back-off (θ=`transitionDampTheta`, same
  option/constants DD-ADI already uses — `rsaGRgammaLo/Hi`,
  `rsaGRreThetaLo`) on γ and Re̅θt **only**, bounds-triggered, matching the
  paper's Algorithm 2 pseudocode exactly. Operates in NK's column-scaled
  space; the physical bound check divides by the column scale. Confirmed
  via smoke test to fire silently (no "exhausted" warnings) — matches the
  paper's description of it as a rarely-firing safety valve when Eq. 59 is
  active.
- Eq. 58 **S_r geometric row-scaling** and **S_a autoscale proxy** —
  implemented 2026-07-16, both gated by their own off-by-default options
  (`transitionRowVolScale`, `transitionResidualAutoscale`). See §5 for
  results — **S_r stalls the NK linear solve badly; not recommended.**
  S_a makes real progress but is noisier than baseline; marginal.
- First-order upwind convection for γ/Re̅θt (and SA via
  `transitionUseApproxSA`).
- NK column scaling (analogue of the S_c part of Eq. 58) — 2026-07-13.
- Bounds/clipping philosophy: ADflow clips transition variables; the paper
  explicitly warns hard clipping "can potentially lead to stalling" and
  prefers Alg. 2 damping. Worth revisiting where the branch clips.

## 5. Eq. 58 S_r/S_a results, and a serious infra bug found along the way (2026-07-16)

**The `.pyf` option-wiring bug.** `src/f2py/adflow.pyf` is a hand-maintained
f2py interface file, not regenerated automatically from the Fortran source.
It was stale — missing `transitionSrcDtRestrict` and its siblings (added
before this session) and everything added this session
(`transitionNK`, `transitionRowVolScale`, `transitionResidualAutoscale`).
f2py's `fortran`-type Python module objects **silently accept and store
arbitrary attribute names** with zero connection to the real Fortran
memory (verified directly: `module.nonexistent_name = 42` succeeds and
reads back 42, no error). So `pyADflow.setOption("transitionRowVolScale",
True)` appeared to work (no error, read-back matched) but the compiled
Fortran code always saw the **declared default** instead, regardless of
what Python set. Confirmed via a Fortran-side debug print during an actual
solve, before and after adding explicit `.pyf` declarations for the three
new variables.

**Consequence for this branch's own work**: `transitionNK`'s Fortran
default is `.true.`, so everything gated only by it (Eq. 59 NK
reactivation, Algorithm 2 in NK) was **genuinely active** the whole time —
those results stand. But `transitionRowVolScale`/
`transitionResidualAutoscale` default to `.false.`, so the "slight
differences" observed in the first pass of smoke testing were **pure
floating-point recompilation noise**, not real effects — the code paths
were never actually exercised until the `.pyf` fix.
`transitionSrcDtRestrict` and its siblings (`transitionSrcDtLimit`,
`srcDtDeactivateIters`, `transitionDampTheta`, `transitionDampMaxIter`)
were also missing from the `.pyf` at the time. **Update 2026-08-12: all
five are now wired (`adflow.pyf` L1129-1133), and `transitionSrcDtEigMode`
has been removed from the codebase entirely — neither is a live bug.** The
currently-known OPEN instance of this bug class is `module anksolver`:
`ANKPhysicalLSTolReTheta` and `omegaMinGamma` are silent no-ops (see
`docs/CORE_BASE/CORE_01_architecture.md`, Turb-ANK options).

**Results once genuinely wired** (`_old/nk_eq59_reactivation_test/run_10iter_step5_Sr_Sa_genuinely_wired_STALL.log`,
`run_10iter_step4_Sa_alone_genuinely_wired.log`, restart from
`r1c_CSANK_entry_rel1e-6_dp.cgns`, `nkswitchtol=1e-5`):

- **S_r geometric row scaling ON**: linear solve stalls almost completely
  (lin res pinned at 0.99–1.00 for 9 consecutive outer iterations,
  totalRes frozen ~21.75, ~300 wasted inner KSP iterations per outer step).
  No NaN, but clearly unusable as implemented. The `vol^(5/3)`/`vol^(2/3)`
  exponents (`setRVec`, `NKSolvers.F90`) were always a best-effort
  interpretation of the paper's SBP-metric-Jacobian `J` in terms of
  ADflow's physical `volRef` (cell volumes here span ~1e-12 to typical
  scale) — this result suggests that interpretation is wrong or at least
  badly conditioned for this mesh. **Stays off by default; not
  recommended for production runs without further investigation.**
- **S_a autoscale proxy ON, S_r OFF**: no stall — lin res varies normally
  (0.2–0.8), totalRes drops 21.75→10.9 over 8 iterations (real progress,
  roughly comparable order to the S_r/S_a-off baseline's 21.75→0.9 over
  the same 8, so **not clearly better**), but step sizes are noisier/more
  often small (0.34, 0.01, 0.01, 0.10) than baseline. Marginal;
  **stays off by default**, not recommended over the baseline without more
  testing.

## 7. Production run E1 fix: NK step pinned at minlambda (2026-07-16)

The first production run (`logs/3_NK_paper_faithful.log`, restart
`r2_after_first_NK_rel1.8e-8_dp.cgns`, `nkswitchtol=8e-8`) pinned `Step` at
exactly `0.01` — `LSCubic`'s hardcoded `minlambda` — for 1000+ consecutive
outer iterations, matching the debugging playbook's **E1** symptom exactly
(`ADFLOW_BASE/ADFLOW_04_debugging_playbook.md`: "very small step sizes...
Fix ladder: relax step control"). Root cause: `LSCubic`'s Armijo
sufficient-decrease coefficient `alpha` (hardcoded `1e-2`, never exposed as
an option) was essentially never satisfied at any λ above the floor, so
every iteration backtracked all the way down and "just took" the minimum
step per the fallback comment in the code.

**Fix, same spirit as CANK's already-relaxed `ANKUnsteadyLSTol`/
`ANKPhysicalLSTol`**: relaxed `LSCubic`'s `alpha` and the turb-residual
pre-limit factor (`turbRes2 > factor*turbRes1`).

- First attempt, `alpha: 1e-2→1e-4`, factor `2.0→5.0`: **SEGV crash** at
  iteration ~23 (`logs/3b_NK_paper_faithful_alpha_relaxed.log`) — but with
  real progress up to that point (totalRes 1.59→0.35 vs the old run's
  1.59→~1.5 over the same span). Isolation test with `transitionNK=False`
  (same alpha) did **not** reproduce the crash within the tested window,
  but also never moved (`Step` pinned at `0.00`) — inconclusive on
  attribution, but consistent with the crash needing the *combination* of
  a large accepted step and no physicality check to catch it. **NK has no
  physicality/bounds check at all** (unlike ANK's
  `physicalityCheckANK`/`physicalityCheckANKTurb`) — an over-permissive
  `alpha` can accept a step that drives density/energy unphysical with
  nothing to catch it, exactly the failure mode
  `ADFLOW_04_debugging_playbook.md` warns about for `ANKPhysicalLSTol`.
- Dialed back to `alpha: 1e-3`, factor `3.0`: **no crash**, verified past
  the previous crash point (100 iterations, `logs/3d_NK_paper_faithful_alpha_1e-3.log`),
  real progress (totalRes 1.59→~0.22, full/near-full steps appearing
  regularly, one self-correcting spike at iter 95 that settled back down).
  **This is the value in production use as of this session.**

**Open item this creates**: NK's missing physicality check (γ/Re̅θt bounds,
or basic ρ/E positivity) is the real gap — `alpha` tuning is a workaround,
not a fix. A faithful per-node Algorithm 2 equivalent for *density/energy*
positivity (not just γ/Re̅θt, which `applyNKAlgorithm2Damping` already
covers) would let `alpha` be relaxed further without the SEGV risk.

**Extended verification** (same session, `alpha=1e-3`/factor `3.0`):
- Deep restart (`r2`, `nkswitchtol=8e-8`,
  `logs/3d_NK_paper_faithful_alpha_1e-3.log`): 2500+ iterations, no crash,
  settles into a slow plateau around totalRes≈0.219 (one self-correcting
  spike at iter 95).
- **Shallow restart / early NK engagement** (`r1c_CSANK_entry_rel1e-6_dp.cgns`,
  `nkswitchtol=1e-5` — matching the paper's stated ~1e-5 engagement point —
  `_old/nk_eq59_reactivation_test/run_early_engage.log`): 35+ iterations, no
  crash, fast initial drop (totalRes 21.75→0.174) then a similar slow
  plateau. Confirms the `alpha` fix isn't restart-depth-specific.
  *(2026-08-12: "reachable without crashing" only — early engagement was
  later shown to stall production runs; see the §2b correction box.)*
- Both plateaus share the same character: a fast drop followed by a slow,
  non-crashing creep with `Step` oscillating between ~0 and ~1 rather than
  monotonic full steps. Worth investigating together with the missing
  physicality check (above) and the PC staleness item (§8.4) — plausibly
  the same root cause (PC/direction quality degrading as the front
  settles) rather than three separate issues.

## 8. Open code items, in impact order (2026-07-16 view; status updated 2026-08-12)

1. ~~Fix the `.pyf` bug for `transitionSrcDtRestrict` and siblings +
   `transitionSrcDtEigMode`~~ — **DONE** (all five wired in `adflow.pyf`
   L1129-1133; eigMode removed from the codebase). The open `.pyf` instance
   is now `anksolver` (`ANKPhysicalLSTolReTheta`/`omegaMinGamma` no-ops —
   see `CORE_BASE/CORE_01_architecture.md`).
2. NK has no physicality/bounds check (§7) — the real fix behind the
   `alpha` workaround. Would let `alpha` relax further without SEGV risk.
3. Investigate why S_r's geometric factor stalls the linear solve — likely
   the `volRef`-vs-paper's-`J` correspondence is wrong, or the exponents
   need recalibrating against how badly this mesh's cell volumes vary.
4. **Stronger NK PC for the front — now the decisive item.** Partly tested
   since (2026-08-07/12, AR5 family): `NKJacobianLag 5` — no effect (PC
   *quality*, not staleness, is the limit); EW-off + slack linear tol —
   falsified; ILU(3) — worse (stalls). No option-level lever remains; the
   deep-NK wall needs actual PC code work. See
   `docs/CORE_BASE/CORE_02_convergence_strategy.md` §"Options analysis".

~~Eq. 59 reactivation in NK~~ — done 2026-07-16.
~~Faithful Algorithm 2~~ — done 2026-07-16 (NK path; DADI already had it).
~~NK step pinned at minlambda (E1)~~ — done 2026-07-16, see §7.
