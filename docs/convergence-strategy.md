# SA-GR convergence strategy (validated 2026-07-15)

Current best practice for converging SA-Gamma-Retheta cases, distilled from
the 2026-07-14→16 campaign (3D plain wing, 175k cells, M=0.2, Tu=0.25%).

**Ready-to-use deliverable** (recipe, per-phase restart CGNS files, proof
logs, phase-entry runner `run_strategy.py`):
`~/Desktop/Run/MDO_PhD/Transition/gama_rethetha/Acelarating convergence/Baseline/3D_Plain_Wing/best_strategie/`

Full evidence: `.../3D_Plain_Wing/claude_attempt/` — master test table in
`TESTS_AND_CONCLUSIONS.md`, narrative in `long_overnight/DECISIONS.md`,
logs in `long_overnight/` (referred to as `RUN/` below). Solver background:
`ADFLOW_BASE/adflow_solvers.md`. Why ADflow differs from the paper's
solver: `adflow-vs-paper-solver.md`.

## The recipe (phase ladder)

| Phase | Activate at (rel totalRes) | Key options | Measured reliable range |
|---|---|---|---|
| ANK segregated | start (freestream) | `ANKUseTurbDADI: True` | -> rel ~1e-5 in ~28 min (flow converges; retheta res parks at ~2e4 — expected). **SANK variant** (`ANKSecondOrdSwitchTol 1e-4`, never couple, 2-leg run): **18 min, -36%** — needs its own leg because secondOrd is one value per run |
| CANK (coupled) | `ANKCoupledSwitchTol: 1e-5` | `ANKADPC: True`, LS below | full 1.00 steps to rel ~1e-7; kills retheta (1.8e4 -> 1.6e3 in 20 iters, ~1 min) |
| CSANK (2nd-order) | `ANKSecondOrdSwitchTol: 1e-6` | same LS | -> rel ~3.5e-8 (one order past CANK; ~40 min, iters get costly) |
| NK | `nkswitchtol` 4.2e-8 (just above CSANK's floor; 1e-6..1e-7 if skipping CSANK) | `NKADPC: True`, `NKSubspaceSize` 200-300 | engaging at CSANK's max: 3.13 -> 0.94 in ONE 3-eval full step; record rel 3.3e-9; **wall below ~5e-9** (lin res -> 0.8) |

Non-negotiable global options:
- `ANKUnsteadyLSTol: 1.5`, `ANKPhysicalLSTol: 0.5` — THE fix for coupled-phase
  step collapse (steps 0.01 -> 1.00). More aggressive (2.0/0.7) gives
  bit-identical results = no benefit; defaults (1.0/0.2) stagnate.
- `solutionPrecision: "double"` — single-precision restart files truncate the
  transition front and poison every restart (~1e-7 noise = the signal at deep
  residuals).
- `ANKADPC`/`NKADPC: True` — the FD-colored PC is unusable for SA-GR Newton
  phases (lin res 0.99); AD-assembled PC is cheaper and stronger.
- `turbResScale: [1e4, 0.1, 1e-4]`, first-order transition advection.

Fastest verified full path: freestream -> rel 1.8e-8 in ~32 min
(28 min segregated + ~1 min CANK + 2 NK iters), `RUN/cank_both.log` after
`RUN/seg_pure.log`. Old reference: 10 h to rel ~1e-5-equivalent (retheta 336).

## Engaging NK too early is the most expensive mistake (2026-08-07)

Measured on the AR5 corrected-foil family (5 levels, 0.46M-7.42M cells,
`11_ar5_corrected_foil/`), launched with the 09 campaign's inherited
`nkswitchtol = 1e-6`. **Every case where NK actually engaged stalled at rel
~1e-6**: `Step = 0.00`, CFL `----`, totalRes rising in the 12th digit, hundreds
of wasted iterations. That is L4 (six Tu-sweep runs plus the refinement leg),
L3, and L0. That threshold engages NK ~1.5 orders above the validated point
(CSANK's floor, rel ~3.5e-8 — see the ladder table).

Restarting the identical state with NK simply **disabled** (`useNKSolver: False`,
letting CSANK finish) closed the gap immediately. Note L2 and L1 were switched
PRE-EMPTIVELY, before NK engaged at all, so their rows show that CSANK alone
reaches 1e-8 — not that NK had stalled on them:

| level | cells | iters after restart | wall | from -> to |
|---|---:|---:|---:|---|
| L4 | 459,452 | 2 | 36 s | 6.81 -> 8.00 orders |
| L3 | 904,134 | 4 | 245 s | 6.89 -> 8.07 orders |
| L2 | 1,783,442 | 17 | - | -> 8.04 orders |
| L1 | 3,667,320 | 27 | - | -> 8.05 orders |

So: **if NK pins at `Step = 0.00`, the first thing to test is that it engaged
too early** — drive CSANK to its own floor instead, and only then hand over.
Do not restart NK at the same threshold.

Operational note: `scancel --signal=USR2` writes the state and ends the current
*solve*, not the job. A runner that issues several staged `CFDSolver(ap)` calls
will simply proceed to the next one; follow with a plain `scancel`.

## Known limits / do NOT

- Do NOT couple early: CANK at rel 1e-2 stagnates (with or without the LS
  relaxations) — the front is too unconverged and the global lambda throttles
  the whole field.
- Do NOT hold CANK past ~1e-7 or CSANK past ~3.5e-8: a front adjustment
  kicks the residual and the phase enters permanent oscillation (never
  recovers its best depth). Hand over before the kick.
- Do NOT restart deep states from single-precision files, and expect ANK CFL
  to re-ramp after any deep restart (floor = ANKCFLMin*(totalR0/totalR)^0.5;
  raise `ANKCFLMin` if stranded).
- `ANKNSubiterTurb` is a DEAD KNOB when `ANKUseTurbDADI=True` (only
  `nSubiterTurb` reaches DADI).

## Open problem

Deep endgame below rel ~1e-8: the NK linear solve saturates (lin res
0.8-0.99, GMRES 200-300 exhausted, 60-200 evals/iter) regardless of
engagement point, JacobianLag, or subspace size. Root causes and code items
(per-node Alg. 2 damping, source-dt reactivation inside NK, stronger PC) in
`adflow-vs-paper-solver.md` §5. Deepest state for PC experiments:
`best_strategie/restarts/r3_deepest_record_rel3.3e-9_dp.cgns`. Note the
wall's depth scales with how settled the field is (from the fully-settled
40k field, NK reached rel 6.4e-11 before the same wall — nk_colscale_test).

## vs the old 40k reference (same wing, old binary, pure segregated SANK)

`RUN/../Atent_1_Converded/out.txt`: 5284 iters / 10.06 h -> rho 1.06e-6,
retheta res 336. Milestones: retheta 1.0e4 @ 52 min, 4.6e3 @ 1.8 h, 1.9e3 @
3.7 h, 638 @ 7.7 h, 336 @ 10 h. The ladder crosses retheta 336 at ~25 min
total and reaches retheta ~0.1 (3000x deeper) by ~1.1 h. Note the reference
is deeper on rho (1.06e-6 vs ~5e-4 at that point) — segregated SANK grinds
the flow forever once turbulence settles; the ladder converges everything
together and passes that rho once NK runs to target.

## Options analysis — no-run exploration (2026-07-15, from logs + `ADFLOW_BASE/adflow_solvers.md`)

- `ANKCFLLimit 1e6`: not a lever. Segregated/CSANK ran AT the limit with
  full steps and lin res ~0.05; every step collapse in the logs coincides
  with a transition-front adjustment, not with CFL growth. Lowering it would
  only slow the good phases.
- `nSubiterTurb` (default 3): the solvers doc recommends 3-7 for RANS ANK.
  The only datapoint (seg_dadi15_sub8, =8, from a plateau) helped transition
  mildly. Plausible shortcut for the 28-min segregated leg; UNTESTED in the
  ladder — try 5-7 if the seg leg ever matters.
- `ANKLinearSolveTol 0.05`: healthy everywhere in the logs (0.04-0.07);
  cank_tight (0.005) already shown not to pay. Leave.
- `NKUseEW True` + weak deep PC = the 100-200-eval stall iterations:
  Eisenstat-Walker demands tighter linear tolerance as convergence improves,
  and with the PC saturating (lin res 0.8) GMRES just burns the full
  subspace. Analysis-based mitigation: `NKUseEW False` + `NKLinearSolveTol 0.3`
  caps the per-iteration waste (accepts poorer Newton steps, but each costs
  ~5x less). **TESTED 2026-08-07 — FALSIFIED.** On the AR5 corrected-foil L0
  (7.42M cells, 64 ranks, at rel 2.5e-8:
  `11_ar5_corrected_foil/results_refinement/L0`, job 1812053), with EW OFF and
  the linear tolerance pinned at a slack 0.3, GMRES still returns **0.998** —
  it cannot reach even that loose target. So EW is not what starves the deep
  NK iterations; the preconditioner alone is. This kills the last
  option-level lever: the deep wall is only addressable by the code items in
  `adflow-vs-paper-solver.md` §5. Same run confirmed the wall is not an
  engagement-point artefact either — NK entered at the correct rel 5e-8 (not
  1e-6) and merely re-attained the depth CSANK had already reached, with the
  linear residual DEGRADING 0.80 -> 0.97 as it went.
- `ANKPCUpdateTol 0.5` / NKJacobianLag: PC freshness. JacLag 5 tested — no
  effect on the deep wall; the PC's *quality*, not staleness, is the limit.
- Verdict: no remaining option is likely to move the deep-NK wall; the real
  fixes are the code items in `adflow-vs-paper-solver.md` §5.

## Index of everything tested (read the log only if you need the details)

All in `RUN/` (= `3D_Plain_Wing/claude_attempt/long_overnight/`) unless
noted; `RUN/DECISIONS.md` has the narrative and
`claude_attempt/TESTS_AND_CONCLUSIONS.md` the same table in Portuguese
with more detail.

| Test (log) | Question | Answer |
|---|---|---|
| archive/ab_matrix/* (7 variants + RESULTS.md) | which solver knobs matter from a plateau | ADPC >> FD PC; tight lin tol & loose LS alone: no effect |
| archive/ab_matrix seg_dadi5/15/sub8 | does more DADI fix segregated retheta stall | no; ANKNSubiterTurb dead knob; seg-from-plateau degrades |
| archive/long_overnight_logs/nk1em6*.log | CSANK@1e-5 bridge; deep restarts | CSANK thrash; restarts stranded at CFL floor + SP poison |
| nk_now.log | NK straight from deep state | works to rel 1e-6 then 200-eval wall |
| archive/long_overnight_logs/nk_pc1.log | bundle of stronger-PC NK knobs | toxic (first step blew transition 125x) — one knob at a time |
| full_dp.log | continuous run, couple at 1e-2, DP | CANK oscillates below rel ~1e-5; DP confirmed |
| seg1em4.log | couple at 1e-4 (no LS) | retheta burst great, rho stagnates (pre-LS) |
| seg_pure.log | pure segregated, never couple | flow converges, retheta parks at 2e4 — segregated can't finish |
| **cank_both.log** | LS 1.5/0.5 + coupled restart from seg state | **breakthrough: full steps, NK quadratic to rel 1.8e-8 in ~3 min** |
| endgame_1em6.log | NK grind below 1e-8 | wall: lin res 0.8, 0.5 order in 90 min (deepest 3.9e-9) |
| testA_cank1em8.log | hold CANK to 1e-8 | no — reliable to ~1.4e-7 then permanent oscillation |
| testB_early1em2.log | couple at 1e-2 WITH LS | still fails — 5 orders behind at equal wall time |
| deep_final.log / deep_final2.log | one-run ladder; NK JacLag5/GMRES300 | ladder works; deep-NK wall unchanged by those knobs |
| csank_1em6.log | CSANK@1e-6 as deep workhorse | best coupled depth: rel 3.5e-8, then front-kick fallback to CANK |
| csank_uls2.log | LS 2.0/0.7 vs 1.5/0.5 | bit-identical — LS is not the deep limiter; keep 1.5/0.5 |
| nk_after_csank.log | NK engaged at CSANK's max depth (rel 4.2e-8) | best deep engagement ever: 3.13 -> 0.94 in ONE 3-eval full step; record depth rel 3.3e-9; same PC wall below ~5e-9 |
| cank1em4_v2.log | couple at 1e-4 WITH LS (user request, 07-16) | loses decisively: at equal wall time (~38 min) it sits at rel 1.25e-5 / rho 0.73 vs the 1e-5 route's rel 7e-7 / rho 1.6e-3 — coupled global-lambda does the flow 10x slower than segregated; **1e-5 coupling confirmed optimal** |
| sank1em4.log | SANK (2nd-order segregated) from rel ~1e-4 | wins the segregated leg: rel 1e-5 in 18.0 min vs 28.3 (ANK); early ANK<->SANK flip-flop costs ~3 min |
| polish_nk.log | segregated polish of the record state + NK re-engage | restart transient (first SANK step) kicks totalRes 0.29 -> 438 — not worth it via restart (would likely work in-run) |
| sank_ladder.log | CANK/CSANK/NK from the SANK state | CANK full steps, no transient; partial validation (stopped at rel 1.8e-7 for reorganization, no anomalies) |
| archive/nk_colscale_test/, totalr0_check/, plateau_restart/ | NK column scaling; totalR0 sanity | scaling fixed NK lin solves; totalR0 consistent (8.91e7). Also: NK from the settled 40k field reached rel 6.4e-11 before the same wall — wall depth scales with how settled the underlying field is |
