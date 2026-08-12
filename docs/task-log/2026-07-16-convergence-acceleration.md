# Convergence acceleration campaign (SA-GR solver strategy) — 2026-07-16

> Path note (2026-08-12): the run tree was renumbered — `best_strategie/` is
> now `03_convergence_strategy/3d_plain_wing/best_strategy/` and
> `claude_attempt/` is `.../3d_plain_wing/_old/campaign_2026-07-14_to_16/`.
> Also: the "regardless of engagement point" clause and the EW mitigation in
> the follow-ups were later falsified/refined — see
> `docs/convergence-strategy.md` (2026-08-07/08, commits ef9fc10d/54c43475).

**Problem:** SA-GR cases only converged via ~40k segregated iterations
(~10 h to rho 1e-6 / retheta res 336, old binary); coupled and Newton
phases stalled, oscillated, or thrashed, and nobody knew which solver
phase/tolerance combination was viable.

**Solution:** 3-day empirical campaign (23 tests, 2026-07-14→16) on the
3D plain wing (175k cells) that produced a validated phase ladder and
isolated three root causes:
- `ANKUnsteadyLSTol 1.5` + `ANKPhysicalLSTol 0.5` cure the coupled-phase
  global-lambda step collapse (the single biggest win).
- `solutionPrecision "double"` is mandatory for restarts (single-precision
  writes truncate the transition front and poison every deep restart).
- AD-assembled PCs (`ANKADPC`/`NKADPC`) are required (FD-colored PC gives
  lin res 0.99 on SA-GR Newton phases).
Final ladder: segregated ANK (28 min; SANK variant 18 min) -> CANK@1e-5
(~1 min, kills retheta 2e4 -> 1.5e3) -> CSANK@1e-6 (~37 min, to rel 3.5e-8)
-> NK@4.2e-8 -> **record rel 3.3e-9 in ~1h20 productive wall** (vs 10 h
reference, 3000x deeper on retheta). No model/Fortran changes — options
and scheduling only (the pre-existing branch changes: NK column scaling,
Eq. 59 srcDt in the coupled path, wall-BC DADI fix, were prerequisites).

**Files created/touched and why:**
- `docs/convergence-strategy.md` — the distilled recipe + measured phase
  limits + no-run options analysis + indexed table of all 23 tests with
  proof-log locations (so nothing is retested from scratch).
- `docs/adflow-vs-paper-solver.md` — why ADflow's solver hierarchy differs
  from Piotrowski & Zingg §IV (per-node Alg. 2 damping, source-dt
  reactivation inside Newton, Eq. 58 scaling) and the code items that
  remain to close the deep-NK gap.
- `docs/ADFLOW_BASE/adflow_solvers.md` — official ADflow solvers doc
  imported into the KB (user-added), indexed.
- Run-side deliverables (case dir `.../3D_Plain_Wing/`): `best_strategie/`
  (STRATEGY.md, per-phase restart CGNS files, proof logs, `run_strategy.py`
  with per-phase entry points) and `claude_attempt/` (full campaign:
  TESTS_AND_CONCLUSIONS.md master table, DECISIONS.md narrative, all logs).

**Verification:** every claim in the strategy is backed by a named log in
`best_strategie/logs/` or `claude_attempt/`; the record run
(`2_CANK_CSANK_NK_to_record.log`) reaches totalRes 0.294 (rel 3.3e-9) from
freestream with the documented switches. No physics validation performed
(per CLAUDE.md rule 7 — user verifies physics).

**Follow-ups (in `../TODO.md` and `adflow-vs-paper-solver.md` §5):**
- Deep-NK linear-solve wall below rel ~5e-9 (lin res 0.8 regardless of
  engagement/JacLag/subspace) — needs code: faithful Alg. 2 per-node
  damping, Eq. 59 reactivation inside NK, stronger PC. Case dir prepared:
  `.../3D_Plain_Wing/ADFLOW_SA_GAMMA_RETHETHA_SOLVER/` (branch
  `sa_gamma_rethetha_paper_solver`).
- `ANKNSubiterTurb` dead-knob with `ANKUseTurbDADI` — wire or document.
- turbResScale used throughout: `[1e4, 0.1, 1e-4]` (differs from the
  default kept in TODO §turbResScale — reconcile).
