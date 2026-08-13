# Design Decisions Log — SA-γ-Re̅θt

> **Not a defining/normative file.** This is a *memory* of what was already
> discussed and decided during the A1-A4 code audits — not a spec, not a
> source of truth, and nothing here overrides the paper (physics, CLAUDE.md
> rule 9) or the code itself. If this log and the current code/paper ever
> disagree, trust the code/paper and treat this file as stale — update it,
> don't defer to it. Its only job is to save a future reader from
> re-litigating a question that was already asked and answered, and to say
> *why* a piece of code looks the way it does.

Condensed record of audit questions raised and resolved during the A1-A4 code
audits (2026-07-06/07). This file is the durable takeaway; the raw
point-in-time audit trail lives in `../_archive/`. The old `findings/`
dir's four files (`A1`/`A2`/`A3`/`A_confirmacao`) were folded into this
single file 2026-07-09 to cut doc sprawl — nothing below lost information,
only the exploratory back-and-forth was trimmed. `D1_transitionRefLength` is
not repeated here — it's fully covered in `architecture.md` (Part 2,
`transitionRefLength` option) and `../ADFLOW_BASE/ADFLOW_08_nondimensionalization.md` §5 exception.

Status: A1-A3 fully closed. A4 (adjoint) closed with follow-up work tracked
in `TODO.md` and `current-task.md`.

---

## Nondimensionalization / paper-vs-code (A1)

### D2 — LM2009 safeguards present in code but not written in the paper (kept, necessary)

- **λ_θ clamp** to [−0.1, +0.1] (`saGammaRetheta.F90:560-561`, smooth
  min/max) and **Tu floor** at 0.027% (`turbUtils.F90:2297`) come from the
  original Langtry–Menter (2009) correlations, not from P&Z (2020) text.
  Third LM2009 limit, Re̅θt ≥ 20, *is* in the paper (Algorithm 2) and matches
  `rsaGRreThetaLo = 20`.
- **Not cosmetic — required for robustness.** Without the λ_θ clamp,
  `exp(−35·λ_θ)` in F1 overflows for λ_θ < −20.3 (Inf → NaN near
  stagnation); with the clamp, F(λ_θ) ∈ [0.532, 1.268] always finite/positive.
- Jacobian (`evalSrcJacBlock`) treats the target as constant, so the clamps
  don't enter it — no residual/Jacobian conflict.
- Verdict: keep both clamps as-is.
- **Update 2026-08-03:** clamp made compile-time switchable
  (`rsaGRclampLambdaTheta`, `paramTurb.F90`, a `parameter` — compile-time,
  not runtime), unclamped tested on S809 and REJECTED (stalls ~1e-4, +38
  counts); default remains ON (`be9d6d1d` → `efed31cf`). The ReThetaT
  correlation also gained a ≥20 floor.

### D3 — Retired distilled-physics doc (historical note)

- A now-deleted "distilled" physics summary had 5 transcription errors found
  across the A1-A2 audits (F_θt branch, F_wake formula, missing √γ in P_γ,
  garbled F1-F3 blend, and Algorithm 2 misread as a single scalar damping
  factor instead of two independent per-variable back-off loops). The code
  was correct in every case; only the distilled doc was wrong.
- Decision (2026-07-07): the distilled doc was deleted. The full paper text
  (`SA_GAMMA_RETHETHA_BASE/Piotrowski_Zingg_2020_...md`) is the sole physics
  reference (CLAUDE.md rule 9).

---

## Convergence strategy vs paper §IV (A2)

### D-A2-1 — No Eq. 59 dt-restriction in coupled ANK/NK paths (accepted, not addressed)

- `grep -n srcLambda src/NKSolver/NKSolvers.F90` → only appears in the
  Turb-ANK KSP path. Fully-coupled ANK/NK rely on the physicality check +
  CFL ramp alone. User decision: leave as-is; not expected to be exercised
  soon. CANK warnings added as a guard instead.
- *Update 2026-08-12: Eq. 59 now also exists in the NK path via the
  `transitionNK` bundle (`applyNKSrcDtDiagonal` etc.).*

### D-A2-2 — Additive (DADI) vs MAX (turbKSP) form of Eq. 59 — both kept, intentional

- DADI: `qq(m,m) += srcLambda/0.9` (additive). turbKSP:
  `max(dtinv, srcLambda/0.9)` (matches paper Eq. 59 literally).
- DD-ADI has no `dtl`/CFL term to MAX against (pseudo-transient is
  `cc/cfl`, not a real Δt) — the additive form is the natural point-implicit
  realization there, same spirit as `qq = max(qq, zero)` in plain SA.
- Both forms only touch the LHS diagonal → converged solution is identical
  either way; worst case the additive form is ~2-2.5× more conservative
  locally (with `alfaTurb=0.8` relaxation stacking an extra 1/α factor).
- Tuning knobs if DADI stalls at the transition front: raise
  `transitionSrcDtLimit`, or move the `+= srcLambda/0.9` after the
  `factor = 1/α` relaxation multiply (1-line change, not done).

### D-A2-3 — Eq. 59 deactivation switch (§IV.B.3): DADI correct, turbKSP fixed

- Paper deactivates the restriction only during the inexact-Newton phase (5
  clean iterations, no backtracking, R_d > 1e-5); reactivates on backtrack or
  residual rise.
- **DADI**: never deactivates — correct, since DADI *is* the
  approximate-Newton globalization phase the paper says stays restricted.
  Runtime-tunable without recompiling (`transitionSrcDtRestrict`,
  `transitionSrcDtLimit`). Open question tracked in `TODO.md`: whether the
  DADI tail slows down visibly with the restriction always on.
- **turbKSP**: had 3 bugs vs the paper, all fixed 2026-07-07 in
  `NKSolvers.F90` (`ANKTurbSolveKSP`, ~lines 3505-3792):
  1. Clean-iteration counter ran from solver start instead of only inside
     the second-order regime (`totalR ≤ ANKSecondOrdSwitchTol·totalR0`).
  2. No reactivation on residual rise — fixed, counter resets when
     `totalR > ANKSecondOrdSwitchTol·totalR0`.
  3. A successful backtrack was still counted as a "clean" iteration — fixed
     with a `backtrackTriggered` flag armed at backtrack entry regardless of
     outcome.
- Operational consequence: default `ANKSecondOrdSwitchTol = 1e-16` means the
  second-order regime is never entered, so turbKSP's restriction never
  deactivates either (conservative). Set to `~1e-5` (paper's phase-switch
  value) to get the paper's acceleration.

### D-A2-4 — `turbResScale` ν̃ scale: 1e4 (ADflow) vs 1e3 (paper §IV.C) — keep ADflow convention

- **Correction 2026-08-12:** the shipped default is `[1e4, 0.1, 1e-4]`
  (~1/state magnitude, see the `pyADflow.py` comment), NOT the
  `[1e4, 10, 1e4]` discussed below — the argument below is historical; the
  reconcile happened in code (the 2026-07-16 campaign value became the
  default).
- Default `[1e4, 10, 1e4]` (γ, Re̅θt scales match the paper exactly; only ν̃
  differs). `turbResScale` is residual row-scaling — a conditioning/tuning
  knob, not part of the converged solution, so "paper wins" (rule 9) doesn't
  apply to it.
- ν̃'s SA-GR residual equation is *exactly* ADflow's plain-SA residual (γ
  only multiplies production) → it should carry SA's historical calibration
  (1e4), not the paper's velocity-based-nondim calibration (1e3). γ and
  Re̅θt are nondim-independent, which is why only ν̃ diverges from the paper.
- Rejected alternative: scaling γ/Re̅θt by ×10 to preserve the paper's
  *ratios* instead — the ×10 gap is specific to ν̃'s nondim, not a
  transferable global ratio. Would trade two paper-matching, already-tested
  values for guesses.
- Future correct calibration (if a coupled run stalls on a turb residual
  line): *measure*, don't transfer ratios — print scaled per-equation
  residual norms on a representative case and tune via runtime
  `turbresscale`. Tracked in `TODO.md`.

### D-A2-5 — Per-variable damping (Algorithm 2), not a single scalar factor — code matches paper

- DADI damps γ and Re̅θt **independently** (two separate back-off `while`
  loops, `saGammaRetheta.F90:2109-2135`, counter `m` reset between them,
  each with its own θ^m); ν̃ is never damped (only `max(w,0)`). This is
  exactly Algorithm 2's real pseudocode (paper lines 622-637) — the
  single-scalar-factor reading was the retired distilled doc's error (D3).
- Fix applied 2026-07-07: the paper's back-off loop is unbounded and
  rejects hard clipping (paper line 620) — the code capped it at 40
  iterations then hard-clipped unconditionally, degenerating into the
  clipper the paper explicitly rejects for large overshoots.
  `transitionDampMaxIter` default raised 40 → 10000 (0.99¹⁰⁰⁰⁰ ≈ 0,
  effectively unbounded); the hard clip is now a last-resort fallback only
  reachable if the loop exhausts, and fires a warning with cell counts.
- Scope: this back-off+clip mechanism is DADI-only. ANK/turbKSP paths
  achieve the same bounds via `physicalityCheckANK{,Turb}` (clamping the
  global Newton step λ), which was already correct — no change needed there.

### N-A2-6 — Blockettes disabled for SA-GR (accepted, deferred)

- `useBlockettes=False` is forced for this model
  (`blocketteResCore` doesn't implement it) — all matrix-free ANK/NK
  mat-vecs use the slower residual path. Performance cost, not correctness.
  User decision: leave for now; implement in `blocketteResCore` when
  performance is evaluated. Tracked in `TODO.md`.
- *Update 2026-08-12: blockette kernels ARE implemented and tested since
  2026-07-24; only the pyADflow force-off remains (`pyADflow.py:~6824`).*

---

## Code-level coherence vs SA / SST (A3)

Scope: SA-GR vs `sa.F90`/`SST.F90` and the multi-equation infrastructure,
not vs the paper (that's A1/A2). Everything checked and found consistent
(dispatch in `turbAPI.F90`, module structure mirroring `sa.F90`, SA source
term reproduced verbatim with γ multiplying only production terms per Eq.
41, ν̃ diffusion identical to `saViscous`, generic advection/unsteady
routines, `tdia3x3` solver, `ResScale`, implicit relaxation, eddy viscosity
dispatch to `saEddyViscosity`, 8-var state vector, wall BC via the standard
`bmt/bvt` machinery, `turbResScale` auto-set) is **not** repeated here since
nothing needed changing. Only the divergences follow.

### W1 — `useft2SA` defaults True in a model named "noft2" (guard added)

- The paper's SA-noft2 base implies `ft2 = 0`, but nothing forces
  `useft2SA=False` for the `SA-noft2-Gamma-Retheta` enum — with all-default
  options the model runs *with* ft2 active, contrary to its name.
- Resolution: a warning was added rather than silently forcing the option
  (keeps user scripts working, flags the mismatch).

### W2 — Wall functions unsupported, no guard (error added)

- `saGammaReThetaSolve` has no wall-function block (unlike `saSolve`) and no
  curve-fit data is initialized for this enum. Physically correct
  (transition needs a resolved BL, y+≈1) but previously ran silently with
  `wallFunctions=True` + SA-GR.
- Resolution: an error is now raised if both options are selected together.

### W3 — Source-Jacobian diagonal clip (resolved 2026-07-07)

- `saSource` clips its diagonal (`qq = max(qq, zero)`) as standard ADflow
  DDADI philosophy — SA-GR's `qq(1,1)` (ν̃) and `qq(2,2)` (γ) lacked this
  clip (only `qq(3,3)` had it, and that one never fires). `qq(2,2)` goes
  negative *routinely* for small γ in the pre-transition zone.
- This is LHS-only stabilization (not a paper equation, rule 9 doesn't
  apply) — converged solution unaffected, only DDADI robustness.
- Fix: `qq(1,1) = max(qq(1,1), zero)` and `qq(2,2) = max(qq(2,2), zero)`
  added, mirroring SA and the pre-existing `qq(3,3)` clip. Both clips are
  inside `#ifndef USE_TAPENADE` — no AD regen needed.

### W4 — Farfield inflow BC aligned to ADflow's standard pattern (resolved 2026-07-07)

- SA-GR had a special-cased farfield-inflow BC (`bvt = 2·wInf`, `bmt = +1`,
  face-value form) instead of the generic ghost-value pattern
  (`bvt = wInf`, `bmt = 0`) used by SA/SST. Neither the paper nor LM2015
  mandates the face-value form (paper doesn't specify discrete BCs at all).
- Fix: removed the SA-GR special case in `bcTurbFarfield`; it now falls
  through to the generic ghost = `wInf` branch, identical to SA/SST. γ_∞=1
  and Re̅θt,∞ (from the Tu∞ correlation) are still imposed correctly.

### W5 — Different initialization defaults (deliberate, kept)

- Interior field γ = 0.02 vs γ_∞ = 1 (suppresses SA production until Fonset
  activates γ — `initializeFlow.F90:2229-2245`). `eddyVisInfRatio` default
  1e-10 for SA-GR vs 0.009 (SA) — freestream quasi-laminar, consistent with
  a transition model. Physics validation of these values is on the user;
  tracked in `TODO.md` ("gamma/Re̅θt initialization value").

### W6 — Stale doc claims (fixed where found)

- `architecture.md` previously said γ's wall BC was Dirichlet (γ=0); it's
  actually zero-gradient/Neumann (`bmt=-1`) — corrected in `architecture.md`.
  A stale file path/line-count reference to a nonexistent
  `saGammaRethetaHelpers.F90` was also corrected there (helpers actually
  live in `turbUtils.F90:2279-2410`).

---

## Provenance

Original audit trails (now retired, folded in here 2026-07-09):
`findings/A1` (nondim, folded into A1 section above + `architecture.md` +
`../ADFLOW_BASE/ADFLOW_08_nondimensionalization.md`), `findings/A2_convergencia.md`,
`findings/A3_coerencia.md`, `findings/A_confirmacao.md`,
`findings/D1_transitionRefLength.md` (folded into `architecture.md` Part 2
+ `../ADFLOW_BASE/ADFLOW_08_nondimensionalization.md` §5). A4 (adjoint) findings live alongside
this file in `adjoint_audit_2026-07-07.md` and `sst_dev_lessons.md`.
