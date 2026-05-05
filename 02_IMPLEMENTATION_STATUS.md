# Implementation Status

> Single source of truth. Read this BEFORE starting any task. Update the
> matching row when a task lands. The ordered task list at the top is the
> roadmap — work it top-to-bottom, one task per Claude Code session.

---

## Ordered Task List (Roadmap)

Do these in order. Each row links to a section in
`03_IMPLEMENTATION_PLAN.md`.

| #  | T-ID    | Title                                            | AD?  | Status |
|----|---------|--------------------------------------------------|------|--------|
| 1  | A1      | Smoke baseline — SA only, transition off         | no   | ✅      |
| 2  | A2      | Wire `transitionDebug` into volume CGNS          | no   | ✅ (48 slots, vol+surf)      |
| 3  | A3      | Smoke run — transition on, γ=1 forced            | no   | ✅      |
| 4  | B1      | Verify timeScale matches ADflow nondim convention| no   | ✅      |
| 5  | B2      | Verify φ_p overflow safety (Algorithm 1)         | yes  | ✅      |
| 6  | B3      | Populate off-diagonal source Jacobian            | yes  | ✅      |
| 7  | B4      | Initialize row/column scaling factors            | no   | ✅      |
| 8  | C1      | First-order upwind option for γ, Re̅θt           | yes  | ✅      |
| 9  | C2      | Source-term dt restriction — DADI (all 3 modes)  | no   | ✅      |
| 10 | C3      | Source-term dt restriction — Turb-ANK CFL cap    | no   | ❌      |
| 11 | C4      | 5-iter deactivation switch                       | no   | ❌      |
| 12 | D1a     | Smoke — DADI decoupled (diagonal only)           | no   | ❌      |
| 13 | D1b     | Smoke — DADI transition (SA decoupled, γ-Re̅θt coupled) | no | ❌  |
| 14 | D2      | Smoke — DADI full 3×3 coupled                    | no   | ❌      |
| 15 | D3      | Smoke — coupled ANK (flow+turb in one Newton)    | no   | ❌      |
| 16 | D4      | Smoke — Turb-ANK KSP (decoupled)                 | no   | ❌      |
| 17 | E1      | Tapenade clean regen — full sweep                | yes  | 🟠      |
| 18 | E2      | Build with AD enabled                            | no   | ❌      |

Legend: ✅ done · 🟡 partial · 🟠 suspect/needs verification · ❌ not started

After E2 lands, the model is implementation-complete. User runs validation
cases (NLF0416, S809) on their own.

---

## Known Bugs / Open Questions

| #   | Where                              | Status                              | Resolved by |
|-----|------------------------------------|-------------------------------------|-------------|
| #1  | `saGammaRetheta.F90:510-559`       | ✅ Fixed — off-diag Jacobian populated (B3) | B3 |
| #2? | `saGammaRetheta.F90:430-431`       | Confrimed — `timeScale` do not need explicit Re factor due to ADflow nondim convention. 
| #3  | feedback loop                      | Symptom (γ pinned at 1) downstream of #1 and possibly #2; should disappear after B1+B3+E1 | (auto) |

---

## Per-Phase Snapshot

### Phase 1 — Plumbing
- T1.1 State-vector indexing (itu1,itu2,itu3) — ✅
- T1.2 Input options — 🟡 (`turbResScale` uninitialized; not exposed in Python)
- T1.3 Halo exchange — ✅
- T1.4 Boundary conditions — ✅
- T1.5 First-order upwind option — ❌ (→ task C1)
- T1.6 Multigrid restriction — SKIP

### Phase 2 — LM2015 Sources
- T2.1 φ_p smooth max/min — ✅ (Algorithm 1 four-branch overflow-safe form)
- T2.2 Re_θt correlation — ✅
- T2.3 F_onset, F_length, Re_θc (LM2015) — 🟡
- T2.4 P_θt source term — 🟠 (timeScale convention unverified, → task B1)
- T2.5 Crossflow D_scf — SKIP (helper exists, user defers)
- T2.6 P_γ, E_γ — ✅
- T2.7 SA coupling (γ × P_SA) — ✅

### Phase 3 — sLM2015 Smooth Variants
- T3.1 smooth F_onset, F_turb (with `useLM2015Fturb` switch) — ✅
- T3.2 smooth F_length, Re_θc — ✅
- T3.3 smooth F(λ_θ) — ✅
- T3.4 Vorticity-limited P_γ, E_γ — ✅
- T3.5 Variant dispatch — ✅

### Phase 4 — Numerical Robustness
- T4.1 Row/column scaling — ✅ (Fortran default=1.0; Python sets model-specific values)
- T4.2 Solution-update damping (Alg. 2) — ✅
- T4.3 3×3 source Jacobian — ✅ (diagonal + off-diagonal)
- T4.4 Source-term dt restriction — 🟡 (DADI diag only; → tasks C2, C3; all 3 TurbDADICoupled modes)
- T4.5 Deactivation switch — ❌ (→ task C4)

### Phase 5 — Coupling (see `07_COUPLING_MODES.md`)
- T5.1a DADI decoupled (diagonal) — ✅ runs (→ smoke D1a)
- T5.1b DADI transition (SA decoupled, γ-Re̅θt coupled) — ✅ implemented (→ smoke D1b)
- T5.1c DADI full 3×3 — ✅ implemented (→ smoke D2)
- T5.2 Turb-ANK KSP — ✅ infrastructure exists (→ smoke D4)
- T5.3 Fully coupled ANK — 🟡 (ANK via Jv works; → smoke D3)
- T5.4 Benchmark all modes — user does at end

### Phase 6 — Adjoint
- T6.1 AD seed zeroing — ✅
- T6.2 Tapenade regen — 🟠 (→ task E1)
- T6.3 Manual diff updates — DEFER
- T6.4 Dot-product verification — DEFER

### Phase 7 — Validation
User runs at end. Not Claude Code's job.

### Phase 8 — Volume Debug Output
Folded into task A2 (no longer separate).

---

## File Inventory (Key Files)

| File                                            | Lines | Status                              |
|-------------------------------------------------|-------|-------------------------------------|
| `src/turbulence/saGammaRetheta.F90`             | 1862  | 🟡 sources live, DADI diag only     |
| `src/turbulence/saGammaRethetaHelpers.F90`      | 367   | ✅ all smooth functions done        |
| `src/modules/constants.F90:128`                 | —     | ✅ enum=8                            |
| `src/modules/paramTurb.F90:32-52`               | —     | ✅ constants match paper             |
| `src/modules/inputParam.F90:293,298`            | —     | 🟠 `turbResScale` uninit             |
| `src/modules/block.F90:662`                     | —     | ✅ `transitionDebug` allocated + filled (48 slots) |
| `src/initFlow/initializeFlow.F90:140,2237`      | —     | ✅ γ_init = 0.02                     |
| `src/turbulence/turbBCRoutines.F90:441,921`     | —     | ✅ wall γ=0, farfield γ=1            |
| `src/NKSolver/NKSolvers.F90:3191,3359`          | —     | ✅ variable bounds                   |
| `src/adjoint/output*/saGammaRetheta_*.f90`      | —     | 🟠 untested after recent changes     |
