# ADFLOW 04 — Debugging Playbook

**Use for:** a run stalls, diverges, or gives bad gradients. Each entry: symptom → likely cause → fix ladder → evidence. Fixes reference options in `05`; the mechanism is in `01`/`02`; the flag↔math map is `03`.

**Prime directive:** the defaults are tuned for **transonic 3-D aircraft, millions–tens of millions of cells**. Before deep debugging, confirm your case is in that regime; if not, expect to re-tune (see "Regime rules" at end).

---

## Fast triage

| Symptom (from solver output) | First move | Detail |
|---|---|---|
| `Step` repeatedly < 0.1 | If coupled → lower `ANKCoupledSwitchTol`; else lower `ANKCFLLimit` or enable 2nd-order before the small steps | E1 |
| `Lin Res` repeatedly > 0.5 | Strengthen PC (`ILUFill`↑, `ASMOverlap`↑) / let auto-CFL cutback act; **do not** enable 2nd-order | E2 |
| Many nonlinear iters, full steps, `Lin Res` fine, but slow nonlinear conv. | Enable 2nd-order and/or coupled at the stall convergence level | E3 |
| Turbulence residual (`resturb`) stalls high | Switch turb solver (turbDADI↔turbKSP) / more `nSubiterTurb` / decouple longer | E4 |
| NK fails/stalls right after switch | Lower `NKSwitchTol` (1e-5→1e-6…1e-8) — but see E5's SA-GR rule first | E5 |
| NK `Step = 0.00`, `CFL ----` (SA-GR) | Premature engagement — disable NK / let CSANK finish, then engage ≈4.2e-8 if at all | E5 |
| Multigrid stalls | Use ANK startup | E6 |
| Overset gradients wrong/inconsistent | `useMatrixFreedrdw=True` | E7 |
| Adjoint won't converge | Strengthen nested PC / AMG global | E8 |
| Gradients ~0.1% off but run completes | You're on FD-Jacobian partials → use AD/Jacobian-free | E9 |
| Coupled ANK: tiny steps near separation | Stay decoupled | E10 |

---

## Detailed entries

### E1 — Very small step sizes (`Step` < 0.1 repeatedly)
- **Cause:** approximate Jacobian (`R1`/`R2`) does not guarantee a descent direction, so the line search over-shrinks; or CFL is at its ceiling and the time term can't stabilize; or moving shocks; or coupled turbulence engaged too early.
- **Fix ladder:** (1) if coupled → lower `ANKCoupledSwitchTol` (converge tighter before coupling). (2) If decoupled and CFL at limit → lower `ANKCFLLimit` from `1e5` (not below a few hundred). (3) Enable exact/2nd-order **just before** the small steps start via `ANKSecondOrdSwitchTol` (better update direction). (4) Relax step control: `ANKUnsteadyLSTol` 1.0→1.5 (allow unsteady-residual rise over a "hill"; risks divergence) or `ANKPhysicalLSTol` 0.2→0.4–0.6 (more aggressive ρ/E; never >1 → risks negative ρ/E). For SA-GR the FD-colored PC is unusable in Newton phases (lin res ~0.99): `ANKADPC: True` / `NKADPC: True` is mandatory and cheaper (convergence-strategy.md).
- **Evidence:** Yildirim §3.6 (descent not guaranteed with `R1`/`R2`; on failure the solver lowers CFL). Strut-braced-wing / nacelle separation cases expected to need more, smaller-step iterations (Yildirim §3.9, §5.6).

### E2 — High linear residual (`Lin Res` > 0.5; >0.9 ≈ stalled)
- **Cause:** linear system too stiff for the preconditioner — weak PC and/or an over-stiff operator (e.g. exact/2nd-order Jacobian, or `R2` on a skewed mesh).
- **Fix ladder:** (1) let the automatic mitigation act — if `Lin Res` exceeds `ANKLinResMax`, ADflow cuts CFL until it drops. (2) Strengthen PC: `ILUFill` 2→3, `ASMOverlap` 1→2–3, more GMRES iters. (3) **Do not** switch to 2nd-order (it makes the linear system *harder*); if you were on 2nd-order, disable it. (4) On **skewed multiblock**, prefer `R1` (default) over `R2`/`R0`. (5) For SA-GR the FD-colored PC is unusable in Newton phases (lin res ~0.99): `ANKADPC: True` / `NKADPC: True` is mandatory and cheaper (convergence-strategy.md).
- **Evidence:** Yildirim Table 6 (CRM WB multiblock): `R2` **fails** (orthogonal-mesh assumption breaks on skew) and `R0` **fails** (added terms increase stiffness → linear solves fail); only `R1` converges. Overset CRM WB (Table 5) preserves orthogonality → `R2` is fine there. ⇒ mesh quality, not a universal ranking, decides R-level.

### E3 — Large nonlinear-iteration count (steps full, linear fine, conv. slow)
- **Cause:** the approximate operator's inherently lower nonlinear rate; you're paying the "more iterations, cheaper each" side of the trade without the payoff.
- **Fix ladder:** record the rel-conv where it stalls; enable exact/2nd-order (`ANKSecondOrdSwitchTol`) and/or coupled (`ANKCoupledSwitchTol`) at a value slightly *above* that stall level.
- **Evidence:** Yildirim central trade-off (Tables 3, 5): approximate `R1`/`R2` need more nonlinear iterations than `R0` but less wall time — until the rate drops too far, when 2nd-order/coupled recovers it. Expected budget: ~100 nonlinear iters / ~2–3k linear for wing-only, ~200 / up to 5k for complex configs, to reach 4–5 orders.

### E4 — Turbulence residual won't converge
- **Cause:** SA residual behaves wildly early ("hill"); wrong turbulence solver for the case size; too few sub-iterations; scaling off.
- **Fix ladder:** (1) switch `ANKUseTurbDADI` (turbDADI ⇄ turbKSP). (2) Increase `nSubiterTurb`/`ANKNSubiterTurb`. (3) Stay **decoupled** until the turb residual clears the hill (rises to ~4 orders, then falls). (4) With turbKSP, set `ANKTurbKSPDebug=True` and read `LIN RES/ITER/REASON/STEP`. (5) `turbResScale`: for SA-GR the 3-element `[1e4, 0.1, 1e-4]` scaling is mandatory (auto-set default; never force a scalar). If the printed turb residual isn't ~4 orders below the flow residuals, that gap (not the flag) is the real symptom.
- **Evidence:** Yildirim §3.8 + Figs 4–6 (turbulence PC goes stale fastest mid-startup; decoupling recovers flow convergence and improves conditioning). Docs: turbDADI for <1M multiblock, turbKSP for >1M / overset.

### E5 — NK fails or stalls after the switch
- **Cause:** switched to NK while the state is still outside the basin of attraction; EW then picks a loose tolerance (state far away).
- **Fix ladder:** lower `NKSwitchTol` (1e-5 → 1e-6, hard cases 1e-8) so ANK does more of the transient. If EW's chosen `Lin Res` sits at the hard cap 0.8 for several iterations → definitively too early → lower `NKSwitchTol`. If NK still won't take, ANK alone can reach 1e-12 on hard cases (see below).
- **SA-GR measured rule (2026-08-12, ef9fc10d / 54c43475):** an NK stall at `Step = 0.00` means PREMATURE engagement — do not re-enter NK at the same threshold; disable NK / let CSANK finish, then engage just above CSANK's floor (~4.2e-8) if at all.
- **Evidence:** Yildirim §5.3 (upwind Roe case): **NK fails** after the switch, but **ANK-only reaches 1e-12** (Table 4) — the approximations make the linear systems solvable. Multiblock CRM WB with `R1`: switching at 1e-5 fails in NK, switching at **1e-6** converges (Table 6, Fig 9). ⇒ the NK switch criterion is genuinely case-dependent; delay it when unsure.

### E6 — Multigrid stalls
- **Cause:** separated flow (MG smoother stalls) or overset mesh (no coarse levels available).
- **Fix:** use ANK startup (`useANKSolver=True`, set `ANKSwitchTol`); it converges both multiblock and overset and handles separation.
- **Evidence:** Docs (Solvers). Yildirim §5.7: on a 1172-airfoil sweep, ANK matches multigrid startup performance while needing **only the finest level** — the reason overset (uncoarsenable) cases rely on ANK. Kenway §3.2: ANK is the default for overset in ADflow.

### E7 — Overset adjoint gradients wrong / inconsistent
- **Cause:** matrix-based `drdw` on overset ignores/incorrectly handles overset connectivity.
- **Fix:** `useMatrixFreedrdw=True` (mandatory on overset). Also verify `frozenTurbulence=False` if accuracy matters.
- **Evidence:** Docs (known bug #204). Kenway §6/§4.6: Jacobian-free (reverse-AD) mat-vec needs no sparsity/connectivity structure → correct and easy on overset; this is a headline advantage of the Jacobian-free adjoint.

### E8 — Adjoint linear system won't converge
- **Cause:** transposed Jacobian from 3-D viscous turbulent flow is ill-conditioned; PC too weak.
- **Fix ladder:** (1) `ILUFill` 2→3, `ASMOverlap` 1→2–3, `outerPreconIts`↑. (2) For large meshes, `globalPreconditioner=multigrid` (AMG). (3) Confirm turbulence is included (`frozenTurbulence=False`) unless deliberately trading accuracy for convergence.
- **Evidence:** Kenway §4.4 (strong PC essential; Table on eigenvalue clustering rationale), §A.3 (PC construction), Fig 3 (nested ASM+ILU+Richardson). Right-preconditioning `(27)`.

### E9 — Gradients ~0.1% off but the run completes
- **Cause:** using the **FD-Jacobian** adjoint (finite-difference partials) → subtractive-cancellation error.
- **Fix:** switch to **Jacobian-free** (or AD-Jacobian) for machine-precision partials; use complex-step only as a reference check.
- **Evidence:** Kenway Tables 3–4 (Jacobian-free/AD match complex-step to 11–12 digits; FD-Jacobian matches to 3–4 digits, ≈0.1% avg error), Table 10. Matters for quasi-Newton (BFGS) optimizers and tightly-converged constraints (Peter & Dwight).

### E10 — Coupled ANK produces tiny steps near separation
- **Cause:** with separation/complex junction flow, the coupled update can't find a direction that reduces the unsteady residual → line search returns the minimum step for many iterations.
- **Fix:** run **decoupled** through those regions (default). Coupling helps only late, after transients settle.
- **Evidence:** Yildirim §5.6 (strut-braced-wing, wing-strut junction with standing shock + separation): the **coupled variant stalls before 1e-5**; decoupled converges without difficulty (Table 7). Explicit conclusion: even for tightly-coupled flow/turbulence, decoupling wins where residuals fluctuate.

---

## Regime rules (when defaults don't apply)

- **Supersonic / moving shocks:** `ANKPhysicalLSTol` 0.2→0.4–0.6 (larger physical changes settle shocks). Consider dissipation-based continuation (`useDissContinuation`, `dissContMagnitude≈⅕–½·M∞`). Evidence: Docs; Seraj2023b.
- **Very low Mach:** ideally use an incompressible code. If not: expect many nonlinear iters; try the 2nd-order switch and a *lower* linear tolerance; set `acousticScaleFactor≈M∞` / `vis4` tuning for low-Mach dissipation. Evidence: Docs; Seraj2023b.
- **Massive separation (up to 90° AoA converged):** more `nSubiterTurb`, start with turbDADI/turbKSP, decoupled. Evidence: Docs; Yildirim §5.6.
- **Actuator/powered regions:** ramp momentum source terms gradually as the solve converges. Evidence: Docs.
- **Complex full-aircraft (tail/nacelle/pylon, heavy overset):** more nonlinear iters + more expensive linear solves are *normal* (inter-block couplings + large/small cell-volume couplings worsen conditioning). Budget ~200 nonlinear iters / up to ~5k linear for 4–5 orders. Evidence: Docs; Yildirim §3.9.
- **2-D / tiny / very large meshes:** minor default adjustments expected; parameters were tuned for millions–tens-of-millions of cells. Evidence: Yildirim §3.9.
- **SA-GR transition model (this branch):** use the validated ANK→CANK→CSANK(→NK) recipe in `../convergence-strategy.md` — it overrides the transonic-SA rules above.

## One-line rules of thumb (memorize these)

- **Scope (2026-08-12):** these are transonic-SA rules; on SA-GR two are inverted (coupling is required; ANK-only floors at ~3.5e-8, not 1e-12) — see `../convergence-strategy.md`.
- Robust default matrix-free operator = **`R1`**; only reach for `R0`/2nd-order when nonlinear convergence stalls with full steps and healthy `Lin Res`.
- **Decouple turbulence before the SA "hill," couple after** (if at all).
- **Delay the NK switch** when unsure; ANK-only can finish hard cases.
- **Overset ⇒ matrix-free everywhere** (primal operator + adjoint `useMatrixFreedrdw=True`).
- **High `Lin Res` = PC problem** (strengthen PC / cut CFL); **tiny `Step` = update-direction/CFL problem** (2nd-order / lower CFL limit / relax LS).
- **Production gradients: Jacobian-free + AD, `frozenTurbulence=False`**; FD-Jacobian only if 0.1% error is acceptable.
