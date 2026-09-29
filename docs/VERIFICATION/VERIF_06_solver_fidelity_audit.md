# VERIF_06 — Solver fidelity vs Piotrowski & Zingg (F-tags)

Code comments in `NKSolvers.F90`, `inputParam.F90` and `pyADflow.py` cite the
tags below (`VERIF_06 Fn`). Each tag is a place where ADflow's solver was
compared with the P&Z solver (MP_03 solver paper, MP_06 thesis §3.1,
Algorithms 2–4). All options default to ADflow's original behaviour unless
stated.

## Findings → options

| Tag | Finding | Option(s) | Default |
|-----|---------|-----------|---------|
| F0 | The thesis's "inexact-Newton" phase `(T + A_2nd)ΔQ = −R` is ADflow's **SANK/CSANK**, not NK. NK is the `T → 0` endpoint. So paper machinery belongs in CSANK, and letting CSANK finish (`useNKSolver False`) is the paper's algorithm, not a workaround | — | — |
| F1 | The ANK CFL floor `ANK_CFLMin` is itself ramped by convergence, so the CFL cutback has no room to act at depth (the thesis uses a fixed `Δt_ref,min`) | `ANKCFLMinCap` (absolute cap on the ramped floor; ≤ 0 = off) | `0.0` |
| F2 | ADflow's unsteady line search (factor 0.7, 12 iters) floors at 0.7¹² = 0.0138, above the 0.01 rejection threshold, so an exhausted backtrack is never treated as failure. Thesis Alg. 4: factor 0.90, ~44 iters, reject-and-retry | `ANKUnsteadyLSFactor`, `ANKUnsteadyLSMaxIter`, `ANKRejectOnLSExhausted` | `0.7`, `12`, `False` |
| F3 | Run configuration, not code: swept-wing runs used a far earlier coupling/NK switch than the validated ladder. See [CORE_02](../CORE_BASE/CORE_02_convergence_strategy.md) | — | — |
| F4 | `ANKPhysicalLSTolReTheta` and `omegaMinGamma` were missing from `adflow.pyf` (silent no-ops). **Fixed**: both are in the `anksolver` block | `ANKPhysicalLSTolReTheta`, `omegaMinGamma` | `0.99`, `0.05` |
| F5 | `transitionSrcDtLimit` is 0.9; the thesis (Eq. 3.14) and MP_03 use 0.8. Parameter only | `transitionSrcDtLimit` | `0.9` |
| F6 | NK has **no ρ/E physicality check** (ANK has `physicalityCheckANK`); its only safety valves are NaN detection and the γ/Re̅θt clamp in `applyNKAlgorithm2Damping`. This is why `LSCubic`'s Armijo alpha is relaxed only to 1e-3 | `NKLSRelax` (alpha 1e-3 / turb pre-limit 3.0; off = upstream 1e-2 / 2.0) | `True` |
| F7 | NK and CANK each held half of the paper's Alg. 2/4: NK had per-node γ/Re̅θt damping, CANK had the physicality check and unsteady LS. Per-node damping now also exists in the coupled path (`applyANKAlgorithm2Damping`) | `ANKAlgorithm2Damping` | `False` |
| F8 | In CANK γ was the variable binding the **global** step in almost every iteration (a few front cells throttle the whole field). `False` removes γ/Re̅θt from the global λ and bounds them per node via `applyANKAlgorithm2Damping` (enabled automatically) | `ANKTransitionGlobalLambda` (`True` = global λ, ADflow behaviour) | `True` |
| F9 | PETSc chooses the MFFD step `h` assuming the residual is exact to machine ε; the SA-GR residual is far noisier (smoothMinMax, correlations, tanh), so `J·a` becomes cancellation noise | `MFFDFunctionError` (≤ 0 keeps PETSc's assumption), `MFFDType` | `0.0`, `"ds"` |
| F10 | Only ρ and ρE are O(1) in the scaled state; momenta and turbulence/transition variables sit 10–70× lower, so one MFFD `h` cannot serve them all. Option scales every variable to unit RMS | `ANKColScaleUnit` | `False` |
| — | Stall diagnostics: reports why the ANK/CANK/NK step controller collapsed the step (physicality, LS, CFL). Diagnostic only | `solverStallDiag`, `solverStallDiagStep` | `False`, `1.0` |

Also noted: `getNKColScale` gates on `transitionNK .and. transitionNKActive`,
`getFullColScale` on `transitionNK` only — they diverge only if the NK
auto-disable latch trips (`transitionNKAutoDisableTol = 0` never trips it).

## Validation: NLF(2)-0415 infinite swept wing

Matrix dissipation (`vis4 0.04`, `Vl = Vn = 0`), crossflow on, 65k-node
cross-section (paper: 81.6k). All seven Re converge (relL2Drop 1.6e-10 to
2.1e-9, ~500–635 iterations). Transition location x_tr/c:

| Re/1e6 | ours  | paper | exp   | vs paper | vs exp  |
|--------|-------|-------|-------|----------|---------|
| 1.796  | 0.680 | 0.710 | 0.780 | −4.2 %   | −12.8 % |
| 2.000  | 0.604 | 0.623 | 0.731 | −3.1 %   | −17.4 % |
| 2.204  | 0.531 | 0.541 | 0.578 | −1.9 %   | −8.1 %  |
| 2.370  | 0.484 | 0.486 | 0.504 | −0.4 %   | −3.9 %  |
| 2.498  | 0.447 | 0.444 | 0.446 | +0.7 %   | +0.1 %  |
| 3.000  | 0.339 | 0.349 | 0.327 | −2.9 %   | +3.6 %  |
| 3.498  | 0.266 | 0.266 | 0.297 | +0.1 %   | −10.5 % |

Mean |error| 1.9 % vs the paper, 8.1 % vs experiment. The largest deviation
from experiment (Re 2.0e6, −17.4 %) is where the paper reports its own ~14 %
error, attributed to the experimental data.

## Lessons

- **Reynolds length is the normal chord (1.0), not √2.** The section is the
  NLF2-0415 airfoil extruded at 45° sweep; Re is based on freestream speed and
  the chord normal to the leading edge. With √2 every run was a factor √2
  below the paper — and the resulting error, growing monotonically with Re,
  looked exactly like weak crossflow. A monotonic error trend points at a
  wrong parameter before a wrong model.
- **Check for a pinned CFL ceiling first.** CFL sitting exactly at
  `ANKCFLLimit` with `Step = 1.00` and a healthy linear residual means the
  ceiling is saturated, not that the solve has settled. Raising `ANKCFLLimit`
  to 1e8 took the swept wing from a 1e-5 plateau to converged in 534
  iterations; F1, F2, F8, F9 and F10 had moved nothing because the binding
  constraint was elsewhere. Tightening the linear tolerance alone never helps
  in that state.
