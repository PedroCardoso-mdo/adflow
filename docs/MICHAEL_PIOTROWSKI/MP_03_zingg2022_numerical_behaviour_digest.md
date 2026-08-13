# MP_03 digest — Piotrowski & Zingg (2022, Special Session), Numerical behaviour in a NK flow solver

Companion to [`MP_03_zingg2022_numerical_behaviour_full.md`](MP_03_zingg2022_numerical_behaviour_full.md).

> **This is the single most relevant paper in the folder for our convergence
> problem.** It is a controlled study of exactly the question we are stuck on:
> how to linearize and time-step SA-sLM2015 in a strong implicit solver so it
> converges deeply. 4 coupling strategies × source-Δt on/off × 3 geometries.

## The headline result

**Fully-coupled linearization + source-term Δt restriction. Both. Neither alone.**

Their measurements, in one table (rel. residual = L2 of the *row-scaled* total
residual normalised by iteration 0; cost in equivalent residual evaluations):

| Case | w/o source-Δt restriction | w/ source-Δt restriction |
|---|---|---|
| **NLF0416** L4 (M 0.10, Re 4e6, α 0°) — easy | fully-coupled **fails to converge**; loose/decoupled OK | all four converge; *fully-coupled is the slowest* |
| **VA-2** L3, α −0.40° (M 0.71, Re 10e6) — no shock | all four cost more | all four similar |
| **VA-2** L3, α **+1.80°** — shock-induced transition | **all four fail** within the compute allocation | **fully-coupled clearly best** |
| **CRM-NLF** L2 (32.8 M nodes, M 0.86, Re 15e6) | **all four stall after ~4 orders**; fully-coupled exits on an unphysical update | **fully-coupled best**; turb-trans-coupled and trans-coupled **stall after 5 orders**; decoupled still converging when it hit the wall clock |

Three things fall out of this that bear directly on us:

1. **Difficulty inverts the ranking.** On the easy 2D case, decoupling wins and
   full coupling is the slowest — which is exactly the observation our
   `convergence-strategy.md` recorded ("do NOT couple early", CANK at 1e-2
   stagnates). But that ranking **reverses** as the case gets hard. On
   CRM-NLF, the loosely-coupled variants are the ones that stall at 5 orders.
   Our AR5 wing is on the hard end of that spectrum, not the easy end.
2. **Stalling after ~4–5 orders with a partially-coupled Jacobian is the
   documented signature of insufficient coupling** — not of a weak
   preconditioner. Worth checking which of our phases is actually running with
   which block structure when the AR5 case parks.
3. **Without the source-Δt restriction, nothing converges on the hard cases** —
   and the fully-coupled strategy is the one that fails *worst* (unphysical
   update, solver exit). A more accurate Jacobian **needs** the Δt restriction
   to be safe. Full coupling and the restriction are a package.

> "The fully-coupled linearization strategy provides a more complete Jacobian
> and produces a more accurate solution update, which significantly improves
> convergence for this more challenging case… **however, the source-term time
> step restriction is required in order to prevent unstable solution updates.**"

## The four linearization strategies (§II.C)

All four are the **same 8×8 coupled block Jacobian with off-diagonal blocks
zeroed**, never a different solver — deliberately, "to isolate the effects of
varying coupling strategies and avoid differences introduced by using a
different solution algorithm, such as the introduction of different
preconditioners":

| Name | Structure |
|---|---|
| fully coupled | full 8×8 (5 mean-flow + ν̃ + γ + Re̅θt) |
| turb-trans coupled | mean-flow decoupled from a 3×3 {ν̃, γ, Re̅θt} block |
| trans coupled | 2×2 {γ, Re̅θt} block, ν̃ linearized independently |
| decoupled | independent linearization for each equation |

Two implementation points we should mirror:

- **The source-term Jacobian used for the eigenvalue must match the solver's
  coupling level.** "the coupling of the source-term Jacobian used to evaluate
  the source-term eigenvalues **mirroring the solver linearization strategy**."
  If our Eq.-59 restriction always evaluates a 3×3 coupled source Jacobian
  while the solver is running a decoupled or 2×2 variant, the restriction is
  inconsistent with the step it is protecting.
- **The reduced-coupling variants are only used in the globalization
  (approximate-Newton) phase.** "To recover the nonlinear convergence
  properties of Newton's method, the loosely-coupled and decoupled
  linearization strategies are **only adopted in the approximate-Newton
  phase**, with the mean-flow, turbulence, and transition model equations
  **fully coupled in the inexact-Newton phase**." Their inexact-Newton phase =
  our NK. **NK must be fully coupled.** Our ladder (segregated ANK → CANK →
  CSANK → NK) matches this shape, which is reassuring.

## Source-term Δt restriction (§II.D)

- Value used throughout: **0.8** (2019 used 0.5).
- Definition: largest positive eigenvalue of the source-term Jacobian × local
  time step.
- **Deactivation rule:** disabled in the inexact-Newton phase after **five
  successive iterations** in which neither the physicality check nor the
  unsteady-residual line search fired; **reinstated the moment either fires
  again.** Our `srcDtDeactivateIters` implements this counter. The trigger set
  matters — it is the *damping/line-search* activity that gates it, not the
  residual level.

## Solver settings they report (useful as sanity references)

- Approximate→inexact-Newton switch at rel **1e-4** (2D) / **1e-5** (3D).
  Our validated CANK coupling point of 1e-5 sits in exactly this band.
- Convergence tolerances: rel **1e-15**, abs **1e-10** — and "for the majority
  of the results presented the **absolute** residual tolerance is met first",
  i.e. they routinely land 12–13 orders down. This is the gap we are measuring
  ourselves against.
- Artificial dissipation (matrix-based, Swanson–Turkel): linear eigenvalue
  limiter **`Vl = 0`** — "found to be **overly dissipative in the laminar
  boundary layer**"; nonlinear limiter `Vn = 0.25` (0.30 for CRM-NLF).
  2nd/4th-difference coefficients: 0.00/0.04 (NLF0416), 2.00/0.04 (VA-2),
  3.00/0.06 (CRM-NLF). Note dissipation is **raised** for the hard 3D case.
- SA variant: **SA-noft2-neg**; QCR2000 on CRM-NLF.
- First-order upwind on turbulence and transition convection (matches our rule 5).
- The L5 (finest) NLF0416 grid converged but "the additional stiffness from the
  small grid spacings significantly affects the **linear solver performance,
  requiring larger Krylov subspace sizes**". Relevant to our NKSubspaceSize
  sizing on fine AR5 levels.

## Cost scaling (calibrates our expectations)

- Free transition costs **1.2–1.6×** the equivalent residual evaluations of the
  fully-turbulent SA run at the same grid level…
- …**on top of** roughly **2× the cost per nonlinear iteration** from the two
  extra transport equations, which the "equivalent residual evaluations" metric
  hides.

## Grid convergence (secondary for us, but load-bearing for validation)

- Free-transition solutions grid-converge similarly to fully-turbulent ones but
  generally need **more streamwise resolution** for the same accuracy —
  the sensitivity is concentrated in the **transition region**, not at the wall.
- Trailing-edge aspect ratio: TE-AR ≈ **20** was chosen for sharp-TE airfoils
  (200 starves streamwise resolution near transition; 2 was also tested).
  LE-AR held at ~200. Transition onset moves **downstream** with refinement.
- Iterative error must sit **2–3 orders below** discretization error (Eça et al.)
  before any grid study means anything — the reason they chase machine zero.

## What I would take from this paper into our AR5 case

1. Confirm the AR5 stall is not a *coupling* stall: our hard case is the CRM-NLF
   analogue, where partial coupling stalls at ~5 orders and only full coupling
   goes deep.
2. Confirm the source-term Δt eigenvalue is evaluated with the same block
   structure as the active solver phase.
3. Confirm the Eq.-59 reactivation trigger is wired to physicality/line-search
   activity (it is, in `NKStep`), and that our threshold is 0.8.
