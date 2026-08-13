# MP_06 digest — Piotrowski, PhD thesis (U. Toronto), Transition prediction for aerodynamic shape optimization

Companion to [`MP_06_piotrowski_phd_thesis_full.md`](MP_06_piotrowski_phd_thesis_full.md) (156 pages).

The superset of every paper in this folder, plus three appendices with material
that appears nowhere else. **Chapter 3 §3.1 is the solver chapter and it is the
most useful text we have for our NK problem** — it gives the actual pseudocode
for the damping, line search, and source-Δt machinery that the papers only
describe.

## Map (read only what you need — rule 9)

| Ch. | Contents | Value to us |
|---|---|---|
| 1 | Intro, transition mechanisms, prediction methods | background |
| 2 | **Transition model equations** — SA-LM2015 → SA-sLM2015 (smoothing, intermittency sources, correlations) → compressibility corrections | equivalent to MP_02 + MP_05 |
| **3.1** | **NK-Schur solver: transition considerations, scaling, update damping, source-Δt** | **highest — see below** |
| 3.2 | Gradient evaluation / adjoint | = MP_04 |
| 4 | Analysis results (NLF0416, S809, NLF2-0415, Sickle, CAST10-2, VA-2, CRM-NLF) | validation targets |
| 5 | Optimization results | = MP_04 |
| 6.2 | **Recommendations** — his own open problems | see §5 below |
| A | `F_onset,scf` validation on the **TU Braunschweig Sickle Wing** | our campaign 07 geometry |
| B | Iterative/grid convergence, linearization strategies | = MP_03 |
| **C** | **Transition length modification** | see §4 below |

---

## 1. §3.1.2 Equation and variable scaling — the exact `S_r` matrix (Eq. 3.7)

**This resolves the open item that has been blocking us.** Our
`transitionRowVolScale` (`S_r`, Eq. 58 geometric row scaling) stalls the NK
linear solve dead — lin res pinned 0.99–1.00 — and
[`../SA_GAMMA_RETHETHA_BASE/SAGR_02_adflow_vs_paper_solver.md`](../SA_GAMMA_RETHETHA_BASE/SAGR_02_adflow_vs_paper_solver.md)
§5 flagged the `volRef`-vs-paper's-`J` correspondence as "a best-effort
interpretation" and probably wrong. The thesis writes the matrix out (verified
against PDF page 61):

The scaled system solved is

```
Sa · Sr · (T + A) · Sc · Sc⁻¹ ΔQ = − Sa · Sr · R(Q)        (Eq. 3.6)
```

with, for node *i*:

```
                ⎡ J_i^(2/3)                                                    ⎤
                ⎢     ⋱                                                        ⎥
        S_r,i = ⎢         J_i^(2/3)                                            ⎥   (5 mean-flow rows)
                ⎢                    ν̃max⁻¹ · J_i^(−1/3)                       ⎥
                ⎢                                 γmax⁻¹ · J_i^(−1/3)          ⎥
                ⎣                                          R̃eθt,max⁻¹ · J_i^(−1/3) ⎦

        S_c,i = diag( 1, 1, 1, 1, 1, ν̃max, γmax, R̃eθt,max )

        ν̃max = 1e3      γmax = 10      R̃eθt,max = 1e4
```

Three things to take from this:

1. **The mean-flow and turbulence/transition rows are scaled by opposite powers
   of `J`** — `J^{+2/3}` versus `J^{−1/3}`. Any implementation that applies the
   same-signed exponent to all rows is not this scaling. Our
   `setRVec` in `NKSolvers.F90` uses `vol^(5/3)` / `vol^(2/3)` — same sign,
   both large, and `J` for an SBP curvilinear discretization is the **metric
   Jacobian**, which scales like an **inverse** volume. If `J ~ 1/vol`, the
   paper's exponents translate to `vol^(−2/3)` on the mean-flow rows and
   `vol^(+1/3)` on the transition rows. With cell volumes spanning ~1e-12 on
   our meshes, that discrepancy is astronomically large — a plausible complete
   explanation for the observed 0.99–1.00 stall, and consistent with "no NaN,
   just no progress."
2. **`S_r` and `S_c` are a matched pair.** `S_c` is the *inverse* variable scale
   (`ν̃max, γmax, R̃eθt,max`), `S_r` carries the *reciprocal* of the same
   constants. Our NK column scaling was implemented separately in 2026-07-13
   and `S_r` bolted on later behind its own option — if the two are not using a
   consistent set of `max` constants, they will not cancel the way Eq. 3.6
   intends. The values are not arbitrary: `ν̃max = 1e3` came out of "a thorough
   algorithm optimization procedure" (Osusky & Zingg); `γmax = 10` and
   `R̃eθt,max = 1e4` were set "following a similar approach".
3. **The residual used for convergence monitoring is the *partially* scaled one**
   — `S_r·R(Q)`, not `Sa·Sr·R` and not the raw residual. Their `Rd^(n)` (their
   Eq. 3.5, the quantity every "12–13 orders" claim in these papers refers to)
   is `||S_r R(Q^n)||₂ / ||S_r R(Q^0)||₂`. If our reported `totalRes` uses a
   different scaling, **our orders and their orders are not the same number**,
   which matters for every comparison we have made to the paper.

**Concrete cheap test:** flip the sign of the `S_r` exponents in `setRVec` and
rerun the existing 10-iteration S_r test from
`r1c_CSANK_entry_rel1e-6_dp.cgns`. This is the §8.3 open code item and we now
have the target values.

## 2. §3.1.3 Solution update damping — Algorithms 2, 3, 4

The whole reason NK survives in their solver. Our NK has **no physicality check
at all** (SAGR_02 §7 open item #2) — this is the missing piece, verbatim.

**Algorithm 2 — physics-based restriction (density & energy):**

```
δ_local ← 1
for each local node j, for each mean-flow variable l ∈ [1,5]:
    if ΔQ_{j,l} < 0:  δ_local = max( ΔQ_{j,l}/(0.90·Q_{j,l}) , 1/δ_local )⁻¹   # local min damping
δ_phys ← MPI_Allreduce(δ_local, MIN)                      # GLOBAL minimum
if δ_phys < 0.01:                                         # poor update
    if Δt_ref > Δt_ref,min:
        Δt_ref = max(0.5·Δt_ref, Δt_ref,min);  δ_phys ← 0;  RETRY the iteration
    else:
        δ_phys = 0.01                                     # clip and continue
```

Rule: **negative updates are limited to 90 % of the current value** of ρ and E
(Eq. 3.9: `|ΔQ_ρ,e| ≤ 0.90·Q_ρ,e`), "a good balance between speed and
robustness". Note the **reject-and-retry-with-halved-Δt_ref** branch — a
failed step is not accepted at a tiny λ; it is thrown away and retried on a
more diagonally dominant system. That is different in kind from our
`LSCubic` behaviour of "just take minlambda", which is the E1 symptom we
worked around by relaxing `alpha`.

**Algorithm 3 — transition-variable bounds by damping, not clipping:**

```
θ_fac ← 0.99
Q_{j,l}^{n+1} = Q_{j,l}^n + δ_phys·ΔQ_{j,l}^n
m = 0;  while (γ^{n+1} > 2 or γ^{n+1} < 1e-10) and θ_fac^m > 0.01:
            m += 1;  γ^{n+1} = γ^n + δ_phys·θ_fac^m·Δγ^n
m = 0;  while (R̃eθt^{n+1} < 20) and θ_fac^m > 0.01:
            m += 1;  R̃eθt^{n+1} = R̃eθt^n + δ_phys·θ_fac^m·ΔR̃eθt^n
```

with the standing warning:

> "Clipping the transition model variables at these bounds was also
> investigated; however, **placing a hard limit on variable bounds can
> potentially lead to stalling of the nonlinear solver.**"

and the reassurance that with the source-Δt restriction plus first-order upwind
convection, "**the upper and lower bounds specified in Algorithm 3 are rarely
reached**" — γ's steady state respects [0.02, 1] implicitly through the source
terms. **If our bounds are firing often, that is a symptom, not a safety net.**
Note the bounds themselves are loose: γ ∈ (1e-10, 2), R̃eθt ≥ 20 — deliberately
wider than the physical range, because they exist to catch divergence, not to
enforce physics.

**Algorithm 4 — unsteady-residual backtracking line search:**

```
if δ_phys > 0.01:
    δ_ls ← δ_phys
    R_unst = || T·δ_ls·ΔQ + R(Q + δ_ls·ΔQ) ||₂                    (Eq. 3.11)
    while R_unst ≥ ||R(Q^n)||₂ and δ_ls > 0.01:
        δ_ls = 0.90·δ_ls;   recompute R_unst
    if δ_ls < 0.01:
        if Δt_ref > Δt_ref,min: halve Δt_ref; Q^{n+1} ← Q^n; RETRY
        else: continue with the backtracked state
```

The acceptance criterion is **decrease of the *unsteady* residual**
`||TΔQ + R(Q+ΔQ)||₂` — the quantity that would decrease if the ODE were
integrated time-accurately — **not** an Armijo condition on `||R||`. This is a
meaningfully different, and much more forgiving, test than `LSCubic`'s
`alpha`-based sufficient-decrease, which we had to relax by hand (SAGR_02 §7)
and which cost us a SEGV when relaxed too far. Their geometric backtracking is
a plain 0.90 factor with a 1 % floor.

They also record what they *rejected*: solution-limited time-stepping (Chisholm
& Zingg; Lian et al.) "can slow convergence for cases with stable updates and
lead to **restrictive time steps for variables approaching zero**, which is
common for the intermittency and eddy-viscosity-like variables."

## 3. §3.1.1 / §3.1.4 Coupling and the source-Δt restriction

Same conclusions as MP_03, stated more bluntly:

> "the results demonstrate that **loosely coupled and decoupled solution
> strategies can result in the solution stalling in the approximate-Newton
> phase**, and that the mean-flow, turbulence, and transition model equations
> **require a tight coupling in order to achieve deep nonlinear convergence**.
> However, a strategy is required to stabilize the implicit treatment of the
> large transition model source terms, where the source terms **can remain
> active throughout the inexact-Newton phase**."

Source-Δt restriction (Eq. 3.13–3.14), with the full 3×3 source Jacobian:

```
A_source = ∂(S_ν̃, S_γ, S_R̃eθt) / ∂(ν̃, γ, R̃eθt)      (3×3, block-diagonal per node)
λ_source = largest positive eigenvalue, via a QR algorithm
constraint:  λ_source · Δt_{j,k,m} ≤ 0.8
```

Deactivation: active throughout the approximate-Newton phase; switched off after
**five successive inexact-Newton iterations with no damping active**
(`δ_phys = δ_ls = 1`); **reactivated** if damping is needed again **or if the
relative residual drop rises back above the approximate→inexact switch
tolerance.** Our NK implementation has the damping-triggered half of this but
deliberately omits the residual-rise leg (module scoping, per SAGR_02 §4) —
worth revisiting now that we know it is not optional in their design.

Also worth knowing: they report the source-Δt restriction **also improves
robustness of plain fully-turbulent SA runs** on complex geometries (TTBW,
HWB) — it is not a transition-only device.

## 4. Appendix C — transition-length modification (reduces streamwise grid cost)

Unpublished elsewhere. To cut the streamwise-resolution requirement identified
in MP_04, γ's production and destruction are **merged into one source term** and
`F_onset` is reshaped to stay active over a longer stretch:

```
S_γ = φ_−300(Ω, M∞ M∞ √Re∞ / 20) · (0.98·F_onset − γ + 0.02)          (C.1)

F_onset   = φ_−300( φ_300(F_onset2, 0), 1 )                           (C.2)
F_onset2  = 0.1575·( F_onset1 − 4/7 )                                 (C.3)
F_onset1  = sqrt( (10·ReS/(2.6·Reθc))² + 3e3·RT·sqrt(μ/(ρ·Ue))^(3/2) )… (C.4 — see PDF)
```

`F_length` and `F_turb` are **removed entirely**. The bounds 0.02 and 1 on γ are
enforced implicitly by the source form itself. Rationale: transition length was
controlled by `F_length`, but shrinking `F_length` moved transition *onset*
downstream instead of lengthening the region; so the length is controlled
through `F_onset`'s sensitivity to the eddy-viscosity ratio instead, with a
normalization (by `Ue^{-1/2}` and `sqrt(ν_t)`) intended to hold the transition
length constant across Mach and Reynolds number.

Status in his own words: **promising but incomplete** — "future work must be
completed to extend this modification to include the stationary crossflow
source terms", and it does **not** fix optimality convergence at high Re.
Filed here as a known option, not a recommendation. Equation C.4 is
mangled in the text extraction — read PDF p. 130 before implementing.

## 5. §6.2 Recommendations — his open problems (i.e. what is *not* solved)

Useful as a reality check on our own expectations:

1. **Turnaround time.** Two extra transport equations "significantly increase
   the memory requirements and computational cost"; free-transition runs "still
   often require significantly more linear and nonlinear iterations to
   converge". His proposed remedies are better domain decomposition and
   *simplifying the model* (one-equation / algebraic transition models) — not
   a better preconditioner.
2. **Validation above Re 15e6** is unavailable, not merely undone.
3. Görtler and attachment-line instabilities are not modelled at all.
4. **Streamwise grid requirements** for optimization — Appendix C is the
   partial answer.
5. **Optimality convergence at high Re is unsolved**, suspected to be gradient
   noise, with the **local pressure-gradient parameter formulation** the named
   suspect.

## Shortlist for us

1. **`S_r` exponents** (§1) — the concrete numbers we were missing. Test by
   sign-flipping `setRVec`. Also check `ν̃max/γmax/R̃eθt,max` consistency between
   our row and column scaling, and check whether our reported residual is the
   `S_r`-scaled one.
2. **Algorithm 2 for ρ/E in NK** (§2) — the missing physicality check, with the
   reject-and-halve-Δt_ref branch. This is SAGR_02 §8 item 2, and here is the
   pseudocode.
3. **Unsteady-residual line search** (§2) as an alternative acceptance test to
   `LSCubic`'s hand-relaxed Armijo `alpha`.
4. Source-Δt reactivation on **residual rise**, not only on damping activity (§3).
