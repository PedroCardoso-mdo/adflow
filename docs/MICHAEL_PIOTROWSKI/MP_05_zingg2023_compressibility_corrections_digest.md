# MP_05 digest — Piotrowski & Zingg (2023, Aeronautical Journal), Compressibility corrections (SA-sLM2015cc)

Companion to [`MP_05_zingg2023_compressibility_corrections_full.md`](MP_05_zingg2023_compressibility_corrections_full.md).

Two deliverables: (a) compressibility corrections `ψ` (Tollmien-Schlichting)
and `ψ_scf` (stationary crossflow), extending the model to 0.71 ≤ M ≤ 0.856;
(b) **a restructuring of how crossflow enters the model**, replacing `D_scf`.

> **Two items here bear directly on our open crossflow stall** (CLAUDE.md rule
> 3: `transitionCrossflow=True` plateaus at ~3e-2 on the tutorial wing).
> See §2 and §3 below — one is an architectural change, one is a
> two-line runtime rule.

## 1. Appendix A is the cleanest full statement of the model

**Practically important:** the 2023 Appendix A writes out the whole
SA-sLM2015 system (A.1–A.39) in **non-dimensional form**, including the
freestream Reynolds number `Re∞` factors in every diffusion, source, and
Reynolds-number definition. Since our biggest recurring class of bug is
non-dimensionalization (CLAUDE.md rule 11), this appendix is a better
cross-check for `Re∞` placement than the 2020 paper's presentation.

Selected forms worth having in front of you (transcript lines ~1253–1430 of the
`_full.md`; **verify against the PDF before implementing**):

```
Re̅θt eq:   ∂R̃eθt/∂t + u_j ∂R̃eθt/∂x_j = Pθt + Dscf + (1/Re∞) ∂/∂x_j[ σθt(ν+νt) ∂R̃eθt/∂x_j ]
           Pθt  = (cθt/t)(Reθt − R̃eθt)(1 − Fθt)
           Fθt  = Fwake·exp(−(d/δ)^4),   Fwake = exp(−ReS×1e-6),   ReS = ρd²S/μ · Re∞
           t    = (500μ)/(ρU²) · (1/Re∞)
           cθt = 0.03, σθt = 2.0, ccrossflow = 0.6

γ eq:      Pγ = ca1·Flength·Fonset·φ_−300(Ω, U∞M∞√Re∞ / (l·20))·√γ·(1 − ce1·γ)
           Eγ = ca2·Fturb·φ_−300(Ω, U∞M∞√Re∞ / (l·20))·γ·(ce2·γ − 1)
           Fonset   = [tanh(6(Fonset,1 − 1.35)) + 1] / 2
           Fonset,1 = sqrt( (ReS/(2.60·Reθc))² + (RT)² )
           Reθc     = 0.67·R̃eθt + 24·sin(R̃eθt/240 + 0.5) + 14
           Flength  = 44 − [44 − (0.50 − 3e-4(R̃eθt − 596))] / (1 + Flength,1)^(1/6)
           Flength,1= exp(−3e-2 (R̃eθt − 460))
           Fturb    = (1 − Fonset)·exp(−RT)
           ce1=1.0, ca1=2.0, σf=1.0, ce2=50, ca2=0.06

SA:        Pν̃ = γ·(cb1/Re∞)·S̃ν̃      ← γ multiplies SA production only (our rule 2)
```

Note the **vorticity limiting** `φ_−300(Ω, U∞M∞√Re∞/(l·20))` sitting inside
both γ source terms — a smooth min of vorticity against a
freestream-referenced ceiling. The paper says this was "re-written in Appendix
A **to clarify the implementation in dimensional solvers**", which is exactly
the audience we are. If our γ sources lack this limiter, that is a real gap,
and it is a *convergence* device, not a physics one.

Also note `Flength` here is the **single smooth expression** (A.34–A.35), not
the 3-Gaussian form of the 2019 paper.

## 2. Crossflow moves from `D_scf` (Re̅θt eq.) to `F_onset,scf` (γ eq.) — §2.2

This is an architectural change and it is the one to weigh against our stall.

**Old (what we implement):** crossflow acts as an extra sink `D_scf` in the
Re̅θt transport equation, pulling `R̃eθt` down in high-helicity regions, which
then feeds transition onset indirectly through `Reθc` inside a single `Fonset`.

**New:** `D_scf` is **removed**. `Fonset` is **split into two onset functions**
— one for 2D mechanisms, one for stationary crossflow — combined by a smooth
max:

```
Fonset,scf   = [tanh(6(Fonset,scf,1 − 1.35)) + 1] / 2
Fonset,scf,1 = sqrt( (ReS / (2.60·Reθc,scf))² + (RT)² )
Reθc,scf     = 0.623 · Reθt,scf                 (calibrated on NLF2-0415)
Fonset      ← φ_300(Fonset, Fonset,scf)          (smooth max, Eq. 17)
```

The stated motive is to let each mechanism carry **its own** compressibility
correction. But the structural consequence matters to us independently of
compressibility: crossflow no longer perturbs the **transported** variable
`R̃eθt` — it acts algebraically on onset. A sink in a transport equation that
fights the `Pθt` production term is precisely the kind of stiff,
counter-directed source pair that produces a residual plateau; moving it out of
the transport equation removes that interaction. Whether that is *our* ~3e-2
plateau is untested, but it is a hypothesis with a concrete implementation
behind it.

The new `Fonset,scf` formulation was validated on the **TU Braunschweig Sickle
Wing** (thesis Appendix A) — the same geometry family as our campaign 07.

## 3. Crossflow initialization dependence — a documented, named pathology

Read this one closely; it describes our symptom almost exactly:

> "Simulations of transonic swept wings revealed that **both crossflow
> source-term approaches, `D_scf` and `F_onset,scf`, produce a transition front
> that is dependent on the initialisation of the flow field.** A solution
> initialised with far-field conditions converges to a transition front
> **upstream** of the converged transition front produced when initialised with
> a converged solution obtained with the crossflow correlation inactive. …this
> behaviour **does not appear in subsonic transition test cases with crossflow**,
> such as the NASA NLF2-0415 Infinite Swept Wing and the TU Braunschweig Sickle
> Wing."

Their fix is a **staged activation**, and it is cheap for us to try:

> "To prevent these inaccurate upstream transition fronts, **the crossflow
> source term is activated after the total residual drops several orders of
> magnitude** without the crossflow source term active. **A relative residual
> drop tolerance of five orders of magnitude** was found to be sufficient…
> This strategy produces the same results as initialising with a converged
> simulation performed without the crossflow correlations active."

So: run to rel 1e-5 with crossflow **off**, then switch it on. That is a runtime
staging rule, implementable as an option (a `transitionCrossflowSwitchTol`
analogous to `ANKCoupledSwitchTol`), or approximated today by a two-leg run with
a restart. Given that our crossflow default was flipped to OFF because of a
plateau reached **from freestream initialization**, this is the first thing I
would test before concluding the crossflow implementation is broken.

Caveat the paper is explicit about: they report this for **transonic** swept
wings, and say it does *not* appear subsonically. Our stalling case is M=0.15
tutorial wing — subsonic. So this is a strong lead, not a diagnosis.

## 4. The corrections themselves

**Tollmien-Schlichting** (applied to `Reθc`, *not* to `Reθt`):

```
Reθc,comp = ψ · Reθc
ψ = (a1·Me² + a2·Me + a3) · exp(b1·λθ,e·Me) · exp(c1 + c2·Tu∞·Me)
a1 = 0.34, a2 = −0.38, a3 = 1.00; b1 = 3.00; c1 = 0.41, c2 = −0.27
```

Why applied to `Reθc` and not `Reθt` — a real implementation trap:

> applying it to `Reθt` "**significantly over-predicted the laminar extent**…
> The rapid acceleration at the upper surface leading edge produced large
> upstream values of the compressibility correction, and therefore `R̃eθt`.
> These large values **convect downstream in the boundary layer and delay
> transition**."

i.e. anything you scale that then gets *transported* contaminates the whole
downstream boundary layer. Scale the algebraic quantity instead. (An earlier
variant `ψ_init` with `a1=0.44, b1=5.00` over-predicted laminar extent on VA-2;
`a1` and `b1` were reduced to get `ψ`.)

**Stationary crossflow** (Malik et al.), applied to the scf onset Reynolds number:

```
Reθt,scf,comp = ψ_scf · Reθt,scf ,     ψ_scf = 1 + ((κ−1)/2)·Me²
```

Both `ψ, ψ_scf → 1` as `Me → 0`, so **neither changes incompressible results** —
they can be implemented and left on without invalidating our subsonic
validation. TS is the stronger correction, especially in favourable gradients.

**Boundary-layer edge quantities**, computed locally by isentropic relations
(no boundary-layer edge search):

```
Me = Ue/ae,   ae = sqrt( κ·p∞/ρ∞ · (p/p∞)^(1−1/κ) )
Ue = sqrt( U∞² + (2κ/(κ−1))·(p∞/ρ∞)·(1 − (p/p∞)^(1−1/κ)) )
λθ,e = (ρθ²/μ)(dUe/ds)·Re∞ ,   dUe/dx_i = −(1/(ρ∞ Ue))·(p/p∞)^(−1/κ)·dp/dx_i
```

**Do not confuse the two pressure-gradient parameters.** `λθ,e` (edge-based) is
used **only** to evaluate the compressibility corrections. The original
velocity-gradient `λθ` is retained for the LM2009 correlation itself, "to be
consistent with the Langtry-Menter model, as it provides better agreement with
experimental data for zero-pressure-gradient flat plate cases". Both are
limited to −0.1 ≤ λ ≤ 0.1; neither is Galilean invariant.

The reason an edge-based `λ` is needed at all: the corrections "are most active
in the **middle of the boundary layer**", where the standard local `λθ` is
invalid — the model normally hides this by using `F_θt` to switch the
correlation off inside the boundary layer.

## 5. Reported cost

> "These corrections and modifications **do not impact the predictive
> capability of the model in the incompressible flow regime and do not have a
> significant impact on its iterative and grid convergence behaviour.**"

So `ψ`/`ψ_scf` are cheap and convergence-neutral. The `D_scf` → `F_onset,scf`
restructuring is the part with real convergence implications, in both
directions.

## Implementation shortlist for us

1. **Test the staged crossflow activation** (crossflow off → rel 1e-5 → on)
   against our tutorial-wing plateau. Cheapest possible test of a documented
   pathology; no code change if done as a two-leg restart.
2. **Check for the vorticity limiter** `φ_−300(Ω, U∞M∞√Re∞/(l·20))` in our γ
   source terms.
3. If (1) does not resolve the plateau, the `D_scf` → `F_onset,scf` restructuring
   is the paper's own answer, and it is the version validated on the sickle wing.
   *(AD-relevant: touches `saGammaRetheta.F90` Source → `TAPENADE NEEDED`.)*
4. `ψ`/`ψ_scf` are only worth implementing when we go transonic — they are inert
   below M≈0.3.
