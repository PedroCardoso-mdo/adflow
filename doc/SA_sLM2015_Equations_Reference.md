# SA-sLM2015 Complete Equation Reference for ADflow Implementation
## Scraped from Piotrowski & Zingg (2020), AIAA J., DOI: 10.2514/1.J059784

---

## HOW TO USE THIS DOCUMENT

- **Equations are grouped by implementation phase**
- For each equation: check how SA (`sa.f90`) or SST implements the analogous term — the discretisation structure (convection loop, diffusion assembly, source accumulation) is identical; only the coefficients change
- The VS Code AI should read `sa.f90` and `blockette.f90` first, then implement each equation block below by cloning the SA pattern

---

## SOLVER STRATEGY

### Decoupled from mean flow, coupled internally
Mean flow (ρ, ρu, ρv, ρw, ρE) is frozen during the turbulence/transition sub-iteration. The three variables (ν̃, γ, Rẽθt) are solved as a **coupled 3-variable system**.

### ANK for the turbulence block
ADflow likely has an ANK (Approximate Newton–Krylov) path for the turbulence equations in `blockette.f90` — check if a matrix-free Krylov iteration already exists there. If it does, **extending it to 3 variables is the right approach** for this coupled system: it handles the stiff source terms better than DADI for transition models and avoids forming the full 3×3 Jacobian explicitly (matrix-free means only the Jacobian-vector product is needed). If only DADI exists, it can be extended to a 3×3 block-diagonal implicit sweep — check which is present and prefer ANK if available.

### Source-term time step limiter (always implement)
Regardless of ANK or DADI: inside the turbulence sub-iteration, the local time step must be limited by the source-term eigenvalue (see §5 below). This is mandatory for stability with transition source terms.

---

## PHASE 1 — Transport Structure (convection + diffusion, zero sources)

### New state variables
```
w(iGamma)    = γ          intermittency             range [1e-10, 1] steady, allow up to 2
w(iReThetat) = Re~_θt     transported Re_theta_t    range [20, ~2000]
```
Add immediately after the SA variable index definition.

---

### Equation 1 — Re~_θt Transport

```
∂Re~_θt/∂t  +  uj ∂Re~_θt/∂xj  =  Pθt  +  Dscf  +  (1/Re) ∂/∂xj [ σθt(ν + νt) ∂Re~_θt/∂xj ]

  σθt = 2.0
```

**SA analogue:** SA transport equation — same LHS structure.  
**Diffusion coeff vs SA:** SA uses `(ν + ν̃)/(σ_SA·Re)`. Here use `σθt·(ν + νt)/Re` where `νt = fv1·ν̃` (physical eddy viscosity, not ν̃ directly). Find where `mut` is assembled in SA and reuse.  
**Phase 1:** set `Pθt = 0`, `Dscf = 0`.

---

### Equation 27 — γ Transport

```
∂γ/∂t  +  uj ∂γ/∂xj  =  Pγ - Eγ  +  (1/Re) ∂/∂xj [ (ν + νt/σf) ∂γ/∂xj ]

  σf = 1.0
```

**SA analogue:** same structure.  
**Diffusion coeff:** `(ν + νt/σf)/Re` with σf = 1.0.  
**Phase 1:** set `Pγ = 0`, `Eγ = 0`.

---

## PHASE 2 — Source Terms

### A. Re~_θt Sources

#### Equation 2 — Production Pθt
```
Pθt = (cθt/t) · (Reθt - Re~_θt) · (1 - Fθt)

  cθt = 0.03
```

#### Equation 7 — Time scale t
```
t = (500·μ) / (ρ·U²) · (1/Re)
```
Note: U = local velocity magnitude.

#### Equation 3 — Fθt (boundary layer sensor)
```
Fθt = Fwake · exp( -(d/δ)⁴ )
```

#### Equation 4 — Boundary layer thickness estimate
```
δBL = (15/2) · θBL
θBL = Re~_θt · μ / (ρ·U) · (1/Re)
δ   = (50·d·Ω/U) · δBL
```
Note: Ω = vorticity magnitude (already computed for SA — reuse).

#### Equation 5 — Fwake (wake sensor)
```
Fwake = exp( -ReS / 1e6 )
ReS   = ρ·d²·S / μ · Re
```
Note: S = strain-rate magnitude (reuse from SA). d = wall distance (reuse from SA).

#### Equations 8–9 — Reθt empirical correlation (target value in freestream)
```
         ⎧ (1173.51 - 589.428·Tu + 0.2196/Tu²) · F(λθ)    if Tu ≤ 1.3
Reθt =  ⎨
         ⎩ 331.50·(Tu - 0.5658)^(-0.671) · F(λθ)          if Tu > 1.3

         ⎧ 1 - [-12.986·λθ - 123.66·λθ² - 405.689·λθ³] · exp(-(Tu/1.5)^1.5)   if λθ ≤ 0
F(λθ) = ⎨
         ⎩ 1 + 0.275·[1 - exp(-35·λθ)] · exp(-Tu/0.5)                           if λθ > 0
```
**Note:** Use Tu_inf (freestream value) uniformly — do not use local Tu. This is an explicit choice from the paper (Medida & Baeder approach).  
**Smooth version of F(λθ):** Equations 54–57 — see §Smooth Functions below.

#### Equations 10–14 — Pressure gradient parameter λθ
```
λθ = (ρ·θBL²/μ) · (dU/ds) · Re

dU/ds = (u/U)·dU/dx + (v/U)·dU/dy + (w/U)·dU/dz

dU/dx = (u·∂u/∂x + v·∂v/∂x + w·∂w/∂x) / U
dU/dy = (u·∂u/∂y + v·∂v/∂y + w·∂w/∂y) / U
dU/dz = (u·∂u/∂z + v·∂v/∂z + w·∂w/∂z) / U

U = sqrt(u² + v² + w²)
```
Clamp: `λθ = max(-0.1, min(0.1, λθ))` to prevent correlation blowup.  
Velocity gradients are already computed for SA viscous terms — reuse.

---

#### Equations 15–22 — Crossflow source Dscf (only when crossflow enabled)
```
Dscf = (cθt/t) · ccrossflow · min(Rescf - Re~_θt, 0) · Fθt

  ccrossflow = 0.6

Rescf = ρ·(U/0.82)·θt/μ · Re  =  -35.088·ln(h·θt) + 319.51 + f(ΔH⁺cf) - f(ΔH⁻cf)

ΔHcf  = Hcrossflow · (1.0 + min(RT, 0.4))      RT = μt/μ

ΔH⁺cf = max(0.1066 - ΔHcf, 0)
f(ΔH⁺cf) = 6200·ΔH⁺cf + 50000·(ΔH⁺cf)²

ΔH⁻cf = max(-(0.1066 - ΔHcf), 0)
f(ΔH⁻cf) = 75·tanh(ΔH⁻cf / 0.0125)
```

#### Equations 23–26 — Helicity / crossflow strength
```
Û = (u, v, w) / U                            unit velocity vector

Ω⃗ = (∂w/∂y - ∂v/∂z,  ∂u/∂z - ∂w/∂x,  ∂v/∂x - ∂u/∂y)    vorticity vector

Ωstreamwise = |Û · Ω⃗|                        streamwise vorticity (helicity)

Hcrossflow  = d · Ωstreamwise / U
```
Protect against U=0: use `U_safe = max(U, 1e-14)`.

---

### B. γ Sources

#### Equation 28 — Pγ (intermittency production) — USE SMOOTH VERSION (Eq. 52)
```
Pγ = ca1 · Flength · Fonset · Ω · sqrt(γ) · (1 - ce1·γ)

  ca1 = 2.0,  ce1 = 1.0,  σf = 1.0
```
**Use smooth Eq. 52 instead** (see below) — vorticity is limited.

#### Equation 34 — Eγ (intermittency destruction) — USE SMOOTH VERSION (Eq. 53)
```
Eγ = ca2 · Fturb · Ω · γ · (ce2·γ - 1)

  ca2 = 0.06,  ce2 = 50
```
**Use smooth Eq. 53 instead** (see below).

#### Equations 29–32 — Fonset (transition trigger) — USE SMOOTH VERSION (Eqs. 46–47)
```
Fonset1 = ReS / (2.6·Reθc)
Fonset2 = min(max(Fonset1, Fonset1⁴), 2)
Fonset3 = max(1 - (RT/2.5)³, 0)
Fonset  = max(sqrt(Fonset2) - Fonset3, 0)
```
**Use smooth Eqs. 46–47 instead** (see below).

#### Equation 35 — Fturb — USE SMOOTH VERSION (Eq. 48)
```
Fturb = exp( -(RT/4)⁴ )
```
**Use smooth Eq. 48 instead**.

#### Equations 37–38 — Flength and Reθc correlations — USE SMOOTH VERSIONS (Eqs. 49–51)
```
         ⎧ 398.189e-1 + (-119.270e-4)·Re~_θt + (-132.567e-6)·Re~_θt²             Re~_θt < 400
         ⎪ 263.404 + (-123.939e-2)·Re~_θt + (194.548e-5)·Re~_θt² + (-101.695e-8)·Re~_θt³   400 ≤ Re~_θt < 596
Flength =⎨ 0.5 - (Re~_θt - 596)·3e-4                                             596 ≤ Re~_θt < 1200
         ⎩ 0.3188                                                                  Re~_θt ≥ 1200

         ⎧ Re~_θt - (396.035e-2 + (-120.656e-4)·Re~_θt + (868.230e-6)·Re~_θt²
         ⎪            + (-696.506e-9)·Re~_θt³ + (174.105e-12)·Re~_θt⁴)           Re~_θt ≤ 1870
Reθc =  ⎨
         ⎩ Re~_θt - (593.11 + 0.482·(Re~_θt - 1870))                              Re~_θt > 1870
```
**Use smooth Eqs. 49–51 instead** — these have kinks at the piecewise boundaries.

---

### C. SA Coupling

#### Equation 41 — Modified SA production
```
P̃_ν̃ = γ · Pν̃
```
In `sa.f90`, find the production term `Pν̃ = cb1·S̃·ν̃/Re`. Multiply by `w(iGamma)` at that cell.  
**This is the only change to `sa.f90` in Phase 2.** Everything else stays in the new transition file.

---

## SMOOTH FUNCTION LIBRARY
### All non-smooth functions above are replaced by these — implement these, not the originals

#### Algorithm 1 — Smooth max/min (p = ±300 throughout)
```
φp(g1, g2, p):
  a = max(g1, g2);   b = min(g1, g2)
  φ_switch = -log(|p| · 1e-15) / |p|

  if p > 0:  (maximum)
    if |a-b| > φ_switch:  return a
    else:                 return a + log(1 + exp(p·(b-a))) / p

  if p < 0:  (minimum)
    if |a-b| > φ_switch:  return b
    else:                 return b + log(1 + exp(p·(a-b))) / p
```
Use `p = +300` for smooth max, `p = -300` for smooth min.

#### Equations 46–47 — Smooth Fonset (replaces Eqs. 29–32)
```
Fonset1 = sqrt( (ReS / (2.6·Reθc))²  +  RT² )

Fonset = ( tanh(6·(Fonset1 - 1.35)) + 1 ) / 2
```
Note: the new Fonset1 includes RT — keeps Fonset active once turbulent BL forms (no relaminarisation, acceptable for transport aircraft).

#### Equation 48 — Smooth Fturb (replaces Eq. 35)
```
Fturb = (1 - Fonset) · exp(-RT)
```

#### Equations 49–50 — Smooth Flength (replaces Eq. 37)
```
Flength1 = exp( -3e-2 · (Re~_θt - 460) )

Flength = 44 - [ 44 - (0.50 - 3e-4·(Re~_θt - 596)) ] / (1 + Flength1)^(1/6)
```

#### Equation 51 — Smooth Reθc (replaces Eq. 38)
```
Reθc = 0.67·Re~_θt  +  24·sin(Re~_θt/240 + 0.5)  +  14
```

#### Equations 52–53 — Smooth Pγ and Eγ with vorticity limiter (replaces Eqs. 28, 34)
```
Pγ = ca1 · Flength · Fonset · φ₋₃₀₀(Ω,  M·sqrt(M·Re)/20) · sqrt(γ) · (1 - ce1·γ)

Eγ = ca2 · Fturb  · φ₋₃₀₀(Ω,  M·sqrt(M·Re)/20) · γ · (ce2·γ - 1)
```
The `φ₋₃₀₀(Ω, M·sqrt(M·Re)/20)` term is a smooth min that limits maximum vorticity.  
M = freestream Mach number, Re = freestream Reynolds number (both are constants in a given run).  
Compute the limiter value once: `Omega_lim = M·sqrt(M·Re)/20`.

#### Equations 54–57 — Smooth F(λθ) (replaces Eq. 9)
```
F1(λθ) = 1 + 0.275·[1 - exp(-35·λθ)] · exp(-Tu/0.5)

F2(λθ) = φ₊₃₀₀(F1(λθ),  1)                              smooth max with 1

F3(λθ) = 1 - [-12.986·λθ - 123.66·λθ² - 405.689·λθ³] · exp(-(Tu/1.5)^1.5)

F(λθ)  = φ₋₃₀₀(F2(λθ),  F3(λθ))                         smooth min
```

---

## SCALING (from paper §IV.B.1, Osusky & Zingg 2013)

Find where SA scaling (`nuTilde_max = 1e3`) is applied — same row/column scaling pattern:

```
Column scale:   Sc(iGamma)    = 10.0
                Sc(iReThetat) = 1e4

Row scale:      Sr(iGamma)    = (1/10.0)  · J^(-1/3)
                Sr(iReThetat) = (1/1e4)   · J^(-1/3)
```
where J = cell volume. Use row-scaled residual for convergence monitoring (not raw residual).

---

## SOLUTION UPDATE DAMPING (Algorithm 2 from paper)

After computing update ΔQ, before applying:
```
theta_fac = 0.99

while γ_new > 2.0  OR  γ_new < 1e-10:
    damp γ update by theta_fac^m

while Re~_θt_new < 20.0:
    damp Re~_θt update by theta_fac^m
```
Find where SA update is damped/clipped in `blockette.f90` or the turbulence update loop and add this alongside it.

---

## SOURCE-TERM TIME STEP LIMITER (paper §IV.B.3, Lian et al. 2010)

The local time step must satisfy:
```
λ_source · Δt ≤ 0.9
```
where λ_source = largest positive eigenvalue of the 3×3 block:
```
Asource = | ∂Sν̃/∂ν̃    ∂Sν̃/∂γ    ∂Sν̃/∂Re~_θt  |
          | ∂Sγ/∂ν̃    ∂Sγ/∂γ    ∂Sγ/∂Re~_θt  |
          | ∂SRe~/∂ν̃  ∂SRe~/∂γ  ∂SRe~/∂Re~_θt |
```
Compute eigenvalues via QR or Gershgorin bound. In **Phase 1** all entries are zero — implement the infrastructure as a no-op stub now, fill in Phase 2.

**ANK→NK switch:** Deactivate source-term time stepping after 5 clean NK iterations; reactivate on backtracking.

---

## BOUNDARY CONDITIONS

| Boundary | γ | Re~_θt |
|----------|---|--------|
| Viscous wall | `1e-10` (Dirichlet) | Zero-gradient (copy interior to ghost) |
| Farfield / inflow | `1.0` | `ReThetat_inf` (from Eq. 8 at Tu_inf, λθ=0) |
| Outflow | Zero-gradient | Zero-gradient |
| Symmetry / periodic | Mirror (same as SA) | Mirror |

---

## INITIALISATION

```
γ        = 1.0            (fully turbulent — SA unchanged until Phase 2 coupling)
Re~_θt   = ReThetat_inf   (from Eq. 8 evaluated at Tu_inf, λθ=0)
```
On restart from old SA file (missing new variables): pad with above values.

---

## PHASE 1 VALIDATION CHECKLIST

1. Unit test all smooth functions against paper Figs 2, 3, 5
2. Flat plate: γ and Re~_θt residuals → machine zero; SA solution bit-for-bit identical to standalone SA
3. NLF0416 α=0°: same checks, smooth convergence histories for new equations
4. Grid refinement (coarse/medium/fine): new variables show smooth grid convergence; lift/drag = SA-only
5. Jacobian finite-difference check on columns iGamma and iReThetat
