# Paper Reference — SA-LM2015 / SA-sLM2015 Transition Model

> **Source**: Piotrowski & Zingg (2020), "Smooth Local Correlation-Based Transition Model
> for the Spalart-Allmaras Turbulence Model", AIAA Journal, Vol. 58, No. 10.
>
> This file contains ALL equations, constants, algorithms, and conventions needed to
> implement the model. It is the single source of truth for "what the paper says."
> When code and paper disagree, paper wins.

---

## 1. Model Overview

Three transport equations coupled to the Spalart-Allmaras (SA) turbulence model:

1. **ν̃** — SA working variable (existing, modified by γ multiplier on production)
2. **γ** — intermittency (new, Eq. 27)
3. **Re̅θt** — transition-onset momentum-thickness Reynolds number (new, Eq. 1)

The intermittency γ multiplies the SA production term:
```
P̃_ν̃ = γ · P_ν̃     (Eq. 41)
```
where P_ν̃ = cb1·(1-ft2)·S̃·ν̃. γ does NOT multiply the SA destruction term.

---

## 2. Re̅θt Transport Equation (Eq. 1)

```
∂(ρ·Re̅θt)/∂t + ∂(ρ·uj·Re̅θt)/∂xj = P_θt + ∂/∂xj [σ_θt·(μ + μ_t)·∂Re̅θt/∂xj]
```

- σ_θt = 2.0 (diffusion coefficient)
- **First-order upwind** recommended for convection (paper Section IV.A)

### 2.1 Source term P_θt (Eq. 2)

```
P_θt = c_θt/t · (Re_θt - Re̅θt) · (1 - F_θt)
```

- c_θt = 0.03 (Eq. 2)
- t = time scale (Eq. 7)
- Re_θt = target value from empirical correlation (Eq. 8)
- F_θt = switching function (Eq. 3)

### 2.2 Time scale t (Eq. 7)

```
t = 500·μ / (ρ·U²) · (1/Re)
```

where:
- μ = dynamic viscosity (dimensional)
- ρ = density
- U = local velocity magnitude
- Re = freestream Reynolds number
- **CRITICAL**: The 1/Re factor is explicit. In nondimensional codes where μ_ref = 1/Re, this becomes `t = 500·μ̃/(ρ·U²)` where μ̃ is the nondimensional viscosity. But if the code stores `rlv = μ/μ_∞` (ratio), then `t = 500·rlv/(ρ·U²·Re)`.

### 2.3 Switching function F_θt (Eq. 3)

```
F_θt = min(max(F_wake · e^(-(y/δ)^4), 1 - ((γ_eff - 1/ce2)/(1 - 1/ce2))^2), 1.0)
```

- γ_eff = max(γ, 0.0) (Eq. 4)
- ce2 = 50.0

### 2.4 Boundary-layer thickness proxies

```
δ_BL = (15/2) · θ_BL          (Eq. 5, Blasius)
θ_BL = Re̅θt · μ / (ρ · U)    (Eq. 5)
δ = (50 · Ω · y / U) · δ_BL   (Eq. 6, where Ω = vorticity magnitude)
```

### 2.5 F_wake (Eq. 3 inner)

```
F_wake = e^(-(Re_ω / (1e5))^2)
```
where Re_ω = ρ·Ω·y²/μ.

---

## 3. Empirical Correlation Re_θt (Eq. 8-9)

```
Re_θt(Tu, λ_θ) = Re_θt0(Tu) · F(λ_θ)   (composite)
```

### 3.1 Re_θt0 (Eq. 8, function of Tu only)

**For Tu ≤ 1.3%:**
```
Re_θt0 = (1173.51 - 589.428·Tu + 0.2196/Tu²)
```

**For Tu > 1.3%:**
```
Re_θt0 = 331.50·(Tu - 0.5658)^(-0.671)
```

Tu is freestream turbulence intensity in percent.

### 3.2 Pressure-gradient function F(λ_θ) — LM2015 (Eq. 10-14, piecewise)

λ_θ = (θ²/ν)·(dU/ds) — Thwaites' pressure-gradient parameter

**For λ_θ ≤ 0 (APG/zero PG):**
```
F(λ_θ) = 1 - (-12.986·λ_θ - 123.66·λ_θ² - 405.689·λ_θ³)·e^(-(Tu/1.5)^1.5)
```

**For λ_θ > 0 (FPG):**
```
F(λ_θ) = 1 + 0.275·(1 - e^(-35·λ_θ))·e^(-Tu/0.5)
```

### 3.3 Pressure-gradient function F(λ_θ) — sLM2015 (Eq. 54-57, smooth)

Replace piecewise with smooth blending using φ_p:

```
F1(λ_θ) = -(12.986·λ_θ + 123.66·λ_θ² + 405.689·λ_θ³)·e^(-(Tu/1.5)^1.5)    (Eq. 55)
F2(λ_θ) = 0.275·(1 - e^(-35·λ_θ))·e^(-Tu/0.5)                                (Eq. 56)
F3(λ_θ) = φ_{-300}(F1, 0) + φ_{+300}(F2, 0)                                   (Eq. 57)
F(λ_θ) = 1 + F3                                                                (Eq. 54)
```

### 3.4 λ_θ computation (Eq. 10-14)

```
θ_BL = Re̅θt · μ / (ρ · U)
λ_θ = (θ_BL² · ρ / μ) · dU/ds
```

where dU/ds is the streamwise velocity gradient (requires velocity gradient projection).

---

## 4. Intermittency γ Transport Equation (Eq. 27)

```
∂(ρ·γ)/∂t + ∂(ρ·uj·γ)/∂xj = P_γ - E_γ + ∂/∂xj [(μ + μ_t/σ_f)·∂γ/∂xj]
```

- σ_f = 1.0 (diffusion coefficient)
- **First-order upwind** recommended for convection

### 4.1 P_γ — Intermittency production (Eq. 28)

**LM2015 form:**
```
P_γ = ca1 · F_length · S · [γ · F_onset]^0.5 · (1 - ce1·γ)
```

**sLM2015 form (Eq. 52, Piotrowski modification — uses vorticity Ω instead of strain S):**
```
P_γ = ca1 · F_length · Ω_lim · F_onset · (1 - ce1·γ)
```

where Ω_lim = φ_{-300}(Ω, M·√(M·Re)/20) — vorticity limited by freestream scaling.

Constants: ca1 = 2.0, ce1 = 1.0

### 4.2 E_γ — Intermittency destruction (Eq. 34)

```
E_γ = ca2 · F_turb · Ω_lim · γ · (ce2·γ - 1)
```

Constants: ca2 = 0.06, ce2 = 50.0

Note: When γ < 1/ce2 = 0.02, the (ce2·γ - 1) term is negative, making E_γ a source
(drives γ toward 0.02). This is the equilibrium point in laminar regions.

### 4.3 F_onset — Transition onset function

**LM2015 form (Eq. 29-32, piecewise):**
```
F_onset1 = Re_V / (2.193 · Re_θc)
F_onset2 = min(max(F_onset1, F_onset1^4), 2.0)
F_onset3 = max(1 - (R_T/2.5)^3, 0)
F_onset = max(F_onset2 - F_onset3, 0)
```

**sLM2015 form (Eq. 46-47, smooth):**
```
F_onset1 = √((Re_S / (2.6·Re_θc))² + R_T²)     (Eq. 46)
F_onset = (tanh(6·(F_onset1 - 1.35)) + 1) / 2    (Eq. 47)
```

where:
- Re_V or Re_S = ρ·S·y²/μ (strain-rate Reynolds number, Eq. 30)
- Re_θc = critical Reynolds number from correlation (Eq. 38 or 51)
- R_T = μ_t/μ (eddy viscosity ratio)

### 4.4 F_turb — Destruction control

**LM2015 form (Eq. 35):**
```
F_turb = e^(-(R_T/4)^4)
```

**sLM2015 form (Eq. 48):**
```
F_turb = (1 - F_onset) · e^(-R_T)
```

**Behavior comparison:**
- LM2015: F_turb → 0 only when R_T > ~4 (more forgiving)
- sLM2015: F_turb → 0 as soon as F_onset activates OR R_T > ~5 (more aggressive)

### 4.5 F_length — Controls P_γ growth rate

**LM2015 form (Eq. 37, piecewise polynomial):**
```
Re̅θt < 400:  F_length = 398.189e-1 - ...  (long polynomial)
400 ≤ Re̅θt < 596: F_length = ...
596 ≤ Re̅θt < 1200: F_length = ...
Re̅θt ≥ 1200: F_length = 0.01
```

**sLM2015 form (Eq. 49-50, generalized logistic):**
```
F_length = F_sub + (F_sup - F_sub) / (1 + e^(-r·(Re̅θt - b)))^(1/ν_gl)
```
with parameters from Eq. 50.

### 4.6 Re_θc — Critical Reynolds number

**LM2015 form (Eq. 38, piecewise):**
```
Re̅θt ≤ 1870: Re_θc = Re̅θt - (396.035e-2 - ...)  (polynomial)
Re̅θt > 1870: Re_θc = Re̅θt - 593.11 - ...
```

**sLM2015 form (Eq. 51, sinusoidal):**
```
Re_θc = 0.67·Re̅θt + 24·sin(Re̅θt/240 + 0.5) + 14
```

---

## 5. Crossflow Source D_scf (Eq. 15-26) — OPTIONAL

For swept wings. Adds a source to the Re̅θt equation:

```
P_θt_total = P_θt + D_scf
```

where D_scf involves:
- Helicity-based crossflow Reynolds number Re_scf (Eq. 22-26)
- Roughness height h
- Streamwise vorticity decomposition

**Skip for 2D cases.** D_scf ≡ 0 when no crossflow exists.

---

## 6. Smooth Max/Min Utility φ_p (Algorithm 1)

The sLM2015 model replaces all piecewise max/min with:

```
φ_p(g1, g2) = (g1 + g2)/2 + log(1 + e^(p·(g2-g1))) / p - log(2)/p
```

- p > 0 → smooth approximation to max(g1, g2)
- p < 0 → smooth approximation to min(g1, g2)
- |p| controls sharpness: |p| = 300 used throughout

**Algorithm 1 — Overflow-safe implementation:**

```fortran
function smooth_max_min(g1, g2, p) result(phi)
    real :: g1, g2, p, phi
    real :: diff, absdiff
    
    diff = p * (g2 - g1)
    absdiff = abs(diff)
    
    if (absdiff < 1e-15) then
        ! Proximity: g1 ≈ g2, return average + correction
        phi = (g1 + g2) / 2.0 + log(2.0) / abs(p) - log(2.0) / p
    else if (diff > 20.0) then
        ! g2 >> g1 (for p>0), max ≈ g2
        phi = g2 - log(2.0) / p
    else if (diff < -20.0) then
        ! g1 >> g2 (for p>0), max ≈ g1
        phi = g1 - log(2.0) / p
    else
        ! General case (no overflow)
        phi = (g1 + g2) / 2.0 + log(1.0 + exp(diff)) / p - log(2.0) / p
    end if
end function
```

**Key**: The `g1 + log(1+exp(p(g2-g1)))/p` rewrite avoids the naïve `log(exp(p·g1)+exp(p·g2))/p` which overflows when p·g > 700.

---

## 7. Numerical Robustness (Section IV)

### 7.1 Source-term time-step restriction (Eq. 59)

```
Δt_local = min(Δt_CFL, 0.9 / max(λ_max, ε))
```

where λ_max is the largest positive eigenvalue of the 3×3 source-term Jacobian:

```
J = ∂[R_ν̃, R_γ, R_θt]ᵀ / ∂[ν̃, γ, Re̅θt]
```

**Diagonal entries** (always needed):
```
J(1,1) = ∂R_ν̃/∂ν̃   (SA source Jacobian, existing)
J(2,2) = ∂R_γ/∂γ     (from P_γ and E_γ differentiation)
J(3,3) = ∂R_θt/∂Re̅θt (from P_θt differentiation)
```

**Off-diagonal entries** (needed for coupled mode):
```
J(1,2) = ∂R_ν̃/∂γ     (from γ·P_ν̃ coupling)
J(2,1) = ∂R_γ/∂ν̃     (from R_T = μ_t/μ in F_turb, E_γ)
J(1,3) = ∂R_ν̃/∂Re̅θt  (usually ~0, indirect via F_onset → γ)
J(3,1) = ∂R_θt/∂ν̃    (usually ~0)
J(2,3) = ∂R_γ/∂Re̅θt   (from Re_θc(Re̅θt) in F_onset)
J(3,2) = ∂R_θt/∂γ     (from F_θt(γ) in P_θt)
```

**Eigenvalue computation**: 3×3 matrix — use LAPACK `dgeev` or closed-form cubic.

**Deactivation switch** (Section IV.B): After 5 successive inexact-Newton iterations without backtracking and R_d > 1e-5, deactivate source-term restriction. Reactivate if backtracking triggers or residual rises.

### 7.2 Solution-update damping (Algorithm 2)

After computing ΔQ from the linear solve:

```
θ_fac = 0.99 (damping factor)

do iter = 1, 40
    if (γ + damping·Δγ < γ_lo) damping *= θ_fac
    if (γ + damping·Δγ > γ_hi) damping *= θ_fac
    if (Re̅θt + damping·ΔRe̅θt < Re̅θt_lo) damping *= θ_fac
end do

ΔQ *= damping
```

Bounds: γ ∈ [1e-10, 2.0], Re̅θt ∈ [20, ∞)

### 7.3 Row/column scaling (Section IV.C)

Scale residuals and state vector to bring all equations to similar magnitude:

```
Scaling factors: (ν̃_max, γ_max, Re̅θt_max) = (1e3, 10, 1e4)
```

Applied to both the linear system (preconditioner) and residual monitoring.

### 7.4 Coupling strategies (Section IV.B)

Three options, selectable at runtime:

1. **Decoupled**: Mean flow → SA → γ,Re̅θt separately (DADI or sub-ANK)
2. **Segregated**: Mean flow → [SA + γ + Re̅θt as one 3×3 block]
3. **Fully coupled**: All 8 equations (5 NS + ν̃ + γ + Re̅θt) in one Newton system

Paper recommends fully coupled for deep convergence. Decoupled is easiest to stabilize.

---

## 8. Boundary Conditions

### Wall
- γ = 0 (Dirichlet — laminar at surface)
- Re̅θt: zero normal gradient (∂Re̅θt/∂n = 0)

### Farfield
- γ = 1 (freestream is "fully turbulent" by convention — transition happens inside BL)
- Re̅θt = Re_θt(Tu_∞, λ_θ=0) from empirical correlation (Eq. 8 with freestream Tu)

### Symmetry / Periodic
- Standard zero-gradient / matched

---

## 9. Constants Summary

| Symbol | Value | Location |
|--------|-------|----------|
| ca1 | 2.0 | P_γ production |
| ca2 | 0.06 | E_γ destruction |
| ce1 | 1.0 | P_γ (1-ce1·γ) term |
| ce2 | 50.0 | E_γ (ce2·γ-1) term |
| c_θt | 0.03 | P_θt coefficient |
| σ_θt | 2.0 | Re̅θt diffusion |
| σ_f | 1.0 | γ diffusion |
| θ_fac | 0.99 | Solution damping |
| p | 300 | φ_p sharpness (sLM2015) |

---

## 10. SA Modification Detail

The ONLY modification to SA is multiplying production by γ:

**Original SA production:**
```
P_ν̃ = cb1·(1-ft2)·S̃·ν̃ + (1/σ)·(cb2·(∇ν̃)² - (cb2+κ²)·(ν̃/d)²·ft2·...)
```

**Modified SA production (Eq. 41):**
```
P̃_ν̃ = γ · [cb1·(1-ft2)·S̃·ν̃]                    ← γ on first production term
      + γ · [(1/σ)·(-cb1·ft2/κ²)·ν̃·fv2·...]     ← γ on fv2 correction (part of production)
      + [(cw1·fw - cb1·ft2/κ²)·(ν̃/d)²]           ← NO γ on destruction
```

γ multiplies ONLY production terms, NOT destruction. This is critical.

---

## 11. Convection Scheme

Paper Section IV.A explicitly recommends **first-order upwind** for both γ and Re̅θt transport equations. This is important because:
- Transition variables have sharp gradients at onset
- Higher-order schemes can cause oscillations in γ near transition
- SA itself uses first-order upwind in many implementations

Implementation: Add runtime option `transitionFirstOrderUpwind` (default True).

---

## 12. Key Figures for Validation

| Figure | Case | What it shows |
|--------|------|---------------|
| Fig. 6 | NLF0416, M=0.1, Re=4e6, Tu=0.15% | Drag polar, lift curve, transition location |
| Fig. 7 | NLF0416 | Convergence history (sLM2015 vs LM2015) |
| Fig. 8-9 | S809, M=0.1, Re=2e6, Tu=0.07% | AoA sweep, LSB detection |
| Fig. 10-11 | NLF2-0415, swept | Crossflow transition |
| Fig. 12-14 | Sickle Wing, 3D | Skin friction, Cp, convergence |
