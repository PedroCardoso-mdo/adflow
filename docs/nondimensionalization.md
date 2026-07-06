# Non-Dimensionalization of the Governing Equations in ADflow

> **Key takeaway:** ADflow uses a **pressure–density ("p-ρ") non-dimensionalization,
> NOT a velocity-based one.** Density and pressure normalize to 1, but the
> free-stream **velocity normalizes to `M·√γ`, not 1**, because the velocity scale
> is `uRef = √(pRef/rhoRef)` (an acoustic-type scale), not the free-stream velocity
> U∞. This is the single most common point of confusion when mapping textbook /
> paper equations (which usually assume U∞ = 1 and μ_ref = 1/Re) onto ADflow.

Authoritative source: subroutine `referenceState`
(`src/initFlow/initializeFlow.F90:10-190`) and the variable comments in
`src/modules/flowVarRefState.F90`. This document just summarizes them.

---

## 1. Reference values (= free-stream dimensional values)

The reference length is **1.0 m** (coordinates are converted to meters), so it
never appears explicitly. For external flow the reference state is simply the
free stream:

| Reference | Set to | Code |
|-----------|--------|------|
| `pRef`   | `pInfDim`   (free-stream pressure, Pa)      | `initializeFlow.F90:57` |
| `rhoRef` | `rhoInfDim` (free-stream density, kg/m³)    | `initializeFlow.F90:59` |
| `TRef`   | `TInfDim`   (free-stream temperature, K)    | `initializeFlow.F90:58` |

> These *could* differ from free-stream for internal-flow cases, but the code
> currently always uses the actual free-stream values.

## 2. Derived scales (with reference length = 1)

| Scale | Formula | Code |
|-------|---------|------|
| Velocity  | `uRef = √(pRef/rhoRef) = √hRef`  | `initializeFlow.F90:76` |
| Enthalpy  | `hRef = pRef/rhoRef`             | `initializeFlow.F90:75` |
| Viscosity | `muRef = √(pRef·rhoRef)`         | `initializeFlow.F90:68` |
| Time      | `timeRef = √(rhoRef/pRef)`       | `initializeFlow.F90:74` |

**Note the velocity scale**: `uRef = √(pRef/rhoRef)`. This is *not* U∞ — it is an
acoustic-type scale (≈ a∞/√γ). Everything downstream follows from this choice.

## 3. Resulting non-dimensional free stream

| Quantity | Value | Code |
|----------|-------|------|
| `rhoInf` | `rhoInfDim/rhoRef = 1`                 | `initializeFlow.F90:82` |
| `pInf`   | `pInfDim/pRef = 1`                     | `initializeFlow.F90:81` |
| `uInf`   | `Mach·√(γ·pInf/rhoInf) = M·√γ` (**not 1**) | `initializeFlow.F90:83` |
| `muInf`  | `muInfDim/muRef`                       | `initializeFlow.F90:85` |
| `nuInf`  | `muInf/rhoInf`                         | `initializeFlow.F90:130` |
| `RGas`   | `RGasDim·rhoRef·TRef/pRef`             | `initializeFlow.F90:84` |

Self-consistency check: the non-dimensional speed of sound is
`a∞ = √(γ·pInf/rhoInf) = √γ`, so `uInf = M·a∞ = M·√γ`. ✔

## 4. Reynolds number

Re is **not** injected by setting `μ_ref = 1/Re`. It appears only through the
non-dimensional viscosity:

```
muInf = muInfDim / √(pRef·rhoRef)      (with L = 1)
Re_nondim = rhoInf·uInf·L / muInf = uInf / muInf
```

State-vector viscosities are stored as **ratios**:
- `rlv(i,j,k) = μ/μ_∞`   (laminar viscosity ratio, dimensionless)
- `rev(i,j,k) = μ_t/μ_∞` (eddy viscosity ratio, dimensionless)

## 5. Why this matters for the SA-γ-Re̅θt transition model

The Piotrowski & Zingg (2020) correlations assume a **velocity-based** nondim
(U∞ = 1, μ_ref = 1/Re). Because ADflow instead uses p-ρ scaling and stores
`rlv = μ/μ_∞`, quantities that mix velocity, viscosity, and Reynolds number must
be rewritten. Concrete example — the transition time scale (see
[`paper-reference.md`](paper-reference.md) §2, and the note at
`paper-reference.md:59`):

```
Paper (velocity-based):   t = 500·μ̃/(ρ·U²)      with μ̃ = 1/Re
ADflow (p-ρ, rlv ratio):  t = 500·rlv/(ρ·U²·Re)
```

**Rule of thumb when porting an equation:** velocities are in units of `uRef`
(so free stream = M·√γ), viscosities are ratios to μ_∞, and any explicit `1/Re`
in the paper has to be reintroduced by hand — it is *not* absorbed into ADflow's
viscosity like it is in a velocity-normalized code.

---

**See also:** [`architecture.md`](architecture.md) (state-vector layout,
`rlv`/`rev`/`timeRef` usage) and [`paper-reference.md`](paper-reference.md)
(equations that consume these non-dimensional quantities).
