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

## 5. The paper → ADflow non-dimensionalization shift

**The paper's equations are already non-dimensional — but by a *different*
scheme than ADflow's**, so terms must be shifted, not copied verbatim.

Piotrowski & Zingg (2020, p. 7) non-dimensionalize by freestream + **sound
speed** (OVERFLOW-style):

```
Re = ρ∞ · a∞ · l / μ∞          (a∞ = freestream sound speed, l = root chord)
```

Because that scheme scales by a∞ and l, the paper's equations carry **explicit
`Re` / `1/Re` factors** wherever a physical Reynolds number is involved — the
time scale (Eq. 7), θ_BL (Eq. 4), Re_S (Eq. 5), λ_θ (Eq. 10), Re_scf (Eq. 17).

ADflow's p-ρ scheme (Sections 1–3) instead absorbs the Reynolds number into the
**non-dimensional viscosity** `muInf = muInfDim/muRef`. Two facts make the
paper's `Re` factors collapse to 1 in the code:

- **`L_ref = 1 m`** — coordinates are metres, so a non-dim length numerically
  equals the physical length in metres.
- **`ρRef·uRef/μRef = 1`** — since `muRef = √(pRef·rhoRef)` and
  `uRef = √(pRef/rhoRef)`.

So ν-based quantities evaluate **directly to physical values**, and **the paper's
explicit `Re`/`1/Re` factors are dropped in code** — the "`Re = 1` implicit"
convention (see the `timeScale` comment in `saGammaRetheta.F90`; `sa.F90` does
the same).

**Timescale example (Eq. 7):**

```
Paper (a∞-based):   t = 500·μ/(ρ·U²) · (1/Re)        ← explicit 1/Re
ADflow code:        timeScale = 500·nu/velMag2        ← NO Re factor
                              = 500·rlv/(ρ·|V|²)        (Re = 1 implicit)
```

> **Correction (2026-07-06):** earlier revisions of this file (and the
> since-retired distilled reference) called the paper "velocity-based
> (U∞ = 1)" and wrote the ADflow timescale with an extra `·Re`. Both were
> wrong: the paper is a∞-based, and the code carries **no** explicit Re factor.

**Rule of thumb when porting a paper equation:** keep velocities in `uRef` units
(freestream = M·√γ) and viscosities as ratios to μ_∞, and **drop** the paper's
explicit `Re`/`1/Re` — it is already absorbed by ADflow's non-dim viscosity.
Lengths then come out physical (metres) and Reynolds-number groupings come out as
physical Reynolds numbers, which is exactly what the empirical correlations expect.

> **Exception (2026-07-07) — vorticity limiter, Eqs. 52–53.** The "drop Re" rule
> only holds where the paper's `Re` is a unit-conversion artifact of a
> scale-invariant equation. In the vorticity limiter `M·√(M·Re)/20` the √Re is a
> **physical calibration scale** (BL wall vorticity ∝ √(ρU³/(μ·l))) — the paper
> itself says so (§IV, "a physical scaling that is independent of the
> nondimensionalization") and ties `l` to the **root chord**. Dropping Re here
> silently substitutes `l = 1 m`, changing the physical cap by √(chord).
> Dimensional test: a threshold on Ω (1/s) built from freestream ρ, U, μ alone
> *needs* a length — nothing else cancels it. Handled by the
> `transitionRefLength` option (auto = AeroProblem `chordRef`; see
> `architecture.md` Part 2). Findings log: `findings/D1_transitionRefLength.md`.
>
> **AD corollary (2026-07-07):** because `uInf`/`muInf` are the code's spelling
> of the paper's freestream M and Re, the cap makes the *residual* depend on
> freestream reference state. Consistency check of the identity:
> `uInf/muInf = ρ∞U∞/μ∞` (since `ρRef·uRef/muRef = 1`), so the code line
> `vortLim = uInf·√(uInf/(muInf·l))/20 = uInf·√Re_U∞/(20·l)` — exactly the
> paper's cap in ADflow vorticity units. For derivatives this means `uInf`,
> `muInf` must be declared active independents of `saGammaRetheta%Source` in
> `Makefile_tapenade`, else Tapenade emits `vortlimd = 0` and `dR/dMach`-type
> partials lose the limiter path in capped cells (state partials dR/dw are
> unaffected). Details + frozen-cap alternative: `adjoint-trace.md` header and
> `audits/adjoint_audit_2026-07-07.md` §3.

## 6. Crossflow (D_scf) term dimension status

The helicity-based crossflow source (Eq. 15–26) was added using the shift above.
Where each new quantity lands dimensionally in the code:

| Quantity | Paper form (Eq.) | Code (`saGammaRetheta.F90`) | Status in ADflow |
|---|---|---|---|
| θt (`thetaBL`) | Re̅θt·μ/(ρU)·(1/Re) (Eq. 4) | `reThetaTilde*nu/velMag` | physical momentum thickness — non-dim by L_ref=1 m ⇒ **value in metres** |
| h (`transitionRoughnessHeight`) | `h/θt` (Eq. 17) | input, default `3.3e-6` | physical roughness length — **must be in mesh units (metres)**; 3.3e-6 = 3.3 µm |
| H_cf (`hcf`) | d·Ω_sw/U (Eq. 26) | `yDist*abs(Û·Ω⃗)/velMag` | **dimensionless** (L_ref cancels; uses pure curl, frame-independent) |
| Re_scf (`reScf`) | correlation (Eq. 17) | `-35.088·ln(h/θt)+319.51+f(ΔH⁺)−f(ΔH⁻)` | **physical Reynolds number** (calibrated correlation) |
| D_scf (`dScf`) | (c_θt/t)·c_cf·min(Re_scf−Re̅θt,0)·F_θt (Eq. 15) | `(rsaGRcthetat/timeScale)*…` | **Re/time**, same units as `P_θt`; no explicit Re |

Since `thetaBL` is the physical θ (metres) and `h` is a physical length, `h/θt`
is the true physical ratio (reference-length-independent) and the correlation
constants apply directly; `min(Re_scf − Re̅θt, 0)` compares two physical Reynolds
numbers. No explicit paper `Re` factor appears in any crossflow term — per §5.

---

**See also:** [`architecture.md`](architecture.md) (state-vector layout,
`rlv`/`rev`/`timeRef` usage, crossflow options) and the full paper,
[`SA_GAMMA_RETHETHA_BASE/Piotrowski_Zingg_2020_SA-sLM2015_clean (1).md`](SA_GAMMA_RETHETHA_BASE/Piotrowski_Zingg_2020_SA-sLM2015_clean%20(1).md)
(equations that consume these non-dimensional quantities).
