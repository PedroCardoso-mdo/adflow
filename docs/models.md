# Transition models

What each model computes as implemented, with the equation numbers of the source papers. The code is the
reference: if this page and the code disagree, the code is right and this page needs fixing.
All quantities are in ADflow's pressure–density non-dimensionalisation (see [developer.md](developer.md));
the papers' explicit Reynolds-number factors do not appear in the code.

Notation: ν̃ SA working variable, γ intermittency, R̅eθt transported transition-onset momentum-thickness
Reynolds number, d wall distance, S strain-rate magnitude, Ω vorticity magnitude, U velocity magnitude,
R_T = ν_t/ν with ν_t = ν̃·f_v1, Tu free-stream turbulence intensity.

## 1. GR: `turbulenceModel="SA-noft2-Gamma-Retheta"`

The smooth γ–R̅eθt model coupled to SA without the ft2 term, called SA-sLM2015 in
Piotrowski and Zingg, "Smooth Local Correlation-Based Transition Model for the Spalart–Allmaras Turbulence
Model", AIAA Journal, 2020, doi:10.2514/1.J059784 ("P&Z" below). Code: `src/turbulence/saGammaRetheta.F90`.

### Equations implemented

| Piece                     | As implemented                                                                                                                                | P&Z equation    |
| ------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------- | --------------- |
| Coupling to SA            | SA production multiplied by γ (clipped to [1e-10, 1]); SA destruction and diffusion unchanged; ft2 off                                        | 41              |
| γ production              | P_γ = c_a1·F_length·F_onset·Ω_lim·√γ·(1 − c_e1·γ)                                                                                             | —               |
| γ destruction             | E_γ = c_a2·F_turb·Ω_lim·γ·(c_e2·γ − 1)                                                                                                        | —               |
| Onset                     | F_onset1 = √[(Re_S/(2.6·Re_θc))² + R_T²], Re_S = ρ·d²·S/μ; F_onset = ½·[tanh(6·(F_onset1 − 1.35)) + 1]                                        | 5 (Re_S)        |
| Relaminarisation switch   | F_turb = (1 − F_onset)·exp(−R_T)                                                                                                              | —               |
| F_length(R̅eθt)            | smooth correlation, routine `flengthCorrelation`                                                                                              | 49–50           |
| Re_θc(R̅eθt)               | smooth correlation, routine `rethetacCorrelation`                                                                                             | 51              |
| Vorticity limiter         | Ω_lim = smooth min(Ω, Ω_cap); Ω_cap = U_∞·√(U_∞/(ν_∞·l))/20, i.e. the paper's M·√(M·Re)/20 with the reference length l explicit               | 52–53           |
| R̅eθt source               | P_θt = (c_θt/t)·(Re_θt − R̅eθt)·(1 − F_θt), time scale t = 500·ν/U²                                                                            | 7               |
| Re_θt(Tu, λ_θ)            | Langtry–Menter correlation with the smooth pressure-gradient function F(λ_θ), routine `reThetaTCorrelation`; λ_θ = (θ²/ν)·dU/ds, θ = R̅eθt·ν/U | 4, 10–11, 54–57 |
| Shielding F_θt            | F_θt = F_wake·exp(−(d/δ)⁴), δ = 50·d·Ω·δ_BL/U, δ_BL = 7.5·θ, F_wake = exp(−Re_S/1e6)                                                          | 4               |
| Crossflow (optional)      | D_scf = (c_θt/t)·c_cf·min(Re_scf − R̅eθt, 0)·F_θt added to the R̅eθt source                                                                     | 15              |
| Crossflow Reynolds number | Re_scf = −35.088·ln(h/θ) + 319.51 + 6200·ΔH⁺ + 50000·(ΔH⁺)² − 75·tanh(ΔH⁻/0.0125); ΔH± = max(±(0.1066 − H_cf·(1 + min(R_T, 0.4))), 0)         | 17              |
| Helicity parameter        | H_cf = d·abs(Û·ω)/U, with the unit velocity vector Û and the vorticity vector ω                                                               | 24–26           |
| Smooth min / max          | φ_p(a, b) with p = ±300 and a proximity switch, routine `smoothMinMax`; replaces every min/max of the original model                          | 43–45           |
| Convection of γ and R̅eθt  | first-order upwind (`transitionFirstOrderUpwind`)                                                                                             | Sec. IV.A       |
| Diffusion                 | σ_f = 1.0 (γ), σ_θt = 2.0 (R̅eθt); routine `Viscous`                                                                                           | —               |

A dash in the last column means that the code carries no equation number for that piece; none is asserted here.

### Constants (`src/modules/paramTurb.F90`)

| Constant                    | Value | Constant                      | Value           |
| --------------------------- | ----- | ----------------------------- | --------------- |
| c_a1                        | 2.0   | c_θt                          | 0.03            |
| c_a2                        | 0.06  | c_cf (crossflow)              | 0.6             |
| c_e1                        | 1.0   | H_cf reference                | 0.1066          |
| c_e2                        | 50    | R_T cap in the crossflow term | 0.4             |
| σ_f                         | 1.0   | vorticity-limiter divisor     | 20              |
| σ_θt                        | 2.0   | smooth min/max exponent p     | ±300            |
| onset: 2.6, k = 6, s = 1.35 | fixed | γ bounds / R̅eθt lower bound   | [1e-10, 2] / 20 |

### Safeguards that are not written in the paper

| Safeguard                                             | Value               | Why it is there                                                                                                                                                     |
| ----------------------------------------------------- | ------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| λ_θ limited (smoothly) to [−0.1, 0.1]                 | compile-time switch | Validity range of the Langtry–Menter correlation. Without it S809 stalls near a residual of 1e-4 and cd moves away from the paper (campaign `12_paper_conditions`). |
| Tu floor in the Re_θt correlation                     | 0.027 %             | Langtry–Menter lower bound; keeps the correlation finite for very low Tu.                                                                                           |
| γ seen by SA production clipped to [1e-10, 1 + 1e-10] | —                   | Divergence guard. The upper bound is padded so that γ = 1 exactly is not a branch point (forward AD and complex step disagreed there).                              |
| Non-negative diagonal of the implicit source Jacobian | —                   | Left-hand side only (not differentiated); the converged solution is unaffected.                                                                                     |

### Boundary, free-stream and initial conditions

| Location                   | ν̃                      | γ                    | R̅eθt                                      |
| -------------------------- | ---------------------- | -------------------- | ----------------------------------------- |
| Free stream / far field    | from `eddyVisInfRatio` | 1                    | Re_θt(Tu_∞, λ_θ = 0) from the correlation |
| Wall                       | 0                      | zero normal gradient | zero normal gradient                      |
| Initial field (cold start) | free-stream value      | 0.02 everywhere      | free-stream value                         |

The interior γ = 0.02 keeps SA production off until the onset function activates it. The value has not
been tuned.

### Free-stream values that were validated

`turbIntensityInf` is a fraction (0.0015 = 0.15 %). `eddyVisInfRatio` is a boundary condition of the model,
not a solver setting: all three models sit on a plateau below about 1e-6 and lose the laminar region at the
ADflow default 0.009 (NLF(1)-0416, α 4°: GR 65.58 → 75.51 counts; campaign `24_minmodel_2d`, sweep 24-D).

| Case                         | `turbIntensityInf` | `eddyVisInfRatio` (GR)      | Note                                                                                               |
| ---------------------------- | ------------------ | --------------------------- | -------------------------------------------------------------------------------------------------- |
| NLF(1)-0416, Re 4e6          | 0.0015             | 2.79e-7 (scalar and matrix) | 2.79e-7 corresponds to χ_∞ = 0.1 of the paper. A matrix series at 5e-7 also agrees with the paper. |
| S809, Re 2e6                 | 0.0007             | 2.79e-7                     | At 1e-7, α 2° falls into a laminar-bubble limit cycle (cd oscillating 58–70 counts).               |
| NACA 0012, Re 2.92e6         | 0.005 and 0.003    | 2.79e-7                     |                                                                                                    |
| NLF(2)-0415 swept, crossflow | 0.002              | 2.79e-7                     |                                                                                                    |
| Sickle wing, crossflow       | 0.0017             | 5e-7                        | Difference from 2.79e-7 not reconciled.                                                            |
| AR5 wing                     | 0.0025             | 0.009 (ADflow default)      | Run before the sensitivity was known; not the recommended value.                                   |
| Flat plate, Re 2.5e6         | 0.0010             | 5e-7                        |                                                                                                    |

A steady GR solution does not always exist: on a strongly laminar optimised airfoil the solver reached
8 orders at Tu ≥ 0.25 % and stalled at 2.8e-6 … 7e-6 for Tu 0.10–0.20 % (campaign `24_minmodel_2d`, 24-H).

### Model options

| Option                       | Default  | Effect                                                                                                                                                                |
| ---------------------------- | -------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `transitionCrossflow`        | `False`  | Adds D_scf. Zero in 2D. Off by default because the term is only validated with matrix dissipation (with scalar dissipation the swept wing shows no crossflow front).  |
| `transitionRoughnessHeight`  | `3.3e-6` | Roughness h in Re_scf, in mesh length units. Rescale it if the mesh chord is not the physical chord (sickle wing at chord 1.0 m: 1.2225e-5).                          |
| `transitionRefLength`        | `-1`     | Length l in Ω_cap. Default: the AeroProblem `chordRef`. The cap scales as 1/√l; it is a physical calibration length (root chord in the paper), not a unit conversion. |
| `transitionFirstOrderUpwind` | `True`   | First-order convection of γ and R̅eθt. Leave on.                                                                                                                       |
| `transitionUseApproxSA`      | `True`   | Approximate-SA simplification of the ν̃ production in the first-order ANK stage. Solver aid; the converged residual is the full one.                                   |
| `transitionArtVisc`          | `0`      | Artificial viscosity on γ. Changes the solution; no validated case. Leave at 0.                                                                                       |

### Known limitations

| Limitation                            | Detail                                                                                                                               |
| ------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------ |
| No multigrid, no wall functions       | `MGCycle="sg"` is required in practice; `useWallFunctions=True` is an error.                                                         |
| Block residual path only              | `useBlockettes` is forced off (performance cost, not a correctness issue).                                                           |
| Crossflow needs matrix dissipation    | Scalar dissipation under-predicts the helicity; see [validation.md](validation.md).                                                  |
| Chordwise mesh resolution             | On the AR5 wing the transition front converges only from 397 points around the section; 181 points give x/c 0.33 instead of 0.62.    |
| Sensitivity to dissipation and to ν̃_∞ | Scalar and matrix dissipation give different drag on the same mesh (S809 α 4°: 66.42 against 63.86 counts on the 193 224-cell mesh). |
| Deep Newton convergence               | Below a relative residual of about 5e-9 the NK linear solve saturates on 3D wings; the limit is the preconditioner.                  |
| Rotating frames                       | The source terms are written in the relative frame, but no rotating case has been run ([developer.md](developer.md)).                |
| F_wake and F_θt                       | Implemented as in the table above. They were not re-checked term by term against the paper for this page.                            |

## 2. SA-γ: `turbulenceModel="SA-noft2-Gamma"`

One transport equation for γ on top of SA, no R̅eθt equation. Code: `src/turbulence/saGamma.F90`; the header
of that file numbers the equations 1–10 and that numbering is used below.

Sources: Menter, Smirnov, Liu and Avancha, "A One-Equation Local Correlation-Based Transition Model",
Flow, Turbulence and Combustion, 2015 (γ equation); Nichols 2019 and Lee, PhD dissertation, University of
Maryland, 2021 (coupling to SA); Lee and Baeder 2021 (constant free-stream Tu); Colonia et al. 2016, via
Jung et al. 2022 (Re_θc constants); P&Z 2020 (smooth min/max). **Status: final form** =
the first onset form (smooth min/max, P&Z F_turb) with artificial viscosity on γ (`sgammaArtVisc`); the tanh onset and the Lee F_turb stay available as options.

### Equations implemented

| No. | As implemented                                                                                                                                     |
| --- | -------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1   | Dγ/Dt = P_γ − E_γ + ∇·[(ν + ν_t/σ_γ)∇γ], σ_γ = 1                                                                                                   |
| 2   | P_γ = F_length·S·γ·(1 − γ)·F_onset, F_length = 100                                                                                                 |
| 3   | E_γ = c_a2·Ω·γ·F_turb·(c_e2·γ − 1), c_a2 = 0.06, c_e2 = 50                                                                                         |
| 4   | F_onset1 = Re_v/(2.2·Re_θc), Re_v = ρ·d²·S/μ; F_onset = φ⁺(φ⁻(F_onset1, 2) − φ⁺(1 − (R_T/3.5)³, 0), 0). Alternative: tanh form (`sgammaOnsetTanh`) |
| 5   | F_turb = (1 − F_onset,normalised)·exp(−R_T). Alternative: exp(−(R_T/2)⁴) (`sgammaFturbLee`)                                                        |
| 7   | Re_θc = 163 + 1002.25·exp(−Tu·F_PG), Tu in percent; F_PG = 1 + φ⁻(14.68·φ⁺(λ, 0), 0.5) + φ⁻(−7.34·φ⁻(λ, 0), 2.0)                                   |
| 8   | λ = 7.57e-3·(d²/ν)·dU/ds + 0.0128, limited smoothly to [−1, 1]                                                                                     |
| 9   | Tu = `turbIntensityInf`, constant in the field                                                                                                     |
| 10  | γ_s = φ⁺((φ⁻(γ, 1) − 1/c_e2)/(1 − 1/c_e2), 0); SA production × γ_s; SA destruction × φ⁺(γ_s, 0.1) (`sgammaCoupleDestruction`)                      |

φ⁺ / φ⁻ are the smooth max / min of P&Z with exponent 300 (`sgammaFPGSmoothP`, `sgammaOnsetSmoothP`).
Convection of γ is first-order upwind; wall, free-stream and initial conditions for ν̃ and γ are those of GR.

### Model options

| Option                                                        | Default       | Effect                                                                                                                                                                                                 |
| ------------------------------------------------------------- | ------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `sgammaCoupleDestruction`                                     | `True`        | Eq. 10 on destruction as well. Used in every reported run.                                                                                                                                             |
| `sgammaOnsetTanh`                                             | `False`       | F_onset = A·½·[tanh(k·(F₁ − s)) + 1] with F₁ = √[(Re_v/(2.2·Re_θc))² + (R_T/c)²].                                                                                                                      |
| `sgammaTanhK`, `sgammaTanhS`, `sgammaTanhAmp`, `sgammaTanhRT` | 6, 1.35, 1, 1 | k, s, A, c of the tanh onset (defaults are the P&Z values). The smoothed variant tested was A 2.5, k 6, s 1.45, c 4 together with `sgammaFturbLee`.                                                    |
| `sgammaFturbLee`                                              | `False`       | Alternative F_turb.                                                                                                                                                                                    |
| `sgammaFPGSmoothP`, `sgammaOnsetSmoothP`                      | 300, 300      | Sharpness of the smooth min/max in F_PG and in the default onset.                                                                                                                                      |
| `sgammaArtVisc`                                               | `0`           | Artificial viscosity on γ (form after Ricci et al. 2026, modified). c = 10 closes 3 of 4 limit-cycle cases and lowers cd by 1.4–3.7 %; c = 100 raises cd by 3.5–4.5 %. Not part of the baseline model. |

### Free-stream values used

`eddyVisInfRatio` 1e-7 (NACA 0012 runs) and 2.79e-7 (NLF(1)-0416 and S809 polars); `turbIntensityInf` as for GR.

### Known limitations

| Limitation            | Detail                                                                                                                        |
| --------------------- | ----------------------------------------------------------------------------------------------------------------------------- |
| Calibration not final | Drag is 4–8 % above GR at low α on the validation airfoils (NLF(1)-0416 α 4°: 69.28 against 65.58 counts).                    |
| Convergence           | 7 of 23 scalar polar points stop at 3.2–6.0 orders; with matrix dissipation S809 α 8° stays at 1.3e-5.                        |
| No R̅eθt history       | Onset is local. Shapes optimised with SA-γ lose much of their gain when re-analysed with GR ([validation.md](validation.md)). |
| No crossflow term     | —                                                                                                                             |

## 3. SA-BCM: `useSABCM=True`

Algebraic intermittency that multiplies the SA production term; no extra equation. Reference: Mura and
Cakmakcioglu, "A Revised One-Equation Transitional Model for External Aerodynamics — Part I: Theory,
Validation and Base Cases", AIAA Paper 2020-2714 (Appendix formulation). Code: `saSource` in
`src/turbulence/sa.F90`; with `useSABCM=False` the routine is upstream SA.

### Equations implemented

| Piece                               | As implemented                                                                                  |
| ----------------------------------- | ----------------------------------------------------------------------------------------------- |
| Critical Reynolds number            | Re_θc = 803.73·(Tu + 0.6067)^(−1.027), Tu in percent (= 100·`turbIntensityInf`)                 |
| Local Reynolds number               | Re_θ = Re_v/2.193, Re_v = ρ·Ω·d²/μ                                                              |
| Term1                               | (Re_θ − Re_θc)/(χ₁·Re_θc)                                                                       |
| Term2                               | f_v1·χ/χ₂ = (ν_t/ν)/χ₂                                                                          |
| Original form (`SABCMSmooth=False`) | γ = 1 − exp(−(√max(Term1, 0) + √Term2))                                                         |
| Smooth form (`SABCMSmooth=True`)    | Term1 → KS-smoothed max(Term1, 0) with parameter ρ_KS; γ = ½·[1 + tanh((Term1 + Term2 − s₀)/w)] |
| Coupling                            | γ multiplies the SA production terms; ft2 = 0; destruction unchanged                            |

The vorticity is in ADflow units and ν̃/ν is a ratio, so no Reynolds-number factor appears.
γ is not stored during the solve; the `intermittency` output variable recomputes it from the state.

### Constants (all are options)

| Option            | Default | Symbol        |
| ----------------- | ------- | ------------- |
| `SABCMChi1`       | 0.002   | χ₁            |
| `SABCMChi2`       | 0.02    | χ₂            |
| `SABCMRho`        | 50      | ρ_KS (smooth) |
| `SABCMTanhCenter` | 0.5     | s₀ (smooth)   |
| `SABCMTanhWidth`  | 0.08    | w (smooth)    |

### Free-stream values used

| Campaign                                 | `eddyVisInfRatio` | Note                                                                                                             |
| ---------------------------------------- | ----------------- | ---------------------------------------------------------------------------------------------------------------- |
| `24_minmodel_2d` (polars, optimisations) | 1e-9              | cd keeps decreasing slowly as ν̃_∞ is lowered (NLF(1)-0416 α 4°: 65.90 → 65.68 → 64.87 counts, the last at 1e-9). |
| SA-BCM article, `15_sabcm_polars`        | 1e-7              |                                                                                                                  |
| ADflow default 0.009, or 0.21            | —                 | Fully turbulent: no transition (NLF(1)-0416 α 4°: 112.7 counts).                                                 |

`doc/options.yaml` recommends 1.4e-10 to 1.1e-9.

### Requirements and limitations

| Item                                     | Detail                                                                                                                                                |
| ---------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
| Requirements                             | RANS, `turbulenceModel="SA"`, `useApproxWallDistance=True` (otherwise the run stops); `useft2SA` is ignored.                                          |
| Local, algebraic onset                   | No transport of onset information. Optimised shapes can exploit this: the gain predicted by SA-BCM is often not confirmed by GR.                      |
| Newton–Krylov                            | In 2D, NK with a fixed step converged to a different state (cd 86 instead of 65.8 counts); the validated recipes use ANK only.                        |
| Dissipation sensitivity                  | One optimised airfoil gave 46.13 counts with `vis4` 0.0156 and 84.64 with 0.01.                                                                       |
| Smooth against original                  | The smooth form is not established as more accurate or better differentiable than the original; see the retraction in [validation.md](validation.md). |
| Matrix dissipation with `epsShear` 0.025 | cd about twice the transitional level (older SA-BCM code, campaign `15_sabcm_polars`); use `epsShear=0`.                                              |
