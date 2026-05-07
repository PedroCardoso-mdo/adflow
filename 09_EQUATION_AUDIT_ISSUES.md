# Equation Audit Issues — SA-sLM2015 Implementation


## Issue 2: F_wake Wrong Formula

**Location**: `src/turbulence/saGammaRetheta.F90:545`

**Paper (Eq. 3 inner)**:
```
F_wake = e^(-(Re_ω / 1e5)^2)
```
where Re_ω = ρ·Ω·y²/μ (vorticity-based Reynolds number).

**Current Code**:
```fortran
fWake_val = exp(-reS_val / 1.0e6_realType)
```

**Problems**:
1. Uses strain-based `reS_val` instead of vorticity-based `Re_ω`
2. Exponent is `-x/1e6` instead of `-(x/1e5)²`

**Fix Plan**:
```fortran
! Compute vorticity-based Reynolds number Re_omega (not strain-based)
reOmega = w(i, j, k, irho) * yDist**2 * vortMag / max(rlv(i, j, k), xminn)

! Correct F_wake formula (Eq. 3)
fWake_val = exp(-(reOmega / 1.0e5_realType)**2)
```

**Note**: Also update `saGammaRethetaHelpers.F90:214-221` if used elsewhere.

---

## Issue 3: P_γ Has Extra √γ Factor

**Location**: `src/turbulence/saGammaRetheta.F90:503-505`

**Paper sLM2015 (Eq. 52)**:
```
P_γ = ca1 · F_length · Ω_lim · F_onset · (1 - ce1·γ)
```

**Paper LM2015 (Eq. 28)**:
```
P_γ = ca1 · F_length · S · [γ · F_onset]^0.5 · (1 - ce1·γ)
```

**Current Code**:
```fortran
pGamma = rsaGRca1 * fLength_val * fOnset * vortMagLim &
         * sqrt(max(gammaLocal, xminn)) &
         * (one - rsaGRce1 * gammaLocal)
```

**Problem**: Code has `sqrt(γ)` but uses `Fonset` directly (not `sqrt(γ·Fonset)`). This is neither LM2015 nor sLM2015 form.

**Fix Plan** (use sLM2015 for consistency):
```fortran
! sLM2015 form (Eq. 52) - no sqrt(gamma), linear in Fonset
pGamma = rsaGRca1 * fLength_val * fOnset * vortMagLim &
         * (one - rsaGRce1 * gammaLocal)
```

**Alternative** (LM2015 form if preferred):
```fortran
! LM2015 form (Eq. 28) - sqrt(gamma * Fonset)
pGamma = rsaGRca1 * fLength_val * strainMag &
         * sqrt(max(gammaLocal * fOnset, xminn)) &
         * (one - rsaGRce1 * gammaLocal)
```

---

## Issue 4: F_turb Uses LM2015 Instead of sLM2015

**Location**: `src/turbulence/saGammaRetheta.F90:499-500`

**Paper sLM2015 (Eq. 48)**:
```
F_turb = (1 - F_onset) · e^(-R_T)
```

**Paper LM2015 (Eq. 35)**:
```
F_turb = e^(-(R_T/4)^4)
```

**Current Code**:
```fortran
!fTurb_val = (one - fOnset) * exp(-rTurb)  ! sLM2015 - commented out
fTurb_val = exp(-(rTurb / 4.0_realType)**4)  ! LM2015 - active
```

**Problem**: Using LM2015 F_turb while other correlations (Fonset, Flength, Re_θc) use sLM2015. This inconsistency may affect convergence behavior.

**Fix Plan** (enable sLM2015):
```fortran
! sLM2015 form (Eq. 48) - couples destruction to onset
fTurb_val = (one - fOnset) * exp(-rTurb)
```

**Note**: The LM2015 form may have been chosen deliberately for stability. Consider making this a runtime option:
```fortran
if (useSmoothFturb) then
    fTurb_val = (one - fOnset) * exp(-rTurb)        ! sLM2015
else
    fTurb_val = exp(-(rTurb / 4.0_realType)**4)     ! LM2015
end if
```

---

## Implementation Order

1. **Issue 2 (F_wake)** — Straightforward fix, low risk
2. **Issue 1 (F_θt)** — Straightforward fix, low risk
3. **Issue 4 (F_turb)** — Consider runtime option, medium risk
4. **Issue 3 (P_γ)** — Most impactful on physics, test both forms

## Testing Strategy

After each fix:
1. Compile: `cd build && make -j`
2. Run Tapenade (all fixes touch AD-relevant code)
3. Smoke test: verify no NaN, residuals decrease
4. Compare transition location on NLF0416 validation case

---

## Verified Equations (No Changes Needed)

- Re_θt correlation (Eqs. 8-9) ✓
- F(λ_θ) pressure gradient function (Eqs. 54-57) ✓
- Fonset (Eqs. 46-47) ✓
- Flength (Eqs. 49-50) ✓
- Re_θc (Eq. 51) ✓
- E_γ destruction (Eq. 34) ✓
- SA γ-production coupling (Eq. 41) ✓
- All constants (ca1, ca2, ce1, ce2, c_θt, σ_θt, σ_f) ✓
- Smooth min/max φ_p (Algorithm 1) ✓
