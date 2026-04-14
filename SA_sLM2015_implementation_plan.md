# SA-sLM2015 Transition Model — Detailed Implementation Plan

## Objective

Add the SA–γ–Reθt (SA-sLM2015) transition model to ADflow as a new 3-equation turbulence-transition system. The model extends SA with two transport equations (intermittency γ and transition onset Reθt~).

**First-version goal (decoupled):** γ is NOT multiplied into SA production. All three equations (SA, γ, Reθt~) are solved in a coupled fashion on the turbulence side but the return is identical to pure SA. This lets us verify the solver infrastructure is correct before enabling physics coupling.

---

## Architecture Overview

```
Pure SA:       Q = [ρ, ρu, ρv, ρw, ρE, ν~]            → nw=6, nt2=6, nwt=1
SA-sLM2015:    Q = [ρ, ρu, ρv, ρw, ρE, ν~, γ, Reθt~]  → nw=8, nt2=8, nwt=3
```

The model follows the v2f pattern (3+ equations) and the SST pattern (multi-equation block solver with DADI).

---

## PHASE 0: Preparation & Infrastructure Audit

### Step 0.1 — Define the new turbulence model ID

**File:** `src/modules/constants.F90` (around L129-135)

Add a new model constant:
```fortran
integer(kind=intType), parameter :: spalartAllmarasLM2015 = 11
```

This sits after `v2f = 10`.

**Verification test:** Grep for `spalartAllmarasLM2015` — should appear only in constants.F90.

---

### Step 0.2 — Define model constants module

**New file:** `src/modules/paramTransition.F90`

Define all SA-sLM2015 specific constants and input parameters:
```fortran
module paramTransition
    use constants, only: realType
    implicit none

    ! Intermittency constants
    real(kind=realType) :: ca1 = 2.0
    real(kind=realType) :: ca2 = 0.06
    real(kind=realType) :: ce1 = 1.0
    real(kind=realType) :: ce2 = 50.0
    real(kind=realType) :: sigmaF = 1.0   ! diffusion coefficient for γ

    ! Reθt~ constants
    real(kind=realType) :: sigmaThetat = 2.0
    real(kind=realType) :: cThetat = 0.03
    real(kind=realType) :: cCrossflow = 0.6

    ! Turbulence intensity at freestream
    real(kind=realType) :: tuInfinity = 0.01  ! default 1% Tu
end module paramTransition
```

**Verification test:** Module compiles without errors.

---

### Step 0.3 — Expose new model in Python

**File:** `adflow/pyADflow.py` (around L5669-5676)

Add `"SA-LM2015"` to the list of allowed turbulenceModel values:
```python
"turbulenceModel": [
    str,
    ["SA", "SA-Edwards", "k-omega Wilcox", "k-omega modified",
     "k-tau", "Menter SST", "v2f", "SA-LM2015"],
],
```

Also add a new input option for freestream turbulence intensity (needed for Reθt):
```python
"turbulenceIntensity": [float, 0.01],  # default 1%
```

**File:** `adflow/pyADflow.py` — mapping turbulence model strings to Fortran constants.

Locate the block that maps string → integer and add:
```python
"sa-lm2015": self.adflow.constants.spalartallmaraslm2015,
```

**Verification test:** `solver.setOption('turbulenceModel', 'SA-LM2015')` does not error.

---

### Step 0.4 — Expose model constant in f2py

**File:** `src/f2py/adflow.pyf` (around L10-36)

Add the new constant to the `constants` module section:
```fortran
integer(kind=inttype) parameter :: spalartallmaraslm2015 = 11_inttype
```

**Verification test:** Python can access `solver.adflow.constants.spalartallmaraslm2015`.

---

## PHASE 1: State Vector & Equation Sizing

### Step 1.1 — Set nw, nt2 for the new model

**File:** `src/inputParam/inputParamRoutines.F90` (around L2123-2210)

Add a new case in `setEquationParameters`:
```fortran
case (spalartAllmarasLM2015)
    nw = 8          ! 5 flow + nu~ + gamma + ReThetat~
    nt2 = 8
    eddyModel = .true.
    ! Note: kPresent = .false. (no k equation, no pressure correction)
```

**What this automatically does (no code changes needed):**
- `nwt = nw - nwf = 3` (3 turbulence equations)
- All generic `nt1:nt2` loops in BCs, solver, output now handle 3 variables
- State vector `w(:,:,:,1:8)` allocates correctly
- Scratch arrays resize to accommodate 3×3 blocks
- NK/ANK solver block sizes adjust via `nTurb = nt2 - nt1 + 1 = 3`

**Verification test:** After setting model, verify `nw=8`, `nt2=8`, `nwt=3` from Python side.

---

### Step 1.2 — Verify scratch array sizing

**File:** `src/modules/constants.F90` (around L52-62)

The current scratch array indices:
```
idvt   = 1   (RHS starts here, sized nwt)
ivort  = 3
iprod  = 3
icd    = 4
if1SST = 5
```

For the SA-sLM2015 model (nwt=3), `dvt` will occupy scratch indices 1,2,3.
We need additional scratch slots for auxiliary quantities:
```
idvt      = 1   (indices 1..3: RHS for nu~, gamma, ReThetat~)
iVortLM   = 4   (vorticity magnitude for transition model)
iProdLM   = 5   (production-related scratch)
iFonset   = 6   (Fonset function)
iFlength  = 7   (Flength function)
```

Either reuse existing scratch allocation (if the total is large enough) or increase the scratch array size. Check `turbWork` allocation in the turbulence module.

**File to check:** `src/modules/turbMod.F90` — look at how `scratch` is dimensioned, and if necessary increase `nScratch` for the new model in the allocation routine.

**Verification test:** Allocation succeeds; no out-of-bounds access.

---

## PHASE 2: Boundary Conditions

BCs are handled generically via `nt1:nt2` loops in `turbBCRoutines.F90`. However, we need model-specific values.

### Step 2.1 — Wall boundary conditions

**File:** `src/turbulence/turbBCRoutines.F90` — `bcTurbWall` subroutine

Currently for SA: `w(halo, itu1) = -w(interior, itu1)` (Dirichlet zero)

For SA-sLM2015, add:
```fortran
case (spalartAllmarasLM2015)
    ! nu~ = 0 at wall (same as SA)
    bvt(i,j,itu1) = zero
    bmt(i,j,itu1,itu1) = -one

    ! gamma = 0 at wall 
    bvt(i,j,itu2) = zero
    bmt(i,j,itu2,itu2) = -one

    ! ReThetat~ = zero gradient at wall (Neumann)
    bvt(i,j,itu3) = zero
    bmt(i,j,itu3,itu3) = one   ! halo = +interior → zero gradient
```

**Verification test:** At wall faces, ν~=0, γ=0, ∂Reθt~/∂n=0.

---

### Step 2.2 — Farfield / Inlet boundary conditions

**File:** `src/turbulence/turbBCRoutines.F90` — `bcTurbFarfield` / `bcTurbInflow`

```fortran
case (spalartAllmarasLM2015)
    ! nu~ from eddyVisInfRatio (same as SA)
    bvt(i,j,itu1) = two * wInf(itu1)
    bmt(i,j,itu1,itu1) = -one

    ! gamma = 1 at farfield
    bvt(i,j,itu2) = two * wInf(itu2)   ! wInf(itu2) = 1.0
    bmt(i,j,itu2,itu2) = -one

    ! ReThetat~ = ReThetat(Tu_inf) at farfield 
    bvt(i,j,itu3) = two * wInf(itu3)   ! wInf(itu3) set from correlation
    bmt(i,j,itu3,itu3) = -one
```

**Verification test:** Farfield values are ν~=SA freestream, γ=1, Reθt~=correlation value.

---

### Step 2.3 — Outlet boundary conditions

**File:** `src/turbulence/turbBCRoutines.F90` — `bcTurbOutflow`

All three variables use zero gradient (extrapolation):
```fortran
case (spalartAllmarasLM2015)
    do l = nt1, nt2
        bvt(i,j,l) = zero
        bmt(i,j,l,l) = one    ! halo = interior → zero gradient
    end do
```

The existing generic code may already handle this. Check if the default behavior for unrecognized models is zero-gradient.

**Verification test:** Outlet has extrapolated values.

---

### Step 2.4 — Symmetry boundary conditions

**File:** `src/turbulence/turbBCRoutines.F90` — `bcTurbSymm`

Zero normal gradient for all turbulence variables — should already be generic via `nt1:nt2` loops.

**Verification test:** Symmetry plane has ∂/∂n = 0 for all three variables.

---

## PHASE 3: Initialization

### Step 3.1 — Set freestream values (wInf)

**File:** `src/initFlow/initializeFlow.F90` (around L131-145)

Add a new case:
```fortran
case (spalartAllmarasLM2015)
    ! nu~ from eddy viscosity ratio (same as SA)
    wInf(itu1) = saNuKnownEddyRatio(eddyVisInfRatio, nuInf)

    ! gamma = 1 at freestream (fully turbulent)
    wInf(itu2) = one

    ! ReThetat~ from empirical correlation at freestream Tu
    wInf(itu3) = reThetatCorrelation(tuInfinity)
```

The `reThetatCorrelation` function implements the empirical formula from implement.txt Section 5.6:
- For Tu ≤ 1.3: `Reθt = (1173.51 - 589.428*Tu + 0.2196/Tu²)`
- For Tu > 1.3: `Reθt = 331.50 * (Tu - 0.5658)^(-0.671)`

(No pressure gradient correction at freestream, so F(λθ) = 1.)

**Verification test:** After init, all cells have ν~=SA freestream value, γ=1.0, Reθt~=empirical value.

---

### Step 3.2 — Restart file reading

**File:** `src/initFlow/variableReading.F90` (around L1257-1280)

Add dispatch for the new model:
```fortran
case (spalartAllmarasLM2015)
    call readTurbSALM2015(nTypeMismatch)
```

**New subroutine:** `readTurbSALM2015`
1. Read ν~ (same as SA: look for `cgnsTurbSANu`)
2. Read γ (look for new CGNS name, e.g., `"Intermittency"`)
3. Read Reθt~ (look for new CGNS name, e.g., `"ReThetatTilde"`)
4. If γ or Reθt~ not found in restart: initialize from freestream values
5. If ν~ not found: fall back to eddy viscosity ratio (same as SA)

**Verification test:** Restart from pure SA solution initializes γ=1 and Reθt~ from freestream.

---

### Step 3.3 — Output variable names

**File:** `src/output/outputMod.F90` (around L330-378)

Add CGNS variable names:
```fortran
case (spalartAllmarasLM2015)
    solNames(itu1) = cgnsTurbSaNu          ! Same as SA
    solNames(itu2) = "Intermittency"       ! γ
    solNames(itu3) = "ReThetatTilde"       ! Reθt~
```

**File:** `src/output/outputMod.F90` (around L3215-3250)

Add CGNS model info routine:
```fortran
case (spalartAllmarasLM2015)
    call writeCGNSSaLM2015Info(cgnsInd, base)
```

**Verification test:** CGNS output has 3 turbulence fields with correct names.

---

## PHASE 4: Eddy Viscosity

### Step 4.1 — Add eddy viscosity dispatch

**File:** `src/turbulence/turbUtils.F90` — `computeEddyViscosity` (around L581-655)

Add the new model to the SA dispatch:
```fortran
case (spalartAllmaras, spalartAllmarasEdwards, spalartAllmarasLM2015)
    call saEddyViscosity(iBeg, iEnd, jBeg, jEnd, kBeg, kEnd)
```

The eddy viscosity is IDENTICAL to pure SA: `νt = ν~ * fv1`. The γ and Reθt~ variables do NOT modify eddy viscosity — they only affect SA production (and even that is disabled in v1).

**Verification test:** `rev` array is identical to pure SA for the same ν~ field.

---

## PHASE 5: Turbulence Residual & Source Terms

This is the core of the implementation. We need a new file with the 3-equation block solver.

### Step 5.1 — Create the main solver file

**New file:** `src/turbulence/saLM2015.F90`

This file contains all model-specific routines. Structure mirrors `SST.F90`:

```
saLM2015_block(resOnly)        ! Main orchestrator (like SST_block)
saLM2015Solve(resOnly)         ! 3×3 block solver (like SSTSolve)
saLM2015Source                 ! Source terms for all 3 equations
saLM2015Viscous                ! Diffusion terms for all 3 equations  
saLM2015ResScale               ! Residual scaling
flengthCorrelation(ReThetat)   ! Flength piecewise polynomial
reThetaCCorrelation(ReThetat)  ! Reθc piecewise polynomial
reThetatCorrelation(Tu, lambdaTheta)  ! Reθt empirical correlation
```

---

### Step 5.2 — saLM2015_block (orchestrator)

```fortran
subroutine saLM2015_block(resOnly)
    ! 1. Apply turbulence BCs
    call bcTurbTreatment

    ! 2. Solve the 3-equation system
    call saLM2015Solve(resOnly)

    ! 3. If not residual-only: update eddy viscosity & reapply BCs
    if (.not. resOnly) then
        call saEddyViscosity(2, il, 2, jl, 2, kl)
        call applyAllTurbBCThisBlock(.true.)
    end if
end subroutine
```

This exactly mirrors `SST_block` in structure.

**Verification test:** Subroutine compiles and is callable from `turbAPI.F90`.

---

### Step 5.3 — saLM2015Solve overview

This is the main solver routine, analogous to `SSTSolve`. Its structure:

```
1. Allocate qq(2:il, 2:jl, 2:kl, 3, 3) — 3×3 Jacobian
2. Set dvt => scratch(:,:,:, idvt:idvt+2) — 3-component RHS
3. Compute source terms → fills dvt(:,:,:,1:3) and qq(:,:,:,1:3,1:3)
4. Add advection → turbAdvection(3, 3, 0, qq)
5. Add unsteady terms → unsteadyTurbTerm(3, 3, 0, qq)
6. Add viscous (diffusion) → saLM2015Viscous, updates dvt and qq
7. Scale residuals → dw(:,:,:,itu1:itu3) = -volRef * dvt * iblank
8. If not resOnly:
   a. Apply wall function mods (if applicable)
   b. DADI solve in j, i, k directions (3×3 block tridiagonal)
   c. Update solution: w(:,:,:,itu1:itu3) += factor * dvt
   d. Clip: nu~ ≥ 0, 0 ≤ γ ≤ 1, ReThetat~ ≥ small positive
```

---

### Step 5.4 — Source terms (DECOUPLED VERSION 1)

**CRITICAL:** In the first version, γ does NOT multiply SA production. The SA source term is computed identically to pure SA. The γ and Reθt~ source terms are computed but their feedback into SA is disabled.

```fortran
subroutine saLM2015Source
    ! This routine fills:
    !   dvt(:,:,:,1) = SA source (IDENTICAL to pure SA saSource)
    !   dvt(:,:,:,2) = γ source
    !   dvt(:,:,:,3) = Reθt~ source
    !   qq(:,:,:,1:3,1:3) = Jacobian block

    ! === EQUATION 1: SA (copy from sa.F90 saSource) ===
    ! Compute vorticity/strain magnitude ss
    ! Compute fv1, fv2, ft2, sst, rr, gg, fwSa
    ! 
    ! dvt(i,j,k,1) = (cb1*(1-ft2)*ss + dist2Inv*(...)) * w(i,j,k,itu1)
    !
    ! In DECOUPLED mode: NO γ multiplier here
    ! In COUPLED mode (future): dvt(i,j,k,1) *= gamma(i,j,k)
    !
    ! qq(i,j,k,1,1) = SA Jacobian (same as pure SA)
    ! qq(i,j,k,1,2) = 0   (no dependence on γ in decoupled mode)
    ! qq(i,j,k,1,3) = 0   (no dependence on Reθt~)

    ! === EQUATION 2: γ (intermittency) ===
    ! Compute:
    !   RT = rev(i,j,k) / rlv(i,j,k)   (eddy viscosity ratio)
    !   Ω = vorticity magnitude
    !   Fonset (from ReS, Reθc correlations)
    !   Flength (from Reθt~ correlation)
    !   Fturb = exp(-(RT/4)^4)
    !
    ! Pγ = ca1 * Flength * Fonset * Ω * sqrt(γ) * (1 - ce1*γ)
    ! Eγ = ca2 * Fturb * Ω * γ * (ce2*γ - 1)
    !
    ! dvt(i,j,k,2) = Pγ - Eγ
    !
    ! qq(i,j,k,2,2) = -∂(Pγ-Eγ)/∂γ  (Jacobian diagonal for γ)
    ! qq(i,j,k,2,1) = 0               (no ν~ dependence in diagonal approx)
    ! qq(i,j,k,2,3) = 0               (Reθt~ dependence through correlations, approximate as 0 for DADI)

    ! === EQUATION 3: Reθt~ ===
    ! Compute:
    !   t = 500*μ / (ρ*U²)            (time scale)
    !   Fθt = Fwake * exp(-(d/δ)^4)   (shielding function)
    !   Reθt from empirical correlation
    !
    ! Pθt = (cθt/t) * (Reθt - Reθt~) * (1 - Fθt)
    !
    ! dvt(i,j,k,3) = Pθt
    !
    ! qq(i,j,k,3,3) = (cθt/t) * (1 - Fθt)    (relaxation Jacobian)
    ! qq(i,j,k,3,1) = 0
    ! qq(i,j,k,3,2) = 0
end subroutine
```

**KEY POINT:** The Jacobian `qq` is block-diagonal in the first version (no cross-coupling). This simplifies the DADI solve and makes debugging much easier. The 3×3 system is essentially three independent scalar solves at this stage.

**Verification test:** With γ=1 everywhere and no feedback, SA residual must match pure SA exactly. Print and compare residual norms.

---

### Step 5.5 — Viscous (diffusion) terms

**File:** `src/turbulence/saLM2015.F90`

Three sets of diffusion terms, one per equation. Each in k, j, i direction sweeps (same loop structure as SA/SST).

```fortran
subroutine saLM2015Viscous
    ! === ν~ diffusion (IDENTICAL to saViscous) ===
    ! Uses cb2, cb3 coefficients
    ! Diffusion coefficient: (ν + (1+cb2)*ν~) / σ
    ! Updates dvt(:,:,:,1) and qq(:,:,:,1,1)

    ! === γ diffusion ===
    ! ∇·[(ν + νt/σf) ∇γ]
    ! Diffusion coefficient: (μ + μt/σf) / ρ
    ! Same discretization stencil as SA but different coefficient
    ! Updates dvt(:,:,:,2) and qq(:,:,:,2,2)

    ! === Reθt~ diffusion ===
    ! ∇·[σθt (ν + νt) ∇Reθt~]
    ! Diffusion coefficient: σθt * (μ + μt) / ρ
    ! Same discretization stencil
    ! Updates dvt(:,:,:,3) and qq(:,:,:,3,3)
end subroutine
```

Each direction sweep (k, j, i) follows the EXACT same metric computation as `saViscous` or `SSTSolve` viscous sections:
1. Compute face metrics `ttm`, `ttp`
2. Compute face-averaged diffusion coefficients
3. Apply central-difference stencil: `c_m * φ_{i-1} - c_0 * φ_i + c_p * φ_{i+1}`
4. Accumulate into `dvt` and `qq`
5. Handle boundary faces (k=2, k=kl) using `bmt` matrices

**Verification test:** With uniform γ=1 and uniform Reθt~, diffusion terms for equations 2 and 3 should be zero.

---

### Step 5.6 — Residual scaling

```fortran
subroutine saLM2015ResScale
    ! Scale residuals into dw
    do k = 2, kl
      do j = 2, jl
        do i = 2, il
          dw(i,j,k,itu1) = -dvt(i,j,k,1) * volRef(i,j,k) * real(iblank(i,j,k))
          dw(i,j,k,itu2) = -dvt(i,j,k,2) * volRef(i,j,k) * real(iblank(i,j,k))
          dw(i,j,k,itu3) = -dvt(i,j,k,3) * volRef(i,j,k) * real(iblank(i,j,k))
        end do
      end do
    end do
end subroutine
```

**Verification test:** `dw(itu1)` matches pure SA `dw(itu1)`.

---

### Step 5.7 — DADI block solver (3×3 tridiagonal)

The DADI solver does alternating-direction implicit sweeps. For SA, it's scalar tridiagonal. For SST, it's 2×2 block tridiagonal (calls `tdia3(2, ...)`). For SA-sLM2015, it's 3×3 block tridiagonal.

**However, for the DECOUPLED first version**, since `qq` is block-diagonal (no cross-coupling between equations), the 3×3 block solve degenerates into three independent scalar tridiagonal solves. This is equivalent to:

```fortran
! Option A (simple): Three scalar DADI sweeps
call tdia3(1, jl, bb1, cc1, dd1, ff1)  ! ν~ sweep
call tdia3(1, jl, bb2, cc2, dd2, ff2)  ! γ sweep
call tdia3(1, jl, bb3, cc3, dd3, ff3)  ! Reθt~ sweep
```

OR

```fortran
! Option B (future-proof): One 3×3 block DADI sweep
call tdia3(3, jl, bb, cc, dd, ff)  ! 3×3 block tridiagonal
```

**Recommendation:** Use Option B (`tdia3(3, ...)`) from the start. When the Jacobian is block-diagonal, the 3×3 tridiagonal solver is mathematically equivalent to 3 scalar solves but we get the infrastructure ready for coupling.

Check if `tdia3` supports arbitrary block size or is hardcoded for 2×2. If hardcoded, we need to generalize it or write a `tdia3_3x3` variant.

**File to check:** `src/utils/` or wherever `tdia3` is defined.

**Verification test:** DADI converges; solution update `dvt` goes to zero as residual converges.

---

### Step 5.8 — Solution update and clipping

After the DADI solve:
```fortran
do k = 2, kl
  do j = 2, jl
    do i = 2, il
      ! ν~ update (same as SA)
      w(i,j,k,itu1) = w(i,j,k,itu1) + factor * dvt(i,j,k,1)
      w(i,j,k,itu1) = max(w(i,j,k,itu1), zero)

      ! γ update (clip to [0, 1])
      w(i,j,k,itu2) = w(i,j,k,itu2) + factor * dvt(i,j,k,2)
      w(i,j,k,itu2) = max(zero, min(w(i,j,k,itu2), one))

      ! Reθt~ update (must be positive)
      w(i,j,k,itu3) = w(i,j,k,itu3) + factor * dvt(i,j,k,3)
      w(i,j,k,itu3) = max(w(i,j,k,itu3), 1.0e-10_realType)
    end do
  end do
end do
```

**Verification test:** γ stays in [0,1], Reθt~ stays positive.

---

## PHASE 6: Solver Integration — All Turbulence Solve Paths

ADflow has **five distinct paths** for solving the turbulence equations. The SA-sLM2015 model must work correctly in ALL of them. This section details each path, what happens automatically, and what requires explicit changes.

### 6.0 — Overview: Five Turbulence Solve Paths

```
┌─────────────────────────────────────────────────────────────────────┐
│                     SOLVER PHASE PROGRESSION                        │
│                                                                     │
│  Phase 1: RK/DADI                 Phase 2: ANK                      │
│  ┌────────────────────┐           ┌─────────────────────────────┐   │
│  │ executeMGCycle()   │  switch   │  ANK decoupled              │   │
│  │ └→ turbSolveDDADI  │ ───────>  │  ├→ Flow:  ANK_KSP         │   │
│  │    └→ sa_block()   │  ANK      │  └→ Turb:  ANK_KSPTurb     │   │
│  │       or           │  switch   │     OR     turbSolveDDADI   │   │
│  │    └→ saLM_block() │  tol      │            (ANK_useTurbDADI)│   │
│  └────────────────────┘           │                             │   │
│                                   │  ANK coupled                │   │
│                                   │  └→ Flow+Turb: single KSP   │   │
│                                   └─────────────────────────────┘   │
│                                              │ NK switch tol        │
│                                   ┌──────────▼──────────────────┐   │
│                                   │  Phase 3: NK                │   │
│                                   │  └→ All vars: single KSP    │   │
│                                   └─────────────────────────────┘   │
└─────────────────────────────────────────────────────────────────────┘
```

**Path 1 — RK/DADI:** The DADI smoother calls `turbSolveDDADI` → `saLM2015_block(.false.)`. Turbulence is solved with the internal DADI tridiagonal solver inside `saLM2015Solve`. This is the simplest path.

**Path 2a — ANK Decoupled + Turb KSP:** Flow solved by the main ANK KSP. Turbulence solved by a **separate** `ANK_KSPTurb` (PETSc GMRES) using a matrix-free Jacobian (`dRdwTurb`) and a block-diagonal preconditioner (`dRdwPreTurb`). Turbulence residuals are computed by calling `blocketteRes(useFlowRes=.False., useTurbRes=.True.)` which dispatches to the model-specific source/diffusion routines.

**Path 2b — ANK Decoupled + Turb DADI:** Same as Path 2a for flow, but turbulence uses `turbSolveDDADI` instead of the KSP (activated by `ANK_useTurbDADI=True`). Calls exactly the same code as Path 1 for turbulence.

**Path 3 — ANK Coupled:** All variables (flow + turbulence) solved in a single PETSc KSP. State vector has `nState = nw = 8`. Residuals for ν~, γ, Reθt~ are packed into the same `rVec` as flow residuals, scaled by `turbResScale`.

**Path 4 — NK:** Full Newton-Krylov. Same structure as ANK coupled: single KSP, all variables.

---

### Step 6.1 — turbAPI.F90 dispatch (Paths 1, 2b)

**File:** `src/turbulence/turbAPI.F90`

`turbSolveDDADI` is the entry point for the DADI turbulence solve. It loops over `nSubIterTurb` sub-iterations, calls model-specific blocks for each domain.

**In `turbSolveDDADI` (L31-83), add:**
```fortran
case (spalartAllmarasLM2015)
    call unsteadyTurbSpectral(itu1, itu3)  ! 3 equations: itu1=6, itu2=7, itu3=8
    call saLM2015_block(.false.)
```

**In `turbResidual` (L97-160), add:**
```fortran
case (spalartAllmarasLM2015)
    call saLM2015_block(.true.)   ! resOnly=.true. → compute residual only, no update
```

**What turbSolveDDADI does:**
```fortran
subroutine turbSolveDDADI
    do iter = 1, nSubIterTurb           ! Sub-iterations (default 3)
        do sps = 1, nTimeIntervalsSpectral
            do nn = 1, nDom
                call setPointers(nn, currentLevel, sps)
                select case (turbModel)
                case (spalartAllmarasLM2015)       ! ← NEW
                    call saLM2015_block(.false.)
                end select
            end do
        end do
        call whalo2(groundLevel, nt1, nt2, ...)    ! Exchange halos for all 3 turb vars
    end do
end subroutine
```

**Verification test:** `turbSolveDDADI` calls `saLM2015_block`; turbulence iterations run.

---

### Step 6.2 — turbAdvection call inside saLM2015Solve

**File:** Inside `saLM2015Solve` (new file `src/turbulence/saLM2015.F90`)

The generic `turbAdvection` routine in `turbUtils.F90` handles arbitrary equation counts:
```fortran
call turbAdvection(3, 3, 0, qq)
! mAdv=3 (system size = 3 equations in this turbulence sub-block)
! nAdv=3 (all 3 equations have advection terms)
! offset=0 (no offset; equations start at dvt(:,:,:,1))
```

Similarly for unsteady terms:
```fortran
call unsteadyTurbTerm(3, 3, 0, qq)
```

These are the SAME generic routines SA and SST use with different sizes.

**Verification test:** Advection contribution appears in all 3 residual components.

---

### Step 6.3 — blockette.F90 residual computation (Paths 2a, 3, 4)

**File:** `src/NKSolver/blockette.F90`

`blocketteRes` computes residuals for the ANK/NK solver. It has flags `useFlowRes` and `useTurbRes` to control which parts are computed.

**For turbulence, it calls (around L815-820):**
```fortran
if (equations == RANSEquations .and. turbRes) then
    select case (turbModel)
    case (spalartAllmaras)
        call saSource
        call saAdvection  
        call saViscous
        call saResScale
    case (menterSST)
        call SSTSource
        ...
    end select
end if
```

**Add the new model dispatch:**
```fortran
case (spalartAllmarasLM2015)
    call saLM2015Source       ! 3-equation source terms
    call turbAdvection(3, 3, 0, qq)
    call saLM2015Viscous      ! 3-equation diffusion
    call saLM2015ResScale     ! Scale into dw(itu1:itu3)
```

**IMPORTANT:** The `blocketteRes` path uses the SAME source/diffusion routines as the DADI path. There is ONE implementation of the physics, called from multiple solve paths.

**This produces `dw(i,j,k,itu1)`, `dw(i,j,k,itu2)`, `dw(i,j,k,itu3)`.** These residuals are then used by whatever solver path is active.

**Verification test:** `dw(itu1)` from blockette matches `dw(itu1)` from DADI path for same state.

---

### Step 6.4 — ANK Decoupled Turbulence KSP (Path 2a)

**File:** `src/NKSolver/NKSolvers.F90`

When `ANK_coupled=.false.` AND `ANK_useTurbDADI=.false.` AND `equations==RANSEquations`, the ANK solver creates a separate turbulence KSP system.

#### 6.4.1 — PETSc vector/matrix creation (L1854-1930)

**What happens automatically (generic via nt1:nt2):**
```fortran
nStateTurb = nt2 - nt1 + 1                     ! = 3 for SA-sLM2015
nDimWTurb = nStateTurb * nCellsLocal * nSps     ! Total turb DOF

! PETSc vectors sized automatically:
call VecCreate(comm, wVecTurb, ierr)             ! State: 3 * nCells
call VecSetSizes(wVecTurb, nDimWTurb, ...)
call VecSetBlockSize(wVecTurb, nStateTurb, ...)  ! Block size = 3

call VecDuplicate(wVecTurb, rVecTurb, ...)       ! Residual
call VecDuplicate(wVecTurb, deltaWTurb, ...)     ! Update
call VecDuplicate(wVecTurb, baseResTurb, ...)    ! Line search base
```

**NO CODE CHANGES NEEDED** — all sizing is derived from `nt1`, `nt2`.

#### 6.4.2 — Preconditioner matrix `dRdwPreTurb` (L1881-1907)

```fortran
! Block-diagonal preconditioner: nStateTurb × nStateTurb blocks per cell
call myMatCreate(dRdwPreTurb, nStateTurb, nDimWTurb, nDimWTurb, ...)
! Block size = 3 for SA-sLM2015
```

**NO CODE CHANGES NEEDED** — block size follows `nStateTurb`.

#### 6.4.3 — Matrix-free Jacobian `dRdwTurb` (L1922-1928)

```fortran
call MatCreateMFFD(comm, nDimWTurb, nDimWTurb, ..., dRdwTurb, ...)
call MatMFFDSetFunction(dRdwTurb, FormFunction_mf_Turb, ctx, ...)
```

`FormFunction_mf_Turb` computes the Jacobian-vector product **by finite differences**:
1. Perturbs `w(nt1:nt2)` → calls `blocketteRes(flowRes=.False., turbRes=.True.)`
2. Gets perturbed residual `dw(nt1:nt2)`
3. PETSc computes `(R(w+h*v) - R(w)) / h`

**NO CODE CHANGES NEEDED** — it calls `blocketteRes` which already dispatches to the right model.

#### 6.4.4 — KSP creation (L1931-1946)

```fortran
call KSPCreate(comm, ANK_KSPTurb, ierr)
if (ANK_useMatrixFree) then
    call KSPSetOperators(ANK_KSPTurb, dRdwTurb, dRdwPreTurb, ...)
else
    call KSPSetOperators(ANK_KSPTurb, dRdwPreTurb, dRdwPreTurb, ...)
end if
```

**NO CODE CHANGES NEEDED.**

#### 6.4.5 — Turbulence sub-iteration loop (ANKTurbSolveKSP, ~L3700-3900)

```fortran
do iTurbIter = 1, ANK_nsubIterTurb       ! default: 3 sub-iterations

    ! 1. Compute turbulence residuals
    call blocketteRes(useFlowRes=.False., useTurbRes=.True.)

    ! 2. Pack residuals into PETSc vector
    call setRVecANKTurb(rVecTurb)
    !   → dw(i,j,k,l) * (1/vol) * turbResScale(l-nt1+1) for l=nt1..nt2

    ! 3. KSP solve: dRdw_turb * Δw_turb = -r_turb
    call KSPSolve(ANK_KSPTurb, rVecTurb, deltaWTurb, ierr)

    ! 4. Apply update: w(nt1:nt2) += lambda * deltaW
    call setwVecANK(wVecTurb, nt1, nt2, deltaWTurb, lambda)

    ! 5. Clip to physical bounds
    do l = nt1, nt2
        w(i,j,k,l) = max(w(i,j,k,l), eps * wInf(l))
    end do

end do
```

**Generic via nt1:nt2** — **NO CODE CHANGES NEEDED in the loop structure.**

#### 6.4.6 — Preconditioner fill: `FormJacobianANKTurb` (L2331-2466)

The preconditioner is block-diagonal with `nStateTurb × nStateTurb` blocks:

```fortran
subroutine FormJacobianANKTurb(...)
    ! blk is (nStateTurb, nStateTurb) — i.e. (3,3) for SA-sLM2015
    blk = zero

    do nn = 1, nDom
        do sps = 1, nTimeIntervalsSpectral
            call setPointers(nn, 1, sps)
            do k = 2, kl; do j = 2, jl; do i = 2, il

                dtinv = 1 / (ANK_CFL * dtl(i,j,k) * volRef(i,j,k))

                do l = nt1, nt2
                    l1 = l - nt1 + 1        ! local index: 1, 2, 3
                    blk(l1, l1) = dtinv * turbResScale(l1) / ANK_turbCFLScale
                end do

                irow = globalCell(i,j,k)
                call MatSetValuesBlocked(dRdwPreTurb, 1, irow, 1, irow, blk, ADD_VALUES, ...)

            end do; end do; end do
        end do
    end do
end subroutine
```

**Generic via `do l = nt1, nt2`** — **NO CODE CHANGES NEEDED.** The 3×3 block diagonal is automatically built.

**Verification test:** PETSc turbulence KSP converges in sub-iterations; ν~ residual matches pure SA.

---

### Step 6.5 — ANK Decoupled + Turb DADI (Path 2b)

**File:** `src/NKSolver/NKSolvers.F90` (L4051-4074)

When `ANK_useTurbDADI=True`, the ANK solver bypasses the turbulence KSP and uses the DADI smoother instead:

```fortran
if ((.not. ANK_coupled) .and. equations == RANSEquations .and. lambda > zero) then
    if (ANK_useTurbDADI) then
        call computeUtau            ! Wall friction velocity
        call turbSolveDDADI         ! ← Same as Path 1 DADI
    else
        call ANKTurbSolveKSP        ! ← Path 2a (KSP)
    end if
end if
```

`turbSolveDDADI` calls `saLM2015_block(.false.)` — the exact same code as Path 1.

**ONLY CHANGE NEEDED:** The dispatch in `turbSolveDDADI` (Step 6.1 above).

**Verification test:** With `ANK_useTurbDADI=True`, turbulence converges via DADI within ANK.

---

### Step 6.6 — ANK Coupled (Path 3)

**File:** `src/NKSolver/NKSolvers.F90`

In coupled mode, ALL variables are solved in a single KSP.

#### 6.6.1 — When does coupled mode activate?

Switching logic (L3669-3706):
```fortran
if (totalR <= ANK_coupledSwitchTol * totalR0 .and. equations == RANSEquations) then
    ANK_coupled = .True.
    ! Destroy turbulence KSP objects
    ! Recreate combined KSP with nState = nw (= 8 for SA-sLM2015)
end if
```

**Default `ANK_coupledSwitchTol = 1e-1`** — coupled mode activates when residual drops by 10x.

#### 6.6.2 — Combined state vector sizing

```fortran
if (ANK_coupled) then
    nState = nw        ! = 8 for SA-sLM2015 (5 flow + 3 turb)
else
    nState = nwf       ! = 5 (flow only)
end if
nDimW = nState * nCellsLocal * nSps    ! Total DOF for KSP
```

**Block Jacobian:** `nw × nw` = `8 × 8` blocks per cell.

#### 6.6.3 — Combined residual assembly: setRVec (L1268-1320)

```fortran
! Pack into single rVec:
do nn = 1, nDom; do sps = ...; do k,j,i = ...
    ovv = 1 / volRef(i,j,k)

    ! FLOW residuals (indices 1..5):
    do l = 1, nwf
        rvec_pointer(ii) = dw(i,j,k,l) * ovv
        ii = ii + 1
    end do

    ! TURBULENCE residuals (indices 6..8):
    do l = nt1, nt2
        tmp = dw(i,j,k,l) * ovv * turbResScale(l - nt1 + 1)
        rvec_pointer(ii) = tmp
        ii = ii + 1
    end do
end do; end do; end do
```

**Generic via `nt1:nt2`** — **NO CODE CHANGES NEEDED.** But `turbResScale` must have 3 entries.

#### 6.6.4 — Coupled preconditioner time-stepping (L2155-2210)

```fortran
if (ANK_coupled) then
    ! Only nt1 is set here — this is a POTENTIAL ISSUE for 3 equations!
    stateToCons(nt1, nt1) = turbResScale(1) / ANK_turbCFLScale
end if
```

**⚠ POSSIBLE CHANGE NEEDED:** For multi-equation turbulence, ALL turbulence variables need their time-step diagonal. Check if the characteristic time-stepping section (L2185+) handles this:

```fortran
if (ANK_coupled) then
    ! The characteristic time-stepping currently assumes single turb variable:
    timeStepBlock(6, 6) = one
    streamToCart(6, 6) = one
    symmToCons(nt1, 6) = one
    consToSymm(6, nt1) = one
end if
```

**This is HARDCODED for index 6 (single turb variable).** For 3 turbulence variables, this needs to become:
```fortran
if (ANK_coupled) then
    do l = nt1, nt2
        llocal = l - nt1 + 6    ! Map to column in streamToCart
        timeStepBlock(llocal, llocal) = one
        streamToCart(llocal, llocal) = one
        symmToCons(l, llocal) = one
        consToSymm(llocal, l) = one
    end do
end if
```

**This is a REQUIRED CODE CHANGE** unless `ANK_charTimeStepType = 'None'` (which skips this block entirely).

#### 6.6.5 — Coupled update application

```fortran
call KSPSolve(ANK_KSP, rVec, deltaW, ierr)
call setwVecANK(wVec, 1, nState, deltaW, lambda)

! Clipping:
do l = nt1, nt2
    w(i,j,k,l) = max(w(i,j,k,l), eps * wInf(l))
end do
```

**Generic.** But ensure `wInf(itu2)=1.0` (γ freestream) and `wInf(itu3)=Reθt` are set (Phase 3).

**Verification test:** With coupled ANK, all 8 variables converge together. In decoupled v1, ν~ matches pure SA.

---

### Step 6.7 — NK Solver (Path 4)

**File:** `src/NKSolver/NKSolvers.F90` (L100-240)

The NK solver is structurally identical to ANK coupled: single KSP, all `nw` variables, full Newton step with line search.

**Activation (ANK→NK switch):**
```fortran
if (residual <= NK_switchTol) then
    ! Switch to Newton-Krylov
end if
```

**NO CODE CHANGES NEEDED** — NK uses the same infrastructure as coupled ANK.

**Verification test:** NK converges to machine precision; ν~ solution identical to pure SA.

---

### Step 6.8 — setRVecANKTurb: Turbulence residual packing (L2935-2973)

```fortran
subroutine setRVecANKTurb(rVecTurb)
    call VecGetArrayF90(rVecTurb, rvec_pointer, ierr)
    ii = 0
    do nn = 1, nDom
        do sps = 1, nTimeIntervalsSpectral
            call setPointers(nn, 1, sps)
            do k = 2, kl; do j = 2, jl; do i = 2, il
                ovv = one / volRef(i, j, k)
                do l = nt1, nt2
                    ii = ii + 1
                    rvec_pointer(ii) = dw(i, j, k, l) * ovv * turbResScale(1)  ! ← NOTE: uses turbResScale(1) for ALL
                end do
            end do; end do; end do
        end do
    end do
    call VecRestoreArrayF90(rVecTurb, rvec_pointer, ierr)
end subroutine
```

**⚠ POTENTIAL ISSUE:** At L2962, the code uses `turbResScale(1)` for ALL turbulence equations instead of `turbResScale(l - nt1 + 1)`. For SA (single equation) this is fine. For 3 equations with different scales, this may need to change to:
```fortran
rvec_pointer(ii) = dw(i, j, k, l) * ovv * turbResScale(l - nt1 + 1)
```

**CHECK:** Compare with `setRVec` (the coupled version at L1303) which correctly uses `turbResScale(l - nt1 + 1)` per equation. If the decoupled version uses `(1)` for all, we may need a fix. Alternatively, if all 3 turbResScale entries are the same, it doesn't matter.

**Verification test:** Print packed residual vector entries; verify scaling is correct per equation.

---

### Step 6.9 — turbResScale setup for SA-sLM2015

**File:** `adflow/pyADflow.py` (L6574-6590)

`turbResScale` is stored as a Fortran array of dimension 4 (hardcoded in `inputParam.F90` L293):
```fortran
real(kind=realType), dimension(4) :: turbResScale
```

**This dimension(4) is sufficient** for 3 turbulence equations (uses indices 1,2,3).

**Python default setup** — add SA-sLM2015 defaults in `_updateTurbResScale`:
```python
def _updateTurbResScale(self):
    if self.getOption("turbresscale") is None:
        turbModel = self.getOption("turbulencemodel")
        if turbModel == "SA":
            self.setOption("turbresscale", 10000.0)
        elif turbModel == "Menter SST":
            self.setOption("turbresscale", [1e3, 1e-6])
        elif turbModel == "SA-LM2015":                        # ← NEW
            self.setOption("turbresscale", [10000.0, 1.0, 0.01])
            # ν~ scale: 10000 (same as SA)
            # γ scale: 1.0 (O(1) variable)
            # Reθt~ scale: 0.01 (O(100-1000) variable, scale down)
        else:
            raise Error(...)
```

**Scale rationale:**
- ν~ ≈ 1e-5 → scale 10000 → scaled residual O(0.1)
- γ ≈ 1.0 → scale 1.0 → scaled residual O(1)
- Reθt~ ≈ 100-1000 → scale 0.01 → scaled residual O(1-10)

**These are initial guesses.** Tuning may be needed based on convergence behavior.

**Verification test:** All 3 turbulence residuals are O(1) at the start of ANK.

---

### Step 6.10 — Coupled ANK preconditioner: stateToCons extension

**File:** `src/NKSolver/NKSolvers.F90` (L2170-2210)

For the coupled preconditioner, the `stateToCons` matrix must include ALL turbulence variables. Currently only `nt1` is set:

**Current code (L2174):**
```fortran
stateToCons(nt1, nt1) = turbResScale(1) / ANK_turbCFLScale
```

**Required change for multi-equation turbulence:**
```fortran
if (ANK_coupled) then
    do l = nt1, nt2
        l1 = l - nt1 + 1
        stateToCons(l, l) = turbResScale(l1) / ANK_turbCFLScale
    end do
end if
```

**This fills entries (6,6), (7,7), (8,8)** of the `stateToCons` diagonal with the correct per-equation scaling.

**⚠ ALSO FIX:** Characteristic time-stepping section (L2185-2200). Replace hardcoded index 6 with loop over `nt1:nt2`:
```fortran
if (ANK_coupled) then
    do l = nt1, nt2
        idx = nwf + (l - nt1 + 1)   ! Column index in the extended matrix
        timeStepBlock(idx, idx) = one
        streamToCart(idx, idx) = one
        symmToCons(l, idx) = one
        consToSymm(idx, l) = one
    end do
end if
```

**Verification test:** Preconditioner diagonal has correct entries for all 3 turb variables.

---

### Step 6.11 — FormJacobianANK for coupled mode (L2476-2600)

There is a SEPARATE `FormJacobianANK` for the coupled (flow+turb) preconditioner.

**At L2594:**
```fortran
dtinv * turbResScale(l - nt1 + 1) / ANK_turbCFLScale
```

**This already uses the per-equation index** — **no change needed** if `turbResScale(1:3)` is set correctly.

**Verification test:** Coupled preconditioner converges.

---

### Step 6.12 — FormFunction_mf (coupled) and FormFunction_mf_Turb (decoupled)

Both matrix-free Jacobian functions work by finite-differencing the residual:
1. Perturb state → compute residual → compare with base residual
2. The residual computation goes through `blocketteRes` → model dispatch

**NO CODE CHANGES NEEDED** — as long as `blocketteRes` dispatches to `saLM2015` routines correctly (Step 6.3).

---

### Step 6.13 — Summary: Required Changes in NKSolvers.F90

| Section | Line | Change | Why |
|---------|------|--------|-----|
| `stateToCons` diagonal | L2174 | Loop `nt1:nt2` instead of single `nt1` | Multi-equation turb in coupled mode |
| Characteristic time‐stepping | L2185-2200 | Loop `nt1:nt2` instead of hardcoded `6` | Multi-equation turb in coupled mode |
| `setRVecANKTurb` | L2962 | Use `turbResScale(l-nt1+1)` not `turbResScale(1)` | Per-equation scaling in decoupled mode |
| PETSc setup | L1854-1946 | **No change** | Generic via nStateTurb |
| Preconditioner fill | L2331-2466 | **No change** | Generic via nt1:nt2 loop |
| Sub-iteration loop | L3700-3900 | **No change** | Generic via nt1:nt2 |
| NK solver | L100-240 | **No change** | Uses same coupled infrastructure |

---

### Step 6.14 — Python solver options for SA-sLM2015

All existing ANK/NK Python options apply without modification:

| Option | Default | Effect on SA-sLM2015 |
|--------|---------|---------------------|
| `ANKUseTurbDADI` | `False` | If True, use DADI for turb instead of KSP |
| `ANKTurbCFLScale` | `1.0` | Separate CFL scaling for turb equations |
| `ANKNSubIterTurb` | `3` | Number of turb sub-iterations per ANK step |
| `ANKCoupledSwitchTol` | `1e-1` | Relative residual threshold for coupled mode |
| `ANKPhysicalLSTolTurb` | `0.1` | Physicality tolerance for turb line search |
| `ANKUseMatrixFree` | `True` | Use matrix-free Jacobian |
| `turbResScale` | `None` → auto | **MUST SET** for SA-sLM2015 (3 values) |
| `NKSwitchTol` | `1e-5` | When to switch from ANK to NK |
| `useANKSolver` | `True` | Enable ANK solver |
| `useNKSolver` | `True` | Enable NK solver |

**New option:**
```python
"turbulenceIntensity": [float, 0.01],   # Freestream Tu for transition model
```

**Verification test:** All solver options work without errors; model runs with default settings.

---

### Step 6.15 — Complete solver execution flow for SA-sLM2015

```
User calls: solver(ap)

1. INITIALIZATION
   ├→ setEquationParameters: nw=8, nt2=8, nwt=3
   ├→ initializeFlow: wInf = [..., nuTilde, 1.0, ReThetat]
   ├→ _updateTurbResScale: [10000, 1.0, 0.01]
   └→ allocate w(1:8), dw(1:8) per cell

2. RK/DADI PHASE (initial iterations)
   ├→ executeMGCycle
   │   ├→ Compute flow residual → dw(1:5)
   │   └→ turbSolveDDADI
   │       └→ saLM2015_block(.false.)
   │           ├→ bcTurbTreatment           # Apply BCs to ν~, γ, Reθt~
   │           ├→ saLM2015Source             # Source: SA + γ + Reθt~
   │           ├→ turbAdvection(3,3,0,qq)    # Generic advection
   │           ├→ unsteadyTurbTerm(3,3,0,qq) # Unsteady terms
   │           ├→ saLM2015Viscous            # Diffusion for all 3
   │           ├→ saLM2015ResScale           # → dw(itu1:itu3)
   │           ├→ solveDDADI_3x3            # 3×3 block DADI
   │           ├→ clip values                # ν~≥0, γ∈[0,1], Reθt~>0
   │           ├→ saEddyViscosity            # Update rev (same as SA)
   │           └→ applyAllTurbBCThisBlock    # Reapply BCs
   │
   └→ whalo2(nt1, nt2)                     # Exchange halo for 3 vars

3. ANK DECOUPLED PHASE (residual < ANK_switchTol)
   ├→ Flow: ANK_KSP solves Δw_flow
   └→ Turb: (if ANK_useTurbDADI)
   │    └→ turbSolveDDADI → same as above
   └→ Turb: (if not ANK_useTurbDADI)
        └→ ANKTurbSolveKSP
            ├→ blocketteRes(flow=F, turb=T) → dw(itu1:itu3)
            ├→ setRVecANKTurb: pack dw * (1/vol) * turbResScale
            ├→ KSPSolve(ANK_KSPTurb, rVec, deltaW)
            ├→ setwVecANK: w += λ * Δw
            └→ clip values

4. ANK COUPLED PHASE (residual < ANK_coupledSwitchTol)
   ├→ nState = nw = 8
   ├→ blocketteRes(flow=T, turb=T) → dw(1:8)
   ├→ setRVec: pack all 8 equations into single rVec
   │   └→ turb residuals scaled by turbResScale(1:3)
   ├→ KSPSolve(ANK_KSP, rVec, deltaW)
   ├→ setwVecANK: update all 8 variables
   └→ clip turbulence values

5. NK PHASE (residual < NK_switchTol)
   ├→ Same as coupled ANK but full Newton step
   └→ Converge to machine precision
```

**Verification test for complete flow:**
1. Start with SA-sLM2015 (decoupled physics)
2. RK/DADI phase runs → residuals decrease
3. ANK decoupled phase runs → faster convergence
4. ANK coupled phase runs → quadratic convergence
5. NK phase (if needed) → machine zero
6. Final ν~ field matches pure SA to roundoff

---

## PHASE 7: Correlation Functions

### Step 7.1 — Flength correlation

Piecewise polynomial function of Reθt~. Can be defined in `saLM2015.F90` or a utility module.

```fortran
function flengthCorrelation(ReThetat) result(Flength)
    ! Implement the piecewise polynomial from the reference paper
    ! Must be continuous and smooth at breakpoints
end function
```

### Step 7.2 — Reθc (critical Reynolds number) correlation

```fortran
function reThetaCCorrelation(ReThetat) result(ReTheataC)
    ! Polynomial correlation from reference paper
end function
```

### Step 7.3 — Reθt (empirical transition Reynolds number) correlation

```fortran
function reThetatCorrelation(Tu, lambdaTheta) result(ReThetat)
    ! For Tu <= 1.3:
    !   ReThetat = (1173.51 - 589.428*Tu + 0.2196/Tu²) * F(lambdaTheta)
    ! For Tu > 1.3:
    !   ReThetat = 331.50 * (Tu - 0.5658)^(-0.671) * F(lambdaTheta)
end function
```

### Step 7.4 — F(λθ) pressure gradient correction

```fortran
function pressureGradientCorrection(lambdaTheta) result(Flambda)
    ! Piecewise function depending on sign of lambdaTheta
end function
```

**Verification test:** Unit test each correlation against published tables / reference data. Evaluate at several input points and compare with expected values from the original paper.

---

## PHASE 8: Build System

### Step 8.1 — Add new files to build

**File:** `src/build/fileList`

Add all new source files:
```
TURB_FILES = ... saLM2015.F90
MODULE_FILES = ... paramTransition.F90
```

**File:** `src/build/fort_depend.py`

Ensure dependency scanning picks up the new module.

**File:** `Makefile`

Verify the new files are compiled and linked.

**Verification test:** `make` succeeds without errors.

---

## PHASE 9: Verification & Testing

### Test 9.1 — Compilation test

```bash
make
```
Must compile without errors or warnings.

---

### Test 9.2 — Decoupled equivalence test (CRITICAL)

**Goal:** SA-sLM2015 (decoupled v1) produces EXACTLY the same flow solution as pure SA.

**Procedure:**
1. Run a flat plate case with pure SA model
2. Run the SAME case with SA-sLM2015 (decoupled: γ not multiplied into SA production)
3. Compare:
   - ν~ field: must be identical (to machine precision)
   - Eddy viscosity field: must be identical
   - Flow solution (ρ, u, v, w, p): must be identical
   - Cf (skin friction): must be identical
   - CD, CL: must be identical
4. Additionally verify:
   - γ converges to its steady-state (should be ~1 everywhere since fully turbulent)
   - Reθt~ converges to its steady-state

**Python test script:**
```python
# Test: SA-sLM2015 decoupled matches pure SA
solver_sa = ADFLOW(options={..., 'turbulenceModel': 'SA'})
solver_sa(ap)
cl_sa = funcs['cl']
cd_sa = funcs['cd']

solver_lm = ADFLOW(options={..., 'turbulenceModel': 'SA-LM2015'})
solver_lm(ap)
cl_lm = funcs['cl']
cd_lm = funcs['cd']

assert abs(cl_sa - cl_lm) < 1e-10
assert abs(cd_sa - cd_lm) < 1e-10
```

---

### Test 9.3 — Convergence test

Run the SA-sLM2015 model and verify:
1. All 3 turbulence residuals decrease monotonically (or at least converge)
2. The ANK solver has reasonable CFL ramp
3. No NaN or Inf in any variable
4. γ stays in [0,1]
5. Reθt~ stays positive and reasonable

---

### Test 9.4 — Boundary condition test

Extract values at boundaries and verify:
- Wall: ν~=0, γ=0, ∂Reθt~/∂n=0
- Farfield: ν~=freestream, γ=1, Reθt~=correlation value
- Outlet: zero gradient
- Symmetry: zero normal gradient

---

### Test 9.5 — Restart test

1. Run SA-sLM2015 for 500 iterations, write restart
2. Read restart file, run 500 more iterations
3. Compare with 1000-iteration run: must match to machine precision


---

## PHASE 7: Adjoint, Tapenade & PC Jacobian Integration — VERIFIED ✅

This phase documents the comprehensive audit of how the SA-GR model connects to ADflow's adjoint infrastructure, Tapenade-generated AD code, and preconditioner Jacobian assembly. All solver and adjoint paths have been verified for correct 3-equation handling.

### 7.0 — Architecture Overview

ADflow has two paths for assembling the preconditioner Jacobian matrix:

```
┌──────────────────────────────────────────────────────────────┐
│  setupStateResidualMatrix (adjointUtils.F90)                 │
│  ├── useAD=True (ANKADPC=True)                              │
│  │   └── block_res_state_d (masterRoutines.F90 L1316)       │
│  │       └── saGRSource_d, turbAdvection_d(3,3,itu1-1,qqGR) │
│  │           saGRViscous_d, saGRResScale_d                   │
│  └── useAD=False (ANKADPC=False, default)                    │
│      └── block_res_state (masterRoutines.F90 L1245)          │
│          └── useBlockettes=False → blockResCore              │
│              └── saGammaRetheta_block(.true.)                │
└──────────────────────────────────────────────────────────────┘
```

**Key design decision**: `useBlockettes` is forced to `False` for SA-GR because `blocketteResCore` (blockette.F90 L299-753) only has SA-specific reimplementations of saSource/saViscous/saAdvection/saResScale. There are NO blockette SA-GR equivalents. When `useBlockettes=True` and SA-GR is active, the select case falls through with no match, producing zero turbulence residuals for γ and Reθt — a silent, catastrophic bug.

### 7.1 — masterRoutines.F90 — All 5 Dispatch Points ✅

All five adjoint/residual paths properly dispatch SA-GR:

| Path | Lines | SA-GR Routines Called | Status |
|------|-------|----------------------|--------|
| Primal (`master`) | L195-213 | allocate qqGR, saGRSource, turbAdvection(3,3,itu1-1,qqGR), saGRViscous, saGRResScale, deallocate qqGR | ✅ |
| Forward AD (`master_d`) | L543-556 | saGRSource_d, turbAdvection_d(3), saGRViscous_d, saGRResScale_d | ✅ |
| Reverse AD (`master_b`) | L797-814 | saGRResScale_b, saGRViscous_b, turbAdvection_b(3), saGRSource_b (reverse order) | ✅ |
| Fast reverse (`master_state_b`) | L1159-1174 | saGRResScale_fast_b, saGRViscous_fast_b, turbAdvection_fast_b(3), saGRSource_fast_b | ✅ |
| PC Jacobian AD (`block_res_state_d`) | L1378-1393 | saGRSource_d, turbAdvection_d(3), saGRViscous_d, saGRResScale_d | ✅ |

**Import statements** are present for all paths:
- Primal: L23 — `use saGammaRetheta, only: saGRSource => Source, ...`
- Forward: L288 — `use saGammaRetheta_d, only: saGRSource_d => Source_d, ...`
- Reverse: L671 — `use saGammaRetheta_b, only: saGRSource_b => Source_b, ...`
- Fast reverse: L1081 — `use saGammaRetheta_fast_b, only: saGRSource_fast_b => Source_fast_b, ...`

### 7.2 — Tapenade AD Files — All 3 Modes ✅

All three AD-generated files exist and contain the correct smooth formulations (tanh Fonset, modified Fturb, vorticity limiting via smoothMinMax):

| Mode | File | Subroutines | Status |
|------|------|-------------|--------|
| Forward (_d) | `src/adjoint/outputForward/saGammaRetheta_d.f90` | source_d (L68), viscous_d (L1218), resscale_d (L2276) + primal copies | ✅ |
| Reverse (_b) | `src/adjoint/outputReverse/saGammaRetheta_b.f90` | source_b (L70), viscous_b (L1331), resscale_b (L2574) + primal copies | ✅ |
| Fast reverse (_fast_b) | `src/adjoint/outputReverseFast/saGammaRetheta_fast_b.f90` | source_fast_b (L65), viscous_fast_b (L1276), resscale_fast_b (L2433) + primal copies | ✅ |

Supporting AD files also handle SA-GR:
- `turbUtils_{d,b,fast_b}.f90` — `computeEddyViscosity` dispatches SA-GR to `saEddyViscosity` ✅
- `turbBCRoutines_{d,b,fast_b}.f90` — Wall/farfield BCs handle SA-GR (4 dispatch points in _d/_b, 2 in _fast_b) ✅
- `initializeFlow_{d,b,fast_b}.f90` — Reference state initialization includes SA-GR ✅

**Note on `qq`/`qqGR` in AD paths**: The DADI diagonal array (`qq`) is declared as `intent(inout)` in `turbAdvection_d/b/fast_b` but is a **dead variable** — Tapenade stripped all actual reads/writes to it. Passing an unallocated `qqGR` module alias is safe because it's never dereferenced. This matches the existing SA pattern.

### 7.3 — AD Preconditioner Path (ANKADPC=True) ✅

When `ANKADPC=True`, `setupStateResidualMatrix` calls `block_res_state_d(nn, sps)` for each domain/spectral instance per coloring step.

**Full call chain:**
```
FormJacobianANKTurb (NKSolvers.F90 L2336)
  └─ setupStateResidualMatrix(useAD=.True., useTurbOnly=.True.)  (adjointUtils.F90)
       └─ block_res_state_d(nn, sps)  (masterRoutines.F90 L1316)
            ├─ computeEddyViscosity_d(.True.)      ← SA-GR uses saEddyViscosity_d
            ├─ BCTurbTreatment_d / applyAllTurbBCthisblock_d
            ├─ applyAllBC_block_d
            ├─ timeStep_block_d
            ├─ case (spalartallmarasnoft2gammaretheta):
            │    call saGRSource_d                  ← Tapenade forward mode
            │    call turbAdvection_d(3,3,itu1-1,qqGR)
            │    call saGRViscous_d                 ← Tapenade forward mode
            │    call saGRResScale_d                ← Tapenade forward mode
            ├─ [mean flow residuals: inviscidCentralFlux_d, viscousFlux_d, ...]
            ├─ sumDwAndFw_d
            └─ resscale_d    ← adjointExtra_d; scales dwd by turbResScale(1:3)
```

**Key properties of this path:**
- `useBlockettes` is **completely irrelevant** — `block_res_state_d` does not check it at all. Even if `useBlockettes=True`, the ANKADPC=True path works perfectly for SA-GR.
- The AD seed `flowDomsd(nn,1,sps)%w(i,j,k,l) = one` is set for one (i,j,k,l) cell at a time per coloring iteration. The resulting `dwd` gives the exact Jacobian column.
- `qqGR` (alias of the module-level `qq` in `saGammaRetheta_d`) is passed to `turbAdvection_d` but is a **dead argument** — Tapenade removed all actual `qq(...)` references from the differentiated `turbAdvection_d`. The unallocated array is never dereferenced.
- `saGRSource_d`, `saGRViscous_d`, `saGRResScale_d` do not use `qq` internally. They operate only on `wd`, `dwd`, and flow-field arrays.
- The 3×3 Jacobian blocks for SA-GR are assembled cell-by-cell through the coloring loop — the same mechanism as SA (1×1) and SST (2×2).

**When to use ANKADPC=True:**
- Provides **exact** first-order PC Jacobian (up to the coloring approximation)
- Costs more compute time than FD (ANKADPC=False) due to larger Tapenade code
- Recommended when FD-based preconditioner convergence is poor (stiff residuals)
- Not needed for convergence-only runs; mainly beneficial for adjoint quality

### 7.4 — FD Preconditioner Path (ANKADPC=False, default) ✅

When `ANKADPC=False` (default), `setupStateResidualMatrix` perturbs each state variable and calls `block_res_state` for finite-difference column assembly:

```
FormJacobianANKTurb (NKSolvers.F90 L2336)
  └─ setupStateResidualMatrix(useAD=.False., useTurbOnly=.True.)  (adjointUtils.F90)
       └─ [perturb w(i,j,k,l) += delta_x, then:]
            block_res_state(nn, sps)  (masterRoutines.F90 L1245)
              ├─ computeEddyViscosity / BCTurbTreatment / applyAllBC
              ├─ useBlockettes?
              │    True  → blocketteResCore  ← ❌ NO SA-GR CASE → zero turb residual
              │    False → blockResCore      ← ✅ dispatches saGammaRetheta_block(.true.)
              ├─ sourceTerms_block
              └─ resScale  ← adjointExtra.F90; scales dw by turbResScale(1:3)
```

**Critical fix applied**: `useBlockettes` is now forced to `False` for SA-GR in `_updateTurbResScale()` (pyADflow.py), **outside** the `turbresscale is None` block, so it applies even when a user explicitly sets turbresscale.

### 7.4 — NKSolvers.F90 — All Generic Paths ✅

The NK/ANK solver infrastructure is **fully parameterized** using `nt1:nt2` loops and `nStateTurb = nt2 - nt1 + 1`. No model-specific dispatches needed:

| Component | Lines | Handles 3 eq? | Mechanism |
|-----------|-------|---------------|-----------|
| `nStateTurb` computation | L1855 | ✅ | `nt2 - nt1 + 1 = 3` |
| `FormJacobianANKTurb` | L2336-2440 | ✅ | `blk(nStateTurb, nStateTurb)`, loops `nt1:nt2` |
| `FormFunction_mf_turb` | L2545-2617 | ✅ | calls `blocketteRes(useFlowRes=.False.)`, loops `nt1:nt2` |
| `setWANK` / `setW` | L2980-3016 / L1340-1376 | ✅ | Generic `lStart:lEnd` / `nt1:nt2` |
| `physicalityCheckANK` | L3018-3223 | ✅ | Has gamma/ReTheta clipping for `nt1+1:nt2` |
| `physicalityCheckANKTurb` | L3225-3340 | ✅ | Same gamma/ReTheta clipping |
| `setRVecANKTurb` | L2962 | ✅ | `turbResScale(l - nt1 + 1)` per equation |
| `stateToCons` diagonal (coupled) | L2175 | ✅ | Loops `nt1:nt2` with per-equation `turbResScale` |
| ANK turb backtracking | L3519-3565 | ✅ | Generic line search |

### 7.5 — Remaining Infrastructure ✅

| Component | Status | Details |
|-----------|--------|---------|
| `turbResScale` array | ✅ | `dimension(4)` in inputParam.F90 — sufficient for 3 eq |
| `resScale` (adjointExtra.F90) | ✅ | Generic: `nTurb = nt2 - nt1 + 1`, loops all turb vars |
| `frozenTurbulence` | ✅ | Generic flag, no model-specific behavior |
| `adjointAPI.F90` | ✅ | No turbulence model dispatches; all generic |
| Eddy viscosity | ✅ | SA-GR uses `saEddyViscosity` (only reads ν~/itu1) |
| `blockResCore` import | ✅ | `use saGammaRetheta, only: saGammaRetheta_block` at blockette.F90 L771 |

### 7.6 — Known Limitation: blocketteResCore

`blocketteResCore` (blockette.F90 L299-753) contains SA-specific reimplementations optimized for cache-blocked ("blockette") memory layout. These are NOT the same as the module-level saSource/saViscous routines — they're entirely rewritten for the blockette local memory scheme. Implementing SA-GR blockette routines is deferred to a future optimization phase. For now, `useBlockettes=False` bypasses the gap with no functional impact.

---

## PHASE 10: Enable Coupling (Future v2)

Once the decoupled version is verified, enable physics coupling:

### Step 10.1 — Multiply γ into SA production

In `saLM2015Source`, equation 1:
```fortran
! BEFORE (decoupled):
! dvt(i,j,k,1) = cb1*(1-ft2)*ss*nutilde + ...

! AFTER (coupled):
! dvt(i,j,k,1) = gamma(i,j,k) * cb1*(1-ft2)*ss*nutilde + ...
```

This is a single-line change.

### Step 10.2 — Add cross-Jacobian terms

The Jacobian `qq` gets off-diagonal entries:
```fortran
qq(i,j,k,1,2) = cb1*(1-ft2)*ss*nutilde  ! ∂R_SA/∂γ = production term
```

### Step 10.3 — Verification

Run a transitional flat plate case and verify:
- Transition location matches published data
- Upstream of transition: γ ≈ 0, flow is laminar
- Downstream of transition: γ → 1, flow is turbulent
- Skin friction coefficient shows laminar-to-turbulent transition

---

## Summary of ALL Files to Modify

| File | Change | Phase |
|------|--------|-------|
| `src/modules/constants.F90` | Add `spalartAllmarasLM2015 = 11` | 0 |
| `src/modules/paramTransition.F90` | **NEW** — model constants | 0 |
| `src/f2py/adflow.pyf` | Add new constant | 0 |
| `adflow/pyADflow.py` | Add model string, turbulenceIntensity option, turbResScale defaults | 0, 6 |
| `src/inputParam/inputParamRoutines.F90` | Add case for nw=8, nt2=8 | 1 |
| `src/turbulence/turbBCRoutines.F90` | Add wall/farfield/inlet/outlet BCs | 2 |
| `src/initFlow/initializeFlow.F90` | Add freestream init for γ, Reθt~ | 3 |
| `src/initFlow/variableReading.F90` | Add restart reading | 3 |
| `src/output/outputMod.F90` | Add CGNS names and model info | 3 |
| `src/turbulence/turbUtils.F90` | Add to eddy viscosity dispatch | 4 |
| `src/turbulence/saLM2015.F90` | **NEW** — main solver file (source, diffusion, DADI) | 5 |
| `src/turbulence/turbAPI.F90` | Add dispatch case in turbSolveDDADI and turbResidual | 6 |
| `src/NKSolver/blockette.F90` | Add model dispatch in turbulence residual section | 6 |
| `src/NKSolver/NKSolvers.F90` | Fix `stateToCons` loop for multi-eq turb (L2174), fix char time-step (L2185-2200), fix `setRVecANKTurb` scaling (L2962) | 6 |
| `src/build/fileList` | Add new source files | 8 |

**Total new files:** 2 (saLM2015.F90, paramTransition.F90)
**Total modified files:** 13

---

## Complete SA Touchpoint Checklist

Every location where pure SA interacts with the code, and what must be done for SA-sLM2015:

| # | SA Touchpoint | File | SA-sLM2015 Action |
|---|---------------|------|-------------------|
| 1 | Model ID constant | constants.F90 | Add `spalartAllmarasLM2015 = 11` |
| 2 | Model string → ID mapping | pyADflow.py | Add `"SA-LM2015"` |
| 3 | f2py constant export | adflow.pyf | Add constant |
| 4 | Equation count (nw, nt2) | inputParamRoutines.F90 | nw=8, nt2=8 |
| 5 | Freestream init (wInf) | initializeFlow.F90 | Add γ=1, Reθt~=corr |
| 6 | Restart reading | variableReading.F90 | Add 3-variable read |
| 7 | Wall BC | turbBCRoutines.F90 | ν~=0, γ=0, ∂Reθt~/∂n=0 |
| 8 | Farfield BC | turbBCRoutines.F90 | ν~=SA, γ=1, Reθt~=corr |
| 9 | Inlet BC | turbBCRoutines.F90 | Same as farfield |
| 10 | Outlet BC | turbBCRoutines.F90 | Zero gradient all |
| 11 | Symmetry BC | turbBCRoutines.F90 | Zero gradient all (likely generic) |
| 12 | Eddy viscosity | turbUtils.F90 | Reuse saEddyViscosity |
| 13 | Source term computation | sa.F90 → saLM2015.F90 | Copy+extend SA source |
| 14 | Diffusion terms | sa.F90 → saLM2015.F90 | Copy SA diffusion + add γ, Reθt~ diffusion |
| 15 | Residual scaling | sa.F90 → saLM2015.F90 | Extend to 3 equations |
| 16 | DADI solver | sa.F90 → saLM2015.F90 | 3×3 block version |
| 17 | Solver dispatch (DDADI) | turbAPI.F90 | Add case |
| 18 | Solver dispatch (residual) | turbAPI.F90 | Add case |
| 19 | Advection | turbUtils.F90 (turbAdvection) | Generic; call with mAdv=3 |
| 20 | Unsteady terms | turbUtils.F90 | Generic; call with itu1:itu3 |
| 21 | Solution clipping | saLM2015.F90 | ν~≥0, γ∈[0,1], Reθt~>0 |
| 22 | CGNS output names | outputMod.F90 | Add 3 variable names |
| 23 | CGNS model info | outputMod.F90 | Add info writer |
| 24 | NK/ANK residual loop | NKSolvers.F90 | Generic (nt1:nt2); verify scaling |
| 25 | NK/ANK CFL scaling | NKSolvers.F90 | Fix stateToCons loop for multi-eq turb |
| 26 | Preconditioner block size | NKSolvers.F90 | Auto via nTurb; verify |
| 27 | ANK decoupled turb residual packing | NKSolvers.F90 (setRVecANKTurb) | Fix per-equation turbResScale |
| 28 | ANK coupled char time-stepping | NKSolvers.F90 (L2185-2200) | Loop nt1:nt2 instead of hardcoded 6 |
| 29 | ANK decoupled turb KSP setup | NKSolvers.F90 (L1854-1946) | Generic; no change |
| 30 | ANK decoupled turb preconditioner fill | NKSolvers.F90 (FormJacobianANKTurb) | Generic; no change |
| 31 | ANK coupled preconditioner fill | NKSolvers.F90 (FormJacobianANK) | Generic; no change |
| 32 | ANK→coupled switch | NKSolvers.F90 (L3669-3706) | Generic; no change |
| 33 | ANK_useTurbDADI path | NKSolvers.F90 (L4051-4074) | Generic; uses turbSolveDDADI |
| 34 | Matrix-free turb Jacobian | NKSolvers.F90 (FormFunction_mf_Turb) | Generic; uses blocketteRes |
| 35 | blockette turb dispatch | blockette.F90 (L815-820) | Add saLM2015 source/diffusion/scale |
| 36 | SA constants usage | paramTurb.F90 | Reuse (no change) |
| 37 | Transition constants | paramTransition.F90 | **NEW** module |
| 38 | Wall distance | wallDistance.F90 | Generic; no change |
| 39 | Build system | fileList, Makefile | Add new files |
| 40 | Python option | pyADflow.py | Add model + Tu + turbResScale |
| 41 | turbResScale defaults | pyADflow.py (_updateTurbResScale) | Add [10000, 1.0, 0.01] |
| 42 | Adjoint (AD) | adjoint/ | Phase 10+ (skip for now) |

---

## Recommended Implementation Order

```
Phase 0: Constants, Python, f2py                    [~1 day]
Phase 1: State vector sizing                        [~0.5 day]
Phase 2: Boundary conditions                        [~1 day]
Phase 3: Initialization + output                    [~1 day]
Phase 4: Eddy viscosity dispatch                    [~0.5 day]
    → CHECKPOINT: Code compiles, model selectable, initializes correctly

Phase 5: Source + diffusion + solver (saLM2015.F90) [~3-5 days]
    → This is the bulk of the work
    → Start with SA source copied verbatim
    → Add γ source (stub: returns 0)
    → Add Reθt~ source (stub: returns 0)
    → Add diffusion for all 3
    → Add DADI solver
    → CHECKPOINT: Model runs, produces pure SA results

Phase 6: Solver integration                         [~1 day]
    → turbAPI dispatch
    → NK/ANK verification
    → CHECKPOINT: Full solver loop works

Phase 7: Correlation functions                      [~1-2 days]
    → Implement Flength, Reθc, Reθt, F(λθ)
    → Unit test against reference values
    → Wire into γ and Reθt~ source terms
    → CHECKPOINT: γ and Reθt~ have nontrivial residuals

Phase 8: Build system                               [~0.5 day]
Phase 9: Testing & verification                     [~2-3 days]
    → Decoupled equivalence test
    → Convergence test
    → BC test
    → Restart test

Phase 10: Enable coupling (future)                  [~1 day]
    → Single-line change to multiply γ into SA production
    → Add cross-Jacobian terms
    → Validate against transition experiments
```

---

## Key Risks & Mitigations

1. **`tdia3` block size limitation:** If hardcoded for ≤2, need to generalize. Mitigation: check `tdia3` source immediately; if limited, write `tdia3_block(n, ...)` variant.

2. **Scratch array overflow:** With 3 turbulence equations, `dvt` uses 3 slots instead of 1-2. Mitigation: verify scratch allocation before writing solver code.

3. **Residual scaling mismatch:** γ~O(1) vs Reθt~~O(1000) vs ν~~O(1e-5) can cause NK/ANK convergence issues. Mitigation: set proper `turbResScale` values; test ANK convergence early.

4. **Correlation function discontinuities:** Piecewise polynomials for Flength and Reθc can have jumps at breakpoints. Mitigation: verify C0/C1 continuity; add smoothing if needed.

5. **Adjoint support:** Not addressed in this plan (Phase 10+). The AD tool (Tapenade) will need to differentiate the new code. Mitigation: write clean Fortran without goto, pointer aliasing, or other AD-unfriendly constructs.

---

## Notes for Implementer

- **Reuse code aggressively.** The SA source, diffusion, and solver are well-tested. Copy them as the starting point for ν~ in saLM2015.F90. Do NOT re-derive — copy and extend.
- **The γ and Reθt~ equations have the same STRUCTURE as SA** (convection + diffusion + source). Only the coefficients and source terms differ. Use the same discretization.
- **Start simple.** Get the infrastructure right with stub source terms (return 0 for γ and Reθt~). Then add physics one equation at a time.
- **Test decoupled first.** The first version must reproduce pure SA exactly. Any discrepancy means a bug in the infrastructure.
- **Wall distance** (`d2wall`) is already available in all turbulence routines via `blockPointers`. No special handling needed.
- **Vorticity and strain** are already computed in the SA source routine. Reuse for γ and Reθt~ source terms — don't recompute.
- **The ANK/NK solver is largely generic.** Most code uses `nt1:nt2` loops that automatically handle 3 equations. The main exceptions are the `stateToCons` diagonal (hardcoded for index 6), the characteristic time-stepping section, and the `setRVecANKTurb` scaling — these need explicit fixes.
- **There are 5 turbulence solve paths** (RK/DADI, ANK decoupled KSP, ANK decoupled DADI, ANK coupled, NK). The physics routines (source, diffusion) are called from ALL paths via the same dispatch. Write the physics ONCE, ensure the dispatch works, and all 5 paths will work.
- **The `turbResScale` array is dimension(4)** in Fortran (hardcoded in `inputParam.F90`). This is sufficient for 3 equations. Values are set from Python via `_updateTurbResScale` and must be tuned for convergence.
- **Clipping matters.** γ must stay in [0,1], Reθt~ must stay positive. The ANK update applies clipping after every step (`w = max(w, eps*wInf)`). Since `wInf(itu2) = 1.0` (γ) and `wInf(itu3) = Reθt(Tu∞)`, the clipping naturally prevents pathological values.


## Solver specifics

Additional stabilization measures (as implemented in the paper)

The paper introduces several solver-level modifications to improve robustness beyond source-term timestep restriction.

A) Equation and variable scaling

The transition variables (γ, Re_theta_t, ν~) have very different magnitudes from the flow variables, which leads to poor conditioning of the linear system.

They solve a scaled linear system:

S_a S_r ( I/Δt + A^(n) ) S_c S_c^{-1} ΔQ^(n) = − S_a S_r R^(n)

Row scaling matrix (per node i):

S_r,i = diag( J_i^(2/3), ..., J_i^(2/3), ν~_max^{-1} J_i^{-1/3}, γ_max^{-1} J_i^{-1/3}, Re_theta_t_max^{-1} J_i^{-1/3} )

Column scaling matrix:

S_c,i = diag( 1, ..., 1, ν~_max, γ_max, Re_theta_t_max )

Typical values:
ν~_max = 10^3
γ_max = 10
Re_theta_t_max = 10^4

This improves conditioning and Krylov convergence.

B) Solution update damping (instead of clipping)

To prevent unstable updates in γ and Re_theta_t, they apply iterative damping of the Newton update.

For each node:

Q^(n+1) = Q^(n) + ΔQ^(n)

If γ is outside bounds:
γ > 2 or γ < 1e−10

then repeatedly apply:

Δγ → (θ)^m Δγ, with θ = 0.99

until bounds are satisfied.

Similarly for Re_theta_t:
Re_theta_t ≥ 20 is enforced via the same damping mechanism.

Important point:

No hard clipping
Hard limits caused solver stalling in their tests

C) Backtracking line search

A line-search/backtracking algorithm is used to stabilize Newton updates.

Purpose:

avoid divergence when Newton step is too aggressive
ensure residual decreases

This is required because transition source terms can cause large nonlinear jumps.

D) Source-term timestep switching

Source-term timestep restriction:

λ_source * Δt ≤ 0.9

is applied only when needed:

Active during approximate-Newton phase
Deactivated after stable convergence (no backtracking for several iterations)
Reactivated if:
backtracking occurs
residual increases

This avoids overly small timesteps in late convergence.

E) Smoothing of non-differentiable functions

All min/max and conditional functions are replaced by smooth approximations:

General form:

φ_p(g1, g2) = g1 + (1/p) log(1 + exp(p (g2 − g1)))

with proximity switch:

|g1 − g2| > −log(|p| * p_switch)/|p| → use exact min/max
otherwise → use smooth function

Typical:
p = ±300
p_switch = 1e−15

Purpose:

remove Jacobian discontinuities
enable Newton convergence
