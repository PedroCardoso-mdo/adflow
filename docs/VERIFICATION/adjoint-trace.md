# SA → SA-γ-R̃e_θt Adjoint Plumbing Trace

*Historical strata: header banners are dated snapshots (2026-05-20 base +
later updates); §§1-12 are reference. Line numbers into `src/` have
drifted — locate by symbol name.*

**Repository**: `/home/mdo/MDOLab_3_v2/adflow_sa_gamma_rethetha`  
**Date**: 2026-05-20  
**Purpose**: Research-only paired inventory of adjoint/AD touchpoints

> **Regeneration status (2026-07-07, audit):** Tapenade WAS rerun after the
> `transitionRefLength` (272ff47c) and farfield-BC (701668a6) primal changes —
> the regenerated `saGammaRetheta_{d,b,fast_b}.f90` and
> `turbBCRoutines_{d,b,fast_b}.f90` sat uncommitted in the working tree at
> the time *(2026-08-12: since committed — Tapenade has been regenerated
> several times since; latest 2026-08-06, commit 2c4ce2c1)* and
> verifiably contain both changes (vortLim `refLenTrans` branch; generic
> farfield ghost form, old `2*wInf` special case removed). Build compiles and
> imports with them.
>
> **Tapenade regeneration DONE (2026-07-07, `./AD_I.sh`):** the `uInf, muInf`
> active declaration in `saGammaRetheta%Source` (fullRoutines head) is now
> baked into the generated files — `vortlimd` carries the `uinfd`/`muinfd`
> contributions (the old `vortlimd = 0.0_8` is gone). The vorticity limiter
> (P&Z Eqs. 52-53) Mach/Re partials are therefore live. State partials
> (dR/dw) were never affected, so `fast_b`/stateOnlyRoutines are unchanged by
> design. Additionally, `turbAdvection` (turbUtils.F90) now forces first-order
> transition advection for SA-GR when `transitionFirstOrderUpwind`, and this
> guard is present in all three regenerated `turbutils` AD files — primal and
> adjoint advection discretizations now agree for any `orderTurb` setting.
> The redundant `orderTurb` swap in `saGammaReTheta_block` was removed
> (2026-07-08): it covered only the decoupled path and mutated a global,
> invisible to the AD sweeps; the in-`turbAdvection` guard is the single
> point of truth. Full finding/fix record: `docs/audits/06_adjoint_wiring.md`.
>
> **Why uInf/muInf are in the residual at all (and the design choice):** the
> cap `vortLim = uInf·√(uInf/(muInf·l))/20` is ADflow's p-ρ spelling of the
> paper's `M·√(M·Re)/20` (Eqs. 52-53) — the model itself defines the cap from
> freestream conditions (laminar-BL wall-vorticity scale `(U∞/l)·√Re`; paper
> §IV calls it "a physical scaling independent of the nondimensionalization").
> Derivation in `nondimensionalization.md` §5 (D1 exception). Since the paper
> also states the limiter is purely numerical (non-predictive), a legitimate
> alternative is to **freeze the cap** ("frozen limiter" practice): revert the
> `Makefile_tapenade` line and accept a small `dR/dMach`-type error confined to
> cells where the cap engages. dR/dw is identical under both options. Current
> decision: differentiate as written (option 1), so partials verification can
> compare against the exact discrete residual.

---

## Tapenade regeneration 2026-08-04 — AD debt paid + a reusable gotcha

Triggered by making `epsAcoustic`/`epsShear` runtime options (see
`architecture.md`, "Matrix-dissipation eigenvalue limiters"): the generated
`fluxes_*.f90` had those values baked in as local `parameter`s and had to be
regenerated. Regenerate with `./AD_I.sh` from the repo root (runs
`ad_forward`, `ad_reverse`, `ad_reverse_fast`, then `make`; it does **not**
`pip install` — do that separately or site-packages stays stale).

**Two things this regeneration surfaced.**

1. **The committed AD files were stale.** They predated `be9d6d1d` /
   `efed31cf` (switchable λ_θ clamp + `ReThetaT ≥ 20` floor), whose commit
   message already said `TAPENADE NEEDED before any adjoint use`. So the
   regen diff legitimately touches `saGammaRetheta_{d,b,fast_b}.f90` and
   `turbUtils_{d,b,fast_b}.f90`, not just `fluxes_*`. Those routines are now
   differentiated as the primal actually computes them.

2. **`use constants, only: …` breaks `_fast_b` for any routine that
   branches.** `autoEditReverseFast.py:222,228` rewrites Tapenade's
   `pushControl1b`/`popControl1b` into bare `myIntPtr`/`myIntStack`
   statements. Those symbols live in `constants`, and Tapenade propagates an
   `only:` list verbatim into the generated code — so a restricted import
   yields `Error: Symbol 'myintptr' at (1) has no IMPLICIT type` at compile
   time. Hit in `reThetaTCorrelation` (`turbUtils.F90:2299`), which gained an
   `if (Tu <= 1.3)` branch with the floor work. Fixed by dropping the
   `only:` there, with an in-source comment.

   **Rule of thumb: any hand-written routine that Tapenade will
   reverse-differentiate and that contains a branch must `use constants`
   unrestricted.** Branchless routines (`flengthCorrelation:2355`,
   `rethetacCorrelation:2379`) are fine with `only:`.

---

## Adjoint-solve diagnostic: psi-history / derivative-convergence reporting (2026-07-24)

Not part of the residual-differentiation inventory below (this touches the
**hand-written adjoint solver**, `adjointAPI.F90`/`adjointUtils.F90` — not
Tapenade-generated, so freely editable per `CLAUDE.md` rule 6). Opt-in via
`storePsiHistory` (bool, default `False`) — off by default, zero behavior
change unless explicitly enabled.

**What it does**: buffers the adjoint solution estimate (`psi`) every
`psiHistoryStep` KSP iterations (default 10, capped at `psiHistoryMax`
snapshots, default 50) during `solveAdjoint`'s `KSPSolve`, via a custom
PETSc `KSPMonitorSet` callback (`MyKSPMonitor`,
`src/adjoint/adjointUtils.F90`) that calls `KSPBuildSolution` mid-solve and
corrects it against the incoming guess (`psi_like1`). After the solve,
Python (`_printPsiHistorySensitivities` /
`_writePsiHistoryJSON`/`_printPsiHistoryTable` in `adflow/pyADflow.py`)
contracts each buffered snapshot through the existing
`computeJacobianVectorProductBwd` path (same machinery the converged-psi
total-derivative call already uses — no new AD code) to show how the total
derivative of each function evolves as the adjoint converges: a
norm-collapsed table to screen (`#`-banner style matching the standard
CL/CD convergence table) and the full per-DV-component values to
`<outputDirectory>/psiHistory_<func>.json` for plotting.

**Storage**: `ADjointPETSc` module gained `psiHistory(nState, snapshot)`,
`psiHistoryIters(snapshot)`, `psiHistoryResid(snapshot)`,
`psiHistoryCount` — allocated/reset once per `solveAdjoint` call (fresh
buffer per objective function, since `psi` means something different for
each). `inputADjoint` gained `storePsiHistory`, `psiHistoryStep`,
`psiHistoryMax`.

**Gotcha that cost ~3 hours to find**: none of those six new module
variables/arrays worked the first build — `storePsiHistory` read `True`
from Python right up to the last line before `solveAdjoint`'s Fortran body,
then read `False` *inside* that same subroutine, with no assignment
anywhere in the Fortran source that could explain it. Root cause: `.pyf`
option wiring is hand-maintained and module-scoped (see
`../architecture.md`'s "Known infra bug: `.pyf` option wiring" section) —
none of the six were listed in `adflow.pyf`'s `module inputadjoint` /
`module adjointpetsc` blocks, so Python was silently writing to phantom
attributes on the wrapper object instead of real Fortran memory. Fixed by
adding all six to their respective `.pyf` blocks. **Any future addition to
`ADjointPETSc`/`inputADjoint` (or any module) needs a matching `.pyf`
entry — it will not error, it will just silently no-op.**

Verification run + plots: `/home/mdo/Desktop/Run/MDO_PhD/Transition/gama_rethetha/adjoint_psihistory_test/` (tutorial-wing mesh, crossflow off, restart from converged state, DVs alpha/span/twist/shape). Finding: total derivatives stabilize (within 1%) after roughly 1 order of magnitude of KSP residual reduction, far before the `adjointL2Convergence` stopping criterion — case-specific (mild subsonic attached flow), not a general rule for how tight `adjointL2Convergence` needs to be on other meshes/regimes.

---

## Executive Inventory

| Category | SA Files/Touchpoints | SA-GR Files/Touchpoints | Status |
|----------|---------------------|------------------------|--------|
| Preprocessor guards | 6 #ifndef + 10 #ifdef in sa.F90 | 5 #ifndef + 10 #ifdef in saGammaRetheta.F90 | Both present |
| Tapenade directives | 10 II-LOOP | 10 II-LOOP | Symmetric |
| Generated AD files | sa_d.f90, sa_b.f90, sa_fast_b.f90 | saGammaRetheta_d.f90, saGammaRetheta_b.f90, saGammaRetheta_fast_b.f90 | All generated |
| Makefile head declarations | saSource, saViscous, saResScale | Source, Viscous, ResScale | All declared |
| masterRoutines wiring | 5 case blocks | 5 case blocks | Fully wired |
| frozenTurbulence | Generic nw/nwf logic | Uses same logic | Shared |
| Empirical correlations | N/A | 3 functions with _d versions | Differentiated |

---

## 1. Preprocessor Guards

### 1.A sa.F90 Guard Inventory

| Line | Guard | Context |
|------|-------|---------|
| 15 | `#ifndef USE_TAPENADE` | Start of qq Jacobian logic in saSource |
| 139 | `#ifdef TAPENADE_REVERSE` | II-LOOP for saSource main loop |
| 304 | `#ifndef USE_TAPENADE` | saViscous qq update |
| 335 | `#ifdef TAPENADE_REVERSE` | II-LOOP |
| 371 | `#ifdef TAPENADE_REVERSE` | II-LOOP |
| 445 | `#ifndef USE_TAPENADE` | Viscous Jacobian updates |
| 468, 478 | `#ifdef TAPENADE_REVERSE` | II-LOOP pair |
| 544, 643, 713 | `#ifndef USE_TAPENADE` | Additional Jacobian blocks |

### 1.B saGammaRetheta.F90 Guard Inventory

| Line | Guard | Context |
|------|-------|---------|
| 258 | `#ifndef USE_TAPENADE` | Source qq(1,1) diagonal update |
| 277 | `#ifdef TAPENADE_REVERSE` | II-LOOP for Source |
| 619 | `#ifndef USE_TAPENADE` | Source qq(2,2), qq(3,3) updates |
| 722 | `#ifdef TAPENADE_REVERSE` | II-LOOP |
| 793 | `#ifdef TAPENADE_REVERSE` | II-LOOP |
| 900 | `#ifndef USE_TAPENADE` | Viscous qq updates |
| 978, 986 | `#ifdef TAPENADE_REVERSE` | II-LOOP pair |
| 1091, 1169, 1177 | `#ifndef USE_TAPENADE` / `#ifdef TAPENADE_REVERSE` | i-direction viscous |
| 1281, 1359, 1385, 1405 | Mixed | j/k-direction viscous |

### 1.C turbUtils.F90 (23 guards total)

Key locations: lines 42, 120, 156, 230, 271, 325, 639, 683, 703, 730, 746, 784, 815, 888, 979, 1069, 1094, 1113, 1206, 1292, 1317, 1336, 1429, 1515, 1540, 1558, 1802, 1827, 2018, 2034, 2043.

### 1.D turbBCRoutines.F90 (4 guards)

Lines 4, 81, 784, 1169: all `#ifndef USE_TAPENADE` for BC Jacobian blocks.

### 1.E turbMod.F90

Line 43: `#ifndef USE_TAPENADE` for module-level allocatables.

### 1.F Pairing Table

| SA Guard (sa.F90) | SA-GR Analog (saGammaRetheta.F90) | Notes |
|-------------------|-----------------------------------|-------|
| L15 qq Jacobian init | L258 qq init | SA computes dfv1...; SA-GR zero-inits only |
| L139 II-LOOP | L277 II-LOOP | Symmetric |
| L304 saViscous qq | L900 Viscous qq | Same role |
| L335/371 II-LOOP | L722/793 II-LOOP | Symmetric |
| L445 Viscous Jac | L1091 i-Viscous Jac | SA-GR split by direction |
| L468/478 II-LOOP | L978/986 II-LOOP | Symmetric |
| L544/643/713 Jac blocks | L1169/1281/1359 | Split into i/j/k directions in SA-GR |

---

## 2. Tapenade In-Source Directives

### SA directives (sa.F90)
```fortran
! Line 140: !$AD II-LOOP
! Line 372: !$AD II-LOOP
! Line 479: !$AD II-LOOP
! Line 578: !$AD II-LOOP
! Line 692: !$AD II-LOOP
```

### SA-GR directives (saGammaRetheta.F90)
```fortran
! Line 278: !$AD II-LOOP
! Line 794: !$AD II-LOOP
! Line 987: !$AD II-LOOP
! Line 1178: !$AD II-LOOP
! Line 1386: !$AD II-LOOP
```

### turbUtils.F90 directives (comprehensive)
```fortran
! Lines: 43, 157, 272, 688, 735, 789, 879 (CHECKPOINT-START),
!        889, 923, 1011, 1102 (CHECKPOINT-END), 1111 (CHECKPOINT-START),
!        1114, 1148, 1234, 1325 (CHECKPOINT-END), 1334 (CHECKPOINT-START),
!        1337, 1371, 1457, 1548 (CHECKPOINT-END), 1807, 2023
```

---

## 3. Manually Differentiated Routines

**Search result**: No routines marked "MANUALLY DIFFERENTIATED" found in `src/`.

All differentiation is handled via Tapenade automatic differentiation.

---

## 4. Tapenade Makefile

### 4.A ALL_RES_FILES list includes both:
```makefile
$(SRC)/turbulence/sa.F90\
$(SRC)/turbulence/saGammaRetheta.F90\
```

### 4.B fullRoutines head declarations

**SA routines**:
```makefile
sa%saSource(w, rlv, vol, si, sj, sk, timeRef, d2wall) > \
           (w, rlv, vol, si, sj, sk, timeREf, scratch) \
sa%saViscous(w, vol, si, sj, sk, rlv, scratch) > \
            (w, vol, si, sj, sk, rlv, scratch) \
sa%saResScale(scratch, dw) > (dw) \
```

**SA-GR routines**:
```makefile
saGammaRetheta%Source(w, rlv, vol, si, sj, sk, timeRef, d2wall) > \
                     (w, rlv, vol, si, sj, sk, timeRef, scratch) \
saGammaRetheta%Viscous(w, vol, si, sj, sk, rlv, scratch) > \
                      (w, vol, si, sj, sk, rlv, scratch) \
saGammaRetheta%ResScale(scratch, dw) > (dw) \
```

### 4.C stateOnlyRoutines (fast_b mode)

**SA**:
```makefile
sa%saSource(w, rlv) > (w, rlv, scratch) \
sa%saViscous(w, rlv, scratch) > (w, rlv, scratch) \
sa%saResScale(scratch) > (dw) \
```

**SA-GR**:
```makefile
saGammaRetheta%Source(w, rlv) > (w, rlv, scratch) \
saGammaRetheta%Viscous(w, rlv, scratch) > (w, rlv, scratch) \
saGammaRetheta%ResScale(scratch) > (dw) \
```

### 4.D Pairing

| SA Declaration | SA-GR Declaration | Active Inputs | Active Outputs |
|----------------|-------------------|---------------|----------------|
| sa%saSource | saGammaRetheta%Source | w, rlv, vol, si, sj, sk, timeRef, d2wall | w, rlv, vol, si, sj, sk, timeRef, scratch |
| sa%saViscous | saGammaRetheta%Viscous | w, vol, si, sj, sk, rlv, scratch | w, vol, si, sj, sk, rlv, scratch |
| sa%saResScale | saGammaRetheta%ResScale | scratch, dw | dw |

---

## 5. Generated AD Output Files

### File sizes

| Directory | SA File | Size | SA-GR File | Size |
|-----------|---------|------|------------|------|
| outputForward | sa_d.f90 | 56KB | saGammaRetheta_d.f90 | 154KB |
| outputReverse | sa_b.f90 | 60KB | saGammaRetheta_b.f90 | 158KB |
| outputReverseFast | sa_fast_b.f90 | 50KB | saGammaRetheta_fast_b.f90 | 148KB |

### Subroutines in saGammaRetheta_d.f90:
- `sagammaretheta_block` (line 63)
- `source_d` (line 126)
- `source` (line 1096)
- `viscous_d` (line 1630)
- `viscous` (line 2331)
- `resscale_d` (line 2688)
- `resscale` (line 2732)
- `sagammarethetasolve` (line 2772)
- `evalsrcjacblock` (line 3455) - primal only
- `computesrclambda` (line 3819) - primal only

---

## 6. Master Routine Wiring

### master (lines 197-212)
```fortran
case (spalartAllmaras)
    allocate (qq(2:il, 2:jl, 2:kl))
    call saSource; call turbAdvection(..., qq); call saViscous; call saResScale
    deallocate (qq)
case (spalartallmarasnoft2gammaretheta)
    allocate (qqGR(2:il, 2:jl, 2:kl, 3, 3))
    call saGRSource; call turbAdvection(3, 3, ..., qqGR); call saGRViscous; call saGRResScale
    deallocate (qqGR)
```

### master_d (lines 544-555)
```fortran
case (spalartAllmaras)
    call saSource_d; call turbAdvection_d(1, 1, ..., qq); call saViscous_d; call saResScale_d
case (spalartallmarasnoft2gammaretheta)
    call saGRSource_d; call turbAdvection_d(3, 3, ..., qqGR); call saGRViscous_d; call saGRResScale_d
```

### master_b (lines 798-812)
```fortran
case (spalartAllmaras)
    call saResScale_b; call saViscous_b; call turbAdvection_b(1, 1, ..., qq); call saSource_b
case (spalartallmarasnoft2gammaretheta)
    call saGRResScale_b; call saGRViscous_b; call turbAdvection_b(3, 3, ..., qqGR); call saGRSource_b
```

### master_fast_b (lines 1160-1171)
```fortran
case (spalartAllmaras)
    call saResScale_fast_b; call saViscous_fast_b; call turbAdvection_fast_b(1, 1, ..., qq); call saSource_fast_b
case (spalartallmarasnoft2gammaretheta)
    call saGRResScale_fast_b; call saGRViscous_fast_b; call turbAdvection_fast_b(3, 3, ..., qqGR); call saGRSource_fast_b
```

### block_res_state_d (lines 1379-1390)
Same pattern as master_d, used for block-level state derivatives.

---

## 7. setupStateResidualMatrix

### Signature (adjointUtils.F90:7):
```fortran
subroutine setupStateResidualMatrix(matrix, useAD, usePC, useTranspose, &
                                    frozenTurb, ...)
```

### lStart/lEnd/nState logic (lines 87-106):
```fortran
if (turbOnly) then
    lStart = nt1; lEnd = nt2; nState = nt2 - nt1 + 1
else
    if (frozenTurb) then
        lStart = 1; lEnd = nwf; nState = nwf
    else
        lStart = 1; lEnd = nw; nState = nw
    end if
end if
```

For SA: `nState = nwf + 1` (5 flow + 1 turb = 6)  
For SA-GR: `nState = nwf + 3` (5 flow + 3 turb = 8)

### Call sites:
- `NKSolvers.F90:407` - NK solver PC
- `NKSolvers.F90:1983` - ANK solver PC
- `NKSolvers.F90:2378` - Turb-only PC

---

## 8. ADPC Flag Plumbing

### Declarations:
- `inputParam.F90:797`: `logical :: frozenTurbulence, viscPC, ADPC`
- `NKSolvers.F90:47`: `logical :: NK_ADPC`
- `NKSolvers.F90:1685`: `logical :: ANK_ADPC`

### Usage:
```fortran
! NKSolvers.F90:394
useAD = NK_ADPC

! NKSolvers.F90:1971, 2366
useAD = ANK_ADPC

! adjointAPI.F90:885
useAD = ADPC
```

---

## 9. frozenTurbulence Semantics

| Line | Context |
|------|---------|
| 106 | `use inputAdjoint, only: frozenTurbulence` |
| 138 | `if (frozenTurbulence) then` - skip turb res in forward |
| 152 | RANS check with frozen turb |
| 189, 215, 228 | Similar pattern for reverse |
| 623, 646 | Pass to `setupStateResidualMatrix` |
| 895 | Same for coarse level |
| 1140, 1165 | PC assembly with frozen turb check |

---

## 10. zeroADSeeds

### Zeroed arrays (adjointUtils.F90:920-1017):
```fortran
flowDomsd(...)%d2wall = zero
flowDomsd(...)%x = zero
flowDomsd(...)%si, sj, sk = zero
flowDomsd(...)%vol = zero
flowDomsd(...)%w, dw, fw = zero
flowDomsd(...)%scratch = zero  ! 5 slots only
flowDomsd(...)%p, gamma, aa = zero
flowDomsd(...)%rlv, rev = zero
flowDomsd(...)%ux..qz = zero
flowDomsd(...)%bmti1..bvtk2 = zero
flowDomsd(...)%BCData(mm)%... = zero
flowDomsd(...)%viscSubface(mm)%tau, q = zero
```

**transitionDebugd**: NOT present. `transitionDebug` is a primal-only diagnostic array with no derivative companion.

---

## 11. Halo Exchange, Surface Integration, getSolution

### whalo2 calls (masterRoutines.F90):
```fortran
! Line 159 (master):
call whalo2(currentLevel, 1_intType, nw, .True., .True., .True.)

! Line 504 (master_d):
call whalo2_d(1, 1, nw, .True., .True., .True.)

! Line 859 (master_b):
call whalo2_b(currentLevel, 1_intType, nw, .True., .True., .True.)

! Line 1205 (master_fast_b):
call whalo2_b(currentLevel, 1_intType, nw, .True., .True., .True.)
```

**nw is generic**: Uses `flowVarRefState%nw` which is 6 for SA, 8 for SA-GR. No hardcoding.

---

## 12. Boundary Condition AD

### turbBCRoutines.F90 guards:
- Line 4: `#ifndef USE_TAPENADE` - module-level
- Line 81: BC jacobian for applyAllTurbBCThisBlock
- Line 784: bcTurbTreatment jacobian
- Line 1169: Additional BC jacobian

### Generated AD files:
- `outputForward/turbBCRoutines_d.f90` (80KB)
- `outputReverse/turbBCRoutines_b.f90` (90KB)
- `outputReverseFast/turbBCRoutines_fast_b.f90` (41KB)

---

## Specific Items Checked

### smoothMinMax (turbUtils.F90:2371-2412)
- **Location**: `src/turbulence/turbUtils.F90:2371`
- **Differentiated versions**: `smoothminmax_d` present in `outputForward/turbUtils_d.f90`
- **Usage in SA-GR**: Imported and used throughout saGammaRetheta_d.f90

### Empirical Correlations
All three are differentiated:
- `reThetaTCorrelation` -> `rethetatcorrelation_d`
- `flengthCorrelation` -> `flengthcorrelation_d`
- `rethetacCorrelation` -> `rethetaccorrelation_d`

### evalSrcJacBlock and computeSrcLambda
- `saGammaRetheta.F90:2062` - `evalSrcJacBlock` defined
- `saGammaRetheta.F90:2303` - `computeSrcLambda` defined
- **Status**: Primal only, not differentiated (correct behavior)

### paramTurb.F90 Constants
All constants are `PARAMETER` - passive to Tapenade.

### transitionDebug Array
- **NO derivative companion** (`transitionDebugd` does not exist)
- Written only in primal-only paths

---

## Interesting Findings

1. **qq vs qqGR allocation**: SA uses `qq(2:il, 2:jl, 2:kl)` (rank-3), SA-GR uses `qqGR(2:il, 2:jl, 2:kl, 3, 3)` (rank-5 for 3x3 coupled Jacobian). Both modules export `qq` via renaming (`qqGR => qq`).

2. **turbAdvection call signature**: SA passes `1_intType, 1_intType` for nEq start/end; SA-GR passes `3_intType, 3_intType`.

3. **saGammaReThetaSolve**: Located after `ResScale` at saGammaRetheta.F90:1639. Contains iterative Newton solve, NOT inside `#ifndef USE_TAPENADE` guard but only called when `resOnly=.false.`, which never happens in AD paths.

4. **evalSrcJacBlock/computeSrcLambda NOT guarded**: These routines (saGammaRetheta.F90:2062, 2303) are visible to Tapenade. They appear in outputForward/saGammaRetheta_d.f90 as primal copies (lines 3455, 3819) but are NOT differentiated (no `_d` suffix). This is correct because they're not called from any differentiated head routine.

5. **Checkpoint directives in turbUtils**: 3 pairs of `!$AD CHECKPOINT-START/END` for memory management in reverse mode.

6. **transitionDebug has no derivative companion**: Declared at block.F90:667 as `real(...), pointer :: transitionDebug`. No `transitionDebugd` exists in block_d or anywhere. It is not zeroed in `zeroADSeeds`. This is intentional as it's diagnostic-only.

7. **No manually differentiated routines found**: Search for "MANUALLY DIFFERENTIATED" returned no results in src/. All AD is via Tapenade.

8. **Guard count verification**:
   - sa.F90: 6 `#ifndef USE_TAPENADE`, 10 `#ifdef TAPENADE_REVERSE`
   - saGammaRetheta.F90: 5 `#ifndef USE_TAPENADE`, 10 `#ifdef TAPENADE_REVERSE`

9. **Generated AD file sizes** (outputForward):
   - sa_d.f90: 56KB
   - saGammaRetheta_d.f90: 154KB (2.75x larger due to 3 transport equations)

10. **Module imports in masterRoutines.F90** (lines 22-24):
    ```fortran
    use sa, only: saSource, saViscous, saResScale, qq
    use saGammaRetheta, only: saGRSource => Source, saGRViscous => Viscous, &
                               saGRResScale => ResScale, qqGR => qq
    ```
    Both SA and SA-GR routines are imported with renaming for consistent calling.

11. **Makefile head declarations verified** (Makefile_tapenade:175-194):
    - SA: `sa%saSource`, `sa%saViscous`, `sa%saResScale`
    - SA-GR: `saGammaRetheta%Source`, `saGammaRetheta%Viscous`, `saGammaRetheta%ResScale`
    - Active variables are symmetric between both models.

12. **zeroADSeeds does NOT zero qq/qqGR derivatives**: The qq arrays are local allocatables in the solve routines, not module-level persistent state. Their derivatives are handled by Tapenade's push/pop mechanism.

---

## Appendix A: Files Touched

### Source files:
- `src/turbulence/sa.F90`
- `src/turbulence/saGammaRetheta.F90`
- `src/turbulence/turbUtils.F90`
- `src/turbulence/turbBCRoutines.F90`
- `src/turbulence/turbMod.F90`
- `src/adjoint/masterRoutines.F90`
- `src/adjoint/adjointUtils.F90`
- `src/adjoint/adjointAPI.F90`
- `src/adjoint/Makefile_tapenade`

### Generated AD:
- `src/adjoint/outputForward/`
- `src/adjoint/outputReverse/`
- `src/adjoint/outputReverseFast/`

---

## Appendix B: Cross-Reference Table

| SA Routine | SA-GR Routine | Forward | Reverse | Fast Reverse |
|------------|---------------|---------|---------|--------------|
| saSource | Source | saSource_d | saSource_b | saSource_fast_b |
| saViscous | Viscous | saViscous_d | saViscous_b | saViscous_fast_b |
| saResScale | ResScale | saResScale_d | saResScale_b | saResScale_fast_b |
| N/A | smoothMinMax | smoothminmax_d | (inlined) | (inlined) |
| N/A | reThetaTCorrelation | rethetatcorrelation_d | (inlined) | N/A |
| N/A | flengthCorrelation | flengthcorrelation_d | (inlined) | N/A |
| N/A | rethetacCorrelation | rethetaccorrelation_d | (inlined) | N/A |
| N/A | evalSrcJacBlock | (primal only) | N/A | N/A |
| N/A | computeSrcLambda | (primal only) | N/A | N/A |

---

## Appendix Z: Open Questions

1. **transitionDebug derivatives**: If transition debug values become optimization targets, this would need addressing.

2. **evalSrcJacBlock/computeSrcLambda**: Primal-only. If source-term dt restriction needs sensitivity, manual differentiation required.

3. **qq vs qqGR module aliasing**: Verify no name collisions in combined builds.

4. **unsteadyTurbTerm**: Called in `saGammaReTheta_block` (line 107) but NOT in any Makefile head declaration. The call is inside a block that is only executed in primal DADI solve (`resOnly=.false.`), not in the coupled master routines. The master routine calls to `unsteadyTurbTerm` are commented out (masterRoutines.F90:201, 547, 801, 1163).

5. **saGammaReTheta_block vs master pathway**: The `saGammaReTheta_block` routine exists for DADI decoupled solve, but the coupled ANK/NK path uses the individual `Source`, `Viscous`, `ResScale` routines via master. These are the ones differentiated. The `_block` routine appears in AD output as primal copy only.

6. **saEddyViscosity in stateOnlyRoutines**: Only `saEddyViscosity` appears in stateOnlyRoutines (Makefile_tapenade:322). The other eddy viscosity routines (`SSTEddyViscosity`, `kwEddyViscosity`, etc.) are not differentiated as they're for other turbulence models.

---

## Sanity-Check Log

- [x] All 5 master routine variants have SA-GR case blocks
- [x] Generated AD files exist for sa and saGammaRetheta in all 3 output directories
- [x] Makefile_tapenade declares both SA and SA-GR head routines
- [x] smoothMinMax has _d version (10 occurrences in turbUtils_d.f90)
- [x] All empirical correlations have _d versions
- [x] evalSrcJacBlock/computeSrcLambda NOT differentiated (correct - primal copies only)
- [x] paramTurb.F90 constants are PARAMETER
- [x] transitionDebug has no derivative companion (block.F90:667)
- [x] whalo2 calls use generic nw
- [x] turbBCRoutines has all 3 AD output files
- [x] sa.F90 guards: 6 `#ifndef USE_TAPENADE`, 10 `#ifdef TAPENADE_REVERSE`
- [x] saGammaRetheta.F90 guards: 5 `#ifndef USE_TAPENADE`, 10 `#ifdef TAPENADE_REVERSE`
- [x] setupStateResidualMatrix lStart/lEnd logic uses nw generically (lines 87-106)
- [x] saEddyViscosity is in stateOnlyRoutines (Makefile_tapenade:322)
- [x] unsteadyTurbTerm calls are commented out in masterRoutines.F90
