# ADflow Architecture & User Notes

> What Claude needs to know about ADflow internals, user constraints, and confirmed facts.

---

## 1. Solver Architecture

### ANK (Approximate Newton-Krylov) — startup solver
- **Coupled mode**: flow + turbulence (all 8 vars) in one PETSc GMRES system
  - Matrix-free Jv captures γ↔ν̃ coupling through residual evaluation
  - Preconditioner: approximate first-order Jacobian (block-ILU)
  - CFL ramps from ANKCFL0 to ANKCFLLimit
- **Decoupled mode**: flow solved first, then turbulence separately
  - Turbulence sub-solve options:
    - **DADI** (DD-ADI block solver): calls `saGammaReThetaSolve`
      - `TurbDADICoupled=True`: should use full 3×3 block → BUG: off-diags zero
      - `TurbDADICoupled=False`: diagonal only
    - **Turb-ANK**: separate ANK for turbulence equations only

### NK (Newton-Krylov) — terminal solver
- Fully coupled: matrix-free Jv with exact AD Jacobian
- Cubic line search, Eisenstat-Walker tolerance
- Activated when residual drops below NKSwitchTol

### Multigrid (RK/D3ADI smoother) — NOT USED by this user
- Skip T1.6 (multigrid restriction)

---

## 2. State-Vector Layout

```
w(i,j,k, 1)   = ρ
w(i,j,k, 2)   = ρu
w(i,j,k, 3)   = ρv
w(i,j,k, 4)   = ρw
w(i,j,k, 5)   = ρE
w(i,j,k, itu1) = ρν̃  (SA working variable)
w(i,j,k, itu2) = γ   (intermittency) — NEW
w(i,j,k, itu3) = Re̅θt (transition onset Re) — NEW
```

All stored as conservative (ρ·φ). Generic nVar extension handles sizing.

---

## 3. Key Module Locations

| What | Where |
|------|-------|
| Turbulence model enum | `src/modules/constants.F90:128` |
| Model constants (ca1,ca2,...) | `src/modules/paramTurb.F90:32-52` |
| Input parameters | `src/modules/inputParam.F90:293,298` |
| Block data (transitionDebug array) | `src/modules/block.F90:662`, `blockPointers.F90:156` |
| Main transition model | `src/turbulence/saGammaRetheta.F90` (1862 lines) |
| Smooth helper functions | `src/turbulence/saGammaRethetaHelpers.F90` (367 lines) |
| Initialization | `src/initFlow/initializeFlow.F90:140-146, 2229-2241` |
| Wall/farfield BCs | `src/turbulence/turbBCRoutines.F90:441-470, 921-983` |
| Dispatch (turbAPI) | `src/turbulence/turbAPI.F90:49,74` |
| ANK/NK variable bounds | `src/NKSolver/NKSolvers.F90:3191,3359` |
| Preconditioner | `src/NKSolver/blockette.F90:815-816` |
| AD forward | `src/adjoint/outputForward/saGammaRetheta_d.f90` |
| AD reverse | `src/adjoint/outputReverse/saGammaRetheta_b.f90` |
| AD reverse fast | `src/adjoint/outputReverseFast/saGammaRetheta_fast_b.f90` |

---

## 4. Key Code Patterns

### Source-term assembly
In `saGammaRetheta.F90`, subroutine `saGammaRetheta_block(calledFromANK)`:
- `calledFromANK = .true.`: compute residual only, don't update w
- `calledFromANK = .false.`: compute residual + run DADI solver + update w

Source routine at line ~300 computes:
1. SA terms (ν̃): term1, term2_prod, term2_dest → with γ multiplier on production
2. γ terms: P_γ, E_γ via F_onset, F_turb, vorticity
3. Re̅θt terms: P_θt via timeScale, F_θt, Re_θt correlation

### DADI solver
`saGammaReThetaSolve` (lines 1251-1861):
- 3×3 block DD-ADI in i,j,k directions
- Uses qq(i,j,k,row,col) matrix from Source routine
- Solution damping (Algorithm 2) at lines 1830-1856
- Row/column scaling at lines 1329-1331 (using turbResScale)

### Variable references
- `rlv(i,j,k)` = μ/μ_∞ (laminar viscosity ratio, dimensionless)
- `rev(i,j,k)` = μ_t/μ_∞ (eddy viscosity ratio, dimensionless)
- `d2Wall(i,j,k)` = wall distance (pre-computed)
- `si/sj/sk(i,j,k,1:3)` = face normals
- `vol(i,j,k)` = cell volume

---

## 5. User Constraints

- "No in-place modifications to SA model" — transition is a modifier, not replacement
- "All coupling strategies selectable at runtime, test all, choose best"
- "SST exists → 2-eq turb infrastructure exists → mirror it"
- "ANK must work (more robust for complex geometries)"
- "Don't care about multigrid"
- "Don't care about crossflow for now" (helper exists if needed later)
- "Adjoint: freeze γ-Re̅θt linearization at this stage"
- Both sLM2015 and LM2015 F_turb forms must be runtime-selectable

---

## 6. Confirmed Answers to Open Questions

| Question | Answer |
|----------|--------|
| State-vector layout | itu1=ν̃, itu2=γ, itu3=Re̅θt, generic nVar extension |
| SST pattern | SST uses itu1=k, itu2=ω; our model follows same pattern |
| ANK modes | Coupled OR decoupled; turb solved with DADI or turb-ANK |
| Wall distance | `d2Wall(i,j,k)` pre-computed, available everywhere |
| Metrics | `si/sj/sk(i,j,k,1:3)`, `vol(i,j,k)` in blockPointers |
| Residual storage | `scratch(i,j,k,idvt+n)` → scaled to `dw(i,j,k,itun)` |
| LAPACK | Available (linked in build system) |
| Tu_∞ | `turbIntensityInf` exists in inputParam.F90:591 |
| Wall BC for γ | γ=0 (Dirichlet), Re̅θt=zero-gradient |
| Roughness | Not implemented yet, helper exists |
