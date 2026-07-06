# ADflow Architecture, Internals & Transition Options

> What Claude needs to know about ADflow internals, user constraints, confirmed
> facts, and every runtime option added for the SA-γ-Re̅θt transition model.
> Physics equations live in [`paper-reference.md`](paper-reference.md); adjoint/AD
> touchpoints in [`adjoint-trace.md`](adjoint-trace.md).

---

# Part 1 — Architecture & Internals

## 1. Solver Architecture

### ANK (Approximate Newton-Krylov) — startup solver
- **Coupled mode**: flow + turbulence (all 8 vars) in one PETSc GMRES system
  - Matrix-free Jv captures γ↔ν̃ coupling through residual evaluation
  - Preconditioner: approximate first-order Jacobian (block-ILU)
  - CFL ramps from ANKCFL0 to ANKCFLLimit
- **Decoupled mode**: flow solved first, then turbulence separately
  - Turbulence sub-solve options:
    - **DADI** (DD-ADI block solver): calls `saGammaReThetaSolve`. Coupling set by
      `TurbDADICoupled` (see Part 2):
      - `"full"`: 3×3 coupled block (default)
      - `"transition"`: SA scalar solve + γ-Re̅θt 2×2 block
      - `"decoupled"`: 3 independent scalar solves
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

> Values are non-dimensional. ADflow uses **pressure–density (p-ρ) scaling**, so
> velocity normalizes to M·√γ (not 1) and viscosities are stored as ratios to
> μ_∞. See [`nondimensionalization.md`](nondimensionalization.md).

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

---

# Part 2 — Transition Options Reference

These options do not exist in upstream ADflow. All were added on this branch.

## Python Options

### New transition-specific options

| Option | Type | Default | What it does |
|---|---|---|---|
| `"transitionFirstOrderUpwind"` | bool | `True` | First-order upwind for γ and Re̅θt convection. More dissipative but more robust. |
| `"transitionSrcDtRestrict"` | bool | `True` | Enable source-term dt restriction (P&Z Eq. 59). Caps λ_source × dt ≤ 0.9. |
| `"transitionSrcDtLimit"` | float | `0.9` | Threshold for source-term dt restriction (λ_source × dt ≤ this value). |
| `"srcDtDeactivateIters"` | int | `5` | Deactivate source-dt restriction after N consecutive ANK iterations without backtracking (P&Z §IV.B.3). Set to 0 to never deactivate. |
| `"TurbDADICoupled"` | str | `"full"` | DADI coupling mode: `"decoupled"` (3 scalar solves), `"transition"` (SA alone + γ-Re̅θt 2×2), `"full"` (3×3 block). |
| `"turbResScale"` | list/None | `None` (auto) | Residual scaling per equation. Auto-set to `[10000, 10, 10000]` for this model. Override to tune convergence balance. |
| `"transitionDampTheta"` | float | `0.99` | Back-off factor for iterative γ/Re̅θt update damping in DD-ADI (P&Z §3). |
| `"transitionDampMaxIter"` | int | `40` | Max back-off iterations for γ/Re̅θt bounds enforcement in DD-ADI. |

### Turb-ANK KSP physicality options (transition-specific)

| Option | Type | Default | What it does |
|---|---|---|---|
| `"ANKPhysicalLSTolReTheta"` | float | `0.99` | Relative physicality tolerance for Re̅θt in Turb-ANK (replaces `ANKPhysicalLSTolTurb` for Re̅θt). |
| `"omegaMinGamma"` | float | `0.05` | Minimum step factor floor for γ. Prevents collapse in laminar regions where γ→0. |

### Existing ADflow options relevant to turbulent solver path

| Option | Type | Default | What it does |
|---|---|---|---|
| `"ANKUseTurbDADI"` | bool | `True` | `True` = DADI for turbulence. `False` = Turb-ANK KSP (Newton-Krylov). |
| `"ANKNSubiterTurb"` | int | `1` | Inner turbulence iterations per outer ANK step. |
| `"ANKTurbCFLScale"` | float | `1.0` | CFL multiplier for turb equations relative to flow. |
| `"ANKTurbKSPDebug"` | bool | `False` | Print linear residual, KSP iters, step size each Turb-ANK iteration. |
| `"ANKPhysicalLSTolTurb"` | float | `0.99` | Physicality line-search tolerance for ν̃ in Turb-ANK (γ uses absolute bounds instead). |

## Examples

### 1. Robust startup (recommended defaults)

```python
solverOptions = {
    # Transition-specific (new)
    "transitionFirstOrderUpwind": True,      # robust convection for γ, Re̅θt
    "transitionSrcDtRestrict": True,         # source limiting ON
    "srcDtDeactivateIters": 5,               # deactivate after 5 clean ANK iters
    "TurbDADICoupled": "full",               # 3×3 coupled DADI
    # turbResScale auto-set to [10000, 10, 10000]

    # Solver path (existing ADflow)
    "ANKUseTurbDADI": True,                  # use DADI for turbulence
}
```

### 2. Accuracy run (restarting from converged solution)

```python
solverOptions = {
    "transitionFirstOrderUpwind": False,  # second-order convection (sharper transition front)
    "transitionSrcDtRestrict": False,     # no source limiting (solution already stable)
    "TurbDADICoupled": "full",
    "ANKUseTurbDADI": True,
}
```

### 3. Debugging convergence — try decoupled DADI

```python
solverOptions = {
    "transitionFirstOrderUpwind": True,
    "transitionSrcDtRestrict": True,
    "TurbDADICoupled": "decoupled",       # simplest: 3 independent scalar solves
    "ANKUseTurbDADI": True,
}
```

### 4. Debugging convergence — try partial coupling

```python
solverOptions = {
    "transitionFirstOrderUpwind": True,
    "transitionSrcDtRestrict": True,
    "TurbDADICoupled": "transition",      # SA alone, γ-Re̅θt coupled as 2×2 block
    "ANKUseTurbDADI": True,
}
```

### 5. Turb-ANK KSP path (Newton-Krylov for turbulence)

```python
solverOptions = {
    "transitionFirstOrderUpwind": True,
    "transitionSrcDtRestrict": True,
    # TurbDADICoupled ignored when ANKUseTurbDADI=False

    "ANKUseTurbDADI": False,              # switch to Turb-ANK KSP
    "ANKNSubiterTurb": 3,                 # more inner Newton iters
    "ANKTurbCFLScale": 0.5,              # lower CFL for turb if unstable
    "ANKTurbKSPDebug": True,             # print convergence info
}
```

### 6. Custom residual scaling

```python
solverOptions = {
    "turbResScale": [5000.0, 1.0, 5000.0],  # [ν̃, γ, Re̅θt] — lower γ scaling
    "TurbDADICoupled": "full",
    "ANKUseTurbDADI": True,
}
```

## Solver Path Summary

```
ANKUseTurbDADI = True ──┬── TurbDADICoupled = "decoupled"   → 3 scalar solves
                        ├── TurbDADICoupled = "transition"  → SA scalar + γ-Re̅θt 2×2 block
                        └── TurbDADICoupled = "full"        → 3×3 coupled block (default)

ANKUseTurbDADI = False ──── Turb-ANK KSP (Newton-Krylov, GMRES)
```

## Source-Term Eigenvalue Control (P&Z 2020, Eq. 59)

The source-term dt restriction prevents unbounded solution updates by limiting:

```
λ_source × Δt ≤ transitionSrcDtLimit  (default 0.9)
```

where `λ_source` is the **largest positive eigenvalue** of the 3×3 source-term Jacobian:

```
            ⎡ ∂S_ν̃/∂ν̃      ∂S_ν̃/∂γ      0           ⎤
A_source =  ⎢ ∂S_γ/∂ν̃      ∂S_γ/∂γ      ∂S_γ/∂Re̅θt   ⎥
            ⎣ 0            0            ∂S_Re̅θt/∂Re̅θt ⎦
```

### Key points

1. **Block-triangular structure**: A13=A31=A32=0 (P&Z §7.1), so eigenvalues are computed exactly without a cubic solver:
   - λ₃ = A33 (Re̅θt diagonal)
   - λ₁,₂ from 2×2 block [A11,A12; A21,A22] via quadratic formula
   - `λ_source = max(0, λ₁, λ₂, λ₃)`

2. **Independent of `TurbDADICoupled` mode** — coupling mode only affects how DADI solves the system, not eigenvalue computation.

3. **Auto-deactivation**: After `srcDtDeactivateIters` consecutive ANK iterations without backtracking, the restriction turns off. Reactivates on backtracking.

### Examples

#### 7. Conservative eigenvalue control (stiff cases)

```python
solverOptions = {
    "transitionSrcDtRestrict": True,
    "transitionSrcDtLimit": 0.7,             # stricter than default 0.9
    "srcDtDeactivateIters": 10,              # wait longer before deactivating
}
```

#### 8. Debug eigenvalue issues

```python
solverOptions = {
    "transitionSrcDtRestrict": True,
    "srcDtDeactivateIters": 0,               # NEVER deactivate — always restrict
    "ANKTurbKSPDebug": True,                 # print iteration info
}
```

#### 9. Disable source-dt restriction entirely

```python
solverOptions = {
    "transitionSrcDtRestrict": False,        # no eigenvalue computation, no restriction
}
```

## γ Physicality Check in Turb-ANK (Redesigned)

In the Turb-ANK KSP path (`ANKUseTurbDADI = False`), γ uses **absolute bound enforcement** instead of a relative tolerance:

- **Full step allowed** if result stays in [gammaLo, gammaHi] (~[1e-10, 2.0])
- **Only reduced** when full step would violate bounds
- **`omegaMinGamma`** (default 0.05) prevents step collapse in laminar regions where γ→0

This differs from ν̃ and Re̅θt which use relative tolerances (`ANKPhysicalLSTolTurb`, `ANKPhysicalLSTolReTheta`).

**Why**: In laminar flow, γ≈0. The old relative check `ratio = γ/update × tol` collapses to near-zero, killing the transition front before it can develop.

## Internal State (not user-settable)

These are managed automatically by the solver when `transitionSrcDtRestrict = True`:

- `srcDtRestrictActive`: starts `True`, flips to `False` after `srcDtDeactivateIters` consecutive no-backtrack ANK iterations. Resets to `True` on backtrack.
- `noBacktrackCount`: counter driving the above.
