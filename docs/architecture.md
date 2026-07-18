# ADflow Architecture, Internals & Transition Options

> What Claude needs to know about ADflow internals, user constraints, confirmed
> facts, and every runtime option added for the SA-γ-Re̅θt transition model.
> Physics equations live in the full paper, [`SA_GAMMA_RETHETHA_BASE/Piotrowski_Zingg_2020_SA-sLM2015_clean (1).md`](SA_GAMMA_RETHETHA_BASE/Piotrowski_Zingg_2020_SA-sLM2015_clean%20(1).md);
> non-dim conventions in [`nondimensionalization.md`](nondimensionalization.md); adjoint/AD
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
| Main transition model | `src/turbulence/saGammaRetheta.F90` (2470 lines) |
| Smooth helper functions (correlations + `smoothMinMax`) | `src/turbulence/turbUtils.F90:2279-2412` (`reThetaTCorrelation`, `flengthCorrelation`, `rethetacCorrelation`, `smoothMinMax`) |
| Initialization | `src/initFlow/initializeFlow.F90:140-146, 2229-2241` |
| Wall/farfield BCs | `src/turbulence/turbBCRoutines.F90:438-465` (farfield: padrão ADflow, sem caso especial), `888-960` (wall, caso SA-GR) |
| Dispatch (turbAPI) | `src/turbulence/turbAPI.F90:49,74` |
| ANK/NK variable bounds | `src/NKSolver/NKSolvers.F90:3191,3359` |
| Preconditioner | `src/NKSolver/blockette.F90:815-816` |
| AD forward | `src/adjoint/outputForward/saGammaRetheta_d.f90` |
| AD reverse | `src/adjoint/outputReverse/saGammaRetheta_b.f90` |
| AD reverse fast | `src/adjoint/outputReverseFast/saGammaRetheta_fast_b.f90` |

---

## 4. Key Code Patterns

### Source-term assembly
In `saGammaRetheta.F90`, subroutine `saGammaReTheta_block(resOnly)` (line 67):
- `resOnly = .true.`: compute residual only, don't update w
- `resOnly = .false.`: compute residual + run DADI solver + update w

Source routine (subroutine `Source`, line 146) computes:
1. SA terms (ν̃): term1, term2_prod, term2_dest → with γ multiplier on production
2. γ terms: P_γ, E_γ via F_onset, F_turb, vorticity
3. Re̅θt terms: P_θt via timeScale, F_θt, Re_θt correlation

### DADI solver
`saGammaReThetaSolve` (lines 1485-2131):
- 3×3 block DD-ADI in i,j,k directions
- Uses qq(i,j,k,row,col) matrix from Source routine
- Solution damping (Algorithm 2, per-variable exponential back-off with
  warned last-resort clip) in the update section after the tri-diagonal
  solves (search `transitionDampMaxIter`)
- Row/column scaling at lines 1582-1584 (using turbResScale)

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
| Wall BC for γ | zero-gradient (Neumann), Re̅θt=zero-gradient — `turbBCRoutines.F90:931` (`bmt=-1`; commit dc1950ef) |
| Roughness | Implemented via crossflow D_scf — `transitionRoughnessHeight` (default 3.3e-6 m) |

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
| `"srcDtDeactivateIters"` | int | `5` | Deactivate source-dt restriction after N consecutive clean (no-backtrack) turbKSP iterations **in the second-order regime** (`totalR ≤ ANKSecondOrdSwitchTol·totalR0`, the inexact-Newton analog of P&Z §IV.B.3). Counter resets when backtracking is triggered (even if it succeeds) or when the residual rises back above the switch tolerance. With the default `ANKSecondOrdSwitchTol = 1e-16` the regime is never entered ⇒ restriction never deactivates; set it to ~`1e-5` (paper's phase-switch value) to enable the acceleration. `0` = restriction inactive in turbKSP from the start (**not** "never deactivate"). DADI ignores this option (restriction always on there). Semantics fixed 2026-07-07 (D-A2-3). |
| `"TurbDADICoupled"` | str | `"full"` | DADI coupling mode: `"decoupled"` (3 scalar solves), `"transition"` (SA alone + γ-Re̅θt 2×2), `"full"` (3×3 block). |
| `"turbResScale"` | list/None | `None` (auto) | Residual scaling per equation. Auto-set to `[10000, 10, 10000]` for this model. Override to tune convergence balance. |
| `"transitionDampTheta"` | float | `0.99` | Back-off factor for per-variable γ/Re̅θt update damping in DD-ADI (P&Z Algorithm 2). |
| `"transitionDampMaxIter"` | int | `10000` | Safety cap on the back-off loop (unbounded in the paper); 10000 ⇒ effectively unbounded (0.99¹⁰⁰⁰⁰ ≈ 0). A hard clip to the bounds remains as **last-resort fallback only**: it can only fire after the loop exhausts, which requires the previous state to already be out of bounds; when it fires, a warning with cell counts prints advising to raise this option or investigate the upstream bound violation. Changed from 40 on 2026-07-07 (D-A2-5). |
| `"transitionCrossflow"` | bool | `True` | Enable the helicity-based crossflow source D_scf (P&Z Eq. 15-26) on the Re̅θt equation. Harmless in 2D (D_scf≡0); enable for swept/3D. |
| `"transitionRoughnessHeight"` | float | `3.3e-6` | Surface roughness height h for the crossflow correlation (Eq. 17), as a physical length in mesh units (metres). 3.3e-6 = 3.3 µm (smooth surface). |
| `"transitionRefLength"` | float | `-1.0` (auto) | Reference length l [mesh units] in the vorticity limiter (P&Z Eqs. 52-53; paper uses root chord — the physical cap scales as 1/√l, a calibration scale, NOT a unit conversion, so the "drop Re" rule of `nondimensionalization.md` does not apply). Negative = auto: uses the AeroProblem `chordRef` (via `inputPhysics%lengthRef`, refreshed at every `setAeroProblem`). Set explicitly to decouple from chordRef; `1.0` recovers the pre-option behavior (l = 1 m). Added 2026-07-07 to close finding D1. |
| `"transitionNK"` | bool | `True` | Master switch for the NK/ANK/turbKSP column-scaling + Eq. 59 bundle (incl. NK reactivation-on-backtrack, Algorithm 2 in NK) — 2026-07-16. Default preserves existing behavior; still additionally gated on `turbModel==SA-Gamma-Retheta` everywhere. |
| `"transitionNKAutoDisableTol"` | float | `0.0` | One-way latch, NK phase only: once the Newton residual norm drops below this fraction of `totalR0` (the **freestream** reference residual from `getFreeStreamResidual`, `solvers.F90:972` — NOT the restart-point residual; e.g. ~8.91e7 on `3D_Plain_Wing`), the `transitionNK` bundle (column scaling, Algorithm 2 damping, Eq. 59 reactivation) is turned off for the rest of the NK phase — i.e. "fall back to native NK." Default `0.0` never trips (unchanged behavior). **Tested 2026-07-18, `nk_switch_crossing_test`, and found unsafe at any point in NK**: tripping it either at NK engagement or ~10 outer iterations later (deep past engagement, residual already down 2+ orders) produces the identical catastrophic blowup (nuturb res → O(1e3), totalRes → O(1e9)) both times. Column scaling is load-bearing for the *entire* NK phase for this model — the 13-orders-of-magnitude state spread (ν̃, γ, Re̅θt) doesn't shrink with the residual, so "native NK" is never safe to fall back to. Kept as a diagnostic knob, not a recommended option. |
| `"transitionRowVolScale"` | bool | `False` | Eq. 58 (P&Z) geometric row-scaling factor on NK's residual rows (`volRef**(5/3)` flow, `volRef**(2/3)` turb, on top of existing `turbResScale`). **Off by default — genuinely tested 2026-07-16 and found to stall NK's linear solve** (lin res pinned ~1.0) on the 3D_Plain_Wing case; see `adflow-vs-paper-solver.md` §5. Not recommended until the volRef-vs-paper's-J correspondence is revisited. |
| `"transitionNKStallStepTol"` / `"transitionNKStallCountTrigger"` / `"transitionNKStallRtolCap"` | float / int / float | `0.1` / `3` / `1.0` | NK stall escape (`nk_switch_crossing_test`, 2026-07-18). Root cause: `getEWTol` (`NKSolvers.F90:2174`, standard PETSc EW) computes the Krylov `rtol` as `(norm/oldNorm)^1.618` — when the Newton step is pinned (`Step`~0), `norm≈oldNorm` so the ratio→1 and `rtol` rises to its 0.8 cap, i.e. EW picks the *loosest* linear solve exactly when stalled. Once `Step` has been below `transitionNKStallStepTol` for `transitionNKStallCountTrigger` consecutive NK iterations (internal counter `nkStallCount`, reset each NK entry), `rtol` is forced down to `transitionNKStallRtolCap` for that iteration instead of trusting EW. `transitionNKStallRtolCap=1.0` (default) disables this (never caps). **Validated 2026-07-18** at `transitionNKStallRtolCap=0.05` on the organic ANK→NK crossing (`run_stall_fix_organic.log` vs `run.log.baseline_stall`): baseline pins at scaledTotalRes~1.7473e-3 and then *drifts back up* to ~1.79e-3 over 300+ iterations with zero net progress; the capped run breaks through the same floor (first drop to 1.522e-3 the instant the cap first engages, iter 35) and keeps inching down via periodic tight-solve "notches" every ~5-6 iterations (reached 1.4772e-3 by iter 308). **Partial fix only** — converts a total freeze into a slow crawl, not fast convergence; the underlying cause of why the Newton direction is poor at this residual level (likely a specific cell near the laminar-turbulent front) is still open. |
| `"transitionResidualAutoscale"` | bool | `False` | Eq. 58 (P&Z) S_a residual-autoscaling proxy — periodically (every NK Jacobian reform) rescales each turbulence variable's row to match the mean-flow block's current residual norm. The paper gives no formula (cites Osusky & Zingg's thesis, unavailable here); this is a same-intent proxy, not verified identical. Off by default — tested 2026-07-16: no stall, real progress, but noisier/smaller steps than baseline; marginal, not clearly better. |

**`transitionRefLength` plumbing & guidance.** No new AeroProblem wiring was
added: `pyADflow.py` already pushes `ap.chordRef` into `inputPhysics%lengthRef`
inside `setAeroProblemData` (~line 3608), and every compute entry point calls
`setAeroProblem` *before* any residual evaluation, so the fallback is always
fresh (pyADflow errors out if `chordRef` is missing; `setDefaultValues` seeds
`lengthRef = 1.0` at init as a safety net). Fortran reads the option in
`saGammaRetheta.F90` (`Source` + `evalSrcJacBlock`): `transitionRefLength > 0`
wins, else `lengthRef`. **Full aircraft:** no single l is "correct" for all
components (cap ∝ 1/√c_local); the paper's own prescription (one global root
chord) is already a compromise — use the MAC or root chord, and remember the
limiter is a numerical safety net whose failure mode (residual oscillation near
LSBs) is visible, not a silent physics error.

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
    "srcDtDeactivateIters": 5,               # deactivate after 5 clean 2nd-order iters
    # (deactivation only engages if ANKSecondOrdSwitchTol is set, e.g. 1e-5)
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

3. **Auto-deactivation** (turbKSP only; P&Z §IV.B.3): after `srcDtDeactivateIters` consecutive clean iterations in the second-order regime (`totalR ≤ ANKSecondOrdSwitchTol·totalR0` — the inexact-Newton analog), the restriction turns off. The counter resets (restriction reactivates) when backtracking is triggered — even if the backtrack succeeds — or when the residual rises back above the switch tolerance. With the default `ANKSecondOrdSwitchTol = 1e-16`, deactivation never engages. DADI has no deactivation: the restriction stays on, matching the paper's approximate-Newton (globalization) phase where it is never deactivated.

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
    # leave ANKSecondOrdSwitchTol at its default (1e-16): the restriction
    # then never deactivates. (Do NOT use srcDtDeactivateIters = 0 for this —
    # 0 means the restriction is inactive in turbKSP from the start.)
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

- `srcDtRestrictActive`: starts `True`, flips to `False` after `srcDtDeactivateIters` consecutive clean turbKSP iterations in the second-order regime. Returns to `True` when backtracking is triggered or when `totalR` rises back above `ANKSecondOrdSwitchTol·totalR0`.
- `noBacktrackCount`: counter driving the above (module variable, `inputParam.F90`; persists across solves — self-corrects via the residual condition on the first iteration of each solve).

## Monitor variable: `"scaledtotalr"` (2026-07-16)

Add `"scaledtotalr"` to `monitorVariables` to print an additional column
showing the Eq. 58 S_r/S_a-scaled residual (`sumAllResidualsScaled`,
`src/utils/utils.F90`) alongside the unchanged `"totalr"` column. Purely
for visibility — it does **not** feed `totalR`/switch tolerances the way
`"totalr"` does (see the `'scaledtotalR'` case in `solvers.F90`, which
computes but never assigns the module-level `totalR`). Reflects whatever
`transitionRowVolScale`/`transitionResidualAutoscale` are currently set to;
identical to `"totalr"` when both are off.

## Known infra bug: `.pyf` option wiring (2026-07-16, not yet fixed for all options)

`src/f2py/adflow.pyf` is a **hand-maintained** f2py interface file, not
auto-regenerated from Fortran source on every build. Any `inputIteration`
module variable not explicitly listed there is invisible to real
Python↔Fortran communication — but f2py's `fortran`-type Python objects
silently accept **arbitrary attribute names** with no backing memory, so
`setOption`/reading the option back gives no error and no warning. This
branch's `.pyf` was stale enough to affect `transitionSrcDtRestrict`,
`transitionSrcDtLimit`, `srcDtDeactivateIters`, `transitionDampTheta`,
`transitionDampMaxIter` (added before this session — **still not fixed**,
out of this session's scope) and `transitionNK`/`transitionRowVolScale`/
`transitionResidualAutoscale` (added this session — **fixed**, entries
added to `module inputiteration` in `adflow.pyf`). `transitionSrcDtEigMode`
has no backing Fortran variable at all (separate, also pre-existing bug).
**Before trusting any Python-side setting of an `inputIteration` option
that isn't already validated working, grep `adflow.pyf`'s `module
inputiteration` block for it first** — if it's missing, the option is a
no-op from Python regardless of what `setOption`/read-back suggest.
