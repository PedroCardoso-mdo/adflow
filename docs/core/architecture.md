# ADflow Architecture, Internals & Transition Options

> What Claude needs to know about ADflow internals, user constraints, confirmed
> facts, and every runtime option added for the SA-γ-Re̅θt transition model.
> Physics equations live in the full paper, [`SA_GAMMA_RETHETHA_BASE/Piotrowski_Zingg_2020_SA-sLM2015_clean (1).md`](../SA_GAMMA_RETHETHA_BASE/Piotrowski_Zingg_2020_SA-sLM2015_clean%20(1).md);
> non-dim conventions in [`nondimensionalization.md`](nondimensionalization.md); adjoint/AD
> touchpoints in [`VERIFICATION/adjoint-trace.md`](../VERIFICATION/adjoint-trace.md).

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
| Input parameters | `src/modules/inputParam.F90` (transition options start ~L327) |
| Block data (transitionDebug array) | `src/modules/block.F90:662`, `blockPointers.F90:156` |
| Main transition model | `src/turbulence/saGammaRetheta.F90` (~2800 lines) |
| Smooth helper functions (correlations + `smoothMinMax`) | `src/turbulence/turbUtils.F90:~2290-2438` (`reThetaTCorrelation` 2290, `flengthCorrelation` 2355, `rethetacCorrelation` 2380, `smoothMinMax` 2397) |
| Initialization | `src/initFlow/initializeFlow.F90:140-146, 2229-2241` |
| Wall/farfield BCs | `src/turbulence/turbBCRoutines.F90:438-465` (farfield: padrão ADflow, sem caso especial), `888-960` (wall, caso SA-GR) |
| Dispatch (turbAPI) | `src/turbulence/turbAPI.F90:49,74` |
| ANK/NK variable bounds | `src/NKSolver/NKSolvers.F90` (`physicalityCheckANKTurb` ~4316, `applyNKAlgorithm2Damping` ~1654) |
| Preconditioner | `src/NKSolver/blockette.F90:815-816` |
| AD forward | `src/adjoint/outputForward/saGammaRetheta_d.f90` |
| AD reverse | `src/adjoint/outputReverse/saGammaRetheta_b.f90` |
| AD reverse fast | `src/adjoint/outputReverseFast/saGammaRetheta_fast_b.f90` |

---

## 4. Key Code Patterns

### Source-term assembly
In `saGammaRetheta.F90`, subroutine `saGammaReTheta_block(resOnly)` (~line 70):
- `resOnly = .true.`: compute residual only, don't update w
- `resOnly = .false.`: compute residual + run DADI solver + update w

Source routine (subroutine `Source`, ~line 212) computes:
1. SA terms (ν̃): term1, term2_prod, term2_dest → with γ multiplier on production
2. γ terms: P_γ, E_γ via F_onset, F_turb, vorticity
3. Re̅θt terms: P_θt via timeScale, F_θt, Re_θt correlation

### DADI solver
`saGammaReThetaSolve` (~lines 1729-2405):
- 3×3 block DD-ADI in i,j,k directions
- Uses qq(i,j,k,row,col) matrix from Source routine
- Solution damping (Algorithm 2, per-variable exponential back-off with
  warned last-resort clip) in the update section after the tri-diagonal
  solves (search `transitionDampMaxIter`)
- Row/column scaling using turbResScale (~lines 1827-1831)

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
- Crossflow transition (D_scf, P&Z Eqs. 15-26): the source code stays live in
  residual/adjoint work, but **`transitionCrossflow` defaults to OFF since
  2026-07-24** (commit `a34441c9`): with it ON, the tutorial-wing mesh at
  M=0.15 stalls at ~3e-2 and never reaches deep convergence — reproduced from
  the original crossflow commit through HEAD; crossflow was only ever
  validated on AR5-type cases. Do not flip it back on without explicit user
  instruction (CLAUDE.md rule 3).
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
| Tu_∞ | `turbIntensityInf` exists in inputParam.F90 (~L702) |
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
| `"transitionUseApproxSA"` | bool | `True` | First-order (approximate) SA convection alongside the transition equations (see `SA_GAMMA_RETHETHA_BASE/adflow-vs-paper-solver.md`). Wired in `.pyf` (L1120). |
| `"transitionSrcDtRestrict"` | bool | `True` | Enable source-term dt restriction (P&Z Eq. 59). Caps λ_source × dt ≤ 0.9. |
| `"transitionSrcDtLimit"` | float | `0.9` | Threshold for source-term dt restriction (λ_source × dt ≤ this value). |
| `"srcDtDeactivateIters"` | int | `5` | Deactivate source-dt restriction after N consecutive clean (no-backtrack) turbKSP iterations **in the second-order regime** (`totalR ≤ ANKSecondOrdSwitchTol·totalR0`, the inexact-Newton analog of P&Z §IV.B.3). Counter resets when backtracking is triggered (even if it succeeds) or when the residual rises back above the switch tolerance. With the default `ANKSecondOrdSwitchTol = 1e-16` the regime is never entered ⇒ restriction never deactivates; set it to ~`1e-5` (paper's phase-switch value) to enable the acceleration. `0` = restriction inactive in turbKSP from the start (**not** "never deactivate"). DADI ignores this option (restriction always on there). Semantics fixed 2026-07-07 (D-A2-3). |
| `"TurbDADICoupled"` | str | `"full"` | DADI coupling mode: `"decoupled"` (3 scalar solves), `"transition"` (SA alone + γ-Re̅θt 2×2), `"full"` (3×3 block). |
| `"turbResScale"` | list/None | `None` (auto) | Residual scaling per equation. Auto-set to **`[10000.0, 0.1, 1.0e-4]`** for this model (`pyADflow.py:~6834` — ≈1/state-magnitude per equation, P&Z §IV.1 row scaling; matches the campaign-validated value in `convergence-strategy.md`). Override only to tune convergence balance. |
| `"transitionDampTheta"` | float | `0.99` | Back-off factor for per-variable γ/Re̅θt update damping in DD-ADI (P&Z Algorithm 2). |
| `"transitionDampMaxIter"` | int | `10000` | Safety cap on the back-off loop (unbounded in the paper); 10000 ⇒ effectively unbounded (0.99¹⁰⁰⁰⁰ ≈ 0). A hard clip to the bounds remains as **last-resort fallback only**: it can only fire after the loop exhausts, which requires the previous state to already be out of bounds; when it fires, a warning with cell counts prints advising to raise this option or investigate the upstream bound violation. Changed from 40 on 2026-07-07 (D-A2-5). |
| `"transitionCrossflow"` | bool | **`False`** | Helicity-based crossflow source D_scf (P&Z Eq. 15-26) on the Re̅θt equation. **Default flipped OFF 2026-07-24** (`a34441c9`): ON it stalls the tutorial-wing case (~3e-2 plateau); only ever validated on AR5-type cases. D_scf≡0 in 2D. ⚠️ The **Fortran default** in `inputParam.F90` is still `.true.` — Python always pushes `False`, but any Fortran-only path that never receives a Python `setOption` runs crossflow ON. |
| `"transitionRoughnessHeight"` | float | `3.3e-6` | Surface roughness height h for the crossflow correlation (Eq. 17), as a physical length in mesh units (metres). 3.3e-6 = 3.3 µm (smooth surface). |
| `"transitionRefLength"` | float | `-1.0` (auto) | Reference length l [mesh units] in the vorticity limiter (P&Z Eqs. 52-53; paper uses root chord — the physical cap scales as 1/√l, a calibration scale, NOT a unit conversion, so the "drop Re" rule of `nondimensionalization.md` does not apply). Negative = auto: uses the AeroProblem `chordRef` (via `inputPhysics%lengthRef`, refreshed at every `setAeroProblem`). Set explicitly to decouple from chordRef; `1.0` recovers the pre-option behavior (l = 1 m). Added 2026-07-07 to close finding D1. |
| `"transitionNK"` | bool | `True` | Master switch for the NK/ANK/turbKSP column-scaling + Eq. 59 bundle (incl. NK reactivation-on-backtrack, Algorithm 2 in NK) — 2026-07-16. Default preserves existing behavior; still additionally gated on `turbModel==SA-Gamma-Retheta` everywhere. |
| `"transitionNKAutoDisableTol"` | float | `0.0` | One-way latch, NK phase only: once the Newton residual norm drops below this fraction of `totalR0` (the **freestream** reference residual from `getFreeStreamResidual`, `solvers.F90:972` — NOT the restart-point residual; e.g. ~8.91e7 on `3D_Plain_Wing`), the `transitionNK` bundle (column scaling, Algorithm 2 damping, Eq. 59 reactivation) is turned off for the rest of the NK phase — i.e. "fall back to native NK." Default `0.0` never trips (unchanged behavior). **Tested 2026-07-18, `nk_switch_crossing_test`, and found unsafe at any point in NK**: tripping it either at NK engagement or ~10 outer iterations later (deep past engagement, residual already down 2+ orders) produces the identical catastrophic blowup (nuturb res → O(1e3), totalRes → O(1e9)) both times. Column scaling is load-bearing for the *entire* NK phase for this model — the 13-orders-of-magnitude state spread (ν̃, γ, Re̅θt) doesn't shrink with the residual, so "native NK" is never safe to fall back to. Kept as a diagnostic knob, not a recommended option. |
| `"transitionRowVolScale"` | bool | `False` | Eq. 58 (P&Z) geometric row-scaling factor on NK's residual rows (`volRef**(5/3)` flow, `volRef**(2/3)` turb, on top of existing `turbResScale`). **Off by default — genuinely tested 2026-07-16 and found to stall NK's linear solve** (lin res pinned ~1.0) on the 3D_Plain_Wing case; see `SA_GAMMA_RETHETHA_BASE/adflow-vs-paper-solver.md` §5. Not recommended until the volRef-vs-paper's-J correspondence is revisited. |
| `"transitionNKStallStepTol"` / `"transitionNKStallCountTrigger"` / `"transitionNKStallRtolCap"` | float / int / float | `0.1` / `3` / `1.0` | Attempted NK stall escape (`nk_switch_crossing_test`, 2026-07-18/19) — **NOT validated, do not rely on it.** Motivation: `getEWTol` (`NKSolvers.F90:2174`, standard PETSc EW) computes the Krylov `rtol` as `(norm/oldNorm)^1.618` — when the Newton step is pinned (`Step`~0), the ratio→1 and `rtol` rises to its 0.8 cap, i.e. EW picks the *loosest* linear solve exactly when stalled. This option forces `rtol` down to `transitionNKStallRtolCap` once `Step` has been below `transitionNKStallStepTol` for `transitionNKStallCountTrigger` consecutive NK iterations. `transitionNKStallRtolCap=1.0` (default) disables this (never caps) — pure diagnostic overhead only, no behavior change, safe to leave on any run. **One run with the cap at 0.05 pushed further than one earlier baseline run** (baseline froze at scaledTotalRes~1.7473e-3 for 300+ iterations and drifted up to ~1.79e-3; the capped run reached 1.4772e-3 by iter 308, later 9.4e-4 by iter ~6100 before re-stalling). **However, a same-day repeat of the unmodified baseline (only new, behaviorally-inert `print` diagnostics added elsewhere in this file) converged cleanly past the same point with no fix applied at all** — i.e. this exact case shows run-to-run variance in outcome from logically-identical code, most likely floating-point-order sensitivity in the matrix-free Jacobian (`-ffast-math`, `mpicc`/`mpifort` codegen) on a system that sits on a numerical knife-edge at this residual level. **Conclusion: the one observed improvement cannot be attributed to this option with confidence — it may equally be the same unexplained run-to-run variance.** Needs a controlled re-test (multiple repeats of both baseline and capped, same binary, before/after) to actually validate. |
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

> ⚠️ **Both options below are currently SILENT NO-OPS from Python** — exactly
> the `.pyf` bug documented at the end of this file (found 2026-08-12): they
> are mapped to the `ank` module in `pyADflow.py:6521-6522`, but the
> `module anksolver` block in `src/f2py/adflow.pyf` does not list
> `ank_physlstolretheta` or `omegamingamma`. The Fortran runs its hard-coded
> defaults (`omegaMinGamma = 0.05`, `NKSolvers.F90:2382`) regardless of any
> `setOption`. Fix = add both to the `.pyf` block; until then treat the
> defaults as the only reachable values.

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

### Matrix-dissipation eigenvalue limiters (2026-08-04)

Swanson & Turkel's Vn / Vl limiters, previously hard-coded Fortran
`parameter`s duplicated across five hand-written routines. Only active when
`discretization = "central plus matrix dissipation"` (the scalar/JST default
does not use them).

| Option | Type | Default | What it does |
|---|---|---|---|
| `"epsAcoustic"` | float | `0.25` | Vn — floor on the two acoustic eigenvalues, `λ = max(λ, Vn·rrad)`. Prevents zero dissipation at sonic lines. |
| `"epsShear"` | float | `0.025` | Vl — floor on the (triple) entropy/vorticity eigenvalue `\|u·n\|`. Prevents zero dissipation where the flow is parallel to the face. |

Why they matter here: inside a boundary layer the flow is parallel to the
wall, so `|u·n| ≈ 0` on wall-normal faces and `rrad ≈ a`. The shear-wave
dissipation there is therefore *entirely* set by Vl and is unrelated to the
local physics. P&Z 2020 §4.1 set **Vl = 0** for exactly this reason
("overly dissipative in the laminar boundary layer"), keeping Vn = 0.25
(0.30 for CRM-NLF). The ADflow default Vl = 0.025 is **not** the paper's
value; it is kept as the default here only for backward compatibility.
Lowering Vl costs robustness near stagnation points.

Wiring (all five hand-written sites read the same module variables):

| File | Routine |
|---|---|
| `src/modules/inputParam.F90:101-120` | `module inputDissipation` — declaration |
| `src/inputParam/inputParamRoutines.F90:~4306-4307` | defaults |
| `src/solver/fluxes.F90:419` | `inviscidDissFluxMatrix` |
| `src/solver/fluxes.F90:4357` | `inviscidDissFluxMatrixApprox` |
| `src/solver/fluxes.F90:5218` | `inviscidDissFluxMatrixCoarse` |
| `src/NKSolver/blockette.F90:2476, 4637` | ANK/NK residual path |
| `src/f2py/adflow.pyf:1024-1028` | `module inputdissipation` block |
| `adflow/pyADflow.py:~5935-5936, ~6206, ~6362-6363` | option defaults + `moduleMap` + option map |

They live in their own `inputDissipation` module rather than
`inputDiscretization` because the Tapenade-generated `fluxes_*.f90` do a
whole-module `use inputdiscretization`, which collided with their own local
`parameter` declarations of the same names.

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
    # turbResScale auto-set to [10000.0, 0.1, 1.0e-4]

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

## Known infra bug: `.pyf` option wiring (2026-07-16, recurring — not module-specific)

`src/f2py/adflow.pyf` is a **hand-maintained** f2py interface file, not
auto-regenerated from Fortran source on every build. Any Fortran module
variable (new option, new diagnostic array, anything) not explicitly listed
in its `.pyf` module block is invisible to real Python↔Fortran
communication — but f2py's `fortran`-type Python objects silently accept
**arbitrary attribute names** with no backing memory, so `setOption`/reading
the value back gives no error and no warning, and reads whatever phantom
Python attribute was last set instead of the (uninitialized/default)
Fortran memory. This has now bitten **three separate modules**, confirming
it's a structural hazard, not a one-off stale-file issue:

- `module inputiteration`: `transitionSrcDtRestrict`, `transitionSrcDtLimit`,
  `srcDtDeactivateIters`, `transitionDampTheta`, `transitionDampMaxIter`
  (**fixed — all five are now in `adflow.pyf` L1129-1133**; the earlier
  "still not fixed" status here was stale, corrected 2026-08-12).
  `transitionNK`/`transitionRowVolScale`/`transitionResidualAutoscale`
  (added 2026-07-16 — **fixed**). `transitionSrcDtEigMode` has since been
  removed from the codebase entirely (no longer a live bug).
- `module inputadjoint`: `storePsiHistory`, `psiHistoryStep`,
  `psiHistoryMax` (added 2026-07-24 for the psi-history adjoint-convergence
  diagnostic, see `VERIFICATION/adjoint-trace.md` — **fixed** same day,
  after ~3 hours of debugging a silent no-op that looked like a Fortran
  module-memory corruption bug before the real cause was found).
- `module adjointpetsc`: `psiHistory`, `psiHistoryIters`, `psiHistoryResid`,
  `psiHistoryCount` (same feature, same date — **fixed**).
- `module inputdissipation`: `epsAcoustic`, `epsShear` (added 2026-08-04 —
  **fixed at the time of writing**; the `.pyf` block was added in the same
  change, and `libadflow.inputdissipation.epsacoustic` was verified readable
  and writable from Python before the feature was declared done).
- `module anksolver`: `ank_physlstolretheta`, `omegamingamma` — **OPEN as of
  2026-08-12**: mapped in `pyADflow.py:6521-6522` but absent from the `.pyf`
  `module anksolver` block, so `ANKPhysicalLSTolReTheta` and `omegaMinGamma`
  are silent no-ops (Fortran hard-codes `omegaMinGamma = 0.05`,
  `NKSolvers.F90:2382`). See the warning in the Turb-ANK options table above.

**Before trusting any Python-side setting of ANY option or diagnostic
array that isn't already validated working — regardless of which Fortran
module it lives in — grep `adflow.pyf` for that module's block first** (`grep
-n "module <modulename>"` then check the variable is listed inside it). If
it's missing, the option/array is a no-op from Python regardless of what
`setOption`/read-back or even a direct `self.adflow.<module>.<var>` read
suggests — that read is reading the phantom Python attribute, not Fortran
memory. Symptom to recognize this by: a value reads back correctly in
Python right after `setOption`, but a Fortran-side `write(*,*)` of the
*same* variable inside a routine that runs later shows the type's zero
value (`.False.`/`0`/`0.0`) — that mismatch is the signature of this bug,
not memory corruption.
