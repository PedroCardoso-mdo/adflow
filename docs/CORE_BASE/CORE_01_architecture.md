# Code reference: transition models on this branch

The single reference for how the transition code is organised and what every
branch-specific option does. Physics equations are in
[`SAGR_01_paper_piotrowski_zingg_2020.md`](../SA_GAMMA_RETHETHA_BASE/SAGR_01_paper_piotrowski_zingg_2020.md)
(SA-GR) and [`docs/papers/`](../papers/) (SA-BCM). The non-dimensional
conventions are in
[`ADFLOW_08_nondimensionalization.md`](../ADFLOW_BASE/ADFLOW_08_nondimensionalization.md),
and the AD wiring is in
[`ADFLOW_09_adjoint_trace.md`](../ADFLOW_BASE/ADFLOW_09_adjoint_trace.md).
Where this file and the code disagree, the code is right.

## 1. Models

| Model  | Selected by                                   | Turbulent state        | Code                          |
|--------|-----------------------------------------------|------------------------|-------------------------------|
| SA-GR  | `turbulenceModel="SA-noft2-Gamma-Retheta"`    | ν̃, γ, Re̅θt (nwt = 3)   | `src/turbulence/saGammaRetheta.F90` |
| SA-sγ  | `turbulenceModel="SA-noft2-Gamma"`            | ν̃, γ (nwt = 2)         | `src/turbulence/saGamma.F90`  |
| SA-BCM | `useSABCM=True` on `turbulenceModel="SA"`     | ν̃ (nwt = 1)            | `useSABCM` block in `src/turbulence/sa.F90` |

- **Enum ids.** The enums are `spalartallmarasnoft2gammaretheta` and
  `spalartallmarasnoft2gamma` (id 9) in `src/modules/constants.F90`. SA-BCM
  has no enum of its own.
- **Transition as a modifier.** In SA-GR and SA-sγ, γ multiplies SA
  production only (P&Z Eq. 41); the SA equation itself is ADflow's.
- **Additive wiring.** Every code path that is not specific to Re̅θt accepts
  both enums (`turbModel == …gammaretheta .or. turbModel == …gamma`). The
  paths that exist only for Re̅θt stay GR-only: Algorithm 2 damping, `cs(3)`,
  the Re̅θt monitor/output/restart names, and `frozenTransition`.

### State vector

```
w(:,:,:,1:5)   ρ, ρu, ρv, ρw, ρE
w(:,:,:,itu1)  ν̃     (SA working variable)
w(:,:,:,itu2)  γ     (SA-GR, SA-sγ)
w(:,:,:,itu3)  Re̅θt  (SA-GR only)
```

- **Units.** All values are p-ρ non-dimensional: velocity scales with M·√γ,
  and `rlv`/`rev` are ratios to μ∞ (CLAUDE.md rule 11).
- **Freestream values.** `initializeFlow` sets γ∞ = 1 and Re̅θt∞ =
  `reThetaTCorrelation(Tu∞, 0)`.
- **Interior values.** The interior field starts at γ = 0.02, which keeps SA
  production suppressed until onset.
- **Restarts without transition fields.** `transitionRestartAlgebraicInit`
  instead maps a converged ν̃ field to γ through the SA-BCM term2
  (`initTransitionAlgebraicWarmStart`).

### Where things live

| What | Where |
|---|---|
| Residual, DADI solve, source Jacobian | `saGammaRetheta.F90`: `saGammaReTheta_block`, `Source`, `Viscous`, `saGammaReThetaSolve`, `evalSrcJacBlock`, `computeSrcLambda` |
| SA-sγ equivalents | `saGamma.F90`: `sgSource`, `sgViscous`, `sgResScale`, 2×2 DADI via `tdia2x2`, `computeSrcLambdaSaGamma` |
| Correlations and smooth min/max | `turbUtils.F90`: `reThetaTCorrelation`, `flengthCorrelation`, `rethetacCorrelation`, `smoothMinMax` |
| Model constants | `src/modules/paramTurb.F90` (`rsaGR*`, `rsaGRclampLambdaTheta`) |
| Options (Fortran) | `src/modules/inputParam.F90`; defaults in `inputParamRoutines.F90`; f2py exposure in `src/f2py/adflow.pyf` |
| Dispatch | `turbAPI.F90`; blockette residual `blockette.F90` (`blocketteRes`) |
| BCs | `turbBCRoutines.F90`: γ and Re̅θt zero-gradient at walls (`bmt = -1`); farfield uses the generic ghost = `wInf` branch |
| ANK/NK hooks | `NKSolvers.F90`: `computeTimeStepBlock`, `FormJacobianANKTurb`, `physicalityCheckANKTurb`, `applyNKAlgorithm2Damping`, `LSCubic`, `getTurbColScale` |
| Diagnostics | the `transitionDebug` block array (`block.F90`), written to the volume CGNS |
| AD output | `src/adjoint/output{Forward,Reverse,ReverseFast}/saGammaRetheta_{d,b,fast_b}.f90` and `saGamma_{d,b,fast_b}.f90` |

## 2. Solver paths

```
ANK (startup)
├─ coupled (ANKCoupledSwitchTol) ── flow + turbulence in one Krylov system
└─ decoupled: flow ANK, then turbulence by
   ├─ ANKUseTurbDADI = True ──┬─ TurbDADICoupled = "decoupled"  → 3 scalar solves
   │                          ├─ TurbDADICoupled = "transition" → SA scalar + γ-Re̅θt 2×2
   │                          └─ TurbDADICoupled = "full"       → 3×3 block (default)
   └─ ANKUseTurbDADI = False ── turb-ANK (turbKSP, GMRES)
NK (terminal, NKSwitchTol) ── fully coupled, matrix-free, LSCubic line search
```

Multigrid is not used with the transition models.

### Source-term time-step restriction (P&Z Eq. 59)

The restriction is λ_source·Δt ≤ `transitionSrcDtLimit`, where λ_source is the
largest positive eigenvalue of the source Jacobian. That Jacobian is
block-triangular (A13 = A31 = A32 = 0), so λ is exact: λ₃ = A33, and λ₁,₂ come
from the 2×2 block. Each path applies it differently:

- **DADI:** additive form, `qq(m,m) += srcLambda/limit`. DD-ADI has no Δt to
  take a max against. It never deactivates, because DADI *is* the
  globalization phase.
- **turbKSP, coupled ANK and NK** (when `transitionNK` is on): max form,
  `max(dtInv, srcLambda/limit)` on the turbulent diagonal.
- **Deactivation (turbKSP):** the restriction switches off after
  `srcDtDeactivateIters` clean iterations inside the second-order regime
  (`totalR ≤ ANKSecondOrdSwitchTol·totalR0`). It re-arms on any backtrack or
  when the residual rises. With the default `ANKSecondOrdSwitchTol = 1e-16`
  it never deactivates.
- **Converged solution:** both forms touch only the LHS diagonal, so the
  converged solution is the same either way.

### Physicality and damping

- **DADI (Algorithm 2):** γ and Re̅θt are damped *independently*. Each has its
  own θ^m back-off loop (`transitionDampTheta`), capped at
  `transitionDampMaxIter`. ν̃ only gets `max(w, 0)`. A hard clip to the bounds
  is a warned last resort that fires only when a loop exhausts.
- **turb-ANK:** γ uses absolute bounds [`rsaGRgammaLo`, `rsaGRgammaHi`]. The
  full step is taken unless it violates them, and `omegaMinGamma` floors the
  step where γ → 0. A relative check would collapse the step in laminar
  regions. ν̃ and Re̅θt use relative tolerances (`ANKPhysicalLSTolTurb`,
  `ANKPhysicalLSTolReTheta`).
- **Coupled path:** the same bounds act through `physicalityCheckANK*`, which
  limits the global λ. `ANKTransitionGlobalLambda` / `ANKAlgorithm2Damping`
  move γ/Re̅θt to per-node damping instead.

## 3. Options (SA-GR / SA-sγ)

The values below were checked against `_getDefaultOptions` in
`adflow/pyADflow.py`.

| Option | Default | Meaning |
|---|---|---|
| `transitionFirstOrderUpwind` | `True` | First-order upwind convection for γ and Re̅θt (paper §IV.A). Keep it on. |
| `transitionUseApproxSA` | `True` | First-order (approximate) SA convection alongside the transition equations. |
| `transitionCrossflow` | `False` | Helicity crossflow source D_scf on Re̅θt (P&Z Eqs. 15-26). It is ≡0 in 2-D. With it on, the tutorial wing plateaus at ~3e-2. The Fortran default is also `.false.`. |
| `transitionRoughnessHeight` | `3.3e-6` | Roughness h used by the crossflow correlation, in mesh units. |
| `transitionRefLength` | `-1.0` | Reference length l of the vorticity limiter (Eqs. 52-53). A value ≤ 0 means use the AeroProblem `chordRef`, which is refreshed on every `setAeroProblem`. This is a calibration scale, not a unit conversion. |
| `transitionSrcDtRestrict` | `True` | Enables the Eq. 59 restriction. |
| `transitionSrcDtLimit` | `0.9` | The limit in λ·Δt ≤ limit. |
| `srcDtDeactivateIters` | `5` | Clean second-order iterations before turbKSP deactivates the restriction. `0` means inactive from the start, not "never". DADI ignores it. |
| `TurbDADICoupled` | `"full"` | DADI coupling: `decoupled` / `transition` / `full`. |
| `transitionDampTheta` | `0.99` | Back-off factor of Algorithm 2 (DADI). |
| `transitionDampMaxIter` | `10000` | Back-off cap; effectively unbounded. The hard clip follows only if it exhausts. |
| `turbResScale` | `None` → auto | Row scaling of about 1/state magnitude. SA-GR `[1e4, 0.1, 1e-4]`, SA-sγ `[1e4, 0.1]`. A tuning knob, not physics. |
| `transitionNK` | `True` | Master switch for column scaling, Eq. 59 and Algorithm 2 in ANK/NK/turbKSP. Column scaling is needed for the whole NK phase, because the state spans about 13 orders of magnitude. |
| `transitionNKAutoDisableTol` | `0.0` | Latch that turns `transitionNK` off below this fraction of `totalR0` in NK. Unsafe at any depth; diagnostic only. |
| `transitionNKStallStepTol` / `…CountTrigger` / `…RtolCap` | `0.1` / `3` / `1.0` | Caps the Eisenstat-Walker `rtol` after repeated tiny NK steps. Inactive at `RtolCap = 1.0`; not validated. |
| `transitionRowVolScale` | `False` | Eq. 58 volume row scaling in NK. Stalls the NK linear solve; leave it off. |
| `transitionResidualAutoscale` | `False` | Eq. 58 S_a proxy. Marginal. |
| `transitionRestartAlgebraicInit` | `False` | Algebraic γ warm start for a restart that has no transition fields. |
| `ANKPhysicalLSTolReTheta` | `0.99` | Relative physicality tolerance for Re̅θt in turb-ANK. |
| `omegaMinGamma` | `0.05` | Step-factor floor for γ in turb-ANK. |
| `ANKNSubiterTurb` | `1` | Inner turbulence iterations. Applies to turbKSP only; DADI makes a single call. |
| `sgamma{VortLimiter, OnsetTanh, FturbLee}` | `False` | SA-sγ formulation switches. |
| `sgammaCoupleDestruction` | `True` | Couples γ into SA destruction (φ⁺(γ, 0.1)). |
| `sgammaFPGSmoothP` / `sgammaCTU1..3` | `300` / `163, 1002.25, 1.0` | φ± sharpness at the F_PG kink, and the Colonia C_TU constants. |

**Deprecated research flags (SA-GR):** `transitionBCMGamma`,
`transitionLocalReTheta` and `transitionReThetaInert` belong to the SA-G-s
study (`22_sa_g_model`). They remain only so those runs can be reproduced, so
build nothing new on them.

**Solver-fidelity options** (tags F1-F10 in the code comments, explained in
[`VERIF_06`](../VERIFICATION/VERIF_06_solver_fidelity_audit.md)):

- `ANKCFLMinCap`
- `ANKUnsteadyLSFactor`, `ANKUnsteadyLSMaxIter`, `ANKRejectOnLSExhausted`
- `ANKAlgorithm2Damping`
- `ANKTransitionGlobalLambda` (`True`)
- `MFFDFunctionError`, `MFFDType`
- `ANKColScaleUnit`
- `solverStallDiag`, `solverStallDiagStep`

**NK line search:** `NKLSRelax` (default `True`) runs `LSCubic` with Armijo
α = 1e-3 and a turbulence blow-up pre-limit factor of 3.0. With `False` it
reverts to upstream's 1e-2 and 2.0. The upstream α was almost never satisfied
on SA-GR. Looser values (α = 1e-4, factor 5.0) diverged, because NK has no
ρ/E physicality check.

**Matrix dissipation:** `epsAcoustic` (Vn, `0.25`) and `epsShear` (Vl,
`0.025`) live in `module inputDissipation`. They apply only with
`discretization="central plus matrix dissipation"`. P&Z use Vl = 0, because
the shear-wave dissipation inside a boundary layer is set entirely by Vl. They
have their own module because the Tapenade `fluxes_*` routines `use` all of
`inputDiscretization`.

**Monitor:** add `"scaledtotalr"` to `monitorVariables` to print the
S_r/S_a-scaled residual. It is for display only and never feeds `totalR` or
the switch tolerances.

### Adjoint options

- **Defaults** (all models): `adjointSolver="LGMRES"` (`adjointLGMRESAugDim=2`),
  `adjointSubspaceSize=400`, `ILUFill=3`, `ASMOverlap=3`. GMRES/LGMRES
  stagnate with a subspace of 200 or less on transition cases.
- **Field split, set automatically:** when `"gamma"` appears in
  `turbulenceModel` (SA-GR and SA-sγ) and the user did not set
  `globalPreconditioner`, the constructor sets it to `"field split"`. The
  split is PCFIELDSPLIT over {flow}{ν̃}{transition}; for plain SA the
  monolithic ASM is better.
  - `adjointFieldSplitType`: `multiplicative` / `additive`.
  - `adjointFieldSplitBlocks`: `3`, or 4 to split γ from Re̅θt.
- `frozenTransition` (`False`, SA-GR only): the adjoint treats γ and Re̅θt as
  constants. Their seeds are zeroed and their rows become identity in
  `master_state_b` / `master_b`. The primal is untouched, and no Tapenade
  regeneration is involved.
- `storePsiHistory` / `psiHistoryStep` / `psiHistoryMax`: write the ψ history
  of each adjoint solve to JSON (a convergence diagnostic).
- **Complex build:** it has no AD preconditioner, so CS re-converges need
  `ANKADPC`/`NKADPC=False` and the decoupled path.

### Input guards (`pyADflow.py`, transition models)

- `useWallFunctions` raises an error, because the transition model needs a
  resolved boundary layer.
- `useft2SA=True` only warns: the model is calibrated on SA-noft2.
- `turbIntensityInf ≤ 0` is rejected.
- `useBlockettes` is forced off. The transition kernels exist in
  `blocketteResCore` and `test_blockette_sagr.py` checks them, but every
  validated solve used the block path.

## 4. SA-BCM (`useSABCM`)

SA-BCM is an algebraic intermittency that multiplies SA production, inside
`saSource` in `sa.F90`. There is no new equation and `nw` is unchanged. The
code is the upstream PR version (mdolab/adflow#410, branch `sabcm-upstream`).

- **Implementation:**
  - The multiplier is the local `gammaBCM` (= 1 when off, which reproduces
    plain SA exactly). It is **not stored**: `saBCMIntermittency` (sa.F90,
    outside Tapenade) recomputes it from the state for the `intermittency`
    volume / isosurface variable.
  - `ft2` is forced to 0; `useft2SA=True` is ignored with a warning.
  - The hand Jacobian adds `dGammaBCM` to the `qq` diagonal, inside
    `#ifndef USE_TAPENADE`.
  - The residual is mirrored in `blockette.F90` (blockettes stay on for BCM).
- **Formulation (source: `docs/papers/`).**
  - The Appendix of AIAA 2020-2714 is the reference formulation.
  - Term1 compares the vorticity Re_θ = ρ|ω|d²/(2.193 μ) against
    Re_θc(Tu) = 803.73(100·Tu+0.6067)^-1.027, scaled by χ₁.
  - Term2 = (fv1·χ)/χ₂, which is the eddy-viscosity ratio, not raw χ.
  - `SABCMSmooth=True` (default): KS smooth max of Term1 + tanh blend.
    `SABCMSmooth=False`: the paper's γ = 1 − exp(−(√max(T1,0) + √T2)),
    no KS (the old `SABCM_Exp=True` KS-smoothed T1 first; not bit-equal).
- **Units:** vorticity is in p-ρ units and ν̃/ν are ratios to μ∞, so no extra
  1/Re appears (see ADFLOW_08).
- **Required:** RANS + SA, and `useApproxWallDistance=True`;
  `inputParamRoutines` terminates otherwise.
- **Hook:** `sa.F90` is guarded by a hook, so each edit asks for approval.
  Approve only for SA-BCM work.

| Option | Default | Meaning | Old name (before 2026-09-30) |
|---|---|---|---|
| `useSABCM` | `False` | Master switch. | `use_SABCM` |
| `SABCMSmooth` | `True` | Smooth (KS + tanh) vs original exp form. | `SABCM_Exp` (inverted) |
| `SABCMChi1` / `SABCMChi2` | `0.002` / `0.02` | χ₁ and χ₂. | `SABCM_Const1` / `SABCM_Const2` |
| `turbIntensityInf` | `0.001` | Tu∞ as a **fraction** (0.005 = 0.5 %). Shared with SA-GR. | `SABCM_TU` (percent) |
| `SABCMTanhCenter` / `SABCMTanhWidth` | `0.5` / `0.08` | Centre and width of the tanh blend. | `SABCM_S0_tanh` / `SABCM_fsmooth` |
| `SABCMRho` | `50.0` | KS sharpness replacing max(Term1, 0) (smooth form only). | `SABCM_maxsmooth` |

`transitionBCMGamma` (SA-GR modifier) reuses `SABCMChi1/Chi2/Rho/TanhCenter/
TanhWidth`. Volume variable `tgamma` is gone: use `intermittency`.
The SA-BCM derivative tests are upstream's `tests/reg_tests/test_sabcm.py`.

## 5. Design rationale (why the code looks like this)

- **Safeguards not written in P&Z.** The λ_θ clamp to [−0.1, 0.1] and the Tu
  floor of 0.027 % (`rsaGRtuFloor`) come from Langtry–Menter 2009. Without the
  clamp, exp(−35 λ_θ) overflows near stagnation. The clamp is a compile-time
  switch (`rsaGRclampLambdaTheta`); with it off, S809 stalls. The Jacobian
  treats the clamped target as constant.
- **Diagonal clips `qq(1,1)`, `qq(2,2)`, `qq(3,3) ≥ 0`.** They mirror SA's
  DDADI clip; `qq(2,2)` goes negative routinely before transition. They are
  LHS-only, inside `#ifndef USE_TAPENADE`.
- **ν̃ row scale 1e4, not the paper's 1e3.** ν̃'s residual is exactly
  ADflow's SA residual, so it keeps SA's scaling; γ and Re̅θt use about
  1/state magnitude.
- **Primal code is written without in-place `x = smoothMinMax(x, …)`.**
  `autoEditReverseFast.py` strips the push/pop that such a line needs, which
  breaks `_fast_b`. Use distinct targets, as in
  `lambdaThetaRaw → lambdaThetaClamped → lambdaThetaLocal`.
- **Option plumbing.** `src/f2py/adflow.pyf` is maintained by hand. A module
  variable missing from its block makes the Python `setOption` a silent no-op:
  the value reads back fine in Python, but Fortran sees the default. Grep the
  module block in `adflow.pyf` before trusting a new option.

## 6. Known limitations

- **Deep-NK wall.** Below a relative residual of about 5e-9, NK's linear
  residual sits at 0.8-0.99 with GMRES exhausted. The limit is the
  preconditioner; no option moves it.
- **No ρ/E physicality check in NK** (ANK has one). This is why the
  line-search settings in `LSCubic` cannot be loosened further.
- **Blockettes are forced off** for SA-GR and SA-sγ (§3). This costs
  performance, not correctness.
- **Crossflow.** With `transitionCrossflow=True` the tutorial-wing case does
  not converge deeply.
- **Rotating frames.** The rotating-frame path (`rotRate ≠ 0`) has never been
  exercised numerically (see
  [`VERIF_02`](../VERIFICATION/VERIF_02_rotating_frame_audit.md)).
- **Initial γ.** The interior initial γ = 0.02 has never been tuned.
