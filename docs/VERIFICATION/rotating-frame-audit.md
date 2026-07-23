# Rotating-frame consistency audit — SA-sLM2015 transition model

**Date:** 2026-07-23 · **Scope:** `src/turbulence/saGammaRetheta.F90` only
(SA model untouched). **Status:** primal edits landed and verified bit-identical
in normal (Ω=0) mode; **Tapenade rerun still required** for rotating-case adjoints
(see §5).

## 1. Motivation

The paper (Piotrowski & Zingg 2020) validates the model only on **inertial-frame**
cases (NLF0416, S809, NLF2-0415 swept wing with *periodic* BCs, Sickle wing). None
rotate. On a **rotating mesh / rotating frame of reference** the model silently
mixed reference frames, giving wrong (or dead) transition behavior on rotors. This
audit makes the whole transition model **frame-consistent in the relative
(rotating) frame** — the frame in which the blade boundary layer is steady and in
which the empirical correlations are defined.

## 2. Verified frame facts (from the code, not assumed)

- **`w(:,:,:,ivx..ivz)` stores ABSOLUTE velocity.** Proof: the rotating no-slip
  wall BC sets `uSlip = Ω×r` (`solverUtils.F90:1491-1496`) — the fluid at a
  rotating wall carries the blade's absolute speed, not zero.
- **`vortx = 2·fact·(wwy−vvz) − 2·omegax = curl(V_abs) − 2Ω = ω_relative`**
  (`saGammaRetheta.F90` transition block; matches `sa.F90:224` and the explicit
  comment `turbUtils.F90:26-32`). This was already correct and consistent with SA.
- Relative velocity: **`V_rel = V_abs − Ω×r`**, with `Ω×r` the local rotational
  grid velocity (`sc` below).

## 3. Findings and fixes (all exact no-ops when Ω=0)

| # | Term (paper Eq.) | Was | Now |
|---|---|---|---|
| 1 | Vorticity cap `vortLim` (52–53) | `uInf·√(uInf/(μ∞ l))/20` — freestream; hover `uInf→0` kills γ prod+dest | `uRefTrans·√(uRefTrans/(μ∞ l))/20`, `uRefTrans=√(uInf²+|Ω×r|²)` (blade-element section speed) |
| 2 | BL velocity scales: θ_BL (4), t (7), δ (4), `û` for λ_θ (10–11) | absolute `velMag=|V_abs|` | relative `velMag=|V_rel|`; `û = V_rel/|V_rel|` |
| 3 | Helicity `H_cf` (24–26) | absolute: undid `−2Ω` (`ω_abs`) · absolute `û` | relative: `ω_rel` (keep `−2Ω`) · relative `û` |

**Not changed (deliberately):**
- Velocity-*gradient* stencils (`uux..wwz`) and strain `S` (`Re_S`, `F_wake`):
  strain is frame-invariant, and in `dUds = û_i û_j ∂u_i/∂x_j` the antisymmetric
  rotation part of `∂(V_rel)/∂x − ∂(V_abs)/∂x` cancels against the symmetric
  `û_i û_j`, so only `û` needs the relative velocity.
- The relative-vs-absolute *vorticity magnitude* choice for SA production (`ss`,
  `vortMag`): ADflow itself flags this as debatable (`turbUtils.F90:30-32`); we keep
  ADflow's default (relative `−2Ω`) and stay consistent — no hard change.

### Why the blade-element section speed for the cap (finding 1)

The cap is `1/20` of the characteristic BL edge wall-vorticity
`Ω_wall ~ U_e/δ ~ √(U_e³/(ν l))`. On a blade the relevant edge velocity `U_e` in
the relative frame is the local section speed; combining freestream and tangential
speed in quadrature gives the standard blade-element resultant
`U_ref = √(uInf² + |Ω×r|²)`:
- Ω=0 → `uInf` (paper, exact no-op).
- hover (uInf=0) → `|Ω×r|` (local blade speed).
- forward flight → `√(uInf² + (Ωr)²)`.

Exact for axial inflow (`Ω×r ⊥ V∞`: hover/climb/propellers/wind turbines), a sound
characteristic magnitude for edgewise flow. **Caveat:** near the rotation axis in
pure hover `U_ref→0` (a small hub region); vastly better than today's *global*
collapse in hover, floor only if a case needs it.

## 4. Implementation

`Ω×r` is computed at the cell center (`sc`) reusing existing machinery — the
already-scaled `omegax/omegay/omegaz = timeRef·rotRate`, `sections(sectionID)%rotCenter`,
and the 8-node `eighth·Σ` cell-center pattern from
`solverUtils.F90:gridVelocitiesFineLevel_block`. Applied in **all three** copies of
the transition logic in `saGammaRetheta.F90`:
1. `Source` — primal residual (always compiled).
2. `Source` — `#ifndef USE_TAPENADE` PC/DADI Jacobian recompute (`qq(3,3)` crossflow).
3. `evalSrcJacBlock` — hand-written source-Jacobian block for `computeSrcLambda`.

For a non-rotating section `rotRate = 0 ⇒ sc = 0` exactly ⇒ `V_rel = V_abs`,
`uRefTrans = uInf`, helicity `+2Ω → 0` — every term is bit-identical.

## 5. Verification

- **Primal residual no-op (bit-exact):** `getResidual` on the converged
  non-rotating AR5 crossflow state (`input_files/ar5_plain_wing_sagr_crossflow_dp.cgns`,
  Mach 0.2, α=0, Tu=0.25%, crossflow ON) before vs. after the edits:
  **max|Δres| = 0.0 across all 1,397,872 states** on both procs (global L2 =
  0.11614292259100134 identical). This case exercises the helicity/crossflow path.
- **Derivative suite:** `./run_sagr_tests.sh real` — **13/13 pass** (Stages 1/2/3
  AD/FD). Since `sc` depends only on `x` and `Ω` (not `w`) and is 0 here, dR/dw,
  dR/dXv and dR/dMach are unchanged for Ω=0, so the (still-stale) generated
  Tapenade files remain consistent for this case.
- **Rotating-frame correctness is NOT yet exercised.** The AR5 case has Ω=0, so
  `sc=0` masks the `sc = Ω×r` formula; the no-op proves the substitutions reduce
  correctly but not that the rotating path is numerically right. A case with
  `rotRate≠0` is needed — this is the user's physics-validation step.

## 6. Adjoint / Tapenade status (frozen-adjoint rule 6)

The hand-written `Source`/`evalSrcJacBlock` changed; the Tapenade-generated files in
`src/adjoint/output{Forward,Reverse,ReverseFast}/` were **not** touched. They are
**stale for Ω≠0** (they still linearize the absolute-frame formulas). For
rotating-case adjoints/optimization the user must **rerun Tapenade** on
`saGammaRetheta.F90` and rebuild. For Ω=0 they are still consistent (verified §5),
so nothing regresses on the existing inertial test suite.

**TAPENADE NEEDED** — regenerate `saGammaRetheta%Source` (and `evalSrcJacBlock` if
in the AD scope) before trusting rotating-case gradients.
