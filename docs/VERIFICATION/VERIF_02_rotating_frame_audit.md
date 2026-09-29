# Rotating-frame consistency of the transition models

The P&Z model is validated only on inertial-frame cases. On a rotating mesh the
transition terms are evaluated in the **relative (rotating) frame**, where the
blade boundary layer is steady and the correlations are defined. Every change
reduces exactly to the inertial formulas when Ω = 0.

## Frame facts (from the code)

- `w(:,:,:,ivx:ivz)` stores the **absolute** velocity: the rotating no-slip
  wall BC sets `uSlip = Ω×r` (`solverUtils.F90`).
- `vortx = curl(V_abs) − 2Ω = ω_rel`, the same convention as `sa.F90`
  (see the comment at the top of `turbUtils.F90`).
- `V_rel = V_abs − Ω×r`, with `Ω×r` computed at the cell centre (`sc`) from
  `omegax/y/z` (already scaled by `timeRef`), `sections(sectionID)%rotCenter`
  and the 8-node `eighth·Σ` pattern of `gridVelocitiesFineLevel_block`.

## Terms evaluated in the relative frame

| Term (P&Z Eq.)                              | Uses |
|---------------------------------------------|------|
| Vorticity cap `vortLim` (52–53)             | `uRefTrans = √(uInf² + |Ω×r|²)` (blade-element section speed) instead of `uInf` |
| θ_BL (4), t (7), δ (4), `û` for λ_θ (10–11) | `|V_rel|`, `û = V_rel/|V_rel|` |
| Helicity `H_cf` (24–26)                     | `ω_rel` and relative `û` |

Deliberately unchanged: velocity-gradient stencils and strain `S` (frame
invariant; in `dU/ds = û_i û_j ∂u_i/∂x_j` the antisymmetric rotation part
cancels), and the SA vorticity magnitude (ADflow's relative `−2Ω` default).

Applied in all copies of the transition logic in `saGammaRetheta.F90`: the
primal `Source`, the `#ifndef USE_TAPENADE` DADI/PC Jacobian recompute, and
`evalSrcJacBlock`. `saGamma.F90` carries the same `uRefTrans` for its
vorticity cap. There is no runtime guard: Ω = 0 ⇒ `sc = 0` ⇒ bit-identical.

`warnRotatingTransition` (rank 0, non-AD) reports once which form is active:
rotation defined in `sections%rotRate` → correction applied; rotation only on
`cgnsDoms%rotRate` (e.g. Python `setRotationRate`) → correction **not**
applied, same limitation as `sa.F90`'s own `−2Ω`.

**Caveats.** Near the axis in pure hover `uRefTrans → 0` (small hub region).
The rotating path (`rotRate ≠ 0`) has never been exercised numerically: the
Ω = 0 no-op is verified bit-exact and the derivative suites pass, but a
rotating case is still needed to validate the physics.
