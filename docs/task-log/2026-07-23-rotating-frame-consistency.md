# Rotating-frame consistency for the SA-sLM2015 transition model — 2026-07-23

**Problem:** The transition model was validated in the paper only on
inertial-frame cases and silently mixed reference frames on a rotating mesh /
rotating frame of reference. ADflow stores **absolute** velocity in `w` (a
rotating no-slip wall carries `u = Ω×r`; `solverUtils.F90:1491-1496`), so the
boundary-layer correlations — which are defined in the frame where the blade BL
is steady — must use the **relative** velocity `V_rel = V_abs − Ω×r`. As written,
`vortLim` used the freestream `uInf` (hover ⇒ `uInf→0` ⇒ γ sources die), and the
BL scales + helicity used the absolute velocity/vorticity.

**Verified frame facts (from code, not assumed):** `w = absolute velocity` (wall
BC above); `vortx = curl(V_abs) − 2Ω = ω_relative` (`turbUtils.F90:26-32`,
matches `sa.F90`).

**Fixes** (`src/turbulence/saGammaRetheta.F90`, all three copies: `Source`
residual, `Source` `#ifndef USE_TAPENADE` PC Jacobian, `evalSrcJacBlock`):
1. `vortLim` reference velocity `uInf → sqrt(uInf² + |Ω×r|²)` — the blade-element
   section speed. Exact for axial inflow (hover/climb/props/turbines), reduces to
   `uInf` when Ω=0. No new option.
2. BL scales (θ_BL, time scale, δ, `û` for λ_θ) → `V_rel`. Gradient stencils and
   strain unchanged (frame-invariant; the antisymmetric rotation part cancels in
   `dUds`).
3. Helicity `H_cf` → relative velocity · relative vorticity (`vortx`, keeping the
   `−2Ω`); the old `+2Ω` undo gave absolute-frame helicity, wrong on a rotor.

`Ω×r` (`sc`) computed at the cell center reusing `omegax/y/z = timeRef·rotRate`,
`sections(sectionID)%rotCenter`, and the 8-node cell-center pattern from
`solverUtils.F90`. No runtime guard: `rotRate` is a fixed input, so `sc` is
computed unconditionally and is 0 (bit-identical no-op) when `rotRate = 0`.

**Where rotation comes from / warning:** `sections%rotRate` (the rate the model
reads, same as `sa.F90`'s `−2Ω`) is set only at partitioning from
`cgnsDoms%rotRate`, populated by (a) the CGNS grid's `RotatingCoordinates`
(`readCGNSGrid.F90`) → propagates to `sections` → correction active; or (b) Python
`CFDSolver.setRotationRate()` (`pyADflow.py:1124` → `updateRotationRate`), which
updates `cgnsDoms` + grid velocities but **not** `sections` → correction not
active (same limitation as SA's own `−2Ω`). `warnRotatingTransition` prints a
once-per-run rank-0 message stating which form is active (globally-replicated
`cgnsDoms`/`sections`, no MPI); it lives in the non-AD driver
`saGammaReTheta_block` and is hidden from Tapenade with `#ifndef USE_TAPENADE`
(its `communication`/`cgnsGrid` deps otherwise broke AD codegen with a spurious
`myid` declaration).

**Verification:**
- `getResidual` on the converged non-rotating AR5 crossflow state
  (`input_files/ar5_plain_wing_sagr_crossflow_dp.cgns`): **BIT-IDENTICAL** to the
  original pre-edit code (max|Δres| = 0 over 1,397,872 states, global L2 =
  0.11614292259100134). This case exercises the helicity/crossflow path.
- Tapenade regenerated (`make -f Makefile_tapenade`), real + complex rebuilt.
- `./run_sagr_tests.sh real` → **13/13**; `./run_sagr_tests.sh cs` → **10/10**
  (forward AD == complex step).

**Not yet done (physics step):** the rotating branch itself is not exercised by a
test — AR5 is non-rotating, so `sc=0` masks the `Ω×r` formula. Numerical
validation on a real `rotRate≠0` case remains open. The abs-vs-relative vorticity
*magnitude* choice for SA production is left at ADflow's default (relative),
flagged as debatable by `turbUtils.F90:30-32`.

**Commits:** `e981a3e8` (frame fixes), `b255f71b` (guard + first Tapenade regen —
guard later removed), `4d28b5ce` (warning + drop guard + Tapenade regen).
See `docs/VERIFICATION/rotating-frame-audit.md` for the full audit.
