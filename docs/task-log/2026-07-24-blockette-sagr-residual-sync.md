# Blockette SA-GR residual synced to the block path (advection sign bug) — 2026-07-24

**Problem:** ADflow builds the RANS residual through two interchangeable paths
inside `blocketteRes` (`src/NKSolver/blockette.F90`), selected by `useBlockettes`:

- `useBlockettes = False` → `blockResCore` → `saGammaReTheta_block`
  (reference residual, `src/turbulence/saGammaRetheta.F90`)
- `useBlockettes = True`  → `blocketteResCore` → inlined
  `saGammaRethetaSource / Advection / Viscous / ResScale` (cache-blocked copies)

The inlined SA-GR kernels had **drifted** from the model (audit-06 F7). pyADflow
force-disables blockettes for SA-GR (`_updateTurbResScale`, `pyADflow.py:6668`),
so the drift was latent, but it meant the blockette path could never be turned on.

**Root cause (the real bug):** the inlined `saGammaRethetaAdvection` used the
first-order upwind difference with the **wrong sign** —
`dw −= uu·(w(upstream) − w(cell))` instead of the validated
`dw −= uu·(w(cell) − w(upstream))` used by `saAdvection` / `turbAdvection` —
plus a different (also sign-flipped) minmod form. Since transition convection is
first-order-upwind-enforced (CLAUDE.md rule 5), this flipped term dominated,
making the blockette turbulence residual disagree with the block path by a
factor ~2 (near sign-flip) on all three turb variables. Meanflow was unaffected.

**Fixes** (`src/NKSolver/blockette.F90`, residual copies only — not Tapenade-
differentiated):
1. `saGammaRethetaAdvection` rewritten as a faithful `nAdv=3` copy of the
   validated `saAdvection` (correct upwind sign + folded minmod), keeping the
   `transitionFirstOrderUpwind` enforcement. **This is what fixed the factor-2.**
2. `saGammaRethetaSource` re-synced to `saGammaRetheta.F90` `Source`: rotating-
   frame relative velocity `V_rel = V_abs − Ω×r` (cell-center `sc`) feeding
   velMag/uₕₐₜ; blade-element `uRefTrans = √(uInf²+|Ω×r|²)` in the vorticity
   limiter; `gammaForSA` clamp margin `one → one+xminn`; distinct lambdaTheta
   clamp targets; helicity on relative velocity·relative vorticity. All reduce
   to bit-identical when `rotRate=0` and `transitionCrossflow` inert.
3. `saGammaRethetaViscous` / `ResScale` were already consistent — no change.

**Verification** — new same-`w` → same-residual test (block vs blockette), now a
permanent part of the SA-GR suite: `tests/reg_tests/test_blockette_sagr.py`
(testflo, `TestBlocketteResidualSAGR`, N_PROCS=2), wired into
`run_sagr_tests.sh` as the `blockette` stage and into `all`/SUMMARY. A
standalone copy also lives at
`.../3D_Plain_Wing/block_vs_blockette_residual/run_block_vs_blockette.py`.
Restarts the converged `mdo_tutorial_sagr_dp.cgns` state, perturbs ~1% (so all
source/advection/diffusion terms are active), evaluates `getResidual()` with
`useBlockettes` False then True about the identical state, asserts per variable
(`rtol=1e-7`). The suite version sets `transitionCrossflow=True` to run the
crossflow `D_scf` Source branch in both paths.
- **Before fix:** meanflow matched to ~1e-10, but nuTilde/gamma/reThetat rel diff
  ≈ 2.0 → FAIL (this caught the bug).
- **After fix:** all 8 variables agree to rel ≤ 1.2e-10 (summation-order
  roundoff) → **PASS**. Run: `OMP_NUM_THREADS=1 mpirun -np 4 --bind-to core`.

**Coverage limitation:** the test case is the non-rotating tutorial wing, so
the rotating-frame Source branches re-synced here are inert (not numerically
exercised); a rotating/swept restart would exercise them. *(Corrected
2026-08-12: this section originally also claimed `transitionCrossflow=False`,
contradicting the paragraph above — verified today that
`test_blockette_sagr.py:73` sets `transitioncrossflow = True`, overriding the
suite's crossflow-OFF base config, so the crossflow `D_scf` branch IS
exercised by this test.)*

**Not done (deliberate):** pyADflow still force-disables blockettes for SA-GR.
Now that the residual operator matches, flipping that on is a separate outward-
facing change needing full NK/ANK *solve* validation — left to the user.
