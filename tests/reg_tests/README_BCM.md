# SA-BCM derivative tests

Test scaffolding for verifying SA-BCM's AD derivatives (Tapenade forward `_d`, reverse `_b`,
reverse-fast `_fast_b`) against finite differences and complex step, brought up to the same
validation standard as the sibling repo's SA-GR harness (`adflow_sa_gamma_rethetha_paper_solver`,
`tests/reg_tests/{reg_sagr.py,run_sagr_tests.sh,dev/}`) -- same tutorial-wing mesh/FFD
(mach=0.15, alpha=1.8), same dot-product / `_b`-vs-`_fast_b` / AD-vs-CS ladder structure, same
"blockette" primal-residual-operator check, same dev-script pattern for raw-output diagnostics.

Two SA-BCM variants exist and are exercised **everywhere** in this suite (test files, dev
scripts): `SABCM_Exp=False` ("smooth", tanh blend, manuscript default) and `SABCM_Exp=True`
("hard", exp-sqrt blend, Mura & Cakmakcioglu original).

| SA-BCM file | Mirrors | Verifies |
|---|---|---|
| `dev/run_bcm_case.py` | (new, SA-BCM has no equivalent before) | interactive raw-log run driver -- start here |
| `dev/diag_blockette_bcm.py` | (new) | block vs blockette residual, raw output -- run before trusting any `useBlockettes=True` run |
| `dev/diag_full_derivatives_bcm.py` | sibling repo's `dev/diag_full_derivatives.py` | full per-DV adjoint-vs-CS table, raw output |
| `generate_bcm_restart.py` | sibling repo's `dev/generate_sagr_restart.py` | converged restart CGNS (one per variant) |
| `test_blockette_bcm.py` | sibling repo's `test_blockette_sagr.py` | registered block==blockette assertion |
| `test_jacVecProdFWD_bcm.py` | `test_jacVecProdFWD.py` | `computeJacobianVectorProductFwd` + meanflow/nuTilde coupling blocks |
| `test_jacVecProdBWDFast_bcm.py` | `test_jacVecProdBWDFast.py` + `reg_test_utils.assert_dot_products_allclose` | `computeJacobianVectorProductBwdFast` + fwd/rev transpose consistency |
| `test_adjoint_bcm.py` | `test_adjoint.py` | total sensitivities from the adjoint solve (twist/span/shape DVGeo) |
| `reg_bcm.py` | `reg_default_options.py` + `reg_aeroproblems.py` + `reg_test_utils.py` + sibling repo's `reg_sagr.py` | fixtures + assert helpers |
| `run_bcm_tests.sh` | sibling repo's `run_sagr_tests.sh` | one entry point for the whole registered suite |

**Block-splitting, not "no coupling machinery".** SA-BCM adds no new transport equation (unlike
a gamma-Re_theta_t model, where `nw` grows and new off-diagonal blocks appear) -- `nw` stays 6,
same as plain SA. But `reg_bcm.py` still ports the sibling repo's `getStateBlocks`/
`maskStateVector`/`assert_coupling_blocks_allclose` machinery (that module's `nw==6` branch
exists specifically for this case): the blocks are just meanflow(5)/nuTilde(1) instead of a
4-way split, and `tTgamma`'s dependence on `rho, rlv, d2wall, chi(nuTilde)` makes the
meanflow<->nuTilde coupling nontrivial and worth locking down explicitly, same spirit as before.

## Key finding: `useBlockettes` is NOT force-disabled for SA-BCM

The sibling SA-GR repo force-disables `useBlockettes` (`pyADflow._updateTurbResScale`) because
its hand-synced `blockette.F90` residual copy was found to have drifted from the
Tapenade-differentiated primary path. **SA-BCM has the identical hand-synced-duplicate
structure in `blockette.F90:1004-1250`, but nothing force-disables it here** (the check only
fires for `turbulencemodel == "SA-noft2-Gamma-Retheta"`; SA-BCM's `turbulencemodel` stays
`"SA"`). `reg_bcm.py`'s `bcmBaseOptions` deliberately keeps `useBlockettes=True` (matching a
previously-used run script), so **today, by default, SA-BCM runs execute the unverified
blockette copy** -- with zero test ever having compared it to `sa.F90`.

Run `dev/diag_blockette_bcm.py` (raw output) and/or `test_blockette_bcm.py` (registered) FIRST,
before trusting any other result in this suite, for both variants.

## How to run (in order)

```bash
# 0. (once) fetch input_files/ if not already populated -- DONE 2026-07-24
cd /home/mdo/MDOLab_3_v2/adflow_sabcm/input_files
bash get-input-files.sh
# or, since the tutorial-wing grid+FFD are the same stock assets the sibling SA-GR repo uses:
#   cp /home/mdo/MDOLab_3_v2/adflow_sa_gamma_rethetha_paper_solver/input_files/{mdo_tutorial_rans_scalar_jst.cgns,mdo_tutorial_ffd.fmt} input_files/

cd tests/reg_tests

# 1. (real build, dev, RAW LOG) interactively stand up the case -- expect this NOT to
#    converge cleanly on the first try; watch stdout and adjust options in dev/run_bcm_case.py
#    between attempts.
mpirun -np 2 --bind-to core /home/mdo/packages_v2/mach/bin/python dev/run_bcm_case.py --variant smooth
mpirun -np 2 --bind-to core /home/mdo/packages_v2/mach/bin/python dev/run_bcm_case.py --variant hard

# 2. (real build, dev, RAW LOG) once a variant runs: is useBlockettes=True actually safe?
mpirun -np 2 --bind-to core /home/mdo/packages_v2/mach/bin/python dev/diag_blockette_bcm.py --variant smooth
mpirun -np 2 --bind-to core /home/mdo/packages_v2/mach/bin/python dev/diag_blockette_bcm.py --variant hard

# 3. (real build) generate the restarts the registered suite linearizes about
mpirun -np 2 --bind-to core /home/mdo/packages_v2/mach/bin/python generate_bcm_restart.py --variant smooth
mpirun -np 2 --bind-to core /home/mdo/packages_v2/mach/bin/python generate_bcm_restart.py --variant hard

# 4. (real build) registered blockette check -- should agree with step 2's raw numbers
./run_bcm_tests.sh blockette

# 5. (real build) TRAIN the reference files
#    !! only after any pending Tapenade rerun + make -- check docs/adjoint-trace.md !!
./run_bcm_tests.sh train

# 6. (real+complex) the full registered ladder: dot products, _b vs _fast_b, AD/FD, CS, adjoint
./run_bcm_tests.sh

# 7. (real+complex build, dev, RAW LOG) full per-DV adjoint-vs-CS table, per variant
mpirun -np 2 --bind-to core /home/mdo/packages_v2/mach/bin/python dev/diag_full_derivatives_bcm.py --variant smooth --mode adjoint --out /tmp/bcm_adj_smooth.json
mpirun -np 2 --bind-to core /home/mdo/packages_v2/mach/bin/python dev/diag_full_derivatives_bcm.py --variant smooth --mode cs --ref /tmp/bcm_adj_smooth.json
mpirun -np 2 --bind-to core /home/mdo/packages_v2/mach/bin/python dev/diag_full_derivatives_bcm.py --variant hard --mode adjoint --out /tmp/bcm_adj_hard.json
mpirun -np 2 --bind-to core /home/mdo/packages_v2/mach/bin/python dev/diag_full_derivatives_bcm.py --variant hard --mode cs --ref /tmp/bcm_adj_hard.json
```

Individual stages, once step 3-5 have been done once:

```bash
./run_bcm_tests.sh real        # Stage 1 (dot products), Stage 2 (_b vs _fast_b), Stage 3 AD/FD
./run_bcm_tests.sh cs          # Stage 3 CS ground truth (complex build)
./run_bcm_tests.sh adjoint     # full total-derivative adjoint vs CS
./run_bcm_tests.sh blockette   # blockette residual == block residual
```

## Verification ladder / how to read failures

0. **`diag_blockette_bcm.py` / `test_blockette_bcm.py` fails** → the hand-synced `blockette.F90`
   copy of `saSource` has drifted from `sa.F90`, same failure mode SA-GR hit before it was
   force-disabled. Every OTHER result in this suite is suspect until this is fixed (or
   `useBlockettes` is switched to `False`, mirroring SA-GR's default -- not done automatically
   here, see the "Key finding" section above).
1. **FD class fails, CS passes** → FD step/kink noise near the tanh blend (manuscript Eq. 8,
   smooth variant only) or the KS max() aggregation (Eq. 10, both variants), not an AD bug.
2. **CS key fails** (e.g. `dFuncs/dalpha`, `Eval Functions Sens: ..._cd shape`) → that specific
   term is wrong in the forward/adjoint code — cross-reference `docs/architecture.md`'s
   paper-symbol ↔ code-flag table to find which manuscript equation feeds it.
3. **`_fast_b` `nuTilde`-row-seed test fails but the full-state seed passes** → suspect
   `autoEditReverseFast.py` push/pop stripping on the `use_SABCM`-gated code before suspecting
   the model — see `docs/adjoint-trace.md`.
4. **Dot products fail** → fwd and rev disagree with each other (bookkeeping); they can also
   both be wrong identically and still pass — hence CS is required, not optional.
5. **`cmplx_test_aero_dvs` fails on `alpha`/`mach` specifically** → points at the vorticity-based
   `Re_theta`/`term1` or eddy-viscosity-ratio `term2` path (both freestream-condition-dependent,
   manuscript Eqs. 4–6).
6. **`cmplx_test_geom_dvs` (`shape`) fails, `cmplx_test_aero_dvs` passes** → points at the
   wall-distance/vorticity geometric-sensitivity path rather than the freestream-condition path.
7. **One variant (smooth/hard) fails and the other passes** → narrows the bug to the
   `SABCM_Exp`-gated blend branch specifically (Eq. 8 tanh vs. the original exp-sqrt) rather than
   the shared `Re_theta`/KS-max machinery both variants go through.

**Never fix a failing comparison by inflating its tolerance.** This harness exists because an
earlier internal manuscript (since replaced in `docs/papers/` by the two published AIAA papers,
see `docs/README.md`'s "Papers were replaced" note) once flagged that a *previous*
implementation's adjoint results needed revalidation. That manuscript's caveat text is no longer
in this repo, but the discipline it implies still applies: a line-by-line paper-vs-code check
(2026-07-24, `docs/adjoint-trace.md`) found the current implementation faithful, and this harness
is what actually proves the derivatives are correct — inflating tolerances to force a pass would
silently re-launder whatever the original bug was instead of catching it.

## NK/ANK convergence: what might help, deferred until the ladder above is green

`NKLSRelax` (off by default, `src/NKSolver/NKSolvers.F90:LSCubic`) was added 2026-07-24 for the
one convergence-tuning idea from an SA-GR NK-mods review that transfers without needing new
transported-state machinery (generic Armijo/turb-blowup relaxation). See `docs/TODO.md` for that
and the other (deferred, bigger-effort) idea — a nuTilde-diagonal source-dt-restriction
prototype. Neither is
implemented here -- tuning solver behavior against unverified derivatives risks chasing a bug
instead of a stiffness problem, which is exactly the ordering mistake the manuscript's own
caveat warns about.
