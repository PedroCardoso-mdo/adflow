# SA-BCM derivative tests

Verifies the SA-BCM AD derivatives with the same ladder, mesh and structure as
the SA-GR suite (see `README_SAGR.md` and
`docs/VERIFICATION/VERIF_00_three_stage_verification.md`). Two variants are
exercised everywhere: `SABCM_Exp=False` ("smooth", tanh blend) and
`SABCM_Exp=True` ("hard", exp-sqrt blend, Mura & Cakmakcioglu original).
Physics source: `docs/papers/AIAA20202714_SABCMPartI.md` (Appendix).

**Status: not runnable yet in this checkout.** The restarts
`input_files/mdo_tutorial_bcm_{smooth,hard}_dp.cgns` are absent and no
`refs/*bcm*.json` has been trained. Run `genw` (flow solve → HPC), then
`train`, before any other stage.

| File                              | Verifies |
|-----------------------------------|----------|
| `reg_bcm.py`                      | case config + block/assert helpers (nw = 6: meanflow(5) / nuTilde(1)) |
| `run_bcm_tests.sh`                | driver: all\|real\|cs\|adjoint\|blockette\|train\|genw |
| `generate_bcm_restart.py`         | converged restart CGNS, one per variant (`genw`) |
| `test_blockette_bcm.py`           | block vs blockette residual (run FIRST) |
| `test_jacVecProdFWD_bcm.py`       | forward mode + meanflow↔nuTilde coupling blocks |
| `test_jacVecProdBWDFast_bcm.py`   | `_fast_b` vs `_b`, fwd↔rev dot products |
| `test_adjoint_bcm.py`             | adjoint totals (twist/span/shape) vs complex step |
| `dev/run_bcm_case.py`             | interactive raw-log run driver, to get a variant converging before `genw` |
| `dev/diag_full_derivatives_bcm.py`| full per-DV adjoint-vs-CS table, raw output |

SA-BCM adds no transport equation, but `tTgamma` depends on `rho, rlv,
d2wall, chi(nuTilde)`, so the meanflow↔nuTilde coupling blocks are still
checked explicitly.

**`useBlockettes` is not forced off for SA-BCM** (the force-off in
`pyADflow._updateTurbResScale` covers only the GR and SA-sγ models, and
`reg_bcm.py` keeps `useBlockettes=True`). BCM runs therefore use the
hand-synced SA-BCM copy in `blockette.F90`, which only `test_blockette_bcm.py`
checks against `sa.F90` — run it first; if it fails, every other result is
suspect.

## Run

```bash
cd tests/reg_tests
./run_bcm_tests.sh genw        # both variants; expect to iterate with dev/run_bcm_case.py first
./run_bcm_tests.sh train       # after any Tapenade rerun + rebuild
./run_bcm_tests.sh             # blockette, real, cs, adjoint
```

Full per-DV table (per variant):

```bash
mpirun -np 2 --bind-to core $PY dev/diag_full_derivatives_bcm.py --variant smooth --mode adjoint --out bcm_adj_smooth.json
mpirun -np 2 --bind-to core $PY dev/diag_full_derivatives_bcm.py --variant smooth --mode cs --ref bcm_adj_smooth.json
```

## Reading failures (BCM-specific)

- `cmplx_test_aero_dvs` fails on `alpha`/`mach` → the vorticity-based
  `Re_θ`/`term1` or eddy-viscosity-ratio `term2` path (freestream-dependent).
- `cmplx_test_geom_dvs` (`shape`) fails, aero passes → wall-distance /
  vorticity geometric sensitivity path.
- One variant fails, the other passes → the `SABCM_Exp`-gated blend branch.
- FD noisy, CS passes → kinks at the tanh blend or the KS max aggregation.
