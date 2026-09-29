# `dev/` — scripts outside the testflo workflow

The derivative suites live one level up (`../run_{sagr,sgamma,bcm}_tests.sh`,
see `docs/VERIFICATION/VERIF_00_three_stage_verification.md`). This folder
holds what testflo cannot do.

| Script                        | Purpose |
|-------------------------------|---------|
| `generate_sagr_restart.py`    | Converges the SA-GR tutorial-wing case from `reg_sagr.sagrBaseOptions` and writes `input_files/mdo_tutorial_sagr_dp.cgns` (there is no download server). Called by `run_sagr_tests.sh genw`. |
| `generate_sgamma_restart.py`  | Same for SA-noft2-Gamma (`reg_sgamma.py` → `mdo_tutorial_sgamma_dp.cgns`). Called by `run_sgamma_tests.sh genw`. |
| `diag_full_derivatives.py`    | Full per-DV adjoint-vs-CS table with solver output visible (see `../README_SAGR.md`). |
| `diag_full_derivatives_bcm.py`| SA-BCM analogue, per variant. |
| `run_bcm_case.py`             | Interactive SA-BCM run driver, to get a variant converging before `genw`. |

## Regenerating the restart state

```bash
# from tests/reg_tests/ -- a real flow solve: run it on the HPC
./run_sagr_tests.sh genw      # or run_sgamma_tests.sh / run_bcm_tests.sh
./run_sagr_tests.sh train     # then retrain the JSON refs
```

The generators read the case from the `reg_<model>.py` config, so the restart
always matches the tests. Each script puts `tests/reg_tests/` on `sys.path`, so
it runs from either directory.
