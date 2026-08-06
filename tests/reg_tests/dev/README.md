# `dev/` — SA-GR scripts outside the standard testflo workflow

The SA-noft2-Gamma-Retheta (SA-GR) derivative tests proper live one level up
in `tests/reg_tests/` and run the ADflow-standard way (testflo + JSON refs):

```
../test_jacVecProdFWD_sagr.py      Stage 3 forward: AD, FD (h-sweep), CS
../test_jacVecProdBWDFast_sagr.py  Stage 1 dot products, Stage 2 _b vs _fast_b
../reg_sagr.py                     shared case config + block helpers
../run_sagr_tests.sh               one entry point for the whole suite
```

This folder holds everything that does **not** fit that workflow but that we
still need to keep — either because it does something testflo cannot, or
because it reproduces an input we cannot download from a server.

## Keep-because-not-reproducible-otherwise

| Script | Why it lives here |
|--------|-------------------|
| `generate_sagr_restart.py` | Generates the converged restart **state (w)** the tests linearize about. Standard ADflow tests download a pre-converged CGNS from a server; we have none, so we regenerate it locally. Uses `reg_sagr.sagrBaseOptions` so it always matches the test case. |

### Regenerating the restart state (w)

```bash
# from tests/reg_tests/
./run_sagr_tests.sh genw            # wraps the call below
# or directly:
mpirun -np 2 <mach-python> dev/generate_sagr_restart.py
```

Then point `reg_sagr.sagrGridFile` / `sagrRestartFile` at the output and
retrain the JSON refs (`./run_sagr_tests.sh train`).

### Regenerating the JSON references

That IS standard testflo and does **not** live here:

```bash
./run_sagr_tests.sh train           # == testflo ... -m "train*"
```

## Diagnostics / one-off investigations (kept for reference)

These reproduce specific analyses; they are not part of the suite and are not
maintained as tests:

| Script | What it probes |
|--------|----------------|
| `check_3way_fwd.py`, `check_3way_fwd_sweep.py` | AD vs FD vs CS at a fixed state; the sweep tunes the FD step (`--turbmodel {sa,sagr}`, `--crossflow`). |
| `sweep_h_fd.py` | FD step-size sweep for a scalar functional derivative. |
| `diag_crosscall_cs.py`, `diag_all16_blocks_cs.py`, `diag_gamma_meanflow_cs.py` | The CS cross-call complex-buffer residue investigation (fixed in `reg_sagr.assert_coupling_blocks_allclose` by re-seating the real state). |
| `debug_cs_ar5_live.py` | Live CS block probing on the AR5 case. |
| `sanity_check_partials_sa.py`, `sanity_check_bwdfast_stage1.py` | Early standalone equivalents of the Stage-1/2 tests, before they were ported to testflo. Superseded by `../test_jacVecProdBWDFast_sagr.py`. |
| `verify_sa_derivatives_fd.py`, `test_adjoint_ar5_sa.py`, `test_adjoint_tutwing_mycode.py` | SA-only (not SA-GR) derivative/adjoint checks used as reference baselines. |

Each script inserts the parent `reg_tests/` dir on `sys.path`, so run them
directly with mpirun from `tests/reg_tests/` (or from here — the shim is
path-independent).
