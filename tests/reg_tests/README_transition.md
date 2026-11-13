# Transition-model regression tests

Tests of the three transition models added to ADflow. They follow the conventions of the other
regression tests (see `tests/README.md`): parameterized `RegTest` classes, one reference file per
case in `refs/`, training through the test handler, run with `testflo`.

| Key      | Model                  | Selected by                                  | States                                 |
|----------|------------------------|----------------------------------------------|----------------------------------------|
| `sagr`   | SA-noft2-Gamma-Retheta | `turbulenceModel = "SA-noft2-Gamma-Retheta"` | 5 mean flow + nuTilde, gamma, reThetat |
| `sgamma` | SA-noft2-Gamma         | `turbulenceModel = "SA-noft2-Gamma"`         | 5 mean flow + nuTilde, gamma           |
| `sabcm`  | SA-BCM                 | `useSABCM = True` with the SA model          | 5 mean flow + nuTilde                  |

## Files

| File                                | Checks                                                                   | Needs a solve |
|-------------------------------------|--------------------------------------------------------------------------|---------------|
| `test_jacVecProdFWD_<model>.py`     | residual and forward Jacobian-vector products vs ref, vs FD, vs CS       | no            |
| `test_jacVecProdBWDFast_<model>.py` | reverse vs reverse-fast products, forward/reverse dot products           | no            |
| `test_adjoint_<model>.py`           | total derivatives: adjoint vs ref, complex step vs adjoint               | yes           |
| `test_blockette_sagr.py`            | blockette residual equals block residual (`sagr` only)                   | no            |
| `reg_transition.py`                 | fixtures: input files, aero problems, options, block-wise assert helpers | -             |
| `dev/generate_<model>_restart.py`   | regenerates the converged restart state                                  | yes           |
| `run_transition_tests.sh`           | optional wrapper around the commands below                               | -             |

There is no `test_adjoint_sabcm.py` and no SA-BCM restart state: the SA-BCM tests evaluate the
residual and its derivatives at states built from the fully turbulent SA solution of the tutorial
wing (see `reg_transition.py`). `sgamma` has no blockette test because the blockette residual does
not implement that model.

### What the fixed-state tests check

All of them read a restart state, evaluate the residual once and compute Jacobian-vector products
about that state. Nothing is converged; each test takes a few seconds on two processors.

- `TestJacVecFwd<MODEL>` (real build): `dR/dw`, `dR/dXv` and `dR/d(alpha, mach, P, T)` products
  against the reference file. For `sagr` and `sgamma` the products are also split by variable
  block (`test_coupling_blocks`, `test_xDvDot_transition_rows`), which checks every block of the
  Jacobian on its own, including the coupling between the SA and the transition equations. For
  `sabcm` the residual norms are recorded as well.
- `TestJacVecFwd<MODEL>FD` (real build, `sagr` and `sgamma`): the same products against finite
  differences with a step-size sweep. The residual comparisons are marked `expectedFailure`
  because the transition residuals are too stiff for FD; complex step is the reference.
- `TestJacVecFwd<MODEL>CS` (complex build, `cmplx_test_*`): the same products by complex step,
  compared with the AD values stored in the reference file by the real build.
- `TestJacVecBWDFast<MODEL>` (real build): reverse-fast products equal the full reverse products,
  for a full seed and for one residual block at a time.
- `TestDotProducts<MODEL>` (real build, `sagr` and `sgamma`; for `sabcm` the dot product is part
  of `TestJacVecBWDFastSABCM`): forward and reverse modes are transposes of each other.
- `TestBlocketteResidualSAGR` (real build): the two residual implementations agree.

Two parameter sets run a model option about the same stored state (a Jacobian-vector product can be
evaluated about any state): `sgamma_artvisc_tut_wing` (`sgammaArtVisc = 10`) and
`sagr_crossflow_tut_wing` (`transitionCrossflow = True`; the term is active on this swept wing).

The adjoint tests have the same two options as extra parameter sets (`sagr_crossflow_tut_wing`,
`sgamma_artvisc_tut_wing`), each about its own converged state
(`dev/generate_sagr_restart.py --variant crossflow`, `dev/generate_sgamma_restart.py --variant artvisc`).
A set appears only when its state file exists in `input_files/`.

### What the adjoint tests check

- `TestAdjoint<MODEL>` (real build): residual norms, total derivatives of `cl`, `cd`, `cmz` and
  `drag` with respect to `alpha`, `mach`, `P`, `T`, twist, span and shape, and the adjoint
  states, against the reference file. One adjoint solve per function.
- `TestCmplxStep<MODEL>` (complex build, `cmplx_test_*`): re-converges the complex flow solution
  with a perturbation on one design variable at a time and compares `imag(f)/h` with the adjoint
  totals in the reference file. The complex build has no AD preconditioner, which limits the
  accuracy of this comparison to about 1e-7. The derivatives with respect to `mach` and the
  derivatives of the dimensional `drag` are printed but not asserted.

## Input files

`input_files/get-input-files.sh` downloads the standard ADflow input files, including the
tutorial-wing grid `mdo_tutorial_rans_scalar_jst.cgns` and FFD `mdo_tutorial_ffd.fmt` used here.

Two restart states are not in that archive and must be placed in `input_files/`:

| File                          | Used by       | Content                                                       |
|-------------------------------|---------------|---------------------------------------------------------------|
| `mdo_tutorial_sagr_dp.cgns`   | `*_sagr.py`   | tutorial-wing grid and converged SA-noft2-Gamma-Retheta state |
| `mdo_tutorial_sgamma_dp.cgns` | `*_sgamma.py` | converged SA-noft2-Gamma state on the tutorial-wing grid      |

Both are converged at Mach 0.15, alpha 1.8 deg, Tu 0.25 % and written in double precision. They
are distributed with the showcase package (`showcase/input_files/`, which also holds a copy of the
standard input files, so no download is needed: copy its content to `input_files/`), or can be regenerated with
`dev/generate_sagr_restart.py` and `dev/generate_sgamma_restart.py` (full flow solves, see
`dev/README.md`). After regenerating a state, retrain the reference files of that model.

## Running

From `tests/reg_tests/`. The tests set their own number of MPI ranks (2); `-n 1` runs one test at
a time.

Real build, fixed-state tests (no solve):

```bash
testflo -n 1 -v test_jacVecProdFWD_sagr.py test_jacVecProdBWDFast_sagr.py test_blockette_sagr.py
testflo -n 1 -v test_jacVecProdFWD_sgamma.py test_jacVecProdBWDFast_sgamma.py
testflo -n 1 -v test_jacVecProdFWD_sabcm.py test_jacVecProdBWDFast_sabcm.py
```

Complex build, fixed-state tests (no solve):

```bash
testflo -n 1 -v test_jacVecProdFWD_sagr.py test_jacVecProdFWD_sgamma.py test_jacVecProdFWD_sabcm.py -m "cmplx_test_*"
```

Adjoint tests (these converge an adjoint or a flow solution: run them on a cluster):

```bash
testflo -n 1 -v test_adjoint_sagr.py test_adjoint_sgamma.py                       # real build
testflo -n 1 -v test_adjoint_sagr.py test_adjoint_sgamma.py -m "cmplx_test_*"     # complex build
```

A single test:

```bash
testflo test_jacVecProdFWD_sagr.py:TestJacVecFwdSAGR_0_sagr_tut_wing.test_coupling_blocks
```

`testflo --dryrun <file>` lists the test names, `-s` shows the solver output.

## Training

The reference files are trained with the real build:

```bash
# fixed-state references (no solve)
testflo -n 1 -v test_jacVecProdFWD_sagr.py test_jacVecProdBWDFast_sagr.py -m "train*"
testflo -n 1 -v test_jacVecProdFWD_sgamma.py test_jacVecProdBWDFast_sgamma.py -m "train*"
testflo -n 1 -v test_jacVecProdFWD_sabcm.py -m "train*"

# adjoint references (adjoint solves: cluster)
testflo -n 1 -v test_adjoint_sagr.py test_adjoint_sgamma.py -m "train*"
```

| Reference file in `refs/`                                                  | Trained by                         |
|----------------------------------------------------------------------------|------------------------------------|
| `jacvecfwd_sagr_tut_wing.json`, `jacvecfwd_sagr_crossflow_tut_wing.json`   | `test_jacVecProdFWD_sagr.py`       |
| `jacvecbwd_sagr_tut_wing.json`, `jacvecbwd_sagr_crossflow_tut_wing.json`   | `test_jacVecProdBWDFast_sagr.py`   |
| `adjoint_sagr_tut_wing.json`                                               | `test_adjoint_sagr.py`             |
| `jacvecfwd_sgamma_tut_wing.json`, `jacvecfwd_sgamma_artvisc_tut_wing.json` | `test_jacVecProdFWD_sgamma.py`     |
| `jacvecbwd_sgamma_tut_wing.json`, `jacvecbwd_sgamma_artvisc_tut_wing.json` | `test_jacVecProdBWDFast_sgamma.py` |
| `adjoint_sgamma_tut_wing.json`                                             | `test_adjoint_sgamma.py`           |
| `funcs_rans_sabcm_{smooth,original}[_transitional]_tut_wing.json`          | `test_jacVecProdFWD_sabcm.py`      |

The complex-step tests are never trained: they compare against the values the real build stored.
Train the adjoint reference before running `TestCmplxStep<MODEL>`.
