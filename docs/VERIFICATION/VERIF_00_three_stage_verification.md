# Derivative verification ladder

The adjoint of every transition model on this branch is checked by the same
ordered ladder, registered as testflo suites in `tests/reg_tests/` and driven
by one shell script per model. All cases linearize about a converged state on
the standard ADflow **tutorial wing** (mach 0.15, α 1.8). The case — mesh,
restart state `w`, AeroProblem, options — lives only in the `reg_<model>.py`
config module.

**How to run it:** use the user skill `adflow-derivative-verification`
(`~/.claude/skills/adflow-derivative-verification/`). Its `preflight.sh
<sagr|sgamma|bcm>` checks libraries, restart and refs deterministically and
prints the fix for anything missing; then run the stages in order.

## Suites

| Model                            | Driver                | Config         | Tests                                                               | Restart `w` (`input_files/`)             |
|----------------------------------|-----------------------|----------------|---------------------------------------------------------------------|------------------------------------------|
| SA-GR (`SA-noft2-Gamma-Retheta`) | `run_sagr_tests.sh`   | `reg_sagr.py`  | `test_{jacVecProdFWD,jacVecProdBWDFast,adjoint,blockette}_sagr.py`  | `mdo_tutorial_sagr_dp.cgns`              |
| SA-sγ (`SA-noft2-Gamma`)         | `run_sgamma_tests.sh` | `reg_sgamma.py`| `test_{jacVecProdFWD,jacVecProdBWDFast,adjoint}_sgamma.py`          | `mdo_tutorial_sgamma_dp.cgns`            |
| SA-BCM (`use_SABCM`, smooth+hard)| `run_bcm_tests.sh`    | `reg_bcm.py`   | `test_{jacVecProdFWD,jacVecProdBWDFast,adjoint,blockette}_bcm.py`   | `mdo_tutorial_bcm_{smooth,hard}_dp.cgns` |

Driver arguments (all three): `all` (default) · `real` · `cs` · `adjoint` ·
`blockette` (not SA-sγ) · `train` · `genw`. Knobs: `PY`, `NP` (default 2),
`NP_CS` (default 1).

- **`genw`** regenerates the restart state. The CGNS is gitignored and has no
  download server, so a fresh checkout **must** run `genw` before anything
  else. It is a flow solve → run it on the HPC.
- **`train`** rewrites `refs/*.json`. Needed after any change to
  differentiated code (Tapenade regen + rebuild first), or when refs are absent.

## The stages, cheapest first

| # | Stage                            | Driver arg   | Build   | What it proves |
|---|----------------------------------|--------------|---------|----------------|
| 0 | Blockette ≡ block residual       | `blockette`  | real    | The cache-blocked residual kernels equal the reference block residual for the same `w` |
| 1 | Dot product fwd ↔ rev            | `real`       | real    | `<ẇ, Jᵀψ> = <ψ, Jẇ>`: `_d` and `_b` agree with each other |
| 2 | `_b` vs `_fast_b`                | `real`       | real    | `autoEditReverseFast.py` did not strip a needed push/pop (row-block seeded) |
| 3a| AD vs trained refs, FD h-sweep   | `real`       | real    | Regression vs refs; FD is context only |
| 3b| **AD vs complex step (partials)**| `cs`         | complex | The linearization is the derivative of the primal: `dR/dw`, `dR/dXv`, `dR/dDV`, all coupling blocks |
| 4 | Total `df/dx`, adjoint           | `adjoint`    | real    | Assembled adjoint totals vs refs (α, mach, twist, span, shape) |
| 5 | Total `df/dx`, complex step      | `adjoint`    | complex | Adjoint totals = CS totals from a re-converged complex solve per DV |

Stages 0–3b evaluate matvecs about the stored `w`; 4–5 converge real and
complex solves and are the expensive part.

Why this order: stage 1 catches fwd/rev bookkeeping cheaply but not a bug
present identically in both modes; stage 2 isolates the fast-reverse
post-processing; only when both agree does stage 3b (h = 1e-40, no
cancellation) decide whether the shared linearization is right. Stages 4–5
can fail for solver reasons (convergence floor) even when the partials are
exact — validated partials ≠ validated gradient.

## Invariants — do not "fix"

- **FD residual tests are `@expectedFailure`** on the transition models. The
  residual spans ~13 orders; element-wise FD cannot meet tolerance at any step
  (cancellation on near-zero cells). `assert_fd_allclose_hsweep` still sweeps
  `h` and reports the best step. CS is the enforced truth. An "unexpected
  success" means drop the decorator — never loosen tolerances to force a pass.
- **CS coupling blocks re-seat the state** (`setStates(real(getStates()))`)
  before each column. Sequential CS calls share a complex work buffer; without
  the re-seat a ~2e-8 imaginary residue appears in the structurally-zero
  mean-flow rows of the next column.
- **The complex build has no AD preconditioner.** The CS classes set
  `ankadpc/nkadpc=False` and `ankcoupledswitchtol=1e-16`; the FD-PC complex
  re-converge stalls near 1e-8. Hence the SA-GR/SA-sγ CS adjoint check uses
  `rtol = atol = 5e-8`, with `mach` and `drag` reported but not asserted (mach
  CS never settles; `drag` duplicates `cd`).
- **Run CS at the same `-np` used for `train`.** A different partition
  linearizes about a slightly different discrete state (~1e-10).

## Rules for differentiated code

- **Never write `x = smoothMinMax(x, …)` in place.** `autoEditReverseFast.py`
  strips the push/pop of the overwritten value, so `_fast_b` goes wrong (it
  once did, with the wrong sign, on the Re̅θt rows) while `_b` and CS agree.
  Use distinct targets (`raw` → `clamped` → `local`) and regenerate Tapenade.
- **Keep residual-path code kind-generic.** The residual files contain no
  `USE_COMPLEX` branches; the complex build is a complexified copy of the same
  source.
- **No bare real literals into `smoothMinMax`** — pass `paramTurb` parameters
  so the complexified kinds match.
- **Real-part casts only inside `#ifdef USE_COMPLEX` step-limiters** in the
  solver (`NKSolvers.F90`), never in residual math.
- An active input missing from a Tapenade head (e.g. `uInf`, `muInf` for the
  vorticity limiter) shows up as a stage-3b failure on the `xDvDot` transition
  rows. Tapenade touchpoints: [ADFLOW_09](../ADFLOW_BASE/ADFLOW_09_adjoint_trace.md).

## Reading a failure

| Symptom                                         | Look at |
|-------------------------------------------------|---------|
| Stage 0 fails                                   | Hand-synced `blockette.F90` kernel drifted from the model file; all other results are suspect |
| Stage 1 fails                                   | Forward/reverse wiring: Tapenade heads, `masterRoutines.F90` dispatch |
| Stage 2 fails on transition-row seeds only      | `autoEditReverseFast.py` stripping; in-place `smoothMinMax` |
| Stage 3b fails, 1–2 pass                        | Differentiated code out of sync with the primal → `TAPENADE NEEDED`, or a missing active input |
| FD fails, CS passes                             | FD kink/cancellation noise, not an AD bug |
| Stages 4–5 fail, 3b passes                      | Primal or adjoint convergence depth, complex re-converge floor |
| SA-BCM: one variant (smooth/hard) fails         | The `SABCM_Exp`-gated blend branch |

## Complex build

`libadflow_cs.so` must be built against `PETSC_ARCH=complex-debug`; building
against `real-debug`, or reusing `.mod` files from a real-arch attempt, gives
spurious `COMPLEX(8)→REAL(8)` interface errors in `fortranPC.F90`,
`adjointAPI.F90` and friends. Use `/build complex`; if the tree is dirty:

```bash
git clean -fdx src_cs/            # -x is required: src_cs is gitignored
PETSC_ARCH=complex-debug make -f Makefile_CS
```

The tests import adflow from site-packages, not `./adflow`: reinstall
(`rm -rf build && pip install . --no-deps`) after every rebuild.

## Not covered by the standing suites

- Crossflow (`D_scf`): the suites run with `transitionCrossflow=False` (the
  model default). Set `transitioncrossflow: True` in `reg_sagr.py`, `genw`,
  `train`, run to exercise it.
- Rotating frames (`rotRate ≠ 0`): see [VERIF_02](VERIF_02_rotating_frame_audit.md).
