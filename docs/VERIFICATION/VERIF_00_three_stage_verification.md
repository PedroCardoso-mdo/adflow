# Three-Stage Low-Level Adjoint Verification Ladder

**Status (plain SA): all three stages PASS**, for plain SA on both the
tutorial-wing and AR5 meshes, as of 2026-07-21.

**Status (SA-GR, `nw=8`): all three stages PASS**, with crossflow
**OFF** in the standing suite (the model default flipped to
`transitionCrossflow=False` on 2026-07-24), as of 2026-07-22. Stage 2 initially
FAILED on the reThetat rows (`_b` vs `_fast_b` diverged, max abs diff ≈
7.5e+2); this was root-caused to an in-place `smoothMinMax` double-clamp of
`lambdaThetaLocal` in the primal and **fixed** (see "Results & provenance"
below and `../task-log/2026-07-22-fastb-retheta-lambdatheta-fix.md`). After
the fix + Tapenade regen, Stage 2's reThetat block matches to 1.4e-9.

**As of 2026-07-22 the ladder is a registered testflo suite**, not the
ad-hoc `dev/` scripts. See "Canonical way to run" immediately below;
condensed results and provenance follow it (the original per-script commands
were trimmed — those scripts now live in `tests/reg_tests/dev/`).

---

## Canonical way to run (start here)

Everything runs from `tests/reg_tests/` through one script,
`run_sagr_tests.sh`. The whole case (mesh, restart state `w`, AeroProblem,
solver options, **crossflow OFF** — `transitioncrossflow: False`, the model
default since 2026-07-24) is centralized in `reg_sagr.py`. The standing case
linearizes about the **tutorial wing** (`mdo_tutorial_sagr_dp.cgns`,
mach 0.15, alpha 1.8).

> Note (2026-08-12): the crossflow adjoint WAS validated once, on
> 2026-07-22, with crossflow ON on the AR5 state (see
> `../task-log/2026-07-22-fastb-retheta-lambdatheta-fix.md` and the pass
> table below) — but it is NOT exercised by the standing suite anymore.

```bash
cd tests/reg_tests/
./run_sagr_tests.sh          # whole ladder: real (Stage 1/2/3-AD-FD) + complex (Stage 3 CS) + adjoint + blockette
./run_sagr_tests.sh real     # Stage 1 (dot products) + Stage 2 (_b vs _fast_b) + Stage 3 AD/FD
./run_sagr_tests.sh cs       # Stage 3 complex-step ground truth (the decisive check)
./run_sagr_tests.sh adjoint  # full total-sensitivity test (test_adjoint_sagr.py)
./run_sagr_tests.sh blockette # blockette SA-GR residual-sync guard (test_blockette_sagr.py)
./run_sagr_tests.sh train    # rebuild the JSON reference files (after changing derivative code)
./run_sagr_tests.sh genw     # regenerate the converged restart state w (rarely needed)
```

**You do NOT need `genw` or `train` to just run the tests** — they read the
existing restart CGNS and JSON refs. `genw` is only for a missing/changed
`w` (the CGNS is gitignored, so a fresh clone lacks it); `train` is only
after the derivative code changes.

**The registered tests (standard testflo, ADflow-conventional):**

| File | Stage(s) | Classes |
|---|---|---|
| `test_jacVecProdBWDFast_sagr.py` | 1 (dot products), 2 (`_b` vs `_fast_b`) | `TestDotProductsSAGR`, `TestJacVecBWDFastSAGR` |
| `test_jacVecProdFWD_sagr.py` | 3 (AD vs ref, FD, CS) | `TestJacVecFwdSAGR` (AD), `TestJacVecFwdSAGRFD` (FD), `TestJacVecFwdSAGRCS` (CS) |
| `reg_sagr.py` | shared config + block helpers | — |
| `refs/jacvec{fwd,bwd}_sagr_tut_wing.json` | trained AD/dot-product references | — |

**Real vs complex build.** Stages 1, 2, and the AD/FD half of Stage 3 use
the real build (`test_*` methods). The decisive CS half of Stage 3 uses the
complex build (`ADFLOW_C`, `cmplx_test_*` methods); the runner selects those
with `-m "cmplx_test_*"`. Both libraries are installed in the same `mach`
env on this branch, so no separate `mach_cs` python is needed. If the
Fortran changed, rebuild + `pip install` first (and rerun Tapenade + rebuild
the complex lib if the change touches differentiated code).

**Two things future-Claude must not "fix":**
- **FD residual tests are `@expectedFailure`.** The SA-GR residual spans
  ~13 orders; the element-wise `assert_allclose` metric cannot be met by FD
  at any step (subtractive cancellation on the near-zero cells), so
  `assert_fd_allclose_hsweep` sweeps `h`, reports the best step, and the
  three FD residual methods are marked expected-fail. This is **by design** —
  CS (no cancellation) is the enforced ground truth. Do **not** loosen the
  tolerance to force a pass (post-mortem SST lesson). An "unexpected success"
  means a better-converged mesh now lets FD pass — then drop the decorator.
- **CS coupling-blocks re-seats the state each column.** Sequential CS
  Jacobian-vector calls share a complex work buffer; a prior large-norm
  column leaves a ~1e-8 imaginary residue that surfaces in the
  structurally-zero mean-flow rows of the next column.
  `assert_coupling_blocks_allclose` calls
  `setStates(real(getStates()))` before each CS column to clear it — this
  is required, not optional; removing it reintroduces a spurious ~2.15e-8
  failure on `dR[meanflow]/dw[gamma]` and `dR[meanflow]/dw[reThetat]`.

**Changing the mesh**: edit
`sagrGridFile` / `sagrRestartFile` (and the AeroProblem conditions if
they differ) in `reg_sagr.py`, then `./run_sagr_tests.sh genw` (if you need
a fresh converged `w`) → `train` (rebuild JSON) → run. Everything routes
through that one config module.

**Generating `w` (the restart state)** is non-standard (no download server),
so its generator lives in `dev/generate_sagr_restart.py` and is documented
in `dev/README.md`. It reads `reg_sagr.sagrBaseOptions`, so it always
matches the test case.

**What this ladder does and does not prove.** It proves the differentiated
code computes the derivative of the primal (partials: `dR/dw`, `dR/dXv`,
`dR/d{aeroDV}`, plus output partials) — validated in all three AD modes
against complex-step (crossflow OFF in the standing suite; crossflow was
validated once on 2026-07-22, see the note above). The *total* sensitivity
`dF/dX` through the assembled adjoint solve is now covered too:
`test_adjoint_sagr.py` is live as the `adjoint` stage and passes (real 4/0;
CS at rtol/atol 5e-8 with mach+drag non-blocking — see
`../task-log/2026-07-23-sagr-full-adjoint-test.md` and
`../task-log/2026-07-24-sagr-cs-tolerance-nonblocking.md`).

---

This is the low-level, raw-API verification campaign — separate from (and
a prerequisite trusted-input for) the higher-level `evalFunctionsSens`-based
campaign (the AR5 dcd/dmach mismatch was re-characterized in
`../task-log/2026-07-24-sagr-cs-tolerance-nonblocking.md`: mach CS does not
settle; cause TBD, non-blocking).

This ladder answers one question at three levels of the AD stack: **does
the differentiated code actually compute the derivative of the primal
code**, independent of any specific functional or DV. It does not check
whether the adjoint *total* sensitivities match a physical/trained
reference — that is `test_adjoint_sagr.py`'s job (the `adjoint` stage). Read
`VERIFICATION/VERIF_01_debugging_derivatives.md` first for the generic
FD→dot-product→CS theory behind why each stage below is structured the
way it is.

## Why three stages, in this order

1. **Reverse ↔ forward consistency (dot-product test).** Cheapest, most
   fundamental check: for the same linearized operator, forward-mode
   (`_d`, computes `dR/dw · wDot`) and reverse-mode (`_b`, computes
   `dR/dw^T · psi`) must satisfy `<wDot, dR/dw^T psi> = <psi, dR/dw wDot>`
   to machine precision. This only checks fwd/rev *bookkeeping* agree with
   each other — it cannot catch a bug present identically in both modes.
2. **Reverse vs. fast-reverse consistency (`_b` vs `_b_fast`).** ADflow
   ships two reverse-mode variants; `_fast_b` is a stripped-down Tapenade
   output (`autoEditReverseFast.py`) used for performance. If this
   diverges from `_b`, `autoEditReverseFast.py`'s push/pop stripping is
   the first suspect (this exact failure mode broke SST's `fast_b`
   upstream — see `../ADFLOW_BASE/ADFLOW_07_sst_dev_lessons.md`).
3. **Forward verification test (3-way AD vs. FD vs. CS).** Only once 1
   and 2 agree does it make sense to ask whether the *shared* fwd/rev
   linearization is actually correct — i.e. matches the true derivative of
   the primal. Complex-step (h=1e-40, no subtractive cancellation) is the
   authoritative ground truth here; single-step real FD is included only
   for context (expected to disagree at ~1e-1 to ~1e-2 relative error due
   to truncation/kinks, not a failure).

## Results & provenance (condensed)

The detailed per-script reproduction (the original one-off `dev/` scripts
`sanity_check_*.py` / `check_3way_fwd*.py`, with their `--mesh` /
`--turbmodel` / `--crossflow` flags) has been trimmed. Those scripts live in
`tests/reg_tests/dev/` and still run; day-to-day use `run_sagr_tests.sh`.

**Key pass numbers** (from the original 2026-07-22 runs, mostly on the AR5
mesh; the standing suite has since moved to the tutorial wing. Pass
criterion is AD == CS to full displayed precision, FD is context only):

| Case | AD | CS | agrees |
|---|---|---|---|
| plain SA, tutorial-wing | 5.1750905670e+10 | 5.1750905670e+10 | 10 sig figs |
| plain SA, AR5 | 3.2957137433e+13 | 3.2957137432e+13 | full |
| SA-GR `nw=8`, AR5, crossflow OFF | 3.3256409775e+13 | 3.3256409775e+13 | 11 sig figs |
| SA-GR `nw=8`, AR5, crossflow ON | 3.3255762835e+13 | 3.3255762835e+13 | full |

Stage 1 (dot products) and Stage 2 (`_b` vs `_fast_b`) pass to `tol=2e-10` /
`atol=1e-16` in every case. FD disagrees at O(1e-1) at a single step but
converges to the AD/CS value as O(h) — see the FD `@expectedFailure` note in
"Canonical way to run" for why the registered suite treats FD that way.

**Stage-2 `_fast_b` reThetat fix (2026-07-22).** The first SA-GR run had
Stage 2 FAIL on the reThetat-seeded block (`_fast_b` wrong by rel_err 1.94,
**wrong sign**) while `_b` and `_d`/CS agreed to ~1e-14 — an in-place
`smoothMinMax` double-clamp of `lambdaThetaLocal` whose Tapenade push/pop was
stripped by `autoEditReverseFast.py`. Fixed by rewriting the primal with
distinct targets (`lambdaThetaRaw`→`lambdaThetaClamped`→`lambdaThetaLocal`) +
Tapenade regen; `_fast_b` then matched to 1.27e-13. Full root-cause in
`../task-log/2026-07-22-fastb-retheta-lambdatheta-fix.md`.

**Relationship to the AR5 CS mismatch.** A *different, higher-level* issue
(re-characterized in `../task-log/2026-07-24-sagr-cs-tolerance-nonblocking.md`:
mach CS does not settle; cause TBD, non-blocking): `evalFunctionsSens`'s
`dcd/dmach` on AR5 (adjoint −0.537 vs CS −231.1). That is **not** a
contradiction of the PASS above — this ladder checks the raw `dR/dw`
Jacobian-vector product (confirmed correct); the `dcd/dmach` gap is
downstream in the *total* sensitivity through a whole-solve re-solve, traced
to an absolute-convergence precision floor on this stalling mesh, not a
linearization bug. Validated partials ≠ validated gradient.

**Complex-build recipe (needed to run the CS stage).** The complex
`libadflow_cs.so` must be built with `PETSC_ARCH=complex-debug`; building
against the default `real-debug` PETSc (or reusing `.mod` files from a prior
real-arch attempt) produces spurious `COMPLEX(8)→REAL(8)` interface errors in
`fortranPC.F90` / `adjointAPI.F90` / `surfaceUtils.F90` / etc. that are
toolchain artifacts, not source problems. Full recipe:

```bash
git clean -fdx src_cs/            # -x is required: src_cs is gitignored
PETSC_ARCH=complex-debug make -f Makefile_CS
```

## What's not yet covered

- **Total sensitivity `dF/dX`** — no longer a gap (updated 2026-08-12):
  `test_adjoint_sagr.py` is live as the `adjoint` stage and passes (real
  4/0; CS at rtol/atol 5e-8 with mach+drag non-blocking — task-log
  2026-07-23/24). The AR5 dcd/dmach mismatch was re-characterized in
  `../task-log/2026-07-24-sagr-cs-tolerance-nonblocking.md` (mach CS does
  not settle; cause TBD, non-blocking).
- Crossflow `D_scf` is **not** exercised by the standing suite (crossflow
  OFF since the 2026-07-24 default flip); it was validated once on
  2026-07-22 with crossflow ON on the AR5 state and remains reproducible
  via the `dev/` scripts' `--crossflow` flag.
