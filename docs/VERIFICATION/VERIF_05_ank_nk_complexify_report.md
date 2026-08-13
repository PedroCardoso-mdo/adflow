# Is the SA-γ-Re̅θt ANK/NK Code Complexified? — Pre-CS-Adjoint Audit

*Historical audit (as of 2026-07-23, branch `sa_gamma_rethetha_paper_solver`,
since merged). Line numbers have drifted — locate by symbol.*

**Purpose:** document complex-step (CS) safety of the transition-model code
**before** running the full adjoint CS verification, so any CS mismatch can be
attributed to the adjoint wiring rather than an un-complexifiable primal.
**Branch:** `sa_gamma_rethetha_paper_solver` · **Date:** 2026-07-23

**Verdict: Yes — the code is complexified and CS-safe.** The complex library
already builds, imports, and the Stage-3 AD = CS check matches to full
precision. Details and scope below.

> Operational caveat found later the same day: the complex build excludes
> the Tapenade AD routines, so ANKADPC/NKADPC must be False in CS runs and
> CS accuracy caps ~1e-8
> (`../task-log/2026-07-23-sagr-full-adjoint-test.md`).

---

## 1. What CS actually exercises (scope)

The complex-step adjoint check perturbs the state `w → w + i·h` and reads the
imaginary part of the **residual** `R(w)`. So the code that MUST be
complex-exact is the **residual-evaluation path**:

- `src/turbulence/saGammaRetheta.F90` (Source / Advection / Viscous / ResScale)
- `src/turbulence/turbUtils.F90` (correlations + `smoothMinMax`)
- `src/turbulence/turbBCRoutines.F90` (transition BCs)
- `src/NKSolver/blockette.F90` (the `saGammaRetheta*` residual assembly branch)

The **solver-iteration** code (column scaling, Eq. 59 source-dt, Algorithm 2
damping, physicality line-search in `NKSolvers.F90`) is **not** differentiated
by CS — it only has to *compile* under `-DUSE_COMPLEX`. It is handled
separately in §3.

---

## 2. Residual path — complexify-clean

**No manual complex handling, and none needed.** The four residual files
contain **zero** `#ifdef USE_COMPLEX` blocks and **zero** raw `real()` /
`dble` / `cmplx` / `aimag` casts. They are written in kind-generic style
(`real(kind=realType)`), which the complexify preprocessor promotes to complex
wholesale — the correct pattern.

Intrinsic inventory in `saGammaRetheta.F90` (all complexify-overloaded, all
CS-safe):

| Intrinsic | Count | CS status |
|-----------|-------|-----------|
| `max` / `min` | 132 / 11 | overloaded (compare real parts) ✓ |
| `sqrt`, `exp`, `log`, `tanh`, `abs` | 30/20/3/9/6 | analytic overloads ✓ |
| `mod` | 10 | **integer index arithmetic only** (`mod(ii,nx)`) — not on complex data ✓ |

The smooth model deliberately routes the non-differentiable `min`/`max` of the
physics through `smoothMinMax` (`turbUtils.F90:2382`), so the branch
discontinuities that would poison a step-differentiation are already smoothed —
CS sees a genuinely differentiable function.

### The two complexify bugs that existed were found and fixed

Commit `ffba7f9c` ("Fix complexify-build bugs blocking the complex ADflow
build", 2026-07-20) — the complex build had **never** been stood up on this
branch before; two real bugs blocked it, both now fixed:

1. **Bare real literals passed to `smoothMinMax`** (e.g. `0.4_realType`) across
   `saGammaRetheta.F90` / `turbUtils.F90` / `blockette.F90`. Fortran will not
   auto-promote a REAL literal actual argument to COMPLEX at a checked module
   interface → compile failure under `-DUSE_COMPLEX`. Fixed by moving them to
   named `parameter`s in `paramTurb.F90` (the `rsaGRpmax`/`rsaGRpmin` pattern).
2. **`blockette.F90`**: a line-continuation right after an `==` comparand broke
   complexify's automatic `== → .ceq.()` parenthesization. Fixed by moving the
   continuation after `.and.`.

Both are exactly the class of issue that only surfaces at complex-build time —
they're now resolved, which is why this audit can assert the path builds.

---

## 3. Solver path (`NKSolvers.F90`) — compiles complex, correctly `#ifdef`-guarded

The solver acceleration code is not CS-differentiated, but it still has to
build in the complex library. It handles this the **stock-ADflow way**: the
21 `USE_COMPLEX` occurrences guard the physicality/step-limiter real-part
extractions:

```fortran
#ifndef USE_COMPLEX
    ratio = abs(wvec_pointer(ii) / (dvec_pointer(ii) + eps)) * ANK_physLSTol
#else
    ratio = abs(real(wvec_pointer(ii)) / real(dvec_pointer(ii) + eps)) * real(ANK_physLSTol)
#endif
```

This is correct: step-length limiting is control logic that must key off the
*real* part only, never propagate the imaginary perturbation. The new γ/Re̅θt
physicality block (`:4110–4165`, `:4317`) was copied from and follows this
exact convention. The Algorithm-2 damping (`applyNKAlgorithm2Damping`) uses
only `realType` and comparisons — no raw casts — so it is complex-clean too.

---

## 4. Evidence it already works

- **Complex build + import:** `ffba7f9c` — "complex ADflow builds and imports
  cleanly"; `test_adjoint.py` harness passes 4/4 real + 2/2 complex on
  `rans_tut_wing`.
- **Fresh complex lib, SA-GR ladder** (`72325c8c`, 2026-07-22) — rebuilt
  `libadflow_cs.so` from scratch (`git clean -fdx src_cs/` +
  `PETSC_ARCH=complex-debug`), reran the SA-GR verification ladder:
  - Stage 2 (reThetat block) `max|_b − _fast_b| = 2.328e-10` — PASS
  - Stage 3 **AD = CS = 3.3256409775e+13** on the fresh complex lib — no stale
    artifacts.

Full procedure/results: `docs/VERIFICATION/VERIF_00_three_stage_verification.md`.

---

## 5. Bottom line for the full-adjoint CS test

- The residual path is complexify-clean and the two historical complex-build
  bugs are fixed → **CS differentiates a genuinely complex-analytic primal.**
- The solver path compiles complex via the standard `#ifdef` guards and is
  outside the differentiated chain.
- A localized 1-point SA-GR CS check already matches AD to full precision.

**Therefore, if the full-adjoint CS test shows a mismatch, suspect the adjoint
seed/wiring or an inactive-derivative declaration — not an un-complexifiable
primal.** Residual-path files carry no `#ifdef` you could have mis-set; the
only complex-conditional logic lives in the non-differentiated solver.

### Pre-flight checklist before running full-adjoint CS
1. Rebuild `libadflow_cs.so` from the **current** primal (`git clean -fdx
   src_cs/`) — avoid a stale complex lib (the trap `72325c8c` guards against).
2. Confirm the CS step `h` is deep enough (≤1e-30 typical) — CS has no
   subtractive cancellation, so use a tiny `h`.
3. Keep the solver-acceleration options irrelevant to the check: CS reads the
   residual, not the iteration path, so `transitionNK` on/off must not change
   the CS derivative (a useful sanity assertion).
