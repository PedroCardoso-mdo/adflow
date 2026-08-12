# ANK / NK Modifications for the SA-γ-Re̅θt Transition Model — Audit Report

*Historical audit (as of 2026-07-23, branch `sa_gamma_rethetha_paper_solver`,
since merged). Line numbers have drifted — locate by symbol.*

**Branch:** `sa_gamma_rethetha_paper_solver` vs `main`
**Files audited:** `src/NKSolver/NKSolvers.F90` (+1335), `src/NKSolver/blockette.F90` (+716)
**Date:** 2026-07-23

This report answers three questions:

1. **What** was changed in ANK/NK for the new model.
2. Is the code **still complexifiable** (complex-step / AD build).
3. Does it reduce to **stock ADflow** when the new options are off, or when
   the transition model isn't the active turbulence model.

---

## 1. What was added

All ANK/NK work is a *convergence-acceleration bundle* layered on top of the
existing solver — it never rewrites the stock Newton/pseudo-transient path.
Three mechanisms from Piotrowski & Zingg (2020):

| Mechanism | Paper | New routines (`NKSolvers.F90`) |
|-----------|-------|--------------------------------|
| **Column scaling** of the 13-order-of-magnitude SA-GR state | §IV.B | `getNKColScale`, `applyNKColumnScaling`, `applyANKColumnScaling`, `applyTurbPCColumnScaling`, `getTurbColScale`, `setW*ANK*Scaled`, `getFullColScale` |
| **Source-dt restriction** (Eq. 59) on the transition rows | §IV.B.3 | `applyNKSrcDtDiagonal`, `NKSrcDtMatMult`, reactivation-on-backtrack in `NKStep` |
| **Algorithm 2** per-node γ/Re̅θt bounds damping + soft-damping diagnostics | §IV.B.2 | `applyNKAlgorithm2Damping` |
| Residual/row autoscale proxies (Eq. 58) | §IV.B | `computeNKResidualAutoscale` (off by default) |
| NK stall-escape (EW rtol cap, auto-disable latch) | — (empirical) | logic in `NKStep`, `getEWTol` |

`blockette.F90` adds a `case (spalartallmarasnoft2gammaretheta)` to the
residual/Jacobian assembly dispatch (`saGammaRethetaSource/Advection/Viscous/
ResScale`). It is a new `select case` branch — other models fall through to
their existing cases untouched.

> Note: the blockette SA-GR kernels audited here were later found
> sign-flipped and re-synced on 2026-07-24
> (`../task-log/2026-07-24-blockette-sagr-residual-sync.md`;
> `test_blockette_sagr.py` now guards them).

---

## 2. Complexifiable? — **Yes.**

- All new declarations use `real(kind=realType)` / `alwaysRealType`, the same
  kinds the complex build promotes. No hard `complex`/`cmplx`/`aimag`/`dble`.
- The only real-part extractions (`real(...)`) sit inside the **physicality
  step-limiter** and are correctly wrapped in `#ifndef USE_COMPLEX / #else /
  #endif` pairs — the identical pattern stock ADflow already uses for its flow
  and SA physicality checks (`NKSolvers.F90:4110–4165`, `:4317`). Step-limiting
  logic is deliberately non-differentiated, so dropping the imaginary part
  there is correct, not a bug.
- The new γ/Re̅θt physicality block follows that same `#ifdef` convention it
  was copied from.
- **Evidence it actually builds/runs complex:** commit `ffba7f9c` ("Fix
  complexify-build bugs blocking the complex ADflow build") and `72325c8c`
  ("Stage-3 CS re-confirmed on freshly rebuilt complex lib") — the complex
  library was rebuilt and the complex-step derivative check passed. See
  `docs/VERIFICATION/three-stage-verification.md`.

**Caveat:** the *soft-damping print diagnostics* (`applyNKAlgorithm2Damping`)
do `print *,` counts — harmless to complex, but they are runtime chatter, not
guarded by a verbosity flag. Not a correctness issue.

---

## 3. Same result as stock ADflow?

### 3a. When the active turbulence model is **not** SA-γ-Re̅θt → **identical.**

Every new code path is gated on
`turbModel == spalartallmarasnoft2gammaretheta`. Verified at **all** call
sites:

| Line | Call | Guard |
|------|------|-------|
| 414 | `computeNKResidualAutoscale` | `turbModel==SAGR .and. transitionNK .and. transitionNKActive .and. transitionResidualAutoscale` |
| 447 | `applyNKColumnScaling` | `turbModel==SAGR .and. transitionNK .and. transitionNKActive` |
| 460 | `applyNKSrcDtDiagonal` | same + `transitionSrcDtRestrict` |
| 752 | `applyNKAlgorithm2Damping` | `turbModel==SAGR .and. transitionNK .and. transitionNKActive` |
| 2598/2730 | `applyANKColumnScaling` | `turbModel==SAGR .and. transitionNK .and. ANK_coupled` |
| 3090 | `applyTurbPCColumnScaling` | `turbModel==SAGR .and. transitionNK` |

Where a routine is called *unconditionally* (`getTurbColScale`/
`setWANKTurbScaled`, ~3273), the scale factor `cs` **collapses to 1** for
non-transition models (`cs = turbResScale` only for SA-GR, `cs = 1` otherwise),
so the arithmetic is a provable no-op — the original `turbResScale` factor is
recovered exactly. This is stated in the in-code comments and confirmed by
inspection.

`blockette.F90`'s added `case` never executes for other models.

➡️ **For any non-SA-GR run (Euler, SA, SST, k-ω…), the branch is
byte-for-byte equivalent to stock ADflow on the ANK/NK path.** This is the
strong invariant and it holds unconditionally. *(Code-inspection claim; the
proposed non-SA-GR regression run was never recorded as executed.)*

### 3b. When SA-γ-Re̅θt **is** active — depends on option defaults.

The bundle is **on by default when the transition model is selected.** To get a
"pure Newton, no acceleration bundle" run for comparison, flip the master
switch. Defaults (`src/modules/inputParam.F90:305–371`):

| Option | Default | Effect at default |
|--------|---------|-------------------|
| `transitionNK` | `.true.` | **master switch** — set `.false.` to disable the entire bundle |
| `transitionSrcDtRestrict` | `.true.` | Eq. 59 source-dt restriction **active** |
| `transitionSrcDtLimit` | `0.9` | — |
| `srcDtDeactivateIters` | `5` | — |
| `transitionResidualAutoscale` | `.false.` | Eq. 58 S_a proxy **off** (unvalidated) |
| `transitionRowVolScale` | `.false.` | Eq. 58 row-vol scale **off** (unvalidated) |
| `transitionNKAutoDisableTol` | `0.0` | auto-disable latch **off** |
| `transitionNKStallRtolCap` | `1.0` | stall rtol-cap **off** |

So with the transition model active and default options, **column scaling +
Eq. 59 + Algorithm 2 damping are engaged** — this is intended (the raw system
is too stiff to converge without them). It is *not* the same iteration path as
"stock ADflow running SA-GR with the bundle off"; setting `transitionNK=.false.`
gives that baseline.

Note `transitionNKActive`/`nkStallCount`/`noBacktrackCount` are **internal
runtime latches**, reset each solve — not user knobs.

---

## Bottom line

- **Complexifiable:** yes — real-part casts are `#ifdef`-guarded exactly like
  stock ADflow; complex build + CS check already re-confirmed on this branch.
- **Reduces to stock ADflow when SA-GR is not the model:** yes,
  unconditionally — all paths gate on `turbModel`, and the one always-called
  scaling path is a mathematical no-op (`cs=1`) for other models.
- **When SA-GR *is* active:** the acceleration bundle is on by default; use
  `transitionNK=.false.` to recover the unaccelerated baseline for A/B
  comparison. The two unvalidated Eq. 58 scalings are already off by default.

### Suggested verification runs
1. **Non-SA-GR regression** — any existing SA/SST reg test: results must match
   `main` refs bit-for-bit (this is the guarantee in 3a).
2. **SA-GR bundle on/off** — same case with `transitionNK=.true.` vs
   `.false.`: same *converged* state, different iteration history.
3. **Complex-step** — rerun the Stage-3 CS check
   (`docs/VERIFICATION/three-stage-verification.md`).
