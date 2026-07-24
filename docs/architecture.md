# Architecture — SA-BCM Differentiable Transition Model

## 1. Solver architecture

SA-BCM is **not** a new turbulence model and **not** a new transport equation. It is a
multiplicative modifier applied inside the existing single-equation SA transport equation
(`turbulenceModel = "SA"`, enum `spalartAllmaras = 2`, `src/modules/constants.F90:128-135` —
there is no separate `SABCM` enum value). Activation is a boolean, `use_SABCM`, layered on top
of plain SA.

**Consequence for state-vector layout:** `nw` is unchanged from plain SA. There is no new
`itu*` index, no new residual block, no new off-diagonal coupling block to isolate in Jacobian
tests. Where the upstream γ-Re̅θt-style test harnesses split the state vector into
`meanflow`/`nuTilde`/`<new vars>` blocks, SA-BCM tests instead directly compare the **modified**
Jacobian terms (those that pick up a `tTgamma` or `dtTgamma` factor) against **plain SA**'s
unmodified terms with `use_SABCM=False` — see `tests/reg_tests/reg_bcm.py`.

### Code path

- `src/turbulence/sa.F90`, `subroutine saSource`:
  - Line 294: `tTgamma = one` (default multiplier — SA-BCM off ⇒ no-op, recovers plain SA exactly).
  - Lines 296–351: `if (use_SABCM) then` — computes vorticity-based `Re_theta` vs.
    `Re_theta_c(SABCM_TU)`, blends via tanh (`SABCM_Exp = .false.`, the manuscript's smoothed
    formulation) or exp-sqrt (`SABCM_Exp = .true.`, the original SA-BC formulation), stores the
    result into the block-pointer field `Tgamma(i,j,k)`.
  - Line 350: when active, `ft2` is forced to zero — SA-BCM replaces `ft2`'s role entirely.
  - Lines 359, 361: SA production term multiplied by `tTgamma`:
    `term1 = tTgamma * rsaCb1 * (one - ft2) * ss`, and the diagonal `term2` similarly.
  - Lines 366–413 (`#ifndef USE_TAPENADE`): hand-written Jacobian path folds a `dtTgamma`
    derivative term into the existing `qq` diagonal — still modifying the existing
    linearization, not adding a new one.
- `src/NKSolver/blockette.F90:1004-1250`: the non-Tapenade ANK/DADI mirror of the same logic —
  must be kept consistent with `sa.F90` by hand (it is not Tapenade-differentiated separately;
  see `adjoint-trace.md`).

### Runtime options

From `adflow/pyADflow.py:5685-5692` (defaults) and `pyADflow.py:6094-6101` (namelist mapping),
`src/modules/inputParam.F90:587-592` (Fortran declarations), `src/f2py/adflow.pyf:1177-1184`
(f2py exposure):

| Python option | Fortran variable | Default | Role |
|---|---|---|---|
| `use_SABCM` | `use_SABCM` (logical) | `False` | Master switch. `False` ⇒ `tTgamma≡1`, byte-identical to plain SA. |
| `SABCM_Exp` | `SABCM_Exp` (logical) | `False` | `True` ("hard") = the paper's own formula verbatim, γ=1-exp(-(√Term1+√Term2)) (AIAA 2020-2714 Appendix / 2020-2706 Eq. 3). `False` ("smooth") = a tanh blend that is **not** in either paper — a deliberate, user-chosen smoothing kept intentionally (not a deviation to "fix"); confirmed by paper-vs-code verification, 2026-07-24. |
| `SABCM_Const1` | `SABCM_Const1` | `0.002` | χ₁ — scales `Re_theta_c` in Term1's denominator (paper: χ₁, both papers agree at 0.002). |
| `SABCM_Const2` | `SABCM_Const2` | `0.02` | χ₂ — Term2 = `(fv1*chi)/SABCM_Const2` = `(mu_T/mu)/χ₂`; **must be the eddy-viscosity ratio `fv1*chi`, not raw `chi`** — this is the trap both papers explicitly flag as "wrong by a factor of 250" if confused with the χ₂=5.0 used in the *unrevised* Term2 derivation. Confirmed correct in code, 2026-07-24. |
| `SABCM_TU` | `SABCM_TU` | `0.5` | Freestream turbulence intensity Tu∞ (%), feeds the shared `Re_theta_c = 803.73(Tu∞+0.6067)^-1.027` correlation (identical in both papers). |
| `SABCM_S0_tanh` | `SABCM_S0_tanh` | `0.5` | S₀ — "smooth" (tanh) variant's blend centering value. Not from either paper (see `SABCM_Exp` row) — part of the deliberate user smoothing, calibrated, do not retune casually (CLAUDE.md rule 5). |
| `SABCM_fsmooth` | `SABCM_fsmooth` | `0.08` | f_smooth — "smooth" (tanh) variant's blend width. Same caveat as `SABCM_S0_tanh`. |
| `SABCM_maxsmooth` | `SABCM_maxsmooth` | `50.0` | KS-function (log-sum-exp) aggregation parameter replacing the paper's `max(Term1, 0)` kink — the other deliberate smoothing, applies to **both** `SABCM_Exp` variants (unlike S0_tanh/fsmooth, which are tanh-variant-only). |
| `NKLSRelax` | `NK_LSRelax` (logical, `NKSolver` module) | `False` | **Not SA-BCM-specific** — generic NK line-search relaxation (`LSCubic`'s Armijo `alpha` 1e-2→1e-3, turb-blowup pre-limit factor 2.0→3.0). Added 2026-07-24 from a sibling repo's NK-convergence-mods review as the one idea that transfers without needing new transported-state machinery. Off by default: NK has no per-node physicality check, so only turn on if a run is observed pinned at `minlambda` for 100+ consecutive iterations. |

**Hard constraint**: `src/inputParam/inputParamRoutines.F90:3413-3423` — SA-BCM with plain SA
and `equations=RANS` requires `useApproxWallDistance=.true.`, else ADflow calls `terminate(...)`.
Always set this option when `use_SABCM=True`.

**Note on citations above**: earlier revisions of this table cited an interim "differentiable
reformulation manuscript" that has since been replaced by the two actual AIAA source papers in
`docs/papers/` (2020-2706, 2020-2714). The physics/constant claims above were re-verified
directly against those papers 2026-07-24 (see `docs/adjoint-trace.md` and
`docs/task-log/2026-07-24-paper-verification-and-nk-relax.md`); any other file in this KB still
saying "manuscript Eq. N" is citing the retired document and should not be trusted for equation
numbers without re-checking against `docs/papers/` directly.

### Paper-symbol ↔ code-flag lookup

| Manuscript symbol | Meaning | Code |
|---|---|---|
| γ | intermittency factor (production multiplier) | `tTgamma`, stored per-cell in `Tgamma(i,j,k)` |
| term1 | natural-transition trigger (manuscript Eq. 4) | `tterm1` in `sa.F90` |
| term2 | bypass/turbulence-sustaining term (manuscript Eq. 6) | `tterm2` in `sa.F90` |
| Rθ, Rv | vorticity Reynolds number | `Re_theta`, computed from vorticity `S`, wall distance `d`, `ρ`, `μ` |
| Rθc | critical Reynolds number correlation | `Re_theta_c`, function of `SABCM_TU` |
| χ1, χ2 | term1/term2 scaling constants | `SABCM_Const1`, `SABCM_Const2` |
| Tu∞ | freestream turbulence intensity | `SABCM_TU` |
| S0 | tanh blend center | `SABCM_S0_tanh` |
| f_smooth | tanh blend width | `SABCM_fsmooth` |
| ρ (KS aggregation) | KS-function sharpness for `max(f,0)` | `SABCM_maxsmooth` |
| k_max | KS anchor, treated as constant under AD | computed inline, `∂k_max/∂f = 0` enforced by construction |

## 2. Diagnostics

Per-cell `γ` (`tTgamma`) is written to the volume CGNS via the `Tgamma` block array
(`src/modules/block.F90:515`, `src/modules/blockPointers.F90:124`) — this is the only sanctioned
SA-BCM diagnostic output channel (CLAUDE.md rule 4).
