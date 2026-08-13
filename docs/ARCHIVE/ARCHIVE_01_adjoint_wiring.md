# 06 — Adjoint wiring & hand-written derivative audit (SA-γ-R̃e_θt)

> **⚠️ STATUS 2026-08-12:** all findings in §5 are RESOLVED or
> premise-changed — F1 (vortlimd) fixed & Tapenade regenerated (verified:
> no `vortlimd = 0.0_8` remains); F5's premise inverted
> (`transitionCrossflow` now defaults False; `inputParam.F90:327`'s Fortran
> default is still `.true.` but Python pushes False); F6 fixed; F7 actually
> fixed 2026-07-24 (blockette resync — not 2026-07-07 as the addendum
> says). `turbResScale` default is now `[1e4, 0.1, 1e-4]`, not
> `[1e4, 10, 1e4]`. AD files are committed. Read this file as a historical
> audit; line numbers have drifted.

**Date:** 2026-07-07. Branch `sa_gamma_rethetha`, working tree as of this audit
(six regenerated AD files uncommitted). Audit-only; no code changed.

> **RESOLUTION ADDENDUM (2026-07-07, same day, post-audit):** all findings
> below were subsequently fixed and the fixes verified (build compiles,
> `libadflow` imports). Tapenade was rerun twice via `./AD_I.sh` (before and
> after the primal edits). Per finding:
> - **F1 — fixed by regeneration.** `vortlimd = 0.0_8` gone; `uinfd`/`muinfd`
>   now propagate into `vortlimd` (`saGammaRetheta_d.f90`, forward and
>   reverse). §3's "TAPENADE NEEDED" status is cleared.
> - **F8 — fixed at the root.** `turbUtils.F90:turbAdvection` now forces
>   `secondOrd = .false.` for SA-GR when `transitionFirstOrderUpwind` (P&Z
>   §IV.A), so DADI, coupled ANK/NK, and every AD master sweep discretize the
>   identical advection residual for any `orderTurb` setting. Verified present
>   in regenerated `turbutils_{d,b,fast_b}`. The old per-caller mechanism —
>   `saGammaReTheta_block` temporarily swapping the global `orderTurb` around
>   its `turbAdvection` call — was subsequently **removed** (2026-07-08): it
>   was redundant with the in-routine guard, covered only the decoupled DADI
>   path, and mutating a global module variable is invisible to the AD master
>   sweeps and not thread-safe. Single point of truth is now the guard inside
>   `turbAdvection` itself. Note the guard applies one `secondOrd` flag to all
>   three SA-GR equations (ν̃ included), same as the old swap did.
> - **F2 — fixed.** `qq(1,2)` drops the `cb1(1−ft2)ss` term in the
>   `approxSA∧transitionUseApproxSA` branch, matching the residual (which sets
>   `term1 = zero` there — the Jacobian was differentiating a term absent from
>   the residual). Root cause: `qq(1,1)` reused the branched intermediates
>   `term1`/`term2` and was consistent; `qq(1,2)` re-expanded from raw
>   primitives and missed the branch. LHS-only (`#ifndef USE_TAPENADE`) —
>   convergence in approx-SA mode, not gradients. `evalSrcJacBlock` A(1,2)
>   deliberately stays unbranched: it linearizes the true source for
>   `computeSrcLambda`'s dt restriction.
> - **F3 — fixed exactly.** `qq(2,3)` now includes `∂E_γ/∂R̃ =
>   −c_a2·e^{−r_T}·dfOnset_dReT·Ω̃·γ(c_e2γ−1)` (entry is now the complete
>   `−∂(P_γ−E_γ)/∂R̃`); `evalSrcJacBlock` A(2,3) FD perturbs the full
>   `pGamma − eGamma`.
> - **F5 — resolved by documentation** (kept zero deliberately): the paper's
>   block treatment doesn't require it and a nonzero A(3,1) would break the
>   triangular structure `computeSrcLambda` relies on; comments no longer
>   claim exactness. Appendix A.7 remains the record of the true partial.
>   Eigenvalue impact assessed (2026-07-08) and accepted: A13 = ∂S_ν̃/∂θ is
>   exactly zero, so A31 enters det(A−λI) only via the product A12·A23·A31
>   (a full ν̃→γ→θ→ν̃ loop) — the λ error is third-order in already-small
>   couplings, and `srcLambda` only needs a conservative max(0, λ) for the
>   dt restriction, not an exact spectrum. Fixing it would mean
>   differentiating the whole crossflow chain plus a per-cell cubic solve.
>   If DADI ever stalls in crossflow-active cells under source-dt limiting,
>   the remedy is Gershgorin-style padding with an FD |A31|, not exactness.
> - **F6 — fixed.** `evalSrcJacBlock` uses `smoothMinMax(vortMag, vortLim,
>   rsaGRpmin)` like the residual. The hard `min` had a derivative kink
>   (1→0) right where the vorticity cap engages — high-shear cells, exactly
>   where the Eq. 59 dt restriction matters — so `srcLambda` was the
>   eigenvalue of a source the solver doesn't converge. dt-restriction only,
>   never gradients. Deliberate asymmetry with F2: no `approxSA` branching
>   here — `evalSrcJacBlock` linearizes the *true* source (protects the
>   physical time step), while `qq` must match the residual DADI iterates.
> - **F7 — fixed.** Blockette SA-GR source kernel synced (refLenTrans branch
>   + full crossflow D_scf) and its advection kernel honors
>   `transitionFirstOrderUpwind`. The `pyADflow` force-off remains as a belt.
>   The kernels were a stale (2026-05-27) copy of the residual missing two
>   later physics changes — dead today only because of the force-off, but a
>   "solver converges A, adjoint linearizes B" landmine (same class as F8)
>   for whoever re-enables blockettes. Synced rather than deleted so that
>   re-enabling starts from correct kernels; blockettes stay off until
>   deliberately validated.
>
> Numerical verification (FD-vs-AD partials tests) is still the pending
> decisive check, as stated below.

**Scope and epistemic status.** This is a *structural and symbolic* audit: I
classified every derivative file by provenance, re-derived every hand-written
partial by hand (Appendix A), checked primal↔AD drift by content comparison,
and traced the adjoint state-vector/coupling wiring with line-level evidence.
**No FD-vs-AD or dot-product tests were run.** Where I say "correct" below I
mean "the code matches my independent derivation / the current primal text";
numerical verification (partials tests) remains the decisive check and is
still pending. Complements `docs/audits/adjoint_audit_2026-07-07.md` (visual
wiring pass) and `docs/VERIFICATION/adjoint-trace.md`. Numerical
verification has since been run — see `docs/VERIFICATION/three-stage-verification.md`.

---

## 1. Provenance: Tapenade-generated vs hand-written, file by file

Every derivative file relevant to the SA-GR residuals carries the literal
header `generated by tapenade (inria, ecuador team) / tapenade 3.16 (develop)
- 22 aug 2023` (that date is the Tapenade *build*, not the generation date).
The uncommitted working-tree modifications to the six SA-GR AD files are
regeneration output, not hand edits — verified in §3.

| Primal residual file | Forward (`_d`) | Reverse (`_b`) | Reverse-fast (`_fast_b`) | Provenance |
|---|---|---|---|---|
| `src/turbulence/saGammaRetheta.F90` (`Source`, `Viscous`, `ResScale`) | `outputForward/saGammaRetheta_d.f90` | `outputReverse/saGammaRetheta_b.f90` | `outputReverseFast/saGammaRetheta_fast_b.f90` | **Tapenade** (all three headers verified) |
| `src/turbulence/turbUtils.F90` (`turbAdvection`, correlations, `smoothMinMax`) | `turbUtils_d.f90` | `turbUtils_b.f90` | `turbUtils_fast_b.f90` | **Tapenade** |
| `src/turbulence/turbBCRoutines.F90` | `turbBCRoutines_d.f90` | `turbBCRoutines_b.f90` | `turbBCRoutines_fast_b.f90` | **Tapenade** |
| `src/initFlow/initializeFlow.F90` (`referenceState`, `wInf(itu2/itu3)`) | `initializeFlow_d.f90` | `initializeFlow_b.f90` | `initializeFlow_fast_b.f90` | **Tapenade** |
| `src/bcdata/BCData.F90` | `BCData_d.f90` | `BCData_b.f90` | temp only (`temp_reverse_fast/`) — same as upstream SA layout | **Tapenade** |
| `src/turbulence/sa.F90` (unmodified baseline) | `sa_d.f90` | `sa_b.f90` | `sa_fast_b.f90` | **Tapenade** |
| `src/adjoint/masterRoutines.F90` (AD driver: `master`, `master_d`, `master_b`, `master_fast_b`, `block_res_state_d`) | — | — | — | **Hand-written wiring** (calls generated routines; contains no derivative math) |

**Hand-written derivative *math* (not part of the adjoint operator):**

| Code | Location | Consumer |
|---|---|---|
| Inline 3×3 source Jacobian `qq(i,j,k,1:3,1:3)` | `saGammaRetheta.F90:660–822` (`#ifndef USE_TAPENADE`, invisible to Tapenade) | DD-ADI implicit LHS (`saGammaReThetaSolve`) |
| `evalSrcJacBlock` (independent 3×3 `A = +∂S/∂Q`) | `saGammaRetheta.F90:2191–2458` | `computeSrcLambda` → Eq. 59 source-dt restriction (ANK) |
| Blockette SA-GR residual kernels | `src/NKSolver/blockette.F90:6903–7559` | Fast residual path — **dead for SA-GR** (see F7) |

Crucial structural point: **the adjoint Jacobian contains no hand-written
math.** All `qq`/`evalSrcJacBlock` code is either behind `#ifndef
USE_TAPENADE` or never called from a differentiated head, and appears in the
AD output files only as undifferentiated primal copies. The hand-written
derivatives affect *solver convergence* (DADI LHS, ANK dt restriction), never
the converged solution or the adjoint gradients. Their errors therefore
manifest as robustness/convergence problems, not wrong gradients. They are
still fully re-derived in Appendix A per the audit brief.

---

## 2. Hand-written partials: verdicts (derivations in Appendix A)

Convention: `qq = −∂S/∂Q` (Source), `A = +∂S/∂Q` (evalSrcJacBlock); they are
sign-mirrors and share the same math. Term-by-term verdicts:

| Entry | Verdict | Notes |
|---|---|---|
| `qq(1,1)` = −∂S_ν/∂ν̃ | **Matches derivation**, with the *same deliberate omission as baseline `sa.F90:325`*: the `term1` production contribution (`term1 + ν̃·∂term1/∂ν̃`) is left out of the diagonal. All retained pieces (dfv1, dfv2, dft2, drr, dgg, dfw chains) re-derive exactly; γ correctly multiplies only the production pieces inside the bracket. Clipped ≥ 0 (LHS-only, matches SA philosophy). | See A.1 |
| `qq(2,2)` = −∂S_γ/∂γ | **Matches derivation exactly** (production `(1.5·ce1·γ−0.5)/√γ` factor and destruction `(2·ce2·γ−1)` factor both correct). Clamp saturation of `gammaLocal` ignored (interior derivative used even when clamped) — LHS-only. Clipped ≥ 0. | A.2 |
| `qq(3,3)` = −∂S_θ/∂R̃ | **Leading term correct** (`cθt/T·(1−FθT)`); deliberately omits ∂FθT/∂R̃ (via θBL→δ) and ∂Reθt_target/∂R̃ (via λθ ∝ R̃²). Crossflow sigmoid derivative `crossflowPhiPrime` re-derives exactly, but the reScf(θBL(R̃)) log-dependence is omitted. All omissions LHS-only. | A.3 |
| `qq(1,2)` = −∂S_ν/∂γ | **Matches derivation** in the normal branch. **Finding F2:** when `approxSA .and. transitionUseApproxSA` the residual sets `term1 = 0` but the Jacobian keeps the `cb1(1−ft2)·ss·ν̃` piece — Jacobian of a different residual in that mode. Clamp saturation of `gammaForSA` outside [0,1] also ignored. | A.4 |
| `qq(2,1)` = −∂S_γ/∂ν̃ | **Matches derivation exactly**: `drTurb/dν̃ = (fv1+χ·dfv1)/ν`, `∂f_onset1/∂r_T = r_T/f_onset1`, `∂F_onset/∂f_onset1 = 12·F_onset(1−F_onset)` (≡ 3·sech²), and the full two-branch chain for `dfTurb/dν̃` are all correct, including signs. | A.5 |
| `qq(2,3)` = −∂S_γ/∂R̃ | dFlength/dR̃ (Eqs. 49–50) and dReθc/dR̃ (Eq. 51) and the F_onset chain **all re-derive exactly**. **Finding F3:** the E_γ dependence on R̃ (through F_turb = (1−F_onset)e^{−r_T}, and F_onset depends on Reθc(R̃)) is silently omitted — only ∂P_γ/∂R̃ is kept, though the comment claims a full "analytical derivative". LHS-only. | A.6 |
| `qq(1,3) = 0` | **Exactly correct** — the SA source has no R̃ dependence. | — |
| `qq(3,1) = 0` | Correct **only when `transitionCrossflow = .false.`**. **Finding F5:** with crossflow on (the *default*: `inputParam.F90:303`), D_scf depends on ν̃ through `crossflowRatio = smoothMinMax(rTurb, 0.4)`, so the true ∂S_θ/∂ν̃ ≠ 0. The "paper §7.1" justification predates the crossflow term. LHS-only. | A.7 |
| `qq(3,2) = 0` | **Exactly correct for this implementation**: its FθT = F_wake·e^{−(y/δ)⁴} carries no γ term (the original LM2015 γ-dependent blanketing branch is not implemented here), so S_θ genuinely has no γ dependence. | — |
| `A(2,3)` (evalSrcJacBlock, one-sided FD of P_γ) | FD consistent with the analytic `qq(2,3)`; same E_γ omission; documented in-code as unused while A(3,2)=0. | A.6 |
| `evalSrcJacBlock` vs residual | **Finding F6:** uses hard `min(vortMag, vortLim)` where the residual uses `smoothMinMax(…, rsaGRpmin)` (`saGammaRetheta.F90:2374` vs `:514`), and has no `approxSA` handling. So the Eq. 59 eigenvalue estimate linearizes a slightly different source than the one actually solved. LHS/dt-restriction only. | — |

---

## 3. Drift check: do the AD files match the current primal?

Method: (a) git history of each primal vs the last AD regeneration; (b) content
fingerprints — distinctive constants and recent-change markers grepped in all
three AD modes; (c) inspection of the uncommitted AD diffs.

**The uncommitted AD modifications are regeneration, not hand edits.** The
working-tree diff of `saGammaRetheta_b.f90` is pure Tapenade idiom
(`pushcontrol1b`, branch flattening of `max`, lowercase) and introduces exactly
the two recent primal changes: the `refLenTrans` branch in the vorticity
limiter (`transitionreflength` appears 17× per file) and the generic
farfield-inflow ghost form in `turbBCRoutines_b.f90` (the old SA-GR
`2*wInf−interior` special case removed), matching primal commits 272ff47c and
701668a6. Fingerprints present and equal across `_d`/`_b`/`_fast_b`: onset
constants (2.6, 1.35), time scale 500.0, δ_BL factor 7.5, F_wake 1.0e6, the
full crossflow set (35.088, 319.51, 6200, 50000, 0.1066, 0.0125), diffusion
coefficients `sigmaF`/`sigmaTheta` (48 hits), the current `smoothMinMax`
(`lambda_switch` proximity-switch form, 10 hits in `turbUtils_d/_b`), and the
correlation derivatives (`rethetatcorrelation_d`, `flengthcorrelation_d`,
`rethetaccorrelation_d`). Post-regeneration primal commits (dac2c78f qq clips,
21f29567) touch only `#ifndef USE_TAPENADE` code, docs, and
`Makefile_tapenade` — no residual-path drift.

**One known, still-open drift — F1: `vortlimd = 0.0_8`.** The working-tree
`saGammaRetheta_d.f90:894` (and `_b`) still hard-zeroes the vorticity-limiter
derivative: `uInf`/`muInf` were passive when these files were generated.
Commit 21f29567 already fixed `Makefile_tapenade` (declares them active in the
`Source` fullRoutines head) but **Tapenade has not been rerun since** —
`uinfd`/`muinfd` have 0 hits in the generated file. Consequence: forward and
reverse partials w.r.t. flow-condition design variables (Mach/Re-type, through
`uInf = M√γ·√(p∞/ρ∞)` and `muInf`) silently drop the limiter path in every
cell where the cap engages. **dR/dw is unaffected** (uInf is constant w.r.t.
w), so the state adjoint, `fast_b`, and preconditioners are fine. STATUS:
TAPENADE NEEDED (user action), then `make`, then commit the six AD files
together with the Makefile change.

**F7: blockette SA-GR kernels are stale — but dead.** `blockette.F90`'s SA-GR
source kernel (last touched 2026-05-27, ed5f2feb) predates both the
`refLenTrans` limiter change (`blockette.F90:7081` still has the old
`uInf*sqrt(uInf/muInf)/20` without the reference length) and the entire
crossflow D_scf term (0 hits of 35.088/0.1066), while `transitionCrossflow`
defaults to `.true.`. This would be a live "solver converges residual A,
adjoint linearizes residual B" bug, except `adflow/pyADflow.py:6668` force-sets
`useBlockettes = False` for `SA-noft2-Gamma-Retheta`, routing residuals through
`saGammaRetheta_block(.true.)` (the true primal, `blockette.F90:821`). Verdict:
dead code today; a latent trap if blockettes are ever re-enabled for this
model. Recommend either updating or deleting the blockette SA-GR kernels.

**F8 (conditional): advection order in the adjoint.** The primal drivers
(`saGammaReTheta_block:101–108`, and its AD-embedded copies) force
`orderTurb = firstOrder` around `turbAdvection` when
`transitionFirstOrderUpwind` (default `.true.`) — note this forces first order
for **all three** turb equations, ν̃ included. The master AD sweeps
(`masterRoutines.F90:205/550/807/1166/1385`) call `turbAdvection*` directly
with **no** order switch, respecting the global `orderTurb`. Today this is
consistent because `orderTurb` defaults to `firstOrder`
(`inputParamRoutines.F90:4235`). If a user ever sets `orderTurb = 'second
order'`, the solver would converge a first-order transition-advection residual
while the adjoint linearizes a second-order one. Cheap hardening: replicate
the switch in the master sweeps, or reject `secondOrder`+SA-GR at option
validation.

---

## 4. Adjoint wiring trace

**State registration and sizes.** `nw = 8`, `nt2 = 8` for the model
(`inputParamRoutines.F90:2152–2154`), with `itu1/itu2/itu3 = 6/7/8`
(`constants.F90:40–43`). All adjoint sizing is generic off `nw`/`nt1:nt2`:
`setupStateResidualMatrix` (`adjointUtils.F90:87–106`) sets
`nState = nw = 8` (full), `nwf = 5` (frozen turb), `nt2−nt1+1 = 3` (turb-only
PC), and allocates/inserts PETSc values in `blk(nState,nState)` dense blocks
(`adjointUtils.F90:109,620–632`) — so the PETSc block size, matrix
preallocation, and `whalo2*` exchanges (`1..nw`, `masterRoutines.F90:159/504/859/1205`)
all pick up 8 automatically with no hardcoded offsets found. Residual monitors
register the three L2 norms (`inputParamRoutines.F90:93–99`). `turbResScale`
is a per-equation 3-vector (default `[1e4, 10, 1e4]`, `pyADflow.py:6673–6675`)
indexed `l−nt1+1` in NKSolvers.

**All five AD sweeps dispatch the model.** `master` (primal replay, allocates
`qqGR(2:il,2:jl,2:kl,3,3)`), `master_d`, `master_b` (exact reverse order:
ResScale_b → Viscous_b → turbAdvection_b → Source_b), `master_fast_b`, and
`block_res_state_d` all carry `case (spalartallmarasnoft2gammaretheta)` blocks
calling the generated SA-GR routines with `turbAdvection*(3,3,itu1−1,qqGR)`
(`masterRoutines.F90:205,550,807,1166,1385`).

**Off-diagonal coupling — both directions, verified in the generated code.**
This was the audit's priority-4 item; all blocks are present because the γ→SA
coupling lives *inside* `saGammaRetheta.F90`'s `Source` (rule 2 of CLAUDE.md)
and Tapenade differentiates that routine w.r.t. the full `w`:

- **d(SA residual)/dγ** (the classically-missed block):
  - Reverse: `saGammaRetheta_b.f90:1169–1191` — the SA-residual seed
    `scratchd(idvt)` back-propagates `term1d/term2d → gammaforsad → x1d →
    wd(i,j,k,itu2)` (clamp branch handled by `popcontrol1b`, so the derivative
    is exactly zero when γ is outside [0,1] — *more* correct than the
    hand-written `qq(1,2)`).
  - Forward: `saGammaRetheta_d.f90:722–751` — `wd(itu2) → x1d → gammaforsad →
    term1d, term2d → scratchd(idvt)`.
  - Reverse-fast: `saGammaRetheta_fast_b.f90:1149,1230` — `wd(itu2) += x1d/x2d`.
- **d(S_γ)/dν̃**: present via the `fv1d → rturbd → fonset1d/fturbd` chains in
  all three modes (r_T = ν̃·fv1/ν is in the residual text).
- **d(S_γ)/dR̃**: present via `rethetac_vald = rethetaccorrelation_d(...)`
  (`saGammaRetheta_d.f90:904`) and `flengthcorrelation_d`; `rethetatilded`
  seeds from `wd(itu3)` (`:842`).
- **d(S_θ)/d(mean flow)**: velocity/density/rlv dependencies (timeScale, λθ,
  FθT, reS) are differentiated — e.g. `wd(itu3)` viscous coupling rows
  `c30` at `_d:2153/2362` and the source chains through `velmag`, `fact`.
- **Viscous cross-coupling**: γ and R̃ diffusion coefficients depend on ν̃
  through `nu_t = ν̃·fv1`; the AD carries this (fv1_m/fv1_p chains), so
  d(γ-viscous)/dν̃ is not dropped.

Because `master`/`master_b` are also what the FD and AD preconditioner
assemblies linearize (`useAD` per `NK_ADPC`/`ANK_ADPC`/`ADPC`), the coupled-ANK
Jacobian and the adjoint PC inherit the same complete coupling.

**BC and freestream derivative paths.** Regenerated `turbBCRoutines_{d,b,fast_b}`
contain the current wall treatment (ν̃ antisymmetric; γ, R̃ Neumann) and the
generic farfield ghost form. `initializeFlow_b/d` contain the SA-GR
`referenceState` case (`initializeFlow_b.f90:141,439`): `winfd(itu2)=0` (γ∞=1
constant) and `winfd(itu3)=0` (correlation of passive Tu∞) — consistent with
the primal, noting Tu∞ is not a design variable.

---

## 5. Findings summary (ranked)

| # | Severity | Finding | Affects |
|---|---|---|---|
| F1 | **High (open, known)** → RESOLVED | `vortlimd = 0.0_8`: current AD files predate the `uInf/muInf` active declaration (Makefile fixed in 21f29567; **Tapenade rerun pending**) | dR/d(Mach, Re-type DVs) wherever the vorticity cap engages. dR/dw unaffected |
| F7 | Medium (latent) | Blockette SA-GR kernels drifted (no refLenTrans, no crossflow D_scf) — currently dead via `pyADflow.py:6668` force-off | Nothing today; residual/adjoint inconsistency if blockettes re-enabled |
| F8 | Medium (conditional) | Master AD sweeps lack the `transitionFirstOrderUpwind` orderTurb switch; safe only because `orderTurb` defaults to first order | Adjoint consistency iff user sets `orderTurb='second order'` |
| F2 | Low (LHS-only) | `qq(1,2)` keeps the `cb1(1−ft2)ss·ν̃` term when `approxSA∧transitionUseApproxSA` zeroes it in the residual | DADI convergence in approx-SA mode only |
| F3 | Low (LHS-only) | `qq(2,3)` omits E_γ's R̃-dependence (via F_onset in F_turb) despite "analytical derivative" comment | DADI convergence |
| F5 | Low (LHS-only) | `qq(3,1)=0` no longer exact with `transitionCrossflow=.true.` (default) — D_scf depends on ν̃ via crossflowRatio | DADI convergence |
| F6 | Low (dt-restriction) | `evalSrcJacBlock` uses hard `min` for vortMagLim vs `smoothMinMax` in the residual; no approxSA handling; same omissions as qq | Eq. 59 eigenvalue estimate / ANK dt restriction |
| — | Info | `qq(1,1)` omits the term1 production diagonal — identical to baseline `sa.F90:325` convention (implicit treatment of diagonal-dominance-adding terms only), then clipped ≥ 0 | By design |

**Conservative bottom line.** The adjoint operator itself (dR/dw in all three
AD modes) is Tapenade-generated from the *current* primal residual text, with
all γ/R̃/ν̃/mean-flow cross-coupling blocks structurally present in both
directions; I found no missing off-diagonal coupling in the differentiated
path. The one confirmed defect in the derivative code that reaches gradients
is F1 (extra-variable partials, fix staged, regeneration pending). All errors
found in hand-written derivative math live strictly in solver LHS / dt
restriction machinery. This audit cannot certify numerical correctness —
run the FD-vs-AD partials tests (dw-dot products per mode, then
`test_jacVecProdFWD`/`test_jacVecProdBWDFast`/`test_adjoint` analogues) before
trusting gradients, with special attention to a Mach-perturbation case where
the vorticity cap is active (targets F1 after regeneration).

---

## Appendix A — Hand derivations (check these, not just the conclusions)

Notation: ν̃ = `w(itu1)`, γ = `w(itu2)`, R̃ = `max(w(itu3), rsaGRreThetaLo)`;
γ_c = min(max(γ,0),1) (SA coupling clamp); γ_ℓ = min(max(γ, γ_lo), γ_hi)
(transition clamp); ν = rlv/ρ; χ = ν̃/ν; d = wall distance, D = 1/d²;
K = 1/κ²; Ω̃ = vortMagLim; c_* are `rsa*`/`rsaGR*` constants. All derivatives
of the clamps are taken in their interior (the code's convention; the AD
handles saturation exactly via branches).

### A.0 SA auxiliary derivatives (shared)

fv1 = χ³/(χ³+c_v1³):
∂fv1/∂χ = [3χ²(χ³+c_v1³) − χ³·3χ²]/(χ³+c_v1³)² = **3χ²c_v1³/(χ³+c_v1³)²**
→ code `dfv1` (w.r.t. χ; the 1/ν of ∂χ/∂ν̃ is applied downstream). ✓

fv2 = 1 − χ/(1+χ·fv1), with ∂χ/∂ν̃ = 1/ν:
∂fv2/∂ν̃ = −[(1/ν)(1+χfv1) − χ((1/ν)fv1 + χ·dfv1·(1/ν))]/(1+χfv1)²
= −(1/ν)[1 + χfv1 − χfv1 − χ²dfv1]/(1+χfv1)²
= **(χ²·dfv1 − 1)/(ν(1+χfv1)²)** → code `dfv2`. ✓

ft2 = c_t3·e^{−c_t4χ²}: ∂ft2/∂ν̃ = ft2·(−c_t4·2χ)(1/ν) = **−2c_t4·χ·ft2/ν**
→ code `dft2`. ✓

S̃ (code `sst`) = ss + ν̃·fv2·K·D, so ∂S̃/∂ν̃ = (fv2 + ν̃·dfv2)·K·D.
r = ν̃KD/S̃ (capped at 10):
∂r/∂ν̃ = KD/S̃ − (ν̃KD/S̃²)(fv2+ν̃dfv2)KD = **(KD/S̃)[1 − r(fv2+ν̃·dfv2)]**
→ code `drr`. ✓ (zero when the r=10 cap is active — code applies the cap
before, and the AD handles the branch; the hand code applies drr regardless of
the cap: harmless because at r=10 the fw curve is flat, and LHS-only anyway.)

g = r + c_w2(r⁶−r): ∂g/∂r = **1 − c_w2 + 6c_w2·r⁵** → code `dgg`. ✓

fw = g·T, T = ((1+c_w3⁶)/(g⁶+c_w3⁶))^{1/6}. ln T = (1/6)[ln(1+c_w3⁶) −
ln(g⁶+c_w3⁶)] ⇒ T′/T = −g⁵/(g⁶+c_w3⁶) ⇒
∂fw/∂g = T(1 − g⁶/(g⁶+c_w3⁶)) = **T·c_w3⁶/(g⁶+c_w3⁶)**
→ code `dfw = (cw36/(gg6+cw36))·termFw·dgg`. ✓

### A.1 qq(1,1) = −∂S_ν/∂ν̃

Residual: S_ν = (term1 + term2·ν̃)·ν̃ with
term1 = γ_c·c_b1(1−ft2)·ss,
term2 = γ_c·P₂ + Dst, P₂ = D·K·c_b1[(1−ft2)fv2 + ft2], Dst = −D·c_w1·fw.

Full derivative:
∂S_ν/∂ν̃ = term1 + ν̃·∂term1/∂ν̃ + 2·term2·ν̃ + ν̃²·∂term2/∂ν̃, where
∂term1/∂ν̃ = −γ_c·c_b1·ss·dft2,
∂term2/∂ν̃ = γ_c·D·K·c_b1[(1−ft2)dfv2 + (1−fv2)dft2] − D·c_w1·dfw.
Expand the bracket: (1−ft2)dfv2 + (1−fv2)dft2 = **dfv2 − ft2·dfv2 − fv2·dft2 +
dft2** — exactly the code's four-term bracket, with γ_c multiplying it and
−c_w1·dfw outside γ_c. ✓ (γ on production only: correct.)

Code stores qq(1,1) = −2·term2·ν̃ − D·ν̃²·[γ_c·c_b1·K(dfv2−ft2dfv2−fv2dft2+dft2)
− c_w1·dfw]: i.e. it **omits** `term1 + ν̃·∂term1/∂ν̃`. This is the identical
omission made by baseline `sa.F90:325–328` (production kept explicit; only
diagonal-dominance-adding pieces treated implicitly), followed by the same
`max(·, 0)` clip. Deliberate convention, LHS-only. ✓-by-convention.

### A.2 qq(2,2) = −∂S_γ/∂γ

S_γ = P_γ − E_γ,
P_γ = c_a1·F_len·F_on·Ω̃·√γ_ℓ·(1−c_e1γ_ℓ) = C_P(γ_ℓ^{1/2} − c_e1γ_ℓ^{3/2}),
C_P = c_a1F_lenF_onΩ̃;
E_γ = c_a2·F_turb·Ω̃·(c_e2γ_ℓ² − γ_ℓ).
F_on, F_turb, Ω̃, F_len are γ-independent (F_on depends on reS, Reθc(R̃), r_T;
F_turb = (1−F_on)e^{−r_T}).

∂P_γ/∂γ = C_P(½γ^{−1/2} − 3/2·c_e1γ^{1/2}) = C_P·(0.5 − 1.5c_e1γ)/√γ.
∂E_γ/∂γ = c_a2F_turbΩ̃(2c_e2γ − 1).
−∂S_γ/∂γ = C_P·**(1.5c_e1γ − 0.5)/√γ** + c_a2F_turbΩ̃·**(2c_e2γ − 1)**
= code exactly (`saGammaRetheta.F90:697–702`). ✓ Then clipped ≥ 0 (triggers in
the pre-transition region where both factors are negative) — LHS-only.
`evalSrcJacBlock` A(2,2) is the exact negation. ✓

### A.3 qq(3,3) = −∂S_θ/∂R̃

S_θ = (c_θt/T)(Reθt_tgt − R̃)(1 − FθT) [+ D_scf],
T = 500ν/|u|² (no R̃), FθT = F_wake·e^{−(y/δ)⁴},
δ = 50·y·Ω·δ_BL/|u|, δ_BL = 7.5θ_BL, θ_BL = R̃ν/|u| ⇒ δ ∝ R̃,
λθ = (θ_BL²/ν)·dU/ds ∝ R̃² ⇒ Reθt_tgt = f(Tu, λθ(R̃)).

Full derivative:
∂S_θ/∂R̃ = (c_θt/T)[ (∂Reθt_tgt/∂R̃ − 1)(1−FθT) − (Reθt_tgt−R̃)·∂FθT/∂R̃ ],
with ∂FθT/∂R̃ = FθT·4(y/δ)⁴·(1/R̃) ≥ 0 (since dδ/dR̃ = δ/R̃, and F_wake has
no R̃).

Code keeps only −(∂/∂R̃)[−(c_θt/T)R̃(1−FθT)]|_{FθT,tgt frozen} =
(c_θt/T)(1−FθT) ≥ 0 — the pure relaxation diagonal. Omitted: the
Reθt_tgt(λθ(R̃)) and FθT(R̃) chains (and, with crossflow, reScf’s
−35.088·log(k_rough/θ_BL) dependence on R̃). All omissions LHS-only; the kept
piece is the dominant, always-stabilizing one. **Approximate by design; not
equal to the true partial.**

Crossflow addition: D_scf = (c_θt/T)·c_cf·Φ_p(reScf−R̃, 0; p_min)·FθT with
smooth-min Φ. For p<0 and x = reScf−R̃ < 0 (so x is the min): Φ = x +
ln(1+e^{p(0−x)})/p ⇒ dΦ/dx = 1 − e^{−px}/(1+e^{−px}) = **1/(1+e^{p(0−x)})** =
code's `crossflowPhiPrime` (guarded to reScf<R̃). ∂D_scf/∂R̃ ≈ −(c_θt/T)c_cf·
FθT·Φ′ (freezing reScf, FθT), so qq(3,3) += +(c_θt/T)c_cf·FθT·Φ′. ✓ matches.

### A.4 qq(1,2) = −∂S_ν/∂γ

Only term1 and the γ_c·P₂ piece of term2 depend on γ:
∂S_ν/∂γ_c = [c_b1(1−ft2)ss]·ν̃ + P₂·ν̃² = **(c_b1(1−ft2)ss + P₂ν̃)·ν̃**,
qq(1,2) = −that. Matches `saGammaRetheta.F90:746–747` (P₂ = `term2_prod`). ✓

Two caveats: (i) ∂γ_c/∂γ = 0 when γ∉[0,1]; code always uses 1 (the AD gets
this right via branch controls). (ii) **F2**: in the
`approxSA∧transitionUseApproxSA` branch the residual has term1 = 0, but this
derivative unconditionally includes the term1 piece — the Jacobian then
corresponds to a residual that isn't the one being solved. LHS-only.

### A.5 qq(2,1) = −∂S_γ/∂ν̃

Only r_T = ν̃fv1/ν depends on ν̃ (reS, Ω̃, γ, F_len do not).
∂r_T/∂ν̃ = (fv1 + ν̃·∂fv1/∂ν̃)/ν = (fv1 + ν̃·dfv1/ν)/ν = **(fv1 + χ·dfv1)/ν**
→ `drTurb_dnu`. ✓

f₁ ≡ F_onset1 = √(R² + r_T²), R = reS/(2.6Reθc):
∂f₁/∂r_T = **r_T/f₁** → `dfOnset1_drT`. ✓

F_on = ½(tanh(6(f₁−1.35)) + 1): with t = tanh(·) = 2F_on−1,
∂F_on/∂f₁ = 3(1−t²) = 3(1−(2F_on−1)²) = **12·F_on(1−F_on)** →
`dfOnset_dfOnset1`. ✓ (the `sech²` form used in A.6 is the same quantity.)

F_turb = (1−F_on)e^{−r_T}:
∂F_turb/∂ν̃ = −e^{−r_T}·∂F_on/∂ν̃ + (1−F_on)e^{−r_T}(−∂r_T/∂ν̃)
= **−e^{−r_T}·dfOnset_dnu − F_turb·drTurb_dnu** → `dfTurb_dnu`. ✓

−∂S_γ/∂ν̃ = −∂P_γ/∂ν̃ + ∂E_γ/∂ν̃
= **−c_a1F_len·dfOnset_dnu·Ω̃√γ_ℓ(1−c_e1γ_ℓ) + c_a2·dfTurb_dnu·Ω̃γ_ℓ(c_e2γ_ℓ−1)**
→ matches `saGammaRetheta.F90:769–774` term for term and sign for sign;
`evalSrcJacBlock` A(2,1) is the exact negation. ✓

### A.6 qq(2,3) = −∂S_γ/∂R̃

Chains through F_len(R̃) and Reθc(R̃) (both correlations of R̃).

dF_len/dR̃ (Eqs. 49–50): F_len = 44 − N/B^{1/6}, N = 44 − (0.5 −
3·10⁻⁴(R̃−596)) = 43.5 + 3·10⁻⁴(R̃−596), B = 1 + F₁, F₁ = e^{−0.03(R̃−460)},
F₁′ = −0.03F₁:
dF_len/dR̃ = −N′/B^{1/6} − N·(−1/6)B^{−7/6}·F₁′
= **−3·10⁻⁴/B^{1/6} − (N/6)·B^{−7/6}·F₁·0.03**
→ matches code (`inner_val` = N). ✓

dReθc/dR̃ (Eq. 51): Reθc = 0.67R̃ + 24sin(R̃/240 + 0.5) + 14 ⇒
**0.67 + (24/240)cos(R̃/240+0.5) = 0.67 + 0.1cos(·)** → code. ✓

∂f₁/∂R̃ = (R·∂R/∂R̃)/f₁ = R·(−R/Reθc)·Reθc′/f₁ =
**−(reS/(2.6Reθc))²·Reθc′/(Reθc·f₁)** → `dfOnset1_dReT`. ✓
∂F_on/∂R̃ = 3(1−tanh²(6(f₁−1.35)))·∂f₁/∂R̃ → `sech2_val` form. ✓

∂P_γ/∂R̃ = c_a1Ω̃√γ_ℓ(1−c_e1γ_ℓ)·[F_on·dF_len + F_len·dF_on] →
`pGamma_common·(fOnset·dFlength_dReT + fLength_val·dfOnset_dReT)`. ✓

**Omission (F3):** E_γ = c_a2(1−F_on)e^{−r_T}Ω̃γ_ℓ(c_e2γ_ℓ−1) also depends on
R̃ through F_on. True qq(2,3) = −∂P_γ/∂R̃ + ∂E_γ/∂R̃ with
∂E_γ/∂R̃ = −c_a2·e^{−r_T}·dfOnset_dReT·Ω̃γ_ℓ(c_e2γ_ℓ−1); the code drops this
term (kept: production only). LHS-only, but the in-code comment overstates
completeness. `evalSrcJacBlock` A(2,3) (one-sided FD of P_γ alone) makes the
same omission, consistently.

### A.7 qq(3,1) = 0 — validity

Without crossflow: T, θ_BL, δ, F_wake, reS, Reθt_tgt contain ρ, u, rlv, y, R̃
but no ν̃ ⇒ ∂S_θ/∂ν̃ = 0 **exactly**. With `transitionCrossflow = .true.`
(default): D_scf depends on ν̃ via crossflowRatio = Φ_p(r_T, 0.4; p_min),
r_T = ν̃fv1/ν, feeding ΔH± → reScf → the smooth-min production. So the true
∂S_θ/∂ν̃ ≠ 0 wherever D_scf is active and r_T is near/below 0.4. Code keeps
qq(3,1)=0 (and `evalSrcJacBlock` likewise, plus its A(2,3)-unused note).
LHS-only (F5). Kept zero after eigenvalue-impact assessment (see F5 in the
resolution addendum): with A13 ≡ 0, A31 perturbs the spectrum only through
A12·A23·A31, negligible for the conservative max(0, λ) dt restriction.
