# Lessons from `mdolab/adflow` `sst_dev` for the SA-γ-Re̅θt AD/adjoint work

**Date:** 2026-07-07
**Sources:** GitHub API compare `main...sst_dev` (31 commits, 82 files, David
Anderegg, Sep 2023 – Jan 2024, merged upstream via PR #331); local analysis
`adflow_sst_dev_derivative_analysis.md` (repo root); direct diff inspection of
`Makefile_tapenade`, `autoEditReverseFast.py`, `masterRoutines.F90`,
`NKSolvers.F90`, `blockette.F90`, `turbUtils.F90`, `turbAPI.F90`, `SST.F90`,
`flowUtils.F90`, `sa.F90`, `adjointUtils.F90`.
**Scope:** analysis only — nothing from `sst_dev` is implemented here.
Adjoint linearization on this branch remains frozen (CLAUDE.md rule 6).

**Baseline note.** This branch forked from a **pre-#331** main: no
`sa_block_residuals` wrappers, `autoEditReverseFast.py` still strips push/pop,
`d2wall` still `(2:il,2:jl,2:kl)`. Some multi-equation solver plumbing that
`sst_dev` introduced upstream (per-equation `turbResScale` indexing,
`nt1:nt2` physicality loops) exists here anyway — added independently by the
SA-GR work. Where the two branches solved the same problem differently, the
table in §5 says so.

---

## 1. Approach taken by `sst_dev` to differentiate a multi-equation model

The strategy was **"make the model look like SA to Tapenade, generalize the
plumbing to `nt1:nt2`, and hand-write only the call-sequence wrappers"**:

1. **Restructure the primal into Tapenade-head-sized routines.**
   `SST.F90` was rewritten from a monolithic `SSTSolve` into
   `SSTSource` / `SSTViscous` / `SSTResScale` / `f1SST`, mirroring the
   `saSource`/`saViscous`/`saResScale` decomposition that already
   differentiated cleanly. Each becomes a head in `Makefile_tapenade`
   (`fullRoutines` and `stateOnlyRoutines` for `fast_b`).

2. **Module-level coupled Jacobian block.** A rank-5 allocatable
   `qq(2:il,2:jl,2:kl,2,2)` holds the 2×2 point-coupled DADI Jacobian,
   passed to the generic `turbAdvection(2,2,offset,qq)` /
   `unsteadyTurbTerm(2,2,offset,qq)`. (SA-GR's `qq(...,3,3)` +
   `turbAdvection(3,3,...)` is the same device one size up.)

3. **Auxiliary fields moved from module arrays into `scratch` slots.**
   Cross-diffusion `kwCD`, blending `f1`, and production `prod` moved from
   `turbMod` pointers into `scratch(:,:,:,icd/if1SST/iprod)` so they are
   *active* through the differentiated heads instead of hidden module state.
   `kwCDterm`, `prodWmag2/Smag2/KatoLaunder` gained explicit index-range
   arguments and II-LOOP forms and became differentiated heads themselves.

4. **Hand-written call-sequence wrappers, generated bodies.** Each model
   module got primal `X_block_residuals(cleanUp)` plus hand-written
   `X_block_residuals_d/_b/_fast_b` that assemble the Tapenade-generated
   pieces in the right (reversed, for `_b`) order, marked
   `Manual Differentiation Warning`. `masterRoutines.F90` and
   `blockette.F90` then dispatch one call per model per sweep. Only the
   *orchestration* is manual; every derivative body is generated.

5. **Widened active-argument lists where the second equation made BCs and
   eddy viscosity state-dependent.** `bcTurbTreatment` and
   `applyAllTurbBCThisBlock` heads gained `bvt*` and `BCData` inlet
   quantities as active in/outputs; `computeEddyViscosity` gained
   `d2wall, vol, si, sj, sk, scratch, timeref` (SST's μ_t depends on
   vorticity/strain and F2, not just `w` and `rlv`).

6. **Wall distance in halos.** `d2wall` grew from `(2:il,...)` to
   `(0:ib,...)` with a new `exchanged2Wall` halo exchange, and the
   wall-distance recompute moved *before* the main loop in
   `master`/`master_d`/`master_b` and `blocketteRes` (mesh-derivative path
   must differentiate through it consistently).

7. **Solver plumbing generalized from "1 turb var" to `nt1:nt2`.**
   Physicality checks loop `l = nt1, nt2` (replacing a skip with a TODO);
   `turbResScale(l - nt1 + 1)` replaces `turbResScale(1)` in
   `setRVecANKTurb` and the ANK time-step block; `approxSA` renamed
   `approxTurb`, with the SST realization "disable the production term" in
   the low-order ANK preconditioner Jacobian.

## 2. What worked (stable / merged)

- **The whole thing eventually merged upstream** (PR #331, Jan 2024) — the
  approach is proven end-to-end for a 2-equation model.
- **Forward tangent mode** validated against complex-step with step sizes
  spanning 1e-40…1e-200 (machine-precision agreement — the differentiated
  model is primal-correct).
- **Reverse ↔ reverse-fast internal consistency** to machine precision via
  dot-product (transpose) tests, including a new `test_jacVecProdBWDFast.py`
  layer for the `fast_b` state-Jacobian routines.
- **The modular Source/Viscous/ResScale decomposition + generic
  `turbAdvection(mAdv,nAdv,offset,qq)`** — no fundamental Tapenade failure;
  the 2×2 coupled block differentiated as routinely as SA's scalar.
- **Per-equation solver plumbing** (`nt1:nt2` loops, indexed
  `turbResScale`) — small, mechanical, uncontroversial.

## 3. What failed or caused problems

- **The k-correction pressure coupling was the recurring villain.**
  `p += (5/3−γ)·ρ·k` (`computePressureSimple`) couples the turbulence state
  into pressure → fluxes → forces. Commits show it disabled for weeks
  ("backwards partials agree … *(k-correction disabled)*", "need to disable
  k-correction") before "activate kpresent again. forward ADJOINTS are
  working!". The term is smooth — the trouble was **bookkeeping** (aliasing,
  loop structure), not physics.
- **Primal aliasing bugs.** The forward vector product *changed the primal
  residuals*; fixed by making `computePressureSimple` write through a temp
  (`pp`) and re-splitting loops so the derivative sweep never mutates primal
  state.
- **Scratch-slot lifetime hazard.** `prod` lives in `scratch(:,:,:,iprod)`
  and is *overwritten by the eddy-viscosity computation*; the hand-written
  reverse wrappers must **recompute the primal production term** before the
  reverse sweep (explicit recompute calls at the top of
  `SST_block_residuals_b/_fast_b`).
- **Unassigned halo/BC values polluted derivatives.**
  `computeEddyViscosity` had to shrink its loop range on boundary faces
  ("saveguard against using values on BC's where they might not be
  assigned"); `kwCDterm` does the same dance.
- **`autoEditReverseFast.py` garbage-stripping broke SST's `fast_b` code.**
  The push/pop deletion patterns that were safe for SA produced non-working
  code for SST and were **commented out** ("does not work with SST"), plus a
  segfault fix and manual garbage-removal automation for `sst_fast_b`.
- **Symmetry-plane normal computation** needed a robustness patch (2/3-point
  vectors instead of corners) to survive differentiation on O-grids.
- **Verification was loosened, not completed.** SST FD-vs-AD tolerances
  blown out to rtol 2.11 (state) and 41 (mesh); adjoint regression covers
  only `cl`; CS coverage drops `span`. The state-derivative 211% is a real
  unresolved signal, most plausibly in the k-coupling path. The transpose
  tests prove *consistency*, not *correctness* — forward and reverse can be
  wrong identically.

## 4. Bearing on the coupling-block issue flagged in our docs

*(The request cut off mid-sentence at "flagged in …". The three flagged
items in this repo that read as "coupling-block" issues are covered below;
if a different one was meant, say which and this section gets refocused.)*

**(a) The turbulence→mean-flow coupling block (the SST killer) is
structurally absent here.** `kPresent` stays `.false.` for
`spalartallmarasnoft2gammaretheta` (`inputParamRoutines.F90:2151-2156`), so
`getCorrectForK()` is false and the pressure never sees a turbulence
variable. γ and Re̅θt reach the mean flow **only** through
μ_t = ρν̃·fv1 — γ does not modulate μ_t (A3 §9) — i.e. exactly the coupling
column SA already has and that ADflow's AD already handles. The single
hardest problem `sst_dev` faced has no analog in SA-GR. Corollary: when our
AD is eventually unfrozen, there is no reason to expect the SST-style
months-long k-coupling fight; the risk profile is SA-like plus source-term
kinks.

**(b) The intra-turbulence coupling blocks (qq off-diagonals, W3/A3) live
on the LHS only — same as SST.** In both branches the coupled block
`qq(m,n)` is the point-implicit DADI Jacobian, not part of the residual;
`sst_dev` never differentiated it (all `qq` logic sits behind
`#ifndef USE_TAPENADE`, as ours does), and the *true* coupling
∂R_ν̃/∂γ etc. is captured automatically because `Source` is differentiated
as one routine with all transport variables active in `w`. Our recent
`qq(1,1)/qq(2,2)` diagonal clipping (commit dac2c78f) therefore has zero AD
footprint — consistent with how SST treated its 2×2 block.

**(c) The coupled-ANK path (D-A2-1, Eq. 59 absent in the 8-eq system).**
`sst_dev` is mildly reassuring here: SST's coupled ANK also has **no
source-eigenvalue time-step restriction** — its only concession to source
stiffness in the coupled Jacobian/PC is `approxTurb` (drop the production
term in the low-order PC) plus the physicality check, and that was enough to
merge. Our branch already has both devices (`transitionUseApproxSA`,
`saGammaRetheta.F90:446`; `physicalityCheckANK{,Turb}` with γ/Re̅θt bounds).
The decision to ignore CANK (D-A2-1 response) is consistent with what
`sst_dev` shipped.

## 5. Comparison table — `sst_dev` vs current SA-γ-Re̅θt

| Problem | `sst_dev` (SST, k+ω) | This branch (SA-GR, ν̃+γ+Re̅θt) | Assessment |
|---|---|---|---|
| Primal structure for Tapenade | Rewrote monolith into `SSTSource/SSTViscous/SSTResScale/f1SST` heads | Built modular from day one (`Source/Viscous/ResScale`), copied from `sa.F90` | Same pattern; `sst_dev` validates it end-to-end |
| Coupled point Jacobian | Module allocatable `qq(...,2,2)`, generic `turbAdvection(2,2,offset,qq)`, behind `#ifndef USE_TAPENADE` | `qq(...,3,3)`, `turbAdvection(3,3,offset,qq)`, same guard | Identical device; LHS-only, no AD footprint |
| AD call-sequence wiring | Hand-written `X_block_residuals_d/_b/_fast_b` wrappers inside the model module ("Manual Differentiation Warning") | Inline `select case` blocks in `masterRoutines.F90` (5 sweeps × 1 case each) | Equivalent function. Wrapper style centralizes ordering in one place per mode (5 inline copies can drift); upstream main now uses wrappers — a merge conflict point later |
| Turb → mean-flow coupling | k enters pressure (`getCorrectForK`), → fluxes → forces; **the** hard problem (aliasing, weeks disabled) | `kPresent=.false.`; only μ_t=ρν̃fv1, γ not in μ_t — SA-identical coupling | **Absent by construction** here; biggest de-risking finding |
| Auxiliary model fields | Moved `kwCD/f1/prod` from `turbMod` arrays into active `scratch` slots; reverse wrappers must **recompute** `prod` (eddy-viscosity overwrites it) | Correlations are functions (`turbUtils.F90:2279-2412`) with `_d` versions; SA-GR uses `scratch` for `dvt` slots | When AD unfreezes: audit every SA-GR `scratch` slot for overwrite between the forward sweep and the reverse sweep — this bit SST twice |
| State-dependent turbulence BCs | Had to widen `bcTurbTreatment`/`applyAllTurbBC` active args (`bvt*`, BCData inlet quantities) | γ/Re̅θt BCs: freestream values from `wInf` (Re̅θt from Tu∞ — input param, not state) via standard `bmt/bvt`; SA-identical activity | No new active BC inputs expected; the widened upstream signature comes free on a future main merge |
| Wall distance in derivatives | `d2wall` extended to halos `(0:ib,...)` + `exchanged2Wall`; recompute moved before main loop | `d2wall(2:il,...)` (pre-#331); SA-GR `Source` uses d2wall at owned cells only, like SA | Not needed for SA-GR per se; becomes relevant only on merging post-#331 main |
| `fast_b` post-processing | `autoEditReverseFast.py` push/pop stripping **broke SST**; patterns disabled upstream | Stripping still **active** here (`autoEditReverseFast.py:23-28`) and SA-GR `fast_b` is built with it | **Watch item.** Same failure class (multi-eq `fast_b` with stripped push/pop) could bite `saGammaRetheta_fast_b`. When AD is regenerated/validated, run the BWDFast dot-product test first; if it fails, suspect the stripping before the model |
| Per-equation solver plumbing | Introduced `nt1:nt2` physicality loops, `turbResScale(l-nt1+1)`, `stateToCons(nt1:nt2)` | Present (`NKSolvers.F90:1303,2178,2794,2988`; 8 `nt1,nt2` loops) — added independently on this branch | Convergent evolution; nothing to port |
| Approx Jacobian for ANK PC | `approxSA`→`approxTurb`; drop production in low-order PC | `transitionUseApproxSA` (default true) — same idea, SA-GR scoped | Equivalent; keep straight that this is a PC approximation, never the adjoint Jacobian |
| Source-term stiffness in coupled path | Nothing beyond approxTurb + physicality check; merged anyway | Same, plus Eq. 59 in DADI/turbKSP (paper device `sst_dev` lacked) | Our coverage is a superset; D-A2-1 "ignore CANK" decision consistent with what shipped |
| Non-smooth model switches | `max`/`min` in μ_t and production limiter — inherent kinks, AD picks active branch | Fonset/Fturb/damping use `smoothMinMax` (differentiated `smoothminmax_d`) | SA-GR is deliberately *smoother* than SST at the same problem points |
| Verification ladder | CS → FD-vs-forward → dot-product → fast dot-product; **tolerances inflated to non-verification** (rtol 2.11 / 41), adjoint regression `cl` only | No AD validation yet (linearization frozen) | Adopt the 4-layer ladder *and* the post-mortem: use norm-relative or masked metrics for mesh derivatives instead of inflating global rtol; validate cd/cm and all DVs with CS, not just lift |

## 6. One-paragraph takeaway

`sst_dev` proves the exact recipe this branch already follows (modular heads
+ generic multi-equation `qq`/advection + guarded LHS Jacobian) carries a
coupled multi-equation turbulence model through Tapenade to a merged,
self-consistent adjoint. Every serious fight it had was **bookkeeping around
the turbulence→mean-flow coupling** (pressure k-correction, aliasing,
scratch lifetimes, halo garbage) — and SA-GR's γ-multiplies-production-only
design means that coupling block does not exist here. The two concrete
watch items to carry forward to the AD-unfreeze task: (1) the
`autoEditReverseFast.py` push/pop stripping that broke SST's `fast_b` is
still active on this branch; (2) audit SA-GR `scratch` slot lifetimes across
the forward/reverse sweep boundary. And when validating: transpose tests
prove consistency, not correctness — budget for CS coverage of drag/moment
and never "fix" a failing FD comparison by inflating its tolerance.
