# SA-γ-Re̅θt adjoint implementation audit (pre-partials-test)

**Date:** 2026-07-07. Visual check of the AD/adjoint wiring against the SST
(`sst_dev`, PR #331) playbook in `audits/sst_dev_lessons.md`. No tests run;
build compiled and imported (`make` → `Module libadflow was successfully
imported`). One fix applied (see §3); everything else verified as-is.

## 1. Verified correct (pass)

| Item | Evidence |
|---|---|
| 5 master-sweep case blocks (`master`, `master_d`, `master_b`, `master_state_b`, `block_res_state_d`) | `masterRoutines.F90:205,550,807,1166,1385`; `_b`/`_fast_b` in exact reverse order (ResScale→Viscous→Advection→Source), `turbAdvection*(3,3,itu1-1,qqGR)` everywhere; allocate/deallocate of `qqGR` only in primal `master` — same discipline as SA and as `sst_dev`'s wrappers |
| Makefile heads (fullRoutines + stateOnlyRoutines) | `Makefile_tapenade:187-193,314-320`; active sets mirror SA's (`w,rlv,vol,si,sj,sk,timeRef,d2wall > …,scratch`) |
| Generated modules export `qq` | rank-5 allocatable at line 63 of all three generated files; `masterRoutines` renames `qqGR => qq` |
| Eddy-viscosity dispatch differentiated | SA-GR shares the SA case *outside* `#ifndef USE_TAPENADE` (`turbUtils.F90` computeEddyViscosity); μ_t needs only `w`,`rlv` — head signature sufficient (unlike SST, which needed d2wall/metrics added) |
| No k-correction coupling | `kPresent=.false.` for the enum (`inputParamRoutines.F90:2152-2156`) → `getCorrectForK()` false → the `sst_dev` pressure-coupling problem class is structurally absent |
| Scratch slots | SA-GR uses `idvt..idvt+2` = slots 1–3 only; primal scratch has 10, derivative scratch 5 (`adjointUtils.F90:797`) — no overflow; no SST-style overwrite between forward and reverse sweeps (SA-GR never reads `iprod/icd/if1SST`, `saEddyViscosity` doesn't touch scratch) |
| Wall BC AD | generated `bcturbtreatment_d` block matches current primal exactly (ν̃ antisym `bmt=+1`; γ, Re̅θt Neumann `bmt=-1`, `bvt(itu3)=0`) |
| Farfield BC AD | regenerated files contain the generic ghost form (`bvt=winf`, `bmt=0`); old `2*wInf` SA-GR special case removed — consistent with primal commit 701668a6 |
| `wInf` seeding | `referencestate_d`: `winfd = 0` wholesale, then `winfd(itu1)` via `sanuknowneddyratio_d`, `winfd(itu2)=0` (γ∞=1), `winfd(itu3)` correctly stays 0 (correlation of passive Tu∞) |
| Correlations + smoothMinMax differentiated | `rethetatcorrelation_d`, `flengthcorrelation_d`, `rethetaccorrelation_d`, `smoothminmax_d` all present in `turbUtils_d.f90` (primal 2026-05-17 < AD 2026-05-18, current) |
| Crossflow (D_scf) in AD | 20 `dscf` hits in `saGammaRetheta_d.f90` (regenerated at e5cd58cd) |
| `transitionRefLength` in AD | present in all three regenerated files (the `docs/adjoint-trace.md` "pending" note was stale — corrected) |
| Primal-only routines excluded | `saGammaReThetaSolve`, `evalSrcJacBlock`, `computeSrcLambda` appear as undifferentiated primal copies only — correct (LHS/DADI machinery) |
| qq Jacobian guards | all `qq` logic (incl. the new dac2c78f diagonal clips) behind `#ifndef USE_TAPENADE` — zero AD footprint, same as SST's 2×2 block |
| Build fileList | `saGammaRetheta_{d,b,fast_b}` all listed (`src/build/fileList:148,169,182`) |
| Solver plumbing | per-equation `turbResScale(l-nt1+1)` (`NKSolvers.F90:1303,2178,2794,2988`), `nt1:nt2` physicality loops (8×) — equivalent to what `sst_dev` had to add |
| Python layer | no SA-GR-specific adjoint blockers in `pyADflow.py` |

## 2. Key state discovery

The six generated AD files (`saGammaRetheta_{d,b,fast_b}`,
`turbBCRoutines_{d,b,fast_b}`) are **modified but uncommitted** — Tapenade was
rerun after the last two primal commits. Content verified consistent with the
current primal. They compile and import. **Commit them together with the
`Makefile_tapenade` change once the next regeneration (§3) is done.**

## 3. Defect found and fixed (needs Tapenade rerun)

**`vortlimd = 0.0_8` — vorticity-limiter derivative hard-zeroed.**
`Source` uses `uInf` and `muInf` in the P&Z Eq. 52-53 cap
(`saGammaRetheta.F90:511`), but neither was declared active in the
`saGammaRetheta%Source` head, so Tapenade treated them passive (0 hits of
`uinfd`/`muinfd` in the generated file, vs. `timerefd` correctly propagated
because `timeRef` *is* in the head). Consequence: wherever the limiter is
active, d(residual)/d(Mach, Re-type DVs) silently loses the `uInf`/`muInf`
path — exactly the kind of localized extra-variable error the `sst_dev`
FD-vs-AD tests would flag on the mesh/DV axis.

- **Fix applied:** `Makefile_tapenade:187` — added `uInf, muInf` to the
  fullRoutines `Source` independents. State-only (`fast_b`) head untouched:
  uInf is constant w.r.t. `w`, so the state Jacobian was never affected.
- **STATUS: TAPENADE NEEDED** (then `make`, then commit generated files).
- No error in dR/dw or the adjoint preconditioner from this; it only affects
  forward/reverse partials w.r.t. flow-condition DVs (mach etc.).

## 4. Watch items (cannot be verified visually — first suspects if partials fail)

1. **`autoEditReverseFast.py` push/pop stripping** (active on this branch;
   upstream disabled it because it *broke SST's* `fast_b`).
   `saGammaRetheta_fast_b.f90` post-strip has 0 push/popreal8 and its 4
   pushcontrol pairs intact. If the BWDFast dot-product test fails, suspect
   the stripping before the model.
2. **`vortLim` and `smoothMinMax` kinks:** derivative is exact but one-sided
   at the cap/blend points — expect FD noise there, not AD error; use CS or
   masked metrics rather than loosening tolerances (SST post-mortem).
3. **`ResScale` volRef path** mirrors SA exactly (`(scratch,dw)>(dw)` head);
   whatever upstream SA does for mesh partials here is inherited unchanged —
   correct-by-regression, not re-derived.
