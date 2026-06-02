# OpenMP Parallelization of SA-γ-Reθ DADI Residual/Jacobian Assembly

## Context

When the turbulence equations are solved with DADI (the `ANKUseTurbDADI=True` path,
or pure DADI), ADflow processes one block at a time with **no OpenMP**, so on
high-core-count machines most cores sit idle during the turbulence solve. The
expensive part of each DADI iteration is assembling the residuals
(`scratch`) and the 3×3 Jacobian (`qq`) in `Source()` and `Viscous()`; the
tridiagonal DADI solve itself is inherently sequential.

**Goal:** Parallelize the residual + Jacobian assembly for the SA-γ-Reθ DADI
solver with OpenMP, with **zero data copying** and **no second copy of the
physics**, so the gains are not eaten by memory traffic on this memory-bound
kernel.

**Decision (Option B):** Add `!$OMP parallel do` directly onto the existing
primal loops in `saGammaRetheta.F90`. We do **not** route DADI through the
blockette path. Rationale — consistency:

- `Source()`/`Viscous()` already write the per-cell `scratch` and `qq` global
  per-block arrays that the DADI solve consumes. There is **nothing to copy**.
- The `qq` Jacobian logic stays in exactly **one** place. Routing DADI through
  blockette would force duplicating the Jacobian derivatives into
  `blockette.F90` (a third copy of the transition physics), which is the main
  consistency hazard we want to avoid.
- `!$OMP` lines are comments, and they go only on the **primal** loop branch
  (`#else`), never the Tapenade branch. The differentiated math is unchanged,
  so **no Tapenade regeneration is required** and the adjoint cannot desync.

**Settings chosen with user:**
- OpenMP is **always compiled in**; thread count is controlled by
  `OMP_NUM_THREADS` (1 thread ⇒ serial). No new Python option.
- Scope is **SA-γ-Reθ only** (`saGammaRetheta.F90`). `sa.F90` is left serial.

---

## Loop structure (already verified)

`Source()` and each directional sweep in `Viscous()` use this shared-body
pattern:

```fortran
#ifdef TAPENADE_REVERSE
        !$AD II-LOOP
        do ii = 0, nx*ny*nz - 1
            i = mod(ii, nx) + 2
            j = mod(ii/nx, ny) + 2
            k = ii/(nx*ny) + 2
#else
            do k = 2, kl
                do j = 2, jl
                    do i = 2, il
#endif
                        ... shared loop body: writes scratch(i,j,k,*) and qq(i,j,k,*,*) ...
#ifdef TAPENADE_REVERSE
                    end do
#else
                end do
            end do
        end do
#endif
```

The OpenMP directive goes **inside the `#else`**, immediately before
`do k = 2, kl`, with the matching `!$OMP end parallel do` immediately after the
three `end do`. The Tapenade branch is untouched.

**Race analysis (confirmed safe):** every loop iteration writes only its own
`scratch(i,j,k,*)`, `qq(i,j,k,*,*)`, and `transitionDebug(i,j,k,*)`; it reads
neighbor values only from `w` and the metrics (read-only, shared). In
`Viscous()` the three directional loops run sequentially and each cell touches
only its own `scratch(i,j,k)` (the `+ neighbor*w` terms read `w`, not
`scratch`), so each directional loop is independently parallelizable. The
imported `turbMod` symbols (`dvt, vort, prod, kwCD, f1`) are not written in the
loop body.

---

## Implementation

**File: `src/turbulence/saGammaRetheta.F90` — the only file changed.**

### 1. `Source()` — one parallel loop
Wrap the primal triple loop (the `#else` branch around line 286, closing around
line 741):

```fortran
#else
            !$OMP parallel do collapse(3) default(none) &
            !$OMP   private(i, j, k, <all per-cell scalar temporaries>) &
            !$OMP   shared(<global per-block arrays + loop-invariant scalars + index constants>)
            do k = 2, kl
                do j = 2, jl
                    do i = 2, il
#endif
```
and after the closing `end do`s in the `#else`:
```fortran
#else
                end do
            end do
        end do
        !$OMP end parallel do
#endif
```

### 2. `Viscous()` — three parallel loops
Apply the same wrapping to each of the three directional sweeps (k-dir ≈ line
910, j-dir ≈ 1101, i-dir ≈ 1291). Each gets its own
`!$OMP parallel do collapse(3) ... / !$OMP end parallel do`. They remain
sequential relative to each other (accumulation order preserved).

### 3. Do NOT touch
- `saGammaReThetaSolve()` — the DADI tridiagonal sweeps (`tdia3x3`, work arrays
  `bb/cc/dd/ff`) are inherently sequential; leave serial.
- `turbAdvection()` / `unsteadyTurbTerm()` in `turbUtils.F90` — shared with
  other models; out of scope for this change.
- The `#ifdef TAPENADE_REVERSE` branches.

### Private/shared classification
Use `default(none)` so the compiler forces every symbol to be classified
(prevents silent races).

- **private:** `i, j, k` and every scalar temporary written per cell — the full
  set declared at `saGammaRetheta.F90:214-250` (e.g. `fv1, fv2, ft2, ss, sst,
  nu, dist2Inv, chi, chi2, chi3, rr, gg, gg6, termFw, fwSa, term1, term2,
  term2_prod, term2_dest, dfv1, dfv2, dft2, drr, dgg, dfw, uux…wwz, div2, fact,
  sxx…syz, vortx/y/z, strainMag2, strainProd, vortProd, vortMag, strainMag,
  nutSA, rTurb, gammaLocal, reThetaTilde, gammaForSA, reS_val, reThetaC_val,
  fLength_val, fTurb_val, fOnset, fOnset1, vortLim, vortMagLim, pGamma, eGamma,
  velMag, velMag2, timeScale, reThetaT_target, thetaBL, deltaBL, delta,
  fWake_val, fThetaT, gammaEff, gammaTerm, pReTheta, yDist, uxhat…dUds,
  lambdaThetaLocal, dudx…dwdz, drTurb_dnu, dfTurb_dnu, dfOnset_dnu,
  dfOnset1_drT, dfOnset_dfOnset1, F1_val, base_val, inner_val, dFlength_dReT,
  dReThetaC_dReT, dfOnset1_dReT, dfOnset_dReT, dPgamma_dReT, pGamma_common,
  sech2_val`, and the per-direction diffusion temporaries in `Viscous()`).
- **shared (read-only in loop):** loop-invariant constants computed before the
  loop (`cv13, kar2Inv, cw36, cb3Inv, omegax, omegay, omegaz`), index constants
  (`idvt, itu1, itu2, itu3`, and the `dbg*` parameters), and all global
  per-block arrays (`w, si, sj, sk, sFaceI/J/K, rlv, vol, volRef, d2wall,
  scratch, qq, transitionDebug, iblank, bmt*` as used in the respective
  routine).
- Named `parameter`s (e.g. `xminn`, `f23`) need no classification.

`ii` (the Tapenade linear index) is declared but unused in the `#else` branch;
it does not appear in the primal pragma.

---

## Verification

1. **Compile:** `cd /home/mdo/MDOLab_3_v2/adflow_sa_gamma_rethetha && make -j`
   (`default(none)` will fail the build on any unclassified symbol — use that as
   the correctness gate.)
2. **Serial equivalence:** run a SA-γ-Reθ case with `OMP_NUM_THREADS=1`; result
   must match the pre-change baseline bit-for-bit (only pragmas added).
3. **Thread equivalence:** same case with `OMP_NUM_THREADS=4`; converged CL and
   final residuals must match the 1-thread run (to solver tolerance). A mismatch
   means a missing `private` — the prime symptom of an OpenMP race.
4. **Smoke:** run N iterations with threads, confirm no NaN.
5. **No Tapenade step required** — confirm only `!$OMP` comment lines and the
   primal `#else` branches changed (`git diff`).

---

## Why not blockette / pointers (recorded for posterity)

- **Copy-based blockette DADI:** would require adding `qq` assembly to
  `blockette.F90` (second copy of Jacobian physics) **and** copying `dw`→`scratch`,
  `qq_blk`→`qq` each iteration — added memory traffic on a memory-bound kernel,
  plus a divergence hazard. Rejected.
- **Pointer-based blockette:** avoids the copy but needs Fortran pointer
  remapping into global arrays under OpenMP; fragile and still duplicates the
  Jacobian code. Rejected.
- **Option B** reuses the one existing implementation, adds no traffic, and
  needs no Tapenade rerun. Chosen.
