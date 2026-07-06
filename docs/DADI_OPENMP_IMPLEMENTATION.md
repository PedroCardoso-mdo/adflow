# OpenMP Parallelization of the SA-γ-Re̅θt DD-ADI Turbulence Solve — Implementation & Rationale

**Branch:** `sa_gamma_rethetha`
**Scope:** the *implementation* of OpenMP in the decoupled DD-ADI turbulence
solve — the code, the design decisions, the points considered, and the reasoning
that led to each one. This document is implementation and justification only: it
does **not** report timings, scaling, or any benchmark.

All references point to the committed source
(`src/turbulence/saGammaRetheta.F90`, `src/turbulence/turbAPI.F90`).

---

## 1. Starting point — the serial structure I had to parallelize

The transition model advances a **3-equation, 3×3 block-coupled** system
(ν̃, γ, Re̅θt) with a γ-weighted SA production term (Eq. 41), using a coupled
diagonally-dominant ADI (DD-ADI) line solver with a per-cell 3×3 block
tridiagonal kernel (`tdia3x3`). As inherited, the entire turbulence phase was
serial.

The per-block driver is `saGammaReTheta_block` (`saGammaRetheta.F90:64`), invoked
from the block loop in `turbSolveDDADI` (`turbAPI.F90:69`). Its phases:

```
saGammaReTheta_block(resOnly):
    bcTurbTreatment                 # boundary Jacobians (shared code)
    allocate qq(2:il,2:jl,2:kl,3,3) # per-cell 3x3 block Jacobian
    Source                          # residual + qq assembly (source terms)
    saGRAdvection                   # first-order upwind advection
    unsteadyTurbTerm                # time term (shared code)
    Viscous                         # viscous residual + qq
    ResScale                        # residual scaling -> dw (monitoring)
    if (.not. resOnly):
        saGammaReThetaSolve         # DD-ADI line-sweep update
        saEddyViscosity / BCs       # shared code
    deallocate qq
```

The targets I owned and could parallelize are `Source`, `saGRAdvection`,
`Viscous`, and `saGammaReThetaSolve`. Everything else is shared/SA code (see D2).

---

## 2. Design decisions and rationale

These are the decisions that shaped the implementation, with the reasoning and the
points I weighed for each.

### D1 — Thread *within* a block, not *across* blocks
The block loop `do nn = 1, nDom` (`turbAPI.F90:69`) calls
`setPointers(nn, …)` (`turbAPI.F90:72`), which **rebinds the shared
`blockPointers` module** to block `nn`. Two blocks cannot be in flight at once
without racing on that global pointer state. So I placed the parallelism on the
**cell/line loops inside each block** and left the block loop serial.
*Point considered:* across-block threading would expose more coarse-grained
parallelism, but it would require making `blockPointers` thread-local — a large,
invasive change to shared infrastructure, rejected.

### D2 — Pragmas only on the primal `#else` branch (no AD regeneration)
Several routines carry a `#ifdef TAPENADE_REVERSE` differentiated branch. Every
OpenMP directive I added sits **only** in the primal `#else` branch, so the
differentiated source is byte-unchanged and Tapenade does **not** need to be
re-run; the adjoint cannot desync. OpenMP directives are comments to the AD tool,
so this is purely a placement discipline.

```
#ifdef TAPENADE_REVERSE
    ! flattened differentiated loop — untouched
#else
    !$OMP parallel do ...        ! OpenMP lives ONLY here
    do k ...; do j ...; do i ...
#endif
```

### D3 — Never touch SA-only or shared code; fork instead
`sa.F90`, the shared `turbAdvection`/`unsteadyTurbTerm`/`tdia3x3` in
`turbUtils.F90`, and the BC routines are off-limits (they are reused by the plain
SA model and by k-ω/SST). Where I needed a threaded version of shared logic —
advection — I did **not** edit the shared routine. I wrote a **specialized local
copy**, `saGRAdvection` (§3.3), fixed to 3 equations and first-order upwind, and
called it from the driver instead of the generic routine. The copy carries an
explicit *KEEP-IN-SYNC* header noting the source it was derived from.
*Point considered:* a fork is a maintenance liability (it can drift from the
shared routine), but it is the only way to thread advection without modifying code
shared with other models.

### D4 — Classify every loop by its data-dependence pattern first
Before adding a single pragma I classified each loop, because the pattern dictates
both *whether* it is safe and *how* to thread it:

```
for each loop nest over cells (i,j,k):
    if  iteration writes ONLY its own cell, reads only read-only shared data
        -> EMBARRASSINGLY PARALLEL          (Source, Viscous, qq-prep, update)
    elif iterations are independent along the OUTER indices, a 1-D solve
         couples the inner index
        -> PARALLEL OVER LINES              (the three ADI sweeps)
    elif iterations accumulate into a shared array across several passes
        -> THREAD EACH PASS, BARRIER BETWEEN PASSES   (advection k->j->i)
    else
        -> LEAVE SERIAL
```

The three resulting patterns are implemented in §3.

### D5 — Per-line tridiagonal buffers must be `private`
In the ADI sweeps each thread solves a different line and `tdia3x3` **writes into**
the working buffers `bb, cc, dd, ff`. If those were shared, threads would clobber
each other's tridiagonal system. They are therefore on the `private` clause — this
is the single most important correctness condition of the whole solve.

### D6 — Keep the three ADI sweeps as separate parallel regions
The j-, i-, and k-sweeps have a true data dependence (each consumes the `scratch`
the previous one produced). I implemented them as **three separate `parallel do`
regions**, so the implicit barrier that closes each region enforces the ADI
ordering. I deliberately did **not** fuse them into one region — fusing would drop
the barrier and corrupt the right-hand-side hand-off between directions.

### D7 — Minimize fork/join where it is free to do so
Within a single routine that has several independent directional loops (advection;
the viscous assembly), I used **one `parallel` region with several `!$OMP do`**
worksharing loops rather than one `parallel do` per loop. This forks the team once
instead of once per direction, while the implicit barrier between the `do`s still
provides the ordering I need.

### D8 — `schedule(static)` and `collapse`
The per-cell work is uniform, so `schedule(static)` gives even, low-overhead
distribution with no scheduling chatter. I used `collapse` to merge the cell loops
into one large iteration space so the thread team has enough work to balance even
when an individual outer loop is short.

---

## 3. The implementation — three patterns

### 3.1 Pattern A — embarrassingly-parallel per-cell assembly
Applies to `Source` (`:329`), `Viscous` (`:882`), the qq-prep
(`saGammaReThetaSolve` `:1630`), and the variable update (`:2159`). Each iteration
touches only its own cell; every loop-body temporary is `private`.

**(a) Source assembly** — demonstrative of the real body (`:329`–`:804`): every
cell recomputes its velocity gradients, the strain/vorticity invariants, the SA
production/destruction, and the γ/Re̅θt source, then stores the residual into
`scratch` and the 3×3 source Jacobian into `qq`. All of these intermediates must be
per-thread, hence the long `private` list.

```
!$OMP parallel do collapse(3) private(i,j,k,
!$OMP     uux..wwz, sxx..syz, vortx..vortz, strainMag, vortMag,   # velocity grads + invariants
!$OMP     nu, chi, fv1, fv2, ft2, ss, sst, rr, gg, fw,            # SA terms
!$OMP     gammaLocal, reThetaTilde, fOnset, fLength_val, ...,     # transition terms
!$OMP     pGamma, eGamma, pReTheta, term1, term2, ... )           # source pieces + jacobian
for k = 2..kl: for j = 2..jl: for i = 2..il:

    # --- reads: ONLY read-only shared arrays at this cell + neighbours ---
    grad(u,v,w)         <- w[i±1, j±1, k±1, ivx:ivz], metrics si/sj/sk, vol
    strainMag, vortMag  <- from the gradients
    (SA terms)          <- w[i,j,k,itu1], rlv[i,j,k], dist[i,j,k]
    (transition terms)  <- w[i,j,k,itu2], w[i,j,k,itu3], correlations()

    # --- writes: ONLY this cell's slice ---
    scratch[i,j,k, idvt+0..2] += S_nu, S_gamma, S_retheta          # residual
    qq[i,j,k, 1..3, 1..3]      = dS/dQ  (5 nonzero entries, P&Z §7.1) # block jacobian
!$OMP end parallel do
```

**(b) qq preparation / scaling** (`:1630`) — same shape, pure per-cell algebra on
the assembled `qq` and right-hand side (decouple unwanted blocks, add the source
dt-restriction to the diagonal, apply symmetric scaling and the CFL factor):

```
!$OMP parallel do collapse(3) schedule(static) private(i,j,k)
for k,j,i:
    if mode == decoupled:  zero the off-diagonal qq[i,j,k, m, n]      # m != n
    if restrict_dt:        qq[i,j,k, m, m] += srcLambda[i,j,k,m] / dtLimit
    qq[i,j,k, m, n] *= s[n]/s[m]        # symmetric scaling, off-diagonals
    qq[i,j,k, *, *] *= factor           # CFL factor, all 9 entries
    scratch[i,j,k, idvt+m] /= s[m]      # scale the RHS
!$OMP end parallel do
```

**(c) Variable update** (`:2159`) — applies the computed increment back to `w`,
with a per-cell exponential back-off so γ/Re̅θt stay in their physical bounds. Note
the inner `while` is *local* to the cell, so the loop is still per-cell independent:

```
!$OMP parallel do collapse(3) schedule(static) private(i,j,k, mm,gammaNew,gammaDelta,dampFactor)
for k,j,i:
    w[i,j,k,itu1] = max(w[i,j,k,itu1] + factor*scaleNu*scratch[i,j,k,idvt+0], 0)
    # gamma: shrink the step by theta^mm until gammaNew lands in [gammaLo, gammaHi]
    gammaDelta = factor*scaleGamma*scratch[i,j,k,idvt+1];  dampFactor = 1
    for mm = 1..maxIter:
        gammaNew = w[i,j,k,itu2] + dampFactor*gammaDelta
        if gammaNew in [gammaLo, gammaHi]: break
        dampFactor *= theta
    w[i,j,k,itu2] = clamp(gammaNew, gammaLo, gammaHi)
    # reTheta: same back-off against the lower bound
    ...
    w[i,j,k,itu3] = max(reThetaNew, reThetaLo)
!$OMP end parallel do
```

**Why it is race-free:** the write set of iteration `(i,j,k)` is exactly
`{scratch[i,j,k], qq[i,j,k], w[i,j,k], dw[i,j,k]}`, disjoint across iterations;
everything *read* (`si/sj/sk`, `vol`, `rlv`, `bmt*`, neighbour `w`) is never
written in the loop. The large `private` lists exist precisely because each cell
recomputes many scratch quantities (and damping counters) that would otherwise be
shared.

### 3.2 Pattern B — the ADI line sweeps (the core new parallelism)
`saGammaReThetaSolve` (`:1512`). The solve is parallel **across lines**: the
j-sweep loops `do k; do i; do j(inner)`, so the outer `(k,i)` pairs are
independent lines; likewise `(k,j)` for i and `(i,j)` for k. I collapse the
**outer two** loops and thread them; the inner solve-direction loop and the
tridiagonal solve stay **serial per line**.

Demonstrative body of one sweep (the j-sweep, `:1699`–`:1847`); the i- and
k-sweeps are identical up to the swept index. `bb`/`dd` are the sub/super-diagonals,
`cc` the 3×3 diagonal blocks, `ff` the right-hand side — **all `private`** so each
thread owns one tridiagonal system (D5):

```
!$OMP parallel do collapse(2) schedule(static) firstprivate(qs)
!$OMP   private(i,j,k, voli,volmi,volpi, xm..zp, xa,ya,za, ttm,ttp,
!$OMP           nu,nuTilde,fv1,nut, cdm,cdp, cdm_gamma,cdp_gamma, cdm_rt,cdp_rt,
!$OMP           c1m,c1p,c2m,c2p,c3m,c3p, uu,um,up, rblank, bb,cc,dd,ff)
for k = 2..kl:                              #  outer two loops -> threads
  for i = 2..il:                            #  each (k,i) = one independent j-line

    # ---- build the tridiagonal system for THIS line ----
    for j = 2..jl:                          #  inner sweep index -> serial
        face metrics ttm,ttp  <- sj[i,j-1,k], sj[i,j,k], vol         # geometry
        cdm,cdp        = SA diffusion coefficients    (nu, nut at j-1/j/j+1)
        cdm_gamma..    = gamma / reTheta diffusion coefficients
        uu             = face-normal velocity - grid velocity qs     # advection speed
        um,up          = upwind split of uu                          # first-order upwind
        bb[*,j]        = -(diffusion_minus) - up                     # sub-diagonal (3 eqns)
        dd[*,j]        = -(diffusion_plus)  + um                     # super-diagonal
        rblank         = iblank[i,j,k]
        cc[1..3,1..3,j]= qq[i,j,k, *, *]   (off-diagonals * rblank)  # 3x3 diagonal block
        ff[1..3,j]     = scratch[i,j,k, idvt+0..2] * rblank          # RHS

    # ---- solve this line, then hand RHS to the next direction ----
    tdia3x3(2, jl, bb, cc, dd, ff)          #  serial block-Thomas, PRIVATE buffers
    for j = 2..jl:
        scratch[i,j,k, idvt+m] = sum_n qq[i,j,k,m,n] * ff[n,j]       # RHS for i-sweep
!$OMP end parallel do
# <-- implicit BARRIER (D6): all j-updates committed before the i-sweep reads scratch

# i-sweep: identical, lines = (k,j), swept index i      -> BARRIER
# k-sweep: identical, lines = (i,j), swept index k, then scratch = ff (final update)
```

**Points considered:**
- The block-Thomas (`tdia3x3`) is inherently sequential along its line; I did not
  attempt a parallel/pipelined solver. The parallelism is the *set of independent
  lines*, which is ample and needs no synchronization.
- `qs` (grid-velocity term) is `firstprivate` so its `qs = zero` default survives
  per-thread when `addGridVelocities` is false.

### 3.3 Pattern C — advection: accumulation across directions
`saGRAdvection` (`:2203`), the local fork from D3. Each of the three directional
passes **accumulates** into the shared `scratch`/`qq`, so the passes cannot run
concurrently with each other — but each pass is itself per-cell parallel. I used
one `parallel` region (D7) with three `!$OMP do`, the barrier between them ordering
the accumulation.

Demonstrative body (`:2240`–`:2414`): one direction shown; `flux_dir(cell)` expands
to the first-order upwind difference plus the implicit diagonal contribution, with
the boundary-face Jacobian `bmt*` folded in only on the domain faces.

```
!$OMP parallel private(i,j,k, ii,jj,kk, qs,voli, xa,ya,za, uu, dwt, impl)  # fork ONCE
  !$OMP do collapse(3)                       # ---- k (zeta) direction ----
  for k,j,i:
      uu = face-normal velocity at k-face - qs            # advection speed
      for ii = 1..3:                                      # the 3 transport eqns
          if uu > 0:                                      # upwind from k-1
              dwt = w[i,j,k, var] - w[i,j,k-1, var]
              scratch[i,j,k, idvt+ii-1] -= uu * dwt
              qq[i,j,k, ii, ii]         += uu             # implicit diagonal
              if k == 2:                                  # domain face -> add bmt jacobian
                  impl[:] = max-clamped bmtk1[i,j, var, :]
                  qq[i,j,k, ii, :] += uu * impl[:]
          else:                                           # upwind from k+1
              dwt = w[i,j,k+1, var] - w[i,j,k, var]
              scratch[i,j,k, idvt+ii-1] -= uu * dwt
              qq[i,j,k, ii, ii]         -= uu
              if k == kl: qq[i,j,k, ii, :] -= uu * (max-clamped bmtk2[...])
  !$OMP end do                               # implicit BARRIER: k-pass done everywhere
  !$OMP do collapse(3)  ... j (eta) direction ...   # same form, sj / bmtj1,bmtj2
  !$OMP end do                               # implicit BARRIER
  !$OMP do collapse(3)  ... i (xi)  direction ...   # same form, si / bmti1,bmti2
  !$OMP end do
!$OMP end parallel                           # join ONCE
```

Within any one `do`, each cell updates only its own `qq[cell]`/`scratch[cell]`, so
there is no race inside a pass; the barriers serialize the three `+=` passes so the
k→j→i accumulation into `qq` is deterministic. `ii,jj,kk` are in the `private` list
because the inner equation loops use them as local indices.

---

## 4. Why the whole scheme is correct (the invariant)

The entire parallelization rests on one invariant: **once the halos are exchanged,
every owned cell is independent.** Each cell writes only its own
`scratch`/`qq`/`w`/`dw`; the flow state `w`, the boundary Jacobians `bmt*`, and the
grid metrics are read-only throughout the solve. The only places that are *not*
purely per-cell — the ADI sweeps (line-coupled) and advection (multi-pass
accumulation) — are handled by line-private buffers (D5) and region barriers
(D6/D7) respectively. No reductions, locks, atomics, or critical sections are
needed anywhere.

The 3×3 transition coupling is **node-local**: a dense 3×3 `qq` per cell plus the
per-line private buffers. It adds arithmetic but no new data sharing and no new
communication — the halo exchange `whalo2(…, nt1, nt2, …)` simply carries a wider
payload (3 variables instead of 1), with the same message pattern. This is what
makes the coupled system a clean OpenMP target.

---

## 5. How correctness was checked (method)

The validation *method* (independent of any specific run): execute the identical
case at one thread and at many threads and require the turbulence solution to be
**bit-for-bit identical**. Because the scheme uses no order-dependent reductions,
a correct parallelization is exactly reproducible; any thread-count-dependent drift
would localize a missing `private` or a misplaced barrier. The production build
(OpenMP enabled, timing macro undefined) was also confirmed to compile and import
cleanly so the instrumentation leaves zero residue. (No quantitative results are
reported here by design.)

---

## 6. Summary of what I implemented

- Classified every DADI loop by data dependence (D4) and applied the matching
  OpenMP pattern: per-cell assembly (Pattern A), line-parallel ADI sweeps
  (Pattern B), barrier-ordered advection accumulation (Pattern C).
- Threaded `Source`, `Viscous`, the qq-preparation, the three ADI sweeps, and the
  variable update inside `saGammaReThetaSolve`; authored a threaded, SA-γ-Re̅θt-local
  advection routine (`saGRAdvection`) to avoid editing shared code (D3).
- Enforced correctness structurally: line-private tridiagonal buffers (D5),
  separate parallel regions to preserve ADI ordering (D6), single-fork
  multi-`do` regions to cut fork/join (D7).
- Confined all parallelism to a single block at a time (D1) and to the primal code
  path only, so no adjoint regeneration is required (D2).
