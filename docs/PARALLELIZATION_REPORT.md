# Parallelization of the SA-γ-Re̅θt Transition Solver vs. Baseline ADflow

**Branch:** `sa_gamma_rethetha`  ·  **Model:** SA-noft2-γ-Re̅θt (Piotrowski & Zingg 2020)
**Scope of this report:** the turbulence-phase timing instrumentation, the OpenMP
parallelization, and the blockette tiling — strictly *what differs from baseline
ADflow and why*. Facts are drawn from code and `git`; inferred rationale is flagged
**[INFERENCE]**.

The transition model adds **2 transport equations** (γ, Re̅θt) on top of the SA
working variable ν̃, giving a **3-equation, 3×3 block-coupled** per-cell system, a
γ-weighted SA production term (paper Eq. 41), solved by a **coupled DD-ADI line
solver** (`tdia3x3`). That structural difference from the scalar SA baseline is the
root cause of every divergence documented below.

---

## 0. Orientation — the two solve paths and where each was instrumented

ADflow's decoupled RANS solve has two distinct turbulence touchpoints, and this
branch treats them as two separate "contexts":

| Context | Routine chain | What it does | Called from |
|---------|---------------|--------------|-------------|
| **A — DADI solve** | `turbSolveDDADI` → `saGammaReTheta_block(resOnly=.false.)` → `Source`/`saGRAdvection`/`Viscous`/`ResScale` → `saGammaReThetaSolve` | builds the 3×3 residual + Jacobian and does the implicit ADI line-sweep update | the turbulence sub-iteration loop |
| **B — residual eval** | `blocketteRes` → `blocketteResCore` (blk=True) **or** `blockResCore` (blk=False) → `saGammaRetheta_block(resOnly=.true.)` | evaluates the residual only (no update), many times, inside ANK's matrix-free Krylov | `FormFunction_mf` / NK / MG |

The whole instrument is built around keeping these two apart, because they
parallelize and scale completely differently.

---

## 1. Turbulence-phase timing — a *separate* instrument from the baseline ANK profiler

### (a) Files & routines
- **`src/modules/turbTiming.F90`** (new module — did not exist in baseline).
- Hook sites: `src/turbulence/turbAPI.F90` (`turbSolveDDADI`),
  `src/turbulence/saGammaRetheta.F90` (`saGammaReTheta_block`),
  `src/NKSolver/blockette.F90` (`blocketteRes` / `blocketteResCore` / `blockResCore`),
  `src/solver/solvers.F90` (reset at `t0Solver`, print in the Solution Timings block).
- Added in commit **`bbff8d42`** *"Parallelize + instrument SA-gamma-Retheta DADI
  turbulence solve"* — *"Phase A: TURB_TIMING-gated wall-clock instrumentation …
  zero residue when the macro is undefined."* The most recent (uncommitted) work
  extends Context-B into the fused blockette path (§3c).

### (b) Real code excerpt — the timing mechanism (`turbTiming.F90:77–98`)

```fortran
subroutine turbTic(id)
    integer(kind=intType), intent(in) :: id
    tTurbStart(id) = mpi_wtime()
end subroutine turbTic

subroutine turbToc(id)
    integer(kind=intType), intent(in) :: id
    tTurbAccum(id) = tTurbAccum(id) + mpi_wtime() - tTurbStart(id)
end subroutine turbToc
```

Every call site is gated, e.g. in `saGammaReTheta_block` (`saGammaRetheta.F90:99–104`):

```fortran
#ifdef TURB_TIMING
        if (.not. resOnly) call turbTic(T_SOURCE)
#endif
        call Source
#ifdef TURB_TIMING
        if (.not. resOnly) call turbToc(T_SOURCE)
#endif
```

### (c) What `TURB_TIMING` gates
Per the module header (`turbTiming.F90:7–10`): *"ALL call sites are wrapped in
`#ifdef TURB_TIMING`, so when TURB_TIMING is not defined nothing here is ever called
(zero overhead). To enable, add `-DTURB_TIMING` to FF90_FLAGS in config.mk."* The
production build (no macro) compiles every hook out to nothing — verified in
`OPENMP_DADI_TIMING_PLAN.md` §STATUS: *"Prod build WITHOUT `-DTURB_TIMING` verified:
compiles clean + imports; timing call-sites compile out … zero runtime overhead.
Build deterministic (md5 stable)."*

### (d) Phases timed separately
15 timer IDs (`turbTiming.F90:19–51`). **Context A** (DADI): `T_SOURCE`, `T_VISCOUS`,
`T_ADV`, `T_UNSTEADY` (these four sum to "residual + qq assembly"), `T_DADI` (the
line-sweep solve), `T_RESSCALE`, `T_TURBHALO` (the `whalo2` turbulence halo),
`T_TURBTOTAL`. **Context B** (residual eval): `T_RESID_FLOW`, `T_RESID_TURB`,
`T_RESID_BOTH` (fused), `T_RESID_HALO`, `T_RESID_TOTAL`, plus two **thread-summed**
accumulators `T_BLKT_TURB_TS` / `T_BLKT_FLOW_TS` (§3c). So *Source, Viscous,
Advection, DADI line-sweep, and turb halo are each timed on their own line*, exactly
as requested.

The print routine reduces each timer two ways — `mpi_max` (critical path) and
`mpi_sum` (→ mean) — so a single high-rank straggler shows up as an
`imb = max/mean` ratio instead of being hidden (`turbTiming.F90:112–121`).

### (e) Why a separate instrument was needed — contrast with baseline
> **Baseline did X:** the unmodified ADflow uses the ANK/NK profiling timers, which
> wrap the *flow* Newton-Krylov loop (`FormFunction`, GMRES, preconditioner) — one
> linear loop with one residual.
>
> **Here we do Y:** a dedicated `turbTiming` module with an explicit A/B phase split.
>
> **Because Z:** the DADI turbulence solve has a *different control structure* than
> the ANK flow loop. It is **not** a Krylov loop — it is a directional **ADI line
> sweep** (`saGammaReThetaSolve`) wrapped in a sub-iteration loop (`turbSolveDDADI`),
> with assembly (`Source`/`Viscous`/advection/unsteady), an implicit tridiagonal
> solve, a residual-scaling step, and a turbulence-only halo. The ANK profiler has
> no bins for any of these. **[INFERENCE]** the A/B split exists specifically because
> the *same* `saGammaReTheta_block` routine is reused for two purposes — a true DADI
> update (`resOnly=.false.`) and a residual-only eval inside ANK
> (`resOnly=.true.`) — and the fine Context-A timers fire **only when
> `resOnly=.false.`** (see the `if (.not. resOnly)` guard in the excerpt above) so a
> residual-eval call cannot pollute the DADI breakdown. This guard *is* the
> separation mechanism; the baseline ANK profiler has no equivalent because nothing
> in the flow path is dual-purpose this way.

### Pseudocode — the instrument

```
# Context A (DADI update): only when resOnly == False
turbSolveDDADI:
  tic(TURBTOTAL)
  for sub-iter in 1..nSubIterTurb:
     saGammaReTheta_block(resOnly=False):
        tic(SOURCE);  Source;        toc(SOURCE)
        tic(ADV);     saGRAdvection; toc(ADV)
        tic(VISCOUS); Viscous;       toc(VISCOUS)   # etc.
        count(DADI); tic(DADI); saGammaReThetaSolve; toc(DADI)
     tic(TURBHALO); whalo2(...); toc(TURBHALO)
  toc(TURBTOTAL)

# Context B (residual eval): bin by which flags are set
blocketteRes:
  resTid = BOTH if (flowRes and turbRes) else FLOW if flowRes else TURB
  count(resTid); tic(resTid); blocketteResCore(...); toc(resTid)

print: MPI_reduce(max) and MPI_reduce(sum)->mean; report imb=max/mean per phase
```

---

## 2. OpenMP — what is NEW in this branch, plus a real correctness fix

### (a) Files, routines, and the commit chain (`git blame` / `git log`)

| Commit | Date | OpenMP work |
|--------|------|-------------|
| `d2b837b6` | 2026-06-02 | *"Fix OpenMP thread-safety for blockette SA and SA-gamma-rethetha"* — the **correctness fix** (§2c) |
| `abefec77` | — | *"Add OpenMP to SA-gamma-rethetha DADI residual/Jacobian assembly"* |
| `057b000b` | — | *"Revert sa.F90 constants to original runtime form"* — undoes one half of `d2b837b6` (§2c) |
| `8c4e6652` | — | *"… implement OpenMP in turbulence DADI assembly"* |
| `bbff8d42` | 2026-06-09 | *"Parallelize + instrument …"* — the bulk: solve sweeps, advection, fork/join reduction (55 pragma lines) |

> **Baseline did X:** stock ADflow's turbulence DADI path (`saSource`, `saViscous`,
> `turbAdvection`, the scalar SA line solve) is **serial**; OpenMP in baseline lives
> only at the block level in `blockette.F90`.
>
> **Here we do Y:** OpenMP is added *inside* the SA-γ-Re̅θt DADI kernels — the 3×3
> assembly (`Source`, `Viscous`), the qq-prep + j/i/k ADI sweeps in
> `saGammaReThetaSolve`, and a brand-new threaded `saGRAdvection`.
>
> **Because Z (the "free cores" argument, `OPENMP_DADI_TIMING_PLAN.md` §Context):**
> in production the run is blockettes+OpenMP on the *flow* residual; during the DADI
> turbulence phase *"the rank's cores are otherwise idle, so any DADI-phase speedup
> is free."* Baseline never threaded this phase because for scalar SA it is a small
> add-on; the 3-equation coupled system makes it heavy enough to be worth threading.

Architectural choice (Decision **D1**, `OPENMP_DADI_TIMING_PLAN.md`): pragmas go
**directly on the existing primal loops**, *not* through the blockette path, because
`Source`/`Viscous` already write the per-block `scratch`/`qq` the solve consumes
(*"nothing to copy"* — the kernel is memory-bound) and routing through blockette
would force a **second copy of the 3×3 transition Jacobian** into `blockette.F90`, a
divergence hazard. Crucially, all pragmas sit in the **primal `#else` branch only**,
so the Tapenade-differentiated math is untouched → **no AD regeneration** (CLAUDE.md
rule #6 honored).

### (b) Real code excerpt — the per-line ADI sweep with private tridiagonal buffers
(`saGammaRetheta.F90:1696–1710`)

```fortran
! Parallel across (k,i) lines: each thread owns a distinct j-line and
! its own tridiagonal buffers bb/cc/dd/ff (private). qs is firstprivate
! so the qs=zero default is preserved when addGridVelocities is false.
!$OMP parallel do collapse(2) schedule(static) firstprivate(qs) &
!$OMP private(i, j, k, voli, volmi, volpi, xm, ym, zm, xp, yp, zp, &
!$OMP xa, ya, za, ttm, ttp, cnud, cam, cap, nutm, nutp, &
!$OMP ... bb, cc, dd, ff)
do k = 2, kl
    do i = 2, il
        do j = 2, jl            ! inner solve-direction loop stays serial
```

The DADI **tridiagonal solve itself stays serial per line** (Decision D4); only the
*outer two* loops over independent lines are collapsed and threaded. `bb/cc/dd/ff`
(the per-line 3×3 tridiagonal working buffers) are `private` — this is the key
correctness requirement, since `tdia3x3` writes into them.

### (c) The OpenMP CORRECTNESS fix — broken vs. fixed

The plan's *"OpenMP implemented and corrected"* refers to commit **`d2b837b6`**
*"Fix OpenMP thread-safety for blockette SA and SA-gamma-rethetha"*. Two distinct
bugs, with opposite resolutions:

**Bug 1 — missing THREADPRIVATE / private (FIXED, retained).** Block-level OpenMP in
`blockette.F90` shared per-thread scratch state. Symptom (Decision D3): *"1-thread
CL=0.0328 vs 4-thread CL=0.0094 (race)."* Fix (`git show d2b837b6`):

```fortran
+    !$OMP THREADPRIVATE(singleHaloStart, doubleHaloStart, nodeStart)
+    !$OMP THREADPRIVATE(dtl, sFaceI, sFaceJ, sFaceK)
...
-        !$OMP parallel do private(i,j,k,l) collapse(2)
+        !$OMP parallel do private(i,j,k,l,ii,jj,kk) collapse(2)
```

The block loop indices `ii,jj,kk` had been **shared** — every thread clobbering the
same block coordinates — and seven halo/face-state scalars lacked `THREADPRIVATE`.

**Bug 2 — the sa.F90 constant episode (REVERTED, do NOT redo).** The *same* commit
also tried to silence a benign race by converting `sa.F90` module scalars
(`cv13`, `kar2Inv`, `cw36`, `cb3Inv`) to compile-time `parameter`s — and **hardcoded
wrong values** (Decision D2): `cw36` became `0.3⁶ = 0.000729` instead of
`2.0⁶ = 64.0` (off by ~88000×), `cb3Inv` became `1/0.622` instead of `1/0.6667`.
This corrupted SA destruction/diffusion — but **only on the blockette path**, which
reads the `sa` module constants; the **DADI path recomputes them at runtime and was
never affected**. The blockette call sites that re-derived the constants were
removed in `d2b837b6`:

```fortran
-        ! Set model constants
-        cv13 = rsaCv1**3
-        kar2Inv = one / (rsaK**2)
-        cw36 = rsaCw3**6
-        cb3Inv = one / rsaCb3
+        ! Model constants cv13, kar2Inv, cw36, cb3Inv are now PARAMETERs in sa module
```

Per CLAUDE.md rule #2 (never touch SA code) this was **reverted in `057b000b`**
*"Revert sa.F90 constants to original runtime form"* — `sa.F90` is byte-identical to
upstream again, and the benign same-value module-scalar races are accepted as
not-bugs.

> **Contrast:** baseline ADflow's scalar SA blockette path is already thread-safe as
> shipped; the race surfaced **here** because the 3-equation coupled assembly enlarged
> the per-thread working set (extra halo/face scalars, the `ii,jj,kk` block indices)
> and *exposed* state that scalar SA never shared at block granularity.

### (d) Parallel structure & race-avoidance — pseudocode

```
# Per-block driver (Context A), once halos are in place:
saGammaReTheta_block:
    Source        # OMP parallel do collapse(3): each cell writes only its
                  #   own scratch(i,j,k,*) and qq(i,j,k,*,*); w/metrics read-only
    saGRAdvection # ONE omp parallel region, three omp do (k,j,i directions);
                  #   implicit barrier after each do => k->j->i qq accumulation safe
    Viscous       # 3 directional omp do fused into ONE parallel region (3->1 fork/join)
    saGammaReThetaSolve:
        qq-prep   # omp parallel do collapse(3): per-cell scaling, independent
        j-sweep   # omp parallel do collapse(2) over (k,i); private bb/cc/dd/ff
        i-sweep   # omp parallel do collapse(2) over (k,j); sweeps run SEQUENTIALLY
        k-sweep   # omp parallel do collapse(2) over (i,j);   (preserves ADI order)
```

Race-avoidance rests on one invariant: **once the halos are exchanged, every owned
cell is independent.** Each cell writes only its own `scratch`/`qq`; `w`, the
boundary Jacobians `bmt*`, and the grid metrics are read-only. Within a sweep each
thread owns a distinct line and its own tridiagonal buffers; the three sweeps run as
**separate parallel regions in sequence**, so the directional ADI ordering is
preserved. Determinism was verified (Phase C / D5): *"omp 1 vs 12 @ 200 cyc,
blk=False → bit-identical CD 0.016210."*

### (e) The 3×3 block coupling — more arithmetic, zero extra communication

The transition coupling lives entirely in the **per-cell 3×3 `qq` Jacobian**
(`saGammaRetheta.F90:1630–1684`):

```fortran
!$OMP parallel do collapse(3) schedule(static) private(i, j, k)
do k = 2, kl; do j = 2, jl; do i = 2, il
    ! Symmetric scaling (§4): qq(m,n) *= s_n / s_m  (off-diagonals only)
    qq(i, j, k, 1, 2) = qq(i, j, k, 1, 2) * s12
    qq(i, j, k, 1, 3) = qq(i, j, k, 1, 3) * s13
    ...
    ! CFL factor scaling (all 9 entries), then RHS scaling
    qq(i, j, k, 1, 1) = factor * qq(i, j, k, 1, 1)
    ...
```

and the line solve calls `tdia3x3` (`saGammaRetheta.F90:1832, 1984, 2136`) instead of
the scalar `tdia` baseline uses. The coupling is **node-local**: it is a dense 3×3
block *at each cell*, factored exactly (`saGammaRetheta.F90:2713`: *"Since A13=A31=A32=0,
the 3×3 factorizes exactly"*). The γ-weighting of SA production
(`saGammaRetheta.F90:495–512`, paper Eq. 41) is also purely local:

```fortran
gammaForSA = min(max(w(i, j, k, itu2), zero), one)
...
term2 = gammaForSA * term2_prod + term2_dest   ! gamma multiplies production ONLY
```

> **Contrast / why no extra comm:** baseline SA exchanges **1** turbulence variable
> per halo; here `whalo2` exchanges **3** (ν̃, γ, Re̅θt) — the same *number* of halo
> messages, just a wider payload. The 3×3 coupling that turns a scalar tridiagonal
> into a block tridiagonal is **entirely within each cell's `qq`** and within each
> line's private `bb/cc/dd/ff`. It therefore adds floating-point work (9 Jacobian
> entries vs. 1, a 3×3 factorization vs. a scalar divide) **but no new
> message-passing** — the halo pattern is unchanged. **[INFERENCE]** this is exactly
> why the coupled system is an attractive OpenMP target: the added cost is arithmetic
> on data already resident per rank, which threads well, rather than communication,
> which does not.

---

## 3. Blockettes in the turbulence path

### (a) Is the blockette tiling applied to the SA-γ-Re̅θt residual? — Yes, but guarded off by default

The blockette SA-GR residual routines were added in commit **`ed5f2feb`**
*"blocket on sa-gamma-rethetha"* (+670 lines in `blockette.F90`): `saGammaRethetaSource`,
`saGammaRethetaAdvection`, `saGammaRethetaViscous`, `saGammaRethetaResScale`
(`blockette.F90:6993, 7227, 7380, …`), dispatched from `blocketteResCore`:

```fortran
case (spalartallmarasnoft2gammaretheta)
    call saGammaRethetaSource
    call saGammaRethetaAdvection
    call saGammaRethetaViscous
    call saGammaRethetaResScale
```

**The flag `ADFLOW_ALLOW_SAGR_BLOCKETTES`** is read in Python, not Fortran
(`adflow/pyADflow.py:6597–6604`):

```python
# SA-GR model is not (fully) validated in blocketteResCore, so by
# default force the blockResCore path for correct residual computation.
# Set ADFLOW_ALLOW_SAGR_BLOCKETTES=1 to bypass this guard ...
_allow_sagr_blk = os.environ.get(
    "ADFLOW_ALLOW_SAGR_BLOCKETTES", "0").lower() in ("1", "true", "yes")
if turbModel == "SA-noft2-Gamma-Retheta" and not _allow_sagr_blk:
    self.setOption("useBlockettes", False)
```

So although `useBlockettes` defaults to `True` in ADflow generally
(`pyADflow.py:5749`), for **this** turbulence model it is **force-disabled** unless
the env var is set. The verdict (`OPENMP_DADI_TIMING_PLAN.md` D2/D5) is that the
blockette SA-GR path is **correct** (CD=0.007010 ≈ blk=False 0.007009) — the guard is
*"conservative, not because the path is wrong."*

> **Contrast:** baseline ADflow runs blockettes on by default for SA. Here the new
> 3-equation residual is gated behind an env var **[INFERENCE]** because it is newer,
> less validated, and — on the only mesh tested so far — a net loss (§3c).

### (b) How the tiling was adapted for the 3-equation coupled system

The blockette tile geometry is **unchanged**: `BS = 8`, working blocks
`bbil = BS+1`, halos `bbib = BS+3` (`blockette.F90:9–12`). What changed is the
**working set per tile**: the SA-GR blockette routines write **3** turbulence
residual slots (ν̃, γ, Re̅θt) and the dense **3×3 `qq`** per cell, versus the scalar
SA blockette which touches **1** turbulence slot and a scalar Jacobian. The flow
state `w`, dissipation `dss`, metrics, etc. are tiled exactly as in baseline; only
the turbulence portion of the per-tile arithmetic and `qq` footprint grew ~3×.

> **[INFERENCE]** the larger turbulence working set per `8³` tile is the mechanism by
> which blockettes are *expected* to pay off on cache-pressured 3D meshes (the tile
> keeps the now-3× turbulence state hot in L2) — but also why the thin 2D case sees
> no benefit (the tile working set never spills cache there to begin with).

### (c) The fused-residual timing problem the 3-equation system created — and the uncommitted fix

Inside `blocketteResCore` the flow flux and the SA-GR turbulence residual are
computed in **one fused OpenMP sweep** (`blockette.F90:673–696`) that *"cannot be
split cheaply or race-free"* (`turbTiming.F90:32–40`). A pure `mpi_wtime()` wall timer
cannot attribute the fused wall-time to "flow" vs "turb." The **uncommitted** change
to `blockette.F90` / `turbTiming.F90` solves this with a **thread-summed OpenMP
reduction** (`blockette.F90:670–696`):

```fortran
!$OMP parallel do private(i,j,k,l,ii,jj,kk,t0sec) collapse(2) reduction(+:tsTurb,tsFlow)
...
#ifdef TURB_TIMING
    t0sec = mpi_wtime()
#endif
    if (equations == RANSEquations .and. turbRes) then
        ... call saGammaRethetaSource/Advection/Viscous/ResScale ...
    end if
#ifdef TURB_TIMING
    tsTurb = tsTurb + mpi_wtime() - t0sec   ! per-thread CPU-seconds, reduced(+)
#endif
```

After the region closes, `turbAddTS(tsTurb, tsFlow)` stores the two thread-summed
totals; `printTurbTiming` forms the ratio `rTurb = TURB_TS/(TURB_TS+FLOW_TS)` and uses
it to split the *fused wall-time* `T_RESID_BOTH` into a turb share and a flow share
(`turbTiming.F90:139–143`). This closes the gap the plan flagged as **MISSING**:
*"Context-B timers don't instrument the blockette path — `T_RESID_FLOW`/`T_RESID_TURB`
read 0.000 under blk=True … hooked on `blockResCore`, not `blocketteResCore`."*

> **Contrast:** the baseline ANK profiler never had to separate a fused flow+turb
> sweep, because in baseline the turbulence is solved *separately* (DADI) and the
> blockette residual is effectively flow-only for scalar SA. The 3-equation fused
> residual is unique to this branch, so a thread-summed-ratio splitter had to be
> invented to attribute its cost.

### (d) git evidence on 2D-vs-3D behaviour — blockettes *hurt* on thin 2D

The full 2×2 matrix (commit **`03c95173`** *"Record full 2x2 (blk x OMP) convergence
matrix + blockette verdict"*, table in `OPENMP_DADI_TIMING_PLAN.md` §STATUS,
RUN_NCYCLES=25000, all reach CD≈0.007009):

| # | blk | omp | wall (s) | CD       |
|---|-----|----:|---------:|----------|
| 1 | F   | 1   | 1082     | 0.007009 |
| 2 | F   | 12  | **1015** | 0.007009 |
| 3 | T   | 1   | 3798     | 0.007010 |
| 4 | T   | 12  | 1856     | 0.007010 |

Verdict (Decision D5): *"blk=False + 12 threads (1015s) is the best config on this
mesh. Blockettes HURT here even with threads (1856 > 1015s): the blockette SA-GR
trajectory needs more outer iters (~6700 vs ~5200) and cache-blocking overhead isn't
amortized on a small 2D case … Blockettes are expected to pay off only on large
cache-pressured 3D meshes."* The open item, flagged in both docs and the new
`PARALLEL_CHARACTERIZATION_PLAN.md`, is the **missing 3D data point** —
*"advection 1.36× and the blockette verdict are both thin-2D-limited."*

---

## 3.5 Connection points to the existing parallel infrastructure — and what I did *not* do

This section is the honest map of *where* the new code plugs into ADflow's existing
parallel machinery, written from the author's side: what is genuinely new parallel
code, what is just a wider argument into infrastructure that was already
equation-count-generic, and what was left un-generalized / faked / disabled. The
equation-count-independent loops are deliberately *not* the focus — the connection
points and the gaps are.

### CP-1 — Variable-count registration: the one integer-config point (extended, trivially)
`src/inputParam/inputParamRoutines.F90:2152–2154`:

```fortran
case (spalartallmarasnoft2gammaretheta)
    nw = 8        ! 5 flow + 3 turbulence
    nt2 = 8       ! (nt1 = 6 set above, so turbulence range = 6..8)
```

This 3-line `case` is the *entire* declaration that "3 turbulence variables exist."
Everything downstream that is range-driven — state allocation `w(:,:,:,1:nw)`, the MG
variable loops, the halo packing — keys off `nt1/nt2/nw` and therefore "just works."
**Admission:** this is not parallel work; it is me telling the existing generic
infrastructure that the turbulence range is now `6..8` instead of `6..6`. I wrote
almost nothing here and inherited almost everything.

### CP-2 — MPI halo exchange: 100% inherited, zero new comm code
`src/turbulence/turbAPI.F90:112`:

```fortran
call whalo2(groundLevel, nt1, nt2, .false., .false., .true.)
```

`whalo2(level, start, end, …)` is range-parametrized (`haloExchange.F90:109`). I did
**not** touch `whalo2`, the pack/unpack, the MPI buffers, or the communication
schedule. Passing `nt1..nt2 = 6..8` makes the existing exchange move **3** turbulence
variables per halo instead of **1** — same messages, wider payload (this is the
"no extra communication" point in §2e, seen from the wiring side). **Admission: there
is no new MPI/parallel-communication code in this branch at all** — the 3-equation
halo is entirely the baseline infrastructure with a wider index range.

### CP-3 — Block-level OpenMP driver in `blockette.F90`: repaired, not authored
The block loop `!$OMP parallel do … collapse(2)` over `(ii,jj,kk)` tiles with its
THREADPRIVATE module-array working set is **baseline infrastructure**. My
contribution was two-fold and both are *extensions of existing infra, not new infra*:
1. The `d2b837b6` thread-safety **fix** — adding the missing
   `THREADPRIVATE(singleHaloStart, …, sFaceK)` and putting `ii,jj,kk` in `private`
   (§2c). I repaired the existing driver so the wider 3-equation working set is
   race-free; I did not write the driver.
2. The SA-GR blockette residual routines (`saGammaRethetaSource/Advection/Viscous/
   ResScale`) run *inside* that existing parallel block loop, writing the global
   per-cell state. They inherit the parallelism — they are dispatched from
   `blocketteResCore`, not separately parallelized.

**Admission:** I authored none of the block-parallel driver; I extended its
THREADPRIVATE set and dropped my residual kernels into it.

### CP-4 — DADI line-sweep OpenMP: this *is* the new parallel code I wrote
`saGammaReThetaSolve` (the qq-prep `collapse(3)` and the j/i/k `collapse(2)` sweeps
with private `bb/cc/dd/ff`) is the one place the branch adds genuinely new
parallelism — baseline DADI is serial. But two honest limits:
- **The tridiagonal solve itself stays serial.** `tdia3x3` is called per line; I
  parallelized only the *outer* loop over independent lines, not the solve. There is
  no parallel/pipelined block-Thomas — the line solver is sequential by design.
- The three sweeps are three separate parallel regions in sequence (fork/join ×3 in
  the solve, on top of assembly), accepted to preserve ADI ordering, not optimized
  away.

### CP-5 — The coupled-Jacobian and block solver: new *serial* numerics, hardcoded to 3
- `qq` is allocated **`allocate(qq(2:il,2:jl,2:kl,3,3))`** (`saGammaRetheta.F90:96`) —
  literally `3,3`, **not** `nwt,nwt`. **Admission: the block size is hardcoded, not
  generalized to the variable count.** A different turbulence-equation count would
  not flow through; the "3" is baked in (and reinforced by `nwt == 3` special-cases
  in `NKSolvers.F90:2416,2613`).
- `tdia3x3` (block tridiagonal solver, `turbUtils.F90`, added with the model in
  `2909a688`) replaces the scalar `tdia`. It is new numerical infrastructure but
  **serial** — it is not part of the parallel story except that independent lines
  call it concurrently.

### CP-6 — Parallel advection: a fork of shared code, not an extension (deliberate divergence)
`saGRAdvection` (`saGammaRetheta.F90`) was **copied** from the shared
`turbUtils.F90:turbAdvection` and specialized to `mAdv=nAdv=3`, first-order only, so
it could be threaded without editing shared code (CLAUDE.md rule #2). The header says
it plainly: *">>> KEEP IN SYNC <<< … this is the same divergence hazard as the sa.F90
incident."* **Admission:** rather than parametrize/parallelize the shared routine, I
duplicated it — a maintenance liability traded for not touching kw/SST code paths.

### CP-7 — Timing reductions: standard, but the blockette split is an approximation
- `printTurbTiming` uses ordinary `mpi_reduce(max)` + `mpi_reduce(sum)` to rank 0
  (`turbTiming.F90:115–125`) — standard, nothing novel.
- **The fused blockette turb/flow split is not a true measurement.** Because flow and
  turb run in one un-splittable OpenMP sweep, the "turb share" is reconstructed from a
  thread-summed *ratio* `rTurb = TURB_TS/(TURB_TS+FLOW_TS)` applied to the fused
  wall-time (§3c). **Admission:** (a) it is a proportional estimate, not a measured
  wall split; (b) it cannot separate the **three turbulence equations from each
  other** — only turbulence-as-a-whole from flow; (c) before the (uncommitted) fix the
  Context-B turb timers read **0.000** on the blockette path entirely.

### CP-8 — The default-disabled path: my blockette parallel turb residual normally never runs
The `pyADflow.py` guard (§3a) force-disables `useBlockettes` for this model unless
`ADFLOW_ALLOW_SAGR_BLOCKETTES=1`. **Admission:** in the default production config the
blockette SA-GR residual — CP-3's kernels and the fused parallel sweep — **does not
execute at all**; the run falls back to the non-blockette `blockResCore`. So the
blockette parallel work, while present and validated correct, is off by default and
was a net loss on the only mesh measured.

### Connection-point summary (author's honest ledger)

| # | Infra touched | New code? | Honest status |
|---|---------------|-----------|---------------|
| CP-1 | var-count (`nw/nt2`) | 3-line case | inherited generic downstream; trivial |
| CP-2 | MPI halo `whalo2` | **none** | wider `(nt1,nt2)` range only; no new comm |
| CP-3 | block OMP driver | fix + drop-in | repaired THREADPRIVATE set; didn't author driver |
| CP-4 | DADI sweep OMP | **yes** | the real new parallelism; solve itself still serial |
| CP-5 | `qq` 3×3 + `tdia3x3` | yes | **hardcoded 3×3**, not `nwt`; solver serial |
| CP-6 | advection | **fork** | copied shared routine; KEEP-IN-SYNC hazard |
| CP-7 | timing reduce | yes | std reduce; blockette split is a **ratio estimate** |
| CP-8 | blockette enable | guard | **off by default**; parallel turb residual normally idle |

**Net:** the only place I added new parallel infrastructure is the DADI line-sweep
(CP-4) and the repair/extension of the existing block driver (CP-3). The MPI layer I
did not touch (CP-2). The coupling is hardcoded to 3 (CP-5), the threaded advection is
a copy not an extension (CP-6), the blockette timing is an estimate (CP-7), and the
blockette parallel residual is disabled by default (CP-8).

---

## 4. Summary table — baseline vs. this branch

| Aspect | Baseline ADflow | SA-γ-Re̅θt branch | Why it differs |
|--------|-----------------|-------------------|----------------|
| Turbulence vars | 1 (ν̃) | 3 (ν̃, γ, Re̅θt) | transition model adds 2 PDEs |
| Per-cell Jacobian | scalar | dense 3×3 `qq`, `tdia3x3` | block coupling (Eq. 41 γ-weighting) |
| Turb-phase profiling | ANK/NK timers (flow Krylov) | new `turbTiming` A/B split, `-DTURB_TIMING` | DADI line-sweep ≠ Krylov loop; dual-purpose `resOnly` routine |
| DADI threading | serial | OpenMP on assembly + sweeps + new `saGRAdvection` | "free cores" during DADI phase; coupled system heavy enough to thread |
| OpenMP correctness | thread-safe as shipped | needed `THREADPRIVATE`+`private(ii,jj,kk)` fix (`d2b837b6`) | wider per-thread working set exposed shared state |
| Extra comm from coupling | n/a | **none** — wider halo payload, same messages | 3×3 coupling is node-local arithmetic |
| Blockettes default | on (SA) | **off** unless `ADFLOW_ALLOW_SAGR_BLOCKETTES=1` | SA-GR blockette path not fully validated; net loss on 2D |
| Fused-residual timing | not needed | thread-summed `rTurb` ratio splitter | flow+turb fused in one un-splittable OMP sweep |

---

## Appendix — provenance

- Commits quoted: `ed5f2feb` (blockette SA-GR), `d2b837b6` (thread-safety fix),
  `057b000b` (sa.F90 revert), `abefec77`/`8c4e6652` (assembly OMP), `bbff8d42`
  (solve+advection OMP + instrumentation), `03c95173` (2×2 matrix), `6446de8a`
  (doc consolidation). Pre-OpenMP baseline tag: `549719f4`.
- Code: `src/modules/turbTiming.F90`, `src/turbulence/saGammaRetheta.F90`,
  `src/NKSolver/blockette.F90`, `src/turbulence/turbAPI.F90`, `adflow/pyADflow.py`.
- Findings register: `docs/OPENMP_DADI_TIMING_PLAN.md`; run matrix:
  `docs/PARALLEL_CHARACTERIZATION_PLAN.md`.
- **[INFERENCE]** items are marked inline; everything else is quoted from code or
  commit messages. Quantitative speedups (3.49× DADI, 4–5× heavy kernels, 1.36×
  advection, 1.07× total) are from this branch's own measurements recorded in
  `OPENMP_DADI_TIMING_PLAN.md`, all on the thin 2D `volumeMesh_L2.cgns`; none is yet
  confirmed on 3D.
