# Plan — Instrument + Aggressively Parallelize the SA-γ-Reθ DADI Turbulence Solve

> **On approval, step 0:** save this plan to `docs/OPENMP_DADI_TIMING_PLAN.md`
> for later reference.

---

## STATUS (updated 2026-06-09)

### ✅ DONE
- **Phase A — timing instrumentation:** `src/modules/turbTiming.F90` + Context-A
  hooks (turbAPI, saGammaRetheta) + Context-B hooks (blockette) + reset/print in
  solvers. All gated behind `#ifdef TURB_TIMING`. Timing box prints correctly.
- **Phase B — parallelization:** Source/Viscous `collapse(3)`; Viscous 3→1
  region; `saGammaReThetaSolve` qq-prep + j/i/k sweeps with private tridiagonal
  buffers; new SA-GR-local `saGRAdvection`. 55 OMP pragma lines, all in primal
  `#else` branches → no Tapenade regen.
- **Phase C — measurement: FULL 2×2 matrix complete + correctness gate PASSED.**
  - Determinism: omp 1 vs 12 @ 200 cyc, blk=False → **bit-identical CD 0.016210**
    (proves race-free/deterministic parallelization).
  - DADI-phase OpenMP speedup @ 12 threads (200 cyc, blk=False): **3.49×**
    subtotal (Source 4.84×, Viscous 4.15×, DADI-solve 3.99×, SA-GR resid 4.39×,
    Advection 1.36× [mesh-limited — light/memory-bound + 3 barriers on thin 2D]).
  - Full-convergence 2×2 matrix (RUN_NCYCLES=25000; all reach CD≈0.007009, no NaN):

    | # | blk | omp | wall (s) | CD       |
    |---|-----|----:|---------:|----------|
    | 1 | F   | 1   | 1082     | 0.007009 |
    | 2 | F   | 12  | **1015** | 0.007009 |
    | 3 | T   | 1   | 3798     | 0.007010 |
    | 4 | T   | 12  | 1856     | 0.007010 |

  - **CONCLUSION: blk=False + 12 threads (1015s) is the best config on this mesh.**
    Blockettes HURT here even with threads (1856 > 1015s): the blockette SA-GR
    trajectory needs more outer iters (~6700 vs ~5200) and cache-blocking overhead
    isn't amortized on a small 2D case. Threads do help within blk=True
    (3798→1856 = 2.05×), just not enough to beat blk=False. Blockettes are
    expected to pay off only on large cache-pressured 3D meshes.
  - blk=True is CORRECT (CD=0.007010 ≈ blk=False 0.007009) → the pyADflow guard
    that force-disables SA-GR blockettes is conservative, not because the path is
    wrong. It can be relaxed, but blk=True is a perf trap without enough threads
    AND on small meshes.
- **Prod build WITHOUT `-DTURB_TIMING` verified:** compiles clean + imports;
  timing call-sites compile out (turbTiming module object still links but is never
  called → zero runtime overhead). Build deterministic (md5 stable).
- Build clean with `-fopenmp -DTURB_TIMING`; venv reinstalled (md5 verified).
- **Committed** validated work as `bbff8d42` on top of `8c4e6652` (revert point).

### ❌ MISSING / TODO (out of this task's scope)
- **3D mesh data point** — advection 1.36× and the blockette verdict are both
  thin-2D-limited; neither is confirmed on a real cache-pressured 3D mesh, where
  blockettes (and advection threading) are expected to actually pay off.
- **Context-B timers don't instrument the blockette path** — `T_RESID_FLOW`/
  `T_RESID_TURB` read 0.000 under blk=True (they're hooked on the non-blockette
  `blockResCore`, not `blocketteResCore`). Measuring per-phase blockette speedup
  would need new hooks in the blockette residual routines.
- **GOTCHA recorded:** `vars_master.py nCycles=1000` caps a plain run at
  itertot~1000 (~iter 110, NOT converged). Must set `RUN_NCYCLES=25000` to
  actually converge.
- **Physics correctness** (CL/CD vs experiment) is still the user's end check —
  these runs only verify self-consistency (all configs agree on CD).

## Rollback points (safe versions if everything goes badly)

- **Pre-OpenMP baseline (the version to `git reset`/`checkout` to if it all goes
  wrong):** `549719f4` — *"Replace FD with analytical derivative for qq(2,3)"*
  (2026-05-29). This is the last commit before ANY OpenMP work. To return:
  `git checkout 549719f4` (or `git reset --hard 549719f4` to discard everything after).
- **First OpenMP commit (where the OpenMP work began):** `d2b837b6` — *"Fix
  OpenMP thread-safety for blockette SA and SA-gamma-rethetha"* (2026-06-02).
- OpenMP commit chain after that: `abefec77` (DADI residual/Jacobian assembly),
  `057b000b` (revert sa.F90 constants — keep SA untouched), `8c4e6652` (DADI
  assembly OpenMP + CLAUDE.md/.gitignore). Current timing work (this plan) is
  built on top of `8c4e6652`.

## Context

OpenMP on the SA-γ-Reθ DADI assembly (`Source`/`Viscous`) is correct (1×1 vs
1×12 converge to identical CD=0.007009) but gives only **1.07× total speedup**.
Two questions drove this plan:

1. **Is adjoint code secretly running in the primal solve?** — **NO.** Traced
   the full primal ANK+DADI path: `solveState → ANKStep → (FormFunction_mf →
   blocketteRes → blocketteResCore) + turbSolveDDADI → saGammaRetheta_block`.
   Every routine is the hand-written primal version. **No `src/adjoint/output*`
   (`_d`/`_b`) code executes during a primal analysis solve** (those only run in
   adjoint drivers via matrix-free Fwd/Bwd products). So that's not the cause.

2. **Why so little speedup?** Amdahl: with `useBlockettes=False` the dominant
   cost — the **flow** residual evaluated many times in ANK's matrix-free Krylov
   (`blockResCore`) — is **serial**, and the DADI turbulence is a small add-on of
   which only `Source`+`Viscous` were parallel. The DADI **solve**, **advection**,
   and **unsteady** terms are still serial.

**User's framing (valid):** in production they run blockettes+OpenMP (flow
parallel) + DADI turbulence. During the DADI phase the rank's cores are
otherwise idle, so **any** DADI-phase speedup is free. Goal: (a) **measure** the
DADI phase in isolation, then (b) **maximize** OpenMP across the whole DADI path
(assembly + solve + advection), minding fork/join overhead on thin cases.

---

## Part 1 — Timing instrumentation (FIRST; measure serial baseline before any new OpenMP)

Goal: see each phase's cost **before** parallelizing — especially **advection on
its own** — so we know the real ROI and ceiling. Mirror the `tic/toc` +
`mpi_wtime()` pattern from `src/overset/oversetUtilities.F90:5-20`.

**New module** `src/modules/turbTiming.F90` — `mpi_wtime()` accumulators
(per rank), in **TWO contexts** so DADI parts and the blockette residual are
differentiated, and it works with blockettes **on or off**:

**Context A — DADI solve** (`turbSolveDDADI` → `saGammaRetheta_block`,
`resOnly=.false.`). Four assembly contributors, each its own line; **their sum =
residual+Jacobian (scratch + qq)** (per user):
- `tSource` — `Source()`            ┐
- `tViscous` — `Viscous()`          │ sum = "residual + qq assembly"
- `tAdv` — `turbAdvection`          │
- `tUnsteady` — `unsteadyTurbTerm`  ┘
- `tDADI` — `saGammaReThetaSolve` (the line-sweep solve)
- `tResScale` — `ResScale` (small)
- `tTurbHalo` — `whalo2` (turbulence halo in turbSolveDDADI)
- `tTurbTotal` — whole `turbSolveDDADI`

**Context B — residual evaluation** (`blocketteRes`, called from ANK
`FormFunction_mf` / MG; `resOnly=.true.`). This is where blockettes matters
(blk=True → parallel blockette flow+turb residual; blk=False → serial block).
Split so the SA-GR turbulence residual is isolated from the flow residual:
- `tResidFlow` — flow residual (inviscid central + dissipation + viscous fluxes)
- `tResidTurb` — **SA-GR turbulence residual**: blk=True → blockette
  `saGammaRethetaSource/Advection/Viscous/ResScale`; blk=False →
  `saGammaRetheta_block(.true.)`. Timed at that call site either way.
- `tResidHalo` — halo exchange in/after the residual path
- `tResidTotal` — whole `blocketteRes` call

This lets you read **blockette speedup on YOUR model** directly from `tResidTurb`
(and `tResidFlow`) across blk=False↔True and OMP 1↔12.

**Separation mechanism:** the fine Context-A timers in `saGammaRetheta_block`
fire **only when `resOnly=.false.`** (the actual DADI solve), so the
`resOnly=.true.` residual-eval calls don't pollute them; the residual-eval cost
is captured by Context B at the `blocketteRes` level. Result: clean split of
"DADI work" vs "blockette/residual work + halos", in both blk modes.

Helpers `turbTic(id)/turbToc(id)` + `printTurbTiming`. Timers are `mpi_wtime()`
at phase boundaries (outside parallel regions) → **negligible overhead, no effect
on blockette behavior or results**.

### Easy to remove / disable (production-friendly) — REQUIRED
Gate **everything** behind a single CPP macro so production builds carry **zero
overhead and zero residue**, and removal is trivial:
```fortran
#ifdef TURB_TIMING
    call turbToc(T_ADV)
#endif
```
- Disabled by default; enabled by adding `-DTURB_TIMING` to `config.mk` (one line)
  for our measurement build.
- When the macro is undefined, the module compiles to empty stubs (or is not
  compiled) and all call sites vanish → nothing in the production binary.
- Removal for good = delete `turbTiming.F90` + grep-delete the `#ifdef TURB_TIMING`
  blocks (all localized at phase boundaries, ~8 sites).

**Hook points (each wrapped in `#ifdef TURB_TIMING`):**
- `turbAPI.F90 turbSolveDDADI`: `tTurbTotal` (around `nSubIterTurb` loop) and
  `tTurbHalo` (the `whalo2` at line 103).
- `saGammaRetheta.F90 saGammaReTheta_block`: separate tic/toc around `Source`,
  `Viscous`, `turbAdvection`, `unsteadyTurbTerm`, `ResScale`, `saGammaReThetaSolve`
  — **only inside the `if (.not. resOnly)` path / guarded by `resOnly==.false.`**
  so they measure the DADI solve only (Context A).
- `NKSolver/blockette.F90 blocketteRes/blocketteResCore/blockResCore`: wrap
  `tResidTotal` (whole call), `tResidFlow` (flux routines), `tResidTurb` (the
  SA-GR residual call site — blockette `saGammaRetheta*` if blk=True, else
  `saGammaRetheta_block(.true.)`), and `tResidHalo` — Context B. Header prints
  the active path (blockette vs block) from `useBlockettes`.
- Reset at `t0Solver = mpi_wtime()` (`solvers.F90:52`); print once at end via the
  Solution Timings block (`solvers.F90` ~1623+), rank 0, `MPI_Reduce` max for comm.

**Phase-A deliverable:** build with `-DTURB_TIMING` on the *current* code
(assembly already OMP, advection+DADI-solve still serial), run 1 thread → get the
**serial per-phase breakdown** (incl. advection %) — this is the baseline the user
wants before adding OpenMP to advection/solve.

This yields the **isolated DADI-phase speedup** (resjac+dadi+comm) at 1 vs 12
threads — the number the "free cores" decision actually needs.

**Output: printed ONCE at the end of the solve** (not per iteration), rank 0,
in/after the existing Solution Timings block. Example format:

```
+--------------------------------------------------------------+
|  ADflow SA-gamma-Retheta timing   ranks=1  OMP=12  BLOCKETTES=ON
|--------------------------------------------------------------|
|  [A] DADI turbulence solve (turbSolveDDADI)                  |
|    Source                      :      6.420 s                |
|    Viscous                     :      5.760 s                |
|    Advection                   :      7.150 s                |
|    Unsteady term               :      0.310 s                |
|      = Residual + qq assembly   :     19.640 s   (sum of 4)  |
|    DADI line-sweep solve        :     22.910 s                |
|    ResScale + other             :      0.310 s                |
|    Turb halo (whalo2)           :      6.050 s                |
|    DADI subtotal                :     48.910 s                |
|                                                              |
|  [B] Residual eval  (path: BLOCKETTE)                        |
|    Flow residual                :    240.100 s                |
|    SA-GR turbulence residual    :     70.100 s                |
|    Residual halo                :     22.400 s                |
|    Residual subtotal            :    332.600 s                |
|--------------------------------------------------------------|
|  (header shows config; comm lines are max-over-ranks)        |
+--------------------------------------------------------------+
```

Header self-labels `ranks / OMP / BLOCKETTES` so each run's output is
unambiguous. Lay two runs side by side to read speedups:
- **DADI OpenMP:** [A] subtotal, blk=False, OMP 1 vs 12.
- **Blockette speedup on your model:** [B] `SA-GR turbulence residual` (and Flow),
  blk=False vs blk=True at the same OMP; and OMP 1 vs 12 within blk=True.

---

## Part 2 — Parallelize the DADI **solve** (biggest untapped win)

`saGammaReThetaSolve` is currently 100% serial but is parallel **across lines**:
- j-sweep `do k; do i; do j(inner)` (line ~1668): outer `(k,i)` lines independent
- i-sweep: outer `(k,j)` independent
- k-sweep: outer `(i,j)` independent
- the qq prep/scaling loop (lines ~1603–1660): per-cell independent

**Changes (all in `saGammaRetheta.F90:saGammaReThetaSolve`):**
- qq prep loop → `!$OMP parallel do collapse(3) private(i,j,k,...)`.
- Each sweep → `!$OMP parallel do collapse(2)` over the **outer two** loops,
  inner solve-direction loop stays serial. **Work arrays `bb, cc, dd, ff` and all
  per-line scalars must be `private`** (per-thread tridiagonal buffers) — this is
  the key correctness requirement (`tdia3x3` writes into them).
- `private` also: `i,j,k`, the recomputed diffusion/advection coefficients, and
  the local geometry temporaries in each sweep.
- Keep `!$OMP end parallel do` after each.

Race-safety: within a sweep each thread owns a distinct line and writes only that
line's `scratch`; sweeps run sequentially (separate regions), preserving the
ADI ordering. `qq`/`scratch` reads are of the thread's own line.

---

## Part 3 — Reduce assembly fork/join overhead

- `Viscous()`: the 3 directional loops currently fork/join 3×. Wrap them in a
  single `!$OMP parallel` with three `!$OMP do` (orphaned) → **3→1** fork/joins.
  (`Source()` already one region.)
- Net assembly fork/joins: 4 → 2.

---

## Part 4 — Parallel advection (SA-GR-local; YES per user)

`turbAdvection`/`unsteadyTurbTerm` live in shared `turbUtils.F90` (used by
kw/SST) — **do not modify** (CLAUDE.md rule #2). Add an **SA-GR-local** parallel
advection routine inside `saGammaRetheta.F90` (specialized to `mAdv=nAdv=3`,
first-order upwind per paper §IV.A) and call it from `saGammaReTheta_block`
instead of the generic `turbAdvection`; parallelize its cell loop. **Sync risk:**
this duplicates ~part of the shared advection — same divergence hazard as the
`sa.F90` incident, so document the source it was derived from.

The user wants this done, but **measure its serial cost first** (Part 1) so the
gain is quantified. Do the parallelization in Phase B alongside the solve.

---

## Part 5 — Overhead awareness (no size guard)

The test case is thin/2D, so OpenMP overhead is visible here — but the goal is
larger 3D meshes where it pays off, so OpenMP stays **always on** (no cell-count
`if()` guard — removed per user). Keep overhead low via:
- `schedule(static)`.
- Avoid `collapse` where the outer loop alone has enough iterations (collapse
  adds index arithmetic).
- Reduce fork/joins (Part 3).

---

## Constraints / non-negotiables

- **Never touch SA-only or shared code** (`sa.F90`, `turbUtils.F90` shared
  routines) — see CLAUDE.md rule #2 / memory `never-touch-sa-model-code`.
- Pragmas only in the **primal** `#else` branches; `#ifdef TAPENADE_REVERSE`
  paths untouched → **no Tapenade regeneration needed** (math unchanged).
- All edits confined to: `saGammaRetheta.F90`, `turbAPI.F90`, new
  `turbTiming.F90`, small hooks in `solvers.F90`. (`saGammaReThetaSolve` is NOT
  Tapenade-differentiated for the solve itself; confirm before relying on that.)

---

## Files to modify

| File | Change |
|------|--------|
| `src/modules/turbTiming.F90` (new) | timing accumulators + tic/toc + print |
| `src/turbulence/turbAPI.F90` | wrap total + `whalo2` comm timers |
| `src/turbulence/saGammaRetheta.F90` | resjac/DADI timers; parallelize solve sweeps + qq prep; Viscous 3→1 region; SA-GR-local parallel advection; overhead guards |
| `src/solver/solvers.F90` | reset + print turb timing in Solution Timings block |

---

## Verification — 1 baseline + 2×2 matrix + build check (1 rank, nkswitchtol=1e-4, L2=1e-6, α=0)

**Phase A (before new OpenMP) — 1 run:** serial breakdown, the "times before OpenMP".

| # | blk | OMP | Why |
|---|-----|-----|-----|
| A0 | False | 1 | Serial per-phase baseline (DADI [A] + residual [B] breakdown) before parallelizing the DADI solve/advection. Go/no-go + the "before" picture. |

**Phase C (after parallelizing) — 2×2 (blk × OMP):** characterize DADI and
blockette of the SA-GR model, with/without OpenMP, with/without blockettes.

| # | blk | OMP | What it isolates |
|---|-----|-----|------------------|
| 1 | False | 1 | Reference (serial DADI + serial block residual). Correctness ref (CL/CD). |
| 2 | False | 12 | **DADI OpenMP speedup** — [A] only ([B] serial both ⇒ clean). |
| 3 | True | 1 | **Blockette restructuring effect** (no threads): [B] block→blockette. |
| 4 | True | 12 | **Full production** — blockette+OMP residual [B] + DADI+OMP [A]. |

Reading the result:
- DADI OpenMP: [A] subtotal #2 vs #1.
- Blockette OpenMP on SA-GR residual: [B] `SA-GR turb residual` #4 vs #3.
- Blockette vs block (no threads): [B] #3 vs #1.
- Full speedup: total #4 vs #1.
- Correctness gate: CL/CD match across all four (CD≈0.007009).

`blk=True` needs `ADFLOW_ALLOW_SAGR_BLOCKETTES=1` (guard bypass; gated edit).

| build | Why |
|-------|-----|
| prod build, no `-DTURB_TIMING` | Confirm timing fully removable (compiles, zero residue). |

Each run: `make` → `pip install .` → confirm `.so` md5 → run. (`./AD_I.sh` not
required — pragmas are comments — run only if user asks.) **Dropped:** MPI hybrid
combos (2×6/4×3/12×1 — scaling nice-to-have only). Smoke every run: no NaN,
NK-switch iteration count unchanged.

---

## Phasing

- **Phase A — timing only (no new OpenMP):** Part 1 behind `-DTURB_TIMING`.
  Build, run 1 thread → **serial per-phase baseline incl. advection %**. This is
  the "times before OpenMP" the user asked for.
- **Phase B — parallelize:** Part 2 (DADI solve sweeps) + Part 3 (fork/join
  reduction) + Part 4 (SA-GR-local parallel advection) + Part 5 (low-overhead).
- **Phase C — re-measure:** same timers at 1/12 threads + hybrid combos →
  quantify per-phase gains; correctness gate (CL/CD match).
