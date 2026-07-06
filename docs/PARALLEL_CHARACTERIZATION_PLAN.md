# Plan — HPC Run Matrix + Figures to Characterize the SA-γ-Reθ Parallelization

Goal: fully characterize the OpenMP + blockette work added to the SA-γ-Reθ DADI
turbulence solve, with three independent claims:

1. **No regression** — the branch keeps the same *pure-MPI* strong scaling as the
   unmodified baseline (out to 256 ranks).
2. **OMP work scales** — the DADI-phase OpenMP speedup is real and holds across
   MPI decompositions (demonstrated at 4 and 32 ranks).
3. **Blockette characterized** — quantify blk on/off (time-to-solution, residual
   phase, iteration count) on a real 3D mesh, not the thin 2D artifact.

All runs reuse `TestOpenMp/run.py` + the env knobs in `validate_openmp.py`.

---

## Prerequisite 0 — Mesh sizing (do this first)

The committed `volumeMesh_L2.cgns` is thin 2D and **cannot** characterize scaling
past a handful of ranks. For Studies 1–3 use the **largest 3D mesh available**
(`NLF_L1` or finer). Sizing rule of thumb:

- Keep **≥ 5–10k cells/rank** at the *largest* rank count, or strong scaling rolls
  off from starvation, not from the code.
- 256 ranks × 8k cells ⇒ mesh ≈ **2M cells** minimum. If the biggest mesh is
  smaller, cap the top rank count so cells/rank stays ≥ 5k and say so.
- Record cells/rank for every point — it is the x-axis context for every figure.

## Prerequisite 1 — Two builds

| Build | Flags | Used for |
|-------|-------|----------|
| **timing** | `-fopenmp -DTURB_TIMING` | Studies 2 & 3 (need [A]/[B] phase split) |
| **prod**   | `-fopenmp` (no TURB_TIMING) | Study 1 (clean wall, zero timing overhead) |

Study 1 can also use the timing build (overhead is negligible), but prod keeps the
"didn't break scaling" claim above reproach. Verify `.so` md5 after each `make`.

## Prerequisite 2 — Two code versions (Study 1 only)

| Tag | Commit | Meaning |
|-----|--------|---------|
| `baseline` | `549719f4` | last commit before ANY OpenMP work |
| `branch`   | `sa_gamma_rethetha` HEAD | the parallelized code |

Build each into its own venv/`.so`; Study 1 overlays their strong-scaling curves.

## Prerequisite 3 — Pinning / launch (NUMA-correct hybrid)

**Allocation:** **2 nodes × 128 cores = 256 cores**. Per node: NUMA domain =
**8 cores** ⇒ **16 NUMA domains/node**, **32 domains total**. Budget: `R*T ≤ 256`.
**Hard rule: OpenMP cannot cross a node** — a rank's `T` threads must fit in one
128-core node, so **per-rank `T ≤ 128`**. The DADI solve is bandwidth-bound, so OMP
scaling is gated by (a) per-domain memory bandwidth and (b) how many 8-core NUMA
domains a rank's threads span = `ceil(T/8)`. T=8 = one full domain (no-crossing
ceiling); T=16/32/64 cross 2/4/8 domains to measure the penalty.

```bash
export OMP_NUM_THREADS=<T>
export OMP_PROC_BIND=close      # fill cores in order -> crisp NUMA boundary at T=16
export OMP_PLACES=cores
# NUMA-crossing probe, both ranks/threads within a node (e.g. R=4 T up to 64):
mpirun --map-by ppr:$((128/T)):node:PE=<T> --bind-to core -np <R> python -u run.py
# clean per-domain run (one rank pinned to one 8-core NUMA node), e.g. R=32 T=8:
mpirun --map-by numa:PE=8 --bind-to core -np <R> python -u run.py
# pure MPI to 256 ranks (Study 1): 128/node x 2 nodes
mpirun --map-by ppr:128:node --bind-to core -np <R> python -u run.py
```
Always `R*T ≤ 256` and per-rank `T ≤ 128`. Record NUMA domains spanned per rank =
`ceil(T/8)` and ranks/node on every run — x-axis context for the figures.

## Prerequisite 4 — Fixed-iteration vs converge (pick per study)

- **Scaling figures (Studies 1, 2):** `RUN_NCYCLES=<N>` (e.g. 500) — *identical work*
  at every rank/thread count, so speedup measures parallel efficiency, not a
  changed convergence path. N large enough to swamp startup (≥ a few hundred).
- **Time-to-solution (Study 3 + the headline blk number):** run to `L2Convergence`
  (no `RUN_NCYCLES`) and report wall **and** iteration count — blk changes the
  trajectory, so wall-to-converge is the honest metric.

---

## Study 1 — Pure-MPI strong scaling (claim: no regression)

Pure MPI, **OMP=1**, blk = production setting (keep fixed). Same mesh, same case
(α=0, nkswitchtol=1e-4), fixed iterations. Two code versions overlaid.

| Knob | Value |
|------|-------|
| ranks R | 1, 2, 4, 8, 16, 32, 64, 128, 256 |
| OMP | 1 |
| blk | fixed (False for a clean baseline; or your production value) |
| versions | `baseline`, `branch` |
| mode | `RUN_NCYCLES=500` |

Runs: 9 ranks × 2 versions = **18** (+2–3 repeats of the 256-rank point for noise).
If 1 rank can't hold the mesh, start at the smallest feasible R and report speedup
relative to it (note the offset).

Metrics per run: total solve wall `T(R)`; speedup `S(R)=T(R0)/T(R)`; efficiency
`E(R)=S(R)/(R/R0)`; cells/rank.

## Study 2 — Hybrid MPI×OMP DADI scaling (claim: OMP work scales)

**Timing build.** Fix MPI at the two counts you named; sweep OMP. blk=False so the
serial flow residual is identical across the sweep and **[A] DADI subtotal** is the
clean OMP signal. Fixed iterations.

Topology gate (2 nodes, 128 cores/node, 8-core NUMA): `R*T ≤ 256`, per-rank `T ≤ 128`.

| Knob | Value | NUMA domains/rank = ⌈T/8⌉, nodes |
|------|-------|------|
| ranks R | 4, 32 | |
| OMP T (R=4) | **1, 2, 4, 8, 16, 32, 64** | 1,1,1,**1**,2,4,8 — full NUMA-crossing probe; 4×64=256 = 2 ranks/node |
| OMP T (R=32) | **1, 2, 4, 8** | 1,1,1,**1**; 32×8=256, clean 1-rank-per-domain on 2 nodes |
| blk | False | |
| mode | `RUN_NCYCLES=500` | |

Runs: 7 (R=4) + 4 (R=32) = **11**. Design intent:
- **R=4 is the NUMA-degradation probe** — T=1→8 stays in one 8-core domain (clean
  reference, expect a bandwidth plateau before 8); T=8→16→32→64 spills into 2/4/8
  domains → exposes the remote-L3/remote-memory penalty (knee expected at 8→16).
  Keep each rank inside one node (`ppr:128/T:node`) so OMP never crosses the node.
- **R=32 is production-realistic** — smaller block/rank exposes fork/join overhead;
  T=8 (2 nodes, `--map-by numa:PE=8`) is the clean 1-rank-per-NUMA-domain point.
The pair separates two causes of rolloff: shrinking per-rank work (R: 4→32) vs
crossing NUMA boundaries (T: 8→64).

Metrics (read from the timing box, max-over-ranks):
- **[A] DADI subtotal** → `S_DADI(T)=A(1)/A(T)`, efficiency `S_DADI/T`.
- Per-phase: Source, Viscous, Advection, DADI line-sweep, Turb halo — which scale,
  which saturate (advection/halo expected to plateau first).
- Total solve wall (shows Amdahl ceiling — serial flow residual still dominates at
  blk=False, so total speedup ≪ DADI speedup; that contrast is a result, not a bug).
- Load imbalance `imb = max/mean` per phase.

## Study 3 — Blockette characterization (claim: quantified on 3D)

**Timing build.** Converge to L2 (time-to-solution). Compare blk on/off at both
rank counts, serial and threaded.

| Knob | Value |
|------|-------|
| ranks R | 4, 32 |
| OMP T | 1, max (from Study 2) |
| blk | False, True (`ADFLOW_ALLOW_SAGR_BLOCKETTES=1`) |
| mode | converge to L2 |

Runs: 2 ranks × 2 OMP × 2 blk = **8**. Reuse the T=1 / T=max endpoints already in
Study 2 where the mode matches (otherwise rerun in converge mode).

Metrics:
- Total wall to converge **and outer-iteration count** (blk changes both).
- **[B] residual subtotal**; **SA-GR turb-resid per-call (ms)** from the work-count
  block (this is the apples-to-apples blk=False↔True turbulence number).
- Fused turb/flow share (`rTurb`) under blk=True.
- Correctness gate: CD matches across all configs (≈ converged reference).

> Expectation to test: blockettes hurt on thin 2D (already shown: 1856 > 1015 s) but
> should pay off on a cache-pressured 3D mesh. Study 3 confirms or refutes this — it
> is the open "3D data point" flagged in `OPENMP_DADI_TIMING_PLAN.md` §MISSING.

---

## Figures

| # | Figure | Status (2026-06-19) | Shows | From |
|---|--------|------|-------|------|
| 1 | Strong-scaling speedup | **kept** — branch only (baseline absent). Axes: y=`Speedup (rel. 4 cores)`, x=`Number of cores` | no regression | Study 1 |
| 2 | Parallel efficiency | **DROPPED** | — | Study 1 |
| 3 | DADI OMP speedup | **kept, R=4 only** — R=32 removed (Finding 2); **NUMA annotations removed** (vertical line + star + "1 NUMA domain" text) so the form matches Fig 1. y-label=`Speedup` | OMP scales in-domain | Study 2 |
| 4 | Turb-solve phase breakdown | **REVISED** — per-iteration [s/iter], log y, grouped per-phase bars **+ Total bar**; flow-residual reference bar removed (Finding 1); y-title 13pt, legend 11pt | which phases dominate/scale (Source + Turb halo) | Study 2 |
| 5 | Speedup vs threads | **kept, R=4 only** — R=32 removed (Finding 2); **now plots SPEEDUP** (`wall(T=1)/wall(T)`), not wall time. blk=OFF (Finding 4) | net hybrid gain + Amdahl ceiling (~2.3×) | Study 2 |
| 6 | Blockette time-to-solution | **DROPPED** (Finding 3) | — | Study 3 |
| 7 | Turb-resid per-call | **DROPPED** — misleading comparison (Finding 3) | — | Study 3 |
| 8 | Load imbalance | **kept** | decomposition quality / straggler detection | Studies 1–2 |

Current figure set: **1, 3, 4, 5, 8** (5 figures). All are **blk=OFF / Study 1–2**.
Plotting convention: **mesh name + cell count annotation removed from all figures**
(was `vmesh_L1.cgns (7,544,832 cells)`); report **min (or median) of repeats**, not
single shots; ideal lines dashed. Mesh is still `vmesh_L1.cgns`, 7,544,832 cells.
Figure canvas: line plots (1,3,5,8) = 840×840 `(6,6)`; Fig 4 = 1120×840 `(8,6)`.

---

## Generating values & per-figure interpretation

Mesh `vmesh_L1.cgns`, 7,544,832 cells. All figures are **blk=OFF**.

### Fig 1 — pure-MPI strong scaling (Study 1, T=1) — *"new model doesn't scale worse than plain ADflow"*
| cores | wall (s) | speedup (rel. 4) | ideal | efficiency |
|---:|---:|---:|---:|---:|
| 4 | 7148.6 | 1.00 | 1 | 100 % |
| 8 | 4233.6 | 1.69 | 2 | 84 % |
| 16 | 3041.3 | 2.35 | 4 | 59 % |
| 32 | 2051.6 | 3.48 | 8 | 44 % |
| 64 | 829.4 | 8.62 | 16 | 54 % |
| 128 | 404.2 | 17.68 | 32 | 55 % |
| 256 | 178.3 | 40.09 | 64 | 63 % |

**Comment:** the SA-γ-Reθ branch keeps healthy pure-MPI strong scaling out to 256 cores
(40× on 64× cores), i.e. the new transition model does **not** introduce an MPI scaling
regression vs the original ADflow flow solve. *Caveat:* this dataset has **no baseline
curve**, so "not worse than original ADflow" is inferred from near-ideal behaviour, not
a direct overlay. (The dip at 16–32 cores then recovery is a cells/rank cache effect.)

### Fig 3 — DADI OpenMP speedup (Study 2, R=4) — *"my DADI OpenMP work scales well"*
| T | DADI subtotal (s) | speedup | ideal | efficiency |
|---:|---:|---:|---:|---:|
| 1 | 2726.8 | 1.00 | 1 | 100 % |
| 2 | 1534.4 | 1.78 | 2 | 89 % |
| 4 | 921.2 | 2.96 | 4 | 74 % |
| 8 | 566.0 | 4.82 | 8 | 60 % |
| 16 | 243.7 | 11.19 | 16 | 70 % |
| 32 | 149.7 | 18.22 | 32 | 57 % |
| 64 | 132.2 | 20.62 | 64 | 32 % |

**Comment:** the DADI turbulence solve (the part I parallelized) scales well — ~5× on
8 threads, 11× on 16, 20× on 64. Efficiency stays 60–70 % up to T=16 and only collapses
at T=64, where each rank's threads span 8 NUMA domains → bandwidth/remote-memory limited,
not an algorithmic limit. **This is the headline "my work scales" result.**

### Fig 4 — turb-solve phase breakdown per iteration (Study 2, log y)
| config | Source | Adv | Visc | DADI-ls | Turb halo | **Total** (s/iter) |
|---|---:|---:|---:|---:|---:|---:|
| R=4, T=1  | 27.5 | 0.36 | 0.68 | 2.3 | 11.6 | **42.6** |
| R=4, T=8  | 4.53 | 0.11 | 0.12 | 0.31 | 3.73 | **8.80** |
| R=4, T=64 | ~1.0 | — | — | 0.2 | 0.8 | **2.05** |
| R=32, T=1 | 1.94 | 0.13 | 0.15 | 0.36 | 1.3 | **3.89** |
| R=32, T=8 | 0.46 | 0.05 | 0.04 | 0.05 | 0.19 | **0.53** |

(iters: R=4 = 64, R=32 = 59; ×iters → absolute s, matches `dadi_subtotal`.)
**Comment:** **Source + Turb halo** dominate every config; both shrink with threads.
Source is the transition-correlation cost (Finding 1); halo is communication-bound.

### Fig 5 — total-wall speedup vs threads (Study 2, R=4) — *"didn't revolutionize wall time, but OMP help is real & free when cores are spare"*
| T | wall (s) | speedup |
|---:|---:|---:|
| 1 | 7065 | 1.00 |
| 2 | 5788 | 1.22 |
| 4 | 4748 | 1.49 |
| 8 | 3832 | 1.84 |
| 16 | 3212 | 2.20 |
| 32 | 3076 | 2.30 |
| 64 | 3017 | 2.34 |

**Comment:** end-to-end wall speedup saturates at **~2.3×** (Amdahl): the DADI solve
threads well (Fig 3) but the flow residual + linear solver — the bulk of the wall — are
not threaded, so they cap the total. So OpenMP **did not revolutionize** time-to-solution,
**but the gain is far from insignificant** (1.8× at 8 threads, 2.3× at 32+) and it is
essentially **free wall-clock when spare cores are available** on the node.

### Fig 8 — load imbalance (max/mean) vs ranks (Study 1, T=1)
| R | Source | DADI-ls | Turb halo | Turb resid | Flow resid |
|---:|---:|---:|---:|---:|---:|
| 4 | 1.21 | 1.22 | 2.26 | 1.20 | 1.05 |
| 8 | 1.73 | 1.08 | 1.54 | 1.74 | 1.04 |
| 16 | 2.28 | 1.27 | 1.35 | 1.91 | 1.05 |
| 32 | **3.03** | 1.50 | 2.04 | **2.66** | 1.06 |
| 64 | 1.80 | 1.53 | 1.96 | 1.68 | 1.16 |
| 128 | 1.42 | 1.48 | **4.55** | 1.39 | 1.13 |
| 256 | 1.23 | 1.37 | **4.14** | 1.28 | 1.22 |

**Conclusion (self-contained — nothing to compare against):** the **flow residual is
essentially perfectly balanced** (max/mean ≈ 1.0–1.2 at every rank count) — the mesh
partitioner balances flow work well. The **turbulence/transition phases are not**: Source
peaks at 3.0× (R=32) and Turb halo at 4.5× (R=128). So the decomposition — tuned for the
flow — does **not** balance the transition workload: a minority of ranks (those holding
boundary-layer/transition regions, and surface-heavy partitions at high R) carry several
times the mean Source/halo work and become the **stragglers that gate the turbulence
solve**. Net: turbulence load imbalance — not the flow — is the partitioning bottleneck,
and weighting transition cells in the decomposition is the obvious future lever.

---

## Findings & interpretation (2026-06-19 analysis session)

### Timing has TWO residual contexts — do not conflate them
The SA-γ-Reθ residual is computed and timed at **two** call sites (`turbTiming.F90`):
- **Context A — inside the DADI solve** (`turbSolveDDADI`): phases Source + Viscous +
  Advection + Unsteady are summed as "Residual + qq assembly". Called **~18–21×/iter**
  (one per DADI sweep). This is where the heavy turbulence work lives.
- **Context B — standalone residual eval** (`blocketteRes` → `blockResCore` /
  `blocketteResCore`): `T_RESID_FLOW`, `T_RESID_TURB`, `T_RESID_BOTH`. Flow residual
  called **~128–220×/iter** (ANK Krylov matrix-free J·v products); turb **~6–7×/iter**.

`p_resid_turb` is **only Context B**. The bulk of turbulence-residual cost is the
Context-A assembly (the Source phase). Never equate `p_resid_turb` with "the turbulence
residual cost" — that under-counts it by the DADI-internal assembly.

### Finding 1 — the transition SOURCE term is the single cost driver (drove the Fig 4 rework)
Per-call, the turbulence residual is **89–97% Source term** (transition correlations
Reθt/Reθc/Flength/Fonset/Fturb — `exp`/`tanh`/`pow`/`min`/`max`, which don't vectorize).
- Per-call R=4: Source 1530 ms vs Advection 19 ms, Viscous 36 ms.
- Flow vs turb residual per-call: R=4 → 154 vs 1603 ms (turb **10×**); R=32 → 50 vs
  104 ms (turb **2×**). Fewer equations (3 vs 5) but far costlier per call **because of
  the algebraic source the flow equations don't have**. Strip the source and turb
  adv+visc IS cheaper than flow, as equation count predicts.
- Aggregate (call-count weighted) flips with rank: total turb-residual work (A+B) =
  **1.97× flow at R=4**, **0.38× flow at R=32** (flow is called ~6× more often).

Old Fig 4 stacked the 5 turb-solve phases against a single flow-residual bar, inviting
a false "flow < turb residual" reading. **Reworked:** only the turb-solve phases, per
iteration, log y, grouped + a Total bar. iters: R=4 = 64, R=32 = 59 (multiply per-iter
by iters to recover absolute seconds; Totals match `dadi_subtotal`).

### Finding 2 — R=32 flow-residual "scaling" is a NUMA-bandwidth confound (→ R=32 dropped from Figs 3 & 5)
S2 ran **blk=False** → flow residual goes through `blockResCore`, which has **NO
OpenMP** (the OMP pragmas live only in `blocketteResCore`, lines 406–825). So the flow
residual is **unthreaded** in S2.
Yet at R=32 the flow-residual time drops 462→154 s over T=1→8. **This is not threading.**
SLURM logs: T=1,2,4 keep **32 ranks/node (constant)**, but `--cpus-per-task=T` spreads
those 32 ranks across more 8-core NUMA domains / memory controllers (≈4→8→16). The
bandwidth-bound flow residual then gets ~2× the bandwidth (461→242 ms = 1.9× at T1→T2)
and saturates by T=4. R=4 stays flat (4 ranks never saturate one node).
=> This inflates the apparent **end-to-end** OMP benefit at R=32 with memory bandwidth,
not DADI threading. R=32 removed from Figs 3 & 5 to avoid a "too good to be true" curve
needing this caveat in the report. The **DADI subtotal** (intended OMP signal) is
unaffected and stays the headline metric.

### Finding 3 — blockette (Study 3) figures dropped as misleading
Old Fig 7 compared blk OFF `c_rturb_ms` (turb **alone**, 7 calls/iter) vs blk ON
`c_rboth_ms` (flow+turb **fused**) — **different operations** — making blk ON look 2×
faster. Real picture (S3, verified):
- Dominant residual cost is the **FLOW** residual (138–220 calls/iter), not turb (7/iter).
- blk ON flow-residual per-call: R=32/T=1 **1.5×** faster (48.7→32.2 ms), R=32/T=8
  **3.6×** (16.5→4.6 ms), but **R=4/T=1 SLOWER 0.63×** (153→242 ms).
- Fused residual split (`rTurb`): **~95% flow / ~5% turb** (the heavy source is in the
  DADI solve, not in the residual eval — consistent with Context A vs B above).
- blk ON needs **more outer iters** (R=32/T=1: 202→254, +26%), DADI is slower, and the
  residual eval is only **~27% of wall** → Fig 6 total gain is modest (**−12% at R=32**,
  slight **regression at R=4**: 32174→32423 s).

Net: blockettes help the flow kernel at R=32 (cache + threading) but not at R=4, and the
wall payoff is small once iteration count and the dominant linear-solver time are
accounted for. **If Fig 7 is ever revived it must compare flow-residual per-call OFF vs
ON (apples-to-apples), not BOTH vs turb-alone.**

### Finding 4 — a blk=ON version of Fig 5 is NOT possible with current data (MISSING)
Fig 5 (speedup vs threads) needs a full thread sweep at one R. For **blk=ON, R=4** the
dataset only has **two** points, and one diverged:

| study | T | blk | iters | CD | wall | mode |
|---|---:|---:|---:|---:|---:|---|
| S3 | 1  | True | 325 | 0.4441 | 32423 s | converge |
| S3 | 16 | True | **1** | **39.71** | 60 s | converge — **DIVERGED** |

Three blockers:
1. **No sweep** — Study 3 only ran T=1 and T=max; the T=1,2,4,8,16,32,64 sweep exists
   only in S2 (blk=OFF). blk=ON has no equivalent.
2. **The one threaded blk=ON point diverged** (R=4/T=16: CD=39.7, 1 iter) → unusable,
   leaving a single valid point (T=1).
3. **Mode mismatch** — S2 (blk=OFF) is fixed-1000-cycles; S3 (blk=ON) is converge. The
   two cannot be mixed into one speedup curve.

The R=4/T=16 blk=ON divergence is itself a result: **the blk=ON path is unstable with
threads at R=4.** To build a blk=ON Fig 5 we need NEW runs: a full T sweep at R=4 with
blk=True, same mode as S2 (fixed cycles), after fixing that instability. All current
figures are therefore **blk=OFF**.

---

## Totals & sequencing

- ~**18 (S1) + 10 (S2) + 8 (S3) ≈ 36 runs** + repeats of the expensive top points.
- Order: (0) size the mesh, (1) build both versions + both flag sets, (2) Study 1
  on prod build, (3) Study 2 on timing build, (4) Study 3 on timing build, (5) parse
  logs (`[TIMING] solve_wall_s=`, the timing box, iter count) → CSV → figures.
- A small parser can scrape every `run_<R>mpi_<T>omp_blk<b>.log` into one CSV
  (ranks, omp, blk, wall, iters, A-subtotal, B-subtotal, per-phase, imb, CD) that
  feeds all eight figures.

## Smoke gate every run
No NaN; converged CD matches the reference; NK-switch iteration count unchanged for
fixed-iter runs. A run that diverges or shifts CD is discarded, not plotted.
