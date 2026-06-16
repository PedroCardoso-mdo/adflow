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

| # | Figure | Axes | Shows | From |
|---|--------|------|-------|------|
| 1 | Strong-scaling speedup | S vs R, log-log, + ideal line; 2 curves (baseline, branch) | **no regression** — curves overlay | Study 1 |
| 2 | Parallel efficiency | E vs R; baseline vs branch | regression seen as efficiency gap | Study 1 |
| 3 | DADI OMP speedup | S_DADI vs T (log2 T), + ideal; R=4 curve to T=32 with a vertical marker at T=8 (1 NUMA domain) | **OMP scales in-domain**, knee at 8→16 = NUMA-crossing penalty | Study 2 |
| 4 | Phase breakdown (stacked bar) | time per phase, T=1 vs T=max, per rank count | *which* phases scale (Source/Visc fast, Adv/halo plateau) | Study 2 |
| 5 | Total wall vs threads | wall vs T at fixed R | net hybrid gain + Amdahl ceiling (flow serial) | Study 2 |
| 6 | Blockette time-to-solution | grouped bars: wall + iters, blk F/T × {4,32} ranks | blk cost/benefit on 3D | Study 3 |
| 7 | Turb-resid per-call | bars: SA-GR turb-resid ms, blk F vs T, T=1 vs max | blockette effect on the turb kernel itself | Study 3 |
| 8 (opt) | Load imbalance | imb=max/mean vs R (and vs T) | decomposition quality / straggler detection | Studies 1–2 |

Plotting convention: always annotate **cells/rank** and **mesh** on each figure;
report **min (or median) of repeats**, not single shots; ideal lines dashed.

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
