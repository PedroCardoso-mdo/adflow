# ADFLOW 02 — Discrete Adjoint & Differentiation (complete extraction: Kenway2019)

**Source:** Kenway, Mader, He, Martins, *Effective adjoint approaches for computational fluid dynamics*, Prog. Aerosp. Sci. (2019), doi:10.1016/j.paerosci.2019.05.002.
**Coverage:** comprehensive technical extraction of §2–§6 + Appendix A. Abstract + narrative dropped; hand facts retained; prose paraphrased; equations and benchmark data reproduced as facts. Eq/§/Table/Fig refs are the paper's. DAFoam (OpenFOAM comparison solver) kept where it clarifies ADflow's choices.

---

## 0. Context facts worth keeping (§1–§2)

- **Adjoint payoff:** cost **independent of #design variables**, scales with **#functions of interest** (objectives+constraints, typically <10 vs hundreds of DVs). Adjoint + gradient-based optimization is the workhorse for 3-D aero shape optimization.
- **Continuous vs discrete:** ADflow uses **discrete** (discretize then linearize). Continuous = low-memory, duality-preserving, but **inaccurate on coarse meshes** (consistent only in the infinite-mesh limit), needs hand-differentiation of the PDE, and the viscous **shear-stress adjoint BC is ill-posed**. Discrete = gradients **numerically consistent with the flow** on any mesh; partials by FD/CS/AD (no hand differentiation needed); can be **machine-exact even through turbulence models**. Discrete's cost = building/storing the exact Jacobian in some form → efficient partials are the whole game.
- **AD reverse mode is analogous to the adjoint** (cost independent of #inputs). **Full-code AD** (differentiate the entire solver) is inefficient / memory-prohibitive on large 3-D; must store all intermediates in the forward pass. **Selective AD** (differentiate only residual/objective routines) is the practical path — prior art: Christianson (reverse accumulation), Griewank–Faure (piggyback), Giles (primal-dual), Xu (Jacobian-trained Krylov-implicit RK, STAMPS), Albring (duality-preserving, SU2), **Mader (selective AD + GMRES, complex-step verified — ADflow's predecessor)**.
- **This paper:** implement **3 adjoint variants** in ADflow + DAFoam, evaluate on **6 factors** (speed, accuracy, memory, scalability, dev effort, extensibility), vs a full-code-AD reference. Propose public **ADODG benchmarks** (Case 3 rect wing, Case 4 CRM wing, Case 5 wing-body-tail).
- **Unsteady/chaotic caveat:** unsteady adjoint costs ≥1 order more (store/recompute linearizations per step); chaotic turbulence makes adjoint derivatives diverge (perturbation growth) — largely unsolved for practical LES/DNS.

## 1. Flow modeling (§3) — ADflow relevant

- **NS `(1)–(3)`:** continuity, momentum, energy. SA `(4)`: `∇·(Uν̃) − (1/σ)∇·[(ν+ν̃)∇ν̃] + C_b2|∇ν̃|² − C_b1 S̃ν̃ + C_w1 f_w (ν̃/d)² = 0`; `ν_t = ν̃·χ³/(χ³+C_v1³)`, `χ = ν̃/ν` `(5)`. Four terms = convective / diffusion / production / near-wall destruction.
- **ADflow:** FV structured multiblock + overset (implicit hole cutting, Kenway [137]); SUmb-based; Fortran90 + Python. Inviscid schemes: **JST scalar** [138], **matrix dissipation** (Turkel–Vatsa [139]), **MUSCL upwind** (van Leer [140] + Roe [141]). Green–Gauss viscous. Compressible → mean-flow coupled. Turbulence: **SA** + **Menter SST**; segregated or fully coupled. Four solvers: **RK, D3ADI, ANK (Yildirim [144] = companion paper), NK [35]**. RK/D3ADI MG-compatible but weak globalization without MG; **ANK is the default globalizer for overset** (robust without MG); NK gives terminal convergence once in the Newton basin. (Enabled the circle→airfoil optimization [145].)
- **DAFoam** (comparison): OpenFOAM `simpleFoam`, **incompressible** (ignore energy), **SIMPLE** [146] + **Rhie–Chow** [147], **segregated**, 2nd-order upwind inviscid + central viscous, Gauss–Seidel + AMG; C++.

## 2. Adjoint equations (§4.1)

- Residuals `R(x,w)=0` `(6)`, `f=f(x,w)` `(7)` (cheap, no iteration). **`w` and `R` include the turbulence variable — no frozen turbulence** (consistency).
- Chain rule `(8)`: `df/dx = ∂f/∂x + (∂f/∂w)(dw/dx)`. `∂f/∂x, ∂f/∂w` cheap; `dw/dx` expensive.
- `dR/dx=0` `(9)` → linear system `(10)`: `(∂R/∂w)(dw/dx) = −∂R/∂x`.
  - **Direct:** solve `(10)` once per DV → cost ∝ #DVs.
  - **Adjoint `(12)`:** `(∂R/∂w)ᵀ ψ = (∂f/∂w)ᵀ` (solve once per function) `(11)`; total `(13)`: `df/dx = ∂f/∂x − ψᵀ(∂R/∂x)`.
- **Four steps:** (1) build `[∂R/∂w]ᵀ`, `[∂f/∂w]ᵀ`; (2) solve `(12)` for `ψ`; (3) build `∂R/∂x`, `∂f/∂x`; (4) assemble `(13)`.

## 3. Partial-derivative computation (§4.2)

Five options; trade-offs:

| Method | Accuracy | Cost / effort | Key facts |
|---|---|---|---|
| **Analytic** (hand) `§4.2.1` | exact | very high effort, error-prone | needs deep solver expertise; usually simplifies turbulence (frozen) → inaccurate; largely abandoned. |
| **Finite difference** `(14) §4.2.2` | truncation + **subtractive cancellation**; step-size sensitive | easy; black-box OK | `(∂y/∂x)_{n,m} ≈ [y_n(x+εe_m) − y_n(x)]/ε + O(ε)`; needs step-size study. |
| **Complex step** `(15) §4.2.3` | machine precision (no cancellation) | 2–4× slower, 2× memory | `(∂y/∂x)_{n,m} = Im[y_n(x+εi e_m)]/ε + O(ε²)`, ε~1e-40; must complexify all floats + logical/abs ops; **no transpose product** → unusable for Jacobian-free GMRES; used for **reference/verification** only. |
| **AD** `(16)–(21) §4.2.4` | machine precision (exact) | tool-dependent | forward + reverse modes; ADflow's choice. |
| **Symbolic** `§4.2.5` | exact | closed-form only; expression blows up | not used for iterative CFD; FEM semi-symbolic (dolfin-adjoint). |

**AD modes:**
- **Forward `(16)–(18)`:** seed input `ẋ = e_n` → `ẏ = [∂y/∂x]ẋ` (one column / regular mat-vec `[∂y/∂x]ψ` when `ẋ=ψ`).
- **Reverse `(19)–(21)`:** seed output `ȳ = e_n` → `x̄ = [∂y/∂x]ᵀȳ` (**transpose** mat-vec `[∂y/∂x]ᵀψ` when `ȳ=ψ`); requires the forward primal run + **stored intermediates** ("tape"). This transpose product is exactly what an efficient Jacobian-free adjoint needs.
- **Implementations:** *source transformation* (Tapenade; Fortran) → faster, lower memory, can match analytic speed [153]; *operator overloading* (CoDiPack/dco; C++) → easy/extensible but **2–4× slower, ≥2× memory** (tapes). ADflow = **source transformation (Tapenade)**; DAFoam/SU2 = overloading.

## 4. Accelerating the Jacobian (§4.3)

Naive `∂R/∂w` = `O(N²)` residual evals, `N > 1e8` → infeasible → exploit sparsity. Even Jacobian-free needs an efficiently-built **PC matrix**.
- **§4.3.1 Single-cell residual:** re-code a reduced one-cell residual; minimal evals but **code duplication** → drift risk between primal and reduced routine; hard for object-oriented codes (OpenFOAM). Used in Mader's original ADjoint [30].
- **§4.3.2 Graph coloring (ADflow's choice):** partition Jacobian columns into **structurally-orthogonal colors** (no two share a nonzero row) → perturb a whole color at once via any forward method (FD / CS / forward-AD) → uses the **original** residual routine (no duplication, guaranteed consistency). Cost ∝ `n_color`. Example: 5-pt stencil on a 5×5 mesh → 26 FD evals naive → **6** with coloring. Coloring cost/complexity depends on mesh topology + stencil.

## 5. Solving the adjoint system (§4.4)

- **Direct (LU):** exact, fixed op count, but prohibitive for large, non-tightly-banded 3-D → not used in CFD here.
- **Fixed-point / stationary `(22)–(23)`:** `R_adj(ψ) = [∂f/∂w]ᵀ − [∂R/∂w]ᵀψ = 0`; `ψ^{n+1} = ψ^n + M·R_adj(ψ^n)`. Reuses the flow solver's iteration → **duality preserving** (same rate as primal), no factorization of `[∂R/∂w]ᵀ`; but uses only last-iterate info → **convergence always linear**.
- **GMRES (ADflow's choice) `(24)`:** Krylov `K_n = span{r₀, Ar₀, …, A^{n-1}r₀}`, `A = [∂R/∂w]ᵀ`, `r₀ = [∂f/∂w]ᵀ − [∂R/∂w]ᵀψ₀ (25)`. Uses info from all iterations → **faster** than fixed-point. NK primal ⇒ primal & transposed-Jacobian eigenvalues coincide ⇒ consistent primal/adjoint convergence.
- **Jacobian-free GMRES `(26)`:** mat-vec by **reverse-mode AD**: `w̄ = [∂R/∂w]ᵀ R̄`, set `R̄ = r₀` → product `[∂R/∂w]ᵀr₀ = Ar₀` **without forming/storing `[∂R/∂w]ᵀ`** → saves time+memory and needs **no sparsity/connectivity info** → extends trivially to any mesh (overset). **But** you still form/store the PC matrix `[∂R/∂w]_PCᵀ` `(27)` (fully-Jacobian-free is slower).
- **Ill-conditioning:** the transposed Jacobian from 3-D viscous turbulent flow is ill-conditioned → **strong PC essential**. Right-preconditioned adjoint `(27)`: `([∂R/∂w]ᵀ [∂R/∂w]_PCᵀ⁻¹)([∂R/∂w]_PCᵀ ψ) = [∂f/∂w]ᵀ`.

## 6. Nested preconditioning (§4.5, Fig 3)

PETSc; PC matrix `[∂R/∂w]_PCᵀ` ≈ `[∂R/∂w]ᵀ` but easily invertible (no general-purpose PC guaranteed → construct per problem). Stack:
- **GMRES** (top-level iterative solver).
- **ASM** global PC (1–2 overlap) with **outer Richardson** iterations.
- **ILU** local PC with **inner Richardson** iterations (Richardson damping = 1.0).
- Richardson iterations emulate a higher ILU fill at **lower memory** → use **low fill (1–2)**. `[∂R/∂w]_PCᵀ` is much sparser than the exact Jacobian → extra Richardson mat-vecs are cheap vs GMRES iters. **Outer (ASM) > inner (ILU)** in effectiveness (ASM gathers cross-block info; ILU is block-local). Defaults: **3 outer, 1 inner** (`innerPreconIts=1` in code; corrected 2026-08-12 — previously stated 2 inner).

## 7. Existing discrete-adjoint solvers (§4.5, Table 1)

| Solver | Org | Lang | Partial deriv | Adjoint solution |
|---|---|---|---|---|
| **ADflow** | U Michigan | Fortran | **Transformation AD** | **Krylov** |
| Cart3D | NASA | Fortran | Hand derived | Fixed point |
| **DAFoam** | U Michigan | C++ | Overloading AD | Krylov |
| elsA | ONERA | C++ | Hand derived | Fixed point |
| FUN3D | NASA | Fortran | Hand derived | Fixed point |
| HYDRA | Rolls-Royce | Fortran | Transformation AD | Fixed point |
| Jetstream | U Toronto | Fortran | Hand derived | Krylov |
| NSU3D | U Wyoming | Fortran | Hand derived | Krylov |
| SimpleFoam (piggy/reverseAcc) | Aachen | C++ | Overloading AD | Fixed point |
| SimpleFoam (checkpointing) | Aachen | C++ | Overloading AD | **Full-code AD** |
| STAMPS | QMUL | Fortran | Transformation AD | Fixed point |
| SU2 | Stanford/TU-KL | C++ | Overloading AD | Fixed point |
| TAU-Code | DLR | C | Hand derived | Krylov |

Fortran → transformation AD (Tapenade, faster/less memory); C++ → overloading (CoDiPack/dco, easy but slower). Full-code AD reference = adjointSimpleShapeCheckpointingFoam (Towara–Naumann). FEM semi-symbolic = dolfin-adjoint (FEniCS). Commercial: CFD-ACE+, Fluent, STAR-CCM+ (limited info).

## 8. The three proposed variants (§4.6)

| Variant | Partials | GMRES mat-vec | Stores `[∂R/∂w]ᵀ`? |
|---|---|---|---|
| **FD Jacobian** | coloring + **finite difference** (all of `[∂R/∂w]ᵀ`, `[∂R/∂w]_PCᵀ`, `∂R/∂x`, `[∂f/∂w]ᵀ`, `∂f/∂x`) | explicit stored Jacobian | yes |
| **AD Jacobian** | coloring + **forward-mode AD** (same set) | explicit stored Jacobian | yes |
| **Jacobian-free** (hybrid) | `[∂R/∂w]_PCᵀ` via **forward-AD + coloring**; `[∂f/∂w]ᵀ`, `∂f/∂x`, mat-vec `[∂R/∂w]ᵀr₀`, and `[∂R/∂x]ᵀψ` via **reverse-AD** | reverse-AD (no stored Jacobian) | **no** |

- All use **graph coloring** (not single-cell) for consistency. `[∂f/∂w]ᵀ` and `∂f/∂x` are **dense** vectors (e.g., drag integrated over all surface faces) → split `f` into `N_D` discrete cell faces, color+sum each `[∂f_i/∂w]ᵀ`.
- Geometric-DV totals go through `df/dx_v` then chain (see A.1), keeping the adjoint independent of the parametrization/warp.

## 9. Results (§5) — reference benchmark data

**WU `(28)`:** `T_WU = nT/T_ref`, `T` = wall time (s), `n` = cores, `T_ref` = TauBench time (`mpirun -np 1 ./TauBench -n 250000 -s 10`); **`T_ref = 2.972 s`** on Stampede2 Skylake (Xeon Platinum 8160, 2.1 GHz, 48 cores/node, 196 GB, 100 GB/s Omni-Path). Memory = peak RSS (GB). (Stampede2-era paper measurement — do not use for Deucalion sizing.)

### 9.1 Case 3 — low-speed rectangular wing (Table 2)
AR 6.12, NACA 0012, **102 912** cells, y⁺≈1.2, 20-chord domain; Re=1e6, M=0.15, SA; obj = C_D, target C_L=0.375; **8 twist DVs**; single Skylake core. Flow/adjoint iters: ADflow **1011/310**, DAFoam **720/143**; flow residual drop 14 (ADflow) / 12 (DAFoam) orders. Reference derivatives: complex-step (ADflow), full-code AD (DAFoam). Flow wall time (ADflow) = 148.6 s.

| Solver / option | Flow WU | Adj WU | Assy | Solve | Adj/flow | Flow GB | Adj GB | GB ratio |
|---|---|---|---|---|---|---|---|---|
| ADflow **Jacobian-free** (transform) | 50.0 | 41.1 | 2.4 | 38.7 | 0.8 | 1.7 | 2.1 | 1.2 |
| ADflow AD-Jacobian (transform) | 50.0 | 56.4 | 19.8 | 36.6 | 1.2 | 1.7 | 3.0 | 1.8 |
| ADflow FD-Jacobian | 50.0 | 45.3 | 9.3 | 36.0 | 0.9 | 1.7 | 3.0 | 1.8 |
| DAFoam Jacobian-free (overload) | 39.7 | 187.4 | 100.3 | 87.1 | 4.7 | 0.3 | 17.2 | 57.3 |
| DAFoam FD-Jacobian | 39.7 | 245.6 | 193.5 | 52.1 | 6.2 | 0.3 | 7.1 | 23.7 |
| adjointSimpleShapeCheckpointingFoam — full-code AD (overload) | 165.2 | 791.4 | — | — | 4.8 | 7.6 | 22.1 | 2.9 |

Findings: **ADflow Jacobian-free best on runtime + memory.** Its adjoint *solve* is ~6% slower than a stored-matrix product (reverse-AD per GMRES iter), but its **assembly is far cheaper** (2.4 vs 19.8/9.3) → total **27% / 9% faster** than AD/FD-Jacobian, **30% less memory**. Caveat: the stored-Jacobian runtime advantage grows as **#functions** grows (assembly amortized). DAFoam Jacobian-free fastest *of DAFoam* but **140% more memory** than DAFoam-FD (overloading tapes). ADflow **4.6× faster** than DAFoam-JF (assembly 2.4 vs 100.3; **162 vs 945 colors**) and **8.2× less memory**. DAFoam-JF **4.2× faster / 22.1% less memory** than full-code AD. Accuracy (Tables 3–4, structure only): Jacobian-free & AD-Jacobian match complex-step to **12 digits** (ADflow) / **10 digits** (DAFoam, vs full-code AD); **FD-Jacobian → 3–4 digits, ≈0.1% avg error**. *(96-digit validation tables not reproduced — verification artifacts.)*

### 9.2 Case 4 — transonic CRM wing (Tables 5–7)
Boeing-777-class, fuselage/tail removed, root at symmetry; C_L=0.5, M=0.85, Re=5e6; **8 twist DVs**; L1 mesh **3 604 480** cells; 48 cores (2 Skylake). Shock on upper surface → large negative sensitivity mid-chord + LE sensitivity at root. Flow converges 13 orders.

**Table 5 (48 cores):**
| ADflow | Flow WU | Adj WU | Assy | Solve | Adj/flow | Flow GB | Adj GB |
|---|---|---|---|---|---|---|---|
| Jacobian-free | 6064.1 | 4781.3 | 110.4 | 4670.9 | 0.8 | 44.6 | 87.2 |
| AD-Jacobian | 6064.1 | 5719.1 | 1176.0 | 4543.1 | 0.9 | 44.6 | 119.5 |
| FD-Jacobian | 6064.1 | 5014.9 | 446.4 | 4568.5 | 0.8 | 44.6 | 119.5 |

→ Jacobian-free **19.5% / 4.9% faster** than AD/FD, **37% less memory**; matches complex-step to **11 digits** (Table 6). 

**Table 7 — parallel scaling** (24 cores/node; efficiency `η = T_WU²⁴/T_WUⁿ`):
| Nodes | Cores | Flow (η%) | JF (η%) | AD-Jac (η%) | FD-Jac (η%) |
|---|---|---|---|---|---|
| 1 | 24 | 5628.3 (100) | 4525.0 (100) | 5455.6 (100) | 4765.9 (100) |
| 2 | 48 | 6064.1 (92.8) | 4781.3 (94.6) | 5719.1 (95.4) | 5014.9 (95.0) |
| 4 | 96 | 7166.4 (78.5) | 6828.0 (66.3) | 7885.9 (69.2) | 6998.3 (68.1) |
| 8 | 192 | 10237.7 (55.0) | 10147.6 (44.6) | 10712.9 (50.9) | 10012.9 (47.6) |
| 16 | 384 | 11581.0 (48.6) | 9208.2 (49.1) | 10006.1 (54.5) | 9407.2 (50.7) |
| 32 | 768 | 14180.6 (39.7) | 12334.0 (36.7) | 13168.2 (41.4) | 12574.1 (37.9) |

→ Adjoint parallel efficiency **>90%** when each proc owns **≥75 000 cells**. Using 48 (vs 24) cores/node lowers walltime but hurts efficiency (1536 cores: 30.1% vs 36.7%) — **memory-bandwidth**-bound, not CPU-bound (Skylake 2.1 GHz). (Stampede2-era paper measurement — do not use for Deucalion sizing.)

### 9.3 Case 5 — wing-body-tail overset (Tables 8–9)
CRM + fuselage + horizontal tail (= DPW4); Re=**4.3e7** (realistic scale); multiblock **overset**, **10 358 373** cells; **9 twist DVs**; 96 cores (4 Skylake). Shock on upper wing → high local sensitivity, weakened near root by wing-fuselage interaction. Only Jacobian-free run.

| ADflow (JF) | Flow WU | Adj WU | Assy | Solve | Adj/flow | Flow GB | Adj GB | GB ratio |
|---|---|---|---|---|---|---|---|---|
| — | 33792.0 | 23577.6 | 460.8 | 23116.8 | 0.7 | 123.7 | 252.9 | 2.0 |

→ **Adjoint faster than flow** (ratio 0.7); matches complex-step to **6 digits** (Table 9). Demonstrates Jacobian-free handling complex overset geometry (no connectivity/sparsity needed).

## 10. Conclusions + verdict (§6, Table 10)

Three variants compared on 6 factors. **Jacobian-free + source-transformation AD = best overall.**

| Option | Speed | Scalability | Memory | Accuracy | Effort | Extensibility |
|---|---|---|---|---|---|---|
| **Jacobian-free (transform)** | ● best | ◐ | ● best | ● exact | ○ high effort | ◐ |
| Jacobian-free (overload) | ○ (~460% slower) | N/A | ○ (~820% more) | ● | ◐ | ◐ |
| AD-Jacobian (transform) | ◐ (~27% slower) | ◐ | ◐ (~30% more) | ● | ○ | ◐ |
| FD-Jacobian | ◐ (~9% slower) | ◐ | ◐ (~30% more) | ○ (~0.1% error) | ○ | ◐ |
| Full-code AD (overload) | ○ (~1930% slower) | N/A | ○ (~1050% more) | ● | ● easy | ● |

Key numbers: Jacobian-free up to **27% / 9%** faster than AD/FD-Jacobian; up to **30% less memory**. Source-transformation JF **4.6× faster / 8.2× less memory** than overloading JF; **19.3× faster / 10.5× less memory** than full-code AD. **FD ≈ 0.1%** gradient error (fine for many uses, risky for tight/quasi-Newton optimization). Structured mesh → **162 colors** for PC; unstructured → **945** (≈5× more assembly cost). Source-transformation costs more **initial** dev effort but extending (new terms/objectives) is then easy.

---

# Appendix A — Solver-specific implementation (the implementation reference)

## A.1 States, residuals, DVs
- **ADflow:** states `w = [ρ, u, v, w, ρE, ν̃]` per cell (primitive set); `R` and `w` length = **`6·n_C`**; **cell-by-cell ordering** (better for block ILU than state-by-state). **Turbulence `ν̃` included in `w`** and differentiated in `∂R/∂w` — critical for accurate totals (no frozen turbulence). *(Branch note 2026-08-12: on this branch SA-GR adds γ and Re̅θt — 8 states per cell, so sizing statements scale accordingly.)*
- **DAFoam:** `w = [u, v, w, p, ϕ, ν̃]` (`ϕ` = cell-face flux, a state due to segregated SIMPLE + Rhie–Chow); length ≈ **`8·n_C`** (mixed cell/face), topology-dependent.
- **Geometric DVs `(A.1)`:** compute `df/dx_v` (volume coords), then `df/dx = df/dx_v · dx_v/dx_s · dx_s/dx` via **pyGeo** (FFD parametrization) + **IDWarp** (analytic inverse-distance mesh warp). Keeps the adjoint independent of parametrization/warp. (ADflow all-3 variants + DAFoam-JF use this; DAFoam-FD instead builds `∂R/∂x`, `∂f/∂x` directly because `∂R/∂x_v` is too dense to color.)

## A.2 Graph coloring
- **Structured multiblock:** analytic per-block colorings for the fixed stencil; halo cells handle cross-block terms. ADflow explicitly stores the full `[∂R/∂w]ᵀ` for FD/AD-Jacobian. **210 colors** for the full state Jacobian, **162** for the first-order PC (`[∂R/∂w]_PCᵀ`); counts based on all state variables tied to a residual (Fig A.10a sparsity).
- **Unstructured / overset:** connectivity is global → NP-hard, no analytic solution → **parallel heuristic coloring** (tentative local colors, resolve conflicts by global exchange). ADflow overset / DAFoam: **~1000 colors** (~3× the max nonzeros/row); state-Jacobian assembly cost ∝ colors, but colors depend only **weakly** on mesh size / core count (Fig A.10b denser, more irregular sparsity).

## A.3 Preconditioner-matrix construction
- **ADflow `[∂R/∂w]_PCᵀ`:** exactly linearize an **alternative simple (level-one) flux** → **≤7 nonzero blocks/row** (greatly reduces PC storage). JST → **second-order dissipation only** (the shock-active term), augmented by `σκ₄`; **shock sensor frozen** for the linearization. Upwind → **first-order reconstruction**. Viscous → velocity/temperature gradients from **adjacent cell-center differences** (per Yildirim). SA → **first-order upwind**. *(Branch note 2026-08-12: the SA-GR transition equations are also Tapenade-differentiated and verified — see `docs/VERIFICATION/three-stage-verification.md`.)* **Overset:** compute the PC on **local blocks only** (ignore overset stencils) via the analytic block coloring.
- **DAFoam `[∂R/∂w]_PCᵀ`:** first-order upwind inviscid (reduce stiffness); **shrink `p`, `ϕ` stencils from level-3 to level-2**; **ignore turbulence production** (reduce stiffness near y⁺=1 viscous layer); scale partials by cell **volume** / **face area** for diagonal dominance; **normalize each state-Jacobian column by its far-field reference**.

## A.4 Differentiation of residuals & functions
- **ADflow master routine** (Alg 1): a standalone re-implementation of the residual + objective, used **only** as an AD source. It recomputes **all** intermediates so derivatives accumulate correctly:
  ```
  Alg 1  master(x,w) -> R,f
    w_int  = F_int(x,w)        # intermediates: p from (ρ,T), μ via Sutherland, …
    w_bc   = F_bc(x,w)         # BC/halo update via MPI, far-field Mach, AoA
    x_geo  = F_geo(x)          # cell volume, face area/normal, wall distance
    R      = F_res(x,w,w_int,w_bc,x_geo)
    f      = F_obj(x,w,w_int,w_bc,x_geo)
  ```
- **Tapenade** produces forward `master_d` (Alg 2) and reverse `master_b` (Alg 3); a high-efficiency `master_state_b` (Alg 4) computes only `∂R/∂w`-type products (reduced routine reused across GMRES iters). In forward mode differentiated routines run in master order; in reverse mode, reverse order.
- **Manual assembly (not "differentiate the whole routine"):** differentiate subroutines, then hand-assemble, for two reasons — (a) forward-pass variables from the converged flow are **reused** across successive transpose mat-vecs (one forward pass serves many GMRES iters), and (b) variables known not to change between subroutine calls can be **kept off the AD tape** → less memory/time. Both matter because reverse-mode runs **every** GMRES iteration.
- **Jacobian-free procedure (ADflow):** (1) set color rows of `ẇ` to 1, call `master_d` per color → assemble `[∂R/∂w]_PCᵀ`; (2) `f̄=1`, `master_b` → `[∂f/∂w]ᵀ`; (3) with PC = `[∂R/∂w]_PCᵀ` and RHS `[∂f/∂w]ᵀ`, seed `r₀` and use `master_state_b` for `[∂R/∂w]ᵀr₀` in GMRES → solve for `ψ`; (4) `R̄=ψ`, `f̄=1`, `master_b` → total `df/dx = ∂f/∂x − ψᵀ ∂R/∂x`.
- **FD/AD-Jacobian (ADflow):** same but compute `[∂R/∂w]ᵀ`, `[∂R/∂w]_PCᵀ` by **FD** (FD-Jacobian) or **forward-AD** (AD-Jacobian) with coloring, store `[∂R/∂w]ᵀ` explicitly, and pass it to GMRES for the mat-vec.
- **Extending:** add a residual term → differentiate only its subroutine with Tapenade and call it in `master_d`/`master_b`. This is why the higher up-front effort of source transformation pays back.
- **DAFoam:** no single master routine (segregated) → separate `R_U, R_p, R_ϕ, R_ν̃, f` functions; **operator overloading (dco/c++)** because OpenFOAM is heavily templated C++ (impractical for source transformation). Its Jacobian-free is like ADflow's except `[∂R/∂w]_PCᵀ` is built by **FD** (not forward-AD); FD-Jacobian computes `∂R/∂x`, `∂f/∂x` by brute-force FD (per A.1).
