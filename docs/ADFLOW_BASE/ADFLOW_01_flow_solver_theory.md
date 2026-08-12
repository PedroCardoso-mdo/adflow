# ADFLOW 01 — Flow Solver & ANK/NK (complete extraction: Yildirim2019)

> **⚠️ SA-GR branch note (2026-08-12):** this is a faithful extraction of
> Yildirim et al. 2019 (transonic aircraft, plain SA). Several of its
> prescriptions are measured NOT to transfer to the SA-γ-Re̅θt model on this
> branch: the NK switch at 1e-4…1e-6 stalls (premature engagement, Step=0.00 —
> see `docs/convergence-strategy.md`, commit ef9fc10d); "decoupled beats
> coupled" is inverted (SA-GR requires CANK→CSANK); L2Convergence code default
> is 1e-8 and 1e-8 is the standing target for transition runs. Hardware/WU
> numbers are Ivybridge/Skylake-era paper data, not Deucalion. Treat everything
> below as paper record, checked against `docs/convergence-strategy.md` before
> use on SA-GR.

**Source:** Yildirim, Kenway, Mader, Martins, *A Jacobian-free approximate Newton–Krylov startup strategy for RANS simulations*, JCP (2019), doi:10.1016/j.jcp.2019.06.018.
**Coverage:** comprehensive technical extraction of §1–§6. Abstract + narrative motivation dropped; hard facts retained. Prose paraphrased; equations and benchmark data reproduced as facts. Eq numbers `(n)` and table/section refs are the paper's.

---

## 0. Context facts worth keeping from §1

- **Why robustness dominates:** high-fidelity aero/aerostructural MDO runs **hundreds** of successive RANS + adjoint solves on meshes up to **tens of millions** of cells; aerostructural cases can take ~3 days on ~1000 cores. A solver failure stalls the whole optimization and needs human intervention, which is costlier than slower-but-robust solves. The optimizer requests shapes a human never would (the "circle → supercritical airfoil" case, He et al. [9]) — RANS may be invalid on intermediate shapes but gradients still give correct trends.
- **NK basics:** Newton solves the nonlinear system; a Krylov method (GMRES) solves each linear system; matrix-free Krylov slashes memory. If the initial guess is outside the **basin of attraction**, NK diverges → a **globalization** method is required.
- **PTC** (pseudo-transient continuation, Kelley–Keyes) is the standard globalization: start with backward Euler (small time step → stable), ramp the step toward Newton. Backtracking/trust-region also used.
- **ANK lineage:** Hicken–Zingg (3-D Euler startup) → Chisholm–Zingg (2-D RANS) → Osusky–Zingg (3-D RANS). "ANK" = globalized NK driven by a **matrix-based** approximate Jacobian, with Jacobian + PC **lagged** to amortize cost.
- **The problem this paper fixes:** lagging a matrix-based Jacobian **3–5 iterations is fine for Euler** but **fails for 3-D RANS** — Osusky–Zingg reported that lagging the matrix-based Jacobian does *not* improve performance and can **diverge**. RANS wants the approximate Jacobian always up to date → huge cost if matrix-based.
- **Two contributions:** (1) use **approximate residual routines in a matrix-free (JFNK) ANK** so the approximate Jacobian is always current, lower-bandwidth, better-conditioned, and cheaper per mat-vec — which lets you **aggressively lag only the preconditioner** without an outdated-Jacobian instability; (2) an **adaptive PC-lagging** algorithm. Neither requires any hand differentiation; both build on the exact/approximate residual routines only. Demonstrated in open-source **ADflow**.

---

## 1. Governing equations (§2.1–§2.2)

- **Navier–Stokes, integral FV form `(1)`:** `∂/∂t ∫_V Q dV + ∮ F·n dA = ∮ D·n dA`.
- **State / fluxes `(2)`:** `Q = [ρ, ρu, ρv, ρw, E]ᵀ`. Convective flux `F` (per direction) has rows `[ρu; ρu²+p; ρuv; ρuw; u(E+p)]` etc.; diffusive flux `D` has rows `[0; τ; …; u·τ − q]` (viscous stresses `τ_ji`, heat flux `q_i`). Laminar viscosity via **Sutherland**; closed by **ideal gas**.
- **Turbulence — Spalart–Allmaras `(3)`:** `∂ν̃/∂t + u_j ∂ν̃/∂x_j = P(ν̃) − M(ν̃) + D(ν̃)` (production/destruction/diffusion). Working variable `ν̃`. Default variant **SA-noft2** (per NASA TMR). Negative-modified-vorticity guard `(4)` (Allmaras): `S̃ = S + S̄` if `S̄ ≥ −c_v2 S`, else `S̃ = S + S(c_v2²S + c_v3 S̄)/((c_v3 − 2c_v2)S − S̄)`, with `S̄ = ν̃ f_v2/(κd)²`, `c_v2 = 0.7`, `c_v3 = 0.9`.
- **Turbulence discretized to first order in space** (advection = 1st-order upwind). They tuned the solver for first-order SA; found little difference vs second-order for MDO. Solver still applies to second-order SA and to non-SA models.

## 2. Discretization & the residual ladder (§2.3)

- Structured **multiblock + overset**; cell-centered ρ, momentum, energy, ν̃. Residual for nonlinear iter n `(5)`: `R0^(n) = R0(Q^(n)) = −∂/∂t ∫_V Q dV |_n`.
- **Inviscid:** JST scalar dissipation. **Viscous gradients:** Green–Gauss. Full residual `R0` = **33-point** stencil (26 nearest + 6 second-level + self).
- **Approximation levels (Table 1):**

| Level | Stencil | Inviscid (NS) | SA advection | SA other |
|---|---|---|---|---|
| `R0` | 33 | JST + scalar dissipation, **up-to-date shock sensor** (entropy) | 2nd-order | 1st-order diffusion; **production kept** |
| `R1` | 27 | **omit 4th-order dissipation**; **freeze shock sensor** (entropy) at `Q^(n)` (drops the 6 second-level neighbors) | 2nd-order | production **omitted** |
| `R2` | 7 | as R1 **+ orthogonal-mesh assumption** (face velocity gradients from face-sharing neighbors only) | 1st-order (orthogonal-mesh) | production **omitted** |

- `R2`'s 7-point (first-order) form is the classic PC/adjoint stencil (Osusky, Kenway, Lyu). SA stencil is always 7-point.
- **Steady state `(6)`:** `R0(Q) = 0`. **Convergence `(7)`:** `η_abs ≥ ‖R0(Q^n)‖₂/‖R0(Q_fs)‖₂`, default `η_abs = 1e-12`.

## 3. Variable/residual scaling (§2.4)

- States normalized by freestream ⇒ entries `O(1)`. Cell volumes span ~**10 orders** (RANS stretching) → without scaling the algorithm preferentially reduces large-volume cells → **scale each cell residual by 1/volume**.
- **Turbulence residual scaled by `1e4`** for adequate magnitude balance vs mean flow. Because scaling is baked into the residual routines, no extra handling is needed in the linear solves. (Strategy à la Chisholm–Zingg.)

## 4. ANK algorithm (§3)

**Per-iteration update `(8)`:** `[ I/Δt^n + (∂R_m/∂Q)^n ] ΔQ^n = −R0(Q^n)`. `I/Δt^n` = diagonal pseudo-time term with the **same CFL in every cell** (global CFL); `m` = matrix-free approximation level.

### 4.1 Linear solver (§3.1)
- **GMRES**, matrix-free → **JFNK** extended to ANK.
- **Mat-vec `(9)`:** `(∂R_m/∂Q)^n v ≈ [R_m(Q^n + εv) − R_m(Q^n)]/ε`.
- **FD step `(10)`:** `ε = e_rel·vᵀQ^n/‖v‖₂²` if `|vᵀQ^n| > u_min‖v‖₁`, else `ε = e_rel·u_min·sign(vᵀQ^n)‖v‖₁/‖v‖₂²`. `u_min=1e-6`, `e_rel=1e-8`.
- **Zero initial guess; no restart across nonlinear iters** — the descent *direction* changes each nonlinear iteration, so a previous linear solution gives no benefit as an initial guess (unlike a fixed-RHS solve).

### 4.2 Approximate Jacobian (§3.2 — contribution #1)
- Matrix-free with `R0`/`R1`/`R2`. **`R0` = true Jacobian**; `R2` = first-order Jacobian (as Osusky's startup); `R1` = intermediate. Selectable at runtime by changing residual options (no differentiation).
- **Trade-off:** little/no approximation → better nonlinear convergence but costlier linear solves (larger bandwidth); heavy approximation → lower nonlinear rate but faster linear solves (more diagonal dominance). **`R1` is the default** (accuracy vs Jacobian stiffness balance).
- Larger stencils let you carry higher-order approximate Jacobians than first-order *without* storing more matrix entries (matrix-free).

### 4.3 Preconditioning (§3.3)
- **Right preconditioning `(11)`:** `[(I/Δt + ∂R_m/∂Q) M⁻¹](MΔQ) = −R0`, with `M = I/Δt + ∂R2/∂Q`. `M⁻¹` = **ASM** (subdomains) → local **RCM** reorder + **ILU**, via **PETSc**. Store `M` and `M⁻¹`.
- `M` built from `R2` ⇒ cheap to store/factorize (7-block rows). `∂R2/∂Q` computed by FD or AD (Kenway) with **coloring** (Lyu; Goldfarb–Toint 7-pt 3-D coloring). *Here, FD only* for ANK PCs.

### 4.4 Adaptive PC lagging (§3.4 — contribution #2)
- Because the matrix-free operator is always current, **only the PC is lagged** ⇒ no outdated-Jacobian instability. Call it **PC lagging**, not Jacobian lagging.
- **Refresh triggers:** (a) the linear solve **hits the max GMRES iterations this iteration AND the previous iteration's linear solve succeeded** (PC has degraded); (b) relative nonlinear convergence since the last PC build reaches `η_pc = 0.5`; (c) safety cap after `n_lag = 10` lagged iterations.
- If the linear solve fails to reach tolerance even with a fresh PC, **terminate at max iterations and take the partial step** (usually still reaches ~0.1 relative). If it stagnates well above that, nonlinear convergence stalls → this is a PC-adequacy problem the lag algorithm can't fix → **restart with a stronger PC / higher iteration limit**.
- **CFL is updated only on PC-refresh iterations** (keeps the time term in `M` consistent with the mat-vec).

### 4.5 Linear-solve tolerance (§3.5)
- **Tolerance `(12)–(13)`:** `η_lin ≥ ‖R_lin^n‖/‖R0^n‖`, `R_lin = [I/Δt + ∂R_m/∂Q]ΔQ + R0`. Default **`η_lin = 0.05`** (largest value that reliably gives a usable update without destabilizing). Literature range **1e-4…1e-1**.
- **Fixed tolerance, not Eisenstat–Walker.** EW is tuned to the nonlinear gain a *full* Jacobian would give; with approximations the realized gain is lower, so EW predicts too-loose a tolerance → inaccurate updates → bad for the robustness-first startup. (EW *is* used in NK.)
- **Max GMRES iterations `n_lin = 40`.** Sized so most solves converge but without wasting effort — another nonlinear step usually beats over-solving. Rare persistent failures → restart with higher limits and/or stronger PCs.

### 4.6 Solution update (§3.6)
- **Update `(14)`:** `Q^{n+1} = Q^n + ω^n ΔQ^n`, `ω^n ∈ [0,1]` (Modisette approach).
- **Physicality check `(15)` / Alg 1:** cap fractional change of **ρ, E** to `θ_phys = 0.2` per cell. For **ν̃**, check only *negative* updates, cap at `θ_phys,ν̃ = 0.99` (lets SA grow fast, forbids negative). **Momentum unchecked** (must change sign). Global `ω_phys` via MPI `min`.
  ```
  Alg 1 (physicality):  ω_local=1
    for each owned cell i:
      for l in {ρ,E}:      ω_local = ( max(|ΔQ_{i,l}/(θ_phys,l Q_{i,l})|, 1/ω_local) )^{-1}
      for l = ν̃, if ΔQ_{i,l}<0: same cap
    ω_phys = MPI_Allreduce(ω_local, MIN)
  ```
- **Backtracking line search `(16)–(17)` / Alg 2:** starting from `ω_phys`, backtrack by `θ_bt = 0.7` until the **unsteady** residual `R_u = ‖ω(I/Δt)ΔQ + R0(Q+ωΔQ)‖₂` satisfies `R_u ≤ ‖R0(Q^n)‖₂`.
  ```
  Alg 2 (backtracking):  ω_u=ω_phys; compute R_u
    if R_u > ‖R0‖:
      while R_u > ‖R0‖ and ω_u>ω_min: ω_u ← ω_u·θ_bt; recompute R_u
    ω^n ← ω_u
  ```
- **Critical:** line search reduces the **unsteady** residual; the printed **steady/total** residual can *increase* while unsteady decreases (time term). Only `R0` guarantees a descent direction; `R1`/`R2` sometimes don't. On a failed reduction the solver **lowers CFL** (weights the time term). If `ω < ω_min = 0.01` above the CFL floor → reject (`ω=0`); next iteration changes CFL + refreshes PC. At the CFL floor, take the (small) step anyway.

### 4.7 Pseudo-time / CFL evolution (§3.7)
Ceze–Fidkowski monotonic residual-difference method (mRDM), from Bücker RDM, with hard CFL limits.
- **mRDM ramp `(18)`:** `CFL^n = CFL^k · α^γ`, `γ = max((‖R0^k‖ − ‖R0^n‖)/‖R0^k‖, 0)` (`k` = last PC/CFL-update iter).
- **SER floor `(19)`** (van Leer–Mulder): `CFL_min^n = (‖R0_fs‖/‖R0^n‖)^β`. Also sets the initial CFL for **warm starts** (previous design's solution as init → starts at higher CFL than a cold start).
- **PC-update CFL logic `(20)`:** if `ω^{n-1} > ω_ramp` → `min(CFL^k α^γ, CFL_max)`; elif `ω_ramp ≥ ω^{n-1} > ω_min` → `max(CFL^k, CFL_min)`; else → `max(θ_red·CFL^k, CFL_min)`. Non-PC iterations hold CFL. `CFL(0)=5`, `CFL_max=1e5`, `α=10`, `β=0.5`, `θ_red=0.5`, `ω_ramp=0.4`, `ω_min=0.01`.

### 4.8 Turbulence handling (§3.8)
- **Decouple by default** (nonlinear block Gauss–Seidel): update flow (frozen `ν̃`) → recompute turb residual with updated flow → update `ν̃` (frozen flow). Two linear systems: **5×5-block flow** + **scalar turbulence**.
- **Why:** recovers the flow's favorable convergence (turb residual can vary wildly without polluting flow); improves conditioning (fixes the ν̃-vs-flow scaling problem); allows **independent PC refresh** (turbulence PC needs it more often). Osusky auto-scales ν̃ by its domain max; Chauhan et al. — optimal solver changes with convergence stage (Gauss–Seidel early, Newton late).
- **Coupled** = standalone Newton on one system, reserved for the terminal stage where the `1e4` scaling suffices.

### 4.9 Full ANK algorithm (Algorithm 3)
```
init n=1; n_pc,f=0; n_pc,t=0; k=0; CFL0=5; ω0=1; R0^0 = R0(Q^0)
while ‖R0^n‖ > η_ANK ‖R0^0‖:
  if mod(n_pc,f, n_lag)=0 or ‖R^{n-1}‖ < η_pc ‖R^k‖:      # flow PC refresh
     CFL^n ← eq(20);  M_f⁻¹ ← factorize[I/Δt + ∂R2/∂Q]_f;  k=n; n_pc,f=0; n_pc,t=0
  ΔQ_f ← solve[I/Δt + ∂R1/∂Q]_f ΔQ_f = −R0_f^n         # default matrix-free op = R1
  ω_phys,f ← Alg1(ρ,E);  ω_f ← Alg2
  if ω_f<ω_min and CFL^n>CFL_min: ω_f=0
  Q_f^{n+1} = Q_f^n + ω_f ΔQ_f
  if n_pc,t=0: M_t⁻¹ ← factorize[I/Δt + ∂R2/∂Q]_t     # turb PC refresh
  R0,t^n ← R0,t(Q_f^{n+1}, Q_t^n)
  ΔQ_t ← solve[I/Δt + ∂R1/∂Q]_t ΔQ_t = −R0,t^n
  ω_phys,t ← Alg1(ν̃);  ω_t ← Alg2
  if ω_t<ω_min: ω_t=0
  Q_t^{n+1} = Q_t^n + ω_t ΔQ_t
  R0^{n+1} ← R0(Q^{n+1}); n_pc,f++ ; n_pc,t++
  if ω_f<ω_min or (η_lin^n>η_lin and η_lin^{n-1}≤η_lin): n_pc,f=0    # flow update broke down
  if ω_t<ω_min or (η_lin,t^n>η_lin and η_lin,t^{n-1}≤η_lin): n_pc,t=0 # turb update broke down
  n++
```
**Coupled variant:** drop the `_t` blocks; solve flow+turbulence as one linear system; apply Alg1/Alg2 to the combined update. Matrix-free op selectable `R0`/`R1`/`R2` (default `R1`; lines 14 & 23).

**Parameter philosophy (§3.9):** many tunables is a feature — users adapt the *solver* to a hard case rather than the case. No single default set is optimal everywhere; robustness-vs-performance trade. Scenarios needing changes: massive separation / heavy shocks → ill-conditioned Jacobians → stronger linear solver (higher iter limit and/or stronger PC); terminal solver failing after ~1e-5 → **lower `η_ANK`**.

### 4.10 Default parameters (Table 2)
| Group | Symbol | Default |
|---|---|---|
| Linear solver | matrix-free op `R_m` | **R1** |
|  | ILU fill | 2 |
|  | ASM overlap | 1 |
|  | `η_lin` | 0.05 |
|  | GMRES iter limit `n_lin` | 40 |
| PC lagging | `η_pc` | 0.5 |
|  | max lag `n_lag` | 10 |
| Step size | `θ_phys,ρ = θ_phys,E` | 0.20 |
|  | `θ_phys,ν̃` | 0.99 |
|  | backtrack `θ_bt` | 0.7 |
| PTC/CFL | mRDM `α` | 10 |
|  | SER `β` | 0.5 |
|  | `CFL(0)` | 5 |
|  | `CFL_max` | 1e5 |
|  | cutback `θ_red` | 0.5 |
|  | ramp threshold `ω_ramp` | 0.4 |
|  | min step `ω_min` | 0.01 |
| Nonlinear | `η_ANK` (switch to NK) | 1e-5 |
|  | `η_abs` (final) | 1e-12 |

Tuned for **transonic 3-D full/partial aircraft, millions–tens of millions of cells**; also expected OK for subsonic and 2-D; much smaller/larger meshes may need adjustment.

## 5. Computational framework (§4)

- **ADflow** (open source): FV, Fortran90 + Python; built from **SUmb** (explicit RK). Added: ANK, Jacobian-free adjoint (Kenway), overset (implicit hole cutting), Python API. Solves Euler / laminar NS / RANS; steady / unsteady / time-spectral; multiblock + overset. Uses **PETSc** for GMRES, preconditioning, mat-free.
- **Solver switching by relative convergence `(21)`:** `η_rel^n = ‖R0^n‖/‖R0_fs‖`.
- **Available solvers:**
  - **RK** — 5-stage 4th-order low-memory Runge–Kutta (the only carryover from SUmb).
  - **D3ADI** — diagonalized diagonally-dominant ADI (Klopfer).
  - **ANK** — this paper.
  - **NK** — fully coupled JFNK; PC formed as in ANK, but matrix-free op is **exact `R0`**; **no time term** (CFL = ∞); linear tolerance via **Eisenstat–Walker**.
- RK/D3ADI are MG smoothers (not applicable to overset). **ANK** works on multiblock + overset → **startup of choice for overset**. **NK** gives the best terminal rate *if inside the basin*, else diverges → needs ANK/MG first. Practice: **ANK to `η_rel=1e-5`, then NK to `1e-12`**. Witherden suggests `1e-4` for the NK switch; a *lower* value gives a smoother transition and avoids switch failures. (Falsified for SA-GR — see banner.)

## 6. Results (§5) — reference benchmark data

**Work unit (WU):** TauBench, `mpirun -np 1 ./TauBench -n 250000 -s 10`, avg of 10 runs. Metrics: **NI** (nonlinear iters), **LI** (cumulative linear iters), **kWU** (thousand WU). `η_rel=1e-5` results = ANK stage; `1e-12` = total (ANK+NK).

### 5.1 Baseline — ONERA M6 (Table 3)
8M-cell multiblock; M=0.84, Re=11.7e6, T_ref=310.93 K, α=3.06°; 120 cores (6 Ivybridge×20), 1 WU=6.24 s. Six tests (R0/R1/R2 × decoupled/coupled).

| # | Level | Turb | NI@1e-5 | LI@1e-5 | kWU@1e-5 | NI@1e-12 | LI@1e-12 | kWU@1e-12 |
|---|---|---|---|---|---|---|---|---|
| 1 | R2 | decoupled | 122 | 2184 | 9.44 | 133 | 2391 | 15.86 |
| 2 | R1 | decoupled | 114 | 1967 | 9.42 | 126 | 2176 | 15.89 |
| 3 | R0 | decoupled | 71 | 2139 | 9.97 | 77 | 2371 | 17.08 |
| 4 | R2 | coupled | — fail — | | | | | |
| 5 | R1 | coupled | 156 | 2723 | 13.30 | 161 | 2905 | 18.78 |
| 6 | R0 | coupled | 95 | 2758 | 14.54 | 102 | 2974 | 21.11 |

Findings: first ~2 orders identical across all (low CFL → approximations inert). `R0` decoupled = fewest NI but **worse walltime** than `R1`/`R2` (costlier linear). **Decoupled beats coupled** (single-configuration result; inverted for SA-GR — see banner); approximate `R1`/`R2` beat exact `R0` on walltime. `R2` coupled diverged.

### 5.2 Lagged-PC eigenvalue study (§5.2, Figs 4–6)
Three convergence stages: **initial** (fs→1e-2, low CFL, strong transients near no-slip surfaces), **intermediate** (1e-2→1e-4, far-field/no-slip interaction, strong turbulence transients), **final** (→1e-5, CFL at `CFL_max`, flow settled). Method: converge to η_target; fix CFL; refresh PC; raise linear iter limit to 300; solve to 0.05 (as baseline) then further to `η_lin=1e-8` for Arnoldi eigen-estimates (not to 1e-16 — FD `(9)` truncation corrupts eigenvalues); lag PC for 10 iters watching spectra at iters 1–5,10.
- **η_target=1e0 (CFL 5):** flow PC spectrum spreads **faster** than turbulence → flow field changing most.
- **η_target=1e-3 (CFL 387.1):** turbulence PC spreads **faster** → turbulence changing most.
- **η_target=1e-5 (CFL 1e5):** both stay clustered (small updates).
→ A single fixed lag would be set by the worst stage; **adaptive lag** + **decoupling** (independent flow/turb PC refresh) win.

### 5.3 Generalization to discretizations — M6 (Table 4)
JST scalar (baseline), JST matrix dissipation (Swanson–Turkel), upwind Roe. For `R1`/`R2`: use **first-order state reconstruction** on faces (Roe flux itself unmodified); same stencils as scalar JST. Default params (tuned for scalar JST) reused.

| Discretization | NI@1e-5 | LI@1e-5 | kWU@1e-5 | NI@1e-12 | LI@1e-12 | kWU@1e-12 |
|---|---|---|---|---|---|---|
| JST scalar | 114 | 1967 | 9.42 | 126 | 2176 | 15.89 |
| JST matrix | 179 | 3570 | 15.56 | 184 | 3943 | 26.83 |
| Upwind (Roe) | 79 | 1973 | 10.59 | — NK fails — | | |
| Upwind, No NK | 79 | 1973 | 11.18 | 187 | 4327 | 23.45 |

→ Upwind **NK fails** after the switch, but **ANK alone reaches 1e-12** — approximations make the linear systems easier ⇒ robustness.

### 5.4 Matrix-free vs matrix-based — CRM WB overset (Table 5)
DPW6 overset, 14M cells (Coder et al.); M=0.85, Re=5e6, T_ref=310.93, α=2.4°; 160 cores (8 Ivybridge×20), 1 WU=6.24 s. 18 tests. **MB** = matrix-based ANK (lag Jacobian **and** PC); **MF** = matrix-free (lag PC only). "Adaptive" = Algorithm 3.

| # | Variant | Lag | NI@1e-5 | LI@1e-5 | kWU@1e-5 | NI@1e-12 | LI@1e-12 | kWU@1e-12 |
|---|---|---|---|---|---|---|---|---|
| 1 | MB-R2 | 1 | 191 | 3242 | 42.23 | 213 | 3908 | 83.88 |
| 2 | MB-R2 | 2 | 412 | 7046 | 59.76 | 428 | 7609 | 95.08 |
| 3 | MB-R2 | 3 | 195 | 3305 | 23.66 | 220 | 3961 | 80.43 |
| 4 | MB-R2 | 4 | 263 | 3855 | 27.47 | 287 | 4548 | 70.33 |
| 5 | MB-R2 | 5 | 270 | 3801 | 25.27 | 293 | 4449 | 65.51 |
| 6 | MB-R2 | 10 | 635 | 7135 | 39.42 | 652 | 7706 | 75.09 |
| 7 | MB-R2 | adaptive | 290 | 4039 | 23.31 | 314 | 4748 | 67.09 |
| 8 | MF-R2 | 1 | 190 | 3232 | 52.42 | 212 | 3895 | 93.82 |
| 9 | MF-R2 | 2 | 199 | 3385 | 59.98 | 223 | 4047 | 100.96 |
| 10 | MF-R2 | 3 | 180 | 3151 | 41.89 | 200 | 3775 | 80.42 |
| 11 | MF-R2 | 4 | 186 | 3239 | 30.62 | 208 | 3884 | 70.89 |
| 12 | MF-R2 | 5 | 240 | 3947 | 34.36 | 260 | 4615 | 75.67 |
| 13 | MF-R2 | 10 | 197 | 3688 | 30.21 | 218 | 4323 | 69.93 |
| 14 | MF-R2 | adaptive | 277 | 4486 | 35.88 | 298 | 5192 | 79.91 |
| 15 | MF-R1 | 10 | 102 | 2971 | 27.74 | 120 | 3563 | 64.52 |
| 16 | MF-R1 | adaptive | 106 | 2850 | 26.96 | 123 | 3456 | 64.64 |
| 17 | MF-R0 | 10 | 71 | 1775 | 16.65 | 87 | 2448 | 58.3 |
| 18 | MF-R0 | adaptive | — fail — | | | | | |

Findings: for large lags, **MF beats MB on NI** (expected); but **MB beats MF on walltime** because a matrix-based mat-vec is **10–20% cheaper** than the FD mat-vec. Tests 1 & 8 (both always-current) → nearly identical NI (small diff = FD truncation accumulation). Approximate formulations need **more NI but similar kWU** (Test 7 vs 16). **Test 18 (MF-R0 adaptive) fails** — adaptive lag refreshes the PC too often → prematurely ramps CFL to unstable values. R0 fails at high CFL; R1 needs fewer NI than R2; R2 oscillates late (line search limits steps; total residual can jump). This overset mesh preserves orthogonality → R2 fine here. **MB not used in production** (destabilizes hard cases; relies on step-limiting to stay stable); MF's physicality/line-search act as a safety net, not a crutch.

### 5.5 Approximation levels — CRM WB multiblock (Table 6)
DPW6 multiblock, **32M** cells (more layer skewness than overset); same flow conditions; 20 Ivybridge nodes. Cap 10 000 cumulative LI. Switch tolerances `η_ANK = 1e-5` and `1e-6`.

| Level | η_ANK | NI@1e-5 | LI | kWU | NI@1e-6 | LI | kWU | NI@1e-12 | LI | kWU |
|---|---|---|---|---|---|---|---|---|---|---|
| MF-R2 | 1e-5 | — fail — | | | | | | | | |
| MF-R1 | 1e-5 | 95 | 2484 | 61.24 | 98 | 2654 | 84.48 | — fail — | | |
| MF-R0 | 1e-5 | — fail — | | | | | | | | |
| MF-R1 | 1e-6 | 95 | 2484 | 62.27 | 119 | 3494 | 82.02 | 135 | 4801 | 258.75 |

→ On **skewed multiblock**: `R2` breaks (orthogonal-mesh assumption invalid), `R0` fails (extra terms → too stiff → linear solves fail); **only `R1` converges**. Premature NK switch (1e-5) fails; **switching at 1e-6 converges**. Justifies `R1` as the default matrix-free operator.

### 5.6 Coupling turbulence — strut-braced wing (Table 7)
PADRI SBW overset, **6.4M** cells (Secco); complex overset → no coarser levels (no MG startup possible); wing-strut junction has a standing shock + shock-induced & junction-induced separation. **SA-noft2-QCR2000** (Spalart QCR). M=0.72, Re=2.3e7, T_ref=228.71, α=1.0°; 48-core Skylake node, 1 WU=2.98 s.

| Turbulence | NI@1e-5 | LI | kWU | NI@1e-12 | LI | kWU |
|---|---|---|---|---|---|---|
| Decoupled | 76 | 915 | 7.72 | 95 | 1356 | 22.7 |
| Coupled | — fail (stalls before 1e-5) — | | | | | |

→ Coupled **stalls**: the complex junction flow prevents finding an update that reduces the unsteady residual → line search returns the minimum step. **Decoupling wins even for tightly-coupled flow where residuals fluctuate.** (Single-configuration result; inverted for SA-GR — see banner.)

### 5.7 vs multigrid — 1172-airfoil sweep (Table 8)
Webfoil subsonic airfoils (Li et al.), 1172 shapes, 35 840 cells each, structured 2-D; M=0.45, Re=6.5e6, T_ref=310.93, α=2.5°; 1 Skylake core each (embarrassingly parallel over 48), 1 WU=2.98 s. Four sets: **3w/4w/5w multigrid** (D3ADI smoother, CFL fixed 5) vs **ANK** (ILU fill reduced to 1, else defaults). Then NK to 1e-12. AM/SD over successful runs; 4 mesh-gen failures per method excluded; other failures = 5000-LI cap (MG: 5000 nonlinear iters).

| Method | Failures | NI@1e-5 (AM/SD) | WU@1e-5 (AM/SD) | NI@1e-12 (AM) | WU@1e-12 (AM) |
|---|---|---|---|---|---|
| 3w | 0 | 604.86 / 129.02 | 78.88 / 13.14 | 616.39 | 93.47 |
| 4w | 5 | 197.42 / 45.94 | 29.62 / 5.74 | 204.65 | 41.28 |
| 5w | 0 | 160.91 / 29.64 | 27.09 / 3.80 | 167.35 | 37.85 |
| ANK | 3 | 35.16 / 10.69 (NI); 381.63 / 110.79 (LI) | 24.54 / 6.24 | 46.79 | 35.66 |

→ >99% of shapes converge for every method. **ANK has the best average WU** and needs **only the finest level**. Note ANK/NK residual routines are **vectorized** (~1.5× speedup); MG residuals are not. Disabling ANK vectorization → 35.67 WU (between 3w and 4w). Net: **ANK ≈ multigrid performance but requires no coarse levels → enables overset**, where coarse levels usually can't be generated.

## 7. Conclusions (§6)

- Two contributions: (1) **approximate residual computations for a matrix-free approximate Jacobian** (eliminates the lagged-matrix-based-Jacobian robustness failure; enables aggressive PC lagging); (2) **adaptive PC-lagging**. Neither needs manual differentiation; both generalize across discretizations.
- **Trade-off confirmed:** approximate formulations need **more NI but cheaper NI** → higher net convergence rate during startup.
- **`R1` chosen as default:** it converged **all** cases (including many during development); `R0`/`R2` each failed on some (stiffness / skewness). Robustness prioritized because failures cost human intervention in optimization loops.
- **Eigenvalue evidence:** PC effectiveness decays as it's lagged, at a rate that depends on convergence stage, and differently for flow vs turbulence → adaptive lag + decoupling justified.
- PCs built from `R2` via FD + coloring → **no manual differentiation**; any residual change is directly reflected in the implicit system → wide discretization generality.
- Demonstrated on M6, CRM WB (overset + multiblock), SBW, and a 100k+-simulation airfoil study.
