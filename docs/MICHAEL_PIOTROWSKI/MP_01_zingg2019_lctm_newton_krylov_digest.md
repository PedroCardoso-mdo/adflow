# MP_01 digest — Piotrowski & Zingg (2019), LCTM in a Newton-Krylov algorithm

Companion to [`MP_01_zingg2019_lctm_newton_krylov_full.md`](MP_01_zingg2019_lctm_newton_krylov_full.md).
**AIAA SciTech 2019 — the origin paper.** This is the pre-smooth model
(designated **LCTM-SA**), i.e. the state *before* SA-sLM2015. Its value to us
is the *diagnosis*: it states exactly which numerical pathologies the strong
implicit solver hits, and which two fixes cured them.

## Why this matters to us

Our deep-NK wall reproduces the two failure modes this paper names. The paper
is explicit that **neither is fixed by the preconditioner** — they are fixed by
(a) restricting the source-term time step and (b) removing kinks from the
residual.

## 1. Two failure modes, measured (§III.C, ONERA M6, 5 M nodes)

Four combinations of source-term time stepping (`stts`) and smoothing (`smth`):

| Condition | `stts=0` | `stts=1, smth=0` | `stts=1, smth=1` |
|---|---|---|---|
| Subsonic (M 0.30, Re 5e6) | **diverges** — large turb/transition updates destabilise; large off-diagonals wreck the approximate-Jacobian PC, linear solver fails | stalls at rel **1e-5** | **machine zero** |
| Transonic (M 0.75, Re 2e7) | does not diverge but **stalls** — negative γ updates force excessive clipping | converges only **3 orders** | **machine zero** |

Read the middle column carefully: **with source-term time stepping already on,
the non-smooth model still stops at 1e-5 / 3 orders.** That is a residual-kink
wall, not a PC wall, and it is the same order of magnitude at which our runs
give up. The cure was smoothing the residual, nothing else.

> "These non-smooth functions, referred to as **simple kinks**, produce
> discontinuities in the Jacobian. This presents a problem for both the
> Newton-Krylov solution algorithm and the discrete-adjoint gradient-based
> optimization algorithm, and **can lead to solution stalling**."

They also note the specific offender: `Fonset`, "through the presence of
several layers of min/max operators" — and that other groups blamed this
stalling on turbulence/transition coupling stiffness and relaxed the coupling
to escape it, whereas P&Z fixed the function itself.

**Directly actionable for us:** `src/turbulence/saGammaRetheta.F90:541, 602,
2556, 2592` are hard `min(max(...))` clips inside the residual. Every cell on a
bound contributes a discontinuous Jacobian row, concentrated exactly along the
transition front. `smoothMinMax` already exists in `turbUtils.F90` and is not
used at those sites.

## 2. Source-term time stepping (§III.A) — the exact rule

Root cause (Lian et al.): implicit Euler applied to a **source** (positive
eigenvalue) is only well-behaved while the non-dimensional source-term time
step is below unity; **above unity the update flips sign**. That is the
mechanism behind loss of positivity of γ / Re̅θt.

Implementation in this paper:

- Form the block-diagonal source-term Jacobian (SA + transition source terms).
- Get its **largest positive eigenvalue** `λ_source` via a **QR algorithm**.
- `if λ_source·Δt ≥ 0.5 then Δt = 0.5 / λ_source`.

Note the **0.5** here; the thesis and the 2022 papers use **0.8** after further
tuning ("an effective balance between speed and robustness"). Our
`transitionSrcDtRestrict` is the same idea — worth checking our threshold and,
more importantly, that we take a genuine largest-positive-eigenvalue of the
coupled 3×3 source Jacobian, not a scalar proxy.

## 3. The smoothing machinery (§III.B) — reusable verbatim

Exponential penalty (log-sum-exp) smooth min/max, with a **switch** so the
expensive branch only runs near the kink:

```
if |g1(x) − g2(x)| > −log(|p|·p_mach)/|p|      then  φ_p = min/max(g1,g2)
else                                                 φ_p = log(Σ exp(p·g_i)) / p
```

- `p = +300` → smooth max, `p = −300` → smooth min. `p_mach ≈ 1e-15`.
- The switch exists for a real reason: near zero the exponentials produce
  **denormalised floats**, which "significantly slow the algorithm". If our
  `smoothMinMax` lacks this guard, it is both slower and less accurate than the
  paper's.

**Smooth `Flength` (Gaussian model)** — replaces the 4-branch piecewise
correlation, whose branch boundaries are kinks:

```
t1 = Σ_{i=1..3} a_i · exp(−((Re̅θt − b_i)/c_i)²)
t2 = 0.5 − (Re̅θt − 596.0)·3e-4
t3 = φ_300(t1, t2)
Flength = φ_300(t3, 0.3188)

a1=38.85,  b1=−13.94,  c1=273.1
a2=13.00,  b2=262.50,  c2=159.5
a3= 3.699, b3=375.90,  c3= 81.92
```

(The 2023 paper later replaces this with a single sigmoid-style expression,
Eq. A.34–A.35 — see [`MP_05`](MP_05_zingg2023_compressibility_corrections_digest.md).
If our code carries the Gaussian form, it is the 2019 vintage.)

## 4. Numerical-dissipation finding (crossflow only)

Relevant to our crossflow stall. On the NLF2-0415 infinite swept wing:

> "The **scalar dissipation model drastically under-predicted the crossflow
> strength and helicity values**… due to the increased dissipation of the
> scalar dissipation model in and above the laminar boundary layer, which
> affected the vorticity profile within the boundary layer. For turbulent
> boundary layers this increased dissipation did not have a significant
> effect."

Matrix-based dissipation (Swanson–Turkel) with `Vl = Vn = 0` and fourth-difference
coefficient 0.04 was used instead. **Crossflow prediction is a
dissipation-scheme-sensitive quantity**, and ADflow's scalar-dissipation
default is on the wrong side of this finding. Worth keeping in mind before
concluding anything about our own `D_scf`.

## 5. Model-form notes (this is the *old* model — do not port blindly)

- Eq. 43: `P̃_ν̃ = γ_eff · P_ν̃` — γ multiplies SA production only. Same as our
  rule 2. `γ_eff = max(γ, γ_sep)`.
- `F_onset2` / `F_onset3` were "modified from their original values" because the
  stock correlations "consistently predicted transition too far upstream" when
  coupled to SA. The 2020/2023 model supersedes this.
- Crossflow enters via `D_scf`, a sink in the Re̅θt equation (Eqs. 19–27). **The
  2023 paper replaces `D_scf` with an `F_onset,scf` term in the γ equation** —
  see MP_05. Our implementation is on the old side of that change.
- `Tu` is taken as the **freestream value throughout the domain** (local Tu is
  unavailable with SA); consistent with Medida & Baeder.

## Bottom line for the current stall

Two things in this paper are testable against our AR5 case without new physics:
**(i)** de-kink the residual at the four clip sites, **(ii)** verify our
source-term Δt restriction really uses the largest positive eigenvalue of the
coupled source Jacobian. The paper's own data says these — not the
preconditioner — are what moved convergence from 1e-5 to machine zero.
