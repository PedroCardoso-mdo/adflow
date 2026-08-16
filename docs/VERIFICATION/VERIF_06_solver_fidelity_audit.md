# VERIF_06 — Solver fidelity audit (SA-GR vs Piotrowski & Zingg)

**Started 2026-08-13.** Restore point: git tag `pre-solver-fidelity-2026-08-13`
(working branch `solver-fidelity-2026-08`).

## VALIDATED (2026-08-16): NLF(2)-0415 swept wing, 1.9% mean error vs the paper

With `ANKCFLLimit 1e8` **and** the corrected Reynolds reference length, all
seven Re converge (relL2Drop 1.6e-10 to 2.1e-9, ~500-635 iterations) with
**matrix dissipation and crossflow on**, and the transition location matches
P&Z's §V.C curve:

| Re/1e6 | ours | paper | exp | vs paper | vs exp |
|---|---|---|---|---|---|
| 1.796 | 0.680 | 0.710 | 0.780 | −4.2% | −12.8% |
| 2.000 | 0.604 | 0.623 | 0.731 | −3.1% | −17.4% |
| 2.204 | 0.531 | 0.541 | 0.578 | −1.9% | −8.1% |
| 2.370 | 0.484 | 0.486 | 0.504 | −0.4% | −3.9% |
| 2.498 | 0.447 | 0.444 | 0.446 | +0.7% | +0.1% |
| 3.000 | 0.339 | 0.349 | 0.327 | −2.9% | +3.6% |
| 3.498 | 0.266 | 0.266 | 0.297 | +0.1% | −10.5% |

**Mean |error| vs paper 1.9%**, vs experiment 8.1%, on a cross-section 20%
coarser than theirs (65k vs 81,608 nodes).

The strongest single piece of evidence is not the mean: our largest deviation
from experiment is at Re 2.0e6 (−17.4%), and the paper reports "an error of
approximately 14%" at exactly that condition, attributing it to experimental
error in Radeztsky's data. **We reproduce not only their curve but their
specific disagreement with experiment, at the same point and magnitude.**

### The second bug: the Reynolds reference length

The runner used `--reynoldsLength 1.41421356`, on the campaign note that Re was
"quoted on the streamwise chord". Wrong. P&Z 2019 §IV.B describe the geometry as
"an NLF2-0415 airfoil **extruded** with a 45 degree sweep angle" and quote Re
"based on free-stream velocity magnitude and chord length" — the chord-1 airfoil
is the section NORMAL to the leading edge, so the reference length is 1.0.

With sqrt(2), asking for Re = 2.0e6 sets the freestream to
`rho*U*sqrt(2)/mu = 2.0e6`, i.e. an actual chord Reynolds number of 1.41e6:
**every run was a factor sqrt(2) below the paper.**

How it was caught, and why it is worth recording: comparing x_tr at the
*nominal* Re gave errors growing 10.8% -> 71.6%, which reads exactly like a
crossflow contribution that is too weak — a plausible, entirely wrong
diagnosis. Comparing instead at Re/sqrt(2), the chord Reynolds number actually
being run, gave −1.2% and +0.9% at the two points that land inside the
digitised range. A monotonically growing error is a signature of a wrong
*parameter*, not a wrong *model*; the transition model was never at fault.

## RESULT (2026-08-15): the binding constraint was the CFL ceiling

**The swept wing converges with matrix dissipation.** Job 1825801, Re 2.0e6,
crossflow on, `vis4 0.04`, `Vl=Vn=0`:

| variant | iters | wall | relL2Drop | converged |
|---|---:|---:|---|---|
| **`cfl8`** (`ANKCFLLimit 1e8`, linear tol left at 0.05) | **534** | **492 s** | 5.47e-09 | **True** |
| `lin4_cfl10` | 325 | 907 s | 6.28e-10 | True |
| `lin3_cfl8` | 1237 | 1609 s | 3.23e-10 | True |
| `cfl10` | 2727 | 3080 s | 2.88e-10 | True |
| previous best (any recipe) | 9282 | — | 8.47e-05 | False |

Four orders deeper and 17x fewer iterations, and the four agree on `cd` to the
fifth digit (0.007071–0.007076) — the recipe, not an artefact of one cell.

**Cause.** CFL sat pinned at `1.00E+06` because that was the `ankCFLLimit` the
runner passed: **saturated at the ceiling, not settling there.** Lifting it
shrinks the pseudo-transient term `T` and CANK becomes the paper's
inexact-Newton phase in fact rather than in name — which is exactly the F0
reframing, arrived at through a knob rather than through new code.

**Two things this retrospectively explains.**

1. Tightening the *linear* tolerance alone (`lin2`/`lin3`/`lin4`) only reached
   ~1e-5, while `cfl8` converged with the linear tolerance at its 0.05
   default. **The linear solve was never the constraint.** The "healthy
   lin res ~0.03" seen throughout was GMRES meeting a slack tolerance and
   exiting, not evidence of a good operator.
2. Why F1, F2, F8, the MFFD `h` and unit column scaling all fired exactly as
   designed and moved nothing. With CFL pinned at the ceiling, `Step` already
   1.00 and the residual flat, none of them *could* help: the binding
   constraint was elsewhere. Each was a correct mechanism aimed at a
   non-binding constraint.

**Standing correction to `CORE_02`:** its "do NOT couple early — CANK at 1e-2
stagnates" rule does not hold here; `ankCoupledSwitchTol 1e-2` is part of the
winning recipe on this mesh. That rule was measured on the 175k plain wing
with the CFL ceiling in place.

**Credit:** the CFL-ceiling question came from the user, as did the linear-
tolerance one that isolated it by failing.

---

## Scope and method

**Deliberately limited** (user instruction 2026-08-13: *"auditoria limitada, já
temos muitos 'maus' resultados do sweep, não gastes tempo e HPC à toa"*). Only
the rows that plausibly explain the observed stall, resolved from **code reading
plus logs that already exist** — no new HPC jobs were spent on this audit.

**Method rule: no verdict without evidence.** Every row quotes our code
(file:line) and the source (paper/thesis + equation or algorithm number). Where
a convention is genuinely ambiguous, the row says *ambiguous* and the resolution
is an experiment, never a guess. Rows not yet audited are listed as OPEN rather
than assumed.

Sources: `docs/MICHAEL_PIOTROWSKI/MP_03_*` (the solver paper),
`MP_06_*` §3.1 (thesis solver chapter, Algorithms 2/3/4 and Eq. 3.6–3.8),
`MP_01_*` (the 2019 origin paper).

## The symptom being explained

`13_nlf0415_swept`, NLF(2)-0415 infinite swept wing, L0 (64.8k cells), the case
that must converge. Two failure signatures across the 2026-08-11/12 rounds:

- **matrix dissipation** (`vis4 0.04`, `Vl=Vn=0` — the scheme P&Z 2019 §IV.B
  requires for crossflow): NK grinds at 4–5.5 orders.
- **scalar `vis4 0.01`**: converges 7–14 min for Re ≤ 2.5e6, but the mid-BL
  helicity is destroyed (`hcf_max` 0.063–0.071 < 0.1066), so the crossflow
  result is invalid; Re 3.0/3.5e6 stagnate in CANK regardless.

---

## F0 — Their "inexact-Newton phase" is still pseudo-transient. Ours is pure Newton.

**Verdict: DIVERGES, architecturally. This is the dominant finding — F1/F2 are
the same disease in the other phase.**

> Raised by the user 2026-08-13: *"os papers que te dei dele focam-se
> essencialmente no NK, esse mecanismo que estás a falar de CFL é no ANK; porque
> é que eles conseguem convergir tudo com NK e eu nem a começar a 1e-8 com as
> alterações deles consigo?"* — correct, and this row is the answer.

### What the authors' NK actually solves

Thesis §3.1, Eq. 3.3 — **both** phases solve the same system:

```
(T^(n) + A^(n)) ΔQ^(n) = −R(Q^(n))
```

> "where `T^(n)` is a diagonal matrix containing the **inverse local time steps**
> and `A^(n)` is the flow Jacobian"

and the second phase is reached not by removing `T`, but by making it small:

> "After the residual drops several orders of magnitude, the inexact-Newton
> method using the full Jacobian is used to converge the system to a residual
> norm of machine zero, where **the inexact-Newton method is recovered by
> aggressively ramping the reference time step, Δt_ref**."

So `T` never leaves their linear system. It is a continuation knob that is
turned down — and, critically, **can be turned back up at any moment and any
depth**. Algorithms 2 and 4 both end with exactly that move:

> `if Δt_ref > Δt_ref,min:  Δt_ref = max(0.5·Δt_ref, Δt_ref,min);  retry`
>
> "This makes the linear system **more diagonally dominant** due to the inverse
> time step in Equation 3.3, which improves convergence of the linear system and
> reduces the magnitude of the solution update."

Note `Δt_ref,min` is a **user-specified** bound, not a convergence-scaled one
(contrast F1).

### What ours solves

`NKSolvers.F90` — the NK KSP is only ever handed two operators:

```fortran
221:            call KSPSetOperators(NK_KSP, dRdw, dRdwPre, ierr)
473:            call KSPSetOperators(NK_KSP, dRdwNKSrcDt, dRdwPre, ierr)
```

- `dRdw` is the raw matrix-free FD Jacobian (`MatCreateMFFD`, `:179-183`) —
  **no `T` term at all**.
- `dRdwNKSrcDt` is that same MFFD wrapped with the Eq. 59 source-Δt diagonal
  (`NKSrcDtMatMult`, `:1576-1651`) — which touches **only the 3 turbulence /
  transition rows**, and switches off after `srcDtDeactivateIters` clean steps.

The one shell that *would* add a pseudo-transient diagonal to **all** rows is
`dRdwPseudo`, whose matmult is:

```fortran
291:        yPtr = yPtr + one / NK_CFL * xPtr
```

It is created at `:188`, given its shell operation at `:192`, destroyed at
`:546` — and **never passed to `KSPSetOperators`**. `NK_CFL` is dutifully
updated at `:644` (`NK_CFL = NK_CFL0 * (totalR0/norm)**1.5`) and consumed by
nothing. Verified by exhaustive grep: those are the only two
`KSPSetOperators(NK_KSP, …)` calls in the file, and the only `dRdwPseudo`
references are create/shell-set/destroy.

### Why this answers the question

**Their NK does 8–11 orders of pseudo-transient continuation; ours is asked to
do the last fraction of an order of pure Newton.** Two different algorithms
sharing a name:

| | P&Z inexact-Newton | ADflow NK |
|---|---|---|
| Engages at | rel **1e-4** (2D) / **1e-5** (3D) | validated ~**4e-8**; 1e-6 falsified as premature (`CORE_02`) |
| Runs to | rel 1e-15 / abs 1e-10 (machine zero) | walls below ~5e-9 |
| Diagonal term | `T = 1/Δt_ref`, all rows, always present | none on mean-flow rows; Eq. 59 on turb rows only, and it deactivates |
| Recovery from a bad step | halve `Δt_ref` and **retry** the solve | shrink λ; at `minlambda = 0.01` "just take it" (documented E1, `SAGR_02` §7) |
| Physicality check | Algorithm 2 on ρ/E | **none** (F6) |

So porting "their NK modifications" (Eq. 59 in NK, Algorithm 2 damping on
γ/Re̅θt) onto a **pure** Newton solve was always going to under-deliver: the
piece those modifications lean on — a restorable `T` — is precisely the piece
ADflow's NK does not have. That is why engaging our NK at their engagement point
(1e-5) stalls, and why even at 4e-8 it walls.

It also unifies this audit: **in both phases the mechanism that restores
diagonal dominance after a bad step is unavailable** — dead in NK (this row),
floored out by the convergence-scaled `ANK_CFLMin` in ANK/CANK (F1).

### REVISION 2026-08-13 — their inexact-Newton is our **CSANK**, not our NK

> User: *"isto tudo não estamos a tornar o NK no ANK? tipo o solver deles que
> chamam Newton já não é +/- o nosso CANK?"* — **yes, and this kills the
> original action below.**

Verified in code (`NKSolvers.F90:5148-5178`): the `*SANK`/`*CSANK` iteration
types are selected when `totalR <= ANK_secondOrdSwitchTol*totalR0`, which is
exactly the regime where `ANK_useDissApprox`, `lumpedDiss` and `approxSA` are
**not** applied — i.e. the true second-order residual — while
`computeTimeStepMat` still contributes the `T` term every iteration.

So:

| Thesis | ADflow |
|---|---|
| approximate-Newton: `(T + A_approx) ΔQ = −R` | **ANK / CANK** |
| inexact-Newton: `(T + A_2nd) ΔQ = −R`, `Δt_ref` ramped up so `T → 0` | **SANK / CSANK** |
| — (the `T = 0` limit) | **NK** |

**ADflow's NK is the `T → 0` endpoint of a continuum that CSANK already
traverses.** Adding a pseudo-transient diagonal to NK, as the original action
proposed, would be re-deriving CSANK — new code for a solver we already have.

**Consequence — the plan changes:** do **not** build `T` into NK. Instead make
**CSANK** do what their inexact-Newton does: engage around rel 1e-5 and run to
convergence. That reframes the target of every remaining item — F1/F2 are
CSANK's step controller, and fixing them fixes the equivalent of their
inexact-Newton phase directly.

This also explains a symptom `CORE_02` already recorded but never accounted
for: *"Do NOT hold … CSANK past ~3.5e-8: a front adjustment kicks the residual
and the phase enters permanent oscillation (never recovers its best depth)."*
**Permanent oscillation with no recovery is precisely the F2 limit cycle** — the
controller cannot classify the collapsed step as a failure, so it never backs
off. It is the same bug, seen one phase earlier.

`useNKSolver: False` + letting CSANK finish — already the recommendation in
`CORE_02` for AR5 L0 — is therefore not a workaround. **It is the paper's actual
algorithm.**

### Action — SUPERSEDED by the revision above

Install a pseudo-transient diagonal on the NK operator and give it a
back-off-and-retry on step rejection, with a **user-specified floor**
(`Δt_ref,min` equivalent), behind an option that defaults to today's behaviour.
The shell (`dRdwPseudo`/`NKMatMult`) and the ramp law (`NK_CFL0·(totalR0/norm)^1.5`)
already exist — this is re-activating and completing existing machinery, not
inventing a mechanism. The matching diagonal must also be added to `dRdwPre`
(e.g. `MatShift`) or the PC no longer preconditions the operator being solved.

**Caveat to respect:** upstream ADflow's NK is deliberately pure Newton and works
that way for fully-turbulent SA, where the problem is far less stiff. The option
must therefore be scoped to the SA-GR path and default off.

---

## F1 — CANK step-controller limit cycle: the CFL cutback never fires

> **CORRECTION 2026-08-13, same day.** This row originally claimed the cutback
> was blocked because the ramped floor `ANK_CFLMin` had reached
> `ANK_CFLLimit`. **That is wrong at the depths where the swept wing actually
> stalls**, and the arithmetic proves it: at the stalled iterations
> `totalRes/totalR0 ≈ 9.4e-6`, so
> `ANK_CFLMin = 1.0·(1/9.4e-6)^0.5 ≈ 3.3e2` — three orders below the 1e6 the
> CFL is pinned at. The cutback had ample room (`max(1e6·0.5, 3.3e2) = 5e5`).
>
> The floor-saturation trap is **real but latent**: it needs
> `(totalR0/totalR)^0.5 ≥ ANK_CFLLimit`, i.e. rel ≈ 1e-12 for a 1e6 limit. It
> plausibly explains the deep-restart strandedness `CORE_02` records ("raise
> `ANKCFLMin` if stranded"), and `ANKCFLMinCap` (commit 6c31217d) remains a
> correct guard for it — but it is **not** this stall's cause. The real
> mechanism is F2, below, and the two rows should be read together.

**Verdict: DIVERGES from the paper — but via F2's mechanism, not the one first
claimed.**

### Evidence — our runs

`logs/nlfsw45_L0_a-4_re3p00_lad.log`, CANK phase, every 200th iteration:

```
iter   type     CFL       Step   LinRes
 217   CANK   6.21E+05    1.00   0.042
 418  *CANK   1.00E+06    1.00   0.017
 618  *CANK   1.00E+06    0.12   0.006
 818   CANK   1.00E+06    0.10   0.050
1018   CANK   1.00E+06    0.22   0.046
1218   CANK   1.00E+06    0.06   0.042
1418   CANK   1.00E+06    1.00   0.041
1618   CANK   1.00E+06    0.50   0.021
1818   CANK   1.00E+06    0.06   0.025
```

and `logs/nlfsw45_L0_a-4_re2p00.log` (matrix dissipation), iterations 6–42:

```
   6   CANK   4.37E+03  0.02  0.047   ...  totalRes 3.2790E+03
  ...
  42   CANK   4.37E+03  0.01  0.041   ...  totalRes 4.2975E+03   <- RISING
```

**The linear solve is healthy throughout — lin res 0.006–0.050.** The
preconditioner is not the limiter here. `CFL` is pinned at its ceiling while
`Step` collapses, and the residual creeps or rises.

### Evidence — our code

`src/NKSolver/NKSolvers.F90:4980-4987`:

```fortran
            ! First of all, update the minimum cfl wrt the overall convergence
            ANK_CFLMin = min(ANK_CFLLimit, ANK_CFLMinBase * (totalR0 / totalR)**ANK_CFLExponent)

            ! Update the CFL number depending on the outcome of the last iteration
            if (lambda < ANK_stepMin * ANK_stepFactor) then
                ! The step was too small, cut back the cfl
                ANK_CFL = max(ANK_CFL * ANK_CFLCutback, ANK_CFLMin)
```

Defaults (`pyADflow.py:6059-6067`): `ANKCFLMin 1.0`, `ANKCFLExponent 0.5`,
`ANKCFLCutback 0.5`, `ANKStepMin 0.01`, `ANKStepFactor 1.0`.

**The trap:** the cutback floor `ANK_CFLMin` is itself *ramped by convergence* —
`ANKCFLMin·(totalR0/totalR)^0.5`. By the time `totalR0/totalR ≈ 1e12`, the floor
has reached `ANK_CFLLimit` (1e6 in these runs). At that point
`max(ANK_CFL*0.5, ANK_CFLMin)` **cannot reduce `ANK_CFL` at all** — the only
recovery mechanism the coupled step controller has is mathematically dead
exactly when it is needed. This matches the logs: CFL locked at 1.00E+06 for
1800+ iterations while `Step` oscillates.

*(`CORE_02` already records the floor formula as a restart annoyance — "raise
`ANKCFLMin` if stranded" — but not that it neutralises the cutback.)*

### Evidence — what the authors do instead

Thesis §3.1.3, Algorithms 2 and 4: a step that needs a damping factor below 1%
is **rejected**, and retried with the reference time step **halved**, against a
**user-specified** lower bound `Δt_ref,min`:

> `if δ_phys < 0.01: if Δt_ref > Δt_ref,min: Δt_ref = max(0.5·Δt_ref, Δt_ref,min);`
> `δ_phys ← 0; return  (retry solution)`
>
> "This makes the linear system **more diagonally dominant** due to the inverse
> time step …, which improves convergence of the linear system and reduces the
> magnitude of the solution update."

Their floor is a fixed user bound; ours auto-ramps to the ceiling. **That is the
divergence.**

### Action
Make the effective cutback floor independent of the convergence ramp when the
step is collapsing, so CFL can actually come back down. Behind an option, with
today's behaviour as the revert path. **This is the highest-priority correction.**

---

## F2 — The line search cannot reach its own rejection threshold

**Verdict: DIVERGES. Promoted 2026-08-13 to the actual cause of the swept-wing
CANK stall** (see the correction in F1).

### The arithmetic

- Unsteady-LS backtracking: starts at `0.7·λ`, then up to **12** further
  factors of 0.7 ⇒ smallest reachable step **`0.7^12 = 0.0138`**.
- Rejection threshold for the CFL cutback: `ANK_stepMin·ANK_stepFactor` =
  **`0.01`** (defaults `0.01 × 1.0`).

**0.0138 > 0.01.** The line search structurally cannot produce a step small
enough to be classified as a failure. The only route to `λ < 0.01` is the
physicality check — and in the coupled branch the turbulence/transition
variables are clipped per-cell to `ratio = one` (see F6), so only ρ or E can do
it. Hence: **the cutback almost never fires, `ANK_CFL` stays pinned at
`ANK_CFLLimit`, every iteration proposes an enormous step, the line search
throttles it back to ~1–3%, and the cycle repeats.** Measured: CFL locked at
1.00E+06 for 1800+ iterations, `Step` 0.01–0.30, lin res healthy 0.006–0.050.

### What the authors do

Thesis Algorithm 4: backtrack by **0.90** (not 0.7) with an explicit **1%
floor**, and on reaching that floor **reject the step and halve `Δt_ref`**.
From λ=1 a factor of 0.90 needs ~44 backtracks to reach 0.01 — so their search
is both *gentler* and allowed to go *further*, and its floor coincides with the
rejection threshold by construction rather than sitting above it.

Two further differences already noted: their acceptance test is against the
*current steady* residual `‖R(Q^n)‖₂`, ours against `unsteadyNorm_old ×
ANK_unstdyLSTol`; and the `ANK_stepMin` floor check at `:5176` is applied to the
*pre*-line-search physicality λ only, never to the accepted λ.

### Action
Make the backtrack factor and budget runtime options so the thesis's 0.90/1%
geometry is reachable, and treat an exhausted budget as a rejection. Both
default to today's 0.7/12 so nothing changes until asked.

`NKSolvers.F90:5217-5248` — the coupled-ANK backtracking loop starts at
`lambda = 0.7*lambda` and runs at most **12** iterations, each `lambda*0.7`:

```fortran
            lambda = 0.7_realType * lambda
            backtrack: do iter = 1, 12
                ...
                if (unsteadyNorm > unsteadyNorm_old * ANK_unstdyLSTol .or. myisnan(unsteadyNorm)) then
                    ...
                    lambda = lambda * 0.7_realType
```

`0.7^12 = 0.0138`. The rejection threshold at `:4984` is
`ANK_stepMin * ANK_stepFactor = 0.01`. **The backtrack floor sits above the
rejection threshold**, so a fully-exhausted line search still reports "success"
at λ ≈ 0.014 and is never counted as a failed step — no CFL cutback, no
`lambda = zero`. The observed `Step 0.01–0.02` is exactly this floor.

Two further differences from thesis Algorithm 4:
- **Acceptance test.** Ours compares the new unsteady norm against
  `unsteadyNorm_old * ANK_unstdyLSTol` (the *previous unsteady* norm, tolerance
  2.0 in these runs). Theirs requires the unsteady residual
  `‖TΔQ + R(Q+ΔQ)‖₂` to fall below `‖R(Q^n)‖₂` — the *current steady* residual.
- **Backtrack factor** 0.7 (ours) vs 0.90 (theirs), and theirs has an explicit
  1% floor that triggers the reject-and-retry of F1.

Note the physicality λ *does* get a floor check (`:5176`,
`if (ANK_CFL .gt. ANK_CFLMin .and. lambda .lt. ANK_stepMin) lambda = zero`) —
but that is applied **before** the line search, not to the post-line-search λ.

### Action
Apply the `ANK_stepMin` rejection to the post-line-search λ as well, so an
exhausted backtrack is treated as the failure it is. Pairs with F1.

---

## F3 — The swept-wing runs did not use the validated ladder

**Verdict: not a code defect — a run-configuration divergence. Free to retest.**

From `logs/nlfsw45_L0_a-4_re2p00.log` modified options:

```
 'ANKCoupledSwitchTol': 0.01,
 'ANKSecondOrdSwitchTol': 0.001,
 'NKSwitchTol': 0.0005,
```

`CORE_02` (validated on the 175k 3d_plain_wing, 2026-07):
- **"Do NOT couple early: CANK at rel 1e-2 stagnates … the front is too
  unconverged and the global lambda throttles the whole field."** Validated
  coupling point: **1e-5**. These runs coupled at **1e-2**.
- **"Engaging NK too early is the most expensive mistake"** — validated NK
  engagement is at CSANK's floor, rel ~4e-8; `1e-6` was explicitly falsified as
  premature. These runs engaged NK at **5e-4**, three orders earlier still.

The later `_lad` round (CANK 1e-3 / CSANK 1e-5 / NK 1e-6) did markedly better —
Re 1.8–2.5e6 converged in 7–14 min — which is consistent with the tolerances
being a first-order effect here, and it is still far from the validated recipe.

### Action
**Before any code change, retest the swept wing with the validated ladder**
(CANK 1e-5, CSANK 1e-6, NK ~4e-8, `ANKUnsteadyLSTol 1.5`, `ANKPhysicalLSTol 0.5`)
under **matrix dissipation**. One job. If it converges, F1/F2 are still real but
may not be on this case's critical path.

---

## F4 — `.pyf` gap: two ANK options are silent no-ops

**Verdict: DIVERGES (live defect) — but NOT the cause of this stall.**

`ANK_physLSTolReTheta` (`NKSolvers.F90:2381`) and `omegaMinGamma` (`:2382`) are
consumed by `physicalityCheckANKTurb` (`:4437`, `:4450`, `:4463`) — the routine
that computes the turbulence/transition λ. The Python side is complete:
defaults at `pyADflow.py:6071-6072`, optionMap at `:6521-6522`. But
**`src/f2py/adflow.pyf` module `anksolver` declares neither** (verified by
reading the whole module block).

Per the documented f2py pathology (`SAGR_02` §5), f2py `fortran`-type module
objects silently accept arbitrary attribute names, so `setOption` succeeds,
reads back the value, and the compiled code keeps its declared default.

**Honest scoping:** the Fortran defaults (0.99, 0.05) equal the Python defaults,
and `run_swept.py` never sets either — so the runs used the intended values and
this did *not* cause the stall. It does mean the knob is unusable, and we may
want it during this campaign.

### Action
Add both to the `.pyf` (2 lines + f2py re-wrap). Low risk, own commit.

---

## F5 — `transitionSrcDtLimit` 0.9 vs the papers' 0.8

**Verdict: DIVERGES (parameter only, no code change needed).**

Ours: `inputParam.F90:393` `transitionSrcDtLimit = 0.9_realType`, mirrored at
`pyADflow.py:5900`. Thesis §3.1.4 Eq. 3.14 and MP_03 §II.D both use
**0.8** ("several source-term time step values were investigated; a value of 0.8
was determined to be an effective balance between speed and robustness").
The 2019 paper used 0.5.

Already a runtime option ⇒ pure A/B, no commit required to test.

---

## F6 — NK has no physicality check (ANK/CANK does)

> **Clarification (user, 2026-08-13: *"acho que temos physicality check no nosso
> solver"* — correct).** We DO have one, in ANK: `physicalityCheckANK`
> (`NKSolvers.F90:4064`) limits ρ and E to a relative change of
> `ANK_physLSTol`, and its `ANK_coupled` branch (`:4155`) additionally covers
> ν̃, γ and Re̅θt *with* the SA-GR bounds. So it **is** active in the phase where
> the swept wing stalls. The gap is NK only.
>
> A consequence worth recording, found while checking this: in the coupled
> branch a turbulence/transition ratio below `ANK_stepMin*ANK_stepFactor` is
> **clipped per-cell and the ratio set to one** (`:4207-4219`, `:4272-4276`), so
> γ/Re̅θt/ν̃ can never throttle the global λ below `ANK_stepMin`. **Only ρ and E
> go straight into the global `min`.** Therefore a collapsed CANK step is either
> ρ/E physicality or the unsteady line search — never the transition variables.
> That is exactly the distinction the STALLDIAG output was added to settle.

**Verdict for NK: ABSENT.** Confirmed at `NKSolvers.F90:860-862` (the code says so
itself) — there is no ρ/E positivity check anywhere in `NKStep`/`LSCubic`/`LSNM`;
the only safety valves are NaN detection and the post-step γ/Re̅θt clamp in
`applyNKAlgorithm2Damping`. This is why `LSCubic`'s Armijo `alpha` was
hand-relaxed 1e-2 → 1e-3 after a 1e-4 attempt SEGV'd.

Thesis Algorithm 2 specifies it: negative ρ/E updates limited to 90% of the
current value, global `MPI_Allreduce(MIN)`. Template already in-tree:
`physicalityCheckANK` (`:4064-4314`).

Already recorded as `SAGR_02` §8 item 2; this audit confirms it unchanged.
**Not on the critical path for the CANK stall** (F1/F2 are), so it stays queued
behind them.

---

## F7 — Which paper items are in NK vs CANK/CSANK: each has half of Alg. 2/4

> User, 2026-08-13: *"temos coisas implementadas no NK que vieram do paper,
> também elas estão implementadas no CANK?"* — the necessary follow-up to F0's
> revision. If CSANK is the paper's inexact-Newton phase, the ported machinery
> has to be there, not only in NK.

Checked by call-site enumeration, not by assumption:

| Paper item | NK | CANK / CSANK |
|---|---|---|
| Eq. 59 source-Δt diagonal | `applyNKSrcDtDiagonal` (`:1514`) + `NKSrcDtMatMult` (`:1576`) | `dtInvSrc` in `computeTimeStepBlock` (`:3029-3045`) |
| Eq. 59 `srcLambda` freeze | `:464` | `:5090` (coupled), `:4678`/`:4772` (turb) |
| Eq. 59 deactivation counter | `:795-799` | `:5330-5338` — **plus the residual-rise leg NK lacks** (F5/`:4947`) |
| Eq. 58 `S_r` row scale / `S_a` autoscale | `setRVec` | **the same `setRVec`** (`:5454`, `:5034`) |
| Eq. 58 `S_c` column scale | `getNKColScale` (`:1476`) | `getFullColScale` (`:3991`) / `getTurbColScale` |
| **Alg. 2 per-node γ/Re̅θt damping** | `applyNKAlgorithm2Damping` (`:760`) | **ABSENT** |
| Alg. 2 physicality on ρ/E | **ABSENT** (F6) | `physicalityCheckANK` (`:5175`) |
| Alg. 4 unsteady-residual line search | ABSENT (cubic Armijo instead) | present — it *is* the unsteady LS |

**Verdict: the two solvers hold mirror-image halves of the paper's
Algorithm-2/4 pair, and neither holds both.** `applyNKAlgorithm2Damping` has
exactly one call site in the whole file (`:760`, inside `NKStep`); there is no
`ANKStep` equivalent.

**Why the gap exists, and why it deserves a re-test.** Per-cell Algorithm-2
damping *was* tried in `physicalityCheckANK` and deliberately reverted — the
comment at `:4241-4247` records that locally damped front cells made the update
inconsistent across the transition front, the γ residual bounced, and it lost to
the global-λ throttle. **But that comparison was run before F2**, i.e. against a
global-λ controller that we now know was in a limit cycle it could not exit. The
baseline it lost to was broken. Re-test after F2 lands, not before.

**Minor inconsistency found in passing:** `getNKColScale` gates on
`transitionNK .and. transitionNKActive`, `getFullColScale` gates on
`transitionNK` only. If the NK auto-disable latch ever trips, NK's column
scaling switches off while CANK's stays on. Inert today
(`transitionNKAutoDisableTol = 0` never trips it) but worth aligning.

## Confirmed-MATCH rows (checked, no action)

| Element | Source | Ours | Verdict |
|---|---|---|---|
| Vorticity limiter in γ sources | MP_05 App. A `φ_−300(Ω, U∞M∞√Re∞/(l·20))` | `saGammaRetheta.F90:660-669`, smooth `smoothMinMax(vortMag, vortLim, rsaGRpmin)`, mirrored in `evalSrcJacBlock:2624-2631` | **match** (also rotating-frame correct) |
| Smooth min/max + proximity switch | MP_01 §III.B Eqs. 46-50, p=±300 | `turbUtils.F90:2397-2440`, `rsaGRpmin/pmax = ∓300`, switch present | **match** |
| γ/Re̅θt damping bounds | Thesis Alg. 3: γ∈(1e-10,2), θ=0.99 | `paramTurb.F90:73-74` `rsaGRgammaLo=1e-10`, `rsaGRgammaHi=2.0`; `transitionDampTheta=0.99` | **match** |
| Eq. 59 largest positive eigenvalue | Thesis Eq. 3.13, QR on the 3×3 | `computeSrcLambda`, `saGammaRetheta.F90:2725-2799` — analytic, exploits block-triangular structure (A13=A31=A32=0), handles the complex-conjugate branch | **match** (analytic shortcut, equivalent) |
| Eq. 59 reactivation on backtrack, coupled path | MP_03 §II.D | `NKSolvers.F90:5330-5338` | **match** |
| First-order upwind on turb/transition convection | MP_03 §II.A | project-wide (CLAUDE.md rule 5) | **match** |

---

## OPEN rows (not audited — do not assume)

- **B2 reported residual.** Whether our `totalRes`/`totalR0` is the `S_r`-scaled
  quantity the thesis monitors (`Rd = ‖S_r R‖₂/‖S_r R⁰‖₂`, Eq. 3.5). `totalR0`
  comes from `getFreeStreamResidual` (`solvers.F90:972`), which does not route
  through `setRVec`. If they differ, our "orders" and the papers' are not the
  same number.
- **B2 `S_r` exponents.** `setRVec:1953-1954` gives net `volRef^(2/3)` (flow) /
  `volRef^(-1/3)` (turb), which matches thesis Eq. 3.7 **iff `J ≡ volRef`** and
  is inverted iff `J ≡ 1/volRef`. **Ambiguous — resolve by experiment, not by
  argument.** Also: `setRVec` and `computeNKResidualAutoscale:1828-1835` hold two
  independent copies of the same literals.
- **B1 coupling structure per phase** vs MP_03's finding that partial coupling
  stalls at 4–5 orders on hard 3D cases. The DADI banner in the swept logs reads
  `fully coupled (3x3 block)`, but what CANK/CSANK/NK each linearize has not been
  traced end-to-end here.
- **B5/B6** deliberately out of scope for this limited pass (already covered in
  `SAGR_02`).

---

## Priority order out of this audit

Revised 2026-08-13 after F0. The user has authorised rewriting the step
controller as an option, **step by step, with a verification test after each
step that shows whether the change produced the intended effect — and if not,
why**. That cadence is part of the plan, not an afterthought: each step below
names its own falsifiable check.

| # | Item | Code? | Verification that it worked |
|---|---|---|---|
| 1 | **F3** — retest swept wing with the validated ladder (CANK 1e-5, CSANK 1e-6, NK ~4e-8, LS 1.5/0.5) under **matrix** dissipation | none | Does it converge at all? Establishes whether F0/F1 are on this case's critical path before any code is written. |
| 2 | **F0a** — pseudo-transient diagonal on the NK operator **and** the PC, behind an option | yes | (a) default-off ⇒ bit-identical; (b) option on with `NKCFL0` huge ⇒ `T→0` ⇒ results ≈ unchanged (proves the term is wired but inert); (c) option on with a modest `NKCFL0` ⇒ NK linear residual should *drop* at fixed state (more diagonally dominant). If (c) fails, the diagonal is not reaching the operator actually solved — check `dRdwPre` got the matching shift. |
| 3 | **F0b** — reject-and-retry: on a collapsed step, halve `Δt_ref` against a **user** floor and re-solve | yes | `Step` should stop pinning: the histogram of `Step` in the NK phase must lose its spike at `minlambda`. If it does not, the rejection threshold is above the line-search floor — the F2 failure mode, one phase over. |
| 4 | **F1 + F2** — the same two fixes in ANK/CANK: decouple the cutback floor from the convergence ramp, and apply `ANK_stepMin` to the post-line-search λ | yes | CFL must be observed *coming back down* in the log while `Step` recovers. Today it never does (evidence in F1). |
| 5 | **F4** — `.pyf` gap (2 lines) | yes | `setOption` to a deliberately absurd value must now visibly change behaviour; today it cannot. |
| 6 | **F5** — `transitionSrcDtLimit` 0.8 vs 0.9 | none | A/B iteration count. |
| 7 | **F6** — NK physicality check (Algorithm 2 on ρ/E) | yes | Allows `LSCubic`'s `alpha` to be tightened back toward 1e-2 without the SEGV that forced 1e-3. |

Items 2–4 are one coherent theme — **restore the ability to put diagonal
dominance back after a bad step** — and are expected to be tested together on
the swept wing once each is individually verified.
