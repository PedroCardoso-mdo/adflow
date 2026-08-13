# VERIF_06 — Solver fidelity audit (SA-GR vs Piotrowski & Zingg)

**Started 2026-08-13.** Restore point: git tag `pre-solver-fidelity-2026-08-13`
(working branch `solver-fidelity-2026-08`).

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

## F1 — CANK step-controller limit cycle: the CFL cutback is a no-op at depth

**Verdict: DIVERGES from the paper. This is the primary finding.**

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

## F2 — The unsteady line search lands just above the rejection threshold

**Verdict: DIVERGES (interacts with F1).**

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

## F6 — NK has no physicality check

**Verdict: ABSENT.** Confirmed at `NKSolvers.F90:860-862` (the code says so
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

1. **F3** — retest the swept wing with the validated ladder under matrix
   dissipation. Zero code, one job, and it may move the whole problem.
2. **F1 + F2** — the CANK step-controller limit cycle. One coherent correction.
3. **F4** — `.pyf` gap (2 lines, unblocks a knob we may want).
4. **F5** — `transitionSrcDtLimit` 0.8 A/B (no code).
5. **F6** — NK physicality check, queued behind the above.
