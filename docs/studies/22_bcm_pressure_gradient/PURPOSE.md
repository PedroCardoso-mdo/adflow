# 22 — SA-BCM with a local pressure-gradient sensor (SABCM_PG)

**Started:** 2026-09-18  **Status: CLOSED 2026-09-19** — sensor 2 implemented, audited and its adjoint verified
(branch `sa-bcm-pg` @ 9b39601e, machV4), but the user's verdict on the results stands: *not a net improvement*
(ranking of the Tu-0.5 % optima fixed, NLF polar degraded, GR C2 unchanged). The C2/C3 BCM-PG optimisation
(job 1936872) was cancelled at SLSQP major 1. `SABCM_PG` stays in the code as an option (default off).

## 1. Why

Cross-model evaluation on the NACA0012 L0 mesh (Tu 0.5 %, cl 0.3, run tree
`08_optimization/2d_tu05_L0/cross_eval/PURPOSE.md`): the SA-BCM optima are *worse than the NACA*
when analysed with SA-γ-R̅eθt (+5.5 % / +14.6 %), while the GR optima keep −25 % / −42 % in SA-BCM.
cp is identical in both models for a given shape; only the cf jump location differs. The BCM
threshold `Re_θc = 803.73 (Tu + 0.6067)^-1.027` depends on Tu∞ only, so the optimiser drives the
laminar run into adverse gradients (nose suction peak, mid-chord dips) that the BCM never
sanctions and the GR model does through `Re_θt(Tu, λ_θ) = Re_θt(Tu)·F(λ_θ)`.

## 2. Model (what was implemented)

Per cell, in `saSource` (`src/turbulence/sa.F90`, `use_SABCM` block; mirrored in
`src/NKSolver/blockette.F90`), when `SABCM_PG = True`:

    n      = nWall(:,i,j,k)                      unit vector nearest-wall-point -> cell centre
    Snn    = n_i n_j du_i/dx_j                    (velocity gradients already in saSource)
    lamL   = -coef * (d^2/nu) * Snn + off         Menter et al. 2015 lambda_thetaL (coef 7.57e-3)
    lam    = smooth clip of gain*lamL to [-lamMax, lamMax]      (smoothMinMax, p = 300, lamMax 0.1)
    F      = smoothMin( smoothMax(F1(lam,Tu), 1), F3(lam,Tu) )  Langtry-Menter Eqs. 54-57 (as in SA-GR)
    Re_theta_c = Re_theta_c^BCM * F

Everything else (Re_θ = Re_v/2.193, KS-max Term1, Term2 = R_T/χ2, γ = ½(1+tanh((T1+T2−S0)/fsmooth)))
is untouched; with `SABCM_PG = False` the code path is bit-identical to `sa-bcm` (the multiplication
sits inside the `if`). No Reynolds factor: `nu = rlv/rho` is non-dimensional with L_ref = 1 m, the
same convention as `Re_v` (the note in `docs/nondimensionalization.md` claiming an explicit 1/Re is
wrong for this purpose).

* `nWall` is produced exactly in `updateWallDistancesQuickly` (`src/wallDistance/wallDistance.F90`)
  from the foot point `xp` and the cell centre `xc`, on the approximate wall-distance path only
  (`checkInputParam` enforces `use_SABCM ∧ useApproxWallDistance` for `SABCM_PG`); it is a new
  block array (`block.F90`, `blockPointers.F90`, `utils.F90`, `adjointUtils.F90` twins of `d2Wall`),
  and a new independent of both Tapenade heads (`Makefile_tapenade`: `updateWallDistancesQuickly`,
  full `saSource`) so the adjoint sees ∂λ/∂n·∂n/∂x.
* Helpers `smoothMinMax` (ported from the SA-GR `turbUtils.F90`) and `bcmFlambda` live in
  `src/turbulence/turbUtils.F90` (outside the `USE_TAPENADE` guard, bare `use constants`).
* Diagnostics in the volume CGNS: `BCM_lambda` (= lamL, before gain/clip) and `BCM_Flambda`
  (`volumeVariables: ["bcmlambda", "bcmflam"]`).
* Options (`pyADflow.py`): `SABCM_PG` (False), `SABCM_PG_coef` (7.57e-3), `SABCM_PG_off` (0.0128),
  `SABCM_PG_gain` (1.0), `SABCM_PG_lamMax` (0.1), `SABCM_PG_p` (300). Runs set off/gain explicitly
  (§3) — the Python defaults are Menter's raw values.
* Side fix: the blockette wrote `Tgamma(i,j,k)` with blockette-local indices (block arrays need the
  `(ii,jj,kk)` offset) — only the CGNS `Tgamma` with `useBlockettes=True` was affected; fixed.

Decision log: the handoff's Ψ/tanh(k, δ) intermittency was **not** adopted — the sa-bcm smooth form
(KS-max + tanh/fsmooth 0.08, studies 17–21) already removes the kinks it targets, and keeping it
means the cross-model change measures one effect only. Menter's own `F_PG` correlation is not
implemented (weak at low Tu; F of Langtry is what the reference model uses).

## 3. Calibration: Falkner–Skan (`falkner_skan_gain.py`, no CFD)

For each β the similarity profile gives, at η* = argmax(η² f'') (the Re_v peak where Term1 fires):
λ_θ exact (2m θ̂²/(m+1)) and λ_θL (with dV/dy = −du/dx from the similarity solution).

| result | value |
|---|---|
| Blasius: η* | 2.09 (FS scaling; ≈ 2.95 in Blasius scaling) |
| Blasius: λ_θL at η*, offset 0 | −0.0162 |
| … with Menter's offset 0.0128 | −0.0034 (F ≈ 0.96 at Tu 0.1 %: **not neutral on a flat plate**) |
| **calibrated offset** (λ_θL = 0 for Blasius at η*) | **0.01623** |
| linear gain λ_θ ≈ gain·(λ_θL + off), LS over \|λ_θ\| ≤ 0.1 | **2.04** (adverse-only 1.37, favourable-only 3.03) |

The map is mildly nonlinear (`falkner_skan_gain.png`): gain 2 is exact for mild gradients
(|λ| ≲ 0.03, the suction-peak regime), over-predicts near separation (λ_θL+off = −0.073 → −0.15,
clipped at −0.1 vs exact −0.068; F3 differs by 0.025 there) and under-predicts strong favourable
gradients (F1 saturates anyway). Run defaults: `off = 0.01623`, `gain = 2.0`; gain swept on the
cross set (§4 C3).

## 3b. Sensor 2 — wall-pressure gradient (2026-09-18, after the low-Tu convergence findings)

The velocity-profile sensor (Menter dV/dy) has a feedback loop: transition → boundary-layer growth →
∂v/∂y at the Re_v peak → F → transition location. On NLF0416/S809 L1 at Tu 0.15/0.07 % it limit-cycles
in ANK even when started from the converged PG-off state (job 1935764) — the PG fixed point is
unstable. A sensor that does not see the velocity profile removes the loop:

    lambda = K (d²/ν) dU_e/ds,   dU_e/ds = −(ŝ·∇p)/(ρ∞ U_e),   U_e² = U∞² + 2(p∞ − p)/ρ∞ (Bernoulli, local p ≈ p_wall)
    ŝ = velocity direction, p from w (constant cp), K = (θ̂/η*)²_Blasius = 0.05064

`falkner_skan_psensor.py`: with this single similarity constant the sensor reproduces the exact
Falkner–Skan λ_θ to +2/−4 % for |λ_θ| ≤ 0.03 (0.70 at separation, where the ±0.1 clip acts anyway);
no offset (U_e' = 0 for Blasius), no gain, no map — a calibration on the similarity family, no CFD
or target data involved. Option `SABCM_PG_sensor = 2` (`SABCM_PG_K`); sensor 1 keeps
`SABCM_PG_map = 2` = the analytic Falkner–Skan map λ_θL → λ_θ (`falkner_skan_map.py`: adverse
−A tanh(−u/B), favourable C(exp(u/D) − 1); A 0.0691, B 0.0326, C 0.0991, D 0.0423, representation
error 7e-4) as the "serious" replacement of the linear gain. Tapenade regenerated for both
(traps met: >132-char lines and long parameter declarations are rejected/mis-wrapped by Tapenade).

## 3c. Loop suavização → verificação (`08_optimization/2d_tu05_L0/pg_loop`, from 2026-09-18 evening)

Per iteration (1 node): adjoint vs complex step on NACA0012 L2 (α, DV 8, DV 12; gate = the PG-off
level, cl 0.008 % / cd 0.26 %) and cold L0 trims of the shapes where the BCM disagrees with GR
(BCM C2 r1 e50, BCM C3 e118, GR C2; GR 79.50 / 71.18 / 34.94, BCM 53.55 / 46.39 / 51.38 counts).
Knobs per iteration are smoothing/formulation only (p, lamMax, ν̃∞, sensor form); no constant is
tuned to these results.

| it | sensor / knobs | adjoint vs CS | BCM C2 r1 | BCM C3 e118 | GR C2 | note |
|---|---|---|---|---|---|---|
| 0 | sensor 2 (K 0.0506, p 300, lamMax 0.1, ν̃∞ 1e-7) | (job died on an int/float option; rerun) | 83.84 · 0.08/0.34 | 88.44 · 0.08/0.14 | 52.12 · 0.75/0.29 | all trims converged in 12–13 solves (no stalls); e118 lower surface over-penalised |
| 0 | (adjoint/CS rerun, job 1936308) | α 1.0 % · DV8 10 % · DV12 35 % (cd); cl ≤ 0.3 % | | | | partials all pass (reg suite); gap = primal branch sensitivity |
| 1 | sensor 3 (θ = Re_θc ν/U_e, no constant) | pending (1936311) | 91.4 · 0.07/0.12 | 89.7 · 0.07/0.11 | 52.7 · 0.75/0.28 | trips both surfaces at 0.07–0.12 c |
| 2 | sensor 2 + fixes of §3e (sdir 2, F(0)=1), gate at the DADI state (job 1936444) | **α 0.004 % · DV12 0.001 % · DV8 0.4 % (cd ≈ 1.6e-4)** — PG-off level | 83.8 · 0.08/0.34 | 88.4 · 0.08/0.14 | 52.1 · 0.75/0.29 | **gate passed**; trims unchanged vs it0 |

**Sensor 3** (2026-09-19): the reference model applies F(λ_θ) with the *model* thickness
θ_t = Re_θt·ν/U (P&Z Eqs. 10–14), bounded and uniform across the layer; sensor 2's physical θ ∝ d
grows with the laminar run and over-reads adverse gradients far downstream (e118 lower surface:
0.14 c vs 0.63 c in GR). Sensor 3 uses θ = Re_θc(Tu)·ν/U_e, λ = θ²U_e′/ν, U_e from the local
pressure — the same definition as the reference, no constant left. Option `SABCM_PG_sensor = 3`.

### 3e. Independent audit of the implementation (workflow, 2026-09-19: 3 lenses + refutation)

Confirmed correct: signs, Green-Gauss stencil scaling, Bernoulli without γ/M factors,
non-dimensionalisation (net exponent 0 in p_ref, ρ_ref, u_ref — also checked by a run with the
reference state scaled ×2/×0.5: BCM_lambda bit-identical), p from w = `computePressureSimple`,
K = 0.05064 (Blasius re-solve), blockette mirror identical, no in-place updates on the differentiated
path, `_b` has zero `pushreal8`. Fixed after the audit (commit 9b39601e): (1) `pInfCorr/rhoInf/uInf`
were input-only in the Tapenade head → `sasource_b` zeroed their seeds at entry (benign in the
present reverse ordering — the PG-off adjoint was bitwise identical to the baseline — but wrong);
now in-out, plus `velDirFreestream`; (2) F(0) was 0.9986 at p = 300 (smoothMax/smoothMin bias) → tanh
blend of the two LM branches, F(0) = 1 exactly; (3) ŝ = local velocity direction flips sign inside
separation-bubble backflow → `SABCM_PG_sdir = 2`: wall tangent oriented by the free stream (the
edge-streamline direction); (4) `cpModel = constant` guard for the pressure sensors.
Reg suite: the adjoint-stage failures (P DV, KeyError) are **pre-existing** — same on the baseline
build; those refs were written by a dev script and only hold `Eval Functions Sens:`.
Adjoint vs CS gap (sensor 2, α 1 %, shape 10–35 %): the two sides sit on **different steady states**
(ANK/NK primal vs the CS DADI path; the CS alpha derivative equals the adjoint *at the DADI state*
to 1e-4 % / 7e-3 %). The loop gate now forms the adjoint at the DADI-continued state (`adjB`) and
runs the CS from it.
Sensor 3 (θ = Re_θc ν/U_e): adjoint vs CS 46–105 % (λ ∝ 1/U_e² near stagnation, saturates), trips
both surfaces at 0.07–0.12 c → dropped.

### 3d. Why the local sensors trip earlier than GR — checked against the GR solution (2026-09-19)

Analytic λ_θ from the surface cp of the BCM C2 r1 shape (`u = U_e/U∞ = √(1−cp)`, λ = Re_θ²/Re · u′/u²):
lower surface x = 0.06–0.10 has a real adverse gradient, λ = −0.12…−0.23 (sensor-3 definition) →
F = 0.61 (clip); favourable at 0.15–0.25 (F 1.10); adverse again at 0.3 (F 0.69). The GR volume of the
same shape (`cross_eval/trim_gr_on_bcm_c2r1_e50/trim_000_vol.cgns`, field `TransitionRetheta`) shows
the transported R̅eθt **outside the layer following exactly this**: 950 (x 0.02–0.04), 560/542/550
(0.06–0.10 ≈ 880·0.61 = 537), 772/966 (0.15/0.20 ≈ 970), 667 (0.30). So the sensor reads the same
edge λ_θ as the reference model — the pressure sensor is *correct* as an implementation of the LM
correlation.

Near the wall (Re_v peak, d ≈ 5e-5…1e-3 c) GR's R̅eθt is lagged/diffused: 600–650 all along the lower
surface (0.68–0.74 of the ZPG value) instead of 540–970 — a mild history effect. The decisive
difference is the **onset criterion level**: GR fires at Re_S ≈ 2.6·1.35·Re_θc(R̅eθt) ≈ 2.46·R̅eθt
(≈ 1470 at R̅eθt 600), the BCM at Re_v ≈ 2.193·Re_θc^BCM·F = 970 at F 0.61, and the γ transport adds
a further delay in GR. With Blasius-like Re_v ≈ 2.19·Re_θ that puts the BCM-PG trip at x ≈ 0.12 c
and GR's at ≈ 0.4 c on this surface — exactly what the trims give (0.12 vs 0.42). The same F(λ)
therefore cannot bring the BCM to GR's x_tr on these shapes: the gap is the base onset criterion
(and its flat-plate calibration), not the gradient sensitivity. Closing it would mean changing the
BCM onset constants to GR's (= becoming the GR model) — not a calibration of the sensor.

What the local sensor does achieve: the ranking of the shapes becomes the GR one (GR C3 < GR C2 <
NACA < BCM optima), i.e. an optimiser driven by BCM-PG no longer exploits the nose suction-peak
loophole. Sensor 1 (Menter, gain 1) happened to land closest to GR's numbers (mean 13 %) because it
under-reads adverse gradients by ~2× (Falkner–Skan map) — a coincidence with GR's lag, and it is not
smooth (§5.1b). Sensor 2 (θ ∝ d, Blasius K) is the smooth, constant-free LM-local model.

## 4. Campaign (Deucalion, machV4 rebuilt from `sa-bcm-pg`)

| step | what | where | status |
|---|---|---|---|
| B1 | PG-off regression: reg suite (blockette/real/cs/adjoint vs existing refs) + NACA L0 cold trim = 68.20 counts | `tests/reg_tests/job_regtests_pg.sh`; `cross_eval_pg` | job 1934904 (reg suite), NACA in the cross jobs |
| B2 | fwd/rev/fast dot products + adjoint-vs-CS with PG on (`bcm_pg_tut_wing` variant) | same job | job 1934904 |
| B3 | adjoint vs complex step, NACA0012 L2, 20 shape DVs + α, PG on | `SA_BCM_ARTICLE/04_cs_naca2d/job_cs_naca_pg.sh` | job 1934905 |
| C2 | TMR flat plate L1/L2, bcm vs bcmpg (ZPG ⇒ identical) | `10_tmr_flatplate/10_bcm_pg/jobs/job_bcm_pg.sh` | job 1934906 |
| C3 | cross shapes trimmed with BCM-PG, gain {2, 1, 3}, vs GR targets | `08_optimization/2d_tu05_L0/cross_eval_pg/submit_cross_pg.sh` | jobs 1934907 (g2), 1934908 (g1), 1934909 (g3) |
| C4 | NLF0416 α 2/4/6 + S809 α 4/6/8, pg vs ref (same ranks) vs experiment/GR | `15_sabcm_polars/jobs/job_polars_pg.sh` | jobs 1934910 (pg), 1934911 (ref11) |
| D | C2/C3 optimisations relaunched with BCM-PG on L0 | `SA_BCM_ARTICLE/2d_opt/dense_L0v2_pg` | pending |

Targets for C3 (cd counts · x_tr up/lo): GR NACA 69.36 · 0.19/0.58 | GR C2 34.94 · 0.76/0.89 |
GR C3 34.02 · 0.77/0.90 | BCM C2 r1 e50 79.50 · 0.16/0.42 | BCM C3 e118 71.18 · 0.17/0.63.
BCM today: 68.20 · 0.22/0.49 | 51.38 · 0.74/0.34 | 39.22 · 0.70/0.78 | 53.55 · 0.33/0.66 | 46.39 · 0.48/0.66.

## 6. Conclusion (2026-09-19) — local sensors closed; line reopened as Plan B (§7)

A local pressure-gradient sensor on the BCM threshold reproduces the LM/GR edge λ_θ exactly, but the
BCM–GR disagreement on the optimised shapes comes from the reference model's onset level and transport
lag, not from missing gradient sensitivity. Applying F(λ) locally therefore over-penalises mild/short
adverse gradients: NLF0416 α4 x_tr 0.33 → 0.21 (GR/experiment 0.34), cd +21 %; NACA 68 → 73 (GR 69);
GR C2 unchanged (+49 %); only the BCM optima are pushed toward the GR values (−33 % → +5 %, −35 % → +24 %).
Closing the gap "seriously" would mean adopting GR's onset criterion and R̅eθt transport — i.e. using the
GR model, which already exists and is the robust optimiser (slide 20). Recommendation: optimise with
SA-γ-R̅eθt; keep SA-BCM for cheap primal evaluation; do not use `SABCM_PG` for predictions.
Side products kept: exact `nWall` (wall-normal) block array with adjoint plumbing, `smoothMinMax` /
`bcmFlambda` helpers, the fast-reverse in-place trap and the reference-scalar head trap (both documented
in `docs/adjoint-trace.md`), and the one-state rule for adjoint-vs-CS checks on this multi-state primal.

## 7. Reopened (2026-09-19 evening): gradient response F=2 (Menter) and Plan B (transported threshold)

User: "descobre como resolver isto, pesquisa e testa até dar". Two more steps.

**F=2 (`SABCM_PG_F=2`, 3dac53f3): Menter-2015 F_PG ratio for adverse + LM F1 for favourable** (§2, Falkner–Skan
inverse map, `falkner_skan_fpg.py`). Direct-prediction tests (all sensor 2, g 2.0, K 0.05064):

| test | BCM | F=1 (LM) | F=2 (Menter) | GR |
|---|---|---|---|---|
| NACA L0 α 0…6 (`pg_calib`, mean \|Δcd/cd_GR\|) | 2.4 % | 6.5 % | **2.1 %** | — |
| NLF0416 L1 α 2/4/6 cd (x_tr,up) | 58.2 (0.38) / 65.5 (0.33) / 85.6 (0.21) | 65.1 (0.31) / 79.4 (0.21) / 102.3* (0.11) | **58.9 (0.38) / 65.9 (0.33) / 88.9 (0.18)** | 57.0 / 64.4 / 79.3 (0.34 exp @ α4) |
| S809 L1 α 4/6/8 cd | 47.6* / 75.0 / 149.1 | 48.7* / 74.6 / 149.9 | 80.1* / 72.5* / 148.3 | 67.5 / 69.8 / 131.8 |
| cross GR C2 / GR C3 / BCM C2 r1 (GR 34.9 / 34.0 / 79.5) | 51.4 / 39.2 / 53.5 | 52.1 / 32.5 / 83.8 | 51.1 / 31.8 / **54.6** | — |
| adjoint vs CS at one state (pg_loop `it3_s2_F2`, α + DV 8/12) | 0.008 % / 0.26 % | 0.00 % / 0.43 % | 0.00 % / 0.43 % | — |

(* not converged, ANK limit cycle — same S809 α4 stall as plain BCM.) F=2 is the least-bad local form: it
no longer degrades NLF/NACA (it is plain BCM there, Menter's F_PG ≈ 1 for mild λ), the derivatives hold —
but it does **not** close the loophole (BCM C2 r1 54.6 vs GR 79.5): Menter's ratio under-penalises the strong
nose gradient that LM over-penalises. Checked against the GR volume (§3d): GR trips at x = 0.16 on that
upper surface where the local λ_θ is mild — R̅eθt there is dragged down from the suction peak by transport.
**No local response fits both** the mild (NACA/NLF) and the strong (BCM optimum) gradients: the missing
ingredient is the lag/history, not the shape of F(λ).

**Plan B (handoff §"plano B"): BCM's γ with GR's transported threshold** — implemented 2026-09-19 as option
`transitionBCMGamma` of the SA-γ-R̅eθt model on branch `transition-models` (`adflow_opt_improvements`,
`e84475f9`; docs `CORE_BASE/CORE_01_architecture.md`): γ_BC algebraic (tanh/KS form, fsmooth 0.08, same
constants) with `Re_θc = Re_θc^BCM(Tu) · R̅eθt / Re_θt(Tu, λ=0)`. ZPG ⇒ exactly BCM; with gradient the
threshold follows the transported R̅eθt (lag + history). The γ equation stays solved but only as
diagnostic. Tests (machV3): `08_optimization/2d_tu05_L0/pg_calib/job_bcmg_naca_sweep.sh` (job 1937059),
`cross_eval_bcmg/submit_cross_bcmg.sh` (1937060), derivative gate `bcmg_gate/job_bcmg_gate.sh` (after
the complex build 1937061). Results: `cross_eval_bcmg/PURPOSE.md` and the 18/09 deck.

## 5. Results (2026-09-18, machV4 @ 50b51a76)

### 5.1 Adjoint / AD
* PG off: every reference-based AD test (fwd, rev, fast-rev, dot products, smooth+hard) passes against
  the 2026-08-25 refs → the regenerated AD is unchanged for the existing model. NACA L0 cold trim is
  in the cross jobs.
* PG on, first build (dcb4c84e): forward AD = complex step (5/5 `cmplx_test_*`), fwd↔rev dot products
  OK, **but `_b` vs `_fast_b` failed** and the NACA L2 adjoint-vs-CS gave up to 180 % errors. Cause:
  `ReThetaCrit = ReThetaCrit*Flam` (in place) — Tapenade guards it with `pushreal8/popreal8` in
  `sa_b.f90`, the fast-reverse autoEdit has no stack, so `flamd` used the overwritten value. Fixed
  with a distinct target (`ReThetaCrit = ReThetaBase*Flam`, commit 50b51a76; no `pushreal8` left in
  `sa_b.f90`); local check: all pg `_b`/`_fast_b`/dot-product tests pass. HPC re-verification: jobs
  1935762 (reg suite) and 1935763 (NACA L2 adjoint vs CS) — see §5.5.
* The AD-vs-FD `wDot` check (`TestJacVecFwdBCMFD.test_wDot`, h=1e-8, rtol 8e-4, atol 0) fails for
  smooth, hard **and the pre-change baseline ae1a7dbb** with identical numbers (18.3 % mismatched,
  max rel 195) — pre-existing fragility of that test, not a PG regression.

### 5.1b Adjoint vs CS after the fix — still wrong, and why (sensor 1)
Rerun 1935763 (fixed fast-reverse): max |Δ| 187 % (cl) / 181 % (cd), median cd 74 % — unchanged. Real
central differences on the same case (job 1935808, cold solves converged to 4e-8, noise pair = 0):
dcd/dα = −0.79 / +0.28 / +0.48 for h = 3e-3 / 1e-3 / 3e-4; shape DVs change sign and order of
magnitude with the step. The primal with the velocity-profile sensor is not a smooth function of the
DVs (transition/profile feedback on top of the BCM branch structure): adjoint and CS measure tangents
of different branches and neither is "the" derivative. This is what motivated sensor 2 and the
smoothing/verification loop (`08_optimization/2d_tu05_L0/pg_loop`).

### 5.2 Cross-model (NACA0012 L0, Tu 0.5 %, cl 0.3; counts · x_tr up/lo; `cross_eval_pg`, jobs 1934907-09)

| shape | GR (target) | BCM | BCM-PG g=1 | BCM-PG g=2 | BCM-PG g=3 |
|---|---|---|---|---|---|
| NACA 0012 | 69.36 · 0.19/0.58 | 68.20 · 0.22/0.49 | 71.24 · 0.18/0.46 | 73.28 · 0.16/0.45 | 74.59 · 0.14/0.44 |
| GR C2 | 34.94 · 0.76/0.89 | 51.38 (+47 %) · 0.74/0.34 | 52.09 · 0.75/0.30 | 52.49 · 0.74/0.29 | 52.80 · 0.74/0.28 |
| GR C3 | 34.02 · 0.77/0.90 | 39.22 (+15 %) · 0.70/0.78 | 33.63 (−1 %) · 0.76/0.83 | 32.50 (−4 %) · 0.76/0.87 | 32.43 (−5 %) |
| BCM C2 r1 e50 | 79.50 · 0.16/0.42 | 53.55 (−33 %) · 0.33/0.66 | 70.63 (−11 %) · 0.12/0.62 | 83.59 (+5 %) · 0.08/0.34 | 85.17 (+7 %) |
| BCM C3 e118 | 71.18 · 0.17/0.63 | 46.39 (−35 %) · 0.48/0.66 | 69.10 (−3 %) · 0.12/0.65 [g1.5: 70.54 (−1 %)] | 74.07 (+4 %) · 0.08/0.56 | [g2.5: 88.91 (+25 %)] |

mean |Δcd/cd_GR| (4 shapes): BCM 24.2 % → PG g1 16.0 %, g2 16.4 %, g3 17.6 %; over 5 shapes (with e118, job 1935767): BCM 26.3 % → PG g1 13.4 %. Ranking vs the own NACA
(what an optimiser sees): GR says BCM C2 r1 is +14.6 % *worse* than NACA and GR C3 −51 %; BCM said
−21 % / −42 %; **BCM-PG g=2 says +14 % / −56 %** — the sensor closes the gap that motivated the work.
GR C2 stays +50 %: all BCM variants trip the lower surface at 0.3 c where GR keeps laminar flow to
0.9 c (`cross_pg_cf.png`); F(λ) cannot raise the threshold enough there (F ≤ 1.10 at Tu 0.5 %) — a
base-threshold/criterion difference (BCM Re_θc(0.5 %) = 725 vs Langtry Re_θt = 880), not a gradient
effect. Some PG trims end at ~1e-6 relative instead of 1e-8 (`fail` flag; forces converged). Sensor 1 is
superseded by sensor 2 (§3b) after the derivative findings (§5.1b).

### 5.3 Polar (NLF0416 / S809 L1, 3 α, `15_sabcm_polars`, jobs 1934910/11, 1935125)
With PG (g=2) the low-Tu airfoils (Tu 0.15/0.07 %) **limit-cycle in ANK** (resrho 0.2–0.7 for 5000+
its, never reach NK) at NLF α6 and S809 α4/α6, while the PG-off twins (same 11 ranks) converge; none
of {SABCM_PG_p 100, gain 1, ANK CFL 1e4, NK @1e-3} fixed it (job 1935125). Where it converged, PG
trips NLF α4 earlier than GR/experiment (x_tr,up 0.19 vs 0.33/0.34; cd 81 vs 65 counts): the sensor
reads λ ≈ −0.005 (raw, → F ≈ 0.89) at the Re_v peak (`analyze_pg_nearwall.py`, x = 0.10 profile),
i.e. a mild adverse gradient that the GR's transported edge-λ_θ does not see. Continuation test
(PG switched on from the converged PG-off state): job 1935764. `nlf a2` twice killed by the
SA-BCM startup NaN (cl/cd monitors) — monitors removed from `run_case_bcm.py`, rerun in 1935764.

### 5.3b Polar with sensor 2 (job 1936549, 11 ranks, rNK; cd counts · x_tr up/lo)

| case | BCM (no PG) | BCM-PG sensor 2 | GR (campaign 12) | exp x_tr,up |
|---|---|---|---|---|
| NLF a2 | 58.2 · 0.38/0.60 | 65.1 · 0.31/0.60 | 57.0 · 0.38/0.60 | — |
| NLF a4 | 65.5 · 0.33/0.61 | 79.4 · 0.21/0.60 | 64.4 · 0.34/0.61 | ≈ 0.34 |
| NLF a6 | 85.6 · 0.21/0.61 | 102.3 · 0.11/0.61 (not conv.) | 79.3 · 0.25/0.61 | — |
| S809 a4 | 47.6 (not conv.) | 48.7 (not conv.) | 67.5 | — |
| S809 a6 | 75.0 · 0.49/0.49 | 74.6 · 0.50/0.49 | 69.8 · 0.50/0.50 | — |
| S809 a8 | 149.1 · 0.96/0.50 | 149.9 · 0.96/0.49 | 131.8 | — |

S809 (transition by laminar separation): the sensor is neutral, as it should be. NLF0416
(Tu 0.15 %, attached-flow transition in a mild adverse gradient): the local F(λ) trips the upper
surface earlier than GR *and* experiment (α4: 0.21 vs 0.34) → cd +21 %. Same mechanism as §3d:
the local LM correlation penalises the mild adverse gradient (F ≈ 0.8), which the reference model's
lagged near-wall R̅eθt and higher onset level do not. Plain BCM already matched experiment here.
Verdict on the user's direct-prediction criterion: **the local sensor is not a net improvement** —
it fixes the ranking of the optimised shapes at Tu 0.5 % (§5.2) but degrades the NLF polar. The
BCM-PG C2/C3 optimisations (job 1936548) will show whether its optima are at least accepted by GR.

### 5.4 Flat plate (TMR, `10_tmr_flatplate/10_bcm_pg`, jobs 1934951, 1935139)
Not usable as the neutrality test: the SA-BCM itself is **bistable** on this plate at Tu 0.1 %
(rNK/CANK recipes converge to the fully turbulent branch, cd_wall 0.00317 = SA; the ANK-only /
CFL-capped recipes stay laminar and never converge). With γ = 1 from Term2 the threshold is
irrelevant, so bcm ≡ bcmpg there trivially. Neutrality on a Blasius layer rests on the
Falkner–Skan analysis (§3) and on the near-wall profiles (λ = offset at the wall, Snn → 0).

### 5.5 Open (updated as jobs land)
- HPC re-verification of the adjoint after the fix (1935762, 1935763).
- PG convergence at low Tu (1935764).
- C2/C3 with BCM-PG sensor 2 (`dense_L0v2_pg`): job 1936548 (32 + 32 ranks, 40 h). Polar NLF/S809 3 α with sensor 2: job 1936549.
