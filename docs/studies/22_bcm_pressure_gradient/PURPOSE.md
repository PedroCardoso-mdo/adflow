# 22 — SA-BCM with a local pressure-gradient sensor (SABCM_PG)

**Started:** 2026-09-18  **Status:** model implemented + Tapenade regenerated (branch `sa-bcm-pg`);
validation/calibration campaign on Deucalion in progress (see §4).

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

## 4. Campaign (Deucalion, machV4 rebuilt from `sa-bcm-pg`)

| step | what | where | status |
|---|---|---|---|
| B1 | PG-off regression: reg suite (blockette/real/cs/adjoint vs existing refs) + NACA L0 cold trim = 68.20 counts | `tests/reg_tests/job_regtests_pg.sh`; `cross_eval_pg` | pending |
| B2 | fwd/rev/fast dot products + adjoint-vs-CS with PG on (`bcm_pg_tut_wing` variant) | same job | pending |
| B3 | adjoint vs complex step, NACA0012 L2, 20 shape DVs + α, PG on | `SA_BCM_ARTICLE/04_cs_naca2d/job_cs_naca_pg.sh` | pending |
| C2 | TMR flat plate L1/L2, bcm vs bcmpg (ZPG ⇒ identical) | `10_tmr_flatplate/10_bcm_pg/jobs/job_bcm_pg.sh` | pending |
| C3 | cross shapes trimmed with BCM-PG, gain {2, 1, 3}, vs GR targets | `08_optimization/2d_tu05_L0/cross_eval_pg/submit_cross_pg.sh` | pending |
| C4 | NLF0416 α 2/4/6 + S809 α 4/6/8, pg vs ref (same ranks) vs experiment/GR | `15_sabcm_polars/jobs/job_polars_pg.sh` | pending |
| D | C2/C3 optimisations relaunched with BCM-PG on L0 | `SA_BCM_ARTICLE/2d_opt/dense_L0v2_pg` | pending |

Targets for C3 (cd counts · x_tr up/lo): GR NACA 69.36 · 0.19/0.58 | GR C2 34.94 · 0.76/0.89 |
GR C3 34.02 · 0.77/0.90 | BCM C2 r1 e50 79.50 · 0.16/0.42 | BCM C3 e118 71.18 · 0.17/0.63.
BCM today: 68.20 · 0.22/0.49 | 51.38 · 0.74/0.34 | 39.22 · 0.70/0.78 | 53.55 · 0.33/0.66 | 46.39 · 0.48/0.66.

## 5. Results

(to be filled)
