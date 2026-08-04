# Paper-validation campaign — S809 / NLF0416 vs Piotrowski & Zingg (2020)

**Date:** 2026-08-03/04.  **Goal:** make this implementation reproduce the
paper's published S809/NLF0416 results.  **Status:** cause isolated to model
*behavior* in the adverse-pressure-gradient / high-cl regime; the drag gap
itself is NOT yet closed. Campaign paused by user with next step defined.

All run artifacts live in the run tree
(`/home/mdo/Desktop/Run/MDO_PhD/Transition/gama_rethetha/06_alpha_sweep/` —
each study folder has a PURPOSE.md); this file is the consolidated record.

## 1. The problem

S809 (M 0.1, Re 2.0e6, Tu 0.07%): drag-bucket corner ~0.5–0.6° early
(cl ≈ 0.70 vs paper 0.79) and cd ~20–30 counts high above the corner
(α=5: ours 87.6 L1 / 98.5 SB3 vs paper ≈ 66; α=6: 113–116 vs ≈ 85).
Bucket floor agrees with experiment; the paper sits 4–6 counts *below*
experiment, so vs the paper we are offset everywhere, worst above the corner.

## 2. Decisive experiment — the paper's own grid

The paper's NLF0416 TCMPS C-grid family is public (1st AIAA Transition
Modeling & Prediction Workshop, Case 2; dims match paper Table 1 verbatim;
grid POC Jim Coder, jcoder@utk.edu — the S809 family is private, request
email drafted). Converted to ADflow CGNS
(`04_mesh_families/nlf0416_tcmps_workshop/convert_tcmps_to_cgns.py`: C-grid
split at wake-cut ends, B2B via cgnsutilities connect; span in y, lift in z
→ liftIndex 3; conditions via mach+reynolds+T since Re 4e6 @ M 0.1 has no
physical altitude).

Result (α-sweep at the paper's exact conditions, vs paper model curve at
matched cl): **+1.1 / +1.2 / +6.7 / +22.3 / +30.2 counts at cl −0.23 /
0.49 / 0.72 / 0.93 / 1.14.**

> With grid, geometry, M, Re, Tu all identical, the S809 gap signature
> reproduces → **not grids, not geometry: implementation behavior under
> adverse pressure gradient.**

## 3. Hypotheses tested and eliminated

| Hypothesis | Test | Verdict |
|---|---|---|
| Mesh resolution / family | mesh family to 381k + exact-paper-grid run | eliminated (family ~few counts; exact grid reproduces gap) |
| Geometry rendering/source | CST refits (earlier work) | real but only ~10 counts; fixed |
| Reynolds offset | Re 1.5/2.0/2.5e6 bracketing | eliminated |
| Freestream Tu level | existing Tu sweep 0.02–0.07% | ~1 count / 0.01% — far too weak |
| JST dissipation level | vis4 ×½ / ×2 on SB3 | ×2 converged: −2.3 counts (wrong size+sign); ×½ destabilizes into limit cycle |
| Matrix dissipation (paper's scheme) | cold, cold+conservative ANK, warm-started | all 3 fail to converge — parked, needs dedicated stabilization |
| λθ ±0.1 clamp not in paper | `rsaGRclampLambdaTheta=.false.` build, S809 α 4–6 | REJECTED: stalls ~1e-4, cd +38 counts wrong direction; clamp is part of correlation validity; reverted (`efed31cf`) |
| C-grid wake cut dislikes α≠0 | TCMPS grid rotated +4°, run at α=0 | ~2 counts vs freestream-α run — eliminated |
| Farfield circulation correction | code inspection | feature doesn't exist in ADflow (dead `vortexCorr` flag); sub-count at 470–1000c anyway |
| Mach 0.10 vs 0.17 (pyHyp sweeps) | pyHyp L-1 at M 0.10 | ~3 counts; the ~14-count α=4 grid-family difference is real but doesn't affect the paper comparison |

## 4. Symbol-level audit vs the paper — all clean

Re-verified line-by-line against the paper (2026-08-03/04): Re̅θt/F(λθ)
correlations (Eqs. 9–14, 54–57 incl. clamp semantics), source terms P_γ/E_γ
(Eqs. 52–53) with vorticity limiter incl. the √γ p-ρ nondim factor,
F_onset (46–47), F_turb (48), F_length (49–50), Re_θc (51), dU/ds
(contraction *and* the 2·fact/2·vol stencil scaling), F_θt/F_wake/δ/θ_BL
(Eqs. 3–5, explicit-Re absorption verified), σ_θt = 2.0 / σ_f = 1.0
placement, timescale 500ν/U². **No transcription error exists.**

## 5. Hard targets extracted from the paper PDF

The markdown KB lacks figure data; the PDF with figures is
`doc/PiotrowskiandZingg2020 (2) (1).pdf`. Digitized at 300 dpi into
`06_alpha_sweep/common/` (run tree):

- `reference_s809_cd_model_fig8e.csv` — cd floor 63–66 counts through α=5,
  corner between 5 and 6, 85 @ α=6.
- `reference_s809_xtr_model.csv` — upper x_tr ≈ 0.53 / **0.49** / 0.29 at
  α 4/5/6 (ours collapses at α≈4–4.5: the defect in physical terms).
- `reference_nlf0416_xtr_model.csv` — upper x_tr(cl): 0.30/0.22 at
  cl 0.95/1.17 (where we are +22/+30 counts).

## 6. Side findings

1. **The bucket corner is non-unique/unsteady in steady RANS:** S809 SB3
   α=5 cold-start converges (98.5 counts, 1e-12); warm-started from the
   converged α=4 solution it *never* converges (300k cycles, residual
   wandering). Caps how precisely any steady code can pin the corner.
2. **Large cycle budgets destabilize the late NK/deep-ANK phase:** NLF
   exact-grid cases at 60k budget reach 6–7 orders cleanly; at 300k, 5/7
   destabilize (one to negative cd). Related to the guarded-NK work in
   `NKSolvers.F90` (checkpoint `03df399e`).

## 7. Next step (defined, not yet run)

One ~2-min HPC job: restart NLF α=4 (and S809 α=4/5) from existing
solutions with nCycles≈1 + `addSlices` cf/cp → our x_tr(α) overlaid on the
§5 curves. Smooth departure ⇒ correlation-level difference; abrupt ⇒
bubble/onset dynamics. NOTE: the runner's `addSlices` hardcodes axis "z";
the TCMPS grid spans y — one-line tweak needed first.

## 8. Code state

Branch `new_conv_strategie_test` @ `efed31cf`: λθ clamp switchable
(default ON = legacy), Re̅θt_eq ≥ 20 floor in `reThetaTCorrelation`
(inert with clamp on). Local mach and Deucalion machV2 builds in sync.
TAPENADE NEEDED before any adjoint use of these changes.
