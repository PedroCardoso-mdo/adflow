# Paper-validation campaign — S809 / NLF0416 vs Piotrowski & Zingg (2020)

**Date:** 2026-08-03/04.  **Goal:** make this implementation reproduce the
paper's published S809/NLF0416 results.  **Status:** cause isolated to model
*behavior* in the adverse-pressure-gradient / high-cl regime; the drag gap
itself is NOT yet closed. Campaign paused by user with next step defined.

All run artifacts live in the run tree
(`/home/mdo/Desktop/Run/MDO_PhD/Transition/gama_rethetha/06_alpha_sweep/` —
each study folder has a PURPOSE.md); this file is the consolidated record.

## 1. The problem

> **⚠️ Superseded by §A1 (2026-08-05):** the SB3-mesh cd numbers in §§1/6
> are from runs that stagnate at 1e-4..1e-3 (`converged: false`) — they
> are not model measurements.

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
| Matrix dissipation (paper's scheme) | cold, cold+conservative ANK, warm-started on SB3; **defaults on L1, 13 alphas (job 1807627, 2026-08-05)** | **RESOLVED 2026-08-06**: it was being run below its stability threshold. Sweeping vis4 inside the matrix scheme (86 runs, jobs 1808897/1808924/1808951/1809388) converges above vis4 ~0.03, working range 0.08–0.1 (blocks of alphas at 12–13 orders). There, cd is insensitive to the coefficient (63/115/173 counts at alpha 0/6/10 across a 3x change) and it is the CLOSEST result to the paper above the corner in the whole campaign (+23.6 counts vs +25.1 for the best scalar). The corner still sits at alpha ~4.5 |
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
   *(⚠️ Superseded by §A1, 2026-08-05: the 98.5-count "converged baseline"
   attribution is wrong — that 1e-12 run was the doubled-vis4 probe; the
   SB3 baseline never converged.)*
2. **Large cycle budgets destabilize the late NK/deep-ANK phase:** NLF
   exact-grid cases at 60k budget reach 6–7 orders cleanly; at 300k, 5/7
   destabilize (one to negative cd). Related to the guarded-NK work in
   `NKSolvers.F90` (checkpoint `03df399e`). *(Re-diagnosed 2026-08-07/08
   as premature NK engagement at nkswitchtol 1e-6 — see
   `../convergence-strategy.md`, commits ef9fc10d/54c43475.)*

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

*(Stale as of 2026-08-12: `new_conv_strategie_test` has since been merged
into `sa_gamma_rethetha` (3b78dd1b), and Tapenade was regenerated
2026-08-04 with AD files committed (2c4ce2c1) — the "TAPENADE NEEDED"
blocker is resolved.)*

---

# ADENDA 2026-08-05 — a SB3 nunca esteve convergida

Campanha de convergência a pedido do utilizador (bucket S809 com dissipação
reduzida). Pastas: `06_alpha_sweep/s809/dissipation_convergence/1_sb3_vis2_0_vis4_0p01/` e
`06_alpha_sweep/s809/dissipation_convergence/` (PURPOSE.md em ambas);
dashboard interativo com as 319 corridas S809 em
`06_alpha_sweep/s809/s809_dashboard.html`.

## A1. O que obriga a rever este relatório

**A malha SB3 — a malha em que assentam as §1–§6 acima — estagna a
1e-4…1e-3 mesmo com a dissipação default.** Todos os α do sweep
`results_geomfix_m010_SB3` têm `converged: false`. A única corrida SB3 no
histórico que chegou a 1e-12 é o probe de **vis4 dobrado**
(`vis4dbl_a5p0`, 96.2 counts) — e é esse o número que a §6 atribui ao
baseline ("SB3 α=5 cold-start converges, 98.5 counts, 1e-12"). A
atribuição está errada: o baseline não convergiu.

Consequência: os desvios ao paper citados na §1 (~20–30 counts acima do
canto) foram medidos em estados estagnados, e a dispersão entre receitas
no mesmo α chega a 40 counts. **Não são medidas do modelo.**

## A2. O que foi tentado para convergir a SB3 (16 receitas, ~136 corridas)

| alavanca | resultado |
|---|---|
| `ANKSecondOrdSwitchTol` 1e-5 → 1e-2 | α=0 passa de 2.6 → **12.6 ordens**; α≥4 fica em ~1e-3 |
| NK forçado a 3e-3, subspace 400 | piso 1e-3 → 1e-5 (a única alavanca que o moveu) |
| `NKLS = non-monotone` | melhor variante do NK (5.8e-6) |
| `turbResScale` do Re̅θt (1e-6, 1e-3, 1e-2) | nada nos dois sentidos; 1e-6 + PC forte → NaN |
| turbulência por KSP, srcDt sempre ligado, CFL limitado | nada (KSP custa 6–8× o tempo) |
| continuação de dissipação (arranque quente e nativa) | o arranque quente **piora** (2e-2) |
| 1e6 ciclos | orçamento não é o limite |

Assinatura em todos: `Step = 0.00` no NK com o **resíduo do Re̅θt** preso
(~8–15, ~46 % do total) enquanto o resíduo de massa já vai em 1e-3…1e-4.
*(Re-diagnosed 2026-08-07/08 as premature NK engagement at nkswitchtol
1e-6 — see `../convergence-strategy.md`, commits ef9fc10d/54c43475.)*

## A3. Onde converge

`05_mesh_independence/s809_v2` (mesh L1): **11.5–12.4 ordens** com vis4 de 0.0156 até 0.005, a
α = 0 e 4 (matriz 2×5) e ao longo de todo o bucket até α = 7. O tail
α ≥ 8 fica em ~5 ordens, como já acontecia no sweep L1 original.

## A4. Dissipação — o que se sabe agora

- **`vis2` é inerte** a M = 0.10: quatro valores (0.25 → 0) dão cd e
  ordens idênticos, no S809 e na flat plate (matriz completa 4×5,
  `10_tmr_flatplate/PURPOSE.md`). A comparação vis2 = 0 vs 0.25 desta
  campanha não estava a medir nada.
- **`vis4` tem limiar ≈ 0.01 nas malhas finas.** Flat plate: L0 converge
  até 0.0078 e falha a 0.005; L1 e L2 falham já a 0.0078. Na S809 L1
  converge-se até 0.005 — o limiar é da malha, não do coeficiente.
- **A sensibilidade do cd ao vis4 era erro de malha**: flat plate L0
  26.5 → 19.3 counts ao longo da gama, L1 20.4 → 20.6, L2 21.1 → 21.0.

## A5. Bucket convergido (L1, job 1805694) vs o paper

| | nosso (vis4 0.0156) | paper |
|---|---|---|
| fundo do bucket | 63.8 counts | 63–66 |
| canto | α ≈ 4.5 | entre 5 e 6 |
| cd a α = 6 | 120.7 | 85 |

Baixar o vis4 corta até 30 counts a α = 10 (207 → 175) e **nada** abaixo
de α = 5. Ou seja: a dissipação explica parte do excesso a alto α e
**nenhuma parte do canto prematuro** — que continua a ser o defeito real,
agora medido em solução convergida.

## A6. O que muda no plano

O passo definido na §7 (x_tr(α) sobreposto às curvas digitalizadas)
mantém-se válido e passa a ser mais informativo, **mas tem de correr na
L1**, não na SB3. Qualquer comparação futura com o paper feita na SB3
mede o piso de convergência dessa malha.

## A7. Dois testes que fecham a hipótese numérica (2026-08-05)

- **Matrix dissipation na L1, tudo no default** (job 1807627): falha nos 13 α
  (1.3–2.9 ordens, cd 1908…20013 counts). Corrige a leitura anterior — o
  problema não era o transiente da SB3. **E é ele próprio corrigido a
  2026-08-06:** varrendo o vis4 dentro do esquema matricial, converge acima de
  ~0.03, com zona útil 0.08–0.1 e o melhor acordo com o paper acima do canto de
  toda a campanha (+23.6 counts). O default 0.0156 estava abaixo do limiar de
  estabilidade do esquema.
- **vis4 acima do default** (job 1807690, 0.05 e 0.08 na L1): converge
  (7.4–11.6 ordens) mas afasta-se do paper em todo o α (+19/+61 e +26/+82
  counts) e **não desloca o canto**, que fica em α ≈ 4.5 com vis4 de 0.005 a
  0.08. Com a direção oposta já testada, o canto prematuro deixa de ter
  explicação dissipativa.
