# A1 — Confirmação SA-GR vs paper (com não-adimensionalização p-ρ)

Data: 2026-07-06. Tarefa A1 de `state.md`: confirmar (não re-auditar) o código
vs Piotrowski & Zingg (2020), com a lente da não-adimensionalização
(`docs/nondimensionalization.md` lido primeiro). **Este ficheiro contém apenas
divergências**; tudo o que não está listado aqui foi conferido e bate certo
com o texto integral do paper
(`docs/SA_GAMMA_RETHETHA_BASE/Piotrowski_Zingg_2020_SA-sLM2015_clean (1).md`).

Âmbito conferido (sem divergência, não repetido abaixo): timeScale (Eq. 7,
Re=1 implícito), θ_BL/δ_BL/δ (Eq. 4), Re_S (Eq. 5), λ_θ e dU/ds (Eqs. 10–14,
contração ûᵢûⱼ∂uᵢ/∂xⱼ idêntica), correlação Re_θt (Eqs. 8–9), F(λ_θ) suave
(Eqs. 54–57 **verbatim**), F_length (Eqs. 49–50), Re_θc (Eq. 51), F_onset
(Eqs. 46–47, constantes 2.6/6/1.35), F_turb (Eq. 48), P_γ (Eq. 52, **com √γ**),
E_γ (Eq. 53, com Ω limitado), F_θt/F_wake (Eqs. 3 e 5, formas simplificadas
SA-LM2015), acoplamento γ só na produção SA (Eq. 41), difusões σ_f=1 divisor /
σ_θt=2 multiplicativo (Eqs. 1, 27), crossflow D_scf (Eqs. 15–26, comprimentos
físicos), φ_p Algorithm 1, bounds do damping (Algorithm 2: γ∈[1e-10,2],
Re̅θt≥20), constantes §9 todas, Tu em percento (`turbIntensityInf*100`,
default 0.001 = fração).

---

## D1 — Limitador de vorticidade usa L_ref = 1 m em vez da corda de raiz (única divergência código↔paper)

- **Código:** `src/turbulence/saGammaRetheta.F90:502` —
  `vortLim = uInf*sqrt(uInf/muInf)/20`.
- **Paper:** Eqs. 52–53 usam `M·√(M·Re)/20` com `Re = ρ∞·a∞·l/μ∞` e, no texto
  (linha ~170 do MD), «the reference length is specified as the root chord».
- **Análise:** convertendo para unidades ADflow (uInf = M√γ,
  ρRef·uRef/μRef = 1, L_ref = 1 m), a expressão do código é *exatamente* a do
  paper com `l = 1 m`. Para uma malha com corda ≠ 1 m, o teto de vorticidade
  difere do paper por um fator `√(corda[m])`.
- **Severidade:** baixa. O próprio paper (§IV, linha ~519) diz que o limitador
  é puramente numérico («does not aim to improve its predictive capability»).
  Mas é a única fórmula portada onde o `Re` do paper NÃO é conversão de
  unidades e sim escala de calibração com `l` explícito — a regra do
  `nondimensionalization.md` §5 («drop Re») não é exata aqui.
- **Prova:** `grep -n "vortLim" src/turbulence/saGammaRetheta.F90` e
  `grep -n "root chord" "docs/SA_GAMMA_RETHETHA_BASE/Piotrowski_Zingg_2020_SA-sLM2015_clean (1).md"`.

## D2 — Salvaguardas LM2009 presentes no código mas não escritas no paper (informativo)

- **λ_θ clamp suave a [−0.1, +0.1]** — `saGammaRetheta.F90:551-552`.
- **Floor de Tu a 0.027%** — `turbUtils.F90:2297` (`Tu_safe`).
- O texto de P&Z (Eqs. 10–14 e 8–9) não enuncia estes limites; eles vêm da
  definição original das correlações de Langtry–Menter (2009) e são prática
  universal. Não considerado erro; registado porque não é derivável do paper.

## D3 — Erros no doc destilado `docs/SA_GAMMA_RETHETHA_BASE/paper-reference.md` (doc↔paper, não código↔paper)

O código está certo nestes quatro pontos; o destilado — declarado «source of
truth» — está errado e pode induzir tarefas futuras em erro. **Não corrigi o
doc** (fora do âmbito A1; decisão do utilizador).

| # | paper-reference.md diz | Paper integral diz | Código |
|---|---|---|---|
| 1 | §2.3 (Eq. 3): F_θt = min(max(F_wake·e^…, 1−((γ_eff−1/ce2)/(1−1/ce2))²), 1) | Eq. 3: `F_θt = F_wake·e^{−(d/δ)⁴}` (forma simplificada SA-LM2015; ramo γ_eff removido de propósito, §II) | igual ao paper (`saGammaRetheta.F90:564`) |
| 2 | §2.5: F_wake = e^{−(Re_ω/1e5)²} com Re_ω = ρΩy²/μ | Eq. 5: `F_wake = e^{−Re_S/1e6}` com Re_S baseado em strain (substituição de Schücker, §II) | igual ao paper (`saGammaRetheta.F90:563`) |
| 3 | §4.1 (Eq. 52 sLM2015): P_γ = ca1·F_length·Ω_lim·F_onset·(1−ce1·γ) — **sem √γ** | Eq. 52 tem `√γ` | igual ao paper (`saGammaRetheta.F90:522-524`) |
| 4 | §3.3 (Eqs. 54–57): F3 = φ₋₃₀₀(F1,0)+φ₊₃₀₀(F2,0); F = 1+F3, com F1 de sinal tal que o ramo APG desaparece | Eqs. 54–57: F1 = 1+0.275(1−e^{−35λ})e^{−Tu/0.5}; F2 = φ₊₃₀₀(F1,1); F3 = 1−(−12.986λ−123.66λ²−405.689λ³)e^{−(Tu/1.5)^1.5}; F = φ₋₃₀₀(F2,F3) | igual ao paper (`turbUtils.F90:2299-2315`) |

- **Prova (ex.):** `grep -n -F "F_{wake}" "docs/SA_GAMMA_RETHETHA_BASE/Piotrowski_Zingg_2020_SA-sLM2015_clean (1).md"` (linha 155) vs `paper-reference.md` §2.5; `sed -n '540,557p'` do paper integral vs §3.3.
