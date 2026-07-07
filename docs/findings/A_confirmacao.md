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

## D1 — Limitador de vorticidade usa L_ref = 1 m em vez da corda de raiz (única divergência código↔paper) _______________________________FEITO________________________________________

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
- **Resolvido** comprimento trocada por default para a corda defenida no aero problem, utilizador pode na mesma defenir um valor manuealmente explicação em D1.md

## D2 — Salvaguardas LM2009 presentes no código mas não escritas no paper (salvaguarda necessária) _______________________________FEITO________________________________________

Validado a fundo em 2026-07-07 (tentativa de refutação clamps × suavizações;
p = ±300 de `paramTurb.F90:44-45`, escala de suavização log(2)/300 ≈ 0.0023).

- **λ_θ clamp suave a [−0.1, +0.1]** — `saGammaRetheta.F90:560-561`
  (smooth max com −0.1, depois smooth min com +0.1). As zonas de transição
  (~1/300) não se sobrepõem (bounds distam 0.2) e ambos os vieses são para
  dentro do intervalo: o resultado está provado ∈ [−0.1, +0.1] sempre.
- **Floor de Tu a 0.027%** — `turbUtils.F90:2297` (`Tu_safe`, smooth max).
  Erra para o lado seguro (`Tu_safe ≥ 0.027` sempre; 0.2196/Tu² nunca explode).
  **Floor efetivo ≈ 0.0293%** devido ao overshoot log(2)/300 do smooth max:
  para Tu ≤ 0.027%, Reθt(Tu) ≈ 1412 em vez dos 1459 do LM2009 com clip duro
  (−3.2%) — diferença sistemática esperada vs. implementações hard-clipped.
- O texto de P&Z (Eqs. 10–14 e 8–9) não enuncia estes limites; eles vêm da
  definição original das correlações de Langtry–Menter (2009). O terceiro
  limite LM2009 (Reθt ≥ 20) **está** no paper (Algorithm 2) e no código
  (`rsaGRreThetaLo = 20`, aplicado em `saGammaRetheta.F90:491`) — por isso
  não consta desta lista.
- **Não são cosméticos — são necessários para robustez**: sem o clamp
  inferior de λ_θ, `exp(−35·λ_θ)` em F1 transborda para λ_θ < −20.3
  (Inf → NaN; perto de estagnação λ_θ bruto atinge essas magnitudes em
  transientes) e F3 fica negativa → Reθt_target negativo. Com o clamp,
  F(λ_θ) ∈ [0.532, 1.268] (cúbica de F3 monótona: discriminante da derivada
  < 0; máximo 0.4677 em λ_θ = −0.1) → target sempre positivo, finito e suave.
  O proximity switch do `smoothMinMax` nunca avalia exp com |arg| ≳ 29.
- Jacobiano DADI (`evalSrcJacBlock`): A(3,3) trata o target como constante,
  logo os clamps nem entram lá — sem conflito resíduo↔Jacobiano; onde o clamp
  satura, a derivada verdadeira é exatamente zero e a aproximação torna-se
  exata. Nota lateral: o Jacobiano usa `min()` duro para o vortLim
  (`saGammaRetheta.F90:2325`) vs. `smoothMinMax` no resíduo (`:514`);
  diferença ≤ log(2)/300·vortLim, afeta só taxa de convergência.

## D3 — Erros no antigo doc destilado de física (doc↔paper, não código↔paper) _______________________________FEITO________________________________________

O código está certo nestes quatro pontos; o destilado — declarado «source of
truth» — está errado e pode induzir tarefas futuras em erro. **Não corrigi o
doc** (fora do âmbito A1; decisão do utilizador).

> **Atualização 2026-07-07:** após um 5.º erro de destilação (Algorithm 2 mal
> transcrito — ver `A2_convergencia.md` §D-A2-5), o utilizador decidiu
> **eliminar o doc destilado**. O paper integral é agora a única referência
> de física; a tabela abaixo fica como registo histórico.

| # | O destilado dizia | Paper integral diz | Código |
|---|---|---|---|
| 1 | §2.3 (Eq. 3): F_θt = min(max(F_wake·e^…, 1−((γ_eff−1/ce2)/(1−1/ce2))²), 1) | Eq. 3: `F_θt = F_wake·e^{−(d/δ)⁴}` (forma simplificada SA-LM2015; ramo γ_eff removido de propósito, §II) | igual ao paper (`saGammaRetheta.F90:564`) |
| 2 | §2.5: F_wake = e^{−(Re_ω/1e5)²} com Re_ω = ρΩy²/μ | Eq. 5: `F_wake = e^{−Re_S/1e6}` com Re_S baseado em strain (substituição de Schücker, §II) | igual ao paper (`saGammaRetheta.F90:563`) |
| 3 | §4.1 (Eq. 52 sLM2015): P_γ = ca1·F_length·Ω_lim·F_onset·(1−ce1·γ) — **sem √γ** | Eq. 52 tem `√γ` | igual ao paper (`saGammaRetheta.F90:522-524`) |
| 4 | §3.3 (Eqs. 54–57): F3 = φ₋₃₀₀(F1,0)+φ₊₃₀₀(F2,0); F = 1+F3, com F1 de sinal tal que o ramo APG desaparece | Eqs. 54–57: F1 = 1+0.275(1−e^{−35λ})e^{−Tu/0.5}; F2 = φ₊₃₀₀(F1,1); F3 = 1−(−12.986λ−123.66λ²−405.689λ³)e^{−(Tu/1.5)^1.5}; F = φ₋₃₀₀(F2,F3) | igual ao paper (`turbUtils.F90:2299-2315`) |

- **Prova (ex.):** `grep -n -F "F_{wake}" "docs/SA_GAMMA_RETHETHA_BASE/Piotrowski_Zingg_2020_SA-sLM2015_clean (1).md"` (linha 155) vs §2.5 do destilado; `sed -n '540,557p'` do paper integral vs §3.3 do destilado.

## D-A2-2 — Forma da restrição Eq. 59: aditiva (DADI) vs MAX (turbKSP) — decidido: manter as duas _______________________________FEITO________________________________________

Data: 2026-07-07. Dúvida levantada em `A2_convergencia.md` §D-A2-2 (a forma
aditiva do DADI é mais restritiva que a Eq. 59 do paper e os dois caminhos
não são numericamente equivalentes — intencional?). **Veredicto: intencional
e correto; nada a alterar.**

- **MAX é a forma do paper e aplica-se onde existe um Δt_CFL real.** O turbKSP
  tem `dtinv = 1/(ANK_CFL·dtl·volRef)` (`NKSolvers.F90:2406`) e aplica
  `max(dtinv, λ/0.9)` (`NKSolvers.F90:2416-2418, 2612-2616`), incluindo o
  switch de desativação do paper via `srcDtRestrictActive`
  (`NKSolvers.F90:2581`). Fiel 1:1 à Eq. 59 e ao §IV.B.3.
- **O DD-ADI não tem Δt_CFL para tomar o MAX.** `dtl` não aparece em nenhum
  solver de `src/turbulence/`; o pseudo-transiente do smoother é definido como
  fração do próprio Jacobiano central (`sa.F90:835-838`: «I/dt = cc/cfl»,
  implementado por `factor = 1/α`). Não existe segundo operando para o `max`
  sem importar `dtl` do escoamento — mudança estrutural sem forma
  correspondente no paper (o paper nem tem caminho DADI).
- **A forma aditiva é o tratamento point-implicit clássico** de fontes
  desestabilizadoras em smoothers de relaxação — mesmo espírito do
  `qq = max(qq, zero)` do SA original (`sa.F90:333`). No DADI,
  `qq(m,m) += λ/0.9` (`saGammaRetheta.F90:1634-1638`) *é* a realização
  natural da Eq. 59, não um desvio. **(Nota 2026-07-07, ver W3 em
  A3_coerencia)** o clip `max(qq, zero)` foi também aplicado diretamente a
  `qq(1,1)` e `qq(2,2)` na fonte SA-GR, antes da soma aditiva — os dois
  mecanismos coexistem: primeiro `max(base,0)`, depois `+ λ/0.9`.
- **Restritividade extra limitada e inofensiva:**
  `max(a,b) ≤ a+b ≤ 2·max(a,b)` ⇒ Δt efetivo no pior caso metade da forma
  MAX, apenas nas células com λ grande. O DADI é fase de arranque —
  conservadorismo aí é barato e desejável.
- **A não-equivalência entre caminhos não é problema:** ambas as formas só
  alteram a diagonal do LHS; o resíduo (RHS) é idêntico ⇒ a solução
  convergida é a mesma nos dois caminhos. DADI e Newton-Krylov já são
  iterações diferentes por construção.
- **Ordem soma-antes-do-factor: também intencional.** Em
  `saGammaRetheta.F90:1635-1637` a soma `λ/0.9` ocorre antes da multiplicação
  por `factor = 1/α` (linha 1650), logo a relaxação implícita amplifica também
  o termo de fonte — efetivo `λ/(0.9·α)`. Coerente com tratar `srcLambda`
  como parte do Jacobiano central (tratamento uniforme com o resto do `qq`);
  satisfaz a Eq. 59 com margem (`λ/(0.9·α) ≥ λ/0.9`). Com o default
  `alfaTurb = 0.8` (`inputParamRoutines.F90:4336`) a margem extra é
  1/α = 1.25 — pior caso combinado ~2.5× menos Δt local que a forma MAX
  literal.
- **Opções de tuning futuro (não defeitos)**, só se o DADI estagnar
  visivelmente na frente de transição: (1) subir `transitionSrcDtLimit`
  (afrouxa sem recompilar); (2) mover a soma para depois do `factor` —
  leitura literal da Eq. 59, mudança de uma linha.
- **Prova:** `grep -rn "dtl" src/turbulence/*.F90` (vazio fora de
  comentários); `sed -n '1629,1660p' src/turbulence/saGammaRetheta.F90`;
  `sed -n '2405,2421p;2580,2620p' src/NKSolver/NKSolvers.F90`;
  `grep -n "alfaTurb = " src/inputParam/inputParamRoutines.F90`.

## D-A2-3 — Switch de desativação da Eq. 59 (§IV.B.3): DADI intencional; turbKSP corrigido _______________________________FEITO________________________________________

Data: 2026-07-07. Dúvida levantada em `A2_convergencia.md` §D-A2-3 (DADI sem
desativação trava a cauda assintótica? turbKSP sem a condição de resíduo
desativa cedo demais?). **Veredicto: DADI correto e fiel ao paper — nada a
alterar; turbKSP divergia do paper em 3 pontos e foi corrigido (código).**

**DADI — «sem desativação» é a leitura fiel do paper, não um desvio:**

- O paper só desativa a restrição na fase **inexact-Newton**; na fase
  approximate-Newton (globalização) fica ativa sem exceção (§IV.B.3). O DADI
  é solver de globalização — não tem fase de Newton nem noção de backtracking
  — logo "sempre ativo" é conforme à filosofia do paper.
- Só altera a diagonal do LHS ⇒ solução convergida idêntica (mesmo argumento
  de D-A2-2). Custo possível: taxa assintótica mais lenta nas células da
  frente de transição em runs DADI-only profundos — o próprio paper avisa que
  os autovalores de fonte «can remain large throughout convergence …
  delaying convergence».
- Controlável em runtime, zero recompilação: `transitionSrcDtRestrict: False`
  desliga; `transitionSrcDtLimit` (default 0.9) afrouxa. **Nada a calibrar:**
  0.9 é o valor calibrado pelo próprio paper (Eq. 59: «0.9 was determined to
  be an effective balance between speed and robustness»); os 5 iters do
  switch idem.
- **Em aberto (→ `docs/TODO.md`):** verificar empiricamente se a cauda do
  DADI fica lenta demais com a restrição sempre ativa.

**turbKSP — 3 divergências vs paper, corrigidas em `NKSolvers.F90` (2026-07-07):**

1. **Contador na fase errada.** Paper: 5 iterações *inexact-Newton* limpas —
   só possíveis depois de R_d < 1e-5 (tolerância do switch de fase). Código:
   contava iterações turbKSP desde o arranque ⇒ desativação prematura com
   resíduo relativo ~1e-1. **Fix:** o contador só incrementa no regime de
   2ª ordem (`totalR ≤ ANK_secondOrdSwitchTol·totalR0`), o análogo interno da
   fase inexact-Newton — mesma condição que a rotina já usava para fluxos de
   2ª ordem + Eisenstat-Walker.
2. **Reativação por subida de resíduo omitida.** O comentário antigo alegava
   «does not track relative residual», mas `totalR`/`totalR0` estavam
   disponíveis na própria rotina. **Fix:** contador reposto a zero sempre que
   `totalR > ANK_secondOrdSwitchTol·totalR0` (cobre "ainda não chegou" e
   "regrediu").
3. **Backtrack bem-sucedido contava como iteração limpa.** `LSFailed` volta a
   `.False.` quando o backtracking é acionado mas reduz a norma (reset antes
   da leitura pelo contador) — contra o paper («reactivated if backtracking
   is required») e contra o próprio comentário do código. **Fix:** novo flag
   `backtrackTriggered`, armado à *entrada* do ramo de backtracking (norma
   unsteady excedida ou NaN), independente do sucesso; é ele que repõe o
   contador.

Implementação mínima, só em `ANKTurbSolveKSP` (`NKSolvers.F90`): declaração
de `backtrackTriggered` (~3505); inicializado/armado junto do line search
(~3672-3684); bloco do contador reescrito (~3778-3792).
`FormJacobianANKTurb`/`FormFunction_mf_turb` inalterados (recalculam a
condição a partir de `noBacktrackCount`). Compila e importa; não é código
diferenciado (sem Tapenade).

**⚠️ Consequência operacional:** com o default `ANKSecondOrdSwitchTol = 1e-16`
o regime de 2ª ordem nunca é atingido ⇒ a restrição **nunca desativa** no
turbKSP (conservador, equivalente a ficar na fase de globalização). Para ter
a aceleração do switch do paper, definir `"ANKSecondOrdSwitchTol": 1e-5`
(valor típico do paper para o switch de fase) ou `1e-4`.

**Erro de doc corrigido de caminho:** `architecture.md` dizia
«`srcDtDeactivateIters = 0` ⇒ never deactivate» — é o oposto (`0 < 0` é
falso ⇒ restrição inativa desde a 1ª iteração no turbKSP). Corrigido na
tabela e no exemplo 8.

Notas residuais (sem ação): `noBacktrackCount` é variável de módulo
(`inputParam.F90:309`) e persiste entre chamadas/fases — com a reativação por
resíduo o estado auto-corrige na 1ª iteração de cada solve. ANK acoplado e NK
não aplicam a restrição de todo — já coberto e fechado em D-A2-1.

- **Prova:** `sed -n '650,658p'` do paper integral (Eq. 59 + §IV.B.3, valor
  0.9 e tolerância 1e-5); `sed -n '3505,3520p;3667,3690p;3775,3795p'
  src/NKSolver/NKSolvers.F90` (código corrigido);
  `grep -n '"ANKSecondOrdSwitchTol"' adflow/pyADflow.py` (default 1e-16,
  linha 5918).

## D-A2-4 — turbResScale de ν̃: 1e4 (ADflow) vs 1e3 (paper §IV.C) — decidido: manter ADflow _______________________________FEITO________________________________________

Data: 2026-07-07. Dúvida levantada em `A2_convergencia.md` §D-A2-4 (default
SA-GR `[1e4, 10, 1e4]` diverge do paper `(1e3, 10, 1e4)` na escala de ν̃ —
seguir ADflow ou paper?). **Veredicto: intencional e correto — convenção
ADflow prevalece; manter `[1e4, 10, 1e4]`.**

- **`turbResScale` não afeta a solução convergida** — é scaling de linha do
  resíduo (condicionamento + comparabilidade de magnitudes); resíduo zero é
  zero com qualquer escala. A regra "paper wins" cobre física, não constantes
  de tuning calibradas para outro código.
- **O valor certo depende da nondimensionalização, que difere entre códigos.**
  O 1e3 do §IV.C foi calibrado para o solver dos autores (velocity-based
  nondim). O ADflow usa p-ρ nondim; a equação de ν̃ no SA-GR é exatamente a
  equação SA do ADflow (γ só multiplica a produção), logo o seu resíduo tem a
  magnitude do SA do ADflow — cuja calibração histórica é 1e4
  (`pyADflow.py:6661`).
- **O padrão nas três escalas confirma decisão deliberada:** γ (10) e Re̅θt
  (1e4) batem com o paper — precisamente as variáveis independentes do esquema
  de nondim (γ adimensional ∈[0,2]; Re̅θt já é um Reynolds). Só ν̃ diverge —
  precisamente a única que carrega a nondim de viscosidade do código. O SST
  logo abaixo (`[1e3, 1e-6]`, `pyADflow.py:6666`) mostra que no ADflow estes
  números são sempre calibração empírica por modelo.
- **Coerência com A3:** com 1e4, a linha de ν̃ do SA-GR fica escalada
  identicamente ao SA puro — mesmas tolerâncias efetivas no coupled ANK/NK e
  históricos de convergência comparáveis SA vs SA-GR.
- **Alternativa "manter ν̃=1e4 e ×10 nas outras" (preservar rácios do paper)
  foi considerada e rejeitada:** o ×10 é específico da nondim do resíduo de
  ν̃, não um shift global transferível — R_γ e R_θt têm fatores próprios.
  Multiplicar γ→100 e θt→1e5 trocaria dois valores calibrados (batem com o
  paper e já correm nos smoke tests / scaling do DADI) por dois palpites, e
  mudaria 10× as tolerâncias efetivas dessas linhas face ao mean-flow.
- **Calibração futura correta (se necessária): medir, não transferir rácios.**
  Imprimir as normas dos resíduos escalados por equação num caso
  representativo e ajustar via opção de runtime `turbresscale` até ficarem
  comparáveis entre si e com o mean-flow. Item registado em `docs/TODO.md`.
- **Prova:** `grep -n -i turbresscale adflow/pyADflow.py` →
  `_updateTurbResScale` (`pyADflow.py:6659-6666`: SA=1e4,
  SA-GR=[1e4,10,1e4], SST=[1e3,1e-6]); paper §IV.C.
