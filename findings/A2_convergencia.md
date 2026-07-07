# A2 — Estratégia de convergência (DDADI + solver tipo-ANK) vs paper §IV

Data: 2026-07-06. Tarefa A2 de `state.md`: confirmar que a estratégia de
convergência implementada faz sentido face a Piotrowski & Zingg (2020) §IV.
**Formato: dúvidas, não veredictos.** Fontes: paper integral
`docs/SA_GAMMA_RETHETHA_BASE/Piotrowski_Zingg_2020_SA-sLM2015_clean (1).md` §IV
(originalmente via um doc destilado, entretanto eliminado — ver D-A2-5),
`docs/architecture.md` Parte 1, `docs/ADFLOW_BASE/ADFLOW_01_flow_solver_theory.md`,
e leitura dirigida de `src/turbulence/saGammaRetheta.F90` e
`src/NKSolver/NKSolvers.F90`.

## Enquadramento — o que bate certo (não repetido nas dúvidas)

- **Estratégia global coerente com o paper.** O paper (§IV.B) oferece 3 níveis
  de acoplamento e recomenda fully-coupled para convergência profunda, sendo o
  decoupled o mais fácil de estabilizar. ADflow mapeia isto em: arranque
  decoupled (DADI com 3 modos via `TurbDADICoupled`, ou Turb-ANK KSP), ANK
  acoplado opcional, e NK terminal fully-coupled. Também coerente com a
  filosofia própria do ADflow (ADFLOW_01: decoupling ganha no arranque mesmo em
  escoamento fortemente acoplado; coupled reservado à fase terminal).
- **Eq. 59, fator 0.9**: `transitionSrcDtLimit` default 0.9
  (`adflow/pyADflow.py:5696`).
- **Alg. 2 (damping)**: presente no DADI (back-off exponencial θ^m por
  variável, secção de update de `saGammaReThetaSolve` em
  `saGammaRetheta.F90`) com os bounds do paper (γ∈[1e-10,2], Re̅θt≥20 — já
  conferidos em A1); desde 2026-07-07 com teto efetivamente ilimitado e clip
  só como último recurso com warning (ver D-A2-5). Nos caminhos ANK, o papel
  do Alg. 2 é
  desempenhado pelo physicality check (`physicalityCheckANK{,Turb}`) com os
  mesmos bounds (`NKSolvers.F90:3219,3228`).
- **Scaling linha/coluna (§IV.C)**: scaling simétrico no DADI
  (`saGammaRetheta.F90:1582-1584,1631-1654`) via `turbResScale`; escalas de
  γ (10) e Re̅θt (1e4) iguais às do paper.
- **Atalho de autovalores A32=0 é legítimo para a forma implementada.** O
  paper lista J(3,2)=∂R_θt/∂γ (via F_θt(γ)), mas a F_θt implementada é a
  simplificada SA-LM2015, `fThetaT = fWake·exp(-(y/δ)^4)`
  (`saGammaRetheta.F90:564`), sem γ — logo A32=0 e o espectro
  {A33} ∪ {bloco 2×2} (`saGammaRetheta.F90:2403`) é exato. As formas LM2015
  selecionáveis em runtime afetam F_turb (linha 2 do Jacobiano), não a linha 3.

---

## Dúvida D-A2-1 — Restrição de Δt por fonte ausente nos caminhos acoplados ____________________________________fEITO__________________________________________

- **Código:** `grep -n srcLambda src/NKSolver/NKSolvers.F90` → só
  `FormJacobianANKTurb` (2417), `FormFunction_mf_turb` (2614) e
  `ANKTurbSolveKSP` (3515-3611). O ANK acoplado (8 eq.) e o NK terminal não
  aplicam a Eq. 59; contam apenas com physicality check + rampa de CFL.
- **Paper:** o solver dos autores é fully-coupled, e a restrição Eq. 59 é
  exatamente o dispositivo de robustez desse modo (§IV.B).
- **Dúvida:** espera-se que o ANK acoplado sobreviva à rigidez dos termos-fonte
  de transição sem a Eq. 59? Se os arranques acoplados divergirem, este é o
  primeiro suspeito. Decisão de desenho a validar por humano. 
- **Resposta** Ponto ignorado, não se espera utilizar CANK warnings foram adicionados ao caso

## Dúvida D-A2-2 — Forma da restrição difere entre caminhos (aditiva vs MAX)  ____________________________________fEITO__________________________________________

- **Código:** DADI usa forma **aditiva** — `qq(m,m) += srcLambda/0.9`, i.e.
  dtinv_eff = dtinv + λ/0.9 (`saGammaRetheta.F90:1620-1629`, comentário
  explícito). O turbKSP usa forma **MAX** sobre o dtinv_CFL existente
  (`NKSolvers.F90:2416-2417,2613-2614`).
- **Paper:** Eq. 59 é `Δt = min(Δt_CFL, 0.9/λ_max)` ⇒ em dtinv, forma MAX.
- **Dúvida:** a forma aditiva do DADI é *mais* restritiva que a do paper
  (soma > max). Conservador para estabilidade, mas pode atrasar a convergência
  do DADI; e os dois caminhos não são numericamente equivalentes entre si.
  Intencional?
- **Resposta** (2026-07-07) Fechado: intencional e correto, manter as duas
  formas. MAX segue o paper onde existe Δt_CFL (turbKSP, incl. switch de
  desativação); a aditiva é o tratamento point-implicit padrão onde não há
  (DD-ADI não tem `dtl` — o pseudo-transiente é `cc/cfl`). Ambas só mexem na
  diagonal do LHS ⇒ solução convergida idêntica; penalidade máxima 2× no Δt
  local (×1.25 adicional da ordem soma-antes-do-factor, com α=0.8). Análise
  completa, prova e opções de tuning em `A_confirmacao.md` §D-A2-2.

## Dúvida D-A2-3 — Switch de desativação (§IV.B.3) implementado parcialmente ____________________________________FEITO__________________________________________

- **Paper:** desativar após 5 iterações inexact-Newton sem backtracking **e**
  R_d > 1e-5; reativar se houver backtracking **ou** o resíduo subir.
- **Código:** (a) DADI — **sem desativação de todo** (comentário
  `saGammaRetheta.F90:1624`: «controlled by transitionSrcDtRestrict only (no
  deactivation)»); (b) `ANKTurbSolveKSP` — conta só iterações sem backtracking
  (`noBacktrackCount`, `NKSolvers.F90:3777-3788`); a condição de resíduo é
  omitida (comentário 3781: «does not track relative residual — only
  backtracking») e a reativação dá-se só por backtrack, não por subida de
  resíduo.
- **Dúvida:** no DADI, faz sentido não desativar (não há noção de backtracking
  ali), mas então a restrição aditiva fica ativa até ao fim — trava a
  convergência assintótica do DADI? No turbKSP, omitir R_d > 1e-5 pode
  desativar cedo demais num arranque que ainda está longe de convergido.
- **Resposta** (2026-07-07) DADI: correto e fiel ao paper — a restrição nunca
  é desativada na fase de globalização (approximate-Newton), e o DADI é
  exatamente essa fase; só altera a diagonal do LHS (solução convergida
  idêntica), custo = cauda mais lenta (verificação empírica →
  `docs/TODO.md`). turbKSP: divergia do paper em 3 pontos (contador corria na
  fase errada, sem reativação por subida de resíduo, backtrack bem-sucedido
  contava como iteração limpa) — **corrigido em `NKSolvers.F90`**
  (2026-07-07); a desativação só atua com `ANKSecondOrdSwitchTol` ~1e-5
  (default 1e-16 = nunca). Análise completa e implementação em
  `A_confirmacao.md` §D-A2-3.

## Dúvida D-A2-4 — Escala de ν̃ no turbResScale: 1e4 vs 1e3 do paper  ____________________________________FEITO__________________________________________

- **Código:** default SA-GR `[10000.0, 10.0, 10000.0]`
  (`adflow/pyADflow.py:6659-6666`, `_updateTurbResScale`).
- **Paper:** §IV.C indica (ν̃, γ, Re̅θt) = (1e3, 10, 1e4).
- **Dúvida:** 1e4 é a convenção histórica do ADflow para SA (default SA puro =
  1e4, e ADFLOW_01 usa 1e4 no coupled), pelo que manter 1e4 é provavelmente
  deliberado e coerente com A3 (coerência com SA). Mas divergem do paper por
  10×; confirmar com humano que a convenção ADflow prevalece aqui
  (nondimensionalização difere, portanto a magnitude "certa" de R_ν̃ também).
- **Resposta** (2026-07-07) Fechado: convenção ADflow prevalece, manter
  `[1e4, 10, 1e4]`. `turbResScale` é tuning numérico (não altera a solução
  convergida); o 1e3 do paper foi calibrado para outra nondim, e a equação de
  ν̃ é a do SA do ADflow (calibração histórica 1e4). γ e Re̅θt batem com o
  paper por serem independentes da nondim. Alternativa de preservar rácios
  (×10 em γ e θt) rejeitada — o ×10 é específico do ν̃. Calibração empírica
  futura registada em `docs/TODO.md`. Análise completa em
  `A_confirmacao.md` §D-A2-4.

## Dúvida D-A2-5 — Damping por variável, não por atualização inteira ____________________________________FEITO__________________________________________

- **Código:** DADI amortece γ e Re̅θt **independentemente** (dois ciclos de
  back-off separados, `saGammaRetheta.F90:2109-2135`) e não amortece ν̃ (só
  `max(w,0)` na linha 2107); após o ciclo há ainda clip duro min/max
  (2122-2123, 2135).
- **Dúvida original:** o antigo doc destilado apresentava o Algorithm 2 como
  **um** fator de damping comum a todo o ΔQ, o que faria do esquema
  por-variável do código uma divergência de algoritmo.
- **Resposta** (2026-07-07) **Fechado: premissa incorreta — o código segue o
  paper.** O pseudocódigo real do Algorithm 2 (paper integral, linhas 622-637)
  amortece **por variável e por nó**: dois ciclos `while` independentes, um
  para γ (bounds [1e-10, 2]) e um para Re̅θt (≥20), com o contador `m`
  **reposto a zero entre eles** (linha 9 do pseudocódigo), cada um com o seu
  θ_fac^m; ν̃ nunca entra no damping. É exatamente a estrutura do código. O
  «fator único» era um erro de destilação (5.º, juntando-se aos 4 de
  `A_confirmacao.md` §D3); o doc destilado foi **eliminado** nesta data e o
  paper integral é agora a única referência de física (índices/routing
  atualizados em `CLAUDE.md`, `docs/README.md`,
  `docs/SA_GAMMA_RETHETHA_BASE/README.md`).
- **Divergência residual corrigida (2026-07-07): clip degradado a fallback.**
  O `while` do paper é ilimitado e sem clip — o paper rejeita explicitamente
  clipping duro (linha 620: «placing a hard limit on variable bounds can
  potentially lead to stalling of the nonlinear solver»). O código limitava
  o back-off a 40 iterações (0.99⁴⁰ ≈ 0.67 ⇒ redução máxima ~33%) e clipava
  incondicionalmente no fim — para overshoots grandes degenerava no clipper
  que o paper descartou. **Alteração:** default `transitionDampMaxIter`
  40 → **10000** (0.99¹⁰⁰⁰⁰ ≈ 0 ⇒ back-off efetivamente ilimitado como no
  paper; opção exposta em Python, `pyADflow.py`); o clip `min/max` passou a
  **último recurso**, só alcançável se o back-off esgotar o teto — o que
  exige que o estado anterior já estivesse fora dos bounds (θ^m·Δ → 0) —
  e nesse caso um warning agregado (contagem de células γ/Re̅θt) aconselha
  subir `transitionDampMaxIter` ou investigar a violação a montante
  (`saGammaRetheta.F90`, secção de update de `saGammaReThetaSolve`). Em
  regime normal o comportamento é exatamente o Algorithm 2. O `max(w,0)`
  do ν̃ mantém-se: comportamento padrão do DADI SA em ADflow (regra 2 —
  não mexer no SA) e o paper também não amortece ν̃.
- **Âmbito por solver:** o back-off+clip era **só do DADI**. Os caminhos
  ANK/turbKSP impõem os mesmos bounds via `physicalityCheckANK{,Turb}`
  (`NKSolvers.F90:3280-3470`): limitam o **passo global** λ (allreduce min)
  para que γ/Re̅θt fiquem nos bounds — sem clip pós-update do estado; o
  papel do Alg. 2 já era desempenhado corretamente aí. Sem alteração.

## Nota N-A2-6 — Blockettes desativados para SA-GR                             ____________________________________FEITO__________________________________________

- `adflow/pyADflow.py:6597-6598` força `useBlockettes=False` para
  `SA-noft2-Gamma-Retheta` («not implemented in blocketteResCore»). Não afeta
  a estratégia de convergência em si, mas todos os produtos matriz-vetor
  matrix-free do ANK/NK usam o caminho de resíduo mais lento. Expectativa de
  custo por iteração maior; registar para quando se avaliar desempenho.
- **Resposta** Ponto ignoradoblockets poderão ser adicionados no fururo
---

**Prova (comandos usados):**
- `grep -n "srcLambda\|evalSrcJac" src/NKSolver/NKSolvers.F90`
- `grep -n "srcDtRestrictActive\|backtrack" src/NKSolver/NKSolvers.F90`
- `sed -n '1610,1660p;2090,2130p' src/turbulence/saGammaRetheta.F90`
- `grep -n "fThetaT" src/turbulence/saGammaRetheta.F90`
- `grep -n "transitionSrcDtLimit\|turbresscale" adflow/pyADflow.py`
