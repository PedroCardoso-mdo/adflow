# A2 — Estratégia de convergência (DDADI + solver tipo-ANK) vs paper §IV

Data: 2026-07-06. Tarefa A2 de `state.md`: confirmar que a estratégia de
convergência implementada faz sentido face a Piotrowski & Zingg (2020) §IV.
**Formato: dúvidas, não veredictos.** Fontes: `docs/SA_GAMMA_RETHETHA_BASE/paper-reference.md`
§7, `docs/architecture.md` Parte 1, `docs/ADFLOW_BASE/ADFLOW_01_flow_solver_theory.md`,
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
- **Alg. 2 (damping)**: presente no DADI (back-off exponencial θ^m,
  `saGammaRetheta.F90:2100-2126`) com os bounds do paper (γ∈[1e-10,2],
  Re̅θt≥20 — já conferidos em A1). Nos caminhos ANK, o papel do Alg. 2 é
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

## Dúvida D-A2-1 — Restrição de Δt por fonte ausente nos caminhos acoplados

- **Código:** `grep -n srcLambda src/NKSolver/NKSolvers.F90` → só
  `FormJacobianANKTurb` (2417), `FormFunction_mf_turb` (2614) e
  `ANKTurbSolveKSP` (3515-3611). O ANK acoplado (8 eq.) e o NK terminal não
  aplicam a Eq. 59; contam apenas com physicality check + rampa de CFL.
- **Paper:** o solver dos autores é fully-coupled, e a restrição Eq. 59 é
  exatamente o dispositivo de robustez desse modo (§IV.B).
- **Dúvida:** espera-se que o ANK acoplado sobreviva à rigidez dos termos-fonte
  de transição sem a Eq. 59? Se os arranques acoplados divergirem, este é o
  primeiro suspeito. Decisão de desenho a validar por humano.

## Dúvida D-A2-2 — Forma da restrição difere entre caminhos (aditiva vs MAX)

- **Código:** DADI usa forma **aditiva** — `qq(m,m) += srcLambda/0.9`, i.e.
  dtinv_eff = dtinv + λ/0.9 (`saGammaRetheta.F90:1620-1629`, comentário
  explícito). O turbKSP usa forma **MAX** sobre o dtinv_CFL existente
  (`NKSolvers.F90:2416-2417,2613-2614`).
- **Paper:** Eq. 59 é `Δt = min(Δt_CFL, 0.9/λ_max)` ⇒ em dtinv, forma MAX.
- **Dúvida:** a forma aditiva do DADI é *mais* restritiva que a do paper
  (soma > max). Conservador para estabilidade, mas pode atrasar a convergência
  do DADI; e os dois caminhos não são numericamente equivalentes entre si.
  Intencional?

## Dúvida D-A2-3 — Switch de desativação (§IV.B.3) implementado parcialmente

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

## Dúvida D-A2-4 — Escala de ν̃ no turbResScale: 1e4 vs 1e3 do paper

- **Código:** default SA-GR `[10000.0, 10.0, 10000.0]`
  (`adflow/pyADflow.py:6600-6602`).
- **Paper:** §IV.C indica (ν̃, γ, Re̅θt) = (1e3, 10, 1e4).
- **Dúvida:** 1e4 é a convenção histórica do ADflow para SA (default SA puro =
  1e4, e ADFLOW_01 usa 1e4 no coupled), pelo que manter 1e4 é provavelmente
  deliberado e coerente com A3 (coerência com SA). Mas divergem do paper por
  10×; confirmar com humano que a convenção ADflow prevalece aqui
  (nondimensionalização difere, portanto a magnitude "certa" de R_ν̃ também).

## Dúvida D-A2-5 — Damping por variável, não por atualização inteira

- **Código:** DADI amortece γ e Re̅θt **independentemente** (dois ciclos de
  back-off separados, `saGammaRetheta.F90:2100-2126`) e não amortece ν̃ (só
  `max(w,0)` na linha 2098); após o ciclo há ainda clip duro min/max
  (2113-2114, 2126).
- **Paper:** Algorithm 2 aplica **um** fator de damping comum a todo o ΔQ
  (reduzido até todas as variáveis caberem nos bounds), sem clip final.
- **Dúvida:** amortecer por variável muda a direção da atualização face ao
  paper (deixa de ser um escalamento uniforme do passo de Newton/relaxação).
  Em DADI (relaxação, não Newton) isso é provavelmente benigno e até menos
  dissipativo; mas é uma divergência de algoritmo a registar, não a corrigir
  por iniciativa própria.

## Nota N-A2-6 — Blockettes desativados para SA-GR

- `adflow/pyADflow.py:6597-6598` força `useBlockettes=False` para
  `SA-noft2-Gamma-Retheta` («not implemented in blocketteResCore»). Não afeta
  a estratégia de convergência em si, mas todos os produtos matriz-vetor
  matrix-free do ANK/NK usam o caminho de resíduo mais lento. Expectativa de
  custo por iteração maior; registar para quando se avaliar desempenho.

---

**Prova (comandos usados):**
- `grep -n "srcLambda\|evalSrcJac" src/NKSolver/NKSolvers.F90`
- `grep -n "srcDtRestrictActive\|backtrack" src/NKSolver/NKSolvers.F90`
- `sed -n '1610,1660p;2090,2130p' src/turbulence/saGammaRetheta.F90`
- `grep -n "fThetaT" src/turbulence/saGammaRetheta.F90`
- `grep -n "transitionSrcDtLimit\|turbresscale" adflow/pyADflow.py`
