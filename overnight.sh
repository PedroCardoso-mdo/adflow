#!/usr/bin/env bash
# overnight.sh -- verificacao SA-gamma-Re_theta, sessao NOVA por tarefa,
# estado duravel em state.md, reinicio ancorado no git. Correr da raiz da repo.
#
#   ./overnight.sh 2>&1 | tee overnight_$(date +%F).log
#
# Parar: Ctrl-C  ou  `touch STOP`.
#
# PORQUE sessao nova por tarefa (nao --continue): contexto limpo. Cada invocacao
# le state.md, faz UMA tarefa, faz commit, actualiza state.md, e sai. Sem deriva
# de contexto acumulado. O estado vive no DISCO (state.md) e no GIT (commits),
# nao na sessao -- por isso sobrevive a rate-limit, crash, ou reboot.
#
# SEGURANCA: --dangerously-skip-permissions da redea solta na working tree.
# Corre numa branch descartavel com tudo commitado.

set -uo pipefail

STATE="state.md"
MAX_TURNS=40
TASK_TIMEOUT="45m"
SLEEP_BETWEEN=300

log(){ echo "[$(date '+%F %T')] $*"; }

if [ ! -f "$STATE" ]; then
cat > "$STATE" <<'EOF'
# Estado do run noturno -- SA-gamma-Re_theta (nTurb=3: nu_tilde + gamma + Re_theta)

Regras deste ficheiro:
- Cada sessao NOVA le este ficheiro primeiro, escolhe a primeira tarefa PENDENTE
  cujas dependencias estao FEITAS, e faz SO essa.
- Marca FEITA **apenas depois** de `git commit` dessa tarefa. Antes disso, EM_CURSO.
- Ao reiniciar, a VERDADE e o `git log`, nao esta lista. Se uma tarefa diz FEITA
  mas nao ha commit dela, trata-a como PENDENTE.
- Duvida / precisa de julgamento humano -> NAO adivinhar: escrever em
  TODO_REVIEW.md, marcar REVER, passar a proxima.
- Teste que nao passa em 2 tentativas -> `git checkout` dessa tarefa, REVER,
  proxima. NUNCA afrouxar tolerancias.

| id  | estado   | tarefa |
|-----|----------|--------|
| A1  | PENDENTE | Confirmar (nao re-auditar) SA-GR vs paper c/ nao-adimensionalizacao (ler docs/nondimensionalization.md 1o: scaling p-rho, vel normaliza a M*sqrt(gamma)). So divergencias em findings/A_confirmacao.md. |
| A2  | PENDENTE | Confirmar estrategia de convergencia (DDADI + solver tipo-ANK) faz sentido vs paper. Duvidas, nao veredictos. |
| A3  | PENDENTE | Confirmar coerencia SA-GR vs SA e SST. O que faz sentido, MANTER. Nao mudar por mudar. |
| B   | PENDENTE | Ler o MD sst_dev. Licoes p/ SA-GR em findings/B_licoes_sst.md: k-correction provavelmente NAO se aplica (SA nao transporta k) -> confirmar no codigo; acoplamento vivo e gamma->fonte-SA sobre kinks (min/max) -> localizar; dot-product prova consistencia, NAO correcao. |
| C   | PENDENTE | Adjoint/Tapenade (usar docs/adjoint-trace.md). Verificar wiring do adjoint e que _d/_b/_fast_b das fontes SA-GR contem o acoplamento gamma->fonte e a correlacao Re_theta. Tabela em findings/C_tapenade.md. NAO diferenciar residuo do adjoint a mao. |
| D1  | PENDENTE | Criar teste transposta (build real), copiando logica de test_jacVecProdFWD.py / test_jacVecProdBWDFast.py: fwd vs reverse e vs reverse-fast, <Ju,v>=<u,J^Tv>. Metrica NORM-RELATIVA (nao rtol por-elemento; MD mostra rtol=41 como anti-padrao). Falha -> localizar alternando acoplamento gamma, nao afrouxar. |
| D2  | PENDENTE | Criar teste forward-AD vs complex-step (REQUER build complexo). Cobrir cl, cd, cmz. Estado apertado; malha c/ mascara near-zero. Se build complexo ausente: PREPARAR tudo, marcar REVER em TODO_REVIEW.md ('D2 precisa build complexo'), PARAR. Nao instalar. |
EOF
log "state.md criado."
fi

for i in $(seq 1 "$MAX_TURNS"); do
  log "===== volta $i / $MAX_TURNS ====="
  [ -f STOP ] && { log "STOP presente; saio."; break; }

  if ! grep -qE '\| PENDENTE \|' "$STATE"; then
    log "sem tarefas PENDENTES em state.md. Terminado."
    break
  fi

  timeout "$TASK_TIMEOUT" claude -p "Es um verificador rigoroso do modelo SA-gamma-Re_theta (Piotrowski & Zingg 2020, SA-sLM2015, nTurb=3) nesta repo ADflow. So tens: esta repo e o MD de analise do branch sst_dev.

ARRANQUE (faz sempre, cada sessao):
1. Le state.md por inteiro. Le docs/README.md e usa a tabela de routing dele para saber que ficheiros abrir.
2. VERIFICA no git log que o estado em state.md e verdade: se uma linha diz FEITA mas nao ha commit correspondente, trata-a como PENDENTE.
3. Escolhe a PRIMEIRA tarefa PENDENTE cujas dependencias (A antes de B antes de C antes de D; D1 antes de D2) estao FEITAS. Faz SO essa tarefa.

EXECUCAO:
- Le apenas os ficheiros que esta tarefa precisa (disciplina de tokens; docs/README.md diz quais).
- Faz o trabalho da tarefa conforme descrito na sua linha em state.md.
- Ao terminar: git add -A && git commit com mensagem descritiva. SO DEPOIS do commit, edita state.md mudando o estado dessa tarefa PENDENTE->FEITA.
- Duvida ou julgamento humano: escreve em TODO_REVIEW.md, muda estado para REVER, NAO facas commit de codigo especulativo.
- Teste que nao passa em 2 tentativas: git checkout dos teus ficheiros dessa tarefa, escreve em TODO_REVIEW.md, estado REVER. NUNCA afrouxes tolerancias nem marques FEITA o que nao passa.
- Nenhum veredicto de correcao que nao consigas provar com um comando.

Faz UMA tarefa e para. Nao encadeies varias numa sessao." \
      --dangerously-skip-permissions
  rc=$?

  if [ "$rc" -ne 0 ]; then
    log "claude saiu rc=$rc (rate-limit, timeout, ou erro). state.md preserva o progresso; re-tento apos sleep."
  fi

  sleep "$SLEEP_BETWEEN"
done

log "loop terminado. Le TODO_REVIEW.md PRIMEIRO, depois git log, depois findings/ e state.md."
