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
| A1  | FEITA    | Confirmar (nao re-auditar) SA-GR vs paper c/ nao-adimensionalizacao (ler docs/nondimensionalization.md 1o: scaling p-rho, vel normaliza a M*sqrt(gamma)). So divergencias em findings/A_confirmacao.md. |
| A2  | FEITA    | Confirmar estrategia de convergencia (DDADI + solver tipo-ANK) faz sentido vs paper. Duvidas, nao veredictos. |
| A3  | PENDENTE | Confirmar coerencia SA-GR vs SA e SST. O que faz sentido, MANTER. Nao mudar por mudar. |
| B   | PENDENTE | Ler o MD sst_dev. Licoes p/ SA-GR em findings/B_licoes_sst.md: k-correction provavelmente NAO se aplica (SA nao transporta k) -> confirmar no codigo; acoplamento vivo e gamma->fonte-SA sobre kinks (min/max) -> localizar; dot-product prova consistencia, NAO correcao. |
| C   | PENDENTE | Adjoint/Tapenade (usar docs/adjoint-trace.md). Verificar wiring do adjoint e que _d/_b/_fast_b das fontes SA-GR contem o acoplamento gamma->fonte e a correlacao Re_theta. Tabela em findings/C_tapenade.md. NAO diferenciar residuo do adjoint a mao. |
| D1  | PENDENTE | Criar teste transposta (build real), copiando logica de test_jacVecProdFWD.py / test_jacVecProdBWDFast.py: fwd vs reverse e vs reverse-fast, <Ju,v>=<u,J^Tv>. Metrica NORM-RELATIVA (nao rtol por-elemento; MD mostra rtol=41 como anti-padrao). Falha -> localizar alternando acoplamento gamma, nao afrouxar. |
| D2  | PENDENTE | Criar teste forward-AD vs complex-step (REQUER build complexo). Cobrir cl, cd, cmz. Estado apertado; malha c/ mascara near-zero. Se build complexo ausente: PREPARAR tudo, marcar REVER em TODO_REVIEW.md ('D2 precisa build complexo'), PARAR. Nao instalar. |
