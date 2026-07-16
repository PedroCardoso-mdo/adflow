# TODO — SA-gamma-Retheta

Lista de itens de tuning/melhoria futuros, decididos mas adiados. Não são
defeitos: cada um foi analisado e fechado em `audits/design-decisions.md` (referência
em cada item). Só atacar quando houver sintoma concreto (run que estagna,
custo excessivo) ou quando o modelo estiver fisicamente validado.

## Tuning numérico

- [ ] **Calibrar escalas `turbResScale`** (de D-A2-4, `audits/design-decisions.md`).
  Default atual `[1e4, 10, 1e4]` — mantido por decisão (convenção ADflow para
  ν̃; γ e Re̅θt do paper). Se um arranque acoplado estagnar numa linha de
  turbulência: **medir, não transferir rácios do paper** — imprimir as normas
  dos resíduos escalados por equação num caso representativo e ajustar via
  opção de runtime `turbresscale` até ficarem comparáveis entre si e com o
  mean-flow. Zero código.
  **Nota 2026-07-16:** a campanha de convergência usou consistentemente
  `[1e4, 0.1, 1e-4]` com sucesso (ver `convergence-strategy.md`) — o item
  passa a ser "reconciliar o default com o valor de campanha", não calibrar
  do zero.

## Solver profundo (da campanha de convergência 2026-07-14→16)

- [ ] **Parede do NK abaixo de rel ~5e-9** — lin res 0.8-0.99 com GMRES
  esgotado, independente do ponto de engate/JacLag/subespaço. Itens de
  código por ordem de impacto em `adflow-vs-paper-solver.md` §5: Alg. 2
  por nó, reativação Eq. 59 dentro do NK, PC mais forte. Case dir e branch
  preparados (`.../3D_Plain_Wing/ADFLOW_SA_GAMMA_RETHETHA_SOLVER/`,
  `sa_gamma_rethetha_paper_solver`). Testbed: `best_strategie/restarts/r3`.
- [ ] **`ANKNSubiterTurb` é knob morto com `ANKUseTurbDADI=True`**
  (NKSolvers.F90: o loop de subiterações só existe no ANKTurbSolveKSP; o
  ramo DADI é uma chamada única) — ligar ou documentar na opção.
- [ ] **Flip-flop ANK<->SANK** no arranque da variante SANK (switch dispara
  a rel ~5e-4 e oscila ~3 min) — histerese ou limiar dedicado poupariam
  esse custo.


- [ ] **Verificar se a cauda do DADI fica lenta demais com a Eq. 59 sempre
  ativa** (de D-A2-3, `audits/design-decisions.md`). No DADI a restrição
  aditiva de fonte nunca desativa (fiel ao paper — fase de globalização), mas
  os autovalores de fonte podem ficar grandes até ao fim ⇒ cauda assintótica
  potencialmente lenta na frente de transição em runs DADI-only profundos.
  Solução convergida não é afetada (só diagonal do LHS). **Teste proposto
  (A/B, zero código):** mesmo caso smoke, DADI-only até L2 profundo, comparar
  o declive da cauda do resíduo turb com (a) `transitionSrcDtRestrict: True`
  (default), (b) `False`, (c) `transitionSrcDtLimit: 2.0–5.0`. Se (a) for
  visivelmente mais lenta: subir `transitionSrcDtLimit` em runtime ou
  desligar a restrição tarde no run; alternativa de 1 linha (soma depois do
  `factor`) já documentada em `audits/design-decisions.md` §D-A2-2.

## Desempenho

- [ ] **Blockettes para SA-GR** (de N-A2-6, `audits/design-decisions.md`).
  `useBlockettes` forçado a False (`pyADflow.py:6656-6657`) — caminho de
  resíduo mais lento em todos os mat-vecs matrix-free do ANK/NK. Implementar
  em `blocketteResCore` quando se avaliar desempenho.
  
  
  
## Fisica


- [ ] **Valor de inicialização de gamma e retheta** para já está a 0.02 diferente da BC, testar o que dá melhor resultado

## Adjoint / partials (de `audits/adjoint_audit_2026-07-07.md`, 2026-07-07)

- [ ] **Rerun Tapenade** para apanhar `uInf, muInf` ativos no head
  `saGammaRetheta%Source` (`Makefile_tapenade`) → `vortlimd` deixa de ser
  hard-zero; depois `make` e commit dos ficheiros gerados (os 6 já
  regenerados + os novos ficam num só estado consistente).
- [ ] **Decisão em aberto (default = diferenciar):** cap do limitador de
  vorticidade diferenciado vs congelado ("frozen limiter"). Se preferir
  congelar: reverter a linha do `Makefile_tapenade` e documentar; dR/dw é
  idêntico nas duas opções. Análise: `adjoint-trace.md` header.
- [ ] **Ao testar partials:** (a) se o dot-product BWDFast falhar, suspeitar
  primeiro do stripping push/pop do `autoEditReverseFast.py` (partiu o
  fast_b do SST upstream — foi desativado lá, continua ativo neste branch);
  (b) perto dos pontos de blend do `smoothMinMax`/vortLim esperar ruído de
  FD, não erro de AD — usar complex-step ou máscara, **nunca** inflar rtol
  global (post-mortem SST: rtol 2.11/41 = não-verificação);
  (c) validar cd/cm e todos os DVs com CS, não só cl (lacuna do SST).
