# TODO — SA-gamma-Retheta

Lista de itens de tuning/melhoria futuros, decididos mas adiados. Não são
defeitos: cada um foi analisado e fechado nos `findings/` (referência em cada
item). Só atacar quando houver sintoma concreto (run que estagna, custo
excessivo) ou quando o modelo estiver fisicamente validado.

## Tuning numérico

- [ ] **Calibrar escalas `turbResScale`** (de D-A2-4, `findings/A_confirmacao.md`).
  Default atual `[1e4, 10, 1e4]` — mantido por decisão (convenção ADflow para
  ν̃; γ e Re̅θt do paper). Se um arranque acoplado estagnar numa linha de
  turbulência: **medir, não transferir rácios do paper** — imprimir as normas
  dos resíduos escalados por equação num caso representativo e ajustar via
  opção de runtime `turbresscale` até ficarem comparáveis entre si e com o
  mean-flow. Zero código.


- [ ] **Verificar se a cauda do DADI fica lenta demais com a Eq. 59 sempre
  ativa** (de D-A2-3, `findings/A_confirmacao.md`). No DADI a restrição
  aditiva de fonte nunca desativa (fiel ao paper — fase de globalização), mas
  os autovalores de fonte podem ficar grandes até ao fim ⇒ cauda assintótica
  potencialmente lenta na frente de transição em runs DADI-only profundos.
  Solução convergida não é afetada (só diagonal do LHS). **Teste proposto
  (A/B, zero código):** mesmo caso smoke, DADI-only até L2 profundo, comparar
  o declive da cauda do resíduo turb com (a) `transitionSrcDtRestrict: True`
  (default), (b) `False`, (c) `transitionSrcDtLimit: 2.0–5.0`. Se (a) for
  visivelmente mais lenta: subir `transitionSrcDtLimit` em runtime ou
  desligar a restrição tarde no run; alternativa de 1 linha (soma depois do
  `factor`) já documentada em `A_confirmacao.md` §D-A2-2.

## Desempenho

- [ ] **Blockettes para SA-GR** (de N-A2-6, `findings/A2_convergencia.md`).
  `useBlockettes` forçado a False (`pyADflow.py:6656-6657`) — caminho de
  resíduo mais lento em todos os mat-vecs matrix-free do ANK/NK. Implementar
  em `blocketteResCore` quando se avaliar desempenho.
