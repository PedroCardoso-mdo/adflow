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

## Test infrastructure

- [x] **Restructure `tests/reg_tests/` to match upstream ADflow's official
  test conventions/fixtures** — DONE 2026-07-22
  (`task-log/2026-07-22-sagr-test-suite-standardization.md`). The ladder is
  now registered testflo tests (`test_jacVecProdFWD_sagr.py`,
  `test_jacVecProdBWDFast_sagr.py` + `reg_sagr.py` + `refs/*.json`), driven by
  `run_sagr_tests.sh` (`all|real|cs|train|genw`). The old one-off scripts moved
  to `tests/reg_tests/dev/` (documented in `dev/README.md`); `w` is produced by
  `dev/generate_sagr_restart.py` (reads `reg_sagr` config), JSON by
  `run_sagr_tests.sh train`. Mesh swap = edit `reg_sagr.py` → `genw` → `train`
  → run. Remaining sub-items below.
- [ ] **Wire the complete-mode `test_adjoint` (total `dF/dX`) in, `@skip`ped**
  with a per-mesh reason: this AR5 state doesn't converge deeply enough for
  total-sensitivity validation. Re-enable on a better-converged mesh. This is
  the one gap between "partials validated" and "gradient validated".
- [ ] **(cosmetic)** delete the now-superseded `dev/sanity_check_*` /
  `dev/check_3way_fwd.py` once their mesh/`--crossflow` flags are no longer
  wanted; rename `_flatplate` refs → `_sagr` (the case is a wing).

## Adjoint / partials (de `audits/adjoint_audit_2026-07-07.md`, 2026-07-07)

- [ ] **Rerun Tapenade** para apanhar `uInf, muInf` ativos no head
  `saGammaRetheta%Source` (`Makefile_tapenade`) → `vortlimd` deixa de ser
  hard-zero; depois `make` e commit dos ficheiros gerados (os 6 já
  regenerados + os novos ficam num só estado consistente).
- [ ] **Decisão em aberto (default = diferenciar):** cap do limitador de
  vorticidade diferenciado vs congelado ("frozen limiter"). Se preferir
  congelar: reverter a linha do `Makefile_tapenade` e documentar; dR/dw é
  idêntico nas duas opções. Análise: `VERIFICATION/adjoint-trace.md` header.
- [ ] **Ao testar partials:** (a) se o dot-product BWDFast falhar, suspeitar
  primeiro do stripping push/pop do `autoEditReverseFast.py` (partiu o
  fast_b do SST upstream — foi desativado lá, continua ativo neste branch);
  (b) perto dos pontos de blend do `smoothMinMax`/vortLim esperar ruído de
  FD, não erro de AD — usar complex-step ou máscara, **nunca** inflar rtol
  global (post-mortem SST: rtol 2.11/41 = não-verificação);
  (c) validar cd/cm e todos os DVs com CS, não só cl (lacuna do SST).
