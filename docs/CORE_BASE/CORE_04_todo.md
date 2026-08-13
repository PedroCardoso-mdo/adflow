# TODO — SA-gamma-Retheta

Lista de itens de tuning/melhoria futuros, decididos mas adiados. Não são
defeitos: cada um foi analisado e fechado em `CORE_03_design_decisions.md` (referência
em cada item). Só atacar quando houver sintoma concreto (run que estagna,
custo excessivo) ou quando o modelo estiver fisicamente validado.

## Tuning numérico

- [x] **Calibrar escalas `turbResScale`** — **FECHADO 2026-08-12**: o default
  auto no código já É o valor de campanha `[1e4, 0.1, 1e-4]`
  (`pyADflow.py:~6834`, ~1/magnitude de estado por equação, P&Z §IV.1). A
  reconciliação pedida na nota de 2026-07-16 aconteceu no código; o registo
  D-A2-4 em `CORE_03_design_decisions.md` ficou histórico (argumentava pelos
  valores antigos `[1e4, 10, 1e4]`).

## Solver profundo (da campanha de convergência 2026-07-14→16)

- [ ] **Parede do NK abaixo de rel ~5e-9** — lin res 0.8-0.99 com GMRES
  esgotado. **Diagnóstico fechado 2026-08-07/12** (commits `ef9fc10d` +
  estudo ILU): EW-off falsificado, JacLag sem efeito, ILU(3) pior — é o PC,
  não há alavanca de opções. Itens de código por ordem de impacto em
  `SA_GAMMA_RETHETHA_BASE/SAGR_02_adflow_vs_paper_solver.md` §8 ("Open code items"):
  PC mais forte (agora o item decisivo), physicality check ρ/E no NK. Case
  dir: `03_convergence_strategy/3d_plain_wing/solver_code_items/`. Testbed:
  `best_strategy/restarts/r3_deepest_record_rel3.3e-9_dp.cgns`.
- [ ] **`ANKNSubiterTurb` é knob morto com `ANKUseTurbDADI=True`**
  (NKSolvers.F90: o loop de subiterações só existe no ANKTurbSolveKSP; o
  ramo DADI é uma chamada única) — ligar ou documentar na opção.
- [ ] **Flip-flop ANK<->SANK** no arranque da variante SANK (switch dispara
  a rel ~5e-4 e oscila ~3 min) — histerese ou limiar dedicado poupariam
  esse custo.


- [ ] **Verificar se a cauda do DADI fica lenta demais com a Eq. 59 sempre
  ativa** (de D-A2-3, `CORE_03_design_decisions.md`). No DADI a restrição
  aditiva de fonte nunca desativa (fiel ao paper — fase de globalização), mas
  os autovalores de fonte podem ficar grandes até ao fim ⇒ cauda assintótica
  potencialmente lenta na frente de transição em runs DADI-only profundos.
  Solução convergida não é afetada (só diagonal do LHS). **Teste proposto
  (A/B, zero código):** mesmo caso smoke, DADI-only até L2 profundo, comparar
  o declive da cauda do resíduo turb com (a) `transitionSrcDtRestrict: True`
  (default), (b) `False`, (c) `transitionSrcDtLimit: 2.0–5.0`. Se (a) for
  visivelmente mais lenta: subir `transitionSrcDtLimit` em runtime ou
  desligar a restrição tarde no run; alternativa de 1 linha (soma depois do
  `factor`) já documentada em `CORE_03_design_decisions.md` §D-A2-2.

## Desempenho

- [ ] **Blockettes para SA-GR** (de N-A2-6, `CORE_03_design_decisions.md`).
  `useBlockettes` forçado a False (`pyADflow.py:~6824`). **Atualização
  2026-07-24:** os kernels blockette SA-GR JÁ estão implementados,
  re-sincronizados e testados (`test_blockette_sagr.py`,
  `task-log/2026-07-24-blockette-sagr-residual-sync.md`) — falta só decidir
  levantar o force-off e medir desempenho.
  
  
  
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
- [x] **Wire the complete-mode `test_adjoint`** — DONE 2026-07-23/24
  (`test_adjoint_sagr.py`, estágio `adjoint` do `run_sagr_tests.sh`; real
  passa, CS a 5e-8 com mach/drag non-blocking —
  `task-log/2026-07-23-sagr-full-adjoint-test.md` + `2026-07-24-…`).
- [ ] **(cosmetic)** delete the now-superseded `dev/sanity_check_*` /
  `dev/check_3way_fwd.py` once their mesh/`--crossflow` flags are no longer
  wanted. ~~rename `_flatplate` refs~~ — done 2026-07-23 (`*_sagr_tut_wing.json`).

## Adjoint / partials (de `../ARCHIVE/ARCHIVE_03_adjoint_audit_2026-07-07.md`, 2026-07-07)

- [x] **Rerun Tapenade** (`uInf, muInf` no head do `Source`) — **FEITO**:
  `Makefile_tapenade:187` tem `uInf, muInf`, `vortlimd = 0.0_8` já não
  existe nos ficheiros gerados, e os AD files estão commitados (última
  regen 2026-08-04, `2c4ce2c1`). Verificado 2026-08-12.
- [ ] **Decisão em aberto (default = diferenciar):** cap do limitador de
  vorticidade diferenciado vs congelado ("frozen limiter"). Se preferir
  congelar: reverter a linha do `Makefile_tapenade` e documentar; dR/dw é
  idêntico nas duas opções. Análise: `ADFLOW_BASE/ADFLOW_09_adjoint_trace.md` header.
- [ ] **Ao testar partials:** (a) se o dot-product BWDFast falhar, suspeitar
  primeiro do stripping push/pop do `autoEditReverseFast.py` (partiu o
  fast_b do SST upstream — foi desativado lá, continua ativo neste branch);
  (b) perto dos pontos de blend do `smoothMinMax`/vortLim esperar ruído de
  FD, não erro de AD — usar complex-step ou máscara, **nunca** inflar rtol
  global (post-mortem SST: rtol 2.11/41 = não-verificação);
  (c) validar cd/cm e todos os DVs com CS, não só cl (lacuna do SST).
