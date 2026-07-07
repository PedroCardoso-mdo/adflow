# A3 — Coerência SA-GR vs SA e SST (código↔código)

Data: 2026-07-06. Tarefa A3 de `state.md`: confirmar que a integração do
SA-γ-Re̅θt no ADflow é coerente com os padrões dos modelos SA e SST já
existentes. **O que faz sentido, manter — nada foi alterado nesta tarefa.**
Este ficheiro lista primeiro o que foi conferido e bate certo, depois as
divergências/dúvidas registadas (sem veredictos onde a decisão é humana).

Nota de âmbito: A1 já confirmou código↔paper; aqui a lente é código↔código
(SA-GR vs `sa.F90` e vs `SST.F90`/infraestrutura multi-equação).

---

## Âmbito conferido — coerente, MANTER

1. **Dispatch (`turbAPI.F90`)** — SA-GR entra nos mesmos dois pontos que os
   outros modelos (`turbSolveDDADI:74-81`, `turbResidual:156-157`), com
   `unsteadyTurbSpectral(itu1, itu3)` análogo ao `(itu1, itu2)` dos 2-eq.
   Prova: `grep -n gammaretheta src/turbulence/turbAPI.F90`.

2. **Estrutura do módulo** — segue o padrão modular do `sa.F90`
   (Source → turbAdvection → unsteadyTurbTerm → Viscous → ResScale → Solve),
   não o monólito do `SST.F90` (que monta tudo dentro de `SSTSolve`). Sensato:
   é o padrão mais amigável ao Tapenade (razão declarada no cabeçalho do
   próprio `sa.F90`). O wrapper `saGammaReTheta_block` (linhas 67-141) tem a
   mesma sequência de `sa_block`/`SST_block`: `bcTurbTreatment` → solve →
   `saEddyViscosity` → `applyAllTurbBCThisBlock(.true.)`.

3. **Parte SA da fonte** — `Source` (saGammaRetheta.F90:299-463) é
   textualmente idêntica a `saSource` (sa.F90:150-302): mesmos gradientes,
   `ss`, fv1/fv2, ft2, `sst`, `rr`/fw, `term1`/`term2`. Diferença única e
   intencional: `gammaForSA = min(max(γ,0),1)` multiplica **só** term1 e
   term2_prod (produção), nunca term2_dest — regra 2 do CLAUDE.md (Eq. 41)
   respeitada. O Jacobiano `qq(1,1)` (linhas 670-674) é o `qq` do SA com γ
   nas derivadas de produção apenas.
   Prova: `diff <(sed -n '299,438p' src/turbulence/saGammaRetheta.F90) <(sed -n '150,289p' src/turbulence/sa.F90)`
   (só whitespace/comentários).

4. **Difusão ν̃** — `Viscous` reproduz `saViscous` exatamente para a equação
   ν̃ (métricas, `cnud` com cb2/cb3, clip `max(cdm+cam, zero)`, tratamento
   implícito das BCs via `bmt`). γ difunde com `(ν + ν_t/σ_f)` e Re̅θt com
   `σ_θt·(ν + ν_t)`, com ν_t = ν̃·fv1 avaliado nas células e mediado nas
   faces — formas do paper (A1) montadas com a mesma discretização métrica
   do SA. As linhas cruzadas do `bmt` (itu1↔itu2↔itu3) entram no tratamento
   implícito das BCs — extensão correta do padrão escalar ao caso multi-eq.
   Prova: comparar saGammaRetheta.F90:875-1047 com sa.F90:382-466.

5. **Advecção/termo instacionário** — usa as rotinas genéricas
   `turbAdvection(3,3,nn,qq)` e `unsteadyTurbTerm(3,3,nn,qq)`; a interface é
   genérica em `mAdv` (turbUtils.F90:825-860: `qq(2:il,2:jl,2:kl,mAdv,mAdv)`),
   tal como os 2-eq usam `2`. First-order upwind para as 3 eqs por opção
   (`transitionFirstOrderUpwind`, regra 5), restaurando `orderTurb` depois.

6. **Solver tridiagonal** — `tdia3x3` (turbUtils.F90:1634) documentado como
   o análogo 3-eq do `tdia3`; mesma estrutura (bloco central cheio,
   off-diagonais escalares).

7. **ResScale** — `ResScale` (1439-1483) = `saResScale` (675-711) estendido
   às 3 variáveis: `dw = -volRef·scratch·iblank`, mesmos comentários.

8. **Relaxação implícita** — `saGammaReThetaSolve` usa o mesmo
   `factor = 1 + (1-alfaTurb)/alfaTurb` do `saSolve` (1597-1599 vs
   sa.F90:845-847) e o mesmo clip do update `w(itu1)=max(w,0)` (linha 2098).

9. **Viscosidade turbilhonar** — `computeEddyViscosity` despacha SA-GR para
   `saEddyViscosity` (turbUtils.F90:637-638): μ_t = ρν̃fv1, sem γ. Coerente
   com o acoplamento só-na-produção (γ não modula μ_t no SA-LM2015).

10. **Vetor de estado / monitores** — `nw=8, nt2=8` (inputParamRoutines:2152-
    2156), extensão genérica do padrão SST (itu1..itu2 → itu1..itu3);
    3 monitores L2 (Nu/Gamma/Retheta, inputParamRoutines:93-97);
    `whalo2(groundLevel, nt1, nt2, ...)` genérico cobre as 3 vars.

11. **BC de parede** (turbBCRoutines.F90:921-983) — ν̃ anti-simétrico
    (`bmt=+1`, halo=−interior), igual ao SA; γ e Re̅θt zero-gradiente
    (`bmt=−1`). Mesma maquinaria `bmt/bvt` do SA e do k-ω/SST.

12. **turbResScale auto** — `[10000, 10, 10000]` (pyADflow.py:6609-6611)
    espelha o `10000.0` escalar do SA para a componente ν̃.

13. **Cobertura dos `select case (turbModel)`** — varrimento de todos os
    ficheiros que tratam `spalartAllmaras` (excl. AD gerado): todos tratam
    também o novo enum, com UMA exceção (turbCurveFits.F90 → W2 abaixo).
    Prova: `for f in $(grep -rl "spalartAllmaras" src/ --include=*.F90 | grep -v adjoint/output | grep -v adjoint/temp); do [ $(grep -ci gammaretheta $f) = 0 ] && echo $f; done`
    → só `src/turbulence/turbCurveFits.F90`.

14. **Blockettes** — SA-GR não está implementado no `blocketteResCore`; o
    wrapper força `useBlockettes=False` para este modelo
    (pyADflow.py:6601-6604). Guarda explícita presente — coerente.

---

## Divergências / dúvidas (registadas, nada alterado)

### W1 — `useft2SA` default True num modelo chamado "noft2" (dúvida principal) _____________________FEITO_____________________________

- O enum é `spalartallmarasnoft2gammaretheta` ("SA-noft2-Gamma-Retheta") e a
  base do paper é SA-noft2, mas a fonte SA-GR mantém o switch `useft2SA`
  (saGammaRetheta.F90:406-410), cujo default é **True** (pyADflow.py:5708) e
  nada o desliga para este modelo (nem no wrapper nem no Fortran).
- Consequência: com opções default, o modelo de transição corre **com** ft2
  activo — contrário ao nome do modelo e à base SA-noft2 do paper. Se o
  utilizador põe `useft2SA=False` nos scripts, não há problema.
- Não alterei: decidir entre (a) forçar ft2=0 para este enum, (b) confiar no
  script do utilizador, é julgamento humano (muda física com defaults).
- Prova: `grep -n useft2 adflow/pyADflow.py src/turbulence/saGammaRetheta.F90`
  e `grep -rn useft2 src/inputParam/`.
-**resposta** Warning adicionadao

### W2 — Wall functions não suportadas e sem guarda _____________________FEITO_____________________________

- `saSolve` tem um bloco inteiro de wall functions (sa.F90:752-833);
  `saGammaReThetaSolve` não tem nenhum (`grep -n wallFunctions
  src/turbulence/saGammaRetheta.F90` → vazio).
- `initCurveFitDataSae` está comentado para o SA-GR
  (inputParamRoutines.F90:2157) e `turbCurveFits.F90` não tem caso para o
  enum (ver ponto 13 acima).
- Fisicamente correto (transição exige BL resolvida, y+≈1), mas não há
  guarda que aborte `wallFunctions=True` + SA-GR — a combinação corre
  silenciosamente sem o enforcement de parede do lado turbulento. Registado
  como dúvida; adicionar guarda é decisão do utilizador.
  -**resposta** Erro adicionado se estas opções forem selecionas juntas

### W3 — Clip de positividade dos diagonais da fonte (RESOLVIDO 2026-07-07) _____________________FEITO_____________________________

- **Antes:** `saSource` clippa o seu diagonal (`qq = max(qq, zero)`,
  sa.F90:333) — filosofia ADflow validada: tratar implicitamente só os
  termos que somam à dominância diagonal; os que desestabilizam vão para o
  RHS explícito. A fonte SA-GR **não** clippava `qq(1,1)` nem `qq(2,2)`;
  só `qq(3,3)` tinha `max(...,zero)` (e esse nunca dispara, pois
  `(1−fThetaT)≥0`).
- **Diagnóstico:**
  - `qq(1,1)` (SA/ν̃) pode ficar negativo pelos termos de destruição
    (`−rsaCw1·dfw`), tal como no SA original.
  - `qq(2,2)` (γ) fica negativo de forma **rotineira** para γ pequeno (zona
    pré-transição): os fatores `(1.5·ce1·γ−0.5)` e `(2·ce2·γ−1)` vão a
    negativo. Este é o caso que **realmente** dispara.
- **Neutro face ao artigo (regra 9):** o `qq` é o Jacobiano implícito no LHS
  do DDADI, não o resíduo. Na convergência o resíduo→0 independentemente do
  LHS, logo a solução estacionária é idêntica — só muda a robustez do
  caminho implícito. Não é equação do paper; é estabilização numérica
  interna do ADflow (a Eq. 59 `srcLambda/limit` é complementar e aditiva).
- **Agora:** `qq(1,1) = max(qq(1,1), zero)` (após o cálculo do diagonal SA)
  e `qq(2,2) = max(qq(2,2), zero)` (após o cálculo do diagonal γ), espelhando
  o SA e o `qq(3,3)` já existente. Ordem final do diagonal SA:
  `max(base,0)` e depois `+ srcLambda/limit` (Eq. 59) quando
  `transitionSrcDtRestrict=True`. Isto fecha o antigo ponto de atenção com
  `transitionSrcDtRestrict=False` — o diagonal SA passa a ter a mesma
  dominância que o solver SA original.
- **AD:** ambos os clips ficam dentro de `#ifndef USE_TAPENADE`; não tocam o
  código gerado (sem `TAPENADE NEEDED`).

### W4 — Farfield inflow do SA-GR alinhado ao padrão ADflow (RESOLVIDO 2026-07-07) _____________________FEITO_____________________________

- **Antes:** SA-GR usava `bvt = 2·wInf`, `bmt = +1` (valor exato `wInf` na
  **face**), enquanto os outros modelos usam ghost = `wInf` diretamente
  (`bvt = wInf`, `bmt = 0`). Não é um desvio "mais forte" — era simplesmente
  o padrão de *inflow* do próprio ADflow (`bcTurbInflow`) aplicado no ramo
  farfield-inflow. As duas formas convergem para o mesmo resultado num
  farfield genuíno (interior → freestream).
- **Nada no paper nem no modelo LM força a forma face-value** (o paper nem
  especifica BCs discretas — usa SBP-SAT). Por decisão do utilizador, o
  SA-GR passa a seguir o padrão testado/validado do ADflow.
- **Agora:** removido o caso especial `if (turbModel == ...gammaretheta)` em
  `bcTurbFarfield`; SA-GR cai no ramo genérico ghost = `wInf` (`bmt = 0`),
  **idêntico a SA/SST**. Também removido o import então inútil de `turbModel`.
  γ_∞ = 1 e Re̅θt,∞ = correlação(Tu∞) continuam a ser impostos corretamente.

### W5 — Inicializações diferentes do padrão (deliberadas, manter) _____________________FEITO_____________________________

- Campo interior γ = 0.02 enquanto γ_∞ = 1 (initializeFlow.F90:2229-2245,
  comentário: suprimir produção SA até Fonset ativar γ). Os outros modelos
  inicializam campo = `wInf`. Verificação vs paper é assunto do utilizador
  (fora do âmbito código↔código).
- `eddyVisInfRatio` default 1e-10 para SA-GR vs 0.009 (SA) / 0.1 (outros)
  (inputParamRoutines.F90:3870-3890) — freestream quase-laminar, coerente
  com um modelo de transição.
-**resposta** de proposito está no to do para verificar se ira causar problemas
  

### W6 — Docs desatualizados (não corrigidos; mesma política do A1-D3)- _____________________FEITO_____________________________

| Doc diz | Código diz |
|---|---|
| `architecture.md` §6: BC de parede de γ é "γ=0 (Dirichlet)" | zero-gradiente/Neumann (turbBCRoutines.F90:931 `bmt=-1`; commit dc1950ef) |
| `architecture.md` §3 + `CLAUDE.md`: helpers em `src/turbulence/saGammaRethetaHelpers.F90` (367 linhas) | ficheiro não existe; correlações + `smoothMinMax` estão em `turbUtils.F90:2279-2410` |
| `architecture.md` §3/§4: `saGammaRetheta.F90` 1862 linhas, solve em 1251-1861 | 2470 linhas; solve em 1485-2131 |

