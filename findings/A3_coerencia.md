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

### W1 — `useft2SA` default True num modelo chamado "noft2" (dúvida principal)

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

### W2 — Wall functions não suportadas e sem guarda

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

### W3 — `qq(1,1)` sem o clip de positividade do SA (deliberado, manter)

- `saSource` clippa o diagonal: `qq = max(qq, zero)` (sa.F90:333). A fonte
  SA-GR **não** clippa `qq(1,1)` (saGammaRetheta.F90:670-674; só `qq(3,3)`
  tem `max(...,zero)` explícito).
- Mitigação existente: com `transitionSrcDtRestrict=True` (default),
  `srcLambda/limit` é somado ao diagonal (Eq. 59, saGammaReThetaSolve:
  1625-1629) — papel estabilizador análogo, alinhado com P&Z (regra 9:
  paper wins). Além disso o update tem bounds (ν̃≥0; γ/Re̅θt damping Alg. 2).
- Ponto de atenção apenas para `transitionSrcDtRestrict=False`: o diagonal
  SA pode ficar menos dominante do que no solver SA original. Manter;
  registado para contexto de debugging futuro.

### W4 — Farfield inflow impõe valor na face, não no ghost (deliberado, manter)

- Outros modelos: ghost = `wInf` diretamente (turbBCRoutines.F90:475-491).
  SA-GR: `bvt = 2·wInf`, `bmt = +1` → ghost = 2·wInf − interior, i.e. valor
  exato `wInf` na **face** (extrapolação linear; turbBCRoutines.F90:441-471,
  comentário no código explica o porquê para γ=1).
- Enforcement mais forte/preciso que o padrão; deliberado e inócuo. Manter.

### W5 — Inicializações diferentes do padrão (deliberadas, manter)

- Campo interior γ = 0.02 enquanto γ_∞ = 1 (initializeFlow.F90:2229-2245,
  comentário: suprimir produção SA até Fonset ativar γ). Os outros modelos
  inicializam campo = `wInf`. Verificação vs paper é assunto do utilizador
  (fora do âmbito código↔código).
- `eddyVisInfRatio` default 1e-10 para SA-GR vs 0.009 (SA) / 0.1 (outros)
  (inputParamRoutines.F90:3870-3890) — freestream quase-laminar, coerente
  com um modelo de transição.

### W6 — Docs desatualizados (não corrigidos; mesma política do A1-D3)- _____________________FEITO_____________________________

| Doc diz | Código diz |
|---|---|
| `architecture.md` §6: BC de parede de γ é "γ=0 (Dirichlet)" | zero-gradiente/Neumann (turbBCRoutines.F90:931 `bmt=-1`; commit dc1950ef) |
| `architecture.md` §3 + `CLAUDE.md`: helpers em `src/turbulence/saGammaRethetaHelpers.F90` (367 linhas) | ficheiro não existe; correlações + `smoothMinMax` estão em `turbUtils.F90:2279-2410` |
| `architecture.md` §3/§4: `saGammaRetheta.F90` 1862 linhas, solve em 1251-1861 | 2470 linhas; solve em 1485-2131 |

### W7 — Cruft menor (não mexido)

- `saGammaReTheta_block` importa `SSTEddyViscosity` que nunca usa
  (saGammaRetheta.F90:73; usa `saEddyViscosity`). Inócuo.
