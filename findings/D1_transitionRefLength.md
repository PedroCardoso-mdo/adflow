# D1 — Limitador de vorticidade: comprimento de referência (RESOLVIDO)

Data: 2026-07-07. Fecho do finding D1 (limitador usava L_ref = 1 m em vez da
corda do paper). Análise completa, decisão de design e implementação da opção
`transitionRefLength`.

## O finding original

- **Código** (antes): `vortLim = uInf·√(uInf/muInf)/20` — equivale à fórmula do
  paper `M·√(M·Re)/20` (Eqs. 52–53) com `l = 1 m`.
- **Paper**: `Re = ρ∞·a∞·l/μ∞` com «the reference length is specified as the
  root chord» (texto integral, linha ~170).
- Divergência: teto físico difere do paper por um fator √(corda[m]).

## Porque NÃO é a regra "drop Re" (a questão-chave)

Objeção natural: «o paper adimensionaliza pela corda, o ADflow por 1 m — não é
tudo consistente?» **Sim para todas as equações de transporte** (invariantes de
escala; o `Re` do paper é artefacto de unidades e cai na conversão — regra §5
do `nondimensionalization.md`). **Não para o limitador**, porque:

1. É um *threshold*, não um termo de equação física. Análise dimensional: um
   teto sobre Ω (1/s) construído só com ρ∞, U∞, μ∞ exige um comprimento, e
   nada na fórmula cancela essa dependência.
2. O próprio paper diz que o √Re aqui é físico, não adimensionalização
   (§IV, linha ~519): «the square root of the freestream Reynolds number …
   is a physical scaling that is independent of the nondimensionalization».
3. Física: vorticidade de parede em CL laminar ∝ √(ρU³/(μ·c)) — o teto do
   paper com l = corda acompanha o nível de vorticidade da geometria real.

Em unidades físicas: paper Ω_lim = (1/20)·√(ρ∞U∞³/(μ∞·c)); código antigo
= idem com c → 1 m. Erro = √(c/1m).

## Avião completo

Nenhum `l` único está "certo" para todos os componentes (asa raiz ~8 m,
empenagem ~2–3 m, slat ~0.3 m; teto local "ideal" ∝ 1/√c_local). A prescrição
do paper (uma corda de raiz global) já é um compromisso, calibrado em perfis
e asas isoladas. A Re de voo o limitador importa menos (transição junto ao
bordo de ataque, LSBs pouco proeminentes). Modo de falha se o teto for alto
demais: reaparecem oscilações de resíduo junto a LSBs — visível, não erro
silencioso de física. Recomendação: MAC ou corda de raiz.

## Implementação (2026-07-07)

| Ficheiro | Alteração |
|---|---|
| `src/modules/inputParam.F90` | `transitionRefLength = -1.0` em `inputIteration` (negativo = auto). |
| `src/turbulence/saGammaRetheta.F90` | `Source` + `evalSrcJacBlock`: `refLenTrans = transitionRefLength` se > 0, senão `lengthRef`; `vortLim = uInf·√(uInf/(muInf·refLenTrans))/20`. |
| `src/f2py/adflow.pyf` | `transitionreflength` exposta. |
| `adflow/pyADflow.py` | Default `[float, -1.0]` + `optionMap` → `["iter", "transitionreflength"]`. |

**Ligação ao AeroProblem — sem plumbing novo**: `setAeroProblemData`
(`pyADflow.py:~3608`) já faz `inputphysics.lengthref = ap.chordRef` em cada
`setAeroProblem`, que corre *antes* de qualquer avaliação de resíduo em todos
os pontos de entrada. `setDefaultValues` semeia `lengthRef = 1.0` no init
(rede de segurança = comportamento antigo). pyADflow dá erro se `chordRef`
faltar no AeroProblem.

Uso: default = corda do AeroProblem (segue trocas de AP automaticamente);
`setOption('transitionRefLength', x)` para fixar; `1.0` = comportamento
pré-opção.

## Estado de verificação

- `make` completo, exit 0; import test OK; `libadflow.inputiteration.transitionreflength = -1.0` confirmado via Python.
- **TAPENADE NEEDED**: ficheiros AD gerados ainda têm a fórmula antiga
  (ver nota em `docs/adjoint-trace.md`).
- Smoke test (N iterações sem NaN) e verificação física: por conta do
  utilizador, como sempre.

## Docs relacionados

- `docs/nondimensionalization.md` §5 — exceção à regra "drop Re" (adicionada 2026-07-07).
- `docs/architecture.md` Parte 2 — entrada da opção + nota de plumbing/guidance.
- `docs/adjoint-trace.md` — nota de regeneração Tapenade pendente.
