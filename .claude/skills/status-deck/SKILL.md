---
name: status-deck
description: Construir/atualizar o PowerPoint de estado do projeto (SA-γ-R̅eθt / SA-BCM) — deck incremental por reunião, com receitas explícitas, proveniência em bibliografia e verificação visual dos slides
---

# Status deck — PowerPoint de reunião

O deck de estado é mostrado uma vez por reunião — e a reunião é **SEMANAL**: um deck por semana, não por sessão de trabalho. Se a semana já tiver um deck intermédio, fundir as secções num único make_deck da semana (o intermédio vai para `absorbed_<data>/` dentro da pasta da semana; ver 2026-09-01). É **incremental e
sistemático**: cada deck novo = deck da reunião anterior COMPLETO + secção
nova no fim. Nunca criar um deck paralelo; nunca perder slides do utilizador
(ex.: LGMRES/adjunto).

## Onde está tudo

- Diretoria: `~/Desktop/Run/MDO_PhD/Transition/gama_rethetha/00_status_decks/`
  — **uma pasta por reunião** (`<data>/`), cada uma com o seu
  `make_deck_<data>.py` e o `SA_gamma_Retheta_status_<data>.pptx` entregue
  (ver o README.md dessa diretoria).
- Script gerador: python-pptx, correr com
  `/home/mdo/packages_v2/mach/bin/python`. Um script por reunião; o script
  novo começa por `Presentation(00_status_decks/<data anterior>/*.pptx)` e
  acrescenta. Nova reunião = nova pasta, copiar o script da anterior como
  ponto de partida e atualizar OLD/OUT.
- As figuras continuam nas campanhas (`15_sabcm_polars/plots/`, etc.) — o
  deck referencia-as por caminho absoluto via `RUN`.
- Helpers no script (copiar do último `make_deck_*.py`): `slide()`,
  `bullets()`, `pic_fit()`, `table()`, `note()`, `kpi()`, `recipe()`,
  `prov()`/`REFS`.

## Regras de conteúdo (pedidos permanentes do utilizador)

1. **Receita de convergência explícita** em todos os slides de resultados —
   linha própria azul bold (`recipe()`) entre o título e as figuras.
2. **Proveniência em bibliografia**: no slide só `Dados: [n]` (canto inferior
   direito, `prov()`); o slide final "Referências de dados [n]" lista job
   Slurm + caminho do funcs.json/log, no HPC ($R) e no PC. `prov()` faz
   `REFS.append()` e o slide de referências itera `REFS`.
   `$R = TransitionModel/SA_GAMMA_RETHETHA` (Deucalion) =
   `~/Desktop/Run/MDO_PhD/Transition/gama_rethetha` (PC).
3. **Comparações têm de ser apples-to-apples**: mesmo nº de ranks, mesma
   malha, mesmo critério de paragem — se uma série mudou de ranks, ou se
   refaz a run ou se declara no slide. Nunca misturar ranks numa tabela sem
   o dizer.
4. **Geometrias/estatísticas separadas por perfil** (min–máx e totais por
   perfil, não agregados); dizer sempre qual é a malha (pyHyp L1, etc.).
5. Números lidos **programaticamente** dos funcs.json/logs sempre que
   possível (menos transcrição manual, proveniência = os próprios ficheiros).
6. cl–Cd e Cl–x_tr são os plots que interessam nas polares; tempo-vs-Cl para
   custo. Tabelas grandes por alfa → preferir plot.
7. Referência "sem nada": nas comparações de custo incluir SA puro; nas de
   warm-start incluir o cold.

## Regras de layout (aprendidas à custa — não repetir erros)

- **Subtítulo de 1 linha** (~100 chars máx a Pt13): a régua azul fica a
  y=1.34" e risca a 2ª linha. Texto longo → `recipe()` (top≈1.44) ou nota.
- Primeira tabela nunca antes de top=1.8 (com recipe: 2.0+).
- Figuras nunca por cima das notas: nota a 6.55–6.85, prov a 7.12.
- `pic_fit()` já preserva aspect ratio — dar-lhe a largura toda (12.1) num
  slide por caso, em vez de espremer 2 figuras largas lado a lado.
- **Apagar slides do deck antigo**: só DEPOIS de adicionar todos os novos
  (remover cedo liberta o partname e o python-pptx duplica slides no zip);
  remover via `prs.part.drop_rel(rId)` + `sldIdLst.remove()` no fim do
  script.
- Trocar uma figura num slide antigo sem tocar no resto: substituição do
  blob da imagem (ver bloco `_swaps` no make_deck).

## Verificação (obrigatória antes de entregar)

1. Rebuild: `mach/bin/python make_deck_<data>.py`.
2. Render: `libreoffice --headless --convert-to pdf` no scratchpad +
   `pdftoppm -png -r 60 -f N -l M` dos slides tocados.
3. **Olhar para os renders** (Read da imagem): sobreposições, texto riscado
   pela régua, tabelas fora do slide, títulos a transbordar.
4. Contar slides e confirmar que não há duplicados
   (`len(prs.slides)` + títulos únicos).
5. Entregar o .pptx com SendUserFile (attach).

## Ao fechar uma reunião

- O deck entregue fica congelado; a próxima reunião cria
  `make_deck_<nova data>.py` que carrega este .pptx.
- Registar nos PURPOSE.md das campanhas os jobs/figuras novos que o deck
  cita, e espelhar scripts/figuras no HPC (rsync).
