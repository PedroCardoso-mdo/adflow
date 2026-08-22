# ADflow RANS-SA ANK+NK: Mapa de Código, Estratégia de Profiling e Diagnóstico do Blockette

**Data:** 2026-03-15  
**Escopo:** RANS SA, solver ANK+NK, blockette ativo, opções Python por defeito  
**Ficheiros-chave de suporte:**
- `doc/parallelization_report.md`
- `doc/parallelization_hotspot_roadmap.md`
- `doc/solver_pass_nk_ank_sa_pc_optimization.md`

---

## Objectivo deste documento

Este documento:
1. Mapeia **todo o caminho de código** desde a chamada Python até ao output, para RANS-SA com ANK+NK e blockette activo.
2. Explica **quantas vezes** cada secção é chamada por iteração não-linear.
3. Fornece uma **estratégia concreta de timers** para medir tempo de computação e comunicação em cada fase.
4. Descreve **ferramentas de profiling externas** (sem modificar o código) para cache misses e bandwidth.
5. Guia o utilizador a **descobrir por si mesmo** onde o tempo é gasto, em vez de concluir à partida.

A ideia central é: mede primeiro, conclui depois.

---

## Parte 1: Mapa Completo do Caminho de Código

### 1.1 Visão de alto nível

Para uma chamada típica `CFDSolver(aeroProblem)` em modo steady RANS-SA com ANK+NK:

```
Python: ADflow.__call__()                    [adflow/pyADflow.py:1185]
  │
  ├─ PRÉ-SOLVER
  │    ├─ setAeroProblem()                   [pyADflow.py:3240]
  │    │    ├─ actualiza parâmetros de voo (Mach, alpha, Re)
  │    │    └─ inicializa flow state se mudança de AP
  │    │
  │    ├─ updateGeometryInfo()               [pyADflow.py ~ L1244]
  │    │    ├─ mesh warp (se DVs geométricos)
  │    │    ├─ shiftCoordAndVolumes          [Fortran: preprocessing]
  │    │    ├─ computeMetrics               [Fortran: preprocessingAPI]
  │    │    └─ updateWallDistances           [Fortran: wallDistance]
  │    │
  │    └─ (checks de falha fatal antes de solver)
  │
  ├─ SOLVER
  │    └─ adflow.solvers.solver()           [src/solver/solvers.F90 ~ L1]
  │         └─ solveState()                [src/solver/solvers.F90:892]
  │              └─ [detalhado em Parte 1.2]
  │
  └─ PÓS-SOLVER
       ├─ writeSolution()                   [pyADflow.py:2664]
       │    ├─ writeCGNS / writeHDF5
       │    └─ writeSurface / writeLift / writeSlices
       │
       └─ (printTiming se printTiming=True)
```

---

### 1.2 Fortran: `solveState()` — estrutura principal

Ficheiro: `src/solver/solvers.F90:892`

```
solveState()
  │
  ├─ [INIT]
  │    ├─ setCycleStrategy()                 — define ciclos MG
  │    ├─ computeResidualNK()               [L:1013] — residual inicial
  │    │    └─ blocketteRes()               — [VER PARTE 1.4]
  │    ├─ timeStep()                        — computa dt local (spectral radii)
  │    └─ convergenceInfo()                 — MPI_Allreduce + print
  │
  ├─ [LOOP PRINCIPAL: do while (approxTotalIts < nMGCycles)]
  │    │
  │    ├─ [FASE RK/DADI — enquanto totalR > ANK_switchTol * totalR0]
  │    │    └─ executeMGCycle()             [src/solver/multiGrid.F90]
  │    │         ├─ coarse grid RK sweeps
  │    │         ├─ residual_block()        [src/solver/residuals.F90]
  │    │         ├─ whalo1/whalo2()         — halo exchange
  │    │         └─ timeStep()
  │    │
  │    ├─ [FASE ANK — enquanto ANK_switchTol > totalR > NK_switchTol]
  │    │    └─ ANKStep(firstANK)            [NKSolvers.F90:3629]
  │    │         └─ [VER PARTE 1.3.A]
  │    │
  │    ├─ [FASE NK — quando totalR <= NK_switchTol * totalR0]
  │    │    └─ NKStep(firstNK)             [NKSolvers.F90:512]
  │    │         └─ [VER PARTE 1.3.B]
  │    │
  │    └─ convergenceInfo()                — MPI_Allreduce + print
  │
  └─ [PÓS-LOOP]
       └─ (dealloc NKLSFuncEvals se necessário)
```

**Nota sobre as fases:** Por defeito com opções standard RANS-SA e ambos ANK+NK activos
(`useANKSolver=True`, `useNKSolver=True`), a sequência é:
1. MG/RK até ~2 ordens de convergência
2. ANK/sANK até ~6 ordens
3. NK para os últimos passos (resolução quadrática)

---

### 1.3.A Fortran: `ANKStep()` — sequência detalhada

Ficheiro: `src/NKSolver/NKSolvers.F90:3629`

```
ANKStep(firstCall)
  │
  ├─ [SETUP — apenas no firstCall]
  │    └─ setupANKSolver()                 [L:1811]
  │         ├─ MatCreateMFFD(...)          — cria operador matrix-free
  │         ├─ KSPCreate(...)
  │         └─ setupANKKSP()              — configura GMRES + PC
  │
  ├─ [INICIAL] blocketteRes(useUpdateIntermed=.True.)  [L:3698]
  │    └─ [computação completa RANS residual + copia intermediários]
  │
  ├─ [PC REBUILD — condicional em mod(ANK_iter, ANK_jacobianLag)==0]
  │    └─ FormJacobianANK()               [L:1935]
  │         ├─ setupStateResidualMatrix(useAD=ANK_ADPC, ...)
  │         │    ├─ setFDReference()                    — 1 blocketteRes()
  │         │    └─ loop sobre (nColors × nVars):
  │         │         └─ blocketteRes() por cor+variável — N_colors × nVars chamadas
  │         │   TOTAL PC: 1 + nColors × nVars = ~20-40 blocketteRes() calls
  │         ├─ computeTimeStepMat()        [L:2041]
  │         └─ MatAXPY(dRdwPre, timeStepMat) — adiciona diagonal pseudo-transient
  │
  ├─ [RESOLUÇÃO LINEAR]
  │    ├─ MatMFFDSetBase(ANK_dRdw, wVec, baseRes)  [L:3908] — base cached
  │    │    └─ 1 call a blocketteRes() para base
  │    │
  │    └─ KSPSolve(ANK_KSP, rVec, deltaW)         [L:3912]
  │         └─ [GMRES — N_ksp_iters iterações]
  │              por iteração GMRES:
  │              ├─ MatMult (MFFD): blocketteRes(R(w+εv))  — 1 blocketteRes()
  │              ├─ MPI_Allreduce (norms/dots GMRES)        — COMUNICAÇÃO GLOBAL
  │              └─ PC apply (ASM/ILU local)                — comunicação overlap
  │
  ├─ [TURBULÊNCIA — se uncoupled]
  │    ├─ (opção A) turbSolveDDADI()       — DADI turb sweeps
  │    └─ (opção B) ANKTurbSolveKSP()     — KSP separado para νt
  │         └─ MFFD turb: blocketteRes(useFlowRes=.False.) por iter KSP
  │
  ├─ [PHYSICALITY CHECK + BACKTRACK]
  │    └─ (se falhar) blocketteRes() adicional por passo de backtrack
  │
  └─ [FINAL]
       └─ blocketteRes(useUpdateIntermed=.True.)   [L:4064]
            — atualiza intermediários em memória com estado final do passo
```

**Opções por defeito relevantes para ANK:**

| Opção Python | Default | Efeito |
|---|---|---|
| `ANKUseCoupled` | `True` | Turb e flow acoplados (1 KSP) |
| `ANKUseMatrixFree` | `True` | Operador MFFD (não monta J explícito) |
| `ANKADPC` | `False` | PC via FD+coloring (não AD) |
| `ANKJacobianLag` | `20` | Rebuild PC a cada 20 iter |
| `ANKSecondOrdSwitchTol` | `1e-16` | Switcha para 2a ordem na tolerância |
| `ANKUseDissApprox` | `True` | Dissipação approx na PC (não no operador) |
| `ANKUseViscApprox` | `True` | Viscoso approx na PC |

---

### 1.3.B Fortran: `NKStep()` — sequência detalhada

Ficheiro: `src/NKSolver/NKSolvers.F90:512`

```
NKStep(firstCall)
  │
  ├─ [SETUP — apenas no firstCall]
  │    └─ setupNKSolver()                  [L:84]
  │         ├─ MatCreateMFFD(...)          — operador NK matrix-free
  │         ├─ MatCreateShell (dRdwPseudo) — shell para Jv + (1/CFL)v
  │         └─ KSPCreate(...)
  │
  ├─ [RESIDUAL INICIAL — apenas no firstCall]
  │    ├─ setwVec(wVec)
  │    ├─ computeResidualNK()             [L:543] — 1 blocketteRes()
  │    └─ setRVec(rVec)                   — popula vector PETSc
  │
  ├─ [PC REBUILD — mod(NK_iter, NK_jacobianLag)==0]
  │    └─ FormJacobianNK()               [L:562]
  │         └─ setupStateResidualMatrix(useAD=NK_ADPC, ...)
  │             TOTAL PC: ~20-40 blocketteRes() calls (FD coloring)
  │
  ├─ [SETUP BASE MFFD]
  │    ├─ formFunction_mf(ctx, wVec, baseRes) [L:630] — 1 blocketteRes()
  │    └─ MatMFFDSetBase(dRdW, wVec, baseRes) — caches base
  │
  ├─ [RESOLUÇÃO LINEAR]
  │    └─ KSPSolve(NK_KSP, rVec, deltaW)     [L:634]
  │         └─ [GMRES — N_ksp_iters iterações]
  │              por iteração GMRES:
  │              ├─ NKMatMult(dRdwPseudo, x, y):    [L:244]
  │              │    ├─ MFFD matvec: blocketteRes(R(w+εx))  — 1 blocketteRes()
  │              │    └─ y += (1/NK_CFL) * x           — adição do termo pseudo-transient
  │              ├─ MPI_Allreduce (norms GMRES)         — COMUNICAÇÃO GLOBAL
  │              └─ PC apply (ASM/ILU)
  │
  ├─ [LINE SEARCH]
  │    └─ LSNone / LSCubic / LSNM
  │         └─ por retrocesso: 1-3+ calls a computeResidualNK() → blocketteRes()
  │
  └─ [VecCopy, update wVec]
```

**Opções por defeito relevantes para NK:**

| Opção Python | Default | Efeito |
|---|---|---|
| `NKJacobianLag` | `20` | Rebuild PC a cada 20 iter NK |
| `NKADPC` | `False` | PC via FD+coloring |
| `NKSubspaceSize` | `60` | Máximo de iterações GMRES |
| `NKLineSearch` | `cubicLineSearch` | Tipo de line search |
| `NKUseEW` | `True` | Tolerância adaptativa Eisenstat-Walker |

---

### 1.4 Fortran: `blocketteRes()` — dispatch da avaliação de residual

Ficheiro: `src/NKSolver/blockette.F90:70`

Esta é a função que faz UMA avaliação completa do residual RANS-SA.
É chamada a partir de todos os contextos acima.

```
blocketteRes(useDissApprox, useViscApprox, useUpdateIntermed, ...)
  │
  ├─ [PRÉ-HALO — loop sobre todos os blocos]
  │    ├─ setPointers(nn, level, sps)
  │    ├─ computePressureSimple(.False.)   [flowUtils.F90:867]
  │    ├─ computeLamViscosity(.False.)     [flowUtils.F90:1201]
  │    ├─ computeEddyViscosity(.False.)    [turbUtils.F90:580]  ← SA: calcula μt
  │    ├─ BCTurbTreatment                 — BCs turb (SA)
  │    ├─ applyAllTurbBCthisblock(.True.) — aplica BCs turb
  │    └─ applyAllBC_block(.True.)        — aplica BCs flow
  │
  ├─ [HALO EXCHANGE]                      ← COMUNICAÇÃO MPI
  │    └─ whalo2(1, lStart, lEnd, .True., .True., .True.)
  │         ├─ pack buffer (local)
  │         ├─ MPI_Isend / MPI_Irecv (non-blocking)
  │         └─ MPI_Waitany (sincronização)
  │
  ├─ [RE-APPLY BCs após halo — apenas se overset presente]
  │
  └─ [LOOP PRINCIPAL: sps → nn]
       ├─ setPointers(nn, level, sps)
       ├─ [se useBlockettes=True]
       │    └─ blocketteResCore(...)      [L:272] ← KERNEL CRÍTICO [VER PARTE 1.5]
       └─ [se useBlockettes=False]
            └─ blockResCore(...)          [L:299] ← fallback sem tiling
```

**Separação compute vs. comunicar em blocketteRes:**
- `computePressureSimple` + `computeLamViscosity` + `computeEddyViscosity`: **computação pura, sem MPI**
- `applyAllBC_block`: **computação local**
- `whalo2`: **MPI puro** — é o único ponto de comunicação directa dentro de cada residual call
- `blocketteResCore`: **computação intensiva** — é onde a aritmética CFD acontece

---

### 1.5 Fortran: `blocketteResCore()` — o kernel interno

Ficheiro: `src/NKSolver/blockette.F90:299`  
Tile size: `BS = 8` → tiles de `9×9×9` (com halos)

```
blocketteResCore(dissApprox, viscApprox, updateIntermed, flowRes, turbRes, storeWall)
  │
  └─ !$OMP PARALLEL DO private(i,j,k,l) collapse(2)        ← único loop OpenMP
       do kk = 2, bkl, BS              ← tile em k
         do jj = 2, bjl, BS            ← tile em j (collapse com k → paralelo)
           do ii = 2, bil, BS          ← tile em i (SERIAL — os threads dividem j×k)
             │
             ├─ [COPY-IN — carregar tile para memória privada de thread]
             │    ├─ SE ii==2 (primeira tile em i): copia TUDO de memória global
             │    │    ├─ w(0:ib, 0:jb, 0:kb, 1:nw)       — estado (7 vars SA)
             │    │    ├─ p(0:ib, ...)                      — pressão
             │    │    ├─ gamma, ss                         — γ e shock sensor
             │    │    ├─ rlv, rev, vol                     — μ_l, μ_t, vol
             │    │    └─ x(0:ie, ...)                      — coordenadas nó
             │    │
             │    └─ SE ii>2 (tile subsequente): reutiliza face final da tile anterior
             │         ├─ w(0:3,...) ← w(BS:BS+3,...) do anterior     — reuse
             │         ├─ rlv,rev,vol,aa,dss idem                     — reuse
             │         ├─ gradientes nodais (ux,uy,...,qz) idem        — reuse
             │         └─ depois preenche i=4..ib da memória global    — load parcial
             │
             ├─ [COPIA porI, porJ, porK, iblank, d2wall, volRef]
             ├─ [COPIA sFaceI/J/K se grid velocities]
             │
             ├─ fw = zero                                   — limpa buffer flux viscoso
             │
             ├─ call metrics                                — face normals/areas (sI, sJ, sK)
             ├─ call initRes(lStart, lEnd)                  — zero dw
             │
             ├─ [SE RANSEquations .and. turbRes]            ← SA turbulência
             │    ├─ call saSource                          — termo fonte SA (destruição/produção)
             │    ├─ call saAdvection                       — advecção de ν̃
             │    ├─ call saViscous                         — fluxo viscoso SA
             │    └─ call saResScale                        — escalamento do residual SA
             │
             ├─ call timeStep(updateIntermed)               — dt_l, spectral radii (radI,J,K)
             │
             ├─ [SE flowRes]
             │    ├─ call inviscidCentralFlux               — JST central (Euler)
             │    │
             │    ├─ [dissipação invíscida]
             │    │    ├─ inviscidDissFluxScalar/Matrix      — se dissScalar/Matrix
             │    │    └─ inviscidUpwindFlux                 — se upwind
             │    │
             │    └─ [SE viscous]
             │         ├─ call computeSpeedOfSoundSquared   — aa = √(γp/ρ)
             │         ├─ [SE não viscApprox]
             │         │    └─ call allNodalGradients        — ux,uy,...,qz por Gauss
             │         └─ call viscousFlux(storeWall)        — fluxos viscosos completos
             │              ou call viscousFluxApprox        — se viscApprox=True
             │
             ├─ call sumDwAndFw                             — acumula fluxos em dw
             │
             ├─ [COPY-OUT — escreve dw de volta em memória global]
             │    └─ bdw(...) ← dw(i,j,k,l)
             │
             └─ [SE updateIntermed — copy-out adicional]
                  ├─ bdtl ← dtl                             — dt local
                  ├─ bradi/j/k ← radi/j/k                  — spectral radii
                  └─ [SE viscous]
                       ├─ baa ← aa                          — speed of sound sq
                       └─ bux..bqz ← ux..qz                — gradientes nodais (12 arrays 3D)
```

**Dados THREADPRIVATE por thread (declarados em `blockette.F90:64-67`):**

| Array | Dimensões | Conteúdo |
|---|---|---|
| `w` | `(0:bbib, 0:bbjb, 0:bbkb, 1:nw)` | estado flow (ρ, ρu, ρv, ρw, ρE, SA) |
| `p, gamma, ss` | `(0:bbib, ...)` | pressão, γ, shock sensor |
| `rlv, rev` | `(1:ie, 1:je, 1:ke)` | μ_lam, μ_turb |
| `vol, volRef, aa` | idem | volume, ref, c² |
| `dw` | `(bbib, bbjb, bbkb, nw)` | residual buffer |
| `fw` | idem | flux viscoso buffer |
| `sI, sJ, sK` | faces | áreas/normais |
| `ux..qz` | `(bbil, bbjl, bbkl)` | 12 arrays grad nodal (u,v,w,q × x,y,z) |
| `x` | nós `(0:ie, 0:je, 0:ke, 3)` | coordenadas nodais |
| `dtl, radI/J/K, dss` | células | dt local, spectral radii, dissipação |

---

### 1.6 Pós-solver e output

```
Python: writeSolution()                   [pyADflow.py:2664]
  ├─ writeCGNS volume solution           — I/O CGNS (pode ser lento em HPC)
  ├─ writeSurface                        — superfície
  ├─ writeLiftDistribution               — distribuição de sustentação
  └─ writeSlices                         — planos de corte
```

O output I/O não é geralmente um bottleneck comparado com o solver, mas pode ser
relevante para casos com muitas iterações de `writeSolEachIter=True`.

---

## Parte 2: Tabela de Chamadas por Iteração

### 2.1 Contagem de `blocketteRes()` por iteração não-linear

Para uma iteração ANK típica com N_ksp iterações GMRES:

| Evento | # blocketteRes calls |
|---|---|
| `blocketteRes` inicial (updateIntermed=True) | 1 |
| PC rebuild (se lag atingido): base FD | 1 |
| PC rebuild: nColors × nStateVars perturbações | ~20 a 50 |
| Base MFFD (`MatMFFDSetBase`) | 1 |
| Por iteração GMRES (MFFD matvec) | N_ksp |
| Turb KSP (se uncoupled): base + N_ksp_turb | 1 + N_ksp_turb |
| Backtrack (se physicality check falhar) | 0–3 |
| Final `blocketteRes` (updateIntermed=True) | 1 |
| **TOTAL (sem PC rebuild, N_ksp=10)** | **~15** |
| **TOTAL (com PC rebuild, N_ksp=10)** | **~40–70** |

Para uma iteração NK típica:

| Evento | # blocketteRes calls |
|---|---|
| Residual inicial (firstCall) | 1 |
| PC rebuild (se lag atingido) | ~20–50 |
| Base MFFD (`formFunction_mf`) | 1 |
| Por iteração GMRES | N_ksp |
| Line search (cubic: 2–4 evals) | 2–4 |
| **TOTAL (sem PC rebuild, N_ksp=15)** | **~20** |
| **TOTAL (com PC rebuild, N_ksp=15)** | **~50–80** |

**Observação:**  
Cada `blocketteRes()` chama internamente:
- 1× `whalo2` (comunicação MPI)
- 1× `blocketteResCore` (computação, paralelizada por OpenMP)
- + loops de setupPressure/viscosity por bloco

---

## Parte 3: Estratégia de Timers no Código Fortran

### 3.1 Princípios

1. Usar `MPI_WTIME()` (disponível via `mpi`, resolução ~microsegundos) para wall time.
2. Medir tempo de **cada secção** separadamente: computação e comunicação.
3. Contar chamadas com contadores inteiros em módulo.
4. Ao fim de cada iteração não-linear, colecionar e imprimir (ou acumular numa tabela) via rank 0.
5. No final da corrida, imprimir tabela de totais.

### 3.2 Módulo de timers proposto

Criar `src/utils/adflowTimers.F90`:

```fortran
module adflowTimers
    use precision
    use mpi
    implicit none

    ! --- Contadores de chamadas ---
    integer(kind=intType) :: cnt_blocketteRes       = 0
    integer(kind=intType) :: cnt_blocketteResCore   = 0
    integer(kind=intType) :: cnt_whalo2             = 0
    integer(kind=intType) :: cnt_whalo1             = 0
    integer(kind=intType) :: cnt_formJacobianANK    = 0
    integer(kind=intType) :: cnt_formJacobianNK     = 0
    integer(kind=intType) :: cnt_KSPSolve_ANK       = 0
    integer(kind=intType) :: cnt_KSPSolve_NK        = 0
    integer(kind=intType) :: cnt_convergenceInfo    = 0
    integer(kind=intType) :: cnt_ANKStep            = 0
    integer(kind=intType) :: cnt_NKStep             = 0

    ! --- Acumuladores de tempo (segundos) ---
    real(kind=realType) :: t_blocketteRes       = 0.0_realType
    real(kind=realType) :: t_blocketteResCore   = 0.0_realType
    real(kind=realType) :: t_whalo2             = 0.0_realType
    real(kind=realType) :: t_whalo1             = 0.0_realType
    real(kind=realType) :: t_formJacobianANK    = 0.0_realType
    real(kind=realType) :: t_formJacobianNK     = 0.0_realType
    real(kind=realType) :: t_KSPSolve_ANK       = 0.0_realType
    real(kind=realType) :: t_KSPSolve_NK        = 0.0_realType
    real(kind=realType) :: t_convergenceInfo    = 0.0_realType
    real(kind=realType) :: t_setupPC            = 0.0_realType
    real(kind=realType) :: t_preSolverSetup     = 0.0_realType

    ! --- Separação compute vs. comunicação em blocketteRes ---
    real(kind=realType) :: t_blocketteRes_compute = 0.0_realType
    real(kind=realType) :: t_blocketteRes_comm    = 0.0_realType

contains

    subroutine printTimerSummary(comm)
        integer, intent(in) :: comm
        integer :: myID, nProcs, ierr
        real(kind=realType) :: global_t(14), local_t(14)

        call MPI_Comm_rank(comm, myID, ierr)
        call MPI_Comm_size(comm, nProcs, ierr)

        local_t(1)  = t_blocketteRes
        local_t(2)  = t_blocketteResCore
        local_t(3)  = t_whalo2
        local_t(4)  = t_whalo1
        local_t(5)  = t_formJacobianANK
        local_t(6)  = t_formJacobianNK
        local_t(7)  = t_KSPSolve_ANK
        local_t(8)  = t_KSPSolve_NK
        local_t(9)  = t_convergenceInfo
        local_t(10) = t_setupPC
        local_t(11) = t_preSolverSetup
        local_t(12) = t_blocketteRes_compute
        local_t(13) = t_blocketteRes_comm
        local_t(14) = 0.0_realType  ! reservado

        call MPI_Reduce(local_t, global_t, 14, adflow_real, MPI_MAX, 0, comm, ierr)

        if (myID == 0) then
            write(*,'(A)') ''
            write(*,'(A)') '+=================================================================+'
            write(*,'(A)') '|              ADflow Timer Summary (wall, max rank)              |'
            write(*,'(A)') '+-------------------------+----------+-----------+-----------------+'
            write(*,'(A)') '| Section                 |  Total(s)| Calls     | Avg/call (ms)   |'
            write(*,'(A)') '+-------------------------+----------+-----------+-----------------+'
            call printRow('blocketteRes',         global_t(1),  cnt_blocketteRes)
            call printRow('blocketteResCore',     global_t(2),  cnt_blocketteResCore)
            call printRow('  └ compute',          global_t(12), cnt_blocketteResCore)
            call printRow('whalo2',               global_t(3),  cnt_whalo2)
            call printRow('whalo1',               global_t(4),  cnt_whalo1)
            call printRow('  └ comm subtotal',    global_t(13), cnt_whalo2)
            call printRow('FormJacobianANK (PC)', global_t(5),  cnt_formJacobianANK)
            call printRow('FormJacobianNK (PC)',  global_t(6),  cnt_formJacobianNK)
            call printRow('KSPSolve ANK',         global_t(7),  cnt_KSPSolve_ANK)
            call printRow('KSPSolve NK',          global_t(8),  cnt_KSPSolve_NK)
            call printRow('convergenceInfo',      global_t(9),  cnt_convergenceInfo)
            call printRow('setupPC (total)',       global_t(10), cnt_formJacobianANK+cnt_formJacobianNK)
            call printRow('preSolverSetup',        global_t(11), 1_intType)
            write(*,'(A)') '+=================================================================+'
        end if
    end subroutine printTimerSummary

    subroutine printRow(label, t_total, n_calls)
        character(len=*), intent(in) :: label
        real(kind=realType), intent(in) :: t_total
        integer(kind=intType), intent(in) :: n_calls
        real(kind=realType) :: avg_ms
        if (n_calls > 0) then
            avg_ms = t_total / real(n_calls, realType) * 1000.0_realType
        else
            avg_ms = 0.0_realType
        end if
        write(*,'(A,A25,A,F10.3,A,I11,A,F17.3,A)') &
            '| ', label, ' |', t_total, ' |', n_calls, ' |', avg_ms, ' |'
    end subroutine printRow

end module adflowTimers
```

### 3.3 Onde inserir os timers

#### Em `blocketteRes()` — `src/NKSolver/blockette.F90:70`

```fortran
subroutine blocketteRes(...)
    use adflowTimers
    real(kind=realType) :: t_start, t_comm_start, t_comm_end, t_compute_start

    t_start = MPI_WTIME()
    cnt_blocketteRes = cnt_blocketteRes + 1

    ! ... [computePressure, computeLamVisc, computeEddyVisc, BCs] ...

    ! Marcar inicio da comunicação
    t_comm_start = MPI_WTIME()
    call whalo2(...)
    t_comm_end = MPI_WTIME()
    cnt_whalo2 = cnt_whalo2 + 1
    t_whalo2   = t_whalo2 + (t_comm_end - t_comm_start)
    t_blocketteRes_comm = t_blocketteRes_comm + (t_comm_end - t_comm_start)

    ! Marcar inicio de compute puro
    t_compute_start = MPI_WTIME()
    ! ... [loop sps, nn, blocketteResCore] ...
    t_blocketteRes_compute = t_blocketteRes_compute + (MPI_WTIME() - t_compute_start)

    t_blocketteRes = t_blocketteRes + (MPI_WTIME() - t_start)
end subroutine blocketteRes
```

#### Em `blocketteResCore()` — `src/NKSolver/blockette.F90:299`

```fortran
subroutine blocketteResCore(...)
    use adflowTimers
    real(kind=realType) :: t_start

    t_start = MPI_WTIME()
    cnt_blocketteResCore = cnt_blocketteResCore + 1

    ! ... [o loop OMP com tiles] ...

    t_blocketteResCore = t_blocketteResCore + (MPI_WTIME() - t_start)
end subroutine blocketteResCore
```

#### Em `ANKStep()` — `NKSolvers.F90:3629`

```fortran
subroutine ANKStep(firstCall)
    use adflowTimers
    real(kind=realType) :: t_start, t_pc_start, t_ksp_start

    cnt_ANKStep = cnt_ANKStep + 1
    t_start = MPI_WTIME()

    ! ...

    if (rebuilding_PC) then
        t_pc_start = MPI_WTIME()
        call FormJacobianANK()
        t_formJacobianANK = t_formJacobianANK + (MPI_WTIME() - t_pc_start)
        cnt_formJacobianANK = cnt_formJacobianANK + 1
        t_setupPC = t_setupPC + (MPI_WTIME() - t_pc_start)
    end if

    ! ...

    t_ksp_start = MPI_WTIME()
    call KSPSolve(ANK_KSP, rVec, deltaW, ierr)
    t_KSPSolve_ANK = t_KSPSolve_ANK + (MPI_WTIME() - t_ksp_start)
    cnt_KSPSolve_ANK = cnt_KSPSolve_ANK + 1
end subroutine ANKStep
```

#### Em `convergenceInfo()` — `src/solver/solvers.F90`

```fortran
! Envolver a chamada em solveState:
t_conv_start = MPI_WTIME()
call convergenceInfo()
t_convergenceInfo = t_convergenceInfo + (MPI_WTIME() - t_conv_start)
cnt_convergenceInfo = cnt_convergenceInfo + 1
```

#### Chamada final — no fim de `solveState()` ou em Python após `adflow.solvers.solver()`:

```fortran
call printTimerSummary(adflow_comm_world)
```

### 3.4 Saída esperada dos timers

Exemplo de output (fictício, para guiar a interpretação):

```
+=================================================================+
|              ADflow Timer Summary (wall, max rank)              |
+-------------------------+----------+-----------+-----------------+
| Section                 |  Total(s)| Calls     | Avg/call (ms)   |
+-------------------------+----------+-----------+-----------------+
| blocketteRes            |   142.30 |      4820 |          29.52  |
| blocketteResCore        |   118.40 |      4820 |          24.56  |
|   └ compute             |   118.40 |      4820 |          24.56  |
| whalo2                  |    23.90 |      4820 |           4.96  |
|   └ comm subtotal       |    23.90 |      4820 |           4.96  |
| FormJacobianANK (PC)    |    38.20 |        24 |        1591.67  |
| FormJacobianNK (PC)     |     9.10 |         6 |        1516.67  |
| KSPSolve ANK            |   112.00 |       210 |         533.33  |
| KSPSolve NK             |    18.50 |        14 |        1321.43  |
| convergenceInfo         |     2.10 |       224 |           9.375 |
| setupPC (total)         |    47.30 |        30 |        1576.67  |
| preSolverSetup          |     0.85 |         1 |         850.00  |
+=================================================================+
```

**Como interpretar os resultados:**

Olha para os números com estas perguntas:
1. Qual é o total acumulado maior? → Esse é o candidato a bottleneck primário.
2. Qual é o avg/call maior? → Se é o PC setup, o lag pode ajudar. Se é o blocketteResCore, é aritmética.
3. Qual a fracção comm/compute? → `t_whalo2 / t_blocketteResCore` diz se és bandwidth-bound de rede ou de memória.
4. Quantas vezes é chamado? → Um valor alto de `cnt_blocketteRes` confirma que é chamado muito.

—

### 3.5 Timer leve alternativo: por iteração ANK/NK

Se não quiseres modificar muito código, podes instrumentar apenas ao nível de
`solveState()` com um print por iteração:

```fortran
! No início do loop principal em solveState():
t_iter_start = MPI_WTIME()

! ... executeMGCycle / ANKStep / NKStep ...

! Depois de convergenceInfo():
if (myID == 0 .and. printIterations) then
    write(*,'(A, I6, A, F8.3, A)') &
        '  [TIMER] iter ', iterTot, ' wall_time=', MPI_WTIME()-t_iter_start, 's'
end if
```

Isto dá-te imediatamente a variação de custo por iteração à medida que o solver muda de fase.

---

## Parte 4: Profiling Sem Modificar o Código (Ferramentas Externas)

### 4.1 Linux `perf` — hotspot por função

**Requisito:** kernel com perf instalado, binário compilado com símbolos (`-g` ou debug).

**Uso básico — encontrar funções mais pesadas:**

```bash
# Compilar com símbolos (se necessário, adicionar -g a config.mk)
# Depois correr com perf:

mpirun -n 4 perf record -F 99 -g -- python run_adflow.py

# Depois de terminar:
perf report --stdio | head -60
```

**Output esperado (exemplo):**

```
# Overhead  Command       Shared Object  Symbol
  38.5%     python       adflow_mdolab  blocketteResCore_
  12.1%     python       adflow_mdolab  inviscidCentralFlux_
   8.3%     python       adflow_mdolab  allNodalGradients_
   6.9%     python       adflow_mdolab  whalo2_
   5.4%     python       adflow_mdolab  viscousFlux_
   ...
```

**Interpretação guiada:** Regista as percentagens de cada função. A ordenação por overhead vai revelar as funções mais pesadas sem qualquer assumção prévia.

Para profiling por thread com OpenMP:

```bash
perf record -F 99 -g --call-graph=dwarf -- python run_adflow.py
perf report --stdio --no-children | head -80
```

---

### 4.2 `perf stat` — contadores de hardware (cache misses, FLOPS)

```bash
# Contadores básicos de cache e branch:
mpirun -n 1 perf stat -e \
  cache-references,cache-misses,\
  L1-dcache-loads,L1-dcache-load-misses,\
  LLC-loads,LLC-load-misses,\
  instructions,cycles \
  -- python run_adflow.py 2>&1 | tee perf_stat_output.txt
```

**O que procurar:**

| Métrica | Bom | Preocupante |
|---|---|---|
| `L1-dcache-load-misses / L1-dcache-loads` | < 5% | > 15% |
| `LLC-load-misses / LLC-loads` | < 5% | > 20% |
| `cache-misses / cache-references` | < 5% | > 30% |
| `instructions / cycles (IPC)` | > 2.0 | < 1.0 |

**Para isolar só o blocketteResCore:**
Podes executar numa corrida muito curta (1-2 iterações ANK) para reduzir o ruído.

---

### 4.3 `perf stat` — bandwidth de memória (eventos DRAM)

```bash
# Eventos específicos para Intel (nomes podem variar por CPU):
perf stat -e \
  uncore_imc/data_reads/,\
  uncore_imc/data_writes/ \
  -- python run_adflow_short.py 2>&1 | tee mem_bw.txt
```

Alternativamente, usar `likwid-perfctr` se disponível:

```bash
# Grupo de bandwidth:
likwid-perfctr -C 0-7 -g MEM_DP -- python run_adflow_short.py

# Grupo de flops:
likwid-perfctr -C 0-7 -g DP_SP -- python run_adflow_short.py
```

**O que calcular:** `bytes_transferred / total_flops` → Arithmetic Intensity (AI).
Se AI < roofline AI do teu processador, o código é **memory-bound**.

---

### 4.4 `valgrind cachegrind` — análise detalhada de cache por instrução

```bash
# ATENÇÃO: corre 10-100x mais lento. Usa malha muito pequena e poucas iterações.
mpirun -n 1 valgrind --tool=cachegrind \
  --I1=32768,8,64 --D1=32768,8,64 --LL=6291456,12,64 \
  python run_adflow_tiny.py

# Analisa resultados:
cg_annotate cachegrind.out.<pid> | head -100
```

**Ajusta os parâmetros I1/D1/LL** para corresponder à cache do teu CPU:

```bash
# Ver tamanhos de cache:
lscpu | grep -i cache
getconf LEVEL1_DCACHE_SIZE  # bytes
getconf LEVEL2_CACHE_SIZE
getconf LEVEL3_CACHE_SIZE
```

---

### 4.5 `valgrind callgrind` — tempo por função (simulação)

```bash
mpirun -n 1 valgrind --tool=callgrind \
  --callgrind-out-file=callgrind.out \
  python run_adflow_tiny.py

# Visualizar com kcachegrind (GUI) ou:
callgrind_annotate callgrind.out | head -60
```

---

### 4.6 Script `run_adflow_short.py` para profiling

Criar um script mínimo para profiling sem gastar tempo desnecessário:

```python
from mpi4py import MPI
from adflow import AERO_SOLVER

options = {
    'gridFile': 'input_files/your_small_mesh.cgns',
    'equationMode': 'steady',
    'equations': 'RANS',
    'turbulenceModel': 'SA',
    'useBlockettes': True,               # garantir blockettes activos
    'ANKUseMatrixFree': True,
    'nCycles': 10,                       # poucas iterações para profiling
    'printTiming': True,
    # Para profiling de cache, prefere 1 rank
}

CFDSolver = AERO_SOLVER(options=options, comm=MPI.COMM_WORLD)
# ... definir aeroProblem ...
CFDSolver(ap)
```

---

## Parte 5: Diagnóstico de Cache no Blockette

### 5.1 Estimativa do working set por thread

Para um tile `BS=8`, o que cada thread precisa carregar antes de calcular:

```
--- Arrays de estado ---
w:       9 × 9 × 9 × 7 vars × 8 bytes = 32,256 bytes  ≈ 31.5 KB
p:       9 × 9 × 9 × 8 bytes           =  4,374 bytes  ≈  4.3 KB
gamma:   9 × 9 × 9 × 8 bytes           =  4,374 bytes  ≈  4.3 KB
ss:      9 × 9 × 9 × 8 bytes           =  4,374 bytes  ≈  4.3 KB

--- Arrays de propriedades ---
rlv/rev: 2 × 10×10×10 × 8 bytes        = 16,000 bytes  ≈ 15.6 KB
vol:     10×10×10 × 8 bytes            =  8,000 bytes  ≈  7.8 KB
x:       10×10×10 × 3 × 8 bytes        = 24,000 bytes  ≈ 23.4 KB

--- Gradientes nodais (if viscous) ---
ux..qz:  12 × 9×9×9 × 8 bytes         = 55,296 bytes  ≈ 54 KB

--- Outros (sI/J/K, dss, aa, dw, fw, etc.) ---
estimativa:                                            ≈ 30 KB

TOTAL APROXIMADO POR THREAD:                          ≈ 175 KB
```

**Comparação com cache sizes típicos:**
- L1 data cache: ~32 KB por core → **working set não cabe em L1**
- L2 cache: ~256–512 KB por core → **working set cabe em L2 (marginal)**
- L3 (shared): vários MB → **sem pressão de L3 se poucos threads**

**Conclusão a testar:** O tile cabe em L2 se o L2 for ≥ 256 KB. Com múltiplos threads a competir pelo L2/L3, pode haver pressure. O `perf stat` vai confirmar.

**Com BS=8 e viscoso+SA:** O working set é dominado pelos 12 arrays de gradientes nodais (~54 KB). Ajustar BS afecta directamente isto.

---

### 5.2 Análise dos padrões de cópia no blockette

O blockette faz múltiplos passes sobre os mesmos dados. Contar os loads/stores:

**Por tile, com `updateIntermed=False` (MFFD matvec, caso mais comum):**

| Operação | Direction | Arrays | Bytes estimados |
|---|---|---|---|
| Copy-in w,p,γ,ss | global→privado | w(full)+3 escalares | ~45 KB |
| Copy-in rlv,rev,vol | global→privado | 3 × 10³ | ~24 KB |
| Copy-in x | global→privado | coordenadas nó | ~24 KB |
| Copy-in iblank,d2wall,volRef | global→privado | 3 × 8³ | ~12 KB |
| Copy-in porI/J/K | global→privado | 3 × faces | ~12 KB |
| `metrics` calc | read x, write sI/J/K | - | - |
| `initRes` | clear dw | - | - |
| SA: saSource,Adv,Visc | read w,rlv,rev, write dw | ~40 KB | |
| `timeStep` | read vol,sI/J/K, write radI/J/K,dtl | ~40 KB | |
| `inviscidCentralFlux` | read w,p,sI/J/K, write fw | ~50 KB | |
| `allNodalGradients` | read w, write ux..qz | ~90 KB | |
| `viscousFlux` | read ux..qz,sI/J/K,rlv, write fw | ~80 KB | |
| `sumDwAndFw` | read fw, write dw | ~20 KB | |
| Copy-out dw | privado→global | ~16 KB | |
| **TOTAL reads/writes estimados** | | | **~452 KB por tile** |

Para que este seja eficiente e não memory-bound desde a RAM:
- Os dados precisam de estar em L2/L3 durante toda a sequência de kernels num tile.
- Se a cache for suficientemente grande (≥ 256 KB L2 por core), a reutilização dentro do tile é boa.
- Se não couber, cada kernel faz um round-trip à RAM.

---

### 5.3 Perguntas para guiar o diagnóstico

Depois de correr `perf stat`, responde a estas perguntas:

1. **Qual a LLC-load-miss rate?**
   - < 5%: working set cabe bem em L3, o código é compute-bound localmente
   - 5–20%: pressão de L3, bandwidth de memória começa a importar
   - > 20%: claramente memory-bandwidth bound

2. **Qual o IPC (instructions per cycle)?**
   - > 2: bom uso de pipeline
   - 1–2: pipeline moderadamente ocupado
   - < 1: stalls dominam (tipicamente memory stalls)

3. **Qual a fracção de tempo em `whalo2` vs `blocketteResCore`?**
   - Se `t_whalo2 / t_blocketteResCore > 0.3`: comunicação MPI é significativa
   - Se `t_whalo2 / t_blocketteResCore < 0.1`: compute domina, não a rede

4. **Qual o breakup dentro de `blocketteResCore`?**
   - Para descobrir, podes adicionar timers finos dentro do loop de tiles:
   ```fortran
   t_metrics_start = MPI_WTIME()
   call metrics
   t_metrics = t_metrics + (MPI_WTIME() - t_metrics_start)
   ! ... idem para cada kernel ...
   ```
   - Os mais pesados num caso RANS-SA viscoso são tipicamente:
     `allNodalGradients`, `viscousFlux`, `inviscidCentralFlux`, `saSource`

5. **Como varia com BS?**
   - Testa BS=4, 6, 8, 10, 12 (modificar `src/NKSolver/blockette.F90:9`)
   - Recompila e corre 10 iterações ANK
   - Plota `t_blocketteResCore` vs BS → o mínimo indica o ponto de melhor uso de cache

---

## Parte 6: Estratégia de Optimização do Blockette (Baseada nos Resultados)

**Esta secção deve ser lida DEPOIS de teres os dados dos timers e do `perf stat`.  
Cada acção abaixo é justificada por um resultado específico de diagnóstico.**

### 6.1 Se: LLC miss rate > 15% → Optimização de cache

**Causa provável:** Working set excede L2, dados são buscados de L3 ou RAM repetidamente.

**Acção A — Reduzir BS:**
```fortran
! src/NKSolver/blockette.F90:9
integer(kind=intType), parameter :: BS = 6  ! (ou 4)
integer(kind=intType), parameter :: bbil = BS + 1, ...
```
Tiles menores → working set menor → maior probabilidade de ficar em L2.
**Trade-off:** tiles menores reduzem reutilização no eixo i → mais copy-in overhead.

**Acção B — Reduzir copy-out condicional:**
O copy-out de `updateIntermed=True` (gradientes nodais, spectral radii) só é necessário
após o passo final. Para as chamadas MFFD intermédias (`updateIntermed=False`) já é
omitido. Verificar se não há caminhos de código que passam `updateIntermed=True`
desnecessariamente.

**Acção C — Fuser loops dentro do tile:**
`metrics` e `inviscidCentralFlux` ambos leem `x` e `sI/sJ/sK`. Actualmente são rotinas
separadas. Fundir os seus loops elimina uma releitura de `sI/J/K`.

### 6.2 Se: IPC < 1.0 → Stalls por dependência de dados ou branches

**Acção A — Verificar vectorização:**
```bash
# Compilar com report de vectorização (Intel):
-qopt-report=5 -qopt-report-phase=vec

# ou gfortran:
-fopt-info-vec-optimized -fopt-info-vec-missed
```
Identificar loops dentro de `inviscidCentralFlux` / `viscousFlux` que não vectorizam.
Causas comuns: halos com índice 0, arrays de dimensão `(0:bbil)` que inibem vectorização.

**Acção B — Remover condicionais de runtime dentro do loop de tiles:**
```fortran
! ANTES: dentro do loop de tiles
if (equations == RANSEquations .and. turbRes) then
    call saSource
    ...
end if
! DEPOIS: dois paths especializados, selecionados fora do loop
```

### 6.3 Se: `t_whalo2 / t_total_blocketteRes > 30%` → Comunicação domina

**Significa:** MPI halo exchange é o bottleneck, não a aritmética.

**Acção A — Verificar sobreposição compute/comm:**
O padrão actual é: compute pressão/visc → halo → compute kern.
Poderia ser: lança halo não-blocking → compute pressão/visc (overlap) → wait halo → compute kern.

```fortran
! Lançar recv/send antes de computePressureSimple:
call MPI_Irecv(...) ; call MPI_Isend(...)   ! lança comunicação
call computePressureSimple(...)             ! overlaps com comm
call MPI_Waitall(...)                       ! aguarda
! ... continua com blocketteResCore ...
```

**Acção B — Reduzir frequência de halo:**
Para o MFFD, cada `blocketteRes` faz 1 `whalo2`. Num passo ANK com N_ksp=15:
15 halos por step não-linear. Este é o custo estrutural do esquema matrix-free.
Reduzir N_ksp (PC melhor) reduz directamente halos.

**Acção C — Menos ranks, mais threads:**
Menos ranks → menos mensagens MPI → menos latência de halo.
Mais threads por rank → `blocketteResCore` paralelo dentro do rank.
O OpenMP loop `collapse(2)` sobre `(kk, jj)` já distribui tiles por threads.
Com BS=8 e bloco de 100×100×100: tiles = 12×12×12 ≈ 1728 tiles por bloco.
Com 8 threads: ~216 tiles por thread — paralelismo suficiente.

### 6.4 Se: `t_formJacobianANK` é grande → PC setup domina

**Acção:** Aumentar `ANKJacobianLag` (ex: de 20 para 40 ou 60).
Risco: mais iterações KSP por step se PC ficar desactualizado.
Diagnóstico: comparar `N_ksp` médio com e sem lag mais alto.

**Acção alternativa:** Activar `ANKADPC=True` — usa AD forward mode para PC (mais exacto
por menos avaliações, mas custo de AD por avaliação é ~2–4x o FD).

### 6.5 Se: `allNodalGradients` é > 20% do tempo do tile → Viscoso caro

**O que acontece:** Para cada tile, `allNodalGradients` computa 12 campos (ux,uy,...,qz)
por Gauss com stencil nodal. É a operação mais cara em fluxo viscoso de alta ordem.

**Acção A:** Verificar se `useViscApprox=True` é aceitável para as iterações ANK:
- Isso bypassa `allNodalGradients` e usa `viscousFluxApprox`
- Por defeito `ANKUseViscApprox=True` — confirmar que está activo na configuração

**Acção B:** Verificar que o storeWall flag não está a forçar operações extra
em `viscousFlux` desnecessariamente no path MFFD.

### 6.6 Adicionar OpenMP a kernels fora do blockette (longo prazo)

Candidatos identificados sem OpenMP actual:

| Ficheiro | Rotina | Cuidados |
|---|---|---|
| `src/solver/flowUtils.F90:867` | `computePressureSimple` | Seguro: loop sobre células sem dependência |
| `src/solver/flowUtils.F90:1201` | `computeLamViscosity` | Seguro: mesmo padrão |
| `src/turbulence/turbUtils.F90:580` | `computeEddyViscosity` | Seguro se SA: cálculo local por célula |
| `src/adjoint/adjointUtils.F90` | `setupStateResidualMatrix` (loop over colors) | Cuidado: cada cor modifica w globalmente |

Para `computePressureSimple` e `computeLamViscosity`, o padrão de adição é:

```fortran
! ANTES (src/solver/flowUtils.F90):
do sps = 1, nTimeIntervalsSpectral
    do nn = 1, nDom
        call setPointers(nn, 1, sps)
        do k = 2, kl
            do j = 2, jl
                do i = 2, il
                    p(i,j,k) = (gamma(i,j,k) - 1) * &
                        (w(i,j,k,irhoE) - 0.5*(...))
                end do
            end do
        end do
    end do
end do

! DEPOIS:
do sps = 1, nTimeIntervalsSpectral
    do nn = 1, nDom
        call setPointers(nn, 1, sps)
        !$OMP parallel do collapse(2) private(i,j,k)
        do k = 2, kl
            do j = 2, jl
                do i = 2, il
                    p(i,j,k) = (gamma(i,j,k) - 1) * &
                        (w(i,j,k,irhoE) - 0.5*(...))
                end do
            end do
        end do
        !$OMP end parallel do
    end do
end do
```

**IMPORTANTE:** Não adicionar OpenMP sem verificar:
1. Ausência de escrita em arrays globais não-privados dentro do loop
2. Ausência de chamadas a MPI dentro de regiões paralelas
3. Teste de convergência antes e depois para confirmar idêntico comportamento

---

## Parte 7: Procedimento de Diagnóstico Passo a Passo

### Passo 1: Instrumentar timers mínimos (1 dia)

1. Criar `src/utils/adflowTimers.F90` com o módulo acima.
2. Adicionar timers em: `blocketteRes`, `blocketteResCore`, `whalo2`, `FormJacobianANK`, `FormJacobianNK`, `KSPSolve` (ANK e NK), `convergenceInfo`.
3. Adicionar chamada a `printTimerSummary` no fim de `solveState`.
4. Recompilar e correr num caso de referência (100 iterações ANK + 20 NK).

### Passo 2: Interpretar os timers (0.5 dia)

- Preencher uma tabela com os tempos medidos.
- Calcular `% do total` para cada secção.
- Calcular `t_whalo2 / t_blocketteRes`.
- Anotar as chamadas mais frequentes e com maior tempo per-call.

### Passo 3: Profiling de hardware (0.5 dia)

```bash
# Corre com perf stat na configuração de 1 rank:
mpirun -n 1 perf stat -e cache-misses,LLC-load-misses,instructions,cycles \
    python run_adflow_10iter.py
```

- Anotar LLC miss rate e IPC.
- Comparar com a estimativa teórica da Parte 5.1.

### Passo 4: Test de sensibilidade BS (0.5 dia)

Para BS in {4, 6, 8, 10, 12}:
1. Editar `blockette.F90:9`, recompilar, correr 10 iter ANK.
2. Registar `t_blocketteResCore` e LLC miss rate.
3. Plotar curva de `t_blocketteResCore` vs BS.

### Passo 5: Tomar decisões de optimização

Com os dados dos Passos 1-4, és capaz de:
- Saber se o bottleneck é compute vs. comm (do Passo 2).
- Saber se é cache ou bandwidth (do Passo 3).
- Saber o BS óptimo (do Passo 4).
- Priorizar acções de optimização da Parte 6 com base em evidências.

---

---

## Parte 8: Decomposição Temporal Exaustiva — O Que Demora em ANKStep e NKStep

Esta parte responde à pergunta central: **dentro de um passo ANK ou NK, qual é a
fracção de tempo gasta em residuais, na álgebra linear do GMRES, no precondicionador,
nas comunicações e em tempo idle?**  
E ainda: **o que escala bem com mais threads OpenMP, e o que escala (mal) com mais
ranks MPI?**

---

### 8.1 As Cinco Categorias de Tempo

Todo o tempo dentro de `ANKStep` e `NKStep` pode ser classificado em cinco categorias:

| # | Categoria | O que inclui | MPI? | OpenMP? |
|---|---|---|---|---|
| **A** | **Avaliação de residual** | `blocketteRes` completo, incluindo `computePressureSimple`, `computeLamViscosity`, `computeEddyViscosity`, BCs, `whalo2`, `blocketteResCore` | halo MPI (vizinhos) | sim (`collapse(2)` sobre tiles) |
| **B** | **Álgebra linear GMRES** (excluindo MatMult) | ortogonalização, normas, rotações de Givens, `VecAXPY`, `VecCopy`, update espaço de Krylov | Allreduce global | não (PETSc serial aqui) |
| **C** | **Precondicionador — apply** | ILU forward/backward substitution, AMG V-cycle (smooth, restrict, prolong) | halo local (vizinhos, ASM overlap); Allreduce nos níveis grossos do AMG | não (PETSc ILU serial) |
| **D** | **Precondicionador — setup** | `FormJacobianANK`/`FormJacobianNK`: FD coloring (N blocketteRes), `MatAssembly`, `MatAXPY`, `setupAMG`/`KSPSetUp` | halo MPI (coloring) + MatAssembly (global) | não (loop de coloring é serial) |
| **E** | **Tempo idle / imbalance** | tempo que cada rank passa à espera nas barreiras implícitas de `MPI_Allreduce` e `whalo2` enquanto ranks mais lentos completam computação | é o custo não-overlap de Allreduce | não aplicável |

---

### 8.2 Linha do Tempo Detalhada de `ANKStep()` com Custos

Ficheiro: `src/NKSolver/NKSolvers.F90:3629`

```
ANKStep(firstCall)
│
├─ [SETUP — apenas no firstCall ou mudança coupled/uncoupled]
│    ├─ destroyANKSolver() + setupANKSolver()    [categoria D]
│    │    ├─ MatCreateMFFD, VecCreate, MatCreate  — trivial (1x)
│    │    └─ setupAMG() se precondType='mg'       — custo moderado (1x)
│    ├─ setwVecANK(wVec, ...)                    [categoria B — VecCopy local]
│    └─ blocketteRes(useUpdateIntermed=True)      [categoria A — 1 residual]
│
├─ [computeTimeStepMat(usePC=True/False)]         [categoria D parcial]
│    ├─ loop cells: computeTimeStepBlock()        — aritmética local
│    ├─ MatSetValuesBlocked()                     — preenchimento local
│    └─ MatAssemblyBegin/End(timeStepMat)         — comunicação global (MatAssembly)
│         → 1 MPI_Allreduce implícito em MatAssembly
│
├─ [FormJacobianANK() — se lag atingido]          [categoria D — CARO]
│    ├─ setupStateResidualMatrix(...)
│    │    ├─ setFDReference()                     → 1 blocketteRes [categoria A]
│    │    └─ loop: nColors × nStateVars:
│    │         → blocketteRes() por cor           → ~20-50 blocketteRes [categoria A]
│    │         (cada um com whalo2 interno)
│    ├─ MatAssemblyBegin/End(dRdwPre)             → MPI global
│    ├─ MatAXPY(dRdwPre, one, timeStepMat)        → local (sparsity)
│    └─ setupAMG() ou KSPSetUp (ASM+ILU)          → ILU factorization local
│         ├─ se ANK_precondType='asm': ILU factor → sem MPI, custo ∝ nnz_local
│         └─ se ANK_precondType='mg':  setupAMG   → construção dos níveis grossos
│              └─ MatPtAP() (Galerkin coarse): MPI local-global
│
├─ [MatAssemblyBegin/End(dRdw) — dummy para MFFD] [categoria B — trivial]
│
├─ [formFunction_mf() — base MFFD]                [categoria A — 1 blocketteRes]
│    └─ blocketteRes(w_atual)                     → guarda R(w) como referência MFFD
│
├─ [KSPSetTolerances(), KSPSetResidualHistory()]   [categoria B — trivial]
│
├─ [KSPSolve(ANK_KSP, rVec, deltaW)]              [ver §8.4 para detalhe completo]
│    └─ GMRES com m iterações (m ≤ ANK_maxIter)
│         por iteração j = 1..m:
│         ├─ PC apply  (right): K^{-1} v          [categoria C]
│         ├─ MatMult (MFFD): R(w+εv')−R(w))/ε    [categoria A — 1 blocketteRes]
│         ├─ VecMDot(new_v, [q1..qj])             [categoria B — 1 Allreduce de tamanho j]
│         ├─ VecAXPY × j (ortogonalização)        [categoria B — local]
│         ├─ VecNorm (normalização)               [categoria B — 1 Allreduce de tamanho 1]
│         └─ Givens rotation (least-squares)      [categoria B — local, negligível]
│
├─ [physicalityCheckANK()]                         [categoria B — local]
│
├─ [VecAXPY(wVec, -lambda, deltaW)]               [categoria B — local]
│
├─ [setWANK(wVec, ...)]                           [categoria B — cópia local]
│
├─ [computeUnsteadyResANK(lambda)]                [categoria A — 1 blocketteRes]
│    └─ chama blocketteRes(useUpdateIntermed=True) internamente
│
├─ [Backtrack line search — se necessário]         [categoria A — até 12 blocketteRes]
│    └─ loop até 12: VecAXPY + setWANK + computeUnsteadyResANK
│
├─ [blocketteRes(useUpdateIntermed=True) — FINAL]  [categoria A — 1 blocketteRes]
│    └─ actualiza intermediários (dtl, radii, gradientes) em memória global
│
└─ [setRVec(rVec), VecNorm(rVec), actualizar contadores]  [categoria B — 1-2 Allreduce]
```

**Contagem total de blocketteRes por ANKStep (categoria A + dentro de D):**

| Contexto | Nº blocketteRes | Custo per blocketteRes | MPI |
|---|---|---|---|
| Initial res + base MFFD + unsteady res + final res | 4 | ~25 ms (RANS) | 1 whalo2 cada |
| FormJacobianANK (FD coloring, se rebuild) | ~20–50 | ~25 ms | 1 whalo2 cada |
| KSPSolve: 1 MFFD matvec por iteração GMRES (m iters) | m (~5–30 típico) | ~25 ms | 1 whalo2 cada |
| Line search backtrack (se necessário) | 0–12 | ~25 ms | 1 whalo2 cada |
| **Total sem rebuild (~m=15)** | **~19** | | **~19 whalo2** |
| **Total com rebuild (~m=15)** | **~55–65** | | **~55–65 whalo2** |

---

### 8.3 Linha do Tempo Detalhada de `NKStep()` com Custos

Ficheiro: `src/NKSolver/NKSolvers.F90:512`

```
NKStep(firstCall)
│
├─ [SETUP — apenas no firstCall]                  [categoria D — 1x]
│    └─ setupNKSolver(): MatCreateMFFD, VecCreate, MatCreate
│
├─ [setwVec(wVec)]                                [categoria B — local]
├─ [computeResidualNK() — apenas firstCall]        [categoria A — 1 blocketteRes]
├─ [setRVec(rVec)]                                [categoria B — local]
│
├─ [VecNorm(rVec, NORM_2, norm)]                  [categoria B — 1 Allreduce]
│
├─ [FormJacobianNK() — se mod(NK_iter, lag)==0]   [categoria D — CARO]
│    └─ setupStateResidualMatrix(...)              → ~20–50 blocketteRes [categoria A]
│         └─ KSPSetUp → ILU factor ou AMG setup [categoria D]
│
├─ [MatAssemblyBegin/End(dRdw) — se sem rebuild]  [categoria B — trivial]
│
├─ [getEWTol() — Eisenstat-Walker]               [categoria B — local]
├─ [KSPSetTolerances()]                           [categoria B — trivial]
│
├─ [formFunction_mf() — base MFFD]                [categoria A — 1 blocketteRes]
├─ [MatMFFDSetBase()]                             [categoria B — trivial]
│
├─ [KSPSolve(NK_KSP, rVec, deltaW)]               [ver §8.4]
│    └─ GMRES com m iterações (m ≤ NK_subspace=60)
│         — mesmo breakdown por iteração que §8.2
│         — mas com subspace até 60 → custo ortogonalização 2× maior que ANK
│
├─ [Line Search: LSNone / LSCubic / LSNM]          [categoria A — 1-10 blocketteRes]
│    └─ LSCubic (default):
│         ├─ MatMult(dRdw, y, w)                  → 1 MFFD = 1 blocketteRes
│         ├─ VecDot(f, w, initslope)              → 1 Allreduce [categoria B]
│         ├─ 1º eval: setW + computeResidualNK    → 1 blocketteRes
│         └─ loop cubic até convergência:
│              → 1 blocketteRes por iteração (típico: 2–4 total)
│
└─ [VecCopy(work, wVec), VecCopy(g, rVec)]        [categoria B — local]
```

**Contagem total de blocketteRes por NKStep:**

| Contexto | Nº blocketteRes | Custo per blocketteRes | MPI |
|---|---|---|---|
| Residual inicial (firstCall) | 1 | ~25 ms | 1 whalo2 |
| FormJacobianNK (se rebuild) | ~20–50 | ~25 ms | 1 whalo2 cada |
| Base MFFD | 1 | ~25 ms | 1 whalo2 |
| KSPSolve: MFFD matvec por iter GMRES (m iters) | m (~10–60) | ~25 ms | 1 whalo2 cada |
| Line search LSCubic (típico 3 evals) | 3–5 | ~25 ms | 1 whalo2 cada |
| **Total sem rebuild (m=20, LS=3)** | **~25** | | **~25 whalo2** |
| **Total com rebuild (m=20, LS=3)** | **~55–80** | | **~55–80 whalo2** |

---

### 8.4 Dentro do `KSPSolve` — GMRES Iteração a Iteração

Este é o nível mais fino de detalhe dentro do solver linear.  
Para o GMRES right-preconditioned (configuração default em ADflow), cada iteração `j`  
faz exactamente esta sequência dentro de `KSPSolve`:

```
Iteração GMRES j:
──────────────────────────────────────────────────────────────────
STEP 1: PC apply (right preconditioning): z_j = K^{-1} p_j
         → Ver §8.5 para detalhe interno do PC apply

STEP 2: MatMult (MFFD): w_j = A z_j = [R(w+ε z_j) − R_base] / ε
         → 1 chamada a blocketteRes (inclui whalo2, blocketteResCore)
         → custo dominante de cada iteração GMRES

STEP 3: Ortogonalização (Classical Gram-Schmidt, sem refinamento)
         → VecMDot(w_j, [q_1, ..., q_j, w_j])  → 1 MPI_Allreduce de tamanho (j+1)
            [PETSc junta todos os j inner products num único Allreduce]
         → VecAXPY × j: w_j ← w_j − Σ H_{ij} q_i              (local, sem MPI)
         → Custo MPI: 1 Allreduce com j+1 doubles por iter j
            → total ao longo de m iters: m Allreduce, tamanhos crescentes 2,3,...,m+1
            → latência total: ~m × latência_Allreduce
            → volume total: 8 bytes × (2+3+...+m+1) = 8 × m(m+3)/2

STEP 4: VecNorm(w_j): normalização do novo vector de Krylov
         → 1 MPI_Allreduce de tamanho 1 (um double)
         → custo: ~1 latência de Allreduce

STEP 5: Givens rotation (Least-squares update)
         → aritmética local 2×2, completamente negligível
──────────────────────────────────────────────────────────────────
```

**Custo acumulado de MPI do GMRES ao longo de m iterações:**

| Fonte | Nº Allreduce | Volume (doubles) | Latência dominante? |
|---|---|---|---|
| Ortogonalização (VecMDot) | m calls | m(m+3)/2 | sim, cresce com m |
| VecNorm por iteração | m calls | m | sim |
| VecNorm final (convergência) | ~1–few | ~1 | sim |
| **TOTAL** | **~2m calls** | **~m²/2 doubles** | **O(m × lat)** |

Para NK com m=60: **~120 Allreduce** só para a álgebra GMRES, antes de contar halos.  
Para ANK com m=200: **~400 Allreduce** por KSPSolve → pode dominar em clusters com muitos ranks.

**Custo de MatMult (MFFD) vs. ortogonalização — exemplo:**

| Parâmetro | Valor típico |
|---|---|
| T_blocketteRes (1 residual) | ~25 ms (malha 2M cells, 8 OMP threads) |
| T_Allreduce (latência, LAN 1 Gbit) | ~0.1–0.5 ms por call |
| T_Allreduce (latência, InfiniBand) | ~0.01–0.05 ms por call |
| MatMult por iter (m=15, ANK) | 15 × 25ms = 375 ms |
| Ortogonalização (m=15 iters) | 30 Allreduce × 0.1ms = 3 ms (IB) ou 15 ms (LAN) |
| PC apply (m=15 iters ASM+ILU) | 15 × T_ilu (ver §8.5) |

**Conclusão imediata:** Para ALT malhas e poucas iterações GMRES, o MatMult domina.  
Para muitas iterações GMRES (NK com m→60) e muitos ranks (P→512+),  
a ortogonalização pode tornar-se comparável ao MatMult.

---

### 8.5 Precondicionador — Custo de Apply: ASM+ILU vs. AMG

#### 8.5.1 ASM+ILU (default: `ANK_precondType = 'asm'`)

O PC apply de ASM+ILU por cada iteração GMRES faz:

```
PCApply (ASM+ILU):
├─ [COMUNICAÇÃO ASM — troca de overlap com ranks vizinhos]
│    ├─ MPI_Isend / MPI_Irecv (não-blocking) para células de overlap
│    │    → volume: overlap × nCell_face × nw × 8 bytes
│    │    → comunicação ponto-a-ponto (não Allreduce) com vizinhos imediatos
│    ├─ MPI_Waitall
│    └─ cost: ~latência + volume_overlap / bandwidth_rede
│        típico com ANK_asmOverlap=1: ~overlap_cells × 7vars × 8bytes por face MPI
│
└─ [ILU SOLVE — factorização local aplicada]
     ├─ Forward substitution: L y = r  (serial, sem MPI, sem OpenMP)
     ├─ Backward substitution: U z = y (serial, sem MPI, sem OpenMP)
     └─ custo ∝ nnz(dRdwPre) × (1 + fill) × n_innerPreConIts
         → escalona com tamanho DO PROBLEMA LOCAL (decresce com mais ranks)
         → NÃO beneficia de mais OpenMP threads (PETSc ILU não usa OMP)
         → ILU quality DEGRADA com mais MPI ranks (interface boundary conditions)
```

**Impacto do fill (`ANKILUFill`):**  
- Fill 0 (ILU(0)): 1 forward/back pass, rápido mas PC fraco → mais iters GMRES  
- Fill 2: mais operações, PC melhor → menos iters GMRES  
- Trade-off: aumentar fill reduz iters mas aumenta custo por apply

#### 8.5.2 AMG (algebraic multigrid custom: `ANK_precondType = 'mg'`, `NK_precondType = 'mg'`)

O PC apply do AMG (V-cycle com `ANK_AMGLevels` níveis) faz:

```
PCApply (AMG V-cycle):
│
├─ [PRÉ-SMOOTH no nível 1 (fino)]: ILU apply local
│    → serial, sem MPI, custo ∝ nnz_fine × fill
│
├─ [RESTRIÇÃO: r_coarse = R × r_fine]
│    → MatMult(R, r) — local (coarseIndices são locais por construção)
│    → sem MPI neste passo
│
├─ [COARSE SOLVE no nível 2..N]     ← recurso (V-cycle recursivo ou K-cycle)
│    ├─ se nível intermédio: repete smooth + restrict + ...
│    └─ se nível mais grosso: ILU direct solve
│         → se problema coarse ainda distribuído: MatMult no coarse requer MPI
│         → comunicação inter-rank nos níveis grossos (coarse grid P2P, não Allreduce)
│
├─ [PROLONGAÇÃO: e_fine = P × e_coarse]
│    → MatMult(P, e) — pode requerer halo nos níveis grossos
│    → custo decresce com nível (malha mais grossa)
│
└─ [PÓS-SMOOTH no nível 1 (fino)]: ILU apply local
     → igual ao pré-smooth
```

**AMG vs. ASM+ILU — trade-offs:**

| Critério | ASM+ILU | AMG |
|---|---|---|
| Custo de setup (FormJacobian) | menor | maior (construção de P, R, A_coarse) |
| Custo de apply por iter GMRES | baixo (1 ILU pass local) | moderado (V-cycle multi-nível) |
| Qualidade do PC | boa para problemas bem condicionados | melhor para problemas stiff |
| Scaling com mais ranks | PC degrada (mais interfaces) | AMG mantém qualidade via coarse grid |
| MPI no apply | apenas overlap vizinhos | também coarse grid (mais comunicação) |
| OpenMP no apply | não (ILU serial) | não (ILU serial nos smoothers) |
| Iters GMRES necessárias | mais (PC mais fraco) | menos (PC mais robusto) |

---

### 8.6 Tempo Idle — Imbalance MPI

O **tempo idle** é o tempo que cada rank passa à espera dentro de operações de sincronização  
enquanto outros ranks ainda estão a fazer computação útil.

#### 8.6.1 Fontes de imbalance em ADflow

```
Fontes de idle time:
│
├─ [whalo2 — halo exchange]
│    ├─ Ranks com blocos maiores (mais células por bloco) levam mais tempo em blocketteResCore
│    ├─ Ranks com blocos mais pequenos terminam cedo e ficam à espera no MPI_Waitall
│    └─ Ocorre em CADA blocketteRes — acumula ao longo das ~50-80 chamadas por passo NK/ANK
│
├─ [MPI_Allreduce — normas GMRES e convergenceInfo]
│    ├─ Todos os ranks chegam ao Allreduce em momentos diferentes
│    ├─ Os ranks rápidos esperam pelos lentos (barrier implícito na maioria das implementações)
│    ├─ Quantas: ~2m Allreduce por KSPSolve + 1-3 por convergenceInfo
│    └─ O idle aqui é PROPORCIONAL ao desvio de tempo entre o rank mais rápido e o mais lento
│
├─ [MatAssemblyEnd — assembly de matrizes PETSc]
│    ├─ PETSc MatAssembly envolve comunicação global (compactação de off-diagonal entries)
│    ├─ Acontece em FormJacobianANK/NK e computeTimeStepMat
│    └─ Idle se ranks terminam a fase MatSetValues em tempos diferentes
│
└─ [convergenceInfo() — MPI_Allreduce + I/O]
     ├─ Chamado uma vez por iteração não-linear — sincroniza todos os ranks
     ├─ Faz Allreduce para normas L2 globais
     └─ Rank 0 imprime — outros esperam (trivial mas existente)
```

#### 8.6.2 Como medir o tempo idle

**Método A — Barrier antes do Allreduce:**

```fortran
! No whalo2 ou antes do VecNorm:
t_barrier_start = MPI_WTIME()
call MPI_Barrier(adflow_comm_world, ierr)   ! barreira explícita
t_idle_barrier = t_idle_barrier + (MPI_WTIME() - t_barrier_start)
! depois o Allreduce/whalo2 real
```

O tempo na barreira é exactamente o tempo que este rank esperou pelos outros.

**Método B — `mpiP` profiler (externo, sem modificar código):**

```bash
# Instalar mpiP (http://mpip.sourceforge.net)
# Compilar ADflow com -lmpiP linkado
mpirun -n 32 python run_adflow.py

# mpiP gera um ficheiro .mpiP com:
# - % tempo em MPI por função (MPI_Allreduce, MPI_Isend, MPI_Waitall, ...)
# - breakdown por call site
# - idle estimado como diferença entre sends e recvs
```

**Método C — `IPM` (Integrated Performance Monitoring) ou `Score-P` + `Vampir`:**

```bash
# Com Score-P: recompilar com Score-P wrappers
# Gera trace completo com timeline de send/recv/wait por rank
# Vampir (GUI) permite ver directamente os gaps de idle entre ranks
```

#### 8.6.3 Magnitude esperada de idle

Sem load balancing optimizado:
- Malhas com blocos multi-zonais de tamanhos diferentes → idle potencialmente 10–30% por whalo2
- Com 512+ ranks: overhead de Allreduce latência pode superar 5–15% do tempo total de GMRES
- convergenceInfo: tipicamente < 1% do tempo total (mas cresce com P)

---

### 8.7 Instrumentação Proposta — Timer Completo por Categoria

Adicionar ao módulo `adflowTimers.F90` (da Parte 3.2) as seguintes variáveis adicionais:

```fortran
! --- Categorias novas ---
! Categoria B: álgebra GMRES pura
real(kind=realType) :: t_gmres_orthog     = 0.0_realType  ! VecMDot + VecAXPY
real(kind=realType) :: t_gmres_vecnorm    = 0.0_realType  ! VecNorm por iteração
real(kind=realType) :: t_gmres_matvec     = 0.0_realType  ! MatMult (MFFD shell)
integer(kind=intType) :: cnt_gmres_matvec = 0

! Categoria C: PC apply
real(kind=realType) :: t_pc_apply_ANK     = 0.0_realType
real(kind=realType) :: t_pc_apply_NK      = 0.0_realType
integer(kind=intType) :: cnt_pc_apply_ANK = 0
integer(kind=intType) :: cnt_pc_apply_NK  = 0

! Categoria D: computeTimeStepMat
real(kind=realType) :: t_computeTimeStepMat = 0.0_realType
integer(kind=intType) :: cnt_computeTimeStepMat = 0

! Categoria E: idle (barreira antes de Allreduce)
real(kind=realType) :: t_idle_allreduce   = 0.0_realType  ! idle em Allreduce
real(kind=realType) :: t_idle_whalo2      = 0.0_realType  ! idle em whalo2

! Categorias de line search
real(kind=realType) :: t_linesearch_ANK   = 0.0_realType
real(kind=realType) :: t_linesearch_NK    = 0.0_realType
integer(kind=intType) :: cnt_ls_blocketteRes_ANK = 0
integer(kind=intType) :: cnt_ls_blocketteRes_NK  = 0
```

**Para medir `t_gmres_matvec` e `t_pc_apply`, é preciso um KSP monitor em PETSc:**

```fortran
! Em setupANKKSP / setupNKKSP, registar um monitor personalizado:
call KSPSetComputeEigenvalues(ANK_KSP, PETSC_TRUE, ierr)  ! opcional

! Alternativa: usar KSPMonitorSet com callback de timer:
call KSPMonitorSet(ANK_KSP, ANKTimerMonitor, PETSC_NULL_FUNCTION, &
                   PETSC_NULL_FUNCTION, ierr)

! O callback ANKTimerMonitor seria:
subroutine ANKTimerMonitor(ksp, it, rnorm, ctx, ierr)
    use adflowTimers
    KSP ksp
    integer(kind=intType), intent(in) :: it
    real(kind=alwaysRealType), intent(in) :: rnorm
    PetscFortranAddr ctx(*)
    integer(kind=intType) :: ierr
    ! Acumular tempo de cada iteração de KSP
    cnt_gmres_ANK_iter = cnt_gmres_ANK_iter + 1
    ! (MPI_WTIME() é lido aqui mas melhor usar MatShellSetContext para passar timers)
end subroutine
```

**Para medir PC apply separado do MatMult:**

A forma mais prática sem introspeccionar PETSc é medir `KSPSolve` total menos  
o tempo acumulado de todas as chamadas a `blocketteRes` durante esse `KSPSolve`:

```
T_KSPSolve_total = T_MatMult_total + T_PC_apply_total + T_orthog_total

T_MatMult_total = Σ T_blocketteRes (apenas os chamados via MFFD durante KSP)
                = cnt_blocketteRes_during_ksp × T_avg_blocketteRes

T_PC_apply_total + T_orthog_total = T_KSPSolve_total − T_MatMult_total
```

Isto é obtido facilmente com os timers já propostos, adicionando:
1. Um contador de `blocketteRes` no início e fim de `KSPSolve`.
2. Subtrair `delta_cnt × T_avg_blocketteRes` do total do KSPSolve.

---

### 8.8 Escalonamento: OpenMP Threads vs. MPI Ranks

A tabela abaixo resume o comportamento esperado de cada categoria quando se  
aumenta o número de OpenMP threads por rank (mantendo ranks fixo) vs. quando se  
aumenta o número de MPI ranks (mantendo threads fixo):

| Categoria | Aumentar OMP threads | Aumentar MPI ranks | Notas |
|---|---|---|---|
| **A — blocketteResCore** (compute) | ✅ escala bem | ✅ escala bem | `collapse(2)` OMP; MPI divide domínio |
| **A — whalo2** (halo comm) | ❌ não ajuda | ⚠️ mais ranks = mais mensagens, menos dados por mensagem | latência cresce com P |
| **A — computePressure/Visc** (pre-halo) | ❌ sem OMP | ✅ escala (domínio menor) | candidatos a OMP futuro |
| **B — VecMDot ortogonalização** | ❌ sem OMP | ⚠️ degrada: mais latência Allreduce com P | custo O(m × log P) |
| **B — VecNorm** | ❌ sem OMP | ⚠️ degrada igualmente | O(log P) por call |
| **B — VecAXPY** | ❌ PETSc pode usar OMP | ✅ volume diminui com P | local, problema menor |
| **C — ILU apply** | ❌ serial no PETSc | ⚠️ degrada (PC mais fraco com mais ranks) | não usa OMP |
| **C — AMG apply (smooth)** | ❌ serial smoothers | ⚠️ degrada nos níveis grossos | comunicação coarse |
| **D — FD coloring loop** | ❌ loop de coloring é serial | ✅ cada blocketteRes menor | candidato a OMP |
| **D — ILU factorization** | ❌ serial | ✅ problema local menor | não usa OMP |
| **D — MatAssembly** | ❌ MPI | ⚠️ overhead global | 1 call por FormJacobian |
| **E — Idle/imbalance** | ⚠️ mais threads = mais variância intra-rank | ⚠️ mais ranks = mais imbalance | difícil de controlar |

**Resumo da estratégia híbrida óptima:**

```
Objectivo: minimizar tempo total = T_compute + T_comm + T_idle

T_compute escala com 1/(N_threads × N_ranks) → máximo paralelismo beneficia

T_comm_halo escalona com N_ranks (mais mensagens) → menos ranks é melhor

T_comm_allreduce escala com log(N_ranks) × N_allreduce → menos ranks é melhor

T_idle cresce com N_ranks (mais oportunidade de imbalance) → menos ranks é melhor

T_ilu_quality DEGRADA com mais N_ranks → menos ranks é melhor para PC

→ SWEET SPOT: poucos MPI ranks (4-8 por nó NUMA), muitos threads por rank
   Ex: 2 ranks × 16 threads vs. 32 ranks × 1 thread: muito melhor com mais threads
   Razão: blocketteResCore paraleliza bem; ILU, Allreduce, halo são penalizados por mais ranks
```

**Experimento recomendado — strong scaling com configuração híbrida:**

```bash
# Fixar trabalho total, variar partições:
# Config 1: 4 ranks × 16 threads  (4 nodes × 1 rank)
OMP_NUM_THREADS=16 mpirun -n 4  python run_adflow.py

# Config 2: 8 ranks × 8 threads
OMP_NUM_THREADS=8  mpirun -n 8  python run_adflow.py

# Config 3: 16 ranks × 4 threads
OMP_NUM_THREADS=4  mpirun -n 16 python run_adflow.py

# Config 4: 64 ranks × 1 thread (puro MPI)
OMP_NUM_THREADS=1  mpirun -n 64 python run_adflow.py

# Medir: t_blocketteResCore, t_whalo2, t_KSPSolve, t_formJacobianANK, t_convergenceInfo
```

---

### 8.9 Tabela Resumo: O Que Demora Mais e Porquê

Esta tabela é uma estimativa orientativa; os valores reais dependem da malha, hardware e configuração.

| Componente | % de T_ANKStep (sem rebuild PC) | % de T_ANKStep (com rebuild PC) | Bottleneck principal |
|---|---|---|---|
| `blocketteResCore` (compute) | ~40–55% | ~25–35% | FLOPs + cache miss L2/L3 |
| `whalo2` (halo comm) | ~10–20% | ~8–12% | latência MPI + bandwidth rede |
| PC apply ILU (dentro KSP) | ~10–20% | ~5–10% | ILU forward/back serial; mais fill = mais lento |
| GMRES ortogonalização | ~5–10% | ~2–5% | latência Allreduce × m²; cresce com P |
| `FormJacobianANK` (rebuild) | 0% (sem rebuild) | ~35–50% | N_colors × blocketteRes (serial) |
| `computeTimeStepMat` | ~1–3% | ~1–3% | local + MatAssembly |
| `convergenceInfo` | < 1% | < 1% | 1 Allreduce por iter |
| Idle/imbalance total | ~5–15% | ~5–15% | load imbalance + Allreduce sync |
| Line search (ANK backtrack) | ~0–10% (raro) | ~0–5% | blocketteRes ×backtrack |
| `computePressure/Visc` (pre-halo) | ~3–6% | ~2–4% | loop células sem OMP |

**O que domina em ANK "normal" (sem rebuild, m~15 GMRES iter):**
→ **blocketteResCore** (~50%) + **whalo2** (~15%) + **PC apply** (~15%) = ~80% do tempo.  
→ Estes três são os candidatos primários a optimização.

**O que domina em ANK com rebuild PC (a cada 20 steps):**
→ **FormJacobianANK** (~40%) domina nos steps de rebuild.  
→ Com lag=20, o custo amortizado do rebuild é ~2–5 pts % do tempo médio por step.

**O que domina em NK (m~30-60 GMRES iter):**
→ **blocketteResCore** (~35%) + **PC apply** (~25%) + **GMRES ortogonalização** (~10-20%).  
→ NK tem subspace maior (60 vs. típico 15–30 em ANK) → ortogonalização pesa mais.

---

## Resumo das Referências de Código

| Rotina | Ficheiro | Linha | Notas |
|---|---|---|---|
| `ADflow.__call__` | `adflow/pyADflow.py` | 1185 | entry point Python |
| `setAeroProblem` | `adflow/pyADflow.py` | 3240 | setup AP |
| `adflow.solvers.solver` | `src/solver/solvers.F90` | ~1 | entry Fortran |
| `solveState` | `src/solver/solvers.F90` | 892 | loop principal |
| `ANKStep` | `src/NKSolver/NKSolvers.F90` | 3629 | passo ANK |
| `NKStep` | `src/NKSolver/NKSolvers.F90` | 512 | passo NK |
| `FormJacobianANK` | `src/NKSolver/NKSolvers.F90` | 1935 | monta PC ANK |
| `FormJacobianNK` | `src/NKSolver/NKSolvers.F90` | 372 | monta PC NK |
| `computeResidualNK` | `src/NKSolver/NKSolvers.F90` | 1084 | wrapper residual |
| `blocketteRes` | `src/NKSolver/blockette.F90` | 70 | dispatch residual |
| `blocketteResCore` | `src/NKSolver/blockette.F90` | 299 | kernel tiled |
| `metrics` | dento de blocketteResCore | - | face normals |
| `saSource` | `src/turbulence/sa.F90` | ~88 | SA source term |
| `saViscous` | `src/turbulence/sa.F90` | ~396 | SA viscous flux |
| `inviscidCentralFlux` | `src/solver/fluxes.F90` | 4 | JST central |
| `allNodalGradients` | `src/solver/flowUtils.F90` | 1676 | Gauss nodal grads |
| `viscousFlux` | `src/solver/fluxes.F90` | 2534 | viscous flux |
| `whalo2` | `src/utils/haloExchange.F90` | - | halo MPI |
| `convergenceInfo` | `src/solver/solvers.F90` | ~571,763 | MPI_Allreduce + print |
| `computePressureSimple` | `src/solver/flowUtils.F90` | 867 | p field |
| `computeLamViscosity` | `src/solver/flowUtils.F90` | 1201 | μ_l field |
| `computeEddyViscosity` | `src/turbulence/turbUtils.F90` | 580 | μ_t SA |
| `BS = 8` | `src/NKSolver/blockette.F90` | 9 | tile size |
| `!$OMP parallel do collapse(2)` | `src/NKSolver/blockette.F90` | 354 | único OpenMP |
