# ANK Timing Flowcharts (Paper Version)

This version is intentionally simplified for publication figures.
It keeps only principal steps, decision points, and timer insertion locations.
Residual evaluations are shown as explicit calls, without internal residual breakdown.

## Diagram A: Nonlinear Loop (ANK-only) with Pre/Post

```mermaid
flowchart TD
    N0[Pre-processing box\ninput/state initialization\nmonitor reset] --> N1[nonlinearIteration loop]
    N1 --> N2[call ANKStep]
    N2 --> N3[convergence/stop check]
    N3 -->|not converged| N1
    N3 -->|converged| N4[Post-processing box\nfinal outputs and report]
```

Notes:
- Only ANK is shown in the nonlinear loop (no RK/NK branches).
- Pre and post boxes are included for paper-level workflow context.

## Diagram B: One ANKStep (Principal Timed Steps)

```mermaid
flowchart TD
    S0[ANKStep start\nSEC1 parent scope] --> S1{Need Jacobian/PC update?}

    S1 -->|yes| S2[computeTimeStepMat\nSEC2 timer inserted]
    S2 --> S3[FormJacobianANK total\nSEC3 timer inserted]
    S3 --> S31[Residual call\nSEC4 timer inserted]
    S3 --> S32[PC setup\nSEC5 timer inserted]
    S32 --> S4[Outside-KSP setup\nSEC36 timer inserted]

    S1 -->|no| S4[Outside-KSP setup\nSEC36 timer inserted]

    S4 --> S5[KSPSolve total\nSEC6 timer inserted]

    S5 --> K1[Branch 1: measurable callbacks\nMatMult branch SEC25 + SEC27/28/29]
    S5 --> K2[Branch 2: PETSc internal remainder\nall non-callback work]

    K1 --> S6[State/local updates\nSEC13 accumulation]
    K2 --> S6[State/local updates\nSEC13 accumulation]

    S6 --> S7[computeUnsteadyResANK\nSEC10 + SEC30/31/32]
    S7 --> S8{Unsteady norm acceptable?}

    S8 -->|no| S9[Line-search/backtracking\nrepeat computeUnsteadyResANK]
    S9 --> S8
    S8 -->|yes| S10[Final residual call\nSEC12 timer inserted]

    S10 --> S11{Turbulence update needed?}
    S11 -->|yes| S12[ANK_turbUpdate total\nSEC24 timer inserted]
    S11 -->|no| S13[ANKStep end]
    S12 --> S13[ANKStep end]
```

Notes:
- Residual is explicitly shown where called, but without internal residual decomposition.
- KSPSolve is intentionally split into only two branches:
  - measurable callback branch,
  - PETSc internal remainder branch.

## Timer Labels Included in This Paper View

- SEC1: ANKStep total (parent scope)
- SEC2: computeTimeStepMat
- SEC3: FormJacobianANK total
- SEC4: FormJac residual call
- SEC5: PC setup
- SEC6: KSPSolve total
- SEC10: computeUnsteadyResANK total
- SEC12: final residual call
- SEC13: local updates accumulation
- SEC24: turbulence update total
- SEC25 + SEC27/28/29: measurable KSPSolve callback branch
- SEC30/31/32: computeUnsteady internal timed blocks
- SEC36: outside-KSPSolve setup
