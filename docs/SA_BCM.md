# SA-BCM Turbulence Model Implementation

## Overview

SA-BCM (Spalart-Allmaras with BCM transition) adds laminar-turbulent transition prediction to the standard SA model via an intermittency function `gamma`.

## Key Files

| File | Purpose |
|------|---------|
| `src/turbulence/sa.F90` | Main implementation in `saSource` (lines 296-351) |
| `src/modules/inputParam.F90` | Parameter definitions |
| `src/NKSolver/blockette.F90` | Newton-Krylov integration (line 1089+) |
| `src/modules/blockPointers.F90` | `Tgamma` pointer (line 124) |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `use_SABCM` | false | Enable SA-BCM model |
| `SABCM_Exp` | true | Use exp-sqrt formulation (false=tanh) |
| `SABCM_Const1` | 0.002 | Transition constant 1 |
| `SABCM_Const2` | 0.06 | Transition constant 2 |
| `SABCM_TU` | 0.01 | Freestream turbulence intensity |
| `SABCM_S0_tanh` | 0.5 | S0 for tanh formulation |
| `SABCM_fsmooth` | 0.08 | Smoothing factor for tanh |
| `SABCM_maxsmooth` | 50.0 | Max smoothing parameter |

## Physics

The model computes intermittency `gamma` (stored as `Tgamma`) that multiplies the SA production term:

1. Calculate critical Reynolds number: `Re_theta_c` from `SABCM_TU`
2. Compute local `Re_theta` from vorticity and wall distance
3. Form transition terms:
   - `tterm1 = (Re_theta - Re_theta_c) / (Re_theta_c * SABCM_Const1)`
   - `tterm2 = fv1*chi / SABCM_Const2`
4. Two formulations:
   - **Exp-sqrt:** `gamma = 1 - exp(-(sqrt(tterm1) + sqrt(tterm2)))`
   - **Tanh:** `gamma = 0.5 * (1 + tanh((tterm1 + tterm2 - S0) / fsmooth))`

## Data Flow

- **Inputs:** velocity gradients, wall distance `d2Wall`, density, viscosity `rlv`
- **Output:** `Tgamma(i,j,k)` in `flowDoms` array

## Parallelization

- No direct MPI calls in `sa.F90`
- Uses ADflow's existing halo exchange via `applyAllTurbBC`
- Block-based parallelism through `flowDoms` structure

## Potential NaN Sources

1. **Line 323:** `Re_vorty = sqrtVort * rho / rlv` - division by `rlv`
2. **Line 327:** Division by `Re_theta_c * SABCM_Const1`
3. **Line 335:** Large exponential from `stransition`
4. **Line 345:** Division by `SABCM_fsmooth`
