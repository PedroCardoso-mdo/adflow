# ADflow SA-BCM Timing Branch

## Project Overview

ADflow CFD solver with SA-BCM turbulence model and ANK timing instrumentation.

## Key Files

| File | Purpose |
|------|---------|
| `src/turbulence/sa.F90` | SA-BCM implementation (lines 296-351) |
| `src/NKSolver/ankProfiling.F90` | Timing instrumentation |
| `src/modules/inputParam.F90` | SA-BCM parameters |
| `adflow/pyADflow.py` | Python API (`use_ANKProfiling` line 5862) |

## Known Issue: NaN with >64 Ranks

When running with >64 MPI ranks, SA-BCM + timing crashes with NaNs.

**Suspect locations in `src/turbulence/sa.F90`:**
- Line 323: `Re_vorty = sqrtVort * rho / rlv` (division by viscosity)
- Line 327: Division by `Re_theta_c * SABCM_Const1`
- Line 335: Large exponential from `stransition`
- Line 345: Division by `SABCM_fsmooth`

**Root cause hypothesis:** Finer domain decomposition exposes cells with zero/invalid viscosity `rlv`.

## Build

```bash
make
```

## Enable Features

```python
solver.setOption('useANKProfiling', True)  # Enable timing
solver.setOption('useSABCM', True)         # Enable SA-BCM
```

## Workflow Instructions

After completing a work section, autonomously evaluate if any docs need updating based on decisions made. If updates are needed, propose the changes and ask user for approval before editing.

Docs to consider: `docs/MEMORY.md`, `docs/PREFERENCES.md`, topic-specific docs.

**After each code change:**
1. Check for errors (compile/lint if applicable)
2. Fix any errors found
3. Commit with descriptive message

## Documentation

See `docs/` for detailed documentation:
- `PROJECT_TREE.md` - Project structure
- `SA_BCM.md` - SA-BCM model details
- `TIMING.md` - Timing instrumentation
- `known_problems/` - Known issues (e.g., `NAN_ISSUE_64RANKS.md`)
- `fixed_problems/` - Resolved issues with root cause, solution, and related commits
