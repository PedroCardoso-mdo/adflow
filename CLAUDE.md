# ADflow SA-BCM Timing Branch

## Project Overview

ADflow CFD solver with SA-BCM turbulence model and ANK timing instrumentation.

## Current Objective

Debug and verify existing OpenMP parallelization in the codebase. OpenMP is already coded but appears not to be working. Goal: get OpenMP working and demonstrate it with measurable speedup. See `docs/OPENMP.md` for detailed analysis.

## Key Files

| File | Purpose |
|------|---------|
| `src/turbulence/sa.F90` | SA-BCM implementation (lines 296-351) |
| `src/NKSolver/ankProfiling.F90` | Timing instrumentation |
| `src/modules/inputParam.F90` | SA-BCM parameters |
| `adflow/pyADflow.py` | Python API (`use_ANKProfiling` line 5862) |

## Build

Use `/build-and-test` skill to:
1. Check if AD needed, run `./AD_sh` if so
2. Run `make`
3. Activate venv: `source /home/mdo/packages_v2/mach/bin/activate`
4. Install: `pip install .`
5. Run test case in `Test_OpenMp/`

Or manually:
```bash
make
source /home/mdo/packages_v2/mach/bin/activate
pip install .
```

Note: `Test_OpenMp/` is excluded from git.

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
- `OPENMP.md` - OpenMP parallelization analysis
- `known_problems/` - Known issues (e.g., `NAN_ISSUE_64RANKS.md`)
- `fixed_problems/` - Resolved issues with root cause, solution, and related commits
