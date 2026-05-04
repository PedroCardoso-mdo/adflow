# CLAUDE.md — ADflow SA-LM2015 Transition Model

## Goal

Finish implementing the Piotrowski & Zingg (2020) γ-Re̅θt transition model in
ADflow. The model is ~60% done; this roadmap drives it to completion. All
four solver paths must remain runtime-selectable and runnable:

1. Decoupled — turbulence solved by DADI (3×3 block coupled)
2. Decoupled — turbulence solved by Turb-ANK
3. Coupled ANK — flow + turbulence in one Newton-Krylov system
4. Coupled NK — terminal Newton-Krylov solver

**Definition of done for any task: code compiles, smoke test runs N
iterations without NaN.** Physics correctness is verified by the user at
the end, not task-by-task.

## Per-Task Workflow

For every task you are asked to do:

1. Read this file (auto-loaded) plus the ONE task block from
   `03_IMPLEMENTATION_PLAN.md` matching the requested T-ID.
2. Read other files (`01_PAPER_REFERENCE.md`, `04_ARCHITECTURE.md`, etc.)
   ONLY if the task block's `Context:` line says to.
3. Make the edits described in the task block.
4. Run `cd build && make -j`. Fix any compile errors. Iterate until clean.
5. Update the corresponding row in `02_IMPLEMENTATION_STATUS.md`.
6. Commit: `git commit -am "T-ID: short msg"`.
7. End your turn with one of these status lines, exactly:
   - `STATUS: TAPENADE NEEDED` — task touched AD-relevant code; user reruns
     Tapenade then `make`, then the smoke script.
   - `STATUS: READY TO RUN` — no AD impact; user runs `make` (if not
     already run) then the smoke script.
   - `STATUS: BLOCKED — <reason>` — could not finish, need user input.
8. Stop. Do not start the next task. The user will `/clear` and start fresh.

## Hard Rules

1. **One task per session.** Read `02_IMPLEMENTATION_STATUS.md` first to
   confirm the task is not already ✅.
2. **Do not modify the SA model directly.** Transition is a modifier — γ
   multiplies SA production only (Eq. 41).
3. **Skip multigrid (T1.6) and crossflow (T2.5).** User defers both.
4. **All transition diagnostics go to the volume CGNS** via the
   `transitionDebug` array. No new ASCII debug files. No per-cell printouts.
5. **First-order upwind** for γ and Re̅θt convection (paper §IV.A).
6. **Adjoint linearization is frozen** at this stage. Touch
   `src/adjoint/output*` only when a task explicitly says to regenerate.
7. **Compile + run = done.** Do not attempt physics verification yourself.
   Do not interpret results. Do not "improve" beyond the task spec.
8. **Do not load files outside the task block's `Context:` list.** Token
   discipline matters — every task block is self-contained for a reason.
9. **Paper wins** when paper and code disagree (`01_PAPER_REFERENCE.md`).

## File Locations (most-edited)

| What            | Where                                                  |
|-----------------|--------------------------------------------------------|
| Main residual   | `src/turbulence/saGammaRetheta.F90`                    |
| Helpers         | `src/turbulence/saGammaRethetaHelpers.F90`             |
| Constants       | `src/modules/paramTurb.F90`                            |
| Input params    | `src/modules/inputParam.F90`                           |
| Block array     | `src/modules/block.F90` (`transitionDebug`)            |
| BCs             | `src/turbulence/turbBCRoutines.F90`                    |
| Init            | `src/initFlow/initializeFlow.F90`                      |
| Output dispatch | locate via `grep -rn "case ('eddy')" src/`             |
| Python wrapper  | `python/pyADflow.py`                                   |
| AD generated    | `src/adjoint/output{Forward,Reverse,ReverseFast}/`     |

## AD-Relevant Code (Triggers `TAPENADE NEEDED`)

If a task edits any of these, the AD copies are stale and must be regenerated:

- `src/turbulence/saGammaRetheta.F90`
- `src/turbulence/saGammaRethetaHelpers.F90`
- `src/turbulence/turbBCRoutines.F90` (transition BC paths)
- Any source term, residual, or Jacobian evaluation reachable from
  `block_res` / `slaveTurbAPI`.

Pure infrastructure changes (input parameters, allocation, output
registration, Python wrapper) do NOT trigger Tapenade.

## Repo Layout

```
adflow/
├── CLAUDE.md                       ← this file, auto-loaded
├── 00_CLAUDE_CODE_SETUP.md         ← how to install/operate Claude Code
├── 01_PAPER_REFERENCE.md           ← paper equations — load on demand
├── 02_IMPLEMENTATION_STATUS.md     ← status + ordered task list — read first
├── 03_IMPLEMENTATION_PLAN.md       ← per-task specs — load ONLY your task
├── 04_ARCHITECTURE.md              ← ADflow internals — load on demand
├── adflow/                         ← ADflow python side
└── src/                            ← Fortran sources
```

## Build Commands

```bash
cd build && make -j                  # Claude Code runs this
# (user reruns Tapenade if STATUS was TAPENADE NEEDED)
# (user reruns make if Tapenade was rerun)
# (user runs python smoke script — see 05_TESTING_AND_DEBUG.md §1)
```
