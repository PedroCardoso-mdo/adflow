# Testing & Diagnostics

> Two things only:
>   1. How to smoke-test (compile + run, no NaN).
>   2. How the volume-output debug array is wired (high-level only — full
>      spec lives in task A2 of `03_IMPLEMENTATION_PLAN.md`).
>
> Physics validation is the user's job, performed at the end after E2.

---

## 1. Smoke-Test Definition

A smoke test is the per-task "did we break anything" check. It is NOT a
physics test.

### 1.1 Compile

```bash
cd build && make -j
```

Pass = exit 0, no errors, no new warnings introduced by the task.

### 1.2 Run (when applicable)

The user provides a Python smoke script — `tests/smoke_transition.py`.
The script:
- Loads a small mesh (NLF0416 coarse is typical).
- Sets options per the task (e.g. for D1-D4, the per-path options block).
- Runs N iterations (N ~ 50-100, set by user).
- Writes a volume CGNS.

Smoke pass criteria — checked by the user:
- The script returns without exception.
- Residuals are finite (no NaN, no Inf) throughout.
- Output CGNS file exists.

Claude Code does NOT run this script. Claude Code's responsibility ends at
`make` exiting 0 and the correct `STATUS:` line.

### 1.3 The Tapenade Loop

When a task ends with `STATUS: TAPENADE NEEDED`:

```
1. Claude Code:  edits source, runs make, ends with STATUS: TAPENADE NEEDED
2. User:         runs the Tapenade regen commands (project-specific)
3. User:         re-runs `cd build && make -j`
4. User:         runs `python tests/smoke_transition.py`
5. User:         confirms run completed; if so, `/clear` and start next task
```

When a task ends with `STATUS: READY TO RUN`:

```
1. Claude Code:  ends with STATUS: READY TO RUN (make already exited 0)
2. User:         runs `python tests/smoke_transition.py`
3. User:         confirms; `/clear`; next task
```

### 1.4 The BLOCKED Path

When a task ends with `STATUS: BLOCKED — <reason>`:

```
1. Claude Code:  could not finish; reason is concrete (file:line, error
                  text, missing dependency, unclear convention)
2. User:         decides:
                   a. resolve manually (e.g. write the missing option, run
                      the missing prereq task), then re-prompt the same
                      task; OR
                   b. refine the task spec in 03_IMPLEMENTATION_PLAN.md
                      and re-prompt
3. /clear and try again
```

---

## 2. Volume Output Debug Array — High-Level

Full spec is in task A2 of `03_IMPLEMENTATION_PLAN.md`. This section
explains what A2 delivers and how to use it.

### 2.1 What A2 Delivers

After A2 lands:

- The block-level array `transitionDebug(:,:,:,:)`
  (`src/modules/block.F90:662`) is allocated with one slot per meaningful
  named quantity in the `Source` routine of `saGammaRetheta.F90`. The
  slot map is enumerated by A2 from the actual source code — not hardcoded
  here — and lives as a comment block at the top of `saGammaRetheta.F90`.
- Each slot is registered as a named field in the volume and surface
  output dispatch (the same dispatch that handles `eddy`, `cf`, `cp`).
- The Python option `storeTransitionDebug` (default `False`) gates filling
  the array. Off by default — zero overhead on production runs.

### 2.2 What's In It

Every variable in `Source` that has physical meaning and a name. That
includes:

- **State-derived**: γ_for_SA, S̃, R_T, vorticity magnitude (raw and
  limited), μ_t.
- **Correlation outputs**: Re_θt target, Re_θc, F_length, F_θt, F_onset,
  F_turb.
- **BL proxies**: θ_BL, δ_BL, λ_θ, Re_V (strain Reynolds number).
- **Assembled sources**: P_γ, E_γ, P_θt, P_SA (raw), P_SA (effective =
  γ·P_SA_raw).
- **Scaling/time**: timeScale.
- **Source Jacobian**: all 9 entries qq(1..3, 1..3) including the off-diags
  that B3 makes non-zero.

What's NOT in it: loop-local scratch variables used to build a single
following expression and then discarded.

### 2.3 How To Use It

In Python:

```python
options = {
    'turbulenceModel':       'SA',
    'useTransitionModel':    True,
    'storeTransitionDebug':  True,
    'volumeVariables': [
        'rho', 'u', 'v', 'w', 'p',           # standard NS
        'eddy', 'nutilde',                    # turbulence state
        'gamma', 'retheta',                   # transition state
        # ... append every name from the slot-map comment block ...
    ],
    'surfaceVariables': ['cp', 'cf', 'gamma', 'fonset'],
}
```

The slot-map comment block at the top of `saGammaRetheta.F90` is the
authoritative list of names available — copy from there.

### 2.4 Inspection

Open the resulting CGNS in Tecplot, ParaView, or PyVista. Every field is
just another cell-centered variable. Example PyVista session:

```python
import pyvista as pv
mesh = pv.read('out_vol.cgns')
print(mesh.array_names)               # see what's available
print(mesh['gamma'].min(),  mesh['gamma'].max())
print(mesh['fonset'].min(), mesh['fonset'].max())
print(mesh['timescale'].mean())       # B1 verification, if relevant
```

For airfoil cases (NLF0416, S809), the surface CGNS is usually enough —
plot `gamma` and `fonset` along chord and compare to paper figures.

---

## 3. Known Bugs / Open Questions (Snapshot)

For full status see `02_IMPLEMENTATION_STATUS.md`.

| # | Site                              | Status                              | Resolved by |
|---|-----------------------------------|-------------------------------------|-------------|
| 1 | `saGammaRetheta.F90:510-559`      | Confirmed — off-diag Jacobian zero  | B3          |
| 2 | `saGammaRetheta.F90:430-431`      | UNVERIFIED — timeScale convention   | B1 (verify) |
| 3 | feedback loop                     | Symptom downstream of #1, maybe #2  | auto after B1+B3 |
