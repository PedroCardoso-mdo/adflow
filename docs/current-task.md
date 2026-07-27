# Current Task

> Overwrite this file's body whenever the active task changes (CLAUDE.md: one task per
> session). When a task finishes, write its summary as a new file in `task-log/` (see
> `task-log/README.md`), add it to that index, and replace this file's body with the next task.

## Task: Root-cause `hard`'s complete-mode (complex-step) state contamination

**Started:** 2026-07-27
**Status:** open, not yet started

### Objective

`hard`'s complete-mode verification (`dev/diag_full_derivatives_bcm.py --mode cs`, all 83 DVs)
fails on **every single one of the 332 derivatives** (4 functionals × 83 DVs) — not just
`cmz`/`drag` as earlier partial-DV sweeps suggested. The smoking gun: the log's
`[L2 reached: ...]` is **identical (4.500e-09) on every one of the 83 DVs**, which should not
happen if each DV's complex primal is genuinely re-converging independently from a fresh
`resetFlow`. This points to state contamination between DVs — likely triggered once `hard`'s
known exp-sqrt kink instability is hit on an early DV — rather than 79 independent physics
failures. See `docs/task-log/2026-07-27-full-verification-sa-smooth-hard-complete.md` for the
full context; `sa` and `smooth` do NOT show this pattern (324/332 and 320/332 rows agree with CS
to <1%), so this is `hard`-specific.

### Context (read only what's listed — CLAUDE.md rule 8)
- `docs/task-log/2026-07-27-full-verification-sa-smooth-hard-complete.md` — how this was found.
- `docs/task-log/2026-07-24-reg-bcm-missing-ank-nk-options.md` — the ANK/NK fix already applied
  (this bug appeared AFTER that fix, with NK genuinely engaged, so it is not the same issue).
- `tests/reg_tests/dev/diag_full_derivatives_bcm.py` — the script, `run_cs`/`resetFlow` logic.
- `tests/reg_tests/dev/logs/cs_bcm_hard_COMPLETO.log` — the raw evidence (search `L2 reached`).

### Suggested approach (not started)
- Instrument/rerun a handful of DVs individually (fresh process each time, e.g. `--shape 0`
  alone vs `--shape 40` alone) to see if the corruption is present from DV #1 or only appears
  after some earlier DV hits the kink.
- Check whether `resetFlow` (or whatever restores the restart state between DVs in
  `diag_full_derivatives_bcm.py`) actually clears the complex perturbation/derivative state, not
  just the real-valued flow state.

### Definition of done
Either the contamination is fixed and `hard`'s complete-mode numbers become DV-dependent and
meaningful again, or the mechanism is understood and documented well enough to know whether it's
a script bug (fixable) or a genuine property of `hard`'s unconverged restart (not fixable without
a converged restart first — see the still-open `NKLSRelax` follow-up in
`../Converging_tuturial_mesh/README.md`).
