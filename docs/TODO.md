# TODO — deferred items

Items here are decided-not-urgent. Each links to its `audits/design-decisions.md` entry if one
exists.

- **Extract `use_SABCM` out of `sa.F90` into its own file.** Considered during initial KB
  setup as a way to make "never touch sa.F90" strictly true again, but rejected as out-of-scope
  extra work not requested by the current task. See `audits/design-decisions.md` §
  "sa.F90 edit-scope rule".
- **`blockette.F90`/`sa.F90` consistency check.** ~~No automated diff exists...~~ DONE:
  `tests/reg_tests/dev/diag_blockette_bcm.py` (raw output) and `test_blockette_bcm.py`
  (registered), both variants, ported from the sibling SA-GR repo's `test_blockette_sagr.py`.
  Not yet run (ready-to-run harness only, see `tests/reg_tests/README_BCM.md`) -- until it's
  run, `useBlockettes=True` (the current default in `reg_bcm.py`, matching the previously-used
  AR5 script) is unverified for SA-BCM, same situation SA-GR was in before this check caught its
  drift.
- **Cherry-pick the divide-by-zero guards from `sa-bcm-timing`** (`e4b0b6ae` + `8e5a4e8f`,
  paired). Deferred until the base harness here reveals whether those guards are actually
  needed for the cases this branch runs. See `adjoint-trace.md`.
- **Remove the stale "external module not seen by Tapenade" comment** in `sa.F90:~331` and
  `blockette.F90:~1220` (`saSource`, right before `k_max = max(stransition, xminn)`) — describes
  a mechanism that isn't what the current inline code does (harmless, just confusing). Blocked in
  the 2026-07-24 session by a `guard_protected_files.sh` hook carried over from a different repo's
  Claude Code settings that blindly matches any `src/turbulence/sa.F90` path regardless of repo;
  needs either a corrected hook scope or an explicit user override to touch.
- **Complex build has never been done in this repo.** `make -f Makefile_CS
  PETSC_ARCH=complex-debug` needs to be run at least once before `run_bcm_tests.sh cs`/`adjoint`
  or any `dev/*` `--mode cs` script can work.

## NK/ANK convergence follow-ups (from an SA-GR NK-mods review, 2026-07-24)

SA-GR's NK/ANK convergence work (column scaling, Eq. 59 source-dt restriction, Algorithm 2
per-node damping) exists to handle two new *transported* state variables (γ, Re̅θt) that span
~13 orders of magnitude and have stiff source Jacobians of their own. SA-BCM has no new
transported state (`nw` stays 6; its `tTgamma` is an algebraic, cell-local multiplier folded
into the existing SA production term), so that machinery does not transfer as written. Two
narrower ideas, NOT tied to having a new state variable, might still be worth prototyping --
**only after** `tests/reg_tests/run_bcm_tests.sh` is green (tuning a solver against unverified
derivatives risks chasing an implementation bug instead of a real stiffness problem, and the
manuscript itself flags an unresolved adjoint bug):

- **NK line-search relaxation** (Armijo `alpha` 1e-2→1e-3, backtrack factor 3.0 in
  `NKSolvers.F90`): a generic Newton-globalization fix for `Step` pinning near a stiff front,
  not specific to new states. ~~Cheapest thing to try~~ **DONE (2026-07-24):** implemented as the
  `NKLSRelax` option (off by default) — see `docs/task-log/2026-07-24-paper-verification-and-nk-relax.md`.
  Not yet tried on an actual pinned run (harness itself hasn't been executed yet).
- **A per-cell source-dt restriction on the *existing* nuTilde equation**, structurally like
  Eq. 59 but landing on nuTilde's diagonal Jacobian term instead of a new variable's: the KS
  smoothmax blend in `tTgamma` (`SABCM_maxsmooth=50.0`) is steep, and its derivative
  (`dtTgamma` in `sa.F90`'s hand-linearization, `#ifndef USE_TAPENADE` block) is folded straight
  into nuTilde's diagonal. Same "locally huge source Jacobian eigenvalue" pathology Eq. 59
  targets, just on an existing row.
