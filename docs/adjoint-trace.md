# Adjoint / AD trace — SA-BCM

## Current sync status (as of commit `85f40a0b`, branch `sa-bcm`)

- `git diff b8a27a25 HEAD -- src/turbulence/sa.F90` → empty.
- `git diff b8a27a25 HEAD -- src/adjoint/outputForward/sa_d.f90 src/adjoint/outputReverse/sa_b.f90 src/adjoint/outputReverseFast/sa_fast_b.f90` → empty.
- Both the hand-written `use_SABCM` block in `sa.F90` and its Tapenade-generated derivatives
  (`sa_d.f90`, `sa_b.f90`, `sa_fast_b.f90`) were introduced together in commit `b8a27a25`
  ("Add SA-BCM version (from 2.12.2 dev 129)"). **No Tapenade rerun is currently needed.**

## What triggers a rerun

Any hand-edit inside the `use_SABCM`-gated block of `sa.F90` (lines 294–413) invalidates
`sa_d.f90`/`sa_b.f90`/`sa_fast_b.f90` and requires Tapenade to be rerun before trusting
`computeJacobianVectorProductFwd`/`Bwd`/`BwdFast` again. This is a user action (per
`Makefile_tapenade`), not something to do by hand-editing the generated files — see CLAUDE.md
rule 6 and `.claude/settings.json`'s `ask` permissions on `src/adjoint/output*/temp_*`.

## Known divide-by-zero risk (not yet on this branch)

Sibling branch `sa-bcm-timing` has two paired commits not present on `sa-bcm`:
- `e4b0b6ae` "Add guards to prevent division by zero in SA-BCM" — hand-edits `sa.F90`/
  `blockette.F90` to add `max(rlv, 1e-20)`, `max(Re_theta_c*SABCM_Const1, 1e-20)`,
  `max(SABCM_fsmooth, 1e-20)` guards.
- `8e5a4e8f` "AD runed" — the corresponding Tapenade rerun updating `sa_d.f90`/`sa_b.f90`/
  `sa_fast_b.f90`.

These guards are **not yet on `sa-bcm`**. If they get cherry-picked/merged here, both commits
must land together — the guard commit alone would desync AD from the hand-written code exactly
like any other un-Tapenade'd edit.

## `blockette.F90` is a manual mirror, not Tapenade-differentiated separately

`src/NKSolver/blockette.F90:1004-1250` duplicates the `saSource` SA-BCM logic for the
non-Tapenade ANK/DADI solve path. It is **not** run through Tapenade — if `sa.F90`'s
`use_SABCM` block changes, `blockette.F90` must be updated by hand to match, and there is no
automated check that the two stay consistent. Any task that edits one should diff against the
other.

## Paper-vs-code verification (2026-07-24)

Line-by-line check of `sa.F90:294-413`/`blockette.F90:1004-1250` against
`docs/papers/AIAA20202714_SABCMPartI.md`'s Appendix (authoritative "copy-for-code" formulation)
and `docs/papers/AIAA20202706_BCMtransitionmodel.md`, including hand-rederiving the
`#ifndef USE_TAPENADE` PC-path derivatives. **Result: faithful implementation, no bugs found.**
Term2 correctly uses the eddy-viscosity ratio `fv1*chi` (= `mu_T/mu`), not raw `chi` — the single
trap both papers explicitly warn is "wrong by a factor of 250" if missed. `Re_v` is local
(2020-2714's formulation), not a boundary-layer max (2020-2706's), matching the LCTM-consistent
choice the reconciliation notes recommend. The only deviations from the papers are the two
already-documented, deliberate smoothings (log-sum-exp smooth-max on Term1's `max(...,0)`, and
the `SABCM_Exp=False` tanh blend as an alternative to the paper's own `SABCM_Exp=True` exp-sqrt
formula) — both are user choices to keep, not defects. Spot-checked `sa_d.f90` and confirmed real
(non-stub) derivative code for the SABCM block, consistent with the "in sync" status above.
`Makefile_tapenade`'s `-head` spec for `saSource` was also confirmed complete: Term1/Term2/tTgamma
are pure local-variable functions of the declared independent vars (`w`, `rlv`, `d2wall`) feeding
the declared dependent var (`scratch`), so Tapenade's activity analysis differentiates them
automatically — no missing `-vars` entries, no rerun needed. Full writeup:
`docs/task-log/2026-07-24-paper-verification-and-nk-relax.md`.

## Verification ladder

See `tests/reg_tests/README_BCM.md` for the full FD/CS/adjoint verification sequence this
branch's task is built around.
