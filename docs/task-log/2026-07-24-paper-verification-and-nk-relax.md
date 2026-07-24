# Paper verification, Tapenade audit, NK line-search relax option — 2026-07-24

**Problem:** Before trusting the newly-built derivative-test harness (`tests/reg_tests/*_bcm.py`,
`dev/run_bcm_case.py`, `dev/diag_blockette_bcm.py`, prior session), the implementation itself
needed a line-by-line check against the two AIAA source papers now in `docs/papers/` (replacing
an earlier, less authoritative "differentiable reformulation" manuscript), and the Tapenade
diff-variable configuration needed confirming as complete for the `use_SABCM` block. Separately,
carry over the NK convergence-transfer analysis from a sibling SA-GR repo review: of the ideas
that could plausibly help SA-BCM's convergence, implement the cheap, generic one behind an
explicit option.

**Solution:**
- Read `docs/papers/AIAA20202714_SABCMPartI.md` (Appendix — the authoritative "copy-for-code"
  formulation) and `docs/papers/AIAA20202706_BCMtransitionmodel.md` side by side against
  `src/turbulence/sa.F90:294-413` and its `src/NKSolver/blockette.F90:1004-1250` mirror.
- **Result: the implementation is faithful to the paper.** Term1, Re_v (local, not
  boundary-layer-max — correctly follows 2020-2714 over 2020-2706 per the reconciliation notes),
  Re_theta_c, and Term2 (`fv1*chi/SABCM_Const2` = `(mu_T/mu)/chi_2`, the eddy-viscosity ratio,
  not raw `chi` — this is the single most dangerous trap the papers flag, and the code gets it
  right) all match the Appendix exactly. Production coupling (γ multiplies both the `ss`-linear
  and the `fv2`-quadratic pieces of SA production, destruction/diffusion untouched) is correct —
  verified algebraically, not just visually, since ADflow's SA production is pre-existing split
  across `term1`/`term2` locals in a way that isn't a literal transcription of the paper's `S̃`.
- **The only deviations are the two already-known, deliberate smoothings, confirmed as such and
  left untouched:** (a) the log-sum-exp smooth-max replacing the paper's `max(Term1, 0)` kink
  (`SABCM_maxsmooth`), and (b) the tanh blend (`SABCM_Exp=False`, "smooth") as an alternative to
  the paper's own `1-exp(-(√T1+√T2))` formula (`SABCM_Exp=True`, "hard" — this variant IS the
  paper's literal formula, not a separate deviation). Both hand-linearization derivatives
  (`dtterm2`, `darg_gamma`/`dtTgamma`, the `#ifndef USE_TAPENADE` PC path) were re-derived by hand
  and confirmed correct, including that Term1 has zero analytic dependence on nuTilde (Re_v
  depends only on vorticity/density/viscosity/wall-distance) — the PC path's omission of a
  Term1 contribution is mathematically exact, not an approximation.
- **Tapenade:** `Makefile_tapenade`'s `-head` spec for `saSource` (`(w, rlv, vol, si, sj, sk,
  timeRef, d2wall) > (..., scratch)`) is unchanged from before SA-BCM and correctly covers the
  chain — Term1/Term2/tTgamma are pure local-variable functions of `w`/`rlv`/`d2wall` feeding
  `scratch`, so Tapenade's whole-routine activity analysis differentiates them automatically; no
  missing `-vars` entries. Spot-checked the *already-generated* `sa_d.f90` and confirmed real,
  non-trivial derivative code exists for `tterm1d`/`tterm2d`/`arg_gammad`/`ttgammad` (not
  stubs) — confirms both that the head-spec is sufficient and that the current
  `sa_d.f90`/`sa_b.f90`/`sa_fast_b.f90` are genuinely correct, not just "in sync by `git diff`".
  **No Tapenade rerun was needed or performed** — no correctness bug was found in `sa.F90` to fix
  (a cosmetic stale-comment cleanup was attempted but blocked by this session's
  `guard_protected_files.sh` hook, which is scoped to a *different* repo's SA-model-protection
  rule and doesn't distinguish repos by path — left as-is, harmless, flagged below).
- **NK convergence:** per the sibling SA-GR repo's NK-mods review, added the one *generic* (not
  transported-state-specific) idea as an explicit, off-by-default option: `NKLSRelax`
  (`src/NKSolver/NKSolvers.F90:LSCubic`) toggles the Armijo `alpha` (1e-2 → 1e-3) and the
  turbulence-blowup pre-limit factor (2.0 → 3.0) together. Off by default — NK has no per-node
  physicality check, so an over-permissive alpha with nothing new added here is a real robustness
  risk; this is a knob to reach for only if a run is observed pinned at `minlambda` for 100+
  consecutive iterations, not a new default. The other NK idea from the review (a per-cell
  source-dt restriction on nuTilde's diagonal, mirroring SA-GR's Eq. 59 mechanism but landing on
  an *existing* variable) was deliberately **not** implemented — it needs new ANK/NK
  Jacobian-assembly plumbing (srcLambda-style eigenvalue estimate, diagonal injection), which is
  a multi-commit effort on the sibling repo and out of scope for "one task per session"; left in
  `TODO.md`.

**Files created/touched and why:**
- `adflow/pyADflow.py` — `"NKLSRelax": [bool, False]` default option
  (`defOpts`, near `NKFixedStep`) + `"nklsrelax": ["nk", "nk_lsrelax"]` namelist mapping (near
  `nkfixedstep`), mirroring the existing NK-option pattern exactly.
- `src/NKSolver/NKSolvers.F90` — `logical :: NK_LSRelax = .False.` module variable (same block as
  `useNKSolver`/`NK_fixedStep`); `LSCubic` now branches `alpha`/`turbBlowupFactor` on it instead
  of hardcoding `1.e-2_realType`/`2.0`.
- `src/f2py/adflow.pyf` — `logical :: nk_lsrelax` added to the `module nksolver` block (same
  pattern as `nk_fixedstep`).
- No `sa.F90`/`blockette.F90`/Tapenade-generated file changes — verification found nothing to fix.

**Verification:**
- Real build (`PETSC_ARCH=real-debug`, `make`) succeeded clean, no Fortran errors from the
  `NKLSRelax` plumbing.
- `pip install . --no-deps` into the mach venv succeeded (`adflow-2.11.0` → `adflow-2.12.2`,
  confirmed via `site-packages/adflow/pyADflow.py` containing `NKLSRelax` — note this **replaced**
  the mach env's previously-active editable install of the sibling SA-GR repo; re-run that repo's
  `/build` skill or `pip install -e .` there before switching back to SA-GR work).
- Did **not** run the actual derivative-test harness (`run_bcm_tests.sh` etc.) — that still needs
  `input_files/` populated and a converged restart (`generate_bcm_restart.py`), both still
  pending per the checklist in `current-task.md`. Complex build (`Makefile_CS`) also not yet done
  in this repo at all (only real has ever been built here) — needed before `run_bcm_tests.sh cs`
  or the CS ground-truth stage can run.

**Follow-ups (added to `../TODO.md`):**
- Populate `input_files/`, run `dev/run_bcm_case.py` → `generate_bcm_restart.py` →
  `run_bcm_tests.sh` (blockette first) — this is what actually proves the "faithful paper
  implementation... passes the derivative test" claim; nothing in this session ran it.
- Build the complex-step tree (`make -f Makefile_CS PETSC_ARCH=complex-debug`) — never done in
  this repo.
- The stale "external module not seen by Tapenade" comment in `sa.F90`/`blockette.F90`
  (harmless, misleading) — needs a session where the cross-repo protected-file hook is scoped
  correctly, or an explicit user override, to touch.
- The nuTilde-diagonal source-dt-restriction NK idea (bigger effort, deferred).
