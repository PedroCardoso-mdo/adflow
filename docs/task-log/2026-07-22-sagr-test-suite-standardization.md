# SA-GR derivative-test suite standardization — 2026-07-22

**Problem:** The SA-GR partial-derivative verification lived in ad-hoc `dev_test/`
scripts (`sanity_check_*`, `check_3way_fwd*`), not testflo-registered, with no
single entry point, no repeatable JSON refs, crossflow only checked in one-off
runs, and two latent issues surfacing once the checks were ported to per-block
CS: (a) a spurious ~2.15e-8 CS failure on the structurally-zero `dR[meanflow]/
dw[gamma]` and `dR[meanflow]/dw[reThetat]` blocks; (b) FD residual comparisons
that can never meet the standard element-wise tolerance on the 13-order residual.

**Solution:**
- **CS cross-call residue fix** (`reg_sagr.assert_coupling_blocks_allclose`):
  re-seat the real state — `setStates(numpy.real(getStates()))` — before each
  CS column. Sequential CS Jacobian-vector calls share a complex work buffer; a
  prior large-norm column (meanflow, O(1e8)) leaves an imaginary machine-eps
  residue (~2.15e-8) in the zero mean-flow rows of the next column. Re-seating
  clears it to exact 0 while preserving the converged linearization point. No
  tolerance loosening. (commit 341559ac)
- **FD h-sweep + expected-fail** (`reg_sagr.assert_fd_allclose_hsweep`,
  `test_jacVecProdFWD_sagr.py`): keep the tight element-wise tolerance, sweep
  `h`, pass if any step meets it, else raise reporting the best step / worst
  tol-ratio / max|Δ| (MPI-reduced). The three FD residual methods are
  `@unittest.expectedFailure` (FD can't pass element-wise on this residual; CS
  is the enforced ground truth). FD is proven correct in the norm sense (clean
  O(h) convergence). (commit 7d63aa27)
- **Crossflow always ON** *(as of 2026-07-22 — the suite later moved to the
  tutorial wing with crossflow False when the model default flipped,
  2026-07-24)* (`reg_sagr.sagrBaseOptions`): `transitioncrossflow=True`,
  grid+restart → the crossflow-converged AR5 volume CGNS. The D_scf block
  `||dR[reThetat]/dw[nuTilde]||` = 2.02e7 is now guarded by every run. Refs
  retrained about the new state. `no_train=True` on the self-contained classes
  (`TestJacVecFwdSAGRFD`, `TestJacVecBWDFastSAGR`). (commit fdf7b4e5)
- **Reorg** (commit after this log): `dev_test/` → `dev/`; add `dev/README.md`;
  add `run_sagr_tests.sh` (`all|real|cs|train|genw`) as the single entry point.

**Files created/touched and why:**
- `tests/reg_tests/reg_sagr.py` — CS state re-seat in the coupling-blocks helper;
  `assert_fd_allclose_hsweep` + `_global_tol_ratio` (MPI-reduced FD sweep);
  crossflow-ON config + crossflow-converged grid/restart.
- `tests/reg_tests/test_jacVecProdFWD_sagr.py` — FD residual asserts routed
  through the sweep; `@expectedFailure` on the 3 FD residual methods;
  `no_train=True` on the FD class.
- `tests/reg_tests/test_jacVecProdBWDFast_sagr.py` — `no_train=True` on the
  self-contained `_b`-vs-`_fast_b` class.
- `tests/reg_tests/refs/jacvec{fwd,bwd}_sagr_flatplate.json` — retrained about
  the crossflow-converged state.
- `tests/reg_tests/run_sagr_tests.sh` — one entry point; per-stage PASS/FAIL
  summary; `train`/`genw` wrappers so the whole workflow is discoverable.
- `tests/reg_tests/dev/` (was `dev_test/`) + `dev/README.md` — non-standard
  scripts: the `w` generator (`generate_sagr_restart.py`, no download server)
  and one-off diagnostics; README says what each is and how to regen `w`/JSON.
- `docs/VERIFICATION/three-stage-verification.md` — new "Canonical way to run"
  section (runner, real-vs-complex, the two do-not-"fix" invariants, mesh-swap
  recipe); old per-script sections banner-marked historical (paths now `dev/`).

**Verification:** `./run_sagr_tests.sh all` → real build 13/13, complex build CS
10/10, with crossflow ON. Retrain (`-m "train*"`) runs clean (self-contained
classes skip via `no_train`). No Fortran changed in this task (test-infra +
docs only); the underlying `_fast_b` primal fix is the separate 2026-07-22
lambdaTheta task.

**Follow-ups (to `../TODO.md`):**
- **RESOLVED (2026-07-23):** Total sensitivity `dF/dX` via the assembled
  adjoint — `test_adjoint_sagr.py` implemented and passing, see
  `2026-07-23-sagr-full-adjoint-test.md`.
- Optionally delete the now-superseded `dev/sanity_check_*` / `check_3way_fwd.py`
  (kept for their mesh/crossflow flags for now).
- Rename `_flatplate` refs → `_sagr` (cosmetic; the case is a wing, not a plate).
