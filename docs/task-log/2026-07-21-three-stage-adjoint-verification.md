# Three-stage low-level adjoint verification (plain SA) — 2026-07-21

**Problem:** Before trusting any adjoint-total derivative check
(`evalFunctionsSens` vs. trained reference / complex-step), the
underlying Tapenade-generated linearization itself needed a raw,
functional-independent verification: does forward-mode agree with
reverse-mode, does reverse-mode agree with its fast variant, and does the
shared forward linearization actually match the true derivative (AD vs.
CS vs. FD)? None of this had been run at the raw `dR/dw` level before —
only higher-level `evalFunctionsSens` totals had been checked
(`../current-task.md`).

**Solution:** Ran a 3-stage ladder against the raw
`computeJacobianVectorProdFwd`/`Bwd` Fortran API, on plain SA (`nw=6`),
both tutorial-wing and AR5 meshes:
1. Reverse↔forward dot-product consistency — PASSED both meshes.
2. Reverse vs. fast-reverse (`_b` vs `_b_fast`) — PASSED (no-op check on
   plain SA; no gamma/reThetat rows to seed).
3. 3-way AD/FD/CS forward check, plus an FD step-size sweep confirming
   first-order convergence to the AD/CS value — PASSED both meshes (AD
   and CS agree to full precision; single-step FD is off as expected,
   sweep confirms `h→0` convergence to the same value).

Full detail (files, exact commands, numeric results) written up in
`../VERIFICATION/three-stage-verification.md` rather than duplicated
here.

**Files created/touched and why:**
- `tests/reg_tests/sanity_check_bwdfast_stage1.py` — Stage 1/2, raw API,
  tutorial-wing, isolates `_b`-vs-`_fast_b` without any reg_sagr
  machinery, to rule that out as a variable before trusting the
  generalized harness below.
- `tests/reg_tests/sanity_check_partials_sa.py` — Stage 1/2, generalizes
  `reg_sagr.py`'s block-splitting dot-product/`_fast_b` helpers to plain
  SA (`nw=6`) so the same harness later serves SA-GR (`nw=8`) unchanged.
  `--mesh {tutorial,ar5}`.
- `tests/reg_tests/check_3way_fwd.py`, `check_3way_fwd_sweep.py` — Stage
  3, AD vs. FD vs. CS on the raw Jacobian-vector product, `--mesh`/
  `--build` selectable since real and complex ADflow can't coexist in one
  process.
- `tests/reg_tests/generate_sagr_restart.py`, `reg_sagr.py` — updated
  `altitude=` → explicit `P=101325.0, T=288.15` (needed so `P`/`T` can be
  DVs per `baseclasses`' constructor-argument-only DV rule) and
  `sagrFFDFile` repointed at the real AR5 FFD (was a stale tutorial-wing
  placeholder); `getStateBlocks` generalized to accept both `nw=8`
  (SA-GR) and `nw=6` (plain SA) layouts so this session's plain-SA
  scripts could reuse the existing block-splitting helpers unchanged.
- `docs/VERIFICATION/` (new directory) — relocated `adjoint-trace.md`
  (from `docs/`) and `debugging_derivatives.md` (from
  `docs/ADFLOW_BASE/`) here, alongside the new
  `three-stage-verification.md`, so all adjoint/AD verification material
  lives under one heading instead of split across the KB.
  `docs/adflow-vs-paper-solver.md` moved to
  `docs/SA_GAMMA_RETHETHA_BASE/` (solver-vs-paper content sits with the
  rest of the paper-grounded material). All cross-references updated in
  `CLAUDE.md`, `docs/README.md`, `docs/TODO.md`, `docs/architecture.md`,
  `docs/nondimensionalization.md`, `docs/audits/adjoint_audit_2026-07-07.md`,
  `docs/audits/00_inventory.md`, `docs/audits/06_adjoint_wiring.md`.

**Verification:** All three stages pass on both tutorial-wing and AR5
meshes for plain SA, per the numeric tables in
`../VERIFICATION/three-stage-verification.md`. No code changes to the
solver itself in this task — verification scripts + doc reorg only.

**Follow-ups:**
- **DONE (same day):** Run this same 3-stage ladder for SA-GR (`nw=8`) —
  see `2026-07-21-three-stage-adjoint-verification-sagr.md`.
- **DONE (2026-07-22):** Restructure `tests/reg_tests/` to match upstream
  ADflow's official test conventions/fixtures, with generator scripts —
  landed as the standardized suite, see
  `2026-07-22-sagr-test-suite-standardization.md`.
