# Finished Task Log — SA-γ-Re̅θt

One file per finished task (a "case"), added here over time as work closes —
not appended into one ever-growing file. When a task in `../current-task.md`
finishes: write `<YYYY-MM-DD>-<slug>.md` in this directory using the
template below, add a row to the index, then reset `../current-task.md` to
the next task.

## Index

| Date | Task | File |
|---|---|---|
| 2026-07-09 | Doc consolidation (`findings/` redistribution) | [2026-07-09-doc-consolidation.md](2026-07-09-doc-consolidation.md) |
| 2026-07-21 | Three-stage low-level adjoint verification (plain SA) | [2026-07-21-three-stage-adjoint-verification.md](2026-07-21-three-stage-adjoint-verification.md) |
| 2026-07-21 | Three-stage low-level adjoint verification (SA-GR, nw=8) — Stages 1&3 PASS, Stage 2 FAILS on reThetat | [2026-07-21-three-stage-adjoint-verification-sagr.md](2026-07-21-three-stage-adjoint-verification-sagr.md) |
| 2026-07-22 | Fix reverse-fast `dR[reThetat]/dw[meanflow]` (lambdaTheta in-place clamp) — Stage 2 SA-GR now PASSES | [2026-07-22-fastb-retheta-lambdatheta-fix.md](2026-07-22-fastb-retheta-lambdatheta-fix.md) |
| 2026-07-22 | Standardize SA-GR test suite: testflo tests + `run_sagr_tests.sh`, CS cross-call fix, FD h-sweep/expected-fail, crossflow always ON, `dev/` folder | [2026-07-22-sagr-test-suite-standardization.md](2026-07-22-sagr-test-suite-standardization.md) |
| 2026-07-23 | Rotating-frame consistency: transition model uses the relative frame (vortLim ref velocity, BL scales, helicity), rotation-form warning; bit-identical no-op when Ω=0, real 13/13 + CS 10/10 | [2026-07-23-rotating-frame-consistency.md](2026-07-23-rotating-frame-consistency.md) |
| 2026-07-23 | SA-GR full total-derivative adjoint test (mirrors `test_adjoint.py`) + `gammaForSA`→[xminn,1+xminn] clamp (Tapenade rerun) + `adjointMaxIter=3000` + complex-build AD-PC overrides + `flatplate`→`tut_wing` rename + `dev/diag_full_derivatives.py`. Partials real 13/0 CS 10/0, adjoint real 4/0; adjoint-CS agrees ~1e-8 (OPEN: complex re-converge depth vs 5e-9 tol) | [2026-07-23-sagr-full-adjoint-test.md](2026-07-23-sagr-full-adjoint-test.md) |
| 2026-07-24 | SA-GR CS adjoint check: tol → 5e-8 (complex-build CS floor); mach + geom `drag` made NON-BLOCKING (reported, not asserted) — alpha + cl/cd/cmz block. Deep-iter study (all logs in `Verification_tuturial_mesh/SaGammaReTheta/`): adjoint reproduces CS to 5–7 sig figs; mach derivative does not settle (rel ~1e-3, OPEN); twist[0]/shape[13,32,33] share an L2=2.648e-14 plateau (OPEN) | [2026-07-24-sagr-cs-tolerance-nonblocking.md](2026-07-24-sagr-cs-tolerance-nonblocking.md) |
| 2026-07-24 | Blockette SA-GR residual synced to block path: fixed sign-flipped first-order upwind in `saGammaRethetaAdvection` (was factor-2 off) + re-synced `saGammaRethetaSource` (rotating frame, `gammaForSA` margin, lambdaTheta, helicity). New block-vs-blockette same-`w` residual test (`test_blockette_sagr.py`, wired into `run_sagr_tests.sh` `blockette` stage): all 8 vars agree to rel ≤ 1.2e-10 (PASS; was ~2.0 FAIL). pyADflow still force-off for SA-GR | [2026-07-24-blockette-sagr-residual-sync.md](2026-07-24-blockette-sagr-residual-sync.md) |

## Template

```
# <task name> — <YYYY-MM-DD>

**Problem:** what was broken/missing/in question, 1-3 sentences.

**Solution:** what was actually changed (files:lines), 1-5 bullets.

**Files created/touched and why:**
- `path/to/file` — purpose in this task, and why it's structured the way
  it is (strategy), not just what it contains.

**Verification:** what was run to confirm done-ness (compiles, smoke test
N iters no NaN, TAPENADE NEEDED?, etc.) per CLAUDE.md's definition of done.

**Follow-ups:** anything punted to `../TODO.md`, with the item name.
```
