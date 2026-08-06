# Three-stage low-level adjoint verification (SA-GR, nw=8) — 2026-07-21

**Problem:** The 3-stage raw-API adjoint verification ladder
(reverse↔forward dot-product consistency, reverse vs. fast-reverse
consistency, 3-way AD/FD/CS forward check) had only ever been run for
plain SA (`nw=6`). SA-GR (`nw=8`, gamma/reThetat rows) had never been run
through it — Stage 2's `_b`-vs-`_fast_b` block check in particular had no
transition content to exercise on plain SA.

**Solution:** Added a `--turbmodel {sa,sagr}` flag to the three plain-SA
ladder scripts (`sagr` pulls `reg_sagr.sagrBaseOptions` +
`ap_sagr_ar5_wing` + `sagrAeroDVs`; `reg_sagr.getStateBlocks` already
handled the `nw=8` layout) and ran all three stages on the AR5 SA-GR
restart. Results:
1. **Stage 1 (dot products): PASS** — fwd == rev to displayed precision
   for w→R, Xv→R, w→F, xV→F, and the combined product.
2. **Stage 2 (`_b` vs `_fast_b`, block-seeded): FAIL on the reThetat
   rows** (max\|_b − _fast_b\| ≈ 7.5e+2, ~44% of elements mismatched);
   meanflow, nuTilde, and gamma seeds pass. This is the first time real
   gamma/reThetat content has been exercised in this stage. Recorded as a
   numeric fact, not diagnosed (CLAUDE.md rule 7); `README_SAGR.md`'s
   "how to read failures" item 3 flags the `autoEditReverseFast.py`
   push/pop stripping as the documented suspect for exactly this
   signature.
3. **Stage 3 (AD/FD/CS forward check): PASS** — AD = CS = 3.3256409775e+13
   (exact to 11 sig figs); one-sided FD sweep converges first-order to the
   AD/CS value (rel_err 8.2e-1 → 4.9e-5 as h: 1e-8 → 1e-12).

Full detail (commands, numeric tables) written up in
`../VERIFICATION/three-stage-verification.md` ("SA-GR (`nw=8`) results").

**Files created/touched and why:**
- `tests/reg_tests/sanity_check_partials_sa.py` — added `--turbmodel`;
  for `sagr` uses `reg_sagr`'s validated AR5 option set + restart. Added a
  per-block `max|_b − _fast_b|` report loop BEFORE the hard assert so
  SA-GR's real gamma/reThetat rows produce a visible table (task asked for
  a per-block result table, not a bare PASS/FAIL).
- `tests/reg_tests/check_3way_fwd.py`, `check_3way_fwd_sweep.py` — added
  the same `--turbmodel {sa,sagr}` flag; guarded the `defaultAeroDVs` add
  so `sagr` (which already adds `sagrAeroDVs`) doesn't double-add. Print
  labels now include `turbmodel`.
- `docs/VERIFICATION/three-stage-verification.md` — added the SA-GR
  results section; updated the top status block; replaced the "SA-GR not
  yet run" bullet under "What's not yet covered" with the open Stage-2
  reThetat failure.

**Verification:** No Fortran changed — verification scripts only, so no
rebuild. All three stages ran on SA-GR/AR5 without NaN in the primal
(restart loaded, `getResidual` once, no reconvergence). Stages 1 and 3
PASS with the numbers above; Stage 2 FAILS on the reThetat rows with the
numbers above. Both real and complex builds (`libadflow.so` /
`libadflow_cs.so`) live in the same `mach` env, so Stage 3's complex run
used `ADFLOW_C` from that env (no separate `mach_cs` needed).

**Follow-ups:**
- **RESOLVED (2026-07-22):** Stage 2 `_b`-vs-`_fast_b` reThetat-row
  divergence root-caused (in-place `smoothMinMax` clamp of
  `lambdaThetaLocal` defeating the push/pop-strip) and fixed in the primal
  + Tapenade regen. Stage 2 SA-GR now PASSES. See
  `2026-07-22-fastb-retheta-lambdatheta-fix.md`.
- SA-GR's higher-level `evalFunctionsSens`-based adjoint/CS validation on
  AR5 (distinct from this raw-API ladder) still not started — tracked in
  `../current-task.md`.
