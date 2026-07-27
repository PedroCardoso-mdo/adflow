# Full 3-layer verification (partials + adjoint + complete/CS) for sa/smooth/hard — 2026-07-27

**Problem:** The harness existed (prior sessions) but had never been run end-to-end for all
three cases (`sa` plain-SA baseline, `smooth`, `hard`) across all three verification layers
(block-partials, adjoint-mode, and the decisive complete-mode complex-step check over all 83
design variables). `reg_bcm.py` also had a real bug — `bcmBaseOptions` never set `useNKSolver`
(pyADflow defaults it to `False`) or most other ANK/NK options, so every prior `hard`-variant run
in this repo had silently run ANK-only, diverging from `dev/run_bcm_case.py`'s option dict
despite a comment claiming they matched.

**Solution:**
- Fixed `reg_bcm.py`: added `useanksolver, ankcfl0, ankchartimesteptype, anksecondordswitchtol,
  ankadpc, ankcfllimit, usenksolver, nkadpc, eddyvisinfratio, infchangecorrection` to
  `bcmBaseOptions`, matching `run_bcm_case.py` exactly. Re-ran `hard --mode cs` after the fix:
  NK genuinely engages now (confirmed via raw iteration log switching ANK/SANK → NK), reproducing
  the known exp-sqrt line-search pinning, not a new artifact of the fix.
- Extended `diag_partials_bcm.py` to support `--variant sa` (previously `smooth`/`hard` only).
- Added `dev/diag_adjoint_isolated_bcm.py`: solves each functional's adjoint on a **fresh**
  AeroProblem instead of the chained `cl,cd,cmz,drag` sweep. Finding: `hard`'s `drag` adjoint
  converges fine on its own (`fail=False`, `ddrag/dmach=2148.3`) — it was never actually being
  attempted in the chained sweep, just zeroed by `skipAfterFailedAdjoint` after `cmz` (the only
  functional that genuinely diverges, even in isolation) failed first.
- Ran all three layers for all three variants, logs saved under `dev/logs/`:
  - **Partials** (`ad/fd/cs/bwd/dot` × `sa/smooth/hard`): Tapenade correct at machine precision
    (1e-12 to 1e-16) for all three, no exceptions. FD-vs-AD at a single h=1e-7 shows the expected
    truncation-error signature (100-200% on the large-magnitude nuTilde/nuTilde block, same order
    across all 3 variants including plain `sa` — confirms it's generic FD noise, not model-specific).
  - **Adjoint mode** (all 83 DVs, real build): `sa`/`smooth` converge cleanly (fail=False).
    `hard` still has `cmz` diverging (fail=True) even with the ANK/NK fix — root cause is
    `hard`'s restart state itself never having converged, not the solver options.
  - **Complete mode** (`--mode cs`, all 83 DVs × 4 functionals = 332 derivatives per variant,
    run in parallel on distinct physical cores — AMD Ryzen 9 7950X, logical CPU N and N+16 are
    the same physical core, confirmed via `taskset -pc`):
    - `sa`: 324/332 rows agree with CS to <1% (essentially exact).
    - `smooth`: 320/332 rows agree with CS to <1% (essentially exact).
    - `hard`: **0/332** — every single row is wrong, not just `cmz`/`drag` as earlier
      partial-DV sweeps had suggested. Root cause candidate: `[L2 reached: ...]` is **identical
      (4.500e-09) on every one of the 83 DVs** in `hard`'s log — not DV-dependent as it should
      be if each DV's complex primal were genuinely re-converging independently. Points to state
      contamination between DVs (possibly `resetFlow` not fully clearing after `hard`'s known
      kink instability is hit on an early DV) rather than 79 independent physics failures. **Not
      yet root-caused further** — open follow-up.
  - Full 332-row side-by-side table (`sa`/`smooth`/`hard`) generated and archived.

**Where the results live:** `/home/mdo/Desktop/Run/MDO_PhD/Transition/SA_BCM/
Verification_tuturial_mesh/` (standalone bundle, rebuilt from scratch this session — scripts,
all 23 logs, refs, the 332-row table, and `PURPOSE.md` fully rewritten with the current state).
The plain-SA baseline numbers were also appended (not replacing anything) to the sibling
SA-Gamma-Retheta repo's `.../SaGammaReTheta/summary_table.txt`, since plain SA is identical
physics regardless of which transition-model repo verifies it.

**Not done:** root-causing `hard`'s complete-mode state-contamination bug (see follow-up above);
`hard` still needs an actual converged restart (`NKLSRelax` untried) before its `cmz` adjoint can
be trusted at all.
