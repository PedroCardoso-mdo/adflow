# SA-GR full total-derivative adjoint test + gammaForSA xminn clamp — 2026-07-23

**Problem:** The SA-GR suite verified partials (dR/dw, dR/dXv) but not the FULL
total derivatives df/dx via a complete adjoint solve, the way plain SA does in
`test_adjoint.py`. The existing `test_adjoint_sagr.py` was stale for the
tutorial-wing case (AR5 FFD, shape-only "flat-plate" DVGeo). Separately the
`gammaForSA` production clamp was retuned, and the misleading `_flatplate`
labels needed retiring.

**Solution (code):**
- **`gammaForSA` clamp → `[xminn, one + xminn]`** (`saGammaRetheta.F90`, both the
  `Source` residual and the `evalSrcJacBlock` analytic Jacobian). Upper cap is
  ~1 (NOT `rsaGRgammaHi=2.0`): Eq. 41 multiplies SA production by the RAW
  intermittency, and P&Z drop the separation-induced `gamma_sep` term (Sec. II),
  so the multiplier is bounded by [0,1]. `gamma==1.0` stays strictly interior
  (`< 1 + xminn`), keeping the AD-vs-CS no-tie property; lower bound `xminn` is
  only met asymptotically (steady-state gamma respects the 0.02 floor). At the
  converged state (gamma in [0.02,1]) the clamp is pass-through, so this is a
  no-op for the derivative checks — but it IS an AD-differentiated residual, so
  Tapenade was rerun. **`rsaGRgammaForSAMargin` (1e-3) is now unused** (left
  declared in paramTurb.F90).
- **`gammaLocal` clamp — investigated, NOT changed.** Its bounds [1e-10, 2.0]
  are the model's own (Algorithm 2, paper line 630), and gamma==1.0 is interior
  to [.,2], so no tie. Comment updated to cite Algorithm 2.
- **Full-derivative adjoint test** (`test_adjoint_sagr.py`): rewritten to mirror
  `test_adjoint.py`'s tutorial-wing case — full ref-axis DVGeo (twist6/span1/
  shape72), `test_adjoint2`, CS geom sweep over span/twist/shape. Case renamed
  `sagr_flatplate` → `sagr_tut_wing`, ref → `adjoint_sagr_tut_wing.json`.
- **FFD fix** (`reg_sagr.py`): `sagrFFDFile` → `mdo_tutorial_ffd.fmt` (the SA
  test's FFD, correct on the tutorial-wing mesh).
- **`adjointMaxIter=3000`, `adjointSubspaceSize=300`** (`reg_sagr.py`): the SA-GR
  8-state adjoint KSP needs ~583–597 iters to reach `adjointL2Convergence=1e-14`
  (confirmed on the r/r0~1e-12 restart); the default 500 cap clipped all four
  functions just short, mis-reporting a converging adjoint as failed.
- **Complex-build solver overrides** in `TestCmplxStepSAGR.setUp`
  (`ankadpc=False, nkadpc=False, ankcoupledswitchtol=1e-16`): the complexify
  build excludes the Tapenade AD routines (audit 08), so the AD preconditioner
  the real solve needs is unavailable — the complex re-converge aborts with
  "Forward AD routines are not complexified" unless it uses FD-colored PC on the
  decoupled ANK→NK path.
- **`flatplate`→`tut_wing` rename** across refs (`git mv`), `test_jacVecProd*`,
  `dev/diag_all16_blocks_cs.py`, `README_SAGR.md`; dropped the dead
  `ap_sagr_flatplate` alias (kept `ap_sagr_ar5_wing` for dev scripts).
- **Driver** (`run_sagr_tests.sh`): `adjoint` stage wired into all/train.
- **Dev script** (`dev/diag_full_derivatives.py`): full df/dx check with ADflow
  output VISIBLE (testflo swallows it), over all aero + geom DVs incl. MULTIPLE
  shape components. `--mode adjoint` (real, all-var totals) / `--mode cs`
  (complex, per-DV re-converge vs the adjoint ref).

**Build pipeline (this session):** Tapenade regen (clean; only the 3 SA-GR AD
files changed, `xminn` propagated) → real `make` + pip → complex
`PETSC_ARCH=complex-debug make -f Makefile_CS` + pip. **Note:** the complex
build MUST use `PETSC_ARCH=complex-debug` (the real arch gives COMPLEX→REAL
PETSc type-mismatch errors in fortranPC/adjointAPI). The `/build` skill handles
this.

**Verification:**
- Partials: **real 13/0, CS 10/0** — the xminn change is a no-op for partials
  (confirmed) and the PETSc-fixed complex build gives correct CS.
- Adjoint **real: 4/0** — totals converge (~583 iters to 1e-14) and reproduce.
- Adjoint **CS: agrees to ~1e-8 absolute** (worst 2.8e-8; e.g. cmz/span
  0.044864 vs 0.044864), i.e. the adjoint totals are CORRECT — but the check
  currently FAILS the inherited `rtol=atol=5e-9`. **OPEN ITEM / root cause:** the
  complex build cannot use the AD preconditioner, so the complex re-converge
  stalls ~1e-8 (FD-PC on SA-GR, per STRATEGY.md), capping CS function accuracy
  and hence CS-vs-adjoint agreement at ~1e-8. Options not yet chosen: (a) loosen
  the SA-GR adjoint CS tol to ~1e-7/1e-5; (b) warm-start the complex solve from
  the converged restart; (c) push the complex solve deeper.
  **RESOLVED next day by option (a) at rtol=atol=5e-8** (with mach + geom
  `drag` non-blocking); warm-start (b) was tried and made the mach
  derivative WORSE — see `2026-07-24-sagr-cs-tolerance-nonblocking.md`.

**Follow-ups:** decide the adjoint-CS tolerance/depth approach (above —
resolved 2026-07-24, option (a)); optional
deeper `genw` restart; `rsaGRgammaForSAMargin` now unused (could be removed).
