# SA-GR complex-step adjoint check: 5e-8 tolerance + non-blocking mach/drag — 2026-07-24

**Problem:** The SA-GR full-adjoint CS regression (`test_adjoint_sagr.py`,
`cmplx_test_aero_dvs` / `cmplx_test_geom_dvs`) failed at the inherited SA
tolerances (aero `atol=5e-10, rtol=1e-8`; geom `5e-9`). A deep verification
study (`.../Verification_tuturial_mesh/SaGammaReTheta/`, all logs archived
there) established WHY, and it is not an adjoint error.

**Findings (from the study, observed — mechanism of the mach/twist floors NOT
concluded, left as open items):**
- Full adjoint-vs-CS sweep over all DVs incl. all 72 shape modes: median rel
  agreement ~1e-7–1e-12; the adjoint reproduces complex-step to 5–7 sig figs.
- Deeper iters (`--ncycles 10000`): the shape "stall" modes drop to rel ~1e-7
  once their real residual reaches L2 ~5e-16. A subset (twist[0],
  shape[13,32,33]) plateau at L2 = 2.648e-14 with rel ~1e-4–1e-5 that does not
  close with more iters.
- **mach:** the CS derivative does NOT settle. Cold-5k stopped at iter 536
  (real plateau ~4e-5, cut short); cold-10k broke the plateau (real → 7e-16)
  but the derivative STAYED at rel ~2.8e-3; warm-start (real → 1e-19) gave rel
  ~1.98 (worse). The real residual and the imaginary/derivative decouple. Cause
  not established — open item (do NOT assume preconditioner without testing).
- **drag vs cd:** `drag = cd·q∞·Sref`, so `d(drag)` and `d(cd)` carry the SAME
  relative error; a `drag` fail with a passing `cd` is purely the dimensional
  magnitude times that shared rel (e.g. twist[0]: cd/drag rel 1.5e-5).

**Decision (user):**
- Set both CS-test tolerances to **`rtol=atol=5e-8`** — matched to the CS floor
  the complex build (no AD preconditioner) can reach. At 5e-8, 15/20 official
  elements pass; the 5 that don't are mach (4) and twist[0]/drag (1).
- **mach → non-blocking:** computed and printed (CS vs ref, rel) with a message,
  but NOT asserted (`cmplx_test_aero_dvs`). Only `alpha` blocks.
- **drag → non-blocking:** in `cmplx_test_geom_dvs`, `drag` is reported (value +
  status) but NOT asserted (cd already blocks the same quantity). `cl/cd/cmz`
  block for span/twist/shape.
- Net: both CS tests go green; mach and drag stay VISIBLE in the output for the
  user to act on. Nothing is silently hidden.

**Files touched:**
- `tests/reg_tests/test_adjoint_sagr.py` — tolerances → 5e-8; alpha blocking +
  mach non-blocking (aero); drag non-blocking (geom); honest comments.
- `tests/reg_tests/dev/diag_full_derivatives.py` — dev diagnostic gained
  `--aero`/`--skip-span`/`--skip-twist`/`--warm` flags, an `L2 reached` table
  column, and standalone-layout auto-detection (runs from the archive folder).
- `tests/reg_tests/README_SAGR.md` — documents the above + the `-np` / `--cpu-set`
  gotchas learned running parallel CS jobs.

**Verification:** `py_compile` clean on both edited files. Behaviour predicted
from the archived 10k logs: `cmplx_test_aero_dvs` and `cmplx_test_geom_dvs` pass
(alpha + cl/cd/cmz within 5e-8); mach and drag print non-blocking FAIL lines.
Full run-logs + `PURPOSE.md` in `.../Verification_tuturial_mesh/SaGammaReTheta/`.

**Follow-ups (open, for the derivative verification — user-owned):**
- mach CS derivative does not converge to a stable value (rel ~1e-3); cause TBD.
- twist[0]/shape[13,32,33] share an L2 = 2.648e-14 plateau with rel ~1e-4; TBD.
- twist[0] `dcd` rel 1.5e-5 persists at deep convergence — isolate adj/CS/FD.
