# reg_bcm.py silently ran without NK (and with drifted ANK settings) — 2026-07-24

**Problem:** All `hard`-variant verification results produced so far via
`dev/diag_full_derivatives_bcm.py` (and by extension anything using `reg_bcm.py`'s
`bcmBaseOptions`) showed the primal getting permanently stuck in ANK, never switching to NK,
even at `totalR` many orders of magnitude below the printed NK switch threshold. Traced to:
`reg_bcm.py`'s `bcmBaseOptions` never set `useNKSolver`, and `pyADflow.py:5788` defaults it to
**`False`**. The file's own comment claimed it "mirrors `dev/run_bcm_case.py`'s option dict" —
that was false; `run_bcm_case.py` sets `useNKSolver=True` explicitly (plus `useANKSolver=True,
ANKCFL0=1.0, ANKCharTimeStepType="VLR", ANKSecondOrdSwitchTol=1e-3, ANKADPC=True,
ANKCFLLimit=3e5, NKADPC=True, eddyVisInfRatio=1e-10, infchangecorrection=True`), none of which
were present in `reg_bcm.py`. Every `hard`/`smooth` adjoint and CS-derivative run in this
session's `Verification_tuturial_mesh` standalone bundle (and the in-repo
`dev/logs/{adjoint,cs}_bcm_*.log`) was therefore computed under ANK-only dynamics that diverge
from the config that generated the restart files themselves — not directly wrong (ANK is still a
valid solver), but not the intended/matching configuration, and it masked whether NK's known
line-search pinning on `hard`'s kink (see prior task-log) shows up in this harness too.

**Fix:** Added the missing keys to `reg_bcm.py`'s `bcmBaseOptions` so it genuinely matches
`run_bcm_case.py`: `useanksolver, ankcfl0, ankchartimesteptype, anksecondordswitchtol, ankadpc,
ankcfllimit, usenksolver, nkadpc, eddyvisinfratio, infchangecorrection`. (The complex-build
override in `diag_full_derivatives_bcm.py` that force-disables `ankadpc`/`nkadpc` for the CS
build, since there's no AD preconditioner there, still applies on top — unaffected by this fix.)

Re-ran `hard --mode cs --shape 0,71` after the fix: NK now genuinely engages (confirmed via raw
iteration log, iter ~55 switches from ANK/SANK to `NK` rows) and reproduces the same line-search
pinning behavior (`Step` collapsing to 0.01–0.10, residual stalling ~3e-1) already characterized
in `run_bcm_case.py`'s cold-start run — consistent with the known exp-sqrt kink, not a new bug.

**Also found (separate, real):** `diag_full_derivatives_bcm.py`'s functional order
(`EVAL_FUNCS = ["cl","cd","cmz","drag"]`) combined with `skipAfterFailedAdjoint=True` (default)
means once `cmz`'s adjoint KSP diverges, `drag`'s `solveAdjoint` is **never actually attempted**
on that AeroProblem — it gets zeroed by the skip logic, not by its own failure
(`pyADflow.py:4090-4096`). Verified via a new isolated-solve script,
`dev/diag_adjoint_isolated_bcm.py` (fresh AP per functional, `evalFunctionsSens(ap, funcsSens,
evalFuncs=[f])`): `drag` solved alone on `hard` **converges** (`fail=False`,
`ddrag/dmach=2148.3`, much closer to the CS value ~2735 than the previously-reported zeroed
"0.0"). `cmz` solved alone still diverges (`fail=True`) — that part of the earlier diagnosis
holds. Net effect: the `Verification_tuturial_mesh` bundle's "cmz/drag both untrustworthy for
hard" conclusion overstates it for `drag` specifically — `drag`'s own adjoint is fine, it was
just never given the chance to run in the chained sweep.

**Not yet done:** propagate the `reg_bcm.py` fix to the three standalone copies under
`/home/mdo/Desktop/Run/MDO_PhD/Transition/SA_BCM/*/reg_bcm.py`, and re-run/re-document the
`hard` adjoint+CS tables with NK genuinely engaged and with `drag` solved via the isolated path
(or `skipAfterFailedAdjoint=False`) instead of the chained one.
