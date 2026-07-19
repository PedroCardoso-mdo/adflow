"""Shared configuration and helpers for the SA-noft2-Gamma-Retheta (SA-GR /
SA-sLM2015) derivative regression tests.

This module plays the role reg_default_options.py / reg_aeroproblems.py /
reg_test_utils.py play for the SA and Euler cases, plus SA-GR-specific
helpers that isolate the transition state variables (gamma, Re_theta_t) and
the SA<->transition coupling blocks flagged in docs/audits/06_adjoint_wiring.md.

State-vector layout (verified against NKSolvers.F90:getStates): the flattened
state/residual vectors are variable-fastest, w(i,j,k,1:nw), so python offset
``iCell*nw + l`` holds variable l (0-based). For SA-GR nw = 8 with
itu1/itu2/itu3 = 6/7/8 (Fortran, 1-based) = nuTilde/gamma/reThetat.

The grid/restart/FFD files below are placeholders until the user supplies a
converged transition case (see docs/audits/08_test_prep.md, checklist).
"""

import os
from collections import OrderedDict

import numpy
from baseclasses import AeroProblem
from baseclasses.testing import getTol

baseDir = os.path.dirname(os.path.abspath(__file__))

# --------------------------------------------------------------------------
# Case inputs — 2026-07-19: switched from the SA tutorial-wing grid to the
# AR5 plain-wing case (`ar5_plain_wing_vol_L3.cgns`) used for the full
# ANK->CANK->NK convergence-strategy investigation
# (`/home/mdo/Desktop/Run/.../3D_Plain_Wing/full_ladder_1e-5_production/`,
# STRATEGY.md ladder). The tutorial-wing case never exercised the NK path at
# all (ANK-only test plan); this case validates derivatives about a real,
# NK-converged (albeit stall-limited, see l2convergence note below) state.
# Grid copied to input_files/ (gitignored, same convention as the downloaded
# SA input files) rather than referencing the external Desktop/Run path.
# Restart written by generate_sagr_restart.py (same convention as the SA
# suite: tests linearize about the restart state via getResidual()).
# --------------------------------------------------------------------------
sagrGridFile = os.path.join(baseDir, "../../input_files/ar5_plain_wing_vol_L3.cgns")
sagrRestartFile = os.path.join(baseDir, "../../input_files/ar5_plain_wing_sagr_dp.cgns")
# NOTE: stale (tutorial-wing FFD) -- AR5 has no FFD yet. Fine for the
# residual/jacVecProd/adjoint-vs-CS tests (they use alpha/mach/P/T DVs, not
# shape); needed before running TestCmplxStepSAGR's geometric "shape" DV.
sagrFFDFile = os.path.join(baseDir, "../../input_files/mdo_tutorial_ffd.fmt")

# AR5 plain-wing AeroProblem — MUST match the conditions
# generate_sagr_restart.py converged (same numbers as the validated
# production ladder run, run_strategy.py / full_ladder_1e-5_production).
ap_sagr_ar5_wing = AeroProblem(
    name="ar5_plain_wing",
    alpha=0.0,
    beta=0.0,
    mach=0.2,
    altitude=0.0,
    areaRef=0.1,
    chordRef=1.0,
    xRef=0.0,
    yRef=0.0,
    zRef=0.0,
    evalFuncs=["cl", "cd", "cmz", "drag"],
)

# backward-compat aliases used by the test files
ap_sagr_tut_wing = ap_sagr_ar5_wing
ap_sagr_flatplate = ap_sagr_ar5_wing

# Aero DVs: mach/P/T drive uInf/muInf and wInf(itu2/itu3) — this is the axis
# that exercises audit-06 finding F1 (vorticity-limiter uInf/muInf path) and
# the farfield transition BC derivatives.
sagrAeroDVs = ["alpha", "mach", "P", "T"]

# Solver options for the SA-GR derivative tests. 2026-07-19: replaced the
# ANK-only tutorial-wing option set with the validated AR5 production ladder
# (ANK -> CANK -> NK, STRATEGY.md / run_strategy.py) so derivatives are
# checked about a real NK-converged state, not just an ANK one. Every
# pre-existing SA-GR-model-specific option (transition* etc.) is kept
# explicit with its rationale; the solver-path options below now mirror
# run_strategy.py exactly. Option semantics: docs/architecture.md Part 2.
sagrBaseOptions = {
    # ---- case ----
    "gridfile": sagrGridFile,
    "restartfile": sagrRestartFile,
    "equationtype": "RANS",
    # the model itself — enables the 8-state system (5 meanflow + nuTilde,
    # gamma, reThetat)
    "turbulencemodel": "SA-noft2-Gamma-Retheta",
    # validated AR5 value (run_strategy.py) -- NOT the pyADflow auto-default
    # [1e4,10,1e4]; conditions the coupled 3x3 DADI block solve directly
    # (saGammaRetheta.F90 scaleNu/scaleGamma/scaleReTheta), do not change
    # without re-validating on this case (see STRATEGY.md).
    "turbResScale": [10000.0, 0.1, 1.0e-4],

    # ---- model requirements enforced/warned by pyADflow ----
    # blockette residual kernels drifted from the SA-GR model (audit-06 F7)
    # and are force-disabled by pyADflow anyway; explicit for documentation
    "useblockettes": False,
    # wall functions replace the low-Re near-wall behaviour the gamma
    # correlations rely on — rejected outright for this model
    "usewallfunctions": False,
    # the model is SA-*noft2*: the ft2 laminar-suppression term must be off
    # (gamma takes over that role via Eq. 41 production multiplication)
    "useft2sa": False,

    # ---- freestream / farfield BC definition (transition-critical) ----
    # Tu_inf [fraction]: sets wInf(itu3) via the Rethetat_inf correlation and
    # enters the gamma correlations; must be > 0 (default 0.001). Matched to
    # the validated AR5 production run (0.25%).
    "turbintensityinf": 0.0025,
    # not set by the validated AR5 run_strategy.py -- leave at the
    # ADflow/SA default rather than the tutorial-wing-tuned 0.009.
    # "eddyvisinfratio": <default>,
    # run_strategy.py DOES set this explicitly (0.2, not the pyADflow
    # default of 1.0 -- confirmed by diffing printed option dumps 2026-07-19,
    # first pass here missed it). PC/solver-side only -- forced to 1.0 inside
    # adjointUtils.F90 before the PC/adjoint is assembled, so this never
    # touches AD/CS derivatives, but it does affect the primal solve.
    "acousticscalefactor": 0.2,

    # ---- transition model options added on this branch ----
    # first-order upwind for gamma/reThetat convection (paper §IV.A, CLAUDE.md
    # rule 5). False = second order (sharper front, adjoint hole: audit-06 F8)
    "transitionfirstorderupwind": True,
    # helicity-based crossflow source D_scf on the Re_thetat equation
    # (P&Z Eqs. 15-26). OFF here (matches the validated AR5 production run):
    # on a real 3D wing the helicity term divides by max(velMag,1e-10) TWICE
    # and blows up in near-wall/near-stagnation cells where velMag is
    # genuinely tiny -- unvalidated per CLAUDE.md rule 3. The tutorial-wing
    # case kept this True only because D_scf is identically zero on its
    # ~2D mesh anyway (free adjoint-block coverage); that reasoning does NOT
    # extend to AR5, a real 3D wing, so dR[reThetat]/d(nuTilde) via D_scf is
    # untested here -- known gap, not a regression.
    "transitioncrossflow": False,
    # surface roughness height h [m] in the crossflow correlation (Eq. 17);
    # 3.3e-6 = smooth surface (default)
    "transitionroughnessheight": 3.3e-6,
    # source-term dt restriction (P&Z Eq. 59): caps lambda_source*dt in
    # DADI/turbKSP. Solver-side only — does not change the residual/adjoint
    "transitionsrcdtrestrict": True,
    # the cap value for the restriction above (lambda_source*dt <= this)
    "transitionsrcdtlimit": 0.9,
    # how lambda_source is estimated from the 3x3 source Jacobian:
    # "eigenvalue" (largest positive eig, default) or "gershgorin" (bound)
    "transitionsrcdteigmode": "eigenvalue",
    # deactivate the source-dt restriction after N clean turbKSP iterations
    # in the second-order regime (only engages if ANKSecondOrdSwitchTol is
    # raised from its 1e-16 default; see architecture.md D-A2-3)
    "srcdtdeactivateiters": 5,
    # per-variable gamma/reThetat update damping in DD-ADI (P&Z Alg. 2):
    # back-off factor and safety cap on the back-off loop (LHS-only)
    "transitiondamptheta": 0.99,
    "transitiondampmaxiter": 10000,
    # drop the production term from the SA row of the low-order PC/DADI
    # Jacobian when approxSA is active (sst_dev "approxTurb" device).
    # PC/LHS-only — never touches the adjoint Jacobian (audit-06 F2 caveat)
    "transitionuseapproxsa": True,
    # reference length l [mesh units] in the vorticity limiter (P&Z
    # Eqs. 52-53). Negative = auto (uses ap.chordRef=1.0 here); this is the
    # uInf/muInf path of audit-06 finding F1
    "transitionreflength": -1.0,
    # DADI coupling mode for the 3 turbulence equations: "full" (3x3 block),
    # "transition" (SA scalar + gamma/reThetat 2x2), "decoupled" (3 scalar)
    "turbdadicoupled": "full",
    #
    # Turb-ANK physicality line search: separate relative tolerance for
    # reThetat (replaces ANKPhysicalLSTolTurb for that variable)
    "ankphysicallstolretheta": 0.99,
    # minimum step-factor floor for gamma in the physicality check —
    # prevents step collapse in laminar regions where gamma -> 0
    "omegamingamma": 0.05,

    # ---- startup smoother (before ANK engages) ----
    "mgcycle": "sg",  # single grid: multigrid is deferred for this model (CLAUDE.md rule 3)
    "smoother": "DADI",
    "cfl": 1.7,
    "cflcoarse": 1.0,
    "resaveraging": "alternate",
    "liftindex": 3,  # AR5 wing convention (run_strategy.py) -- was 2 for the tutorial wing
    # NOT set by run_strategy.py -- true pyADflow default (1), not the
    # tutorial-wing-tuned 3 (confirmed by diffing printed option dumps
    # 2026-07-19; leftover from the old tutorial-wing config, missed on
    # first pass here)
    "nsubiter": 1,
    # DADI turbulence sub-iterations per flow iteration
    "nsubiterturb": 3,
    # capped to exactly reproduce full_ladder_1e-5_production/run_full.log
    # at outer iter 282 / Iter Tot 3594 (totalRes=9.9223314015988036E-04),
    # the specific point requested for this restart
    "ncycles": 3594,
    "monitorvariables": ["cpu", "resrho", "resturb", "cl", "cd", "totalr", "scaledtotalr"],

    # ---- solver path: full ladder, ANK -> CANK -> NK (STRATEGY.md /
    # run_strategy.py) -- 2026-07-19: replaces the tutorial-wing test plan's
    # ANK-only path so derivatives are validated about a real NK state.
    "usenksolver": True,
    "useanksolver": True,
    # NOT set by run_strategy.py -- true pyADflow default (1000), not the
    # tutorial-wing-tuned 1.0 (confirmed by diffing printed option dumps
    # 2026-07-19). Both trivially engage ANK at iter 1, but 1000 also lets
    # ANK reactivate if totalR ever spikes back above totalR0 during a
    # partial blowup/recovery event; 1.0 would not allow that fallback.
    "ankswitchtol": 1000.0,
    # decoupled path: turbulence handled by DADI inside ANK (path 1 of the
    # roadmap); set False to use Turb-ANK KSP instead (path 2)
    "ankuseturbdadi": True,
    # phase 2 (CANK): validated switch point, NEVER earlier (STRATEGY.md:
    # "coupling before 1e-5 stagnates, with or without LS")
    "ankcoupledswitchtol": 1e-5,
    # phase 3 (CSANK): validated switch point. In practice the early NK
    # trigger below (nkswitchtol 1e-5) crosses before this does, so CSANK is
    # typically skipped entirely (see STRATEGY.md "Gatilho NK cedo").
    "anksecondordswitchtol": 1e-6,
    # essential -- without ANKUnsteadyLSTol/ANKPhysicalLSTol at these values
    # the coupled phase stagnates at step ~0.01 (STRATEGY.md)
    "ankunsteadylstol": 2.0,
    "ankphysicallstol": 0.8,
    "ankphysicallstolturb": 0.99,
    # inner Newton iterations for the turbulence sub-solve inside ANK
    "anknsubiterturb": 1,
    # CFL ramp for the pseudo-transient continuation
    "ankcfl0": 5.0,
    "ankcfllimit": 1e6,  # matched to STRATEGY.md
    # use the approximate (first-order) SA Jacobian in the ANK PC — this is
    # what transitionUseApproxSA modifies for the SA-GR source terms
    "ankuseapproxsa": False,
    # essential: FD-colored PC is unusable on SA-GR (STRATEGY.md)
    "ankadpc": True,
    "nkadpc": True,
    "nksubspacesize": 300,
    # local pseudo-time-step type for ANK's CFL calc (solver-side only, not
    # part of R(w) -- excluded from AD/CS); matched to STRATEGY.md
    "ankchartimesteptype": "VLR",
    # phase 4 (NK): early-engage trigger validated 2026-07-18 (Eq.59/
    # Algorithm 2 reactivation inside NK, NKSolvers.F90) -- 3 orders earlier
    # than the pre-2026-07-18 4.2e-8/8e-8 value
    "nkswitchtol": 1e-5,
    # Eq. 59 (P&Z) source-dt restriction on the coupled path (already set
    # above too, kept here for the NK-phase reader)
    "transitionsrcdtrestrict": True,
    # 2026-07-16: Eq.59 reactivation + Algorithm 2 damping *inside* NK
    # (default True; explicit for the record) -- this is what makes the
    # early nkswitchtol above survivable
    "transitionnk": True,
    # full double precision so the tests linearize about exactly this state
    "solutionprecision": "double",

    # ---- convergence / adjoint ----
    # Back to a tight/effectively-disabled value -- attempted 1e-3 first
    # (2026-07-19) as an auto-stop, but L2Convergence (totalR < L2Conv *
    # totalR0, solvers.F90:1755) is checked every iteration regardless of
    # solver phase, and 1e-3 is LOOSER than ANKCoupledSwitchTol/nkswitchtol
    # (1e-5): the run satisfied it during plain ANK descent and stopped
    # after only 20 iterations, never reaching CANK/NK at all. There is no
    # single L2Convergence value that both survives past the 1e-5 switch
    # and auto-stops before the indefinite post-NK stall (see
    # full_ladder_1e-5_production/PURPOSE.md) -- so convergence depth is
    # controlled by `ncycles` above instead (capped just past the point
    # that run was already flat and sub-1e-3).
    "l2convergence": 1e-14,
    "l2convergencecoarse": 1e-4,
    "adjointl2convergence": 1e-14,
    # the adjoint must include the transition equations — never freeze
    "frozenturbulence": False,

    # ---- output: everything needed for derivative validation ----
    "writevolumesolution": True,
    "writesurfacesolution": True,
    # transitionDebug fields (CLAUDE.md rule 4: all transition diagnostics
    # go to the volume CGNS) plus the standard resrho/resturb/intermittency/
    # fonset fields already used in the validated production run
    "volumevariables": [
        "resrho", "resturb", "intermittency", "fonset",
        "thetabl", "deltabl", "transdelta", "transvelmag", "duds",
        "nutsa", "rethetatilde", "transvortlim",
    ],
    "surfacevariables": ["cp", "mach", "cf", "vx", "vy", "vz", "cfx", "cfy", "cfz"],
}


# --------------------------------------------------------------------------
# State-vector block helpers
# --------------------------------------------------------------------------
def getStateBlocks(CFDSolver):
    """Return (nw, blocks) where blocks is an ordered dict mapping block name
    to the list of 0-based variable offsets inside each cell's nw-slot chunk.

    Reads nw/nt1/nt2 from the compiled module so it works identically on the
    real and complex builds.
    """
    nw = int(CFDSolver.adflow.flowvarrefstate.nw)
    nt1 = int(CFDSolver.adflow.flowvarrefstate.nt1)  # Fortran 1-based
    nt2 = int(CFDSolver.adflow.flowvarrefstate.nt2)

    if nw != 8 or (nt2 - nt1) != 2:
        raise ValueError(
            "SA-GR tests expect the 8-state SA-noft2-Gamma-Retheta layout "
            "(nw=8, three turbulence variables); got nw=%d, nt1=%d, nt2=%d" % (nw, nt1, nt2)
        )

    blocks = OrderedDict()
    blocks["meanflow"] = list(range(0, nt1 - 1))  # rho, rho*u, rho*v, rho*w, rho*E
    blocks["nuTilde"] = [nt1 - 1]  # itu1
    blocks["gamma"] = [nt1]  # itu2
    blocks["reThetat"] = [nt1 + 1]  # itu3
    return nw, blocks


def maskStateVector(vec, nw, offsets):
    """Zero every entry of a flattened (variable-fastest) state-shaped vector
    except the given 0-based per-cell variable offsets."""
    masked = numpy.zeros_like(vec)
    for off in offsets:
        masked[off::nw] = vec[off::nw]
    return masked


# --------------------------------------------------------------------------
# Assert helpers (mirroring reg_test_utils.py style: record norms in the
# BaseRegTest handler; the same keys are recorded by the AD (real-build) and
# CS (complex-build) test classes, so training with the real build and then
# running the cmplx_* tests compares AD against complex step).
# --------------------------------------------------------------------------
def assert_coupling_blocks_allclose(handler, CFDSolver, seed=314, mode=None, h=None, **kwargs):
    """Column-by-column forward products dR/dw * e[block], with the result
    norms split by residual row block.

    Covers, block by block, the coupling partials from audit 06 §4:
      - dR[nuTilde]/dw[gamma]      (the classically-missed gamma->SA block)
      - dR[gamma]/dw[nuTilde]      (r_T = nuTilde*fv1/nu chain)
      - dR[gamma]/dw[reThetat]     (Flength/Rethetac correlation chain)
      - dR[reThetat]/dw[meanflow]  (timescale/lambda_theta/F_thetat chain)
      - dR[nuTilde]/dw[reThetat]   (structurally zero -- regression-guards it)
      - dR[reThetat]/dw[nuTilde]   (nonzero only through crossflow D_scf, F5)
    Wall/farfield BC halo chains are inside these products because
    applyAllTurbBC runs within the residual evaluation.
    """
    rtol, atol = getTol(**kwargs)
    nw, blocks = getStateBlocks(CFDSolver)

    wDotFull = CFDSolver.getStatePerturbation(seed)

    extraArgs = {}
    if mode is not None:
        extraArgs = {"mode": mode, "h": h}

    for colName, colOffsets in blocks.items():
        wDot = maskStateVector(wDotFull, nw, colOffsets)
        resDot = CFDSolver.computeJacobianVectorProductFwd(wDot=wDot, residualDeriv=True, **extraArgs)

        for rowName, rowOffsets in blocks.items():
            rowDot = maskStateVector(resDot, nw, rowOffsets)
            key = "||dR[%s]/dw[%s] * eDot||" % (rowName, colName)
            handler.root_print(key)
            handler.par_add_norm(key, rowDot, rtol=rtol, atol=atol)


def assert_transition_xdvdot_allclose(handler, CFDSolver, ap, seed=1.0, mode=None, h=None, **kwargs):
    """Forward products dR/d(aeroDV), with the result norms split by residual
    row block, for each aero DV on the AeroProblem.

    mach/P/T reach the transition residuals through uInf/muInf (the P&Z
    Eq. 52-53 vorticity limiter -- audit-06 finding F1) and through the
    farfield ghost state wInf(itu2/itu3) (turbBCRoutines farfield form), so
    the gamma/reThetat rows of these products are the targeted check for the
    transition-model BC and flow-condition derivatives.
    """
    rtol, atol = getTol(**kwargs)
    nw, blocks = getStateBlocks(CFDSolver)

    extraArgs = {}
    if mode is not None:
        extraArgs = {"mode": mode, "h": h}

    for aeroDV in ap.DVs.values():
        key = aeroDV.key
        xDvDot = {key: seed}
        resDot = CFDSolver.computeJacobianVectorProductFwd(xDvDot=xDvDot, residualDeriv=True, **extraArgs)

        for rowName, rowOffsets in blocks.items():
            rowDot = maskStateVector(resDot, nw, rowOffsets)
            name = "||dR[%s]/d%s||" % (rowName, key)
            handler.root_print(name)
            handler.par_add_norm(name, rowDot, rtol=rtol, atol=atol)


def assert_bwdfast_blocks_allclose(CFDSolver, seed=314, atol=1e-16):
    """Row-block-seeded consistency between the full reverse mode (_b) and the
    reverse-fast state-only mode (_fast_b, computeJacobianVectorProductBwdFast).

    Seeding one residual row block at a time checks each transposed coupling
    column separately; a failure isolated to the gamma/reThetat seeds is the
    signature of the autoEditReverseFast.py push/pop stripping hazard flagged
    in docs/audits/sst_dev_lessons.md (watch item 1).
    """
    nw, blocks = getStateBlocks(CFDSolver)
    dwBarFull = CFDSolver.getStatePerturbation(seed)

    for rowName, rowOffsets in blocks.items():
        dwBar = maskStateVector(dwBarFull, nw, rowOffsets)

        wBar = CFDSolver.computeJacobianVectorProductBwd(resBar=dwBar, wDeriv=True)
        wBarFast = CFDSolver.computeJacobianVectorProductBwdFast(resBar=dwBar)

        numpy.testing.assert_allclose(
            wBar, wBarFast, atol=atol, err_msg="BWD vs BWDFast, resBar seeded on %s rows" % rowName
        )


def assert_coupling_dot_products_allclose(handler, CFDSolver, seed=314, **kwargs):
    """Blockwise transpose (dot-product) tests: for selected (column, row)
    pairs, seed wDot on one variable block and dwBar on another and require
    dwBar^T (dR/dw wDot) == wDot^T (dR/dw^T dwBar). Each pair isolates one
    off-diagonal Jacobian block in both AD modes simultaneously.
    """
    rtol, atol = getTol(**kwargs)
    nw, blocks = getStateBlocks(CFDSolver)

    wDotFull = CFDSolver.getStatePerturbation(seed)
    dwBarFull = CFDSolver.getStatePerturbation(seed + 1)

    # (column block seeded in wDot, row block seeded in dwBar)
    pairs = [
        ("gamma", "nuTilde"),  # dR[nuTilde]/d(gamma): Eq. 41 production coupling
        ("nuTilde", "gamma"),  # dR[gamma]/d(nuTilde): r_T chains
        ("reThetat", "gamma"),  # dR[gamma]/d(reThetat): correlation chains
        ("nuTilde", "reThetat"),  # dR[reThetat]/d(nuTilde): crossflow D_scf (F5)
        ("meanflow", "reThetat"),  # dR[reThetat]/d(meanflow): timescale/F_thetat
        ("gamma", "gamma"),  # diagonal sanity
        ("reThetat", "reThetat"),  # diagonal sanity
    ]

    for colName, rowName in pairs:
        wDot = maskStateVector(wDotFull, nw, blocks[colName])
        dwBar = maskStateVector(dwBarFull, nw, blocks[rowName])

        dwDot = CFDSolver.computeJacobianVectorProductFwd(wDot=wDot, residualDeriv=True)
        wBar = CFDSolver.computeJacobianVectorProductBwd(resBar=dwBar, wDeriv=True)

        dotLocal1 = numpy.sum(dwDot * dwBar)
        dotLocal2 = numpy.sum(wDot * wBar)

        key = "Dot product test for w[%s] -> R[%s]" % (colName, rowName)
        handler.root_print(key)
        handler.par_add_sum(key, dotLocal1, rtol=rtol, atol=atol)
        handler.par_add_sum(key, dotLocal2, rtol=rtol, atol=atol, compare=True)
