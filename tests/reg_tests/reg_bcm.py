"""Shared configuration and helpers for the SA-BCM (differentiable-reformulation) derivative
regression tests.

This module plays the role reg_default_options.py / reg_aeroproblems.py / reg_test_utils.py
play for the plain SA and Euler cases, plus SA-BCM-specific helpers.

Unlike a transition model that adds new transport equations (e.g. gamma-Re_theta_t, where nw
grows and off-diagonal coupling blocks between the mean flow / nuTilde / new equations need to
be isolated), SA-BCM modifies the *existing* SA nuTilde production term in place
(src/turbulence/sa.F90:294-413, use_SABCM-gated). nw is unchanged from plain SA -- there is no
new state variable and no new residual block.

Case/mesh: deliberately matched to the sibling SA-GR harness
(adflow_sa_gamma_rethetha_paper_solver/tests/reg_tests/reg_sagr.py) for apples-to-apples
turnaround/comparability -- SAME tutorial-wing grid+FFD, SAME flow conditions
(mach=0.15, alpha=1.8), NOT the manuscript's NLF-0416 physics-validation case (that stays a
separate, later concern; this harness's job is derivative consistency, not aerodynamic
accuracy -- exact agreement with the manuscript is not required here).

Two SA-BCM variants exist and are both exercised throughout this suite, selected by the
`SABCM_Exp` flag (manuscript Eq. 8 tanh blend vs. the original exp-sqrt blend):
  - "smooth" (SABCM_Exp=False, default): the manuscript's smoothed tanh formulation
  - "hard"   (SABCM_Exp=True): the original (Mura & Cakmakcioglu) exp-sqrt formulation

Restart files are written by generate_bcm_restart.py, one per variant (the two blends converge
to different states). useBlockettes is left True (not force-disabled the way SA-GR's is) --
SA-BCM's blockette.F90 mirror of saSource has never been checked against sa.F90, unlike SA-GR's
(which WAS force-disabled after test_blockette_sagr.py caught drift). See
dev/diag_blockette_bcm.py / test_blockette_bcm.py, which exist specifically to answer that
question rather than hide it behind a default change.
"""

import os
from collections import OrderedDict

import numpy
from baseclasses import AeroProblem
from baseclasses.testing import getTol

baseDir = os.path.dirname(os.path.abspath(__file__))

# --------------------------------------------------------------------------
# Case inputs -- SA-GR's tutorial wing, wholesale (grid, FFD, flow conditions).
# See module docstring for why: apples-to-apples comparability with the
# validated SA-GR harness, not the manuscript's NLF-0416 case.
# --------------------------------------------------------------------------
bcmGridFile = os.path.join(baseDir, "../../input_files/mdo_tutorial_rans_scalar_jst.cgns")
# One restart per SABCM_Exp variant -- written by generate_bcm_restart.py.
bcmRestartFileSmooth = os.path.join(baseDir, "../../input_files/mdo_tutorial_bcm_smooth_dp.cgns")
bcmRestartFileHard = os.path.join(baseDir, "../../input_files/mdo_tutorial_bcm_hard_dp.cgns")
# backward-compat default (smooth variant) for any code that still imports the singular name
bcmRestartFile = bcmRestartFileSmooth
# same stock ADflow tutorial-wing FFD the SA-GR/SA adjoint tests use -- wraps this mesh
bcmFFDFile = os.path.join(baseDir, "../../input_files/mdo_tutorial_ffd.fmt")

# Tutorial-wing AeroProblem -- MUST match the conditions generate_bcm_restart.py converged
# (mach=0.15, alpha=1.8; same as ap_sagr_tut_wing in the sibling SA-GR repo).
ap_bcm_tut_wing = AeroProblem(
    name="mdo_tutorial_bcm",
    alpha=1.8,
    beta=0.0,
    mach=0.15,
    P=20000.0,
    T=220.0,
    R=287.87,
    areaRef=45.5,
    chordRef=3.25,
    xRef=0.0,
    yRef=0.0,
    zRef=0.0,
    evalFuncs=["cl", "cd", "cmz", "drag"],
)

# backward-compat alias -- older BCM scaffolding (reg_test files written before the tutorial-wing
# switch) imported ap_bcm_case by name.
ap_bcm_case = ap_bcm_tut_wing

# Aero DVs that exercise SA-BCM's modified terms specifically: alpha/mach feed the vorticity S
# and freestream conditions that set Re_theta (term1) and the eddy-viscosity ratio (term2); P/T
# feed rlv/rev through the freestream viscosity.
bcmAeroDVs = ["alpha", "mach", "P", "T"]

# --------------------------------------------------------------------------
# Solver options -- every SA-BCM-relevant option explicit and commented.
# Deliberately mirrors dev/run_bcm_case.py's option dict so the interactive
# driver and the registered suite never silently diverge in defaults.
# --------------------------------------------------------------------------
bcmBaseOptions = {
    "gridfile": bcmGridFile,
    "restartfile": bcmRestartFileSmooth,
    "equationtype": "RANS",
    "turbulencemodel": "SA",  # SA-BCM is not a separate turbulenceModel string, see architecture.md
    "useapproxwalldistance": True,  # required by inputParamRoutines.F90:3413-3423 when use_SABCM=True
    "useft2sa": False,  # ft2 off (SA-BCM overrides ft2=0 internally anyway; matches the AR5 script)
    # blockette.F90 has a hand-synced mirror of saSource (like SA-GR's did before it was force-
    # disabled) -- deliberately left True here so any drift shows up via diag_blockette_bcm.py /
    # test_blockette_bcm.py instead of being hidden by a default change. See module docstring.
    "useblockettes": True,
    "use_sabcm": True,
    # sabcm_exp is NOT set here -- bcmBaseOptionsSmooth/Hard below set it explicitly per variant.
    "sabcm_const1": 0.002,  # chi_1, term1 denominator scale (manuscript Eq. 4)
    "sabcm_const2": 0.02,  # chi_2, term2 denominator scale (manuscript Eq. 6) -- repo default, NOT
    # the 2.0 used on the old AR5 script (100x different; 0.02 chosen deliberately here)
    "sabcm_tu": 0.5,  # freestream turbulence intensity Tu_inf (%), feeds Re_theta_c (Eq. 5)
    "sabcm_s0_tanh": 0.5,  # S0, tanh blend center (Eq. 8)
    "sabcm_fsmooth": 0.08,  # f_smooth, tanh blend width (Eq. 8) -- calibrated, do not retune casually
    "sabcm_maxsmooth": 50.0,  # rho, KS aggregation parameter replacing max(f,0) (Eq. 10)
    "frozenturbulence": False,  # adjoint must include the SA-BCM-modified production terms
    "acousticscalefactor": 0.15,  # matches the tutorial-wing/M=0.15 case (same as sagrBaseOptions)
    "mgcycle": "sg",
    "smoother": "DADI",
    "infchangecorrection": True,
    "eddyvisinfratio": 1e-10,
    # ANK/NK solver options -- pyADflow's own default for useNKSolver is False (pyADflow.py:5788),
    # so without setting it explicitly here every run in this file was ANK-only, never reaching
    # NK, silently diverging from dev/run_bcm_case.py's option dict (the one that generated the
    # restarts). 2026-07-24: added to actually mirror run_bcm_case.py, see module docstring.
    "useanksolver": True,
    "ankcfl0": 1.0,
    "ankchartimesteptype": "VLR",
    "anksecondordswitchtol": 1e-3,
    "ankadpc": True,
    "ankcfllimit": 3.0e5,
    "usenksolver": True,
    "nkadpc": True,
    "ncycles": 10000,
    "l2convergence": 1e-14,
    "adjointl2convergence": 1e-14,
    "adjointmaxiter": 3000,
    "solutionprecision": "double",
    "writevolumesolution": True,
    "writesurfacesolution": True,
}

# Per-variant option dicts -- every downstream test/script iterates BOTH.
bcmBaseOptionsSmooth = dict(bcmBaseOptions)
bcmBaseOptionsSmooth["sabcm_exp"] = False  # tanh blend (manuscript default)
bcmBaseOptionsSmooth["restartfile"] = bcmRestartFileSmooth

bcmBaseOptionsHard = dict(bcmBaseOptions)
bcmBaseOptionsHard["sabcm_exp"] = True  # exp-sqrt blend (Mura & Cakmakcioglu original)
bcmBaseOptionsHard["restartfile"] = bcmRestartFileHard

# Same case, plain SA (use_SABCM off) -- the "before" reference for the direct term comparison
# in assert_bcm_vs_plain_sa_wdot_allclose. Uses the smooth-variant restart (arbitrary -- plain SA
# doesn't care which SABCM_Exp the restart was generated under, only that use_sabcm=False here).
bcmPlainSAOptions = dict(bcmBaseOptionsSmooth)
bcmPlainSAOptions["use_sabcm"] = False


# --------------------------------------------------------------------------
# State-vector block helpers -- ported verbatim from the SA-GR harness's
# reg_sagr.py. SA-BCM is nw=6 (no gamma/reThetat), which that module's
# getStateBlocks already documents and supports unmodified (the nw==6 branch
# exists there specifically so these helpers also run against a plain-SA-
# shaped state, which is exactly what SA-BCM is).
# --------------------------------------------------------------------------
def getStateBlocks(CFDSolver):
    """Return (nw, blocks): blocks maps block name -> 0-based per-cell variable
    offsets. SA-BCM only ever hits the nw==6 (meanflow/nuTilde) branch -- there
    is no gamma/reThetat layout to support here."""
    nw = int(CFDSolver.adflow.flowvarrefstate.nw)
    nt1 = int(CFDSolver.adflow.flowvarrefstate.nt1)  # Fortran 1-based
    nt2 = int(CFDSolver.adflow.flowvarrefstate.nt2)

    if nw == 6 and (nt2 - nt1) == 0:
        blocks = OrderedDict()
        blocks["meanflow"] = list(range(0, nt1 - 1))  # rho, rho*u, rho*v, rho*w, rho*E
        blocks["nuTilde"] = [nt1 - 1]  # itu1, the only turbulence state SA-BCM touches
    else:
        raise ValueError(
            "expected the 6-state plain-SA layout (nw=6, one turbulence variable) -- "
            "SA-BCM adds no new state variable; got nw=%d, nt1=%d, nt2=%d" % (nw, nt1, nt2)
        )
    return nw, blocks


def maskStateVector(vec, nw, offsets):
    """Zero every entry of a flattened (variable-fastest) state-shaped vector
    except the given 0-based per-cell variable offsets."""
    masked = numpy.zeros_like(vec)
    for off in offsets:
        masked[off::nw] = vec[off::nw]
    return masked


def assert_coupling_blocks_allclose(handler, CFDSolver, seed=314, mode=None, h=None, **kwargs):
    """Column-by-column forward products dR/dw * e[block], with the result norms
    split by residual row block (meanflow/nuTilde x meanflow/nuTilde -- 4 combos).

    Even without a new state variable, SA-BCM's tTgamma multiplier depends on
    rho, rlv, d2wall, and chi(nuTilde), so dR[meanflow]/dw[nuTilde] and
    dR[nuTilde]/dw[meanflow] are both nontrivial and worth locking down
    explicitly, same spirit as the SA-GR coupling-block check.
    """
    rtol, atol = getTol(**kwargs)
    nw, blocks = getStateBlocks(CFDSolver)

    wDotFull = CFDSolver.getStatePerturbation(seed)

    extraArgs = {}
    if mode is not None:
        extraArgs = {"mode": mode, "h": h}

    for colName, colOffsets in blocks.items():
        wDot = maskStateVector(wDotFull, nw, colOffsets)
        if mode == "CS":
            # see reg_sagr.py's identical note: re-seat the real state between CS
            # columns to avoid imaginary-buffer contamination from the prior column
            CFDSolver.setStates(numpy.real(CFDSolver.getStates()))
        resDot = CFDSolver.computeJacobianVectorProductFwd(wDot=wDot, residualDeriv=True, **extraArgs)

        for rowName, rowOffsets in blocks.items():
            rowDot = maskStateVector(resDot, nw, rowOffsets)
            key = "||dR[%s]/dw[%s] * eDot||" % (rowName, colName)
            handler.root_print(key)
            handler.par_add_norm(key, rowDot, rtol=rtol, atol=atol)


def assert_bwdfast_blocks_allclose(CFDSolver, seed=314, atol=1e-16):
    """Row-block-seeded consistency between the full reverse mode (_b) and the
    reverse-fast state-only mode (_fast_b), split meanflow vs. nuTilde."""
    nw, blocks = getStateBlocks(CFDSolver)
    dwBarFull = CFDSolver.getStatePerturbation(seed)

    for rowName, rowOffsets in blocks.items():
        dwBar = maskStateVector(dwBarFull, nw, rowOffsets)

        wBar = CFDSolver.computeJacobianVectorProductBwd(resBar=dwBar, wDeriv=True)
        wBarFast = CFDSolver.computeJacobianVectorProductBwdFast(resBar=dwBar)

        numpy.testing.assert_allclose(
            wBar, wBarFast, atol=atol, err_msg="BWD vs BWDFast, resBar seeded on %s rows" % rowName
        )


# --------------------------------------------------------------------------
# Direct term-comparison helpers (SA-BCM modified vs. plain SA unmodified) --
# unchanged from the original reg_bcm.py.
# --------------------------------------------------------------------------
def assert_bcm_vs_plain_sa_wdot_allclose(handler, CFDSolverBCM, CFDSolverSA, seed=314, mode=None, h=None, **kwargs):
    """Forward products dR/dw * wDot for the SA-BCM case and the matched plain-SA case, recorded
    side by side under distinct keys. This does not assert they're equal (they shouldn't be --
    SA-BCM changes the production term) -- it records both so a reviewer can see exactly which
    residual rows differ once use_SABCM is turned on, and by how much.
    """
    rtol, atol = getTol(**kwargs)
    extraArgs = {"mode": mode, "h": h} if mode is not None else {}

    wDotBCM = CFDSolverBCM.getStatePerturbation(seed)
    resDotBCM = CFDSolverBCM.computeJacobianVectorProductFwd(wDot=wDotBCM, residualDeriv=True, **extraArgs)
    handler.root_print("||dR/dw * wDot|| (use_SABCM=True)")
    handler.par_add_norm("||dR/dw * wDot|| (use_SABCM=True)", resDotBCM, rtol=rtol, atol=atol)

    wDotSA = CFDSolverSA.getStatePerturbation(seed)
    resDotSA = CFDSolverSA.computeJacobianVectorProductFwd(wDot=wDotSA, residualDeriv=True, **extraArgs)
    handler.root_print("||dR/dw * wDot|| (use_SABCM=False, plain SA)")
    handler.par_add_norm("||dR/dw * wDot|| (use_SABCM=False, plain SA)", resDotSA, rtol=rtol, atol=atol)


def assert_bcm_xdvdot_allclose(handler, CFDSolver, ap, seed=1.0, mode=None, h=None, **kwargs):
    """Forward products dR/d(aeroDV) for the SA-BCM-relevant freestream DVs (bcmAeroDVs) --
    targets exactly the path that feeds Re_theta/Re_theta_c/term2 (see module docstring)."""
    rtol, atol = getTol(**kwargs)
    extraArgs = {"mode": mode, "h": h} if mode is not None else {}

    for aeroDV in ap.DVs.values():
        key = aeroDV.key
        if key not in bcmAeroDVs:
            continue
        resDot = CFDSolver.computeJacobianVectorProductFwd(xDvDot={key: seed}, residualDeriv=True, **extraArgs)
        name = "||dR/d%s||" % key
        handler.root_print(name)
        handler.par_add_norm(name, resDot, rtol=rtol, atol=atol)
