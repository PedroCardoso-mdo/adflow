"""Shared fixtures and helpers for the transition-model regression tests.

Plays the role of reg_default_options.py, reg_aeroproblems.py and reg_test_utils.py for the three
transition models:

    sagr    turbulenceModel = "SA-noft2-Gamma-Retheta"   states: 5 mean flow + nuTilde, gamma, reThetat
    sgamma  turbulenceModel = "SA-noft2-Gamma"           states: 5 mean flow + nuTilde, gamma
    sabcm   useSABCM = True on the standard SA model     states: 5 mean flow + nuTilde

The flattened state and residual vectors are variable-fastest, w(i,j,k,1:nw): python offset
``iCell*nw + l`` holds variable l (0-based).
"""

# built-ins
import os
import copy
from collections import OrderedDict

# external
import numpy
from baseclasses import AeroProblem
from baseclasses.testing import getTol

from reg_aeroproblems import ap_tutorial_wing

baseDir = os.path.dirname(os.path.abspath(__file__))
inputDir = os.path.join(baseDir, "../../input_files")

# --------------------------------------------------------------------------
# Input files
# --------------------------------------------------------------------------
# Standard tutorial-wing RANS grid and FFD (from input_files/get-input-files.sh)
tutWingGridFile = os.path.join(inputDir, "mdo_tutorial_rans_scalar_jst.cgns")
tutWingFFDFile = os.path.join(inputDir, "mdo_tutorial_ffd.fmt")
# Converged double-precision states written by dev/generate_<model>_restart.py. The derivative tests
# linearize about these states, so the options and flow conditions below must match the ones the
# states were converged with.
sagrRestartFile = os.path.join(inputDir, "mdo_tutorial_sagr_dp.cgns")
sgammaRestartFile = os.path.join(inputDir, "mdo_tutorial_sgamma_dp.cgns")

# --------------------------------------------------------------------------
# Aero problems
# --------------------------------------------------------------------------
# Tutorial wing at Mach 0.15 (Re_chord ~ 4.2e6): natural transition sits around mid-chord, so the
# laminar, transitional and turbulent branches of the models are all active.
_tutWingTransitionFlow = dict(
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

ap_sagr_tut_wing = AeroProblem(name="mdo_tutorial_sagr", **_tutWingTransitionFlow)
ap_sgamma_tut_wing = AeroProblem(name="mdo_tutorial_sgamma", **_tutWingTransitionFlow)

# mach, P and T reach the transition residuals through uInf and muInf (vorticity limiter) and
# through the farfield values of gamma and reThetat.
transitionAeroDVs = ["alpha", "mach", "P", "T"]

# --------------------------------------------------------------------------
# Solver options
# --------------------------------------------------------------------------
# Options common to the two transport-equation models. The solver path is only used by the tests
# that converge a solution (test_adjoint_*, dev/generate_*_restart.py); the fixed-state tests only
# evaluate residuals and Jacobian-vector products.
_gammaModelOptions = {
    "equationtype": "RANS",
    # ---- model requirements ----
    # pyADflow disables the blockette residual for the transition models; explicit here
    "useblockettes": False,
    # the gamma correlations need the resolved near-wall profile
    "usewallfunctions": False,
    # SA-noft2: gamma replaces the ft2 laminar-suppression term
    "useft2sa": False,
    # ---- freestream ----
    # Tu_inf (fraction): sets the farfield reThetat and enters the onset correlations
    "turbintensityinf": 0.0025,
    # only scales the dissipation of the primal solve; value the restart states were converged with
    "acousticscalefactor": 0.15,
    # ---- transition model ----
    # first-order upwind convection of gamma and reThetat (Piotrowski & Zingg 2020, section IV.A)
    "transitionfirstorderupwind": True,
    # crossflow source on the reThetat equation; off (the default). The sagr_crossflow_tut_wing
    # parameter sets switch it on (the term is active on this swept wing).
    "transitioncrossflow": False,
    # surface roughness height [m] in the crossflow correlation
    "transitionroughnessheight": 3.3e-6,
    # source-term time-step restriction (P&Z Eq. 59) and its limit; solver side only
    "transitionsrcdtrestrict": True,
    "transitionsrcdtlimit": 0.9,
    "srcdtdeactivateiters": 5,
    # update damping of gamma and reThetat in DADI (P&Z Algorithm 2); solver side only
    "transitiondamptheta": 0.99,
    "transitiondampmaxiter": 10000,
    # drop the production term from the SA row of the approximate Jacobian; preconditioner only
    "transitionuseapproxsa": True,
    # reference length of the vorticity limiter (P&Z Eqs. 52-53); negative = chordRef
    "transitionreflength": -1.0,
    # DADI coupling of the turbulence equations: one fully coupled block
    "turbdadicoupled": "full",
    # keep the source time-step restriction and the update damping active inside NK
    "transitionnk": True,
    # ---- startup smoother ----
    # multigrid is not supported for the transition models
    "mgcycle": "sg",
    "smoother": "DADI",
    "cfl": 1.7,
    "cflcoarse": 1.0,
    "resaveraging": "alternate",
    "nsubiterturb": 3,
    "ncycles": 10000,
    "monitorvariables": ["cpu", "resrho", "resturb", "cl", "cd", "totalr", "scaledtotalr"],
    # ---- ANK / NK ----
    "usenksolver": True,
    "useanksolver": True,
    "ankswitchtol": 1000.0,
    # turbulence solved by DADI inside the uncoupled ANK phase
    "ankuseturbdadi": True,
    # relaxed line-search and physicality limits; with the defaults the coupled phase stalls
    "ankunsteadylstol": 2.0,
    "ankphysicallstol": 0.8,
    "ankphysicallstolturb": 0.99,
    "anknsubiterturb": 1,
    "ankcfl0": 5.0,
    "ankuseapproxsa": False,
    # the AD preconditioner is required; the FD-coloured one does not converge these models
    "ankadpc": True,
    "nkadpc": True,
    "nksubspacesize": 300,
    "ankchartimesteptype": "VLR",
    # ---- convergence / adjoint ----
    # states are written and read in double precision so the tests linearize about the exact state
    "solutionprecision": "double",
    "l2convergence": 1e-14,
    "l2convergencecoarse": 1e-4,
    "adjointl2convergence": 1e-14,
    # the transposed Jacobian is stiff: the default 500 iterations do not reach 1e-14
    "adjointmaxiter": 3000,
    # the adjoint includes the transition equations
    "frozenturbulence": False,
    # ---- output ----
    "writevolumesolution": True,
    "writesurfacesolution": True,
    "surfacevariables": ["cp", "mach", "cf", "vx", "vy", "vz", "cfx", "cfy", "cfz"],
}

_transitionVolumeVariables = [
    "resrho",
    "resturb",
    "intermittency",
    "fonset",
    "thetabl",
    "deltabl",
    "transdelta",
    "transvelmag",
    "duds",
    "nutsa",
    "transvortlim",
]

# SA-noft2-Gamma-Retheta: ANK -> coupled ANK at 1e-5 -> NK at 1e-5
sagrOptions = {
    **copy.deepcopy(_gammaModelOptions),
    # the restart file holds the grid and the converged solution
    "gridfile": sagrRestartFile,
    "restartfile": sagrRestartFile,
    "turbulencemodel": "SA-noft2-Gamma-Retheta",
    # residual scaling of nuTilde, gamma, reThetat in the coupled DADI block
    "turbresscale": [10000.0, 0.1, 1.0e-4],
    # coupling earlier than 1e-5 stalls
    "ankcoupledswitchtol": 1e-5,
    "anksecondordswitchtol": 1e-6,
    "ankcfllimit": 1e6,
    "nkswitchtol": 1e-5,
    "volumevariables": _transitionVolumeVariables + ["rethetatilde"],
}

# SA-noft2-Gamma: ANK -> second-order uncoupled ANK at 1e-2 -> NK at 1e-8
# Same state with the crossflow transition term (stationary crossflow, helicity based) switched on.
# Jacobian-vector products can be evaluated about any state, so no separate restart is needed.
sagrCrossflowOptions = {**copy.deepcopy(sagrOptions), "transitioncrossflow": True}

# Total derivatives with crossflow need a state converged with crossflow on. It is written by
# dev/generate_sagr_restart.py --variant crossflow; test_adjoint_sagr.py adds the parameter set only
# when the file exists.
sagrCrossflowRestartFile = os.path.join(inputDir, "mdo_tutorial_sagr_crossflow_dp.cgns")
sagrCrossflowAdjointOptions = {
    **copy.deepcopy(sagrOptions),
    "transitioncrossflow": True,
    "gridfile": sagrCrossflowRestartFile,
    "restartfile": sagrCrossflowRestartFile,
}

sgammaOptions = {
    **copy.deepcopy(_gammaModelOptions),
    "gridfile": tutWingGridFile,
    "restartfile": sgammaRestartFile,
    "turbulencemodel": "SA-noft2-Gamma",
    # residual scaling of nuTilde, gamma
    "turbresscale": [10000.0, 0.1],
    # never couple: the coupled ANK step stalls for this model
    "ankcoupledswitchtol": 1e-20,
    "anksecondordswitchtol": 1e-2,
    "ankcfllimit": 1e8,
    # NK only for the last orders
    "nkswitchtol": 1e-8,
    "volumevariables": _transitionVolumeVariables,
}

# Same state with the gradient-driven artificial viscosity on gamma switched on. Residuals and
# Jacobian-vector products can be evaluated about any state, so no separate restart is needed.
sgammaArtViscOptions = {**copy.deepcopy(sgammaOptions), "sgammaartvisc": 10.0}
# Total derivatives with the artificial viscosity need a state converged with it. It is written by
# dev/generate_sgamma_restart.py --variant artvisc; test_adjoint_sgamma.py adds the parameter set only
# when the file exists.
sgammaArtViscRestartFile = os.path.join(inputDir, "mdo_tutorial_sgamma_artvisc_dp.cgns")
sgammaArtViscAdjointOptions = {
    **copy.deepcopy(sgammaArtViscOptions),
    "gridfile": sgammaArtViscRestartFile,
    "restartfile": sgammaArtViscRestartFile,
}

# SA-BCM (algebraic intermittency on the SA model). There is no converged SA-BCM restart: the
# residuals and Jacobian-vector products are checked at two states built from the converged fully
# turbulent SA solution of the tutorial wing (flow conditions of ap_tutorial_wing).
#   "turbulent":    the SA solution itself. Term2 saturates the intermittency in the boundary
#                   layer, which exercises Term2 and the intermittency blend.
#   "transitional": the SA solution with nuTilde scaled by 0.01 and Tu = 2 %. The boundary layer is
#                   then laminar or transitional and Term1 (Re_theta against its critical value)
#                   is active as well.
sabcmOptions = {
    "gridfile": tutWingGridFile,
    "restartfile": tutWingGridFile,
    "equationType": "RANS",
    "useSABCM": True,
    "turbIntensityInf": 0.005,
}

sabcm_test_params = [
    {
        "name": "rans_sabcm_smooth_tut_wing",
        "options": {**sabcmOptions, "SABCMSmooth": True},
        "ref_file": "funcs_rans_sabcm_smooth_tut_wing.json",
        "aero_prob": copy.deepcopy(ap_tutorial_wing),
        "nuTildeScale": 1.0,
    },
    {
        "name": "rans_sabcm_original_tut_wing",
        "options": {**sabcmOptions, "SABCMSmooth": False},
        "ref_file": "funcs_rans_sabcm_original_tut_wing.json",
        "aero_prob": copy.deepcopy(ap_tutorial_wing),
        "nuTildeScale": 1.0,
    },
    {
        "name": "rans_sabcm_smooth_transitional_tut_wing",
        "options": {**sabcmOptions, "SABCMSmooth": True, "turbIntensityInf": 0.02},
        "ref_file": "funcs_rans_sabcm_smooth_transitional_tut_wing.json",
        "aero_prob": copy.deepcopy(ap_tutorial_wing),
        "nuTildeScale": 0.01,
    },
    {
        "name": "rans_sabcm_original_transitional_tut_wing",
        "options": {**sabcmOptions, "SABCMSmooth": False, "turbIntensityInf": 0.02},
        "ref_file": "funcs_rans_sabcm_original_transitional_tut_wing.json",
        "aero_prob": copy.deepcopy(ap_tutorial_wing),
        "nuTildeScale": 0.01,
    },
]


def setSABCMTestState(CFDSolver, ap, nuTildeScale):
    """Propagate the restart state and scale nuTilde to build the SA-BCM test state."""
    CFDSolver.getResidual(ap)
    if nuTildeScale != 1.0:
        nw = CFDSolver.adflow.flowvarrefstate.nw
        w = CFDSolver.getStates()
        w.reshape(-1, nw)[:, -1] *= nuTildeScale
        CFDSolver.setStates(w)
        CFDSolver.getResidual(ap)


# --------------------------------------------------------------------------
# State-vector block helpers
# --------------------------------------------------------------------------
def getStateBlocks(CFDSolver):
    """Return (nw, blocks). blocks is an ordered dict mapping a block name (meanflow, nuTilde,
    gamma, reThetat) to the list of 0-based variable offsets inside each cell's nw-slot chunk.

    nw, nt1 and nt2 are read from the compiled module, so this works on the real and the complex
    build and for the three layouts: nw = 8 (SA-noft2-Gamma-Retheta), 7 (SA-noft2-Gamma), 6 (SA).
    """
    nw = int(CFDSolver.adflow.flowvarrefstate.nw)
    nt1 = int(CFDSolver.adflow.flowvarrefstate.nt1)  # Fortran 1-based
    nt2 = int(CFDSolver.adflow.flowvarrefstate.nt2)

    turbNames = {6: ["nuTilde"], 7: ["nuTilde", "gamma"], 8: ["nuTilde", "gamma", "reThetat"]}
    if nw not in turbNames or (nt2 - nt1 + 1) != len(turbNames[nw]):
        raise ValueError(
            "expected the SA (nw=6), SA-noft2-Gamma (nw=7) or SA-noft2-Gamma-Retheta (nw=8) layout; "
            "got nw=%d, nt1=%d, nt2=%d" % (nw, nt1, nt2)
        )

    blocks = OrderedDict()
    blocks["meanflow"] = list(range(0, nt1 - 1))  # rho, rho*u, rho*v, rho*w, rho*E
    for i, name in enumerate(turbNames[nw]):
        blocks[name] = [nt1 - 1 + i]
    return nw, blocks


def maskStateVector(vec, nw, offsets):
    """Zero every entry of a flattened state-shaped vector except the given 0-based per-cell
    variable offsets."""
    masked = numpy.zeros_like(vec)
    for off in offsets:
        masked[off::nw] = vec[off::nw]
    return masked


# (column block seeded in wDot, row block seeded in dwBar) pairs of the blockwise dot-product test
couplingDotProductPairs = {
    "sagr": [
        ("gamma", "nuTilde"),  # dR[nuTilde]/d(gamma): gamma multiplies the SA production
        ("nuTilde", "gamma"),  # dR[gamma]/d(nuTilde): turbulent viscosity ratio R_T
        ("reThetat", "gamma"),  # dR[gamma]/d(reThetat): Flength and Re_theta_c correlations
        ("nuTilde", "reThetat"),  # dR[reThetat]/d(nuTilde): nonzero only with crossflow
        ("meanflow", "reThetat"),  # dR[reThetat]/d(meanflow): time scale and F_theta_t
        ("gamma", "gamma"),
        ("reThetat", "reThetat"),
    ],
    "sgamma": [
        ("gamma", "nuTilde"),  # dR[nuTilde]/d(gamma): gamma multiplies the SA production
        ("nuTilde", "gamma"),  # dR[gamma]/d(nuTilde): turbulent viscosity ratio R_T
        ("meanflow", "gamma"),  # dR[gamma]/d(meanflow): Re_v, lambda_theta, Re_theta_c
        ("gamma", "gamma"),
    ],
}


# --------------------------------------------------------------------------
# Assert helpers. Like the ones in reg_test_utils.py they record norms in the BaseRegTest handler.
# The AD (real build) and CS (complex build) test classes record the same keys, so training with
# the real build and then running the cmplx_* tests compares AD against complex step.
# --------------------------------------------------------------------------
def assert_coupling_blocks_allclose(handler, CFDSolver, seed=314, mode=None, h=None, **kwargs):
    """Forward products dR/dw * e[column block], with the result norms split by residual row block.

    This checks each block of the Jacobian separately, in particular the coupling between the SA
    and the transition equations, which a whole-vector norm would not resolve. The boundary
    condition chains are included because the turbulence BCs are applied inside the residual.
    """
    rtol, atol = getTol(**kwargs)
    nw, blocks = getStateBlocks(CFDSolver)

    wDotFull = CFDSolver.getStatePerturbation(seed)

    extraArgs = {}
    if mode is not None:
        extraArgs = {"mode": mode, "h": h}

    failed = []
    for colName, colOffsets in blocks.items():
        wDot = maskStateVector(wDotFull, nw, colOffsets)
        if mode == "CS":
            # Sequential CS products share a complex work buffer: a large-norm column leaves an
            # imaginary residue of machine-epsilon size that shows up in the structurally zero rows
            # of the next column. Re-seating the real state removes it.
            CFDSolver.setStates(numpy.real(CFDSolver.getStates()))
        resDot = CFDSolver.computeJacobianVectorProductFwd(wDot=wDot, residualDeriv=True, **extraArgs)

        for rowName, rowOffsets in blocks.items():
            rowDot = maskStateVector(resDot, nw, rowOffsets)
            key = "||dR[%s]/dw[%s] * eDot||" % (rowName, colName)
            handler.root_print(key)
            # check every block before failing, so one mismatch does not hide the others
            try:
                handler.par_add_norm(key, rowDot, rtol=rtol, atol=atol)
            except AssertionError as e:
                handler.root_print("BLOCK FAIL " + key + " :: " + " ".join(str(e).split()))
                failed.append(key)
    if failed:
        raise AssertionError("coupling blocks failed: %s" % failed)


def assert_transition_xdvdot_allclose(handler, CFDSolver, ap, seed=1.0, mode=None, h=None, **kwargs):
    """Forward products dR/d(aeroDV) for each aero DV of the AeroProblem, with the result norms
    split by residual row block.

    The gamma and reThetat rows are the targeted check of the flow-condition derivatives of the
    transition models (uInf and muInf in the vorticity limiter, farfield transition BC).
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
    """Consistency between the full reverse mode and the reverse-fast (state-only) mode, seeding
    one residual row block at a time so each transposed coupling column is checked separately."""
    nw, blocks = getStateBlocks(CFDSolver)
    dwBarFull = CFDSolver.getStatePerturbation(seed)

    for rowName, rowOffsets in blocks.items():
        dwBar = maskStateVector(dwBarFull, nw, rowOffsets)

        wBar = CFDSolver.computeJacobianVectorProductBwd(resBar=dwBar, wDeriv=True)
        wBarFast = CFDSolver.computeJacobianVectorProductBwdFast(resBar=dwBar)

        numpy.testing.assert_allclose(
            wBar, wBarFast, atol=atol, err_msg="BWD vs BWDFast, resBar seeded on %s rows" % rowName
        )


def assert_coupling_dot_products_allclose(handler, CFDSolver, pairs, seed=314, **kwargs):
    """Blockwise dot-product tests: for each (column, row) pair, seed wDot on one variable block
    and dwBar on another and require dwBar^T (dR/dw wDot) == wDot^T (dR/dw^T dwBar). Each pair
    isolates one Jacobian block in both AD modes."""
    rtol, atol = getTol(**kwargs)
    nw, blocks = getStateBlocks(CFDSolver)

    wDotFull = CFDSolver.getStatePerturbation(seed)
    dwBarFull = CFDSolver.getStatePerturbation(seed + 1)

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


# --------------------------------------------------------------------------
# FD-vs-AD comparison with a step-size sweep
# --------------------------------------------------------------------------
# The transition residuals have one-sided switches (limiters, blending functions) and span about
# 13 orders of magnitude, so no single FD step suits every cell. The check keeps a tight tolerance
# and sweeps the step: it passes if any step meets the tolerance and otherwise reports the best
# step. The complex-step classes are the step-free reference.
FD_H_SWEEP = [3e-5, 1e-5, 3e-6, 1e-6, 3e-7, 1e-7, 3e-8, 1e-8, 3e-9, 1e-9, 1e-10]


def _global_tol_ratio(fd, ref, rtol, atol, comm):
    """Worst element ratio |fd-ref| / (atol + rtol|ref|) over all MPI ranks; <= 1.0 means the
    array meets numpy.allclose(rtol, atol). Also returns max|fd-ref|."""
    diff = numpy.abs(fd - ref)
    allowed = atol + rtol * numpy.abs(ref)
    with numpy.errstate(divide="ignore", invalid="ignore"):
        ratio = numpy.where(allowed > 0, diff / allowed, numpy.where(diff > 0, numpy.inf, 0.0))
    local_ratio = float(numpy.max(ratio)) if ratio.size else 0.0
    local_absdiff = float(numpy.max(diff)) if diff.size else 0.0
    if comm is not None:
        from mpi4py import MPI

        local_ratio = comm.allreduce(local_ratio, op=MPI.MAX)
        local_absdiff = comm.allreduce(local_absdiff, op=MPI.MAX)
    return local_ratio, local_absdiff


def assert_fd_allclose_hsweep(
    computeFD, refDot, rtol=0.0, atol=0.0, hPrimary=1e-8, hSweep=None, comm=None, root_print=None, err_msg=""
):
    """AD-vs-FD array comparison with a step-size sweep.

    computeFD(h) returns the FD directional derivative to compare against the AD array refDot.
    hPrimary is tried first, then the steps in hSweep; the check passes as soon as one step meets
    (rtol, atol). If none does, an AssertionError reports the best step, its worst tolerance ratio
    and max|delta|.
    """
    if hSweep is None:
        hSweep = FD_H_SWEEP
    if root_print is None:
        root_print = print

    ratio, absdiff = _global_tol_ratio(computeFD(hPrimary), refDot, rtol, atol, comm)
    if ratio <= 1.0:
        return

    best_h, best_ratio, best_absdiff = hPrimary, ratio, absdiff
    for h in hSweep:
        if h == hPrimary:
            continue
        ratio, absdiff = _global_tol_ratio(computeFD(h), refDot, rtol, atol, comm)
        if ratio < best_ratio:
            best_h, best_ratio, best_absdiff = h, ratio, absdiff
        if ratio <= 1.0:
            root_print("  [h-sweep] %s: PASS at h=%.0e (failed at primary h=%.0e)" % (err_msg, h, hPrimary))
            return

    raise AssertionError(
        "%s: FD-vs-AD failed at every step of the sweep. Best h=%.0e -> worst tol-ratio=%.2f (>1), "
        "max|delta|=%.3e. Use the complex-step test for a step-free check."
        % (err_msg, best_h, best_ratio, best_absdiff)
    )
