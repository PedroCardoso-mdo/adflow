"""Validated ADflow option sets for the three transition models of this fork.

    GR      turbulenceModel = "SA-noft2-Gamma-Retheta"   (main model; crossflow by transitionCrossflow)
    SA-gamma  turbulenceModel = "SA-noft2-Gamma"
    SA-BCM  turbulenceModel = "SA" + useSABCM = True     (SABCMSmooth True = smooth form, False = original)
    SA      plain SA (no ft2), the fully turbulent reference

Every dict below is copied from the job script / option dump of a run that actually
converged; the comment next to each block names that run. Where a document and the
job script disagreed, the job script won (the cases are listed in README.md).

How the pieces combine (see ``options()`` at the bottom):

    COMMON  +  DISSIPATION[...]  +  model block  +  RECIPES[...]  +  case-specific keys

Notation used in the comments:
    ANK    approximate Newton-Krylov, turbulence segregated (DADI)
    CANK   ANK with flow and turbulence in one coupled Krylov system   (below ANKCoupledSwitchTol)
    CSANK  coupled ANK with the exact (second-order) matrix-free Jacobian (below ANKSecondOrdSwitchTol)
    NK     Newton-Krylov                                               (below NKSwitchTol)
All switch tolerances are relative to the free-stream residual.

The campaign folders quoted as provenance (e.g. 24_minmodel_2d, 13_nlf0415_swept) belong to the
author's run tree, which is not in git; the reduced results ship with the showcase package.
In the comments below, `L/` is the root of that run tree and `H/` its twin on the cluster.
Inventories with every number: showcase/_manifest_2d.md, showcase/_manifest_3d.md.
"""

# --------------------------------------------------------------------------------------
# 1. Common block (identical for every model and case)
# --------------------------------------------------------------------------------------
COMMON = {
    "equationType": "RANS",
    "useft2SA": False,                # "SA-noft2": the ft2 trip term is off for every model
    "useWallFunctions": False,
    "L2Convergence": 1e-8,            # 8 orders on totalR: the standing target of the validated runs
    "printIterations": True,
    "printTiming": True,
    "writeSurfaceSolution": True,
    "writeVolumeSolution": False,     # the example scripts switch it on with --writeVolume
    "writeSolutionEachIter": False,   # never True: it writes every iteration
    "solutionPrecision": "double",    # single precision truncates the transition front on restart
    "surfaceVariables": ["cp", "cf", "cfx", "cfy", "cfz", "yplus", "mach", "blank"],
}

# --------------------------------------------------------------------------------------
# 2. Dissipation schemes
# --------------------------------------------------------------------------------------
DISSIPATION = {
    # Scalar JST with the fourth-difference coefficient lowered to 0.01.
    # Campaign 24 (L/24_minmodel_2d/config/setup.py:70-82): every 2D article run.
    "scalar": {
        "discretization": "central plus scalar dissipation",
        "vis2": 0.25,
        "vis4": 0.01,
        "dissipationScalingExponent": 0.67,
        "acousticScaleFactor": 1.0,   # = Mach was tried and dropped (3x more cycles, limit cycles)
    },
    # Matrix dissipation as in Piotrowski & Zingg 2020: dis4 0.04, Vl = 0 (epsShear), Vn = 0.25 (epsAcoustic).
    # Campaign 12 job 1827893 (jobs/job_2dfoils_final.sh:51-56), campaign 13 job 1826778, campaign 24-C.
    # Needed for crossflow: scalar dissipation under-predicts the helicity (swept wing: no crossflow front).
    "matrix": {
        "discretization": "central plus matrix dissipation",
        "vis2": 0.25,
        "vis4": 0.04,
        "epsShear": 0.0,
        "epsAcoustic": 0.25,
        "dissipationScalingExponent": 0.67,
        "acousticScaleFactor": 1.0,
    },
    # ADflow default scalar JST (vis4 0.0156): the AR5 refinement study (campaign 11) and the
    # SA-BCM article ran with it. Kept only to reproduce those numbers.
    "scalar_default": {
        "discretization": "central plus scalar dissipation",
        "vis2": 0.25,
        "vis4": 0.0156,
    },
}

# --------------------------------------------------------------------------------------
# 3. Model blocks
# --------------------------------------------------------------------------------------
# Free-stream nu_tilde (eddyVisInfRatio) is a boundary condition of each model, not a solver knob.
# L/24_minmodel_2d/config/setup.py:104-112 (sweep 24-D: all models sit on a plateau below ~1e-6;
# the ADflow default 9e-3 destroys the laminar region: NLF0416 alpha 4, GR +17 % cd).
NUINF = {
    "gr": 2.79e-7,     # chi_inf = 0.1 of Piotrowski & Zingg; at 1e-7 S809 falls into a bubble limit cycle
    "sgamma": 1e-7,    # campaign-24 setup; the L0 polars of the figures ran at 2.79e-7 (pass --nuinf)
    "bcm": 1e-9,       # smooth SA-BCM
    "bcmo": 1e-9,      # original SA-BCM
    "sa": 9e-3,        # fully turbulent: ADflow default
}

# Options of the transition models that pyADflow now sets by itself when turbulenceModel contains
# "Gamma" and the user does not pass them (adflow/pyADflow.py, "Model-gated defaults").
# Written out here so that the recipe is complete on its own.
TRANSITION_PYADFLOW_DEFAULTS = {
    "globalPreconditioner": "field split",   # adjoint PC: 2D -29 %, 3D 5.1x vs additive Schwarz
    "ILUFill": 3,
    "ASMOverlap": 3,
    "adjointSolver": "LGMRES",
    "adjointSubspaceSize": 400,
    "NKLSRelax": True,                       # relaxed NK line search (only acts in the NK phase)
}

# Shared by GR and SA-gamma (L/24_minmodel_2d/config/recipes.py:49-56, L/13_nlf0415_swept/run_swept.py:295-301).
_TRANSITION = {
    "transitionFirstOrderUpwind": True,   # paper IV.A: first-order upwind for gamma and Re_theta_t
    "transitionUseApproxSA": True,
    "transitionSrcDtRestrict": True,
    "srcDtDeactivateIters": 5,
    "TurbDADICoupled": "full",            # 3x3 (GR) / 2x2 (SA-gamma) block DADI
    "transitionCrossflow": False,         # the crossflow cases switch it on themselves
}

MODEL_LABEL = {
    "gr": "SA-noft2-Gamma-Retheta (GR)",
    "sgamma": "SA-noft2-Gamma (SA-gamma)",
    "bcm": "SA-BCM (smooth)",
    "bcmo": "SA-BCM (original)",
    "sa": "SA (fully turbulent)",
}
MODELS = tuple(MODEL_LABEL)


def modelBlock(model, tuPercent, nuinf=None):
    """Options that select the model and its boundary conditions. Tu in percent."""
    if model == "gr":
        b = dict(TRANSITION_PYADFLOW_DEFAULTS, **_TRANSITION)
        b.update(turbulenceModel="SA-noft2-Gamma-Retheta", volumeVariables=["intermittency", "fonset", "eddyratio"])
    elif model == "sgamma":
        b = dict(TRANSITION_PYADFLOW_DEFAULTS, **_TRANSITION)
        b.update(
            turbulenceModel="SA-noft2-Gamma",
            sgammaCoupleDestruction=True,     # L/24_minmodel_2d/config/setup.py:115
            volumeVariables=["intermittency", "fonset", "eddyratio"],
        )
    elif model in ("bcm", "bcmo", "sa"):
        b = {"turbulenceModel": "SA", "useSABCM": model != "sa", "volumeVariables": ["eddyratio"]}
        if model != "sa":
            # Constants of the SA-BCM article (L/24_minmodel_2d/config/recipes.py:176-179); all are the
            # pyADflow defaults. With SABCMSmooth False (original model) Rho / Tanh* are unused.
            b.update(
                SABCMSmooth=(model == "bcm"),
                SABCMChi1=0.002,
                SABCMChi2=0.02,
                SABCMTanhCenter=0.5,
                SABCMTanhWidth=0.08,
                SABCMRho=50.0,
            )
    else:
        raise ValueError(f"model must be one of {MODELS}, got {model!r}")
    b["turbIntensityInf"] = tuPercent / 100.0     # fraction; SA ignores it
    b["eddyVisInfRatio"] = NUINF[model] if nuinf is None else float(nuinf)
    return b


# --------------------------------------------------------------------------------------
# 4. Solver recipes
# --------------------------------------------------------------------------------------
# Staged-ANK base shared by GR and SA-gamma in 2D (L/24_minmodel_2d/config/recipes.py:59-77,
# first proven in L/08_optimization/2d_tu05_L0/solver_setup.py:68-99, job 1927080: 126/126 solves).
_ANK_BASE = {
    "smoother": "DADI",
    "MGCycle": "sg",                  # multigrid is not supported for the transition models
    "useANKSolver": True,
    "ANKCharTimeStepType": "VLR",
    "ANKUseTurbDADI": True,           # segregated turbulence in the first (ANK) phase
    "ANKUnsteadyLSTol": 2.0,
    "ANKPhysicalLSTol": 0.8,
    "ANKADPC": True,
    "NKADPC": True,
    "ANKCFLLimit": 1.0e8,             # the old 1e6 ceiling pinned the CFL and stalled CANK (job 1825801)
    "useNKSolver": False,             # no NK in the 2D article runs (L/.../recipes.py:212-215)
    # adjoint
    "ADPC": True,
    "adjointL2Convergence": 1e-9,
    "adjointMaxIter": 10000,
    "adjointMonitorStep": 10,
    "innerPreconIts": 4,
    "outerPreconIts": 10,
    "skipAfterFailedAdjoint": False,
}
_MON_FULL = ["resrho", "cl", "cd", "cdp", "cdv", "resturb", "totalr"]

RECIPES = {}

# ---- GR, 2D, scalar dissipation -------------------------------------------------------
# ANK -> CANK at 1e-5 -> CSANK at 1e-6, no NK.
# Campaign 24 "opt" recipe (L/24_minmodel_2d/config/recipes.py:83-88). Runs: NACA 0012 mesh study
# jobs 1962375-79 (6/6 levels at 1e-8, L0: 299 s on 32 ranks), single-point optimisations (job 1988112).
RECIPES["gr_2d_scalar"] = dict(
    _ANK_BASE,
    ANKCoupledSwitchTol=1e-5,
    ANKSecondOrdSwitchTol=1e-6,
    turbResScale=[1e4, 0.1, 1e-4],    # nu_tilde, gamma, Re_theta_t (the state spans 13 orders)
    nCycles=60000,                    # counts "Iter Tot"; 30000 capped 3/7 evaluations (job 1961106)
    monitorVariables=_MON_FULL,
)

# GR polars on L0, scalar: ANK -> CANK at 1e-2 -> CSANK at 1e-5, no NK.
# Campaign 24-B, jobs 1961834-39 and 1961894 (result.json: purpose "polar_old", useNKSolver False):
# NLF0416 13/13 and S809 10/10 at 1e-8, 32 ranks, 570-900 s per point.
# S809 alpha 1 and 2 needed ANKSecondOrdSwitchTol 1e-4 (fig03_polares/RUNS.md).
RECIPES["gr_2d_scalar_polar"] = dict(
    _ANK_BASE,
    ANKCoupledSwitchTol=1e-2,
    ANKSecondOrdSwitchTol=1e-5,
    turbResScale=[1e4, 0.1, 1e-4],
    nCycles=40000,
    monitorVariables=_MON_FULL,
)

# ---- GR, 2D, matrix dissipation -------------------------------------------------------
# ANK -> CANK at 1e-2 -> CSANK at 1e-4, no NK.
# Campaign 12-A jobs 1827893 / 1828015 (jobs/job_2dfoils_final.sh:51-56; 21/21 points at ~1e-8) and
# campaign 24-C jobs 2201434 / 2201491 (job_s809_gr_mtx_L0.sh:24-25; NLF 13/13, S809 10/10).
# With the stock tolerances (1e-5 / 1e-6) CSANK never engages and S809 stalls at ~2e-5.
RECIPES["gr_2d_matrix"] = dict(
    _ANK_BASE,
    ANKCoupledSwitchTol=1e-2,
    ANKSecondOrdSwitchTol=1e-4,
    turbResScale=[1e4, 0.1, 1e-4],
    nCycles=80000,
    monitorVariables=_MON_FULL,
)

# GR, 2D, matrix dissipation, optimisation-type solves (cold start on a deformed, laminar shape).
# = gr_2d_scalar with ONE change found on 2026-11-12: ANKSecondOrdSwitchTol 1e-6 -> 1e-4, so the
# exact matrix-free Jacobian is on before the 1-2e-5 level where the laminar-bubble stall sits.
# Jobs 2219449 (lanes B_so1e4_cl05 / D_so1e4_cl07) and 2219697: the multipoint shape that never
# converged from free stream (stall 1-2e-5, jobs 2215275, 2217571) reaches 1e-8; cd 47.34 counts at
# Tu 0.3 %, cl 0.5 (29 min on 8 ranks). In use: L/24_minmodel_2d/07_multipoint/job_mp7_mtx_gr_so.sh.
# Note the order: second-order ANK starts at 1e-4, the coupling at 1e-5.
RECIPES["gr_2d_matrix_opt"] = dict(
    _ANK_BASE,
    ANKCoupledSwitchTol=1e-5,
    ANKSecondOrdSwitchTol=1e-4,
    ANKMaxIter=80,                    # -7..-8 % wall time (recipe_ab, job 2197279)
    turbResScale=[1e4, 0.1, 1e-4],
    nCycles=30000,
    adjointL2Convergence=1e-8,        # job_mp7_mtx*.sh: --diag adjointL2Convergence=1e-8
    monitorVariables=_MON_FULL,
)

# ---- GR, 3D with crossflow: NLF(2)-0415 infinite swept wing ----------------------------
# ANK -> CANK at 1e-2 -> CSANK at 1e-5 -> NK at 1e-6. Job 1826778
# (L/13_nlf0415_swept/jobs/job_swept_L0_reynolds_fix.sh:84-95, run_swept.py:295-324):
# 7/7 Reynolds numbers at 1e-8 in 495-634 iterations, 9 ranks, 36-72 min each.
# The fix that made it converge was ANKCFLLimit 1e8 (job 1825801: 534 its instead of 9282 unconverged).
RECIPES["gr_3d_crossflow_swept"] = {
    "smoother": "DADI",
    "MGCycle": "sg",
    "useANKSolver": True,
    "ANKCharTimeStepType": "VLR",
    "ANKUseTurbDADI": True,
    "ANKCoupledSwitchTol": 1e-2,
    "ANKSecondOrdSwitchTol": 1e-5,
    "ANKUnsteadyLSTol": 2.0,
    "ANKPhysicalLSTol": 0.8,
    "ANKADPC": True,
    "ANKCFLLimit": 1.0e8,
    "useNKSolver": True,
    "NKSwitchTol": 1e-6,
    "NKADPC": True,
    "NKSubspaceSize": 400,
    "NKViscPC": True,
    "NKUseEW": False,
    "NKLinearSolveTol": 0.01,
    "NKLS": "cubic",
    "NKLSRelax": False,               # as run: these binaries predate NKLSRelax (original line search)
    "turbResScale": [1e4, 0.1, 1e-4],
    "nCycles": 60000,
    "monitorVariables": ["resrho", "cl", "cd", "resturb", "totalr"],
}

# Same wing with scalar dissipation (vis4 0.01): the counter-example. It converges ~8 orders in
# 7-14 min but shows NO crossflow transition (x_tr ~0.83 at every Re).
# L/13_nlf0415_swept/jobs/job_swept_L0_resweep_ladder.sh:33-48: CANK 1e-3, CSANK 1e-5, NK 1e-6,
# ANKCFLLimit at the runner default 1e6. (Those runs used reynoldsLength sqrt(2), i.e. a lower Re.)
RECIPES["gr_3d_crossflow_swept_scalar"] = dict(
    RECIPES["gr_3d_crossflow_swept"],
    ANKCoupledSwitchTol=1e-3,
    ANKCFLLimit=1.0e6,
    nCycles=150000,
)

# ---- GR, 3D with crossflow: sickle wing on the paper authors' grid (run v14) ------------
# ANK -> CANK at 1e-3 -> CSANK at 1e-6, no NK, per-node damping, CFL floor cap 1e3.
# L/07_sickle_wing/mesh_paper_authors/job_sickle_paper_v14_cfon.sh:118-144 (64 ranks, 17.7 h,
# 6426 iterations, 6.25 orders: the 1e-8 target was missed by ~3x).
# The per-node ("Algorithm 2") damping has no switch of its own any more: it runs exactly when
# ANKTransitionGlobalLambda is False.
# Earlier project notes quote "CANK 1e-2 -> CSANK 1e-4" for this wing; no finished sickle run used that.
RECIPES["gr_3d_crossflow_sickle"] = {
    "smoother": "DADI",
    "MGCycle": "sg",
    "useANKSolver": True,
    "ANKCharTimeStepType": "VLR",
    "ANKUseTurbDADI": True,
    "ANKCoupledSwitchTol": 1e-3,
    "ANKSecondOrdSwitchTol": 1e-6,
    "ANKUnsteadyLSTol": 2.0,
    "ANKPhysicalLSTol": 0.8,
    "ANKADPC": True,
    "NKADPC": True,
    "ANKCFLLimit": 1.0e8,
    "ANKCFLMinCap": 1.0e3,
    "ANKTransitionGlobalLambda": False,
    "useNKSolver": False,
    "turbResScale": [1e4, 0.1, 1e-4],
    "nCycles": 100000,
    "monitorVariables": ["resrho", "cl", "cd", "resturb", "totalr"],
}

# ---- GR, 3D wing without crossflow: AR5 wing + winglet (campaign 11) --------------------
# First leg: ANK -> CANK at 1e-4 -> CSANK at 3e-5 -> NK at 1e-6. Second leg: restart from the
# last volume with useNKSolver False (CSANK closes the last order in 2-27 iterations).
# L/11_ar5_corrected_foil/batch_jobs/ar5cf.sh:77-80 + common/run_ar5cf.py:211-252;
# jobs 1810888-92 (leg 1) and 1810916-25 (leg 2): L4..L1 at 1e-8, L0 at 7.6 orders.
# ILU(3) stalls NK on this family and ILU(1) costs 18 % more iterations: fill 2 on both.
RECIPES["gr_3d_wing_ar5"] = {
    "smoother": "DADI",
    "MGCycle": "sg",
    "useANKSolver": True,
    "ANKCharTimeStepType": "VLR",
    "ANKUseTurbDADI": True,
    "ANKCoupledSwitchTol": 1e-4,
    "ANKSecondOrdSwitchTol": 3e-5,
    "ANKUnsteadyLSTol": 2.0,
    "ANKPhysicalLSTol": 0.8,
    "ANKADPC": True,
    "ANKPCILUFill": 2,
    "ANKCFLLimit": 1.0e8,
    "useNKSolver": True,              # leg 2: False, with restartFile
    "NKSwitchTol": 1e-6,
    "NKADPC": True,
    "NKPCILUFill": 2,
    "NKLSRelax": False,               # as run: these binaries predate NKLSRelax (original line search)
    "turbResScale": [1e4, 0.1, 1e-4],
    "blockSplitting": True,           # 5 blocks: must be split to fill the ranks
    "loadImbalance": 0.1,
    "nCycles": 60000,
    "monitorVariables": ["resrho", "cl", "cd", "resturb", "totalr"],
}

# ---- AR5 wing, matrix dissipation, 397-wrap mesh c376 (13.68 M cells): the three models ----
# GR: ANK -> CANK at 1e-2 -> CSANK at 1e-4, no NK, CFL limit 1e8, source time-step restriction
# kept for 200 iterations, target 1e-7 (no mesh of this family reaches 1e-8 unattended).
# H/08_optimization/3d_c376/{job_c376_mtx_cross.sh,run_c376_mtx_cross.py}, job 1849168 (2 nodes,
# 128 ranks, restart from a GR solution): 7.7 orders with crossflow off, 8.7 with crossflow on.
RECIPES["gr_3d_wing_matrix"] = dict(
    RECIPES["gr_3d_wing_ar5"],
    ANKCoupledSwitchTol=1e-2,
    ANKSecondOrdSwitchTol=1e-4,
    useNKSolver=False,
    srcDtDeactivateIters=200,
    L2Convergence=1e-7,
    nCycles=60000,
)

# SA-gamma on the same mesh: NOT a validated recipe (the model had never run on this wing when
# this was written). The GR 3D base with the 2D SA-gamma switching: second order at 1e-2, never
# coupled (CANK pins the step for this model in 2D, job 1942269).
RECIPES["sgamma_3d_wing_matrix"] = dict(
    RECIPES["gr_3d_wing_matrix"],
    ANKCoupledSwitchTol=1e-20,
    ANKSecondOrdSwitchTol=1e-2,
    turbResScale=[1e4, 0.1],
    globalPreconditioner="additive Schwarz",
    adjointSolver="GMRES",
)

# GR 3D wing optimisation on the 397-wrap mesh c376 (13.68 M cells) -- an ATTEMPT, not finished.
# H/08_optimization/3d_c376/{vars_master.py,solver_setup.py}. ANK -> CANK at 1e-2 -> CSANK at 1e-3,
# no NK, 7 orders (no mesh of this family reaches 8 unattended).
# The file defaults are 1e-4 / 3e-5 and the first job (1830863) ran them; every later job
# (1844851, 1931137, 1939946) was submitted with ANKCOUP=1e-2 ANKSO=1e-3 (option dump in the logs).
# Adjoint: additive Schwarz + GMRES, because field split does not fit 4 nodes on this mesh
# (jobs 1928936/37, 1929286/87 OOM); tolerance 1e-3 read off a psi-history study (job 1825800).
RECIPES["gr_3d_wing_ar5_opt"] = dict(
    RECIPES["gr_3d_wing_ar5"],
    ANKCoupledSwitchTol=1e-2,
    ANKSecondOrdSwitchTol=1e-3,
    useNKSolver=False,
    NKSwitchTol=1e-9,                 # unreachable on purpose
    L2Convergence=1e-7,
    nCycles=30000,
    monitorVariables=["resrho", "cl", "cd", "cdp", "cdv", "yplus", "resturb", "totalr"],
    ADPC=True,
    adjointL2Convergence=1e-3,
    adjointMaxIter=3000,
    adjointMonitorStep=10,
    adjointSolver="GMRES",
    globalPreconditioner="additive Schwarz",
    adjointSubspaceSize=400,
    ILUFill=3,
    ASMOverlap=3,
    innerPreconIts=4,
    outerPreconIts=10,
    skipAfterFailedAdjoint=False,
)

# ---- SA-gamma, 2D -----------------------------------------------------------------------
# ANK -> second-order ANK at 1e-2, never coupled (CANK pins the step, job 1942269), no NK.
# L/24_minmodel_2d/config/recipes.py:100-105. Runs: NACA mesh study jobs 1962375-79 (6/6 at 1e-8),
# polars jobs 1961969-74 (16-17 of 23 at 1e-8; the others at 3.2-6 orders), matrix polars 2201621/29.
# Adjoint: additive Schwarz + GMRES (field split needed 7-10x the GR iterations for this model).
RECIPES["sgamma_2d"] = dict(
    _ANK_BASE,
    ANKCoupledSwitchTol=1e-20,
    ANKSecondOrdSwitchTol=1e-2,
    turbResScale=[1e4, 0.1],          # nu_tilde, gamma
    nCycles=40000,
    monitorVariables=_MON_FULL,
    globalPreconditioner="additive Schwarz",
    adjointSolver="GMRES",
)

# SA-gamma tuned (24-I, jobs 2197279 / 2199671: -44 % wall time on the NACA at Tu 0.3 %, cl 0.5):
# couple at 1e-5, second order at 1e-4, 3 outer adjoint PC iterations (adjoint cd 457 -> 209 s).
# The matrix optimisations (job 2203553, job_mp7_mtx.sh:36-37) also add the model option
# sgammaArtVisc = 5.0 (gamma artificial viscosity); it changes cd by -1.4..-3.7 % and is NOT part
# of this dict -- pass it explicitly if you want that run.
RECIPES["sgamma_2d_tuned"] = dict(
    RECIPES["sgamma_2d"],
    ANKCoupledSwitchTol=1e-5,
    ANKSecondOrdSwitchTol=1e-4,
    outerPreconIts=3,
    adjointL2Convergence=1e-8,        # as run: job 2203553 and job_mp7_mtx*.sh
)

# ---- SA-BCM (smooth and original) and plain SA, 2D ---------------------------------------
# ANK (AD preconditioner) -> second-order ANK at 1e-3, no NK (NK stalls or lands on a wrong state:
# cd 86 instead of 65.8), 3 turbulence sub-iterations.
# SA-BCM article recipe (SA_BCM_ARTICLE/2d_opt, 242 solves) as restated in
# L/24_minmodel_2d/config/recipes.py:109-115. Runs: mesh study job 1962377, polars 1961957-62 (smooth,
# 23/23), 2168686 (original), matrix polars 2201620/28 (original, 23/23), SA polars 1962365-70.
# Monitors are residuals only: with cl/cd in the monitor the old SA-BCM code produced a start-up NaN.
RECIPES["bcm_2d"] = {
    "smoother": "DADI",
    "MGCycle": "sg",
    "useANKSolver": True,
    "ANKADPC": True,
    "ANKSecondOrdSwitchTol": 1e-3,
    "useNKSolver": False,
    "nSubiterTurb": 3,
    "useBlockettes": False,
    "infChangeCorrection": True,
    "L2ConvergenceCoarse": 1e-4,
    "nCycles": 80000,                 # 20000 (article, 32 ranks) hit the cap at 64 ranks (job 1960817)
    "monitorVariables": ["resrho", "resturb"],
    # adjoint
    "ADPC": True,
    "adjointL2Convergence": 1e-9,
    "adjointMaxIter": 10000,
    "adjointSolver": "GMRES",
    "globalPreconditioner": "additive Schwarz",
    "adjointSubspaceSize": 400,
    "ILUFill": 3,
    "ASMOverlap": 3,
    "outerPreconIts": 3,
    "skipAfterFailedAdjoint": False,
}

# SA-BCM tuned = bcm_2d + the GR ANK base (VLR time step, line-search tolerances, CFL limit 1e8)
# and second order at 1e-4. 24-I job 2198926 (-24 % wall time, scalar); used for every SA-BCM and SA
# matrix optimisation (job 2203553, job_mp7_mtx.sh:38-41).
RECIPES["bcm_2d_tuned"] = dict(
    RECIPES["bcm_2d"],
    ANKCharTimeStepType="VLR",
    ANKPhysicalLSTol=0.8,
    ANKUnsteadyLSTol=2.0,
    ANKCFLLimit=1.0e8,
    ANKSecondOrdSwitchTol=1e-4,
    nCycles=60000,
    adjointL2Convergence=1e-8,        # as run: job 2203553 and job_mp7_mtx*.sh
)

RECIPES["sa_2d"] = dict(RECIPES["bcm_2d"])           # plain SA = SA-BCM recipe with useSABCM False
RECIPES["sa_2d_tuned"] = dict(RECIPES["bcm_2d_tuned"])

# ---- SA-BCM and SA, 3D wing (SA-BCM article, AR5 wing, 943 104 cells) --------------------
# showcase/bcm_article/cluster_package/06_opt_3d/02_c2_opt_tu025/{run.py,vars_master.py}
# (jobs 1904652 optimisation, 1904793 baseline; 64 ranks; baseline 163 iterations / 312 s).
# ANK from CFL 1 -> second-order ANK at 1e-3; target 1e-6, so NK (switch 1e-7) is never reached.
# The article ran acousticScaleFactor = Mach (0.133) and the ADflow default vis4: both are set
# by the case script, not here. Old option names: use_SABCM -> useSABCM, SABCM_TU [%] -> turbIntensityInf.
RECIPES["bcm_3d_wing"] = {
    "smoother": "DADI",
    "MGCycle": "sg",
    "useANKSolver": True,
    "ANKCharTimeStepType": "VLR",
    "ANKCFL0": 1.0,
    "ANKSecondOrdSwitchTol": 1e-3,
    "ANKADPC": True,
    "ANKCFLLimit": 9.0e6,
    "useNKSolver": True,
    "NKSwitchTol": 1e-7,
    "NKADPC": True,
    "useBlockettes": True,
    "infChangeCorrection": True,
    "L2Convergence": 1e-6,
    "nCycles": 15000,
    "monitorVariables": ["resrho", "resturb", "totalr"],
    # adjoint
    "ADPC": True,
    "adjointL2Convergence": 1e-4,
    "adjointMaxIter": 3000,
    "adjointMonitorStep": 10,
    "adjointSubspaceSize": 200,
    "ILUFill": 4,
    "ASMOverlap": 3,
    "innerPreconIts": 3,
    "outerPreconIts": 6,
    "skipAfterFailedAdjoint": True,
}

# Plain SA on the same wing (article baseline 08_base_sa and C1 07_c1_sa): the SA-BCM recipe with a
# lower CFL ceiling and a tighter adjoint. The article ran SA with eddyVisInfRatio 1e-10 (run.py:41),
# not the SA value 9e-3: pass --nuinf to reproduce it.
RECIPES["sa_3d_wing"] = dict(
    RECIPES["bcm_3d_wing"],
    ANKCFLLimit=3.0e5,
    NKSwitchTol=1e-6,
    nCycles=10000,
    adjointL2Convergence=1e-8,
    adjointMaxIter=650,
)

# One-line description of every recipe (printed by the example scripts)
# SA-BCM (and SA) on the same mesh: coupled turbulence at 1e-2, second order at 1e-4, no NK,
# CFL limit 1e8, target 1e-7. L/15_sabcm_polars/jobs/job_c376_bcm_mtx.sh, job 1849321 (2 nodes,
# 128 ranks, cold start, 9.2 h: stopped at 6.5 orders). That job also set acousticScaleFactor 0.15
# and eddyVisInfRatio 5e-7; the case script sets them.
RECIPES["bcm_3d_wing_matrix"] = dict(
    RECIPES["bcm_3d_wing"],
    ANKCFL0=5.0,
    ANKCoupledSwitchTol=1e-2,
    ANKSecondOrdSwitchTol=1e-4,
    ANKCFLLimit=1.0e8,
    ANKPhysicalLSTol=0.8,
    ANKUnsteadyLSTol=2.0,
    useNKSolver=False,
    L2Convergence=1e-7,
    nCycles=60000,
)

DESCRIPTION = {
    "gr_2d_scalar": "GR 2D scalar: ANK -> CANK 1e-5 -> CSANK 1e-6, no NK (campaign 24, jobs 1962375-79)",
    "gr_2d_scalar_polar": "GR 2D scalar polars: ANK -> CANK 1e-2 -> CSANK 1e-5, no NK (24-B, jobs 1961834-39)",
    "gr_2d_matrix": "GR 2D matrix: ANK -> CANK 1e-2 -> CSANK 1e-4, no NK (12-A job 1827893, 24-C job 2201434)",
    "gr_2d_matrix_opt": "GR 2D matrix optimisation solves: second order 1e-4, CANK 1e-5 (jobs 2219449, 2219697)",
    "gr_3d_crossflow_swept": "GR swept wing: ANK -> CANK 1e-2 -> CSANK 1e-5 -> NK 1e-6 (job 1826778)",
    "gr_3d_crossflow_swept_scalar": "GR swept wing, scalar ladder: CANK 1e-3 -> CSANK 1e-5 -> NK 1e-6 (resweep_ladder)",
    "gr_3d_crossflow_sickle": "GR sickle wing v14: ANK -> CANK 1e-3 -> CSANK 1e-6, no NK, per-node damping",
    "gr_3d_wing_ar5": "GR AR5 wing: ANK -> CANK 1e-4 -> CSANK 3e-5 -> NK 1e-6, then no-NK restart (jobs 1810888-925)",
    "gr_3d_wing_ar5_opt": "GR AR5 c376 optimisation: ANK -> CANK 1e-2 -> CSANK 1e-3, no NK, 1e-7 (jobs 1844851, 1931137)",
    "sgamma_2d": "SA-gamma 2D: ANK -> second order 1e-2, never coupled, no NK (jobs 1962375-79, 1961969-74)",
    "sgamma_2d_tuned": "SA-gamma 2D tuned: couple 1e-5, second order 1e-4 (jobs 2197279, 2199671)",
    "bcm_2d": "SA-BCM 2D: ANK -> second order 1e-3, no NK, 3 turbulence sub-iterations (jobs 1962377, 1961957-62)",
    "bcm_2d_tuned": "SA-BCM 2D tuned: VLR + line-search tolerances + CFL 1e8, second order 1e-4 (job 2198926)",
    "sa_2d": "SA 2D: the SA-BCM recipe with the model off (jobs 1962365-70)",
    "sa_2d_tuned": "SA 2D tuned: the tuned SA-BCM recipe with the model off (job_mp7_mtx.sh)",
    "bcm_3d_wing": "SA-BCM AR5 wing (article): ANK from CFL 1 -> second order 1e-3, target 1e-6 (jobs 1904652, 1904793)",
    "sa_3d_wing": "SA AR5 wing (article baseline / C1): as bcm_3d_wing, CFL limit 3e5",
    "gr_3d_wing_matrix": "GR AR5 13.7 M mesh, matrix: ANK -> CANK 1e-2 -> CSANK 1e-4, no NK, 1e-7 (job 1849168)",
    "sgamma_3d_wing_matrix": "SA-gamma AR5 13.7 M mesh, matrix: second order 1e-2, never coupled, no NK (not validated in 3D)",
    "bcm_3d_wing_matrix": "SA-BCM / SA AR5 13.7 M mesh, matrix: CANK 1e-2 -> CSANK 1e-4, no NK, 1e-7 (job 1849321)",
}
assert set(DESCRIPTION) == set(RECIPES)


def recipe2D(model, dissipation="scalar", purpose="analysis"):
    """Name of the validated recipe for a 2D solve.

    purpose: "analysis" (NACA 0012 fixed alpha / trim), "polar" (NLF0416, S809 polars) or
    "opt" (cold-started design points of an optimisation).
    """
    if model == "gr":
        if dissipation == "matrix":
            return "gr_2d_matrix_opt" if purpose == "opt" else "gr_2d_matrix"
        return "gr_2d_scalar_polar" if purpose == "polar" else "gr_2d_scalar"
    if model == "sgamma":
        return "sgamma_2d_tuned" if (purpose == "opt" and dissipation == "matrix") else "sgamma_2d"
    if model in ("bcm", "bcmo"):
        return "bcm_2d_tuned" if (purpose == "opt" and dissipation == "matrix") else "bcm_2d"
    if model == "sa":
        return "sa_2d_tuned" if (purpose == "opt" and dissipation == "matrix") else "sa_2d"
    raise ValueError(model)


def options(model, recipe, dissipation, tuPercent, gridFile, outputDir, nuinf=None, **extra):
    """Assemble one ADflow option dict: COMMON + dissipation + model + recipe + extra (later wins)."""
    opts = dict(COMMON)
    opts.update(DISSIPATION[dissipation])
    opts.update(modelBlock(model, tuPercent, nuinf))
    opts.update(RECIPES[recipe])
    opts.update(gridFile=gridFile, outputDirectory=outputDir)
    opts.update(extra)
    return opts
