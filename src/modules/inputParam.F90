module inputDiscretization
    !
    !       Input parameters which are related to the discretization of
    !       the governing equations, i.e. scheme parameters, time accuracy
    !       (in case of an unsteady computation) and preconditioning info.
    !
    use constants, only: intType, realType
    implicit none
    save
    !
    !       Definition of the discretization input parameters.
    !
    ! spaceDiscr:             Fine grid discretization.
    ! spaceDiscrCoarse:       Coarse grid discretization.
    ! orderTurb:              Order of the discretization of the advective
    !                         terms of the turbulent transport equations.
    !                         Possibilities are 1st and 2nd order.
    ! riemann:                Fine grid riemann solver, upwind schemes only.
    ! riemannCoarse:          Idem, but on the coarse grids.
    ! limiter:                Limiter, upwind schemes only.
    ! precond:                Preconditioner.
    ! eulerWallBCTreatment:   Wall boundary condition treatment for inviscid
    !                         simulations.
    ! viscWallBCTreatment:    Wall boundary condition treatment for viscous
    !                         simulations.
    ! outflowTreatment:       Treatment of the outflow boundaries. Either
    !                         constantExtrapol or linExtrapol.
    ! nonMatchTreatment:      Treatment of the non-matching block
    !                         boundaries. Either NonConservative or
    !                         Conservative.
    ! vis2:                   Coefficient of the second order dissipation.
    ! vis4:                   Coefficient of the fourth order dissipation.
    ! vis2Coarse:             Coefficient of the second order dissipation
    !                         on the coarser grids in the mg cycle. On the
    !                         coarser grids a first order scheme is used.
    ! adis:                   Exponent for directional scaling of the
    !                         dissipation. adis == 0: no directional scaling,
    !                                      adis == 1: isotropic dissipation.
    ! kappaCoef:              Coefficient in the upwind reconstruction
    !                         schemes, both linear and nonlinear.
    ! vortexCorr:             Whether or not a vortex correction must be
    !                         applied. Steady flow only.
    ! dirScaling:             Whether or not directional scaling must be
    !                         applied.
    ! hScalingInlet:          Whether or not the outgoing Riemann invariant
    !                         must be scaled for a subsonic inlet. May be
    !                         needed for stability when strong total
    !                         temperature gradients are present.
    ! radiiNeededFine:        Whether or not the spectral radii are needed
    !                         to compute the fluxes of the fine grid.
    ! radiiNeededCoarse:      Idem for the coarse grid.
    ! lumpedDiss :            logical factor for determining whether or not
    !                         lumped dissipation is used for preconditioner
    ! approxSA:               Determines if the approximate source terms form
    !                         the SA model is used.
    ! sigma      :            Scaling parameter for dissipation lumping in
    !                         approximateprecondtioner
    ! useApproxWallDistance : logical to determine if the user wants to
    !                         use the fast approximate wall distance
    !                         computations. Typically only used for
    !                         repeated calls when the wall distance would
    !                         not have changed significantly
    ! updateWallAssociations: Logical to determine if the full wall distance
    !                         assocation is to be performed on the next
    !                         wall distance calculation. This is only
    !                         significant when useApproxWallDistance is
    !                         set to True. This allows the user to
    !                         reassociate the face a cell is associated
    !                         with.
    ! lowspeedpreconditoner:  Whether or not to use low-speed precondioner

    integer(kind=intType) :: spaceDiscr, spaceDiscrCoarse
    integer(kind=intType) :: orderTurb, limiter
    integer(kind=intType) :: riemann, riemannCoarse, precond
    integer(kind=intType) :: eulerWallBCTreatment, viscWallBCTreatment, outflowTreatment
    integer(kind=intType) :: nonMatchTreatment

    real(kind=realType) :: vis2, vis4, vis2Coarse, adis
    real(kind=realType) :: acousticScaleFactor
    real(kind=realType) :: kappaCoef
    logical :: lumpedDiss
    logical :: approxSA
    real(kind=realType) :: sigma
    logical :: useBlockettes

#ifndef USE_TAPENADE
    real(kind=realType) :: vis2b, vis4b, vis2Coarseb, adisb
    real(kind=realType) :: kappaCoefb
    real(kind=realType) :: sigmab
#endif
    logical :: vortexCorr, dirScaling, hScalingInlet
    logical :: radiiNeededFine, radiiNeededCoarse

    logical :: useApproxWallDistance
    logical :: updateWallAssociations
    logical :: lowSpeedPreconditioner
end module inputDiscretization

!      ==================================================================

module inputDissipation
    !
    !       Input parameters of the artificial dissipation scheme that
    !       cannot live in inputDiscretization, because the Tapenade
    !       generated flux routines use-associate that module in full and
    !       still declare local parameters with these names.
    !
    !       epsAcoustic:  Swanson-Turkel eigenvalue limiter coefficient for
    !                     the acoustic eigenvalues of the matrix
    !                     dissipation scheme.
    !       epsShear:     Swanson-Turkel eigenvalue limiter coefficient for
    !                     the shear (convective) eigenvalues of the matrix
    !                     dissipation scheme.
    !
    use precision, only: realType
    implicit none
    save

    real(kind=realType) :: epsAcoustic, epsShear
end module inputDissipation


!      ==================================================================

module inputIO
    !
    !       Input parameters which are related to io issues, like file
    !       names and corresponding info.
    !
    use constants
    implicit none
    save
    !
    !       Definition of the IO input parameters.
    !
    ! paramFile:           Parameter file, command line argument.
    ! firstWrite:          Whether or not this is the first time a
    !                      solution is written. Needed when different
    !                      file formats are used for reading and
    !                      writing.
    ! gridFile:            Grid file.
    ! newGridFile:         File to which the changed grid is
    !                      written. Needed for moving and/or
    !                      deforming geometries.
    ! restartFiles:        Restart solution files; for cgns this
    !                      could be the same as the grid file, but
    !                      not necesarrily.
    ! solFile:             Solution file; for cgns this could be the
    !                      same as the grid or restart file, but not
    !                      necesarrily.
    ! surfaceSolFile:      Surface solution file.
    ! sliceSolFile:        File name of a slice of a surface solution. TEMPORARY
    ! liftDistributionFile:File name of a lift file. TEMPORARY
    ! cpFile:              File which contains the curve fits for cp.
    ! precisionGrid:       Precision of the grid file to be written.
    !                      Possibilities are precisionSingle and
    !                      precisionDouble.
    ! precisionSol:        Idem for the solution file(s).

    ! precisionSurfGrid:   Precision of the grid in the surface file
    ! precisionSurfSol:    Precision of the solution in the surface file
    ! storeRindLayer:      Whether or not to store 1 layer of rind
    !                      (halo) cells in the solution file.
    ! checkRestartSol:     Whether or not the solution in the restart
    !                      file must be checked for correct
    !                      nondimensionalization.
    ! autoParameterUpdate: Whether or not the parameter file must be
    !                      updated automatically. After a restart file
    !                      is written, such that a restart can be made
    !                      without editing the parameter file.
    ! writeCoorMeter:      Whether or not the coordinates in the
    !                      solution files must be written in meters.
    !                      If not, the original units are used.
    ! storeConvInnerIter:  Whether or not to store the convergence of
    !                      the inner iterations for unsteady mode.
    !                      On systems with a limited amount of memory
    !                      the storage of this info could be a
    !                      bottleneck for memory.

    integer(kind=intType) :: precisionGrid, precisionSol
    integer(kind=intType) :: precisionSurfGrid, precisionSurfSol
    character(len=maxStringLen) :: paramFile, gridFile
    character(len=maxStringLen) :: newGridFile
    character(len=maxStringLen) :: solFile
    character(len=maxstringlen), dimension(:), allocatable :: restartFiles
    character(len=maxStringLen) :: surfaceSolFile, cpFile, sliceSolFile, liftDistributionFile

    logical :: storeRindLayer, checkRestartSol
    logical :: autoParameterUpdate, writeCoorMeter
    logical :: storeConvInnerIter
    logical :: firstWrite = .true.
    logical :: viscousSurfaceVelocities = .True.

    ! Extra file names (set from python) that specify the name of
    ! the volume, surface, lift and slice files written from an
    ! interrupt.
    character(len=maxStringLen) :: forcedSurfaceFile, forcedVolumeFile
    character(len=maxStringLen) :: forcedLiftFile, forcedSliceFile
    character(len=maxStringLen) :: convSolFileBasename
    ! logical to control the us of the transition model
    logical :: laminarToTurbulent

end module inputIO

!      ==================================================================

module inputIteration
    !
    !       Input parameters which are related to the iteration process,
    !       i.e. multigrid parameters, cfl numbers, smoothers and
    !       convergence.
    !
    use constants
    implicit none
    save
    !
    !       Definition of the iteration input parameters.
    !
    ! nCycles:          Maximum number of multigrid cycles.
    ! nCyclesCoarse:    Idem, but on the coarse grids in full multigrid.
    ! nSaveVolume:      Number of fine grid cycles after which a volume
    !                   solution file is written.
    ! nSaveSurface:     Number of fine grid cycles after which a
    !                   surface solution file is written.
    ! nsgStartup:       Number of single grid iterations, before
    !                   switching to multigrid. Could be useful for
    !                   supersonic problems with strong shocks.
    ! nSubIterTurb:     Number of turbulent subiterations when using
    !                   a decoupled approach for the turbulence.
    ! nUpdateBleeds:    Number of iterations after which the bleed
    !                   boundary conditions must be updated.
    ! smoother:         Smoother to be used.
    ! nRKStages:        Number of stages in the runge kutta scheme.
    ! nSubiterations:   Maximum number of subiterations used in
    !                   DADI.
    ! turbTreatment:    Treatment of the turbulent transport equations;
    !                   either decoupled or coupled.
    ! turbSmoother:     Smoother to use in case a decoupled solver
    !                   is to be used.
    ! turbRelax:        What kind of turbulent relaxation to use.
    !                   Either turbRelaxExplicit or
    !                   turbRelaxImplicit.
    ! resAveraging:     What kind of residual averaging to use.
    ! freezeTurbSource: Whether or not the turbulent source terms must
    !                   be frozen on the coarser grid levels; only if
    !                   a coupled solver is to be used.
    ! mgBoundCorr:      Treatment of the boundary halo's for the
    !                   multigrid corrections. Either dirichlet0,
    !                   set the corrections to zero, or neumann.
    ! mgStartlevel:     Grid level on which the multigrid must be
    !                   started in the full mg cycle. In case a restart
    !                   is specified this info is overruled and the
    !                   start level is the finest grid.
    ! nMGSteps:         Number of steps in the array cycleStrategy.
    ! nMGLevels:        Number of levels in the multigrid. This info
    !                   is derived from the cycle strategy.
    ! cycleStrategy:    Array which describes the mg cycle.
    ! cfl:              Cfl number on the fine grid.
    ! cflCoarse:        Idem, but on the coarse grids.
    ! cfllimit          Limit used to determine how much residuals are smoothed
    ! alfaTurb:         Relaxation factor in turbulent dd-adi smoother.
    ! betaTurb:         Relaxation factor in vf dd-adi smoother.
    ! relaxBleeds:      Relaxation coefficient for the update
    !                   of the bleed boundary condition.
    ! smoop:            Coefficient in the implicit smoothing.
    ! fcoll:            Relaxation factor for the restricted residuals.
    ! L2Conv:           Relative L2 norm of the density residuals for
    !                   which the computation is assumed converged.
    ! L2ConvCoarse:     Idem, but on the coarse grids during full mg.
    ! etaRk:            Coefficients in the runge kutta scheme. The
    !                   values depend on the number of stages specified.
    ! cdisRk:           Dissipative coefficients in the runge kutta
    !                   scheme. The values depend on the number of
    !                   stages specified.
    ! printIterations:  If True, iterations are printed to stdout
    ! turbresscale:     Scaling factor for turbulent residual. Necessary for
    !                   NKsolver with RANS. Only tested on SA.
    ! meshMaxSkewness   If one cell has a highe skewness than this, the Solver
    !                   errors out.
    ! iterType : String used for specifying which type of iteration was taken
    !
    ! Definition of the string, which stores the multigrid cycling
    ! strategy.
    !

    integer(kind=intType) :: nCycles, nCyclesCoarse
    integer(kind=intType) :: nSaveVolume, nSaveSurface
    integer(kind=intType) :: nsgStartup, smoother, nRKStages
    integer(kind=intType) :: nSubiterations
    integer(kind=intType) :: nSubIterTurb, nUpdateBleeds
    integer(kind=intType) :: resAveraging
    real(kind=realType) :: CFLLimit
    integer(kind=intType) :: turbTreatment, turbSmoother, turbRelax
    integer(kind=intType) :: mgBoundCorr, mgStartlevel
    integer(kind=intType) :: nMGSteps, nMGLevels
    real(kind=realType) :: timeLimit
    integer(kind=intType), allocatable, dimension(:) :: cycleStrategy
    integer(kind=intType) :: miniterNum
    real(kind=realType) :: cfl, cflCoarse, fcoll, smoop
    real(kind=realType) :: alfaTurb, betaTurb
    real(kind=realType) :: L2Conv, L2ConvCoarse
    real(kind=realType) :: L2ConvRel
    real(kind=realType) :: maxL2DeviationFactor
    real(kind=realType) :: relaxBleeds
    real(kind=realtype) :: epscoefconv
    integer(kind=inttype) :: convcheckwindowsize
    real(kind=realType), allocatable, dimension(:) :: etaRK, cdisRK
    character(len=maxStringLen) :: mgDescription
    logical :: rkReset
    logical :: useLinResMonitor
    logical :: freezeTurbSource
    logical :: printIterations
    logical :: printWarnings
    logical :: printNegativeVolumes
    logical :: printBadlySkewedCells
    logical :: printBCWarnings
    real(kind=realType), dimension(4) :: turbResScale = (/1.0_realType, &
                                                          1.0_realType, &
                                                          1.0_realType, &
                                                          1.0_realType/)
    real(kind=realType) :: meshMaxSkewness
    logical :: useSkewnessCheck
    logical :: useDissContinuation
    real(kind=realType) :: dissContMagnitude, dissContMidpoint, dissContSharpness
    integer(kind=intType) :: TurbDADICoupled = 2
    logical :: transitionFirstOrderUpwind = .true.
    logical :: transitionCrossflow = .true.
    real(kind=realType) :: transitionRoughnessHeight = 3.3e-6_realType
    ! Master switch for the NK/ANK/turbKSP convergence-acceleration bundle
    ! (column scaling + Eq. 59 source-dt restriction, incl. NK reactivation-
    ! on-backtrack): default True keeps today's behavior for SA-Gamma-Retheta
    ! unchanged (all of it is still additionally gated on
    ! turbModel==spalartallmarasnoft2gammaretheta / nwt==3, since the
    ! underlying srcLambda eigenvalue machinery is model-specific). Exposed
    ! as its own switch so it can be turned off independently, or reused
    ! later for another stiff turbulence model without touching the
    ! turbModel gates themselves.
    logical :: transitionNK = .true.
    ! One-way runtime latch for the NK phase only: once the Newton residual
    ! (relative to totalR0) drops below transitionNKAutoDisableTol, the whole
    ! transitionNK bundle (PC/column scaling, Eq. 59 reactivation, Algorithm 2
    ! damping) is turned off for the remainder of the NK phase, i.e. NK falls
    ! back to ADflow's native turbModel-agnostic behavior. Added 2026-07-18 to
    ! probe an NK stall (Step pinned near minlambda for 100s of iterations
    ! around scaledTotalRes~1.7e-3, nk_switch_crossing_test) in case the
    ! bundle's PC scaling is contributing indirectly to the stalled steps.
    ! Default 0.0 disables the latch: NK phase always follows the static
    ! transitionNK option, unchanged behavior. Not user-facing state beyond
    ! this trigger tolerance -- transitionNKActive itself is internal
    ! (like srcDtRestrictActive/noBacktrackCount), reset each solve by
    ! setupNKSolver.
    real(kind=realType) :: transitionNKAutoDisableTol = 0.0_realType
    logical :: transitionNKActive = .true.
    ! Stall escape (2026-07-18, nk_switch_crossing_test): once Step
    ! (stepMonitor) has been below transitionNKStallStepTol for
    ! transitionNKStallCountTrigger consecutive NK iterations, force the
    ! Eisenstat-Walker linear-solve rtol down to transitionNKStallRtolCap
    ! for that iteration. Rationale: getEWTol's rtol = (norm/oldNorm)^1.618
    ! -- when stalled, norm~oldNorm so the ratio~1 and rtol rises to its
    ! 0.8 cap, i.e. EW asks for the LOOSEST linear solve exactly when the
    ! stall needs the tightest one (a poor direction from an under-solved
    ! linear system is the likely cause of the tiny step in the first
    ! place). transitionNKStallRtolCap=1.0 (default) disables this (never
    ! caps below whatever EW already picked); nkStallCount is internal
    ! (like noBacktrackCount), reset each NK entry in NKStep.
    real(kind=realType) :: transitionNKStallStepTol = 0.1_realType
    integer(kind=intType) :: transitionNKStallCountTrigger = 3
    real(kind=realType) :: transitionNKStallRtolCap = 1.0_realType
    integer(kind=intType) :: nkStallCount = 0
    ! Eq. 58 (P&Z 2020) geometric row-scaling factor, on top of the existing
    ! turbResScale 1/max row scaling: multiplies NK's residual rows (setRVec)
    ! by volRef**(5/3) for flow rows, volRef**(2/3) for turb rows, bringing
    ! ADflow's uniform 1/volRef row normalization toward the paper's
    ! per-block J^{2/3}/J^{-1/3} exponents. NOT verified to match the
    ! paper's SBP metric Jacobian J exactly (ADflow's volRef is physical
    ! cell volume, not necessarily the same normalization) -- off by
    ! default pending validation, same caution level as
    ! transitionResidualAutoscale.
    logical :: transitionRowVolScale = .false.
    ! Eq. 58 (P&Z 2020) residual-autoscaling (S_a) proxy: the paper gives no
    ! formula (cites Osusky & Zingg's thesis, unavailable here). This is an
    ! adaptive per-equation-block row rescaling from monitored residual-norm
    ! ratios -- same intent, NOT verified identical to their method. Off by
    ! default.
    logical :: transitionResidualAutoscale = .false.
    ! Per-block autoscale factor [flow, nuTilde, gamma, reTheta] computed by
    ! NKSolver's computeNKResidualAutoscale (lagged, same cadence as
    ! NK_jacobianLag). Lives here (not in NKSolver) so both NKSolver and
    ! utils (generic monitor/residual printing) can read it without a
    ! circular module dependency.
    real(kind=realType) :: nkAutoScaleFac(4) = one
    logical :: transitionSrcDtRestrict = .true.
    real(kind=realType) :: transitionSrcDtLimit = 0.9_realType
    ! srcDtDeactivateIters: deactivate after N clean ANK turb iters (P&Z §IV.B.3)
    integer(kind=intType) :: srcDtDeactivateIters = 5
    integer(kind=intType) :: noBacktrackCount = 0
    ! Damping parameters for γ/Reθt iterative update (P&Z Algorithm 2).
    ! MaxIter is a safety cap on the (unbounded in the paper) back-off loop;
    ! 10000 makes it effectively unbounded (0.99^10000 ~ 0), so the hard
    ! clip after it is a warned last-resort fallback, not the working
    ! bounds mechanism.
    real(kind=realType) :: transitionDampTheta = 0.99_realType
    integer(kind=intType) :: transitionDampMaxIter = 10000
    ! Use approxSA simplification in SA-gamma-rethetha (default true for stability)
    logical :: transitionUseApproxSA = .true.
    ! Reference length l [grid units] in the vorticity limiter (P&Z Eqs. 52-53);
    ! the paper uses the root chord. Negative => auto: use lengthRef
    ! (the AeroProblem chordRef).
    real(kind=realType) :: transitionRefLength = -1.0_realType

    ! ------------------------------------------------------------------
    ! Stall diagnostics (VERIF_06). Purely diagnostic: when enabled the
    ! ANK/CANK and NK step controllers report WHY the step collapsed --
    ! which limiter bound the global step (rho/E physicality, turbulence
    ! bound, unsteady line search, cubic line search), where the binding
    ! cell is, and whether the CFL cutback was floored out by ANK_CFLMin.
    ! No effect on the solution path; default off so production logs are
    ! unchanged.
    !
    ! solverStallDiagStep: only report on iterations whose accepted step
    ! is below this value (1.0 => report every iteration).
    logical :: solverStallDiag = .false.
    real(kind=realType) :: solverStallDiagStep = 1.0_realType

    ! ------------------------------------------------------------------
    ! ANK/CANK CFL floor cap (VERIF_06 F1).
    !
    ! ADflow ramps the CFL floor with convergence,
    !     ANK_CFLMin = min(ANK_CFLLimit, ANK_CFLMinBase*(totalR0/totalR)**ANK_CFLExponent),
    ! and every CFL update is clipped from below by it. At depth the floor
    ! reaches ANK_CFLLimit, at which point the cutback
    !     ANK_CFL = max(ANK_CFL*ANK_CFLCutback, ANK_CFLMin)
    ! cannot reduce ANK_CFL at all -- the coupled step controller loses its
    ! only recovery mechanism exactly when a collapsing step needs it, and
    ! the solver limit-cycles (huge CFL -> bad step -> backtrack -> repeat).
    !
    ! Piotrowski's thesis (§3.1.3, Algorithms 2 and 4) instead halves the
    ! reference time step against a USER-SPECIFIED bound, dt_ref,min. This
    ! option supplies that bound: an absolute ceiling on how far the ramped
    ! floor may rise, so the cutback always has room to act.
    !
    ! <= 0 (default) => disabled, floor ramps exactly as before.
    real(kind=realType) :: ankCFLMinCap = 0.0_realType

    ! ------------------------------------------------------------------
    ! ANK/CANK unsteady line-search geometry (VERIF_06 F2).
    !
    ! ADflow backtracks by 0.7 for at most 12 iterations, so the smallest
    ! step the line search can produce is 0.7**12 = 0.0138 -- ABOVE its own
    ! rejection threshold ANK_stepMin*ANK_stepFactor (0.01 by default).
    ! The search therefore cannot classify a collapsed step as a failure,
    ! the CFL cutback never fires, the CFL stays pinned at ANK_CFLLimit,
    ! and the solver limit-cycles: huge proposed step -> throttled to ~1-3%
    ! -> repeat. Measured on the NLF(2)-0415 swept wing: CFL locked at 1e6
    ! for 1800+ iterations with a healthy linear residual.
    !
    ! Piotrowski's thesis (Algorithm 4) backtracks by 0.90 with an explicit
    ! 1% floor and rejects on reaching it -- gentler, allowed to go further,
    ! and its floor coincides with the rejection threshold by construction.
    ! From lambda = 1 that needs ~44 backtracks to reach 0.01.
    !
    ! Defaults reproduce ADflow's current 0.7 / 12 geometry exactly.
    real(kind=realType) :: ankUnsteadyLSFactor = 0.7_realType
    integer(kind=intType) :: ankUnsteadyLSMaxIter = 12

    ! Treat an exhausted backtracking budget as a step rejection (thesis
    ! Algorithm 4), so the CFL cutback actually fires instead of the solver
    ! accepting the floor step forever. Default .false. = previous behaviour.
    logical :: ankRejectOnLSExhausted = .false.

    ! Algorithm 2 (P&Z §IV.B.2) per-node gamma/Re-theta-t damping in the
    ! COUPLED ANK path (VERIF_06 F7). It previously existed only in NK, but
    ! the thesis's inexact-Newton phase is CSANK, so the paper's algorithm
    ! needs it here. Default .false. = previous behaviour.
    !
    ! Note: per-cell damping was tried in physicalityCheckANK and reverted
    ! (see the comment there), but that comparison predates the F2 line-search
    ! fix and so was measured against a global-lambda controller stuck in a
    ! limit cycle. Re-test with F2 active before drawing a conclusion.
    logical :: ankAlgorithm2Damping = .false.

end module inputIteration

module inputCostFunctions
    use constants
    real(kind=realtype) :: sepSensorOffset = zero
    real(kind=realtype) :: sepSensorSharpness = 10.0_realType
    real(kind=realtype) :: cavSensorOffset
    real(kind=realtype) :: cavSensorSharpness
    integer(kind=inttype) :: cavExponent
    logical :: computeCavitation

end module inputCostFunctions

!      ==================================================================

module inputMotion
    !
    !       Input parameters which are related to the rigid body motion of
    !       the entire mesh, i.e. translation and rotation.
    !       These parameters can only be specified for an external flow
    !       computation.
    !
    use precision
    implicit none
    save
    ! rotPoint(3): Rotation point of the rigid body rotation.

    real(kind=realType), dimension(3) :: rotPoint
    real(kind=realType), dimension(3) :: rotPointd

    ! degreePolXRot: Degree of the x-rotation polynomial.
    ! degreePolYRot: Degree of the y-rotation polynomial.
    ! degreePolZRot: Degree of the z-rotation polynomial.

    integer(kind=intType) :: degreePolXRot
    integer(kind=intType) :: degreePolYRot
    integer(kind=intType) :: degreePolZRot

    ! coefPolXRot(0:): coefficients of the x-rotation polynomial.
    ! coefPolYRot(0:): coefficients of the y-rotation polynomial.
    ! coefPolZRot(0:): coefficients of the z-rotation polynomial.

    real(kind=realType), dimension(:), allocatable :: coefPolXRot
    real(kind=realType), dimension(:), allocatable :: coefPolYRot
    real(kind=realType), dimension(:), allocatable :: coefPolZRot

    ! degreeFourXRot: Degree of the x-rotation fourier series.
    ! degreeFourYRot: Degree of the y-rotation fourier series.
    ! degreeFourZRot: Degree of the z-rotation fourier series.

    integer(kind=intType) :: degreeFourXRot
    integer(kind=intType) :: degreeFourYRot
    integer(kind=intType) :: degreeFourZRot

    ! omegaFourXRot: Fourier frequency of the x-rotation; the
    !                   period of the motion is 2*pi/omega.
    ! omegaFourYRot: Fourier frequency of the y-rotation.
    ! omegaFourZRot: Fourier frequency of the z-rotation.

    real(kind=realType) :: omegaFourXRot, omegaFourXRotb
    real(kind=realType) :: omegaFourYRot, omegaFourYRotb
    real(kind=realType) :: omegaFourZRot, omegaFourZRotb

    ! cosCoefFourXRot(0:): cosine coefficients of the
    !                      x-rotation fourier series.
    ! cosCoefFourYRot(0:): cosine coefficients of the
    !                      y-rotation fourier series.
    ! cosCoefFourZRot(0:): cosine coefficients of the
    !                      z-rotation fourier series.

    real(kind=realType), dimension(:), allocatable :: cosCoefFourXRot
    real(kind=realType), dimension(:), allocatable :: cosCoefFourYRot
    real(kind=realType), dimension(:), allocatable :: cosCoefFourZRot

    ! sinCoefFourXRot(1:): sine coefficients of the
    !                      x-rotation fourier series.
    ! sinCoefFourYRot(1:): sine coefficients of the
    !                      y-rotation fourier series.
    ! sinCoefFourZRot(1:): sine coefficients of the
    !                      z-rotation fourier series.

    real(kind=realType), dimension(:), allocatable :: sinCoefFourXRot
    real(kind=realType), dimension(:), allocatable :: sinCoefFourYRot
    real(kind=realType), dimension(:), allocatable :: sinCoefFourZRot

    ! degreePolAlpha: Degree of the Alpha polynomial.

    integer(kind=intType) :: degreePolAlpha

    ! coefPolAlpha(0:): coefficients of the Alpha polynomial.

    real(kind=realType), dimension(:), allocatable :: coefPolAlpha
    real(kind=realType), dimension(:), allocatable :: coefPolAlphab

    ! degreeFourAlpha: Degree of the Alpha fourier series.

    integer(kind=intType) :: degreeFourAlpha

    ! omegaFourAlpha: Fourier frequency of the Alpha; the
    !                   period of the motion is 2*pi/omega.

    real(kind=realType) :: omegaFourAlpha, omegafouralphab

    ! cosCoefFourAlpha(0:): cosine coefficients of the
    !                      x-rotation fourier series.

    real(kind=realType), dimension(:), allocatable :: cosCoefFourAlpha
    real(kind=realType), dimension(:), allocatable :: cosCoefFourAlphab

    ! sinCoefFourAlpha(1:): sine coefficients of the
    !                      Alpha fourier series.

    real(kind=realType), dimension(:), allocatable :: sinCoefFourAlpha
    real(kind=realType), dimension(:), allocatable :: sinCoefFourAlphab

    ! degreePolXRot: Degree of the Beta polynomial.

    integer(kind=intType) :: degreePolBeta

    ! coefPolXRot(0:): coefficients of the Beta polynomial.

    real(kind=realType), dimension(:), allocatable :: coefPolBeta
    real(kind=realType), dimension(:), allocatable :: coefPolBetab

    ! degreeFourBeta: Degree of the Beta fourier series.

    integer(kind=intType) :: degreeFourBeta

    ! omegaFourBeta: Fourier frequency of the Beta; the
    !                   period of the motion is 2*pi/omega.

    real(kind=realType) :: omegaFourBeta, omegafourbetab

    ! cosCoefFourBeta(0:): cosine coefficients of the
    !                      Beta fourier series.

    real(kind=realType), dimension(:), allocatable :: cosCoefFourBeta
    real(kind=realType), dimension(:), allocatable :: cosCoefFourBetab

    ! sinCoefFourBeta(1:): sine coefficients of the
    !                      Beta fourier series.

    real(kind=realType), dimension(:), allocatable :: sinCoefFourBeta
    real(kind=realType), dimension(:), allocatable :: sinCoefFourBetab

    ! degreePolMach: Degree of the Mach polynomial.

    integer(kind=intType) :: degreePolMach

    ! coefPolMach(0:): coefficients of the Mach polynomial.

    real(kind=realType), dimension(:), allocatable :: coefPolMach
    real(kind=realType), dimension(:), allocatable :: coefPolMachb

    ! degreeFourMach: Degree of the Mach fourier series.

    integer(kind=intType) :: degreeFourMach

    ! omegaFourMach: Fourier frequency of the Mach Number; the
    !                   period of the motion is 2*pi/omega.

    real(kind=realType) :: omegaFourMach, omegafourmachb

    ! cosCoefFourMach(0:): cosine coefficients of the
    !                      Mach Number fourier series.

    real(kind=realType), dimension(:), allocatable :: cosCoefFourMach
    real(kind=realType), dimension(:), allocatable :: cosCoefFourMachb

    ! sinCoefFourMach(1:): sine coefficients of the
    !                      Mach Number fourier series.

    real(kind=realType), dimension(:), allocatable :: sinCoefFourMach
    real(kind=realType), dimension(:), allocatable :: sinCoefFourMachb

    ! gridMotionSpecified: Whether or not a rigid body motion of
    !                      the grid has been specified.

    logical :: gridMotionSpecified

end module inputMotion

!      ==================================================================

module inputParallel
    !
    !       Input parameters which are related to the parallelization.
    !
    use precision
    implicit none
    save

    ! loadImbalance: Allowable load imbalance
    ! splitBlocks:   Whether or not blocks can be split to improve
    !                the load balance.
    ! loadBalanceIter: The number of refinment iterations to run to try
    !                   to get better load balancing.
    real(realType) :: loadImbalance
    logical :: splitBlocks
    integer(kind=inttype) :: loadBalanceIter, partitionlikenproc
end module inputParallel

!      ==================================================================

module inputPhysics
    !
    !       Input parameters which are related to the physics of the flow,
    !       like governing equations, mode of the equations, turbulence
    !       model and free stream conditions.
    !
    use precision
    implicit none
    save

    !       Definition of the physics input parameters.
    !
    ! equations:           Governing equations to be solved.
    ! equationMode:        Mode of the equations, steady, unsteady
    !                      or timeSpectral.
    ! flowType:            Type of flow, internal or external.
    ! cpModel:             Which cp model, constant or function of
    !                      temperature via curve fits.
    ! turbModel:           Turbulence model.
    ! turbProd:            Which production term to use in the transport
    !                      turbulence equations, strain, vorticity or
    !                      kato-launder.
    ! rvfN:                Determines the version of v2f turbulence model.
    ! rvfB:                Whether or not to solve v2f with an
    !                      upper bound.
    ! useQCR:              Determines if the QCR term is applied to the shear tensor computation
    !                      when considering turbulence model effects
    ! useRotationSA:       Determines if we will use rotation correction (SA model only)
    ! useft2SA:            Determines if we will use the ft2 term (SA model only)
    ! wallFunctions:       Whether or not to use wall functions.
    ! wallDistanceNeeded:  Whether or not the wall distance is needed
    !                      for the turbulence model in a RANS problem.
    ! Mach:                Free stream Mach number.
    ! MachCoef:            Mach number used to compute coefficients;
    !                      only relevant for translating geometries.
    ! MachGrid:            Mach number of the Mesh. Used in stability
    !                      derivative calculations. Specified as the
    !                      negative of the desired freestream Mach number.
    !                      When this option is set, set Mach = 0.0...
    ! velDirFreestream(3): Direction of the free-stream velocity.
    !                      Internally this vector is scaled to a unit
    !                      vector, so there is no need to specify a
    !                      unit vector. Specifying this vector solves
    !                      the problem of angle of attack and yaw angle
    !                      definition as well as the direction of the
    !                      axis (e.g. y- or z-axis in spanwise direction).
    ! liftDirection(3):    Direction vector for the lift.
    ! dragDirection(3):    Direction vector for the drag.
    ! gammaConstant:       Constant specific heat ratio.
    ! RGasDim:             Gas constant in S.I. units.
    ! Prandtl:             Prandtl number.
    ! PrandtlTurb:         Turbulent prandtl number.
    ! pklim:               Limiter for the production of k, the production
    !                      is limited to pklim times the destruction.
    ! wallOffset:          Offset from the wall when wall functions
    !                      are used.
    ! eddyVisInfRatio:     Free stream value of the eddy viscosity.
    ! turbIntensityInf:    Free stream value of the turbulent intensity.
    ! surfaceRef:          Reference area for the force and moments
    !                      computation.
    ! lengthRef:           Reference length for the moments computation.
    ! pointRef(3):         Moment reference point.
    ! pointRefEC(3):       Elastic center. Bending moment refernce point
    ! SSuthDim:            Sutherlands law temperature (SI Units)
    ! muSuthDim:           Reference viscosity at reference temperature for Sutherlands law (SI Units)
    ! TSuthDim:            Reference temperature for Sutherlands law (SI Units)
    ! momentAxis(3,2)      Axis about which to calculate a moment, provided as 2 points in 3-D
    ! cavitationnumber     Negative Cp value that triggers the traditional
    !                      step-function based cavitation sensor.
    ! cpmin_rho            The rho parameter used with the KS-based cavitation sensor.
    ! cpmin_family         The cpmin for a given surface family that does not use
    !                      KS-aggregation, but rather an exact min computation.

    integer(kind=intType) :: equations, equationMode, flowType
    integer(kind=intType) :: turbModel, cpModel, turbProd
    integer(kind=intType) :: rvfN
    logical :: rvfB
        logical :: useQCR, useRotationSA, useft2SA

    logical :: wallFunctions, wallDistanceNeeded

    real(kind=realType) :: alpha, beta
    integer(kind=intType) :: liftIndex
    real(kind=realType) :: Mach, MachCoef, MachGrid
    real(kind=realType) :: Reynolds, ReynoldsLength
    real(kind=realType) :: gammaConstant, RGasDim
    real(kind=realType) :: Prandtl, PrandtlTurb, pklim, wallOffset, wallDistCutoff
    real(kind=realType) :: eddyVisInfRatio, turbIntensityInf
    real(kind=realType) :: surfaceRef, lengthRef
    real(kind=realType), dimension(3) :: velDirFreestream
    real(kind=realType), dimension(3) :: liftDirection
    real(kind=realType), dimension(3) :: dragDirection
    real(kind=realType), dimension(3) :: pointRef
    real(kind=realType), dimension(3, 2) :: momentAxis
    real(kind=realType) :: SSuthDim, muSuthDim, TSuthDim
    real(kind=realType) :: cavitationnumber
    real(kind=realType) :: cpmin_rho
    real(kind=realType), dimension(:), allocatable :: cpmin_family

#ifndef USE_TAPENADE
    real(kind=realType) :: alphad, betad
    real(kind=realType), dimension(3) :: velDirFreestreamd, velDirFreeStreamb
    real(kind=realType), dimension(3) :: liftDirectiond, liftDirectionb
    real(kind=realType), dimension(3) :: dragDirectiond, dragDirectionb
    real(kind=realType), dimension(3) :: pointRefd, pointRefb
    real(kind=realType), dimension(3, 2) :: momentAxisd, momentAxisb
    real(kind=realType) :: Machd, MachCoefd, MachGridd
    real(kind=realType) :: reynoldsd, reynoldslengthd
    real(kind=realType) :: gammaconstantd
    real(kind=realType) :: surfaceRefd, lengthRefd
    real(kind=realType) :: rgasdimd
    real(kind=realType) :: Prandtlb, PrandtlTurbb
#endif

    real(kind=realType), dimension(3) :: pointRefEC

    ! Return forces as tractions instead of forces:
    logical :: forcesAsTractions

end module inputPhysics

!      ==================================================================

module inputTimeSpectral
    !
    !       Input parameters for time spectral problems.
    !
    use precision
    implicit none
    save

    ! nTimeIntervalsSpectral: Number of time instances used.

    integer(kind=intType) :: nTimeIntervalsSpectral

    ! dscalar(:,:,:): Matrix for the time derivatices of scalar
    !                 quantities; different for every section to
    !                 allow for different periodic angles.
    !                 The second and third dimension equal the
    !                 number of time intervals.
    ! dvector(:,:,:): Matrices for the time derivatives of vector
    !                 quantities; different for every section to
    !                 allow for different periodic angles and for
    !                 sector periodicity.
    !                 The second and third dimension equal 3 times
    !                 the number of time intervals.

    real(kind=realType), dimension(:, :, :), allocatable :: dscalar
    real(kind=realType), dimension(:, :, :), allocatable :: dvector

    ! writeUnsteadyRestartSpectral: Whether or not a restart file
    !                               must be written, which is
    !                               capable to do a restart in
    !                               unsteady mode.
    ! dtUnsteadyRestartSpectral:    The corresponding time step.

    real(kind=realType) :: dtUnsteadyRestartSpectral
    logical :: writeUnsteadyRestartSpectral

    ! writeUnsteadyVolSpectral:  Whether or not the corresponding
    !                            unsteady volume solution files
    !                            must be written after the
    !                            computation.
    ! writeUnsteadySurfSpectral: Idem for the surface solution
    !                            files.
    ! nUnsteadySolSpectral:      The corresponding number of
    !                            unsteady solutions to be created.

    integer(kind=intType) :: nUnsteadySolSpectral
    logical :: writeUnsteadyVolSpectral
    logical :: writeUnsteadySurfSpectral

    ! rotMatrixSpectral(:,3,3):  The corresponding rotation matrices
    !                            for the velocity. No rotation
    !                            point is needed, because only the
    !                            velocities need to be transformed.
    !                            The matrix stored is the one used
    !                            when the upper bound of the mode
    !                            number is exceeded; for the lower
    !                            bound the inverse (== transpose)
    !                            must be used. The 1st dimension
    !                            is the number of sections.

    real(kind=realType), dimension(:, :, :), allocatable :: &
        rotMatrixSpectral
    logical :: useTSInterpolatedGridVelocity

    real(kind=realType) :: omegaFourier

end module inputTimeSpectral

!      ==================================================================

module inputUnsteady
    !
    !       Input parameters for unsteady problems.
    !
    use constants
    implicit none
    save

    ! timeIntegrationScheme: Time integration scheme to be used for
    !                        unsteady problems. Possibilities are
    !                        Backward difference schemes, explicit
    !                        RungeKutta schemes and implicit
    !                        RungeKutta schemes.

    integer(kind=intType) :: timeIntegrationScheme

    ! timeAccuracy:     Accuracy of the time integrator for unsteady
    !                   problems. Possibilities are 1st, 2nd and 3rd
    !                   order accurate schemes.
    ! nTimeStepsCoarse: Number of time steps on the coarse mesh;
    !                   only relevant for periodic problems for
    !                   which a full mg can be used.
    ! nTimeStepsFine:   Number of time steps on the fine mesh.
    ! deltaT:           Physical time step in seconds.

    integer(kind=intType) :: timeAccuracy
    integer(kind=intType) :: nTimeStepsCoarse, nTimeStepsFine

    real(kind=realType) :: deltaT

    ! nRKStagesUnsteady:   Number of stages used in the Runge-Kutta
    !                      schemes for a time accurate computation.
    ! betaRKUnsteady(:,:): Matrix with the Runge-Kutta coefficients
    !                      for the residuals.
    ! gammaRKUnsteady(:):  Vector with the time portion of the
    !                      Runge-Kutta stages.

    integer(kind=intType) :: nRKStagesUnsteady

    real(kind=realType), dimension(:, :), allocatable :: betaRKUnsteady
    real(kind=realType), dimension(:), allocatable :: gammaRKUnsteady

    ! nOldGridRead: Number of old grid levels read from the grid
    !               files. Needed only for a consistent restart
    !               on the deforming meshes.

    integer(kind=intType) :: nOldGridRead

    ! useALE: Use the deforming mesh ale formuation.
    logical :: useALE

    ! updateWallDistanceUnsteady: Whether or not to update the wall
    !                             distance in unsteady mode. For a
    !                             RANS simulation on a changing grid
    !                             this should be done if the
    !                             turbulence model requires the wall
    !                             distance. However, the user may
    !                             overrule this if he thinks it is
    !                             not necessary.

    logical :: updateWallDistanceUnsteady

end module inputUnsteady

module inputADjoint
    !
    !       Definition of some parameters ADjoint.
    !       The actual values of this parameters are arbitrary;
    !       in the code always the symbolic names are (should be) used.
    !
    use constants
    implicit none
    save
    !
    !       Definition of the adjoint input parameters.
    !

    ! Monitor      : Whether or not to enable the monitor for the KSP
    !                contexts.
    ! ApproxPC     : Whether or not to use the approximate jacobian
    !                preconditioner
    ! ADPC         : Whether or not to use AD for preconditioning
    ! viscPC       : Whether or not to keep cross derivative terms
    !                in viscous preconditioner.
    ! FrozenTurbulence: Whether to use frozen turbulence assumption
    ! useDiagTSPC   : Whether or not the off time instance terms are
    !                 included in the TS preconditioner.
    logical :: setMonitor, ApproxPC, useDiagTSPC
    logical :: frozenTurbulence, viscPC, ADPC

    ! ADjointSolverType: Type of linear solver for the ADjoint
    ! PreCondType      : Type of Preconditioner to use
    ! Matrix Ordering  : Type of matrix ordering to use
    ! LocalPCType      : Type of preconditioner to use on subdomains
    character(maxStringLen) :: ADjointSolverType
    character(maxStringLen) :: GMRESOrthogType
    character(maxStringLen) :: PreCondType
    character(maxStringLen) :: matrixOrdering
    character(maxStringLen) :: adjointPCSide
    character(maxStringLen) :: LocalPCType

    ! FillLevel     : Number of levels of fill for the ILU local PC
    ! Overlap       : Amount of overlap in the ASM PC
    integer(kind=intType) :: fillLevel, overlap
    integer(kind=intType) :: fillLevelCoarse, overlapCoarse

    ! adjRelTol     : Relative tolerance
    ! adjAbsTol     : Absolute tolerance
    ! adjDivTol     : Relative tolerance increase to divergence
    ! adjMaxIter    : Maximum number of iterations
    ! adjRestart    : Maximum number of steps before restart
    !                 It has a high impact on the required memory!
    ! adjMonStep    : Convergence monitor step

    real(kind=alwaysRealType) :: adjRelTol
    real(kind=alwaysRealType) :: adjAbsTol
    real(kind=alwaysRealType) :: adjRelTolRel
    real(kind=alwaysRealType) :: adjDivTol
    real(kind=realType) :: adjMaxL2Dev
    integer(kind=intType) :: adjMaxIter
    integer(kind=intType) :: adjRestart
    integer(kind=intType) :: adjMonStep

    ! storePsiHistory : Whether to buffer intermediate adjoint solution
    !                   estimates during the KSP solve, so their total
    !                   derivatives can be reported after convergence
    !                   (see MyKSPMonitor / solveAdjoint).
    ! psiHistoryStep  : Buffer a snapshot every this many KSP iterations.
    ! psiHistoryMax   : Maximum number of snapshots to buffer.
    logical :: storePsiHistory
    integer(kind=intType) :: psiHistoryStep
    integer(kind=intType) :: psiHistoryMax

    ! outerPCIts : Number of iterations to run for on (global) preconditioner
    ! intterPCIts : Number of iterations to run on local preconditioner
    integer(kind=intType) :: outerPreConIts
    integer(kind=intType) :: innerPreConIts, innerPreConItsCoarse
    integer(kind=intType) :: adjAMGLevels, adjAMGNSmooth

    logical :: printTiming
    integer(kind=intType) :: subKSPSubspaceSize
    integer(kind=intType) :: applyAdjointPCSubSpaceSize

    ! firstRun     :  Whether this is the first run of the TGT debugger
    ! verifyState  :  Whether to verify state
    ! verifySpatial:  Whether to verify spatial
    ! verifyExtra  :  Whether to verify extra
    logical :: firstRun
    logical :: verifyState
    logical :: verifySpatial
    logical :: verifyExtra

    ! Logicals for specifiying if we are using matrix-free forms of
    ! drdw
    logical :: useMatrixFreedRdw

end module inputADjoint

module inputTSStabDeriv
    !
    !       Definition of some parameters for Time Spectral stability
    !       derivatives.
    !       The actual values of this parameters are arbitrary;
    !       in the code always the symbolic names are (should be) used.
    !

    ! TSStability : Whether or not the TS stability derivatives should
    !               be computed
    logical :: TSStability, TSAlphaMode, TSBetaMode, TSpMode, &
               TSqMode, TSrMode, TSAltitudeMode, TSMachMode
    ! TSAlphaFollowing : Whether or not alpha follows the body in p,q,r mode
    logical :: TSAlphaFollowing

    ! useWindAxis : whether to rotate around the wind axis or the body
    !               axis...
    logical :: useWindAxis

end module inputTSStabDeriv

module inputOverset
    use constants
    implicit none
    save
    !
    !       Definition of parameters for the overset implementation
    !
    logical :: useoversetLoadBalance = .True.
    real(kind=realType) :: overlapFactor = 0.9
    real(kind=realType) :: nearWallDist = 0.1
    real(kind=realType) :: oversetProjTol = 1e-12
    real(kind=realType) :: backgroundVolScale = 1.0
    logical :: debugZipper = .False.
    integer(kind=intType) :: oversetUpdateMode
    real(kind=realType) :: selfZipCutoff
    ! nRefine: number of connectivity loops to run
    integer(kind=intType) :: nRefine
    integer(kind=intType) :: nFloodIter
    logical :: useZipperMesh
    logical :: useOversetWallScaling
    logical :: recomputeOverlapMatrix
    logical :: oversetDebugPrint
end module inputOverset
