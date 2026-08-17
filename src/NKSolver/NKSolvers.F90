module NKSolver

    use constants

    ! MPI comes from constants, so we need to avoid MPIF_H in PETSc
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    ! PETSc Matrices:
    ! dRdw: This is the actual matrix-free matrix computed with FD
    ! dRdwPre: The preconditoner matrix for NK method. This matrix is stored.
    ! dRdwPseudo: Shell matrix used with the pseudo-transient
    !             continuation method.
    ! dRdwNKSrcDt: Shell matrix wrapping dRdw with the SA-Gamma-Retheta
    !              Eq. 59 source-dt-restriction diagonal added on the
    !              transition rows, so GMRES solves the same restricted
    !              system the preconditioner assumes (applyNKSrcDtDiagonal)
    !              instead of the exact, unregularized Newton system.
    !              Used as the KSP operator only for that transition
    !              model's NK path (see FormJacobianNK).

    Mat dRdw, dRdwPre, dRdwPseudo, dRdwNKSrcDt

    ! PETSc Vectors:
    ! wVec: PETsc version of ADflow 'w'
    ! rVec: PETSc version of ADflow 'dw', but divided by volume
    ! deltaW: Update to the wVec from linear solution
    ! diagV: Diagonal lumping term

    Vec wVec, rVec, deltaW, work, g, baseRes

    ! NK_KSP: The ksp object for solving the newton udpate
    KSP NK_KSP

    PetscFortranAddr ctx(1)

    ! Options for NK Slver
    logical :: useNKSolver
    integer(kind=intType) :: NK_jacobianLag
    integer(kind=intType) :: NK_subspace
    integer(kind=intType) :: NK_asmOverlap
    integer(kind=intType) :: NK_asmOverlapCoarse
    integer(kind=intType) :: NK_iluFill
    integer(kind=intType) :: NK_iluFillCoarse
    integer(kind=intType) :: NK_innerPreConIts
    integer(kind=intType) :: NK_innerPreConItsCoarse
    integer(kind=intType) :: NK_outerPreConIts
    integer(kind=intType) :: NK_AMGLevels
    integer(kind=intType) :: NK_AMGNSmooth
    integer(kind=intType) :: NK_LS
    character(len=maxStringLen) :: NK_precondType
    logical :: NK_useEW
    logical :: NK_ADPC
    logical :: NK_viscPC
    real(kind=realType) :: NK_CFL0
    real(kind=realType) :: NK_switchTol
    real(kind=realType) :: NK_rtolInit
    real(kind=realType) :: NK_divTol = 10
    real(kind=realType) :: NK_fixedStep

    ! Misc variables
    logical :: NK_solverSetup = .False.
    integer(kind=intType) :: NK_iter

    ! Eisenstat-Walker Parameters
    integer(kind=intType) :: ew_version
    real(kind=realType) :: ew_rtol_0
    real(kind=realType) :: ew_rtol_max
    real(kind=realType) :: ew_gamma
    real(kind=realType) :: ew_alpha
    real(kind=realType) :: ew_alpha2
    real(kind=realType) :: ew_threshold
    real(kind=alwaysRealType) :: rtolLast, oldNorm

    ! Misc Parameters
    logical :: freeStreamResSet = .False.
    real(kind=realType) :: NK_CFL

    ! Variables for non-monotone line search
    real(kind=realType), dimension(:), allocatable :: NKLSFuncEvals
    integer(kind=intType) :: Mmax = 5
    integer(kind=intType) :: iter_k
    integer(kind=intType) :: iter_m

    ! Parameter for external preconditioner
    integer(kind=intType) :: applyPCSubSpaceSize

    ! Guarded fixed-step (SA-GR LSNone only): lowest accepted residual norm
    ! this NK phase, the reference for the global rise cap. Re-armed on every
    ! fresh NK entry in NKStep(firstCall).
    real(kind=alwaysRealType) :: nkGuardBestNorm = huge(1.0_alwaysRealType)

contains

    subroutine setupNKsolver

        ! Setup the PETSc objects for the Newton-Krylov solver.
        ! destroyNKsolver can be used to destroy the objects created in this function.

        use constants
        use stencils, only: visc_pc_stencil, euler_pc_stencil, N_visc_pc, N_euler_pc
        use communication, only: adflow_comm_world
        use inputTimeSpectral, only: nTimeIntervalsSpectral
        use inputIteration, only: useLinResMonitor
        use flowVarRefState, only: nw, viscous
        use InputAdjoint, only: viscPC
        use ADjointVars, only: nCellsLocal
        use utils, only: EChk
        use adjointUtils, only: myMatCreate, statePreAllocation
        use amg, only: setupAMG
        implicit none

        ! Working Variables
        integer(kind=intType) :: ierr, nDimw
        integer(kind=intType), dimension(:), allocatable :: nnzDiagonal, nnzOffDiag
        integer(kind=intType) :: n_stencil
        integer(kind=intType), dimension(:, :), pointer :: stencil
        integer(kind=intType) :: level

        ! We don't have memory for the approximate and exact Newton solvers kicking around at the same time.
        ! destroyANKSolver() is called in solveState in solvers.F90

        if (.not. NK_solverSetup) then
            nDimW = nw * nCellsLocal(1_intTYpe) * nTimeIntervalsSpectral

            call VecCreate(ADFLOW_COMM_WORLD, wVec, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call VecSetSizes(wVec, nDimW, PETSC_DECIDE, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call VecSetBlockSize(wVec, nw, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call VecSetType(wVec, VECMPI, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            !  Create duplicates for residual and delta
            call VecDuplicate(wVec, rVec, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call VecDuplicate(wVec, deltaW, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call VecDuplicate(wVec, baseRes, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Create the two additional work vectors for the line search:
            call VecDuplicate(wVec, g, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call VecDuplicate(wVec, work, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Create Pre-Conditioning Matrix
            allocate (nnzDiagonal(nCellsLocal(1_intType) * nTimeIntervalsSpectral), &
                      nnzOffDiag(nCellsLocal(1_intType) * nTimeIntervalsSpectral))

            if (viscous .and. NK_viscPC) then
                stencil => visc_pc_stencil
                n_stencil = N_visc_pc
            else
                stencil => euler_pc_stencil
                n_stencil = N_euler_pc
            end if

            level = 1
            call statePreAllocation(nnzDiagonal, nnzOffDiag, nDimW / nw, stencil, n_stencil, &
                                    level, .False.)
            call myMatCreate(dRdwPre, nw, nDimW, nDimW, nnzDiagonal, nnzOffDiag, &
                             __FILE__, __LINE__)

            call matSetOption(dRdwPre, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE, ierr)
            call EChk(ierr, __FILE__, __LINE__)
            deallocate (nnzDiagonal, nnzOffDiag)

            ! Setup Matrix-Free dRdw matrix and its function
            call MatCreateMFFD(adflow_comm_world, nDimW, nDimW, &
                               PETSC_DETERMINE, PETSC_DETERMINE, dRdw, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call MatMFFDSetFunction(dRdw, FormFunction_mf, ctx, ierr)
            call EChk(ierr, __FILE__, __LINE__)
            call configureMFFD(dRdw)

            ! Setup a matrix free matrix for drdw
            call MatCreateShell(ADFLOW_COMM_WORLD, nDimW, nDimW, PETSC_DETERMINE, &
                                PETSC_DETERMINE, ctx, dRdwPseudo, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Set the shell operation for doing matrix vector multiplies
            call MatShellSetOperation(dRdwPseudo, MATOP_MULT, NKMatMult, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Setup a matrix free matrix wrapping dRdw with the SA-Gamma-Retheta
            ! Eq. 59 source-dt-restriction diagonal (see NKSrcDtMatMult).
            call MatCreateShell(ADFLOW_COMM_WORLD, nDimW, nDimW, PETSC_DETERMINE, &
                                PETSC_DETERMINE, ctx, dRdwNKSrcDt, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call MatShellSetOperation(dRdwNKSrcDt, MATOP_MULT, NKSrcDtMatMult, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Set the mat_row_oriented option to false so that dense
            ! subblocks can be passed in in fortran column-oriented format
            call MatSetOption(dRdWPre, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call MatSetOption(dRdW, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            if (NK_precondType == 'mg') then
                call setupAMG(drdwpre, nDimW / nw, nw, NK_AMGLevels, NK_AMGNSmooth)
            end if

            !  Create the linear solver context
            call KSPCreate(ADFLOW_COMM_WORLD, NK_KSP, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Set operators for the solver
            call KSPSetOperators(NK_KSP, dRdw, dRdwPre, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            if (useLinResMonitor) then
#if PETSC_VERSION_GE(3,8,0)
                ! This could be wrong. There is no petsc_null_context???
                call KSPMonitorSet(NK_KSP, linearResidualMonitor, PETSC_NULL_FUNCTION, &
                                   PETSC_NULL_FUNCTION, ierr)
#else
                call KSPMonitorSet(NK_KSP, linearResidualMonitor, PETSC_NULL_OBJECT, &
                                   PETSC_NULL_FUNCTION, ierr)
#endif
                call EChk(ierr, __FILE__, __LINE__)
            end if

            NK_solverSetup = .True.
            NK_iter = 0
        end if

    end subroutine setupNKsolver

    subroutine linearResidualMonitor(myKSP, n, rnorm, dummy, ierr)
        use communication, only: myid
        implicit none
        !
        !     Subroutine arguments.
        !
        ! myKsp - Iterative context
        ! n     - Iteration number
        ! rnorm - 2-norm (preconditioned) residual value
        ! dummy - Optional user-defined monitor context (unused here)
        ! ierr  - Return error code

        KSP myKSP
        integer(kind=intType) :: n, dummy, ierr
        real(kind=alwaysRealType) :: rnorm

        ! Write the residual norm to stdout every adjMonStep iterations.
        if (myid == 0) then
            print *, n, rnorm
        end if
        ierr = 0
    end subroutine LinearResidualMonitor

    subroutine NKMatMult(A, vecX, vecY, ierr)

        ! PETSc user-defied call back function for computing the product of
        ! dRdw with a vector. Here we just call the much more broadly
        ! useful routine computeMatrixFreeProductFwd()

        use constants
        use utils, only: EChk
        implicit none

        ! PETSc Arguments
        Mat A
        Vec vecX, vecY
        integer(kind=intType) :: ierr, i, j, k, l, nn, sps, ii
        real(kind=realType) :: dt
        real(kind=realType), pointer :: yPtr(:), xPtr(:)

        ! Frist run the underlying matrix-free mult
        call matMult(dRdw, vecX, vecY, ierr)

        call VecGetArrayF90(vecY, yPtr, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecGetArrayReadF90(vecX, xPtr, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        yPtr = yPtr + one / NK_CFL * xPtr

        call VecRestorearrayF90(vecY, yPtr, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecRestorearrayReadF90(vecX, xPtr, ierr)
        call EChk(ierr, __FILE__, __LINE__)

    end subroutine NKMatMult

    subroutine getFreeStreamResidual(rhoRes, totalRRes)

        use constants
        use blockPointers, only: nDom, ib, jb, kb, w
        use inputTimeSpectral, only: nTimeIntervalsSpectral
        use flowVarRefState, only: nw, winf
        use utils, only: setPointers
        implicit none

        real(kind=realType), intent(out) :: rhoRes, totalRRes
        real(kind=realType), dimension(:), allocatable :: tmp
        integer(kind=intType) :: nDimW, nDimP, counter
        integer(kind=intType) :: nn, sps, i, j, k, l, n

        ! Get the residual cooresponding to the free-stream on the fine
        ! grid-level --- This saves the current values in W, P, rlv and rev
        ! and restores them when finished.

        call getInfoSize(n)
        allocate (tmp(n))
        call getInfo(tmp, n)

        ! Set the w-variables to the ones of the uniform flow field.
        spectralLoop4b: do sps = 1, nTimeIntervalsSpectral
            domains4b: do nn = 1, nDom
                call setPointers(nn, 1, sps)
                do l = 1, nw
                    do k = 0, kb
                        do j = 0, jb
                            do i = 0, ib
                                w(i, j, k, l) = winf(l)
                            end do
                        end do
                    end do
                end do
            end do domains4b
        end do spectralLoop4b

        ! Evaluate the residual now
        call computeResidualNK(useUpdateIntermed=.True.)
        call getCurrentResidual(rhoRes, totalRRes)

        ! Put everything back
        call setInfo(tmp, n)

        deallocate (tmp)

        ! propogate the old values throught the code.
        ! This is not needed for euler and shouldn't be needed for
        ! viscous equations either, but becuase of an issue else where it is.
        ! see https://github.com/mdolab/adflow/pull/46 for the discussion.
        call computeResidualNK(useUpdateIntermed=.True.)

    end subroutine getFreeStreamResidual

    subroutine getCurrentResidual(rhoRes, totalRRes)

        use constants
        use communication, only: adflow_comm_world
        use block, only: nCellGlobal
        use blockPointers, only: nDom
        use inputTimeSpectral, only: nTimeIntervalsSpectral
        use iteration, only: currentLevel
        use monitor, only: monLoc, monGlob, nMonSum
        use utils, only: setPointers, sumResiduals, sumAllResiduals
        implicit none

        ! Compute rhoRes and totalR. The actual residual must have already
        ! been evaluated

        real(kind=realType), intent(out) :: rhoRes, totalRRes
        integer(kind=intType) :: sps, nn, ierr

        monLoc = zero
        do sps = 1, nTimeIntervalsSpectral
            do nn = 1, nDom
                call setPointers(nn, currentLevel, sps)
                call sumResiduals(1, 1) ! Sum 1st state res into first mon location
                call sumAllResiduals(2) ! Sum into second mon location
            end do
        end do

        ! This is the same calc as in convergence info, just for rehoRes and
        ! totalR only.
        call mpi_allreduce(monLoc, monGlob, nMonSum, adflow_real, &
                           mpi_sum, ADflow_comm_world, ierr)

        rhoRes = sqrt(monGlob(1) / nCellGlobal(currentLevel))
        totalRRes = sqrt(monGlob(2))

    end subroutine getCurrentResidual

    subroutine FormJacobianNK

        use constants
        use inputADjoint, only: viscPC
        use inputPhysics, only: turbModel
        use inputIteration, only: transitionSrcDtRestrict, noBacktrackCount, srcDtDeactivateIters, &
                                  transitionNK, transitionNKActive, transitionResidualAutoscale
        use communication, only: myID
        use utils, only: EChk
        use adjointUtils, only: setupStateResidualMatrix, setupStandardKSP, setupStandardMultigrid
        use paramTurb, only: srcLambdaModeFull
        use saGammaRetheta, only: computeSrcLambda
        implicit none

        ! Local Variables
        character(len=maxStringLen) :: preConSide, localPCType, kspObjectType, globalPCType, localOrdering
        integer(kind=intType) :: ierr
        logical :: useAD, usePC, useTranspose, useObjective, tmp
        integer(kind=intType) :: i, j, k, l, ii, nn, sps
        logical :: useCoarseMats

        ! Eq. 58 S_a proxy: refresh the per-block autoscale factor once per
        ! Jacobian reform (same lagged cadence as everything else here),
        ! before any residual/PC work below uses it (via setRVec).
        if (turbModel == spalartallmarasnoft2gammaretheta .and. transitionNK .and. transitionNKActive .and. &
            transitionResidualAutoscale) then
            call computeNKResidualAutoscale()
        end if

        ! Dummy assembly begin/end calls for the matrix-free Matrx
        call MatAssemblyBegin(dRdw, MAT_FINAL_ASSEMBLY, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        call MatAssemblyEnd(dRdw, MAT_FINAL_ASSEMBLY, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Assemble the approximate PC (fine level, level 1)
        useAD = NK_ADPC
        usePC = .True.
        useTranspose = .False.
        useObjective = .False.
        tmp = viscPC ! Save what is in viscPC and set to the NK variable
        viscPC = NK_viscPC

        if (NK_precondType == 'mg') then
            useCoarseMats = .True.
        else
            useCoarseMats = .False.
        end if

        call setupStateResidualMatrix(dRdwPre, useAD, usePC, useTranspose, &
                                      useObjective, .False., 1_intType, useCoarseMats=useCoarseMats)
        ! Reset saved value
        viscPC = tmp

        ! For SA-Gamma-Retheta the NK system is solved in column-scaled
        ! variables (see getNKColScale); make the assembled PC consistent
        ! with the scaled MFFD operator. The MG coarse levels are NOT
        ! scaled — use the (default) ASM preconditioner with this model.
        if (turbModel == spalartallmarasnoft2gammaretheta .and. transitionNK .and. transitionNKActive) then
            call applyNKColumnScaling(dRdwPre)
            if (NK_precondType == 'mg' .and. myID == 0) then
                print *, 'Warning: NK MG coarse levels are not column-scaled for ', &
                    'SA-Gamma-Retheta; use NKGlobalPreconditioner=additive Schwarz'
            end if

            ! Source-dt restriction (P&Z Eq. 59) reactivation-on-backtrack:
            ! see NKStep for the noBacktrackCount update. Added after column
            ! scaling and pre-multiplied by turbResScale, same convention as
            ! ANKStep's timeStepMat (see the comment there) so it is already
            ! consistent with the scaled PC.
            if (transitionSrcDtRestrict .and. (noBacktrackCount < srcDtDeactivateIters)) then
                call computeSrcLambda(srcLambdaModeFull)
                call applyNKSrcDtDiagonal(dRdwPre)
            end if

            ! Solve the same Eq. 59 source-dt-restricted system the PC
            ! above assumes, instead of the exact unregularized Newton
            ! system: dRdwNKSrcDt wraps dRdw and adds the identical
            ! diagonal to the true operator (see NKSrcDtMatMult), gated
            ! dynamically the same way on every call.
            call KSPSetOperators(NK_KSP, dRdwNKSrcDt, dRdwPre, ierr)
            call EChk(ierr, __FILE__, __LINE__)
        end if

        ! Set up KSP Options
        preConSide = 'right'
        localPCType = 'ilu'
        kspObjectType = 'gmres'
        globalPCType = 'asm'
        localOrdering = 'rcm'

        ! Set up the KSP using the same code as used for the adjoint
        if (NK_precondType == 'asm') then
            call setupStandardKSP(NK_KSP, kspObjectType, NK_subSpace, &
                                  preConSide, globalPCType, NK_asmOverlap, NK_outerPreConIts, localPCType, &
                                  localOrdering, NK_iluFill, NK_innerPreConIts)
        else
            call setupStandardMultigrid(NK_KSP, kspObjectType, NK_subSpace, &
                                        preConSide, NK_asmOverlap, NK_outerPreConIts, &
                                        localOrdering, NK_iluFill, NK_innerPreConIts, &
                                        NK_asmOverlapCoarse, NK_iluFillCoarse, NK_innerPreConItsCoarse)
        end if

        ! Don't do iterative refinement
        call KSPGMRESSetCGSRefinementType(NK_KSP, KSP_GMRES_CGS_REFINE_NEVER, ierr)
        call EChk(ierr, __FILE__, __LINE__)

    end subroutine FormJacobianNK

    subroutine FormFunction_mf(ctx, wVec, rVec, ierr)

        ! This is basically a copy of FormFunction, however it has a
        ! different calling sequence from PETSc. It performs the identical
        ! function. This is used for linear solve application for the
        ! aerostructural system pre-conditioner

        use constants
        implicit none

        ! PETSc Variables
        PetscFortranAddr ctx(*)
        Vec wVec, rVec
        integer(kind=intType) :: ierr

        ! This is just a shell routine that runs the more broadly useful
        ! computeResidualNK subroutine

        call setW(wVec)
        call computeResidualNK(useUpdateIntermed=.False.)
        call setRVec(rVec)
        ! We don't check an error here, so just pass back zero
        ierr = 0

    end subroutine FormFunction_mf

    subroutine destroyNKsolver

        ! Destroy all the PETSc objects for the Newton-Krylov solver.

        use constants
        use utils, only: EChk
        use amg, only: destroyAMG
        implicit none
        integer(kind=intType) :: ierr

        if (NK_solverSetup) then

            call MatDestroy(dRdw, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call MatDestroy(dRdwPre, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call MatDestroy(dRdwPseudo, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call MatDestroy(dRdwNKSrcDt, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call VecDestroy(wVec, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call VecDestroy(rVec, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call VecDestroy(deltaW, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call VecDestroy(baseRes, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call VecDestroy(g, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call VecDestroy(work, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call KSPDestroy(NK_KSP, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call destroyAMG()

            NK_solverSetup = .False.
        end if

    end subroutine destroyNKsolver

    subroutine NKStep(firstCall)

        use constants
        use flowVarRefState, only: nw
        use inputPhysics, only: equations, turbModel
        use flowVarRefState, only: nw, nwf
        use inputIteration, only: solverStallDiag, solverStallDiagStep, &
                                  L2conv, transitionSrcDtRestrict, noBacktrackCount, transitionNK, &
                                  transitionNKAutoDisableTol, transitionNKActive, transitionNKStallStepTol, &
                                  transitionNKStallCountTrigger, transitionNKStallRtolCap, nkStallCount
        use iteration, only: approxTotalIts, totalR0, stepMonitor, LinResMonitor, iterType
        use utils, only: EChk
        use communication, only: myid
        use killSignals, only: routineFailed
        implicit none

        ! Input Variables
        logical, intent(in) :: firstCall

        ! Working Variables
        integer(kind=intType) :: iter, ierr, kspIterations
        integer(kind=intType) :: maxNonLinearIts, nfevals, maxIt
        real(kind=alwaysRealType) :: norm, rtol, atol
        real(kind=alwaysrealType) :: fnorm, ynorm, gnorm
        logical :: flag
        real(kind=alwaysRealType) :: resHist(NK_subspace + 1)

        if (firstCall) then
            call setupNKSolver()

            ! Fresh entry into the NK phase: re-arm the transitionNKAutoDisableTol
            ! latch (see inputParam.F90) so it re-evaluates from scratch even if
            ! a previous NK excursion already tripped it off.
            transitionNKActive = .true.
            nkStallCount = 0
            nkGuardBestNorm = huge(1.0_alwaysRealType)

            ! Copy the adflow 'w' into the petsc wVec
            call setwVec(wVec)

            ! Evaluate the residual before we start and put the residual in
            ! 'g', which is what would be the case after a linesearch.
            call computeResidualNK(useUpdateIntermed=.False.)
            call setRVec(rVec)
            iter_k = 1
            iter_m = 0
        else
            NK_iter = NK_iter + 1

            ! Increment counter for the nonmonotne line serach
            iter_k = iter_k + 1
        end if

        ! Compute the norm of rVec for use in EW Criteria
        call VecNorm(rVec, NORM_2, norm, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! One-way latch: once past transitionNKAutoDisableTol (relative to
        ! totalR0), drop the transitionNK bundle for the rest of this NK
        ! phase. transitionNKAutoDisableTol <= 0 (default) never trips this.
        if (transitionNKAutoDisableTol > zero .and. (norm / totalR0) < transitionNKAutoDisableTol) then
            transitionNKActive = .false.
        end if

        ! Determine if if we need to form the Preconditioner
        if (mod(NK_iter, NK_jacobianLag) == 0) then
            NK_CFL = NK_CFL0 * (totalR0 / norm)**1.5
            iterType = "     *NK"
            call FormJacobianNK()
        else

            call MatAssemblyBegin(dRdw, MAT_FINAL_ASSEMBLY, ierr)
            call EChk(ierr, __FILE__, __LINE__)
            call MatAssemblyEnd(dRdw, MAT_FINAL_ASSEMBLY, ierr)
            call EChk(ierr, __FILE__, __LINE__)
            iterType = "      NK"
        end if

        if (NK_iter == 0 .or. .not. NK_useEW) then
            rtol = NK_rtolInit
        else
            call getEWTol(norm, oldNorm, rtolLast, rtol)
        end if

        ! Stall escape: EW picks a loose rtol exactly when stalled (norm
        ! barely changing -> ratio~1 -> rtol->0.8 cap, see getEWTol above).
        ! Force it tighter once Step has been pinned for several iterations
        ! in a row (see inputParam.F90 for the rationale).
        if (turbModel == spalartallmarasnoft2gammaretheta .and. transitionNK .and. transitionNKActive .and. &
            transitionNKStallRtolCap < one .and. nkStallCount >= transitionNKStallCountTrigger) then
            rtol = min(rtol, transitionNKStallRtolCap)
        end if

        ! Save the old rtol and norm for the next iteration
        oldNorm = norm
        rtolLast = rtol

        ! Set all tolerances for linear solver.

        ! The 0.01 multiplier for atol requires some explanation:
        ! The linear residual is roughly the same magnitude
        ! as the nonlinear one at the start of the linear solution,
        ! assuming the initial guess does not have a large effect
        ! on the linear residual. KSPSolve exits when either the
        ! rtol or atol is satisfied, which means that the atol
        ! only comes into play near the end of the nonlinear solution.
        ! We use atol because in the final Newton step when we are
        ! close to the L2 target, we won't need to solve the linear
        ! system tightly. For example, if we are one order of magnitude
        ! away from the nonlinear solver target, then there is
        ! no point in solving the linear system to 8 orders of
        ! magnitude convergence in the linear residual. Instead,
        ! we can stop early using atol. However, in very rare situations,
        ! it can happen that the nonlinear residual is *just* above
        ! the convergence criteria, while the linear residual is
        ! *just* below. What happens is that the linear solver hits
        ! the atol limit immediately, doesn't do anything, and then
        ! the nonlinear convergence check can't do anything either.
        ! By multiplying by 0.01, we make sure that the linear solver
        ! actually has to do *something* and not just exit immediately.

#ifndef USE_COMPLEX
        ! in the real mode, we set the atol slightly lower than the target L2 convergence
        ! as explained in the comment block above
        atol = totalR0 * L2Conv * 0.01_realType
#else
        ! in complex mode, we want to tightly solve the linear system every time
        ! because even though the real residuals converge, complex ones might (and do) lag
        ! this approach makes sure that even with a converged real system, the linear solver
        ! still converges the linear system tightly and this helps with the complex system convergence
        atol = totalR0 * L2Conv * 1e-6_realType
#endif

        maxIt = NK_subspace

        call KSPSetTolerances(NK_KSP, real(rtol), &
                              real(atol), real(NK_divTol), maxIt, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call KSPSetResidualHistory(NK_KSP, resHist, maxIt + 1, PETSC_TRUE, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! set the BaseVector of the matrix-free matrix
        call formFunction_mf(ctx, wVec, baseRes, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        call MatMFFDSetBase(dRdW, wVec, baseRes, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Actually do the Linear Krylov Solve
        call KSPSolve(NK_KSP, rVec, deltaW, ierr)

        ! DON'T just check the error. We want to catch error code 72
        ! which is a floating point error. This is ok, we just reset and
        ! keep going
        if (ierr == 72) then
            ! The convergence check will get the nan
        else
            call EChk(ierr, __FILE__, __LINE__)
        end if

        nfevals = 0
        if (NK_LS == noLineSearch) then
            call LSNone(wVec, rVec, g, deltaW, work, nfevals, flag, stepMonitor)
        else if (NK_LS == cubicLineSearch) then
            call LSCubic(wVec, rVec, g, deltaW, work, fnorm, ynorm, gnorm, &
                         nfevals, flag, stepMonitor)
        else if (NK_LS == nonMonotoneLineSearch) then
            iter_m = min(iter_m + 1, mMax)
            call LSNM(wVec, rVec, g, deltaW, work, fnorm, ynorm, gnorm, &
                      nfevals, flag, stepMonitor)
        end if

        if (.not. flag) then
            routineFailed = .True.
        end if

        ! Algorithm 2 (P&Z 2020 SS IV.B.2): per-node bounds-triggered damping
        ! of gamma/Re-theta-t on the accepted step, before it becomes the new
        ! state. See applyNKAlgorithm2Damping for the full rationale.
        if (turbModel == spalartallmarasnoft2gammaretheta .and. transitionNK .and. transitionNKActive) then
            call applyNKAlgorithm2Damping(wVec, work)
        end if

        ! Copy the work vector to wVec. This is our new state vector
        call VecCopy(work, wVec, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Use the result from the line sesarch for the residual
        call vecCopy(g, rVec, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Update the approximate iteration counter. The +nFevals is for the
        ! iterations taken during the linesearch

        call KSPGetIterationNumber(NK_KSP, kspIterations, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        linResMonitor = resHist(kspIterations + 1) / resHist(1)

        approxTotalIts = approxTotalIts + nfEvals + kspIterations

        ! ============== Source-dt reactivation-on-backtrack (P&Z Eq. 59 / SS4.B.3) ==============
        ! NK has no pseudo-transient term of its own (pure Newton on the steady
        ! residual), so a rejected/backtracked step just repeats the same
        ! direction with no regularization -- unlike the paper, which switches
        ! the source-dt restriction back on when a Newton step needs
        ! backtracking. Mirror ANKStep's coupled-path counter (same
        ! noBacktrackCount/srcDtDeactivateIters globals, see FormJacobianNK):
        ! a backtracked step (stepMonitor < 1) or an outright line-search
        ! failure resets the counter, reactivating the restriction on the PC
        ! for the next srcDtDeactivateIters clean Newton steps. Unlike ANK,
        ! there is no totalR-vs-secondOrdSwitchTol leg here: NK only engages
        ! once the residual is already well past that regime, so a residual
        ! rise is already caught by the backtrack check.
        if (turbModel == spalartallmarasnoft2gammaretheta .and. transitionNK .and. transitionNKActive .and. &
            transitionSrcDtRestrict) then
            if ((.not. flag) .or. stepMonitor < one) then
                noBacktrackCount = 0
            else
                noBacktrackCount = noBacktrackCount + 1
            end if
        end if

        ! Stall detector feeding the rtol cap above: count consecutive
        ! pinned-step iterations.
        if (turbModel == spalartallmarasnoft2gammaretheta .and. transitionNK .and. transitionNKActive) then
            if (stepMonitor < transitionNKStallStepTol) then
                nkStallCount = nkStallCount + 1
            else
                nkStallCount = 0
            end if
        end if

        ! ============== Stall diagnostics (VERIF_06) ==============
        ! NK has no physicality check and no pseudo-transient term (F0/F6), so
        ! the only things that can throttle the step are the cubic line search
        ! and the turbulence-residual pre-limit. Report what actually happened:
        !   lam      : accepted step (minlambda = 0.01 means the search ran out
        !              of room and the step was taken anyway -- the E1 symptom)
        !   nfeval   : residual evaluations spent in the line search
        !   linRes   : linear residual achieved by GMRES this iteration
        !   stall    : consecutive pinned-step iterations
        ! A run showing lam=1.0E-02 with a healthy linRes is step-limited, not
        ! preconditioner-limited -- which is the distinction the whole audit
        ! turns on.
        if (solverStallDiag .and. myid == 0 .and. stepMonitor < solverStallDiagStep) then
            write (*, "(a,i6,a,es9.2,a,i4,a,es9.2,a,l1,a,i4)") &
                " STALLDIAG NK  iter=", NK_iter, &
                "  lam=", stepMonitor, &
                "  nfeval=", nfevals, &
                "  linRes=", linResMonitor, &
                "  LSok=", flag, &
                "  stall=", nkStallCount
        end if

    end subroutine NKStep

    subroutine LSCubic(x, f, g, y, w, fnorm, ynorm, gnorm, nfevals, flag, lambda)

        use constants
        use utils, only: EChk
        use genericISNAN, only: myisnan
        use communication, only: myid
        use initializeFlow, only: setUniformFlow
        use iteration, only: totalR0
        implicit none

        ! Input/Output
        Vec x, f, g, y, w
        !x    - current iterate
        !f    - residual evaluated at x
        !y    - search direction
        !w    - work vector -> On output, new iterate
        !g    - residual evaluated at new iterate y

        real(kind=alwaysrealType) :: fnorm, gnorm, ynorm
        real(kind=realType) :: alpha
        logical :: flag
        integer(kind=intType) :: nfevals
        !   Note that for line search purposes we work with with the related
        !   minimization problem:
        !      min  z(x):  R^n -> R,
        !   where z(x) = .5 * fnorm*fnorm, and fnorm = || f ||_2.
        !

        real(kind=realType) :: initslope, lambdaprev, gnormprev, a, b, d, t1, t2
        real(kind=alwaysRealType) :: minlambda, lambda, lambdatemp
        real(kind=alwaysRealType) :: rellength
        integer(kind=intType) :: ierr, iter
        real(kind=alwaysRealType) :: turbRes1, turbRes2, flowRes1, flowRes2, totalRes1, totalRes2
        logical :: hadANan
        ! Call to get the split norms
        call setRVec(g, flowRes1, turbRes1, totalRes1)

        ! Set some defaults:
        ! alpha relaxed 1e-2 -> 1e-3 (2026-07-16): the production run
        ! (best_strategie/logs/3_NK_paper_faithful.log) pinned step at
        ! minlambda for 1000+ consecutive iterations -- the Armijo test
        ! was essentially never satisfied at any lambda above the floor.
        ! A smaller alpha accepts a step with a weaker guarantee of
        ! decrease, trading some robustness for actually taking full-ish
        ! steps (same idea as relaxing ANKUnsteadyLSTol/ANKPhysicalLSTol
        ! for CANK -- see docs/ADFLOW_BASE/ADFLOW_04_debugging_playbook.md
        ! E1, which explicitly warns "more aggressive rho/E ... risks
        ! negative rho/E"). 1e-4 was tried first and SEGV'd after ~23
        ! iterations (see logs/3b_NK_paper_faithful_alpha_relaxed.log) --
        ! NK has no physicality check at all (unlike ANK), so an
        ! over-permissive alpha can accept a step that drives density/
        ! energy unphysical with nothing to catch it. 1e-3 is the
        ! compromise; revert to 1e-2 if this still diverges/crashes.
        alpha = 1.e-3_realType
        minlambda = .01
        nfevals = 0
        flag = .True.
        lambda = 1.0_realType
        ! Compute the two norms we need:
        call VecNorm(y, NORM_2, ynorm, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecNorm(f, NORM_2, fnorm, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call MatMult(dRdw, y, w, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        nfevals = nfevals + 1

        call VecDot(f, w, initslope, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        if (initslope > zero) then
            initslope = -initslope
        end if

        if (initslope == 0.0_realType) then
            initslope = -1.0_realType
        end if
#ifdef USE_COMPLEX
        call VecWAXPY(w, cmplx(-lambda, 0.0), y, x, ierr)
        call EChk(ierr, __FILE__, __LINE__)
#else
        call VecWAXPY(w, -lambda, y, x, ierr)
        call EChk(ierr, __FILE__, __LINE__)
#endif

        ! Compute Function:
        call setW(w)
        call computeResidualNK(useUpdateIntermed=.True.)
        call setRVec(g, flowRes2, turbRes2, gnorm)

        nfevals = nfevals + 1

        ! Before we get to the actual line search we do two additional
        ! checks:

        ! 1. If the full step has a nan, we backtrack until we get a valid
        ! step. We then lower the NK switch tol such that the solver is
        ! forced back up to ANK or DADI/RK to keep going a bit further.
        !
        ! 2. If the turbulence residual goes up by large-ish factor (2.0),
        ! we pre-limit the step. The reason for this is that a unit step
        ! might lower the total residual, but the turb res could go up an
        ! order of magnitude or more.

        ! Turb-blowup pre-limit threshold relaxed 2.0 -> 5.0 (2026-07-16,
        ! same investigation as the alpha change above): a factor of 2.0
        ! was tripping on essentially every full step in the production
        ! run, forcing the special backtrack path (which itself floors
        ! near minlambda) instead of letting the normal Armijo-based
        ! cubic backtrack (now much more permissive via the lower alpha)
        ! handle it. Dialed back from 5.0 to 3.0 after alpha=1e-4+5.0
        ! together SEGV'd (see alpha comment above).
        ! Retried at 5.0 with alpha still at 1e-3 (2026-07-18,
        ! nk_switch_crossing_test): confirmed unsafe even at alpha=1e-3 --
        ! the very next NK iteration after entry let a step through that
        ! blew nuturb res up to O(1e3) and totalRes to O(1e9), and the
        ! solve never recovered (diverged in ANK for 1000+ iters
        ! afterward). Reverted to 3.0. Do not raise this again without a
        ! physicality check to back it up; the stall investigation should
        ! go through transitionNKAutoDisableTol instead (see inputParam.F90).
        hadANan = .False.
        if (myisnan(gnorm) .or. turbRes2 > 3.0 * turbRes1) then
            ! Special testing for nans

            if (myisnan(gnorm)) then
                hadANan = .True.
                call setUniformFlow()
                lambda = 0.5
            else
                ! Large turb jump
                lambda = lambda * (turbRes1 / turbRes2)
                lambda = max(lambda, 0.1)
            end if

            backtrack: do iter = 1, 10
                ! Compute new x value:
#ifdef USE_COMPLEX
                call VecWAXPY(w, cmplx(-lambda, 0.0), y, x, ierr)
                call EChk(ierr, __FILE__, __LINE__)
#else
                call VecWAXPY(w, -lambda, y, x, ierr)
                call EChk(ierr, __FILE__, __LINE__)
#endif

                ! Compute Function
                call setW(w)
                call computeResidualNK(useUpdateIntermed=.True.)
                call setRVec(g, flowRes2, turbRes2, gnorm)

                nfevals = nfevals + 1

                if (myisnan(gnorm)) then
                    ! Just reset the flow, adjust the step back and keep
                    ! going
                    call setUniformFlow()
                    lambda = lambda*.5
                else

                    ! Sufficient reduction! Whoo! This is great we're done!
                    if (0.5_realType * gnorm * gnorm <= 0.5_realType * fnorm * fnorm + alpha * initslope) then
                        exit
                    end if

                    ! If we're less than min lambda, just take it. This could
                    ! let the residual go up slightly. That's ok.
                    if (lambda < minlambda) then
                        exit
                    end if

                    ! Otherwise, cut back the lambda
                    lambda = lambda * 0.5

                end if
            end do backtrack

            if (hadANan) then
                ! Adjust the NK switch tolerance such that the ANK or DADI
                ! goes a little further.
                nk_switchtol = 0.8 * (gnorm / totalR0)
            end if

            ! All finished with this "pre" line search.
            return
        end if

        ! Sufficient reduction from the basic step. This is the return for
        ! a unit step. This is what we want.
        if (0.5_realType * gnorm * gnorm <= 0.5_realType * fnorm * fnorm + alpha * initslope) then
            goto 100
        end if

        ! Fit points with quadratic
        lambda = 1.0_realType
        lambdatemp = -initslope / (gnorm * gnorm - fnorm * fnorm - 2.0_realType * initslope)
        lambdaprev = lambda
        gnormprev = gnorm
        if (lambdatemp > 0.5_realType * lambda) then
            lambdatemp = 0.5_realType * lambda
        end if

        if (lambdatemp <= .1_realType * lambda) then
            lambda = .1_realType * lambda
        else
            lambda = lambdatemp
        end if

#ifdef USE_COMPLEX
        call VecWAXPY(w, -cmplx(lambda, 0.0), y, x, ierr)
        call EChk(ierr, __FILE__, __LINE__)
#else
        call VecWAXPY(w, -lambda, y, x, ierr)
        call EChk(ierr, __FILE__, __LINE__)
#endif

        ! Compute new function again:
        call setW(w)
        call computeResidualNK(useUpdateIntermed=.True.)
        call setRVec(g)

        nfevals = nfevals + 1

        call VecNorm(g, NORM_2, gnorm, ierr)
        if (ierr == PETSC_ERR_FP) then
            flag = .False.
            return
        end if
        call EChk(ierr, __FILE__, __LINE__)

        ! Sufficient reduction
        if (0.5_realType * gnorm * gnorm <= 0.5_realType * fnorm * fnorm + lambda * alpha * initslope) then
            goto 100
        end if

        ! Fit points with cubic
        cubic_loop: do while (.True.)

            if (lambda <= minlambda) then
                exit cubic_loop
            end if
            t1 = 0.5_realType * (gnorm * gnorm - fnorm * fnorm) - lambda * initslope
            t2 = 0.5_realType * (gnormprev * gnormprev - fnorm * fnorm) - lambdaprev * initslope

            a = (t1 / (lambda * lambda) - t2 / (lambdaprev * lambdaprev)) / (lambda - lambdaprev)
            b = (-lambdaprev * t1 / (lambda * lambda) + lambda * t2 / (lambdaprev * lambdaprev)) / (lambda - lambdaprev)
            d = b * b - three * a * initslope
            if (d < 0.0_realType) then
                d = 0.0_realType
            end if

            if (a == 0.0_realType) then
                lambdatemp = -initslope / (2.0_realType * b)
            else
                lambdatemp = (-b + sqrt(d)) / (3.0_realType * a)
            end if

            lambdaprev = lambda
            gnormprev = gnorm

            if (lambdatemp > 0.5_realType * lambda) then
                lambdatemp = 0.5_realType * lambda
            end if
            if (lambdatemp <= .1_realType * lambda) then
                lambda = .1_realType * lambda
            else
                lambda = lambdatemp
            end if

            if (myisnan(lambda)) then
                flag = .False.
                exit cubic_loop
            end if

#ifdef USE_COMPLEX
            call VecWAXPY(w, cmplx(-lambda, 0.0), y, x, ierr)
            call EChk(ierr, __FILE__, __LINE__)
#else
            call VecWAXPY(w, -lambda, y, x, ierr)
            call EChk(ierr, __FILE__, __LINE__)
#endif
            ! Compute new function again:
            call setW(w)
            call computeResidualNK(useUpdateIntermed=.True.)
            call setRVec(g)
            nfevals = nfevals + 1

            call VecNorm(g, NORM_2, gnorm, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Is reduction enough?
            if (0.5_realType * gnorm * gnorm <= 0.5_realType * fnorm * fnorm + lambda * alpha * initslope) then
                exit cubic_loop
            end if
        end do cubic_loop

100     continue

    end subroutine LSCubic

    subroutine LSNone(x, f, g, y, w, nfevals, flag, step)

        use constants
        use utils, only: EChk
        use communication
        use genericISNAN, only: myisnan
        use inputPhysics, only: turbModel
        use inputIteration, only: transitionNK, transitionNKActive
        implicit none

        ! Input/Output
        Vec x, f, g, y, w
        !x    - current iterate
        !f    - residual evaluated at x
        !y    - search direction
        !w    - work vector -> On output, new iterate
        !g    - residual evaluated at new iterate y

        integer(kind=intType) :: nfevals
        integer(kind=intType) :: ierr, try
        logical :: flag, useGuard
        real(kind=alwaysRealType) :: step
        real(kind=realType) :: tmp
        real(kind=alwaysRealType) :: fnorm, gnorm

        ! Guard constants (SA-GR only, see below). growFac bounds the
        ! per-step residual growth; riseCap bounds the total drift above the
        ! best norm seen this NK phase -- the fine-S809 campaign showed the
        ! unguarded fixed step failing by exactly this slow drift (0.11 -> 9
        ! over ~110 iterations at ~2%/step, each step individually harmless,
        ! ending in NaN), which a per-step check alone cannot catch. A
        ! rejected step halves and retries, and the resulting step < 1
        ! feeds the existing Eq. 59 srcDt-reactivation-on-backtrack logic in
        ! NKStep -- the paper's own globalization -- which a constant full
        ! step can never trigger.
        real(kind=alwaysRealType), parameter :: growFac = 1.5_alwaysRealType
        real(kind=alwaysRealType), parameter :: riseCap = 10.0_alwaysRealType
        integer(kind=intType), parameter :: maxHalvings = 6

        flag = .True.
        nfevals = 0
        step = nk_fixedStep

        ! Stock behavior (accept the fixed step blindly) for everything
        ! except the SA-GR transition model with transitionNK active.
        useGuard = (turbModel == spalartallmarasnoft2gammaretheta .and. &
                    transitionNK .and. transitionNKActive)

        if (.not. useGuard) then
            tmp = -step
            call VecWAXPY(w, tmp, y, x, ierr)
            call EChk(ierr, __FILE__, __LINE__)
            call setW(w)
            call computeResidualNK(useUpdateIntermed=.True.)
            call setRVec(g)
            nfevals = nfevals + 1
            return
        end if

        call VecNorm(f, NORM_2, fnorm, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        do try = 0, maxHalvings
            tmp = -step
            call VecWAXPY(w, tmp, y, x, ierr)
            call EChk(ierr, __FILE__, __LINE__)
            call setW(w)
            call computeResidualNK(useUpdateIntermed=.True.)
            call setRVec(g)
            nfevals = nfevals + 1

            call VecNorm(g, NORM_2, gnorm, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            if (.not. myisnan(gnorm)) then
                if (gnorm <= growFac * fnorm .and. &
                    gnorm <= riseCap * nkGuardBestNorm) exit
                ! Last resort: accept the smallest-step (least-damage)
                ! non-NaN state rather than fail the whole solve.
                if (try == maxHalvings) exit
            else if (try == maxHalvings) then
                ! NaN survives even the smallest step: restore the current
                ! iterate (zero step) so the state stays finite, and report
                ! line-search failure so the caller can fall back.
                call VecCopy(x, w, ierr)
                call EChk(ierr, __FILE__, __LINE__)
                call setW(w)
                call computeResidualNK(useUpdateIntermed=.True.)
                call setRVec(g)
                step = zero
                flag = .False.
                if (myid == 0) then
                    print *, 'LSNone guard: NaN persists at step ', &
                        nk_fixedStep / real(2**maxHalvings, alwaysRealType), &
                        '; step rejected.'
                end if
                return
            end if
            step = step * half
        end do

        nkGuardBestNorm = min(nkGuardBestNorm, gnorm)
    end subroutine LSNone

    subroutine LSNM(x, f, g, y, w, fnorm, ynorm, gnorm, nfevals, flag, step)

        use constants
        use utils, only: EChk
        implicit none

        ! Input/Output
        Vec x, f, g, y, w
        !x    - current iterate
        !f    - residual evaluated at x
        !y    - search direction
        !w    - work vector -> On output, new iterate
        !g    - residual evaluated at new iterate y

        real(kind=alwaysRealType) :: fnorm, gnorm, ynorm
        real(kind=realType) :: alpha
        real(kind=alwaysRealType) :: step
        logical :: flag
        integer(kind=intType) :: nfevals
        !   Note that for line search purposes we work with with the related
        !   minimization problem:
        !      min  z(x):  R^n -> R,
        !   where z(x) = .5 * fnorm*fnorm, and fnorm = || f ||_2.
        !
        real(kind=realType) :: initslope, gamma, sigma, max_val
        integer(kind=intType) :: ierr, iter, j

        ! Set some defaults:
        gamma = 1e-3_realType
        sigma = 0.5_realType

        nfevals = 0
        flag = .True.

        ! Compute the two norms we need:
        call VecNorm(y, NORM_2, ynorm, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecNorm(f, NORM_2, fnorm, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        NKLSFuncEvals(iter_k) = 0.5_realType * fnorm * fnorm

        call MatMult(dRdw, y, w, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        nfevals = nfevals + 1

        call VecDot(f, w, initslope, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        if (initslope > 0.0_realType) then
            initslope = -initslope
        end if

        if (initslope == 0.0_realType) then
            initslope = -1.0_realType
        end if

        alpha = 1.0 ! Initial step length:
        backtrack: do iter = 1, 10

            ! Compute new x value:
            call VecWAXPY(w, -alpha, y, x, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Compute Function @ new x (w is the work vector
            call setW(w)
            call computeResidualNK(useUpdateIntermed=.True.)
            call setRVec(g)
            nfevals = nfevals + 1

            ! Compute the norm at the new trial location
            call VecNorm(g, NORM_2, gnorm, ierr)
            if (ierr == PETSC_ERR_FP) then ! Error code 72 floating point error
                ! Just apply the step limit and keep going (back to the loop start)
                alpha = alpha * sigma
            else
                call EChk(ierr, __FILE__, __LINE__)

                max_val = NKLSFuncEvals(iter_k) + alpha * gamma * initSlope

                ! Loop over the previous, m function values and find the max:
                do j = iter_k - 1, iter_k - iter_m + 1, -1
                    max_val = max(max_val, NKLSFuncEvals(j) + alpha * gamma * initSlope)
                end do

                ! Sufficient reduction
                if (0.5_realType * gnorm * gnorm <= max_val) then
                    exit backtrack
                else
                    alpha = alpha * sigma
                end if
            end if
        end do backtrack
        step = alpha
    end subroutine LSNM

    subroutine computeResidualNK(useUpdateIntermed)

        use constants
        use blockette, only: blocketteRes
        implicit none

        logical, intent(in), optional :: useUpdateIntermed
        logical :: updateIntermed

        ! Only update the time step if explicitly requested
        updateIntermed = .false.

        if (present(useUpdateIntermed)) then
            updateIntermed = useUpdateIntermed
        end if

        ! Shell function to maintain backward compatibility with code using computeResidualNK
        call blocketteRes(useUpdateIntermed=updateIntermed)

    end subroutine computeResidualNK

    subroutine applyPC(in_vec, out_vec, ndof)

        ! Apply the NK PC to the in_vec. This subroutine is ONLY used as a
        ! preconditioner for a global Aero-Structural Newton-Krylov Method

        use constants
        use utils, only: EChk

        implicit none

        ! Input/Output
        integer(kind=intType) :: ndof
        real(kind=realType), dimension(ndof), intent(in) :: in_vec
        real(kind=realTYpe), dimension(ndof), intent(inout) :: out_vec

        ! Working Variables
        integer(kind=intType) :: ierr

        ! Setup the NKsolver if not already done so
        if (.not. NK_solverSetup) then
            call setupNKSolver
        end if

        ! We possibly need to re-form the jacobian
        if (mod(NK_iter, NK_jacobianLag) == 0) then
            call FormJacobianNK()
        end if

        ! Place the two arrays into two vectos. We reuse 'work' and 'g'.
        call VecPlaceArray(work, in_vec, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecPlaceArray(g, out_vec, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Set the base vec
        call setwVec(wVec)

        ! Set the base vec
        call setwVec(wVec)
        call formFunction_mf(ctx, wVec, baseRes, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        call MatMFFDSetBase(dRdW, wVec, baseRes, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        ! This needs to be a bit better...
        call KSPSetTolerances(NK_KSP, 1e-8, 1e-16, 10.0, &
                              applyPCSubSpaceSize, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Actually do the Linear Krylov Solve
        call KSPSolve(NK_KSP, work, g, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Reset the array pointers:
        call VecResetArray(work, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecResetArray(g, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        NK_iter = NK_iter + 1

    end subroutine applyPC

    subroutine applyAdjointPC(in_vec, out_vec, ndof)

        ! Apply the Adjoint PC to the in_vec. This subroutine is ONLY used as a
        ! preconditioner for a global Aero-Structural Krylov Method

        use constants
        use ADjointPETSc, only: adjointKSP, KSP_NORM_NONE, PETSC_DEFAULT_REAL, &
                                psi_like1, psi_like2
        use inputAdjoint, only: applyAdjointPCSubSpaceSize
        use utils, only: EChk
        implicit none

        ! Input/Output
        integer(kind=intType) :: ndof
        real(kind=realType), dimension(ndof), intent(in) :: in_vec
        real(kind=realTYpe), dimension(ndof), intent(inout) :: out_vec

        ! Working Variables
        integer(kind=intType) :: ierr

        ! Hijack adjoint and adjointRes with in_vec and out_vec
        call VecPlaceArray(psi_like1, in_vec, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecPlaceArray(psi_like2, out_vec, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Set KSP_NORM Type to none. Implictly turns off convergence
        ! check. Since we just want to run a fixed number of iterations this
        ! is fine. The should be set regardless of the KSPType.

        call KSPSetNormType(adjointKSP, KSP_NORM_NONE, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! This needs to be a bit better...
        call KSPSetTolerances(adjointKSP, PETSC_DEFAULT_REAL, &
                              PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, &
                              applyAdjointPCSubSpaceSize, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Actually do the Linear Krylov Solve
        call KSPSolve(adjointKSP, psi_like1, psi_like2, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Reset the array pointers:
        call VecResetArray(psi_like1, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecResetArray(psi_like2, ierr)
        call EChk(ierr, __FILE__, __LINE__)

    end subroutine applyAdjointPC


    subroutine configureMFFD(mat)
        ! VERIF_06 F9: apply the MFFD differencing controls. PETSc's default is
        ! Dennis-Schnabel with an assumed function error of machine epsilon,
        ! which is wrong for a residual carrying smoothMinMax blends,
        ! empirical correlations, tanh and clipping -- h then comes out far too
        ! small and J*a is dominated by cancellation. Defaults leave PETSc
        ! exactly as before.
        use constants
        use inputIteration, only: mffdFunctionError, mffdType
        use utils, only: EChk
        implicit none
        Mat mat
        integer(kind=intType) :: ierr

        if (trim(mffdType) == 'wp') then
            call MatMFFDSetType(mat, MATMFFD_WP, ierr)
            call EChk(ierr, __FILE__, __LINE__)
        end if

        if (mffdFunctionError > zero) then
            call MatMFFDSetFunctionError(mat, mffdFunctionError, ierr)
            call EChk(ierr, __FILE__, __LINE__)
        end if

    end subroutine configureMFFD

    subroutine getNKColScale(cs)
        ! Column scale for the NK state vector: one for the mean-flow
        ! entries, turbResScale for the turbulence entries. For models other
        ! than SA-Gamma-Retheta every factor is one, so wVec is the raw state
        ! and all NK behavior is unchanged. For SA-Gamma-Retheta the state
        ! components span ~13 orders of magnitude and the single MFFD
        ! differencing step makes the nuTilde/gamma Jacobian columns FD noise
        ! unless the vector is scaled to O(1) per component.
        use constants
        use inputPhysics, only: turbModel
        use inputIteration, only: turbResScale, transitionNK, transitionNKActive
        use flowVarRefState, only: nw, nt1, nt2
        implicit none
        real(kind=realType), intent(out) :: cs(1:nw)
        integer(kind=intType) :: l

        cs = one
        if (turbModel == spalartallmarasnoft2gammaretheta .and. transitionNK .and. transitionNKActive) then
            do l = nt1, nt2
                cs(l) = turbResScale(l - nt1 + 1)
            end do
        end if
    end subroutine getNKColScale

    subroutine applyNKColumnScaling(matrix)
        ! Right-multiply the assembled NK preconditioner by diag(1/cs) so it
        ! is consistent with the column-scaled state vector. Only called for
        ! SA-Gamma-Retheta.
        use constants
        use flowVarRefState, only: nw
        use utils, only: EChk
        implicit none

        Mat matrix
        Vec colVec
        integer(kind=intType) :: ierr, jj, lState
        real(kind=realType), pointer :: cvec_pointer(:)
        real(kind=realType) :: cs(1:nw)

        call getNKColScale(cs)

        call VecDuplicate(wVec, colVec, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecGetArrayF90(colVec, cvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        do jj = 0, size(cvec_pointer) / nw - 1
            do lState = 1, nw
                cvec_pointer(jj * nw + lState) = one / cs(lState)
            end do
        end do
        call VecRestoreArrayF90(colVec, cvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call MatDiagonalScale(matrix, PETSC_NULL_VEC, colVec, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecDestroy(colVec, ierr)
        call EChk(ierr, __FILE__, __LINE__)

    end subroutine applyNKColumnScaling

    subroutine applyNKSrcDtDiagonal(matrix)
        ! Add the P&Z Eq. 59 source-term dt-restriction diagonal
        ! (srcLambda/transitionSrcDtLimit) to the transition rows of the
        ! assembled, column-scaled NK preconditioner. Row values are
        ! pre-multiplied by turbResScale, the same convention used by
        ! ANKStep's timeStepMat (see the comment at its MatAXPY into
        ! dRdwPre) so the added diagonal is already consistent with the
        ! scaled PC without needing further scaling here. Only called for
        ! SA-Gamma-Retheta, and only when the caller (FormJacobianNK) has
        ! determined the restriction is active and called computeSrcLambda.
        use constants
        use blockPointers, only: nDom, il, jl, kl, srcLambda
        use inputTimeSpectral, only: nTimeIntervalsSpectral
        use inputIteration, only: transitionSrcDtLimit, turbResScale
        use flowVarRefState, only: nw, nt1, nt2
        use utils, only: EChk, setPointers
        implicit none

        Mat matrix
        Vec srcDtVec
        integer(kind=intType) :: ierr, nn, sps, i, j, k, l, l1, ii
        real(kind=realType), pointer :: svec_pointer(:)

        call VecDuplicate(wVec, srcDtVec, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecZeroEntries(srcDtVec, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecGetArrayF90(srcDtVec, svec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        ii = 1
        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                call setPointers(nn, 1_intType, sps)
                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il
                            do l = 1, nw
                                if (l >= nt1 .and. l <= nt2) then
                                    l1 = l - nt1 + 1
                                    svec_pointer(ii) = (srcLambda(i, j, k, l1) / transitionSrcDtLimit) &
                                                       * turbResScale(l1)
                                end if
                                ii = ii + 1
                            end do
                        end do
                    end do
                end do
            end do
        end do
        call VecRestoreArrayF90(srcDtVec, svec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call MatDiagonalSet(matrix, srcDtVec, ADD_VALUES, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecDestroy(srcDtVec, ierr)
        call EChk(ierr, __FILE__, __LINE__)

    end subroutine applyNKSrcDtDiagonal

    subroutine NKSrcDtMatMult(A, vecX, vecY, ierr)
        ! Shell matrix-vector product: dRdw*x with the P&Z Eq. 59
        ! source-dt-restriction diagonal (srcLambda/transitionSrcDtLimit)
        ! added on the transition rows, mirroring applyNKSrcDtDiagonal but
        ! applied directly to the Krylov vector instead of the assembled
        ! PC matrix.
        !
        ! Previously this restriction was only added to dRdwPre (the PC),
        ! never to the actual operator GMRES solves (dRdw). That meant NK
        ! was still solving the exact, unregularized Newton system at
        ! every iteration; the PC only made that exact system easier to
        ! solve; it did not change the solution GMRES converged to. At the
        ! stiff transition-front cells this produced enormous raw updates
        ! that only got bounded post-hoc by applyNKAlgorithm2Damping's
        ! hard clip -- a nonsmooth correction invisible to the next
        ! Jacobian/residual evaluation, which reproduced a residual spike
        ! at the same cells and collapsed the line-search step. Adding the
        ! same diagonal here makes the true operator consistent with the
        ! PC, so the *solution* of the linear system is itself restricted,
        ! not just easier to compute.
        !
        ! Same row/column-scaled space as dRdw already operates in (see
        ! setW/setRVec and the applyNKSrcDtDiagonal comment) -- no
        ! additional scaling needed. Gated by the same dynamic condition
        ! used for the PC (see FormJacobianNK) so the two stay matched on
        ! every call, not just at the last Jacobian reform.
        use constants
        use blockPointers, only: nDom, il, jl, kl, srcLambda
        use inputTimeSpectral, only: nTimeIntervalsSpectral
        use inputIteration, only: transitionSrcDtRestrict, transitionSrcDtLimit, &
                                  turbResScale, noBacktrackCount, srcDtDeactivateIters
        use flowVarRefState, only: nw, nt1, nt2
        use utils, only: EChk, setPointers
        implicit none

        Mat A
        Vec vecX, vecY
        integer(kind=intType) :: ierr, nn, sps, i, j, k, l, l1, ii
        real(kind=realType), pointer :: yPtr(:), xPtr(:)

        call matMult(dRdw, vecX, vecY, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        if (transitionSrcDtRestrict .and. (noBacktrackCount < srcDtDeactivateIters)) then
            call VecGetArrayF90(vecY, yPtr, ierr)
            call EChk(ierr, __FILE__, __LINE__)
            call VecGetArrayReadF90(vecX, xPtr, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ii = 1
            do nn = 1, nDom
                do sps = 1, nTimeIntervalsSpectral
                    call setPointers(nn, 1_intType, sps)
                    do k = 2, kl
                        do j = 2, jl
                            do i = 2, il
                                do l = 1, nw
                                    if (l >= nt1 .and. l <= nt2) then
                                        l1 = l - nt1 + 1
                                        yPtr(ii) = yPtr(ii) + (srcLambda(i, j, k, l1) / transitionSrcDtLimit) &
                                                   * turbResScale(l1) * xPtr(ii)
                                    end if
                                    ii = ii + 1
                                end do
                            end do
                        end do
                    end do
                end do
            end do

            call VecRestoreArrayF90(vecY, yPtr, ierr)
            call EChk(ierr, __FILE__, __LINE__)
            call VecRestoreArrayReadF90(vecX, xPtr, ierr)
            call EChk(ierr, __FILE__, __LINE__)
        end if

    end subroutine NKSrcDtMatMult

    subroutine applyNKAlgorithm2Damping(x, work)
        ! Paper Algorithm 2 (P&Z 2020 SS IV.B.2): per-node, per-variable
        ! exponential back-off damping of the gamma and Re-theta-t updates
        ! only (NOT nu-tilde). Applied after the global NK line search
        ! accepts `work` as the new state and before it is copied into
        ! wVec (see NKStep). Mirrors the DD-ADI damping loop in
        ! saGammaReThetaSolve (saGammaRetheta.F90) exactly: same bounds
        ! (rsaGRgammaLo/Hi, rsaGRreThetaLo), same
        ! transitionDampTheta/transitionDampMaxIter options -- reused, not
        ! reinvented.
        !
        ! Operates directly in the column-scaled space NK already uses
        ! (see getNKColScale): the per-node update x->work is a linear
        ! interpolation, so damping it by a scalar factor is scale-
        ! invariant; only the *bounds check* needs the physical value,
        ! obtained by dividing by the column scale.
        !
        ! Per the paper this is a rarely-firing safety valve (Eq. 59
        ! keeps updates bounded in practice). Like DD-ADI's damping loop,
        ! no residual re-evaluation is done after damping here -- the
        ! next NKStep's residual evaluation picks up the corrected state
        ! naturally. Same trade-off already accepted for DD-ADI.
        use constants
        use paramTurb, only: rsaGRgammaLo, rsaGRgammaHi, rsaGRreThetaLo
        use inputIteration, only: transitionDampTheta, transitionDampMaxIter
        use flowVarRefState, only: nw, nt1
        use communication, only: myid
        use utils, only: EChk
        implicit none

        Vec x, work
        integer(kind=intType) :: ierr, jj, nCells, mm
        integer(kind=intType) :: nDampCapGamma, nDampCapReTheta
        integer(kind=intType) :: nSoftDampGamma, nSoftDampReTheta
        real(kind=realType) :: minDampFactorGamma, minDampFactorReTheta
        real(kind=realType), pointer :: xPtr(:), workPtr(:)
        real(kind=realType) :: cs(1:nw)
        real(kind=realType) :: xOld, deltaScaled, dampFactor, candScaled, candPhys
        integer(kind=intType) :: gammaOff, rethetaOff

        call getNKColScale(cs)
        gammaOff = nt1 + 1
        rethetaOff = nt1 + 2

        call VecGetArrayReadF90(x, xPtr, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        call VecGetArrayF90(work, workPtr, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        nDampCapGamma = 0
        nDampCapReTheta = 0
        ! Diagnostic-only counters (2026-07-19, nk_switch_crossing_test stall
        ! investigation): the hard clip above (nDampCapGamma/ReTheta) never
        ! fires in practice, but the *soft* exponential back-off runs
        ! silently every iteration with no visibility -- it happens AFTER
        ! the line search, so a cell sitting at its gamma/Re-theta-t bound
        ! could have its accepted update crushed back to near-zero net
        ! change regardless of how good the search direction was, and
        ! neither the printed Step nor Lin Res columns would show it. These
        ! track how many cells needed ANY back-off (dampFactor < 1) and the
        ! worst-case factor, printed alongside the existing clip warning.
        nSoftDampGamma = 0
        nSoftDampReTheta = 0
        minDampFactorGamma = one
        minDampFactorReTheta = one
        nCells = size(workPtr) / nw

        do jj = 0, nCells - 1
            ! Gamma: exponential back-off until in [gammaLo, gammaHi]
            xOld = xPtr(jj * nw + gammaOff)
            deltaScaled = workPtr(jj * nw + gammaOff) - xOld
            dampFactor = one
            candScaled = xOld + dampFactor * deltaScaled
            do mm = 1, transitionDampMaxIter
                candPhys = candScaled / cs(gammaOff)
                if (candPhys >= rsaGRgammaLo .and. candPhys <= rsaGRgammaHi) exit
                dampFactor = dampFactor * transitionDampTheta
                candScaled = xOld + dampFactor * deltaScaled
            end do
            candPhys = candScaled / cs(gammaOff)
            if (candPhys < rsaGRgammaLo .or. candPhys > rsaGRgammaHi) then
                nDampCapGamma = nDampCapGamma + 1
                candPhys = min(max(candPhys, rsaGRgammaLo), rsaGRgammaHi)
                candScaled = candPhys * cs(gammaOff)
            end if
            if (dampFactor < one) then
                nSoftDampGamma = nSoftDampGamma + 1
                minDampFactorGamma = min(minDampFactorGamma, dampFactor)
            end if
            workPtr(jj * nw + gammaOff) = candScaled

            ! Re-theta-t: exponential back-off until >= rsaGRreThetaLo (lower bound only)
            xOld = xPtr(jj * nw + rethetaOff)
            deltaScaled = workPtr(jj * nw + rethetaOff) - xOld
            dampFactor = one
            candScaled = xOld + dampFactor * deltaScaled
            do mm = 1, transitionDampMaxIter
                candPhys = candScaled / cs(rethetaOff)
                if (candPhys >= rsaGRreThetaLo) exit
                dampFactor = dampFactor * transitionDampTheta
                candScaled = xOld + dampFactor * deltaScaled
            end do
            candPhys = candScaled / cs(rethetaOff)
            if (candPhys < rsaGRreThetaLo) then
                nDampCapReTheta = nDampCapReTheta + 1
                candScaled = rsaGRreThetaLo * cs(rethetaOff)
            end if
            if (dampFactor < one) then
                nSoftDampReTheta = nSoftDampReTheta + 1
                minDampFactorReTheta = min(minDampFactorReTheta, dampFactor)
            end if
            workPtr(jj * nw + rethetaOff) = candScaled
        end do

        call VecRestoreArrayReadF90(x, xPtr, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        call VecRestoreArrayF90(work, workPtr, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        if (nDampCapGamma + nDampCapReTheta > 0) then
            print *, 'Warning: NK Algorithm 2 damping exhausted transitionDampMaxIter (', &
                transitionDampMaxIter, ') in ', nDampCapGamma, ' gamma / ', &
                nDampCapReTheta, ' reTheta cells on proc ', myid, '; values clipped.'
        end if

        if (nSoftDampGamma + nSoftDampReTheta > 0) then
            print *, 'NK Algorithm 2 soft damping: ', nSoftDampGamma, ' gamma cells (worst factor ', &
                minDampFactorGamma, '), ', nSoftDampReTheta, ' reTheta cells (worst factor ', &
                minDampFactorReTheta, ') on proc ', myid
        end if

    end subroutine applyNKAlgorithm2Damping


    subroutine computeNKResidualAutoscale()
        ! Eq. 58 (P&Z 2020) residual-autoscaling (S_a) proxy. The paper
        ! gives no formula for S_a (cites Osusky & Zingg's thesis,
        ! unavailable here) -- this is a same-intent proxy, NOT verified
        ! identical to their method: periodically (called once per
        ! FormJacobianNK, same cadence as NK_jacobianLag) measure the
        ! row-scaled (S_r-scaled, if transitionRowVolScale is on) residual
        ! norm of the mean-flow block and of each turbulence variable
        ! separately, then set nkAutoScaleFac so every turbulence
        ! variable's block is rescaled to match the current mean-flow
        ! block's magnitude -- preventing one equation's contribution
        ! from vanishing or dominating the combined norm GMRES sees.
        ! Gated by transitionResidualAutoscale (default False, own switch
        ! -- see inputParam.F90).
        use constants
        use blockPointers, only: nDom, volRef, il, jl, kl, dw
        use inputTimeSpectral, only: nTimeIntervalsSpectral
        use flowVarRefState, only: nwf, nt1, nt2
        use inputIteration, only: turbResScale, transitionRowVolScale, transitionNK, transitionNKActive, nkAutoScaleFac
        use inputPhysics, only: turbModel
        use utils, only: setPointers, EChk
        use communication, only: adflow_comm_world
        implicit none

        integer(kind=intType) :: ierr, nn, sps, i, j, k, l
        real(kind=realType) :: ovv, flowRowFac, turbRowFac, tmp
        real(kind=realType) :: sumLocal(4), sumGlobal(4)
        real(kind=realType), parameter :: floorNorm = 1.0e-30_realType
        logical :: useRowVolScale

        sumLocal = zero
        useRowVolScale = transitionNK .and. transitionNKActive .and. transitionRowVolScale .and. &
                        turbModel == spalartallmarasnoft2gammaretheta

        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                call setPointers(nn, 1_intType, sps)
                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il
                            ovv = one / volRef(i, j, k)
                            if (useRowVolScale) then
                                flowRowFac = volRef(i, j, k)**(5.0_realType / 3.0_realType)
                                turbRowFac = volRef(i, j, k)**(2.0_realType / 3.0_realType)
                            else
                                flowRowFac = one
                                turbRowFac = one
                            end if
                            do l = 1, nwf
                                tmp = dw(i, j, k, l) * ovv * flowRowFac
                                sumLocal(1) = sumLocal(1) + tmp**2
                            end do
                            do l = nt1, nt2
                                tmp = dw(i, j, k, l) * ovv * turbResScale(l - nt1 + 1) * turbRowFac
                                sumLocal(1 + (l - nt1 + 1)) = sumLocal(1 + (l - nt1 + 1)) + tmp**2
                            end do
                        end do
                    end do
                end do
            end do
        end do

        call mpi_allreduce(sumLocal, sumGlobal, 4, adflow_real, mpi_sum, adflow_comm_world, ierr)
        sumGlobal = sqrt(sumGlobal)

        nkAutoScaleFac(1) = one
        do l = 2, 4
            if (sumGlobal(l) > floorNorm) then
                nkAutoScaleFac(l) = sumGlobal(1) / sumGlobal(l)
            else
                nkAutoScaleFac(l) = one
            end if
        end do

    end subroutine computeNKResidualAutoscale

    subroutine setWVec(wVec)

        ! Set the current state in w into the PETSc Vector, column-scaled
        ! (see getNKColScale; scale is one except for SA-Gamma-Retheta)

        use constants
        use blockPointers, only: nDom, il, jl, kl, w
        use inputtimespectral, only: ntimeIntervalsSpectral
        use flowvarrefstate, only: nw
        use utils, only: setPointers, EChk
        implicit none

        Vec wVec
        integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii
        real(kind=realType), pointer :: wvec_pointer(:)
        real(kind=realType) :: cs(1:nw)

        call getNKColScale(cs)

        call VecGetArrayF90(wVec, wvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        ii = 1
        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                call setPointers(nn, 1_intType, sps)
                ! Copy off w to wVec
                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il
                            do l = 1, nw
                                wvec_pointer(ii) = w(i, j, k, l) * cs(l)
                                ii = ii + 1
                            end do
                        end do
                    end do
                end do
            end do
        end do

        call VecRestoreArrayF90(wVec, wvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

    end subroutine setWVec

    subroutine setRVec(rVec, flowRes, turbRes, totalRes)

        ! Set the current residual in dw into the PETSc Vector
        use constants
        use blockPointers, only: nDom, volRef, il, jl, kl, dw
        use inputtimespectral, only: nTimeIntervalsSpectral
        use flowvarrefstate, only: nw, nwf, nt1, nt2
        use inputPhysics, only: turbModel
        use inputIteration, only: turbResScale, transitionNK, transitionNKActive, transitionRowVolScale, &
                                  transitionResidualAutoscale, nkAutoScaleFac
        use utils, only: setPointers, EChk
        use communication, only: adflow_comm_world
        implicit none

        Vec rVec
        integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii
        real(kind=realType), pointer :: rvec_pointer(:)
        real(Kind=realType) :: ovv, flowRowFac, turbRowFac
        real(kind=alwaysRealType), intent(out), optional :: flowRes, turbRes, totalRes
        real(kind=realType) :: tmp, tmp2(2), flowResLocal, turbResLocal
        logical :: useRowVolScale, useAutoscale

        flowResLocal = zero
        turbResLocal = zero
        useRowVolScale = transitionNK .and. transitionNKActive .and. transitionRowVolScale .and. &
                         turbModel == spalartallmarasnoft2gammaretheta
        useAutoscale = transitionNK .and. transitionNKActive .and. transitionResidualAutoscale .and. &
                       turbModel == spalartallmarasnoft2gammaretheta

        call VecGetArrayF90(rVec, rvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        ii = 1

        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                call setPointers(nn, 1_intType, sps)
                ! Copy off dw/vol to rVec
                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il
                            ovv = 1 / volRef(i, j, k)
                            ! Eq. 58 (P&Z) geometric row scaling, additional to the
                            ! existing 1/volRef: brings flow rows toward J^{2/3} and
                            ! turb rows toward J^{-1/3} (see transitionRowVolScale
                            ! comment, inputParam.F90 -- unverified exponent match).
                            if (useRowVolScale) then
                                flowRowFac = volRef(i, j, k)**(5.0_realType / 3.0_realType)
                                turbRowFac = volRef(i, j, k)**(2.0_realType / 3.0_realType)
                            else
                                flowRowFac = one
                                turbRowFac = one
                            end if
                            do l = 1, nwf
                                tmp = dw(i, j, k, l) * ovv * flowRowFac
                                if (useAutoscale) tmp = tmp * nkAutoScaleFac(1)
                                rvec_pointer(ii) = tmp
                                ii = ii + 1
                                flowResLocal = flowResLocal + tmp**2
                            end do
                            do l = nt1, nt2
                                tmp = dw(i, j, k, l) * ovv * turbResScale(l - nt1 + 1) * turbRowFac
                                if (useAutoscale) tmp = tmp * nkAutoScaleFac(1 + (l - nt1 + 1))
                                rvec_pointer(ii) = tmp
                                ii = ii + 1
                                turbResLocal = turbResLocal + tmp**2
                            end do
                        end do
                    end do
                end do
            end do
        end do

        call VecRestoreArrayF90(rVec, rvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        if (present(flowRes) .and. present(turbRes) .and. present(totalRes)) then
            call mpi_allreduce((/flowResLocal, turbResLocal/), tmp2, 2, adflow_real, &
                               mpi_sum, ADflow_comm_world, ierr)
            flowRes = sqrt(tmp2(1))
            totalRes = sqrt(tmp2(1) + tmp2(2))
            if (tmp2(2) > zero) then
                turbRes = sqrt(tmp2(2))
            else
                turbRes = zero
            end if
        end if

    end subroutine setRVec

    subroutine setW(wVec)

        use constants
        use blockPointers, only: nDom, il, jl, kl, w
        use inputTimeSpectral, only: nTimeIntervalsSpectral
        use flowVarRefState, only: nw, nwf, nt1, nt2, winf
        use utils, only: setPointers, EChk

        implicit none

        Vec wVec
        integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii
        real(kind=realType), pointer :: wvec_pointer(:)
        real(kind=realType) :: cs(1:nw)

        call getNKColScale(cs)

        call VecGetArrayReadF90(wVec, wvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ii = 1
        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                call setPointers(nn, 1_intType, sps)

                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il
                            do l = 1, nwf
                                w(i, j, k, l) = wvec_pointer(ii)
                                ii = ii + 1
                            end do
                            ! Clip the turb to prevent negative turb SA
                            ! values. This is similar to the pressure
                            ! clip. Need to check this for other Turb models.
                            do l = nt1, nt2
                                w(i, j, k, l) = max(1e-6 * winf(l), wvec_pointer(ii) / cs(l))
                                ii = ii + 1
                            end do
                        end do
                    end do
                end do
            end do
        end do

        call VecRestoreArrayReadF90(wVec, wvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

    end subroutine setW

    subroutine getStates(states, ndimw)
        ! Return the state vector, w to Python

        use constants
        use blockPointers, only: il, jl, kl, nDom, w
        use inputTimeSpectral, only: nTimeIntervalsSpectral
        use flowvarrefstate, only: nw
        use utils, only: setPointers

        implicit none

        integer(kind=intType), intent(in) :: ndimw
        real(kind=realType), dimension(ndimw), intent(out) :: states(ndimw)

        ! Local Variables
        integer(kind=intType) :: nn, i, j, k, l, counter, sps

        counter = 0
        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                call setPointers(nn, 1, sps)
                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il
                            do l = 1, nw
                                counter = counter + 1
                                states(counter) = w(i, j, k, l)
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end subroutine getStates

    subroutine getRes(res, ndimw)

        ! Compute the residual and return result to Python
        use constants
        use blockPointers, only: il, jl, kl, nDom, dw, volRef
        use inputTimeSpectral, only: nTimeIntervalsSpectral
        use flowvarrefstate, only: nw
        use utils, only: setPointers

        implicit none

        integer(kind=intType), intent(in) :: ndimw
        real(kind=realType), dimension(ndimw), intent(inout) :: res(ndimw)

        ! Local Variables
        integer(kind=intType) :: nn, i, j, k, l, counter, sps
        real(kind=realType) :: ovv

        call computeResidualNK(useUpdateIntermed=.True.)
        counter = 0
        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                call setPointers(nn, 1, sps)
                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il
                            ovv = one / volRef(i, j, k)
                            do l = 1, nw
                                counter = counter + 1
                                res(counter) = dw(i, j, k, l) * ovv
                            end do
                        end do
                    end do
                end do
            end do
        end do

    end subroutine getRes

    subroutine setStates(states, ndimw)

        ! Take in externallly generated states and set them in ADflow
        use constants
        use blockPointers, only: il, jl, kl, nDom, w
        use inputTimeSpectral, only: nTimeIntervalsSpectral
        use flowvarrefstate, only: nw
        use utils, only: setPointers

        implicit none

        integer(kind=intType), intent(in) :: ndimw
        real(kind=realType), dimension(ndimw), intent(in) :: states(ndimw)

        ! Local Variables
        integer(kind=intType) :: nn, i, j, k, l, counter, sps

        counter = 0
        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                call setPointers(nn, 1, sps)
                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il
                            do l = 1, nw
                                counter = counter + 1
                                w(i, j, k, l) = states(counter)
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end subroutine setStates

    subroutine getInfoSize(iSize)
        use constants
        use blockPointers, only: ib, jb, kb, nDom
        use inputTimeSpectral, only: nTimeIntervalsSpectral
        use flowvarrefstate, only: nw, viscous, eddymodel
        use utils, only: setPointers

        implicit none
        integer(kind=intType), intent(out) :: iSize
        integer(kind=intType) :: nn, sps, nc
        ! Determine the size of a flat array needed to store w, P, ( and
        ! rlv, rev if necessary) with full double halos.
        iSize = 0
        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                call setPointers(nn, 1_intType, sps)
                nc = (kb + 1) * (jb + 1) * (ib + 1)
                iSize = iSize + nc * (nw + 1) ! plus 1 for the P
                if (viscous) then
                    iSize = iSize + nc
                end if
                if (eddyModel) then
                    iSize = iSize + nc
                end if
            end do
        end do
    end subroutine getInfoSize

    subroutine setInfo(info, iSize)

        use constants
        use blockPointers, only: w, p, ib, jb, kb, rlv, rev, nDom
        use inputTimeSpectral, only: nTimeIntervalsSpectral
        use flowvarrefstate, only: nw, viscous, eddymodel
        use utils, only: setPointers
        implicit none

        integer(kind=intType), intent(in) :: iSize
        real(kind=realType), intent(in), dimension(iSize) :: info
        integer(kind=intType) :: nn, counter, i, j, k, l, sps
        ! Determine the size of a flat array needed to store w, P, ( and
        ! rlv, rev if necessary) with full double halos.
        counter = 0
        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                call setPointers(nn, 1, sps)
                do k = 0, kb
                    do j = 0, jb
                        do i = 0, ib
                            do l = 1, nw
                                counter = counter + 1
                                w(i, j, k, l) = info(counter)
                            end do

                            counter = counter + 1
                            P(i, j, k) = info(counter)

                            if (viscous) then
                                counter = counter + 1
                                rlv(i, j, k) = info(counter)
                            end if

                            if (eddyModel) then
                                counter = counter + 1
                                rev(i, j, k) = info(counter)
                            end if
                        end do
                    end do
                end do
            end do
        end do
    end subroutine setInfo

    subroutine getInfo(info, iSize)

        use constants
        use blockPointers, only: w, p, ib, jb, kb, rlv, rev, nDom
        use inputTimeSpectral, only: nTimeIntervalsSpectral
        use flowvarrefstate, only: nw, viscous, eddymodel
        use utils, only: setPointers

        implicit none

        integer(kind=intType), intent(in) :: iSize
        real(kind=realType), intent(out), dimension(iSize) :: info
        integer(kind=intType) :: nn, counter, i, j, k, l, sps
        ! Determine the size of a flat array needed to store w, P, ( and
        ! rlv, rev if necessary) with full double halos.
        counter = 0
        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                call setPointers(nn, 1, sps)
                do k = 0, kb
                    do j = 0, jb
                        do i = 0, ib
                            do l = 1, nw
                                counter = counter + 1
                                info(counter) = w(i, j, k, l)
                            end do

                            counter = counter + 1
                            info(counter) = P(i, j, k)

                            if (viscous) then
                                counter = counter + 1
                                info(counter) = rlv(i, j, k)
                            end if

                            if (eddyModel) then
                                counter = counter + 1
                                info(counter) = rev(i, j, k)
                            end if
                        end do
                    end do
                end do
            end do
        end do
    end subroutine getInfo

    subroutine getEWTol(norm, old_norm, rtol_last, rtol)

        use constants
        implicit none

        ! There are the default EW Parameters from PETSc. They seem to work well
        !version:           2
        !rtol_0:  0.300000000000000
        !rtol_max:  0.900000000000000
        !gamma:   1.00000000000000
        !alpha:   1.61803398874989
        !alpha2:   1.61803398874989
        !threshold:  0.100000000000000

        real(kind=alwaysrealType), intent(in) :: norm, old_norm, rtol_last
        real(kind=alwaysrealType), intent(out) :: rtol
        real(kind=alwaysrealType) :: rtol_max, gamma, alpha, alpha2, threshold, stol

        rtol_max = 0.8_realType
        gamma = 1.0_realType
        alpha = (1.0_realType + sqrt(five)) / 2.0_realType
        alpha2 = (1.0_realType + sqrt(five)) / 2.0_realType
        threshold = 0.10_realType
        ! We use version 2:
        rtol = gamma * (norm / old_norm)**alpha
        stol = gamma * rtol_last**alpha

        if (stol > threshold) then
            rtol = max(rtol, stol)
        end if

        ! Safeguard: avoid rtol greater than one
        rtol = min(rtol, rtol_max)

    end subroutine getEWTol
end module NKSolver

module ANKSolver

    use constants
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    Mat dRdw, dRdwPre, timeStepMat
    Vec wVec, rVec, deltaW, baseRes
    KSP ANK_KSP

    ! Turb KSP related PETSc objects
    Mat dRdwTurb, dRdwPreTurb
    Vec wVecTurb, rVecTurb, deltaWTurb, baseResTurb
    KSP ANK_KSPTurb

    PetscFortranAddr ctx(1)

    ! Options for ANK Solver
    logical :: useANKSolver
    integer(kind=intType) :: ANK_jacobianLag
    integer(kind=intType) :: ANK_subSpace
    integer(kind=intType) :: ANK_maxIter
    integer(kind=intType) :: ANK_asmOverlap
    integer(kind=intType) :: ANK_asmOverlapCoarse
    integer(kind=intType) :: ANK_iluFill
    integer(kind=intType) :: ANK_iluFillCoarse
    integer(kind=intType) :: ANK_innerPreConIts
    integer(kind=intType) :: ANK_innerPreConItsCoarse
    integer(kind=intType) :: ANK_outerPreConIts
    integer(kind=intType) :: ANK_AMGLevels
    integer(kind=intType) :: ANK_AMGNSmooth
    character(len=maxStringLen) :: ANK_precondType
    real(kind=realType) :: ANK_rtol
    real(kind=realType) :: ANK_atol_buffer
    real(kind=realType) :: ANK_linResMax
    real(kind=realType) :: ANK_switchTol
    real(kind=realType) :: ANK_divTol = 10
    logical :: ANK_useTurbDADI
    logical :: ANK_useApproxSA
    real(kind=realType) :: ANK_turbcflscale
    logical :: ANK_useFullVisc
    logical :: ANK_ADPC
    logical :: ANK_turbDebug
    logical :: ANK_useMatrixFree
    character(len=maxStringLen) :: ANK_charTimeStepType
    integer(kind=intType) :: ANK_nsubIterTurb

    ! Misc variables
    real(kind=realType) :: ANK_CFL, ANK_CFL0, ANK_CFLLimit, ANK_CFLFactor, ANK_CFLCutback
    real(kind=realType) :: ANK_CFLMin0, ANK_CFLMin, ANK_CFLMinBase, ANK_CFLExponent
    real(kind=realType) :: ANK_stepMin, ANK_StepFactor, ANK_constCFLStep
    real(kind=realType) :: ANK_secondOrdSwitchTol, ANK_coupledSwitchTol
    real(kind=realType) :: ANK_physLSTol, ANK_unstdyLSTol
    real(kind=realType) :: ANK_pcUpdateTol, ANK_pcUpdateTol2
    real(kind=realType) :: ANK_pcUpdateCutoff
    real(kind=realType) :: lambda
    ! ------------------------------------------------------------------
    ! P&Z stepping mode (ANKPZStepping): when ON, ADflow's own CFL
    ! controller, physicality check and unsteady line search are fully
    ! replaced by the Piotrowski/Zingg (Diablo) system, expressed in CFL
    ! units instead of dt_ref:
    !   phase 1 (approximate-Newton, totalR > ANK_secondOrdSwitchTol*totalR0):
    !     CFL = CFL0 * b^n           (Osusky & Zingg Eq. after (30): dt_ref = a b^n)
    !     + source-term dt restriction (thesis Eq. 3.14) forced ON
    !   phase 2 (inexact-Newton analog, CSANK):
    !     CFL = max(alpha * Rd^-beta, CFL_prev)   (Mulder & van Leer SER)
    !     alpha anchored at the phase switch for continuity
    !   rejection (thesis Algs. 2 & 4): failed update => CFL = max(CFL/2, CFLmin)
    ! ------------------------------------------------------------------
    logical :: ANK_pzStepping = .False.
    real(kind=realType) :: ANK_pzCFL0 = 1.0_realType     ! 'a' of a*b^n, in CFL units
    real(kind=realType) :: ANK_pzGrowth = 1.3_realType   ! 'b' of a*b^n (thesis: 1.3)
    real(kind=realType) :: ANK_pzBeta = 1.75_realType    ! SER exponent beta (thesis: [1.5, 2.0])
    real(kind=realType) :: ANK_pzCFLMin = 0.1_realType   ! dt_ref,min analogue, in CFL units
    real(kind=realType) :: ANK_pzCFLMax = 1.0e10_realType ! numerical ceiling for CFL -> infinity
    ! PZ internal state (not options)
    real(kind=realType) :: pz_CFLRef = -one        ! running reference CFL
    real(kind=realType) :: pz_alphaSER = -one      ! SER alpha; < 0 = not anchored yet
    logical :: pz_prevAccepted = .False.           ! last outer step was accepted
    logical :: ANK_solverSetup = .False.
    logical :: ANK_CFLReset
    integer(kind=intTYpe) :: ANK_iter
    integer(kind=intType) :: nState
    real(kind=alwaysRealType) :: totalR_old, totalR_pcUpdate ! for recording the previous residual
    real(kind=alwaysRealType) :: rtolLast, linResOld ! for recording the previous relativel tolerance for Eisenstat-Walker
    logical :: ANK_useDissApprox

    ! Turb KSP related modifications
    logical :: ANK_coupled = .False.
    logical :: ANK_turbSetup = .False.
    integer(kind=intType) :: ANK_iterTurb, nStateTurb
    real(kind=realType) :: lambdaTurb, ANK_physLSTolTurb
    real(kind=realType) :: ANK_physLSTolReTheta = 0.99_realType
    real(kind=realType) :: omegaMinGamma = 0.05_realType
    real(kind=alwaysRealType) :: linResOldTurb

    ! ------------------------------------------------------------------
    ! Stall diagnostics (VERIF_06), gated by inputIteration%solverStallDiag.
    ! These record WHICH limiter bound the global step on the last
    ! iteration, so a stalling log says why instead of just showing a small
    ! Step. Purely diagnostic -- never read back into the solution path.
    !
    ! stallBindVar codes (the variable whose cell produced the binding ratio
    ! in physicalityCheckANK/Turb):
    !   0 = nothing bound it (lambda stayed at its incoming value)
    !   1 = density        2 = total energy
    !   3 = nuTilde        4 = gamma           5 = reTheta
    integer(kind=intType) :: stallBindVar = 0
    integer(kind=intType) :: stallBindLoc(4) = 0     ! (nn, i, j, k) of that cell
    real(kind=alwaysRealType) :: stallLamPhys = one  ! lambda after physicality
    real(kind=alwaysRealType) :: stallLamLS = one    ! lambda after the line search
    integer(kind=intType) :: stallNBacktrack = 0     ! unsteady-LS backtracks used
    logical :: stallLSFailed = .false.               ! LS exhausted its budget
    logical :: stallCFLFloored = .false.             ! cutback requested but ANK_CFLMin blocked it
    real(kind=alwaysRealType) :: stallCFLBefore = zero, stallCFLAfter = zero
    ! Algorithm 2 activity on the last coupled step: how many cells needed ANY
    ! back-off, and the worst factor applied. A front that is being crushed
    ! every iteration shows up here and nowhere else -- the damping happens
    ! AFTER the line search, so neither Step nor Lin Res reveals it.
    integer(kind=intType) :: stallSoftDampG = 0, stallSoftDampR = 0
    real(kind=alwaysRealType) :: stallMinDampG = one, stallMinDampR = one
    ! Oscillation detector: previous total residual and a run length of
    ! consecutive iterations in which it ROSE. A solver that is oscillating
    ! rather than stalling shows a nonzero, repeatedly-resetting run length.
    ! Unit column scale (VERIF_06 F10): per-variable factors making the
    ! scaled state RMS one. <=0 entries mean 'not measured yet', in which
    ! case getFullColScale falls back to the legacy scaling.
    real(kind=realType) :: colScaleUnitFac(8) = -one
    real(kind=alwaysRealType) :: stallResPrev = -one
    integer(kind=intType) :: stallRiseCount = 0, stallRiseMax = 0

contains

    function stallVarName(code) result(nameOut)
        ! Human-readable name for a stallBindVar code (see above).
        use constants
        implicit none
        integer(kind=intType), intent(in) :: code
        character(len=8) :: nameOut

        select case (code)
        case (1)
            nameOut = 'rho'
        case (2)
            nameOut = 'energy'
        case (3)
            nameOut = 'nuTilde'
        case (4)
            nameOut = 'gamma'
        case (5)
            nameOut = 'reTheta'
        case default
            nameOut = 'none'
        end select

    end function stallVarName

    subroutine setupANKsolver

        ! Setup the PETSc objects for the approximate Newton-Krylov solver.
        ! destroyANKsolver can be used to destroy the objects created in this function.

        use constants
        use stencils, only: euler_PC_stencil, N_euler_PC
        use communication, only: adflow_comm_world, myid
        use inputTimeSpectral, only: nTimeIntervalsSpectral
        use inputIteration, only: useLinResMonitor
        use inputPhysics, only: equations
        use flowVarRefState, only: nw, viscous, nwf, nt1, nt2
        use ADjointVars, only: nCellsLocal
        use NKSolver, only: destroyNKSolver, linearResidualMonitor, configureMFFD
        use utils, only: EChk
        use adjointUtils, only: myMatCreate, statePreAllocation
        use amg, only: setupAMG
        implicit none

        ! Working Variables
        integer(kind=intType) :: ierr, nDimw, nDimWTurb
        integer(kind=intType), dimension(:), allocatable :: nnzDiagonal, nnzOffDiag
        integer(kind=intType) :: n_stencil
        integer(kind=intType), dimension(:, :), pointer :: stencil
        integer(kind=intType) :: level

        ! Make sure we don't have memory for the approximate and exact
        ! Newton solvers kicking around at the same time.
        call destroyNKSolver()

        if (.not. ANK_solverSetup) then

            ! Determine if we are in coupled mode
            if (ANK_coupled) then
                nState = nw
            else
                nState = nwf
            end if

            nDimW = nState * nCellsLocal(1_intTYpe) * nTimeIntervalsSpectral

            call VecCreate(ADFLOW_COMM_WORLD, wVec, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call VecSetSizes(wVec, nDimW, PETSC_DECIDE, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call VecSetBlockSize(wVec, nState, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call VecSetType(wVec, VECMPI, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            !  Create duplicates for residual and delta
            call VecDuplicate(wVec, rVec, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call VecDuplicate(wVec, deltaW, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call VecDuplicate(wVec, baseRes, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Create Pre-Conditioning Matrix
            allocate (nnzDiagonal(nCellsLocal(1_intType) * nTimeIntervalsSpectral), &
                      nnzOffDiag(nCellsLocal(1_intType) * nTimeIntervalsSpectral))

            stencil => euler_pc_stencil
            n_stencil = N_euler_pc

            level = 1
            call statePreAllocation(nnzDiagonal, nnzOffDiag, nDimW / nState, stencil, n_stencil, &
                                    level, .False.)
            call myMatCreate(dRdwPre, nState, nDimW, nDimW, nnzDiagonal, nnzOffDiag, &
                             __FILE__, __LINE__)
            call matSetOption(dRdwPre, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call myMatCreate(timeStepMat, nState, nDimW, nDimW, nnzDiagonal, nnzOffDiag, &
                             __FILE__, __LINE__)
            call matSetOption(timeStepMat, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE, ierr)
            call EChk(ierr, __FILE__, __LINE__)
            deallocate (nnzDiagonal, nnzOffDiag)

            ! Set the mat_row_oriented option to false so that dense
            ! subblocks can be passed in in fortran column-oriented format
            call MatSetOption(dRdWPre, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Setup Matrix-Free dRdw matrix and its function
            call MatCreateMFFD(ADFLOW_COMM_WORLD, nDimW, nDimW, &
                               PETSC_DETERMINE, PETSC_DETERMINE, dRdw, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call MatMFFDSetFunction(dRdw, FormFunction_mf, ctx, ierr)
            call EChk(ierr, __FILE__, __LINE__)
            call configureMFFD(dRdw)

            call MatSetOption(dRdW, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            if (ANK_precondType == 'mg') then
                call setupAMG(drdwpre, nDimW / nState, nState, ANK_AMGLevels, ANK_AMGNSmooth)
            end if

            !  Create the linear solver context
            call KSPCreate(ADFLOW_COMM_WORLD, ANK_KSP, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Set operators for the solver
            if (ANK_useMatrixFree) then
                ! Matrix free drdw
                call KSPSetOperators(ANK_KSP, dRdw, dRdwPre, ierr)
            else
                ! Matrix based drdw = drdwpre
                call KSPSetOperators(ANK_KSP, dRdwPre, dRdwPre, ierr)
            end if
            call EChk(ierr, __FILE__, __LINE__)

            if (useLinResMonitor) then

#if PETSC_VERSION_GE(3,8,0)
                ! This is probably wrong. NO petsc_null_context
                call KSPMonitorSet(ANK_KSP, LinearResidualMonitor, PETSC_NULL_FUNCTION, &
                                   PETSC_NULL_FUNCTION, ierr)
#else
                call KSPMonitorSet(ANK_KSP, LinearResidualMonitor, PETSC_NULL_OBJECT, &
                                   PETSC_NULL_FUNCTION, ierr)

#endif
                call EChk(ierr, __FILE__, __LINE__)
            end if

            ANK_solverSetup = .True.
            ANK_iter = 0
            ANK_useDissApprox = .False.

            ! Check if we need to set up the Turb KSP
            if ((.not. ANK_coupled) .and. (.not. ANK_useTurbDADI) .and. equations == RANSEquations) then
                nStateTurb = nt2 - nt1 + 1

                nDimWTurb = nStateTurb * nCellsLocal(1_intTYpe) * nTimeIntervalsSpectral

                call VecCreate(ADFLOW_COMM_WORLD, wVecTurb, ierr)
                call EChk(ierr, __FILE__, __LINE__)

                call VecSetSizes(wVecTurb, nDimWTurb, PETSC_DECIDE, ierr)
                call EChk(ierr, __FILE__, __LINE__)

                call VecSetBlockSize(wVecTurb, nStateTurb, ierr)
                call EChk(ierr, __FILE__, __LINE__)

                call VecSetType(wVecTurb, VECMPI, ierr)
                call EChk(ierr, __FILE__, __LINE__)

                !  Create duplicates for residual and delta
                call VecDuplicate(wVecTurb, rVecTurb, ierr)
                call EChk(ierr, __FILE__, __LINE__)

                call VecDuplicate(wVecTurb, deltaWTurb, ierr)
                call EChk(ierr, __FILE__, __LINE__)

                call VecDuplicate(wVecTurb, baseResTurb, ierr)
                call EChk(ierr, __FILE__, __LINE__)

                ! Create Pre-Conditioning Matrix
                allocate (nnzDiagonal(nCellsLocal(1_intType) * nTimeIntervalsSpectral), &
                          nnzOffDiag(nCellsLocal(1_intType) * nTimeIntervalsSpectral))

                stencil => euler_pc_stencil
                n_stencil = N_euler_pc

                level = 1
                call statePreAllocation(nnzDiagonal, nnzOffDiag, nDimWTurb / nStateTurb, stencil, n_stencil, &
                                        level, .False.)
                call myMatCreate(dRdwPreTurb, nStateTurb, nDimWTurb, nDimWTurb, nnzDiagonal, nnzOffDiag, &
                                 __FILE__, __LINE__)

                call matSetOption(dRdwPreTurb, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE, ierr)
                call EChk(ierr, __FILE__, __LINE__)
                deallocate (nnzDiagonal, nnzOffDiag)

                ! Set the mat_row_oriented option to false so that dense
                ! subblocks can be passed in in fortran column-oriented format
                call MatSetOption(dRdWPreTurb, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
                call EChk(ierr, __FILE__, __LINE__)

                ! Setup Matrix-Free dRdw matrix and its function
                call MatCreateMFFD(ADFLOW_COMM_WORLD, nDimWTurb, nDimWTurb, &
                                   PETSC_DETERMINE, PETSC_DETERMINE, dRdwTurb, ierr)
                call EChk(ierr, __FILE__, __LINE__)

                call MatMFFDSetFunction(dRdwTurb, FormFunction_mf_Turb, ctx, ierr)
                call EChk(ierr, __FILE__, __LINE__)
                call configureMFFD(dRdwTurb)

                call MatSetOption(dRdWTurb, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
                call EChk(ierr, __FILE__, __LINE__)

                !  Create the linear solver context
                call KSPCreate(ADFLOW_COMM_WORLD, ANK_KSPTurb, ierr)
                call EChk(ierr, __FILE__, __LINE__)

                ! Set operators for the solver
                if (ANK_useMatrixFree) then
                    ! Matrix free
                    call KSPSetOperators(ANK_KSPTurb, dRdwTurb, dRdwPreTurb, ierr)
                else
                    ! Matrix based
                    call KSPSetOperators(ANK_KSPTurb, dRdwPreTurb, dRdwPreTurb, ierr)
                end if
                call EChk(ierr, __FILE__, __LINE__)

                ANK_turbSetup = .True.
                ANK_iterTurb = 0
            end if
        end if

    end subroutine setupANKsolver

    subroutine FormJacobianANK

        use constants
        use flowVarRefState, only: nw, nwf, nt1, nt2
        use blockPointers, only: nDom, volRef, il, jl, kl, w, dw, dtl, globalCell, iblank
        use inputTimeSpectral, only: nTimeIntervalsSpectral
        use inputIteration, only: turbResScale, transitionNK
        use inputADjoint, only: viscPC
        use inputDiscretization, only: approxSA
        use inputPhysics, only: turbModel
        use iteration, only: totalR0, totalR
        use utils, only: EChk, setPointers
        use adjointUtils, only: setupStateResidualMatrix, setupStandardKSP, setupStandardMultigrid
        use communication
        use amg, only: setupShellPC, destroyShellPC, applyShellPC, coarseIndices, A
        implicit none

        ! Local Variables
        character(len=maxStringLen) :: preConSide, localPCType, kspObjectType, globalPCType, localOrdering
        integer(kind=intType) :: ierr
        logical :: useAD, usePC, useTranspose, useObjective, tmp, frozenTurb
        real(kind=realType) :: dtinv, rho
        integer(kind=intType) :: i, j, k, l, ii, irow, nn, sps, outerPreConIts, subspace, lvl
        integer(kind=intType), dimension(2:10) :: coarseRows
        real(kind=realType), dimension(:, :), allocatable :: blk
        logical :: useCoarseMats
        PC shellPC

        if (ANK_precondType == 'mg') then
            useCoarseMats = .True.
        else
            useCoarseMats = .False.
        end if

        ! Assemble the approximate PC (fine level, level 1)
        useAD = ANK_ADPC
        frozenTurb = (.not. ANK_coupled)
        usePC = .True.
        useTranspose = .False.
        useObjective = .False.
        tmp = viscPC ! Save what is in viscPC and set to False
        viscPC = .False.

        if (totalR > ANK_secondOrdSwitchTol * totalR0) &
            approxSA = .True.

        ! Create the preconditoner matrix
        call setupStateResidualMatrix(dRdwPre, useAD, usePC, useTranspose, &
                                      useObjective, frozenTurb, 1_intType, useCoarseMats=useCoarseMats)

        ! Reset saved value
        viscPC = tmp
        approxSA = .False.

        ! Begin PETSc matrix assembly
        call MatAssemblyBegin(dRdwPre, MAT_FINAL_ASSEMBLY, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Complete the matrix assembly
        call MatAssemblyEnd(dRdwPre, MAT_FINAL_ASSEMBLY, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Column-scale the coupled PC to match the column-scaled state
        ! vector (timeStepMat is scaled at its own assembly, so it is added
        ! afterwards already consistent). Skipped for non-transition models.
        if (turbModel == spalartallmarasnoft2gammaretheta .and. transitionNK .and. ANK_coupled) then
            call applyANKColumnScaling(dRdwPre)
            if (ANK_precondType == 'mg' .and. myid == 0) then
                print *, 'Warning: ANK multigrid coarse levels are not ', &
                    'column-scaled for SA-Gamma-Retheta; use ', &
                    'ANKGlobalPreconditioner=additive Schwarz.'
            end if
        end if

        ! Add the contribution from the time step matrix
        call MatAXPY(dRdwPre, one, timeStepMat, SUBSET_NONZERO_PATTERN, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        if (useCoarseMats) then
            do lvl = 2, ANK_AMGLevels
                call MatAssemblyBegin(A(lvl), MAT_FINAL_ASSEMBLY, ierr)
                call EChk(ierr, __FILE__, __LINE__)
                call MatAssemblyEnd(A(lvl), MAT_FINAL_ASSEMBLY, ierr)
                call EChk(ierr, __FILE__, __LINE__)
            end do
        end if

        ! Set up KSP options
        preConSide = 'right'
        localPCType = 'ilu'
        kspObjectType = 'gmres'
        globalPCType = 'asm'
        localOrdering = 'rcm'
        outerPreConIts = ank_outerPreconIts

        if (ANK_subspace < 0) then
            subspace = ANK_maxIter
        else
            subspace = ANK_subspace
        end if

        ! Set up the KSP using the same code as used for the adjoint
        if (ANK_precondType == 'asm') then
            call setupStandardKSP(ANK_KSP, kspObjectType, subSpace, &
                                  preConSide, globalPCType, ANK_asmOverlap, outerPreConIts, localPCType, &
                                  localOrdering, ANK_iluFill, ANK_innerPreConIts)
        else if (ANK_precondType == 'mg') then
            call setupStandardMultigrid(ANK_KSP, kspObjectType, subSpace, &
                                        preConSide, ANK_asmOverlap, outerPreConIts, &
                                        localOrdering, ANK_iluFill, ANK_innerPreConIts, &
                                        ANK_asmOverlapCoarse, ANK_iluFillCoarse, ANK_innerPreConItsCoarse)
        end if

        ! Don't do iterative refinement
        call KSPGMRESSetCGSRefinementType(ANK_KSP, KSP_GMRES_CGS_REFINE_NEVER, ierr)
        call EChk(ierr, __FILE__, __LINE__)

    end subroutine FormJacobianANK

    subroutine computeTimeStepMat(usePC)
        ! Loops through all cells and computes the time step terms
        ! The terms are stored in the PETSc matrix timeStepMat

        use constants
        use blockPointers, only: nDom, il, jl, kl, globalCell
        use inputTimeSpectral, only: nTimeIntervalsSpectral
        use inputIteration, only: turbResScale, transitionNK
        use inputPhysics, only: turbModel
        use utils, only: EChk, setPointers
        use amg, only: coarseIndices, A
        implicit none

        ! Input variables
        logical, intent(in) :: usePC

        ! Local Variables
        character(len=maxStringLen) :: preConSide, localPCType, kspObjectType, globalPCType, localOrdering
        integer(kind=intType) :: ierr
        integer(kind=intType) :: i, j, k, irow, nn, sps, lvl
        integer(kind=intType), dimension(2:10) :: coarseRows
        real(kind=realType), dimension(nState, nState) :: timeStepBlock
        logical :: useCoarseMats

        if (ANK_precondType == 'mg') then
            useCoarseMats = .True.
        else
            useCoarseMats = .False.
        end if

        ! Zero out the time step matrix before we start
        call MatZeroEntries(timeStepMat, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                call setPointers(nn, 1_intType, sps)
                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il

                            ! Compute the block for this cell
                            call computeTimeStepBlock(i, j, k, timeStepBlock)

                            ! Get the global cell index
                            irow = globalCell(i, j, k)

                            ! Add the contribution to the PETSc matrix
                            call MatSetValuesBlocked(timeStepMat, 1, irow, 1, irow, timeStepBlock, ADD_VALUES, ierr)
                            call EChk(ierr, __FILE__, __LINE__)

                            ! Extension for setting coarse grids
                            ! We only do this when we are updating the ANK PC
                            if (useCoarseMats .and. usePC) then
                                do lvl = 2, ANK_AMGLevels
                                    coarseRows(lvl) = coarseIndices(nn, lvl - 1)%arr(i, j, k)
                                    call MatSetValuesBlocked(A(lvl), 1, coarseRows(lvl), 1, coarseRows(lvl), &
                                                             timeStepBlock, ADD_VALUES, ierr)
                                    call EChk(ierr, __FILE__, __LINE__)
                                end do
                            end if

                        end do
                    end do
                end do
            end do
        end do

        ! PETSc Matrix Assembly
        call MatAssemblyBegin(timeStepMat, MAT_FINAL_ASSEMBLY, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        call MatAssemblyEnd(timeStepMat, MAT_FINAL_ASSEMBLY, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Column-scale to match the column-scaled state vector (both the
        ! matrix-free MatMultAdd in FormFunction_mf and the PC MatAXPY
        ! consume this matrix). No-op arithmetic for non-transition models
        ! (skipped entirely there).
        if (turbModel == spalartallmarasnoft2gammaretheta .and. transitionNK .and. ANK_coupled) then
            call applyANKColumnScaling(timeStepMat)
        end if

    end subroutine computeTimeStepMat

    subroutine computeTimeStepBlock(i, j, k, timeStepBlock)
        ! Computes the time step block matrix for a given cell i, j, k.

        use constants
        use inputPhysics, only: machInf => mach, turbModel
        use blockPointers, only: volRef, w, dtl, gamma, p, aa, srcLambda
        use flowVarRefState, only: viscous, nt1, nt2
        use inputIteration, only: turbResScale, transitionSrcDtRestrict, &
                                  transitionSrcDtLimit, noBacktrackCount, srcDtDeactivateIters, transitionNK
        use communication
        implicit none

        ! Input variables
        integer(kind=intType), intent(in) :: i, j, k

        ! Output variables
        real(kind=realType), dimension(nState, nState), intent(out) :: timeStepBlock

        ! Local variables
        integer(kind=intType) :: l, l1
        real(kind=realType) :: dtInvSrc
        real(kind=realType) :: blendFactor, dtInv, rho, velX, velY, velZ
        real(kind=realType) :: speed, speedOfSound, mach, machSqr, machSqrTrunc, alpha, beta, tau, gammaMinusOne
        real(kind=realType) :: speedXY, sinTheta, cosTheta, sinAlpha, cosAlpha
        real(kind=realType), dimension(nState, nState) :: streamToCart, symmToCons, consToSymm, stateToCons

        ! Zero the block matrices
        timeStepBlock = zero
        stateToCons = zero
        streamToCart = zero
        symmToCons = zero
        consToSymm = zero

        ! Save density and velocity components for convenience
        rho = w(i, j, k, iRho)
        velX = w(i, j, k, ivx)
        velY = w(i, j, k, ivy)
        velZ = w(i, j, k, ivz)

        ! Calculate one over time step for this cell.
        ! Multiply dtl by cell volume to get the time step required for a CFL of one,
        ! then multiply with the actual CFL number in the solver.
        dtInv = one / (ANK_CFL * dtl(i, j, k) * volRef(i, j, k))

        ! We need to convert the velocity updates to momentum updates to get the desired effect from time steps.
        ! To do this, we form a "pseudo" Jacobian dU/du for this cell, where U is the vector of conservative
        ! variables, and u is the vector of state variables. Only the velocity entries are modified because ADflow
        ! saves density, velocities, and total energy in the state vector w(:,:,:,:).

        stateToCons(iRho, iRho) = one
        stateToCons(ivx, iRho) = velX
        stateToCons(ivx, ivx) = rho
        stateToCons(ivy, iRho) = velY
        stateToCons(ivy, ivy) = rho
        stateToCons(ivz, iRho) = velZ
        stateToCons(ivz, ivz) = rho
        stateToCons(iRhoE, iRhoE) = one

        if (ANK_coupled) then
            ! The turbulence variable can get a different CFL number, so we scale it by ANK_turbCFLScale.
            ! In addition, turbResScale is required because the turbulent residuals are scaled with it.
            do l = nt1, nt2
                stateToCons(l, l) = turbResScale(l - nt1 + 1) / ANK_turbCFLScale
            end do
        end if

        if (ANK_charTimeStepType == 'None') then

            timeStepBlock = stateToCons * dtInv

        else

            ! The equations referenced in this block of the code are from the paper:
            ! "Improving the Performance of a Compressible RANS Solver for Low and High Mach Number Flows" (Seraj2022c)

            ! Characteristic time-stepping is not applied to the turbulence equation
            if (ANK_coupled) then
                do l = nt1, nt2
                    timeStepBlock(l, l) = one
                    streamToCart(l, l) = one
                    symmToCons(l, l) = one
                    consToSymm(l, l) = one
                end do
            end if

            ! Compute the speed of sound squared for inviscid flow
            if (.not. viscous) then
                aa(i, j, k) = gamma(i, j, k) * p(i, j, k) / rho
            end if

            ! Compute the Mach number in each cell and some other repeated terms
            speed = SQRT(velX**2 + velY**2 + velZ**2)
            speedOfSound = SQRT(aa(i, j, k))
            mach = speed / speedOfSound
            machSqr = mach**2
            gammaMinusOne = gamma(i, j, k) - one

            ! Define the state transformation matrix from Euler symmetrizing variables to conservative variables (Eq. 22)
            symmToCons(iRho, 1) = rho / speedOfSound
            symmToCons(iRho, 5) = -one / aa(i, j, k)
            symmToCons(ivx, 1) = rho * velX / speedOfSound
            symmToCons(ivx, 2) = rho
            symmToCons(ivx, 5) = -velX / aa(i, j, k)
            symmToCons(ivy, 1) = rho * velY / speedOfSound
            symmToCons(ivy, 3) = rho
            symmToCons(ivy, 5) = -velY / aa(i, j, k)
            symmToCons(ivz, 1) = rho * velZ / speedOfSound
            symmToCons(ivz, 4) = rho
            symmToCons(ivz, 5) = -velZ / aa(i, j, k)
            symmToCons(iRhoE, 1) = rho * speedOfSound * (machSqr / 2 + 1 / gammaMinusOne)
            symmToCons(iRhoE, 2) = rho * velX
            symmToCons(iRhoE, 3) = rho * velY
            symmToCons(iRhoE, 4) = rho * velZ
            symmToCons(iRhoE, 5) = -machSqr / 2

            ! Define the inverse of the state transformation matrix above
            consToSymm(1, iRho) = gammaMinusOne / 2 * speedOfSound * machSqr / rho
            consToSymm(1, ivx) = -gammaMinusOne * velX / (rho * speedOfSound)
            consToSymm(1, ivy) = -gammaMinusOne * velY / (rho * speedOfSound)
            consToSymm(1, ivz) = -gammaMinusOne * velZ / (rho * speedOfSound)
            consToSymm(1, iRhoE) = gammaMinusOne / (rho * speedOfSound)
            consToSymm(2, iRho) = -velX / rho
            consToSymm(2, ivx) = one / rho
            consToSymm(3, iRho) = -velY / rho
            consToSymm(3, ivy) = one / rho
            consToSymm(4, iRho) = -velZ / rho
            consToSymm(4, ivz) = one / rho
            consToSymm(5, iRho) = aa(i, j, k) * (gammaMinusOne / 2 * machSqr - one)
            consToSymm(5, ivx) = -gammaMinusOne * velX
            consToSymm(5, ivy) = -gammaMinusOne * velY
            consToSymm(5, ivz) = -gammaMinusOne * velZ
            consToSymm(5, iRhoE) = gammaMinusOne

            ! Compute the CFL-based blending factor
            blendFactor = ANK_CFL / ANK_CFLLimit

            if (ANK_charTimeStepType == 'VLR') then

                ! Truncate the squared Mach number (Eq. 25)
                machSqrTrunc = MAX(machSqr, 1e-4 * machInf**2)

                ! Define the beta and tau terms in the VLR matrix (Eq. 24)
                if (mach < one) then
                    beta = SQRT(one - machSqrTrunc)
                    tau = beta
                else
                    beta = SQRT(machSqrTrunc - one)
                    tau = SQRT(one - one / machSqrTrunc) + 1e-4
                end if

                ! Define the VLR matrix in Euler symmetrizing variables and in the streamwise coordinate frame
                ! This is Eq. 23 combined with the blending in Eq. 16
                timeStepBlock(1, 1) = blendFactor * (beta**2 + tau) / (machSqrTrunc * tau) + (one - blendFactor) * one
                timeStepBlock(1, 2) = blendFactor * one / mach
                timeStepBlock(2, 1) = blendFactor * one / mach
                timeStepBlock(2, 2) = one
                timeStepBlock(3, 3) = blendFactor * one / tau + (one - blendFactor) * one
                timeStepBlock(4, 4) = blendFactor * one / tau + (one - blendFactor) * one
                timeStepBlock(5, 5) = one

                ! Define repeated terms in the rotation matrix
                speedXY = SQRT(velX**2 + velY**2)
                sinTheta = velY / speedXY
                cosTheta = velX / speedXY
                sinAlpha = velZ / speed
                cosAlpha = speedXY / speed

                ! Define the rotation matrix from streamwise coordinates to Cartesian coordinates (Eq. 27)
                streamToCart(1, 1) = one
                streamToCart(2, 2) = cosAlpha * cosTheta
                streamToCart(2, 3) = -sinTheta
                streamToCart(2, 4) = -sinAlpha * cosTheta
                streamToCart(3, 2) = cosAlpha * sinTheta
                streamToCart(3, 3) = cosTheta
                streamToCart(3, 4) = -sinAlpha * sinTheta
                streamToCart(4, 2) = sinAlpha
                streamToCart(4, 4) = cosAlpha
                streamToCart(5, 5) = one

                ! Transform the VLR matrix to conservative variables and Cartesian coordinates (Eq. 26)
                timeStepBlock = MATMUL(streamToCart, timeStepBlock)
                timeStepBlock = MATMUL(timeStepBlock, TRANSPOSE(streamToCart))
                timeStepBlock = MATMUL(symmToCons, timeStepBlock)
                timeStepBlock = MATMUL(timeStepBlock, consToSymm)

                ! Construct the preconditioned time step matrix
                timeStepBlock = MATMUL(timeStepBlock, stateToCons)
                timeStepBlock = timeStepBlock * dtInv

            else if (ANK_charTimeStepType == 'Turkel') then

                ! Truncate the squared Mach number (Eq. 19)
                machSqrTrunc = MIN(one, MAX(machSqr, 1e-4 * machInf**2))

                ! Set alpha to turn off preconditioning for locally supersonic flow (Eq. 20)
                alpha = one - machSqrTrunc**10

                ! Define the Turkel matrix in Euler symmetrizing variables
                ! This is Eq. 18 combined with the blending in Eq. 16
                timeStepBlock(1, 1) = blendFactor * one / machSqrTrunc + (one - blendFactor) * one
                timeStepBlock(2, 1) = blendFactor * alpha * velX / speedOfSound / machSqrTrunc
                timeStepBlock(3, 1) = blendFactor * alpha * velY / speedOfSound / machSqrTrunc
                timeStepBlock(4, 1) = blendFactor * alpha * velZ / speedOfSound / machSqrTrunc
                timeStepBlock(2, 2) = one
                timeStepBlock(3, 3) = one
                timeStepBlock(4, 4) = one
                timeStepBlock(5, 5) = one

                ! Transform the Turkel matrix to conservative variables (Eq. 21)
                timeStepBlock = MATMUL(symmToCons, timeStepBlock)
                timeStepBlock = MATMUL(timeStepBlock, consToSymm)

                ! Construct the preconditioned time step matrix
                timeStepBlock = MATMUL(timeStepBlock, stateToCons)
                timeStepBlock = timeStepBlock * dtInv

            end if

        end if

        ! Source-term time-step restriction for the coupled SA-Gamma-Retheta
        ! solve (P&Z Eq. 59): keep lambda_source * dt <= transitionSrcDtLimit
        ! on the transition rows. MAX form, same as the turbKSP path in
        ! FormJacobianANKTurb. srcLambda is frozen at the base state in
        ! ANKStep before this matrix is formed. The turb rows are purely
        ! diagonal here (characteristic time stepping is not applied to
        ! them), so overriding the diagonal after the transforms is exact.
        if (ANK_coupled .and. turbModel == spalartallmarasnoft2gammaretheta .and. transitionNK) then
            if (transitionSrcDtRestrict .and. (noBacktrackCount < srcDtDeactivateIters)) then
                do l = nt1, nt2
                    l1 = l - nt1 + 1
                    dtInvSrc = srcLambda(i, j, k, l1) / transitionSrcDtLimit
                    timeStepBlock(l, l) = max(dtInv, dtInvSrc) * turbResScale(l1) / ANK_turbCFLScale
                end do
            end if
        end if

    end subroutine computeTimeStepBlock

    subroutine FormJacobianANKTurb

        use constants
        use flowVarRefState, only: nw, nwf, nt1, nt2, nwt
        use blockPointers, only: nDom, volRef, il, jl, kl, w, dw, dtl, globalCell, srcLambda
        use inputTimeSpectral, only: nTimeIntervalsSpectral
        use inputIteration, only: turbResScale, transitionSrcDtRestrict, transitionSrcDtLimit, &
                                   noBacktrackCount, srcDtDeactivateIters, transitionNK
        use inputADjoint, only: viscPC
        use inputDiscretization, only: approxSA
        use inputPhysics, only: turbModel
        use iteration, only: totalR0, totalR
        use utils, only: EChk, setPointers
        use adjointUtils, only: setupStateResidualMatrix, setupStandardKSP
        use communication
        implicit none

        ! Local Variables
        character(len=maxStringLen) :: preConSide, localPCType, kspObjectType, globalPCType, localOrdering
        integer(kind=intType) :: ierr
        logical :: useAD, usePC, useTranspose, useObjective, tmp, frozenTurb, srcDtRestrictActive
        real(kind=realType) :: dtinv, rho, dtinv_src
        integer(kind=intType) :: i, j, k, l, l1, ii, irow, nn, sps, outerPreConIts, subspace
        real(kind=realType), dimension(:, :), allocatable :: blk

        ! Derived condition for source dt restriction (ANK turbKSP path)
        srcDtRestrictActive = transitionNK .and. transitionSrcDtRestrict .and. (noBacktrackCount < srcDtDeactivateIters)

        ! Assemble the approximate PC (fine level, level 1)
        useAD = ANK_ADPC
        frozenTurb = .False.
        usePC = .True.
        useTranspose = .False.
        useObjective = .False.
        tmp = viscPC ! Save what is in viscPC and set to False
        viscPC = .False.

        if (totalR > ANK_secondOrdSwitchTol * totalR0) &
            approxSA = .True.

        ! Create the preconditoner matrix
        call setupStateResidualMatrix(dRdwPreTurb, useAD, usePC, useTranspose, &
                                      useObjective, frozenTurb, 1_intType, .True.)

        ! Reset saved value
        viscPC = tmp
        approxSA = .False.

        ! Add the contribution from the time step term

        ! Generic block to use while setting values
        allocate (blk(nStateTurb, nStateTurb))

        ! Zero the block once, since the previous entries will be overwritten
        ! for each cell, and zero entries will remain zero.
        blk = zero

        ! For the coupled solver, CFL number for the turbulent variable needs scaling
        ! because the residuals are scaled, and additional scaling of the time step
        ! for the turbulence variable might be required.
        ii = 1
        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                call setPointers(nn, 1_intType, sps)
                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il

                            ! Base dtinv from CFL
                            dtinv = one / (ANK_CFL * dtl(i, j, k) * volRef(i, j, k))

                            do l = nt1, nt2
                                ! l1 is just l that starts with 1 on the turb variables
                                l1 = l - nt1 + 1

                                ! Source dt restriction (P&Z Eq. 59): MAX form.
                                ! turbKSP has a dtinv_CFL from the CFL ramp; Eq. 59 restricts it.
                                ! Contrast with DD-ADI (saGammaRetheta.F90) which uses ADDITIVE
                                ! form because there is no separate dtinv_CFL in the smoother.
                                if (transitionSrcDtRestrict .and. srcDtRestrictActive .and. nwt == 3) then
                                    dtinv_src = srcLambda(i, j, k, l1) / transitionSrcDtLimit
                                    blk(l1, l1) = max(dtinv, dtinv_src) * turbResScale(l1) / ANK_turbCFLScale
                                else
                                    blk(l1, l1) = dtinv * turbResScale(l1) / ANK_turbCFLScale
                                end if
                            end do

                            ! get the global cell index
                            irow = globalCell(i, j, k)

                            ! Add the contribution to the matrix in PETSc
                            call setBlock()
                        end do
                    end do
                end do
            end do
        end do

        ! De-allocate the generic block
        deallocate (blk)

        ! Being PETSc matrix assembly
        call MatAssemblyBegin(dRdwPreTurb, MAT_FINAL_ASSEMBLY, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Complete the matrix assembly
        call MatAssemblyEnd(dRdwPreTurb, MAT_FINAL_ASSEMBLY, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Column-scale the PC to match the column-scaled state vector used
        ! by the matrix-free operator (see getTurbColScale). The assembled
        ! matrix rows already carry turbResScale; right-multiplying by
        ! diag(1/cs) makes it a consistent preconditioner for
        ! S_row * dRdw * diag(1/cs). Skipped entirely for models other
        ! than SA-Gamma-Retheta (cs = 1 there).
        if (turbModel == spalartallmarasnoft2gammaretheta .and. transitionNK) then
            call applyTurbPCColumnScaling()
        end if

        ! Set up KSP options
        preConSide = 'right'
        localPCType = 'ilu'
        kspObjectType = 'gmres'
        globalPCType = 'asm'
        localOrdering = 'rcm'
        outerPreConIts = 1

        ! Set up the KSP using the same code as used for the adjoint
        if (ank_subspace < 0) then
            subspace = ANK_maxIter
        else
            subspace = ANK_subspace
        end if

        call setupStandardKSP(ANK_KSPTurb, kspObjectType, subSpace, &
                              preConSide, globalPCType, ANK_asmOverlap, outerPreConIts, localPCType, &
                              localOrdering, ANK_iluFill, ANK_innerPreConIts)

        ! Don't do iterative refinement
        call KSPGMRESSetCGSRefinementType(ANK_KSPTurb, KSP_GMRES_CGS_REFINE_NEVER, ierr)
        call EChk(ierr, __FILE__, __LINE__)

    contains
        subroutine setBlock()
            ! This subroutine is used to set the diagonal time stepping terms
            ! for the Jacobians in ANK. It is only used to set diagonal blocks

            implicit none

            call MatSetValuesBlocked(dRdwPreTurb, 1, irow, 1, irow, blk, &
                                     ADD_VALUES, ierr)
            call EChk(ierr, __FILE__, __LINE__)

        end subroutine setBlock

        subroutine applyTurbPCColumnScaling()
            ! Right-multiply the assembled turb PC by diag(1/cs), where cs
            ! is the column scale of the state vector (getTurbColScale).
            implicit none

            Vec colVec
            integer(kind=intType) :: iCell, lState, jj
            real(kind=realType), pointer :: cvec_pointer(:)
            real(kind=realType) :: cs(nStateTurb)

            call getTurbColScale(cs, nStateTurb)

            call VecDuplicate(wVecTurb, colVec, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call VecGetArrayF90(colVec, cvec_pointer, ierr)
            call EChk(ierr, __FILE__, __LINE__)
            do jj = 0, size(cvec_pointer) / nStateTurb - 1
                do lState = 1, nStateTurb
                    cvec_pointer(jj * nStateTurb + lState) = one / cs(lState)
                end do
            end do
            call VecRestoreArrayF90(colVec, cvec_pointer, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call MatDiagonalScale(dRdwPreTurb, PETSC_NULL_VEC, colVec, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call VecDestroy(colVec, ierr)
            call EChk(ierr, __FILE__, __LINE__)

        end subroutine applyTurbPCColumnScaling
    end subroutine FormJacobianANKTurb

    subroutine FormFunction_mf(ctx, inVec, rVec, ierr)

        ! This is the function used for the matrix-free matrix-vector products
        ! for the GMRES solver used in ANK

        use constants
        use blockPointers, only: nDom, volRef, il, jl, kl, dw, dtl
        use inputtimespectral, only: nTimeIntervalsSpectral
        use inputIteration, only: turbResScale
        use flowvarrefstate, only: nwf, nt1, nt2
        use NKSolver, only: setRvec
        use utils, only: setPointers, EChk
        use blockette, only: blocketteRes
        implicit none

        ! PETSc Variables
        PetscFortranAddr ctx(*)
        Vec inVec, rVec
        real(kind=realType) :: dtinv, rho
        integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii, iiRho
        real(kind=realType), pointer :: rvec_pointer(:)
        real(kind=realType), pointer :: invec_pointer(:)
        real(kind=realType), pointer :: wvec_pointer(:)
        logical :: useViscApprox

        ! get the input vector (column-scaled state, see getFullColScale)
        call setWANKScaled(inVec, 1, nState)

        ! determine if we want the approximate viscous fluxes
        useViscApprox = (.not. ANK_useFullVisc) .and. ANK_useDissApprox

        ! Determine if we want the turb residuals
        call blocketteRes(useDissApprox=ANK_useDissApprox, useViscApprox=useViscApprox, &
                          useTurbRes=ANK_coupled, useStoreWall=.False.)

        ! Copy the residuals to rVec in petsc
        if (ANK_coupled) then
            call setRVec(rVec)
        else
            call setRVecANK(rVec)
        end if

        ! rVec contains the full steady residual vector
        call VecGetArrayF90(rVec, rvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! inVec contains the perturbed state vector
        call VecGetArrayReadF90(inVec, invec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Also read the wVec to access the un-perturbed state vector.
        call VecGetArrayReadF90(wVec, wvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Multiply the perturbed state by the time step terms and add it to the residual
        call MatMultAdd(timeStepMat, inVec, rVec, rVec, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecRestoreArrayF90(rVec, rvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecRestoreArrayReadF90(wVec, wvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecRestoreArrayReadF90(inVec, invec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! We don't check an error here, so just pass back zero
        ierr = 0

    end subroutine FormFunction_mf

    subroutine FormFunction_mf_turb(ctx, inVec, rVec, ierr)

        ! This is the function used for the matrix-free matrix-vector products
        ! for the GMRES solver used in ANK

        use constants
        use blockPointers, only: nDom, volRef, il, jl, kl, dw, dtl, srcLambda
        use inputtimespectral, only: nTimeIntervalsSpectral
        use inputIteration, only: turbResScale, transitionSrcDtRestrict, transitionSrcDtLimit, &
                                   noBacktrackCount, srcDtDeactivateIters, transitionNK
        use flowvarrefstate, only: nwf, nt1, nt2, nwt
        use NKSolver, only: setRvec
        use utils, only: setPointers, EChk
        use blockette, only: blocketteRes
        implicit none

        ! PETSc Variables
        PetscFortranAddr ctx(*)
        Vec inVec, rVec
        real(kind=realType) :: dtinv, rho, dtinv_src
        integer(kind=intType) :: ierr, nn, sps, i, j, k, l, l1, ii, iiRho
        logical :: srcDtRestrictActive
        real(kind=realType), pointer :: rvec_pointer(:)
        real(kind=realType), pointer :: invec_pointer(:)
        real(kind=realType) :: cs(nt2 - nt1 + 1), rcf(nt2 - nt1 + 1)

        ! Derived condition for source dt restriction (ANK turbKSP path)
        srcDtRestrictActive = transitionNK .and. transitionSrcDtRestrict .and. (noBacktrackCount < srcDtDeactivateIters)

        ! Row scale (turbResScale) combined with the column scale of the
        ! state vector: the pseudo-time diagonal of the scaled system is
        ! dtinv * S_row / cs. For SA-Gamma-Retheta cs = turbResScale so the
        ! factor collapses to one; for all other models cs = 1 and the
        ! original turbResScale factor is recovered.
        call getTurbColScale(cs, nt2 - nt1 + 1)
        rcf = turbResScale(1:nt2 - nt1 + 1) / cs

        ! get the input vector (column-scaled turbulence state)
        call setWANKTurbScaled(inVec)

        call blocketteRes(useFlowRes=.False., useStoreWall=.False.)
        call setRVecANKTurb(rVec)

        ! Add the contribution from the time stepping term

        call VecGetArrayF90(rVec, rvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! inVec contains the perturbed state vector
        call VecGetArrayReadF90(inVec, invec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Include time step for turbulence
        ii = 1
        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                call setPointers(nn, 1_intType, sps)
                ! read the density residuals and set local CFL
                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il
                            ! needs to be modified
                            dtinv = one / (ANK_CFL * dtl(i, j, k) * volRef(i, j, k))

                            do l = nt1, nt2
                                l1 = l - nt1 + 1
                                ! Source dt restriction (P&Z Eq. 59): MAX form - must match PC above.
                                if (transitionSrcDtRestrict .and. srcDtRestrictActive .and. nwt == 3) then
                                    dtinv_src = srcLambda(i, j, k, l1) / transitionSrcDtLimit
                                    rvec_pointer(ii) = rvec_pointer(ii) + invec_pointer(ii) * &
                                                       max(dtinv, dtinv_src) * rcf(l1) / ANK_turbCFLScale
                                else
                                    rvec_pointer(ii) = rvec_pointer(ii) + invec_pointer(ii) * &
                                                       dtinv * rcf(l1) / ANK_turbCFLScale
                                end if
                                ii = ii + 1
                            end do
                        end do
                    end do
                end do
            end do
        end do

        call VecRestoreArrayF90(rVec, rvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecRestoreArrayReadF90(inVec, invec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! We don't check an error here, so just pass back zero
        ierr = 0

    end subroutine FormFunction_mf_turb

    subroutine computeUnsteadyResANK(omega)

        ! This routine calculates the unsteady residual in a given iteration.
        ! It needs the following variables/vectors:
        !
        !   omega:      This is the step size taken in the last update to the state
        !   deltaW:     Vector that contains the full update given from the
        !               Newton/Euler iteration.
        !   w(:,:,:,:): Should contain the updated state with the given step size
        !               lambdaLS and given update deltaW
        !   ANK_CFL:    The CFL number used for this non-linear iteration
        !   dtl:        Array containing time step values giving a CFL number of 1
        !               on each cell.
        !
        ! The routine calculates the unsteady residual and leaves the result in
        ! rVec, which was previously used to keep the steady residual only. This
        ! is done because the norm of this vector can easily be calculated with
        ! PETSc, however, after the line search, the rVec vector needs to be
        ! updated to contain only the steady state residuals. This can be done with
        ! setRVecANK/setRVec, with a dw(:,:,:,:) that is also up to date.

        use constants
        use blockPointers, only: nDom, volRef, il, jl, kl, w, dw, dtl
        use inputtimespectral, only: nTimeIntervalsSpectral
        use inputIteration, only: turbResScale
        use flowvarrefstate, only: nwf, nt1, nt2
        use NKSolver, only: setRvec
        use utils, only: setPointers, EChk
        use blockette, only: blocketteRes
        implicit none

        real(kind=realType), intent(in) :: omega

        real(kind=realType) :: dtinv, rho, uu, vv, ww
        integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii, iiRho
        real(kind=realType), pointer :: rvec_pointer(:)
        real(kind=realType), pointer :: dvec_pointer(:)

        real(kind=realType), dimension(nState, nState) :: timeStepMatrix
        real(kind=realType), dimension(nState) :: wPrev

        ! Allocate a PETSc vector like deltaW for intermediate computations
        Vec unsteadyVec

        call VecDuplicate(deltaW, unsteadyVec, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Calculate the steady residuals
        call blocketteRes(useTurbRes=ANK_coupled)
        if (ANK_coupled) then
            ! coupled rVec needs flow AND (turbResScale-scaled) turb rows;
            ! setRVecANK packs flow rows only and would misalign the vector
            call setRvec(rVec)
        else
            call setRVecANK(rVec)
        end if

        ! rVec contains the full steady residual vector
        call VecGetArrayF90(rVec, rvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! deltaW contains the full state update vector
        call VecGetArrayReadF90(deltaW, dvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! TODO AY: check if this routine is fine with complex mode...
        ! dtl and volume can both have complex values in them

        ! Multiply the delta by the time step terms
        call MatMult(timeStepMat, deltaW, unsteadyVec, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Add unsteady term to the steady residual
        call VecAXPY(rVec, -omega, unsteadyVec, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Deallocate the intermediate vector
        call VecDestroy(unsteadyVec, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecRestoreArrayF90(rVec, rvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecRestoreArrayReadF90(deltaW, dvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! We don't check an error here, so just pass back zero
        ierr = 0

    end subroutine computeUnsteadyResANK

    subroutine computeUnsteadyResANKTurb(omega)

        ! This routine calculates the unsteady residual in a given iteration.
        ! It needs the following variables/vectors:
        !
        !   omega:      This is the step size taken in the last update to the state
        !   deltaWTurb: Vector that contains the full update given from the
        !               Newton/Euler iteration.
        !   w(:,:,:,:): Should contain the updated state with the given step size
        !               lambdaLS and given update deltaW
        !   ANK_CFL:    The CFL number used for this non-linear iteration
        !   dtl:        Array containing time step values giving a CFL number of 1
        !               on each cell.
        !
        ! The routine calculates the unsteady residual and leaves the result in
        ! rVecTurb, which was previously used to keep the steady residual only. This
        ! is done because the norm of this vector can easily be calculated with
        ! PETSc, however, after the line search, the rVec vector needs to be
        ! updated to contain only the steady state residuals. This can be done with
        ! setRVecANK/setRVec, with a dw(:,:,:,:) that is also up to date.

        use constants
        use blockPointers, only: nDom, volRef, il, jl, kl, w, dw, dtl
        use inputtimespectral, only: nTimeIntervalsSpectral
        use inputIteration, only: turbResScale
        use flowvarrefstate, only: nwf, nt1, nt2
        use NKSolver, only: setRvec
        use utils, only: setPointers, EChk
        use blockette, only: blocketteRes
        implicit none

        real(kind=realType), intent(in) :: omega

        real(kind=realType) :: dtinv, rho, uu, vv, ww
        integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii, iiRho
        real(kind=realType), pointer :: rvec_pointer(:)
        real(kind=realType), pointer :: dvec_pointer(:)
        real(kind=realType) :: cs(nt2 - nt1 + 1), rcf(nt2 - nt1 + 1)

        ! TODO AY: check if this routine is fine in complex mode...
        ! dtl and volume can both have complex values in them

        ! deltaWTurb holds the column-scaled update; the pseudo-time term of
        ! the scaled system carries turbResScale / cs (see FormFunction_mf_turb).
        call getTurbColScale(cs, nt2 - nt1 + 1)
        rcf = turbResScale(1:nt2 - nt1 + 1) / cs

        ! Calculate the steady residuals
        call blocketteRes(useFlowRes=.False.)
        call setRVecANKTurb(rVecTurb)

        ! Add the contribution from the time stepping term

        call VecGetArrayF90(rVecTurb, rvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! deltaW contains the full update to the state
        call VecGetArrayReadF90(deltaWTurb, dvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Include time step for turbulence
        ii = 1
        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                call setPointers(nn, 1_intType, sps)
                ! read the density residuals and set local CFL
                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il
                            dtinv = one / (ANK_CFL * dtl(i, j, k) * volRef(i, j, k))

                            do l = nt1, nt2
                                ! turbulence variable needs additional scaling, and it may
                                ! get a different CFL number
                                rvec_pointer(ii) = rvec_pointer(ii) - omega * dvec_pointer(ii) * &
                                                   dtinv * rcf(l - nt1 + 1) / ANK_turbCFLScale
                                ii = ii + 1
                            end do
                        end do
                    end do
                end do
            end do
        end do

        call VecRestoreArrayF90(rVecTurb, rvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecRestoreArrayReadF90(deltaWTurb, dvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! We don't check an error here, so just pass back zero
        ierr = 0

    end subroutine computeUnsteadyResANKTurb

    subroutine destroyANKsolver

        ! Destroy all the PETSc objects for the approximate Newton-Krylov solver.

        use constants
        use utils, only: EChk
        use amg, only: destroyAMG
        implicit none
        integer(kind=intType) :: ierr

        if (ANK_SolverSetup) then

            call MatDestroy(dRdw, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call MatDestroy(dRdwPre, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call MatDestroy(timeStepMat, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call VecDestroy(wVec, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call VecDestroy(rVec, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call VecDestroy(deltaW, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call VecDestroy(baseRes, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call KSPDestroy(ANK_KSP, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call destroyAMG()

            ANK_SolverSetup = .False.

            if (ANK_turbSetup) then

                call MatDestroy(dRdwTurb, ierr)
                call EChk(ierr, __FILE__, __LINE__)

                call MatDestroy(dRdwPreTurb, ierr)
                call EChk(ierr, __FILE__, __LINE__)

                call VecDestroy(wVecTurb, ierr)
                call EChk(ierr, __FILE__, __LINE__)

                call VecDestroy(rVecTurb, ierr)
                call EChk(ierr, __FILE__, __LINE__)

                call VecDestroy(deltaWTurb, ierr)
                call EChk(ierr, __FILE__, __LINE__)

                call VecDestroy(baseResTurb, ierr)
                call EChk(ierr, __FILE__, __LINE__)

                call KSPDestroy(ANK_KSPTurb, ierr)
                call EChk(ierr, __FILE__, __LINE__)

                ANK_turbSetup = .False.
            end if
        end if
    end subroutine destroyANKsolver

    subroutine setWVecANK(wVec, lStart, lEnd)
        ! Set the current FLOW variables in the PETSc Vector

        use constants
        use blockPointers, only: nDom, il, jl, kl, w
        use inputtimespectral, only: ntimeIntervalsSpectral
        use utils, only: setPointers, EChk
        implicit none

        Vec wVec
        integer(kind=intType), intent(in) :: lStart, lEnd
        integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii
        real(kind=realType), pointer :: wvec_pointer(:)

        call VecGetArrayF90(wVec, wvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        ii = 0
        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                call setPointers(nn, 1_intType, sps)
                ! Copy off w to wVec
                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il
                            do l = lStart, lEnd
                                ii = ii + 1
                                wvec_pointer(ii) = w(i, j, k, l)
                            end do
                        end do
                    end do
                end do
            end do
        end do

        call VecRestoreArrayF90(wVec, wvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

    end subroutine setWVecANK

    subroutine setRVecANK(rVec)

        ! Set the current FLOW residual in dw into the PETSc Vector
        use constants
        use blockPointers, only: nDom, volRef, il, jl, kl, dw
        use inputtimespectral, only: nTimeIntervalsSpectral
        use flowvarrefstate, only: nwf, nt1, nt2
        use inputIteration, only: turbResScale
        use utils, only: setPointers, EChk
        implicit none
        Vec rVec
        integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii
        real(kind=realType), pointer :: rvec_pointer(:)
        real(Kind=realType) :: ovv
        call VecGetArrayF90(rVec, rvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        ii = 0
        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                call setPointers(nn, 1_intType, sps)
                ! Copy off dw/vol to rVec
                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il
                            ovv = one / volRef(i, j, k)
                            do l = 1, nwf
                                ii = ii + 1
                                rvec_pointer(ii) = dw(i, j, k, l) * ovv
                            end do
                        end do
                    end do
                end do
            end do
        end do

        call VecRestoreArrayF90(rVec, rvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

    end subroutine setRVecANK

    subroutine setRVecANKTurb(rVecTurb)

        ! Set the current Turb residual in dw into the PETSc Vector
        use constants
        use blockPointers, only: nDom, volRef, il, jl, kl, dw
        use inputtimespectral, only: nTimeIntervalsSpectral
        use flowvarrefstate, only: nt1, nt2
        use inputIteration, only: turbResScale
        use utils, only: setPointers, EChk
        implicit none
        Vec rVecTurb
        integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii
        real(kind=realType), pointer :: rvec_pointer(:)
        real(Kind=realType) :: ovv
        call VecGetArrayF90(rVecTurb, rvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        ii = 0
        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                call setPointers(nn, 1_intType, sps)
                ! Copy off dw/vol to rVec
                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il
                            ovv = one / volRef(i, j, k)
                            do l = nt1, nt2
                                ii = ii + 1
                                rvec_pointer(ii) = dw(i, j, k, l) * ovv * turbResScale(l - nt1 + 1)
                            end do
                        end do
                    end do
                end do
            end do
        end do

        call VecRestoreArrayF90(rVecTurb, rvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

    end subroutine setRVecANKTurb

    subroutine setWANK(wVec, lStart, lEnd)
        ! Get the updated solution from the PETSc Vector

        use constants
        use blockPointers, only: nDom, vol, il, jl, kl, w
        use inputtimespectral, only: nTimeIntervalsSpectral
        use utils, only: setPointers, EChk
        implicit none

        Vec wVec
        integer(kind=intType), intent(in) :: lStart, lEnd
        integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii
        real(kind=realType), pointer :: wvec_pointer(:)
        call VecGetArrayReadF90(wVec, wvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ii = 0
        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                call setPointers(nn, 1_intType, sps)

                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il
                            do l = lStart, lEnd
                                ii = ii + 1
                                w(i, j, k, l) = wvec_pointer(ii)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        call VecRestoreArrayReadF90(wVec, wvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

    end subroutine setWANK

    subroutine getTurbColScale(cs, nState)
        ! Column scale for the turbulence KSP state vector.
        !
        ! For SA-Gamma-Retheta the turbulence states span ~13 orders of
        ! magnitude in one PETSc vector (nuTilde ~1e-10..1e-4, gamma
        ! ~0.02..1, ReThetaTilde ~1e2..1e4). PETSc's matrix-free product
        ! computes a single differencing step from whole-vector norms
        ! (dominated by ReTheta), so without column scaling the
        ! nuTilde/gamma Jacobian columns are finite-difference noise
        ! whenever the pseudo-time diagonal does not dominate (high
        ! CFL / second-order phase). The state is scaled by turbResScale
        ! (~1/state magnitude): this is the column-scaling measure of
        ! Piotrowski & Zingg Sec. IV.1 (gamma_max, ReThetaT_max).
        !
        ! For every other turbulence model the factor is one, so the
        ! validated SA (and SST/kw) turbKSP behavior is unchanged.
        use constants
        use inputPhysics, only: turbModel
        use inputIteration, only: turbResScale, transitionNK
        implicit none
        integer(kind=intType), intent(in) :: nState
        real(kind=realType), intent(out) :: cs(nState)

        cs = one
        if (turbModel == spalartallmarasnoft2gammaretheta .and. transitionNK) then
            cs(1:nState) = turbResScale(1:nState)
        end if
    end subroutine getTurbColScale

    subroutine setWVecANKTurbScaled(wVec)
        ! Pack the turbulence state into the PETSc vector in COLUMN-SCALED
        ! form: wVec = w * cs. Counterpart of setWANKTurbScaled.
        use constants
        use blockPointers, only: nDom, il, jl, kl, w
        use inputtimespectral, only: ntimeIntervalsSpectral
        use flowVarRefState, only: nt1, nt2
        use utils, only: setPointers, EChk
        implicit none

        Vec wVec
        integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii
        real(kind=realType), pointer :: wvec_pointer(:)
        real(kind=realType) :: cs(nt2 - nt1 + 1)

        call getTurbColScale(cs, nt2 - nt1 + 1)

        call VecGetArrayF90(wVec, wvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        ii = 0
        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                call setPointers(nn, 1_intType, sps)
                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il
                            do l = nt1, nt2
                                ii = ii + 1
                                wvec_pointer(ii) = w(i, j, k, l) * cs(l - nt1 + 1)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        call VecRestoreArrayF90(wVec, wvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

    end subroutine setWVecANKTurbScaled

    subroutine setWANKTurbScaled(wVec)
        ! Unpack the COLUMN-SCALED turbulence state from the PETSc vector:
        ! w = wVec / cs. Counterpart of setWVecANKTurbScaled.
        use constants
        use blockPointers, only: nDom, il, jl, kl, w
        use inputtimespectral, only: ntimeIntervalsSpectral
        use flowVarRefState, only: nt1, nt2
        use utils, only: setPointers, EChk
        implicit none

        Vec wVec
        integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii
        real(kind=realType), pointer :: wvec_pointer(:)
        real(kind=realType) :: cs(nt2 - nt1 + 1), csInv(nt2 - nt1 + 1)

        call getTurbColScale(cs, nt2 - nt1 + 1)
        csInv = one / cs

        call VecGetArrayReadF90(wVec, wvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        ii = 0
        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                call setPointers(nn, 1_intType, sps)
                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il
                            do l = nt1, nt2
                                ii = ii + 1
                                w(i, j, k, l) = wvec_pointer(ii) * csInv(l - nt1 + 1)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        call VecRestoreArrayReadF90(wVec, wvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

    end subroutine setWANKTurbScaled

    subroutine getFullColScale(fac, lStart, lEnd)
        ! Column scale for the (possibly coupled) ANK state vector: one for
        ! the mean-flow entries, getTurbColScale for the turbulence entries.
        ! For models other than SA-Gamma-Retheta every factor is one.
        use constants
        use inputPhysics, only: turbModel
        use inputIteration, only: turbResScale, transitionNK, ankColScaleUnit
        use flowVarRefState, only: nt1, nt2
        implicit none
        integer(kind=intType), intent(in) :: lStart, lEnd
        real(kind=realType), intent(out) :: fac(lStart:lEnd)
        integer(kind=intType) :: l

        fac = one
        if (turbModel == spalartallmarasnoft2gammaretheta .and. transitionNK) then
            do l = max(lStart, nt1), min(lEnd, nt2)
                fac(l) = turbResScale(l - nt1 + 1)
            end do
        end if

        ! VERIF_06 F10: replace the above with factors that put EVERY variable
        ! at a unit RMS, so the single MFFD differencing step h is equally
        ! appropriate for all of them. Only applied once the factors have
        ! actually been measured (they are refreshed at PC-reform cadence);
        ! until then the legacy scaling above stands, so the first iterations
        ! behave exactly as before.
        if (ankColScaleUnit) then
            do l = lStart, min(lEnd, size(colScaleUnitFac))
                if (colScaleUnitFac(l) > zero) fac(l) = colScaleUnitFac(l)
            end do
        end if
    end subroutine getFullColScale

    subroutine applyANKAlgorithm2Damping(wv, dw, lam)
        ! Paper Algorithm 2 (P&Z 2020 §IV.B.2) for the COUPLED ANK path.
        !
        ! VERIF_06 F0/F7: the thesis's inexact-Newton phase is (T + A_2nd)dQ =
        ! -R, which is ADflow's CSANK -- not ADflow's NK. All the other ported
        ! machinery (Eq. 58 S_r/S_a/S_c, Eq. 59 diagonal + freeze + deactivation)
        ! is already present in both paths, but the per-node gamma/Re-theta-t
        ! damping existed ONLY in NK (single call site in NKStep). This supplies
        ! the missing half so CANK/CSANK runs the paper's algorithm.
        !
        ! Identical back-off to applyNKAlgorithm2Damping and to the DD-ADI loop
        ! in saGammaReThetaSolve: same bounds (rsaGRgammaLo/Hi, rsaGRreThetaLo),
        ! same transitionDampTheta/transitionDampMaxIter options.
        !
        ! Difference from the NK version, and the reason for the extra argument:
        ! NK carries the pre-step state in a separate Vec, whereas ANKStep has
        ! already applied VecAXPY(wv, -lam, dw). The pre-step value is therefore
        ! reconstructed per node as wv + lam*dw. When lam == 0 (rejected step)
        ! the delta is zero and this is a no-op, as it should be.
        !
        ! Operates in the column-scaled space (getFullColScale): damping a
        ! linear interpolation by a scalar is scale-invariant, so only the
        ! bounds check needs the physical value.
        use constants
        use paramTurb, only: rsaGRgammaLo, rsaGRgammaHi, rsaGRreThetaLo
        use inputIteration, only: transitionDampTheta, transitionDampMaxIter, solverStallDiag
        use flowVarRefState, only: nw, nt1
        use communication, only: myid
        use utils, only: EChk
        implicit none

        Vec wv, dw
        real(kind=realType), intent(in) :: lam
        integer(kind=intType) :: ierr, jj, nCells, mm
        integer(kind=intType) :: nDampCapGamma, nDampCapReTheta
        integer(kind=intType) :: nSoftDampGamma, nSoftDampReTheta
        real(kind=realType) :: minDampFactorGamma, minDampFactorReTheta
        real(kind=realType), pointer :: wPtr(:), dPtr(:)
        real(kind=realType) :: cs(1:nw)
        real(kind=realType) :: xOld, deltaScaled, dampFactor, candScaled, candPhys
        integer(kind=intType) :: gammaOff, rethetaOff

        call getFullColScale(cs, 1_intType, nw)
        gammaOff = nt1 + 1
        rethetaOff = nt1 + 2

        call VecGetArrayF90(wv, wPtr, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        call VecGetArrayReadF90(dw, dPtr, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        nDampCapGamma = 0
        nDampCapReTheta = 0
        nSoftDampGamma = 0
        nSoftDampReTheta = 0
        minDampFactorGamma = one
        minDampFactorReTheta = one
        nCells = size(wPtr) / nw

        do jj = 0, nCells - 1
            ! Gamma: exponential back-off until in [gammaLo, gammaHi]
            deltaScaled = -lam * dPtr(jj * nw + gammaOff)
            xOld = wPtr(jj * nw + gammaOff) - deltaScaled
            dampFactor = one
            candScaled = xOld + dampFactor * deltaScaled
            do mm = 1, transitionDampMaxIter
                candPhys = candScaled / cs(gammaOff)
                if (candPhys >= rsaGRgammaLo .and. candPhys <= rsaGRgammaHi) exit
                dampFactor = dampFactor * transitionDampTheta
                candScaled = xOld + dampFactor * deltaScaled
            end do
            candPhys = candScaled / cs(gammaOff)
            if (candPhys < rsaGRgammaLo .or. candPhys > rsaGRgammaHi) then
                nDampCapGamma = nDampCapGamma + 1
                candPhys = min(max(candPhys, rsaGRgammaLo), rsaGRgammaHi)
                candScaled = candPhys * cs(gammaOff)
            end if
            if (dampFactor < one) then
                nSoftDampGamma = nSoftDampGamma + 1
                minDampFactorGamma = min(minDampFactorGamma, dampFactor)
            end if
            wPtr(jj * nw + gammaOff) = candScaled

            ! Re-theta-t: back-off until >= rsaGRreThetaLo (lower bound only)
            deltaScaled = -lam * dPtr(jj * nw + rethetaOff)
            xOld = wPtr(jj * nw + rethetaOff) - deltaScaled
            dampFactor = one
            candScaled = xOld + dampFactor * deltaScaled
            do mm = 1, transitionDampMaxIter
                candPhys = candScaled / cs(rethetaOff)
                if (candPhys >= rsaGRreThetaLo) exit
                dampFactor = dampFactor * transitionDampTheta
                candScaled = xOld + dampFactor * deltaScaled
            end do
            candPhys = candScaled / cs(rethetaOff)
            if (candPhys < rsaGRreThetaLo) then
                nDampCapReTheta = nDampCapReTheta + 1
                candScaled = rsaGRreThetaLo * cs(rethetaOff)
            end if
            if (dampFactor < one) then
                nSoftDampReTheta = nSoftDampReTheta + 1
                minDampFactorReTheta = min(minDampFactorReTheta, dampFactor)
            end if
            wPtr(jj * nw + rethetaOff) = candScaled
        end do

        call VecRestoreArrayF90(wv, wPtr, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        call VecRestoreArrayReadF90(dw, dPtr, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        if (nDampCapGamma + nDampCapReTheta > 0) then
            print *, 'Warning: ANK Algorithm 2 damping exhausted transitionDampMaxIter (', &
                transitionDampMaxIter, ') in ', nDampCapGamma, ' gamma / ', &
                nDampCapReTheta, ' reTheta cells on proc ', myid, '; values clipped.'
        end if

        ! Hand the soft-damping activity to the stall diagnostics rather than
        ! printing per-proc lines: it is reported on the single STALLDIAG ANK
        ! line so stall/oscillation causes stay readable in one place.
        stallSoftDampG = nSoftDampGamma
        stallSoftDampR = nSoftDampReTheta
        stallMinDampG = minDampFactorGamma
        stallMinDampR = minDampFactorReTheta

    end subroutine applyANKAlgorithm2Damping

    subroutine reportScaledStateStats(wv, dw)
        ! VERIF_06 F9 metric: per-variable statistics of the COLUMN-SCALED
        ! state and update.
        !
        ! Why this matters. ANK/CANK and NK all build their Jacobian-vector
        ! products with PETSc MFFD (MatCreateMFFD; MatMFFDSetType and
        ! MatMFFDSetFunctionError are never called, so it is Dennis-Schnabel
        ! with an assumed function error of machine epsilon). MFFD uses ONE
        ! scalar differencing step h for the whole vector:
        !     J*a ~ [F(u + h*a) - F(u)] / h
        ! so the perturbation seen by component i is h*a_i. If the scaled
        ! components span orders of magnitude, no single h can be small
        ! relative to the large components AND large relative to the noise
        ! floor of the small ones -- the columns belonging to the small
        ! components come back as noise, and the Newton DIRECTION is wrong in
        ! exactly those equations even though the linear solve converges
        ! beautifully. That is the observed signature: Step = 1.00, linear
        ! residual 0.02-0.04, and no progress.
        !
        ! ADflow reuses turbResScale for BOTH the row (residual) scaling in
        ! setRVec and the column (variable) scaling here, whereas the thesis
        ! keeps S_r and S_c as separate reciprocal roles (Eq. 3.6-3.8). Values
        ! chosen to normalise residuals need not normalise variables.
        !
        ! This routine reports what the scaling actually achieves, so the fix
        ! is chosen from data rather than from argument.
        use constants
        use flowVarRefState, only: nw, nt1, nt2
        use inputIteration, only: ankColScaleUnit
        use communication, only: ADflow_comm_world, myid
        use utils, only: EChk
        implicit none

        Vec wv, dw
        integer(kind=intType) :: ierr, jj, nCells, l
        real(kind=realType), pointer :: wPtr(:), dPtr(:)
        real(kind=realType) :: cs(1:nw)
        real(kind=alwaysRealType) :: sw(nw), sd(nw), mw(nw)
        real(kind=alwaysRealType) :: swG(nw), sdG(nw), mwG(nw)
        integer(kind=intType) :: nTot, nTotG
        character(len=8) :: vname
        real(kind=alwaysRealType) :: rmsL

        call getFullColScale(cs, 1_intType, nw)
        call VecGetArrayReadF90(wv, wPtr, ierr); call EChk(ierr, __FILE__, __LINE__)
        call VecGetArrayReadF90(dw, dPtr, ierr); call EChk(ierr, __FILE__, __LINE__)
        ! NOTE: both arrays are restored before the scaling update below, which
        ! repacks wv -- PETSc forbids writing a vector whose array is checked out.

        sw = zero; sd = zero; mw = zero
        nCells = size(wPtr) / nw
        do jj = 0, nCells - 1
            do l = 1, nw
                sw(l) = sw(l) + wPtr(jj * nw + l)**2
                sd(l) = sd(l) + dPtr(jj * nw + l)**2
                mw(l) = max(mw(l), abs(wPtr(jj * nw + l)))
            end do
        end do
        nTot = nCells

        call VecRestoreArrayReadF90(wv, wPtr, ierr); call EChk(ierr, __FILE__, __LINE__)
        call VecRestoreArrayReadF90(dw, dPtr, ierr); call EChk(ierr, __FILE__, __LINE__)

        call mpi_allreduce(sw, swG, nw, MPI_DOUBLE, mpi_sum, ADflow_comm_world, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        call mpi_allreduce(sd, sdG, nw, MPI_DOUBLE, mpi_sum, ADflow_comm_world, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        call mpi_allreduce(mw, mwG, nw, MPI_DOUBLE, mpi_max, ADflow_comm_world, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        call mpi_allreduce(nTot, nTotG, 1_intType, MPI_INTEGER, mpi_sum, ADflow_comm_world, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! VERIF_06 F10: fold the measured RMS back into the column scale so
        ! every variable ends up at unit RMS. rmsW is the RMS of the ALREADY
        ! scaled state, so the update is new = old / rmsW, which drives the
        ! scaled RMS to one. Done here, at PC-reform cadence, so the assembled
        ! preconditioner and the operator never disagree about the scaling.
        if (ankColScaleUnit) then
            do l = 1, min(nw, size(colScaleUnitFac))
                rmsL = sqrt(swG(l) / real(max(nTotG, 1), alwaysRealType))
                ! Clamp the per-update change. A single step can only move a
                ! factor by 10x either way, which keeps a bad measurement from
                ! running the scaling away (the first implementation had no
                ! clamp and no repack, and drove nuTilde's factor to 3e13).
                rmsL = min(max(rmsL, 0.1_alwaysRealType), 10.0_alwaysRealType)
                colScaleUnitFac(l) = cs(l) / rmsL
            end do
            ! The factors have changed, so wVec -- packed with the OLD ones --
            ! is now inconsistent with them. Repack from ADflow's w, which is
            ! the authoritative unscaled state, before anything downstream
            ! (FormJacobianANK, the MFFD base) uses either.
            call setWVecANKScaled(wv, 1_intType, nw)
        end if

        if (myid == 0) then
            do l = 1, nw
                if (l < nt1) then
                    write (vname, "(a,i1)") 'flow', l
                else if (l == nt1) then
                    vname = 'nuTilde'
                else if (l == nt1 + 1) then
                    vname = 'gamma'
                else
                    vname = 'reTheta'
                end if
                ! rmsW/maxW are the SCALED magnitudes a single MFFD h has to
                ! serve simultaneously; spread across rows is the diagnosis.
                write (*, "(a,a8,a,es10.3,a,es10.3,a,es10.3,a,es10.3)") &
                    " SCALEDIAG ", vname, &
                    "  colScale=", cs(l), &
                    "  rmsW=", sqrt(swG(l) / real(max(nTotG, 1), alwaysRealType)), &
                    "  maxW=", mwG(l), &
                    "  rmsDelta=", sqrt(sdG(l) / real(max(nTotG, 1), alwaysRealType))
            end do
        end if

    end subroutine reportScaledStateStats

    subroutine setWVecANKScaled(wVec, lStart, lEnd)
        ! Pack the state into the PETSc vector in COLUMN-SCALED form:
        ! wVec = w * fac (fac = 1 for flow entries and for non-transition
        ! models, so the validated behavior is unchanged there).
        use constants
        use blockPointers, only: nDom, il, jl, kl, w
        use inputtimespectral, only: ntimeIntervalsSpectral
        use utils, only: setPointers, EChk
        implicit none

        Vec wVec
        integer(kind=intType), intent(in) :: lStart, lEnd
        integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii
        real(kind=realType), pointer :: wvec_pointer(:)
        real(kind=realType) :: fac(lStart:lEnd)

        call getFullColScale(fac, lStart, lEnd)

        call VecGetArrayF90(wVec, wvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        ii = 0
        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                call setPointers(nn, 1_intType, sps)
                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il
                            do l = lStart, lEnd
                                ii = ii + 1
                                wvec_pointer(ii) = w(i, j, k, l) * fac(l)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        call VecRestoreArrayF90(wVec, wvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

    end subroutine setWVecANKScaled

    subroutine setWANKScaled(wVec, lStart, lEnd)
        ! Unpack the COLUMN-SCALED state from the PETSc vector: w = wVec / fac.
        ! Counterpart of setWVecANKScaled.
        use constants
        use blockPointers, only: nDom, il, jl, kl, w
        use inputtimespectral, only: ntimeIntervalsSpectral
        use utils, only: setPointers, EChk
        implicit none

        Vec wVec
        integer(kind=intType), intent(in) :: lStart, lEnd
        integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii
        real(kind=realType), pointer :: wvec_pointer(:)
        real(kind=realType) :: fac(lStart:lEnd), facInv(lStart:lEnd)

        call getFullColScale(fac, lStart, lEnd)
        facInv = one / fac

        call VecGetArrayReadF90(wVec, wvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        ii = 0
        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                call setPointers(nn, 1_intType, sps)
                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il
                            do l = lStart, lEnd
                                ii = ii + 1
                                w(i, j, k, l) = wvec_pointer(ii) * facInv(l)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        call VecRestoreArrayReadF90(wVec, wvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

    end subroutine setWANKScaled

    subroutine applyANKColumnScaling(matrix)
        ! Right-multiply an assembled coupled ANK matrix (dRdwPre or
        ! timeStepMat) by diag(1/fac) so it is consistent with the
        ! column-scaled state vector. Only called for SA-Gamma-Retheta.
        use constants
        use flowVarRefState, only: nt1, nt2
        use utils, only: EChk
        implicit none

        Mat matrix
        Vec colVec
        integer(kind=intType) :: ierr, jj, lState
        real(kind=realType), pointer :: cvec_pointer(:)
        real(kind=realType) :: fac(1:nState)

        call getFullColScale(fac, 1_intType, nState)

        call VecDuplicate(wVec, colVec, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecGetArrayF90(colVec, cvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        do jj = 0, size(cvec_pointer) / nState - 1
            do lState = 1, nState
                cvec_pointer(jj * nState + lState) = one / fac(lState)
            end do
        end do
        call VecRestoreArrayF90(colVec, cvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call MatDiagonalScale(matrix, PETSC_NULL_VEC, colVec, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecDestroy(colVec, ierr)
        call EChk(ierr, __FILE__, __LINE__)

    end subroutine applyANKColumnScaling

    subroutine physicalityCheckANK(lambdaP)

        use constants
        use blockPointers, only: ndom, il, jl, kl
        use flowVarRefState, only: nw, nwf, nt1, nt2
        use inputPhysics, only: turbModel
        use inputIteration, only: turbResScale, solverStallDiag, ankTransitionGlobalLambda
        use paramTurb, only: rsaGRgammaLo, rsaGRgammaHi, rsaGRreThetaLo
        use inputtimespectral, only: nTimeIntervalsSpectral
        use utils, only: setPointers, EChk
        use genericISNAN, only: myisnan
        use communication, only: ADflow_comm_world, myid
        implicit none

        ! input variable
        real(kind=realType), intent(inout) :: lambdaP

        ! local variables
        integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii
        real(kind=realType), pointer :: wvec_pointer(:)
        real(kind=realType), pointer :: dvec_pointer(:)
        real(kind=alwaysRealType) :: lambdaL ! L is for local
        real(kind=alwaysRealType) :: lambdaP_recv ! to receive the global step
        real(kind=alwaysRealType) :: ratio
        real(kind=alwaysRealType) :: wval, dval, ratioBound

        ! Stall diagnostics (VERIF_06): track which variable/cell produced the
        ! binding ratio. Only maintained when solverStallDiag is on.
        integer(kind=intType) :: bindVar, bindLoc(4), minRank
        real(kind=alwaysRealType) :: diagIn(2), diagOut(2)

        ! Determine the maximum step size that would yield
        ! a maximum relative change of ANK_physLSTol in density, and total energy.
        ! We also check for turbulence, but only limit the step
        ! for updates that decrease the value of the turbulence working variable.

        ! Initialize the local step size as lambdaP which is an i/o variable
        lambdaL = real(lambdaP)

        ! First we need to read both the update and the state
        ! from PETSc because the w in ADFlow currently contains
        ! the state that is perturbed during the matrix-free
        ! operations.

        ! wVec contains the state vector
        call VecGetArrayF90(wVec, wvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! deltaW contains the full update
        call VecGetArrayF90(deltaW, dvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Stall diagnostics: nothing has bound the step yet.
        bindVar = 0
        bindLoc = 0

        ! in decoupled, we just have the flow variables
        if (.not. ANK_coupled) then
            ii = 1
            do nn = 1, nDom
                do sps = 1, nTimeIntervalsSpectral
                    call setPointers(nn, 1_intType, sps)
                    do k = 2, kl
                        do j = 2, jl
                            do i = 2, il
                                ! multiply the ratios by ANK_physLSTol to check if the change in a
                                ! variable is greater than ANK_physLSTol of the variable itself.

                                ! check density
#ifndef USE_COMPLEX
                                ! to have the real mode slightly more efficient, we do not check if variables are real
                                ratio = abs(wvec_pointer(ii) / (dvec_pointer(ii) + eps)) * ANK_physLSTol
#else
                                ! We dont care what happens to the complex part of the update because
                                ! that is a linear system. So again check the real update for the physical
                                ! line search. Towards the end of the simulation, real part gets smaller and
                                ! and smaller updates, so this routine will always give a step of 1 which is what
                                ! we want for the complex parts.
                                ratio = abs(real(wvec_pointer(ii)) / real(dvec_pointer(ii) + eps)) * real(ANK_physLSTol)
#endif
                                if (solverStallDiag .and. ratio < lambdaL) then
                                    bindVar = 1; bindLoc = (/nn, i, j, k/)
                                end if
                                lambdaL = min(lambdaL, ratio)

                                ! increment by 4 because we want to skip momentum variables
                                ii = ii + 4

                                ! check energy
#ifndef USE_COMPLEX
                                ! see the comment above for the difference between real and complex versions
                                ratio = abs(wvec_pointer(ii) / (dvec_pointer(ii) + eps)) * ANK_physLSTol
#else
                                ratio = abs(real(wvec_pointer(ii)) / real(dvec_pointer(ii) + eps)) * real(ANK_physLSTol)
#endif
                                if (solverStallDiag .and. ratio < lambdaL) then
                                    bindVar = 2; bindLoc = (/nn, i, j, k/)
                                end if
                                lambdaL = min(lambdaL, ratio)
                                ii = ii + 1
                            end do
                        end do
                    end do
                end do
            end do
            ! in coupled, we also have the turbulence variables
        else
            ii = 1
            do nn = 1, nDom
                do sps = 1, nTimeIntervalsSpectral
                    call setPointers(nn, 1_intType, sps)
                    do k = 2, kl
                        do j = 2, jl
                            do i = 2, il

                                ! multiply the ratios by ANK_physLSTol to check if the change in a
                                ! variable is greater than ANK_physLSTol of the variable itself.

                                ! check density
#ifndef USE_COMPLEX
                                ! to have the real mode slightly more efficient, we do not check if variables are real
                                ratio = abs(wvec_pointer(ii) / (dvec_pointer(ii) + eps)) * ANK_physLSTol
#else
                                ! We dont care what happens to the complex part of the update because
                                ! that is a linear system. So again check the real update for the physical
                                ! line search. Towards the end of the simulation, real part gets smaller and
                                ! and smaller updates, so this routine will always give a step of 1 which is what
                                ! we want for the complex parts.
                                ratio = abs(real(wvec_pointer(ii)) / real(dvec_pointer(ii) + eps)) * real(ANK_physLSTol)
#endif
                                if (solverStallDiag .and. ratio < lambdaL) then
                                    bindVar = 1; bindLoc = (/nn, i, j, k/)
                                end if
                                lambdaL = min(lambdaL, ratio)

                                ! increment by 4 because we want to skip momentum variables
                                ii = ii + 4

                                ! check energy
#ifndef USE_COMPLEX
                                ! see the comment above for the difference between real and complex versions
                                ratio = abs(wvec_pointer(ii) / (dvec_pointer(ii) + eps)) * ANK_physLSTol
#else
                                ratio = abs(real(wvec_pointer(ii)) / real(dvec_pointer(ii) + eps)) * real(ANK_physLSTol)
#endif
                                if (solverStallDiag .and. ratio < lambdaL) then
                                    bindVar = 2; bindLoc = (/nn, i, j, k/)
                                end if
                                lambdaL = min(lambdaL, ratio)
                                ii = ii + 1

                                ! if coupled ank is used, nstate = nw and this loop is executed
                                ! if no turbulence variables, this loop will be automatically skipped
                                ! check turbulence variable
#ifndef USE_COMPLEX
                                ratio = (wvec_pointer(ii) / (dvec_pointer(ii) + eps)) * ANK_physLSTolTurb
#else
                                ratio = (real(wvec_pointer(ii)) &
                                         / real(dvec_pointer(ii) + eps)) * real(ANK_physLSTolTurb)
#endif
                                ! if the ratio is less than min step, the update is either
                                ! in the positive direction, therefore we do not clip it,
                                ! or the update is very limiting, so we just clip the
                                ! individual update for this cell.
                                if (ratio .lt. ANK_stepFactor * ANK_stepMin) then
                                    ! The update was very limiting, so just clip this
                                    ! individual update and dont change the overall
                                    ! step size. To select the new update, instead of
                                    ! clipping to zero, we clip to 1 percent of the original.
                                    if (ratio .gt. zero) &
                                        dvec_pointer(ii) = wvec_pointer(ii) * ANK_physLSTolTurb

                                    ! Either case, set the ratio to one. Positive updates
                                    ! do not limit the step, negative updates below minimum
                                    ! step were already clipped.
                                    ratio = one
                                end if
                                if (solverStallDiag .and. ratio < lambdaL) then
                                    bindVar = 3; bindLoc = (/nn, i, j, k/)
                                end if
                                lambdaL = min(lambdaL, ratio)
                                ii = ii + 1

                                ! Physicality checks for gamma and ReTheta
                                do l = nt1 + 1, nt2
#ifndef USE_COMPLEX
                                    ratio = (wvec_pointer(ii) / (dvec_pointer(ii) + eps)) * ANK_physLSTolTurb
                                    wval = wvec_pointer(ii)
                                    dval = dvec_pointer(ii)
#else
                                    ratio = (real(wvec_pointer(ii)) &
                                             / real(dvec_pointer(ii) + eps)) * real(ANK_physLSTolTurb)
                                    wval = real(wvec_pointer(ii))
                                    dval = real(dvec_pointer(ii))
#endif

                                    ! For SA-gamma-ReTheta, also enforce variable bounds by
                                    ! limiting the global step with the local admissible step.
                                    ! wVec/deltaW hold the COLUMN-SCALED state/update, so the
                                    ! absolute bounds are scaled to match (relative ratios are
                                    ! scale-invariant).
                                    ! NOTE (2026-07-14): per-cell Algorithm-2 damping (scale only
                                    ! the offending cell's update, never limit the global step)
                                    ! was tried and REVERTED: full global steps with locally
                                    ! damped front cells made the update inconsistent across the
                                    ! transition front — gamma residual bounced and wall-time
                                    ! progress was worse than the global-lambda throttle
                                    ! (paper_mimic run 2 vs run 1).
                                    if (turbModel == spalartallmarasnoft2gammaretheta) then
                                        ratioBound = one

                                        if (l == nt1 + 1) then
                                            ! gamma lower and upper bounds
                                            if (dval > zero) then
                                                ratioBound = min(ratioBound, &
                                                    (wval - rsaGRgammaLo * turbResScale(2)) / (dval + eps))
                                            else if (dval < zero) then
                                                ratioBound = min(ratioBound, &
                                                    (rsaGRgammaHi * turbResScale(2) - wval) / (-dval + eps))
                                            end if
                                        else if (l == nt1 + 2) then
                                            ! ReTheta lower bound
                                            if (dval > zero) then
                                                ratioBound = min(ratioBound, &
                                                    (wval - rsaGRreThetaLo * turbResScale(3)) / (dval + eps))
                                            end if
                                        end if

                                        ratioBound = max(ratioBound, zero)
                                        ratio = min(ratio, ratioBound)
                                    end if

                                    if (ratio .lt. ANK_stepFactor * ANK_stepMin) then
                                        if (ratio .gt. zero) &
                                            dvec_pointer(ii) = wvec_pointer(ii) * ANK_physLSTolTurb
                                        ratio = one
                                    end if
                                    if (solverStallDiag .and. ratio < lambdaL) then
                                        bindVar = 3 + (l - nt1); bindLoc = (/nn, i, j, k/)
                                    end if
                                    ! VERIF_06 F8: gamma/Re-theta-t take part in
                                    ! the GLOBAL step limit only when
                                    ! ankTransitionGlobalLambda is on. Measured
                                    ! on the swept wing, gamma was the binding
                                    ! variable in 7143 of ~7200 coupled
                                    ! iterations -- a few front cells throttling
                                    ! the whole field, which is exactly what the
                                    ! paper's per-node Algorithm 2 avoids. With
                                    ! the switch off the per-cell clipping above
                                    ! still runs; the bounds are then enforced by
                                    ! applyANKAlgorithm2Damping instead.
                                    if (ankTransitionGlobalLambda) lambdaL = min(lambdaL, ratio)
                                    ii = ii + 1
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end if

        ! Restore the pointers to PETSc vectors

        call VecRestoreArrayF90(wVec, wvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecRestoreArrayF90(deltaW, dvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Make sure that we did not get any NaN's in the process
        if (myisnan(lambdaL)) then
            lambdaL = zero
        end if

        ! Finally, communicate the step size across processes and return
        ! mpi allreduce is not defined for complex numbers with the min operation
        ! so we will use the lambdaP_recv variable to receive
        call mpi_allreduce(lambdaL, lambdaP_recv, 1_intType, MPI_DOUBLE, &
                           mpi_min, ADflow_comm_world, ierr)
        call EChk(ierr, __FILE__, __LINE__)

#ifndef USE_COMPLEX
        lambdaP = lambdaP_recv
#else
        ! finally, as a safety check, purge the complex part of lambda
        lambdaP = cmplx(lambdaP_recv, 0.0_realType)
#endif

        ! Stall diagnostics: find which rank owns the globally binding cell
        ! (MINLOC on the local lambdas) and broadcast its variable/location so
        ! every rank -- and hence the root's log line -- reports the same cell.
        if (solverStallDiag) then
            diagIn(1) = lambdaL
            diagIn(2) = real(myid, alwaysRealType)
            call mpi_allreduce(diagIn, diagOut, 1_intType, MPI_2DOUBLE_PRECISION, &
                               MPI_MINLOC, ADflow_comm_world, ierr)
            call EChk(ierr, __FILE__, __LINE__)
            minRank = int(diagOut(2), intType)
            call mpi_bcast(bindVar, 1_intType, MPI_INTEGER, minRank, ADflow_comm_world, ierr)
            call EChk(ierr, __FILE__, __LINE__)
            call mpi_bcast(bindLoc, 4_intType, MPI_INTEGER, minRank, ADflow_comm_world, ierr)
            call EChk(ierr, __FILE__, __LINE__)
            stallBindVar = bindVar
            stallBindLoc = bindLoc
            stallLamPhys = lambdaP_recv
        end if

    end subroutine physicalityCheckANK

    subroutine pzPhysicsRestriction(lambdaP)

        ! Thesis Algorithm 2 (physics-based restriction), PZ mode only:
        ! delta_phys = min over all cells with a NEGATIVE density or energy
        ! update of 0.90*Q/|dQ| (Eq. 3.10). The update applied later is
        ! -lambdaP*deltaW, so a negative physical update corresponds to
        ! deltaW > 0. Operates on the column-scaled wVec/deltaW; the ratio is
        ! scale-invariant because rho and E share one positive scale factor.
        ! The caller handles delta < 1% (reject + halve reference CFL, or clip
        ! at the CFL floor).

        use constants
        use utils, only: EChk
        use genericISNAN, only: myisnan
        use communication, only: ADflow_comm_world
        implicit none

        real(kind=realType), intent(inout) :: lambdaP

        integer(kind=intType) :: ierr, ii, vecSize
        real(kind=realType), pointer :: wvec_pointer(:)
        real(kind=realType), pointer :: dvec_pointer(:)
        real(kind=alwaysRealType) :: deltaL, deltaG, wval, dval

        call VecGetArrayF90(wVec, wvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        call VecGetArrayF90(deltaW, dvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        deltaL = real(lambdaP)

        ! The local state vector is nState entries per cell:
        ! [rho, mx, my, mz, E, (turb...)] — check rho (offset 0) and E (offset 4).
        vecSize = size(wvec_pointer)
        do ii = 1, vecSize, nState
            ! density
#ifndef USE_COMPLEX
            wval = wvec_pointer(ii)
            dval = dvec_pointer(ii)
#else
            wval = real(wvec_pointer(ii))
            dval = real(dvec_pointer(ii))
#endif
            if (dval > eps .and. wval > zero) then
                deltaL = min(deltaL, 0.90_alwaysRealType * wval / dval)
            end if

            ! total energy
#ifndef USE_COMPLEX
            wval = wvec_pointer(ii + 4)
            dval = dvec_pointer(ii + 4)
#else
            wval = real(wvec_pointer(ii + 4))
            dval = real(dvec_pointer(ii + 4))
#endif
            if (dval > eps .and. wval > zero) then
                deltaL = min(deltaL, 0.90_alwaysRealType * wval / dval)
            end if
        end do

        call VecRestoreArrayF90(wVec, wvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        call VecRestoreArrayF90(deltaW, dvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        if (myisnan(deltaL)) deltaL = zero

        call mpi_allreduce(deltaL, deltaG, 1_intType, MPI_DOUBLE, &
                           mpi_min, ADflow_comm_world, ierr)
        call EChk(ierr, __FILE__, __LINE__)

#ifndef USE_COMPLEX
        lambdaP = deltaG
#else
        lambdaP = cmplx(deltaG, 0.0_realType)
#endif

    end subroutine pzPhysicsRestriction

    subroutine physicalityCheckANKTurb(lambdaP)

        use constants
        use blockPointers, only: ndom, il, jl, kl
        use flowVarRefState, only: nw, nwf, nt1, nt2
        use inputPhysics, only: turbModel
        use paramTurb, only: rsaGRgammaLo, rsaGRgammaHi, rsaGRreThetaLo
        use inputtimespectral, only: nTimeIntervalsSpectral
        use utils, only: setPointers, EChk
        use genericISNAN, only: myisnan
        use communication, only: ADflow_comm_world
        implicit none

        ! input variable
        real(kind=realType), intent(inout) :: lambdaP

        ! local variables
        integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii
        real(kind=realType), pointer :: wvec_pointer(:)
        real(kind=realType), pointer :: dvec_pointer(:)
        real(kind=alwaysRealType) :: lambdaL ! L is for local
        real(kind=alwaysRealType) :: lambdaP_recv ! to receive the global step
        real(kind=alwaysRealType) :: ratio
        real(kind=alwaysRealType) :: wval, dval, ratioBound, gammaFull
        real(kind=realType) :: cs(nt2 - nt1 + 1), gLoS, gHiS, rLoS

        ! Determine the maximum step size that would yield
        ! a maximum change of 10% in density, total energy,
        ! and turbulence variable after a KSP solve.

        ! wVecTurb/deltaWTurb hold the COLUMN-SCALED state/update
        ! (see getTurbColScale). Relative ratios (w/d) are scale-invariant,
        ! but the absolute gamma/retheta bounds must be scaled to match.
        if (turbModel == spalartallmarasnoft2gammaretheta) then
            call getTurbColScale(cs, nt2 - nt1 + 1)
            gLoS = rsaGRgammaLo * cs(2)
            gHiS = rsaGRgammaHi * cs(2)
            rLoS = rsaGRreThetaLo * cs(3)
        end if

        ! Initialize the local step size as ANK_stepFactor
        ! because the initial step is likely to be equal to this.
        lambdaL = real(lambdaP)

        ! First we need to read both the update and the state
        ! from PETSc because the w in ADFlow currently contains
        ! the state that is perturbed during the matrix-free
        ! operations.

        ! wVec contains the state vector
        call VecGetArrayF90(wVecTurb, wvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! deltaW contains the full update
        call VecGetArrayF90(deltaWTurb, dvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ii = 1
        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                call setPointers(nn, 1_intType, sps)
                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il
                            ! multiply the ratios by 10 to check if the change in a
                            ! variable is greater than 10% of the variable itself.

                            ! needs to be modified
                            ! if coupled ank is used, nstate = nw and this loop is executed
                            ! if no turbulence variables, this loop will be automatically skipped
                            ! check turbulence variable
#ifndef USE_COMPLEX
                            ratio = (wvec_pointer(ii) / (dvec_pointer(ii) + eps)) * ANK_physLSTolTurb
#else
                            ratio = (real(wvec_pointer(ii)) / real(dvec_pointer(ii) + eps)) * real(ANK_physLSTolTurb)
#endif
                            ! if the ratio is less than min step, the update is either
                            ! in the positive direction, therefore we do not clip it,
                            ! or the update is very limiting, so we just clip the
                            ! individual update for this cell.
                            if (ratio .lt. ANK_stepFactor * ANK_stepMin) then
                                ! The update was very limiting, so just clip this
                                ! individual update and dont change the overall
                                ! step size. To select the new update, instead of
                                ! clipping to zero, we clip to 1 percent of the original.
                                if (ratio .gt. zero) &
                                    dvec_pointer(ii) = wvec_pointer(ii) * ANK_physLSTolTurb

                                ! Either case, set the ratio to one. Positive updates
                                ! do not limit the step, negative updates below minimum
                                ! step were already clipped.
                                ratio = one
                            end if
                            lambdaL = min(lambdaL, ratio)
                            ii = ii + 1

                            ! Physicality checks for gamma and ReTheta
                            do l = nt1 + 1, nt2
#ifndef USE_COMPLEX
                                wval = wvec_pointer(ii)
                                dval = dvec_pointer(ii)
#else
                                wval = real(wvec_pointer(ii))
                                dval = real(dvec_pointer(ii))
#endif

                                if (turbModel == spalartallmarasnoft2gammaretheta) then
                                    if (l == nt1 + 1) then
                                        ! GAMMA: absolute bound enforcement, no relative check.
                                        ! Full step allowed if result stays in [gammaLo, gammaHi].
                                        ratio = one
                                        gammaFull = wval - dval

                                        if (gammaFull > gHiS) then
                                            ratio = (wval - gHiS) / (dval + eps)
                                        else if (gammaFull < gLoS) then
                                            ratio = (wval - gLoS) / (dval + eps)
                                        end if

                                        ratio = max(ratio, zero)

                                        if (ratio < omegaMinGamma) then
                                            ! Clip individual update to stay in bounds
                                            if (dval > zero) then
                                                dvec_pointer(ii) = wval - gLoS
                                            else if (dval < zero) then
                                                dvec_pointer(ii) = wval - gHiS
                                            end if
                                            ratio = one
                                        end if

                                    else if (l == nt1 + 2) then
                                        ! RETHETA: relative check with own tolerance + lower bound.
#ifndef USE_COMPLEX
                                        ratio = (wval / (dval + eps)) * ANK_physLSTolReTheta
#else
                                        ratio = (wval / real(dval + eps)) * real(ANK_physLSTolReTheta)
#endif
                                        ! Lower bound enforcement
                                        if (dval > zero) then
                                            ratioBound = (wval - rLoS) / (dval + eps)
                                            ratioBound = max(ratioBound, zero)
                                            ratio = min(ratio, ratioBound)
                                        end if

                                        if (ratio < ANK_stepFactor * ANK_stepMin) then
                                            if (ratio > zero) &
                                                dvec_pointer(ii) = wval * ANK_physLSTolReTheta
                                            ratio = one
                                        end if
                                    end if
                                else
                                    ! Non-transition: use standard relative check
#ifndef USE_COMPLEX
                                    ratio = (wval / (dval + eps)) * ANK_physLSTolTurb
#else
                                    ratio = (wval / real(dval + eps)) * real(ANK_physLSTolTurb)
#endif
                                    if (ratio < ANK_stepFactor * ANK_stepMin) then
                                        if (ratio > zero) &
                                            dvec_pointer(ii) = wval * ANK_physLSTolTurb
                                        ratio = one
                                    end if
                                end if

                                lambdaL = min(lambdaL, ratio)
                                ii = ii + 1
                            end do
                        end do
                    end do
                end do
            end do
        end do

        ! Restore the pointers to PETSc vectors

        call VecRestoreArrayF90(wVecTurb, wvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecRestoreArrayF90(deltaWTurb, dvec_pointer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Make sure that we did not get any NaN's in the process
        if (myisnan(lambdaL)) then
            lambdaL = zero
        end if

        ! Finally, communicate the step size across processes and return
        ! mpi allreduce is not defined for complex numbers with the min operation
        ! so we will use the lambdaP_recv variable to receive
        call mpi_allreduce(lambdaL, lambdaP_recv, 1_intType, MPI_DOUBLE, &
                           mpi_min, ADflow_comm_world, ierr)
        call EChk(ierr, __FILE__, __LINE__)

#ifndef USE_COMPLEX
        lambdaP = lambdaP_recv
#else
        ! finally, as a safety check, purge the complex part of lambda
        lambdaP = cmplx(lambdaP_recv, 0.0_realType)
#endif

    end subroutine physicalityCheckANKTurb

    subroutine ANKTurbSolveKSP

        ! This routine solves the turbulence model equation using
        ! a similar approach to the main ank solver.

        use constants
        use blockPointers, only: nDom, flowDoms
        use inputIteration, only: L2conv, transitionSrcDtRestrict, noBacktrackCount, srcDtDeactivateIters, transitionNK
        use paramTurb, only: srcLambdaModeFull
        use saGammaReTheta, only: computeSrcLambda
        use inputTimeSpectral, only: nTimeIntervalsSpectral
        use inputDiscretization, only: approxSA, orderturb
        use iteration, only: approxTotalIts, totalR0, totalR, currentLevel
        use utils, only: EChk, setPointers
        use genericISNAN, only: myisnan
        use solverUtils, only: computeUTau
        use NKSolver, only: getEwTol
        use BCRoutines, only: applyAllBC, applyAllBC_block
        use haloExchange, only: whalo2
        use oversetData, only: oversetPresent
        use flowVarRefState, only: nw, nwf, nt1, nt2, kPresent, pInfCorr
        use communication
        use blockette, only: blocketteRes
        implicit none

        ! Working Variables
        integer(kind=intType) :: ierr, maxIt, kspIterations, nn, sps, reason, nHist, iter, feval, orderturbsave
        integer(kind=intType) :: i, j, k, n
        real(kind=realType) :: atol, val, v2, factK, gm1
        real(kind=alwaysRealType) :: rtol, totalR_dummy, linearRes, norm
        real(kind=alwaysRealType) :: resHist(ANK_maxIter + 1)
        real(kind=alwaysRealType) :: unsteadyNorm, unsteadyNorm_old
        real(kind=alwaysRealType) :: linResMonitorTurb, totalRTurb
        logical :: correctForK, LSFailed, srcDtRestrictActive, backtrackTriggered

        ! Derived condition for source dt restriction (ANK turbKSP path)
        ! DD-ADI uses transitionSrcDtRestrict alone; ANK adds deactivation after clean iters.
        srcDtRestrictActive = transitionNK .and. transitionSrcDtRestrict .and. (noBacktrackCount < srcDtDeactivateIters)

        ! Calculate the residuals and set rVecTurb before the first iteration
        call blocketteRes(useFlowRes=.False., useStoreWall=.False.)
        call setRVecANKTurb(rVecTurb)

        ! Freeze srcLambda at base state before ANK iterations
        if (transitionSrcDtRestrict .and. srcDtRestrictActive) then
            call computeSrcLambda(srcLambdaModeFull)
        end if

        do n = 1, ANK_nsubIterTurb

            ! Compute the norm of rVecTurb, which is identical to the
            ! norm of the unsteady residual vector.
            call VecNorm(rVecTurb, NORM_2, totalRTurb, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Determine if we need to form the Preconditioner
            if (mod(ANK_iterTurb, ANK_jacobianLag) == 0) then

                ! Actually form the preconditioner and factorize it.
                if (myid .eq. 0 .and. ANK_turbDebug) &
                    write (*, *) "Re-doing turb PC"
                call FormJacobianANKTurb()
                ANK_iterTurb = 0
            end if

            ! Increment the iteration counter
            ANK_iterTurb = ANK_iterTurb + 1

            ! Start with trying to take the full step set by the user.
#ifndef USE_COMPLEX
            lambdaTurb = ANK_StepFactor
#else
            ! make sure we zero out the complex part of the step size
            lambdaTurb = cmplx(real(ANK_StepFactor), 0.0_realType)
#endif

            ! Dummy matrix assembly for the matrix-free matrix
            call MatAssemblyBegin(dRdwTurb, MAT_FINAL_ASSEMBLY, ierr)
            call EChk(ierr, __FILE__, __LINE__)
            call MatAssemblyEnd(dRdwTurb, MAT_FINAL_ASSEMBLY, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            if (totalR > ANK_secondOrdSwitchTol * totalR0) then
                ! Save if second order turbulence is used, we will only use 1st order during ANK (only matters for the coupled solver)
                approxSA = .True.
                orderturbsave = orderturb
                orderturb = firstOrder

                ! Determine the relative convergence for the KSP solver
                rtol = ANK_rtol ! Just use the input relative tolerance for approximate fluxes
            else
                ! If the second order fluxes are used, Eisenstat-Walker algorithm to determine relateive
                ! convergence tolerance helps with performance.
                totalR_dummy = totalR
                call getEWTol(totalR_dummy, totalR_old, rtolLast, rtol)

                ! Use the ANK rtol if E-W algorithm is not picking anything lower
                rtol = min(ANK_rtol, rtol)
            end if

            ! also check if we are using approxSA always
            if (ANK_useApproxSA) &
                approxSA = .True.

            ! Record the total residual and relative convergence for next iteration
            totalR_old = totalR
            rtolLast = rtol

            ! Set all tolerances for linear solve:
#ifndef USE_COMPLEX
            ! in the real mode, we set the atol slightly lower than the target L2 convergence
            ! the reasoning for this is detailed in the NKStep subroutine
            atol = totalR0 * L2Conv * 0.01_realType
#else
            ! in complex mode, we want to tightly solve the linear system every time
            ! again, see the NKStep subroutine for the explanation
            atol = totalR0 * L2Conv * 1e-6_realType
#endif

            ! Set the iteration limit to maxIt, determined by which fluxes are used.
            ! This is because ANK step require 0.1 convergence for stability during initial stages.
            ! Due to an outdated preconditioner, the KSP solve might take more iterations.
            ! If this happens, the preconditioner is re-computed and because of this,
            ! ANK iterations usually don't take more than 2 times number of ANK_subSpace size iterations
            call KSPSetTolerances(ANK_KSPTurb, rtol, &
                                  real(atol), real(ANK_divTol), ank_maxIter, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call KSPSetResidualHistory(ANK_KSPTurb, resHist, ank_maxIter + 1, PETSC_TRUE, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Set the BaseVector of the matrix-free matrix:
            call formFunction_mf_turb(ctx, wVecTurb, baseResTurb, ierr)
            call EChk(ierr, __FILE__, __LINE__)
            call MatMFFDSetBase(dRdWTurb, wVecTurb, baseResTurb, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Refresh srcLambda from updated base state before KSPSolve
            if (transitionSrcDtRestrict .and. srcDtRestrictActive) then
                call computeSrcLambda(srcLambdaModeFull)
            end if

            ! Actually do the Linear Krylov Solve
            call KSPSolve(ANK_KSPTurb, rVecTurb, deltaWTurb, ierr)

            ! DON'T just check the error. We want to catch error code 72
            ! which is a floating point error. This is ok, we just reset and
            ! keep going
            if (ierr == 72) then
                ! The convergence check will get the nan
            else
                call EChk(ierr, __FILE__, __LINE__)
            end if

            ! Get the number of iterations from the KSP solver
            call KSPGetIterationNumber(ANK_KSPTurb, kspIterations, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call KSPGetConvergedReason(ANK_KSPTurb, reason, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Return previously changed variables back to normal, VERY IMPORTANT
            if (totalR > ANK_secondOrdSwitchTol * totalR0) then
                ! Replace the second order turbulence option
                orderturb = orderturbsave
                approxSA = .False.
            end if

            ! put back the approxsa flag if we were using it
            if (ANK_useApproxSA) &
                approxSA = .False.

            ! Compute the maximum step that will limit the change
            ! in SA variable to some user defined fraction.
            call physicalityCheckANKTurb(lambdaTurb)
            !if (myid .eq. 0) write(*,*)"physicality check lambda: ",lambdaTurb
            !lambdaTurb = max(ANK_stepMin, lambdaTurb)

            ! Take the uodate after the physicality check.
            call VecAXPY(wVecTurb, -lambdaTurb, deltaWTurb, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Set the updated state variables
            call setWANKTurbScaled(wVecTurb)

            ! Compute the unsteady residuals. The actual residuals
            ! also get calculated in the process, and are stored in
            ! dw. Make sure to call setRVec/setRVecANK after this
            ! routine because rVec contains the unsteady residuals,
            ! and we need the steady residuals for the next iteration.
            call computeUnsteadyResANKTurb(lambdaTurb)

            ! Count the number of of residual evaluations outside the KSP solve
            feval = 1_intType

            ! Check if the norm of the rVec is bad:
            call VecNorm(rVecTurb, NORM_2, unsteadyNorm, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! initialize this outside the ls
            LSFailed = .False.
            backtrackTriggered = .False.

            if ((unsteadyNorm > totalRTurb * ANK_unstdyLSTol .or. myisnan(unsteadyNorm))) then
                ! The unsteady residual is too high or we have a NAN. Do a
                ! backtracking line search until we get a residual that is lower.

                LSFailed = .True.
                backtrackTriggered = .True.

                ! Restore the starting (old) w value by adding lamda*deltaW
                call VecAXPY(wVecTurb, lambdaTurb, deltaWTurb, ierr)
                call EChk(ierr, __FILE__, __LINE__)

                ! Set the initial new lambda. This is working off the
                ! potentially already physically limited step.
                lambdaTurb = 0.7_realType * lambdaTurb

                backtrack: do iter = 1, 12

                    ! Apply the new step
                    call VecAXPY(wVecTurb, -lambdaTurb, deltaWTurb, ierr)
                    call EChk(ierr, __FILE__, __LINE__)

                    ! Set and recompute
                    call setWANKTurbScaled(wVecTurb)

                    ! Compute the unsteady residuals with the current step
                    call computeUnsteadyResANKTurb(lambdaTurb)
                    feval = feval + 1

                    call VecNorm(rVecTurb, NORM_2, unsteadyNorm, ierr)
                    call EChk(ierr, __FILE__, __LINE__)

                    if (unsteadyNorm > totalRTurb * ANK_unstdyLSTol .or. myisnan(unsteadyNorm)) then

                        ! Restore back to the original wVec
                        call VecAXPY(wVecTurb, lambdaTurb, deltaWTurb, ierr)
                        call EChk(ierr, __FILE__, __LINE__)

                        ! Haven't backed off enough yet....keep going
                        lambdaTurb = lambdaTurb * 0.7_realType
                    else
                        ! We have succefssfully reduced the norm
                        LSFailed = .False.
                        exit
                    end if
                end do backtrack

                if (LSFailed .or. myisnan(unsteadyNorm)) then
                    ! the line search wasn't much help.

                    if (ANK_CFL > ANK_CFLMin) then
                        ! the cfl number is not already at the lower limit.  We
                        ! can cut the CFL back and try again. Set lambda to zero
                        ! to indicate we never took a step.
                        lambdaTurb = zero
                    else
                        ! cfl is as low as it goes, try taking the step
                        ! anyway. We can't do  anything else
                        call VecAXPY(wVecTurb, -lambdaTurb, deltaWTurb, ierr)
                        call EChk(ierr, __FILE__, __LINE__)
                    end if

                    ! Set the state vec and compute the new residual
                    call setWANKTurbScaled(wVecTurb)
                    call blocketteRes(useFlowRes=.False., &
                                      useStoreWall=.False.)
                    feval = feval + 1
                end if
            end if

            call setRvecANKTurb(rVecTurb)

            linResMonitorTurb = resHist(kspIterations + 1) / resHist(1)

            if ((linResMonitorTurb .ge. ANK_rtol .and. &
                 totalR > ANK_secondOrdSwitchTol * totalR0 .and. &
                 linResOldTurb .le. ANK_rtol) &
                !.or. LSFailed) then
                !            .or. lambdaTurb .le. ANK_stepMin) then
                .or. (lambdaTurb .eq. zero)) then

                ! We should reform the PC since it took longer than we want,
                ! or we need to adjust the CFL because the last update was bad,
                ! or convergence since the last PC update was good enough and we
                ! would benefit from re-calculating the PC.
                ANK_iterTurb = 0
            end if

            ! update the linear residual for next iteration
            linResOldTurb = linResMonitorTurb

            ! Update step monitor
            ! stepMonitor = lambda

            ! Update the approximate iteration counter. The +1 is for the
            ! residual evaluations.
            ! approxTotalIts = approxTotalIts + feval + kspIterations

            ! Print some info about the turbulence ksp
            if (myid == 0 .and. ANK_turbDebug) then
                Write (*, *) "LIN RES, ITER, INITRES, REASON, STEP", linResMonitorTurb, kspIterations, &
                    reshist(1), reason, lambdaTurb
            end if

            ! ============== Source-dt deactivation switch (P&Z §IV.B.3) ==============
            ! Clean iterations count only in the inexact-Newton analog phase, i.e.
            ! the second-order regime (totalR <= ANK_secondOrdSwitchTol * totalR0).
            ! After srcDtDeactivateIters clean iterations, deactivate source dt
            ! restriction. Reset the counter (reactivate) when backtracking is
            ! triggered — even if the backtrack succeeds — or when the relative
            ! residual rises back above the phase-switch tolerance.
            if (transitionNK .and. transitionSrcDtRestrict) then
                if (backtrackTriggered .or. totalR > ANK_secondOrdSwitchTol * totalR0) then
                    noBacktrackCount = 0
                else
                    noBacktrackCount = noBacktrackCount + 1
                end if
                srcDtRestrictActive = transitionNK .and. transitionSrcDtRestrict .and. (noBacktrackCount < srcDtDeactivateIters)
            end if

        end do

    end subroutine ANKTurbSolveKSP

    subroutine ANKStep(firstCall)

        use constants
        use blockPointers, only: nDom, flowDoms, shockSensor, ib, jb, kb, p, w, gamma
        use inputPhysics, only: equations, turbModel
        use inputIteration, only: L2conv, transitionSrcDtRestrict, noBacktrackCount, srcDtDeactivateIters, transitionNK, &
                                  solverStallDiag, solverStallDiagStep, ankCFLMinCap, &
                                  ankUnsteadyLSFactor, ankUnsteadyLSMaxIter, ankRejectOnLSExhausted, &
                                  ankAlgorithm2Damping, ankTransitionGlobalLambda, ankColScaleUnit
        use paramTurb, only: srcLambdaModeFull
        use saGammaReTheta, only: computeSrcLambda
        use inputTimeSpectral, only: nTimeIntervalsSpectral
        use inputDiscretization, only: lumpedDiss, approxSA, orderturb
        use iteration, only: approxTotalIts, totalR0, totalR, stepMonitor, linResMonitor, currentLevel, iterType
        use utils, only: EChk, setPointers
        use genericISNAN, only: myisnan
        use turbAPI, only: turbSolveDDADI
        use solverUtils, only: computeUTau
        use adjointUtils, only: referenceShockSensor
        use NKSolver, only: setRVec, getEwTol
        use initializeFlow, only: setUniformFlow
        use BCRoutines, only: applyAllBC, applyAllBC_block
        use haloExchange, only: whalo2
        use oversetData, only: oversetPresent
        use flowVarRefState, only: nw, nwf, nt1, nt2, kPresent, pInfCorr
        use flowUtils, only: computeLamViscosity
        use turbUtils, only: computeEddyViscosity
        use communication
        use blockette, only: blocketteRes
        implicit none

        ! Input Variables
        logical, intent(in) :: firstCall

        ! Working Variables
        integer(kind=intType) :: ierr, maxIt, kspIterations, nn, sps, reason, nHist, iter, feval, orderturbsave
        integer(kind=intType) :: i, j, k
        real(kind=realType) :: atol, val, v2, factK, gm1
        real(kind=alwaysRealType) :: rtol, totalR_dummy, linearRes, norm
        real(kind=alwaysRealType) :: resHist(ANK_maxIter + 1)
        real(kind=alwaysRealType) :: unsteadyNorm, unsteadyNorm_old, rel_pcUpdateTol
        real(kind=alwaysRealType) :: pzRd, lsTolEff
        real(kind=realType) :: lsFactorEff
        integer(kind=intType) :: lsMaxIterEff
        logical :: correctForK, LSFailed, backtrackTriggeredANK

        ! Enter this check if this is the first ANK step OR we are switching to the coupled ANK solver
        if (firstCall .or. &
            ((totalR .le. ANK_coupledSwitchTol * totalR0) .and. (.not. ANK_coupled) &
             .and. (equations .eq. RANSEquations))) then

            ! If this is a first call, we need to change the coupled switch
            ! to the correct value.
            if (firstCall) then

                ! Check if we are above or below the coupled switch tolerance.
                ! PZ mode: the thesis's main solver is fully coupled from
                ! iteration 1 (set ANKCoupledSwitchTol >= 1 for that), but
                ! Appendix B also runs decoupled/loosely-coupled variants in
                ! the approximate-Newton phase (fully coupled in phase 2), so
                ! PZ mode respects ANKCoupledSwitchTol like standard ADflow.
                if (totalR .le. ANK_coupledSwitchTol * totalR0 &
                    .and. equations .eq. RANSEquations) then
                    ANK_coupled = .True.
                else
                    ANK_coupled = .False.
                end if

                ! This is not a first call, and the only option left is that,
                ! we may be switching from uncoupled to coupled
            else
                ANK_coupled = .True.
            end if

            ! If we are in here, destroy the solver regardless,
            ! and set up with the correct coupling mode.
            call destroyANKSolver()
            call setupANKSolver()

            ! Copy the adflow 'w' into the petsc wVec
            call setWVecANKScaled(wVec, 1, nstate)

            ! Evaluate the residual before we start
            call blocketteRes(useUpdateIntermed=.True.)
            if (ANK_coupled) then
                call setRvec(rVec)
            else
                call setRVecANK(rVec)
            end if

            ! Check if we are using the turb KSP
            if ((.not. ANK_coupled) .and. (.not. ANK_useTurbDADI) .and. equations == RANSEquations) then
                call setWVecANKTurbScaled(wVecTurb)
                call setRVecANKTurb(rVecTurb)
            end if

            if (firstCall) then
                ! Start with the selected fraction of the ANK_StepFactor
                lambda = ANK_StepFactor
                lambdaTurb = ANK_stepFactor

                ! Initialize some variables
                totalR_old = totalR ! Record the old residual for the first iteration
                rtolLast = ANK_rtol ! Set the previous relative convergence tolerance for the first iteration
                if (ANK_CFLReset) then
                    ! we are asked to reset the cfl to the initial value on every first ANK iteration.
                    ! if this is set to false, the CFL must be set from the python layer using either the
                    ! initial CFL or the last CFL for the aeroproblem
                    ANK_CFL = ANK_CFL0
                end if
                ANK_CFLMinBase = ANK_CFLMin0
                totalR_pcUpdate = totalR ! only update the residual at last PC calculation for the first iteration
                linResOld = zero
                linResOldTurb = zero
                ANK_iter = 0

                ! PZ mode initialization: the whole P&Z stabilization stack is
                ! one switch — no hybrid with the ADflow controller. Force the
                ! pieces of their system that already exist on this branch ON,
                ! and start the reference CFL at its phase-1 seed.
                if (ANK_pzStepping) then
                    transitionSrcDtRestrict = .True.  ! thesis Eq. 3.14 (source-term dt restriction)
                    ankAlgorithm2Damping = .True.     ! thesis Alg. 3 (per-node gamma/ReTheta damping)
                    pz_CFLRef = ANK_pzCFL0
                    pz_alphaSER = -one
                    pz_prevAccepted = .False.
                    if (myid == 0) then
                        print *, 'PZ stepping mode ON: CFL = ', real(ANK_pzCFL0), &
                            ' * ', real(ANK_pzGrowth), '^n (phase 1), SER beta = ', &
                            real(ANK_pzBeta), ' (phase 2); ADflow CFL controller disabled.'
                    end if
                end if
            end if
        else
            ANK_iter = ANK_iter + 1
        end if

        ! figure out if we want to change the ANKPCUpdateTol:
        ! If we are converged below the L2 target specified by ANK_pcUpdateCutoff
        ! we use a different tolerance for the relative PC update.
        if (totalR / totalR0 .lt. ANK_pcUpdateCutoff) then
            rel_pcUpdateTol = ANK_pcUpdateTol2
        else
            ! if we are not that far down converged, use the option directly
            rel_pcUpdateTol = ANK_pcUpdateTol
        end if

        ! Compute the norm of rVec, which is identical to the
        ! norm of the unsteady residual vector.
        call VecNorm(rVec, NORM_2, unsteadyNorm_old, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! For the coupled SA-Gamma-Retheta solve, freeze the source-term
        ! eigenvalues at the base state before the time-step matrix is
        ! formed (P&Z Eq. 59; both computeTimeStepMat branches below
        ! consume srcLambda through computeTimeStepBlock).
        if (ANK_coupled .and. turbModel == spalartallmarasnoft2gammaretheta .and. transitionNK .and. &
            transitionSrcDtRestrict .and. (noBacktrackCount < srcDtDeactivateIters)) then
            call computeSrcLambda(srcLambdaModeFull)
        end if

        ! ============== PZ reference-CFL law (replaces the ADflow controller) ==============
        ! Updated EVERY outer iteration (the Diablo dt_ref is), not only at PC
        ! reforms. Growth is only applied after an ACCEPTED step; a rejected
        ! step instead halves pz_CFLRef at the bottom of this routine.
        if (ANK_pzStepping) then
            if (totalR > ANK_secondOrdSwitchTol * totalR0) then
                ! Phase 1 (approximate-Newton): geometric ramp, CFL = CFL0 * b^n
                if (pz_prevAccepted) pz_CFLRef = min(pz_CFLRef * ANK_pzGrowth, ANK_pzCFLMax)
                ! Not in phase 2 (yet, or anymore): drop the SER anchor so it is
                ! re-anchored (continuously) on the next phase-2 entry.
                pz_alphaSER = -one
            else
                ! Phase 2 (inexact-Newton analog): Mulder & van Leer SER,
                ! CFL^(n) = max(alpha * Rd^-beta, CFL^(n-1)), monotone while accepted.
                pzRd = totalR / totalR0
                if (pz_alphaSER < zero) pz_alphaSER = pz_CFLRef * pzRd**ANK_pzBeta
                if (pz_prevAccepted) then
                    pz_CFLRef = min(max(pz_alphaSER * pzRd**(-ANK_pzBeta), pz_CFLRef), ANK_pzCFLMax)
                end if
            end if
            ANK_CFL = pz_CFLRef
            ! dt_ref,min analogue: the single floor the rejection logic tests.
            ANK_CFLMin = ANK_pzCFLMin
        end if

        ! Determine if if we need to form the Preconditioner
        if (mod(ANK_iter, ANK_jacobianLag) == 0 .or. totalR / totalR_pcUpdate < rel_pcUpdateTol) then

            if (.not. ANK_pzStepping) then
            ! First of all, update the minimum cfl wrt the overall convergence
            ANK_CFLMin = min(ANK_CFLLimit, ANK_CFLMinBase * (totalR0 / totalR)**ANK_CFLExponent)

            ! VERIF_06 F1: cap how far that ramp may carry the floor. Without a
            ! cap the floor reaches ANK_CFLLimit at depth and the cutback below
            ! becomes a no-op, so the controller can never back off a bad step.
            ! This is the analogue of the thesis's user-specified dt_ref,min
            ! (§3.1.3, Algorithms 2 and 4). <= 0 => disabled (previous behaviour).
            if (ankCFLMinCap > zero) ANK_CFLMin = min(ANK_CFLMin, ankCFLMinCap)

            ! Update the CFL number depending on the outcome of the last iteration
            stallCFLBefore = ANK_CFL
            stallCFLFloored = .False.
            if (lambda < ANK_stepMin * ANK_stepFactor) then

                ! The step was too small, cut back the cfl
                ! VERIF_06 F1: ANK_CFLMin is itself ramped by convergence
                ! (ANK_CFLMinBase*(totalR0/totalR)**ANK_CFLExponent), so at depth
                ! it can reach ANK_CFLLimit and this max() blocks the cutback
                ! entirely -- the controller loses its only recovery mechanism.
                ! Record when that happens so a stalling log says so.
                if (solverStallDiag .and. ANK_CFL * ANK_CFLCutback < ANK_CFLMin) &
                    stallCFLFloored = .True.
                ANK_CFL = max(ANK_CFL * ANK_CFLCutback, ANK_CFLMin)

            else if (totalR < totalR_pcUpdate .and. lambda .ge. ANK_constCFLStep * ANK_stepFactor) then

                ! total residuals have decreased since the last cfl
                ! change, or the step was large enough, we can ramp
                ! the cfl up
                ANK_CFL = max(min(ANK_CFL * ANK_CFLFactor** &
                                  ((totalR_pcUpdate - totalR) / totalR_pcUpdate), ANK_CFLLimit), ANK_CFLMin)

            else

                ! The step was not small, but it was not large enough
                ! to ramp up the cfl, so we keep it constant. Just
                ! make sure that the cfl does not go below the
                ! minimum value.
                ANK_CFL = max(ANK_CFL, ANK_CFLMin)

            end if
            end if ! .not. ANK_pzStepping (PZ mode: CFL set by the PZ law above)

            ! F9 metric: report what the column scaling actually achieves,
            ! at PC-reform cadence so the cost is negligible.
            if ((solverStallDiag .or. ankColScaleUnit) .and. ANK_coupled) &
                call reportScaledStateStats(wVec, deltaW)

            ! Record the total residuals when the PC is calculated.
            totalR_pcUpdate = totalR

            ! Update the time step terms before forming the PC because the states and CFL may have changed
            ! This call also updates the time step terms for the AMG PC if required
            call computeTimeStepMat(usePC=.True.)

            ! Actually form the preconditioner and factorize it.
            call FormJacobianANK()

            if (totalR .le. ANK_secondOrdSwitchTol * totalR0) then
                if (ANK_coupled) then
                    iterType = "  *CSANK"
                else
                    iterType = "   *SANK"
                end if
            else
                if (ANK_coupled) then
                    iterType = "   *CANK"
                else
                    iterType = "    *ANK"
                end if
            end if
            ANK_iter = 0

            ! Also update the turb PC bec. the CFL has changed
            ANK_iterTurb = 0
        else
            ! Update the time step terms because the states may have changed from the last step
            ! This call does not update the time step terms for the AMG PC
            call computeTimeStepMat(usePC=.False.)

            if (totalR .le. ANK_secondOrdSwitchTol * totalR0) then
                if (ANK_coupled) then
                    iterType = "   CSANK"
                else
                    iterType = "    SANK"
                end if
            else
                if (ANK_coupled) then
                    iterType = "    CANK"
                else
                    iterType = "     ANK"
                end if
            end if
        end if

        ! Start with trying to take the full step set by the user.
#ifndef USE_COMPLEX
        lambda = ANK_StepFactor
#else
        ! make sure we zero out the complex part of the step size
        lambda = cmplx(real(ANK_StepFactor), 0.0_realType)
#endif

        ! Dummy matrix assembly for the matrix-free matrix
        call MatAssemblyBegin(dRdw, MAT_FINAL_ASSEMBLY, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        call MatAssemblyEnd(dRdw, MAT_FINAL_ASSEMBLY, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! ============== Flow Update =============

        ! For the approximate solver, we need the approximate flux routines
        ! We set the variables required for approximate fluxes here and they will be used
        ! for the matrix-free matrix-vector product routines when the KSP solver calls it
        ! Very important to set the variables back to their original values after each
        ! KSP solve because we want actual flux functions when calculating residuals
        if (totalR > ANK_secondOrdSwitchTol * totalR0) then
            ! Setting lumped dissipation to true gives approximate fluxes
            ANK_useDissApprox = .True.
            lumpedDiss = .True.
            approxSA = .True.

            ! Save the turbulence order, we will only use 1st order during ANK (only matters for the coupled solver)
            orderturbsave = orderturb
            orderturb = firstOrder

            ! Calculate the shock sensor here because the approximate routines do not
            call referenceShockSensor()

            ! Determine the relative convergence for the KSP solver
            rtol = ANK_rtol ! Just use the input relative tolerance for approximate fluxes

        else
            ! If the second order fluxes are used, Eisenstat-Walker algorithm to determine relateive
            ! convergence tolerance helps with performance.
            totalR_dummy = totalR
            call getEWTol(totalR_dummy, totalR_old, rtolLast, rtol)

            ! Use the ANK rtol if E-W algorithm is not picking anything lower
            rtol = min(ANK_rtol, rtol)
        end if

        ! also check if we are using approxSA always
        if (ANK_useApproxSA) &
            approxSA = .True.

        ! Record the total residual and relative convergence for next iteration
        totalR_old = totalR
        rtolLast = rtol

        ! Set all tolerances for linear solve:
#ifndef USE_COMPLEX
        ! in the real mode, we set the atol slightly lower than the target L2 convergence
        ! the reasoning for this is detailed in the NKStep subroutine
        atol = totalR0 * L2Conv * ANK_atol_buffer
#else
        ! in complex mode, we want to tightly solve the linear system every time
        ! again, see the NKStep subroutine for the explanation
        atol = totalR0 * L2Conv * 1e-6_realType
#endif

        ! Set the iteration limit to maxIt, determined by which fluxes are used.
        ! This is because ANK step require 0.1 convergence for stability during initial stages.
        ! Due to an outdated preconditioner, the KSP solve might take more iterations.
        ! If this happens, the preconditioner is re-computed and because of this,
        ! ANK iterations usually don't take more than 2 times number of ANK_subSpace size iterations
        call KSPSetTolerances(ANK_KSP, rtol, &
                              real(atol), real(ANK_divTol), ank_maxIter, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call KSPSetResidualHistory(ANK_KSP, resHist, ank_maxIter + 1, PETSC_TRUE, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Set the BaseVector of the matrix-free matrix:
        call formFunction_mf(ctx, wVec, baseRes, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        call MatMFFDSetBase(dRdW, wVec, baseRes, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Actually do the Linear Krylov Solve
        call KSPSolve(ANK_KSP, rVec, deltaW, ierr)

        ! DON'T just check the error. We want to catch error code 72
        ! which is a floating point error. This is ok, we just reset and
        ! keep going
        if (ierr == 72) then
            ! The convergence check will get the nan
        else
            call EChk(ierr, __FILE__, __LINE__)
        end if

        ! Get the number of iterations from the KSP solver
        call KSPGetIterationNumber(ANK_KSP, kspIterations, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call KSPGetConvergedReason(ANK_KSP, reason, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Return previously changed variables back to normal, VERY IMPORTANT
        if (totalR > ANK_secondOrdSwitchTol * totalR0) then
            ! Set ANK_useDissApprox back to False to go back to using actual flux routines
            ANK_useDissApprox = .False.
            lumpedDiss = .False.
            approxSA = .False.

            ! Replace turbulence order
            orderturb = orderturbsave

        end if

        ! put back the approxsa flag if we were using it
        if (ANK_useApproxSA) &
            approxSA = .False.

        ! Compute the maximum step that will limit the change in pressure
        ! and energy to some user defined fraction.
        if (ANK_pzStepping) then
            ! Thesis Algorithm 2 (physics-based restriction): global
            ! delta_phys = min over cells with a NEGATIVE rho/E update of
            ! 0.90*Q/|dQ|. delta < 1% => reject the update and retry with the
            ! reference CFL halved (done at the bottom of this routine); at the
            ! CFL floor, clip to 1% and continue.
            call pzPhysicsRestriction(lambda)
            if (lambda < 0.01_realType) then
                if (ANK_CFL > ANK_CFLMin) then
                    lambda = zero ! reject; pz_CFLRef is halved below
                else
                    lambda = 0.01_realType
                end if
            end if
        else
            call physicalityCheckANK(lambda)
            if (ANK_CFL .gt. ANK_CFLMin .and. lambda .lt. ANK_stepMin) &
                lambda = zero
        end if

        ! Take the uodate after the physicality check.
        call VecAXPY(wVec, -lambda, deltaW, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Set the updated state variables
        call setWANKScaled(wVec, 1, nState)

        ! Compute the unsteady residuals. The actual residuals
        ! also get calculated in the process, and are stored in
        ! dw. Make sure to call setRVec/setRVecANK after this
        ! routine because rVec contains the unsteady residuals,
        ! and we need the steady residuals for the next iteration.
        call computeUnsteadyResANK(lambda)

        ! Count the number of of residual evaluations outside the KSP solve
        feval = 1_intType

        ! Check if the norm of the rVec is bad:
        call VecNorm(rVec, NORM_2, unsteadyNorm, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! initialize this outside the ls
        LSFailed = .False.
        backtrackTriggeredANK = .False.
        ! Stall diagnostics: lambda entering the line search is the
        ! physicality-limited one; reset the per-iteration LS counters.
        stallNBacktrack = 0
        stallLSFailed = .False.

        ! Effective line-search parameters. PZ mode = thesis Algorithm 4:
        ! accept when R_unst < ||R(Q^n)||_2 (tol 1), back off by 0.90, floor
        ! at delta_ls = 0.01 by VALUE (checked inside the loop), rejection
        ! handled by the CFL-halving branch below.
        if (ANK_pzStepping) then
            lsTolEff = one
            lsFactorEff = 0.90_realType
            lsMaxIterEff = 64
        else
            lsTolEff = ANK_unstdyLSTol
            lsFactorEff = ankUnsteadyLSFactor
            lsMaxIterEff = ankUnsteadyLSMaxIter
        end if

        if ((unsteadyNorm > unsteadyNorm_old * lsTolEff .or. myisnan(unsteadyNorm))) then
            ! The unsteady residual is too high or we have a NAN. Do a
            ! backtracking line search until we get a residual that is lower.

            LSFailed = .True.
            backtrackTriggeredANK = .True.

            ! Restore the starting (old) w value by adding lamda*deltaW
            call VecAXPY(wVec, lambda, deltaW, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Set the initial new lambda. This is working off the
            ! potentially already physically limited step.
            ! VERIF_06 F2: factor and budget are options. Defaults 0.7/12
            ! reproduce ADflow exactly, whose floor 0.7**12 = 0.0138 sits ABOVE
            ! the 0.01 rejection threshold and so can never trigger a cutback.
            ! The thesis's Algorithm 4 geometry is 0.90 with ~44 iterations.
            lambda = lsFactorEff * lambda

            backtrack: do iter = 1, lsMaxIterEff

                ! Apply the new step
                call VecAXPY(wVec, -lambda, deltaW, ierr)
                call EChk(ierr, __FILE__, __LINE__)

                ! Set and recompute
                call setWANKScaled(wVec, 1, nState)

                ! Compute the unsteady residuals with the current step
                call computeUnsteadyResANK(lambda)
                feval = feval + 1

                call VecNorm(rVec, NORM_2, unsteadyNorm, ierr)
                call EChk(ierr, __FILE__, __LINE__)

                if (unsteadyNorm > unsteadyNorm_old * lsTolEff .or. myisnan(unsteadyNorm)) then

                    ! Restore back to the original wVec
                    call VecAXPY(wVec, lambda, deltaW, ierr)
                    call EChk(ierr, __FILE__, __LINE__)

                    ! Haven't backed off enough yet....keep going
                    lambda = lambda * lsFactorEff
                    stallNBacktrack = iter

                    ! PZ mode / thesis Algorithm 4: the floor is on the VALUE of
                    ! delta_ls, not on an iteration budget. Below 1% the search
                    ! has failed; leave LSFailed set and let the rejection
                    ! branch below halve the reference CFL (or take the floored
                    ! step if the CFL is already at its floor).
                    if (ANK_pzStepping .and. lambda < 0.01_realType) exit backtrack
                else
                    ! We have succefssfully reduced the norm
                    LSFailed = .False.
                    stallNBacktrack = iter

                    ! VERIF_06 F2 / thesis Algorithm 4: a step accepted only
                    ! because the backtracking budget ran out is a failed step,
                    ! not a successful one. Without this the controller can
                    ! never classify the collapse as a failure, so the CFL
                    ! cutback never fires and the solver limit-cycles. Undo the
                    ! step and fall through to the rejection branch below, which
                    ! already handles restoring the state and cutting the CFL.
                    if ((.not. ANK_pzStepping) .and. ankRejectOnLSExhausted .and. &
                        lambda < ANK_stepMin * ANK_stepFactor) then
                        call VecAXPY(wVec, lambda, deltaW, ierr)
                        call EChk(ierr, __FILE__, __LINE__)
                        LSFailed = .True.
                    end if
                    exit
                end if
            end do backtrack
            stallLSFailed = LSFailed

            if (LSFailed .or. myisnan(unsteadyNorm)) then
                ! the line search wasn't much help.

                if (ANK_CFL > ANK_CFLMin) then
                    ! the cfl number is not already at the lower limit.  We
                    ! can cut the CFL back and try again. Set lambda to zero
                    ! to indicate we never took a step.
                    lambda = zero
                else
                    ! cfl is as low as it goes, try taking the step
                    ! anyway. We can't do  anything else
                    call VecAXPY(wVec, -lambda, deltaW, ierr)
                    call EChk(ierr, __FILE__, __LINE__)
                end if

                ! Set the state vec and compute the new residual
                call setWANKScaled(wVec, 1, nState)
                if (.not. ANK_coupled) then
                    call blocketteRes(useTurbRes=.False., useStoreWall=.False.)
                else
                    call blocketteRes()
                end if
                feval = feval + 1
            else
            end if
        end if

        ! ============== Turb Update =============
        if ((.not. ANK_coupled) .and. equations == RANSEquations .and. lambda > zero) then

            if (ANK_useTurbDADI) then
                ! actually do the turbulence update
                call computeUtau
                call turbSolveDDADI
            else
                call ANKTurbSolveKSP
            end if
        end if

        ! ===== Algorithm 2 per-node damping, coupled path (VERIF_06 F0/F7) =====
        ! The thesis's inexact-Newton phase is CSANK, not NK, so the per-node
        ! gamma/Re-theta-t back-off has to live here too -- it previously had a
        ! single call site, inside NKStep. Placed after the line search has
        ! settled and before the residual evaluation below, mirroring NK's
        ! placement (after the search accepts, before the state is final), and
        ! covering every exit path of the search including a rejected step
        ! (lambda = 0 makes it a no-op).
        ! Removing gamma/Re-theta-t from the global step limit leaves Algorithm 2
        ! as the ONLY thing enforcing their bounds, so it must be active.
        if ((.not. ankTransitionGlobalLambda) .and. (.not. ankAlgorithm2Damping) &
            .and. ANK_coupled .and. transitionNK .and. myid == 0 .and. firstCall) then
            print *, 'Warning: ANKTransitionGlobalLambda=False requires ', &
                'ANKAlgorithm2Damping; enabling it so gamma/ReTheta stay bounded.'
        end if
        if ((ankAlgorithm2Damping .or. (.not. ankTransitionGlobalLambda)) &
            .and. ANK_coupled .and. transitionNK .and. &
            turbModel == spalartallmarasnoft2gammaretheta) then
            call applyANKAlgorithm2Damping(wVec, deltaW, lambda)
            call setWANKScaled(wVec, 1, nState)
        end if

        ! We need to now compute the residual for the next iteration.  We
        ! also need the to update the update the time step and the
        ! viscWall pointer stuff

        call blocketteRes(useUpdateIntermed=.True.)

        feval = feval + 1
        if (ANK_coupled) then
            call setRvec(rVec)
        else
            call setRVecANK(rVec)
        end if

        linResMonitor = resHist(kspIterations + 1) / resHist(1)

        if ((linResMonitor .ge. ANK_rtol .and. &
             totalR > ANK_secondOrdSwitchTol * totalR0 .and. &
             linResOld .le. ANK_rtol) &
            !.or. LSFailed) then
            !.or. lambda .le. ANK_stepMin) then
            .or. (lambda .eq. zero)) then
            ! We should reform the PC since it took longer than we want,
            ! or we need to adjust the CFL because the last update was bad,
            ! or convergence since the last PC update was good enough and we
            ! would benefit from re-calculating the PC.
            ANK_iter = -1
        end if

        ! update the linear residual for next iteration
        linResOld = linResMonitor

        ! Update step monitor
        stepMonitor = lambda

        ! ============== Stall diagnostics (VERIF_06) ==============
        ! Report WHY the step collapsed. Printed on the root proc only, and
        ! only for iterations whose accepted step is below solverStallDiagStep
        ! (default 1.0 => every iteration once the feature is enabled).
        !   lamPhys  : step surviving the physicality check
        !   bind     : variable+cell whose ratio produced lamPhys
        !   lamLS    : step after the unsteady-residual line search
        !   bt       : backtracks used (12 = budget exhausted)
        !   CFL/min  : current CFL and the floor the cutback is clipped to
        !   FLOORED  : a cutback was requested but ANK_CFLMin blocked it (F1)
        ! Oscillation detector: track consecutive iterations in which the total
        ! residual ROSE. Stalling (residual flat, step pinned) and oscillating
        ! (residual sawtoothing while the step swings) look identical in the
        ! standard columns but need different fixes, so separate them here.
        if (solverStallDiag) then
            if (stallResPrev > zero .and. totalR > stallResPrev) then
                stallRiseCount = stallRiseCount + 1
                stallRiseMax = max(stallRiseMax, stallRiseCount)
            else
                stallRiseCount = 0
            end if
            stallResPrev = totalR
        end if

        if (solverStallDiag .and. myid == 0 .and. lambda < solverStallDiagStep) then
            stallLamLS = lambda
            stallCFLAfter = ANK_CFL
            ! One line per iteration carrying every cause we can distinguish:
            !   lamPhys/bind/blk/ijk : physicality limit and the exact cell
            !   lamLS/bt/LSfail      : line-search outcome (bt at its budget
            !                          means the search ran out of room)
            !   CFL/CFLmin/FLOORED   : whether the controller could back off
            !   dampG/dampR/wf       : Algorithm 2 cells damped + worst factor
            !                          (a front being crushed every iteration)
            !   rise                 : consecutive residual increases (oscillation)
            write (*, "(a,es9.2,a,i6,a,es9.2,a,a,a,i3,a,i4,i4,i4,a,es9.2,a,i3,a,l1,a,es9.2,a,es9.2,a,l1,a,i7,i7,a,es9.2,a,i4)") &
                " STALLDIAG ANK rel=", totalR / totalR0, &
                "  pcIter=", ANK_iter, &
                "  lamPhys=", stallLamPhys, &
                "  bind=", trim(stallVarName(stallBindVar)), &
                "  blk=", stallBindLoc(1), &
                "  ijk=", stallBindLoc(2), stallBindLoc(3), stallBindLoc(4), &
                "  lamLS=", stallLamLS, &
                "  bt=", stallNBacktrack, &
                "  LSfail=", stallLSFailed, &
                "  CFL=", stallCFLAfter, &
                "  CFLmin=", ANK_CFLMin, &
                "  FLOORED=", stallCFLFloored, &
                "  damp=", stallSoftDampG, stallSoftDampR, &
                "  wf=", min(stallMinDampG, stallMinDampR), &
                "  rise=", stallRiseCount
        end if

        ! ============== Source-dt deactivation switch (P&Z §IV.B.3), coupled path ==============
        ! Mirror of the logic in ANKTurbSolveKSP: count clean iterations in the
        ! second-order (inexact-Newton analog) regime; reset the counter —
        ! reactivating the source-dt restriction — when backtracking is
        ! triggered, the step is rejected, or the relative residual rises back
        ! above the phase-switch tolerance. In segregated mode the counter is
        ! owned by ANKTurbSolveKSP, so only update it here when coupled.
        if (ANK_coupled .and. turbModel == spalartallmarasnoft2gammaretheta .and. transitionNK .and. &
            transitionSrcDtRestrict) then
            if (backtrackTriggeredANK .or. lambda == zero .or. &
                totalR > ANK_secondOrdSwitchTol * totalR0) then
                noBacktrackCount = 0
            else
                noBacktrackCount = noBacktrackCount + 1
            end if
        end if

        ! ============== PZ rejection bookkeeping ==============
        ! Thesis Algorithms 2 and 4: a rejected update (lambda = 0) halves the
        ! reference CFL (dt_ref = max(0.5 dt_ref, dt_ref,min)) and re-anchors
        ! the SER law so phase 2 continues from the halved value. An accepted
        ! step arms the growth law for the next iteration.
        if (ANK_pzStepping) then
            if (lambda == zero) then
                pz_CFLRef = max(half * pz_CFLRef, ANK_pzCFLMin)
                pz_alphaSER = -one
                pz_prevAccepted = .False.
                ANK_iter = -1 ! reform the PC with the halved CFL
            else
                pz_prevAccepted = .True.
            end if
        else
            ! Check if the linear solutions are failing.
            ! If the lin res is above .5 or so, the solver
            ! might stall, so we might be better off just
            ! reducing the CFL and keep going. We Modify
            ! the CFLMin by altering CFLMinBase.
            ! (ADflow controller only — the PZ system has no linear-residual
            ! CFL response; its rejection logic covers this.)
            if (linResMonitor .gt. ANK_linResMax) then
                ! This will adjust MinBase such that we can halve the cfl
                ! based on the current CFL.
                ANK_CFLMinBase = ANK_CFLCutback * ANK_CFL * ((totalR / totalR0)**ANK_CFLExponent)
                ! flags to refresh the Jacobian and cut back the CFL
                ANK_iter = -1
                lambda = zero
            end if
        end if

        ! Update the approximate iteration counter. The +1 is for the
        ! residual evaluations.
        approxTotalIts = approxTotalIts + feval + kspIterations

    end subroutine ANKStep
end module ANKSolver
