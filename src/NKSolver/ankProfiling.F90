module ankProfiling

    use constants
    use communication, only: adflow_comm_world, myID, nProc

    implicit none
    save

    integer(kind=intType), parameter :: ANK_CTX_NONE = 0_intType
    integer(kind=intType), parameter :: ANK_CTX_FORMJAC = 1_intType
    integer(kind=intType), parameter :: ANK_CTX_MATMULT = 2_intType
    integer(kind=intType), parameter :: ANK_CTX_UNSTEADY = 3_intType
    integer(kind=intType), parameter :: ANK_CTX_LINESEARCH = 4_intType
    integer(kind=intType), parameter :: ANK_CTX_FINALRES = 5_intType
    integer(kind=intType), parameter :: ANK_CTX_MATMULT_IN_KSPSOLVE = 6_intType
    integer(kind=intType), parameter :: ANK_CTX_MATMULT_OUTSIDE_KSPSOLVE = 7_intType
    integer(kind=intType), parameter :: ANK_CTX_TURBUPDATE = 8_intType

    integer(kind=intType), parameter :: ANK_SEC_ANKSTEP_TOTAL = 1_intType
    integer(kind=intType), parameter :: ANK_SEC_COMPUTETIMESTEPMAT = 2_intType
    integer(kind=intType), parameter :: ANK_SEC_FORMJAC_TOTAL = 3_intType
    integer(kind=intType), parameter :: ANK_SEC_FORMJAC_RES_TOTAL = 4_intType
    integer(kind=intType), parameter :: ANK_SEC_PCSETUP_TOTAL = 5_intType
    integer(kind=intType), parameter :: ANK_SEC_KSPSOLVE_TOTAL = 6_intType
    integer(kind=intType), parameter :: ANK_SEC_MATMULT_TOTAL = 7_intType
    integer(kind=intType), parameter :: ANK_SEC_PCAPPLY_TOTAL = 8_intType
    integer(kind=intType), parameter :: ANK_SEC_GMRES_ORTHOG_TOTAL = 9_intType
    integer(kind=intType), parameter :: ANK_SEC_COMPUTEUNSTEADY_TOTAL = 10_intType
    integer(kind=intType), parameter :: ANK_SEC_LINESEARCH_TOTAL = 11_intType
    integer(kind=intType), parameter :: ANK_SEC_FINALRES_TOTAL = 12_intType
    integer(kind=intType), parameter :: ANK_SEC_LOCALOPS_TOTAL = 13_intType
    integer(kind=intType), parameter :: ANK_SEC_BLOCKETTERES_TOTAL_IN_ANK = 14_intType
    integer(kind=intType), parameter :: ANK_SEC_BLOCKETTERES_COMPUTE_IN_ANK = 15_intType
    integer(kind=intType), parameter :: ANK_SEC_BLOCKETTERES_COMM_IN_ANK = 16_intType
    integer(kind=intType), parameter :: ANK_SEC_BLOCKETTERESCORE_TOTAL_IN_ANK = 17_intType
    integer(kind=intType), parameter :: ANK_SEC_WHALO2_TOTAL_IN_ANK = 18_intType
    integer(kind=intType), parameter :: ANK_SEC_BLOCKETTERES_IN_FORMJAC = 19_intType
    integer(kind=intType), parameter :: ANK_SEC_BLOCKETTERES_IN_MATMULT = 20_intType
    integer(kind=intType), parameter :: ANK_SEC_BLOCKETTERES_IN_UNSTEADY = 21_intType
    integer(kind=intType), parameter :: ANK_SEC_BLOCKETTERES_IN_LINESEARCH = 22_intType
    integer(kind=intType), parameter :: ANK_SEC_BLOCKETTERES_IN_FINALRES = 23_intType
    integer(kind=intType), parameter :: ANK_SEC_TURB_UPDATE_TOTAL = 24_intType
    integer(kind=intType), parameter :: ANK_SEC_MATMULT_IN_KSPSOLVE = 25_intType
    integer(kind=intType), parameter :: ANK_SEC_MATMULT_OUTSIDE_KSPSOLVE = 26_intType
    integer(kind=intType), parameter :: ANK_SEC_MATMULT_BLOCKETTERES_TOTAL = 27_intType
    integer(kind=intType), parameter :: ANK_SEC_MATMULT_TIMESTEPMATADD_TOTAL = 28_intType
    integer(kind=intType), parameter :: ANK_SEC_MATMULT_VECOPS_TOTAL = 29_intType
    integer(kind=intType), parameter :: ANK_SEC_COMPUTEUNSTEADY_BLOCKETTERES_TOTAL = 30_intType
    integer(kind=intType), parameter :: ANK_SEC_COMPUTEUNSTEADY_TIMESTEPMATMULT_TOTAL = 31_intType
    integer(kind=intType), parameter :: ANK_SEC_COMPUTEUNSTEADY_VECOPS_TOTAL = 32_intType
    integer(kind=intType), parameter :: ANK_SEC_BLOCKETTERES_IN_MATMULT_IN_KSPSOLVE = 33_intType
    integer(kind=intType), parameter :: ANK_SEC_BLOCKETTERES_IN_MATMULT_OUTSIDE_KSPSOLVE = 34_intType
    integer(kind=intType), parameter :: ANK_SEC_LOCALOPS_MATASSEMBLY_TOTAL = 35_intType
    integer(kind=intType), parameter :: ANK_SEC_LOCALOPS_OUTSIDE_KSPSOLVE_TOTAL = 36_intType
    integer(kind=intType), parameter :: ANK_SEC_LOCALOPS_PHYSICALITY_TOTAL = 37_intType
    integer(kind=intType), parameter :: ANK_SEC_LOCALOPS_STATEUPDATE_TOTAL = 38_intType
    integer(kind=intType), parameter :: ANK_SEC_LOCALOPS_VECNORM_TOTAL = 39_intType
    integer(kind=intType), parameter :: ANK_SEC_TURBUPDATE_KSP_TOTAL = 40_intType
    integer(kind=intType), parameter :: ANK_SEC_TURBUPDATE_DDADI_TOTAL = 41_intType
    integer(kind=intType), parameter :: ANK_SEC_TURBUPDATE_SA_BLOCKETTERES_TOTAL = 42_intType
    integer(kind=intType), parameter :: ANK_SEC_BLOCKETTERES_IN_TURBUPDATE = 43_intType
    integer(kind=intType), parameter :: ANK_SEC_TURBUPDATE_COMPUTE_OTHER_TOTAL = 44_intType
    integer(kind=intType), parameter :: ANK_SEC_APPLYSHELLPC = 45_intType
    integer(kind=intType), parameter :: ANK_SEC_MGPRECON = 46_intType
    integer(kind=intType), parameter :: ANK_SEC_KSPSOLVE_KSPLEVELS1 = 47_intType
    integer(kind=intType), parameter :: ANK_SEC_BLOCKETTERES_PREP_IN_ANK = 48_intType
    integer(kind=intType), parameter :: ANK_SEC_BLOCKETTERES_COPY_IN_ANK = 49_intType
    integer(kind=intType), parameter :: ANK_SEC_WHALO2_PACK_IN_ANK = 50_intType
    integer(kind=intType), parameter :: ANK_SEC_WHALO2_POST_IN_ANK = 51_intType
    integer(kind=intType), parameter :: ANK_SEC_WHALO2_LOCAL_IN_ANK = 52_intType
    integer(kind=intType), parameter :: ANK_SEC_WHALO2_UNPACK_IN_ANK = 53_intType
    integer(kind=intType), parameter :: ANK_SEC_WHALO2_RECV_IN_ANK = 54_intType
    integer(kind=intType), parameter :: ANK_SEC_WHALO2_SEND_IN_ANK = 55_intType
    integer(kind=intType), parameter :: ANK_N_SECTIONS = 55_intType

    integer(kind=intType), parameter :: ANK_CNT_N_ANKSTEP = 1_intType
    integer(kind=intType), parameter :: ANK_CNT_N_FORMJAC = 2_intType
    integer(kind=intType), parameter :: ANK_CNT_N_PCSETUP = 3_intType
    integer(kind=intType), parameter :: ANK_CNT_N_KSPSOLVE = 4_intType
    integer(kind=intType), parameter :: ANK_CNT_N_GMRES_ITER = 5_intType
    integer(kind=intType), parameter :: ANK_CNT_N_MATMULT = 6_intType
    integer(kind=intType), parameter :: ANK_CNT_N_PCAPPLY = 7_intType
    integer(kind=intType), parameter :: ANK_CNT_N_GMRES_ORTHOG = 8_intType
    integer(kind=intType), parameter :: ANK_CNT_N_BLOCKETTERES_TOTAL = 9_intType
    integer(kind=intType), parameter :: ANK_CNT_N_BLOCKETTERES_FORMJAC = 10_intType
    integer(kind=intType), parameter :: ANK_CNT_N_BLOCKETTERES_MATMULT = 11_intType
    integer(kind=intType), parameter :: ANK_CNT_N_BLOCKETTERES_UNSTEADY = 12_intType
    integer(kind=intType), parameter :: ANK_CNT_N_BLOCKETTERES_LINESEARCH = 13_intType
    integer(kind=intType), parameter :: ANK_CNT_N_BLOCKETTERES_FINALRES = 14_intType
    integer(kind=intType), parameter :: ANK_CNT_N_WHALO2 = 15_intType
    integer(kind=intType), parameter :: ANK_N_COUNTERS = 15_intType

    character(len=48), dimension(ANK_N_SECTIONS), parameter :: sectionNames = (/ character(len=48) :: &
        'ANKStep_total', &
        'computeTimeStepMat', &
        'FormJacobianANK_total', &
        '-->FormJacobianANK_residuals_total', &
        '-->PCSetup_total', &
        'KSPSolve_ANK_total', &
        '-->MatMult_ANK_total', &
        '-->PCApply_ANK_total', &
        '-->GMRES_orthog_total', &
        'computeUnsteadyResANK_total', &
        '-->ANKLineSearch_total', &
        'ANK_finalResidual_total', &
        'ANK_localOps_total', &
        'blocketteRes_total_in_ANK', &
        'blocketteRes_compute_in_ANK', &
        'blocketteRes_comm_in_ANK', &
        'blocketteResCore_total_in_ANK', &
        'whalo2_total_in_ANK', &
        'blocketteRes_in_FormJacobianANK', &
        'blocketteRes_in_MatMult_ANK', &
        'blocketteRes_in_computeUnsteadyResANK', &
        'blocketteRes_in_ANKLineSearch', &
        'blocketteRes_in_ANK_finalResidual', &
        'ANK_turbUpdate_total', &
        '-->MatMult_in_KSPSolve_ANK', &
        '-->MatMult_outside_KSPSolve_ANK', &
        '-->MatMult_blocketteRes_total', &
        '-->MatMult_timeStepMatAdd_total', &
        '-->MatMult_vecOps_total', &
        '-->computeUnsteadyResANK_blocketteRes_total', &
        '-->computeUnsteadyResANK_timeStepMatMult_total', &
        '-->computeUnsteadyResANK_vecOps_total', &
        'blocketteRes_in_MatMult_in_KSPSolve', &
        'blocketteRes_in_MFFD_base', &
        '-->matrixAssembly_total', &
        '-->outside_KSPSolve_total', &
        '-->physicality_total', &
        '-->stateUpdate_total', &
        '-->vecNorm_total', &
        '-->ANK_turbUpdate_KSP_total', &
        '-->ANK_turbUpdate_DDADI_total', &
        '-->ANK_turbUpdate_comm_SA_residual', &
        'blocketteRes_in_ANK_turbUpdate', &
        '-->ANK_turbUpdate_compute_other', &
        '-->applyShellPC', &
        '-->MGPreCon_in_applyShellPC', &
        '-->KSPSolve_kspLevels1_in_applyShellPC', &
        'blocketteRes_prep_in_ANK', &
        'blocketteRes_copy_in_ANK', &
        'whalo2_pack_in_ANK', &
        'whalo2_post_in_ANK', &
        'whalo2_local_in_ANK', &
        'whalo2_unpack_in_ANK', &
        'whalo2_recv_wait_in_ANK', &
        'whalo2_send_wait_in_ANK' /)

    character(len=48), dimension(ANK_N_COUNTERS), parameter :: counterNames = (/ character(len=48) :: &
        'n_ANKStep', &
        'n_FormJacobianANK', &
        'n_PCSetup_ANK', &
        'n_KSPSolve_ANK', &
        'n_GMRES_iter_ANK', &
        'n_MatMult_ANK', &
        'n_PCApply_ANK', &
        'n_GMRES_orthog_ANK', &
        'n_blocketteRes_total_in_ANK', &
        'n_blocketteRes_in_FormJacobianANK', &
        'n_blocketteRes_in_MatMult_ANK', &
        'n_blocketteRes_in_computeUnsteadyResANK', &
        'n_blocketteRes_in_ANKLineSearch', &
        'n_blocketteRes_in_ANK_finalResidual', &
        'n_whalo2_in_ANK' /)

    logical :: ankProfileEnabled = .false.
    logical :: ankStepActive = .false.
    integer(kind=intType) :: ankContext = ANK_CTX_NONE

    real(kind=alwaysRealType), dimension(ANK_N_SECTIONS) :: totalLocal = 0.0_alwaysRealType
    real(kind=alwaysRealType), dimension(ANK_N_SECTIONS) :: computeLocal = 0.0_alwaysRealType
    real(kind=alwaysRealType), dimension(ANK_N_SECTIONS) :: commLocal = 0.0_alwaysRealType
    integer(kind=intType), dimension(ANK_N_SECTIONS) :: callsLocal = 0_intType

    integer(kind=intType), dimension(ANK_N_COUNTERS) :: countersLocal = 0_intType

contains

    function ankNow() result(t)
        implicit none
        real(kind=alwaysRealType) :: t
        t = mpi_wtime()
    end function ankNow

    subroutine ankProfSetEnabled(flag)
        implicit none
        logical, intent(in) :: flag
        ankProfileEnabled = flag
    end subroutine ankProfSetEnabled

    subroutine ankProfReset()
        implicit none
        totalLocal = 0.0_alwaysRealType
        computeLocal = 0.0_alwaysRealType
        commLocal = 0.0_alwaysRealType
        callsLocal = 0_intType
        countersLocal = 0_intType
        ankContext = ANK_CTX_NONE
        ankStepActive = .false.
    end subroutine ankProfReset

    subroutine ankProfEnterStep()
        implicit none
        ankStepActive = .true.
        call ankProfIncrementCounter(ANK_CNT_N_ANKSTEP, 1_intType)
    end subroutine ankProfEnterStep

    subroutine ankProfExitStep()
        implicit none
        ankStepActive = .false.
        ankContext = ANK_CTX_NONE
    end subroutine ankProfExitStep

    logical function ankProfIsEnabled()
        implicit none
        ankProfIsEnabled = ankProfileEnabled
    end function ankProfIsEnabled

    logical function ankProfIsActive()
        implicit none
        ankProfIsActive = ankProfileEnabled .and. ankStepActive
    end function ankProfIsActive

    subroutine ankProfSetContext(context)
        implicit none
        integer(kind=intType), intent(in) :: context
        ankContext = context
    end subroutine ankProfSetContext

    integer(kind=intType) function ankProfGetContext()
        implicit none
        ankProfGetContext = ankContext
    end function ankProfGetContext

    subroutine ankProfIncrementCounter(counterID, amount)
        implicit none
        integer(kind=intType), intent(in) :: counterID, amount

        if (.not. ankProfileEnabled) return
        if (counterID < 1 .or. counterID > ANK_N_COUNTERS) return

!$omp critical(ankprof_counter_accum)
        countersLocal(counterID) = countersLocal(counterID) + amount
!$omp end critical(ankprof_counter_accum)
    end subroutine ankProfIncrementCounter

    subroutine ankProfAddSection(sectionID, dtTotal, dtCompute, dtComm)
        implicit none
        integer(kind=intType), intent(in) :: sectionID
        real(kind=alwaysRealType), intent(in) :: dtTotal, dtCompute, dtComm

        if (.not. ankProfileEnabled) return
        if (sectionID < 1 .or. sectionID > ANK_N_SECTIONS) return

!$omp critical(ankprof_section_accum)
        callsLocal(sectionID) = callsLocal(sectionID) + 1_intType
        totalLocal(sectionID) = totalLocal(sectionID) + dtTotal
        computeLocal(sectionID) = computeLocal(sectionID) + dtCompute
        commLocal(sectionID) = commLocal(sectionID) + dtComm
!$omp end critical(ankprof_section_accum)
    end subroutine ankProfAddSection

    subroutine ankProfAddSectionNoCall(sectionID, dtTotal, dtCompute, dtComm)
        implicit none
        integer(kind=intType), intent(in) :: sectionID
        real(kind=alwaysRealType), intent(in) :: dtTotal, dtCompute, dtComm

        if (.not. ankProfileEnabled) return
        if (sectionID < 1 .or. sectionID > ANK_N_SECTIONS) return

!$omp critical(ankprof_section_accum)
        totalLocal(sectionID) = totalLocal(sectionID) + dtTotal
        computeLocal(sectionID) = computeLocal(sectionID) + dtCompute
        commLocal(sectionID) = commLocal(sectionID) + dtComm
!$omp end critical(ankprof_section_accum)
    end subroutine ankProfAddSectionNoCall

    subroutine ankProfGetSectionAccum(sectionID, dtTotal, dtCompute, dtComm, nCalls)
        implicit none
        integer(kind=intType), intent(in) :: sectionID
        real(kind=alwaysRealType), intent(out) :: dtTotal, dtCompute, dtComm
        integer(kind=intType), intent(out) :: nCalls

        dtTotal = 0.0_alwaysRealType
        dtCompute = 0.0_alwaysRealType
        dtComm = 0.0_alwaysRealType
        nCalls = 0_intType

        if (.not. ankProfileEnabled) return
        if (sectionID < 1 .or. sectionID > ANK_N_SECTIONS) return

!$omp critical(ankprof_section_accum)
        dtTotal = totalLocal(sectionID)
        dtCompute = computeLocal(sectionID)
        dtComm = commLocal(sectionID)
        nCalls = callsLocal(sectionID)
!$omp end critical(ankprof_section_accum)
    end subroutine ankProfGetSectionAccum

    subroutine ankProfReport(useMatrixFree, precondType)
        implicit none
        logical, intent(in) :: useMatrixFree
        character(len=*), intent(in) :: precondType

        integer(kind=intType) :: ierr, i
        integer(kind=intType), dimension(ANK_N_SECTIONS) :: callsGlobal
        integer(kind=intType), dimension(ANK_N_COUNTERS) :: countersGlobal
        real(kind=alwaysRealType), dimension(ANK_N_SECTIONS) :: totalMax, totalSum
        real(kind=alwaysRealType), dimension(ANK_N_SECTIONS) :: computeMax, computeSum
        real(kind=alwaysRealType), dimension(ANK_N_SECTIONS) :: commMax, commSum
        real(kind=alwaysRealType), dimension(ANK_N_SECTIONS) :: totalAvg, computeAvg, commAvg
        real(kind=alwaysRealType) :: nProcReal

        integer(kind=intType), dimension(7) :: parentSections
        integer(kind=intType), dimension(2) :: formJacChildren
        integer(kind=intType), dimension(1) :: kspChildren
        integer(kind=intType), dimension(1) :: localOpsChildren
        integer(kind=intType), dimension(3) :: unsteadyChildren
        integer(kind=intType), dimension(3) :: kspMatvecChildren
        integer(kind=intType), dimension(3) :: residualTotals
        integer(kind=intType), dimension(2) :: residualComputeSubtotals
        integer(kind=intType), dimension(7) :: residualByContext
        real(kind=alwaysRealType) :: parentTotalAvg, parentComputeAvg, parentCommAvg
        real(kind=alwaysRealType) :: parentTotalDiff, parentComputeDiff, parentCommDiff

        if (.not. ankProfileEnabled) return
        if (countersLocal(ANK_CNT_N_ANKSTEP) <= 0_intType) return

        call mpi_allreduce(callsLocal, callsGlobal, ANK_N_SECTIONS, adflow_integer, mpi_sum, adflow_comm_world, ierr)
        call mpi_allreduce(totalLocal, totalMax, ANK_N_SECTIONS, MPI_DOUBLE, mpi_max, adflow_comm_world, ierr)
        call mpi_allreduce(totalLocal, totalSum, ANK_N_SECTIONS, MPI_DOUBLE, mpi_sum, adflow_comm_world, ierr)
        call mpi_allreduce(computeLocal, computeMax, ANK_N_SECTIONS, MPI_DOUBLE, mpi_max, adflow_comm_world, ierr)
        call mpi_allreduce(computeLocal, computeSum, ANK_N_SECTIONS, MPI_DOUBLE, mpi_sum, adflow_comm_world, ierr)
        call mpi_allreduce(commLocal, commMax, ANK_N_SECTIONS, MPI_DOUBLE, mpi_max, adflow_comm_world, ierr)
        call mpi_allreduce(commLocal, commSum, ANK_N_SECTIONS, MPI_DOUBLE, mpi_sum, adflow_comm_world, ierr)

        call mpi_allreduce(countersLocal, countersGlobal, ANK_N_COUNTERS, adflow_integer, mpi_sum, adflow_comm_world, ierr)

        nProcReal = real(nProc, alwaysRealType)
        totalAvg = totalSum / nProcReal
        computeAvg = computeSum / nProcReal
        commAvg = commSum / nProcReal

        if (myID /= 0) return

        parentSections = (/ANK_SEC_COMPUTETIMESTEPMAT, ANK_SEC_FORMJAC_TOTAL, ANK_SEC_KSPSOLVE_TOTAL, &
            ANK_SEC_LOCALOPS_TOTAL, ANK_SEC_COMPUTEUNSTEADY_TOTAL, ANK_SEC_FINALRES_TOTAL, ANK_SEC_TURB_UPDATE_TOTAL/)

        formJacChildren = (/ANK_SEC_FORMJAC_RES_TOTAL, ANK_SEC_PCSETUP_TOTAL/)
        kspChildren = (/ANK_SEC_MATMULT_IN_KSPSOLVE/)
        localOpsChildren = (/ANK_SEC_LOCALOPS_OUTSIDE_KSPSOLVE_TOTAL/)
        unsteadyChildren = (/ANK_SEC_COMPUTEUNSTEADY_BLOCKETTERES_TOTAL, ANK_SEC_COMPUTEUNSTEADY_TIMESTEPMATMULT_TOTAL, &
            ANK_SEC_COMPUTEUNSTEADY_VECOPS_TOTAL/)
        kspMatvecChildren = (/ANK_SEC_MATMULT_BLOCKETTERES_TOTAL, ANK_SEC_MATMULT_TIMESTEPMATADD_TOTAL, &
            ANK_SEC_MATMULT_VECOPS_TOTAL/)

        parentTotalAvg = sum(totalAvg(parentSections))
        parentComputeAvg = sum(computeAvg(parentSections))
        parentCommAvg = sum(commAvg(parentSections))

        parentTotalDiff = totalAvg(ANK_SEC_ANKSTEP_TOTAL) - parentTotalAvg
        parentComputeDiff = computeAvg(ANK_SEC_ANKSTEP_TOTAL) - parentComputeAvg
        parentCommDiff = commAvg(ANK_SEC_ANKSTEP_TOTAL) - parentCommAvg

        residualTotals = (/ANK_SEC_BLOCKETTERES_TOTAL_IN_ANK, &
            ANK_SEC_BLOCKETTERES_COMPUTE_IN_ANK, ANK_SEC_BLOCKETTERES_COMM_IN_ANK/)

        residualComputeSubtotals = (/ANK_SEC_BLOCKETTERES_PREP_IN_ANK, ANK_SEC_BLOCKETTERES_COPY_IN_ANK/)

        residualByContext = (/ANK_SEC_BLOCKETTERES_IN_FORMJAC, ANK_SEC_BLOCKETTERES_IN_MATMULT_IN_KSPSOLVE, &
            ANK_SEC_BLOCKETTERES_IN_MATMULT_OUTSIDE_KSPSOLVE, ANK_SEC_BLOCKETTERES_IN_UNSTEADY, &
            ANK_SEC_BLOCKETTERES_IN_LINESEARCH, ANK_SEC_BLOCKETTERES_IN_FINALRES, ANK_SEC_BLOCKETTERES_IN_TURBUPDATE/)

        write (*, '(A)') ' '
        write (*, '(A)') '================ ANK Profiling Report (MPI max/avg) ================'
        write (*, '(A)') 'Wall-clock interpretation: use max_* as the primary parallel metric.'
        write (*, '(A)') 'Comm columns only include directly measured communication timings.'
        write (*, '(A,L1,A,A,A)') 'ANK_useMatrixFree=', useMatrixFree, '  ANK_precondType=', trim(precondType), ' '
        write (*, '(A)') ' '

        write (*, '(A)') '1. ANKStep hierarchy (non-arrow rows are additive parents)'
        call writeHeader()
        call writeSectionRow(ANK_SEC_ANKSTEP_TOTAL)

        call writeSectionRow(ANK_SEC_COMPUTETIMESTEPMAT)

        call writeSectionRow(ANK_SEC_FORMJAC_TOTAL)
        call writeSectionRows(formJacChildren)

        call writeSectionRow(ANK_SEC_KSPSOLVE_TOTAL)
        call writeSectionRows(kspChildren)

        call writeSectionRow(ANK_SEC_LOCALOPS_TOTAL)
        call writeSectionRows(localOpsChildren)

        call writeSectionRow(ANK_SEC_COMPUTEUNSTEADY_TOTAL)
        call writeSectionRows(unsteadyChildren)

        call writeSectionRow(ANK_SEC_FINALRES_TOTAL)

        call writeSectionRow(ANK_SEC_TURB_UPDATE_TOTAL)

        write (*, '(A)') ' '
        write (*, '(A)') '1b. Parent closure check (avg columns, strict additive check)'
        write (*, '(A,1X,ES12.4,1X,ES12.4,1X,ES12.4)') 'sum(parent avg_total/avg_compute/avg_comm):', &
            parentTotalAvg, parentComputeAvg, parentCommAvg
        write (*, '(A,1X,ES12.4,1X,ES12.4,1X,ES12.4)') 'ANKStep avg_total/avg_compute/avg_comm:', &
            totalAvg(ANK_SEC_ANKSTEP_TOTAL), computeAvg(ANK_SEC_ANKSTEP_TOTAL), commAvg(ANK_SEC_ANKSTEP_TOTAL)
        write (*, '(A,1X,ES12.4,1X,ES12.4,1X,ES12.4)') 'difference (ANKStep - parent_sum):', &
            parentTotalDiff, parentComputeDiff, parentCommDiff

        write (*, '(A)') ' '
        write (*, '(A)') '2. Residual totals inside ANKStep'
        call writeHeader()
        call writeSectionRows(residualTotals)

        write (*, '(A)') ' '
        write (*, '(A)') '2b. Residual compute sub-breakdown (nested inside blocketteRes_compute_in_ANK)'
        call writeHeader()
        call writeSectionRows(residualComputeSubtotals)
        write (*, '(A48,1X,A8,1X,ES12.4,1X,ES12.4,1X,ES12.4,1X,ES12.4,1X,ES12.4,1X,ES12.4)') &
            'blocketteRes_compute_other', '       -', &
            totalMax(ANK_SEC_BLOCKETTERES_COMPUTE_IN_ANK) - totalMax(ANK_SEC_BLOCKETTERES_PREP_IN_ANK) - &
                totalMax(ANK_SEC_BLOCKETTERES_COPY_IN_ANK), &
            totalAvg(ANK_SEC_BLOCKETTERES_COMPUTE_IN_ANK) - totalAvg(ANK_SEC_BLOCKETTERES_PREP_IN_ANK) - &
                totalAvg(ANK_SEC_BLOCKETTERES_COPY_IN_ANK), &
            computeMax(ANK_SEC_BLOCKETTERES_COMPUTE_IN_ANK) - computeMax(ANK_SEC_BLOCKETTERES_PREP_IN_ANK) - &
                computeMax(ANK_SEC_BLOCKETTERES_COPY_IN_ANK), &
            computeAvg(ANK_SEC_BLOCKETTERES_COMPUTE_IN_ANK) - computeAvg(ANK_SEC_BLOCKETTERES_PREP_IN_ANK) - &
                computeAvg(ANK_SEC_BLOCKETTERES_COPY_IN_ANK), &
            commMax(ANK_SEC_BLOCKETTERES_COMPUTE_IN_ANK) - commMax(ANK_SEC_BLOCKETTERES_PREP_IN_ANK) - &
                commMax(ANK_SEC_BLOCKETTERES_COPY_IN_ANK), &
            commAvg(ANK_SEC_BLOCKETTERES_COMPUTE_IN_ANK) - commAvg(ANK_SEC_BLOCKETTERES_PREP_IN_ANK) - &
                commAvg(ANK_SEC_BLOCKETTERES_COPY_IN_ANK)

        write (*, '(A)') ' '
        write (*, '(A)') '3. Residual totals by calling context'
        call writeHeader()
        call writeSectionRows(residualByContext)
        write (*, '(A48,1X,A8,1X,ES12.4,1X,ES12.4,1X,ES12.4,1X,ES12.4,1X,ES12.4,1X,ES12.4)') &
            'blocketteRes_unaccounted', '       -', &
            totalMax(ANK_SEC_BLOCKETTERES_TOTAL_IN_ANK) - sum(totalMax(residualByContext)), &
            totalAvg(ANK_SEC_BLOCKETTERES_TOTAL_IN_ANK) - sum(totalAvg(residualByContext)), &
            computeMax(ANK_SEC_BLOCKETTERES_TOTAL_IN_ANK) - sum(computeMax(residualByContext)), &
            computeAvg(ANK_SEC_BLOCKETTERES_TOTAL_IN_ANK) - sum(computeAvg(residualByContext)), &
            commMax(ANK_SEC_BLOCKETTERES_TOTAL_IN_ANK) - sum(commMax(residualByContext)), &
            commAvg(ANK_SEC_BLOCKETTERES_TOTAL_IN_ANK) - sum(commAvg(residualByContext))

        write (*, '(A)') ' '
        write (*, '(A)') '4. KSPSolve callback-only subtimers'
        write (*, '(A)') '   (KSPSetTolerances/MatMFFDSetBase are setup calls and are not inside KSPSolve_total)'
        call writeHeader()
        call writeSectionRow(ANK_SEC_KSPSOLVE_TOTAL)
        call writeSectionRow(ANK_SEC_MATMULT_IN_KSPSOLVE)
        call writeSectionRows(kspMatvecChildren)
        write (*, '(A48,1X,A8,1X,ES12.4,1X,ES12.4,1X,ES12.4,1X,ES12.4,1X,ES12.4,1X,ES12.4)') &
            'KSPSolve_remainder', '       -', &
            totalMax(ANK_SEC_KSPSOLVE_TOTAL) - totalMax(ANK_SEC_MATMULT_IN_KSPSOLVE), &
            totalAvg(ANK_SEC_KSPSOLVE_TOTAL) - totalAvg(ANK_SEC_MATMULT_IN_KSPSOLVE), &
            computeMax(ANK_SEC_KSPSOLVE_TOTAL) - computeMax(ANK_SEC_MATMULT_IN_KSPSOLVE), &
            computeAvg(ANK_SEC_KSPSOLVE_TOTAL) - computeAvg(ANK_SEC_MATMULT_IN_KSPSOLVE), &
            commMax(ANK_SEC_KSPSOLVE_TOTAL) - commMax(ANK_SEC_MATMULT_IN_KSPSOLVE), &
            commAvg(ANK_SEC_KSPSOLVE_TOTAL) - commAvg(ANK_SEC_MATMULT_IN_KSPSOLVE)

        write (*, '(A)') ' '
        write (*, '(A)') '5. whalo2 compact decomposition (direct timers + residual other)'
        write (*, '(A)') '   buffer = pack + post + unpack; wait = recv_wait + send_wait'
        call writeHeader()
        call writeSectionRow(ANK_SEC_WHALO2_TOTAL_IN_ANK)
        write (*, '(A48,1X,A8,1X,ES12.4,1X,ES12.4,1X,ES12.4,1X,ES12.4,1X,ES12.4,1X,ES12.4)') &
            'whalo2_buffer_ops', '       -', &
            totalMax(ANK_SEC_WHALO2_PACK_IN_ANK) + totalMax(ANK_SEC_WHALO2_POST_IN_ANK) + &
                totalMax(ANK_SEC_WHALO2_UNPACK_IN_ANK), &
            totalAvg(ANK_SEC_WHALO2_PACK_IN_ANK) + totalAvg(ANK_SEC_WHALO2_POST_IN_ANK) + &
                totalAvg(ANK_SEC_WHALO2_UNPACK_IN_ANK), &
            computeMax(ANK_SEC_WHALO2_PACK_IN_ANK) + computeMax(ANK_SEC_WHALO2_POST_IN_ANK) + &
                computeMax(ANK_SEC_WHALO2_UNPACK_IN_ANK), &
            computeAvg(ANK_SEC_WHALO2_PACK_IN_ANK) + computeAvg(ANK_SEC_WHALO2_POST_IN_ANK) + &
                computeAvg(ANK_SEC_WHALO2_UNPACK_IN_ANK), &
            commMax(ANK_SEC_WHALO2_PACK_IN_ANK) + commMax(ANK_SEC_WHALO2_POST_IN_ANK) + &
                commMax(ANK_SEC_WHALO2_UNPACK_IN_ANK), &
            commAvg(ANK_SEC_WHALO2_PACK_IN_ANK) + commAvg(ANK_SEC_WHALO2_POST_IN_ANK) + &
                commAvg(ANK_SEC_WHALO2_UNPACK_IN_ANK)
        write (*, '(A48,1X,A8,1X,ES12.4,1X,ES12.4,1X,ES12.4,1X,ES12.4,1X,ES12.4,1X,ES12.4)') &
            'whalo2_local_exchange', '       -', &
            totalMax(ANK_SEC_WHALO2_LOCAL_IN_ANK), totalAvg(ANK_SEC_WHALO2_LOCAL_IN_ANK), &
            computeMax(ANK_SEC_WHALO2_LOCAL_IN_ANK), computeAvg(ANK_SEC_WHALO2_LOCAL_IN_ANK), &
            commMax(ANK_SEC_WHALO2_LOCAL_IN_ANK), commAvg(ANK_SEC_WHALO2_LOCAL_IN_ANK)
        write (*, '(A48,1X,A8,1X,ES12.4,1X,ES12.4,1X,ES12.4,1X,ES12.4,1X,ES12.4,1X,ES12.4)') &
            'whalo2_wait', '       -', &
            totalMax(ANK_SEC_WHALO2_RECV_IN_ANK) + totalMax(ANK_SEC_WHALO2_SEND_IN_ANK), &
            totalAvg(ANK_SEC_WHALO2_RECV_IN_ANK) + totalAvg(ANK_SEC_WHALO2_SEND_IN_ANK), &
            computeMax(ANK_SEC_WHALO2_RECV_IN_ANK) + computeMax(ANK_SEC_WHALO2_SEND_IN_ANK), &
            computeAvg(ANK_SEC_WHALO2_RECV_IN_ANK) + computeAvg(ANK_SEC_WHALO2_SEND_IN_ANK), &
            commMax(ANK_SEC_WHALO2_RECV_IN_ANK) + commMax(ANK_SEC_WHALO2_SEND_IN_ANK), &
            commAvg(ANK_SEC_WHALO2_RECV_IN_ANK) + commAvg(ANK_SEC_WHALO2_SEND_IN_ANK)
        write (*, '(A48,1X,A8,1X,ES12.4,1X,ES12.4,1X,ES12.4,1X,ES12.4,1X,ES12.4,1X,ES12.4)') &
            'whalo2_other', '       -', &
            totalMax(ANK_SEC_WHALO2_TOTAL_IN_ANK) - (totalMax(ANK_SEC_WHALO2_PACK_IN_ANK) + &
                totalMax(ANK_SEC_WHALO2_POST_IN_ANK) + totalMax(ANK_SEC_WHALO2_UNPACK_IN_ANK) + &
                totalMax(ANK_SEC_WHALO2_LOCAL_IN_ANK) + totalMax(ANK_SEC_WHALO2_RECV_IN_ANK) + &
                totalMax(ANK_SEC_WHALO2_SEND_IN_ANK)), &
            totalAvg(ANK_SEC_WHALO2_TOTAL_IN_ANK) - (totalAvg(ANK_SEC_WHALO2_PACK_IN_ANK) + &
                totalAvg(ANK_SEC_WHALO2_POST_IN_ANK) + totalAvg(ANK_SEC_WHALO2_UNPACK_IN_ANK) + &
                totalAvg(ANK_SEC_WHALO2_LOCAL_IN_ANK) + totalAvg(ANK_SEC_WHALO2_RECV_IN_ANK) + &
                totalAvg(ANK_SEC_WHALO2_SEND_IN_ANK)), &
            computeMax(ANK_SEC_WHALO2_TOTAL_IN_ANK) - (computeMax(ANK_SEC_WHALO2_PACK_IN_ANK) + &
                computeMax(ANK_SEC_WHALO2_POST_IN_ANK) + computeMax(ANK_SEC_WHALO2_UNPACK_IN_ANK) + &
                computeMax(ANK_SEC_WHALO2_LOCAL_IN_ANK) + computeMax(ANK_SEC_WHALO2_RECV_IN_ANK) + &
                computeMax(ANK_SEC_WHALO2_SEND_IN_ANK)), &
            computeAvg(ANK_SEC_WHALO2_TOTAL_IN_ANK) - (computeAvg(ANK_SEC_WHALO2_PACK_IN_ANK) + &
                computeAvg(ANK_SEC_WHALO2_POST_IN_ANK) + computeAvg(ANK_SEC_WHALO2_UNPACK_IN_ANK) + &
                computeAvg(ANK_SEC_WHALO2_LOCAL_IN_ANK) + computeAvg(ANK_SEC_WHALO2_RECV_IN_ANK) + &
                computeAvg(ANK_SEC_WHALO2_SEND_IN_ANK)), &
            commMax(ANK_SEC_WHALO2_TOTAL_IN_ANK) - (commMax(ANK_SEC_WHALO2_PACK_IN_ANK) + &
                commMax(ANK_SEC_WHALO2_POST_IN_ANK) + commMax(ANK_SEC_WHALO2_UNPACK_IN_ANK) + &
                commMax(ANK_SEC_WHALO2_LOCAL_IN_ANK) + commMax(ANK_SEC_WHALO2_RECV_IN_ANK) + &
                commMax(ANK_SEC_WHALO2_SEND_IN_ANK)), &
            commAvg(ANK_SEC_WHALO2_TOTAL_IN_ANK) - (commAvg(ANK_SEC_WHALO2_PACK_IN_ANK) + &
                commAvg(ANK_SEC_WHALO2_POST_IN_ANK) + commAvg(ANK_SEC_WHALO2_UNPACK_IN_ANK) + &
                commAvg(ANK_SEC_WHALO2_LOCAL_IN_ANK) + commAvg(ANK_SEC_WHALO2_RECV_IN_ANK) + &
                commAvg(ANK_SEC_WHALO2_SEND_IN_ANK))

        write (*, '(A)') '====================================================================='

    contains

        subroutine writeHeader()
            implicit none

            write (*, '(A)') 'name                                           calls    max_total     avg_total'// &
                '    max_compute   avg_compute   max_comm      avg_comm'
        end subroutine writeHeader

        subroutine writeSectionRows(sectionList)
            implicit none
            integer(kind=intType), intent(in), dimension(:) :: sectionList
            integer(kind=intType) :: ii, sid

            do ii = 1, size(sectionList)
                sid = sectionList(ii)
                call writeSectionRow(sid)
            end do
        end subroutine writeSectionRows

        subroutine writeSectionRow(sid)
            implicit none
            integer(kind=intType), intent(in) :: sid

            write (*, '(A48,1X,I8,1X,ES12.4,1X,ES12.4,1X,ES12.4,1X,ES12.4,1X,ES12.4,1X,ES12.4)') &
                trim(sectionNames(sid)), callsGlobal(sid), totalMax(sid), totalAvg(sid), &
                computeMax(sid), computeAvg(sid), commMax(sid), commAvg(sid)
        end subroutine writeSectionRow

    end subroutine ankProfReport

end module ankProfiling
