module turbTiming
    !
    !  Lightweight wall-clock timers for the SA-gamma-Retheta DADI turbulence
    !  solve and the residual evaluation, used to characterize where time is
    !  spent and the OpenMP / blockette speedups.
    !
    !  ALL call sites are wrapped in #ifdef TURB_TIMING, so when TURB_TIMING is
    !  not defined nothing here is ever called (zero overhead). To enable, add
    !  -DTURB_TIMING to FF90_FLAGS in config.mk. To remove for production:
    !  delete this file and grep-delete the "#ifdef TURB_TIMING" blocks.
    !
    use constants
    implicit none
    save

    ! ---- Timer IDs ------------------------------------------------------
    ! Context A: DADI turbulence solve (turbSolveDDADI -> saGammaReTheta_block,
    !            resOnly = .false.). The first four sum to residual+qq assembly.
    integer(kind=intType), parameter :: T_SOURCE     = 1
    integer(kind=intType), parameter :: T_VISCOUS    = 2
    integer(kind=intType), parameter :: T_ADV        = 3
    integer(kind=intType), parameter :: T_UNSTEADY   = 4
    integer(kind=intType), parameter :: T_DADI       = 5   ! line-sweep solve
    integer(kind=intType), parameter :: T_RESSCALE   = 6
    integer(kind=intType), parameter :: T_TURBHALO   = 7   ! whalo2 in DDADI
    integer(kind=intType), parameter :: T_TURBTOTAL  = 8   ! whole turbSolveDDADI
    ! Context B: residual evaluation (blocketteRes -> blockResCore / blocketteResCore)
    integer(kind=intType), parameter :: T_RESID_FLOW  = 9   ! flow-only resid evals
    integer(kind=intType), parameter :: T_RESID_TURB  = 10  ! SA-GR turb-only resid evals
    integer(kind=intType), parameter :: T_RESID_HALO  = 11
    integer(kind=intType), parameter :: T_RESID_TOTAL = 12
    ! Fused flow+turb residual evals. In the non-blockette path (blockResCore) the
    ! flow and turb sections are timed separately even within a fused call, so this
    ! stays ~0 there. In the blockette path (blocketteResCore) flow and turb are
    ! computed in one fused OpenMP sweep that cannot be split cheaply or race-free,
    ! so a fused call (flowRes .and. turbRes) is timed as a whole and lands here
    ! rather than being mis-attributed to FLOW/TURB. Pure-flow / pure-turb calls
    ! still land in FLOW / TURB in both paths and ARE directly comparable across
    ! blk=False <-> blk=True (that is the apples-to-apples number for Study 3).
    integer(kind=intType), parameter :: T_RESID_BOTH  = 13

    ! Thread-summed (CPU-seconds, NOT wall) accumulators for the turb-assembly vs
    ! flow-flux sections INSIDE the fused blockette OpenMP sweep. Summed per-thread
    ! via an OpenMP reduction in blocketteResCore, so they cannot be compared as
    ! wall time directly. Their RATIO r = TURB_TS / (TURB_TS + FLOW_TS) is
    ! thread-count-robust and is used in printTurbTiming to split the fused
    ! wall-time (T_RESID_BOTH) into a turbulent and a flow share. This is the only
    ! way to isolate the turbulent part of a fused blockette residual eval, since
    ! the two sections run interleaved in one cache-blocked sweep.
    integer(kind=intType), parameter :: T_BLKT_TURB_TS = 14
    integer(kind=intType), parameter :: T_BLKT_FLOW_TS = 15

    integer(kind=intType), parameter :: nTurbTimers = 15

    ! ---- Global work counters (for offline per-iteration normalization) -
    ! Headline OpenMP metric = (DADI subtotal) / nDADISweep.
    ! Blockette per-eval metric = T_RESID_TURB / nResidTurb (pure-turb evals).
    integer(kind=intType), parameter :: C_DADI  = 1   ! DADI line-sweep solves (sub-iters)
    integer(kind=intType), parameter :: C_RTURB = 2   ! pure-turb residual evals
    integer(kind=intType), parameter :: C_RFLOW = 3   ! pure-flow residual evals
    integer(kind=intType), parameter :: C_RBOTH = 4   ! fused flow+turb residual evals
    integer(kind=intType), parameter :: nTurbCounters = 4

    real(kind=realType), dimension(nTurbTimers) :: tTurbStart = zero
    real(kind=realType), dimension(nTurbTimers) :: tTurbAccum = zero
    integer(kind=intType), dimension(nTurbCounters) :: tTurbCount = 0

contains

    subroutine turbTimingReset
        ! Zero all accumulators and counters (call where t0Solver is set).
        tTurbAccum = zero
        tTurbStart = zero
        tTurbCount = 0
    end subroutine turbTimingReset

    subroutine turbTic(id)
        integer(kind=intType), intent(in) :: id
        tTurbStart(id) = mpi_wtime()
    end subroutine turbTic

    subroutine turbToc(id)
        integer(kind=intType), intent(in) :: id
        tTurbAccum(id) = tTurbAccum(id) + mpi_wtime() - tTurbStart(id)
    end subroutine turbToc

    subroutine turbCount(id)
        integer(kind=intType), intent(in) :: id
        tTurbCount(id) = tTurbCount(id) + 1
    end subroutine turbCount

    subroutine turbAddTS(turbSec, flowSec)
        ! Add the thread-summed turb/flow section times from one fused blockette
        ! sweep. Called once, serially, after the OpenMP parallel region closes.
        real(kind=realType), intent(in) :: turbSec, flowSec
        tTurbAccum(T_BLKT_TURB_TS) = tTurbAccum(T_BLKT_TURB_TS) + turbSec
        tTurbAccum(T_BLKT_FLOW_TS) = tTurbAccum(T_BLKT_FLOW_TS) + flowSec
    end subroutine turbAddTS

    subroutine printTurbTiming
        !$ use omp_lib, only: omp_get_max_threads
        use communication, only: ADflow_comm_world, myID, nProc
        use inputDiscretization, only: useBlockettes
        use iteration, only: t0Solver
        implicit none
        real(kind=realType), dimension(nTurbTimers) :: tMax, tSum
        integer(kind=intType), dimension(nTurbCounters) :: cMax
        real(kind=realType) :: assembly, dadiSub, residSub, oth, wallLoc, wallMax
        real(kind=realType) :: rTurb, fusedTurb, fusedFlow
        integer(kind=intType) :: ierr, nth

        ! Reductions to rank 0: time as max-over-ranks (the critical path) and as
        ! sum-over-ranks (-> mean) so a single straggler at high rank count is
        ! visible via the max/mean imbalance ratio rather than hidden by max-only.
        call mpi_reduce(tTurbAccum, tMax, nTurbTimers, adflow_real, &
                        mpi_max, 0, ADflow_comm_world, ierr)
        call mpi_reduce(tTurbAccum, tSum, nTurbTimers, adflow_real, &
                        mpi_sum, 0, ADflow_comm_world, ierr)
        ! Counts: every rank marches in lockstep on iteration count, so max == mean.
        call mpi_reduce(tTurbCount, cMax, nTurbCounters, adflow_integer, &
                        mpi_max, 0, ADflow_comm_world, ierr)
        ! Total solve wall (max over ranks).
        wallLoc = mpi_wtime() - t0Solver
        call mpi_reduce(wallLoc, wallMax, 1, adflow_real, &
                        mpi_max, 0, ADflow_comm_world, ierr)

        if (myID /= 0) return

        nth = 1
        !$ nth = omp_get_max_threads()

        assembly = tMax(T_SOURCE) + tMax(T_VISCOUS) + tMax(T_ADV) + tMax(T_UNSTEADY)
        oth      = tMax(T_RESSCALE)
        dadiSub  = assembly + tMax(T_DADI) + oth + tMax(T_TURBHALO)
        residSub = tMax(T_RESID_FLOW) + tMax(T_RESID_TURB) + tMax(T_RESID_BOTH) + tMax(T_RESID_HALO)

        ! Split the fused blockette wall-time into turb/flow by the global
        ! thread-summed ratio (see T_BLKT_TURB_TS comment). Uses sum-over-ranks.
        rTurb = zero
        if (tSum(T_BLKT_TURB_TS) + tSum(T_BLKT_FLOW_TS) > 1.0e-12_realType) &
            rTurb = tSum(T_BLKT_TURB_TS) / (tSum(T_BLKT_TURB_TS) + tSum(T_BLKT_FLOW_TS))
        fusedTurb = rTurb * tMax(T_RESID_BOTH)
        fusedFlow = (one - rTurb) * tMax(T_RESID_BOTH)

        write (*, '(A)') '+----------------------------------------------------------------------+'
        write (*, '(A,I0,A,I0,A,A)') '|  ADflow SA-gamma-Retheta timing   ranks=', nProc, &
            '  OMP=', nth, '  BLOCKETTES=', trim(blkLabel(useBlockettes))
        write (*, '(A,F12.3,A)') '|  total solve wall              : ', wallMax, ' s'
        write (*, '(A)') '|----------------------------------------------------------------------|'
        write (*, '(A)') '|  phase                            max(s)     mean(s)   imb           |'
        write (*, '(A)') '|  [A] DADI turbulence solve (turbSolveDDADI)                          |'
        call line('Source',                 T_SOURCE,   tMax, tSum)
        call line('Viscous',                T_VISCOUS,  tMax, tSum)
        call line('Advection',              T_ADV,      tMax, tSum)
        call line('Unsteady term',          T_UNSTEADY, tMax, tSum)
        write (*, '(A,F10.3,A)') '|      = Residual + qq assembly   ', assembly, '                       |'
        call line('DADI line-sweep solve',  T_DADI,     tMax, tSum)
        call line('ResScale + other',       T_RESSCALE, tMax, tSum)
        call line('Turb halo (whalo2)',     T_TURBHALO, tMax, tSum)
        write (*, '(A,F10.3,A)') '|    DADI subtotal                ', dadiSub, '                       |'
        write (*, '(A)') '|                                                                      |'
        write (*, '(A,A,A)') '|  [B] Residual eval  (path: ', trim(blkLabel(useBlockettes)), ')                                  |'
        call line('Flow residual (pure)',   T_RESID_FLOW, tMax, tSum)
        call line('SA-GR turb resid (pure)', T_RESID_TURB, tMax, tSum)
        call line('Flow+turb fused',        T_RESID_BOTH, tMax, tSum)
        write (*, '(A,F6.1,A,F10.3,A)') '|      -> turb  share (', 100.0_realType * rTurb, &
            ' %) wall ', fusedTurb, '             |'
        write (*, '(A,F6.1,A,F10.3,A)') '|      -> flow  share (', 100.0_realType * (one - rTurb), &
            ' %) wall ', fusedFlow, '             |'
        call line('Residual halo',          T_RESID_HALO, tMax, tSum)
        write (*, '(A,F10.3,A)') '|    Residual subtotal            ', residSub, '                       |'
        write (*, '(A)') '|----------------------------------------------------------------------|'
        write (*, '(A)') '|  work counts (max over ranks)        per-call (max-time/count)        |'
        call cline('DADI sweeps',           cMax(C_DADI),  tMax(T_DADI))
        call cline('turb-resid evals(pure)', cMax(C_RTURB), tMax(T_RESID_TURB))
        call cline('flow-resid evals(pure)', cMax(C_RFLOW), tMax(T_RESID_FLOW))
        call cline('fused resid evals',     cMax(C_RBOTH), tMax(T_RESID_BOTH))
        write (*, '(A)') '|----------------------------------------------------------------------|'
        write (*, '(A)') '|  imb = max/mean (load imbalance); per-call in ms; -DTURB_TIMING       |'
        write (*, '(A)') '+----------------------------------------------------------------------+'

    contains

        subroutine line(label, id, tmx, tsm)
            character(len=*), intent(in) :: label
            integer(kind=intType), intent(in) :: id
            real(kind=realType), dimension(:), intent(in) :: tmx, tsm
            real(kind=realType) :: mean, imb
            mean = tsm(id) / real(nProc, realType)
            imb = zero
            if (mean > 1.0e-12_realType) imb = tmx(id) / mean
            write (*, '(A,A24,F10.3,F10.3,F8.2,A)') '|    ', label, tmx(id), mean, imb, '   |'
        end subroutine line

        subroutine cline(label, cnt, tmx)
            character(len=*), intent(in) :: label
            integer(kind=intType), intent(in) :: cnt
            real(kind=realType), intent(in) :: tmx
            real(kind=realType) :: perCallMs
            perCallMs = zero
            if (cnt > 0) perCallMs = 1000.0_realType * tmx / real(cnt, realType)
            write (*, '(A,A24,I10,A,F12.4,A)') '|    ', label, cnt, '   ', perCallMs, ' ms   |'
        end subroutine cline

    end subroutine printTurbTiming

    function blkLabel(on) result(s)
        logical, intent(in) :: on
        character(len=3) :: s
        s = 'OFF'
        if (on) s = 'ON '
    end function blkLabel

end module turbTiming
