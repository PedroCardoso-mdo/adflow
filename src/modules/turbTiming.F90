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
    ! Context B: residual evaluation (blocketteRes / blockResCore, resOnly=.true.)
    integer(kind=intType), parameter :: T_RESID_FLOW  = 9
    integer(kind=intType), parameter :: T_RESID_TURB  = 10  ! SA-GR turb residual
    integer(kind=intType), parameter :: T_RESID_HALO  = 11
    integer(kind=intType), parameter :: T_RESID_TOTAL = 12

    integer(kind=intType), parameter :: nTurbTimers = 12

    real(kind=realType), dimension(nTurbTimers) :: tTurbStart = zero
    real(kind=realType), dimension(nTurbTimers) :: tTurbAccum = zero

contains

    subroutine turbTimingReset
        ! Zero all accumulators (call where t0Solver is set).
        tTurbAccum = zero
        tTurbStart = zero
    end subroutine turbTimingReset

    subroutine turbTic(id)
        integer(kind=intType), intent(in) :: id
        tTurbStart(id) = mpi_wtime()
    end subroutine turbTic

    subroutine turbToc(id)
        integer(kind=intType), intent(in) :: id
        tTurbAccum(id) = tTurbAccum(id) + mpi_wtime() - tTurbStart(id)
    end subroutine turbToc

    subroutine printTurbTiming
        !$ use omp_lib, only: omp_get_max_threads
        use communication, only: ADflow_comm_world, myID, nProc
        use inputDiscretization, only: useBlockettes
        implicit none
        real(kind=realType), dimension(nTurbTimers) :: tMax
        real(kind=realType) :: assembly, dadiSub, residSub, oth
        integer(kind=intType) :: ierr, nth
        character(len=3) :: blkStr

        ! Max-over-ranks for every timer (comm in particular is meaningful as max)
        call mpi_reduce(tTurbAccum, tMax, nTurbTimers, adflow_real, &
                        mpi_max, 0, ADflow_comm_world, ierr)

        if (myID /= 0) return

        nth = 1
        !$ nth = omp_get_max_threads()
        blkStr = 'OFF'
        if (useBlockettes) blkStr = 'ON '

        assembly = tMax(T_SOURCE) + tMax(T_VISCOUS) + tMax(T_ADV) + tMax(T_UNSTEADY)
        oth      = tMax(T_RESSCALE)
        dadiSub  = assembly + tMax(T_DADI) + oth + tMax(T_TURBHALO)
        residSub = tMax(T_RESID_FLOW) + tMax(T_RESID_TURB) + tMax(T_RESID_HALO)

        write (*, '(A)') '+--------------------------------------------------------------+'
        write (*, '(A,I0,A,I0,A,A)') '|  ADflow SA-gamma-Retheta timing   ranks=', nProc, &
            '  OMP=', nth, '  BLOCKETTES=', blkStr
        write (*, '(A)') '|--------------------------------------------------------------|'
        write (*, '(A)') '|  [A] DADI turbulence solve (turbSolveDDADI)                  |'
        write (*, '(A,F12.3,A)') '|    Source                      : ', tMax(T_SOURCE),   ' s'
        write (*, '(A,F12.3,A)') '|    Viscous                     : ', tMax(T_VISCOUS),  ' s'
        write (*, '(A,F12.3,A)') '|    Advection                   : ', tMax(T_ADV),      ' s'
        write (*, '(A,F12.3,A)') '|    Unsteady term               : ', tMax(T_UNSTEADY), ' s'
        write (*, '(A,F12.3,A)') '|      = Residual + qq assembly   : ', assembly, ' s  (sum of 4)'
        write (*, '(A,F12.3,A)') '|    DADI line-sweep solve        : ', tMax(T_DADI),     ' s'
        write (*, '(A,F12.3,A)') '|    ResScale + other             : ', oth,              ' s'
        write (*, '(A,F12.3,A)') '|    Turb halo (whalo2)           : ', tMax(T_TURBHALO), ' s'
        write (*, '(A,F12.3,A)') '|    DADI subtotal                : ', dadiSub,          ' s'
        write (*, '(A)') '|                                                              |'
        write (*, '(A,A,A)') '|  [B] Residual eval  (path: ', blkStr, ')                       |'
        write (*, '(A,F12.3,A)') '|    Flow residual                : ', tMax(T_RESID_FLOW), ' s'
        write (*, '(A,F12.3,A)') '|    SA-GR turbulence residual    : ', tMax(T_RESID_TURB), ' s'
        write (*, '(A,F12.3,A)') '|    Residual halo                : ', tMax(T_RESID_HALO), ' s'
        write (*, '(A,F12.3,A)') '|    Residual subtotal            : ', residSub,           ' s'
        write (*, '(A)') '|--------------------------------------------------------------|'
        write (*, '(A)') '|  (max-over-ranks; enable with -DTURB_TIMING)                 |'
        write (*, '(A)') '+--------------------------------------------------------------+'

    end subroutine printTurbTiming

end module turbTiming
