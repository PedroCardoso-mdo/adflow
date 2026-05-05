











module ALEUtils
contains

!
!      ******************************************************************
!      *                                                                *
!      * File:          slipVelocities.f90                              *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-12-2004                                      *
!      * Last modified: 06-28-2005                                      *
!      *                                                                *
!      ******************************************************************
!
    subroutine slipVelocitiesFineLevelALE_block(useOldCoor, t, sps)
        !
        !      ******************************************************************
        !      *                                                                *
        !      * slipVelocitiesFineLevel computes the slip velocities for       *
        !      * viscous subfaces on all viscous boundaries on groundLevel for  *
        !      * the given spectral solution. If useOldCoor is .true. the       *
        !      * velocities are determined using the unsteady time integrator;  *
        !      * otherwise the analytic form is used.                           *
        !      *                                                                *
        !      * Calculates the surface normal and normal velocity on BC using  *
        !      * FIRST order BDF.                                               *
        !      *                                                                *
        !      ******************************************************************
        !
        use constants
        use inputTimeSpectral
        use blockPointers
        use cgnsGrid
        use flowVarRefState
        use inputMotion
        use inputUnsteady
        use iteration
        use inputPhysics
        use inputTSStabDeriv
        use monitor
        use communication
        use utils, only: setCoefTimeIntegrator
        implicit none
        !
        !      Subroutine arguments.
        !
        integer(kind=intType), intent(in) :: sps
        logical, intent(in) :: useOldCoor

        real(kind=realType), dimension(*), intent(in) :: t
        !
        !      Local variables.
        !
        integer(kind=intType) :: nn, mm, i, j, level

        real(kind=realType) :: oneOver4dt
        real(kind=realType) :: velxGrid, velyGrid, velzGrid, ainf
        real(kind=realType) :: velxGrid0, velyGrid0, velzGrid0

        real(kind=realType), dimension(3) :: xc, xxc
        real(kind=realType), dimension(3) :: rotCenter, rotRate

        real(kind=realType), dimension(3) :: rotationPoint
        real(kind=realType), dimension(3, 3) :: rotationMatrix, &
                                                derivRotationMatrix

        real(kind=realType) :: tNew, tOld

        real(kind=realType), dimension(:, :, :), pointer :: uSlip
        real(kind=realType), dimension(:, :, :), pointer :: xFace
        real(kind=realType), dimension(:, :, :, :), pointer :: xFaceOld

        real(kind=realType) :: intervalMach, alphaTS, alphaIncrement, &
                               betaTS, betaIncrement
        real(kind=realType), dimension(3) :: velDir
        real(kind=realType), dimension(3) :: refDirection

        !Function Definitions

        real(kind=realType) :: TSAlpha, TSBeta, TSMach
        !
        !      ******************************************************************
        !      *                                                                *
        !      * Begin execution                                                *
        !      *                                                                *
        !      ******************************************************************
        !
        ! Determine the situation we are having here.

        ! *******************************
        ! REMOVED the rigid body rotation part for simplicity
        ! *******************************

        ! The velocities must be determined via a finite difference
        ! formula using the coordinates of the old levels.

        ! Set the coefficients for the time integrator and store the
        ! inverse of the physical nonDimensional time step, divided
        ! by 4, a bit easier.

        call setCoefTimeIntegrator
        oneOver4dt = fourth * timeRef / deltaT

        ! Loop over the number of viscous subfaces.

        bocoLoop1: do mm = 1, nViscBocos

            ! Set the pointer for uSlip to make the code more
            ! readable.

            uSlip => BCData(mm)%uSlip

            ! Determine the grid face on which the subface is located
            ! and set some variables accordingly.

            select case (BCFaceID(mm))

            case (iMin)
                xFace => x(1, :, :, :); xFaceOld => xOld(:, 1, :, :, :)

            case (iMax)
                xFace => x(il, :, :, :); xFaceOld => xOld(:, il, :, :, :)

            case (jMin)
                xFace => x(:, 1, :, :); xFaceOld => xOld(:, :, 1, :, :)

            case (jMax)
                xFace => x(:, jl, :, :); xFaceOld => xOld(:, :, jl, :, :)

            case (kMin)
                xFace => x(:, :, 1, :); xFaceOld => xOld(:, :, :, 1, :)

            case (kMax)
                xFace => x(:, :, kl, :); xFaceOld => xOld(:, :, :, kl, :)

            end select

            ! Some boundary faces have a different rotation speed than
            ! the corresponding block. This happens e.g. in the tip gap
            ! region of turboMachinary problems where the casing does
            ! not rotate. As the coordinate difference corresponds to
            ! the rotation rate of the block, a correction must be
            ! computed. Therefore compute the difference in rotation
            ! rate and store the rotation center a bit easier. Note that
            ! the rotation center of subface is taken, because if there
            ! is a difference in rotation rate this info for the subface
            ! must always be specified.

            j = nbkGlobal
            i = cgnsSubface(mm)

            ! Loop over the quadrilateral faces of the viscous subface.
            ! Note that due to the usage of the pointers xFace and
            ! xFaceOld an offset of +1 must be used in the coordinate
            ! arrays, because x and xOld originally start at 0 for the
            ! i, j and k indices.

            do j = BCData(mm)%jcBeg, BCData(mm)%jcEnd
                do i = BCData(mm)%icBeg, BCData(mm)%icEnd

                    ! Determine the coordinates of the centroid of the
                    ! face, multiplied by 4.

                    uSlip(i, j, 1) = (xFace(i + 1, j + 1, 1) + xFace(i + 1, j, 1) &
                                      + xFace(i, j + 1, 1) + xFace(i, j, 1))
                    uSlip(i, j, 2) = (xFace(i + 1, j + 1, 2) + xFace(i + 1, j, 2) &
                                      + xFace(i, j + 1, 2) + xFace(i, j, 2))
                    uSlip(i, j, 3) = (xFace(i + 1, j + 1, 3) + xFace(i + 1, j, 3) &
                                      + xFace(i, j + 1, 3) + xFace(i, j, 3))

                    ! Loop over the older time levels and take their
                    ! contribution into account.

                    level = 1 ! There was a loop over all old levels
                    uSlip(i, j, 1) = uSlip(i, j, 1) &
                                     + (xFaceOld(level, i + 1, j + 1, 1) &
                                        + xFaceOld(level, i + 1, j, 1) &
                                        + xFaceOld(level, i, j + 1, 1) &
                                        + xFaceOld(level, i, j, 1)) &
                                     * (-1.0_realType)
                    uSlip(i, j, 2) = uSlip(i, j, 2) &
                                     + (xFaceOld(level, i + 1, j + 1, 2) &
                                        + xFaceOld(level, i + 1, j, 2) &
                                        + xFaceOld(level, i, j + 1, 2) &
                                        + xFaceOld(level, i, j, 2)) &
                                     * (-1.0_realType)
                    uSlip(i, j, 3) = uSlip(i, j, 3) &
                                     + (xFaceOld(level, i + 1, j + 1, 3) &
                                        + xFaceOld(level, i + 1, j, 3) &
                                        + xFaceOld(level, i, j + 1, 3) &
                                        + xFaceOld(level, i, j, 3)) &
                                     * (-1.0_realType)

                    ! Divide by 4 times the time step to obtain the
                    ! correct velocity.

                    uSlip(i, j, 1) = uSlip(i, j, 1) * oneOver4dt
                    uSlip(i, j, 2) = uSlip(i, j, 2) * oneOver4dt
                    uSlip(i, j, 3) = uSlip(i, j, 3) * oneOver4dt
                end do
            end do

        end do bocoLoop1

    end subroutine slipVelocitiesFineLevelALE_block

! ===========================================================
    subroutine interpLevelALE_block
        !
        !      ******************************************************************
        !      *                                                                *
        !      * interpLevelALE_block interpolates geometric data over the      *
        !      * latest time step.                                              *
        !      *                                                                *
        !      ******************************************************************
        !
        use blockPointers
        use iteration
        use inputUnsteady
        use inputPhysics
        implicit none
        !
        !      Local variables.
        !
        integer(kind=intType) :: i, j, k, l, nn, mm, kk

        if (.not. useALE .or. equationMode .ne. unsteady) then
            return
        end if

        ! --------------------------------
        ! First store then clear current data
        ! --------------------------------
        clearI: do k = 1, ke
            do j = 1, je
                do i = 0, ie
                    sFaceIALE(0, i, j, k) = sFaceI(i, j, k)
                    sIALE(0, i, j, k, 1) = sI(i, j, k, 1)
                    sIALE(0, i, j, k, 2) = sI(i, j, k, 2)
                    sIALE(0, i, j, k, 3) = sI(i, j, k, 3)
                    sFaceI(i, j, k) = zero
                    sI(i, j, k, 1) = zero
                    sI(i, j, k, 2) = zero
                    sI(i, j, k, 3) = zero
                end do
            end do
        end do clearI

        clearJ: do k = 1, ke
            do j = 0, je
                do i = 1, ie
                    sFaceJALE(0, i, j, k) = sFaceJ(i, j, k)
                    sJALE(0, i, j, k, 1) = sJ(i, j, k, 1)
                    sJALE(0, i, j, k, 2) = sJ(i, j, k, 2)
                    sJALE(0, i, j, k, 3) = sJ(i, j, k, 3)
                    sFaceJ(i, j, k) = zero
                    sJ(i, j, k, 1) = zero
                    sJ(i, j, k, 2) = zero
                    sJ(i, j, k, 3) = zero
                end do
            end do
        end do clearJ

        clearK: do k = 0, ke
            do j = 1, je
                do i = 1, ie
                    sFaceKALE(0, i, j, k) = sFaceK(i, j, k)
                    sKALE(0, i, j, k, 1) = sK(i, j, k, 1)
                    sKALE(0, i, j, k, 2) = sK(i, j, k, 2)
                    sKALE(0, i, j, k, 3) = sK(i, j, k, 3)
                    sFaceK(i, j, k) = zero
                    sK(i, j, k, 1) = zero
                    sK(i, j, k, 2) = zero
                    sK(i, j, k, 3) = zero
                end do
            end do
        end do clearK

        ALEloop: do l = 1, nALEsteps
            ! --------------------------------
            ! Then average surface normal and normal velocity from array of old variables
            ! This eq. 10a and 10b, found paper by C.Farhat http://dx.doi.org/10.1016/S0021-9991(03)00311-5
            ! --------------------------------
            updateI: do k = 1, ke
                do j = 1, je
                    do i = 0, ie
                        sFaceI(i, j, k) = sFaceI(i, j, k) + coefTimeALE(l) * sFaceIALE(l, i, j, k)
                        sI(i, j, k, 1) = sI(i, j, k, 1) + coefTimeALE(l) * sIALE(l, i, j, k, 1)
                        sI(i, j, k, 2) = sI(i, j, k, 2) + coefTimeALE(l) * sIALE(l, i, j, k, 2)
                        sI(i, j, k, 3) = sI(i, j, k, 3) + coefTimeALE(l) * sIALE(l, i, j, k, 3)
                    end do
                end do
            end do updateI

            updateJ: do k = 1, ke
                do j = 0, je
                    do i = 1, ie
                        sFaceJ(i, j, k) = sFaceJ(i, j, k) + coefTimeALE(l) * sFaceJALE(l, i, j, k)
                        sJ(i, j, k, 1) = sJ(i, j, k, 1) + coefTimeALE(l) * sJALE(l, i, j, k, 1)
                        sJ(i, j, k, 2) = sJ(i, j, k, 2) + coefTimeALE(l) * sJALE(l, i, j, k, 2)
                        sJ(i, j, k, 3) = sJ(i, j, k, 3) + coefTimeALE(l) * sJALE(l, i, j, k, 3)
                    end do
                end do
            end do updateJ

            updateK: do k = 0, ke
                do j = 1, je
                    do i = 1, ie
                        sFaceK(i, j, k) = sFaceK(i, j, k) + coefTimeALE(l) * sFaceKALE(l, i, j, k)
                        sK(i, j, k, 1) = sK(i, j, k, 1) + coefTimeALE(l) * sKALE(l, i, j, k, 1)
                        sK(i, j, k, 2) = sK(i, j, k, 2) + coefTimeALE(l) * sKALE(l, i, j, k, 2)
                        sK(i, j, k, 3) = sK(i, j, k, 3) + coefTimeALE(l) * sKALE(l, i, j, k, 3)
                    end do
                end do
            end do updateK
        end do ALEloop

    end subroutine interpLevelALE_block

! ===========================================================
    subroutine recoverLevelALE_block
        !
        !      ******************************************************************
        !      *                                                                *
        !      * recoverLevelALE_block recovers current geometric data from     *
        !      * temporary interpolation                                        *
        !      *                                                                *
        !      ******************************************************************
        !
        use blockPointers
        use inputUnsteady
        use inputPhysics
        implicit none
        !
        !      Local variables.
        !
        integer(kind=intType) :: i, j, k, nn, mm, kk

        if (.not. useALE .or. equationMode .ne. unsteady) then
            return
        end if

        recoverI: do k = 1, ke
            do j = 1, je
                do i = 0, ie
                    sFaceI(i, j, k) = sFaceIALE(0, i, j, k)
                    sI(i, j, k, 1) = sIALE(0, i, j, k, 1)
                    sI(i, j, k, 2) = sIALE(0, i, j, k, 2)
                    sI(i, j, k, 3) = sIALE(0, i, j, k, 3)
                end do
            end do
        end do recoverI

        recoverJ: do k = 1, ke
            do j = 0, je
                do i = 1, ie
                    sFaceJ(i, j, k) = sFaceJALE(0, i, j, k)
                    sJ(i, j, k, 1) = sJALE(0, i, j, k, 1)
                    sJ(i, j, k, 2) = sJALE(0, i, j, k, 2)
                    sJ(i, j, k, 3) = sJALE(0, i, j, k, 3)
                end do
            end do
        end do recoverJ

        recoverK: do k = 0, ke
            do j = 1, je
                do i = 1, ie
                    sFaceK(i, j, k) = sFaceKALE(0, i, j, k)
                    sK(i, j, k, 1) = sKALE(0, i, j, k, 1)
                    sK(i, j, k, 2) = sKALE(0, i, j, k, 2)
                    sK(i, j, k, 3) = sKALE(0, i, j, k, 3)
                end do
            end do
        end do recoverK

    end subroutine recoverLevelALE_block

! ===========================================================
    subroutine interpLevelALEBC_block
        !
        !      ******************************************************************
        !      *                                                                *
        !      * interpLevelALEBC_block interpolates geometric data on boundary *
        !      * over the latest time step.                                     *
        !      *                                                                *
        !      ******************************************************************
        !
        use blockPointers
        use iteration
        use inputUnsteady
        use inputPhysics
        implicit none
        !
        !      Local variables.
        !
        integer(kind=intType) :: i, j, k, l, nn, mm, kk

        if (.not. useALE .or. equationMode .ne. unsteady) then
            return
        end if

        ! --------------------------------
        ! First store then clear current data
        ! --------------------------------
        clearNM: do mm = 1, nBocos
            do j = BCData(mm)%jcBeg, BCData(mm)%jcEnd
                do i = BCData(mm)%icBeg, BCData(mm)%icEnd
                    BCData(mm)%normALE(0, i, j, 1) = BCData(mm)%norm(i, j, 1)
                    BCData(mm)%normALE(0, i, j, 2) = BCData(mm)%norm(i, j, 2)
                    BCData(mm)%normALE(0, i, j, 3) = BCData(mm)%norm(i, j, 3)
                    BCData(mm)%norm(i, j, 1) = zero
                    BCData(mm)%norm(i, j, 2) = zero
                    BCData(mm)%norm(i, j, 3) = zero
                end do
            end do
        end do clearNM

        clearRF: do mm = 1, nBocos
            testAssoc1: if (associated(BCData(mm)%rFace)) then
                do j = BCData(mm)%jcBeg, BCData(mm)%jcEnd
                    do i = BCData(mm)%icBeg, BCData(mm)%icEnd
                        BCData(mm)%rFaceALE(0, i, j) = BCData(mm)%rFace(i, j)
                        BCData(mm)%rFace(i, j) = zero
                    end do
                end do
            end if testAssoc1
        end do clearRF

        clearUS: do mm = 1, nViscBocos
            do j = BCData(mm)%jcBeg, BCData(mm)%jcEnd
                do i = BCData(mm)%icBeg, BCData(mm)%icEnd
                    BCData(mm)%uSlipALE(0, i, j, 1) = BCData(mm)%uSlip(i, j, 1)
                    BCData(mm)%uSlipALE(0, i, j, 2) = BCData(mm)%uSlip(i, j, 2)
                    BCData(mm)%uSlipALE(0, i, j, 3) = BCData(mm)%uSlip(i, j, 3)
                    BCData(mm)%uSlip(i, j, 1) = zero
                    BCData(mm)%uSlip(i, j, 2) = zero
                    BCData(mm)%uSlip(i, j, 3) = zero
                end do
            end do
        end do clearUS

        ALEloop: do l = 1, nALEsteps
            ! --------------------------------
            ! Then average surface normal and normal velocity from array of old variables
            ! --------------------------------
            updateNM: do mm = 1, nBocos
                do j = BCData(mm)%jcBeg, BCData(mm)%jcEnd
                    do i = BCData(mm)%icBeg, BCData(mm)%icEnd
                        BCData(mm)%norm(i, j, 1) = BCData(mm)%norm(i, j, 1) &
                                                   + coefTimeALE(l) * BCData(mm)%normALE(l, i, j, 1)
                        BCData(mm)%norm(i, j, 2) = BCData(mm)%norm(i, j, 2) &
                                                   + coefTimeALE(l) * BCData(mm)%normALE(l, i, j, 2)
                        BCData(mm)%norm(i, j, 3) = BCData(mm)%norm(i, j, 3) &
                                                   + coefTimeALE(l) * BCData(mm)%normALE(l, i, j, 3)
                    end do
                end do
            end do updateNM

            updateRF: do mm = 1, nBocos
                testAssoc2: if (associated(BCData(mm)%rFace)) then
                    do j = BCData(mm)%jcBeg, BCData(mm)%jcEnd
                        do i = BCData(mm)%icBeg, BCData(mm)%icEnd
                            BCData(mm)%rFace(i, j) = BCData(mm)%rFace(i, j) &
                                                     + coefTimeALE(l) * BCData(mm)%rFaceALE(0, i, j)
                        end do
                    end do
                end if testAssoc2
            end do updateRF

            updateUS: do mm = 1, nViscBocos
                do j = BCData(mm)%jcBeg, BCData(mm)%jcEnd
                    do i = BCData(mm)%icBeg, BCData(mm)%icEnd
                        BCData(mm)%uSlip(i, j, 1) = BCData(mm)%uSlip(i, j, 1) &
                                                    + coefTimeALE(l) * BCData(mm)%uSlipALE(l, i, j, 1)
                        BCData(mm)%uSlip(i, j, 2) = BCData(mm)%uSlip(i, j, 2) &
                                                    + coefTimeALE(l) * BCData(mm)%uSlipALE(l, i, j, 2)
                        BCData(mm)%uSlip(i, j, 3) = BCData(mm)%uSlip(i, j, 3) &
                                                    + coefTimeALE(l) * BCData(mm)%uSlipALE(l, i, j, 3)
                    end do
                end do
            end do updateUS
        end do ALEloop

    end subroutine interpLevelALEBC_block

! ===========================================================
    subroutine recoverLevelALEBC_block
        !
        !      ******************************************************************
        !      *                                                                *
        !      * recoverLevelALEBC_block recovers current geometric data on     *
        !      * boundary from temporary interpolation                          *
        !      *                                                                *
        !      ******************************************************************
        !
        use blockPointers
        use inputUnsteady
        use inputPhysics
        implicit none
        !
        !      Local variables.
        !
        integer(kind=intType) :: i, j, k, nn, mm, kk

        if (.not. useALE .or. equationMode .ne. unsteady) then
            return
        end if

        recoverNM: do mm = 1, nBocos
            do j = BCData(mm)%jcBeg, BCData(mm)%jcEnd
                do i = BCData(mm)%icBeg, BCData(mm)%icEnd
                    BCData(mm)%norm(i, j, 1) = BCData(mm)%normALE(0, i, j, 1)
                    BCData(mm)%norm(i, j, 2) = BCData(mm)%normALE(0, i, j, 2)
                    BCData(mm)%norm(i, j, 3) = BCData(mm)%normALE(0, i, j, 3)
                end do
            end do
        end do recoverNM

        recoverRF: do mm = 1, nBocos
            testAssoc: if (associated(BCData(mm)%rFace)) then
                do j = BCData(mm)%jcBeg, BCData(mm)%jcEnd
                    do i = BCData(mm)%icBeg, BCData(mm)%icEnd
                        BCData(mm)%rFace(i, j) = BCData(mm)%rFaceALE(0, i, j)
                    end do
                end do
            end if testAssoc
        end do recoverRF

        recoverUS: do mm = 1, nViscBocos
            do j = BCData(mm)%jcBeg, BCData(mm)%jcEnd
                do i = BCData(mm)%icBeg, BCData(mm)%icEnd
                    BCData(mm)%uSlip(i, j, 1) = BCData(mm)%uSlipALE(0, i, j, 1)
                    BCData(mm)%uSlip(i, j, 2) = BCData(mm)%uSlipALE(0, i, j, 2)
                    BCData(mm)%uSlip(i, j, 3) = BCData(mm)%uSlipALE(0, i, j, 3)
                end do
            end do
        end do recoverUS

    end subroutine recoverLevelALEBC_block

    ! ----------------------------------------------------------------------
    !                                                                      |
    !                    No Tapenade Routine below this line               |
    !                                                                      |
    ! ----------------------------------------------------------------------

end module ALEUtils
