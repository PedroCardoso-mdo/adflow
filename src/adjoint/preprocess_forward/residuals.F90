










module residuals
contains

    subroutine residual_block
        !
        !       residual computes the residual of the mean flow equations on
        !       the current MG level.
        !
        use blockPointers
        use cgnsGrid
        use flowVarRefState
        use inputIteration
        use inputDiscretization
        use inputTimeSpectral
        use inputUnsteady ! Added by HDN
        use iteration
        use inputAdjoint
        use flowUtils, only: computeSpeedOfSoundSquared, allNodalGradients
        use fluxes
        implicit none
        !
        !      Local variables.
        !
        integer(kind=intType) :: discr
        integer(kind=intType) :: i, j, k, l
        integer(kind=intType) :: iale, jale, kale, lale, male ! For loops of ALE
        real(kind=realType), parameter :: K1 = 1.05_realType
        ! The line below is only used for the low-speed preconditioner part of this routine
        real(kind=realType), parameter :: K2 = 0.6_realType ! Random given number

        real(kind=realType), parameter :: M0 = 0.2_realType ! Mach number preconditioner activation
        real(kind=realType), parameter :: alpha = 0_realType
        real(kind=realType), parameter :: delta = 0_realType
        !real(kind=realType), parameter :: hinf = 2_realType ! Test phase
        real(kind=realType), parameter :: Cpres = 4.18_realType ! Test phase
        real(kind=realType), parameter :: TEMP = 297.15_realType

        !
        !     Local variables
        !
        real(kind=realType) :: K3, h, velXrho, velYrho, velZrho, SoS, hinf
        real(kind=realType) :: resM, A11, A12, A13, A14, A15, A21, A22, A23, A24, A25, A31, A32, A33, A34, A35
        real(kind=realType) :: A41, A42, A43, A44, A45, A51, A52, A53, A54, A55, B11, B12, B13, B14, B15
        real(kind=realType) :: B21, B22, B23, B24, B25, B31, B32, B33, B34, B35
        real(kind=realType) :: B41, B42, B43, B44, B45, B51, B52, B53, B54, B55
        real(kind=realType) :: rhoHdash, betaMr2
        real(kind=realType) :: G, q
        real(kind=realType) :: b1, b2, b3, b4, b5
        real(kind=realType) :: dwo(nwf)
        logical :: fineGrid

        ! Set the value of rFil, which controls the fraction of the old
        ! dissipation residual to be used. This is only for the runge-kutta
        ! schemes; for other smoothers rFil is simply set to 1.0.
        ! Note the index rkStage+1 for cdisRK. The reason is that the
        ! residual computation is performed before rkStage is incremented.

        if (smoother == RungeKutta) then
            rFil = cdisRK(rkStage + 1)
        else
            rFil = one
        end if

        ! Set the value of the discretization, depending on the grid level,
        ! and the logical fineGrid, which indicates whether or not this
        ! is the finest grid level of the current mg cycle.

        discr = spaceDiscrCoarse
        if (currentLevel == 1) discr = spaceDiscr

        fineGrid = .false.
        if (currentLevel == groundLevel) fineGrid = .true.

        ! ===========================================================
        !
        ! Assuming ALE has nothing to do with MG
        ! The geometric data will be interpolated if in MD mode
        !
        ! ===========================================================
        ! ===========================================================
        !
        ! The fluxes are calculated as usual
        !
        ! ===========================================================

        call inviscidCentralFlux

        select case (discr)

        case (dissScalar) ! Standard scalar dissipation scheme.

            if (fineGrid) then
                if (.not. lumpedDiss) then
                    call inviscidDissFluxScalar
                else
                    call inviscidDissFluxScalarApprox
                end if
            else
            end if

            !===========================================================

        case (dissMatrix) ! Matrix dissipation scheme.

            if (fineGrid) then
                if (.not. lumpedDiss) then
                    call inviscidDissFluxMatrix
                else
                    call inviscidDissFluxMatrixApprox
                end if
            else
            end if

            !===========================================================

        case (upwind) ! Dissipation via an upwind scheme.

            call inviscidUpwindFlux(fineGrid)

        end select

        !-------------------------------------------------------
        ! Lastly, recover the old s[I,J,K], sFace[I,J,K]
        ! This shall be done before difussive and source terms
        ! are computed.
        !-------------------------------------------------------

        if (viscous) then
            ! Only compute viscous fluxes if rFil > 0
            if (abs(rFil) > thresholdReal) then
                ! not lumpedDiss means it isn't the PC...call the vicousFlux
                if (.not. lumpedDiss) then
                    call computeSpeedOfSoundSquared
                    call allNodalGradients
                    call viscousFlux
                else
                    ! This is a PC calc...only include viscous fluxes if viscPC
                    ! is used
                    ! if full visc is true, also need full viscous terms, even if
                    ! lumpedDiss is true
                    call computeSpeedOfSoundSquared
                    if (viscPC) then
                        call allNodalGradients
                        call viscousFlux
                    else
                        call viscousFluxApprox
                    end if
                end if
            end if
        end if

        !===========================================================

        ! Add the dissipative and possibly viscous fluxes to the
        ! Euler fluxes. Loop over the owned cells and add fw to dw.
        ! Also multiply by iblank so that no updates occur in holes
        if (lowspeedpreconditioner) then
            do k = 2, kl
                do j = 2, jl
                    do i = 2, il
                        !    Compute speed of sound
                        SoS = sqrt(gamma(i, j, k) * p(i, j, k) / w(i, j, k, irho))

                        ! Compute velocities without rho from state vector
                        !      (w is pointer.. see type blockType setup in block.F90)
                        !      w(0:ib,0:jb,0:kb,1:nw) is allocated in block.F90
                        !      these are per definition nw=[rho,u,v,w,rhoeE]
                        !      so the velocity is simply just taken out below...
                        !      we do not have to divide with rho since it is already
                        !      without rho...
                        velXrho = w(i, j, k, ivx) ! ivx: l. 60 in constants.F90
                        velYrho = w(i, j, k, ivy)
                        velZrho = w(i, j, k, ivz)

                        q = (velXrho**2 + velYrho**2 + velZrho**2)

                        resM = sqrt(q) / SoS
                        ! resM above is used as M_a (thesis) and M (paper 2015)
                        ! and is the Free Stream Mach number

                        ! see routine setup above:
                        ! l. 30: real(kind=realType), parameter :: K1 = 1.05_realType
                        ! Random given number for K2:
                        ! l. 31: real(kind=realType), parameter :: K2 = 0.6_realType
                        ! Mach number preconditioner activation for K3:
                        ! l. 32: real(kind=realType), parameter :: M0 = 0.2_realType
                        !
                        !    Compute K3
                        ! eq. 2.7 in Garg 2015. K1, M0 and resM are scalars
                        !
                        ! unfortunately, Garg has switched the K1 and K3 here in the
                        ! code. In both paper and thesis it is K3 that is used to det-
                        ! ermine K1 below

                        !
                        !    Compute K3

                        K3 = K1 * (1 + ((1 - K1 * M0**2) * resM**2) / (K1 * M0**4))
                        !    Compute BetaMr2
                        ! betaMr2 -> eq. 7 in Garg 2015
                        ! (use eq. 2.6 in thesis thesis since paper has an error)
                        ! where a==SoS
                        !
                        ! again, K1 and K3 are switched compared with paper/thesis
                        !    Compute BetaMr2
                        betaMr2 = min(max(K3 * (velXrho**2 + velYrho**2 &
                                                + velZrho**2), ((K2) * (wInf(ivx)**2 &
                                                                        + wInf(ivy)**2 + wInf(ivz)**2))), SoS**2)

                        ! above, the wInf is the free stream velocity
                        !
                        ! Should this first line's first element have SoS^4 or SoS^2

                        A11 = (betaMr2) * (1 / SoS**4)
                        A12 = zero
                        A13 = zero
                        A14 = zero
                        A15 = (-betaMr2) / SoS**4

                        A21 = one * velXrho / SoS**2
                        A22 = one * w(i, j, k, irho)
                        A23 = zero
                        A24 = zero
                        A25 = one * (-velXrho) / SoS**2

                        A31 = one * velYrho / SoS**2
                        A32 = zero
                        A33 = one * w(i, j, k, irho)
                        A34 = zero
                        A35 = one * (-velYrho) / SoS**2

                        A41 = one * velZrho / SoS**2
                        A42 = zero
                        A43 = zero
                        A44 = one * w(i, j, k, irho)
                        A45 = zero + one * (-velZrho) / SoS**2

                        ! mham: seems he fixed the above line an irregular way?

                        A51 = one * ((1 / (gamma(i, j, k) - 1)) + (resM**2) / 2)
                        A52 = one * w(i, j, k, irho) * velXrho
                        A53 = one * w(i, j, k, irho) * velYrho
                        A54 = one * w(i, j, k, irho) * velzrho
                        A55 = one * ((-(resM**2)) / 2)

                        B11 = A11 * (gamma(i, j, k) - 1) * q / 2 + A12 * (-velXrho) &
                              / w(i, j, k, irho) + A13 * (-velYrho) / w(i, j, k, irho) + &
                              A14 * (-velZrho) / w(i, j, k, irho) &
                              + A15 * (((gamma(i, j, k) - 1) * q / 2) - SoS**2)
                        B12 = A11 * (1 - gamma(i, j, k)) * velXrho + A12 * 1 / w(i, j, k, irho) &
                              + A15 * (1 - gamma(i, j, k)) * velXrho
                        B13 = A11 * (1 - gamma(i, j, k)) * velYrho + A13 &
                              / w(i, j, k, irho) + A15 * (1 - gamma(i, j, k)) * velYrho
                        B14 = A11 * (1 - gamma(i, j, k)) * velZrho &
                              + A14 / w(i, j, k, irho) + A15 * (1 - gamma(i, j, k)) * velZrho
                        B15 = A11 * (gamma(i, j, k) - 1) + A15 * (gamma(i, j, k) - 1)

                        B21 = A21 * (gamma(i, j, k) - 1) * q / 2 + A22 * (-velXrho) &
                              / w(i, j, k, irho) + A23 * (-velYrho) / w(i, j, k, irho) + A24 * (-velZrho) &
                              / w(i, j, k, irho) + A25 * (((gamma(i, j, k) - 1) * q / 2) - SoS**2)
                        B22 = A21 * (1 - gamma(i, j, k)) * velXrho + A22 &
                              / w(i, j, k, irho) + A25 * (1 - gamma(i, j, k)) * velXrho
                        B23 = A21 * (1 - gamma(i, j, k)) * velYrho &
                              + A23 * 1 / w(i, j, k, irho) + A25 * (1 - gamma(i, j, k)) * velYrho
                        B24 = A21 * (1 - gamma(i, j, k)) * velZrho &
                              + A24 * 1 / w(i, j, k, irho) + A25 * (1 - gamma(i, j, k)) * velZrho
                        B25 = A21 * (gamma(i, j, k) - 1) + A25 * (gamma(i, j, k) - 1)

                        B31 = A31 * (gamma(i, j, k) - 1) * q / 2 + A32 * (-velXrho) &
                              / w(i, j, k, irho) + A33 * (-velYrho) / w(i, j, k, irho) + &
                              A34 * (-velZrho) / w(i, j, k, irho) &
                              + A35 * (((gamma(i, j, k) - 1) * q / 2) - SoS**2)
                        B32 = A31 * (1 - gamma(i, j, k)) * velXrho + A32 &
                              / w(i, j, k, irho) + A35 * (1 - gamma(i, j, k)) * velXrho
                        B33 = A31 * (1 - gamma(i, j, k)) * velYrho &
                              + A33 * 1 / w(i, j, k, irho) + A35 * (1 - gamma(i, j, k)) * velYrho
                        B34 = A31 * (1 - gamma(i, j, k)) * velZrho &
                              + A34 * 1 / w(i, j, k, irho) + A35 * (1 - gamma(i, j, k)) * velZrho
                        B35 = A31 * (gamma(i, j, k) - 1) + A35 * (gamma(i, j, k) - 1)

                        B41 = A41 * (gamma(i, j, k) - 1) * q / 2 + A42 * (-velXrho) &
                              / w(i, j, k, irho) + A43 * (-velYrho) / w(i, j, k, irho) + A44 * (-velZrho) &
                              / w(i, j, k, irho) + A45 * (((gamma(i, j, k) - 1) * q / 2) - SoS**2)
                        B42 = A41 * (1 - gamma(i, j, k)) * velXrho + A42 &
                              / w(i, j, k, irho) + A45 * (1 - gamma(i, j, k)) * velXrho
                        B43 = A41 * (1 - gamma(i, j, k)) * velYrho &
                              + A43 * 1 / w(i, j, k, irho) + A45 * (1 - gamma(i, j, k)) * velYrho
                        B44 = A41 * (1 - gamma(i, j, k)) * velZrho &
                              + A44 * 1 / w(i, j, k, irho) + A45 * (1 - gamma(i, j, k)) * velZrho
                        B45 = A41 * (gamma(i, j, k) - 1) + A45 * (gamma(i, j, k) - 1)

                        B51 = A51 * (gamma(i, j, k) - 1) * q / 2 + A52 * (-velXrho) &
                              / w(i, j, k, irho) + A53 * (-velYrho) / w(i, j, k, irho) + A54 * (-velZrho) &
                              / w(i, j, k, irho) + A55 * (((gamma(i, j, k) - 1) * q / 2) - SoS**2)
                        B52 = A51 * (1 - gamma(i, j, k)) * velXrho + A52 &
                              / w(i, j, k, irho) + A55 * (1 - gamma(i, j, k)) * velXrho
                        B53 = A51 * (1 - gamma(i, j, k)) * velYrho &
                              + A53 * 1 / w(i, j, k, irho) + A55 * (1 - gamma(i, j, k)) * velYrho
                        B54 = A51 * (1 - gamma(i, j, k)) * velZrho &
                              + A54 * 1 / w(i, j, k, irho) + A55 * (1 - gamma(i, j, k)) * velZrho
                        B55 = A51 * (gamma(i, j, k) - 1) + A55 * (gamma(i, j, k) - 1)

                        ! dwo is the orginal redisual
                        do l = 1, nwf
                            dwo(l) = (dw(i, j, k, l) + fw(i, j, k, l)) * max(real(iblank(i, j, k), realType), zero)
                        end do

                        dw(i, j, k, 1) = B11 * dwo(1) + B12 * dwo(2) + B13 * dwo(3) + B14 * dwo(4) + B15 * dwo(5)
                        dw(i, j, k, 2) = B21 * dwo(1) + B22 * dwo(2) + B23 * dwo(3) + B24 * dwo(4) + B25 * dwo(5)
                        dw(i, j, k, 3) = B31 * dwo(1) + B32 * dwo(2) + B33 * dwo(3) + B34 * dwo(4) + B35 * dwo(5)
                        dw(i, j, k, 4) = B41 * dwo(1) + B42 * dwo(2) + B43 * dwo(3) + B44 * dwo(4) + B45 * dwo(5)
                        dw(i, j, k, 5) = B51 * dwo(1) + B52 * dwo(2) + B53 * dwo(3) + B54 * dwo(4) + B55 * dwo(5)

                    end do
                end do
            end do ! end of lowspeedpreconditioners three cells loops

        else ! else.. i.e. if we do not have preconditioner turned on...
            do l = 1, nwf
                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il
                            dw(i, j, k, l) = (dw(i, j, k, l) + fw(i, j, k, l)) &
                                             * max(real(iblank(i, j, k), realType), zero)
                        end do
                    end do
                end do
            end do
        end if

    end subroutine residual_block

    subroutine sourceTerms_block(nn, res, iRegion, pLocal)

        ! Apply the source terms for the given block. Assume that the
        ! block pointers are already set.
        use constants
        use actuatorRegionData
        use blockPointers, only: vol, dw, w
        use flowVarRefState, only: pRef, uRef, LRef
        use communication
        use iteration, only: ordersConverged
        implicit none

        ! Input
        integer(kind=intType), intent(in) :: nn, iRegion
        logical, intent(in) :: res
        real(kind=realType), intent(inout) :: pLocal

        ! Working
        integer(kind=intType) :: i, j, k, ii, iStart, iEnd
        real(kind=realType) :: Ftmp(3), Vx, Vy, Vz, F_fact(3), Q_fact, Qtmp, reDim, factor, oStart, oEnd

        reDim = pRef * uRef

        ! Compute the relaxation factor based on the ordersConverged

        ! How far we are into the ramp:
        if (ordersConverged < actuatorRegions(iRegion)%relaxStart) then
            factor = zero
        else if (ordersConverged > actuatorRegions(iRegion)%relaxEnd) then
            factor = one
        else ! In between
            oStart = actuatorRegions(iRegion)%relaxStart
            oEnd = actuatorRegions(iRegion)%relaxEnd
            factor = (ordersConverged - oStart) / (oEnd - oStart)
        end if

        ! Compute the constant force factor
        F_fact = factor * actuatorRegions(iRegion)%force / actuatorRegions(iRegion)%volume / pRef

        ! Heat factor. This is heat added per unit volume per unit time
        Q_fact = factor * actuatorRegions(iRegion)%heat / actuatorRegions(iRegion)%volume / (pRef * uRef * LRef * LRef)

        ! Loop over the ranges for this block
        iStart = actuatorRegions(iRegion)%blkPtr(nn - 1) + 1
        iEnd = actuatorRegions(iRegion)%blkPtr(nn)

        !$AD II-LOOP
        do ii = iStart, iEnd

            ! Extract the cell ID.
            i = actuatorRegions(iRegion)%cellIDs(1, ii)
            j = actuatorRegions(iRegion)%cellIDs(2, ii)
            k = actuatorRegions(iRegion)%cellIDs(3, ii)

            ! This actually gets the force
            FTmp = vol(i, j, k) * F_fact

            Vx = w(i, j, k, iVx)
            Vy = w(i, j, k, iVy)
            Vz = w(i, j, k, iVz)

            ! this gets the heat addition rate
            QTmp = vol(i, j, k) * Q_fact

            if (res) then
                ! Momentum residuals
                dw(i, j, k, imx:imz) = dw(i, j, k, imx:imz) - Ftmp

                ! energy residuals
                dw(i, j, k, iRhoE) = dw(i, j, k, iRhoE) - &
                                     Ftmp(1) * Vx - Ftmp(2) * Vy - Ftmp(3) * Vz - Qtmp
            else
                ! Add in the local power contribution:
                pLocal = pLocal + (Vx * Ftmp(1) + Vy * FTmp(2) + Vz * Ftmp(3)) * reDim
            end if
        end do

    end subroutine sourceTerms_block

    subroutine initres_block(varStart, varEnd, nn, sps)
        !
        !       initres initializes the given range of the residual. Either to
        !       zero, steady computation, or to an unsteady term for the time
        !       spectral and unsteady modes. For the coarser grid levels the
        !       residual forcing term is taken into account.
        !
        use blockPointers
        use flowVarRefState
        use inputIteration
        use inputPhysics
        use inputTimeSpectral
        use inputUnsteady
        use iteration

        implicit none
        !
        !      Subroutine arguments.
        !
        integer(kind=intType), intent(in) :: varStart, varEnd, nn, sps
        !
        !      Local variables.
        !
        integer(kind=intType) :: mm, ll, ii, jj, i, j, k, l, m
        real(kind=realType) :: oneOverDt, tmp

        real(kind=realType), dimension(:, :, :, :), pointer :: ww, wsp, wsp1
        real(kind=realType), dimension(:, :, :), pointer :: volsp

        ! Return immediately of no variables are in the range.

        if (varEnd < varStart) return

        ! Determine the equation mode and act accordingly.

        select case (equationMode)
        case (steady)

            ! Steady state computation.
            ! Determine the currently active multigrid level.

            steadyLevelTest: if (currentLevel == groundLevel) then

                ! Ground level of the multigrid cycle. Initialize the
                ! owned residuals to zero.

                do l = varStart, varEnd
                    do k = 2, kl
                        do j = 2, jl
                            do i = 2, il
                                dw(i, j, k, l) = zero
                            end do
                        end do
                    end do
                end do

                else steadyLevelTest

                ! Coarse grid level. Initialize the owned cells to the
                ! residual forcing terms.

                do l = varStart, varEnd
                    do k = 2, kl
                        do j = 2, jl
                            do i = 2, il
                                dw(i, j, k, l) = wr(i, j, k, l)
                            end do
                        end do
                    end do
                end do
            end if steadyLevelTest

        end select

        ! Set the residual in the halo cells to zero. This is just
        ! to avoid possible problems. Their values do not matter.

        do l = varStart, varEnd
            do k = 0, kb
                do j = 0, jb
                    dw(0, j, k, l) = zero
                    dw(1, j, k, l) = zero
                    dw(ie, j, k, l) = zero
                    dw(ib, j, k, l) = zero
                end do
            end do

            do k = 0, kb
                do i = 2, il
                    dw(i, 0, k, l) = zero
                    dw(i, 1, k, l) = zero
                    dw(i, je, k, l) = zero
                    dw(i, jb, k, l) = zero
                end do
            end do

            do j = 2, jl
                do i = 2, il
                    dw(i, j, 0, l) = zero
                    dw(i, j, 1, l) = zero
                    dw(i, j, ke, l) = zero
                    dw(i, j, kb, l) = zero
                end do
            end do
        end do

    end subroutine initres_block

    ! ----------------------------------------------------------------------
    !                                                                      |
    !                    No Tapenade Routine below this line               |
    !                                                                      |
    ! ----------------------------------------------------------------------

end module residuals
