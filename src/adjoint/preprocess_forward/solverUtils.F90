










module solverUtils
contains
    subroutine timeStep_block(onlyRadii)
        !
        !       timeStep computes the time step, or more precisely the time
        !       step divided by the volume per unit CFL, in the owned cells.
        !       However, for the artificial dissipation schemes, the spectral
        !       radIi in the halo's are needed. Therefore the loop is taken
        !       over the the first level of halo cells. The spectral radIi are
        !       stored and possibly modified for high aspect ratio cells.
        !
        use constants
        use blockPointers, only: ie, je, ke, il, jl, kl, w, p, rlv, rev, &
                                 radi, radj, radk, si, sj, sk, sFaceI, sfaceJ, sfaceK, dtl, gamma, vol, &
                                 addGridVelocities, sectionID
        use flowVarRefState, only: timeRef, eddyModel, gammaInf, pInfCorr, &
                                   viscous, rhoInf
        use inputDiscretization, only: adis, dirScaling, radiiNeededCoarse, &
                                       radiiNeededFine, precond, acousticScaleFactor
        use inputPhysics, only: equationMode
        use iteration, only: groundLevel, currentLevel
        use section, only: sections
        use inputTimeSpectral, only: nTimeIntervalsSpectral
        use utils, only: terminate
        implicit none
        !
        !      Subroutine argument.
        !
        logical, intent(in) :: onlyRadii
        !
        !      Local parameters.
        !
        real(kind=realType), parameter :: b = 2.0_realType
        !
        !      Local variables.
        !
        integer(kind=intType) :: i, j, k, ii

        real(kind=realType) :: plim, rlim, clim2
        real(kind=realType) :: uux, uuy, uuz, cc2, qsi, qsj, qsk, sx, sy, sz, rmu
        real(kind=realType) :: ri, rj, rk, rij, rjk, rki
        real(kind=realType) :: vsi, vsj, vsk, rfl, dpi, dpj, dpk
        real(kind=realType) :: sFace, tmp

        logical :: radiiNeeded, doScaling

        ! Determine whether or not the spectral radii are needed for the
        ! flux computation.

        radiiNeeded = radiiNeededCoarse
        if (currentLevel <= groundLevel) radiiNeeded = radiiNeededFine

        ! Return immediately if only the spectral radii must be computed
        ! and these are not needed for the flux computation.

        if (onlyRadii .and. (.not. radiiNeeded)) return

        ! Set the value of plim. To be fully consistent this must have
        ! the dimension of a pressure. Therefore a fraction of pInfCorr
        ! is used. Idem for rlim; compute clim2 as well.

        plim = 0.001_realType * pInfCorr
        rlim = 0.001_realType * rhoInf
        clim2 = 0.000001_realType * gammaInf * pInfCorr / rhoInf

        doScaling = (dirScaling .and. currentLevel <= groundLevel)

        ! Initialize sFace to zero. This value will be used if the
        ! block is not moving.

        sFace = zero
        !
        !           Inviscid contribution, depending on the preconditioner.
        !           Compute the cell centered values of the spectral radii.
        !
        select case (precond)

        case (noPrecond)

            ! No preconditioner. Simply the standard spectral radius.
            ! Loop over the cells, including the first level halo.

                do k = 1, ke
                    do j = 1, je
                        do i = 1, ie
                            ! Compute the velocities and speed of sound squared.

                            uux = w(i, j, k, ivx)
                            uuy = w(i, j, k, ivy)
                            uuz = w(i, j, k, ivz)
                            cc2 = gamma(i, j, k) * p(i, j, k) / w(i, j, k, irho)
                            cc2 = max(cc2, clim2)

                            ! Set the dot product of the grid velocity and the
                            ! normal in i-direction for a moving face. To avoid
                            ! a number of multiplications by 0.5 simply the sum
                            ! is taken.

                            if (addGridVelocities) &
                                sFace = sFaceI(i - 1, j, k) + sFaceI(i, j, k)

                            ! Spectral radius in i-direction.

                            sx = si(i - 1, j, k, 1) + si(i, j, k, 1)
                            sy = si(i - 1, j, k, 2) + si(i, j, k, 2)
                            sz = si(i - 1, j, k, 3) + si(i, j, k, 3)

                            qsi = uux * sx + uuy * sy + uuz * sz - sFace

                            ri = half * (abs(qsi) &
                                         + acousticScaleFactor * sqrt(cc2 * (sx**2 + sy**2 + sz**2)))

                            ! The grid velocity in j-direction.

                            if (addGridVelocities) &
                                sFace = sFaceJ(i, j - 1, k) + sFaceJ(i, j, k)

                            ! Spectral radius in j-direction.

                            sx = sj(i, j - 1, k, 1) + sj(i, j, k, 1)
                            sy = sj(i, j - 1, k, 2) + sj(i, j, k, 2)
                            sz = sj(i, j - 1, k, 3) + sj(i, j, k, 3)

                            qsj = uux * sx + uuy * sy + uuz * sz - sFace

                            rj = half * (abs(qsj) &
                                         + acousticScaleFactor * sqrt(cc2 * (sx**2 + sy**2 + sz**2)))

                            ! The grid velocity in k-direction.

                            if (addGridVelocities) &
                                sFace = sFaceK(i, j, k - 1) + sFaceK(i, j, k)

                            ! Spectral radius in k-direction.

                            sx = sk(i, j, k - 1, 1) + sk(i, j, k, 1)
                            sy = sk(i, j, k - 1, 2) + sk(i, j, k, 2)
                            sz = sk(i, j, k - 1, 3) + sk(i, j, k, 3)

                            qsk = uux * sx + uuy * sy + uuz * sz - sFace

                            rk = half * (abs(qsk) &
                                         + acousticScaleFactor * sqrt(cc2 * (sx**2 + sy**2 + sz**2)))

                            ! Compute the inviscid contribution to the time step.

                            if (.not. onlyRadii) dtl(i, j, k) = ri + rj + rk

                            !
                            !           Adapt the spectral radii if directional scaling must be
                            !           applied.
                            !
                            if (doScaling) then

                                ! Avoid division by zero by clipping radi, radJ and
                                ! radK.

                                ri = max(ri, eps)
                                rj = max(rj, eps)
                                rk = max(rk, eps)

                                ! Compute the scaling in the three coordinate
                                ! directions.

                                rij = (ri / rj)**adis
                                rjk = (rj / rk)**adis
                                rki = (rk / ri)**adis

                                ! Create the scaled versions of the aspect ratios.
                                ! Note that the multiplication is done with radi, radJ
                                ! and radK, such that the influence of the clipping
                                ! is negligible.

                                radi(i, j, k) = ri * (one + one / rij + rki)
                                radJ(i, j, k) = rj * (one + one / rjk + rij)
                                radK(i, j, k) = rk * (one + one / rki + rjk)
                            else
                                radi(i, j, k) = ri
                                radj(i, j, k) = rj
                                radk(i, j, k) = rk
                            end if
                    end do
                end do
            end do

        case (Turkel)
            call terminate("timeStep", "Turkel preconditioner not implemented yet")

        case (ChoiMerkle)
            call terminate("timeStep", &
                           "choi merkle preconditioner not implemented yet")
        end select

        ! The rest of this file can be skipped if only the spectral
        ! radii need to be computed.
        testRadiiOnly: if (.not. onlyRadii) then

            ! The viscous contribution, if needed.

            viscousTerm: if (viscous) then

                ! Loop over the owned cell centers.

                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il

                            ! Compute the effective viscosity coefficient. The
                            ! factor 0.5 is a combination of two things. In the
                            ! standard central discretization of a second
                            ! derivative there is a factor 2 multiplying the
                            ! central node. However in the code below not the
                            ! average but the sum of the left and the right face
                            ! is taken and squared. This leads to a factor 4.
                            ! Combining both effects leads to 0.5. Furthermore,
                            ! it is divided by the volume and density to obtain
                            ! the correct dimensions and multiplied by the
                            ! non-dimensional factor factVis.

                            rmu = rlv(i, j, k)
                            if (eddyModel) rmu = rmu + rev(i, j, k)
                            rmu = half * rmu / (w(i, j, k, irho) * vol(i, j, k))

                            ! Add the viscous contribution in i-direction to the
                            ! (inverse) of the time step.

                            sx = si(i, j, k, 1) + si(i - 1, j, k, 1)
                            sy = si(i, j, k, 2) + si(i - 1, j, k, 2)
                            sz = si(i, j, k, 3) + si(i - 1, j, k, 3)

                            vsi = rmu * (sx * sx + sy * sy + sz * sz)
                            dtl(i, j, k) = dtl(i, j, k) + vsi

                            ! Add the viscous contribution in j-direction to the
                            ! (inverse) of the time step.

                            sx = sj(i, j, k, 1) + sj(i, j - 1, k, 1)
                            sy = sj(i, j, k, 2) + sj(i, j - 1, k, 2)
                            sz = sj(i, j, k, 3) + sj(i, j - 1, k, 3)

                            vsj = rmu * (sx * sx + sy * sy + sz * sz)
                            dtl(i, j, k) = dtl(i, j, k) + vsj

                            ! Add the viscous contribution in k-direction to the
                            ! (inverse) of the time step.

                            sx = sk(i, j, k, 1) + sk(i, j, k - 1, 1)
                            sy = sk(i, j, k, 2) + sk(i, j, k - 1, 2)
                            sz = sk(i, j, k, 3) + sk(i, j, k - 1, 3)

                            vsk = rmu * (sx * sx + sy * sy + sz * sz)
                            dtl(i, j, k) = dtl(i, j, k) + vsk

                        end do
                    end do
                end do

            end if viscousTerm

            ! For the spectral mode an additional term term must be
            ! taken into account, which corresponds to the contribution
            ! of the highest frequency.

            if (equationMode == timeSpectral) then

                tmp = nTimeIntervalsSpectral * pi * timeRef &
                      / sections(sectionID)%timePeriod

                ! Loop over the owned cell centers and add the term.

                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il
                            dtl(i, j, k) = dtl(i, j, k) + tmp * vol(i, j, k)
                        end do
                    end do
                end do

            end if

            ! Currently the inverse of dt/vol is stored in dtl. Invert
            ! this value such that the time step per unit cfl number is
            ! stored and correct in cases of high gradients.

            do k = 2, kl
                do j = 2, jl
                    do i = 2, il
                        dpi = abs(p(i + 1, j, k) - two * p(i, j, k) + p(i - 1, j, k)) &
                              / (p(i + 1, j, k) + two * p(i, j, k) + p(i - 1, j, k) + plim)
                        dpj = abs(p(i, j + 1, k) - two * p(i, j, k) + p(i, j - 1, k)) &
                              / (p(i, j + 1, k) + two * p(i, j, k) + p(i, j - 1, k) + plim)
                        dpk = abs(p(i, j, k + 1) - two * p(i, j, k) + p(i, j, k - 1)) &
                              / (p(i, j, k + 1) + two * p(i, j, k) + p(i, j, k - 1) + plim)
                        rfl = one / (one + b * (dpi + dpj + dpk))

                        dtl(i, j, k) = rfl / dtl(i, j, k)
                    end do
                end do
            end do

        end if testRadiiOnly

    end subroutine timeStep_block

    subroutine gridVelocitiesFineLevel_block(useOldCoor, t, sps, nn)
        !
        !       gridVelocitiesFineLevel computes the grid velocities for
        !       the cell centers and the normal grid velocities for the faces
        !       of moving blocks for the currently finest grid, i.e.
        !       groundLevel. The velocities are computed at time t for
        !       spectral mode sps. If useOldCoor is .true. the velocities
        !       are determined using the unsteady time integrator in
        !       combination with the old coordinates; otherwise the analytic
        !       form is used.
        !
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
        use flowUtils, only: derivativeRotMatrixRigid, getDirVector
        use utils, only: setCoefTimeIntegrator, tsAlpha, tsBeta, tsMach, terminate, &
                         rotMatrixRigidBody, getDirAngle

        implicit none
        !
        !      Subroutine arguments.
        !
        integer(kind=intType), intent(in) :: sps, nn
        logical, intent(in) :: useOldCoor

        real(kind=realType), dimension(*), intent(in) :: t
        !
        !      Local variables.
        !
        integer(kind=intType) :: mm
        integer(kind=intType) :: i, j, k, ii, iie, jje, kke

        real(kind=realType) :: oneOver4dt, oneOver8dt
        real(kind=realType) :: velxGrid, velyGrid, velzGrid, ainf
        real(kind=realType) :: velxGrid0, velyGrid0, velzGrid0

        real(kind=realType), dimension(3) :: sc, xc, xxc
        real(kind=realType), dimension(3) :: rotCenter, rotRate

        real(kind=realType), dimension(3) :: rotationPoint
        real(kind=realType), dimension(3, 3) :: rotationMatrix, &
                                                derivRotationMatrix

        real(kind=realType) :: tNew, tOld
        real(kind=realType), dimension(:, :), pointer :: sFace

        real(kind=realType), dimension(:, :, :), pointer :: xx, ss
        real(kind=realType), dimension(:, :, :, :), pointer :: xxOld

        real(kind=realType) :: intervalMach, alphaTS, alphaIncrement, &
                               betaTS, betaIncrement
        real(kind=realType), dimension(3) :: velDir
        real(kind=realType), dimension(3) :: refDirection

        ! Compute the mesh velocity from the given mesh Mach number.

        ! vel{x,y,z}Grid0 is the ACTUAL velocity you want at the
        ! geometry.
        aInf = sqrt(gammaInf * pInf / rhoInf)
        velxGrid0 = (aInf * machgrid) * (-velDirFreestream(1))
        velyGrid0 = (aInf * machgrid) * (-velDirFreestream(2))
        velzGrid0 = (aInf * machgrid) * (-velDirFreestream(3))

        ! Compute the derivative of the rotation matrix and the rotation
        ! point; needed for velocity due to the rigid body rotation of
        ! the entire grid. It is assumed that the rigid body motion of
        ! the grid is only specified if there is only 1 section present.

        call derivativeRotMatrixRigid(derivRotationMatrix, rotationPoint, t(1))

        !compute the rotation matrix to update the velocities for the time
        !spectral stability derivative case...


        testMoving: if (blockIsMoving) then
            ! Determine the situation we are having here.

            testUseOldCoor: if (useOldCoor) then
                else testUseOldCoor
                !
                !             The velocities must be determined analytically.
                !
                ! Store the rotation center and determine the
                ! nonDimensional rotation rate of this block. As the
                ! reference length is 1 timeRef == 1/uRef and at the end
                ! the nonDimensional velocity is computed.

                j = nbkGlobal

                rotCenter = cgnsDoms(j)%rotCenter
                rotRate = timeRef * cgnsDoms(j)%rotRate

                velXgrid = velXGrid0
                velYgrid = velYGrid0
                velZgrid = velZGrid0
                !
                !             Grid velocities of the cell centers, including the
                !             1st level halo cells.
                !
                ! Loop over the cells, including the 1st level halo's.

                do k = 1, ke
                    do j = 1, je
                        do i = 1, ie

                            ! Determine the coordinates of the cell center,
                            ! which are stored in xc.

                            xc(1) = eighth * (flowDoms(nn, groundLevel, sps)%x(i - 1, j - 1, k - 1, 1) &
                                              + flowDoms(nn, groundLevel, sps)%x(i, j - 1, k - 1, 1) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, j, k - 1, 1) &
                                              + flowDoms(nn, groundLevel, sps)%x(i, j, k - 1, 1) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, j - 1, k, 1) &
                                              + flowDoms(nn, groundLevel, sps)%x(i, j - 1, k, 1) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, j, k, 1) &
                                              + flowDoms(nn, groundLevel, sps)%x(i, j, k, 1))
                            xc(2) = eighth * (flowDoms(nn, groundLevel, sps)%x(i - 1, j - 1, k - 1, 2) &
                                              + flowDoms(nn, groundLevel, sps)%x(i, j - 1, k - 1, 2) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, j, k - 1, 2) &
                                              + flowDoms(nn, groundLevel, sps)%x(i, j, k - 1, 2) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, j - 1, k, 2) &
                                              + flowDoms(nn, groundLevel, sps)%x(i, j - 1, k, 2) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, j, k, 2) &
                                              + flowDoms(nn, groundLevel, sps)%x(i, j, k, 2))
                            xc(3) = eighth * (flowDoms(nn, groundLevel, sps)%x(i - 1, j - 1, k - 1, 3) &
                                              + flowDoms(nn, groundLevel, sps)%x(i, j - 1, k - 1, 3) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, j, k - 1, 3) &
                                              + flowDoms(nn, groundLevel, sps)%x(i, j, k - 1, 3) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, j - 1, k, 3) &
                                              + flowDoms(nn, groundLevel, sps)%x(i, j - 1, k, 3) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, j, k, 3) &
                                              + flowDoms(nn, groundLevel, sps)%x(i, j, k, 3))

                            ! Determine the coordinates relative to the
                            ! center of rotation.

                            xxc(1) = xc(1) - rotCenter(1)
                            xxc(2) = xc(2) - rotCenter(2)
                            xxc(3) = xc(3) - rotCenter(3)

                            ! Determine the rotation speed of the cell center,
                            ! which is omega*r.

                            sc(1) = rotRate(2) * xxc(3) - rotRate(3) * xxc(2)
                            sc(2) = rotRate(3) * xxc(1) - rotRate(1) * xxc(3)
                            sc(3) = rotRate(1) * xxc(2) - rotRate(2) * xxc(1)

                            ! Determine the coordinates relative to the
                            ! rigid body rotation point.

                            xxc(1) = xc(1) - rotationPoint(1)
                            xxc(2) = xc(2) - rotationPoint(2)
                            xxc(3) = xc(3) - rotationPoint(3)

                            ! Determine the total velocity of the cell center.
                            ! This is a combination of rotation speed of this
                            ! block and the entire rigid body rotation.

                            s(i, j, k, 1) = sc(1) + velxGrid &
                                            + derivRotationMatrix(1, 1) * xxc(1) &
                                            + derivRotationMatrix(1, 2) * xxc(2) &
                                            + derivRotationMatrix(1, 3) * xxc(3)
                            s(i, j, k, 2) = sc(2) + velyGrid &
                                            + derivRotationMatrix(2, 1) * xxc(1) &
                                            + derivRotationMatrix(2, 2) * xxc(2) &
                                            + derivRotationMatrix(2, 3) * xxc(3)
                            s(i, j, k, 3) = sc(3) + velzGrid &
                                            + derivRotationMatrix(3, 1) * xxc(1) &
                                            + derivRotationMatrix(3, 2) * xxc(2) &
                                            + derivRotationMatrix(3, 3) * xxc(3)
                        end do
                    end do
                end do
                !
                !             Normal grid velocities of the faces.
                !
                ! Loop over the three directions.

                ! The original code is elegant but the Tapenade has a difficult time
                ! to understand it. Thus, we unfold it and make it easier for the
                ! Tapenade.

                ! i-direction
                do k = 1, ke
                    do j = 1, je
                        do i = 0, ie

                            ! Determine the coordinates of the face center,
                            ! which are stored in xc.

                            xc(1) = fourth * (flowDoms(nn, groundLevel, sps)%x(i, j - 1, k - 1, 1) + &
                                              flowDoms(nn, groundLevel, sps)%x(i, j, k - 1, 1) &
                                              + flowDoms(nn, groundLevel, sps)%x(i, j - 1, k, 1) + &
                                              flowDoms(nn, groundLevel, sps)%x(i, j, k, 1))
                            xc(2) = fourth * (flowDoms(nn, groundLevel, sps)%x(i, j - 1, k - 1, 2) + &
                                              flowDoms(nn, groundLevel, sps)%x(i, j, k - 1, 2) &
                                              + flowDoms(nn, groundLevel, sps)%x(i, j - 1, k, 2) + &
                                              flowDoms(nn, groundLevel, sps)%x(i, j, k, 2))
                            xc(3) = fourth * (flowDoms(nn, groundLevel, sps)%x(i, j - 1, k - 1, 3) + &
                                              flowDoms(nn, groundLevel, sps)%x(i, j, k - 1, 3) &
                                              + flowDoms(nn, groundLevel, sps)%x(i, j - 1, k, 3) + &
                                              flowDoms(nn, groundLevel, sps)%x(i, j, k, 3))

                            call cellFaceVelocities(xc, rotCenter, rotRate, &
                                                    velxGrid, velyGrid, velzGrid, derivRotationMatrix, sc)

                            ! Store the dot product of grid velocity sc and
                            ! the normal ss in sFace.

                            sFaceI(i, j, k) = sc(1) * si(i, j, k, 1) + sc(2) * si(i, j, k, 2) &
                                              + sc(3) * si(i, j, k, 3)
                        end do
                    end do
                end do

                ! j-direction
                do k = 1, ke
                    do j = 0, je
                        do i = 1, ie

                            ! Determine the coordinates of the face center,
                            ! which are stored in xc.

                            xc(1) = fourth * (flowDoms(nn, groundLevel, sps)%x(i - 1, j, k, 1) + &
                                              flowDoms(nn, groundLevel, sps)%x(i, j, k - 1, 1) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, j, k - 1, 1) + &
                                              flowDoms(nn, groundLevel, sps)%x(i, j, k, 1))
                            xc(2) = fourth * (flowDoms(nn, groundLevel, sps)%x(i - 1, j, k, 2) + &
                                              flowDoms(nn, groundLevel, sps)%x(i, j, k - 1, 2) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, j, k - 1, 2) + &
                                              flowDoms(nn, groundLevel, sps)%x(i, j, k, 2))
                            xc(3) = fourth * (flowDoms(nn, groundLevel, sps)%x(i - 1, j, k, 3) + &
                                              flowDoms(nn, groundLevel, sps)%x(i, j, k - 1, 3) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, j, k - 1, 3) + &
                                              flowDoms(nn, groundLevel, sps)%x(i, j, k, 3))

                            call cellFaceVelocities(xc, rotCenter, rotRate, &
                                                    velxGrid, velyGrid, velzGrid, derivRotationMatrix, sc)

                            ! Store the dot product of grid velocity sc and
                            ! the normal ss in sFace.

                            sFaceJ(i, j, k) = sc(1) * sj(i, j, k, 1) + sc(2) * sj(i, j, k, 2) &
                                              + sc(3) * sj(i, j, k, 3)
                        end do
                    end do
                end do

                ! k-direction
                do k = 0, ke
                    do j = 1, je
                        do i = 1, ie

                            ! Determine the coordinates of the face center,
                            ! which are stored in xc.

                            xc(1) = fourth * (flowDoms(nn, groundLevel, sps)%x(i, j - 1, k, 1) + &
                                              flowDoms(nn, groundLevel, sps)%x(i - 1, j, k, 1) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, j - 1, k, 1) + &
                                              flowDoms(nn, groundLevel, sps)%x(i, j, k, 1))
                            xc(2) = fourth * (flowDoms(nn, groundLevel, sps)%x(i, j - 1, k, 2) + &
                                              flowDoms(nn, groundLevel, sps)%x(i - 1, j, k, 2) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, j - 1, k, 2) + &
                                              flowDoms(nn, groundLevel, sps)%x(i, j, k, 2))
                            xc(3) = fourth * (flowDoms(nn, groundLevel, sps)%x(i, j - 1, k, 3) + &
                                              flowDoms(nn, groundLevel, sps)%x(i - 1, j, k, 3) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, j - 1, k, 3) + &
                                              flowDoms(nn, groundLevel, sps)%x(i, j, k, 3))

                            call cellFaceVelocities(xc, rotCenter, rotRate, &
                                                    velxGrid, velyGrid, velzGrid, derivRotationMatrix, sc)

                            ! Store the dot product of grid velocity sc and
                            ! the normal ss in sFace.

                            sFaceK(i, j, k) = sc(1) * sk(i, j, k, 1) + sc(2) * sk(i, j, k, 2) &
                                              + sc(3) * sk(i, j, k, 3)
                        end do
                    end do
                end do

            end if testUseOldCoor
        end if testMoving

    end subroutine gridVelocitiesFineLevel_block

    subroutine cellFaceVelocities(xc, rotCenter, rotRate, velxGrid, velyGrid, velzGrid, derivRotationMatrix, sc)
        !
        !  Returns the cell face velocities for a given face center
        !
        use constants

        implicit none
        !
        !      Subroutine arguments.
        !
        real(kind=realType), dimension(3), intent(in) :: xc, rotCenter, rotRate
        real(kind=realType), intent(in) :: velxGrid, velyGrid, velzGrid
        real(kind=realType), dimension(3, 3), intent(in) :: derivRotationMatrix
        real(kind=realType), dimension(3), intent(out) :: sc
        !
        !      Local variables.
        !
        real(kind=realType), dimension(3) :: rotationPoint, xxc

        ! Determine the coordinates relative to the
        ! center of rotation.

        xxc(1) = xc(1) - rotCenter(1)
        xxc(2) = xc(2) - rotCenter(2)
        xxc(3) = xc(3) - rotCenter(3)

        ! Determine the rotation speed of the face center,
        ! which is omega*r.

        sc(1) = rotRate(2) * xxc(3) - rotRate(3) * xxc(2)
        sc(2) = rotRate(3) * xxc(1) - rotRate(1) * xxc(3)
        sc(3) = rotRate(1) * xxc(2) - rotRate(2) * xxc(1)

        ! Determine the coordinates relative to the
        ! rigid body rotation point.

        xxc(1) = xc(1) - rotationPoint(1)
        xxc(2) = xc(2) - rotationPoint(2)
        xxc(3) = xc(3) - rotationPoint(3)

        ! Determine the total velocity of the cell face.
        ! This is a combination of rotation speed of this
        ! block and the entire rigid body rotation.

        sc(1) = sc(1) + velxGrid &
                + derivRotationMatrix(1, 1) * xxc(1) &
                + derivRotationMatrix(1, 2) * xxc(2) &
                + derivRotationMatrix(1, 3) * xxc(3)
        sc(2) = sc(2) + velyGrid &
                + derivRotationMatrix(2, 1) * xxc(1) &
                + derivRotationMatrix(2, 2) * xxc(2) &
                + derivRotationMatrix(2, 3) * xxc(3)
        sc(3) = sc(3) + velzGrid &
                + derivRotationMatrix(3, 1) * xxc(1) &
                + derivRotationMatrix(3, 2) * xxc(2) &
                + derivRotationMatrix(3, 3) * xxc(3)

    end subroutine cellFaceVelocities


    subroutine slipVelocitiesFineLevel_block(useOldCoor, t, sps, nn)
        !
        !       slipVelocitiesFineLevel computes the slip velocities for
        !       viscous subfaces on all viscous boundaries on groundLevel for
        !       the given spectral solution. If useOldCoor is .true. the
        !       velocities are determined using the unsteady time integrator;
        !       otherwise the analytic form is used.
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
        use flowUtils, only: derivativeRotMatrixRigid, getDirVector
        use utils, only: tsAlpha, tsBeta, tsMach, terminate, rotMatrixRigidBody, &
                         setCoefTimeIntegrator, getDirAngle
        implicit none
        !
        !      Subroutine arguments.
        !
        integer(kind=intType), intent(in) :: sps, nn
        logical, intent(in) :: useOldCoor

        real(kind=realType), dimension(*), intent(in) :: t
        !
        !      Local variables.
        !
        integer(kind=intType) :: mm, i, j, level, ii

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

        ! Determine the situation we are having here.

        testUseOldCoor: if (useOldCoor) then

            continue
        else

            ! The velocities must be determined analytically.

            ! Compute the mesh velocity from the given mesh Mach number.

            !  aInf = sqrt(gammaInf*pInf/rhoInf)
            !  velxGrid = aInf*MachGrid(1)
            !  velyGrid = aInf*MachGrid(2)
            !  velzGrid = aInf*MachGrid(3)

            aInf = sqrt(gammaInf * pInf / rhoInf)
            velxGrid0 = (aInf * machgrid) * (-velDirFreestream(1))
            velyGrid0 = (aInf * machgrid) * (-velDirFreestream(2))
            velzGrid0 = (aInf * machgrid) * (-velDirFreestream(3))

            ! Compute the derivative of the rotation matrix and the rotation
            ! point; needed for velocity due to the rigid body rotation of
            ! the entire grid. It is assumed that the rigid body motion of
            ! the grid is only specified if there is only 1 section present.

            call derivativeRotMatrixRigid(derivRotationMatrix, rotationPoint, &
                                          t(1))

            !compute the rotation matrix to update the velocities for the time
            !spectral stability derivative case...


            ! Loop over the number of viscous subfaces.

            bocoLoop2: do mm = 1, nViscBocos

                ! Store the rotation center and the rotation rate
                ! for this subface.

                ii = cgnsSubface(mm)

                rotCenter = cgnsDoms(nbkGlobal)%bocoInfo(ii)%rotCenter
                rotRate = timeRef * cgnsDoms(nbkGlobal)%bocoInfo(ii)%rotRate

                ! useWindAxis should go back here!
                velXgrid = velXGrid0
                velYgrid = velYGrid0
                velZgrid = velZGrid0

                ! Loop over the quadrilateral faces of the viscous
                ! subface.

                ! The new procedure is less elegant as the previous one.
                ! But the new stands up to Tapenade.
                if (BCFaceID(mm) == iMin) then

                    do j = BCData(mm)%jcBeg, BCData(mm)%jcEnd
                        do i = BCData(mm)%icBeg, BCData(mm)%icEnd

                            ! Compute the coordinates of the centroid of the face.
                            ! Normally this would be an average of i-1 and i, but
                            ! due to the usage of the pointer xFace and the fact
                            ! that x starts at index 0 this is shifted 1 index.

                            xc(1) = fourth * (flowDoms(nn, groundLevel, sps)%x(1, i, j, 1) &
                                              + flowDoms(nn, groundLevel, sps)%x(1, i, j - 1, 1) &
                                              + flowDoms(nn, groundLevel, sps)%x(1, i - 1, j, 1) &
                                              + flowDoms(nn, groundLevel, sps)%x(1, i - 1, j - 1, 1))
                            xc(2) = fourth * (flowDoms(nn, groundLevel, sps)%x(1, i, j, 2) &
                                              + flowDoms(nn, groundLevel, sps)%x(1, i, j - 1, 2) &
                                              + flowDoms(nn, groundLevel, sps)%x(1, i - 1, j, 2) &
                                              + flowDoms(nn, groundLevel, sps)%x(1, i - 1, j - 1, 2))
                            xc(3) = fourth * (flowDoms(nn, groundLevel, sps)%x(1, i, j, 3) &
                                              + flowDoms(nn, groundLevel, sps)%x(1, i, j - 1, 3) &
                                              + flowDoms(nn, groundLevel, sps)%x(1, i - 1, j, 3) &
                                              + flowDoms(nn, groundLevel, sps)%x(1, i - 1, j - 1, 3))

                            ! Determine the coordinates relative to the center
                            ! of rotation.

                            xxc(1) = xc(1) - rotCenter(1)
                            xxc(2) = xc(2) - rotCenter(2)
                            xxc(3) = xc(3) - rotCenter(3)

                            ! Compute the velocity, which is the cross product
                            ! of rotRate and xc.

                            BCData(mm)%uSlip(i, j, 1) = rotRate(2) * xxc(3) - rotRate(3) * xxc(2)
                            BCData(mm)%uSlip(i, j, 2) = rotRate(3) * xxc(1) - rotRate(1) * xxc(3)
                            BCData(mm)%uSlip(i, j, 3) = rotRate(1) * xxc(2) - rotRate(2) * xxc(1)

                            ! Determine the coordinates relative to the
                            ! rigid body rotation point.

                            xxc(1) = xc(1) - rotationPoint(1)
                            xxc(2) = xc(2) - rotationPoint(2)
                            xxc(3) = xc(3) - rotationPoint(3)

                            ! Determine the total velocity of the cell center.
                            ! This is a combination of rotation speed of this
                            ! block and the entire rigid body rotation.

                            BCData(mm)%uSlip(i, j, 1) = BCData(mm)%uSlip(i, j, 1) + velxGrid &
                                                        + derivRotationMatrix(1, 1) * xxc(1) &
                                                        + derivRotationMatrix(1, 2) * xxc(2) &
                                                        + derivRotationMatrix(1, 3) * xxc(3)
                            BCData(mm)%uSlip(i, j, 2) = BCData(mm)%uSlip(i, j, 2) + velyGrid &
                                                        + derivRotationMatrix(2, 1) * xxc(1) &
                                                        + derivRotationMatrix(2, 2) * xxc(2) &
                                                        + derivRotationMatrix(2, 3) * xxc(3)
                            BCData(mm)%uSlip(i, j, 3) = BCData(mm)%uSlip(i, j, 3) + velzGrid &
                                                        + derivRotationMatrix(3, 1) * xxc(1) &
                                                        + derivRotationMatrix(3, 2) * xxc(2) &
                                                        + derivRotationMatrix(3, 3) * xxc(3)
                        end do
                    end do

                else if (BCFaceID(mm) == iMax) then

                    do j = BCData(mm)%jcBeg, BCData(mm)%jcEnd
                        do i = BCData(mm)%icBeg, BCData(mm)%icEnd

                            ! Compute the coordinates of the centroid of the face.
                            ! Normally this would be an average of i-1 and i, but
                            ! due to the usage of the pointer xFace and the fact
                            ! that x starts at index 0 this is shifted 1 index.

                            xc(1) = fourth * (flowDoms(nn, groundLevel, sps)%x(il, i, j, 1) &
                                              + flowDoms(nn, groundLevel, sps)%x(il, i, j - 1, 1) &
                                              + flowDoms(nn, groundLevel, sps)%x(il, i - 1, j, 1) &
                                              + flowDoms(nn, groundLevel, sps)%x(il, i - 1, j - 1, 1))
                            xc(2) = fourth * (flowDoms(nn, groundLevel, sps)%x(il, i, j, 2) &
                                              + flowDoms(nn, groundLevel, sps)%x(il, i, j - 1, 2) &
                                              + flowDoms(nn, groundLevel, sps)%x(il, i - 1, j, 2) &
                                              + flowDoms(nn, groundLevel, sps)%x(il, i - 1, j - 1, 2))
                            xc(3) = fourth * (flowDoms(nn, groundLevel, sps)%x(il, i, j, 3) &
                                              + flowDoms(nn, groundLevel, sps)%x(il, i, j - 1, 3) &
                                              + flowDoms(nn, groundLevel, sps)%x(il, i - 1, j, 3) &
                                              + flowDoms(nn, groundLevel, sps)%x(il, i - 1, j - 1, 3))

                            ! Determine the coordinates relative to the center
                            ! of rotation.

                            xxc(1) = xc(1) - rotCenter(1)
                            xxc(2) = xc(2) - rotCenter(2)
                            xxc(3) = xc(3) - rotCenter(3)

                            ! Compute the velocity, which is the cross product
                            ! of rotRate and xc.

                            BCData(mm)%uSlip(i, j, 1) = rotRate(2) * xxc(3) - rotRate(3) * xxc(2)
                            BCData(mm)%uSlip(i, j, 2) = rotRate(3) * xxc(1) - rotRate(1) * xxc(3)
                            BCData(mm)%uSlip(i, j, 3) = rotRate(1) * xxc(2) - rotRate(2) * xxc(1)

                            ! Determine the coordinates relative to the
                            ! rigid body rotation point.

                            xxc(1) = xc(1) - rotationPoint(1)
                            xxc(2) = xc(2) - rotationPoint(2)
                            xxc(3) = xc(3) - rotationPoint(3)

                            ! Determine the total velocity of the cell center.
                            ! This is a combination of rotation speed of this
                            ! block and the entire rigid body rotation.

                            BCData(mm)%uSlip(i, j, 1) = BCData(mm)%uSlip(i, j, 1) + velxGrid &
                                                        + derivRotationMatrix(1, 1) * xxc(1) &
                                                        + derivRotationMatrix(1, 2) * xxc(2) &
                                                        + derivRotationMatrix(1, 3) * xxc(3)
                            BCData(mm)%uSlip(i, j, 2) = BCData(mm)%uSlip(i, j, 2) + velyGrid &
                                                        + derivRotationMatrix(2, 1) * xxc(1) &
                                                        + derivRotationMatrix(2, 2) * xxc(2) &
                                                        + derivRotationMatrix(2, 3) * xxc(3)
                            BCData(mm)%uSlip(i, j, 3) = BCData(mm)%uSlip(i, j, 3) + velzGrid &
                                                        + derivRotationMatrix(3, 1) * xxc(1) &
                                                        + derivRotationMatrix(3, 2) * xxc(2) &
                                                        + derivRotationMatrix(3, 3) * xxc(3)
                        end do
                    end do

                else if (BCFaceID(mm) == jMin) then

                    do j = BCData(mm)%jcBeg, BCData(mm)%jcEnd
                        do i = BCData(mm)%icBeg, BCData(mm)%icEnd

                            ! Compute the coordinates of the centroid of the face.
                            ! Normally this would be an average of i-1 and i, but
                            ! due to the usage of the pointer xFace and the fact
                            ! that x starts at index 0 this is shifted 1 index.

                            xc(1) = fourth * (flowDoms(nn, groundLevel, sps)%x(i, 1, j, 1) &
                                              + flowDoms(nn, groundLevel, sps)%x(i, 1, j - 1, 1) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, 1, j, 1) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, 1, j - 1, 1))
                            xc(2) = fourth * (flowDoms(nn, groundLevel, sps)%x(i, 1, j, 2) &
                                              + flowDoms(nn, groundLevel, sps)%x(i, 1, j - 1, 2) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, 1, j, 2) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, 1, j - 1, 2))
                            xc(3) = fourth * (flowDoms(nn, groundLevel, sps)%x(i, 1, j, 3) &
                                              + flowDoms(nn, groundLevel, sps)%x(i, 1, j - 1, 3) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, 1, j, 3) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, 1, j - 1, 3))

                            ! Determine the coordinates relative to the center
                            ! of rotation.

                            xxc(1) = xc(1) - rotCenter(1)
                            xxc(2) = xc(2) - rotCenter(2)
                            xxc(3) = xc(3) - rotCenter(3)

                            ! Compute the velocity, which is the cross product
                            ! of rotRate and xc.

                            BCData(mm)%uSlip(i, j, 1) = rotRate(2) * xxc(3) - rotRate(3) * xxc(2)
                            BCData(mm)%uSlip(i, j, 2) = rotRate(3) * xxc(1) - rotRate(1) * xxc(3)
                            BCData(mm)%uSlip(i, j, 3) = rotRate(1) * xxc(2) - rotRate(2) * xxc(1)

                            ! Determine the coordinates relative to the
                            ! rigid body rotation point.

                            xxc(1) = xc(1) - rotationPoint(1)
                            xxc(2) = xc(2) - rotationPoint(2)
                            xxc(3) = xc(3) - rotationPoint(3)

                            ! Determine the total velocity of the cell center.
                            ! This is a combination of rotation speed of this
                            ! block and the entire rigid body rotation.

                            BCData(mm)%uSlip(i, j, 1) = BCData(mm)%uSlip(i, j, 1) + velxGrid &
                                                        + derivRotationMatrix(1, 1) * xxc(1) &
                                                        + derivRotationMatrix(1, 2) * xxc(2) &
                                                        + derivRotationMatrix(1, 3) * xxc(3)
                            BCData(mm)%uSlip(i, j, 2) = BCData(mm)%uSlip(i, j, 2) + velyGrid &
                                                        + derivRotationMatrix(2, 1) * xxc(1) &
                                                        + derivRotationMatrix(2, 2) * xxc(2) &
                                                        + derivRotationMatrix(2, 3) * xxc(3)
                            BCData(mm)%uSlip(i, j, 3) = BCData(mm)%uSlip(i, j, 3) + velzGrid &
                                                        + derivRotationMatrix(3, 1) * xxc(1) &
                                                        + derivRotationMatrix(3, 2) * xxc(2) &
                                                        + derivRotationMatrix(3, 3) * xxc(3)
                        end do
                    end do

                else if (BCFaceID(mm) == jMax) then

                    do j = BCData(mm)%jcBeg, BCData(mm)%jcEnd
                        do i = BCData(mm)%icBeg, BCData(mm)%icEnd

                            ! Compute the coordinates of the centroid of the face.
                            ! Normally this would be an average of i-1 and i, but
                            ! due to the usage of the pointer xFace and the fact
                            ! that x starts at index 0 this is shifted 1 index.

                            xc(1) = fourth * (flowDoms(nn, groundLevel, sps)%x(i, jl, j, 1) &
                                              + flowDoms(nn, groundLevel, sps)%x(i, jl, j - 1, 1) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, jl, j, 1) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, jl, j - 1, 1))
                            xc(2) = fourth * (flowDoms(nn, groundLevel, sps)%x(i, jl, j, 2) &
                                              + flowDoms(nn, groundLevel, sps)%x(i, jl, j - 1, 2) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, jl, j, 2) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, jl, j - 1, 2))
                            xc(3) = fourth * (flowDoms(nn, groundLevel, sps)%x(i, jl, j, 3) &
                                              + flowDoms(nn, groundLevel, sps)%x(i, jl, j - 1, 3) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, jl, j, 3) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, jl, j - 1, 3))

                            ! Determine the coordinates relative to the center
                            ! of rotation.

                            xxc(1) = xc(1) - rotCenter(1)
                            xxc(2) = xc(2) - rotCenter(2)
                            xxc(3) = xc(3) - rotCenter(3)

                            ! Compute the velocity, which is the cross product
                            ! of rotRate and xc.

                            BCData(mm)%uSlip(i, j, 1) = rotRate(2) * xxc(3) - rotRate(3) * xxc(2)
                            BCData(mm)%uSlip(i, j, 2) = rotRate(3) * xxc(1) - rotRate(1) * xxc(3)
                            BCData(mm)%uSlip(i, j, 3) = rotRate(1) * xxc(2) - rotRate(2) * xxc(1)

                            ! Determine the coordinates relative to the
                            ! rigid body rotation point.

                            xxc(1) = xc(1) - rotationPoint(1)
                            xxc(2) = xc(2) - rotationPoint(2)
                            xxc(3) = xc(3) - rotationPoint(3)

                            ! Determine the total velocity of the cell center.
                            ! This is a combination of rotation speed of this
                            ! block and the entire rigid body rotation.

                            BCData(mm)%uSlip(i, j, 1) = BCData(mm)%uSlip(i, j, 1) + velxGrid &
                                                        + derivRotationMatrix(1, 1) * xxc(1) &
                                                        + derivRotationMatrix(1, 2) * xxc(2) &
                                                        + derivRotationMatrix(1, 3) * xxc(3)
                            BCData(mm)%uSlip(i, j, 2) = BCData(mm)%uSlip(i, j, 2) + velyGrid &
                                                        + derivRotationMatrix(2, 1) * xxc(1) &
                                                        + derivRotationMatrix(2, 2) * xxc(2) &
                                                        + derivRotationMatrix(2, 3) * xxc(3)
                            BCData(mm)%uSlip(i, j, 3) = BCData(mm)%uSlip(i, j, 3) + velzGrid &
                                                        + derivRotationMatrix(3, 1) * xxc(1) &
                                                        + derivRotationMatrix(3, 2) * xxc(2) &
                                                        + derivRotationMatrix(3, 3) * xxc(3)
                        end do
                    end do

                else if (BCFaceID(mm) == kMin) then

                    do j = BCData(mm)%jcBeg, BCData(mm)%jcEnd
                        do i = BCData(mm)%icBeg, BCData(mm)%icEnd

                            ! Compute the coordinates of the centroid of the face.
                            ! Normally this would be an average of i-1 and i, but
                            ! due to the usage of the pointer xFace and the fact
                            ! that x starts at index 0 this is shifted 1 index.

                            xc(1) = fourth * (flowDoms(nn, groundLevel, sps)%x(i, j, 1, 1) &
                                              + flowDoms(nn, groundLevel, sps)%x(i, j - 1, 1, 1) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, j, 1, 1) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, j - 1, 1, 1))
                            xc(2) = fourth * (flowDoms(nn, groundLevel, sps)%x(i, j, 1, 2) &
                                              + flowDoms(nn, groundLevel, sps)%x(i, j - 1, 1, 2) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, j, 1, 2) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, j - 1, 1, 2))
                            xc(3) = fourth * (flowDoms(nn, groundLevel, sps)%x(i, j, 1, 3) &
                                              + flowDoms(nn, groundLevel, sps)%x(i, j - 1, 1, 3) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, j, 1, 3) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, j - 1, 1, 3))

                            ! Determine the coordinates relative to the center
                            ! of rotation.

                            xxc(1) = xc(1) - rotCenter(1)
                            xxc(2) = xc(2) - rotCenter(2)
                            xxc(3) = xc(3) - rotCenter(3)

                            ! Compute the velocity, which is the cross product
                            ! of rotRate and xc.

                            BCData(mm)%uSlip(i, j, 1) = rotRate(2) * xxc(3) - rotRate(3) * xxc(2)
                            BCData(mm)%uSlip(i, j, 2) = rotRate(3) * xxc(1) - rotRate(1) * xxc(3)
                            BCData(mm)%uSlip(i, j, 3) = rotRate(1) * xxc(2) - rotRate(2) * xxc(1)

                            ! Determine the coordinates relative to the
                            ! rigid body rotation point.

                            xxc(1) = xc(1) - rotationPoint(1)
                            xxc(2) = xc(2) - rotationPoint(2)
                            xxc(3) = xc(3) - rotationPoint(3)

                            ! Determine the total velocity of the cell center.
                            ! This is a combination of rotation speed of this
                            ! block and the entire rigid body rotation.

                            BCData(mm)%uSlip(i, j, 1) = BCData(mm)%uSlip(i, j, 1) + velxGrid &
                                                        + derivRotationMatrix(1, 1) * xxc(1) &
                                                        + derivRotationMatrix(1, 2) * xxc(2) &
                                                        + derivRotationMatrix(1, 3) * xxc(3)
                            BCData(mm)%uSlip(i, j, 2) = BCData(mm)%uSlip(i, j, 2) + velyGrid &
                                                        + derivRotationMatrix(2, 1) * xxc(1) &
                                                        + derivRotationMatrix(2, 2) * xxc(2) &
                                                        + derivRotationMatrix(2, 3) * xxc(3)
                            BCData(mm)%uSlip(i, j, 3) = BCData(mm)%uSlip(i, j, 3) + velzGrid &
                                                        + derivRotationMatrix(3, 1) * xxc(1) &
                                                        + derivRotationMatrix(3, 2) * xxc(2) &
                                                        + derivRotationMatrix(3, 3) * xxc(3)
                        end do
                    end do

                else if (BCFaceID(mm) == kMax) then

                    do j = BCData(mm)%jcBeg, BCData(mm)%jcEnd
                        do i = BCData(mm)%icBeg, BCData(mm)%icEnd

                            ! Compute the coordinates of the centroid of the face.
                            ! Normally this would be an average of i-1 and i, but
                            ! due to the usage of the pointer xFace and the fact
                            ! that x starts at index 0 this is shifted 1 index.

                            xc(1) = fourth * (flowDoms(nn, groundLevel, sps)%x(i, j, kl, 1) &
                                              + flowDoms(nn, groundLevel, sps)%x(i, j - 1, kl, 1) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, j, kl, 1) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, j - 1, kl, 1))
                            xc(2) = fourth * (flowDoms(nn, groundLevel, sps)%x(i, j, kl, 2) &
                                              + flowDoms(nn, groundLevel, sps)%x(i, j - 1, kl, 2) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, j, kl, 2) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, j - 1, kl, 2))
                            xc(3) = fourth * (flowDoms(nn, groundLevel, sps)%x(i, j, kl, 3) &
                                              + flowDoms(nn, groundLevel, sps)%x(i, j - 1, kl, 3) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, j, kl, 3) &
                                              + flowDoms(nn, groundLevel, sps)%x(i - 1, j - 1, kl, 3))

                            ! Determine the coordinates relative to the center
                            ! of rotation.

                            xxc(1) = xc(1) - rotCenter(1)
                            xxc(2) = xc(2) - rotCenter(2)
                            xxc(3) = xc(3) - rotCenter(3)

                            ! Compute the velocity, which is the cross product
                            ! of rotRate and xc.

                            BCData(mm)%uSlip(i, j, 1) = rotRate(2) * xxc(3) - rotRate(3) * xxc(2)
                            BCData(mm)%uSlip(i, j, 2) = rotRate(3) * xxc(1) - rotRate(1) * xxc(3)
                            BCData(mm)%uSlip(i, j, 3) = rotRate(1) * xxc(2) - rotRate(2) * xxc(1)

                            ! Determine the coordinates relative to the
                            ! rigid body rotation point.

                            xxc(1) = xc(1) - rotationPoint(1)
                            xxc(2) = xc(2) - rotationPoint(2)
                            xxc(3) = xc(3) - rotationPoint(3)

                            ! Determine the total velocity of the cell center.
                            ! This is a combination of rotation speed of this
                            ! block and the entire rigid body rotation.

                            BCData(mm)%uSlip(i, j, 1) = BCData(mm)%uSlip(i, j, 1) + velxGrid &
                                                        + derivRotationMatrix(1, 1) * xxc(1) &
                                                        + derivRotationMatrix(1, 2) * xxc(2) &
                                                        + derivRotationMatrix(1, 3) * xxc(3)
                            BCData(mm)%uSlip(i, j, 2) = BCData(mm)%uSlip(i, j, 2) + velyGrid &
                                                        + derivRotationMatrix(2, 1) * xxc(1) &
                                                        + derivRotationMatrix(2, 2) * xxc(2) &
                                                        + derivRotationMatrix(2, 3) * xxc(3)
                            BCData(mm)%uSlip(i, j, 3) = BCData(mm)%uSlip(i, j, 3) + velzGrid &
                                                        + derivRotationMatrix(3, 1) * xxc(1) &
                                                        + derivRotationMatrix(3, 2) * xxc(2) &
                                                        + derivRotationMatrix(3, 3) * xxc(3)
                        end do
                    end do

                end if

            end do bocoLoop2

        end if testUseOldCoor

    end subroutine slipVelocitiesFineLevel_block


    subroutine normalVelocities_block(sps)
        !
        !       normalVelocitiesAllLevels computes the normal grid
        !       velocities of some boundary faces of the moving blocks for
        !       spectral mode sps. All grid levels from ground level to the
        !       coarsest level are considered.
        !
        use constants
        use blockPointers, only: il, jl, kl, addGridVelocities, nBocos, BCData, &
                                 sfaceI, sfaceJ, sfaceK, bcFaceID, si, sj, sk
        !use iteration
        implicit none
        !
        !      Subroutine arguments.
        !
        integer(kind=intType), intent(in) :: sps
        !
        !      Local variables.
        !
        integer(kind=intType) :: mm
        integer(kind=intType) :: i, j

        real(kind=realType) :: weight, mult

        real(kind=realType), dimension(:, :), pointer :: sFace
        real(kind=realType), dimension(:, :, :), pointer :: ss

        ! Check for a moving block. As it is possible that in a
        ! multidisicplinary environment additional grid velocities
        ! are set, the test should be done on addGridVelocities
        ! and not on blockIsMoving.

        testMoving: if (addGridVelocities) then
            !
            !             Determine the normal grid velocities of the boundaries.
            !             As these values are based on the unit normal. A division
            !             by the length of the normal is needed.
            !             Furthermore the boundary unit normals are per definition
            !             outward pointing, while on the iMin, jMin and kMin
            !             boundaries the face normals are inward pointing. This
            !             is taken into account by the factor mult.
            !
            ! Loop over the boundary subfaces.

            bocoLoop: do mm = 1, nBocos

                ! Check whether rFace is allocated.

                testAssoc: if (associated(BCData(mm)%rFace)) then

                    ! Determine the block face on which the subface is
                    ! located and set some variables accordingly.

                    ! The new procedure is less elegant as the previous one.
                    ! But the new stands up to Tapenade.
                    if (BCFaceID(mm) == iMin) then

                        mult = -one

                        do j = BCData(mm)%jcBeg, BCData(mm)%jcEnd
                            do i = BCData(mm)%icBeg, BCData(mm)%icEnd

                                ! Compute the inverse of the length of the normal
                                ! vector and possibly correct for inward pointing.

                                weight = sqrt(si(1, i, j, 1)**2 + si(1, i, j, 2)**2 &
                                              + si(1, i, j, 3)**2)
                                if (weight > zero) weight = mult / weight

                                ! Compute the normal velocity based on the outward
                                ! pointing unit normal.

                                BCData(mm)%rFace(i, j) = weight * sFaceI(1, i, j)
                            end do
                        end do

                    else if (BCFaceID(mm) == iMax) then

                        mult = one

                        do j = BCData(mm)%jcBeg, BCData(mm)%jcEnd
                            do i = BCData(mm)%icBeg, BCData(mm)%icEnd

                                ! Compute the inverse of the length of the normal
                                ! vector and possibly correct for inward pointing.

                                weight = sqrt(si(il, i, j, 1)**2 + si(il, i, j, 2)**2 &
                                              + si(il, i, j, 3)**2)
                                if (weight > zero) weight = mult / weight

                                ! Compute the normal velocity based on the outward
                                ! pointing unit normal.

                                BCData(mm)%rFace(i, j) = weight * sFaceI(il, i, j)
                            end do
                        end do

                    else if (BCFaceID(mm) == jMin) then

                        mult = -one

                        do j = BCData(mm)%jcBeg, BCData(mm)%jcEnd
                            do i = BCData(mm)%icBeg, BCData(mm)%icEnd

                                ! Compute the inverse of the length of the normal
                                ! vector and possibly correct for inward pointing.

                                weight = sqrt(sj(i, 1, j, 1)**2 + sj(i, 1, j, 2)**2 &
                                              + sj(i, 1, j, 3)**2)
                                if (weight > zero) weight = mult / weight

                                ! Compute the normal velocity based on the outward
                                ! pointing unit normal.

                                BCData(mm)%rFace(i, j) = weight * sFaceJ(i, 1, j)
                            end do
                        end do

                    else if (BCFaceID(mm) == jMax) then

                        mult = one

                        do j = BCData(mm)%jcBeg, BCData(mm)%jcEnd
                            do i = BCData(mm)%icBeg, BCData(mm)%icEnd

                                ! Compute the inverse of the length of the normal
                                ! vector and possibly correct for inward pointing.

                                weight = sqrt(sj(i, jl, j, 1)**2 + sj(i, jl, j, 2)**2 &
                                              + sj(i, jl, j, 3)**2)
                                if (weight > zero) weight = mult / weight

                                ! Compute the normal velocity based on the outward
                                ! pointing unit normal.

                                BCData(mm)%rFace(i, j) = weight * sFaceJ(i, jl, j)
                            end do
                        end do

                    else if (BCFaceID(mm) == kMin) then

                        mult = -one

                        do j = BCData(mm)%jcBeg, BCData(mm)%jcEnd
                            do i = BCData(mm)%icBeg, BCData(mm)%icEnd

                                ! Compute the inverse of the length of the normal
                                ! vector and possibly correct for inward pointing.

                                weight = sqrt(sk(i, j, 1, 1)**2 + sk(i, j, 1, 2)**2 &
                                              + sk(i, j, 1, 3)**2)
                                if (weight > zero) weight = mult / weight

                                ! Compute the normal velocity based on the outward
                                ! pointing unit normal.

                                BCData(mm)%rFace(i, j) = weight * sFaceK(i, j, 1)
                            end do
                        end do

                    else if (BCFaceID(mm) == kMax) then

                        mult = one

                        do j = BCData(mm)%jcBeg, BCData(mm)%jcEnd
                            do i = BCData(mm)%icBeg, BCData(mm)%icEnd

                                ! Compute the inverse of the length of the normal
                                ! vector and possibly correct for inward pointing.

                                weight = sqrt(sk(i, j, kl, 1)**2 + sk(i, j, kl, 2)**2 &
                                              + sk(i, j, kl, 3)**2)
                                if (weight > zero) weight = mult / weight

                                ! Compute the normal velocity based on the outward
                                ! pointing unit normal.

                                BCData(mm)%rFace(i, j) = weight * sFaceK(i, j, kl)
                            end do
                        end do

                    end if

                end if testAssoc
            end do bocoLoop

            else testMoving

            ! Block is not moving. Loop over the boundary faces and set
            ! the normal grid velocity to zero if allocated.

            do mm = 1, nBocos
                if (associated(BCData(mm)%rFace)) &
                    BCData(mm)%rFace = zero
            end do

        end if testMoving

    end subroutine normalVelocities_block

    ! ----------------------------------------------------------------------
    !                                                                      |
    !                    No Tapenade Routine below this line               |
    !                                                                      |
    ! ----------------------------------------------------------------------

end module solverUtils
