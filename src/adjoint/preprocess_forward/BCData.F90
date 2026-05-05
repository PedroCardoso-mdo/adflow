










module BCData
    use constants
    use BCDataMod

contains
    ! ---------------------------------------------------------------
    ! Routines that set the appropriate variable names for BCs with
    ! BCdata.

    subroutine setBCVarNamesIsothermalWall
        use cgnsNames
        use constants
        implicit none
        nbcVar = nbcVarIsothermalWall
        bcVarNames(1) = cgnsTemp

    end subroutine setBCVarNamesIsothermalWall

    subroutine setBCVarNamesSubsonicInflow
        use constants
        use cgnsNames
        use inputPhysics, only: equations
        use flowVarRefState, only: nwt
        implicit none
        !
        !      Local variables.
        !
        logical :: varAllowed

        nbcVar = nbcVarSubsonicInflow
        if (equations == RANSEquations) then
            nbcVar = nbcVar + nwt
        end if

        bcVarNames(1) = cgnsPtot
        bcVarNames(2) = cgnsTtot
        bcVarNames(3) = cgnsRhotot
        bcVarNames(4) = cgnsVelAnglex
        bcVarNames(5) = cgnsVelAngley
        bcVarNames(6) = cgnsVelAnglez
        bcVarNames(7) = cgnsVelVecx
        bcVarNames(8) = cgnsVelVecy
        bcVarNames(9) = cgnsVelVecz
        bcVarNames(10) = cgnsVelVecr
        bcVarNames(11) = cgnsVelVectheta
        bcVarNames(12) = cgnsDensity
        bcVarNames(13) = cgnsVelx
        bcVarNames(14) = cgnsVely
        bcVarNames(15) = cgnsVelz
        bcVarNames(16) = cgnsVelr
        bcVarNames(17) = cgnsVeltheta

        call setBcVarNamesTurb(17_intType)

    end subroutine setBCVarNamesSubsonicInflow

    subroutine setBCVarNamesSubsonicOutflow
        use cgnsNames
        use constants
        use flowVarRefState, only: nwt

        nbcVar = nbcVarSubsonicOutflow

        bcVarNames(1) = cgnsPressure

    end subroutine setBCVarNamesSubsonicOutflow

    subroutine setBCVarNamesSupersonicInflow
        use constants
        use cgnsNames
        use inputPhysics, only: equations
        use flowVarRefState, only: nwt

        nbcVar = nbcVarSupersonicInflow
        if (equations == RANSEquations) then
            nbcVar = nbcVar + nwt
        end if

        bcVarNames(1) = cgnsDensity
        bcVarNames(2) = cgnsPressure
        bcVarNames(3) = cgnsVelx
        bcVarNames(4) = cgnsVely
        bcVarNames(5) = cgnsVelz
        bcVarNames(6) = cgnsVelr
        bcVarNames(7) = cgnsVeltheta

        call setBCVarNamesTurb(7_intType)

    end subroutine setBCVarNamesSupersonicInflow

    subroutine setBCVarNamesTurb(offset)
        !
        !       setBCVarNamesTurb sets the names for the turbulence
        !       variables to be determined. This depends on the turbulence
        !       model. If not the RANS equations are solved an immediate
        !       return is made.
        !
        use constants
        use cgnsNames
        use inputPhysics, only: equations, turbModel
        implicit none
        !
        !      Subroutine arguments.
        !
        integer(kind=intType), intent(in) :: offset

        ! Return immediately if not the RANS equations are solved.

        if (equations /= RANSEquations) return

        ! Determine the turbulence model and set the names accordingly.

        select case (turbModel)
        case (spalartAllmaras, spalartAllmarasEdwards)
            bcVarNames(offset + 1) = cgnsTurbSaNu

        case (spalartallmarasnoft2gammaretheta)
            bcVarNames(offset + 1) = cgnsTurbSaNu
            bcVarNames(offset + 2) = cgnsTurbGamma
            bcVarNames(offset + 3) = cgnsTurbRetheta

        case (komegaWilcox, komegaModified, menterSST)
            bcVarNames(offset + 1) = cgnsTurbK
            bcVarNames(offset + 2) = cgnsTurbOmega

        case (ktau)
            bcVarNames(offset + 1) = cgnsTurbK
            bcVarNames(offset + 2) = cgnsTurbTau

        case (v2f)
            bcVarNames(offset + 1) = cgnsTurbK
            bcVarNames(offset + 2) = cgnsTurbEpsilon
            bcVarNames(offset + 3) = cgnsTurbV2
            bcVarNames(offset + 4) = cgnsTurbF

        end select

    end subroutine setBCVarNamesTurb
    ! ---------------------------------------------------------------
    ! --------------------------------------
    !                Utilities
    ! --------------------------------------

    subroutine computeHtot(tt, ht)
        !
        !       computeHtot computes the total enthalpy from the given total
        !       temperature. The total enthalpy is the integral of cp, which
        !       is a very simple expression for constant cp. For a variable cp
        !       it is a bit more work.
        !
        use constants
        use cpCurveFits
        use communication, only: myid
        use inputPhysics, only: cpModel, gammaConstant, rGasDim
        use flowVarRefState, only: PinfDim
        implicit none
        !
        !      Subroutine arguments.
        !
        real(kind=realType), intent(in) :: tt
        real(kind=realType), intent(out) :: ht
        !
        !      Local variables.
        !
        integer(kind=intType) :: ii, nn, mm, start

        real(kind=realType) :: t2

        ! Determine the cp model used in the computation.

        select case (cpModel)

        case (cpConstant)

            ! Constant cp. The total enthalpy is simply cp*tt.

            ht = gammaConstant * RGasDim * tt / (gammaConstant - one)

            !        ================================================================
        end select

    end subroutine computeHtot

    subroutine unitVectorsCylSystem(boco)
        !
        !       unitVectorsCylSystem determines the unit vectors of the
        !       local coordinate systen of the boundary face defined by the
        !       data in BCDataMod. In that local system the axial direction
        !       is rotation axis.
        !
        use constants
        use blockPointers, only: BCFaceID, BCData, x, si, sj, sk, il, jl, kl, &
                                 sectionID
        use section, only: sections
        implicit none
        !
        !      Subroutine arguments.
        !
        integer(kind=intType), intent(in) :: boco
        !
        !      Local variables.
        !
        integer(kind=intType) :: i, j
        real(kind=realType) :: factInlet, var

        real(kind=realType), dimension(3) :: dir

        real(kind=realType), dimension(:, :, :), pointer :: ss

        ! Set the pointers for coordinates and normals of the block
        ! face on which this subface is located. Set factInlet
        ! such that factInlet*normals points into the domain.

        select case (BCFaceID(boco))
        case (iMin)
            xf => x(1, :, :, :); ss => si(1, :, :, :); factInlet = one
        case (iMax)
            xf => x(il, :, :, :); ss => si(il, :, :, :); factInlet = -one
        case (jMin)
            xf => x(:, 1, :, :); ss => sj(:, 1, :, :); factInlet = one
        case (jMax)
            xf => x(:, jl, :, :); ss => sj(:, jl, :, :); factInlet = -one
        case (kMin)
            xf => x(:, :, 1, :); ss => sk(:, :, 1, :); factInlet = one
        case (kMax)
            xf => x(:, :, kl, :); ss => sk(:, :, kl, :); factInlet = -one
        end select

        ! Loop over the physical range of the subface to store the sum of
        ! the normals. Note that jBeg, jEnd, iBeg, iEnd cannot be used
        ! here, because they may include the halo faces. Instead the
        ! nodal range is used, which defines the original subface. The
        ! offset of +1 in the start index is there because you need
        ! the face id's.

        dir(1) = zero; dir(2) = zero; dir(3) = zero

        do j = (BCData(boco)%jnBeg + 1), BCData(boco)%jnEnd
            do i = (BCData(boco)%inBeg + 1), BCData(boco)%inEnd
                dir(1) = dir(1) + ss(i, j, 1)
                dir(2) = dir(2) + ss(i, j, 2)
                dir(3) = dir(3) + ss(i, j, 3)
            end do
        end do

        ! Multiply by factInlet to make sure that the normal
        ! is inward pointing.

        dir(1) = dir(1) * factInlet
        dir(2) = dir(2) * factInlet
        dir(3) = dir(3) * factInlet

        ! Determine three unit vectors, which define the local cartesian
        ! coordinate system of the rotation axis. First the axial
        ! direction. If the axis cannot be determined from rotation info,
        ! it is assumed to be the x-axis.

        axis = sections(sectionId)%rotAxis
        var = axis(1)**2 + axis(2)**2 + axis(3)**2
        if (var < half) then

            ! No rotation axis specified. Assume the x-axis
            ! and set the logical axAssumed to .True.

            axis(1) = one; axis(2) = zero; axis(3) = zero
            axAssumed = .true.
        end if

        ! The axial axis must be such that it points into the
        ! computational domain. If the dot product with dir is
        ! negative the direction of axis should be reversed.

        var = axis(1) * dir(1) + axis(2) * dir(2) + axis(3) * dir(3)
        if (var < zero) then
            axis(1) = -axis(1); axis(2) = -axis(2); axis(3) = -axis(3)
        end if

        ! Two unit vectors define the radial plane. These vectors are
        ! defined up to a constants. Just pick a direction for the second
        ! and create a unit vector normal to axis.

        if (abs(axis(2)) < 0.707107_realType) then
            radVec1(1) = zero; radVec1(2) = one; radVec1(3) = zero
        else
            radVec1(1) = zero; radVec1(2) = zero; radVec1(3) = one
        end if

        var = radVec1(1) * axis(1) + radVec1(2) * axis(2) &
              + radVec1(3) * axis(3)
        radVec1(1) = radVec1(1) - var * axis(1)
        radVec1(2) = radVec1(2) - var * axis(2)
        radVec1(3) = radVec1(3) - var * axis(3)

        var = one / sqrt(radVec1(1)**2 + radVec1(2)**2 &
                         + radVec1(3)**2)
        radVec1(1) = radVec1(1) * var
        radVec1(2) = radVec1(2) * var
        radVec1(3) = radVec1(3) * var

        ! The second vector of the radial plane is obtained
        ! by taking the cross product of axis and radVec1.

        radVec2(1) = axis(2) * radVec1(3) - axis(3) * radVec1(2)
        radVec2(2) = axis(3) * radVec1(1) - axis(1) * radVec1(3)
        radVec2(3) = axis(1) * radVec1(2) - axis(2) * radVec1(1)

    end subroutine unitVectorsCylSystem

    ! ---------------------------------------------------------------
    ! Routines that set the actual BCdata values from the CGNS data set
    ! information.
    ! ---------------------------------------------------------------

    subroutine BCDataIsothermalWall(boco, bcVarArray, iBeg, iEnd, jBeg, jEnd)
        !
        !       BCDataIsothermalWall tries to extract the wall temperature
        !       for the currently active boundary face, which is an isothermal
        !       viscous wall.
        !
        use constants
        use cgnsNames
        use blockPointers, only: BCFaceID, BCData, nBKGlobal
        use utils, only: terminate, siTemperature
        use flowVarRefState, only: Tref
        implicit none
        !
        !      Subroutine arguments.
        !
        integer(kind=intType) :: boco
        integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
        real(kind=realType), dimension(iBeg:iEnd, jBeg:jEnd, nbcVarMax) :: bcVarArray
        !
        !      Local variables.
        !
        integer :: ierr

        integer(kind=intType) :: i, j

        real(kind=realType) :: mult, trans

        character(len=maxStringLen) :: errorMessage

        ! Write an error message and terminate if it was not
        ! possible to determine the temperature.


        ! Convert to si-units and store the temperature in TNS_Wall.

        call siTemperature(temp(1), mult, trans)

        do j = jBeg, jEnd
            do i = iBeg, iEnd
                BCData(boco)%TNS_Wall(i, j) = (mult * bcVarArray(i, j, 1) + trans) / Tref
            end do
        end do

    end subroutine BCDataIsothermalWall

    subroutine BCDataSubsonicInflow(boco, bcVarArray, iBeg, iEnd, jBeg, jEnd, allTurbPresent)
        !
        !       BCDataSubsonicInflow tries to extract the prescribed data
        !       for the currently active boundary face, which is a subsonic
        !       inflow. Either total conditions and velocity direction or the
        !       velocity and density can be prescribed. In the latter case the
        !       mass flow is prescribed, which is okay as long as the flow is
        !       not choked.
        !
        use constants
        use cgnsNames
        use blockPointers, only: nbkGlobal, sectionID, BCFaceID, BCData
        use flowVarRefState, only: Tref, Pref, Href, rhoRef, muRef, nwt, wInf
        use inputPhysics, only: equations
        use utils, only: siDensity, siVelocity, siPressure, siAngle, &
                         siTemperature, terminate
        implicit none
        !
        !      Subroutine arguments.
        !
        integer(kind=intType), intent(in) :: boco
        integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
        real(kind=realType), dimension(iBeg:iEnd, jBeg:jEnd, nbcVarMax) :: bcVarArray
        logical, intent(inout) :: allTurbPresent
        !
        !      Local variables.
        !
        integer :: ierr, nn

        logical :: ptPresent, ttPresent, rhotPresent
        logical :: axPresent, ayPresent, azPresent
        logical :: xdirPresent, ydirPresent, zdirPresent
        logical :: rdirPresent, tdirPresent
        logical :: velxPresent, velyPresent, velzPresent
        logical :: rhoPresent, velrPresent, veltPresent
        logical :: totPresent, velPresent, dirPresent

        character(len=maxStringLen) :: errorMessage

        ! Store the logicals, which indicate succes or failure
        ! a bit more readable.

        ptPresent = bcVarPresent(1)
        ttPresent = bcVarPresent(2)
        rhotPresent = bcVarPresent(3)
        axPresent = bcVarPresent(4)
        ayPresent = bcVarPresent(5)
        azPresent = bcVarPresent(6)
        xdirPresent = bcVarPresent(7)
        ydirPresent = bcVarPresent(8)
        zdirPresent = bcVarPresent(9)
        rdirPresent = bcVarPresent(10)
        tdirPresent = bcVarPresent(11)
        rhoPresent = bcVarPresent(12)
        velxPresent = bcVarPresent(13)
        velyPresent = bcVarPresent(14)
        velzPresent = bcVarPresent(15)
        velrPresent = bcVarPresent(16)
        veltPresent = bcVarPresent(17)

        ! Check if the total conditions are present.

        nn = 0
        if (ptPresent) nn = nn + 1
        if (ttPresent) nn = nn + 1
        if (rhotPresent) nn = nn + 1

        totPresent = .false.
        if (nn >= 2) totPresent = .true.

        ! Check if a velocity direction is present.

        dirPresent = .false.
        if (xdirPresent .and. rdirPresent) dirPresent = .true.
        if ((axPresent .or. xdirPresent) .and. &
            (ayPresent .or. ydirPresent) .and. &
            (azPresent .or. zdirPresent)) dirPresent = .true.

        ! Check if a velocity vector is present.

        velPresent = .false.
        if (velxPresent .and. velrPresent) velPresent = .true.
        if (velxPresent .and. velyPresent .and. velzPresent) &
            velPresent = .true.

        ! Determine the situation we have here.

        if (totPresent .and. dirPresent) then

            ! Total conditions and velocity direction are prescribed.
            ! Determine the values for the faces of the subface.

            call totalSubsonicInlet

        else

            ! Not enough data is prescribed. Print an error message
            ! and exit.


        end if

        ! Set the turbulence variables and check if all of them are
        ! prescribed. If not set allTurbPresent to .false.

        allTurbPresent = setBcVarTurb(17_intType, boco, bcVarArray, &
                                      iBeg, iEnd, jBeg, jEnd, BCData(boco)%turbInlet)

        !=================================================================

    contains

        !===============================================================

        subroutine totalSubsonicInlet
            !
            !         TotalSubsonicInlet converts the prescribed total
            !         conditions and velocity direction into a useable format.
            !
            use constants
            use communication, only: adflow_comm_world
            use inputPhysics, only: RGasDim
            use section, only: sections
            implicit none
            !
            !        Local variables.
            !
            integer(kind=intType) :: i, j, nn

            real(kind=realType) :: rhot, mult, trans, Hdim, Tdim
            real(kind=realType) :: ax, r1, r2, var, wax, wrad, wtheta

            real(kind=realType), dimension(3) :: xc, dir

            integer :: ierr

            ! Set the subsonic inlet treatment to totalConditions.

            BCData(boco)%subsonicInletTreatment = totalConditions

            ! If the total pressure is present, convert it to SI-units and
            ! store it.

            if (ptPresent) then
                call siPressure(mass(1), length(1), time(1), mult, trans)

                do j = jBeg, jEnd
                    do i = iBeg, iEnd
                        BCData(boco)%ptInlet(i, j) = (mult * bcVarArray(i, j, 1) &
                                                      + trans) / Pref
                    end do
                end do
            end if

            ! If the total temperature is present, convert it to SI-units
            ! and store it.

            if (ttPresent) then
                call siTemperature(temp(2), mult, trans)

                do j = jBeg, jEnd
                    do i = iBeg, iEnd
                        BCData(boco)%ttInlet(i, j) = (mult * bcVarArray(i, j, 2) &
                                                      + trans) / Tref
                    end do
                end do
            end if

            ! Check if the total density is present. If so, it may be used
            ! to determine the total temperature or pressure if one of these
            ! variables was not specified.

            if (rhotPresent) then
                call siDensity(mass(3), length(3), mult, trans)

                if (ptPresent .and. (.not. ttPresent)) then

                    ! Total pressure is present but total temperature is not.
                    ! Convert the total density to SI-units and use the perfect
                    ! gas law to obtain the total temperature.

                    do j = jBeg, jEnd
                        do i = iBeg, iEnd
                            rhot = mult * bcVarArray(i, j, 3) + trans
                            BCData(boco)%ttInlet(i, j) = &
                                (BCData(boco)%ptInlet(i, j) * pRef / (RGasDim * rhot)) / Tref

                        end do
                    end do

                else if (ttPresent .and. (.not. ptPresent)) then

                    ! Total temperature is present but total pressure is not.
                    ! Convert the total density to SI-units and use the perfect
                    ! gas law to obtain the total pressure.

                    do j = jBeg, jEnd
                        do i = iBeg, iEnd
                            rhot = mult * bcVarArray(i, j, 3) + trans

                            BCData(boco)%ptInlet(i, j) = (RGasDim * rhot &
                                                          * BCData(boco)%ttInlet(i, j) * Tref) / Pref
                        end do
                    end do

                end if
            end if

            ! Determine the velocity direction. There are multiple
            ! possibilities to specify this direction.

            radialTest: if (rdirPresent) then

                ! Radial direction specified, i.e. a cylindrical coordinate
                ! system is used for the velocity direction.

                ! Determine the unit vectors, which define the cylindrical
                ! coordinate system aligned with the rotation axis.

                call unitVectorsCylSystem(boco)

                ! Initialize wtheta to zero. This value will be used if no
                ! theta velocity component was specified.

                wtheta = zero

                ! Loop over the faces of the subface.

                do j = jBeg, jEnd
                    do i = iBeg, iEnd

                        ! Determine the coordinates of the face center relative to
                        ! the rotation point of this section. Normally this is an
                        ! average of i-1, i, j-1, j, but due to the usage of the
                        ! pointer xf and the fact that x originally starts at 0,
                        ! an offset of 1 is introduced and thus the average should
                        ! be taken of i, i+1, j and j+1.

                        xc(1) = fourth * (xf(i, j, 1) + xf(i + 1, j, 1) &
                                          + xf(i, j + 1, 1) + xf(i + 1, j + 1, 1)) &
                                - sections(sectionId)%rotCenter(1)
                        xc(2) = fourth * (xf(i, j, 2) + xf(i + 1, j, 2) &
                                          + xf(i, j + 1, 2) + xf(i + 1, j + 1, 2)) &
                                - sections(sectionId)%rotCenter(2)
                        xc(3) = fourth * (xf(i, j, 3) + xf(i + 1, j, 3) &
                                          + xf(i, j + 1, 3) + xf(i + 1, j + 1, 3)) &
                                - sections(sectionId)%rotCenter(3)

                        ! Determine the coordinates in the local cartesian frame,
                        ! i.e. the frame determined by axis, radVec1 and radVec2.

                        ax = xc(1) * axis(1) + xc(2) * axis(2) &
                             + xc(3) * axis(3)
                        r1 = xc(1) * radVec1(1) + xc(2) * radVec1(2) &
                             + xc(3) * radVec1(3)
                        r2 = xc(1) * radVec2(1) + xc(2) * radVec2(2) &
                             + xc(3) * radVec2(3)

                        ! Determine the weights of the unit vectors in the local
                        ! cylindrical system.

                        wax = bcVarArray(i, j, 7)
                        wrad = bcVarArray(i, j, 10)
                        if (tdirPresent) wtheta = bcVarArray(i, j, 11)

                        ! Determine the direction in the local cartesian frame,
                        ! determined by axis, radVec1 and radVec2.

                        var = one / sqrt(max(eps, (r1 * r1 + r2 * r2)))
                        dir(1) = wax
                        dir(2) = var * (wrad * r1 - wtheta * r2)
                        dir(3) = var * (wrad * r2 + wtheta * r1)

                        ! Transform this direction to the global cartesian frame.

                        BCData(boco)%flowXdirInlet(i, j) = dir(1) * axis(1) &
                                                           + dir(2) * radVec1(1) &
                                                           + dir(3) * radVec2(1)

                        BCData(boco)%flowYdirInlet(i, j) = dir(1) * axis(2) &
                                                           + dir(2) * radVec1(2) &
                                                           + dir(3) * radVec2(2)

                        BCData(boco)%flowZdirInlet(i, j) = dir(1) * axis(3) &
                                                           + dir(2) * radVec1(3) &
                                                           + dir(3) * radVec2(3)
                    end do
                end do

                else radialTest

                ! Cartesian direction specified. Either the angle or the
                ! direction should be present.

                ! X-direction.

                if (axPresent) then

                    ! Angle specified. Convert it to SI-units and determine
                    ! the corresponding direction.

                    call siAngle(angle(4), mult, trans)

                    do j = jBeg, jEnd
                        do i = iBeg, iEnd
                            BCData(boco)%flowXdirInlet(i, j) = &
                                cos(mult * bcVarArray(i, j, 4) + trans)
                        end do
                    end do

                else

                    ! Direction specified. Simply copy it.

                    do j = jBeg, jEnd
                        do i = iBeg, iEnd
                            BCData(boco)%flowXdirInlet(i, j) = bcVarArray(i, j, 7)
                        end do
                    end do

                end if

                ! Y-direction.

                if (ayPresent) then

                    ! Angle specified. Convert it to SI-units and determine
                    ! the corresponding direction.

                    call siAngle(angle(5), mult, trans)

                    do j = jBeg, jEnd
                        do i = iBeg, iEnd
                            BCData(boco)%flowYdirInlet(i, j) = &
                                cos(mult * bcVarArray(i, j, 5) + trans)
                        end do
                    end do

                else

                    ! Direction specified. Simply copy it.

                    do j = jBeg, jEnd
                        do i = iBeg, iEnd
                            BCData(boco)%flowYdirInlet(i, j) = bcVarArray(i, j, 8)
                        end do
                    end do

                end if

                ! Z-direction.

                if (azPresent) then

                    ! Angle specified. Convert it to SI-units and determine
                    ! the corresponding direction.

                    call siAngle(angle(6), mult, trans)

                    do j = jBeg, jEnd
                        do i = iBeg, iEnd
                            BCData(boco)%flowZdirInlet(i, j) = &
                                cos(mult * bcVarArray(i, j, 6) + trans)
                        end do
                    end do

                else

                    ! Direction specified. Simply copy it.

                    do j = jBeg, jEnd
                        do i = iBeg, iEnd
                            BCData(boco)%flowZdirInlet(i, j) = bcVarArray(i, j, 9)
                        end do
                    end do

                end if

            end if radialTest

            ! Loop over the faces of the subface to compute some
            ! additional info.

            do j = jBeg, jEnd
                do i = iBeg, iEnd

                    ! Compute the total enthalpy from the given
                    ! total temperature.
                    TDim = BCData(boco)%ttInlet(i, j) * Tref
                    call computeHtot(TDim, Hdim)
                    BCData(boco)%htInlet(i, j) = Hdim / Href

                    ! Determine the unit vector of the flow direction.

                    dir(1) = BCData(boco)%flowXdirInlet(i, j)
                    dir(2) = BCData(boco)%flowYdirInlet(i, j)
                    dir(3) = BCData(boco)%flowZdirInlet(i, j)

                    var = one / max(eps, sqrt(dir(1)**2 + dir(2)**2 + dir(3)**2))

                    BCData(boco)%flowXdirInlet(i, j) = var * dir(1)
                    BCData(boco)%flowYdirInlet(i, j) = var * dir(2)
                    BCData(boco)%flowZdirInlet(i, j) = var * dir(3)

                end do
            end do

            ! Check if the prescribed direction is an inflow. No halo's
            ! should be included here and therefore the nodal range
            ! (with an offset) must be used.

            nn = 0
            do j = (BCData(boco)%jnbeg + 1), BCData(boco)%jnend
                do i = (BCData(boco)%inbeg + 1), BCData(boco)%inend

                    var = BCData(boco)%flowXdirInlet(i, j) &
                          * BCData(boco)%norm(i, j, 1) &
                          + BCData(boco)%flowYdirInlet(i, j) &
                          * BCData(boco)%norm(i, j, 2) &
                          + BCData(boco)%flowZdirInlet(i, j) &
                          * BCData(boco)%norm(i, j, 3)

                    if (var > zero) nn = nn + 1

                end do
            end do
        end subroutine totalSubsonicInlet

    end subroutine BCDataSubsonicInflow

    subroutine BCDataSubsonicOutflow(boco, bcVarArray, iBeg, iEnd, jBeg, jEnd)
        !
        !       BCDataSubsonicOutflow tries to extract the static pressure
        !       for the currently active boundary face, which is a subsonic
        !       outflow boundary.
        !
        use constants
        use cgnsNames
        use blockPointers, only: BCData, nbkGlobal, BCFaceID
        use utils, only: terminate, siPressure
        use flowVarRefState, only: pRef
        implicit none
        !
        !      Subroutine arguments.
        !
        integer(kind=intType) :: boco
        integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
        real(kind=realType), dimension(iBeg:iEnd, jBeg:jEnd, nbcVarMax) :: bcVarArray
        !
        !      Local variables.
        !
        integer :: ierr

        integer(kind=intType) :: i, j

        real(kind=realType) :: mult, trans

        character(len=maxStringLen) :: errorMessage

        ! Write an error message and terminate if it was not
        ! possible to determine the static pressure.


        ! Convert to SI-units and store the pressure in ps.

        call siPressure(mass(1), length(1), time(1), mult, trans)
        do j = jBeg, jEnd
            do i = iBeg, iEnd
                BCData(boco)%ps(i, j) = (mult * bcVarArray(i, j, 1) + trans) / Pref
            end do
        end do

    end subroutine BCDataSubsonicOutflow

    subroutine BCDataSupersonicInflow(boco, bcVarArray, iBeg, iEnd, jBeg, jEnd, &
                                      allFlowPresent, allTurbPresent)
        !
        !       BCDataSupersonicInflow tries to extract the primitive state
        !       vector for the currently active boundary face, which is a
        !       supersonic inflow.
        !
        use constants
        use cgnsNames
        use blockPointers, only: BCData, nbkGlobal, BCFaceID, sectionID
        use flowVarRefState, only: nwt, pInfCorr, wInf, uRef, rhoRef, pRef, muRef
        use inputPhysics, onlY: equations, flowType, velDirFreeStream
        use utils, only: siDensity, siPressure, siVelocity, siTemperature, terminate
        implicit none
        !
        !      Subroutine arguments.
        !
        integer(kind=intType), intent(in) :: boco
        integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
        real(kind=realType), dimension(iBeg:iEnd, jBeg:jEnd, nbcVarMax) :: bcVarArray
        logical, intent(inout) :: allFlowPresent
        logical, intent(inout) :: allTurbPresent
        !
        !      Local variables.
        !
        integer :: ierr

        integer(kind=intType) :: i, j, nn

        real(kind=realType) :: var

        character(len=maxStringLen) :: errorMessage

        logical :: rhoPresent, pPresent, velPresent
        logical :: velxPresent, velyPresent, velzPresent
        logical :: velrPresent, veltPresent

        ! Store the logicals, which indicate success or failure
        ! a bit more readable.

        rhoPresent = bcVarPresent(1)
        pPresent = bcVarPresent(2)
        velxPresent = bcVarPresent(3)
        velyPresent = bcVarPresent(4)
        velzPresent = bcVarPresent(5)
        velrPresent = bcVarPresent(6)
        veltPresent = bcVarPresent(7)

        ! Check if a velocity vector is present.

        velPresent = .false.
        if (velxPresent .and. velrPresent) velPresent = .true.
        if (velxPresent .and. velyPresent .and. velzPresent) &
            velPresent = .true.

        ! Check if rho, p and the velocity vector are present.

        testPresent: if (rhoPresent .and. pPresent .and. velPresent) then

            ! All the variables needed are prescribed. Set them.

            call prescribedSupersonicInlet

            else testPresent

            ! Not all variables are present. Check what type of flow
            ! is to be solved.

            select case (flowType)

            case (internalFlow)

                ! Internal flow. Data at the inlet must be specified;
                ! no free stream data can be taken.


                !=============================================================

            case (externalFlow)

                ! External flow. Free stream data is used.

                do j = jBeg, jEnd
                    do i = iBeg, iEnd
                        BCData(boco)%rho(i, j) = wInf(iRho)
                        BCData(boco)%velx(i, j) = wInf(ivx)
                        BCData(boco)%vely(i, j) = wInf(ivy)
                        BCData(boco)%velz(i, j) = wInf(ivz)
                        BCData(boco)%ps(i, j) = PinfCorr
                    end do
                end do

                ! Set the turbulence values
                allTurbPresent = setBCVarTurb(7_intType, boco, bcVarArray, &
                                              iBeg, iEnd, jBeg, jEnd, BCData(boco)%turbInlet)

                ! Set allFlowPresent to .false.

                allFlowPresent = .false.

            end select

        end if testPresent


    contains

        subroutine prescribedSupersonicInlet
            !
            !         prescribedSupersonicInlet sets the variables for this
            !         supersonic inlet to prescribed values.
            !
            use section, only: sections
            implicit none
            !
            !        Local variables.
            !
            integer(kind=intType) :: i, j

            real(kind=realType) :: mult, trans
            real(kind=realType) :: ax, r1, r2, var, vax, vrad, vtheta

            real(kind=realType), dimension(3) :: xc, vloc
            real(kind=realType), dimension(3) :: multVel, transVel

            ! Set the density. Take the conversion factor to SI-units
            ! into account.

            call siDensity(mass(1), length(1), mult, trans)

            do j = jBeg, jEnd
                do i = iBeg, iEnd
                    BCData(boco)%rho(i, j) = (mult * bcVarArray(i, j, 1) + trans) / rhoRef
                end do
            end do

            ! Set the pressure. Take the conversion factor to SI-units
            ! into account.

            call siPressure(mass(1), length(2), time(2), mult, trans)

            do j = jBeg, jEnd
                do i = iBeg, iEnd
                    BCData(boco)%ps(i, j) = (mult * bcVarArray(i, j, 2) + trans) / pRef
                end do
            end do

            ! Check the situation we are having here for the velocity.

            testRadial: if (velrPresent) then

                ! Radial velocity component prescribed. This must be converted
                ! to cartesian components.

                ! Determine the unit vectors, which define the cylindrical
                ! coordinate system aligned with the rotation axis.

                call unitVectorsCylSystem(boco)

                ! Determine the conversion factor to SI-units for the three
                ! components. Note that a test must be made whether the theta
                ! component is present.

                call siVelocity(length(3), time(3), multVel(1), transVel(1))
                call siVelocity(length(6), time(6), multVel(2), transVel(2))

                if (veltPresent) &
                    call siVelocity(length(7), time(7), multVel(3), transVel(3))

                ! Initialize vtheta to zero. This value will be used
                ! if no theta velocity component was specified.

                vtheta = zero

                ! Loop over the faces of the subface.

                do j = jBeg, jEnd
                    do i = iBeg, iEnd

                        ! Determine the coordinates of the face center relative to
                        ! the rotation point of this section. Normally this is an
                        ! average of i-1, i, j-1, j, but due to the usage of the
                        ! pointer xf and the fact that x originally starts at 0,
                        ! an offset of 1 is introduced and thus the average should
                        ! be taken of i, i+1, j and j+1.

                        xc(1) = fourth * (xf(i, j, 1) + xf(i + 1, j, 1) &
                                          + xf(i, j + 1, 1) + xf(i + 1, j + 1, 1)) &
                                - sections(sectionID)%rotCenter(1)
                        xc(2) = fourth * (xf(i, j, 2) + xf(i + 1, j, 2) &
                                          + xf(i, j + 1, 2) + xf(i + 1, j + 1, 2)) &
                                - sections(sectionID)%rotCenter(2)
                        xc(3) = fourth * (xf(i, j, 3) + xf(i + 1, j, 3) &
                                          + xf(i, j + 1, 3) + xf(i + 1, j + 1, 3)) &
                                - sections(sectionID)%rotCenter(3)

                        ! Determine the coordinates in the local cartesian frame,
                        ! i.e. the frame determined by axis, radVec1 and radVec2.

                        ax = xc(1) * axis(1) + xc(2) * axis(2) &
                             + xc(3) * axis(3)
                        r1 = xc(1) * radVec1(1) + xc(2) * radVec1(2) &
                             + xc(3) * radVec1(3)
                        r2 = xc(1) * radVec2(1) + xc(2) * radVec2(2) &
                             + xc(3) * radVec2(3)

                        ! Determine the velocity components in the local
                        ! cylindrical system. Take the conversion to si units
                        ! into account.

                        vax = multVel(1) * bcVarArray(i, j, 3) + transVel(1)
                        vrad = multVel(2) * bcVarArray(i, j, 6) + transVel(2)
                        if (veltPresent) &
                            vtheta = multVel(3) * bcVarArray(i, j, 7) + transVel(3)

                        ! Determine the velocities in the local cartesian
                        ! frame determined by axis, radVec1 and radVec2.

                        var = one / sqrt(max(eps, (r1 * r1 + r2 * r2)))
                        vloc(1) = vax
                        vloc(2) = var * (vrad * r1 - vtheta * r2)
                        vloc(3) = var * (vrad * r2 + vtheta * r1)

                        ! Transform vloc to the global cartesian frame and
                        ! store the values.

                        BCData(boco)%velx(i, j) = (vloc(1) * axis(1) &
                                                   + vloc(2) * radVec1(1) &
                                                   + vloc(3) * radVec2(1)) / uRef

                        BCData(boco)%vely(i, j) = (vloc(1) * axis(2) &
                                                   + vloc(2) * radVec1(2) &
                                                   + vloc(3) * radVec2(2)) / uRef

                        BCData(boco)%velz(i, j) = (vloc(1) * axis(3) &
                                                   + vloc(2) * radVec1(3) &
                                                   + vloc(3) * radVec2(3)) / uRef
                    end do
                end do

                else testRadial

                ! Cartesian components prescribed.

                ! Determine the conversion factor to SI-units for the three
                ! components.

                call siVelocity(length(3), time(3), multVel(1), transVel(1))
                call siVelocity(length(4), time(4), multVel(2), transVel(2))
                call siVelocity(length(5), time(5), multVel(3), transVel(3))

                ! Set the velocities.

                do j = jBeg, jEnd
                    do i = iBeg, iEnd
                        BCData(boco)%velx(i, j) = (multVel(1) * bcVarArray(i, j, 3) &
                                                   + transVel(1)) / uRef
                        BCData(boco)%vely(i, j) = (multVel(2) * bcVarArray(i, j, 4) &
                                                   + transVel(2)) / uRef
                        BCData(boco)%velz(i, j) = (multVel(3) * bcVarArray(i, j, 5) &
                                                   + transVel(3)) / uRef
                    end do
                end do

            end if testRadial

            ! Set the turbulence variables and check if all of them are
            ! prescribed. If not set allTurbPresent to .false.

            allTurbPresent = setBCVarTurb(7_intType, boco, bcVarArray, &
                                          iBeg, iEnd, jBeg, jEnd, BCData(boco)%turbInlet)

        end subroutine prescribedSupersonicInlet

    end subroutine BCDataSupersonicInflow

    !=================================================================

    logical function setBCVarTurb(offset, boco, bcVarArray, &
                                  iBeg, iEnd, jBeg, jEnd, turbInlet)
        !
        !       SetBCVarTurb sets the array for the turbulent halo data
        !       for inlet boundaries. This function returns .true. If all
        !       turbulence variables could be interpolated and .false.
        !       otherwise.
        !
        use constants
        use flowVarRefState, only: nt1, nt2, muRef, Pref, rhoRef, wInf
        use inputPhysics, only: equations, turbModel
        use utils, only: terminate, siTurb

        implicit none
        !
        !      Subroutine arguments.
        !
        integer(kind=intType), intent(in) :: offset, boco, iBeg, iEnd, jBeg, jEnd
        real(kind=realType), dimension(iBeg:iEnd, jBeg:jEnd, nbcVarMax) :: bcVarArray
        real(kind=realType), dimension(:, :, :), pointer :: turbInlet
        !
        !      Local variables.
        !
        integer(kind=intType) :: nn, mm, i, j
        real(kind=realType) :: mult, trans, nuRef
        real(kind=realType), dimension(nt1:nt2) :: ref

        ! Initialize setBCVarTurb to .true. And return immediately
        ! if not the rans equations are solved.

        setBCVarTurb = .true.
        if (equations /= RANSEquations) return

        ! Set the reference values depending on the turbulence model.

        nuRef = muRef / rhoRef
        select case (turbModel)

        case (spalartAllmaras, spalartAllmarasEdwards)
            ref(itu1) = nuRef

        case (spalartallmarasnoft2gammaretheta)
            ref(itu1) = nuRef
            ref(itu2) = one
            ref(itu3) = one

        case (komegaWilcox, komegaModified, menterSST)
            ref(itu1) = pRef / rhoRef
            ref(itu2) = ref(itu1) / nuRef

        case (ktau)
            ref(itu1) = pRef / rhoRef
            ref(itu2) = nuRef / ref(itu1)

        case (v2f)
            ref(itu1) = pRef / rhoRef
            ref(itu4) = ref(itu1) / nuRef
            ref(itu2) = ref(itu1) * ref(itu4)
            ref(itu3) = ref(itu1)

        end select

        ! Loop over the number of turbulent variables. mm is the counter
        ! in the arrays bcVarArray and bcVarPresent.

        mm = offset
        turbLoop: do nn = nt1, nt2
            mm = mm + 1

            ! Check if the variable is present. If so, use the
            ! interpolated data.

            if (bcVarPresent(mm)) then

                ! Conversion to SI units if possible.

                call siTurb(mass(mm), length(mm), time(mm), temp(mm), &
                            bcVarNames(mm), mult, trans)

                ! Set the turbulent variables.

                do j = jBeg, jEnd
                    do i = iBeg, iEnd
                        turbInlet(i, j, nn) = (mult * bcVarArray(i, j, mm) + trans) / ref(nn)
                    end do
                end do

            else

                ! Turbulent variable not present. Use the free stream data.
                do j = jBeg, jEnd
                    do i = iBeg, iEnd
                        turbInlet(i, j, nn) = wInf(nn)
                    end do
                end do

                ! Set the logical value to false to indicate that indeed not
                ! all the values were present
                setBCVarTurb = .false.

            end if
        end do turbLoop
    end function setBCVarTurb

end module BCData
