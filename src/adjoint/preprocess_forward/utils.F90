










module utils
    implicit none

contains

    function char2str(charArray, n)
        use constants
        !
        ! some gymnastics to cast a char array to string
        !
        implicit none
        !
        !      Function arguments.
        !
        character, dimension(maxCGNSNameLen), intent(in) :: charArray
        integer(kind=intType), intent(in) :: n
        !
        !      Function type
        !
        character(len=n) :: char2str
        !
        !      Local variables.
        !
        integer(kind=intType) :: i
        do i = 1, n
            char2str(i:i) = charArray(i)
        end do

    end function char2str

    function TSbeta(degreePolBeta, coefPolBeta, &
                    degreeFourBeta, omegaFourBeta, &
                    cosCoefFourBeta, sinCoefFourBeta, t)
        !
        !       TSbeta computes the angle of attack for a given Time interval
        !       in a time spectral solution.
        !
        use constants
        use inputPhysics, only: equationMode
        implicit none
        !
        !      Function type
        !
        real(kind=realType) :: TSbeta
        !
        !      Function arguments.
        !
        integer(kind=intType), intent(in) :: degreePolBeta
        integer(kind=intType), intent(in) :: degreeFourBeta

        real(kind=realType), intent(in) :: omegaFourBeta, t

        real(kind=realType), dimension(0:*), intent(in) :: coefPolBeta
        real(kind=realType), dimension(0:*), intent(in) :: cosCoefFourBeta
        real(kind=realType), dimension(*), intent(in) :: sinCoefFourBeta
        !
        !      Local variables.
        !
        integer(kind=intType) :: nn

        real(kind=realType) :: beta, val

        ! Return immediately if this is a steady computation.

        if (equationMode == steady) then
            TSBeta = zero
            return
        end if

        ! Compute the polynomial contribution. If no polynomial was
        ! specified, the value of index 0 is set to zero automatically.

        beta = coefPolBeta(0)
        do nn = 1, degreePolBeta
            beta = beta + coefPolBeta(nn) * (t**nn)
        end do

        ! Compute the fourier contribution. Again the cosine coefficient
        ! of index 0 is defaulted to zero if not specified.

        beta = beta + cosCoefFourBeta(0)
        do nn = 1, degreeFourBeta
            val = nn * omegaFourBeta * t
            beta = beta + cosCoefFourbeta(nn) * cos(val) &
                   + sinCoefFourbeta(nn) * sin(val)
        end do

        ! Set TSBeta to phi.

        TSBeta = beta

    end function TSbeta

    function TSbetadot(degreePolBeta, coefPolBeta, &
                       degreeFourBeta, omegaFourBeta, &
                       cosCoefFourBeta, sinCoefFourBeta, t)
        !
        !       TSbeta computes the angle of attack for a given Time interval
        !       in a time spectral solution.
        !
        use constants
        use inputPhysics, only: equationMode
        implicit none
        !
        !      Function type
        !
        real(kind=realType) :: TSbetadot
        !
        !      Function arguments.
        !
        integer(kind=intType), intent(in) :: degreePolBeta
        integer(kind=intType), intent(in) :: degreeFourBeta

        real(kind=realType), intent(in) :: omegaFourBeta, t

        real(kind=realType), dimension(0:*), intent(in) :: coefPolBeta
        real(kind=realType), dimension(0:*), intent(in) :: cosCoefFourBeta
        real(kind=realType), dimension(*), intent(in) :: sinCoefFourBeta
        !
        !      Local variables.
        !
        integer(kind=intType) :: nn

        real(kind=realType) :: betadot, val

        ! Return immediately if this is a steady computation.

        if (equationMode == steady) then
            TSBetadot = zero
            return
        end if

        ! Compute the polynomial contribution. If no polynomial was
        ! specified, the value of index 0 is set to zero automatically.

        betadot = zero
        do nn = 1, degreePolBeta
            betadot = betadot + nn * coefPolBeta(nn) * (t**(nn - 1))
        end do

        ! Compute the fourier contribution. Again the cosine coefficient
        ! of index 0 is defaulted to zero if not specified.

        do nn = 1, degreeFourBeta
            val = nn * omegaFourBeta
            betadot = betadot - val * cosCoefFourbeta(nn) * sin(val * t) &
                      + val * sinCoefFourbeta(nn) * cos(val * t)
        end do

        ! Set TSBeta to phi.

        TSBetadot = betadot

    end function TSbetadot

    function TSMach(degreePolMach, coefPolMach, &
                    degreeFourMach, omegaFourMach, &
                    cosCoefFourMach, sinCoefFourMach, t)
        !
        !       TSMach computes the Mach Number for a given time interval
        !       in a time spectral solution.
        !
        use constants
        use inputPhysics, only: equationMode
        implicit none
        !
        !      Function type
        !
        real(kind=realType) :: TSmach
        !
        !      Function arguments.
        !
        integer(kind=intType), intent(in) :: degreePolMach
        integer(kind=intType), intent(in) :: degreeFourMach

        real(kind=realType), intent(in) :: omegaFourMach, t

        real(kind=realType), dimension(0:*), intent(in) :: coefPolMach
        real(kind=realType), dimension(0:*), intent(in) :: cosCoefFourMach
        real(kind=realType), dimension(*), intent(in) :: sinCoefFourMach
        !
        !      Local variables.
        !
        integer(kind=intType) :: nn

        real(kind=realType) :: intervalMach, val

        ! Return immediately if this is a steady computation.

        if (equationMode == steady) then
            TSMach = zero
            return
        end if

        ! Compute the polynomial contribution. If no polynomial was
        ! specified, the value of index 0 is set to zero automatically.

        intervalMach = coefPolMach(0)
        do nn = 1, degreePolMach
            intervalMach = intervalMach + coefPolMach(nn) * (t**nn)
        end do

        ! Compute the fourier contribution. Again the cosine coefficient
        ! of index 0 is defaulted to zero if not specified.

        intervalMach = intervalMach + cosCoefFourMach(0)
        do nn = 1, degreeFourMach
            val = nn * omegaFourMach * t
            intervalMach = intervalMach + cosCoefFourmach(nn) * cos(val) &
                           + sinCoefFourmach(nn) * sin(val)
        end do
        print *, 'inTSMach', intervalMach, nn, val, t
        ! Set TSMach to phi.

        TSMach = intervalMach

    end function TSmach

    function TSMachdot(degreePolMach, coefPolMach, &
                       degreeFourMach, omegaFourMach, &
                       cosCoefFourMach, sinCoefFourMach, t)
        !
        !       TSmach computes the angle of attack for a given Time interval
        !       in a time spectral solution.
        !
        use constants
        use inputPhysics, only: equationMode
        implicit none
        !
        !      Function type
        !
        real(kind=realType) :: TSmachdot
        !
        !      Function arguments.
        !
        integer(kind=intType), intent(in) :: degreePolMach
        integer(kind=intType), intent(in) :: degreeFourMach

        real(kind=realType), intent(in) :: omegaFourMach, t

        real(kind=realType), dimension(0:*), intent(in) :: coefPolMach
        real(kind=realType), dimension(0:*), intent(in) :: cosCoefFourMach
        real(kind=realType), dimension(*), intent(in) :: sinCoefFourMach
        !
        !      Local variables.
        !
        integer(kind=intType) :: nn

        real(kind=realType) :: machdot, val

        ! Return immediately if this is a steady computation.

        if (equationMode == steady) then
            TSMachdot = zero
            return
        end if

        ! Compute the polynomial contribution. If no polynomial was
        ! specified, the value of index 0 is set to zero automatically.

        machdot = zero
        do nn = 1, degreePolMach
            machdot = machdot + nn * coefPolMach(nn) * (t**(nn - 1))
        end do

        ! Compute the fourier contribution. Again the cosine coefficient
        ! of index 0 is defaulted to zero if not specified.

        do nn = 1, degreeFourMach
            val = nn * omegaFourMach
            machdot = machdot - val * cosCoefFourmach(nn) * sin(val * t) &
                      + val * sinCoefFourmach(nn) * cos(val * t)
        end do

        ! Set TSMach to phi.

        TSMachdot = machdot

    end function TSmachdot

    function TSalpha(degreePolAlpha, coefPolAlpha, &
                     degreeFourAlpha, omegaFourAlpha, &
                     cosCoefFourAlpha, sinCoefFourAlpha, t)
        !
        !       TSalpha computes the angle of attack for a given Time interval
        !       in a time spectral solution.
        !
        use constants
        use inputPhysics, only: equationMode
        implicit none
        !
        !      Function type
        !
        real(kind=realType) :: TSalpha
        !
        !      Function arguments.
        !
        integer(kind=intType), intent(in) :: degreePolAlpha
        integer(kind=intType), intent(in) :: degreeFourAlpha

        real(kind=realType), intent(in) :: omegaFourAlpha, t

        real(kind=realType), dimension(0:*), intent(in) :: coefPolAlpha
        real(kind=realType), dimension(0:*), intent(in) :: cosCoefFourAlpha
        real(kind=realType), dimension(*), intent(in) :: sinCoefFourAlpha
        !
        !      Local variables.
        !
        integer(kind=intType) :: nn

        real(kind=realType) :: alpha, val

        ! Return immediately if this is a steady computation.

        if (equationMode == steady) then
            TSAlpha = zero
            return
        end if

        ! Compute the polynomial contribution. If no polynomial was
        ! specified, the value of index 0 is set to zero automatically.
        alpha = coefPolAlpha(0)
        do nn = 1, degreePolAlpha
            alpha = alpha + coefPolAlpha(nn) * (t**nn)
        end do

        ! Compute the fourier contribution. Again the cosine coefficient
        ! of index 0 is defaulted to zero if not specified.

        alpha = alpha + cosCoefFourAlpha(0)
        do nn = 1, degreeFourAlpha
            val = nn * omegaFourAlpha * t
            alpha = alpha + cosCoefFouralpha(nn) * cos(val) &
                    + sinCoefFouralpha(nn) * sin(val)
        end do
        !print *,'inTSalpha',alpha,nn,val,t
        ! Set TSAlpha to phi.

        TSAlpha = alpha

    end function TSalpha

    function TSalphadot(degreePolAlpha, coefPolAlpha, &
                        degreeFourAlpha, omegaFourAlpha, &
                        cosCoefFourAlpha, sinCoefFourAlpha, t)
        !
        !       TSalpha computes the angle of attack for a given Time interval
        !       in a time spectral solution.
        !
        use constants
        use inputPhysics, only: equationMode
        implicit none
        !
        !      Function type
        !
        real(kind=realType) :: TSalphadot
        !
        !      Function arguments.
        !
        integer(kind=intType), intent(in) :: degreePolAlpha
        integer(kind=intType), intent(in) :: degreeFourAlpha

        real(kind=realType), intent(in) :: omegaFourAlpha, t

        real(kind=realType), dimension(0:*), intent(in) :: coefPolAlpha
        real(kind=realType), dimension(0:*), intent(in) :: cosCoefFourAlpha
        real(kind=realType), dimension(*), intent(in) :: sinCoefFourAlpha
        !
        !      Local variables.
        !
        integer(kind=intType) :: nn

        real(kind=realType) :: alphadot, val

        ! Return immediately if this is a steady computation.

        if (equationMode == steady) then
            TSAlphadot = zero
            return
        end if

        ! Compute the polynomial contribution. If no polynomial was
        ! specified, the value of index 0 is set to zero automatically.

        alphadot = zero
        do nn = 1, degreePolAlpha
            alphadot = alphadot + nn * coefPolAlpha(nn) * (t**(nn - 1))
        end do

        ! Compute the fourier contribution. Again the cosine coefficient
        ! of index 0 is defaulted to zero if not specified.

        do nn = 1, degreeFourAlpha
            val = nn * omegaFourAlpha
            alphadot = alphadot - val * cosCoefFouralpha(nn) * sin(val * t) &
                       + val * sinCoefFouralpha(nn) * cos(val * t)
        end do

        ! Set TSAlpha to phi.

        TSAlphadot = alphadot

    end function TSalphadot

    function derivativeRigidRotAngle(degreePolRot, &
                                     coefPolRot, &
                                     degreeFourRot, &
                                     omegaFourRot, &
                                     cosCoefFourRot, &
                                     sinCoefFourRot, t)
        !
        !       derivativeRigidRotAngle computes the time derivative of the
        !       rigid body rotation angle at the given time for the given
        !       arguments. The angle is described by a combination of a
        !       polynomial and fourier series.
        !
        use constants
        use inputPhysics, only: equationMode
        use flowVarRefState, only: timeRef
        implicit none
        !
        !      Function type
        !
        real(kind=realType) :: derivativeRigidRotAngle
        !
        !      Function arguments.
        !
        integer(kind=intType), intent(in) :: degreePolRot
        integer(kind=intType), intent(in) :: degreeFourRot

        real(kind=realType), intent(in) :: omegaFourRot, t

        real(kind=realType), dimension(0:*), intent(in) :: coefPolRot
        real(kind=realType), dimension(0:*), intent(in) :: cosCoefFourRot
        real(kind=realType), dimension(*), intent(in) :: sinCoefFourRot
        !
        !      Local variables.
        !
        integer(kind=intType) :: nn

        real(kind=realType) :: dPhi, val

        ! Return immediately if this is a steady computation.

        if (equationMode == steady) then
            derivativeRigidRotAngle = zero
            return
        end if

        ! Compute the polynomial contribution.

        dPhi = zero
        do nn = 1, degreePolRot
            dPhi = dPhi + nn * coefPolRot(nn) * (t**(nn - 1))
        end do

        ! Compute the fourier contribution.

        do nn = 1, degreeFourRot
            val = nn * omegaFourRot
            dPhi = dPhi - val * cosCoefFourRot(nn) * sin(val * t)
            dPhi = dPhi + val * sinCoefFourRot(nn) * cos(val * t)
        end do

        ! Set derivativeRigidRotAngle to dPhi. Multiply by timeRef
        ! to obtain the correct non-dimensional value.

        derivativeRigidRotAngle = timeRef * dPhi

    end function derivativeRigidRotAngle

    function myDim(x, y)

        use constants

        real(kind=realType) x, y
        real(kind=realType) :: myDim

        myDim = x - y
        if (myDim < 0.0) then
            myDim = 0.0
        end if

    end function myDim

    function getCorrectForK()

        use constants
        use flowVarRefState, only: kPresent
        use iteration, only: currentLevel, groundLevel
        implicit none

        logical :: getCorrectForK

        if (kPresent .and. currentLevel <= groundLevel) then
            getCorrectForK = .true.
        else
            getCorrectForK = .false.
        end if
    end function getCorrectForK
    subroutine terminate(routineName, errorMessage)
        !
        !       terminate writes an error message to standard output and
        !       terminates the execution of the program.
        !
        use constants
        use communication, only: adflow_comm_world, myid
        implicit none
        !
        !      Subroutine arguments
        !
        character(len=*), intent(in) :: routineName
        character(len=*), intent(in) :: errorMessage

    end subroutine terminate

    subroutine rotMatrixRigidBody(tNew, tOld, rotationMatrix, &
                                  rotationPoint)
        !
        !       rotMatrixRigidBody determines the rotation matrix and the
        !       rotation point to determine the coordinates of the new time
        !       level starting from the coordinates of the old time level.
        !
        use constants
        use inputMotion
        use flowVarRefState, only: Lref
        implicit none
        !
        !      Subroutine arguments.
        !
        real(kind=realType), intent(in) :: tNew, tOld

        real(kind=realType), dimension(3), intent(out) :: rotationPoint
        real(kind=realType), dimension(3, 3), intent(out) :: rotationMatrix
        !
        !      Local variables.
        !
        integer(kind=intType) :: i, j

        real(kind=realType) :: phi
        real(kind=realType) :: cosX, cosY, cosZ, sinX, sinY, sinZ

        real(kind=realType), dimension(3, 3) :: mNew, mOld

        ! Determine the rotation angle around the x-axis for the new
        ! time level and the corresponding values of the sine and cosine.

        phi = rigidRotAngle(degreePolXRot, coefPolXRot, &
                            degreeFourXRot, omegaFourXRot, &
                            cosCoefFourXRot, sinCoefFourXRot, tNew)
        sinX = sin(phi)
        cosX = cos(phi)

        ! Idem for the y-axis.

        phi = rigidRotAngle(degreePolYRot, coefPolYRot, &
                            degreeFourYRot, omegaFourYRot, &
                            cosCoefFourYRot, sinCoefFourYRot, tNew)
        sinY = sin(phi)
        cosY = cos(phi)

        ! Idem for the z-axis.

        phi = rigidRotAngle(degreePolZRot, coefPolZRot, &
                            degreeFourZRot, omegaFourZRot, &
                            cosCoefFourZRot, sinCoefFourZRot, tNew)
        sinZ = sin(phi)
        cosZ = cos(phi)

        ! Construct the transformation matrix at the new time level.
        ! It is assumed that the sequence of rotation is first around the
        ! x-axis then around the y-axis and finally around the z-axis.

        mNew(1, 1) = cosY * cosZ
        mNew(2, 1) = cosY * sinZ
        mNew(3, 1) = -sinY

        mNew(1, 2) = sinX * sinY * cosZ - cosX * sinZ
        mNew(2, 2) = sinX * sinY * sinZ + cosX * cosZ
        mNew(3, 2) = sinX * cosY

        mNew(1, 3) = cosX * sinY * cosZ + sinX * sinZ
        mNew(2, 3) = cosX * sinY * sinZ - sinX * cosZ
        mNew(3, 3) = cosX * cosY

        ! Determine the rotation angle around the x-axis for the old
        ! time level and the corresponding values of the sine and cosine.

        phi = rigidRotAngle(degreePolXRot, coefPolXRot, &
                            degreeFourXRot, omegaFourXRot, &
                            cosCoefFourXRot, sinCoefFourXRot, tOld)
        sinX = sin(phi)
        cosX = cos(phi)

        ! Idem for the y-axis.

        phi = rigidRotAngle(degreePolYRot, coefPolYRot, &
                            degreeFourYRot, omegaFourYRot, &
                            cosCoefFourYRot, sinCoefFourYRot, tOld)
        sinY = sin(phi)
        cosY = cos(phi)

        ! Idem for the z-axis.

        phi = rigidRotAngle(degreePolZRot, coefPolZRot, &
                            degreeFourZRot, omegaFourZRot, &
                            cosCoefFourZRot, sinCoefFourZRot, tOld)
        sinZ = sin(phi)
        cosZ = cos(phi)

        ! Construct the transformation matrix at the old time level.

        mOld(1, 1) = cosY * cosZ
        mOld(2, 1) = cosY * sinZ
        mOld(3, 1) = -sinY

        mOld(1, 2) = sinX * sinY * cosZ - cosX * sinZ
        mOld(2, 2) = sinX * sinY * sinZ + cosX * cosZ
        mOld(3, 2) = sinX * cosY

        mOld(1, 3) = cosX * sinY * cosZ + sinX * sinZ
        mOld(2, 3) = cosX * sinY * sinZ - sinX * cosZ
        mOld(3, 3) = cosX * cosY

        ! Construct the transformation matrix between the new and the
        ! old time level. This is mNew*inverse(mOld). However the
        ! inverse of mOld is the transpose.

        do j = 1, 3
            do i = 1, 3
                rotationMatrix(i, j) = mNew(i, 1) * mOld(j, 1) &
                                       + mNew(i, 2) * mOld(j, 2) &
                                       + mNew(i, 3) * mOld(j, 3)
            end do
        end do

        ! Determine the rotation point at the old time level; it is
        ! possible that this value changes due to translation of the grid.

        !  aInf = sqrt(gammaInf*pInf/rhoInf)

        !  rotationPoint(1) = LRef*rotPoint(1) &
        !                   + MachGrid(1)*aInf*tOld/timeRef
        !  rotationPoint(2) = LRef*rotPoint(2) &
        !                   + MachGrid(2)*aInf*tOld/timeRef
        !  rotationPoint(3) = LRef*rotPoint(3) &
        !                   + MachGrid(3)*aInf*tOld/timeRef

        rotationPoint(1) = LRef * rotPoint(1)
        rotationPoint(2) = LRef * rotPoint(2)
        rotationPoint(3) = LRef * rotPoint(3)

    end subroutine rotMatrixRigidBody

    function secondDerivativeRigidRotAngle(degreePolRot, &
                                           coefPolRot, &
                                           degreeFourRot, &
                                           omegaFourRot, &
                                           cosCoefFourRot, &
                                           sinCoefFourRot, t)
        !
        !       2ndderivativeRigidRotAngle computes the 2nd time derivative of
        !       the rigid body rotation angle at the given time for the given
        !       arguments. The angle is described by a combination of a
        !       polynomial and fourier series.
        !
        use constants
        use flowVarRefState, only: timeRef
        use inputPhysics, only: equationMode
        implicit none
        !
        !      Function type
        !
        real(kind=realType) :: secondDerivativeRigidRotAngle
        !
        !      Function arguments.
        !
        integer(kind=intType), intent(in) :: degreePolRot
        integer(kind=intType), intent(in) :: degreeFourRot

        real(kind=realType), intent(in) :: omegaFourRot, t

        real(kind=realType), dimension(0:*), intent(in) :: coefPolRot
        real(kind=realType), dimension(0:*), intent(in) :: cosCoefFourRot
        real(kind=realType), dimension(*), intent(in) :: sinCoefFourRot
        !
        !      Local variables.
        !
        integer(kind=intType) :: nn

        real(kind=realType) :: dPhi, val

        ! Return immediately if this is a steady computation.

        if (equationMode == steady) then
            secondDerivativeRigidRotAngle = zero
            return
        end if

        ! Compute the polynomial contribution.

        dPhi = zero
        do nn = 2, degreePolRot
            dPhi = dPhi + (nn - 1) * nn * coefPolRot(nn) * (t**(nn - 2))
        end do

        ! Compute the fourier contribution.

        do nn = 1, degreeFourRot
            val = nn * omegaFourRot
            dPhi = dPhi - val**2 * sinCoefFourRot(nn) * sin(val * t)
            dPhi = dPhi - val**2 * cosCoefFourRot(nn) * cos(val * t)
        end do

        ! Set derivativeRigidRotAngle to dPhi. Multiply by timeRef
        ! to obtain the correct non-dimensional value.

        secondDerivativeRigidRotAngle = timeRef**2 * dPhi

    end function secondDerivativeRigidRotAngle

    function rigidRotAngle(degreePolRot, coefPolRot, &
                           degreeFourRot, omegaFourRot, &
                           cosCoefFourRot, sinCoefFourRot, t)
        !
        !       rigidRotAngle computes the rigid body rotation angle at the
        !       given time for the given arguments. The angle is described by
        !       a combination of a polynomial and fourier series.
        !
        use constants
        use inputPhysics, only: equationMode
        implicit none
        !
        !      Function type
        !
        real(kind=realType) :: rigidRotAngle
        !
        !      Function arguments.
        !
        integer(kind=intType), intent(in) :: degreePolRot
        integer(kind=intType), intent(in) :: degreeFourRot

        real(kind=realType), intent(in) :: omegaFourRot, t

        real(kind=realType), dimension(0:*), intent(in) :: coefPolRot
        real(kind=realType), dimension(0:*), intent(in) :: cosCoefFourRot
        real(kind=realType), dimension(*), intent(in) :: sinCoefFourRot
        !
        !      Local variables.
        !
        integer(kind=intType) :: nn

        real(kind=realType) :: phi, val

        ! Return immediately if this is a steady computation.

        if (equationMode == steady) then
            rigidRotAngle = zero
            return
        end if

        ! Compute the polynomial contribution. If no polynomial was
        ! specified, the value of index 0 is set to zero automatically.

        phi = coefPolRot(0)
        do nn = 1, degreePolRot
            phi = phi + coefPolRot(nn) * (t**nn)
        end do

        ! Compute the fourier contribution. Again the cosine coefficient
        ! of index 0 is defaulted to zero if not specified.

        phi = phi + cosCoefFourRot(0)
        do nn = 1, degreeFourRot
            val = nn * omegaFourRot * t
            phi = phi + cosCoefFourRot(nn) * cos(val) &
                  + sinCoefFourRot(nn) * sin(val)
        end do

        ! Set rigidRotAngle to phi.

        rigidRotAngle = phi

    end function rigidRotAngle

    subroutine setBCPointers(nn, spatialPointers)
        !
        !       setBCPointers sets the pointers needed for the boundary
        !       condition treatment on a general face, such that the boundary
        !       routines are only implemented once instead of 6 times.
        !
        use constants
        use blockPointers, only: w, p, rlv, rev, gamma, x, d2wall, &
                                 si, sj, sk, s, globalCell, BCData, nx, il, ie, ib, &
                                 ny, jl, je, jb, nz, kl, ke, kb, BCFaceID, &
                                 addgridvelocities, sFaceI, sFaceJ, sFaceK, addGridVelocities
        use BCPointers, only: ww0, ww1, ww2, ww3, pp0, pp1, pp2, pp3, &
                              rlv0, rlv1, rlv2, rlv3, rev0, rev1, rev2, rev3, &
                              gamma0, gamma1, gamma2, gamma3, gcp, xx, ss, ssi, ssj, ssk, dd2wall, &
                              sFace, iStart, iEnd, jStart, jEnd, iSize, jSize
        use inputPhysics, only: cpModel, equations
        implicit none

        ! Subroutine arguments.
        integer(kind=intType), intent(in) :: nn
        logical, intent(in) :: spatialPointers

        ! Determine the sizes of each face and point to just the range we
        ! need on each face.
        iStart = BCData(nn)%icBeg
        iEnd = BCData(nn)%icEnd
        jStart = BCData(nn)%jcBeg
        jEnd = BCData(nn)%jcEnd

        ! Set the size of the subface
        isize = iEnd - iStart + 1
        jsize = jEnd - jStart + 1

        ! Determine the face id on which the subface is located and set
        ! the pointers accordinly.

        select case (BCFaceID(nn))

            !---------------------------------------------------------------------------
        case (iMin)

            ww3 => w(3, 1:, 1:, :)
            ww2 => w(2, 1:, 1:, :)
            ww1 => w(1, 1:, 1:, :)
            ww0 => w(0, 1:, 1:, :)

            pp3 => p(3, 1:, 1:)
            pp2 => p(2, 1:, 1:)
            pp1 => p(1, 1:, 1:)
            pp0 => p(0, 1:, 1:)

            rlv3 => rlv(3, 1:, 1:)
            rlv2 => rlv(2, 1:, 1:)
            rlv1 => rlv(1, 1:, 1:)
            rlv0 => rlv(0, 1:, 1:)

            rev3 => rev(3, 1:, 1:)
            rev2 => rev(2, 1:, 1:)
            rev1 => rev(1, 1:, 1:)
            rev0 => rev(0, 1:, 1:)

            gamma3 => gamma(3, 1:, 1:)
            gamma2 => gamma(2, 1:, 1:)
            gamma1 => gamma(1, 1:, 1:)
            gamma0 => gamma(0, 1:, 1:)

            gcp => globalCell(2, 1:, 1:)
            !---------------------------------------------------------------------------

        case (iMax)

            ww3 => w(nx, 1:, 1:, :)
            ww2 => w(il, 1:, 1:, :)
            ww1 => w(ie, 1:, 1:, :)
            ww0 => w(ib, 1:, 1:, :)

            pp3 => p(nx, 1:, 1:)
            pp2 => p(il, 1:, 1:)
            pp1 => p(ie, 1:, 1:)
            pp0 => p(ib, 1:, 1:)

            rlv3 => rlv(nx, 1:, 1:)
            rlv2 => rlv(il, 1:, 1:)
            rlv1 => rlv(ie, 1:, 1:)
            rlv0 => rlv(ib, 1:, 1:)

            rev3 => rev(nx, 1:, 1:)
            rev2 => rev(il, 1:, 1:)
            rev1 => rev(ie, 1:, 1:)
            rev0 => rev(ib, 1:, 1:)

            gamma3 => gamma(nx, 1:, 1:)
            gamma2 => gamma(il, 1:, 1:)
            gamma1 => gamma(ie, 1:, 1:)
            gamma0 => gamma(ib, 1:, 1:)

            gcp => globalCell(il, 1:, 1:)
            !---------------------------------------------------------------------------

        case (jMin)

            ww3 => w(1:, 3, 1:, :)
            ww2 => w(1:, 2, 1:, :)
            ww1 => w(1:, 1, 1:, :)
            ww0 => w(1:, 0, 1:, :)

            pp3 => p(1:, 3, 1:)
            pp2 => p(1:, 2, 1:)
            pp1 => p(1:, 1, 1:)
            pp0 => p(1:, 0, 1:)

            rlv3 => rlv(1:, 3, 1:)
            rlv2 => rlv(1:, 2, 1:)
            rlv1 => rlv(1:, 1, 1:)
            rlv0 => rlv(1:, 0, 1:)

            rev3 => rev(1:, 3, 1:)
            rev2 => rev(1:, 2, 1:)
            rev1 => rev(1:, 1, 1:)
            rev0 => rev(1:, 0, 1:)

            gamma3 => gamma(1:, 3, 1:)
            gamma2 => gamma(1:, 2, 1:)
            gamma1 => gamma(1:, 1, 1:)
            gamma0 => gamma(1:, 0, 1:)

            gcp => globalCell(1:, 2, 1:)
            !---------------------------------------------------------------------------

        case (jMax)

            ww3 => w(1:, ny, 1:, :)
            ww2 => w(1:, jl, 1:, :)
            ww1 => w(1:, je, 1:, :)
            ww0 => w(1:, jb, 1:, :)

            pp3 => p(1:, ny, 1:)
            pp2 => p(1:, jl, 1:)
            pp1 => p(1:, je, 1:)
            pp0 => p(1:, jb, 1:)

            rlv3 => rlv(1:, ny, 1:)
            rlv2 => rlv(1:, jl, 1:)
            rlv1 => rlv(1:, je, 1:)
            rlv0 => rlv(1:, jb, 1:)

            rev3 => rev(1:, ny, 1:)
            rev2 => rev(1:, jl, 1:)
            rev1 => rev(1:, je, 1:)
            rev0 => rev(1:, jb, 1:)

            gamma3 => gamma(1:, ny, 1:)
            gamma2 => gamma(1:, jl, 1:)
            gamma1 => gamma(1:, je, 1:)
            gamma0 => gamma(1:, jb, 1:)

            gcp => globalCell(1:, jl, 1:)
            !---------------------------------------------------------------------------

        case (kMin)

            ww3 => w(1:, 1:, 3, :)
            ww2 => w(1:, 1:, 2, :)
            ww1 => w(1:, 1:, 1, :)
            ww0 => w(1:, 1:, 0, :)

            pp3 => p(1:, 1:, 3)
            pp2 => p(1:, 1:, 2)
            pp1 => p(1:, 1:, 1)
            pp0 => p(1:, 1:, 0)

            rlv3 => rlv(1:, 1:, 3)
            rlv2 => rlv(1:, 1:, 2)
            rlv1 => rlv(1:, 1:, 1)
            rlv0 => rlv(1:, 1:, 0)

            rev3 => rev(1:, 1:, 3)
            rev2 => rev(1:, 1:, 2)
            rev1 => rev(1:, 1:, 1)
            rev0 => rev(1:, 1:, 0)

            gamma3 => gamma(1:, 1:, 3)
            gamma2 => gamma(1:, 1:, 2)
            gamma1 => gamma(1:, 1:, 1)
            gamma0 => gamma(1:, 1:, 0)

            gcp => globalCell(1:, 1:, 2)
            !---------------------------------------------------------------------------

        case (kMax)

            ww3 => w(1:, 1:, nz, :)
            ww2 => w(1:, 1:, kl, :)
            ww1 => w(1:, 1:, ke, :)
            ww0 => w(1:, 1:, kb, :)

            pp3 => p(1:, 1:, nz)
            pp2 => p(1:, 1:, kl)
            pp1 => p(1:, 1:, ke)
            pp0 => p(1:, 1:, kb)

            rlv3 => rlv(1:, 1:, nz)
            rlv2 => rlv(1:, 1:, kl)
            rlv1 => rlv(1:, 1:, ke)
            rlv0 => rlv(1:, 1:, kb)

            rev3 => rev(1:, 1:, nz)
            rev2 => rev(1:, 1:, kl)
            rev1 => rev(1:, 1:, ke)
            rev0 => rev(1:, 1:, kb)

            gamma3 => gamma(1:, 1:, nz)
            gamma2 => gamma(1:, 1:, kl)
            gamma1 => gamma(1:, 1:, ke)
            gamma0 => gamma(1:, 1:, kb)

            gcp => globalCell(1:, 1:, kl)
        end select

        if (spatialPointers) then
            select case (BCFaceID(nn))
            case (iMin)
                xx => x(1, :, :, :)
                ssi => si(1, :, :, :)
                ssj => sj(2, :, :, :)
                ssk => sk(2, :, :, :)
                ss => s(2, :, :, :)
            case (iMax)
                xx => x(il, :, :, :)
                ssi => si(il, :, :, :)
                ssj => sj(il, :, :, :)
                ssk => sk(il, :, :, :)
                ss => s(il, :, :, :)
            case (jMin)
                xx => x(:, 1, :, :)
                ssi => sj(:, 1, :, :)
                ssj => si(:, 2, :, :)
                ssk => sk(:, 2, :, :)
                ss => s(:, 2, :, :)
            case (jMax)
                xx => x(:, jl, :, :)
                ssi => sj(:, jl, :, :)
                ssj => si(:, jl, :, :)
                ssk => sk(:, jl, :, :)
                ss => s(:, jl, :, :)
            case (kMin)
                xx => x(:, :, 1, :)
                ssi => sk(:, :, 1, :)
                ssj => si(:, :, 2, :)
                ssk => sj(:, :, 2, :)
                ss => s(:, :, 2, :)
            case (kMax)
                xx => x(:, :, kl, :)
                ssi => sk(:, :, kl, :)
                ssj => si(:, :, kl, :)
                ssk => sj(:, :, kl, :)
                ss => s(:, :, kl, :)
            end select

            if (addGridVelocities) then
                select case (BCFaceID(nn))
                case (iMin)
                    sFace => sFaceI(1, :, :)
                case (iMax)
                    sFace => sFaceI(il, :, :)
                case (jMin)
                    sFace => sFaceJ(:, 1, :)
                case (jMax)
                    sFace => sFaceJ(:, jl, :)
                case (kMin)
                    sFace => sFaceK(:, :, 1)
                case (kMax)
                    sFace => sFaceK(:, :, kl)
                end select
            end if

            if (equations == RANSEquations) then
                select case (BCFaceID(nn))
                case (iMin)
                    dd2Wall => d2Wall(2, :, :)
                case (iMax)
                    dd2Wall => d2Wall(il, :, :)
                case (jMin)
                    dd2Wall => d2Wall(:, 2, :)
                case (jMax)
                    dd2Wall => d2Wall(:, jl, :)
                case (kMin)
                    dd2Wall => d2Wall(:, :, 2)
                case (kMax)
                    dd2Wall => d2Wall(:, :, kl)
                end select
            end if
        end if
    end subroutine setBCPointers

    subroutine computeRootBendingMoment(cf, cm, bendingMoment)

        !                                                      *
        ! Compute a normalized bending moment coefficient from *
        ! the force and moment coefficient. At the moment this *
        ! Routine only works for a half body. Additional logic *
        ! would be needed for a full body.                     *
        !                                                      *

        use constants
        use inputPhysics, only: lengthRef, pointRef, pointRefEC, liftIndex
        implicit none

        !input/output variables
        real(kind=realType), intent(in), dimension(3) :: cf, cm
        real(kind=realType), intent(out) :: bendingMoment

        !Subroutine Variables
        real(kind=realType) :: elasticMomentx, elasticMomenty, elasticMomentz
        bendingMoment = zero
        if (liftIndex == 2) then
            !z out wing sum momentx,momentz
            elasticMomentx = cm(1) + cf(2) * (pointRefEC(3) - &
                                              pointRef(3)) / lengthref - cf(3) * &
                             (pointRefEC(2) - pointRef(2)) / lengthref
            elasticMomentz = cm(3) - cf(2) * (pointRefEC(1) - &
                                              pointref(1)) / lengthref + cf(1) * &
                             (pointRefEC(2) - pointRef(2)) / lengthref
            bendingMoment = sqrt(elasticMomentx**2 + elasticMomentz**2)
        elseif (liftIndex == 3) then
            !y out wing sum momentx,momenty
            elasticMomentx = cm(1) + cf(3) * (pointrefEC(2) - &
                                              pointRef(2)) / lengthref + &
                             cf(3) * (pointrefEC(3) - pointref(3)) / lengthref
            elasticMomenty = cm(2) + cf(3) * (pointRefEC(1) - &
                                              pointRef(1)) / lengthref + &
                             cf(1) * (pointrefEC(3) - pointRef(3)) / lengthref
            bendingMoment = sqrt(elasticMomentx**2 + elasticMomenty**2)
        end if

    end subroutine computeRootBendingMoment

    subroutine computeLeastSquaresRegression(y, x, npts, m, b)
        !
        !       Computes the slope of best fit for a set of x,y data of length
        !       npts
        !
        use constants
        implicit none
        !Subroutine arguments
        integer(kind=intType) :: npts
        real(kind=realType), dimension(npts) :: x, y
        real(kind=realType) :: m, b

        !local variables
        real(kind=realType) :: sumx, sumy, sumx2, sumxy
        integer(kind=intType) :: i

        !begin execution
        sumx = 0.0
        sumy = 0.0
        sumx2 = 0.0
        sumxy = 0.0
        do i = 1, npts

            sumx = sumx + x(i)
            sumy = sumy + y(i)
            sumx2 = sumx2 + x(i) * x(i)
            sumxy = sumxy + x(i) * y(i)
        end do

        m = ((npts * sumxy) - (sumy * sumx)) / ((npts * sumx2) - (sumx)**2)
        b = (sumy * sumx2 - (sumx * sumxy)) / ((npts * sumx2) - (sumx)**2)

    end subroutine computeLeastSquaresRegression

    subroutine computeTSDerivatives(force, moment, coef0, dcdalpha, &
                                    dcdalphadot, dcdq, dcdqdot)
        !
        !      Computes the stability derivatives based on the time spectral
        !      solution of a given mesh. Takes in the force coefficients at
        !      all time instantces and computes the agregate parameters
        !
        use constants
        use communication
        use inputPhysics
        use inputTimeSpectral
        use inputTSStabDeriv
        use flowvarrefstate
        use monitor
        use section
        use inputMotion
        implicit none

        !
        !     Subroutine arguments.
        !
        real(kind=realType), dimension(3, nTimeIntervalsSpectral) :: force, moment
        real(kind=realType), dimension(8) :: dcdq, dcdqdot
        real(kind=realType), dimension(8) :: dcdalpha, dcdalphadot
        real(kind=realType), dimension(8) :: Coef0

        ! Working Variables
        real(kind=realType), dimension(nTimeIntervalsSpectral, 8) :: baseCoef
        real(kind=realType), dimension(8) :: coef0dot
        real(kind=realType), dimension(nTimeIntervalsSpectral, 8) :: ResBaseCoef
        real(kind=realType), dimension(nTimeIntervalsSpectral) :: intervalAlpha, intervalAlphadot
        real(kind=realType), dimension(nTimeIntervalsSpectral) :: intervalMach, intervalMachdot
        real(kind=realType), dimension(nSections) :: t
        integer(kind=intType) :: i, sps, nn
        !speed of sound: for normalization of q derivatives
        real(kind=realType) :: a
        real(kind=realType) :: fact, factMoment
        ! Functions
        real(kind=realType), dimension(nTimeIntervalsSpectral) :: dPhix, dPhiy, dphiz
        real(kind=realType), dimension(nTimeIntervalsSpectral) :: dPhixdot, dPhiydot, dphizdot
        real(kind=realType) :: derivativeRigidRotAngle, secondDerivativeRigidRotAngle

        fact = two / (gammaInf * pInf * MachCoef**2 &
                      * surfaceRef * LRef**2)
        factMoment = fact / (lengthRef * LRef)

        if (TSqMode) then

            print *, 'TS Q Mode code needs to be updated in computeTSDerivatives!'
            stop

            ! !q is pitch
            ! do sps =1,nTimeIntervalsSpectral
            !    !compute the time of this intervavc
            !    t = timeUnsteadyRestart

            !    if(equationMode == timeSpectral) then
            !       do nn=1,nSections
            !          t(nn) = t(nn) + (sps-1)*sections(nn)%timePeriod &
            !               /         (nTimeIntervalsSpectral*1.0)
            !       enddo
            !    endif

            !    ! Compute the time derivative of the rotation angles around the
            !    ! z-axis. i.e. compute q

            !    dphiZ(sps) = derivativeRigidRotAngle(degreePolZRot,   &
            !         coefPolZRot,     &
            !         degreeFourZRot,  &
            !         omegaFourZRot,   &
            !         cosCoefFourZRot, &
            !         sinCoefFourZRot, t)

            !    ! add in q_dot computation
            !    dphiZdot(sps) = secondDerivativeRigidRotAngle(degreePolZRot,   &
            !         coefPolZRot,     &
            !         degreeFourZRot,  &
            !         omegaFourZRot,   &
            !         cosCoefFourZRot, &
            !         sinCoefFourZRot, t)
            ! end do

            ! !now compute dCl/dq
            ! do i =1,8
            !    call computeLeastSquaresRegression(BaseCoef(:,i),dphiz,nTimeIntervalsSpectral,dcdq(i),coef0(i))
            ! end do

            ! ! now subtract off estimated cl,cmz and use remainder to compute
            ! ! clqdot and cmzqdot.
            ! do i = 1,8
            !    do sps = 1,nTimeIntervalsSpectral
            !       ResBaseCoef(sps,i) = BaseCoef(sps,i)-(dcdq(i)*dphiz(sps)+Coef0(i))
            !    enddo
            ! enddo

            ! !now normalize the results...
            ! a  = sqrt(gammaInf*pInfDim/rhoInfDim)
            ! dcdq = dcdq*timeRef*2*(machGrid*a)/lengthRef

            ! !now compute dCl/dpdot
            ! do i = 1,8
            !    call computeLeastSquaresRegression(ResBaseCoef(:,i),dphizdot,nTimeIntervalsSpectral,dcdqdot(i),Coef0dot(i))
            ! enddo

        elseif (TSAlphaMode) then

            do sps = 1, nTimeIntervalsSpectral

                !compute the time of this interval
                t = timeUnsteadyRestart

                if (equationMode == timeSpectral) then
                    do nn = 1, nSections
                        t(nn) = t(nn) + (sps - 1) * sections(nn)%timePeriod &
                                / (nTimeIntervalsSpectral * 1.0)
                    end do
                end if

                intervalAlpha(sps) = TSAlpha(degreePolAlpha, coefPolAlpha, &
                                             degreeFourAlpha, omegaFourAlpha, &
                                             cosCoefFourAlpha, sinCoefFourAlpha, t(1))

                intervalAlphadot(sps) = TSAlphadot(degreePolAlpha, coefPolAlpha, &
                                                   degreeFourAlpha, omegaFourAlpha, &
                                                   cosCoefFourAlpha, sinCoefFourAlpha, t(1))

                ! THIS CALL IS WRONG!!!!
                !call getDirAngle(velDirFreestream,liftDirection,liftIndex,alpha+intervalAlpha(sps), beta)

                BaseCoef(sps, 1) = fact * ( &
                                   force(1, sps) * liftDirection(1) + &
                                   force(2, sps) * liftDirection(2) + &
                                   force(3, sps) * liftDIrection(3))
                BaseCoef(sps, 2) = fact * ( &
                                   force(1, sps) * dragDirection(1) + &
                                   force(2, sps) * dragDirection(2) + &
                                   force(3, sps) * dragDIrection(3))
                BaseCoef(sps, 3) = force(1, sps) * fact
                BaseCoef(sps, 4) = force(2, sps) * fact
                BaseCoef(sps, 5) = force(3, sps) * fact
                BaseCoef(sps, 6) = moment(1, sps) * factMoment
                BaseCoef(sps, 7) = moment(2, sps) * factMoment
                BaseCoef(sps, 8) = moment(3, sps) * factMoment
            end do

            !now compute dCl/dalpha
            do i = 1, 8
                call computeLeastSquaresRegression(BaseCoef(:, i), &
                                                   intervalAlpha, nTimeIntervalsSpectral, dcdAlpha(i), coef0(i))
            end do

            ! now subtract off estimated cl,cmz and use remainder to compute
            ! clalphadot and cmzalphadot.
            do i = 1, 8
                do sps = 1, nTimeIntervalsSpectral
                    ResBaseCoef(sps, i) = BaseCoef(sps, i) - (dcdalpha(i) * intervalAlpha(sps) + Coef0(i))
                end do
            end do

            !now compute dCi/dalphadot
            do i = 1, 8
                call computeLeastSquaresRegression(ResBaseCoef(:, i), &
                                                   intervalAlphadot, nTimeIntervalsSpectral, &
                                                   dcdalphadot(i), Coef0dot(i))
            end do

            a = sqrt(gammaInf * pInfDim / rhoInfDim)
            dcdalphadot = dcdalphadot * 2 * (machGrid * a) / lengthRef

        else
            call terminate('computeTSDerivatives', 'Not a valid stability motion')
        end if

    end subroutine computeTSDerivatives

    subroutine getDirAngle(freeStreamAxis, liftAxis, liftIndex, alpha, beta)
        !
        !      Convert the wind axes to angle of attack and side slip angle.
        !      The direction angles alpha and beta are computed given the
        !      components of the wind direction vector (freeStreamAxis), the
        !      lift direction vector (liftAxis) and assuming that the
        !      body direction (xb,yb,zb) is in the default ijk coordinate
        !      system. The rotations are determined by first determining
        !      whether the lift is primarily in the j or k direction and then
        !      determining the angles accordingly.
        !      direction vector:
        !        1) Rotation about the zb or yb -axis: alpha clockwise (CW)
        !           (xb,yb,zb) -> (x1,y1,z1)
        !        2) Rotation about the yl or z1 -axis: beta counter-clockwise
        !           (CCW) (x1,y1,z1) -> (xw,yw,zw)
        !         input arguments:
        !            freeStreamAxis = wind vector in body axes
        !            liftAxis       = lift direction vector in body axis
        !         output arguments:
        !            alpha    = angle of attack in radians
        !            beta     = side slip angle in radians
        !
        use constants

        implicit none
        !
        !     Subroutine arguments.
        !
        !      real(kind=realType), intent(in)  :: xw, yw, zw
        real(kind=realType), dimension(3), intent(in) :: freeStreamAxis
        real(kind=realType), dimension(3), intent(in) :: liftAxis
        real(kind=realType), intent(out) :: alpha, beta
        integer(kind=intType), intent(out) :: liftIndex
        !
        !     Local variables.
        !
        real(kind=realType) :: rnorm
        integer(kind=intType) :: flowIndex, i
        real(kind=realType), dimension(3) :: freeStreamAxisNorm
        integer(kind=intType) :: temp

        ! Assume domoniate flow is x

        flowIndex = 1

        ! Determine the dominant lift direction
        if (abs(liftAxis(1)) > abs(liftAxis(2)) .and. &
            abs(liftAxis(1)) > abs(liftAxis(3))) then
            temp = 1
        else if (abs(liftAxis(2)) > abs(liftAxis(1)) .and. &
                 abs(liftAxis(2)) > abs(liftAxis(3))) then
            temp = 2
        else
            temp = 3
        end if

        liftIndex = temp

        ! Normalize the freeStreamDirection vector.
        rnorm = sqrt(freeStreamAxis(1)**2 + freeStreamAxis(2)**2 + freeStreamAxis(3)**2)
        do i = 1, 3
            freeStreamAxisNorm(i) = freeStreamAxis(i) / rnorm
        end do

        if (liftIndex == 2) then
            ! different coordinate system for aerosurf
            ! Wing is in z- direction
            ! Compute angle of attack alpha.

            alpha = asin(freeStreamAxisNorm(2))

            ! Compute side-slip angle beta.

            beta = -atan2(freeStreamAxisNorm(3), freeStreamAxisNorm(1))

        elseif (liftIndex == 3) then
            ! Wing is in y- direction

            ! Compute angle of attack alpha.

            alpha = asin(freeStreamAxisNorm(3))

            ! Compute side-slip angle beta.

            beta = atan2(freeStreamAxisNorm(2), freeStreamAxisNorm(1))
        else
            call terminate('getDirAngle', 'Invalid Lift Direction')
        end if
    end subroutine getDirAngle

    subroutine stabilityDerivativeDriver
        !
        !      Runs the Time spectral stability derivative routines from the
        !      main program file
        !
        use precision
        implicit none
        !
        !     Local variables.
        !
        real(kind=realType), dimension(8) :: dcdalpha, dcdalphadot, dcdbeta, &
                                             dcdbetadot, dcdMach, dcdMachdot
        real(kind=realType), dimension(8) :: dcdp, dcdpdot, dcdq, dcdqdot, dcdr, dcdrdot
        real(kind=realType), dimension(8) :: Coef0, Coef0dot

        !call computeTSDerivatives(coef0,dcdalpha,dcdalphadot,dcdq,dcdqdot)

    end subroutine stabilityDerivativeDriver
    subroutine setCoefTimeIntegrator
        !
        !       setCoefTimeIntegrator determines the coefficients of the
        !       time integration scheme in unsteady mode. Normally these are
        !       equal to the coefficients corresponding to the specified
        !       accuracy. However during the initial phase there are not
        !       enough states in the past and the accuracy is reduced.
        !
        use constants
        use inputUnsteady
        use inputPhysics
        use iteration
        use monitor
        implicit none
        !
        !      Local variables.
        !
        integer(kind=intType) :: nn, nLevelsSet

        ! Determine which time integrator must be used.

        ! Modified by HDN
        select case (timeAccuracy)
        case (firstOrder)

            ! 1st order. No need to check the number of available
            ! states in the past. Set the two coefficients and
            ! nLevelsSet to 2.

            coefTime(0) = 1.0_realType
            coefTime(1) = -1.0_realType

            if (useALE .and. equationMode .eq. unsteady) then
                coefTimeALE(1) = 1.0_realType
                coefMeshALE(1, 1) = half
                coefMeshALE(1, 2) = half
            end if

            nLevelsSet = 2

            !--------------------------------------------------

        case (secondOrder)

            ! Second order time integrator. Determine the amount of
            ! available states and set the coefficients accordingly.
            select case (nOldSolAvail)

            case (1_intType)
                coefTime(0) = 1.0_realType
                coefTime(1) = -1.0_realType

                if (useALE .and. equationMode .eq. unsteady) then
                    coefTimeALE(1) = half
                    coefTimeALE(2) = half
                    coefTimeALE(3) = zero
                    coefTimeALE(4) = zero

                    coefMeshALE(1, 1) = half
                    coefMeshALE(1, 2) = half
                    coefMeshALE(2, 1) = half
                    coefMeshALE(2, 2) = half
                end if

                nLevelsSet = 2

            case default   ! 2 or bigger.
                coefTime(0) = 1.5_realType
                coefTime(1) = -2.0_realType
                coefTime(2) = 0.5_realType

                if (useALE .and. equationMode .eq. unsteady) then
                    coefTimeALE(1) = threefourth
                    coefTimeALE(2) = threefourth
                    coefTimeALE(3) = -fourth
                    coefTimeALE(4) = -fourth

                    coefMeshALE(1, 1) = half * (1.0_realType + 1.0_realType / sqrtthree)
                    coefMeshALE(1, 2) = half * (1.0_realType - 1.0_realType / sqrtthree)
                    coefMeshALE(2, 1) = coefMeshALE(1, 2)
                    coefMeshALE(2, 2) = coefMeshALE(1, 1)
                end if

                nLevelsSet = 3

            end select

            !--------------------------------------------------

        case (thirdOrder)

            ! Third order time integrator.  Determine the amount of
            ! available states and set the coefficients accordingly.

            select case (nOldSolAvail)

            case (1_intType)
                coefTime(0) = 1.0_realType
                coefTime(1) = -1.0_realType

                if (useALE .and. equationMode .eq. unsteady) then
                    coefTimeALE(1) = 1.0_realType
                    coefMeshALE(1, 1) = half
                    coefMeshALE(1, 2) = half
                end if

                nLevelsSet = 2

            case (2_intType)
                coefTime(0) = 1.5_realType
                coefTime(1) = -2.0_realType
                coefTime(2) = 0.5_realType

                if (useALE .and. equationMode .eq. unsteady) then
                    coefTimeALE(1) = threefourth
                    coefTimeALE(2) = -fourth
                    coefMeshALE(1, 1) = half * (1.0_realType + 1.0_realType / sqrtthree)
                    coefMeshALE(1, 2) = half * (1.0_realType - 1.0_realType / sqrtthree)
                    coefMeshALE(2, 1) = coefMeshALE(1, 2)
                    coefMeshALE(2, 2) = coefMeshALE(1, 1)
                end if

                nLevelsSet = 3

            case default   ! 3 or bigger.
                coefTime(0) = 11.0_realType / 6.0_realType
                coefTime(1) = -3.0_realType
                coefTime(2) = 1.5_realType
                coefTime(3) = -1.0_realType / 3.0_realType

                ! These numbers are NOT correct
                ! DO NOT use 3rd order ALE for now
                if (useALE .and. equationMode .eq. unsteady) then
                    print *, 'Third-order ALE not implemented yet.'
                    coefTimeALE(1) = threefourth
                    coefTimeALE(2) = threefourth
                    coefTimeALE(3) = -fourth
                    coefTimeALE(4) = -fourth
                    coefMeshALE(1, 1) = half * (1.0_realType + 1.0_realType / sqrtthree)
                    coefMeshALE(1, 2) = half * (1.0_realType - 1.0_realType / sqrtthree)
                    coefMeshALE(2, 1) = coefMeshALE(1, 2)
                    coefMeshALE(2, 2) = coefMeshALE(1, 1)
                    coefMeshALE(3, 1) = coefMeshALE(1, 2)
                    coefMeshALE(3, 2) = coefMeshALE(1, 1)
                end if

                nLevelsSet = 4

            end select

        end select

        ! Set the rest of the coefficients to 0 if not enough states
        ! in the past are available.

        do nn = nLevelsSet, nOldLevels
            coefTime(nn) = zero
        end do

    end subroutine setCoefTimeIntegrator

    function myNorm2(x)
        use constants
        implicit none
        real(kind=realType), dimension(3), intent(in) :: x
        real(kind=realType) :: myNorm2
        myNorm2 = sqrt(x(1)**2 + x(2)**2 + x(3)**2)
    end function myNorm2

    function isWallType(bType)

        use constants
        implicit none
        integer(kind=intType) :: bType
        logical :: isWallType

        isWallType = .False.
        if (bType == NSWallAdiabatic .or. &
            bType == NSWallIsoThermal .or. &
            bType == EulerWall) then
            isWallType = .True.
        end if

    end function isWallType

    subroutine cross_prod(a, b, c)

        use precision

        ! Inputs
        real(kind=realType), dimension(3), intent(in) :: a, b

        ! Outputs
        real(kind=realType), dimension(3), intent(out) :: c

        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)

    end subroutine cross_prod

    subroutine siAngle(angle, mult, trans)

        use constants
        use su_cgns, only: Radian, Degree
        implicit none
        !
        !      Subroutine arguments.
        !
        integer, intent(in) :: angle
        real(kind=realType), intent(out) :: mult, trans

        ! Determine the situation we are having here.

        if (angle == Radian) then

            ! Angle is already given in radIans. No need for a conversion.

            mult = one
            trans = zero

        else if (angle == Degree) then

            ! Angle is given in degrees. A multiplication must be performed.

            mult = pi / 180.0_realType
            trans = zero

        else

            call terminate("siAngle", &
                           "No idea how to convert this to SI units")

        end if

    end subroutine siAngle

    subroutine siDensity(mass, len, mult, trans)
        !
        !       siDensity computes the conversion from the given density
        !       unit, which can be constructed from mass and length, to the
        !       SI-unit kg/m^3. The conversion will look like:
        !       density in kg/m^3 = mult*(density in NCU) + trans.
        !       NCU means non-christian units, i.e. everything that is not SI.
        !
        use constants
        use su_cgns, only: Kilogram, meter
        implicit none
        !
        !      Subroutine arguments.
        !
        integer, intent(in) :: mass, len
        real(kind=realType), intent(out) :: mult, trans

        ! Determine the situation we are having here.

        if (mass == Kilogram .and. len == Meter) then

            ! Density is given in kg/m^3, i.e. no need for a conversion.

            mult = one
            trans = zero

        else

            call terminate("siDensity", &
                           "No idea how to convert this to SI units")

        end if

    end subroutine siDensity

    subroutine siLen(len, mult, trans)
        !
        !       siLen computes the conversion from the given length unit to
        !       the SI-unit meter. The conversion will look like:
        !       length in meter = mult*(length in NCU) + trans.
        !       NCU means non-christian units, i.e. everything that is not SI.
        !
        use constants
        use su_cgns, only: Meter, Centimeter, millimeter, Foot, Inch
        implicit none
        !
        !      Subroutine arguments.
        !
        integer, intent(in) :: len
        real(kind=realType), intent(out) :: mult, trans

        ! Determine the situation we are having here.

        select case (len)

        case (Meter)
            mult = one; trans = zero

        case (CenTimeter)
            mult = 0.01_realType; trans = zero

        case (Millimeter)
            mult = 0.001_realType; trans = zero

        case (Foot)
            mult = 0.3048_realType; trans = zero

        case (Inch)
            mult = 0.0254_realType; trans = zero

        case default
            call terminate("siLen", &
                           "No idea how to convert this to SI units")

        end select

    end subroutine siLen

    subroutine siPressure(mass, len, time, mult, trans)
        !
        !       siPressure computes the conversion from the given pressure
        !       unit, which can be constructed from mass, length and time, to
        !       the SI-unit Pa. The conversion will look like:
        !       pressure in Pa = mult*(pressure in NCU) + trans.
        !       NCU means non-christian units, i.e. everything that is not SI.
        !
        use constants
        use su_cgns, only: Kilogram, Meter, Second
        implicit none
        !
        !      Subroutine arguments.
        !
        integer, intent(in) :: mass, len, time
        real(kind=realType), intent(out) :: mult, trans

        ! Determine the situation we are having here.

        if (mass == Kilogram .and. len == Meter .and. time == Second) then

            ! Pressure is given in Pa, i.e. no need for a conversion.

            mult = one
            trans = zero

        else

            call terminate("siPressure", &
                           "No idea how to convert this to SI units")

        end if

    end subroutine siPressure

    subroutine siTemperature(temp, mult, trans)
        !
        !       siTemperature computes the conversion from the given
        !       temperature unit to the SI-unit kelvin. The conversion will
        !       look like:
        !       temperature in K = mult*(temperature in NCU) + trans.
        !       NCU means non-christian units, i.e. everything that is not SI.
        !
        use constants
        use su_cgns, only: Kelvin, Celsius, Rankine, Fahrenheit
        implicit none
        !
        !      Subroutine arguments.
        !
        integer, intent(in) :: temp
        real(kind=realType), intent(out) :: mult, trans

        ! Determine the situation we are having here.

        select case (temp)

        case (Kelvin)

            ! Temperature is already given in Kelvin. No need to convert.

            mult = one
            trans = zero

        case (Celsius)      ! is it Celcius or Celsius?

            ! Temperature is in Celsius. Only an offset must be applied.

            mult = one
            trans = 273.16_realType

        case (Rankine)

            ! Temperature is in Rankine. Only a multiplication needs to
            ! be performed.

            mult = 5.0_realType / 9.0_realType
            trans = zero

        case (Fahrenheit)

            ! Temperature is in Fahrenheit. Both a multiplication and an
            ! offset must be applied.

            mult = 5.0_realType / 9.0_realType
            trans = 255.382

        case default

            ! Unknown temperature unit.

            call terminate("siTemperature", &
                           "No idea how to convert this to SI units")

        end select

    end subroutine siTemperature
    subroutine siTurb(mass, len, time, temp, turbName, mult, trans)
        !
        !       siTurb computes the conversion from the given turbulence
        !       unit, which can be constructed from mass, len, time and temp,
        !       to the SI-unit for the given variable. The conversion will
        !       look like: var in SI = mult*(var in NCU) + trans.
        !       NCU means non-christian units, i.e. everything that is not SI.
        !
        use constants
        use su_cgns, only: Kilogram, Meter, Second, Kelvin
        implicit none
        !
        !      Subroutine arguments.
        !
        integer, intent(in) :: mass, len, time, temp
        character(len=*), intent(in) :: turbName
        real(kind=realType), intent(out) :: mult, trans

        ! Determine the situation we are having here.

        if (mass == Kilogram .and. len == Meter .and. &
            time == Second .and. temp == Kelvin) then

            ! Everthing is already in SI units. No conversion needed.

            mult = one
            trans = zero

        else

            call terminate("siTurb", &
                           "No idea how to convert this to SI units")

        end if

    end subroutine siTurb

    subroutine siVelocity(length, time, mult, trans)
        !
        !       siVelocity computes the conversion from the given velocity
        !       unit, which can be constructed from length and time, to the
        !       SI-unit m/s. The conversion will look like:
        !       velocity in m/s = mult*(velocity in ncu) + trans.
        !       Ncu means non-christian units, i.e. everything that is not SI.
        !
        use constants
        use su_cgns, only: Meter, CentiMeter, Millimeter, Foot, Inch, Second
        implicit none
        !
        !      Subroutine arguments.
        !
        integer, intent(in) :: length, time
        real(kind=realType), intent(out) :: mult, trans

        ! Determine the situation we are having here.
        ! First the length.

        select case (length)

        case (Meter)
            mult = one; trans = zero

        case (CenTimeter)
            mult = 0.01_realType; trans = zero

        case (Millimeter)
            mult = 0.001_realType; trans = zero

        case (Foot)
            mult = 0.3048_realType; trans = zero

        case (Inch)
            mult = 0.0254_realType; trans = zero

        case default
            call terminate("siVelocity", &
                           "No idea how to convert this length to SI units")

        end select

        ! And the time.

        select case (time)

        case (Second)
            mult = mult

        case default
            call terminate("siVelocity", &
                           "No idea how to convert this time to SI units")

        end select

    end subroutine siVelocity

    ! ----------------------------------------------------------------------
    !                                                                      |
    !                    No Tapenade Routine below this line               |
    !                                                                      |
    ! ----------------------------------------------------------------------

end module utils
