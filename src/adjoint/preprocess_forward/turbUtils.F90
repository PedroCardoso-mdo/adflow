










module turbUtils

contains

    subroutine prodKatoLaunder
        !
        !       prodKatoLaunder computes the turbulent production term using
        !       the Kato-Launder formulation.
        !
        use constants
        use blockPointers, only: nx, ny, nz, il, jl, kl, w, si, sj, sk, vol, sectionID, scratch
        use flowVarRefState, only: timeRef
        use section, only: sections
        use turbMod, only: prod
        implicit none
        !
        !      Local variables.
        !
        integer(kind=intType) :: i, j, k, ii

        real(kind=realType) :: uux, uuy, uuz, vvx, vvy, vvz, wwx, wwy, wwz
        real(kind=realType) :: qxx, qyy, qzz, qxy, qxz, qyz, sijsij
        real(kind=realType) :: oxy, oxz, oyz, oijoij
        real(kind=realType) :: fact, omegax, omegay, omegaz

        ! Determine the non-dimensional wheel speed of this block.
        ! The vorticity term, which appears in Kato-Launder is of course
        ! not frame invariant. To approximate frame invariance the wheel
        ! speed should be substracted from oxy, oxz and oyz, which results
        ! in the vorticity in the rotating frame. However some people
        ! claim that the absolute vorticity should be used to obtain the
        ! best results. In that omega should be set to zero.

        omegax = timeRef * sections(sectionID)%rotRate(1)
        omegay = timeRef * sections(sectionID)%rotRate(2)
        omegaz = timeRef * sections(sectionID)%rotRate(3)

        ! Loop over the cell centers of the given block. It may be more
        ! efficient to loop over the faces and to scatter the gradient,
        ! but in that case the gradients for u, v and w must be stored.
        ! In the current approach no extra memory is needed.
            do k = 2, kl
                do j = 2, jl
                    do i = 2, il

                        ! Compute the gradient of u in the cell center. Use is made
                        ! of the fact that the surrounding normals sum up to zero,
                        ! such that the cell i,j,k does not give a contribution.
                        ! The gradient is scaled by a factor 2*vol.

                        uux = w(i + 1, j, k, ivx) * si(i, j, k, 1) - w(i - 1, j, k, ivx) * si(i - 1, j, k, 1) &
                              + w(i, j + 1, k, ivx) * sj(i, j, k, 1) - w(i, j - 1, k, ivx) * sj(i, j - 1, k, 1) &
                              + w(i, j, k + 1, ivx) * sk(i, j, k, 1) - w(i, j, k - 1, ivx) * sk(i, j, k - 1, 1)
                        uuy = w(i + 1, j, k, ivx) * si(i, j, k, 2) - w(i - 1, j, k, ivx) * si(i - 1, j, k, 2) &
                              + w(i, j + 1, k, ivx) * sj(i, j, k, 2) - w(i, j - 1, k, ivx) * sj(i, j - 1, k, 2) &
                              + w(i, j, k + 1, ivx) * sk(i, j, k, 2) - w(i, j, k - 1, ivx) * sk(i, j, k - 1, 2)
                        uuz = w(i + 1, j, k, ivx) * si(i, j, k, 3) - w(i - 1, j, k, ivx) * si(i - 1, j, k, 3) &
                              + w(i, j + 1, k, ivx) * sj(i, j, k, 3) - w(i, j - 1, k, ivx) * sj(i, j - 1, k, 3) &
                              + w(i, j, k + 1, ivx) * sk(i, j, k, 3) - w(i, j, k - 1, ivx) * sk(i, j, k - 1, 3)

                        ! Idem for the gradient of v.

                        vvx = w(i + 1, j, k, ivy) * si(i, j, k, 1) - w(i - 1, j, k, ivy) * si(i - 1, j, k, 1) &
                              + w(i, j + 1, k, ivy) * sj(i, j, k, 1) - w(i, j - 1, k, ivy) * sj(i, j - 1, k, 1) &
                              + w(i, j, k + 1, ivy) * sk(i, j, k, 1) - w(i, j, k - 1, ivy) * sk(i, j, k - 1, 1)
                        vvy = w(i + 1, j, k, ivy) * si(i, j, k, 2) - w(i - 1, j, k, ivy) * si(i - 1, j, k, 2) &
                              + w(i, j + 1, k, ivy) * sj(i, j, k, 2) - w(i, j - 1, k, ivy) * sj(i, j - 1, k, 2) &
                              + w(i, j, k + 1, ivy) * sk(i, j, k, 2) - w(i, j, k - 1, ivy) * sk(i, j, k - 1, 2)
                        vvz = w(i + 1, j, k, ivy) * si(i, j, k, 3) - w(i - 1, j, k, ivy) * si(i - 1, j, k, 3) &
                              + w(i, j + 1, k, ivy) * sj(i, j, k, 3) - w(i, j - 1, k, ivy) * sj(i, j - 1, k, 3) &
                              + w(i, j, k + 1, ivy) * sk(i, j, k, 3) - w(i, j, k - 1, ivy) * sk(i, j, k - 1, 3)

                        ! And for the gradient of w.

                        wwx = w(i + 1, j, k, ivz) * si(i, j, k, 1) - w(i - 1, j, k, ivz) * si(i - 1, j, k, 1) &
                              + w(i, j + 1, k, ivz) * sj(i, j, k, 1) - w(i, j - 1, k, ivz) * sj(i, j - 1, k, 1) &
                              + w(i, j, k + 1, ivz) * sk(i, j, k, 1) - w(i, j, k - 1, ivz) * sk(i, j, k - 1, 1)
                        wwy = w(i + 1, j, k, ivz) * si(i, j, k, 2) - w(i - 1, j, k, ivz) * si(i - 1, j, k, 2) &
                              + w(i, j + 1, k, ivz) * sj(i, j, k, 2) - w(i, j - 1, k, ivz) * sj(i, j - 1, k, 2) &
                              + w(i, j, k + 1, ivz) * sk(i, j, k, 2) - w(i, j, k - 1, ivz) * sk(i, j, k - 1, 2)
                        wwz = w(i + 1, j, k, ivz) * si(i, j, k, 3) - w(i - 1, j, k, ivz) * si(i - 1, j, k, 3) &
                              + w(i, j + 1, k, ivz) * sj(i, j, k, 3) - w(i, j - 1, k, ivz) * sj(i, j - 1, k, 3) &
                              + w(i, j, k + 1, ivz) * sk(i, j, k, 3) - w(i, j, k - 1, ivz) * sk(i, j, k - 1, 3)

                        ! Compute the strain and vorticity terms. The multiplication
                        ! is present to obtain the correct gradients. Note that
                        ! the wheel speed is substracted from the vorticity terms.

                        fact = half / vol(i, j, k)

                        qxx = fact * uux
                        qyy = fact * vvy
                        qzz = fact * wwz

                        qxy = fact * half * (uuy + vvx)
                        qxz = fact * half * (uuz + wwx)
                        qyz = fact * half * (vvz + wwy)

                        oxy = fact * half * (vvx - uuy) - omegaz
                        oxz = fact * half * (uuz - wwx) - omegay
                        oyz = fact * half * (wwy - vvz) - omegax

                        ! Compute the summation of the strain and vorticity tensors.

                        sijsij = two * (qxy**2 + qxz**2 + qyz**2) &
                                 + qxx**2 + qyy**2 + qzz**2
                        oijoij = two * (oxy**2 + oxz**2 + oyz**2)

                        ! Compute the production term.

                        scratch(i, j, k, iprod) = two * sqrt(sijsij * oijoij)
                end do
            end do
        end do
    end subroutine prodKatoLaunder

    subroutine prodSmag2
        !
        !       prodSmag2 computes the term:
        !              2*sij*sij - 2/3 div(u)**2 with  sij=0.5*(duidxj+dujdxi)
        !       which is used for the turbulence equations.
        !       It is assumed that the pointer prod, stored in turbMod, is
        !       already set to the correct entry.
        !
        use constants
        use blockPointers, only: nx, ny, nz, il, jl, kl, w, si, sj, sk, vol, sectionID, scratch
        implicit none
        !
        !      Local parameter
        !
        real(kind=realType), parameter :: f23 = two * third
        !
        !      Local variables.
        !
        integer(kind=intType) :: i, j, k, ii
        real(kind=realType) :: uux, uuy, uuz, vvx, vvy, vvz, wwx, wwy, wwz
        real(kind=realType) :: div2, fact, sxx, syy, szz, sxy, sxz, syz

        ! Loop over the cell centers of the given block. It may be more
        ! efficient to loop over the faces and to scatter the gradient,
        ! but in that case the gradients for u, v and w must be stored.
        ! In the current approach no extra memory is needed.

            do k = 2, kl
                do j = 2, jl
                    do i = 2, il

                        ! Compute the gradient of u in the cell center. Use is made
                        ! of the fact that the surrounding normals sum up to zero,
                        ! such that the cell i,j,k does not give a contribution.
                        ! The gradient is scaled by the factor 2*vol.

                        uux = w(i + 1, j, k, ivx) * si(i, j, k, 1) - w(i - 1, j, k, ivx) * si(i - 1, j, k, 1) &
                              + w(i, j + 1, k, ivx) * sj(i, j, k, 1) - w(i, j - 1, k, ivx) * sj(i, j - 1, k, 1) &
                              + w(i, j, k + 1, ivx) * sk(i, j, k, 1) - w(i, j, k - 1, ivx) * sk(i, j, k - 1, 1)
                        uuy = w(i + 1, j, k, ivx) * si(i, j, k, 2) - w(i - 1, j, k, ivx) * si(i - 1, j, k, 2) &
                              + w(i, j + 1, k, ivx) * sj(i, j, k, 2) - w(i, j - 1, k, ivx) * sj(i, j - 1, k, 2) &
                              + w(i, j, k + 1, ivx) * sk(i, j, k, 2) - w(i, j, k - 1, ivx) * sk(i, j, k - 1, 2)
                        uuz = w(i + 1, j, k, ivx) * si(i, j, k, 3) - w(i - 1, j, k, ivx) * si(i - 1, j, k, 3) &
                              + w(i, j + 1, k, ivx) * sj(i, j, k, 3) - w(i, j - 1, k, ivx) * sj(i, j - 1, k, 3) &
                              + w(i, j, k + 1, ivx) * sk(i, j, k, 3) - w(i, j, k - 1, ivx) * sk(i, j, k - 1, 3)

                        ! Idem for the gradient of v.

                        vvx = w(i + 1, j, k, ivy) * si(i, j, k, 1) - w(i - 1, j, k, ivy) * si(i - 1, j, k, 1) &
                              + w(i, j + 1, k, ivy) * sj(i, j, k, 1) - w(i, j - 1, k, ivy) * sj(i, j - 1, k, 1) &
                              + w(i, j, k + 1, ivy) * sk(i, j, k, 1) - w(i, j, k - 1, ivy) * sk(i, j, k - 1, 1)
                        vvy = w(i + 1, j, k, ivy) * si(i, j, k, 2) - w(i - 1, j, k, ivy) * si(i - 1, j, k, 2) &
                              + w(i, j + 1, k, ivy) * sj(i, j, k, 2) - w(i, j - 1, k, ivy) * sj(i, j - 1, k, 2) &
                              + w(i, j, k + 1, ivy) * sk(i, j, k, 2) - w(i, j, k - 1, ivy) * sk(i, j, k - 1, 2)
                        vvz = w(i + 1, j, k, ivy) * si(i, j, k, 3) - w(i - 1, j, k, ivy) * si(i - 1, j, k, 3) &
                              + w(i, j + 1, k, ivy) * sj(i, j, k, 3) - w(i, j - 1, k, ivy) * sj(i, j - 1, k, 3) &
                              + w(i, j, k + 1, ivy) * sk(i, j, k, 3) - w(i, j, k - 1, ivy) * sk(i, j, k - 1, 3)

                        ! And for the gradient of w.

                        wwx = w(i + 1, j, k, ivz) * si(i, j, k, 1) - w(i - 1, j, k, ivz) * si(i - 1, j, k, 1) &
                              + w(i, j + 1, k, ivz) * sj(i, j, k, 1) - w(i, j - 1, k, ivz) * sj(i, j - 1, k, 1) &
                              + w(i, j, k + 1, ivz) * sk(i, j, k, 1) - w(i, j, k - 1, ivz) * sk(i, j, k - 1, 1)
                        wwy = w(i + 1, j, k, ivz) * si(i, j, k, 2) - w(i - 1, j, k, ivz) * si(i - 1, j, k, 2) &
                              + w(i, j + 1, k, ivz) * sj(i, j, k, 2) - w(i, j - 1, k, ivz) * sj(i, j - 1, k, 2) &
                              + w(i, j, k + 1, ivz) * sk(i, j, k, 2) - w(i, j, k - 1, ivz) * sk(i, j, k - 1, 2)
                        wwz = w(i + 1, j, k, ivz) * si(i, j, k, 3) - w(i - 1, j, k, ivz) * si(i - 1, j, k, 3) &
                              + w(i, j + 1, k, ivz) * sj(i, j, k, 3) - w(i, j - 1, k, ivz) * sj(i, j - 1, k, 3) &
                              + w(i, j, k + 1, ivz) * sk(i, j, k, 3) - w(i, j, k - 1, ivz) * sk(i, j, k - 1, 3)

                        ! Compute the components of the stress tensor.
                        ! The combination of the current scaling of the velocity
                        ! gradients (2*vol) and the definition of the stress tensor,
                        ! leads to the factor 1/(4*vol).

                        fact = fourth / vol(i, j, k)

                        sxx = two * fact * uux
                        syy = two * fact * vvy
                        szz = two * fact * wwz

                        sxy = fact * (uuy + vvx)
                        sxz = fact * (uuz + wwx)
                        syz = fact * (vvz + wwy)

                        ! Compute 2/3 * divergence of velocity squared

                        div2 = f23 * (sxx + syy + szz)**2

                        ! Store the square of strain as the production term.

                        scratch(i, j, k, iprod) = two * (two * (sxy**2 + sxz**2 + syz**2) &
                                                         + sxx**2 + syy**2 + szz**2) - div2
                end do
            end do
        end do
    end subroutine prodSmag2

    subroutine prodWmag2
        !
        !       prodWmag2 computes the term:
        !          2*oij*oij  with oij=0.5*(duidxj - dujdxi).
        !       This is equal to the magnitude squared of the vorticity.
        !       It is assumed that the pointer vort, stored in turbMod, is
        !       already set to the correct entry.
        !
        use constants
        use blockPointers, only: nx, ny, nz, il, jl, kl, w, si, sj, sk, vol, sectionID, scratch
        use flowVarRefState, only: timeRef
        use section, only: sections
        implicit none
        !
        !      Local variables.
        !
        integer :: i, j, k, ii

        real(kind=realType) :: uuy, uuz, vvx, vvz, wwx, wwy
        real(kind=realType) :: fact, vortx, vorty, vortz
        real(kind=realType) :: omegax, omegay, omegaz

        ! Determine the non-dimensional wheel speed of this block.

        omegax = timeRef * sections(sectionID)%rotRate(1)
        omegay = timeRef * sections(sectionID)%rotRate(2)
        omegaz = timeRef * sections(sectionID)%rotRate(3)

        ! Loop over the cell centers of the given block. It may be more
        ! efficient to loop over the faces and to scatter the gradient,
        ! but in that case the gradients for u, v and w must be stored.
        ! In the current approach no extra memory is needed.
            do k = 2, kl
                do j = 2, jl
                    do i = 2, il

                        ! Compute the necessary derivatives of u in the cell center.
                        ! Use is made of the fact that the surrounding normals sum up
                        ! to zero, such that the cell i,j,k does not give a
                        ! contribution. The gradient is scaled by a factor 2*vol.

                        uuy = w(i + 1, j, k, ivx) * si(i, j, k, 2) - w(i - 1, j, k, ivx) * si(i - 1, j, k, 2) &
                              + w(i, j + 1, k, ivx) * sj(i, j, k, 2) - w(i, j - 1, k, ivx) * sj(i, j - 1, k, 2) &
                              + w(i, j, k + 1, ivx) * sk(i, j, k, 2) - w(i, j, k - 1, ivx) * sk(i, j, k - 1, 2)
                        uuz = w(i + 1, j, k, ivx) * si(i, j, k, 3) - w(i - 1, j, k, ivx) * si(i - 1, j, k, 3) &
                              + w(i, j + 1, k, ivx) * sj(i, j, k, 3) - w(i, j - 1, k, ivx) * sj(i, j - 1, k, 3) &
                              + w(i, j, k + 1, ivx) * sk(i, j, k, 3) - w(i, j, k - 1, ivx) * sk(i, j, k - 1, 3)

                        ! Idem for the gradient of v.

                        vvx = w(i + 1, j, k, ivy) * si(i, j, k, 1) - w(i - 1, j, k, ivy) * si(i - 1, j, k, 1) &
                              + w(i, j + 1, k, ivy) * sj(i, j, k, 1) - w(i, j - 1, k, ivy) * sj(i, j - 1, k, 1) &
                              + w(i, j, k + 1, ivy) * sk(i, j, k, 1) - w(i, j, k - 1, ivy) * sk(i, j, k - 1, 1)
                        vvz = w(i + 1, j, k, ivy) * si(i, j, k, 3) - w(i - 1, j, k, ivy) * si(i - 1, j, k, 3) &
                              + w(i, j + 1, k, ivy) * sj(i, j, k, 3) - w(i, j - 1, k, ivy) * sj(i, j - 1, k, 3) &
                              + w(i, j, k + 1, ivy) * sk(i, j, k, 3) - w(i, j, k - 1, ivy) * sk(i, j, k - 1, 3)

                        ! And for the gradient of w.

                        wwx = w(i + 1, j, k, ivz) * si(i, j, k, 1) - w(i - 1, j, k, ivz) * si(i - 1, j, k, 1) &
                              + w(i, j + 1, k, ivz) * sj(i, j, k, 1) - w(i, j - 1, k, ivz) * sj(i, j - 1, k, 1) &
                              + w(i, j, k + 1, ivz) * sk(i, j, k, 1) - w(i, j, k - 1, ivz) * sk(i, j, k - 1, 1)
                        wwy = w(i + 1, j, k, ivz) * si(i, j, k, 2) - w(i - 1, j, k, ivz) * si(i - 1, j, k, 2) &
                              + w(i, j + 1, k, ivz) * sj(i, j, k, 2) - w(i, j - 1, k, ivz) * sj(i, j - 1, k, 2) &
                              + w(i, j, k + 1, ivz) * sk(i, j, k, 2) - w(i, j, k - 1, ivz) * sk(i, j, k - 1, 2)

                        ! Compute the three components of the vorticity vector.
                        ! Substract the part coming from the rotating frame.

                        fact = half / vol(i, j, k)

                        vortx = fact * (wwy - vvz) - two * omegax
                        vorty = fact * (uuz - wwx) - two * omegay
                        vortz = fact * (vvx - uuy) - two * omegaz

                        ! Compute the magnitude squared of the vorticity.

                        scratch(i, j, k, ivort) = vortx**2 + vorty**2 + vortz**2
                end do
            end do
        end do
    end subroutine prodWmag2
    function saNuKnownEddyRatio(eddyRatio, nuLam)
        !
        !       saNuKnownEddyRatio computes the Spalart-Allmaras transport
        !       variable nu for the given eddy viscosity ratio.
        !
        use constants
        use paramTurb
        implicit none
        !
        !      Function type.
        !
        real(kind=realType) :: saNuKnownEddyRatio
        !
        !      Function arguments.
        !
        real(kind=realType), intent(in) :: eddyRatio, nuLam
        !
        !      Local variables.
        !
        real(kind=realType) :: cv13, chi, chi2, chi3, chi4, f, df, dchi

        ! Take care of the exceptional cases.

        if (eddyRatio <= zero) then
            saNuKnownEddyRatio = zero
            return
        end if

        ! Set the value of cv1^3, which is the constant appearing in the
        ! sa function fv1 to compute the eddy viscosity

        cv13 = rsaCv1**3

        ! Determine the value of chi, which is given by the quartic
        ! polynomial chi^4 - ratio*(chi^3 + cv1^3) = 0.
        ! First determine the start value, depending on the eddyRatio.

        if (eddyRatio < 1.e-4_realType) then
            chi = 0.5_realType
        else if (eddyRatio < 1.0_realType) then
            chi = 5.0_realType
        else if (eddyRatio < 10.0_realType) then
            chi = 10.0_realType
        else
            chi = eddyRatio
        end if

        ! The actual newton algorithm.

        do
            ! Compute the function value and the derivative.

            chi2 = chi * chi
            chi3 = chi * chi2
            chi4 = chi * chi3

            f = chi4 - eddyRatio * (chi3 + cv13)
            df = four * chi3 - three * eddyRatio * chi2

            ! Compute the negative update and the new value of chi.

            dchi = f / df
            chi = chi - dchi

            ! Condition to exit the loop.

            if (abs(dchi / chi) <= thresholdReal) exit
        end do

        ! Chi is the ratio of the spalart allmaras transport variable and
        ! the laminar viscosity. So multiply chi with the laminar viscosity
        ! to obtain the correct value.

        saNuKnownEddyRatio = nuLam * chi

    end function saNuKnownEddyRatio

    subroutine unsteadyTurbTerm(mAdv, nAdv, offset, qq)
        !
        !       unsteadyTurbTerm discretizes the time derivative of the
        !       turbulence transport equations and add it to the residual.
        !       As the time derivative is the same for all turbulence models,
        !       this generic routine can be used; both the discretization of
        !       the time derivative and its contribution to the central
        !       jacobian are computed by this routine.
        !       Only nAdv equations are treated, while the actual system has
        !       size mAdv. The reason is that some equations for some
        !       turbulence equations do not have a time derivative, e.g. the
        !       f equation in the v2-f model. The argument offset indicates
        !       the offset in the w vector where this subsystem starts. As a
        !       consequence it is assumed that the indices of the current
        !       subsystem are contiguous, e.g. if a 2*2 system is solved the
        !       Last index in w is offset+1 and offset+2 respectively.
        !
        use blockPointers
        use flowVarRefState
        use inputPhysics
        use inputTimeSpectral
        use inputUnsteady
        use iteration
        use section
        use turbMod

        implicit none
        !
        !      Subroutine arguments.
        !
        integer(kind=intType), intent(in) :: mAdv, nAdv, offset

        real(kind=realType), dimension(2:il, 2:jl, 2:kl, mAdv, mAdv), &
            intent(inout) :: qq
        !
        !      Local variables.
        !
        integer(kind=intType) :: i, j, k, ii, jj, nn

        real(kind=realType) :: oneOverDt, tmp

        ! Determine the equation mode.

        select case (equationMode)

        case (steady)

            ! Steady computation. No time derivative present.

            return

            !===============================================================

        case (unsteady)

            ! The time deritvative term depends on the integration
            ! scheme used.

            select case (timeIntegrationScheme)

            case (BDF)

                ! Backward difference formula is used as time
                ! integration scheme.

                ! Store the inverse of the physical nonDimensional
                ! time step a bit easier.

                oneOverDt = timeRef / deltaT

                ! Loop over the number of turbulent transport equations.

                nAdvLoopUnsteady: do ii = 1, nAdv

                    ! Store the index of the current turbulent variable in jj.

                    jj = ii + offset

                    ! Loop over the owned cells of this block to compute the
                    ! time derivative.

                    do k = 2, kl
                        do j = 2, jl
                            do i = 2, il

                                ! Initialize tmp to the value of the current
                                ! level multiplied by the corresponding coefficient
                                ! in the time integration scheme.

                                tmp = coefTime(0) * w(i, j, k, jj)

                                ! Loop over the old time levels and add the
                                ! corresponding contribution to tmp.

                                do nn = 1, noldLevels
                                    tmp = tmp + coefTime(nn) * wold(nn, i, j, k, jj)
                                end do

                                ! Update the residual. Note that in the turbulent
                                ! routines the residual is defined with an opposite
                                ! sign compared to the residual of the flow equations.
                                ! Therefore the time derivative must be substracted
                                ! from dvt.

                                scratch(i, j, k, idvt + ii - 1) = scratch(i, j, k, idvt + ii - 1) - oneOverDt * tmp

                                ! Update the central jacobian.

                                qq(i, j, k, ii, ii) = qq(i, j, k, ii, ii) &
                                                      + coefTime(0) * oneOverDt
                            end do
                        end do
                    end do

                end do nAdvLoopUnsteady

                !===========================================================

            case (explicitRK)

                ! Explicit time integration scheme. The time derivative
                ! is handled differently.

                return

            end select

            !===============================================================

        case (timeSpectral)

            ! Time spectral method.

            ! Loop over the number of turbulent transport equations.

            nAdvLoopSpectral: do ii = 1, nAdv

                ! Store the index of the current turbulent variable in jj.

                jj = ii + offset

                ! The time derivative has been computed earlier in
                ! unsteadyTurbSpectral and stored in entry jj of scratch.
                ! Substract this value for all owned cells. It must be
                ! substracted, because in the turbulent routines the
                ! residual is defined with an opposite sign compared to
                ! the residual of the flow equations.
                ! Also add a term to the diagonal matrix, which corresponds
                ! to to the contribution of the highest frequency. This is
                ! equivalent to an explicit treatment of the time derivative
                ! and may need to be changed.

                tmp = nTimeIntervalsSpectral * pi * timeRef &
                      / sections(sectionID)%timePeriod

                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il
                            scratch(i, j, k, idvt + ii - 1) = scratch(i, j, k, idvt + ii - 1) - dw(i, j, k, jj)
                            qq(i, j, k, ii, ii) = qq(i, j, k, ii, ii) + tmp
                        end do
                    end do
                end do

            end do nAdvLoopSpectral

        end select

    end subroutine unsteadyTurbTerm

    subroutine computeEddyViscosity(includeHalos)
        !
        !       computeEddyViscosity computes the eddy viscosity in the
        !       owned cell centers of the given block. It is assumed that the
        !       pointes already point to the correct block before entering
        !       this subroutine.
        !
        use constants
        use flowVarRefState
        use inputPhysics
        use iteration
        use blockPointers
        implicit none

        ! Input Parameter
        logical, intent(in) :: includeHalos

        !
        !      Local variables.
        !
        logical :: returnImmediately
        integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd

        ! Check if an immediate return can be made.

        if (eddyModel) then
            if ((currentLevel <= groundLevel)) then
                returnImmediately = .false.
            else
                returnImmediately = .true.
            end if
        else
            returnImmediately = .true.
        end if

        if (returnImmediately) return

        ! Determine the turbulence model and call the appropriate
        ! routine to compute the eddy viscosity.
        if (includeHalos) then
            iBeg = 1
            iEnd = ie
            jBeg = 1
            jEnd = je
            kBeg = 1
            kEnd = ke
        else
            iBeg = 2
            iEnd = il
            jBeg = 2
            jEnd = jl
            kBeg = 2
            kEnd = kl
        end if

        select case (turbModel)

        case (spalartAllmaras, spalartAllmarasEdwards, spalartallmarasnoft2gammaretheta)
            call saEddyViscosity(iBeg, iEnd, jBeg, jEnd, kBeg, kEnd)
        end select

    end subroutine computeEddyViscosity

    subroutine saEddyViscosity(iBeg, iEnd, jBeg, jEnd, kBeg, kEnd)
        !
        !       saEddyViscosity computes the eddy-viscosity according to the
        !       Spalart-Allmaras model for the block given in blockPointers.
        !       This routine for both the original version as well as the
        !       modified version according to Edwards.
        !
        use constants
        use blockPointers
        use constants
        use paramTurb
        implicit none
        ! Input variables
        integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd

        !
        !      Local variables.
        !
        integer(kind=intType) :: i, j, k, ii, iSize, jSize, kSize
        real(kind=realType) :: chi, chi3, fv1, rnuSA, cv13

        ! Store the cv1^3; cv1 is a constant of the Spalart-Allmaras model.

        cv13 = rsaCv1**3

        ! Loop over the cells of this block and compute the eddy viscosity.
        ! Do not include halo's.
            do k = kBeg, kEnd
                do j = jBeg, jEnd
                    do i = iBeg, iEnd
                        rnuSA = w(i, j, k, itu1) * w(i, j, k, irho)
                        chi = rnuSA / rlv(i, j, k)
                        chi3 = chi**3
                        fv1 = chi3 / (chi3 + cv13)
                        rev(i, j, k) = fv1 * rnuSA
                end do
            end do
        end do
    end subroutine saEddyViscosity

    subroutine kwEddyViscosity(iBeg, iEnd, jBeg, jEnd, kBeg, kEnd)
        !
        !       kwEddyViscosity computes the eddy viscosity according to the
        !       k-omega models (both the original Wilcox as well as the
        !       modified version) for the block given in blockPointers.
        !
        use constants
        use blockPointers
        implicit none
        ! Input variables
        integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd
        !
        !      Local variables.
        !
        integer(kind=intType) :: i, j, k, ii, iSize, jSize, kSize

        ! Loop over the cells of this block and compute the eddy viscosity.
        ! Do not include halo's.
            do k = kBeg, kEnd
                do j = jBeg, jEnd
                    do i = iBeg, iEnd
                        rev(i, j, k) = abs(w(i, j, k, irho) * w(i, j, k, itu1) / w(i, j, k, itu2))
                end do
            end do
        end do

    end subroutine kwEddyViscosity

    subroutine SSTEddyViscosity(iBeg, iEnd, jBeg, jEnd, kBeg, kEnd)
        !
        !       SSTEddyViscosity computes the eddy viscosity according to
        !       menter's SST variant of the k-omega turbulence model for the
        !       block given in blockPointers.
        !
        use constants
        use blockPointers
        use paramTurb
        use turbMod
        implicit none
        ! Input variables
        integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd
        !
        !      Local variables.
        !
        integer(kind=intType) :: i, j, k, ii, iSize, jSize, kSize
        real(kind=realType) :: t1, t2, arg2, f2, vortMag

        ! Compute the vorticity squared in the cell centers. The reason
        ! for computing the vorticity squared is that a routine exists
        ! for it; for the actual eddy viscosity computation the vorticity
        ! itself is needed.

        call prodWmag2

        ! Loop over the cells of this block and compute the eddy viscosity.
        ! Do not include halo's.
            do k = kBeg, kEnd
                do j = jBeg, jEnd
                    do i = iBeg, iEnd
                        ! Compute the value of the function f2, which occurs in the
                        ! eddy-viscosity computation.

                        t1 = two * sqrt(w(i, j, k, itu1)) &
                             / (0.09_realType * w(i, j, k, itu2) * d2Wall(i, j, k))
                        t2 = 500.0_realType * rlv(i, j, k) &
                             / (w(i, j, k, irho) * w(i, j, k, itu2) * d2Wall(i, j, k)**2)

                        arg2 = max(t1, t2)
                        f2 = tanh(arg2**2)

                        ! And compute the eddy viscosity.

                        vortMag = sqrt(scratch(i, j, k, iprod))
                        rev(i, j, k) = w(i, j, k, irho) * rSSTA1 * w(i, j, k, itu1) &
                                       / max(rSSTA1 * w(i, j, k, itu2), f2 * vortMag)
                end do
            end do
        end do

    end subroutine SSTEddyViscosity

    subroutine turbAdvection(mAdv, nAdv, offset, qq)
        !
        !       turbAdvection discretizes the advection part of the turbulent
        !       transport equations. As the advection part is the same for all
        !       models, this generic routine can be used. Both the
        !       discretization and the central jacobian are computed in this
        !       subroutine. The former can either be 1st or 2nd order
        !       accurate; the latter is always based on the 1st order upwind
        !       discretization. When the discretization must be second order
        !       accurate, the fully upwind (kappa = -1) scheme in combination
        !       with the minmod limiter is used.
        !       Only nAdv equations are treated, while the actual system has
        !       size mAdv. The reason is that some equations for some
        !       turbulence equations do not have an advection part, e.g. the
        !       f equation in the v2-f model. The argument offset indicates
        !       the offset in the w vector where this subsystem starts. As a
        !       consequence it is assumed that the indices of the current
        !       subsystem are contiguous, e.g. if a 2*2 system is solved the
        !       Last index in w is offset+1 and offset+2 respectively.
        !
        use constants
        use blockPointers, only: nx, ny, nz, il, jl, kl, vol, sfaceI, sfaceJ, sfaceK, &
                                 w, si, sj, sk, addGridVelocities, bmti1, bmti2, bmtj1, bmtj2, &
                                 bmtk1, bmtk2, scratch
        use inputDiscretization, only: orderTurb
        use iteration, only: groundLevel
        use turbMod, only: secondOrd
        implicit none
        !
        !      Subroutine arguments.
        !
        integer(kind=intType), intent(in) :: nAdv, mAdv, offset

        real(kind=realType), dimension(2:il, 2:jl, 2:kl, mAdv, mAdv), &
            intent(inout) :: qq
        !
        !      Local variables.
        !
        integer(kind=intType) :: i, j, k, ii, jj, kk, iii

        real(kind=realType) :: qs, voli, xa, ya, za
        real(kind=realType) :: uu, dwt, dwtm1, dwtp1, dwti, dwtj, dwtk

        real(kind=realType), dimension(mAdv) :: impl

        ! Determine whether or not a second order discretization for the
        ! advective terms must be used.
        secondOrd = .false.
        if (groundLevel == 1_intType .and. &
            orderTurb == secondOrder) secondOrd = .true.

        ! Initialize the grid velocity to zero. This value will be used
        ! if the block is not moving.
        continue
        !$AD CHECKPOINT-START
        qs = zero
        !
        !       Upwind discretization of the convective term in k (zeta)
        !       direction. Either the 1st order upwind or the second order
        !       fully upwind interpolation scheme, kappa = -1, is used in
        !       combination with the minmod limiter.
        !       The possible grid velocity must be taken into account.
        !
            do k = 2, kl
                do j = 2, jl
                    do i = 2, il
                        ! Compute the grid velocity if present.
                        ! It is taken as the average of k and k-1,

                        voli = half / vol(i, j, k)
                        if (addGridVelocities) &
                            qs = (sFaceK(i, j, k) + sFaceK(i, j, k - 1)) * voli

                        ! Compute the normal velocity, where the normal direction
                        ! is taken as the average of faces k and k-1.

                        xa = (sk(i, j, k, 1) + sk(i, j, k - 1, 1)) * voli
                        ya = (sk(i, j, k, 2) + sk(i, j, k - 1, 2)) * voli
                        za = (sk(i, j, k, 3) + sk(i, j, k - 1, 3)) * voli

                        uu = xa * w(i, j, k, ivx) + ya * w(i, j, k, ivy) + za * w(i, j, k, ivz) - qs
                        ! This term has unit: velocity/length

                        ! Determine the situation we are having here, i.e. positive
                        ! or negative normal velocity.

                        velKdir: if (uu > zero) then

                            ! Velocity has a component in positive k-direction.
                            ! Loop over the number of advection equations.
                            !$AD II-LOOP
                            do ii = 1, nAdv

                                ! Set the value of jj such that it corresponds to the
                                ! turbulent entry in w.

                                jj = ii + offset

                                ! Check whether a first or a second order discretization
                                ! must be used.

                                if (secondOrd) then

                                    ! Second order; store the three differences for the
                                    ! discretization of the derivative in k-direction.

                                    dwtm1 = w(i, j, k - 1, jj) - w(i, j, k - 2, jj)
                                    dwt = w(i, j, k, jj) - w(i, j, k - 1, jj)
                                    dwtp1 = w(i, j, k + 1, jj) - w(i, j, k, jj)

                                    ! Construct the derivative in this cell center. This
                                    ! is the first order upwind derivative with two
                                    ! nonlinear corrections.

                                    dwtk = dwt

                                    if (dwt * dwtp1 > zero) then
                                        if (abs(dwt) < abs(dwtp1)) then
                                            dwtk = dwtk + half * dwt
                                        else
                                            dwtk = dwtk + half * dwtp1
                                        end if
                                    end if

                                    if (dwt * dwtm1 > zero) then
                                        if (abs(dwt) < abs(dwtm1)) then
                                            dwtk = dwtk - half * dwt
                                        else
                                            dwtk = dwtk - half * dwtm1
                                        end if
                                    end if

                                else

                                    ! 1st order upwind scheme.

                                    dwtk = w(i, j, k, jj) - w(i, j, k - 1, jj)

                                end if

                                ! Update the residual. The convective term must be
                                ! substracted, because it appears on the other side of
                                ! the equation as the source and viscous terms.
                                ! uu*dwtk = (V.dot.face_normal)*delta(nuTilde)/delta(x)

                                scratch(i, j, k, idvt + ii - 1) = scratch(i, j, k, idvt + ii - 1) - uu * dwtk

                            end do

                            else velKdir

                            ! Velocity has a component in negative k-direction.
                            ! Loop over the number of advection equations
                            !$AD II-LOOP
                            do ii = 1, nAdv

                                ! Set the value of jj such that it corresponds to the
                                ! turbulent entry in w.

                                jj = ii + offset

                                ! Check whether a first or a second order discretization
                                ! must be used.

                                if (secondOrd) then

                                    ! Store the three differences for the discretization of
                                    ! the derivative in k-direction.

                                    dwtm1 = w(i, j, k, jj) - w(i, j, k - 1, jj)
                                    dwt = w(i, j, k + 1, jj) - w(i, j, k, jj)
                                    dwtp1 = w(i, j, k + 2, jj) - w(i, j, k + 1, jj)

                                    ! Construct the derivative in this cell center. This is
                                    ! the first order upwind derivative with two nonlinear
                                    ! corrections.

                                    dwtk = dwt

                                    if (dwt * dwtp1 > zero) then
                                        if (abs(dwt) < abs(dwtp1)) then
                                            dwtk = dwtk - half * dwt
                                        else
                                            dwtk = dwtk - half * dwtp1
                                        end if
                                    end if

                                    if (dwt * dwtm1 > zero) then
                                        if (abs(dwt) < abs(dwtm1)) then
                                            dwtk = dwtk + half * dwt
                                        else
                                            dwtk = dwtk + half * dwtm1
                                        end if
                                    end if

                                else

                                    ! 1st order upwind scheme.

                                    dwtk = w(i, j, k + 1, jj) - w(i, j, k, jj)

                                end if

                                ! Update the residual. The convective term must be
                                ! substracted, because it appears on the other side
                                ! of the equation as the source and viscous terms.

                                scratch(i, j, k, idvt + ii - 1) = scratch(i, j, k, idvt + ii - 1) - uu * dwtk

                                ! Update the central jacobian. First the term which is
                                ! always present, i.e. -uu.
                            end do

                        end if velKdir
                end do
            end do
        end do
        continue
        !$AD CHECKPOINT-END
        !
        !       Upwind discretization of the convective term in j (eta)
        !       direction. Either the 1st order upwind or the second order
        !       fully upwind interpolation scheme, kappa = -1, is used in
        !       combination with the minmod limiter.
        !       The possible grid velocity must be taken into account.
        !
        continue
        !$AD CHECKPOINT-START
        qs = zero
            do k = 2, kl
                do j = 2, jl
                    do i = 2, il

                        ! Compute the grid velocity if present.
                        ! It is taken as the average of j and j-1,

                        voli = half / vol(i, j, k)
                        if (addGridVelocities) &
                            qs = (sFaceJ(i, j, k) + sFaceJ(i, j - 1, k)) * voli

                        ! Compute the normal velocity, where the normal direction
                        ! is taken as the average of faces j and j-1.

                        xa = (sj(i, j, k, 1) + sj(i, j - 1, k, 1)) * voli
                        ya = (sj(i, j, k, 2) + sj(i, j - 1, k, 2)) * voli
                        za = (sj(i, j, k, 3) + sj(i, j - 1, k, 3)) * voli

                        uu = xa * w(i, j, k, ivx) + ya * w(i, j, k, ivy) + za * w(i, j, k, ivz) - qs

                        ! Determine the situation we are having here, i.e. positive
                        ! or negative normal velocity.

                        velJdir: if (uu > zero) then

                            ! Velocity has a component in positive j-direction.
                            ! Loop over the number of advection equations.
                            !$AD II-LOOP
                            do ii = 1, nAdv

                                ! Set the value of jj such that it corresponds to the
                                ! turbulent entry in w.

                                jj = ii + offset

                                ! Check whether a first or a second order discretization
                                ! must be used.

                                if (secondOrd) then

                                    ! Second order; store the three differences for the
                                    ! discretization of the derivative in j-direction.

                                    dwtm1 = w(i, j - 1, k, jj) - w(i, j - 2, k, jj)
                                    dwt = w(i, j, k, jj) - w(i, j - 1, k, jj)
                                    dwtp1 = w(i, j + 1, k, jj) - w(i, j, k, jj)

                                    ! Construct the derivative in this cell center. This is
                                    ! the first order upwind derivative with two nonlinear
                                    ! corrections.

                                    dwtj = dwt

                                    if (dwt * dwtp1 > zero) then
                                        if (abs(dwt) < abs(dwtp1)) then
                                            dwtj = dwtj + half * dwt
                                        else
                                            dwtj = dwtj + half * dwtp1
                                        end if
                                    end if

                                    if (dwt * dwtm1 > zero) then
                                        if (abs(dwt) < abs(dwtm1)) then
                                            dwtj = dwtj - half * dwt
                                        else
                                            dwtj = dwtj - half * dwtm1
                                        end if
                                    end if

                                else

                                    ! 1st order upwind scheme.

                                    dwtj = w(i, j, k, jj) - w(i, j - 1, k, jj)

                                end if

                                ! Update the residual. The convective term must be
                                ! substracted, because it appears on the other side of
                                ! the equation as the source and viscous terms.

                                scratch(i, j, k, idvt + ii - 1) = scratch(i, j, k, idvt + ii - 1) - uu * dwtj

                                ! Update the central jacobian. First the term which is
                                ! always present, i.e. uu.
                            end do

                            else velJdir

                            ! Velocity has a component in negative j-direction.
                            ! Loop over the number of advection equations.
                            !$AD II-LOOP
                            do ii = 1, nAdv

                                ! Set the value of jj such that it corresponds to the
                                ! turbulent entry in w.

                                jj = ii + offset

                                ! Check whether a first or a second order discretization
                                ! must be used.

                                if (secondOrd) then

                                    ! Store the three differences for the discretization of
                                    ! the derivative in j-direction.

                                    dwtm1 = w(i, j, k, jj) - w(i, j - 1, k, jj)
                                    dwt = w(i, j + 1, k, jj) - w(i, j, k, jj)
                                    dwtp1 = w(i, j + 2, k, jj) - w(i, j + 1, k, jj)

                                    ! Construct the derivative in this cell center. This is
                                    ! the first order upwind derivative with two nonlinear
                                    ! corrections.

                                    dwtj = dwt

                                    if (dwt * dwtp1 > zero) then
                                        if (abs(dwt) < abs(dwtp1)) then
                                            dwtj = dwtj - half * dwt
                                        else
                                            dwtj = dwtj - half * dwtp1
                                        end if
                                    end if

                                    if (dwt * dwtm1 > zero) then
                                        if (abs(dwt) < abs(dwtm1)) then
                                            dwtj = dwtj + half * dwt
                                        else
                                            dwtj = dwtj + half * dwtm1
                                        end if
                                    end if

                                else

                                    ! 1st order upwind scheme.

                                    dwtj = w(i, j + 1, k, jj) - w(i, j, k, jj)

                                end if

                                ! Update the residual. The convective term must be
                                ! substracted, because it appears on the other side
                                ! of the equation as the source and viscous terms.

                                scratch(i, j, k, idvt + ii - 1) = scratch(i, j, k, idvt + ii - 1) - uu * dwtj

                                ! Update the central jacobian. First the term which is
                                ! always present, i.e. -uu.
                            end do

                        end if velJdir
                end do
            end do
        end do
        continue
        !$AD CHECKPOINT-END
        !
        !       Upwind discretization of the convective term in i (xi)
        !       direction. Either the 1st order upwind or the second order
        !       fully upwind interpolation scheme, kappa = -1, is used in
        !       combination with the minmod limiter.
        !       The possible grid velocity must be taken into account.
        !
        continue
        !$AD CHECKPOINT-START
        qs = zero
            do k = 2, kl
                do j = 2, jl
                    do i = 2, il

                        ! Compute the grid velocity if present.
                        ! It is taken as the average of i and i-1,

                        voli = half / vol(i, j, k)
                        if (addGridVelocities) &
                            qs = (sFaceI(i, j, k) + sFaceI(i - 1, j, k)) * voli

                        ! Compute the normal velocity, where the normal direction
                        ! is taken as the average of faces i and i-1.

                        xa = (si(i, j, k, 1) + si(i - 1, j, k, 1)) * voli
                        ya = (si(i, j, k, 2) + si(i - 1, j, k, 2)) * voli
                        za = (si(i, j, k, 3) + si(i - 1, j, k, 3)) * voli

                        uu = xa * w(i, j, k, ivx) + ya * w(i, j, k, ivy) + za * w(i, j, k, ivz) - qs

                        ! Determine the situation we are having here, i.e. positive
                        ! or negative normal velocity.

                        velIdir: if (uu > zero) then

                            ! Velocity has a component in positive i-direction.
                            ! Loop over the number of advection equations.
                            !$AD II-LOOP
                            do ii = 1, nAdv

                                ! Set the value of jj such that it corresponds to the
                                ! turbulent entry in w.

                                jj = ii + offset

                                ! Check whether a first or a second order discretization
                                ! must be used.

                                if (secondOrd) then

                                    ! Second order; store the three differences for the
                                    ! discretization of the derivative in i-direction.

                                    dwtm1 = w(i - 1, j, k, jj) - w(i - 2, j, k, jj)
                                    dwt = w(i, j, k, jj) - w(i - 1, j, k, jj)
                                    dwtp1 = w(i + 1, j, k, jj) - w(i, j, k, jj)

                                    ! Construct the derivative in this cell center. This is
                                    ! the first order upwind derivative with two nonlinear
                                    ! corrections.

                                    dwti = dwt

                                    if (dwt * dwtp1 > zero) then
                                        if (abs(dwt) < abs(dwtp1)) then
                                            dwti = dwti + half * dwt
                                        else
                                            dwti = dwti + half * dwtp1
                                        end if
                                    end if

                                    if (dwt * dwtm1 > zero) then
                                        if (abs(dwt) < abs(dwtm1)) then
                                            dwti = dwti - half * dwt
                                        else
                                            dwti = dwti - half * dwtm1
                                        end if
                                    end if

                                else

                                    ! 1st order upwind scheme.

                                    dwti = w(i, j, k, jj) - w(i - 1, j, k, jj)

                                end if

                                ! Update the residual. The convective term must be
                                ! substracted, because it appears on the other side of
                                ! the equation as the source and viscous terms.

                                scratch(i, j, k, idvt + ii - 1) = scratch(i, j, k, idvt + ii - 1) - uu * dwti

                                ! Update the central jacobian. First the term which is
                                ! always present, i.e. uu.
                            end do

                            else velIdir

                            ! Velocity has a component in negative i-direction.
                            ! Loop over the number of advection equations.
                            !$AD II-LOOP
                            do ii = 1, nAdv

                                ! Set the value of jj such that it corresponds to the
                                ! turbulent entry in w.

                                jj = ii + offset

                                ! Check whether a first or a second order discretization
                                ! must be used.

                                if (secondOrd) then

                                    ! Second order; store the three differences for the
                                    ! discretization of the derivative in i-direction.

                                    dwtm1 = w(i, j, k, jj) - w(i - 1, j, k, jj)
                                    dwt = w(i + 1, j, k, jj) - w(i, j, k, jj)
                                    dwtp1 = w(i + 2, j, k, jj) - w(i + 1, j, k, jj)

                                    ! Construct the derivative in this cell center. This is
                                    ! the first order upwind derivative with two nonlinear
                                    ! corrections.

                                    dwti = dwt

                                    if (dwt * dwtp1 > zero) then
                                        if (abs(dwt) < abs(dwtp1)) then
                                            dwti = dwti - half * dwt
                                        else
                                            dwti = dwti - half * dwtp1
                                        end if
                                    end if

                                    if (dwt * dwtm1 > zero) then
                                        if (abs(dwt) < abs(dwtm1)) then
                                            dwti = dwti + half * dwt
                                        else
                                            dwti = dwti + half * dwtm1
                                        end if
                                    end if

                                else

                                    ! 1st order upwind scheme.

                                    dwti = w(i + 1, j, k, jj) - w(i, j, k, jj)

                                end if

                                ! Update the residual. The convective term must be
                                ! substracted, because it appears on the other side
                                ! of the equation as the source and viscous terms.

                                scratch(i, j, k, idvt + ii - 1) = scratch(i, j, k, idvt + ii - 1) - uu * dwti

                                ! Update the central jacobian. First the term which is
                                ! always present, i.e. -uu.
                            end do

                        end if velIdir
                end do
            end do
        end do

        !$AD CHECKPOINT-END
        continue
    end subroutine turbAdvection

    ! ----------------------------------------------------------------------
    !                                                                      |
    !                    No Tapenade Routine below this line               |
    !                                                                      |
    ! ----------------------------------------------------------------------


    function reThetaTCorrelation(Tu, lambdaTheta) result(reThetaT)
        !
        !       Compute the critical momentum-thickness Reynolds number
        !       Re_theta_t from the Langtry-Menter correlation.
        !       Uses smooth F(lambda_theta) from Eqs. 54-57.
        !
        !       Tu          : freestream turbulence intensity in percent
        !       lambdaTheta : pressure-gradient parameter (0 for uniform inflow)
        !
        use constants, only: realType, one
        implicit none

        real(kind=realType), intent(in) :: Tu, lambdaTheta
        real(kind=realType) :: reThetaT

        real(kind=realType) :: Flambda, F1val, F2val, F3val, Tu_safe

        Tu_safe = max(Tu, 0.027_realType)

        ! --- Smooth F(lambda_theta) Eqs. 54-57 ---
        ! Eq. 54: F1 = 1 + 0.275*(1 - exp(-35*lam))*exp(-Tu/0.5)
        F1val = one + 0.275_realType &
                * (one - exp(-35.0_realType * lambdaTheta)) &
                * exp(-Tu_safe / 0.5_realType)

        ! Eq. 56: F2 = smoothMax(F1, 1)
        F2val = smoothMinMax(F1val, one, 300.0_realType)

        ! Eq. 55: F3 = 1 - (-12.986*lam - 123.66*lam^2 - 405.689*lam^3)*exp(-(Tu/1.5)^1.5)
        F3val = one - (-12.986_realType * lambdaTheta &
                       - 123.66_realType * lambdaTheta**2 &
                       - 405.689_realType * lambdaTheta**3) &
                * exp(-(Tu_safe / 1.5_realType)**1.5_realType)

        ! Eq. 57: F(lambda) = smoothMin(F2, F3)
        Flambda = smoothMinMax(F2val, F3val, -300.0_realType)

        ! --- Re_theta_t(Tu) * F(lambda_theta) ---
        if (Tu_safe <= 1.3_realType) then
            reThetaT = (1173.51_realType &
                        - 589.428_realType * Tu_safe &
                        + 0.2196_realType / (Tu_safe**2)) * Flambda
        else
            reThetaT = 331.50_realType &
                       * (Tu_safe - 0.5658_realType)**(-0.671_realType) * Flambda
        end if

    end function reThetaTCorrelation

    function flengthCorrelation(reThetaTilde) result(Flength)
        !
        !       Smooth Flength correlation (Eqs. 49-50).
        !       Flength1 = exp(-3e-2 * (ReThetaTilde - 460))
        !       Flength  = 44 - (44 - (0.5 - 3e-4*(ReThetaTilde-596))) / (1+Flength1)^(1/6)
        !
        use constants, only: realType, one
        implicit none

        real(kind=realType), intent(in) :: reThetaTilde
        real(kind=realType) :: Flength

        real(kind=realType) :: Flength1, base

        ! Eq. 49
        Flength1 = exp(-3.0e-2_realType * (reThetaTilde - 460.0_realType))

        ! Eq. 50
        base = one + Flength1
        Flength = 44.0_realType - (44.0_realType &
                  - (0.5_realType - 3.0e-4_realType * (reThetaTilde - 596.0_realType))) &
                  / base**(one / 6.0_realType)

    end function flengthCorrelation

    function rethetacCorrelation(reThetaTilde) result(reThetaC)
        !
        !       Smooth Reθc correlation (Eq. 51).
        !       Reθc = 0.67*ReThetaTilde + 24*sin(ReThetaTilde/240 + 0.5) + 14
        !
        use constants, only: realType
        implicit none

        real(kind=realType), intent(in) :: reThetaTilde
        real(kind=realType) :: reThetaC

        reThetaC = 0.67_realType * reThetaTilde &
                 + 24.0_realType * sin(reThetaTilde / 240.0_realType + 0.5_realType) &
                 + 14.0_realType

    end function rethetacCorrelation

    function smoothMinMax(g1, g2, p) result(phi)
        !
        !       Smooth approximation of min (p<0) or max (p>0).
        !
        !       phi_p(g1,g2) = g1 + (1/p)*ln(1 + exp(p*(g2-g1)))
        !
        !       p = +300 : smooth max
        !       p = -300 : smooth min
        !
        !       Overflow-safe: uses asymptotic limits when the
        !       exponent is large.
        !
        use constants, only: realType, one
        implicit none

        real(kind=realType), intent(in) :: g1, g2, p
        real(kind=realType) :: phi

        real(kind=realType) :: arg

        arg = p * (g2 - g1)

        if (arg > 500.0_realType) then
            ! exp(arg) >> 1: ln(1+exp(arg)) ~ arg
            phi = g2
        else if (arg < -500.0_realType) then
            ! exp(arg) ~ 0: ln(1+exp(arg)) ~ 0
            phi = g1
        else
            phi = g1 + log(one + exp(arg)) / p
        end if

    end function smoothMinMax

end module turbUtils
