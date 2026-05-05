










! This module contains the source code related to the SA turbulence
! model. It is slightly more modularized than the original which makes
! performing reverse mode AD simplier.

module sa

    use constants
    real(kind=realType) :: cv13, kar2Inv, cw36, cb3Inv
    real(kind=realType), dimension(:, :, :), allocatable :: qq
    real(kind=realType), dimension(:, :, :), pointer :: ddw, ww, ddvt
    real(kind=realType), dimension(:, :), pointer :: rrlv
    real(kind=realType), dimension(:, :), pointer :: dd2Wall

contains

    subroutine saSource
        !
        !  Source terms.
        !  Determine the source term and its derivative w.r.t. nuTilde
        !  for all internal cells of the block.
        !  Remember that the SA field variable nuTilde = w(i,j,k,itu1)

        use blockPointers
        use constants
        use paramTurb
        use section
        use inputPhysics
        use inputDiscretization, only: approxSA
        use flowVarRefState
        implicit none

        ! Local parameters
        real(kind=realType), parameter :: f23 = two * third

        ! Local variables.
        integer(kind=intType) :: i, j, k, nn, ii
        real(kind=realType) :: fv1, fv2, ft2
        real(kind=realType) :: ss, sst, nu, dist2Inv, chi, chi2, chi3
        real(kind=realType) :: rr, gg, gg6, termFw, fwSa, term1, term2
        real(kind=realType) :: dfv1, dfv2, dft2, drr, dgg, dfw
        real(kind=realType) :: uux, uuy, uuz, vvx, vvy, vvz, wwx, wwy, wwz
        real(kind=realType) :: div2, fact, sxx, syy, szz, sxy, sxz, syz
        real(kind=realType) :: vortx, vorty, vortz
        real(kind=realType) :: omegax, omegay, omegaz
        real(kind=realType) :: strainMag2, strainProd, vortProd
        real(kind=realType), parameter :: xminn = 1.e-10_realType

        ! Set model constants
        cv13 = rsaCv1**3
        kar2Inv = one / (rsaK**2)
        cw36 = rsaCw3**6
        cb3Inv = one / rsaCb3

        ! Determine the non-dimensional wheel speed of this block.

        omegax = timeRef * sections(sectionID)%rotRate(1)
        omegay = timeRef * sections(sectionID)%rotRate(2)
        omegaz = timeRef * sections(sectionID)%rotRate(3)

        ! Create switches to production term depending on the variable that
        ! should be used
        if (turbProd .eq. katoLaunder) then
            print *, 'katoLaunder production term not supported for SA'
            stop
        end if

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

                        if (turbProd .eq. strain) then

                            sxx = two * fact * uux
                            syy = two * fact * vvy
                            szz = two * fact * wwz

                            sxy = fact * (uuy + vvx)
                            sxz = fact * (uuz + wwx)
                            syz = fact * (vvz + wwy)

                            ! Compute 2/3 * divergence of velocity squared

                            div2 = f23 * (sxx + syy + szz)**2

                            ! Compute strain production term

                            strainMag2 = two * (sxy**2 + sxz**2 + syz**2) &
                                         + sxx**2 + syy**2 + szz**2

                            strainProd = two * strainMag2 - div2

                            ss = sqrt(strainProd)

                        else if (turbProd .eq. vorticity) then

                            ! Compute the three components of the vorticity vector.
                            ! Substract the part coming from the rotating frame.

                            vortx = two * fact * (wwy - vvz) - two * omegax
                            vorty = two * fact * (uuz - wwx) - two * omegay
                            vortz = two * fact * (vvx - uuy) - two * omegaz

                            ! Compute the vorticity production term

                            vortProd = vortx**2 + vorty**2 + vortz**2

                            ! First take the square root of the production term to
                            ! obtain the correct production term for spalart-allmaras.
                            ! We do this to avoid if statements.

                            ss = sqrt(vortProd)

                        end if

                        ! Compute the laminar kinematic viscosity, the inverse of
                        ! wall distance squared, the ratio chi (ratio of nuTilde
                        ! and nu) and the functions fv1 and fv2. The latter corrects
                        ! the production term near a viscous wall.

                        nu = rlv(i, j, k) / w(i, j, k, irho)
                        dist2Inv = one / (d2Wall(i, j, k)**2)
                        chi = w(i, j, k, itu1) / nu
                        chi2 = chi * chi
                        chi3 = chi * chi2
                        fv1 = chi3 / (chi3 + cv13)
                        fv2 = one - chi / (one + chi * fv1)

                        ! The function ft2, which is designed to keep a laminar
                        ! solution laminar. When running in fully turbulent mode
                        ! this function should be set to 0.0.

                        if (useft2SA) then
                            ft2 = rsaCt3 * exp(-rsaCt4 * chi2)
                        else
                            ft2 = zero
                        end if

                        ! Correct the production term to account for the influence
                        ! of the wall.

                        sst = ss + w(i, j, k, itu1) * fv2 * kar2Inv * dist2Inv

                        ! Add rotation term (useRotationSA defined in inputParams.F90)

                        if (useRotationSA) then
                            sst = sst + rsaCrot * min(zero, sqrt(two * strainMag2))
                        end if

                        ! Make sure that this term remains positive
                        ! (the function fv2 is negative between chi = 1 and 18.4,
                        ! which can cause sst to go negative, which is undesirable).

                        sst = max(sst, xminn)

                        ! Compute the function fw. The argument rr is cut off at 10
                        ! to avoid numerical problems. This is ok, because the
                        ! asymptotical value of fw is then already reached.

                        rr = w(i, j, k, itu1) * kar2Inv * dist2Inv / sst
                        rr = min(rr, 10.0_realType)
                        gg = rr + rsaCw2 * (rr**6 - rr)
                        gg6 = gg**6
                        termFw = ((one + cw36) / (gg6 + cw36))**sixth
                        fwSa = gg * termFw

                        ! Compute the source term; some terms are saved for the
                        ! linearization. The source term is stored in dvt.

                        if (approxSA) then
                            term1 = zero
                        else
                            term1 = rsaCb1 * (one - ft2) * ss
                        end if
                        term2 = dist2Inv * (kar2Inv * rsaCb1 * ((one - ft2) * fv2 + ft2) &
                                            - rsaCw1 * fwSa)

                        scratch(i, j, k, idvt) = (term1 + term2 * w(i, j, k, itu1)) * w(i, j, k, itu1)

                end do
            end do
        end do
    end subroutine saSource

    subroutine saViscous
        !
        !  Viscous term.
        !  Determine the viscous contribution to the residual
        !  for all internal cells of the block.

        use blockPointers
        use paramTurb
        implicit none
        ! Local variables.
        integer(kind=intType) :: i, j, k, nn, ii
        real(kind=realType) :: nu
        real(kind=realType) :: fv1, fv2, ft2
        real(kind=realType) :: voli, volmi, volpi, xm, ym, zm, xp, yp, zp
        real(kind=realType) :: xa, ya, za, ttm, ttp, cnud, cam, cap
        real(kind=realType) :: nutm, nutp, num, nup, cdm, cdp
        real(kind=realType) :: c1m, c1p, c10, b1, c1, d1, qs

        ! Set model constants
        cv13 = rsaCv1**3
        kar2Inv = one / (rsaK**2)
        cw36 = rsaCw3**6
        cb3Inv = one / rsaCb3

        !
        !       Viscous terms in k-direction.
        !
            do k = 2, kl
                do j = 2, jl
                    do i = 2, il
                        ! Compute the metrics in zeta-direction, i.e. along the
                        ! line k = constant.

                        voli = one / vol(i, j, k)
                        volmi = two / (vol(i, j, k) + vol(i, j, k - 1))
                        volpi = two / (vol(i, j, k) + vol(i, j, k + 1))

                        xm = sk(i, j, k - 1, 1) * volmi
                        ym = sk(i, j, k - 1, 2) * volmi
                        zm = sk(i, j, k - 1, 3) * volmi
                        xp = sk(i, j, k, 1) * volpi
                        yp = sk(i, j, k, 2) * volpi
                        zp = sk(i, j, k, 3) * volpi

                        xa = half * (sk(i, j, k, 1) + sk(i, j, k - 1, 1)) * voli
                        ya = half * (sk(i, j, k, 2) + sk(i, j, k - 1, 2)) * voli
                        za = half * (sk(i, j, k, 3) + sk(i, j, k - 1, 3)) * voli
                        ttm = xm * xa + ym * ya + zm * za
                        ttp = xp * xa + yp * ya + zp * za

                        ! ttm and ttp ~ 1/deltaX^2

                        ! Computation of the viscous terms in zeta-direction; note
                        ! that cross-derivatives are neglected, i.e. the mesh is
                        ! assumed to be orthogonal.
                        ! Furthermore, the grad(nu)**2 has been rewritten as
                        ! div(nu grad(nu)) - nu div(grad nu) to enhance stability.
                        ! The second derivative in zeta-direction is constructed as
                        ! the central difference of the first order derivatives, i.e.
                        ! d^2/dzeta^2 = d/dzeta (d/dzeta k+1/2 - d/dzeta k-1/2).
                        ! In this way the metric can be taken into account.

                        ! Compute the diffusion coefficients multiplying the nodes
                        ! k+1, k and k-1 in the second derivative. Make sure that
                        ! these coefficients are nonnegative.

                        cnud = -rsaCb2 * w(i, j, k, itu1) * cb3Inv
                        cam = ttm * cnud
                        cap = ttp * cnud

                        ! Compute nuTilde at the faces

                        nutm = half * (w(i, j, k - 1, itu1) + w(i, j, k, itu1))
                        nutp = half * (w(i, j, k + 1, itu1) + w(i, j, k, itu1))

                        ! Compute nu at the faces

                        nu = rlv(i, j, k) / w(i, j, k, irho)
                        num = half * (rlv(i, j, k - 1) / w(i, j, k - 1, irho) + nu)
                        nup = half * (rlv(i, j, k + 1) / w(i, j, k + 1, irho) + nu)

                        cdm = (num + (one + rsaCb2) * nutm) * ttm * cb3Inv
                        cdp = (nup + (one + rsaCb2) * nutp) * ttp * cb3Inv

                        c1m = max(cdm + cam, zero)
                        c1p = max(cdp + cap, zero)
                        c10 = c1m + c1p

                        ! Update the residual for this cell and store the possible
                        ! coefficients for the matrix in b1, c1 and d1.

                        scratch(i, j, k, idvt) = scratch(i, j, k, idvt) + c1m * w(i, j, k - 1, itu1) &
                                                 - c10 * w(i, j, k, itu1) + c1p * w(i, j, k + 1, itu1)
                end do
            end do
        end do
        !
        !       Viscous terms in j-direction.
        !
            do k = 2, kl
                do j = 2, jl
                    do i = 2, il
                        ! Compute the metrics in eta-direction, i.e. along the
                        ! line j = constant.

                        voli = one / vol(i, j, k)
                        volmi = two / (vol(i, j, k) + vol(i, j - 1, k))
                        volpi = two / (vol(i, j, k) + vol(i, j + 1, k))

                        xm = sj(i, j - 1, k, 1) * volmi
                        ym = sj(i, j - 1, k, 2) * volmi
                        zm = sj(i, j - 1, k, 3) * volmi
                        xp = sj(i, j, k, 1) * volpi
                        yp = sj(i, j, k, 2) * volpi
                        zp = sj(i, j, k, 3) * volpi

                        xa = half * (sj(i, j, k, 1) + sj(i, j - 1, k, 1)) * voli
                        ya = half * (sj(i, j, k, 2) + sj(i, j - 1, k, 2)) * voli
                        za = half * (sj(i, j, k, 3) + sj(i, j - 1, k, 3)) * voli
                        ttm = xm * xa + ym * ya + zm * za
                        ttp = xp * xa + yp * ya + zp * za

                        ! Computation of the viscous terms in eta-direction; note
                        ! that cross-derivatives are neglected, i.e. the mesh is
                        ! assumed to be orthogonal.
                        ! Furthermore, the grad(nu)**2 has been rewritten as
                        ! div(nu grad(nu)) - nu div(grad nu) to enhance stability.
                        ! The second derivative in eta-direction is constructed as
                        ! the central difference of the first order derivatives, i.e.
                        ! d^2/deta^2 = d/deta (d/deta j+1/2 - d/deta j-1/2).
                        ! In this way the metric can be taken into account.

                        ! Compute the diffusion coefficients multiplying the nodes
                        ! j+1, j and j-1 in the second derivative. Make sure that
                        ! these coefficients are nonnegative.

                        cnud = -rsaCb2 * w(i, j, k, itu1) * cb3Inv
                        cam = ttm * cnud
                        cap = ttp * cnud

                        nutm = half * (w(i, j - 1, k, itu1) + w(i, j, k, itu1))
                        nutp = half * (w(i, j + 1, k, itu1) + w(i, j, k, itu1))
                        nu = rlv(i, j, k) / w(i, j, k, irho)
                        num = half * (rlv(i, j - 1, k) / w(i, j - 1, k, irho) + nu)
                        nup = half * (rlv(i, j + 1, k) / w(i, j + 1, k, irho) + nu)
                        cdm = (num + (one + rsaCb2) * nutm) * ttm * cb3Inv
                        cdp = (nup + (one + rsaCb2) * nutp) * ttp * cb3Inv

                        c1m = max(cdm + cam, zero)
                        c1p = max(cdp + cap, zero)
                        c10 = c1m + c1p

                        ! Update the residual for this cell and store the possible
                        ! coefficients for the matrix in b1, c1 and d1.

                        scratch(i, j, k, idvt) = scratch(i, j, k, idvt) + c1m * w(i, j - 1, k, itu1) &
                                                 - c10 * w(i, j, k, itu1) + c1p * w(i, j + 1, k, itu1)
                end do
            end do
        end do
        !
        !       Viscous terms in i-direction.
        !
            do k = 2, kl
                do j = 2, jl
                    do i = 2, il
                        ! Compute the metrics in xi-direction, i.e. along the
                        ! line i = constant.

                        voli = one / vol(i, j, k)
                        volmi = two / (vol(i, j, k) + vol(i - 1, j, k))
                        volpi = two / (vol(i, j, k) + vol(i + 1, j, k))

                        xm = si(i - 1, j, k, 1) * volmi
                        ym = si(i - 1, j, k, 2) * volmi
                        zm = si(i - 1, j, k, 3) * volmi
                        xp = si(i, j, k, 1) * volpi
                        yp = si(i, j, k, 2) * volpi
                        zp = si(i, j, k, 3) * volpi

                        xa = half * (si(i, j, k, 1) + si(i - 1, j, k, 1)) * voli
                        ya = half * (si(i, j, k, 2) + si(i - 1, j, k, 2)) * voli
                        za = half * (si(i, j, k, 3) + si(i - 1, j, k, 3)) * voli
                        ttm = xm * xa + ym * ya + zm * za
                        ttp = xp * xa + yp * ya + zp * za

                        ! Computation of the viscous terms in xi-direction; note
                        ! that cross-derivatives are neglected, i.e. the mesh is
                        ! assumed to be orthogonal.
                        ! Furthermore, the grad(nu)**2 has been rewritten as
                        ! div(nu grad(nu)) - nu div(grad nu) to enhance stability.
                        ! The second derivative in xi-direction is constructed as
                        ! the central difference of the first order derivatives, i.e.
                        ! d^2/dxi^2 = d/dxi (d/dxi i+1/2 - d/dxi i-1/2).
                        ! In this way the metric can be taken into account.

                        ! Compute the diffusion coefficients multiplying the nodes
                        ! i+1, i and i-1 in the second derivative. Make sure that
                        ! these coefficients are nonnegative.

                        cnud = -rsaCb2 * w(i, j, k, itu1) * cb3Inv
                        cam = ttm * cnud
                        cap = ttp * cnud

                        nutm = half * (w(i - 1, j, k, itu1) + w(i, j, k, itu1))
                        nutp = half * (w(i + 1, j, k, itu1) + w(i, j, k, itu1))
                        nu = rlv(i, j, k) / w(i, j, k, irho)
                        num = half * (rlv(i - 1, j, k) / w(i - 1, j, k, irho) + nu)
                        nup = half * (rlv(i + 1, j, k) / w(i + 1, j, k, irho) + nu)
                        cdm = (num + (one + rsaCb2) * nutm) * ttm * cb3Inv
                        cdp = (nup + (one + rsaCb2) * nutp) * ttp * cb3Inv

                        c1m = max(cdm + cam, zero)
                        c1p = max(cdp + cap, zero)
                        c10 = c1m + c1p

                        ! Update the residual for this cell and store the possible
                        ! coefficients for the matrix in b1, c1 and d1.

                        scratch(i, j, k, idvt) = scratch(i, j, k, idvt) + c1m * w(i - 1, j, k, itu1) &
                                                 - c10 * w(i, j, k, itu1) + c1p * w(i + 1, j, k, itu1)
                end do
            end do
        end do
    end subroutine saViscous

    subroutine saResScale

        !
        !  Multiply the residual by the volume and store this in dw; this
        ! * is done for monitoring reasons only. The multiplication with the
        ! * volume is present to be consistent with the flow residuals; also
        !  the negative value is taken, again to be consistent with the
        ! * flow equations. Also multiply by iblank so that no updates occur
        !  in holes or the overset boundary.
        use blockPointers
        implicit none

        ! Local variables
        integer(kind=intType) :: i, j, k, ii
        real(kind=realType) :: rblank

            do k = 2, kl
                do j = 2, jl
                    do i = 2, il
                        rblank = max(real(iblank(i, j, k), realType), zero)
                        dw(i, j, k, itu1) = -volRef(i, j, k) * scratch(i, j, k, idvt) * rblank
                end do
            end do
        end do
    end subroutine saResScale

end module sa
