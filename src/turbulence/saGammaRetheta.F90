module saGammaReTheta

    ! This module contains the source code related to the SST turbulence
    ! model. It is slightly more modularized than the original which makes
    ! performing reverse mode AD simplier.

    use constants, only: realType

    real(kind=realType), dimension(:, :, :, :, :), allocatable :: qq

contains

    subroutine saGammaReTheta_block(resOnly)

        use constants
        use blockPointers, only: il, jl, kl
        use inputTimeSpectral
        use iteration
        use turbUtils, only: SSTEddyViscosity, turbAdvection, unsteadyTurbTerm, saEddyViscosity
        use turbBCRoutines, only: bcTurbTreatment, applyAllTurbBCThisBlock
        implicit none

        !
        !      Subroutine argument.
        !
        logical, intent(in) :: resOnly
        !
        !      Local variables.
        !
        integer(kind=intType) :: nn, sps

        ! Set the arrays for the boundary condition treatment.

        call bcTurbTreatment

        ! Alloc central jacobian memory
        allocate (qq(2:il,2:jl,2:kl,3,3))

        ! Source Terms
        call Source

        ! Advection Term
        nn = itu1 - 1
        call turbAdvection(3_intType, 3_intType, nn, qq)

        call unsteadyTurbTerm(3_intType, 3_intType, nn, qq)

        ! Viscous Terms
        call Viscous

        ! Perform the residual scaling
        call ResScale



        ! The eddy viscosity and the boundary conditions are only
        ! applied if an actual update has been computed in saGammaReThetaSolve.

        if (.not. resOnly) then

            ! Solve the transport equations for v, gamma, Retheta.

            call saGammaReThetaSolve(resOnly)

            ! Compute the corresponding eddy viscosity.

            call saEddyViscosity(2, il, 2, jl, 2, kl)

            ! Set the halo values for the turbulent variables.
            ! We are on the finest mesh, so the second layer of halo
            ! cells must be computed as well.

            call applyAllTurbBCThisBlock(.true.)
        end if
        deallocate (qq)

    end subroutine saGammaReTheta_block



    
    subroutine Source
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
        use turbMod, only: dvt, vort, prod, kwCD, f1
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
        real(kind=realType) :: cv13, kar2Inv, cw36, cb3Inv
        real(kind=realType), parameter :: gammaStatic = one
        real(kind=realType), parameter :: reThetaStatic = one
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

#ifdef TAPENADE_REVERSE
        !$AD II-LOOP
        do ii = 0, nx * ny * nz - 1
            i = mod(ii, nx) + 2
            j = mod(ii / nx, ny) + 2
            k = ii / (nx * ny) + 2
#else
            do k = 2, kl
                do j = 2, jl
                    do i = 2, il
#endif
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
                        ! linearization. The source term is stored in scratch.

                        if (approxSA) then
                            term1 = zero
                        else
                            term1 = rsaCb1 * (one - ft2) * ss
                        end if
                        term2 = dist2Inv * (kar2Inv * rsaCb1 * ((one - ft2) * fv2 + ft2) &
                                            - rsaCw1 * fwSa)

                        scratch(i, j, k, idvt) = (term1 + term2 * w(i, j, k, itu1)) * w(i, j, k, itu1)

                        ! Placeholder equations for SA-gamma-ReTheta integration.
                        ! Gamma and ReTheta are kept frozen at static values and
                        ! contribute zero to the source residual for now.

                        scratch(i, j, k, idvt+1) = scratch(i, j, k, idvt+1) + zero * (w(i, j, k, itu2) - gammaStatic)
                        scratch(i, j, k, idvt+2) = scratch(i, j, k, idvt+2) + zero * (w(i, j, k, itu3) - reThetaStatic)

#ifndef USE_TAPENADE
                        ! Compute some derivatives w.r.t. nuTilde. These will occur
                        ! in the left hand side, i.e. the matrix for the implicit
                        ! treatment.

                        dfv1 = three * chi2 * cv13 / ((chi3 + cv13)**2)
                        dfv2 = (chi2 * dfv1 - one) / (nu * ((one + chi * fv1)**2))
                        dft2 = -two * rsaCt4 * chi * ft2 / nu

                        drr = (one - rr * (fv2 + w(i, j, k, itu1) * dfv2)) &
                              * kar2Inv * dist2Inv / sst
                        dgg = (one - rsaCw2 + six * rsaCw2 * (rr**5)) * drr
                        dfw = (cw36 / (gg6 + cw36)) * termFw * dgg

                        ! Compute the source term jacobian. Note that the part
                        ! containing term1 is treated explicitly. The reason is that
                        ! implicit treatment of this part leads to a decrease of the
                        ! diagonal dominance of the jacobian and it thus decreases
                        ! the stability. You may want to play around and try to
                        ! take this term into account in the jacobian.
                        ! Note that -dsource/dnu is stored.
                        qq(i, j, k, 1, 1) = -two * term2 * w(i, j, k, itu1) &
                                      - dist2Inv * w(i, j, k, itu1) * w(i, j, k, itu1) &
                                      * (rsaCb1 * kar2Inv * (dfv2 - ft2 * dfv2 - fv2 * dft2 + dft2) &
                                         - rsaCw1 * dfw)

                        ! A couple of terms in qq may lead to a negative
                        ! contribution. Clip qq to zero, if the total is negative.

                        qq(i, j, k, 1, 1) = max(qq(i, j, k, 1, 1), zero)
                        qq(i, j, k, 2, 2) = zero
                        qq(i, j, k, 3, 3) = zero
#endif
#ifdef TAPENADE_REVERSE
                    end do
#else
                end do
            end do
        end do
#endif
    end subroutine Source

    subroutine Viscous
        !
        !  Viscous term for SA variable.
        !  Determine the viscous contribution to the residual
        !  for all internal cells of the block.

        use constants
        use blockPointers
        use paramTurb
        implicit none

        integer(kind=intType) :: i, j, k, ii

        ! viscosity variables
        real(kind=realType) :: nu, nu_m, nu_p
        real(kind=realType) :: nut, nut_m, nut_p
        real(kind=realType) :: nu_tm, nu_tp
        real(kind=realType) :: nuTilde, nuTilde_m, nuTilde_p

        ! SA auxiliary functions
        real(kind=realType) :: chi, chi3
        real(kind=realType) :: chi_m, chi3_m
        real(kind=realType) :: chi_p, chi3_p
        real(kind=realType) :: fv1, fv1_m, fv1_p

        ! geometry terms
        real(kind=realType) :: voli, volmi, volpi
        real(kind=realType) :: xm, ym, zm
        real(kind=realType) :: xp, yp, zp
        real(kind=realType) :: xa, ya, za
        real(kind=realType) :: ttm, ttp

        ! diffusion coefficients
        real(kind=realType) :: num, nup
        real(kind=realType) :: cdm, cdp
        real(kind=realType) :: cdm_gamma, cdp_gamma
        real(kind=realType) :: cdm_rt, cdp_rt

        ! SA nonlinear correction
        real(kind=realType) :: cnud, cam, cap
        real(kind=realType) :: nutm, nutp

        ! matrix coefficients
        real(kind=realType) :: c1m, c1p, c10
        real(kind=realType) :: c2m, c2p, c20
        real(kind=realType) :: c3m, c3p, c30
        real(kind=realType) :: b1, c1, d1
        real(kind=realType) :: b2, c2, d2
        real(kind=realType) :: b3, c3, d3

        ! constants
        real(kind=realType) :: cb3Inv, cv13

        

        cb3Inv = one / rsaCb3
        cv13 = rsaCv1**3



        ! Viscous terms in k-direction.
#ifdef TAPENADE_REVERSE
        !$AD II-LOOP
        do ii = 0, nx * ny * nz - 1
            i = mod(ii, nx) + 2
            j = mod(ii / nx, ny) + 2
            k = ii / (nx * ny) + 2
#else
            do k = 2, kl
                do j = 2, jl
                    do i = 2, il
#endif
                        !Mesh Geometry contribution
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

                        cnud = -rsaCb2 * w(i, j, k, itu1) * cb3Inv
                        cam = ttm * cnud
                        cap = ttp * cnud

                        nutm = half * (w(i, j, k - 1, itu1) + w(i, j, k, itu1))
                        nutp = half * (w(i, j, k + 1, itu1) + w(i, j, k, itu1))


                        ! ----- viscosities in cell (i,j,k)
                        nu  = rlv(i,j,k) / w(i,j,k,irho)
                        nuTilde = w(i,j,k,itu1)

                        chi  = nuTilde / nu
                        chi3 = chi*chi*chi
                        fv1  = chi3 / (chi3 + cv13)

                        nut = nuTilde * fv1


                        ! ----- k-1 cell
                        nu_m = rlv(i,j,k-1) / w(i,j,k-1,irho)
                        nuTilde_m = w(i,j,k-1,itu1)

                        chi_m  = nuTilde_m / nu_m
                        chi3_m = chi_m*chi_m*chi_m
                        fv1_m  = chi3_m / (chi3_m + cv13)

                        nut_m = nuTilde_m * fv1_m


                        ! ----- k+1 cell
                        nu_p = rlv(i,j,k+1) / w(i,j,k+1,irho)
                        nuTilde_p = w(i,j,k+1,itu1)

                        chi_p  = nuTilde_p / nu_p
                        chi3_p = chi_p*chi_p*chi_p
                        fv1_p  = chi3_p / (chi3_p + cv13)

                        nut_p = nuTilde_p * fv1_p


                        num = half * (nu_m + nu)
                        nup = half * (nu_p + nu)

                        nu_tm = half * (nut_m + nut)
                        nu_tp = half * (nut_p + nut)

                        !SA diffusion contribution
                        cdm = (num + (one + rsaCb2) * nutm) * ttm * cb3Inv
                        cdp = (nup + (one + rsaCb2) * nutp) * ttp * cb3Inv

                        !Gamma diffusion contribution
                        cdm_gamma = (num + nu_tm/sigmaF) * ttm
                        cdp_gamma = (nup + nu_tp/sigmaF) * ttp

                        !Rethata diffusion contribution
                        cdm_rt = sigmaTheta * (num + nu_tm) * ttm
                        cdp_rt = sigmaTheta * (nup + nu_tp) * ttp


                        c1m = max(cdm + cam, zero)
                        c1p = max(cdp + cap, zero)
                        c10 = c1m + c1p

                        c2m = max(cdm_gamma, zero)
                        c2p = max(cdp_gamma, zero)
                        c20 = c2m + c2p

                        c3m = max(cdm_rt, zero)
                        c3p = max(cdp_rt, zero)
                        c30 = c3m + c3p

                        scratch(i, j, k, idvt) = scratch(i, j, k, idvt) + c1m * w(i, j, k - 1, itu1) &
                                                 - c10 * w(i, j, k, itu1) + c1p * w(i, j, k + 1, itu1)
                        scratch(i, j, k, idvt + 1) = scratch(i, j, k, idvt + 1) + c2m * w(i, j, k - 1, itu2) &
                                                     - c20 * w(i, j, k, itu2) + c2p * w(i, j, k + 1, itu2)
                        scratch(i, j, k, idvt + 2) = scratch(i, j, k, idvt + 2) + c3m * w(i, j, k - 1, itu3) &
                                                     - c30 * w(i, j, k, itu3) + c3p * w(i, j, k + 1, itu3)
#ifndef USE_TAPENADE
                        b1 = -c1m
                        c1 = c10
                        d1 = -c1p

                        ! Update the central jacobian. For nonboundary cells this
                        ! is simply c1. For boundary cells this is slightly more
                        ! complicated, because the boundary conditions are treated
                        ! implicitly and the off-diagonal terms b1 and d1 must be
                        ! taken into account.
                        ! The boundary conditions are only treated implicitly if
                        ! the diagonal dominance of the matrix is increased.

                        if (k == 2) then
                            qq(i,j,k,1,1) = qq(i,j,k,1,1) + c1 &
                                            - b1 * max(bmtk1(i,j,itu1,itu1), zero)

                        else if (k == kl) then
                            qq(i,j,k,1,1) = qq(i,j,k,1,1) + c1 &
                                            - d1 * max(bmtk2(i,j,itu1,itu1), zero)

                        else
                            qq(i,j,k,1,1) = qq(i,j,k,1,1) + c1
                        end if

                        !====================================================
                        ! GAMMA EQUATION
                        !====================================================

                        b2 = -c2m
                        c2 =  c20
                        d2 = -c2p

                        if (k == 2) then
                            qq(i,j,k,2,2) = qq(i,j,k,2,2) + c2 &
                                            - b2 * max(bmtk1(i,j,itu2,itu2), zero)

                        else if (k == kl) then
                            qq(i,j,k,2,2) = qq(i,j,k,2,2) + c2 &
                                            - d2 * max(bmtk2(i,j,itu2,itu2), zero)

                        else
                            qq(i,j,k,2,2) = qq(i,j,k,2,2) + c2
                        end if


                        !====================================================
                        ! RETHETA EQUATION
                        !====================================================

                        b3 = -c3m
                        c3 =  c30
                        d3 = -c3p

                        if (k == 2) then
                            qq(i,j,k,3,3) = qq(i,j,k,3,3) + c3 &
                                            - b3 * max(bmtk1(i,j,itu3,itu3), zero)

                        else if (k == kl) then
                            qq(i,j,k,3,3) = qq(i,j,k,3,3) + c3 &
                                            - d3 * max(bmtk2(i,j,itu3,itu3), zero)

                        else
                            qq(i,j,k,3,3) = qq(i,j,k,3,3) + c3
                        end if
#endif
#ifdef TAPENADE_REVERSE
                    end do
#else
                end do
            end do
        end do
#endif
        ! Viscous terms in j-direction.
#ifdef TAPENADE_REVERSE
        !$AD II-LOOP
        do ii = 0, nx * ny * nz - 1
            i = mod(ii, nx) + 2
            j = mod(ii / nx, ny) + 2
            k = ii / (nx * ny) + 2
#else
            do k = 2, kl
                do j = 2, jl
                    do i = 2, il
#endif
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

                        cnud = -rsaCb2 * w(i, j, k, itu1) * cb3Inv
                        cam = ttm * cnud
                        cap = ttp * cnud

                        nutm = half * (w(i, j - 1, k, itu1) + w(i, j, k, itu1))
                        nutp = half * (w(i, j + 1, k, itu1) + w(i, j, k, itu1))

                       ! ----- viscosities in cell (i,j,k)
                        nu  = rlv(i,j,k) / w(i,j,k,irho)
                        nuTilde = w(i,j,k,itu1)

                        chi  = nuTilde / nu
                        chi3 = chi*chi*chi
                        fv1  = chi3 / (chi3 + cv13)

                        nut = nuTilde * fv1


                        ! ----- j-1 cell
                        nu_m = rlv(i,j-1,k) / w(i,j-1,k,irho)
                        nuTilde_m = w(i,j-1,k,itu1)

                        chi_m  = nuTilde_m / nu_m
                        chi3_m = chi_m*chi_m*chi_m
                        fv1_m  = chi3_m / (chi3_m + cv13)

                        nut_m = nuTilde_m * fv1_m


                        ! ----- j+1 cell
                        nu_p = rlv(i,j+1,k) / w(i,j+1,k,irho)
                        nuTilde_p = w(i,j+1,k,itu1)

                        chi_p  = nuTilde_p / nu_p
                        chi3_p = chi_p*chi_p*chi_p
                        fv1_p  = chi3_p / (chi3_p + cv13)

                        nut_p = nuTilde_p * fv1_p


                        num = half * (nu_m + nu)
                        nup = half * (nu_p + nu)

                        nu_tm = half * (nut_m + nut)
                        nu_tp = half * (nut_p + nut)

                        !SA diffusion contribution
                        cdm = (num + (one + rsaCb2) * nutm) * ttm * cb3Inv
                        cdp = (nup + (one + rsaCb2) * nutp) * ttp * cb3Inv

                        !Gamma diffusion contribution
                        cdm_gamma = (num + nu_tm/sigmaF) * ttm
                        cdp_gamma = (nup + nu_tp/sigmaF) * ttp

                        !Rethata diffusion contribution
                        cdm_rt = sigmaTheta * (num + nu_tm) * ttm
                        cdp_rt = sigmaTheta * (nup + nu_tp) * ttp


                        c1m = max(cdm + cam, zero)
                        c1p = max(cdp + cap, zero)
                        c10 = c1m + c1p

                        c2m = max(cdm_gamma, zero)
                        c2p = max(cdp_gamma, zero)
                        c20 = c2m + c2p

                        c3m = max(cdm_rt, zero)
                        c3p = max(cdp_rt, zero)
                        c30 = c3m + c3p

                        scratch(i, j, k, idvt) = scratch(i, j, k, idvt) + c1m * w(i, j-1, k, itu1) &
                                                 - c10 * w(i, j, k, itu1) + c1p * w(i, j+1, k, itu1)
                        scratch(i, j, k, idvt + 1) = scratch(i, j, k, idvt + 1) + c2m * w(i, j-1, k, itu2) &
                                                     - c20 * w(i, j, k, itu2) + c2p * w(i, j+1, k, itu2)
                        scratch(i, j, k, idvt + 2) = scratch(i, j, k, idvt + 2) + c3m * w(i, j-1, k, itu3) &
                                                     - c30 * w(i, j, k, itu3) + c3p * w(i, j+1, k, itu3)
#ifdef TAPENADE_REVERSE
                        b1 = -c1m
                        c1 = c10
                        d1 = -c1p

                        ! Update the central jacobian. For nonboundary cells this
                        ! is simply c1. For boundary cells this is slightly more
                        ! complicated, because the boundary conditions are treated
                        ! implicitly and the off-diagonal terms b1 and d1 must be
                        ! taken into account.
                        ! The boundary conditions are only treated implicitly if
                        ! the diagonal dominance of the matrix is increased.

                        if (j == 2) then
                            qq(i,j,k,1,1) = qq(i,j,k,1,1) + c1 &
                                            - b1 * max(bmtj1(i,j,itu1,itu1), zero)

                        else if (j == jl) then
                            qq(i,j,k,1,1) = qq(i,j,k,1,1) + c1 &
                                            - d1 * max(bmtj2(i,j,itu1,itu1), zero)

                        else
                            qq(i,j,k,1,1) = qq(i,j,k,1,1) + c1
                        end if

                        !====================================================
                        ! GAMMA EQUATION
                        !====================================================

                        b2 = -c2m
                        c2 =  c20
                        d2 = -c2p

                        if (j == 2) then
                            qq(i,j,k,2,2) = qq(i,j,k,2,2) + c2 &
                                            - b2 * max(bmtj1(i,j,itu2,itu2), zero)

                        else if (j == jl) then
                            qq(i,j,k,2,2) = qq(i,j,k,2,2) + c2 &
                                            - d2 * max(bmtj2(i,j,itu2,itu2), zero)

                        else
                            qq(i,j,k,2,2) = qq(i,j,k,2,2) + c2
                        end if


                        !====================================================
                        ! RETHETA EQUATION
                        !====================================================

                        b3 = -c3m
                        c3 =  c30
                        d3 = -c3p

                        if (j == 2) then
                            qq(i,j,k,3,3) = qq(i,j,k,3,3) + c3 &
                                            - b3 * max(bmtj1(i,j,itu3,itu3), zero)

                        else if (j == jl) then
                            qq(i,j,k,3,3) = qq(i,j,k,3,3) + c3 &
                                            - d3 * max(bmtj2(i,j,itu3,itu3), zero)

                        else
                            qq(i,j,k,3,3) = qq(i,j,k,3,3) + c3
                        end if
#endif
#ifdef TAPENADE_REVERSE
                    end do
#else
                end do
            end do
        end do
#endif
        ! Viscous terms in i-direction.
#ifdef TAPENADE_REVERSE
        !$AD II-LOOP
        do ii = 0, nx * ny * nz - 1
            i = mod(ii, nx) + 2
            j = mod(ii / nx, ny) + 2
            k = ii / (nx * ny) + 2
#else
            do k = 2, kl
                do j = 2, jl
                    do i = 2, il
#endif
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

                        cnud = -rsaCb2 * w(i, j, k, itu1) * cb3Inv
                        cam = ttm * cnud
                        cap = ttp * cnud

                        nutm = half * (w(i - 1, j, k, itu1) + w(i, j, k, itu1))
                        nutp = half * (w(i + 1, j, k, itu1) + w(i, j, k, itu1))
                        ! ----- viscosities in cell (i,j,k)
                        nu  = rlv(i,j,k) / w(i,j,k,irho)
                        nuTilde = w(i,j,k,itu1)

                        chi  = nuTilde / nu
                        chi3 = chi*chi*chi
                        fv1  = chi3 / (chi3 + cv13)

                        nut = nuTilde * fv1


                        ! ----- i-1 cell
                        nu_m = rlv(i-1,j,k) / w(i-1,j,k,irho)
                        nuTilde_m = w(i-1,j,k,itu1)

                        chi_m  = nuTilde_m / nu_m
                        chi3_m = chi_m*chi_m*chi_m
                        fv1_m  = chi3_m / (chi3_m + cv13)

                        nut_m = nuTilde_m * fv1_m


                        ! ----- i+1 cell
                        nu_p = rlv(i+1,j,k) / w(i+1,j,k,irho)
                        nuTilde_p = w(i+1,j,k,itu1)

                        chi_p  = nuTilde_p / nu_p
                        chi3_p = chi_p*chi_p*chi_p
                        fv1_p  = chi3_p / (chi3_p + cv13)

                        nut_p = nuTilde_p * fv1_p


                        num = half * (nu_m + nu)
                        nup = half * (nu_p + nu)

                        nu_tm = half * (nut_m + nut)
                        nu_tp = half * (nut_p + nut)

                        !SA diffusion contribution
                        cdm = (num + (one + rsaCb2) * nutm) * ttm * cb3Inv
                        cdp = (nup + (one + rsaCb2) * nutp) * ttp * cb3Inv

                        !Gamma diffusion contribution
                        cdm_gamma = (num + nu_tm/sigmaF) * ttm
                        cdp_gamma = (nup + nu_tp/sigmaF) * ttp

                        !Rethata diffusion contribution
                        cdm_rt = sigmaTheta * (num + nu_tm) * ttm
                        cdp_rt = sigmaTheta * (nup + nu_tp) * ttp


                        c1m = max(cdm + cam, zero)
                        c1p = max(cdp + cap, zero)
                        c10 = c1m + c1p

                        c2m = max(cdm_gamma, zero)
                        c2p = max(cdp_gamma, zero)
                        c20 = c2m + c2p

                        c3m = max(cdm_rt, zero)
                        c3p = max(cdp_rt, zero)
                        c30 = c3m + c3p

                        scratch(i, j, k, idvt) = scratch(i, j, k, idvt) + c1m * w(i-1, j, k, itu1) &
                                                 - c10 * w(i, j, k, itu1) + c1p * w(i+1, j, k, itu1)
                        scratch(i, j, k, idvt + 1) = scratch(i, j, k, idvt + 1) + c2m * w(i-1, j, k, itu2) &
                                                     - c20 * w(i, j, k, itu2) + c2p * w(i+1, j, k, itu2)
                        scratch(i, j, k, idvt + 2) = scratch(i, j, k, idvt + 2) + c3m * w(i-1, j, k, itu3) &
                                                     - c30 * w(i, j, k, itu3) + c3p * w(i+1, j, k, itu3)
#ifdef TAPENADE_REVERSE
                        b1 = -c1m
                        c1 = c10
                        d1 = -c1p

                        ! Update the central jacobian. For nonboundary cells this
                        ! is simply c1. For boundary cells this is slightly more
                        ! complicated, because the boundary conditions are treated
                        ! implicitly and the off-diagonal terms b1 and d1 must be
                        ! taken into account.
                        ! The boundary conditions are only treated implicitly if
                        ! the diagonal dominance of the matrix is increased.

                        if (i == 2) then
                            qq(i,j,k,1,1) = qq(i,j,k,1,1) + c1 &
                                            - b1 * max(bmti1(i,j,itu1,itu1), zero)

                        else if (i == il) then
                            qq(i,j,k,1,1) = qq(i,j,k,1,1) + c1 &
                                            - d1 * max(bmti2(i,j,itu1,itu1), zero)

                        else
                            qq(i,j,k,1,1) = qq(i,j,k,1,1) + c1
                        end if

                        !====================================================
                        ! GAMMA EQUATION
                        !====================================================

                        b2 = -c2m
                        c2 =  c20
                        d2 = -c2p

                        if (i == 2) then
                            qq(i,j,k,2,2) = qq(i,j,k,2,2) + c2 &
                                            - b2 * max(bmti1(i,j,itu2,itu2), zero)

                        else if (i == il) then
                            qq(i,j,k,2,2) = qq(i,j,k,2,2) + c2 &
                                            - d2 * max(bmti2(i,j,itu2,itu2), zero)

                        else
                            qq(i,j,k,2,2) = qq(i,j,k,2,2) + c2
                        end if


                        !====================================================
                        ! RETHETA EQUATION
                        !====================================================

                        b3 = -c3m
                        c3 =  c30
                        d3 = -c3p

                        if (i == 2) then
                            qq(i,j,k,3,3) = qq(i,j,k,3,3) + c3 &
                                            - b3 * max(bmti1(i,j,itu3,itu3), zero)

                        else if (i == il) then
                            qq(i,j,k,3,3) = qq(i,j,k,3,3) + c3 &
                                            - d3 * max(bmti2(i,j,itu3,itu3), zero)

                        else
                            qq(i,j,k,3,3) = qq(i,j,k,3,3) + c3
                        end if
#endif
#ifdef TAPENADE_REVERSE
                    end do
#else
                end do
            end do
        end do
#endif
    end subroutine Viscous

    subroutine ResScale

        !
        !  Multiply the residual by the volume and store this in dw; this
        ! * is done for monitoring reasons only. The multiplication with the
        ! * volume is present to be consistent with the flow residuals; also
        !  the negative value is taken, again to be consistent with the
        ! * flow equations. Also multiply by iblank so that no updates occur
        !  in holes or the overset boundary.
        use constants
        use blockPointers
        implicit none

        ! Local variables
        integer(kind=intType) :: i, j, k, ii
        real(kind=realType) :: rblank

#ifdef TAPENADE_REVERSE
        !$AD II-LOOP
        do ii = 0, nx * ny * nz - 1
            i = mod(ii, nx) + 2
            j = mod(ii / nx, ny) + 2
            k = ii / (nx * ny) + 2
#else
            do k = 2, kl
                do j = 2, jl
                    do i = 2, il
#endif
                        rblank = max(real(iblank(i, j, k), realType), zero)
                        ! SA variable ν~
                        dw(i,j,k,itu1) = -volRef(i,j,k) * scratch(i, j, k, idvt) * rblank

                        ! Gamma
                        dw(i,j,k,itu2) = -volRef(i,j,k) * scratch(i,j,k,idvt+1) * rblank

                        ! ReTheta_t
                        dw(i,j,k,itu3) = -volRef(i,j,k) * scratch(i,j,k,idvt+2) * rblank
#ifdef TAPENADE_REVERSE
                    end do
#else
                end do
            end do
        end do
#endif
    end subroutine ResScale

    subroutine saGammaReThetaSolve(resOnly)
        !
        !       Residual-first SA-gamma-ReTheta transport solve.
        !       The residual (source + advection + unsteady + viscous)
        !       is assembled for 3 equations and updates are applied
        !       explicitly without the SST ADI solver path.
        !
        use blockPointers
        use constants
        use flowVarRefState
        use inputIteration
        use inputPhysics
        use paramTurb
        implicit none
        !
        !      Subroutine arguments.
        !
        logical, intent(in) :: resOnly
        !
        !      Local variables.
        !
        integer(kind=intType) :: i, j, k, nn
        real(kind=realType) :: rblank, factor
        real(kind=realType), dimension(:, :, :, :), pointer :: dvt

        ! Set a couple of pointers to the correct entries in dw to
        ! make the code more readable.

        dvt => scratch(1:, 1:, 1:, idvt:)


        if (resOnly) return

        ! Explicit update only; implicit ADI solver path intentionally removed.

        factor = one
        if (turbRelax == turbRelaxExplicit) factor = alfaTurb

        do k = 2, kl
            do j = 2, jl
                do i = 2, il
                    rblank = real(iblank(i, j, k), realType)

                    w(i, j, k, itu1) = w(i, j, k, itu1) + factor * dvt(i, j, k, 1) * rblank
                    w(i, j, k, itu1) = max(w(i, j, k, itu1), zero)

                    w(i, j, k, itu2) = w(i, j, k, itu2) + factor * dvt(i, j, k, 2) * rblank
                    w(i, j, k, itu2) = max(zero, min(one, w(i, j, k, itu2)))

                    w(i, j, k, itu3) = w(i, j, k, itu3) + factor * dvt(i, j, k, 3) * rblank
                    w(i, j, k, itu3) = max(w(i, j, k, itu3), 1.e-5_realType * wInf(itu3))
                end do
            end do
        end do

    end subroutine saGammaReThetaSolve

end module saGammaReTheta