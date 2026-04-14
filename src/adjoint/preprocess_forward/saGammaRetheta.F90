










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
        use turbUtils, only: reThetaTCorrelation, flengthCorrelation, rethetacCorrelation, smoothMinMax
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

        ! Gamma-ReTheta source term variables
        real(kind=realType) :: vortMag, strainMag
        real(kind=realType) :: nutSA, rTurb, gammaLocal, reThetaTilde
        real(kind=realType) :: reS_val, reThetaC_val, fLength_val, fTurb_val
        real(kind=realType) :: fOnset, fOnset1
        real(kind=realType) :: vortLim, vortMagLim
        real(kind=realType) :: pGamma, eGamma
        real(kind=realType) :: velMag, velMag2, timeScale, reThetaT_target
        real(kind=realType) :: thetaBL, deltaBL, fWake_val, fThetaT
        real(kind=realType) :: pReTheta, yDist
        real(kind=realType) :: uxhat, uyhat, uzhat, dUds, lambdaThetaLocal



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
                        ! linearization. The source term is stored in scratch.

                        if (approxSA) then
                            term1 = zero
                        else
                            term1 = rsaCb1 * (one - ft2) * ss
                        end if
                        term2 = dist2Inv * (kar2Inv * rsaCb1 * ((one - ft2) * fv2 + ft2) &
                                            - rsaCw1 * fwSa)

                        scratch(i, j, k, idvt) = (term1 + term2 * w(i, j, k, itu1)) * w(i, j, k, itu1)

                        ! ========================================================
                        ! Gamma and ReTheta source terms (Langtry-Menter)
                        ! ========================================================

                        ! --- Compute vorticity and strain magnitudes ---
                        vortx = two * fact * (wwy - vvz) - two * omegax
                        vorty = two * fact * (uuz - wwx) - two * omegay
                        vortz = two * fact * (vvx - uuy) - two * omegaz
                        vortMag = sqrt(max(vortx**2 + vorty**2 + vortz**2, xminn))

                        sxx = two * fact * uux
                        syy = two * fact * vvy
                        szz = two * fact * wwz
                        sxy = fact * (uuy + vvx)
                        sxz = fact * (uuz + wwx)
                        syz = fact * (vvz + wwy)
                        strainMag = sqrt(max(two * (sxy**2 + sxz**2 + syz**2) &
                                        + sxx**2 + syy**2 + szz**2, xminn))

                        ! --- Local variables ---
                        nutSA = w(i, j, k, itu1) * fv1
                        rTurb = nutSA / max(nu, xminn)
                        gammaLocal = max(w(i, j, k, itu2), zero)
                        reThetaTilde = max(w(i, j, k, itu3), xminn)
                        yDist = max(d2Wall(i, j, k), xminn)
                        velMag2 = w(i, j, k, ivx)**2 + w(i, j, k, ivy)**2 &
                                  + w(i, j, k, ivz)**2
                        velMag = sqrt(max(velMag2, xminn))

                        ! --- Vorticity limiting ---
                        vortLim = MachCoef * sqrt(max(MachCoef * Reynolds, xminn)) &
                                  / 20.0_realType
                        vortMagLim = smoothMinMax(vortMag, vortLim, rsaGRvortLimP)

                        ! --- Fonset (smooth tanh-based transition onset) ---
                        reS_val = w(i, j, k, irho) * yDist**2 * strainMag &
                                  / max(rlv(i, j, k), xminn)
                        reThetaC_val = rethetacCorrelation(reThetaTilde)
                        reThetaC_val = max(reThetaC_val, xminn)
                        fOnset1 = sqrt((reS_val / (rsaGRfonsetC * reThetaC_val))**2 &
                                       + rTurb**2)
                        fOnset = (tanh(rsaGRfonsetK &
                                       * (fOnset1 - rsaGRfonsetS)) + one) * half

                        ! --- Flength and Fturb (modified) ---
                        fLength_val = flengthCorrelation(reThetaTilde)
                        fTurb_val = (one - fOnset) * exp(-rTurb)

                        ! --- Gamma production and destruction ---
                        pGamma = rsaGRca1 * fLength_val * fOnset * vortMagLim &
                                 * sqrt(max(gammaLocal, xminn)) &
                                 * (one - rsaGRce1 * gammaLocal)
                        eGamma = rsaGRca2 * fTurb_val * vortMagLim * gammaLocal &
                                 * (rsaGRce2 * gammaLocal - one)

                        scratch(i, j, k, idvt + 1) = pGamma - eGamma

                        ! --- ReTheta production (relaxation toward correlation) ---
                        timeScale = 500.0_realType * rlv(i, j, k) &
                                    / max(w(i, j, k, irho) * velMag2, xminn)

                        ! Compute thetaBL first (needed for lambdaTheta)
                        thetaBL = reThetaTilde * rlv(i, j, k) &
                                  / max(w(i, j, k, irho) * velMag, xminn)

                        ! Compute local lambdaTheta = (thetaBL^2 / nu) * dU/ds
                        uxhat = w(i, j, k, ivx) / max(velMag, xminn)
                        uyhat = w(i, j, k, ivy) / max(velMag, xminn)
                        uzhat = w(i, j, k, ivz) / max(velMag, xminn)
                        dUds = two * fact &
                             * (uxhat * (uxhat * uux + uyhat * uuy + uzhat * uuz) &
                              + uyhat * (uxhat * vvx + uyhat * vvy + uzhat * vvz) &
                              + uzhat * (uxhat * wwx + uyhat * wwy + uzhat * wwz))
                        lambdaThetaLocal = (thetaBL**2 / max(nu, xminn)) * dUds
                        lambdaThetaLocal = max(lambdaThetaLocal, -0.1_realType)
                        lambdaThetaLocal = min(lambdaThetaLocal, 0.1_realType)

                        reThetaT_target = reThetaTCorrelation( &
                            turbIntensityInf * 100.0_realType, lambdaThetaLocal)

                        ! Ftheta_t shielding: shields BL interior, allows
                        ! freestream to drive ReTheta toward correlation value
                        deltaBL = 375.0_realType * vortMag * yDist * thetaBL &
                                  / max(velMag, xminn)
                        deltaBL = max(deltaBL, xminn)
                        fWake_val = exp(-reS_val / 1.0e5_realType)
                        fThetaT = min(fWake_val &
                                  * exp(-(yDist / deltaBL)**4), one)

                        pReTheta = rsaGRcthetat / max(timeScale, xminn) &
                                   * (reThetaT_target - reThetaTilde) &
                                   * (one - fThetaT)

                        scratch(i, j, k, idvt + 2) = pReTheta

                end do
            end do
        end do
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
            do k = 2, kl
                do j = 2, jl
                    do i = 2, il
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
                end do
            end do
        end do
        ! Viscous terms in j-direction.
            do k = 2, kl
                do j = 2, jl
                    do i = 2, il
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
                end do
            end do
        end do
        ! Viscous terms in i-direction.
            do k = 2, kl
                do j = 2, jl
                    do i = 2, il
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
                end do
            end do
        end do
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

            do k = 2, kl
                do j = 2, jl
                    do i = 2, il
                        rblank = max(real(iblank(i, j, k), realType), zero)
                        ! SA variable ν~
                        dw(i,j,k,itu1) = -volRef(i,j,k) * scratch(i, j, k, idvt) * rblank

                        ! Gamma
                        dw(i,j,k,itu2) = -volRef(i,j,k) * scratch(i, j, k, idvt + 1) * rblank

                        ! ReTheta
                        dw(i,j,k,itu3) = -volRef(i,j,k) * scratch(i, j, k, idvt + 2) * rblank
                end do
            end do
        end do
    end subroutine ResScale

    subroutine saGammaReThetaSolve(resOnly)
        !
        !       Point-implicit SA-gamma-ReTheta transport solve.
        !       The residual (source + advection + unsteady + viscous)
        !       is assembled for 3 equations and updates are applied
        !       using a point-implicit scheme (dividing by the diagonal
        !       of the Jacobian qq), analogous to the SA DADI approach.
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
        integer(kind=intType) :: i, j, k, m
        real(kind=realType) :: rblank, factor, delta, theta, wNew

        if (resOnly) return

        ! For implicit relaxation, scale the diagonal by the CFL factor
        ! (same as SA's saSolve: factor = 1 + (1-alfa)/alfa).

        factor = one
        if (turbRelax == turbRelaxImplicit) &
            factor = one + (one - alfaTurb) / alfaTurb

        do k = 2, kl
            do j = 2, jl
                do i = 2, il
                    rblank = real(iblank(i, j, k), realType)

                    ! Scale the diagonal Jacobian for CFL-based damping
                    qq(i, j, k, 1, 1) = factor * qq(i, j, k, 1, 1)
                    qq(i, j, k, 2, 2) = factor * qq(i, j, k, 2, 2)
                    qq(i, j, k, 3, 3) = factor * qq(i, j, k, 3, 3)

                    ! SA equation: point-implicit update (divide by diagonal)
                    if (qq(i, j, k, 1, 1) > zero) then
                        w(i, j, k, itu1) = w(i, j, k, itu1) &
                            + scratch(i, j, k, idvt) * rblank / qq(i, j, k, 1, 1)
                    end if
                    w(i, j, k, itu1) = max(w(i, j, k, itu1), zero)

                    ! Gamma equation: point-implicit update with damping
                    if (qq(i, j, k, 2, 2) > zero) then
                        delta = scratch(i, j, k, idvt + 1) * rblank / qq(i, j, k, 2, 2)
                        theta = one
                        do m = 1, 200
                            wNew = w(i, j, k, itu2) + theta * delta
                            if (wNew > rsaGRgammaLo .and. wNew < rsaGRgammaHi) exit
                            theta = rsaGRdampTheta * theta
                        end do
                        w(i, j, k, itu2) = w(i, j, k, itu2) + theta * delta
                    end if

                    ! ReTheta equation: point-implicit update with damping
                    if (qq(i, j, k, 3, 3) > zero) then
                        delta = scratch(i, j, k, idvt + 2) * rblank / qq(i, j, k, 3, 3)
                        theta = one
                        do m = 1, 200
                            wNew = w(i, j, k, itu3) + theta * delta
                            if (wNew > rsaGRreThetaLo) exit
                            theta = rsaGRdampTheta * theta
                        end do
                        w(i, j, k, itu3) = w(i, j, k, itu3) + theta * delta
                    end if
                end do
            end do
        end do

    end subroutine saGammaReThetaSolve

end module saGammaReTheta
