module saGammaReTheta

    ! This module contains the source code for SA-sLM2015 turbulence
    ! model. It is slightly more modularized than the original which makes
    ! performing reverse mode AD simplier.

    ! transitionDebug slot map (nTransitionDebug = 48; see paramTurb.F90).
    ! Do NOT renumber existing slots — append new ones at the end.
    !  1  fonset            fOnset
    !  2  fonset1           fOnset1
    !  3  flength           fLength_val
    !  4  rturb             rTurb
    !  5  rethetatarget     reThetaT_target
    !  6  res               reS_val
    !  7  rethetac          reThetaC_val
    !  8  resovercrit       reS_val / (2.6 * reThetaC_val)
    !  9  strainmag         strainMag
    ! 10  fthetat           fThetaT
    ! 11  fwake             fWake_val
    ! 12  dudx              dudx
    ! 13  dudy              dudy
    ! 14  dudz              dudz
    ! 15  dvdx              dvdx
    ! 16  dvdy              dvdy
    ! 17  dvdz              dvdz
    ! 18  dwdx              dwdx
    ! 19  dwdy              dwdy
    ! 20  dwdz              dwdz
    ! 21  transgamma        w(i,j,k,itu2)  [γ state; slot 21 = cgnsTurbGamma]
    ! 22  transwalldist     yDist
    ! 23  transrho          w(i,j,k,irho)
    ! 24  transmu           rlv(i,j,k)
    ! 25  gammaprod         pGamma
    ! 26  gammadest         eGamma
    ! 27  transtimescale    timeScale
    ! 28  translambdatheta  lambdaThetaLocal
    ! 29  transpretheta     pReTheta
    ! 30  gammaforsa        gammaForSA  [min(max(γ,0),1) seen by SA production]
    ! 31  gammalocal        gammaLocal  [min(max(γ,γ_lo),1) in transition eqs]
    ! 32  transvortmag      vortMag     [raw |ω| before limiter]
    ! 33  transvortmaglim   vortMagLim  [limited |ω| used in sources]
    ! 34  fturb             fTurb_val   [γ destruction switch exp(-(RT/4)^4)]
    ! 35  sastrainrate      ss          [S, raw strain rate]
    ! 36  samodstrainrate   sst         [S̃, SA-modified strain with wall correction]
    ! 37  ft2               ft2         [SA ft2 function]
    ! 38  thetabl           thetaBL     [BL momentum thickness]
    ! 39  deltabl           deltaBL     [7.5·θ_BL]
    ! 40  transdelta        delta       [BL thickness for fThetaT shielding]
    ! 41  transvelmag       velMag      [|U|]
    ! 42  duds              dUds        [dU/ds streamwise gradient]
    ! 43  nutsa             nutSA       [w(itu1)·fv1, SA eddy viscosity]
    ! 44  rethetatilde      reThetaTilde [max(w(itu3),1)]
    ! 45  transvortlim      vortLim     [vorticity limiter threshold]
    ! 46  qq11              qq(i,j,k,1,1) [SA diagonal Jacobian]
    ! 47  qq22              qq(i,j,k,2,2) [γ diagonal Jacobian]
    ! 48  qq33              qq(i,j,k,3,3) [Re̅θt diagonal Jacobian]

    use constants, only: realType, zero

    real(kind=realType), dimension(:, :, :, :, :), allocatable :: qq

contains

    subroutine saGammaReTheta_block(resOnly)

        use constants
        use blockPointers, only: il, jl, kl
        use inputTimeSpectral
        use iteration
        use turbUtils, only: SSTEddyViscosity, turbAdvection, unsteadyTurbTerm, saEddyViscosity
        use turbBCRoutines, only: bcTurbTreatment, applyAllTurbBCThisBlock
        use inputIteration, only: transitionFirstOrderUpwind
        use inputDiscretization, only: orderTurb
        implicit none

        !
        !      Subroutine argument.
        !
        logical, intent(in) :: resOnly
        !
        !      Local variables.
        !
        integer(kind=intType) :: nn, sps
        integer(kind=intType) :: orderTurbSave

        ! Set the arrays for the boundary condition treatment.

        call bcTurbTreatment

        ! Alloc central jacobian memory
        allocate (qq(2:il,2:jl,2:kl,3,3))

        ! Source Terms
        call Source

        ! Advection Term
        nn = itu1 - 1
        if (transitionFirstOrderUpwind) then
            orderTurbSave = orderTurb
            orderTurb = firstOrder
        end if
        call turbAdvection(3_intType, 3_intType, nn, qq)
        if (transitionFirstOrderUpwind) then
            orderTurb = orderTurbSave
        end if

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
        integer(kind=intType), parameter :: dbgFonset = 1_intType
        integer(kind=intType), parameter :: dbgFonset1 = 2_intType
        integer(kind=intType), parameter :: dbgFlength = 3_intType
        integer(kind=intType), parameter :: dbgRturb = 4_intType
        integer(kind=intType), parameter :: dbgReThetaTarget = 5_intType
        integer(kind=intType), parameter :: dbgReS = 6_intType
        integer(kind=intType), parameter :: dbgReThetaC = 7_intType
        integer(kind=intType), parameter :: dbgReSOverCrit = 8_intType
        integer(kind=intType), parameter :: dbgStrainMag = 9_intType
        integer(kind=intType), parameter :: dbgFthetaT = 10_intType
        integer(kind=intType), parameter :: dbgFwake = 11_intType
        integer(kind=intType), parameter :: dbgDudx = 12_intType
        integer(kind=intType), parameter :: dbgDudy = 13_intType
        integer(kind=intType), parameter :: dbgDudz = 14_intType
        integer(kind=intType), parameter :: dbgDvdx = 15_intType
        integer(kind=intType), parameter :: dbgDvdy = 16_intType
        integer(kind=intType), parameter :: dbgDvdz = 17_intType
        integer(kind=intType), parameter :: dbgDwdx = 18_intType
        integer(kind=intType), parameter :: dbgDwdy = 19_intType
        integer(kind=intType), parameter :: dbgDwdz = 20_intType
        integer(kind=intType), parameter :: dbgGamma = 21_intType
        integer(kind=intType), parameter :: dbgWallDist = 22_intType
        integer(kind=intType), parameter :: dbgRho = 23_intType
        integer(kind=intType), parameter :: dbgMu = 24_intType
        integer(kind=intType), parameter :: dbgGammaProd = 25_intType
        integer(kind=intType), parameter :: dbgGammaDest = 26_intType
        integer(kind=intType), parameter :: dbgTimeScale = 27_intType
        integer(kind=intType), parameter :: dbgLambdaTheta = 28_intType
        integer(kind=intType), parameter :: dbgPReTheta = 29_intType
        integer(kind=intType), parameter :: dbgGammaForSA = 30_intType
        integer(kind=intType), parameter :: dbgGammaLocal = 31_intType
        integer(kind=intType), parameter :: dbgVortMag = 32_intType
        integer(kind=intType), parameter :: dbgVortMagLim = 33_intType
        integer(kind=intType), parameter :: dbgFturb = 34_intType
        integer(kind=intType), parameter :: dbgSS = 35_intType
        integer(kind=intType), parameter :: dbgSST = 36_intType
        integer(kind=intType), parameter :: dbgFt2 = 37_intType
        integer(kind=intType), parameter :: dbgThetaBL = 38_intType
        integer(kind=intType), parameter :: dbgDeltaBL = 39_intType
        integer(kind=intType), parameter :: dbgDelta = 40_intType
        integer(kind=intType), parameter :: dbgVelMag = 41_intType
        integer(kind=intType), parameter :: dbgDUds = 42_intType
        integer(kind=intType), parameter :: dbgNutSA = 43_intType
        integer(kind=intType), parameter :: dbgReThetaTilde = 44_intType
        integer(kind=intType), parameter :: dbgVortLim = 45_intType
        integer(kind=intType), parameter :: dbgQQ11 = 46_intType
        integer(kind=intType), parameter :: dbgQQ22 = 47_intType
        integer(kind=intType), parameter :: dbgQQ33 = 48_intType

        ! Local variables.
        integer(kind=intType) :: i, j, k, nn, ii
        real(kind=realType) :: fv1, fv2, ft2
        real(kind=realType) :: ss, sst, nu, dist2Inv, chi, chi2, chi3
        real(kind=realType) :: rr, gg, gg6, termFw, fwSa, term1, term2
        real(kind=realType) :: term2_prod, term2_dest
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
        real(kind=realType) :: gammaForSA
        real(kind=realType) :: reS_val, reThetaC_val, fLength_val, fTurb_val
        real(kind=realType) :: fOnset, fOnset1
        real(kind=realType) :: vortLim, vortMagLim
        real(kind=realType) :: pGamma, eGamma
        real(kind=realType) :: velMag, velMag2, timeScale, reThetaT_target
        real(kind=realType) :: thetaBL, deltaBL, delta, fWake_val, fThetaT
        real(kind=realType) :: gammaEff, gammaTerm
        real(kind=realType) :: pReTheta, yDist
        real(kind=realType) :: uxhat, uyhat, uzhat, dUds, lambdaThetaLocal
        real(kind=realType) :: dudx, dudy, dudz, dvdx, dvdy, dvdz
        real(kind=realType) :: dwdx, dwdy, dwdz
        real(kind=realType) :: epsRT, reThetaTilde_p, reThetaC_p
        real(kind=realType) :: fOnset1_p, fOnset_p, fLength_p, pGamma_p
        real(kind=realType) :: drTurb_dnu, dfTurb_dnu, dfOnset_dnu
        real(kind=realType) :: dfOnset1_drT, dfOnset_dfOnset1



        ! Set model constants
        cv13 = rsaCv1**3
        kar2Inv = one / (rsaK**2)
        cw36 = rsaCw3**6
        cb3Inv = one / rsaCb3

#ifndef USE_TAPENADE
        ! Initialize the full 3x3 turbulence Jacobian block to avoid
        ! uninitialized off-diagonal entries contaminating ANK/DADI.
        qq = zero
#endif

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

                        ! Clamp gamma to the physical range for SA production coupling.
                        gammaForSA = min(max(w(i, j, k, itu2), zero), one)

                        if (approxSA) then
                            term1 = zero
                        else
                            term1 = gammaForSA * rsaCb1 * (one - ft2) * ss
                        end if

                        ! Split term2 into production and destruction parts.
                        ! Production: near-wall correction from Cb1*fv2/(kappa^2*d^2)
                        ! Destruction: -Cw1*fw/d^2
                        ! Gamma multiplies only production.
                        term2_prod = dist2Inv * kar2Inv * rsaCb1 &
                                     * ((one - ft2) * fv2 + ft2)
                        term2_dest = -dist2Inv * rsaCw1 * fwSa

                        ! Effective term2 with gamma on production only
                        term2 = gammaForSA * term2_prod + term2_dest

                        scratch(i, j, k, idvt) = (term1 + term2 * w(i, j, k, itu1)) * w(i, j, k, itu1)

                        ! ========================================================
                        ! Gamma and ReTheta source terms (sLangtry-Menter)
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

                        strainMag2 = two*(sxy**2 + sxz**2 + syz**2) + sxx**2 + syy**2 + szz**2
                        strainMag  = sqrt(max(two*strainMag2, xminn))

                        ! --- Local variables ---
                        !v_t= ν̃ · fv1 is the SA eddy viscosity
                        nutSA = w(i, j, k, itu1) * fv1
                        !rTurb= ν_t/ν
                        rTurb = nutSA / nu
                        gammaLocal = min(max(w(i, j, k, itu2), rsaGRgammaLo), rsaGRgammaHi)
                        reThetaTilde = max(w(i, j, k, itu3), rsaGRreThetaLo)
                        yDist = d2Wall(i, j, k)
                        velMag2 = w(i, j, k, ivx)**2 + w(i, j, k, ivy)**2 &
                                  + w(i, j, k, ivz)**2
                        velMag = sqrt(max(velMag2, xminn))

                        ! --- Vorticity limiting ---
                        ! ADflow nondim of paper Eqs. 52–53. Paper writes M·√(M·Re)/20
                        ! using a∞ as velocity scale; ADflow uses √(p/ρ) as velocity
                        ! scale and √(p*ρ) for dynamic viscosity.
                        ! Rotating frame not adress here!!!!! uInf has no meaning on it.
                        vortLim = uInf * sqrt(max(uInf / max(muInf, xminn), xminn)) &
                                / 20.0_realType

                        vortMagLim = smoothMinMax(vortMag, vortLim, rsaGRpmin)

                        ! --- Fonset (smooth tanh-based transition onset) ---
                        reS_val = w(i, j, k, irho) * yDist**2 * strainMag &
                                  / rlv(i, j, k)
                        reThetaC_val = rethetacCorrelation(reThetaTilde)                       
                        fOnset1 = sqrt((reS_val / (2.6_realType * reThetaC_val))**2 &
                                       + rTurb**2)
                        fOnset = (tanh(6.0_realType * (fOnset1 - 1.35_realType)) + one) * half

                        ! --- Flength and Fturb (modified) ---
                        fLength_val = flengthCorrelation(reThetaTilde)
                        fTurb_val = (one - fOnset) * exp(-rTurb)
                        !Check here if needed 
                        !fTurb_val = exp(-(rTurb / 4.0_realType)**4)

                        ! --- Gamma production and destruction ---
                        pGamma = rsaGRca1 * fLength_val * fOnset * vortMagLim &
                                 * sqrt(gammaLocal) &
                                 * (one - rsaGRce1 * gammaLocal)
                        eGamma = rsaGRca2 * fTurb_val * vortMagLim * gammaLocal &
                                 * (rsaGRce2 * gammaLocal - one)

                        scratch(i, j, k, idvt + 1) = pGamma - eGamma
                        

                        ! --- ReTheta production (relaxation toward correlation) ---
                        ! NOTE: No explicit Reynolds factor here.
                        ! Nondim form: rlv = mu/muRef with L_ref=1m so Re=1 implicitly
                        ! (see initializeFlow.F90:62-66). No explicit Re factor needed,
                        ! consistent with nu = rlv/rho in sa.F90:245.
                        timeScale = 500.0_realType * nu / max(velMag2, xminn)

                        ! Compute thetaBL first (needed for lambdaTheta)
                        thetaBL = reThetaTilde * nu &
                                  / max(velMag, xminn)

                        ! Compute local lambdaTheta = (thetaBL^2 / nu) * dU/ds
                        uxhat = w(i, j, k, ivx) / max(velMag, xminn)
                        uyhat = w(i, j, k, ivy) / max(velMag, xminn)
                        uzhat = w(i, j, k, ivz) / max(velMag, xminn)
                        dUds = two * fact &
                             * (uxhat * (uxhat * uux + uyhat * uuy + uzhat * uuz) &
                              + uyhat * (uxhat * vvx + uyhat * vvy + uzhat * vvz) &
                              + uzhat * (uxhat * wwx + uyhat * wwy + uzhat * wwz))
                        lambdaThetaLocal = (thetaBL**2 / nu) * dUds
                        lambdaThetaLocal = smoothMinMax(lambdaThetaLocal, -0.1_realType, rsaGRpmax)
                        lambdaThetaLocal = smoothMinMax(lambdaThetaLocal, 0.1_realType, rsaGRpmin)

                        reThetaT_target = reThetaTCorrelation( &
                            turbIntensityInf * 100.0_realType, lambdaThetaLocal)

                        ! Ftheta_t shielding: shields BL interior, allows
                        ! freestream to drive ReTheta toward correlation value
                        deltaBL = 7.5_realType * thetaBL
                        delta = 50.0_realType * yDist * vortMag * deltaBL &
                            / max(velMag, xminn)
                        delta = max(delta, xminn)
                        fWake_val = exp(-reS_val / 1.0e6_realType)
                        fThetaT    = fWake_val * exp(-(yDist/delta)**4)

                        pReTheta = rsaGRcthetat / max(timeScale, xminn) &
                                   * (reThetaT_target - reThetaTilde) &
                                   * (one - fThetaT)

                        scratch(i, j, k, idvt + 2) = pReTheta

                        if (associated(transitionDebug)) then
                            dudx = two * fact * uux
                            dudy = two * fact * uuy
                            dudz = two * fact * uuz
                            dvdx = two * fact * vvx
                            dvdy = two * fact * vvy
                            dvdz = two * fact * vvz
                            dwdx = two * fact * wwx
                            dwdy = two * fact * wwy
                            dwdz = two * fact * wwz

                            transitionDebug(i, j, k, dbgFonset) = fOnset
                            transitionDebug(i, j, k, dbgFonset1) = fOnset1
                            transitionDebug(i, j, k, dbgFlength) = fLength_val
                            transitionDebug(i, j, k, dbgRturb) = rTurb
                            transitionDebug(i, j, k, dbgReThetaTarget) = reThetaT_target
                            transitionDebug(i, j, k, dbgReS) = reS_val
                            transitionDebug(i, j, k, dbgReThetaC) = reThetaC_val
                            transitionDebug(i, j, k, dbgReSOverCrit) = reS_val / (2.6_realType * reThetaC_val)
                            transitionDebug(i, j, k, dbgStrainMag) = strainMag
                            transitionDebug(i, j, k, dbgFthetaT) = fThetaT
                            transitionDebug(i, j, k, dbgFwake) = fWake_val
                            transitionDebug(i, j, k, dbgDudx) = dudx
                            transitionDebug(i, j, k, dbgDudy) = dudy
                            transitionDebug(i, j, k, dbgDudz) = dudz
                            transitionDebug(i, j, k, dbgDvdx) = dvdx
                            transitionDebug(i, j, k, dbgDvdy) = dvdy
                            transitionDebug(i, j, k, dbgDvdz) = dvdz
                            transitionDebug(i, j, k, dbgDwdx) = dwdx
                            transitionDebug(i, j, k, dbgDwdy) = dwdy
                            transitionDebug(i, j, k, dbgDwdz) = dwdz
                            transitionDebug(i, j, k, dbgGamma) = w(i, j, k, itu2)
                            transitionDebug(i, j, k, dbgWallDist) = yDist
                            transitionDebug(i, j, k, dbgRho) = w(i, j, k, irho)
                            transitionDebug(i, j, k, dbgMu) = rlv(i, j, k)
                            transitionDebug(i, j, k, dbgGammaProd) = pGamma
                            transitionDebug(i, j, k, dbgGammaDest) = eGamma
                            transitionDebug(i, j, k, dbgTimeScale) = timeScale
                            transitionDebug(i, j, k, dbgLambdaTheta) = lambdaThetaLocal
                            transitionDebug(i, j, k, dbgPReTheta) = pReTheta
                            transitionDebug(i, j, k, dbgGammaForSA) = gammaForSA
                            transitionDebug(i, j, k, dbgGammaLocal) = gammaLocal
                            transitionDebug(i, j, k, dbgVortMag) = vortMag
                            transitionDebug(i, j, k, dbgVortMagLim) = vortMagLim
                            transitionDebug(i, j, k, dbgFturb) = fTurb_val
                            transitionDebug(i, j, k, dbgSS) = ss
                            transitionDebug(i, j, k, dbgSST) = sst
                            transitionDebug(i, j, k, dbgFt2) = ft2
                            transitionDebug(i, j, k, dbgThetaBL) = thetaBL
                            transitionDebug(i, j, k, dbgDeltaBL) = deltaBL
                            transitionDebug(i, j, k, dbgDelta) = delta
                            transitionDebug(i, j, k, dbgVelMag) = velMag
                            transitionDebug(i, j, k, dbgDUds) = dUds
                            transitionDebug(i, j, k, dbgNutSA) = nutSA
                            transitionDebug(i, j, k, dbgReThetaTilde) = reThetaTilde
                            transitionDebug(i, j, k, dbgVortLim) = vortLim
                        end if

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

                        ! Compute the source term jacobian.
                        ! term2 already contains gammaForSA on production part.
                        ! The derivative chain also needs gamma on production
                        ! derivatives but NOT on destruction derivatives.
                        ! Note that -dsource/dnu is stored.
                        qq(i, j, k, 1, 1) = -two * term2 * w(i, j, k, itu1) &
                                      - dist2Inv * w(i, j, k, itu1) * w(i, j, k, itu1) &
                                      * (gammaForSA * rsaCb1 * kar2Inv &
                                         * (dfv2 - ft2 * dfv2 - fv2 * dft2 + dft2) &
                                         - rsaCw1 * dfw)

                        ! Full gamma source Jacobian: -d(pGamma - eGamma)/dgamma
                        ! Includes both production and destruction linearization.
                        qq(i, j, k, 2, 2) = rsaGRca1 * fLength_val * fOnset &
                            * vortMagLim &
                            * (1.5_realType * rsaGRce1 * gammaLocal - half) &
                            / sqrt(max(gammaLocal, xminn)) &
                            + rsaGRca2 * fTurb_val * vortMagLim &
                            * (two * rsaGRce2 * gammaLocal - one)

                        ! ReTheta Jacobian diagonal: -dP_theta/dReTheta_tilde
                        ! Always >= 0 (relaxation), so restriction is a no-op.
                        qq(i, j, k, 3, 3) = max( &
                            rsaGRcthetat / max(timeScale, xminn) &
                            * (one - fThetaT), zero)

                        ! --- Off-diagonal source Jacobian entries ---

                        ! qq(1,2) = -dS_nu/dgamma: SA production depends on gamma
                        qq(i, j, k, 1, 2) = -(rsaCb1 * (one - ft2) * ss &
                            + term2_prod * w(i, j, k, itu1)) * w(i, j, k, itu1)

                        ! qq(1,3) = -dS_nu/dReThetaTilde: ~0 (paper §7.1)
                        qq(i, j, k, 1, 3) = zero

                        ! qq(2,1) = -dS_gamma/dnu_tilde: both pGamma and
                        ! eGamma depend on nu_tilde through rTurb.
                        ! rTurb = nu_tilde * fv1 / nu
                        drTurb_dnu = (fv1 + chi * dfv1) / nu

                        ! dfOnset/dnu = dfOnset/dfOnset1 * dfOnset1/drTurb * drTurb/dnu
                        ! fOnset1 = sqrt((reS/(2.6*reThetaC))^2 + rTurb^2)
                        dfOnset1_drT = rTurb / max(fOnset1, xminn)
                        ! fOnset = 0.5*(tanh(6*(fOnset1-1.35))+1)
                        dfOnset_dfOnset1 = 12.0_realType * fOnset * (one - fOnset)
                        dfOnset_dnu = dfOnset_dfOnset1 * dfOnset1_drT * drTurb_dnu

                        ! dfTurb/dnu: fTurb = (1-fOnset)*exp(-rTurb)
                        ! Full chain rule: d/dnu = d/dfOnset*dfOnset/dnu + d/drTurb*drTurb/dnu
                        dfTurb_dnu = -exp(-rTurb) * dfOnset_dnu - fTurb_val * drTurb_dnu

                        ! -d(pGamma - eGamma)/dnu_tilde
                        qq(i, j, k, 2, 1) = &
                            -(rsaGRca1 * fLength_val * dfOnset_dnu * vortMagLim &
                              * sqrt(max(gammaLocal, xminn)) &
                              * (one - rsaGRce1 * gammaLocal)) &
                            + rsaGRca2 * dfTurb_dnu * vortMagLim * gammaLocal &
                              * (rsaGRce2 * gammaLocal - one)

                        ! qq(2,3) = -dS_gamma/dReThetaTilde: P_gamma depends on
                        ! fOnset(reThetaC(reThetaTilde)) and fLength(reThetaTilde).
                        ! Use one-sided finite difference on P_gamma.
                        epsRT = max(1.0e-4_realType * reThetaTilde, 1.0e-2_realType)
                        reThetaTilde_p = reThetaTilde + epsRT
                        reThetaC_p = rethetacCorrelation(reThetaTilde_p)
                        reThetaC_p = max(reThetaC_p, xminn)
                        fOnset1_p = sqrt((reS_val &
                            / (2.6_realType * reThetaC_p))**2 + rTurb**2)
                        fOnset_p = (tanh(6.0_realType &
                            * (fOnset1_p - 1.35_realType)) + one) * half
                        fLength_p = flengthCorrelation(reThetaTilde_p)
                        pGamma_p = rsaGRca1 * fLength_p * fOnset_p * vortMagLim &
                            * sqrt(max(gammaLocal, xminn)) &
                            * (one - rsaGRce1 * gammaLocal)
                        qq(i, j, k, 2, 3) = -(pGamma_p - pGamma) / epsRT

                        ! qq(3,1) = -dS_retheta/dnu_tilde: ~0 (paper §7.1)
                        qq(i, j, k, 3, 1) = zero

                        ! qq(3,2) = -dS_retheta/dgamma: fThetaT does not
                        ! depend on gamma in the current implementation.
                        qq(i, j, k, 3, 2) = zero

                        if (associated(transitionDebug)) then
                            transitionDebug(i, j, k, dbgQQ11) = qq(i, j, k, 1, 1)
                            transitionDebug(i, j, k, dbgQQ22) = qq(i, j, k, 2, 2)
                            transitionDebug(i, j, k, dbgQQ33) = qq(i, j, k, 3, 3)
                        end if
#endif
#ifdef TAPENADE_REVERSE
                    end do
#else
                end do
            end do
        end do
#endif

        call computeSrcLambda()

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
                            qq(i,j,k,1,2) = qq(i,j,k,1,2) - b1 * bmtk1(i,j,itu1,itu2)
                            qq(i,j,k,1,3) = qq(i,j,k,1,3) - b1 * bmtk1(i,j,itu1,itu3)

                        else if (k == kl) then
                            qq(i,j,k,1,1) = qq(i,j,k,1,1) + c1 &
                                            - d1 * max(bmtk2(i,j,itu1,itu1), zero)
                            qq(i,j,k,1,2) = qq(i,j,k,1,2) - d1 * bmtk2(i,j,itu1,itu2)
                            qq(i,j,k,1,3) = qq(i,j,k,1,3) - d1 * bmtk2(i,j,itu1,itu3)

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
                            qq(i,j,k,2,1) = qq(i,j,k,2,1) - b2 * bmtk1(i,j,itu2,itu1)
                            qq(i,j,k,2,2) = qq(i,j,k,2,2) + c2 &
                                            - b2 * max(bmtk1(i,j,itu2,itu2), zero)
                            qq(i,j,k,2,3) = qq(i,j,k,2,3) - b2 * bmtk1(i,j,itu2,itu3)

                        else if (k == kl) then
                            qq(i,j,k,2,1) = qq(i,j,k,2,1) - d2 * bmtk2(i,j,itu2,itu1)
                            qq(i,j,k,2,2) = qq(i,j,k,2,2) + c2 &
                                            - d2 * max(bmtk2(i,j,itu2,itu2), zero)
                            qq(i,j,k,2,3) = qq(i,j,k,2,3) - d2 * bmtk2(i,j,itu2,itu3)

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
                            qq(i,j,k,3,1) = qq(i,j,k,3,1) - b3 * bmtk1(i,j,itu3,itu1)
                            qq(i,j,k,3,2) = qq(i,j,k,3,2) - b3 * bmtk1(i,j,itu3,itu2)
                            qq(i,j,k,3,3) = qq(i,j,k,3,3) + c3 &
                                            - b3 * max(bmtk1(i,j,itu3,itu3), zero)

                        else if (k == kl) then
                            qq(i,j,k,3,1) = qq(i,j,k,3,1) - d3 * bmtk2(i,j,itu3,itu1)
                            qq(i,j,k,3,2) = qq(i,j,k,3,2) - d3 * bmtk2(i,j,itu3,itu2)
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

                        if (j == 2) then
                            qq(i,j,k,1,1) = qq(i,j,k,1,1) + c1 &
                                            - b1 * max(bmtj1(i,k,itu1,itu1), zero)
                            qq(i,j,k,1,2) = qq(i,j,k,1,2) - b1 * bmtj1(i,k,itu1,itu2)
                            qq(i,j,k,1,3) = qq(i,j,k,1,3) - b1 * bmtj1(i,k,itu1,itu3)

                        else if (j == jl) then
                            qq(i,j,k,1,1) = qq(i,j,k,1,1) + c1 &
                                            - d1 * max(bmtj2(i,k,itu1,itu1), zero)
                            qq(i,j,k,1,2) = qq(i,j,k,1,2) - d1 * bmtj2(i,k,itu1,itu2)
                            qq(i,j,k,1,3) = qq(i,j,k,1,3) - d1 * bmtj2(i,k,itu1,itu3)

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
                            qq(i,j,k,2,1) = qq(i,j,k,2,1) - b2 * bmtj1(i,k,itu2,itu1)
                            qq(i,j,k,2,2) = qq(i,j,k,2,2) + c2 &
                                            - b2 * max(bmtj1(i,k,itu2,itu2), zero)
                            qq(i,j,k,2,3) = qq(i,j,k,2,3) - b2 * bmtj1(i,k,itu2,itu3)

                        else if (j == jl) then
                            qq(i,j,k,2,1) = qq(i,j,k,2,1) - d2 * bmtj2(i,k,itu2,itu1)
                            qq(i,j,k,2,2) = qq(i,j,k,2,2) + c2 &
                                            - d2 * max(bmtj2(i,k,itu2,itu2), zero)
                            qq(i,j,k,2,3) = qq(i,j,k,2,3) - d2 * bmtj2(i,k,itu2,itu3)

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
                            qq(i,j,k,3,1) = qq(i,j,k,3,1) - b3 * bmtj1(i,k,itu3,itu1)
                            qq(i,j,k,3,2) = qq(i,j,k,3,2) - b3 * bmtj1(i,k,itu3,itu2)
                            qq(i,j,k,3,3) = qq(i,j,k,3,3) + c3 &
                                            - b3 * max(bmtj1(i,k,itu3,itu3), zero)

                        else if (j == jl) then
                            qq(i,j,k,3,1) = qq(i,j,k,3,1) - d3 * bmtj2(i,k,itu3,itu1)
                            qq(i,j,k,3,2) = qq(i,j,k,3,2) - d3 * bmtj2(i,k,itu3,itu2)
                            qq(i,j,k,3,3) = qq(i,j,k,3,3) + c3 &
                                            - d3 * max(bmtj2(i,k,itu3,itu3), zero)

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

                        if (i == 2) then
                            qq(i,j,k,1,1) = qq(i,j,k,1,1) + c1 &
                                            - b1 * max(bmti1(j,k,itu1,itu1), zero)
                            qq(i,j,k,1,2) = qq(i,j,k,1,2) - b1 * bmti1(j,k,itu1,itu2)
                            qq(i,j,k,1,3) = qq(i,j,k,1,3) - b1 * bmti1(j,k,itu1,itu3)

                        else if (i == il) then
                            qq(i,j,k,1,1) = qq(i,j,k,1,1) + c1 &
                                            - d1 * max(bmti2(j,k,itu1,itu1), zero)
                            qq(i,j,k,1,2) = qq(i,j,k,1,2) - d1 * bmti2(j,k,itu1,itu2)
                            qq(i,j,k,1,3) = qq(i,j,k,1,3) - d1 * bmti2(j,k,itu1,itu3)

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
                            qq(i,j,k,2,1) = qq(i,j,k,2,1) - b2 * bmti1(j,k,itu2,itu1)
                            qq(i,j,k,2,2) = qq(i,j,k,2,2) + c2 &
                                            - b2 * max(bmti1(j,k,itu2,itu2), zero)
                            qq(i,j,k,2,3) = qq(i,j,k,2,3) - b2 * bmti1(j,k,itu2,itu3)

                        else if (i == il) then
                            qq(i,j,k,2,1) = qq(i,j,k,2,1) - d2 * bmti2(j,k,itu2,itu1)
                            qq(i,j,k,2,2) = qq(i,j,k,2,2) + c2 &
                                            - d2 * max(bmti2(j,k,itu2,itu2), zero)
                            qq(i,j,k,2,3) = qq(i,j,k,2,3) - d2 * bmti2(j,k,itu2,itu3)

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
                            qq(i,j,k,3,1) = qq(i,j,k,3,1) - b3 * bmti1(j,k,itu3,itu1)
                            qq(i,j,k,3,2) = qq(i,j,k,3,2) - b3 * bmti1(j,k,itu3,itu2)
                            qq(i,j,k,3,3) = qq(i,j,k,3,3) + c3 &
                                            - b3 * max(bmti1(j,k,itu3,itu3), zero)

                        else if (i == il) then
                            qq(i,j,k,3,1) = qq(i,j,k,3,1) - d3 * bmti2(j,k,itu3,itu1)
                            qq(i,j,k,3,2) = qq(i,j,k,3,2) - d3 * bmti2(j,k,itu3,itu2)
                            qq(i,j,k,3,3) = qq(i,j,k,3,3) + c3 &
                                            - d3 * max(bmti2(j,k,itu3,itu3), zero)

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
                        dw(i,j,k,itu2) = -volRef(i,j,k) * scratch(i, j, k, idvt + 1) * rblank

                        ! ReTheta
                        dw(i,j,k,itu3) = -volRef(i,j,k) * scratch(i, j, k, idvt + 2) * rblank
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
        !       Coupled DD-ADI SA-gamma-ReTheta transport solve.
        !       Solves the 3-equation system using a diagonally-dominant
        !       alternating-direction-implicit scheme with full 3x3 block
        !       coupling, following the same pattern as the SST solver.
        !
        use blockPointers
        use constants
        use flowVarRefState
        use inputIteration
        use inputPhysics
        use paramTurb
        use turbUtils, only: tdia3x3
        implicit none
        !
        !      Subroutine arguments.
        !
        logical, intent(in) :: resOnly
        !
        !      Local variables.
        !
        integer(kind=intType) :: i, j, k

        ! viscosity variables
        real(kind=realType) :: nu, nu_m, nu_p
        real(kind=realType) :: nut, nut_m, nut_p
        real(kind=realType) :: nu_tm, nu_tp
        real(kind=realType) :: nuTilde, nuTilde_m, nuTilde_p
        real(kind=realType) :: chi, chi3, fv1
        real(kind=realType) :: chi_m, chi3_m, fv1_m
        real(kind=realType) :: chi_p, chi3_p, fv1_p

        ! geometry terms
        real(kind=realType) :: voli, volmi, volpi
        real(kind=realType) :: xm, ym, zm, xp, yp, zp, xa, ya, za
        real(kind=realType) :: ttm, ttp

        ! diffusion coefficients
        real(kind=realType) :: num_v, nup_v
        real(kind=realType) :: cdm, cdp
        real(kind=realType) :: cdm_gamma, cdp_gamma
        real(kind=realType) :: cdm_rt, cdp_rt

        ! SA nonlinear correction
        real(kind=realType) :: cnud, cam, cap
        real(kind=realType) :: nutm, nutp

        ! clipped diffusion coefficients
        real(kind=realType) :: c1m, c1p, c2m, c2p, c3m, c3p

        ! advection
        real(kind=realType) :: qs, uu, um, up

        ! constants
        real(kind=realType) :: cb3Inv, cv13

        ! misc
        real(kind=realType) :: rblank, factor
        real(kind=realType) :: gammaNew, gammaDelta, dampFactor
        integer(kind=intType) :: mm

        ! Source dt restriction (Eq. 59)
        real(kind=realType) :: dt_inv

        ! Scaling values from existing turbulence residual scaling options
        real(kind=realType) :: scaleNu, scaleGamma, scaleReTheta

        ! Scaling ratios (precomputed from scaling values)
        real(kind=realType) :: s12, s13, s21, s23, s31, s32

        ! ADI work arrays
        real(kind=realType), dimension(3, 2:max(il, jl, kl)) :: bb, dd, ff
        real(kind=realType), dimension(3, 3, 2:max(il, jl, kl)) :: cc

        logical, save :: printedCoupling = .false.

        if (resOnly) return

        if (.not. printedCoupling) then
            printedCoupling = .true.
            select case (TurbDADICoupled)
            case (0)
                print *, 'SA-gamma-ReTheta DADI coupling: fully decoupled (diagonal only)'
            case (1)
                print *, 'SA-gamma-ReTheta DADI coupling: SA decoupled, gamma-ReTheta coupled'
            case (2)
                print *, 'SA-gamma-ReTheta DADI coupling: fully coupled (3x3 block)'
            case default
                print *, 'SA-gamma-ReTheta DADI coupling: unknown mode', TurbDADICoupled
            end select
        end if

        cb3Inv = one / rsaCb3
        cv13 = rsaCv1**3

        ! Use existing per-equation scaling controls. Protect against zero.
        scaleNu = max(abs(turbResScale(1)), one)
        scaleGamma = max(abs(turbResScale(2)), one)
        scaleReTheta = max(abs(turbResScale(3)), one)

        ! Precompute scaling ratios: s_col / s_row for off-diagonal entries
        s12 = scaleGamma   / scaleNu       ! eq1, var2
        s13 = scaleReTheta / scaleNu       ! eq1, var3
        s21 = scaleNu      / scaleGamma    ! eq2, var1
        s23 = scaleReTheta / scaleGamma    ! eq2, var3
        s31 = scaleNu      / scaleReTheta  ! eq3, var1
        s32 = scaleGamma   / scaleReTheta  ! eq3, var2

        ! Prepare qq: decoupled mode, scaling, and CFL factor.
        ! For implicit relaxation: factor = 1 + (1-alfa)/alfa = 1/alfa.

        factor = one
        if (turbRelax == turbRelaxImplicit) &
            factor = one + (one - alfaTurb) / alfaTurb

        do k = 2, kl
            do j = 2, jl
                do i = 2, il
                    ! TurbDADICoupled: 0=decoupled, 1=transition-only, 2=full
                    if (TurbDADICoupled == 0) then
                        qq(i, j, k, 1, 2) = zero
                        qq(i, j, k, 1, 3) = zero
                        qq(i, j, k, 2, 1) = zero
                        qq(i, j, k, 2, 3) = zero
                        qq(i, j, k, 3, 1) = zero
                        qq(i, j, k, 3, 2) = zero
                    else if (TurbDADICoupled == 1) then
                        ! SA decoupled, gamma-retheta coupled
                        qq(i, j, k, 1, 2) = zero
                        qq(i, j, k, 1, 3) = zero
                        qq(i, j, k, 2, 1) = zero
                        qq(i, j, k, 3, 1) = zero
                    end if

                    ! Source dt restriction (Eq. 59): additive I/Δt inflation
                    if (transitionSrcDtRestrict .and. srcDtRestrictActive) then
                        dt_inv = srcLambda(i,j,k) / transitionSrcDtLimit
                        qq(i,j,k,1,1) = qq(i,j,k,1,1) + dt_inv
                        qq(i,j,k,2,2) = qq(i,j,k,2,2) + dt_inv
                        qq(i,j,k,3,3) = qq(i,j,k,3,3) + dt_inv
                    end if

                    ! Symmetric scaling (§4): qq(m,n) *= s_n / s_m
                    ! Diagonal entries unchanged; only off-diag scaled.
                    qq(i, j, k, 1, 2) = qq(i, j, k, 1, 2) * s12
                    qq(i, j, k, 1, 3) = qq(i, j, k, 1, 3) * s13
                    qq(i, j, k, 2, 1) = qq(i, j, k, 2, 1) * s21
                    qq(i, j, k, 2, 3) = qq(i, j, k, 2, 3) * s23
                    qq(i, j, k, 3, 1) = qq(i, j, k, 3, 1) * s31
                    qq(i, j, k, 3, 2) = qq(i, j, k, 3, 2) * s32

                    ! CFL factor scaling (all 9 entries)
                    qq(i, j, k, 1, 1) = factor * qq(i, j, k, 1, 1)
                    qq(i, j, k, 1, 2) = factor * qq(i, j, k, 1, 2)
                    qq(i, j, k, 1, 3) = factor * qq(i, j, k, 1, 3)
                    qq(i, j, k, 2, 1) = factor * qq(i, j, k, 2, 1)
                    qq(i, j, k, 2, 2) = factor * qq(i, j, k, 2, 2)
                    qq(i, j, k, 2, 3) = factor * qq(i, j, k, 2, 3)
                    qq(i, j, k, 3, 1) = factor * qq(i, j, k, 3, 1)
                    qq(i, j, k, 3, 2) = factor * qq(i, j, k, 3, 2)
                    qq(i, j, k, 3, 3) = factor * qq(i, j, k, 3, 3)

                    ! Scale the RHS: scratch(m) /= s_m
                    scratch(i, j, k, idvt)     = scratch(i, j, k, idvt)     / scaleNu
                    scratch(i, j, k, idvt + 1) = scratch(i, j, k, idvt + 1) / scaleGamma
                    scratch(i, j, k, idvt + 2) = scratch(i, j, k, idvt + 2) / scaleReTheta
                end do
            end do
        end do

        ! Initialize grid velocity to zero.
        qs = zero
        !
        !       DD-ADI step in j-direction. As we solve in j-direction,
        !       the j-loop is the innermost loop.
        !
        do k = 2, kl
            do i = 2, il
                do j = 2, jl

                    ! Recompute viscous diffusion coefficients in j-direction.

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

                    nu = rlv(i, j, k) / w(i, j, k, irho)
                    nuTilde = w(i, j, k, itu1)
                    chi = nuTilde / nu; chi3 = chi * chi * chi
                    fv1 = chi3 / (chi3 + cv13)
                    nut = nuTilde * fv1

                    nu_m = rlv(i, j - 1, k) / w(i, j - 1, k, irho)
                    nuTilde_m = w(i, j - 1, k, itu1)
                    chi_m = nuTilde_m / nu_m; chi3_m = chi_m * chi_m * chi_m
                    fv1_m = chi3_m / (chi3_m + cv13)
                    nut_m = nuTilde_m * fv1_m

                    nu_p = rlv(i, j + 1, k) / w(i, j + 1, k, irho)
                    nuTilde_p = w(i, j + 1, k, itu1)
                    chi_p = nuTilde_p / nu_p; chi3_p = chi_p * chi_p * chi_p
                    fv1_p = chi3_p / (chi3_p + cv13)
                    nut_p = nuTilde_p * fv1_p

                    num_v = half * (nu_m + nu)
                    nup_v = half * (nu_p + nu)
                    nu_tm = half * (nut_m + nut)
                    nu_tp = half * (nut_p + nut)

                    ! SA diffusion
                    cdm = (num_v + (one + rsaCb2) * nutm) * ttm * cb3Inv
                    cdp = (nup_v + (one + rsaCb2) * nutp) * ttp * cb3Inv

                    ! Gamma diffusion
                    cdm_gamma = (num_v + nu_tm / sigmaF) * ttm
                    cdp_gamma = (nup_v + nu_tp / sigmaF) * ttp

                    ! ReTheta diffusion
                    cdm_rt = sigmaTheta * (num_v + nu_tm) * ttm
                    cdp_rt = sigmaTheta * (nup_v + nu_tp) * ttp

                    c1m = max(cdm + cam, zero)
                    c1p = max(cdp + cap, zero)
                    c2m = max(cdm_gamma, zero)
                    c2p = max(cdp_gamma, zero)
                    c3m = max(cdm_rt, zero)
                    c3p = max(cdp_rt, zero)

                    bb(1, j) = -c1m
                    dd(1, j) = -c1p
                    bb(2, j) = -c2m
                    dd(2, j) = -c2p
                    bb(3, j) = -c3m
                    dd(3, j) = -c3p

                    ! Add advection off-diagonal terms in j-direction.

                    if (addGridVelocities) &
                        qs = half * (sFaceJ(i, j, k) + sFaceJ(i, j - 1, k)) * voli

                    uu = xa * w(i, j, k, ivx) + ya * w(i, j, k, ivy) &
                         + za * w(i, j, k, ivz) - qs
                    um = zero; up = zero
                    if (uu < zero) um = uu
                    if (uu > zero) up = uu

                    bb(1, j) = bb(1, j) - up
                    dd(1, j) = dd(1, j) + um
                    bb(2, j) = bb(2, j) - up
                    dd(2, j) = dd(2, j) + um
                    bb(3, j) = bb(3, j) - up
                    dd(3, j) = dd(3, j) + um

                    ! Store central jacobian and rhs in cc and ff.
                    ! Multiply off-diagonal qq entries and rhs by iblank
                    ! so the update for iblank=0 cells is zero.

                    rblank = real(iblank(i, j, k), realType)

                    cc(1, 1, j) = qq(i, j, k, 1, 1)
                    cc(1, 2, j) = qq(i, j, k, 1, 2) * rblank
                    cc(1, 3, j) = qq(i, j, k, 1, 3) * rblank
                    cc(2, 1, j) = qq(i, j, k, 2, 1) * rblank
                    cc(2, 2, j) = qq(i, j, k, 2, 2)
                    cc(2, 3, j) = qq(i, j, k, 2, 3) * rblank
                    cc(3, 1, j) = qq(i, j, k, 3, 1) * rblank
                    cc(3, 2, j) = qq(i, j, k, 3, 2) * rblank
                    cc(3, 3, j) = qq(i, j, k, 3, 3)

                    ff(1, j) = scratch(i, j, k, idvt) * rblank
                    ff(2, j) = scratch(i, j, k, idvt + 1) * rblank
                    ff(3, j) = scratch(i, j, k, idvt + 2) * rblank

                    bb(:, j) = bb(:, j) * rblank
                    dd(:, j) = dd(:, j) * rblank

                end do

                ! Solve the tri-diagonal system in j-direction.

                call tdia3x3(2_intType, jl, bb, cc, dd, ff)

                ! Determine the new rhs for the next direction.

                do j = 2, jl
                    scratch(i, j, k, idvt) = qq(i, j, k, 1, 1) * ff(1, j) &
                        + qq(i, j, k, 1, 2) * ff(2, j) + qq(i, j, k, 1, 3) * ff(3, j)
                    scratch(i, j, k, idvt + 1) = qq(i, j, k, 2, 1) * ff(1, j) &
                        + qq(i, j, k, 2, 2) * ff(2, j) + qq(i, j, k, 2, 3) * ff(3, j)
                    scratch(i, j, k, idvt + 2) = qq(i, j, k, 3, 1) * ff(1, j) &
                        + qq(i, j, k, 3, 2) * ff(2, j) + qq(i, j, k, 3, 3) * ff(3, j)
                end do

            end do
        end do
        !
        !       DD-ADI step in i-direction. As we solve in i-direction,
        !       the i-loop is the innermost loop.
        !
        do k = 2, kl
            do j = 2, jl
                do i = 2, il

                    ! Recompute viscous diffusion coefficients in i-direction.

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

                    nu = rlv(i, j, k) / w(i, j, k, irho)
                    nuTilde = w(i, j, k, itu1)
                    chi = nuTilde / nu; chi3 = chi * chi * chi
                    fv1 = chi3 / (chi3 + cv13)
                    nut = nuTilde * fv1

                    nu_m = rlv(i - 1, j, k) / w(i - 1, j, k, irho)
                    nuTilde_m = w(i - 1, j, k, itu1)
                    chi_m = nuTilde_m / nu_m; chi3_m = chi_m * chi_m * chi_m
                    fv1_m = chi3_m / (chi3_m + cv13)
                    nut_m = nuTilde_m * fv1_m

                    nu_p = rlv(i + 1, j, k) / w(i + 1, j, k, irho)
                    nuTilde_p = w(i + 1, j, k, itu1)
                    chi_p = nuTilde_p / nu_p; chi3_p = chi_p * chi_p * chi_p
                    fv1_p = chi3_p / (chi3_p + cv13)
                    nut_p = nuTilde_p * fv1_p

                    num_v = half * (nu_m + nu)
                    nup_v = half * (nu_p + nu)
                    nu_tm = half * (nut_m + nut)
                    nu_tp = half * (nut_p + nut)

                    ! SA diffusion
                    cdm = (num_v + (one + rsaCb2) * nutm) * ttm * cb3Inv
                    cdp = (nup_v + (one + rsaCb2) * nutp) * ttp * cb3Inv

                    ! Gamma diffusion
                    cdm_gamma = (num_v + nu_tm / sigmaF) * ttm
                    cdp_gamma = (nup_v + nu_tp / sigmaF) * ttp

                    ! ReTheta diffusion
                    cdm_rt = sigmaTheta * (num_v + nu_tm) * ttm
                    cdp_rt = sigmaTheta * (nup_v + nu_tp) * ttp

                    c1m = max(cdm + cam, zero)
                    c1p = max(cdp + cap, zero)
                    c2m = max(cdm_gamma, zero)
                    c2p = max(cdp_gamma, zero)
                    c3m = max(cdm_rt, zero)
                    c3p = max(cdp_rt, zero)

                    bb(1, i) = -c1m
                    dd(1, i) = -c1p
                    bb(2, i) = -c2m
                    dd(2, i) = -c2p
                    bb(3, i) = -c3m
                    dd(3, i) = -c3p

                    ! Add advection off-diagonal terms in i-direction.

                    if (addGridVelocities) &
                        qs = half * (sFaceI(i, j, k) + sFaceI(i - 1, j, k)) * voli

                    uu = xa * w(i, j, k, ivx) + ya * w(i, j, k, ivy) &
                         + za * w(i, j, k, ivz) - qs
                    um = zero; up = zero
                    if (uu < zero) um = uu
                    if (uu > zero) up = uu

                    bb(1, i) = bb(1, i) - up
                    dd(1, i) = dd(1, i) + um
                    bb(2, i) = bb(2, i) - up
                    dd(2, i) = dd(2, i) + um
                    bb(3, i) = bb(3, i) - up
                    dd(3, i) = dd(3, i) + um

                    ! Store central jacobian and rhs in cc and ff.

                    rblank = real(iblank(i, j, k), realType)

                    cc(1, 1, i) = qq(i, j, k, 1, 1)
                    cc(1, 2, i) = qq(i, j, k, 1, 2) * rblank
                    cc(1, 3, i) = qq(i, j, k, 1, 3) * rblank
                    cc(2, 1, i) = qq(i, j, k, 2, 1) * rblank
                    cc(2, 2, i) = qq(i, j, k, 2, 2)
                    cc(2, 3, i) = qq(i, j, k, 2, 3) * rblank
                    cc(3, 1, i) = qq(i, j, k, 3, 1) * rblank
                    cc(3, 2, i) = qq(i, j, k, 3, 2) * rblank
                    cc(3, 3, i) = qq(i, j, k, 3, 3)

                    ff(1, i) = scratch(i, j, k, idvt) * rblank
                    ff(2, i) = scratch(i, j, k, idvt + 1) * rblank
                    ff(3, i) = scratch(i, j, k, idvt + 2) * rblank

                    bb(:, i) = bb(:, i) * rblank
                    dd(:, i) = dd(:, i) * rblank

                end do

                ! Solve the tri-diagonal system in i-direction.

                call tdia3x3(2_intType, il, bb, cc, dd, ff)

                ! Determine the new rhs for the next direction.

                do i = 2, il
                    scratch(i, j, k, idvt) = qq(i, j, k, 1, 1) * ff(1, i) &
                        + qq(i, j, k, 1, 2) * ff(2, i) + qq(i, j, k, 1, 3) * ff(3, i)
                    scratch(i, j, k, idvt + 1) = qq(i, j, k, 2, 1) * ff(1, i) &
                        + qq(i, j, k, 2, 2) * ff(2, i) + qq(i, j, k, 2, 3) * ff(3, i)
                    scratch(i, j, k, idvt + 2) = qq(i, j, k, 3, 1) * ff(1, i) &
                        + qq(i, j, k, 3, 2) * ff(2, i) + qq(i, j, k, 3, 3) * ff(3, i)
                end do

            end do
        end do
        !
        !       DD-ADI step in k-direction. As we solve in k-direction,
        !       the k-loop is the innermost loop.
        !
        do j = 2, jl
            do i = 2, il
                do k = 2, kl

                    ! Recompute viscous diffusion coefficients in k-direction.

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

                    nu = rlv(i, j, k) / w(i, j, k, irho)
                    nuTilde = w(i, j, k, itu1)
                    chi = nuTilde / nu; chi3 = chi * chi * chi
                    fv1 = chi3 / (chi3 + cv13)
                    nut = nuTilde * fv1

                    nu_m = rlv(i, j, k - 1) / w(i, j, k - 1, irho)
                    nuTilde_m = w(i, j, k - 1, itu1)
                    chi_m = nuTilde_m / nu_m; chi3_m = chi_m * chi_m * chi_m
                    fv1_m = chi3_m / (chi3_m + cv13)
                    nut_m = nuTilde_m * fv1_m

                    nu_p = rlv(i, j, k + 1) / w(i, j, k + 1, irho)
                    nuTilde_p = w(i, j, k + 1, itu1)
                    chi_p = nuTilde_p / nu_p; chi3_p = chi_p * chi_p * chi_p
                    fv1_p = chi3_p / (chi3_p + cv13)
                    nut_p = nuTilde_p * fv1_p

                    num_v = half * (nu_m + nu)
                    nup_v = half * (nu_p + nu)
                    nu_tm = half * (nut_m + nut)
                    nu_tp = half * (nut_p + nut)

                    ! SA diffusion
                    cdm = (num_v + (one + rsaCb2) * nutm) * ttm * cb3Inv
                    cdp = (nup_v + (one + rsaCb2) * nutp) * ttp * cb3Inv

                    ! Gamma diffusion
                    cdm_gamma = (num_v + nu_tm / sigmaF) * ttm
                    cdp_gamma = (nup_v + nu_tp / sigmaF) * ttp

                    ! ReTheta diffusion
                    cdm_rt = sigmaTheta * (num_v + nu_tm) * ttm
                    cdp_rt = sigmaTheta * (nup_v + nu_tp) * ttp

                    c1m = max(cdm + cam, zero)
                    c1p = max(cdp + cap, zero)
                    c2m = max(cdm_gamma, zero)
                    c2p = max(cdp_gamma, zero)
                    c3m = max(cdm_rt, zero)
                    c3p = max(cdp_rt, zero)

                    bb(1, k) = -c1m
                    dd(1, k) = -c1p
                    bb(2, k) = -c2m
                    dd(2, k) = -c2p
                    bb(3, k) = -c3m
                    dd(3, k) = -c3p

                    ! Add advection off-diagonal terms in k-direction.

                    if (addGridVelocities) &
                        qs = half * (sFaceK(i, j, k) + sFaceK(i, j, k - 1)) * voli

                    uu = xa * w(i, j, k, ivx) + ya * w(i, j, k, ivy) &
                         + za * w(i, j, k, ivz) - qs
                    um = zero; up = zero
                    if (uu < zero) um = uu
                    if (uu > zero) up = uu

                    bb(1, k) = bb(1, k) - up
                    dd(1, k) = dd(1, k) + um
                    bb(2, k) = bb(2, k) - up
                    dd(2, k) = dd(2, k) + um
                    bb(3, k) = bb(3, k) - up
                    dd(3, k) = dd(3, k) + um

                    ! Store central jacobian and rhs in cc and ff.

                    rblank = real(iblank(i, j, k), realType)

                    cc(1, 1, k) = qq(i, j, k, 1, 1)
                    cc(1, 2, k) = qq(i, j, k, 1, 2) * rblank
                    cc(1, 3, k) = qq(i, j, k, 1, 3) * rblank
                    cc(2, 1, k) = qq(i, j, k, 2, 1) * rblank
                    cc(2, 2, k) = qq(i, j, k, 2, 2)
                    cc(2, 3, k) = qq(i, j, k, 2, 3) * rblank
                    cc(3, 1, k) = qq(i, j, k, 3, 1) * rblank
                    cc(3, 2, k) = qq(i, j, k, 3, 2) * rblank
                    cc(3, 3, k) = qq(i, j, k, 3, 3)

                    ff(1, k) = scratch(i, j, k, idvt) * rblank
                    ff(2, k) = scratch(i, j, k, idvt + 1) * rblank
                    ff(3, k) = scratch(i, j, k, idvt + 2) * rblank

                    bb(:, k) = bb(:, k) * rblank
                    dd(:, k) = dd(:, k) * rblank

                end do

                ! Solve the tri-diagonal system in k-direction.

                call tdia3x3(2_intType, kl, bb, cc, dd, ff)

                ! Store the final update in scratch.

                do k = 2, kl
                    scratch(i, j, k, idvt) = ff(1, k)
                    scratch(i, j, k, idvt + 1) = ff(2, k)
                    scratch(i, j, k, idvt + 2) = ff(3, k)
                end do

            end do
        end do
        !
        !       Update the turbulent variables. For explicit relaxation the
        !       update must be relaxed; for implicit relaxation this has been
        !       done via the time step.
        !       Unscale the update: ΔQ(m) = s_m * ΔQ_scaled(m).
        !
        factor = one
        if (turbRelax == turbRelaxExplicit) factor = alfaTurb

        do k = 2, kl
            do j = 2, jl
                do i = 2, il
                    w(i, j, k, itu1) = w(i, j, k, itu1) &
                                       + factor * scaleNu &
                                       * scratch(i, j, k, idvt)
                    w(i, j, k, itu1) = max(w(i, j, k, itu1), zero)

                    ! Gamma update with exponential back-off damping (§3).
                    ! If the raw update overshoots [gammaLo, gammaHi],
                    ! reduce the step by theta^m until it stays in range.
                    gammaDelta = factor * scaleGamma &
                                 * scratch(i, j, k, idvt + 1)
                    gammaNew = w(i, j, k, itu2) + gammaDelta
                    dampFactor = one
                    do mm = 1, 40
                        if (gammaNew >= rsaGRgammaLo .and. &
                            gammaNew <= rsaGRgammaHi) exit
                        dampFactor = dampFactor * rsaGRdampTheta
                        gammaNew = w(i, j, k, itu2) + dampFactor * gammaDelta
                    end do
                    w(i, j, k, itu2) = min(max(gammaNew, rsaGRgammaLo), &
                                           rsaGRgammaHi)

                    ! ReTheta update with exponential back-off damping (§3).
                    gammaDelta = factor * scaleReTheta &
                                 * scratch(i, j, k, idvt + 2)
                    gammaNew = w(i, j, k, itu3) + gammaDelta
                    dampFactor = one
                    do mm = 1, 40
                        if (gammaNew >= rsaGRreThetaLo) exit
                        dampFactor = dampFactor * rsaGRdampTheta
                        gammaNew = w(i, j, k, itu3) + dampFactor * gammaDelta
                    end do
                    w(i, j, k, itu3) = max(gammaNew, rsaGRreThetaLo)
                end do
            end do
        end do

    end subroutine saGammaReThetaSolve

    subroutine computeSrcLambda()
        ! Compute srcLambda = max(0, λ_max) where λ_max is the largest positive
        ! eigenvalue of A_source = -qq (P&Z 2020 Eq. 59).
        ! Note: qq stores -∂S/∂Q, so A_source = -qq.
        !
        ! Mode 0: Signed Gershgorin upper bound (AD-safe)
        !   srcLambda = max(0, max_i[ A_ii + Σ_{j≠i} |A_ij| ])
        ! Mode 1: Exact 3x3 eigenvalue via cubic formula

        use constants
        use blockPointers, only: il, jl, kl, srcLambda
        use inputIteration, only: transitionSrcDtEigMode
        implicit none

        integer(kind=intType) :: i, j, k
        real(kind=realType) :: A11, A12, A13, A21, A22, A23, A31, A32, A33
        real(kind=realType) :: g1, g2, g3, lambdaMax
        real(kind=realType) :: p, qc, r, a, b, c, disc, phi, t, sqrtP
        real(kind=realType), parameter :: oneThird = one / three
        real(kind=realType), parameter :: piVal = 3.14159265358979323846_realType

        do k = 2, kl
            do j = 2, jl
                do i = 2, il
                    ! A_source = -qq (qq stores -∂S/∂Q)
                    A11 = -qq(i,j,k,1,1)
                    A12 = -qq(i,j,k,1,2)
                    A13 = -qq(i,j,k,1,3)
                    A21 = -qq(i,j,k,2,1)
                    A22 = -qq(i,j,k,2,2)
                    A23 = -qq(i,j,k,2,3)
                    A31 = -qq(i,j,k,3,1)
                    A32 = -qq(i,j,k,3,2)
                    A33 = -qq(i,j,k,3,3)

                    if (transitionSrcDtEigMode == 0) then
                        ! Mode 0: Signed Gershgorin upper bound
                        g1 = A11 + abs(A12) + abs(A13)
                        g2 = A22 + abs(A21) + abs(A23)
                        g3 = A33 + abs(A31) + abs(A32)
                        lambdaMax = max(g1, g2, g3)
                    else
                        ! Mode 1: Exact eigenvalues via cubic formula
                        ! Characteristic polynomial: λ³ - aλ² + bλ - c = 0
                        ! where a = tr(A), b = (tr(A)² - tr(A²))/2, c = det(A)
                        a = A11 + A22 + A33
                        b = A11*A22 + A22*A33 + A33*A11 &
                          - A12*A21 - A23*A32 - A31*A13
                        c = A11*(A22*A33 - A23*A32) &
                          - A12*(A21*A33 - A23*A31) &
                          + A13*(A21*A32 - A22*A31)

                        ! Depressed cubic: t³ + pt + qc = 0, λ = t + a/3
                        p = b - a*a*oneThird
                        qc = two*a*a*a/27.0_realType - a*b*oneThird + c
                        disc = qc*qc/four + p*p*p/27.0_realType

                        if (disc <= zero) then
                            ! Three real roots (trigonometric solution)
                            sqrtP = sqrt(-p*oneThird)
                            if (sqrtP > 1.0e-30_realType) then
                                phi = acos(max(-one, min(one, -qc/(two*sqrtP**3))))
                                t = two * sqrtP * cos(phi * oneThird)
                            else
                                t = zero
                            end if
                            lambdaMax = t + a * oneThird
                            t = two * sqrtP * cos((phi + two*piVal) * oneThird)
                            lambdaMax = max(lambdaMax, t + a * oneThird)
                            t = two * sqrtP * cos((phi + four*piVal) * oneThird)
                            lambdaMax = max(lambdaMax, t + a * oneThird)
                        else
                            ! One real root, two complex conjugates
                            r = sqrt(disc)
                            t = sign(abs(-qc*half + r)**oneThird, -qc*half + r) &
                              + sign(abs(-qc*half - r)**oneThird, -qc*half - r)
                            lambdaMax = t + a * oneThird
                        end if

                        ! Fallback to Gershgorin if eigenvalue computation produces NaN
                        if (lambdaMax /= lambdaMax) then
                            g1 = A11 + abs(A12) + abs(A13)
                            g2 = A22 + abs(A21) + abs(A23)
                            g3 = A33 + abs(A31) + abs(A32)
                            lambdaMax = max(g1, g2, g3)
                        end if
                    end if

                    srcLambda(i,j,k) = max(zero, lambdaMax)
                end do
            end do
        end do

    end subroutine computeSrcLambda

end module saGammaReTheta