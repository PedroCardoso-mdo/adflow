module saGamma

    ! SA-noft2-Gamma ("SA-sgammama"): two-equation transition model
    ! (nuTilde = itu1, gamma = itu2, NO ReThetaTilde).
    !
    ! Model: the one-equation gamma transition model of Menter et al. (2015)
    ! coupled to Spalart-Allmaras as in Nichols (2019) / Lee (2021, UMD
    ! dissertation, Ch. 3) with constant freestream Tu (Lee & Baeder 2021) and
    ! the Colonia et al. (2016) Re_theta_c constants (via Jung et al. 2022),
    ! smoothed with the Piotrowski & Zingg (2020) machinery (smoothMinMax
    ! phi_{+-p} on every min/max kink). Equations and the decision record:
    ! ~/Desktop/Run/MDO_PhD/Transition/sa_sgammama/00_proposal/RESUMO_FINAL_SA_sgammama.md
    !
    !   1  Dgammama/Dt = P_g - E_g + div[(nu + nu_t/sigma_g) grad gamma],  sigma_g = 1
    !   2  P_g = Flength * S * g (1 - g) * F_onset,  Flength = 100
    !   3  E_g = c_a2 * Omega * g * F_turb * (c_e2 g - 1),  c_a2 = 0.06, c_e2 = 50
    !   4  F_onset1 = Re_v/(2.2 Re_theta_c); F2 = phi-(F1, 2); F3 = phi+(1-(R_T/3.5)^3, 0);
    !      F_onset = phi+(F2 - F3, 0)          [option sgammamaOnsetTanh: P&Z tanh form]
    !   5  F_turb = (1 - F_onset) exp(-R_T)   [option sgammamaFturbLee: exp(-(R_T/2)^4)]
    !   6  Flength = 100
    !   7  Re_theta_c = C_TU1 + C_TU2 exp(-C_TU3 Tu F_PG),
    !      F_PG = 1 + phi-(14.68 phi+(lam,0), 0.5) + phi-(-7.34 phi-(lam,0), 2.0)
    !   8  lam = 7.57e-3 (d^2/nu) dU/ds + 0.0128, clipped to [-1, 1] (smooth)
    !   9  Tu = turbIntensityInf (constant)
    !  10  gammama_s = phi+((phi-(g,1) - 1/c_e2)/(1 - 1/c_e2), 0);
    !      P_nu *= gammama_s;  D_nu *= phi+(gammama_s, 0.1)  [option sgammamaCoupleDestruction]
    !
    ! The module has the same public interface as saGammaRetheta
    ! (saGamma_block, saGammaSolve, sgSource, sgViscous, sgResScale, qq,
    ! computeSrcLambdaSaGamma) so the solver/adjoint dispatch is symmetric.
    ! Blocks marked LOCKSTEP are copies of saGammaRetheta.F90 (minus the
    ! ReThetaTilde equation) and must be kept in step with it.
    !
    ! transitionDebug slot map (reuses the GR slot numbers where the meaning
    ! is the same, so the volume/surface output plumbing is unchanged):
    !  1  fonset      F_onset          2  fonset1     F_onset1
    !  3  flength     Flength (100)    4  rturb       R_T
    !  5  rethetatarget  F_PG          6  res         Re_v
    !  7  rethetac    Re_theta_c       8  resovercrit Re_v/(2.2 Re_theta_c)
    !  9  strainmag   S               10  fthetat     gammama_s (SA production multiplier)
    ! 11  fwake       SA destruction multiplier phi+(gammama_s, 0.1)
    ! 12-20 velocity gradients       21  transgammama  gamma
    ! 22  transwalldist d            23  transrho    rho
    ! 24  transmu     mu             25  gammamaprod   P_g
    ! 26  gammamadest   E_g            28  translambdatheta  lambda_thetaL (clipped)
    ! 30  gammamaforsa  gammama_s        31  gammamalocal  clipped gamma in the sources
    ! 32  transvortmag Omega         33  transvortmaglim  S/Omega after limiter
    ! 34  fturb       F_turb         35  sastrainrate ss   36  samodstrainrate sst
    ! 37  ft2         ft2            41  transvelmag |U|   42  duds  dU/ds
    ! 43  nutsa       nu_t           45  transvortlim vortLim
    ! 46  qq11        qq(1,1)        47  qq22        qq(2,2)

    use constants, only: realType, zero

    real(kind=realType), dimension(:, :, :, :, :), allocatable :: qq

    ! Model constants (Menter 2015 / Lee 2021 Eq. 3.7, 3.12, 3.15; Nichols coupling)
    real(kind=realType), parameter :: sgFlength = 100.0_realType
    real(kind=realType), parameter :: sgFonsetC = 2.2_realType
    real(kind=realType), parameter :: sgFonsetRT = 3.5_realType
    real(kind=realType), parameter :: sgCPG1 = 14.68_realType
    real(kind=realType), parameter :: sgCPG2 = -7.34_realType
    real(kind=realType), parameter :: sgCPG1lim = 1.5_realType
    real(kind=realType), parameter :: sgCPG2lim = 3.0_realType
    real(kind=realType), parameter :: sgLamCoef = 7.57e-3_realType
    real(kind=realType), parameter :: sgLamOffset = 0.0128_realType
    real(kind=realType), parameter :: sgLamClip = 1.0_realType
    real(kind=realType), parameter :: sgDestFloor = 0.1_realType
    ! P&Z tanh onset (plan B, option sgammaOnsetTanh)
    real(kind=realType), parameter :: sgTanhK = 6.0_realType
    real(kind=realType), parameter :: sgTanhS = 1.35_realType

contains

    subroutine saGamma_block(resOnly)

        use constants
        use blockPointers, only: il, jl, kl
        use inputTimeSpectral
        use iteration
        use turbUtils, only: turbAdvection, unsteadyTurbTerm, saEddyViscosity
        use turbBCRoutines, only: bcTurbTreatment, applyAllTurbBCThisBlock
        implicit none

        logical, intent(in) :: resOnly
        integer(kind=intType) :: nn

        call bcTurbTreatment

        allocate (qq(2:il, 2:jl, 2:kl, 2, 2))

        call sgSource
        nn = itu1 - 1
        call turbAdvection(2_intType, 2_intType, nn, qq)
        call unsteadyTurbTerm(2_intType, 2_intType, nn, qq)

        call sgViscous
        call sgResScale
        if (.not. resOnly) then
            call saGammaSolve(resOnly)
            call saEddyViscosity(2, il, 2, jl, 2, kl)
            call applyAllTurbBCThisBlock(.true.)
        end if
        deallocate (qq)

    end subroutine saGamma_block

    subroutine cellSources(nuTilde, gam, nu, dist2Inv, ss, strainMag, vortMag, &
                           yDist, dUds, rho, mu, tuPct, vortLim, ft2, approxTerm1, &
                           useLimiter, onsetTanh, fturbLee, coupleDest, pSmooth, &
                           cTU1, cTU2, cTU3, &
                           sNu, sGamma, gammaS, destMult, term2prod, term2dest, &
                           fv1, fv2, sst, rr, gg, gg6, termFw, fwSa, &
                           reVort, lamL, fPG, reThetaC, fOnset1, fOnset, fTurb, &
                           pGamma, eGamma, gammaLocal, sEff, omEff, rTurb, nutSA)
        !
        ! Per-cell SA and gam sources of SA-sgamma from the state
        ! (nuTilde, gam) and the cell kinematics. Used by Source (residual,
        ! differentiated) and, under #ifndef USE_TAPENADE, for the finite-
        ! difference LHS Jacobian and by evalSrcJacBlockSaGamma (Eq. 59).
        ! Every smoothMinMax has its own target variable (no in-place
        ! update on the differentiated path -- fast-reverse rule).
        !
        use constants
        use paramTurb, only: rsaCv1, rsaK, rsaCw3, rsaCb1, rsaCw1, rsaCw2, rsaCt3, rsaCt4, &
                             rsaGRca2, rsaGRce2, rsaGRgammaLo, rsaGRgammaHi, rsaGRpmin, rsaGRpmax
        use turbUtils, only: smoothMinMax
        implicit none

        real(kind=realType), intent(in) :: nuTilde, gam, nu, dist2Inv, ss, strainMag, vortMag
        real(kind=realType), intent(in) :: yDist, dUds, rho, mu, tuPct, vortLim, ft2
        logical, intent(in) :: approxTerm1, useLimiter, onsetTanh, fturbLee, coupleDest
        real(kind=realType), intent(in) :: pSmooth, cTU1, cTU2, cTU3
        real(kind=realType), intent(out) :: sNu, sGamma, gammaS, destMult, term2prod, term2dest
        real(kind=realType), intent(out) :: fv1, fv2, sst, rr, gg, gg6, termFw, fwSa
        real(kind=realType), intent(out) :: reVort, lamL, fPG, reThetaC, fOnset1, fOnset, fTurb
        real(kind=realType), intent(out) :: pGamma, eGamma, gammaLocal, sEff, omEff, rTurb, nutSA

        real(kind=realType), parameter :: xminn = 1.e-10_realType
        real(kind=realType) :: cv13, kar2Inv, cw36, chi, chi2, chi3
        real(kind=realType) :: term1, term2
        real(kind=realType) :: gClip1, gS0, lamRaw, lamLo, lamPlus, lamMinus
        real(kind=realType) :: fPGplusRaw, fPGplus, fPGminusRaw, fPGminus
        real(kind=realType) :: fOnset2, fOnset3raw, fOnset3, fOnsetRaw, fOnset1t

        cv13 = rsaCv1**3
        kar2Inv = one / (rsaK**2)
        cw36 = rsaCw3**6

        ! ---------------- SA part (LOCKSTEP saGammaRetheta.F90 Source) ----------------
        chi = nuTilde / nu
        chi2 = chi * chi
        chi3 = chi * chi2
        fv1 = chi3 / (chi3 + cv13)
        fv2 = one - chi / (one + chi * fv1)

        sst = ss + nuTilde * fv2 * kar2Inv * dist2Inv
        sst = max(sst, xminn)

        rr = nuTilde * kar2Inv * dist2Inv / sst
        rr = min(rr, 10.0_realType)
        gg = rr + rsaCw2 * (rr**6 - rr)
        gg6 = gg**6
        termFw = ((one + cw36) / (gg6 + cw36))**sixth
        fwSa = gg * termFw

        ! Eq. 10: rescaled intermittency gamma_s (Nichols/Lee 3.18), smooth
        gClip1 = smoothMinMax(gam, one, rsaGRpmin)
        gS0 = (gClip1 - one / rsaGRce2) / (one - one / rsaGRce2)
        gammaS = smoothMinMax(gS0, zero, rsaGRpmax)
        if (coupleDest) then
            destMult = smoothMinMax(gammaS, sgDestFloor, rsaGRpmax)
        else
            destMult = one
        end if

        if (approxTerm1) then
            term1 = zero
        else
            term1 = gammaS * rsaCb1 * (one - ft2) * ss
        end if
        term2prod = dist2Inv * kar2Inv * rsaCb1 * ((one - ft2) * fv2 + ft2)
        term2dest = -dist2Inv * rsaCw1 * fwSa
        term2 = gammaS * term2prod + destMult * term2dest
        sNu = (term1 + term2 * nuTilde) * nuTilde

        ! ---------------- gam part ----------------
        nutSA = nuTilde * fv1
        rTurb = nutSA / nu
        gammaLocal = min(max(gam, rsaGRgammaLo), rsaGRgammaHi)

        ! Eq. 8: local pressure-gradient parameter (continuity form), smooth clip [-1, 1]
        lamRaw = sgLamCoef * yDist**2 / nu * dUds + sgLamOffset
        lamLo = smoothMinMax(lamRaw, -sgLamClip, rsaGRpmax)
        lamL = smoothMinMax(lamLo, sgLamClip, rsaGRpmin)

        ! Eq. 7: F_PG (V with caps) and Re_theta_c
        lamPlus = smoothMinMax(lamL, zero, pSmooth)
        lamMinus = smoothMinMax(lamL, zero, -pSmooth)
        fPGplusRaw = sgCPG1 * lamPlus
        fPGplus = smoothMinMax(fPGplusRaw, sgCPG1lim - one, -pSmooth)
        fPGminusRaw = sgCPG2 * lamMinus
        fPGminus = smoothMinMax(fPGminusRaw, sgCPG2lim - one, -pSmooth)
        fPG = one + fPGplus + fPGminus
        reThetaC = cTU1 + cTU2 * exp(-cTU3 * tuPct * fPG)

        ! Eq. 4: onset
        reVort = rho * yDist**2 * strainMag / mu
        if (onsetTanh) then
            fOnset1t = sqrt((reVort / (sgFonsetC * reThetaC))**2 + rTurb**2)
            fOnset1 = fOnset1t
            fOnset = (tanh(sgTanhK * (fOnset1t - sgTanhS)) + one) * half
        else
            fOnset1 = reVort / (sgFonsetC * reThetaC)
            fOnset2 = smoothMinMax(fOnset1, two, rsaGRpmin)
            fOnset3raw = one - (rTurb / sgFonsetRT)**3
            fOnset3 = smoothMinMax(fOnset3raw, zero, rsaGRpmax)
            fOnsetRaw = fOnset2 - fOnset3
            fOnset = smoothMinMax(fOnsetRaw, zero, rsaGRpmax)
        end if

        ! Eq. 5: F_turb
        if (fturbLee) then
            fTurb = exp(-(rTurb / two)**4)
        else
            fTurb = (one - fOnset) * exp(-rTurb)
        end if

        ! Eqs. 2-3: production with S, destruction with Omega (optional P&Z floor)
        if (useLimiter) then
            sEff = smoothMinMax(strainMag, vortLim, rsaGRpmin)
            omEff = smoothMinMax(vortMag, vortLim, rsaGRpmin)
        else
            sEff = strainMag
            omEff = vortMag
        end if
        pGamma = sgFlength * sEff * gammaLocal * (one - gammaLocal) * fOnset
        eGamma = rsaGRca2 * omEff * gammaLocal * fTurb * (rsaGRce2 * gammaLocal - one)
        sGamma = pGamma - eGamma

    end subroutine cellSources

    subroutine cellKinematics(i, j, k, ss, strainMag, vortMag, dUds, velMag, vortLim, &
                              uux, uuy, uuz, vvx, vvy, vvz, wwx, wwy, wwz, fact)
        !
        ! Cell-centre velocity gradients, SA production magnitude ss (strain or
        ! vorticity per turbProd), strain and vorticity magnitudes, streamwise
        ! velocity gradient dU/ds (relative frame) and the P&Z vorticity-limiter
        ! threshold. LOCKSTEP saGammaRetheta.F90 Source (gradient + rotating-
        ! frame blocks), minus everything that needs ReThetaTilde.
        !
        use blockPointers
        use constants
        use section
        use inputPhysics
        use flowVarRefState
        use inputIteration, only: transitionRefLength
        implicit none

        integer(kind=intType), intent(in) :: i, j, k
        real(kind=realType), intent(out) :: ss, strainMag, vortMag, dUds, velMag, vortLim
        real(kind=realType), intent(out) :: uux, uuy, uuz, vvx, vvy, vvz, wwx, wwy, wwz, fact

        real(kind=realType), parameter :: f23 = two * third
        real(kind=realType), parameter :: xminn = 1.e-10_realType
        real(kind=realType) :: div2, sxx, syy, szz, sxy, sxz, syz
        real(kind=realType) :: vortx, vorty, vortz, omegax, omegay, omegaz
        real(kind=realType) :: strainMag2, strainProd, vortProd
        real(kind=realType) :: xc(3), xxc(3), sc(3)
        real(kind=realType) :: velRelx, velRely, velRelz, velMag2, uRefTrans, refLenTrans
        real(kind=realType) :: uxhat, uyhat, uzhat

        omegax = timeRef * sections(sectionID)%rotRate(1)
        omegay = timeRef * sections(sectionID)%rotRate(2)
        omegaz = timeRef * sections(sectionID)%rotRate(3)

        uux = w(i + 1, j, k, ivx) * si(i, j, k, 1) - w(i - 1, j, k, ivx) * si(i - 1, j, k, 1) &
              + w(i, j + 1, k, ivx) * sj(i, j, k, 1) - w(i, j - 1, k, ivx) * sj(i, j - 1, k, 1) &
              + w(i, j, k + 1, ivx) * sk(i, j, k, 1) - w(i, j, k - 1, ivx) * sk(i, j, k - 1, 1)
        uuy = w(i + 1, j, k, ivx) * si(i, j, k, 2) - w(i - 1, j, k, ivx) * si(i - 1, j, k, 2) &
              + w(i, j + 1, k, ivx) * sj(i, j, k, 2) - w(i, j - 1, k, ivx) * sj(i, j - 1, k, 2) &
              + w(i, j, k + 1, ivx) * sk(i, j, k, 2) - w(i, j, k - 1, ivx) * sk(i, j, k - 1, 2)
        uuz = w(i + 1, j, k, ivx) * si(i, j, k, 3) - w(i - 1, j, k, ivx) * si(i - 1, j, k, 3) &
              + w(i, j + 1, k, ivx) * sj(i, j, k, 3) - w(i, j - 1, k, ivx) * sj(i, j - 1, k, 3) &
              + w(i, j, k + 1, ivx) * sk(i, j, k, 3) - w(i, j, k - 1, ivx) * sk(i, j, k - 1, 3)

        vvx = w(i + 1, j, k, ivy) * si(i, j, k, 1) - w(i - 1, j, k, ivy) * si(i - 1, j, k, 1) &
              + w(i, j + 1, k, ivy) * sj(i, j, k, 1) - w(i, j - 1, k, ivy) * sj(i, j - 1, k, 1) &
              + w(i, j, k + 1, ivy) * sk(i, j, k, 1) - w(i, j, k - 1, ivy) * sk(i, j, k - 1, 1)
        vvy = w(i + 1, j, k, ivy) * si(i, j, k, 2) - w(i - 1, j, k, ivy) * si(i - 1, j, k, 2) &
              + w(i, j + 1, k, ivy) * sj(i, j, k, 2) - w(i, j - 1, k, ivy) * sj(i, j - 1, k, 2) &
              + w(i, j, k + 1, ivy) * sk(i, j, k, 2) - w(i, j, k - 1, ivy) * sk(i, j, k - 1, 2)
        vvz = w(i + 1, j, k, ivy) * si(i, j, k, 3) - w(i - 1, j, k, ivy) * si(i - 1, j, k, 3) &
              + w(i, j + 1, k, ivy) * sj(i, j, k, 3) - w(i, j - 1, k, ivy) * sj(i, j - 1, k, 3) &
              + w(i, j, k + 1, ivy) * sk(i, j, k, 3) - w(i, j, k - 1, ivy) * sk(i, j, k - 1, 3)

        wwx = w(i + 1, j, k, ivz) * si(i, j, k, 1) - w(i - 1, j, k, ivz) * si(i - 1, j, k, 1) &
              + w(i, j + 1, k, ivz) * sj(i, j, k, 1) - w(i, j - 1, k, ivz) * sj(i, j - 1, k, 1) &
              + w(i, j, k + 1, ivz) * sk(i, j, k, 1) - w(i, j, k - 1, ivz) * sk(i, j, k - 1, 1)
        wwy = w(i + 1, j, k, ivz) * si(i, j, k, 2) - w(i - 1, j, k, ivz) * si(i - 1, j, k, 2) &
              + w(i, j + 1, k, ivz) * sj(i, j, k, 2) - w(i, j - 1, k, ivz) * sj(i, j - 1, k, 2) &
              + w(i, j, k + 1, ivz) * sk(i, j, k, 2) - w(i, j, k - 1, ivz) * sk(i, j, k - 1, 2)
        wwz = w(i + 1, j, k, ivz) * si(i, j, k, 3) - w(i - 1, j, k, ivz) * si(i - 1, j, k, 3) &
              + w(i, j + 1, k, ivz) * sj(i, j, k, 3) - w(i, j - 1, k, ivz) * sj(i, j - 1, k, 3) &
              + w(i, j, k + 1, ivz) * sk(i, j, k, 3) - w(i, j, k - 1, ivz) * sk(i, j, k - 1, 3)

        fact = fourth / vol(i, j, k)

        sxx = two * fact * uux
        syy = two * fact * vvy
        szz = two * fact * wwz
        sxy = fact * (uuy + vvx)
        sxz = fact * (uuz + wwx)
        syz = fact * (vvz + wwy)
        strainMag2 = two * (sxy**2 + sxz**2 + syz**2) + sxx**2 + syy**2 + szz**2

        vortx = two * fact * (wwy - vvz) - two * omegax
        vorty = two * fact * (uuz - wwx) - two * omegay
        vortz = two * fact * (vvx - uuy) - two * omegaz
        vortProd = vortx**2 + vorty**2 + vortz**2

        ! SA production magnitude per turbProd (strain or vorticity), as sa.F90
        if (turbProd .eq. strain) then
            div2 = f23 * (sxx + syy + szz)**2
            strainProd = two * strainMag2 - div2
            ss = sqrt(strainProd)
        else
            ss = sqrt(vortProd)
        end if

        vortMag = sqrt(max(vortProd, xminn))
        strainMag = sqrt(max(two * strainMag2, xminn))

        ! Relative (rotating-frame) velocity and streamwise gradient dU/ds
        xc(1) = eighth * (x(i - 1, j - 1, k - 1, 1) + x(i, j - 1, k - 1, 1) &
                          + x(i - 1, j, k - 1, 1) + x(i, j, k - 1, 1) + x(i - 1, j - 1, k, 1) &
                          + x(i, j - 1, k, 1) + x(i - 1, j, k, 1) + x(i, j, k, 1))
        xc(2) = eighth * (x(i - 1, j - 1, k - 1, 2) + x(i, j - 1, k - 1, 2) &
                          + x(i - 1, j, k - 1, 2) + x(i, j, k - 1, 2) + x(i - 1, j - 1, k, 2) &
                          + x(i, j - 1, k, 2) + x(i - 1, j, k, 2) + x(i, j, k, 2))
        xc(3) = eighth * (x(i - 1, j - 1, k - 1, 3) + x(i, j - 1, k - 1, 3) &
                          + x(i - 1, j, k - 1, 3) + x(i, j, k - 1, 3) + x(i - 1, j - 1, k, 3) &
                          + x(i, j - 1, k, 3) + x(i - 1, j, k, 3) + x(i, j, k, 3))
        xxc(1) = xc(1) - sections(sectionID)%rotCenter(1)
        xxc(2) = xc(2) - sections(sectionID)%rotCenter(2)
        xxc(3) = xc(3) - sections(sectionID)%rotCenter(3)
        sc(1) = omegay * xxc(3) - omegaz * xxc(2)
        sc(2) = omegaz * xxc(1) - omegax * xxc(3)
        sc(3) = omegax * xxc(2) - omegay * xxc(1)
        velRelx = w(i, j, k, ivx) - sc(1)
        velRely = w(i, j, k, ivy) - sc(2)
        velRelz = w(i, j, k, ivz) - sc(3)
        velMag2 = velRelx**2 + velRely**2 + velRelz**2
        velMag = sqrt(max(velMag2, xminn))

        uxhat = velRelx / max(velMag, xminn)
        uyhat = velRely / max(velMag, xminn)
        uzhat = velRelz / max(velMag, xminn)
        dUds = two * fact &
               * (uxhat * (uxhat * uux + uyhat * uuy + uzhat * uuz) &
                  + uyhat * (uxhat * vvx + uyhat * vvy + uzhat * vvz) &
                  + uzhat * (uxhat * wwx + uyhat * wwy + uzhat * wwz))

        ! P&Z vorticity-limiter threshold (Eqs. 52-53), only used with sgammaVortLimiter
        if (transitionRefLength > zero) then
            refLenTrans = transitionRefLength
        else
            refLenTrans = lengthRef
        end if
        uRefTrans = sqrt(uInf**2 + sc(1)**2 + sc(2)**2 + sc(3)**2)
        vortLim = uRefTrans * sqrt(max(uRefTrans / max(muInf * refLenTrans, xminn), xminn)) &
                  / 20.0_realType

    end subroutine cellKinematics

    subroutine sgSource
        !
        !  SA-sgamma source terms for all internal cells of the block, plus
        !  the 2x2 source Jacobian qq (LHS only, finite differences of
        !  cellSources, hidden from Tapenade).
        !
        use blockPointers
        use constants
        use paramTurb
        use inputPhysics
        use inputDiscretization, only: approxSA
        use flowVarRefState
        use inputIteration, only: transitionUseApproxSA, sgammaVortLimiter, sgammaOnsetTanh, &
                                  sgammaFturbLee, sgammaCoupleDestruction, sgammaFPGSmoothP, &
                                  sgammaCTU1, sgammaCTU2, sgammaCTU3
        implicit none

        integer(kind=intType), parameter :: dbgFonset = 1_intType
        integer(kind=intType), parameter :: dbgFonset1 = 2_intType
        integer(kind=intType), parameter :: dbgFlength = 3_intType
        integer(kind=intType), parameter :: dbgRturb = 4_intType
        integer(kind=intType), parameter :: dbgFPG = 5_intType
        integer(kind=intType), parameter :: dbgReV = 6_intType
        integer(kind=intType), parameter :: dbgReThetaC = 7_intType
        integer(kind=intType), parameter :: dbgReSOverCrit = 8_intType
        integer(kind=intType), parameter :: dbgStrainMag = 9_intType
        integer(kind=intType), parameter :: dbgGammaS = 10_intType
        integer(kind=intType), parameter :: dbgDestMult = 11_intType
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
        integer(kind=intType), parameter :: dbgLambdaTheta = 28_intType
        integer(kind=intType), parameter :: dbgGammaForSA = 30_intType
        integer(kind=intType), parameter :: dbgGammaLocal = 31_intType
        integer(kind=intType), parameter :: dbgVortMag = 32_intType
        integer(kind=intType), parameter :: dbgVortMagLim = 33_intType
        integer(kind=intType), parameter :: dbgFturb = 34_intType
        integer(kind=intType), parameter :: dbgSS = 35_intType
        integer(kind=intType), parameter :: dbgSST = 36_intType
        integer(kind=intType), parameter :: dbgFt2 = 37_intType
        integer(kind=intType), parameter :: dbgVelMag = 41_intType
        integer(kind=intType), parameter :: dbgDUds = 42_intType
        integer(kind=intType), parameter :: dbgNutSA = 43_intType
        integer(kind=intType), parameter :: dbgVortLim = 45_intType
        integer(kind=intType), parameter :: dbgQQ11 = 46_intType
        integer(kind=intType), parameter :: dbgQQ22 = 47_intType

        real(kind=realType), parameter :: xminn = 1.e-10_realType

        integer(kind=intType) :: i, j, k, ii
        real(kind=realType) :: ss, strainMag, vortMag, dUds, velMag, vortLim
        real(kind=realType) :: uux, uuy, uuz, vvx, vvy, vvz, wwx, wwy, wwz, fact
        real(kind=realType) :: nu, dist2Inv, chi, chi2, ft2, tuPct, yDist
        real(kind=realType) :: nuTilde, gam
        real(kind=realType) :: sNu, sGamma, gammaS, destMult, term2prod, term2dest
        real(kind=realType) :: fv1, fv2, sst, rr, gg, gg6, termFw, fwSa
        real(kind=realType) :: reVort, lamL, fPG, reThetaC, fOnset1, fOnset, fTurb
        real(kind=realType) :: pGamma, eGamma, gammaLocal, sEff, omEff, rTurb, nutSA
        logical :: approxTerm1, useLimiter, onsetTanh, fturbLee, coupleDest
        real(kind=realType) :: pSmooth, cTU1, cTU2, cTU3
        ! FD Jacobian scratch (LHS only)
        real(kind=realType) :: epsNu, epsG, sNuP, sGamP
        real(kind=realType) :: d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14
        real(kind=realType) :: d15, d16, d17, d18, d19, d20, d21, d22, d23, d24, d25, d26

        ! Options -> locals (copied once, before the differentiated loop)
        approxTerm1 = approxSA .and. transitionUseApproxSA
        useLimiter = sgammaVortLimiter
        onsetTanh = sgammaOnsetTanh
        fturbLee = sgammaFturbLee
        coupleDest = sgammaCoupleDestruction
        pSmooth = sgammaFPGSmoothP
        cTU1 = sgammaCTU1
        cTU2 = sgammaCTU2
        cTU3 = sgammaCTU3
        tuPct = turbIntensityInf * 100.0_realType

#ifndef USE_TAPENADE
        qq = zero
#endif

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
                        call cellKinematics(i, j, k, ss, strainMag, vortMag, dUds, velMag, vortLim, &
                                            uux, uuy, uuz, vvx, vvy, vvz, wwx, wwy, wwz, fact)

                        nu = rlv(i, j, k) / w(i, j, k, irho)
                        dist2Inv = one / (d2Wall(i, j, k)**2)
                        yDist = d2Wall(i, j, k)
                        nuTilde = w(i, j, k, itu1)
                        gam = w(i, j, k, itu2)
                        chi = nuTilde / nu
                        chi2 = chi * chi
                        if (useft2SA) then
                            ft2 = rsaCt3 * exp(-rsaCt4 * chi2)
                        else
                            ft2 = zero
                        end if

                        call cellSources(nuTilde, gam, nu, dist2Inv, ss, strainMag, vortMag, &
                                         yDist, dUds, w(i, j, k, irho), rlv(i, j, k), tuPct, vortLim, &
                                         ft2, approxTerm1, useLimiter, onsetTanh, fturbLee, coupleDest, &
                                         pSmooth, cTU1, cTU2, cTU3, &
                                         sNu, sGamma, gammaS, destMult, term2prod, term2dest, &
                                         fv1, fv2, sst, rr, gg, gg6, termFw, fwSa, &
                                         reVort, lamL, fPG, reThetaC, fOnset1, fOnset, fTurb, &
                                         pGamma, eGamma, gammaLocal, sEff, omEff, rTurb, nutSA)

                        scratch(i, j, k, idvt) = sNu
                        scratch(i, j, k, idvt + 1) = sGamma

                        if (associated(transitionDebug)) then
                            transitionDebug(i, j, k, dbgFonset) = fOnset
                            transitionDebug(i, j, k, dbgFonset1) = fOnset1
                            transitionDebug(i, j, k, dbgFlength) = sgFlength
                            transitionDebug(i, j, k, dbgRturb) = rTurb
                            transitionDebug(i, j, k, dbgFPG) = fPG
                            transitionDebug(i, j, k, dbgReV) = reVort
                            transitionDebug(i, j, k, dbgReThetaC) = reThetaC
                            transitionDebug(i, j, k, dbgReSOverCrit) = reVort / (sgFonsetC * reThetaC)
                            transitionDebug(i, j, k, dbgStrainMag) = strainMag
                            transitionDebug(i, j, k, dbgGammaS) = gammaS
                            transitionDebug(i, j, k, dbgDestMult) = destMult
                            transitionDebug(i, j, k, dbgDudx) = two * fact * uux
                            transitionDebug(i, j, k, dbgDudy) = two * fact * uuy
                            transitionDebug(i, j, k, dbgDudz) = two * fact * uuz
                            transitionDebug(i, j, k, dbgDvdx) = two * fact * vvx
                            transitionDebug(i, j, k, dbgDvdy) = two * fact * vvy
                            transitionDebug(i, j, k, dbgDvdz) = two * fact * vvz
                            transitionDebug(i, j, k, dbgDwdx) = two * fact * wwx
                            transitionDebug(i, j, k, dbgDwdy) = two * fact * wwy
                            transitionDebug(i, j, k, dbgDwdz) = two * fact * wwz
                            transitionDebug(i, j, k, dbgGamma) = gam
                            transitionDebug(i, j, k, dbgWallDist) = yDist
                            transitionDebug(i, j, k, dbgRho) = w(i, j, k, irho)
                            transitionDebug(i, j, k, dbgMu) = rlv(i, j, k)
                            transitionDebug(i, j, k, dbgGammaProd) = pGamma
                            transitionDebug(i, j, k, dbgGammaDest) = eGamma
                            transitionDebug(i, j, k, dbgLambdaTheta) = lamL
                            transitionDebug(i, j, k, dbgGammaForSA) = gammaS
                            transitionDebug(i, j, k, dbgGammaLocal) = gammaLocal
                            transitionDebug(i, j, k, dbgVortMag) = vortMag
                            transitionDebug(i, j, k, dbgVortMagLim) = omEff
                            transitionDebug(i, j, k, dbgFturb) = fTurb
                            transitionDebug(i, j, k, dbgSS) = ss
                            transitionDebug(i, j, k, dbgSST) = sst
                            transitionDebug(i, j, k, dbgFt2) = ft2
                            transitionDebug(i, j, k, dbgVelMag) = velMag
                            transitionDebug(i, j, k, dbgDUds) = dUds
                            transitionDebug(i, j, k, dbgNutSA) = nutSA
                            transitionDebug(i, j, k, dbgVortLim) = vortLim
                        end if

#ifndef USE_TAPENADE
                        ! ---- 2x2 source Jacobian, -dS/dQ, by one-sided finite
                        ! differences of cellSources (LHS only; the residual is
                        ! the analytic call above). Diagonals clipped >= 0 as in
                        ! sa.F90 / saGammaRetheta.F90 (implicit treatment only for
                        ! terms that add diagonal dominance).
                        epsNu = max(1.0e-6_realType * abs(nuTilde), 1.0e-12_realType * nu)
                        epsG = 1.0e-6_realType
                        call cellSources(nuTilde + epsNu, gam, nu, dist2Inv, ss, strainMag, vortMag, &
                                         yDist, dUds, w(i, j, k, irho), rlv(i, j, k), tuPct, vortLim, &
                                         ft2, approxTerm1, useLimiter, onsetTanh, fturbLee, coupleDest, &
                                         pSmooth, cTU1, cTU2, cTU3, &
                                         sNuP, sGamP, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, &
                                         d12, d13, d14, d15, d16, d17, d18, d19, d20, d21, d22, d23, d24, d25, d26)
                        qq(i, j, k, 1, 1) = -(sNuP - sNu) / epsNu
                        qq(i, j, k, 2, 1) = -(sGamP - sGamma) / epsNu
                        call cellSources(nuTilde, gam + epsG, nu, dist2Inv, ss, strainMag, vortMag, &
                                         yDist, dUds, w(i, j, k, irho), rlv(i, j, k), tuPct, vortLim, &
                                         ft2, approxTerm1, useLimiter, onsetTanh, fturbLee, coupleDest, &
                                         pSmooth, cTU1, cTU2, cTU3, &
                                         sNuP, sGamP, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, &
                                         d12, d13, d14, d15, d16, d17, d18, d19, d20, d21, d22, d23, d24, d25, d26)
                        qq(i, j, k, 1, 2) = -(sNuP - sNu) / epsG
                        qq(i, j, k, 2, 2) = -(sGamP - sGamma) / epsG

                        qq(i, j, k, 1, 1) = max(qq(i, j, k, 1, 1), zero)
                        qq(i, j, k, 2, 2) = max(qq(i, j, k, 2, 2), zero)

                        if (associated(transitionDebug)) then
                            transitionDebug(i, j, k, dbgQQ11) = qq(i, j, k, 1, 1)
                            transitionDebug(i, j, k, dbgQQ22) = qq(i, j, k, 2, 2)
                        end if
#endif
#ifdef TAPENADE_REVERSE
                    end do
#else
                end do
            end do
        end do
#endif

    end subroutine sgSource
    subroutine sgViscous
        !
        !  Viscous (diffusion) terms of nuTilde and gam. LOCKSTEP
        !  saGammaRetheta.F90 Viscous minus the ReThetaTilde equation.
        !
        use constants
        use blockPointers
        use paramTurb
        implicit none

        integer(kind=intType) :: i, j, k, ii
        real(kind=realType) :: nu, nu_m, nu_p, nut, nut_m, nut_p, nu_tm, nu_tp
        real(kind=realType) :: nuTilde, nuTilde_m, nuTilde_p
        real(kind=realType) :: chi, chi3, chi_m, chi3_m, chi_p, chi3_p, fv1, fv1_m, fv1_p
        real(kind=realType) :: voli, volmi, volpi, xm, ym, zm, xp, yp, zp, xa, ya, za, ttm, ttp
        real(kind=realType) :: num, nup, cdm, cdp, cdm_gamma, cdp_gamma
        real(kind=realType) :: cnud, cam, cap, nutm, nutp
        real(kind=realType) :: c1m, c1p, c10, c2m, c2p, c20
        real(kind=realType) :: b1, c1, d1, b2, c2, d2
        real(kind=realType) :: cb3Inv, cv13

        cb3Inv = one / rsaCb3
        cv13 = rsaCv1**3

        ! ---------------- k-direction ----------------
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

                        num = half * (nu_m + nu)
                        nup = half * (nu_p + nu)
                        nu_tm = half * (nut_m + nut)
                        nu_tp = half * (nut_p + nut)

                        cdm = (num + (one + rsaCb2) * nutm) * ttm * cb3Inv
                        cdp = (nup + (one + rsaCb2) * nutp) * ttp * cb3Inv
                        cdm_gamma = (num + nu_tm / sigmaF) * ttm
                        cdp_gamma = (nup + nu_tp / sigmaF) * ttp

                        c1m = max(cdm + cam, zero)
                        c1p = max(cdp + cap, zero)
                        c10 = c1m + c1p
                        c2m = max(cdm_gamma, zero)
                        c2p = max(cdp_gamma, zero)
                        c20 = c2m + c2p

                        scratch(i, j, k, idvt) = scratch(i, j, k, idvt) + c1m * w(i, j, k - 1, itu1) &
                                                 - c10 * w(i, j, k, itu1) + c1p * w(i, j, k + 1, itu1)
                        scratch(i, j, k, idvt + 1) = scratch(i, j, k, idvt + 1) + c2m * w(i, j, k - 1, itu2) &
                                                     - c20 * w(i, j, k, itu2) + c2p * w(i, j, k + 1, itu2)
#ifndef USE_TAPENADE
                        b1 = -c1m; c1 = c10; d1 = -c1p
                        if (k == 2) then
                            qq(i, j, k, 1, 1) = qq(i, j, k, 1, 1) + c1 - b1 * max(bmtk1(i, j, itu1, itu1), zero)
                            qq(i, j, k, 1, 2) = qq(i, j, k, 1, 2) - b1 * bmtk1(i, j, itu1, itu2)
                        else if (k == kl) then
                            qq(i, j, k, 1, 1) = qq(i, j, k, 1, 1) + c1 - d1 * max(bmtk2(i, j, itu1, itu1), zero)
                            qq(i, j, k, 1, 2) = qq(i, j, k, 1, 2) - d1 * bmtk2(i, j, itu1, itu2)
                        else
                            qq(i, j, k, 1, 1) = qq(i, j, k, 1, 1) + c1
                        end if
                        b2 = -c2m; c2 = c20; d2 = -c2p
                        if (k == 2) then
                            qq(i, j, k, 2, 1) = qq(i, j, k, 2, 1) - b2 * bmtk1(i, j, itu2, itu1)
                            qq(i, j, k, 2, 2) = qq(i, j, k, 2, 2) + c2 - b2 * bmtk1(i, j, itu2, itu2)
                        else if (k == kl) then
                            qq(i, j, k, 2, 1) = qq(i, j, k, 2, 1) - d2 * bmtk2(i, j, itu2, itu1)
                            qq(i, j, k, 2, 2) = qq(i, j, k, 2, 2) + c2 - d2 * bmtk2(i, j, itu2, itu2)
                        else
                            qq(i, j, k, 2, 2) = qq(i, j, k, 2, 2) + c2
                        end if
#endif
#ifdef TAPENADE_REVERSE
                    end do
#else
                end do
            end do
        end do
#endif

        ! ---------------- j-direction ----------------
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

                        num = half * (nu_m + nu)
                        nup = half * (nu_p + nu)
                        nu_tm = half * (nut_m + nut)
                        nu_tp = half * (nut_p + nut)

                        cdm = (num + (one + rsaCb2) * nutm) * ttm * cb3Inv
                        cdp = (nup + (one + rsaCb2) * nutp) * ttp * cb3Inv
                        cdm_gamma = (num + nu_tm / sigmaF) * ttm
                        cdp_gamma = (nup + nu_tp / sigmaF) * ttp

                        c1m = max(cdm + cam, zero)
                        c1p = max(cdp + cap, zero)
                        c10 = c1m + c1p
                        c2m = max(cdm_gamma, zero)
                        c2p = max(cdp_gamma, zero)
                        c20 = c2m + c2p

                        scratch(i, j, k, idvt) = scratch(i, j, k, idvt) + c1m * w(i, j - 1, k, itu1) &
                                                 - c10 * w(i, j, k, itu1) + c1p * w(i, j + 1, k, itu1)
                        scratch(i, j, k, idvt + 1) = scratch(i, j, k, idvt + 1) + c2m * w(i, j - 1, k, itu2) &
                                                     - c20 * w(i, j, k, itu2) + c2p * w(i, j + 1, k, itu2)
#ifndef USE_TAPENADE
                        b1 = -c1m; c1 = c10; d1 = -c1p
                        if (j == 2) then
                            qq(i, j, k, 1, 1) = qq(i, j, k, 1, 1) + c1 - b1 * max(bmtj1(i, k, itu1, itu1), zero)
                            qq(i, j, k, 1, 2) = qq(i, j, k, 1, 2) - b1 * bmtj1(i, k, itu1, itu2)
                        else if (j == jl) then
                            qq(i, j, k, 1, 1) = qq(i, j, k, 1, 1) + c1 - d1 * max(bmtj2(i, k, itu1, itu1), zero)
                            qq(i, j, k, 1, 2) = qq(i, j, k, 1, 2) - d1 * bmtj2(i, k, itu1, itu2)
                        else
                            qq(i, j, k, 1, 1) = qq(i, j, k, 1, 1) + c1
                        end if
                        b2 = -c2m; c2 = c20; d2 = -c2p
                        if (j == 2) then
                            qq(i, j, k, 2, 1) = qq(i, j, k, 2, 1) - b2 * bmtj1(i, k, itu2, itu1)
                            qq(i, j, k, 2, 2) = qq(i, j, k, 2, 2) + c2 - b2 * bmtj1(i, k, itu2, itu2)
                        else if (j == jl) then
                            qq(i, j, k, 2, 1) = qq(i, j, k, 2, 1) - d2 * bmtj2(i, k, itu2, itu1)
                            qq(i, j, k, 2, 2) = qq(i, j, k, 2, 2) + c2 - d2 * bmtj2(i, k, itu2, itu2)
                        else
                            qq(i, j, k, 2, 2) = qq(i, j, k, 2, 2) + c2
                        end if
#endif
#ifdef TAPENADE_REVERSE
                    end do
#else
                end do
            end do
        end do
#endif

        ! ---------------- i-direction ----------------
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

                        num = half * (nu_m + nu)
                        nup = half * (nu_p + nu)
                        nu_tm = half * (nut_m + nut)
                        nu_tp = half * (nut_p + nut)

                        cdm = (num + (one + rsaCb2) * nutm) * ttm * cb3Inv
                        cdp = (nup + (one + rsaCb2) * nutp) * ttp * cb3Inv
                        cdm_gamma = (num + nu_tm / sigmaF) * ttm
                        cdp_gamma = (nup + nu_tp / sigmaF) * ttp

                        c1m = max(cdm + cam, zero)
                        c1p = max(cdp + cap, zero)
                        c10 = c1m + c1p
                        c2m = max(cdm_gamma, zero)
                        c2p = max(cdp_gamma, zero)
                        c20 = c2m + c2p

                        scratch(i, j, k, idvt) = scratch(i, j, k, idvt) + c1m * w(i - 1, j, k, itu1) &
                                                 - c10 * w(i, j, k, itu1) + c1p * w(i + 1, j, k, itu1)
                        scratch(i, j, k, idvt + 1) = scratch(i, j, k, idvt + 1) + c2m * w(i - 1, j, k, itu2) &
                                                     - c20 * w(i, j, k, itu2) + c2p * w(i + 1, j, k, itu2)
#ifndef USE_TAPENADE
                        b1 = -c1m; c1 = c10; d1 = -c1p
                        if (i == 2) then
                            qq(i, j, k, 1, 1) = qq(i, j, k, 1, 1) + c1 - b1 * max(bmti1(j, k, itu1, itu1), zero)
                            qq(i, j, k, 1, 2) = qq(i, j, k, 1, 2) - b1 * bmti1(j, k, itu1, itu2)
                        else if (i == il) then
                            qq(i, j, k, 1, 1) = qq(i, j, k, 1, 1) + c1 - d1 * max(bmti2(j, k, itu1, itu1), zero)
                            qq(i, j, k, 1, 2) = qq(i, j, k, 1, 2) - d1 * bmti2(j, k, itu1, itu2)
                        else
                            qq(i, j, k, 1, 1) = qq(i, j, k, 1, 1) + c1
                        end if
                        b2 = -c2m; c2 = c20; d2 = -c2p
                        if (i == 2) then
                            qq(i, j, k, 2, 1) = qq(i, j, k, 2, 1) - b2 * bmti1(j, k, itu2, itu1)
                            qq(i, j, k, 2, 2) = qq(i, j, k, 2, 2) + c2 - b2 * bmti1(j, k, itu2, itu2)
                        else if (i == il) then
                            qq(i, j, k, 2, 1) = qq(i, j, k, 2, 1) - d2 * bmti2(j, k, itu2, itu1)
                            qq(i, j, k, 2, 2) = qq(i, j, k, 2, 2) + c2 - d2 * bmti2(j, k, itu2, itu2)
                        else
                            qq(i, j, k, 2, 2) = qq(i, j, k, 2, 2) + c2
                        end if
#endif
#ifdef TAPENADE_REVERSE
                    end do
#else
                end do
            end do
        end do
#endif

    end subroutine sgViscous
    subroutine sgResScale
        use constants
        use blockPointers
        implicit none

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
                        dw(i, j, k, itu1) = -volRef(i, j, k) * scratch(i, j, k, idvt) * rblank
                        dw(i, j, k, itu2) = -volRef(i, j, k) * scratch(i, j, k, idvt + 1) * rblank
#ifdef TAPENADE_REVERSE
                    end do
#else
                end do
            end do
        end do
#endif
    end subroutine sgResScale
    subroutine saGammaSolve(resOnly)
        !
        !  Coupled DD-ADI solve of the 2x2 (nuTilde, gam) system. LOCKSTEP
        !  saGammaReThetaSolve minus the ReThetaTilde row/column, tdia2x2.
        !
        use blockPointers
        use constants
        use flowVarRefState
        use inputIteration
        use inputPhysics
        use paramTurb
        use turbUtils, only: tdia2x2
        implicit none

        logical, intent(in) :: resOnly

        integer(kind=intType) :: i, j, k, mm
        real(kind=realType) :: nu, nu_m, nu_p, nut, nut_m, nut_p, nu_tm, nu_tp
        real(kind=realType) :: nuTilde, nuTilde_m, nuTilde_p
        real(kind=realType) :: chi, chi3, fv1, chi_m, chi3_m, fv1_m, chi_p, chi3_p, fv1_p
        real(kind=realType) :: voli, volmi, volpi, xm, ym, zm, xp, yp, zp, xa, ya, za, ttm, ttp
        real(kind=realType) :: num_v, nup_v, cdm, cdp, cdm_gamma, cdp_gamma
        real(kind=realType) :: cnud, cam, cap, nutm, nutp
        real(kind=realType) :: c1m, c1p, c2m, c2p
        real(kind=realType) :: qs, uu, um, up
        real(kind=realType) :: cb3Inv, cv13, rblank, factor
        real(kind=realType) :: gammaNew, gammaDelta, dampFactor
        integer(kind=intType) :: nDampCapGamma
        real(kind=realType) :: scaleNu, scaleGamma, s12, s21
        real(kind=realType), dimension(2, 2:max(il, jl, kl)) :: bb, dd, ff
        real(kind=realType), dimension(2, 2, 2:max(il, jl, kl)) :: cc
        logical, save :: printedCoupling = .false.

        if (resOnly) return

        if (.not. printedCoupling) then
            printedCoupling = .true.
            if (TurbDADICoupled == 0) then
                print *, 'SA-sgamma DADI coupling: decoupled (diagonal only)'
            else
                print *, 'SA-sgamma DADI coupling: fully coupled (2x2 block)'
            end if
        end if

        cb3Inv = one / rsaCb3
        cv13 = rsaCv1**3

        scaleNu = max(abs(turbResScale(1)), 1.0e-12_realType)
        scaleGamma = max(abs(turbResScale(2)), 1.0e-12_realType)
        s12 = scaleGamma / scaleNu
        s21 = scaleNu / scaleGamma

        factor = one
        if (turbRelax == turbRelaxImplicit) factor = one + (one - alfaTurb) / alfaTurb

        do k = 2, kl
            do j = 2, jl
                do i = 2, il
                    if (TurbDADICoupled == 0) then
                        qq(i, j, k, 1, 2) = zero
                        qq(i, j, k, 2, 1) = zero
                    end if
                    ! Source dt restriction (P&Z Eq. 59), additive form as in the GR DADI
                    if (transitionSrcDtRestrict) then
                        qq(i, j, k, 1, 1) = qq(i, j, k, 1, 1) + srcLambda(i, j, k, 1) / transitionSrcDtLimit
                        qq(i, j, k, 2, 2) = qq(i, j, k, 2, 2) + srcLambda(i, j, k, 2) / transitionSrcDtLimit
                    end if
                    qq(i, j, k, 1, 2) = qq(i, j, k, 1, 2) * s12
                    qq(i, j, k, 2, 1) = qq(i, j, k, 2, 1) * s21
                    qq(i, j, k, 1, 1) = factor * qq(i, j, k, 1, 1)
                    qq(i, j, k, 1, 2) = factor * qq(i, j, k, 1, 2)
                    qq(i, j, k, 2, 1) = factor * qq(i, j, k, 2, 1)
                    qq(i, j, k, 2, 2) = factor * qq(i, j, k, 2, 2)
                    scratch(i, j, k, idvt) = scratch(i, j, k, idvt) / scaleNu
                    scratch(i, j, k, idvt + 1) = scratch(i, j, k, idvt + 1) / scaleGamma
                end do
            end do
        end do

        qs = zero

        ! ---------------- j-direction ----------------
        do k = 2, kl
            do i = 2, il
                do j = 2, jl
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

                    cdm = (num_v + (one + rsaCb2) * nutm) * ttm * cb3Inv
                    cdp = (nup_v + (one + rsaCb2) * nutp) * ttp * cb3Inv
                    cdm_gamma = (num_v + nu_tm / sigmaF) * ttm
                    cdp_gamma = (nup_v + nu_tp / sigmaF) * ttp

                    c1m = max(cdm + cam, zero)
                    c1p = max(cdp + cap, zero)
                    c2m = max(cdm_gamma, zero)
                    c2p = max(cdp_gamma, zero)

                    bb(1, j) = -c1m; dd(1, j) = -c1p
                    bb(2, j) = -c2m; dd(2, j) = -c2p

                    if (addGridVelocities) qs = half * (sFaceJ(i, j, k) + sFaceJ(i, j - 1, k)) * voli
                    uu = xa * w(i, j, k, ivx) + ya * w(i, j, k, ivy) + za * w(i, j, k, ivz) - qs
                    um = zero; up = zero
                    if (uu < zero) um = uu
                    if (uu > zero) up = uu
                    bb(1, j) = bb(1, j) - up; dd(1, j) = dd(1, j) + um
                    bb(2, j) = bb(2, j) - up; dd(2, j) = dd(2, j) + um

                    rblank = real(iblank(i, j, k), realType)
                    cc(1, 1, j) = qq(i, j, k, 1, 1)
                    cc(1, 2, j) = qq(i, j, k, 1, 2) * rblank
                    cc(2, 1, j) = qq(i, j, k, 2, 1) * rblank
                    cc(2, 2, j) = qq(i, j, k, 2, 2)
                    ff(1, j) = scratch(i, j, k, idvt) * rblank
                    ff(2, j) = scratch(i, j, k, idvt + 1) * rblank
                    bb(:, j) = bb(:, j) * rblank
                    dd(:, j) = dd(:, j) * rblank
                end do

                call tdia2x2(2_intType, jl, bb, cc, dd, ff)

                do j = 2, jl
                    scratch(i, j, k, idvt) = qq(i, j, k, 1, 1) * ff(1, j) + qq(i, j, k, 1, 2) * ff(2, j)
                    scratch(i, j, k, idvt + 1) = qq(i, j, k, 2, 1) * ff(1, j) + qq(i, j, k, 2, 2) * ff(2, j)
                end do
            end do
        end do

        ! ---------------- i-direction ----------------
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

                    cdm = (num_v + (one + rsaCb2) * nutm) * ttm * cb3Inv
                    cdp = (nup_v + (one + rsaCb2) * nutp) * ttp * cb3Inv
                    cdm_gamma = (num_v + nu_tm / sigmaF) * ttm
                    cdp_gamma = (nup_v + nu_tp / sigmaF) * ttp

                    c1m = max(cdm + cam, zero)
                    c1p = max(cdp + cap, zero)
                    c2m = max(cdm_gamma, zero)
                    c2p = max(cdp_gamma, zero)

                    bb(1, i) = -c1m; dd(1, i) = -c1p
                    bb(2, i) = -c2m; dd(2, i) = -c2p

                    if (addGridVelocities) qs = half * (sFaceI(i, j, k) + sFaceI(i - 1, j, k)) * voli
                    uu = xa * w(i, j, k, ivx) + ya * w(i, j, k, ivy) + za * w(i, j, k, ivz) - qs
                    um = zero; up = zero
                    if (uu < zero) um = uu
                    if (uu > zero) up = uu
                    bb(1, i) = bb(1, i) - up; dd(1, i) = dd(1, i) + um
                    bb(2, i) = bb(2, i) - up; dd(2, i) = dd(2, i) + um

                    rblank = real(iblank(i, j, k), realType)
                    cc(1, 1, i) = qq(i, j, k, 1, 1)
                    cc(1, 2, i) = qq(i, j, k, 1, 2) * rblank
                    cc(2, 1, i) = qq(i, j, k, 2, 1) * rblank
                    cc(2, 2, i) = qq(i, j, k, 2, 2)
                    ff(1, i) = scratch(i, j, k, idvt) * rblank
                    ff(2, i) = scratch(i, j, k, idvt + 1) * rblank
                    bb(:, i) = bb(:, i) * rblank
                    dd(:, i) = dd(:, i) * rblank
                end do

                call tdia2x2(2_intType, il, bb, cc, dd, ff)

                do i = 2, il
                    scratch(i, j, k, idvt) = qq(i, j, k, 1, 1) * ff(1, i) + qq(i, j, k, 1, 2) * ff(2, i)
                    scratch(i, j, k, idvt + 1) = qq(i, j, k, 2, 1) * ff(1, i) + qq(i, j, k, 2, 2) * ff(2, i)
                end do
            end do
        end do

        ! ---------------- k-direction ----------------
        do j = 2, jl
            do i = 2, il
                do k = 2, kl
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

                    cdm = (num_v + (one + rsaCb2) * nutm) * ttm * cb3Inv
                    cdp = (nup_v + (one + rsaCb2) * nutp) * ttp * cb3Inv
                    cdm_gamma = (num_v + nu_tm / sigmaF) * ttm
                    cdp_gamma = (nup_v + nu_tp / sigmaF) * ttp

                    c1m = max(cdm + cam, zero)
                    c1p = max(cdp + cap, zero)
                    c2m = max(cdm_gamma, zero)
                    c2p = max(cdp_gamma, zero)

                    bb(1, k) = -c1m; dd(1, k) = -c1p
                    bb(2, k) = -c2m; dd(2, k) = -c2p

                    if (addGridVelocities) qs = half * (sFaceK(i, j, k) + sFaceK(i, j, k - 1)) * voli
                    uu = xa * w(i, j, k, ivx) + ya * w(i, j, k, ivy) + za * w(i, j, k, ivz) - qs
                    um = zero; up = zero
                    if (uu < zero) um = uu
                    if (uu > zero) up = uu
                    bb(1, k) = bb(1, k) - up; dd(1, k) = dd(1, k) + um
                    bb(2, k) = bb(2, k) - up; dd(2, k) = dd(2, k) + um

                    rblank = real(iblank(i, j, k), realType)
                    cc(1, 1, k) = qq(i, j, k, 1, 1)
                    cc(1, 2, k) = qq(i, j, k, 1, 2) * rblank
                    cc(2, 1, k) = qq(i, j, k, 2, 1) * rblank
                    cc(2, 2, k) = qq(i, j, k, 2, 2)
                    ff(1, k) = scratch(i, j, k, idvt) * rblank
                    ff(2, k) = scratch(i, j, k, idvt + 1) * rblank
                    bb(:, k) = bb(:, k) * rblank
                    dd(:, k) = dd(:, k) * rblank
                end do

                call tdia2x2(2_intType, kl, bb, cc, dd, ff)

                do k = 2, kl
                    scratch(i, j, k, idvt) = ff(1, k)
                    scratch(i, j, k, idvt + 1) = ff(2, k)
                end do
            end do
        end do

        ! ---------------- update (Algorithm 2 damping on gam) ----------------
        factor = one
        if (turbRelax == turbRelaxExplicit) factor = alfaTurb
        nDampCapGamma = 0

        do k = 2, kl
            do j = 2, jl
                do i = 2, il
                    w(i, j, k, itu1) = w(i, j, k, itu1) + factor * scaleNu * scratch(i, j, k, idvt)
                    w(i, j, k, itu1) = max(w(i, j, k, itu1), zero)

                    gammaDelta = factor * scaleGamma * scratch(i, j, k, idvt + 1)
                    gammaNew = w(i, j, k, itu2) + gammaDelta
                    dampFactor = one
                    do mm = 1, transitionDampMaxIter
                        if (gammaNew >= rsaGRgammaLo .and. gammaNew <= rsaGRgammaHi) exit
                        dampFactor = dampFactor * transitionDampTheta
                        gammaNew = w(i, j, k, itu2) + dampFactor * gammaDelta
                    end do
                    if (gammaNew < rsaGRgammaLo .or. gammaNew > rsaGRgammaHi) then
                        nDampCapGamma = nDampCapGamma + 1
                        gammaNew = min(max(gammaNew, rsaGRgammaLo), rsaGRgammaHi)
                    end if
                    w(i, j, k, itu2) = gammaNew
                end do
            end do
        end do

        if (nDampCapGamma > 0) then
            print *, 'Warning: SA-sgamma update damping exhausted transitionDampMaxIter (', &
                transitionDampMaxIter, ') in ', nDampCapGamma, ' gam cells; values clipped.'
        end if

    end subroutine saGammaSolve

    subroutine evalSrcJacBlockSaGamma(i, j, k, A)
        !
        ! Source Jacobian A = dS/dQ (2x2) for cell (i,j,k) by finite
        ! differences of cellSources. Independent of qq. Used by
        ! computeSrcLambdaSaGamma (P&Z Eq. 59).
        !
        use blockPointers
        use constants
        use paramTurb
        use inputPhysics
        use flowVarRefState
        use inputDiscretization, only: approxSA
        use inputIteration, only: transitionUseApproxSA, sgammaVortLimiter, sgammaOnsetTanh, &
                                  sgammaFturbLee, sgammaCoupleDestruction, sgammaFPGSmoothP, &
                                  sgammaCTU1, sgammaCTU2, sgammaCTU3
        implicit none

        integer(kind=intType), intent(in) :: i, j, k
        real(kind=realType), intent(out) :: A(2, 2)

        real(kind=realType) :: ss, strainMag, vortMag, dUds, velMag, vortLim
        real(kind=realType) :: uux, uuy, uuz, vvx, vvy, vvz, wwx, wwy, wwz, fact
        real(kind=realType) :: nu, dist2Inv, chi, chi2, ft2, tuPct, yDist, nuTilde, gam
        real(kind=realType) :: sNu, sGamma, sNuP, sGamP, epsNu, epsG
        real(kind=realType) :: d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14
        real(kind=realType) :: d15, d16, d17, d18, d19, d20, d21, d22, d23, d24, d25, d26
        logical :: approxTerm1

        approxTerm1 = approxSA .and. transitionUseApproxSA
        tuPct = turbIntensityInf * 100.0_realType

        call cellKinematics(i, j, k, ss, strainMag, vortMag, dUds, velMag, vortLim, &
                            uux, uuy, uuz, vvx, vvy, vvz, wwx, wwy, wwz, fact)
        nu = rlv(i, j, k) / w(i, j, k, irho)
        dist2Inv = one / (d2Wall(i, j, k)**2)
        yDist = d2Wall(i, j, k)
        nuTilde = w(i, j, k, itu1)
        gam = w(i, j, k, itu2)
        chi = nuTilde / nu
        chi2 = chi * chi
        if (useft2SA) then
            ft2 = rsaCt3 * exp(-rsaCt4 * chi2)
        else
            ft2 = zero
        end if

        call cellSources(nuTilde, gam, nu, dist2Inv, ss, strainMag, vortMag, yDist, dUds, &
                         w(i, j, k, irho), rlv(i, j, k), tuPct, vortLim, ft2, approxTerm1, &
                         sgammaVortLimiter, sgammaOnsetTanh, sgammaFturbLee, sgammaCoupleDestruction, &
                         sgammaFPGSmoothP, sgammaCTU1, sgammaCTU2, sgammaCTU3, &
                         sNu, sGamma, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, &
                         d12, d13, d14, d15, d16, d17, d18, d19, d20, d21, d22, d23, d24, d25, d26)
        epsNu = max(1.0e-6_realType * abs(nuTilde), 1.0e-12_realType * nu)
        epsG = 1.0e-6_realType
        call cellSources(nuTilde + epsNu, gam, nu, dist2Inv, ss, strainMag, vortMag, yDist, dUds, &
                         w(i, j, k, irho), rlv(i, j, k), tuPct, vortLim, ft2, approxTerm1, &
                         sgammaVortLimiter, sgammaOnsetTanh, sgammaFturbLee, sgammaCoupleDestruction, &
                         sgammaFPGSmoothP, sgammaCTU1, sgammaCTU2, sgammaCTU3, &
                         sNuP, sGamP, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, &
                         d12, d13, d14, d15, d16, d17, d18, d19, d20, d21, d22, d23, d24, d25, d26)
        A(1, 1) = (sNuP - sNu) / epsNu
        A(2, 1) = (sGamP - sGamma) / epsNu
        call cellSources(nuTilde, gam + epsG, nu, dist2Inv, ss, strainMag, vortMag, yDist, dUds, &
                         w(i, j, k, irho), rlv(i, j, k), tuPct, vortLim, ft2, approxTerm1, &
                         sgammaVortLimiter, sgammaOnsetTanh, sgammaFturbLee, sgammaCoupleDestruction, &
                         sgammaFPGSmoothP, sgammaCTU1, sgammaCTU2, sgammaCTU3, &
                         sNuP, sGamP, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, &
                         d12, d13, d14, d15, d16, d17, d18, d19, d20, d21, d22, d23, d24, d25, d26)
        A(1, 2) = (sNuP - sNu) / epsG
        A(2, 2) = (sGamP - sGamma) / epsG

    end subroutine evalSrcJacBlockSaGamma

    subroutine computeSrcLambdaSaGamma(mode)
        !
        ! srcLambda(i,j,k,1:2) for the SA-sgamma 2x2 source Jacobian (P&Z Eq. 59).
        ! Same modes as saGammaRetheta%computeSrcLambda; srcLambda(:,:,:,3)
        ! (allocated for both models) is set to zero and never read for nwt=2.
        !
        use constants
        use blockPointers, only: il, jl, kl, srcLambda
        use paramTurb, only: srcLambdaModeDecoupled
        implicit none

        integer(kind=intType), intent(in) :: mode
        integer(kind=intType) :: i, j, k
        real(kind=realType) :: A(2, 2), tr2, disc2, lambda2x2

        do k = 2, kl
            do j = 2, jl
                do i = 2, il
                    call evalSrcJacBlockSaGamma(i, j, k, A)
                    if (mode == srcLambdaModeDecoupled) then
                        srcLambda(i, j, k, 1) = max(zero, A(1, 1))
                        srcLambda(i, j, k, 2) = max(zero, A(2, 2))
                    else
                        tr2 = A(1, 1) + A(2, 2)
                        disc2 = ((A(1, 1) - A(2, 2)) * half)**2 + A(1, 2) * A(2, 1)
                        if (disc2 >= zero) then
                            lambda2x2 = tr2 * half + sqrt(disc2)
                        else
                            lambda2x2 = tr2 * half
                        end if
                        lambda2x2 = max(zero, lambda2x2)
                        srcLambda(i, j, k, 1) = lambda2x2
                        srcLambda(i, j, k, 2) = lambda2x2
                    end if
                    srcLambda(i, j, k, 3) = zero
                end do
            end do
        end do

    end subroutine computeSrcLambdaSaGamma

end module saGamma
