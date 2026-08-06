module paramTurb
!
!       Module that contains the constants for the turbulence models
!       as well as some global variables/parameters for the turbulent
!       routines.
!
    use constants, only: realType, intType
    implicit none
    save
!
!       Spalart-Allmaras constants.
!
    real(kind=realType), parameter :: rsaK = 0.41_realType
    real(kind=realType), parameter :: rsaCb1 = 0.1355_realType
    real(kind=realType), parameter :: rsaCb2 = 0.622_realType
    real(kind=realType), parameter :: rsaCb3 = 0.66666666667_realType
    real(kind=realType), parameter :: rsaCv1 = 7.1_realType
    real(kind=realType), parameter :: rsaCw1 = rsaCb1 / (rsaK * rsaK) &
                                      + (1.+rsaCb2) / rsaCb3
    real(kind=realType), parameter :: rsaCw2 = 0.3_realType
    real(kind=realType), parameter :: rsaCw3 = 2.0_realType
    real(kind=realType), parameter :: rsaCt1 = 1.0_realType
    real(kind=realType), parameter :: rsaCt2 = 2.0_realType
    real(kind=realType), parameter :: rsaCt3 = 1.2_realType
    real(kind=realType), parameter :: rsaCt4 = 0.5_realType
    real(kind=realType), parameter :: rsaCrot = 2.0_realType
!
!       SA-noft2-Gamma-Retheta constants.
!
    real(kind=realType), parameter :: sigmaTheta = 2.0_realType
    real(kind=realType), parameter :: sigmaF     = 1.0_realType
    real(kind=realType), parameter :: rsaGRca1   = 2.0_realType
    real(kind=realType), parameter :: rsaGRca2   = 0.06_realType
    real(kind=realType), parameter :: rsaGRce1   = 1.0_realType
    real(kind=realType), parameter :: rsaGRce2   = 50.0_realType
    real(kind=realType), parameter :: rsaGRcthetat = 0.03_realType
    real(kind=realType), parameter :: rsaGRccrossflow = 0.6_realType
!
!       SA-noft2-Gamma-Retheta correlation constants (paper Eqs. 54-57, 15-26).
!       Named module-level parameters, NOT bare literals, at smoothMinMax()
!       call sites: a bare `0.4_realType`-style literal argument fails to
!       complexify (the complex build's smoothMinMax has an explicit
!       COMPLEX interface, and Fortran does not auto-promote a REAL literal
!       actual argument to COMPLEX at a checked interface -- it only
!       auto-promotes on assignment/parameter-initialization, which is why
!       `one`/`zero`/`rsaGRpmax` work here but a bare literal doesn't).
!
    real(kind=realType), parameter :: rsaGRtuFloor = 0.027_realType
    real(kind=realType), parameter :: rsaGRlambdaThetaMin = -0.1_realType
    real(kind=realType), parameter :: rsaGRlambdaThetaMax = 0.1_realType
    ! 2026-08-03 test (NLF0416 exact-grid showed +22/+30 counts vs paper
    ! above the drag bucket): the +-0.1 lambda_theta clamp is a legacy
    ! LM2009 prescription that is NOT part of the P&Z smooth model -- their
    ! F(lambda_theta) (Eqs. 54-57) carries its own smooth saturations.
    ! TESTED 2026-08-03 and REJECTED: unclamped lambda_theta stalls the
    ! solver (~1e-4 plateau on S809 SB3 alpha 4-6 at 300k cycles) and moves
    ! cd AWAY from the paper (+106..129 counts vs clamped 67..116 vs paper
    ! ~63-85). The +-0.1 bound is part of the correlation's validity domain
    ! (separated-flow dU/ds drives lambda_theta far outside it near the
    ! bubble). Keep .true.; the switch remains for future experiments.
    logical, parameter :: rsaGRclampLambdaTheta = .true.
    real(kind=realType), parameter :: rsaGRcrossflowRatioCap = 0.4_realType
    real(kind=realType), parameter :: rsaGRhcfRef = 0.1066_realType
!
!       SA-noft2-Gamma-Retheta smooth-Fonset constants.
!
    real(kind=realType), parameter :: rsaGRfonsetC  = 2.6_realType
    real(kind=realType), parameter :: rsaGRfonsetK  = 6.0_realType
    real(kind=realType), parameter :: rsaGRfonsetS  = 1.35_realType
    real(kind=realType), parameter :: rsaGRpmin = -300.0_realType
    real(kind=realType), parameter :: rsaGRpmax  =  300.0_realType
    ! rsaGRdampTheta moved to inputIteration as transitionDampTheta (configurable)
    real(kind=realType), parameter :: rsaGRgammaLo  = 1.e-10_realType
    real(kind=realType), parameter :: rsaGRgammaHi  = 2.0_realType
    real(kind=realType), parameter :: rsaGRreThetaLo = 20.0_realType
    ! Safety margin for the gammaForSA clamp (saGammaRetheta.F90) that keeps
    ! gamma's own natural saturation values (0 = fully laminar, 1 = fully
    ! turbulent) strictly INSIDE the unclamped pass-through zone, not sitting
    ! bit-exact on the clamp boundary. A raw min(max(gamma,0),1) ties exactly
    ! at those physically-routine values, and Tapenade's forward-mode tangent
    ! picks the wrong branch there vs. complex-step ground truth (confirmed
    ! 2026-07-23: every AD-vs-CS mismatch cell in dR[nuTilde]/dw[gamma] had
    ! gamma==1.0 bit-exact, CS=0 correctly, AD nonzero incorrectly). Padding
    ! is >> the observed NK overshoot (~2e-16, i.e. 1 ULP) but still tiny
    ! relative to gamma's O(1) physical range, so this is purely a
    ! failure/divergence safeguard, not a change to normal-operation physics.
    real(kind=realType), parameter :: rsaGRgammaForSAMargin = 1.e-3_realType
!
!       SA-Gamma-Retheta stabilization: source dt restriction and scaling.
!
    real(kind=realType), parameter :: rsaGRsrcDtLimit = 0.9_realType

!       srcLambda eigenvalue modes (distinct from TurbDADICoupled).
!
    integer(kind=intType), parameter :: srcLambdaModeDecoupled  = 0_intType
    integer(kind=intType), parameter :: srcLambdaModeTransition = 1_intType
    integer(kind=intType), parameter :: srcLambdaModeFull       = 2_intType

!       transitionDebug slot count (must match comment block in saGammaRetheta.F90).
!       Do NOT renumber existing slots — append new ones at the end.
    integer(kind=intType), parameter :: nTransitionDebug = 51_intType

!
!       K-omega constants.
!
    real(kind=realType), parameter :: rkwK = 0.41_realType
    real(kind=realType), parameter :: rkwSigk1 = 0.5_realType
    real(kind=realType), parameter :: rkwSigw1 = 0.5_realType
    real(kind=realType), parameter :: rkwSigd1 = 0.5_realType
    real(kind=realType), parameter :: rkwBeta1 = 0.0750_realType
    real(kind=realType), parameter :: rkwBetas = 0.09_realType
!
!       K-omega SST constants.
!
    real(kind=realType), parameter :: rSSTK = 0.41_realType
    real(kind=realType), parameter :: rSSTA1 = 0.31_realType
    real(kind=realType), parameter :: rSSTBetas = 0.09_realType

    real(kind=realType), parameter :: rSSTSigk1 = 0.85_realType
    real(kind=realType), parameter :: rSSTSigw1 = 0.5_realType
    real(kind=realType), parameter :: rSSTBeta1 = 0.0750_realType

    real(kind=realType), parameter :: rSSTSigk2 = 1.0_realType
    real(kind=realType), parameter :: rSSTSigw2 = 0.856_realType
    real(kind=realType), parameter :: rSSTBeta2 = 0.0828_realType
!
!       K-tau constants.
!
    real(kind=realType), parameter :: rktK = 0.41_realType
    real(kind=realType), parameter :: rktSigk1 = 0.5_realType
    real(kind=realType), parameter :: rktSigt1 = 0.5_realType
    real(kind=realType), parameter :: rktSigd1 = 0.5_realType
    real(kind=realType), parameter :: rktBeta1 = 0.0750_realType
    real(kind=realType), parameter :: rktBetas = 0.09_realType
!
!       V2-f constants.
!
    real(kind=realType), parameter :: rvfC1 = 1.4_realType
    real(kind=realType), parameter :: rvfC2 = 0.3_realType
    real(kind=realType), parameter :: rvfBeta = 1.9_realType
    real(kind=realType), parameter :: rvfSigk1 = 1.0_realType
    real(kind=realType), parameter :: rvfSige1 = 0.7692307692_realType
    real(kind=realType), parameter :: rvfSigv1 = 1.00_realType
    real(kind=realType), parameter :: rvfCn = 70.0_realType

    real(kind=realType), parameter :: rvfN1Cmu = 0.190_realType
    real(kind=realType), parameter :: rvfN1A = 1.300_realType
    real(kind=realType), parameter :: rvfN1B = 0.250_realType
    real(kind=realType), parameter :: rvfN1Cl = 0.300_realType
    real(kind=realType), parameter :: rvfN6Cmu = 0.220_realType
    real(kind=realType), parameter :: rvfN6A = 1.400_realType
    real(kind=realType), parameter :: rvfN6B = 0.045_realType
    real(kind=realType), parameter :: rvfN6Cl = 0.230_realType

    real(kind=realType) :: rvfLimitK, rvfLimitE, rvfCl
    real(kind=realType) :: rvfCmu
!
!       Variables to store the parameters for the wall functions fits.
!       As these variables depend on the turbulence model they are set
!       during runtime. Allocatables are used, because the number of
!       fits could be different for the different models.
!       The curve is divided in a number of intervals and is
!       constructed such that both the function and the derivatives
!       are continuous. Consequently cubic polynomials are used.
!
    ! nFit:               Number of intervals of the curve.
    ! ypT(0:nFit):        y+ values at the interval boundaries.
    ! reT(0:nFit):        Reynolds number at the interval
    !                     boundaries, where the Reynolds number is
    !                     defined with the local velocity and the
    !                     wall distance.
    ! up0(nFit):          Coefficient 0 in the fit for the
    !                     nondimensional tangential velocity as a
    !                     function of the Reynolds number.
    ! up1(nFit):          Idem for coefficient 1.
    ! up2(nFit):          Idem for coefficient 2.
    ! up3(nFit):          Idem for coefficient 3.
    ! tup0(nFit,nt1:nt2): Coefficient 0 in the fit for the
    !                     nondimensional turbulence variables as a
    !                     function of y+.
    ! tup1(nFit,nt1:nt2): Idem for coefficient 1.
    ! tup2(nFit,nt1:nt2): Idem for coefficient 2.
    ! tup3(nFit,nt1:nt2): Idem for coefficient 3.
    ! tuLogFit(nt1:nt2):  Whether or not the logarithm of the variable
    !                     has been fitted.

    integer(kind=intType) :: nFit

    real(kind=realType), dimension(:), allocatable :: ypT, reT
    real(kind=realType), dimension(:), allocatable :: up0, up1
    real(kind=realType), dimension(:), allocatable :: up2, up3

    real(kind=realType), dimension(:, :), allocatable :: tup0, tup1
    real(kind=realType), dimension(:, :), allocatable :: tup2, tup3
#ifndef USE_TAPENADE
    real(kind=realType), dimension(:), allocatable :: ypTb, reTb
    real(kind=realType), dimension(:), allocatable :: up0b, up1b
    real(kind=realType), dimension(:), allocatable :: up2b, up3b
#endif

    logical, dimension(:), allocatable :: tuLogFit

end module paramTurb
