










module cgnsNames
!
!       Parametrized cgns names of the variables used in this code.
!
    use constants, only: maxCGNSNameLen
    implicit none
    save
!
!       Time history values.
!
    character(len=maxCGNSNameLen), parameter :: &
        cgnsTimeValue = "TimeValues"
!
!       coordinate names.
!
    character(len=maxCGNSNameLen), parameter :: &
        cgnsCoorX = "CoordinateX"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsCoorY = "CoordinateY"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsCoorZ = "CoordinateZ"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsCoorR = "CoordinateR"
!
!       Variable names.
!
    character(len=maxCGNSNameLen), parameter :: &
        cgnsDensity = "Density"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsMomX = "MomentumX"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsMomY = "MomentumY"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsMomZ = "MomentumZ"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsEnergy = "EnergyStagnationDensity"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsTurbSANu = "TurbulentSANuTilde"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsTurbK = "TurbulentEnergyKinetic"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsTurbOmega = "TurbulentDissipationRate"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsTurbTau = "TurbulentInvDissipationRate"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsTurbEpsilon = "TurbulentDissipation"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsTurbV2 = "TurbulentScalarV2"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsTurbF = "TurbulentScalarF"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsTurbGamma = "TransitionGamma"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsTurbRetheta = "TransitionRetheta"

    character(len=maxCGNSNameLen), parameter :: &
        cgnsVelX = "VelocityX"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsVelY = "VelocityY"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsVelZ = "VelocityZ"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsRelVelX = "RelativeVelocityX"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsRelVelY = "RelativeVelocityY"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsRelVelZ = "RelativeVelocityZ"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsVelr = "VelocityR"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsVelTheta = "VelocityTheta"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsPressure = "Pressure"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsTemp = "Temperature"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsCp = "CoefPressure"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsMach = "Mach"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsRelMach = "RelativeMach"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsMachTurb = "MachTurbulent"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsViscMol = "ViscosityMolecular"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsViscKin = "ViscosityKinematic"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsEddy = "ViscosityEddy"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsEddyRatio = "ViscosityEddyRatio"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsWallDist = "TurbulentDistance"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsVortMagn = "VorticityMagnitude"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsVortX = "VorticityX"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsVorty = "VorticityY"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsVortZ = "VorticityZ"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsPtotLoss = "RelativePressureStagnationLoss"

    character(len=maxCGNSNameLen), parameter :: &
        cgnsRhoTot = "DensityStagnation"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsPTot = "PressureStagnation"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsTTot = "TemperatureStagnation"

    character(len=maxCGNSNameLen), parameter :: &
        cgnsSkinFMag = "SkinFrictionMagnitude"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsSkinFX = "SkinFrictionX"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsSkinFY = "SkinFrictionY"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsSkinFZ = "SkinFrictionZ"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsForceInDragDir = "ForceInDragDir"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsForceInLiftDir = "ForceInLiftDir"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsStanton = "StantonNumber"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsYPlus = "YPlus"

    character(len=maxCGNSNameLen), parameter :: &
        cgnsVelocity = "Mach_Velocity"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsSoundSpeed = "Mach_VelocitySound"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsLength = "LengthReference"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsReyn = "Reynolds"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsReynLen = "Reynolds_Length"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsHeatRatio = "SpecificHeatRatio"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsPrandtl = "Prandtl"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsPrandtlTurb = "PrandtlTurbulent"

    character(len=maxCGNSNameLen), parameter :: &
        cgnsVelAngleX = "VelocityAngleX"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsVelAngleY = "VelocityAngleY"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsVelAngleZ = "VelocityAngleZ"

    character(len=maxCGNSNameLen), parameter :: &
        cgnsVelVecX = "VelocityUnitVectorX"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsVelVecY = "VelocityUnitVectorY"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsVelVecZ = "VelocityUnitVectorZ"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsVelVecR = "VelocityUnitVectorR"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsVelVecTheta = "VelocityUnitVectorTheta"

    character(len=maxCGNSNameLen), parameter :: &
        cgnsMassFlow = "MassFlow"

    character(len=maxCGNSNameLen), parameter :: &
        cgnsShock = "Shock"

    character(len=maxCGNSNameLen), parameter :: &
        cgnsFilteredShock = "FilteredShock"

    character(len=maxCGNSNameLen), parameter :: &
        cgnsGC = "globalCell"

    character(len=maxCGNSNameLen), parameter :: &
        cgnsStatus = "status"

    character(len=maxCGNSNameLen), parameter :: &
        cgnsIntermittency = "Intermittency"

    character(len=maxCGNSNameLen), parameter :: &
        cgnsFonset = "TransitionFonset"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsFlength = "TransitionFlength"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsRturb = "TransitionRturb"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsReThetaTarget = "TransitionReThetaTarget"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsFonset1 = "TransitionFonset1"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsReS = "TransitionReS"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsReThetaC = "TransitionReThetaC"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsReSOverCrit = "TransitionReSOverCrit"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsStrainMag = "StrainMagnitude"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsDudx = "VelocityGradientDudx"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsDudy = "VelocityGradientDudy"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsDudz = "VelocityGradientDudz"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsDvdx = "VelocityGradientDvdx"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsDvdy = "VelocityGradientDvdy"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsDvdz = "VelocityGradientDvdz"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsDwdx = "VelocityGradientDwdx"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsDwdy = "VelocityGradientDwdy"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsDwdz = "VelocityGradientDwdz"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsFthetaT = "TransitionFthetaT"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsFwake = "TransitionFwake"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsGammaProd = "TransitionGammaProduction"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsGammaDest = "TransitionGammaDestruction"
    ! slots 22-24, 27-29 (slot 21 = "TransitionGamma" already via cgnsTurbGamma)
    character(len=maxCGNSNameLen), parameter :: &
        cgnsTransWallDist = "TransitionWallDist"
    ! slots 30-48
    character(len=maxCGNSNameLen), parameter :: &
        cgnsGammaForSA = "TransitionGammaForSA"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsGammaLocal = "TransitionGammaLocal"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsTransVortMag = "TransitionVortMag"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsTransVortMagLim = "TransitionVortMagLim"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsTransFturb = "TransitionFturb"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsSAStrainRate = "TransitionSAStrainRate"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsSAModStrainRate = "TransitionSAModStrainRate"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsTransFt2 = "TransitionFt2"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsTransThetaBL = "TransitionThetaBL"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsTransDeltaBL = "TransitionDeltaBL"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsTransDelta = "TransitionDelta"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsTransVelMag = "TransitionVelMag"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsTransDUds = "TransitionDUds"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsTransNutSA = "TransitionNutSA"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsTransReThetaTilde = "TransitionReThetaTilde"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsTransVortLim = "TransitionVortLim"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsTransQQ11 = "TransitionQQ11"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsTransQQ22 = "TransitionQQ22"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsTransQQ33 = "TransitionQQ33"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsTransRho = "TransitionRho"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsTransMu = "TransitionMu"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsTransTimeScale = "TransitionTimeScale"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsTransLambdaTheta = "TransitionLambdaTheta"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsTransPReTheta = "TransitionPReTheta"

!
!       Residual names.
!
    character(len=maxCGNSNameLen), parameter :: &
        cgnsResRho = "ResDensity"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsResMomX = "ResMomentumX"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsResMomY = "ResMomentumY"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsResMomZ = "ResMomentumZ"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsResRhoE = "ResEnergyStagnationDensity"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsResNu = "ResTurbulentSANuTilde"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsResGamma = "ResTransitionGamma"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsResRetheta = "ResTransitionRetheta"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsResK = "ResTurbulentEnergyKinetic"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsResOmega = "ResTurbulentDissipationRate"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsResTau = "ResTurbulentInvDissipationRate"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsResEpsilon = "ResTurbulentDissipation"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsResV2 = "ResTurbulentScalarV2"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsResF = "ResTurbulentScalarF"
!
!       Residual L2 norm names.
!
    character(len=maxCGNSNameLen), parameter :: &
        cgnsL2ResRho = "RSDMassRMS"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsL2ResMomX = "RSDMomentumXRMS"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsL2ResMomy = "RSDMomentumYRMS"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsL2ResMomZ = "RSDMomentumZRMS"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsL2ResRhoE = "RSDEnergyStagnationDensityRMS"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsL2ResNu = "RSDTurbulentSANuTildeRMS"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsL2ResGamma = "RSDTransitionGammaRMS"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsL2ResRetheta = "RSDTransitionRethetaRMS"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsL2ResK = "RSDTurbulentEnergyKineticRMS"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsL2ResOmega = "RSDTurbulentDissRateRMS"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsL2ResTau = "RSDTurbulentInvDissRateRMS"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsL2resEpsilon = "RSDTurbulentDissRMS"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsL2resV2 = "RSDTurbulentScalarV2RMS"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsL2ResF = "RSDTurbulentScalarFRMS"
!
!       Force and moment coefficients names.
!
    character(len=maxCGNSNameLen), parameter :: &
        cgnsCL = "CoefLift"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsCLp = "CoefPressureLift"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsCLv = "CoefViscousLift"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsCD = "CoefDrag"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsCDp = "CoefPressureDrag"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsCDv = "CoefViscousDrag"

    character(len=maxCGNSNameLen), parameter :: &
        cgnsCFx = "CoefForceX"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsCFy = "CoefForceY"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsCFz = "CoefForceZ"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsCMx = "CoefMomentX"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsCMy = "CoefMomentY"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsCMz = "CoefMomentZ"
!
!       Names of the "maximum" variables.
!
    character(len=maxCGNSNameLen), parameter :: &
        cgnsHDiffMax = "MaxDiffHAndHinf"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsMachMax = "MaxMach"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsYPlusMax = "MaxYplus"
    character(len=maxCGNSNameLen), parameter :: &
        cgnsEddyMax = "MaxRatioEddyAndLaminarViscosity"
!
!       Names of the blanking paramter.
!
    character(len=maxCGNSNameLen), parameter :: cgnsBlank = "Iblank"
!
!       Names of the "lift" force, separation sensor and cavitation
!
    character(len=maxCGNSNameLen), parameter :: cgnsSepSensor = "SepSensor"
    character(len=maxCGNSNameLen), parameter :: cgnsCavitation = "Cavitation"
    character(len=maxCGNSNameLen), parameter :: cgnsAxisMoment = "AxisMoment"
!
!       Names for the convergence history and time history.
!
    character(len=maxCGNSNameLen), parameter :: &
        ConvHistory = "ConvergenceHistory"
    character(len=maxCGNSNameLen), parameter :: &
        TimeHistory = "TimeHistory"

end module cgnsNames
