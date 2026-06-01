# Build and Test ADflow

Build ADflow with optional AD, install to venv, and run OpenMP test case.

## Steps

1. **Check if AD is needed**: Check if any AD-trigger files have changed (see list below)
2. **Run AD if needed**: `./AD_sh` (if automatic differentiation required)
3. **Build**: `make`
4. **Activate venv**: `source /home/mdo/packages_v2/mach/bin/activate`
5. **Install**: `pip install .`
6. **Run test**: 
   ```bash
   cd Test_OpenMp
   OMP_NUM_THREADS=4 mpirun -np 1 python analysis.py
   ```

## AD Trigger Files

If any of these files have changed since last AD run, `./AD_sh` must be executed:

```
src/adjoint/adjointExtra.F90
src/modules/surfaceFamilies.F90
src/modules/monitor.f90
src/modules/diffSizes.f90
src/modules/block.F90
src/modules/inputParam.F90
src/modules/constants.F90
src/modules/precision_tapenade.f90
src/modules/iteration.f90
src/modules/section.f90
src/modules/communication.F90
src/modules/paramTurb.F90
src/modules/cgnsGrid.F90
src/modules/CpCurveFits.f90
src/modules/blockPointers.F90
src/modules/BCPointers.F90
src/modules/flowVarRefState.F90
src/modules/wallDistanceData.F90
src/modules/actuatorRegionData.F90
src/modules/overset.F90
src/modules/cgnsNames.f90
src/modules/su_cgns.F90
src/modules/BCDataMod.F90
src/solver/BCRoutines.F90
src/solver/fluxes.F90
src/solver/solverUtils.F90
src/solver/residuals.F90
src/solver/surfaceIntegrations.F90
src/solver/zipperIntegrations.F90
src/solver/ALEUtils.F90
src/solver/actuatorRegion.F90
src/initFlow/initializeFlow.F90
src/turbulence/turbUtils.F90
src/turbulence/turbBCRoutines.F90
src/turbulence/turbMod.F90
src/turbulence/sa.F90
src/utils/flowUtils.F90
src/utils/utils.F90
src/utils/sorting.F90
src/overset/oversetUtilities.F90
src/wallDistance/wallDistance.F90
src/bcdata/BCData.F90
```

## Usage

Invoke with `/build-and-test` or ask Claude to build and run the OpenMP test.

## Notes

- Test_OpenMp/ is excluded from git
- Virtual environment: `/home/mdo/packages_v2/mach/bin/activate`
- Test with 1 MPI rank and 4 OpenMP threads to isolate OpenMP performance
- AD trigger files extracted from `src/adjoint/Makefile_tapenade` (ALL_RES_FILES)
