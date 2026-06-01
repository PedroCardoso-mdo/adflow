# Build and Test ADflow

Build ADflow with optional AD, install to venv, and run OpenMP test case.

## Steps

1. **Check if AD is needed**: Look for changes in AD-dependent files
2. **Run AD if needed**: `./AD_sh` (if automatic differentiation required)
3. **Build**: `make`
4. **Activate venv**: `source /home/mdo/packages_v2/mach/bin/activate`
5. **Install**: `pip install .`
6. **Run test**: Execute python case in `Test_OpenMp/`

## Usage

Invoke with `/build-and-test` or ask Claude to build and run the OpenMP test.

## Notes

- Test_OpenMp/ is excluded from git
- Virtual environment: `/home/mdo/packages_v2/mach/bin/activate`
