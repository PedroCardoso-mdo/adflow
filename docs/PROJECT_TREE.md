# ADflow Project Tree (SA-BCM Timing Branch)

ADflow is a CFD (Computational Fluid Dynamics) solver from MDOLab.

```
adflow_sabcm_timing/
├── adflow/              # Python package (main API/wrapper)
├── config/              # Build configuration files
├── doc/                 # Documentation
├── examples/            # Example cases/scripts
├── input_files/         # Input mesh/config files for tests
├── src/                 # Fortran source code (core solver)
│   ├── adjoint/         # Adjoint solver for sensitivities
│   ├── ADT/             # Alternating Digital Tree (search)
│   ├── bcdata/          # Boundary condition data handling
│   ├── f2py/            # Fortran-to-Python interface
│   ├── initFlow/        # Flow initialization routines
│   ├── inputParam/      # Input parameter parsing
│   ├── metis-4.0/       # METIS graph partitioning library
│   ├── modules/         # Fortran modules (shared data)
│   ├── NKSolver/        # Newton-Krylov nonlinear solver
│   ├── output/          # Output/post-processing
│   ├── overset/         # Overset mesh handling
│   ├── partitioning/    # Domain decomposition
│   ├── preprocessing/   # Mesh preprocessing
│   ├── solver/          # Main flow solver routines
│   ├── turbulence/      # Turbulence models (SA, SA-BCM)
│   ├── utils/           # Utility functions
│   ├── wallDistance/    # Wall distance computation
│   └── warping/         # Mesh warping
├── src_cs/              # Complex-step derivative source
├── tests/               # Test suite
└── setup.py             # Python package setup
```

## Key Directories

| Directory | Purpose |
|-----------|---------|
| `src/turbulence/` | SA, SA-BCM turbulence models |
| `src/NKSolver/` | Newton-Krylov solver + timing profiling |
| `src/modules/` | Shared Fortran modules and parameters |
| `src/solver/` | Main flow solver routines |
