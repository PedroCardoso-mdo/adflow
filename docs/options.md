# Options added by this fork

Every option below exists in `adflow/pyADflow.py` of this fork and not in upstream ADflow.
The list is generated from the difference between the two default-option dictionaries
(`_getDefaultOptions`), so it is complete: **51 new options** (24 model, 24 solver, 3 adjoint),
plus three upstream options that gained a choice. Upstream options are documented in
`doc/options.yaml`, which also carries an entry for every option listed here.

"Models" says which model reads the option: GR = `SA-noft2-Gamma-Retheta`, SA-γ = `SA-noft2-Gamma`,
SA-BCM = `useSABCM`, all = independent of the turbulence model. Model equations: [models.md](models.md).
Validated values: [recipes.md](recipes.md). Options that existed during development and were removed:
[tried_and_rejected.md](tried_and_rejected.md).

## Model options (24)

| Option                       | Default   | Models           | Meaning                                                                                                                                                                                                                                   |
| ---------------------------- | --------- | ---------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `transitionFirstOrderUpwind` | `True`    | GR, SA-γ         | First-order upwind convection for γ and R̅eθt whatever `turbulenceOrder` says (P&Z Sec. IV.A). Every validated run has it on.                                                                                                              |
| `transitionCrossflow`        | `False`   | GR               | Adds the helicity-based stationary-crossflow sink D_scf to the R̅eθt equation (P&Z Eqs. 15–26). Identically zero in 2D. Validated only with matrix dissipation.                                                                            |
| `transitionRoughnessHeight`  | `3.3e-06` | GR               | Surface roughness height h of the crossflow correlation (P&Z Eq. 17), in mesh length units (3.3e-6 = 3.3 µm on a metre-scaled mesh). Only read with `transitionCrossflow`.                                                                |
| `transitionArtVisc`          | `0.0`     | GR               | Coefficient c of a gradient-driven artificial viscosity on the γ equation, ν_art = c·V·h·Δγ per face (velocity magnitude, face-normal cell size, γ jump across the face). 0 = off. Changes the solution; experimental, no validated case. |
| `transitionUseApproxSA`      | `True`    | GR, SA-γ         | Lets the ν̃ equation of the transition models use ADflow's approximate-SA simplification (production term dropped) while the first-order ANK stage requests it.                                                                            |
| `sgammaOnsetTanh`            | `False`   | SA-γ             | Use the P&Z tanh onset function instead of the smoothed Menter/Lee min–max onset.                                                                                                                                                         |
| `sgammaFturbLee`             | `False`   | SA-γ             | F_turb = exp(−(R_T/2)⁴) (Lee) instead of (1 − F_onset)·exp(−R_T).                                                                                                                                                                         |
| `sgammaCoupleDestruction`    | `True`    | SA-γ             | Also scale the SA destruction term by φ⁺(γ_s, 0.1) (Nichols/Lee coupling). False: γ acts on production only.                                                                                                                              |
| `sgammaFPGSmoothP`           | `300.0`   | SA-γ             | Sharpness p of the smooth min/max operators inside the pressure-gradient function F_PG.                                                                                                                                                   |
| `sgammaTanhK`                | `6.0`     | SA-γ             | Slope k of the tanh onset, F_onset = A·[tanh(k·(F₁ − s)) + 1]/2. Only with `sgammaOnsetTanh`.                                                                                                                                             |
| `sgammaTanhS`                | `1.35`    | SA-γ             | Centre s of the tanh onset. Only with `sgammaOnsetTanh`.                                                                                                                                                                                  |
| `sgammaTanhAmp`              | `1.0`     | SA-γ             | Amplitude A of the tanh onset (the default min–max onset peaks at 2). Only with `sgammaOnsetTanh`.                                                                                                                                        |
| `sgammaTanhRT`               | `1.0`     | SA-γ             | Scale c of R_T in the tanh onset argument, F₁ = √[(Re_v/(2.2·Re_θc))² + (R_T/c)²]. Only with `sgammaOnsetTanh`.                                                                                                                           |
| `sgammaOnsetSmoothP`         | `300.0`   | SA-γ             | Sharpness p of the three smooth min/max operators of the default onset function. Unused with `sgammaOnsetTanh`.                                                                                                                           |
| `sgammaArtVisc`              | `0.0`     | SA-γ             | Coefficient c of a gradient-driven artificial viscosity on γ, ν_art = c·V·h·Δγ per face (velocity magnitude, face-normal cell size, γ jump across the face). 0 = off. Changes cd by −1.4…−3.7 % at c = 5–10.                              |
| `transitionRefLength`        | `-1.0`    | GR               | Reference length l of the vorticity limiter (P&Z Eqs. 52–53), in mesh length units. ≤ 0: use the AeroProblem `chordRef`.                                                                                                                  |
| `useSABCM`                   | `False`   | SA-BCM           | Switches on the SA-BCM algebraic intermittency (needs `turbulenceModel="SA"`, RANS, `useApproxWallDistance=True`). ft2 is forced off.                                                                                                     |
| `SABCMSmooth`                | `True`    | SA-BCM           | True: differentiable form (KS-smoothed Term1, tanh intermittency). False: original form γ = 1 − exp(−(√Term1 + √Term2)).                                                                                                                  |
| `SABCMChi1`                  | `0.002`   | SA-BCM           | Model constant χ₁ (Term1).                                                                                                                                                                                                                |
| `SABCMChi2`                  | `0.02`    | SA-BCM           | Model constant χ₂ (Term2).                                                                                                                                                                                                                |
| `SABCMRho`                   | `50.0`    | SA-BCM           | KS aggregation parameter that replaces max(Term1, 0). Smooth form only.                                                                                                                                                                   |
| `SABCMTanhCenter`            | `0.5`     | SA-BCM           | Centre of the tanh intermittency blend. Smooth form only.                                                                                                                                                                                 |
| `SABCMTanhWidth`             | `0.08`    | SA-BCM           | Width of the tanh intermittency blend. Smooth form only.                                                                                                                                                                                  |
| `turbIntensityInf`           | `0.001`   | GR, SA-γ, SA-BCM | Free-stream turbulence intensity as a fraction (0.001 = 0.1 %). Input of the onset correlations of all three models; sets R̅eθt∞ for GR.                                                                                                   |

## Solver options (23)

`epsAcoustic` and `epsShear` are discretisation coefficients; they are listed here because they are
not tied to a model.

| Option                      | Default  | Models   | Meaning                                                                                                                                                                                                             |
| --------------------------- | -------- | -------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `transitionNK`              | `True`   | GR, SA-γ | Master switch of the transition-specific Newton–Krylov machinery in turbulence KSP, coupled ANK and NK: column scaling of the turbulence unknowns, source-term time-step restriction and the guarded fixed NK step. |
| `transitionSrcDtRestrict`   | `True`   | GR, SA-γ | Source-term time-step restriction λ_source·Δt ≤ limit (P&Z Eq. 59) on the turbulence diagonal.                                                                                                                      |
| `transitionSrcDtLimit`      | `0.9`    | GR, SA-γ | The limit in λ_source·Δt ≤ limit.                                                                                                                                                                                   |
| `srcDtDeactivateIters`      | `5`      | GR, SA-γ | Number of consecutive steps without backtracking after which the restriction is released in the Krylov paths (re-armed on a backtrack). 0 = released from the start. The DADI path ignores it.                      |
| `transitionDampTheta`       | `0.99`   | GR, SA-γ | Back-off factor θ of the per-node update damping of γ and R̅eθt (P&Z Algorithm 2).                                                                                                                                   |
| `transitionDampMaxIter`     | `10000`  | GR, SA-γ | Cap on the back-off loop of the per-node damping; a hard clip to the bounds follows only if it is exhausted.                                                                                                        |
| `ANKCFLMinCap`              | `0.0`    | all      | Absolute ceiling on the convergence-ramped ANK CFL floor, so that the CFL cutback keeps room to act deep in the convergence. ≤ 0 = off (upstream behaviour).                                                        |
| `ANKUnsteadyLSFactor`       | `0.7`    | all      | Backtracking factor of the ANK unsteady line search (upstream hard-codes 0.7).                                                                                                                                      |
| `ANKUnsteadyLSMaxIter`      | `12`     | all      | Maximum number of backtracks of the ANK unsteady line search (upstream hard-codes 12).                                                                                                                              |
| `ANKRejectOnLSExhausted`    | `False`  | all      | Treat an exhausted ANK line search as a rejected step, so that the CFL cutback fires. False = upstream behaviour.                                                                                                   |
| `ANKTransitionGlobalLambda` | `True`   | GR, SA-γ | True: γ and R̅eθt take part in the global physicality step limit of coupled ANK. False: they are excluded from it and bounded per node (P&Z Algorithm 2) instead.                                                    |
| `MFFDFunctionError`         | `0.0`    | all      | Relative residual-evaluation error handed to PETSc for the matrix-free differencing step of ANK/NK. ≤ 0 = PETSc's machine-epsilon assumption (upstream behaviour).                                                  |
| `MFFDType`                  | `"ds"`   | all      | Matrix-free differencing formula: `ds` (PETSc default) or `wp` (Walker–Pernice).                                                                                                                                    |
| `ANKColScaleUnit`           | `False`  | all      | Scale every unknown of the ANK system to unit RMS before the matrix-free product. Factors are refreshed when the preconditioner is rebuilt.                                                                         |
| `epsAcoustic`               | `0.25`   | all      | Discretisation: eigenvalue-limiter coefficient of the acoustic eigenvalues in the matrix dissipation scheme (V_n). The default is the value hard-coded upstream.                                                    |
| `epsShear`                  | `0.025`  | all      | Discretisation: eigenvalue-limiter coefficient of the convective (shear) eigenvalues in the matrix dissipation scheme (V_l). The default is the value hard-coded upstream; the transition validations use 0.        |
| `NKLSRelax`                 | `False`  | all      | Relaxed NK cubic line search: Armijo constant 1e-3 and turbulence pre-limit factor 3.0 instead of upstream's 1e-2 and 2.0. pyADflow sets it True for the two γ models.                                              |
| `TurbDADICoupled`           | `"full"` | GR, SA-γ | Coupling of the turbulence DADI solve: `decoupled` (scalar solves), `transition` (ν̃ scalar + γ–R̅eθt block), `full` (3×3 block for GR, 2×2 for SA-γ).                                                                |
| `turbL2Convergence`         | `0.0`    | all      | Per-equation convergence gate: the solve also requires every turbulence and transition residual to be below this factor times its own peak. 0 = off (total residual only, upstream behaviour).                      |
| `transitionGammaFloor`      | `0.0`    | GR, SA-γ | Lower bound of γ in the ANK update limiters only (the residual keeps 1e-10). Releases cells stuck on the bound with a constant γ residual; 1e-4 to 1e-2. 0 = off.                                                   |
| `ANKStallSwitchIters`       | `0`      | all      | Stall-triggered phase advance of ANK: if the residual is not halved within this many iterations (below 1e-3 relative, steps ≥ 0.5), go coupled, then second order, without waiting for the thresholds. 0 = off.     |
| `ANKStallSwitchCharOff`     | `False`  | all      | Third stage of the stall switch: a stall in second-order coupled ANK turns the characteristic time step off for the rest of the solve (floor of laminar-bubble cases at low Mach).                                  |
| `ANKCharMachFloor`          | `1e-4`   | all      | K in the cut-off M² ≥ K·M∞² of the local Mach number in the characteristic time step. The default is the value hard-coded upstream.                                                                                 |
| `ANKCharPresCutoff`         | `0.0`    | all      | ε of a pressure cut-off of the local Mach number in the characteristic time step, M² ≥ ε²·Δp/(ρa²). 0 = off (upstream behaviour).                                                                                   |

## Adjoint options (3)

| Option                    | Default                                                    | Models   | Meaning                                                                                                             |
| ------------------------- | ---------------------------------------------------------- | -------- | ------------------------------------------------------------------------------------------------------------------- |
| `adjointLGMRESAugDim`     | `2`                                                        | all      | Number of error-approximation vectors LGMRES carries across restarts. Only with `adjointSolver="LGMRES"`.           |
| `adjointFieldSplitType`   | `"multiplicative"` (choices: `multiplicative`, `additive`) | GR, SA-γ | Composition of the field-split preconditioner blocks. Only with `globalPreconditioner="field split"`.               |
| `adjointFieldSplitBlocks` | `3`                                                        | GR, SA-γ | 3: blocks {flow}, {ν̃}, {γ, R̅eθt}. 4: γ and R̅eθt in separate blocks. Only with `globalPreconditioner="field split"`. |

## Upstream options with a new choice

| Option                 | New choice                 | Meaning                                                                           |
| ---------------------- | -------------------------- | --------------------------------------------------------------------------------- |
| `turbulenceModel`      | `"SA-noft2-Gamma-Retheta"` | GR model: SA (no ft2) + γ + R̅eθt transport equations.                             |
| `turbulenceModel`      | `"SA-noft2-Gamma"`         | SA-γ model: SA (no ft2) + one γ transport equation.                               |
| `adjointSolver`        | `"LGMRES"`                 | PETSc LGMRES: restarted GMRES augmented with `adjointLGMRESAugDim` error vectors. |
| `globalPreconditioner` | `"field split"`            | PETSc PCFIELDSPLIT over the physical blocks, each block with ASM + ILU.           |

## Defaults that depend on the model

Upstream defaults are unchanged. When `turbulenceModel` contains "Gamma" (GR and SA-γ) pyADflow replaces
the following defaults, and only for options the user did not pass:

| Option                 | Upstream default     | Default for GR and SA-γ | Reason                                                                                         |
| ---------------------- | -------------------- | ----------------------- | ---------------------------------------------------------------------------------------------- |
| `globalPreconditioner` | `"additive Schwarz"` | `"field split"`         | adjoint time −29 % in 2D, 5.0–5.7× faster on a 1.48 M-cell wing (campaign `02_adjoint_checks`) |
| `ILUFill`              | `2`                  | `3`                     | adjoint preconditioner strength used in those tests                                            |
| `ASMOverlap`           | `1`                  | `3`                     | as above                                                                                       |
| `adjointSolver`        | `"GMRES"`            | `"LGMRES"`              | does not stagnate at small subspaces where GMRES does                                          |
| `adjointSubspaceSize`  | `100`                | `400`                   | field split stagnates at 200 or less                                                           |
| `NKLSRelax`            | `False`              | `True`                  | NK line search otherwise stays at its minimum step                                             |

Also set by pyADflow for GR and SA-γ:

| Setting                       | Value                                    | Note                                                                                                                 |
| ----------------------------- | ---------------------------------------- | -------------------------------------------------------------------------------------------------------------------- |
| `turbResScale` (if not given) | GR `[1e4, 0.1, 1e-4]`, SA-γ `[1e4, 0.1]` | about 1 / state magnitude for ν̃, γ, R̅eθt                                                                             |
| `useBlockettes`               | forced `False`                           | the blockette kernels exist and are tested against the block residual, but every validated solve used the block path |
| `useWallFunctions=True`       | error                                    | the models need a resolved boundary layer                                                                            |
| `useft2SA=True`               | warning                                  | the models are calibrated without ft2                                                                                |
| `turbIntensityInf ≤ 0`        | warning                                  | an internal floor of 0.027 % takes over                                                                              |

SA-BCM keeps every upstream default. `useSABCM` without RANS + SA or with `useApproxWallDistance=False`
stops with an error; `useft2SA=True` is ignored with a message.

## Output and monitoring keywords

These are values of existing options (`volumeVariables`, `surfaceVariables`, `isoVariables`,
`monitorVariables`), not new options.

| Keyword                                                                                                                                                                                           | Where              | Content                                                                                |
| ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------ | -------------------------------------------------------------------------------------- |
| `intermittency`                                                                                                                                                                                   | volume, isosurface | γ of the transported models, the recomputed algebraic γ for SA-BCM, 1 otherwise        |
| `fonset`, `fonset1`, `flength`, `fturb`, `rturb`                                                                                                                                                  | volume, surface    | onset and length functions, R_T                                                        |
| `rethetac`, `rethetatilde`, `resovercrit`, `fthetat`, `fwake`, `translambdatheta`, `transpretheta`, `transtimescale`, `thetabl`, `deltabl`, `transdelta`                                          | volume, surface    | R̅eθt-equation quantities (GR; SA-γ reuses some slots, see the header of `saGamma.F90`) |
| `gammaprod`, `gammadest`, `gammaforsa`, `gammalocal`                                                                                                                                              | volume, surface    | γ source terms and the clipped γ values                                                |
| `dscf`, `rescf`, `hcf`                                                                                                                                                                            | volume, surface    | crossflow source, crossflow Reynolds number, helicity parameter                        |
| `transvortmag`, `transvortmaglim`, `transvortlim`, `strainmag`, `sastrainrate`, `samodstrainrate`, `ft2`, `nutsa`, `transvelmag`, `duds`, `dudx` … `dwdz`, `transwalldist`, `transrho`, `transmu` | volume, surface    | kinematics and SA quantities seen by the transition sources                            |
| `qq11`, `qq22`, `qq33`                                                                                                                                                                            | volume, surface    | diagonal of the implicit turbulence source Jacobian                                    |
| `scaledtotalr`                                                                                                                                                                                    | monitor            | total residual with the turbulence rows scaled; display only, not used by any switch   |

All transition diagnostics come from one per-block array filled in the source routine; there is no
text-file debug output.

## Removed options

36 options of the development branch are not in this code. Passing one raises the usual
pyADflow "option not valid" error. They are listed with the reason in
[tried_and_rejected.md](tried_and_rejected.md).
