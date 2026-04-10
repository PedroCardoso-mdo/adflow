# ADflow implementation brief: smooth local correlation-based transition model for SA

Based on: Piotrowski and Zingg, *Smooth Local Correlation-Based Transition Model for the Spalart-Allmaras Turbulence Model*, AIAA Journal, 2020.

## Goal of this implementation stage

Implement the infrastructure for the SA-sLM2015 transition model in ADflow, but stop before turning on the transition source terms and before coupling intermittency into the SA production term. The objective of this stage is to get the two new transported variables integrated cleanly into the solver, Jacobian assembly path, halo exchange, boundary-condition path, restart/output path, and diagnostics, using the already-working SA-gamma-Retheta and SST-based transition implementations as the template.

This stage should therefore produce a solver variant with:

- two extra transported variables:
  - intermittency `gamma`
  - transition onset momentum-thickness Reynolds number `ReThetaTilde`
- advection and diffusion terms implemented
- all metric, gradient, halo, BC, restart, residual, linearization, and output plumbing implemented
- correlation helper functions implemented and callable
- source terms disabled or hard-set to zero
- SA production coupling disabled, so SA remains unchanged for this first validation stage

The point is not to re-invent the model. Copy as much structure as possible from the existing `SA-gamma-Retheta` and `SST transition` paths. Only the model-specific formulas should differ.

---

## High-level model structure

The full SA-LM2015 / SA-sLM2015 model is a 3-equation system:

- SA transported eddy viscosity `nuTilde`
- intermittency `gamma`
- transported transition onset momentum-thickness Reynolds number `ReThetaTilde`

In the final model, `gamma` multiplies the SA production term. In this first stage, do **not** activate that coupling yet.

The paper modifies the LM2015 `gamma-ReThetaTilde` model so it can be coupled to SA instead of SST, while preserving a fully local formulation. The smooth version mainly replaces min/max and piecewise kinks with smooth approximations and reformulates some stiff trigger functions to behave better in Newton-Krylov. That is the version to target for the final code path.

For this first implementation stage, focus on the transport-equation framework and helper functions. Keep source terms off.

---

## Equations to implement structurally now

### 1. Transport equation for `ReThetaTilde`

Use the non-conservative form consistent with the SA implementation:

`d(ReThetaTilde)/dt + u_j d(ReThetaTilde)/dx_j = P_theta_t + D_scf + (1/Re) d/dx_j [ sigma_theta_t (nu + nu_t) d(ReThetaTilde)/dx_j ]`

Constants:

- `c_theta_t = 0.03`
- `sigma_theta_t = 2.0`

For this stage:

- implement advection
- implement diffusion with coefficient `sigma_theta_t * (nu + nu_t)`
- keep `P_theta_t = 0`
- keep `D_scf = 0`

But still implement the helper functions used by the full model so the path is ready.

### 2. Transport equation for `gamma`

`d(gamma)/dt + u_j d(gamma)/dx_j = P_gamma - E_gamma + (1/Re) d/dx_j [ (nu + nu_t/sigma_f) d(gamma)/dx_j ]`

Constants:

- `sigma_f = 1.0`

For this stage:

- implement advection
- implement diffusion with coefficient `(nu + nu_t / sigma_f)`
- keep `P_gamma = 0`
- keep `E_gamma = 0`

### 3. SA coupling to leave off for now

The final coupling is:

`P_nuTilde_modified = gamma * P_nuTilde`

Do **not** activate this yet. Keep SA exactly as it already works.

---

## Full model formulas that should already exist as helper functions

Even though source terms stay off in this stage, implement the functions now. This avoids redesign later and lets you unit test them independently.

## Core local quantities

### Velocity magnitude

`U = sqrt(u^2 + v^2 + w^2)`

Protect against division by zero with the same kind of small floor already used elsewhere in ADflow.

### Pressure-gradient parameter

`lambda_theta = (rho * theta^2 / mu) * dU/ds * Re`

with

`dU/ds = (u/U) dU/dx + (v/U) dU/dy + (w/U) dU/dz`

and

`dU/dx = (u du/dx + v dv/dx + w dw/dx) / U`

same idea for `dU/dy`, `dU/dz`.

### Strain-rate Reynolds number

`Re_S = (rho * d^2 * S / mu) * Re`

### Eddy-viscosity ratio

`R_T = mu_t / mu`

For SA this should use the same `mu_t` path already used elsewhere in the code.

### Boundary-layer proxies used in `F_theta_t`

`theta_BL = ReThetaTilde * mu / (rho * U) * (1/Re)`

`delta_BL = (15/2) * theta_BL`

`delta = (50 * d * Omega / U) * delta_BL`

### Wake factor used by SA-coupled form

`F_wake = exp( -Re_S / 1e6 )`

### F_theta_t

`F_theta_t = F_wake * exp( -(d/delta)^4 )`

### Re_theta correlation

Use free-stream turbulence intensity `Tu_inf` throughout the domain, not a local turbulence intensity.

If `Tu <= 1.3`:

`Re_theta = [1173.51 - 589.428 Tu + 0.2196 / Tu^2] * F(lambda_theta)`

If `Tu > 1.3`:

`Re_theta = 331.50 * (Tu - 0.5658)^(-0.671) * F(lambda_theta)`

### Original piecewise `F(lambda_theta)`

If `lambda_theta <= 0`:

`F = 1 - (-12.986 lambda_theta - 123.66 lambda_theta^2 - 405.689 lambda_theta^3) * exp(-(Tu/1.5)^1.5)`

If `lambda_theta > 0`:

`F = 1 + 0.275 * [1 - exp(-35 lambda_theta)] * exp(-Tu/0.5)`

### Crossflow helicity quantities

`Omega_streamwise = | U_hat · omega_vec |`

`H_crossflow = d * Omega_streamwise / U`

where `U_hat = V / |V|` and `omega_vec = curl(V)`.

### Original crossflow sink helper

`DeltaH_crossflow = H_crossflow * (1 + min(R_T, 0.4))`

`DeltaH_plus  = max(0.1066 - DeltaH_crossflow, 0)`

`DeltaH_minus = max(-(0.1066 - DeltaH_crossflow), 0)`

`f(DeltaH_plus)  = 6200 DeltaH_plus + 50000 DeltaH_plus^2`

`f(DeltaH_minus) = 75 tanh(DeltaH_minus / 0.0125)`

`Re_scf = -35.088 ln(h / theta_t) + 319.51 + f(DeltaH_plus) - f(DeltaH_minus)`

and the full sink term is:

`D_scf = (c_theta_t / t) * c_crossflow * min(Re_scf - ReThetaTilde, 0) * F_theta_t`

with `c_crossflow = 0.6`

This stays disabled in this stage.

---

## Smooth SA-sLM2015 helper functions to target for the final path

Implement these helper functions now as reusable utilities, but do not use them in active source terms yet unless needed by debug outputs.

### Smooth min/max operator

Use the exponential penalty form with proximity switch from the paper.

For two arguments `g1`, `g2`:

- if far from the kink, use ordinary `min` or `max`
- only near the kink use the smooth form

Recommended paper settings:

- `pswitch = 1e-15`
- smoothing parameter magnitude `|p| = 300`

This utility should be shared by all smooth approximations.

### Smooth `Fonset`

Final smooth form:

`Fonset1 = sqrt( (Re_S / (2.6 * Re_theta_c))^2 + R_T^2 )`

`Fonset = [tanh( 6 * (Fonset1 - 1.35) ) + 1] / 2`

### Smooth `Fturb`

Final smooth form:

`Fturb = (1 - Fonset) * exp(-R_T)`

### Smooth `Flength`

`Flength1 = exp( -3e-2 * (ReThetaTilde - 460) )`

`Flength = 44 - [44 - (0.50 - 3e-4 * (ReThetaTilde - 596))] / (1 + Flength1)^(1/6)`

### Smooth `Re_theta_c`

`Re_theta_c = 0.67 * ReThetaTilde + 24 * sin(ReThetaTilde / 240 + 0.5) + 14`

### Smooth `F(lambda_theta)`

`F1 = 1 + 0.275 * (1 - exp(-35 lambda_theta)) * exp(-Tu/0.5)`

`F2 = smooth_max(F1, 1)` with `p = +300`

`F3 = 1 - (-12.986 lambda_theta - 123.66 lambda_theta^2 - 405.689 lambda_theta^3) * exp(-(Tu/1.5)^1.5)`

`F = smooth_min(F2, F3)` with `p = -300`

---

## What to copy directly from the existing code

Do not design the plumbing from scratch. Mirror the existing working implementation.

### Copy from the existing SA-gamma-Retheta path

Use it as the primary template for:

- extra transported variables in the state vector
- indexing conventions
- restart read/write
- halo exchange of extra scalars
- residual storage and norm reporting
- scalar gradient reconstruction
- scalar diffusion assembly
- Jacobian seed and derivative path
- preconditioner / approximate Jacobian treatment
- wall, farfield, symmetry, overset, and block-interface BC hooks
- visualization and monitor outputs
- option parsing and runtime flags

### Copy from the SST transition path only where needed

Use it mainly for:

- correlation helper logic already present there
- any existing `gamma` and `ReTheta` boundary-condition treatment
- any debug output that already exists for transition triggers
- any special treatment of diffusion coefficients or clipping strategies for the transition variables

The shortest path is:

1. clone the already-working `gamma-ReTheta` transport infrastructure
2. swap in the SA-LM2015 / SA-sLM2015 helper formulas
3. keep source terms off
4. keep SA coupling off

---

## Variable scaling and solver guidance

This matters because ADflow uses ANK/NK and the paper is explicit that transition variables can ruin scaling and Jacobian conditioning.

### Column scaling targets from the paper

Use the paper values as the default normalization scales for the new columns:

- `nuTilde_max = 1e3`
- `gamma_max = 10`
- `ReThetaTilde_max = 1e4`

### Row scaling

Follow the same row-scaling structure already used by SA in the solver. The paper uses equation scaling that includes geometric factors and separate variable normalization. Keep the implementation consistent with the existing ADflow scaled system. Do not add a special one-off scaling path for transition if the current framework already supports turbulence/transition scalars generically.

### ANK/NK implications

Even with sources off, adding two scalars changes the Jacobian size, coupling blocks, and conditioning.

Recommendations:

- keep the first implementation maximally conservative
- reuse the same approximate Jacobian structure used by the existing transition model path
- do not try to improve the model and the solver simultaneously
- if the existing code supports first-order scalar discretization during early ANK, keep that behaviour for the new equations too
- keep the same pseudo-transient continuation machinery already used for SA / transition systems
- do not invent a new segregated strategy for this stage

When source terms are later activated, the paper strongly recommends keeping the system fully coupled and controlling stiffness through source-term-aware time stepping instead of deleting Jacobian couplings.

---

## Boundary conditions

The paper section provided here does **not** give a dedicated BC table for `gamma` and `ReThetaTilde`. So do not guess a novel BC treatment.

For this stage, the correct implementation strategy is:

- copy the boundary-condition treatment already used in ADflow for the existing `gamma-ReTheta` model
- keep the same SAT / ghost / interface handling pattern as the other transported scalars
- keep block-interface and halo treatment identical to the existing turbulence-scalar path

Practical guidance:

### At solid walls

Use the same wall BC definitions already used by the existing `gamma-ReTheta` implementation in ADflow. Do not improvise based only on the paper.

### At farfield / freestream

Again, copy the existing transition-model BC path. Typically `ReThetaTilde` needs a freestream-prescribed value from the correlation using `Tu_inf`, while `gamma` should follow the already-established freestream treatment in the code. Reuse exactly what the working model does.

### Symmetry / inviscid / interfaces

Treat both variables like existing transported scalars in ADflow.

Because your goal is code integration, not model redesign, matching the current ADflow BC conventions matters more than theoretical purity here.

---

## What the Jacobian path should include now

Even with source terms off, do not leave the equations half-implicit.

Implement derivatives for:

- advection of `gamma` and `ReThetaTilde`
- diffusion of `gamma` and `ReThetaTilde`
- dependence of diffusion coefficient on `nu_t`
- dependence on primitive/state variables exactly the same way existing scalar equations do

If helper functions are already coded now, their derivatives can remain unused until sources are enabled. But the transport-equation linearization must be complete.

Do not add temporary finite-difference Jacobians for these equations if the existing path is analytic or AD-based. Follow the existing implementation style.

---

## Recommended implementation order inside ADflow

1. **State vector and metadata**
   - add `gamma` and `ReThetaTilde` everywhere the state size is declared
   - add names, indices, restart support, I/O labels, residual labels

2. **Initialization**
   - initialize using the same conventions as the working transition model path
   - safe defaults for the first stage:
     - `gamma` near laminar value used by existing code path
     - `ReThetaTilde` from freestream correlation based on `Tu_inf`

3. **Halo and block communication**
   - extend halo pack/unpack and interface exchange to include the two new variables

4. **Residual assembly**
   - add advection and diffusion operators for both scalars
   - set all source contributions to zero

5. **Jacobian / preconditioner**
   - include the two equations in the same linearization framework used by the current transition model path

6. **Boundary conditions**
   - copy from the existing `gamma-ReTheta` implementation

7. **Output / debug hooks**
   - volume and surface outputs for the new fields and helper diagnostics

8. **Verification with sources off**
   - run transport-only tests before activating any physics trigger

---

## Debug outputs to add now

Add these as volume outputs if cheap enough, even if some are not yet used in active residual terms:

- `gamma`
- `ReThetaTilde`
- `mu_t`
- `R_T`
- `Re_S`
- `U`
- `lambda_theta`
- `F_lambdaTheta`
- `F_theta_t`
- `F_wake`
- `H_crossflow`
- `Omega_streamwise`
- `Fonset`
- `Fturb`
- `Flength`
- `Re_theta_c`
- `Re_theta`

At this stage, these are mainly to confirm that helper quantities are computed consistently and stay finite.

Also add residual norms for the two new equations separately.

---

## Minimum validation for this stage

Before sources are turned on, validate only code integration and numerics.

### Required checks

1. Solver runs with the enlarged state vector.
2. Restarts work.
3. Halo exchange and block interfaces work in parallel.
4. Residuals for the new equations converge for a passive transport-diffusion case.
5. Jacobian-based ANK/NK still behaves normally.
6. No NaNs from helper functions in laminar zones, stagnation regions, or near-zero velocity regions.
7. Output fields are sensible and smooth.

### Good initial tests

- a simple attached 2D case with uniform initialization
- a coarse 3D wing case just to exercise parallel block coupling
- a manufactured-style passive scalar transport test if easy to add

Do not evaluate transition prediction yet. This stage is only about plumbing and solver stability.

---

## Numerical cautions

- Protect all divisions by `U`, `mu`, `theta_t`, and any near-zero quantity.
- Reuse existing ADflow floors, not ad hoc new ones, unless absolutely needed.
- Keep helper functions side-effect free so they can be reused later in residual and Jacobian code.
- Where the paper uses smooth min/max operators, keep both the raw and smooth versions available behind a clean interface.
- Do not scatter constants through the code. Centralize them like other turbulence-model constants.
- Avoid branching differences between residual and Jacobian paths.
- The final model is intended for deep convergence in fully-coupled Newton-Krylov. So code structure now should already be compatible with exact same-state residual/Jacobian evaluation later.

---

## What to postpone deliberately

Do **not** do these yet:

- activate `P_gamma`
- activate `E_gamma`
- activate `P_theta_t`
- activate `D_scf`
- multiply SA production by `gamma`
- add source-term time-step restriction
- add transition-specific update damping
- tune model constants
- judge physical transition locations

All of that belongs to the next phase, once the transport-only integration is validated.

---

## Short implementation brief for the coding agent

Implement SA-sLM2015 in ADflow by cloning the existing working `gamma-ReTheta` transport infrastructure and adapting it to the SA-based model of Piotrowski and Zingg (2020). Add two transported variables, `gamma` and `ReThetaTilde`, through the full solver stack: state vector, restart, residual assembly, Jacobian, preconditioner path, communication, BCs, output, and diagnostics. Implement the advection and diffusion operators only. Implement all helper correlations and smooth helper functions from the paper, but keep all transition source terms disabled and keep SA production uncoupled for this first stage. Boundary conditions should be copied from the existing ADflow `gamma-ReTheta` implementation, not re-invented from scratch. The goal is a stable transport-only integration that behaves correctly in ANK/NK and is ready for later activation of the source terms and SA coupling.

---

## Reference

Piotrowski, M. G. H., and Zingg, D. W., “Smooth Local Correlation-Based Transition Model for the Spalart-Allmaras Turbulence Model,” *AIAA Journal*, 2020, doi:10.2514/1.J059784.
