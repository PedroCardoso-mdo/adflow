# Turbulence Coupling Modes — Annotation

This file documents the turbulence coupling options available in the
SA-gamma-ReTheta transition model. All modes are runtime-selectable
from Python.

---

## DADI Coupling (`TurbDADICoupled`)

Controls how the 3x3 source Jacobian block is used in the DD-ADI solver.

| Python value   | Fortran int | Description                                           |
|----------------|-------------|-------------------------------------------------------|
| `"decoupled"`  | 0           | Diagonal only. SA, gamma, ReTheta solved independently. All off-diagonal qq entries zeroed. |
| `"transition"` | 1           | SA decoupled from transition. gamma-ReTheta coupled via qq(2,3) and qq(3,2). SA rows/cols (qq(1,2), qq(1,3), qq(2,1), qq(3,1)) zeroed. |
| `"full"`       | 2 (default) | Full 3x3 block coupling. All off-diagonal entries active. |

### Off-diagonal source Jacobian entries (set in `Source`, used by DADI)

| Entry    | Derivative             | Dependency chain                                      | Method   |
|----------|------------------------|-------------------------------------------------------|----------|
| qq(1,2)  | -dS_nu/dgamma          | gammaForSA in SA production (term1, term2_prod)       | Analytic |
| qq(1,3)  | -dS_nu/dReThetaTilde   | None — SA source independent of ReTheta               | Zero     |
| qq(2,1)  | -dS_gamma/dnu_tilde    | rTurb(nu_tilde) in fOnset (pGamma) and fTurb (eGamma) | Analytic |
| qq(2,3)  | -dS_gamma/dReThetaTilde | reThetaC and fLength correlations in pGamma           | FD       |
| qq(3,1)  | -dS_retheta/dnu_tilde  | None — pReTheta independent of nu_tilde               | Zero     |
| qq(3,2)  | -dS_retheta/dgamma     | None — fThetaT has no gamma branch in current code    | Zero     |

### Active off-diags per mode

| Mode         | qq(1,2) | qq(1,3) | qq(2,1) | qq(2,3) | qq(3,1) | qq(3,2) |
|--------------|---------|---------|---------|---------|---------|---------|
| `decoupled`  | 0       | 0       | 0       | 0       | 0       | 0       |
| `transition` | 0       | 0       | 0       | active  | 0       | 0       |
| `full`       | active  | 0       | active  | active  | 0       | active  |

Note: qq(1,3) and qq(3,1) are always zero because the SA source does not
depend on ReTheta and vice versa.

---

## Turb Sub-Solver Selection (`ANKUseTurbDADI`)

Controls which sub-solver handles turbulence in the decoupled ANK path.

| Python value | Description                                                        |
|--------------|--------------------------------------------------------------------|
| `True`       | Use DD-ADI (with coupling level set by `TurbDADICoupled`)          |
| `False`      | Use Turb-ANK KSP (Newton-Krylov sub-solver for turbulence block)  |

Only applies when flow and turbulence are decoupled (`ANKCoupledSwitchTol`
not yet triggered).

---

## Flow-Turb Coupling (`ANKCoupledSwitchTol`)

Controls whether flow and turbulence are solved in one Newton system.

| Value        | Effect                                                             |
|--------------|--------------------------------------------------------------------|
| `0.0`        | Never switch to coupled — flow and turb always decoupled           |
| `1e-16`      | Switch to coupled only at deep convergence (default)               |
| `1.0`        | Always coupled from the start                                      |

When coupled, all 8 equations (5 NS + nu_tilde + gamma + ReTheta) are in
one Krylov system. The `TurbDADICoupled` and `ANKUseTurbDADI` options do
not apply in coupled mode.

---

## Summary: The Four Solver Paths

| Path | Flow-Turb  | Turb Solver    | DADI Coupling      | Python config                                        |
|------|-----------|----------------|---------------------|------------------------------------------------------|
| 1    | Decoupled | DADI diagonal  | `"decoupled"`       | `ANKUseTurbDADI=True, TurbDADICoupled="decoupled"`  |
| 1b   | Decoupled | DADI transition| `"transition"`      | `ANKUseTurbDADI=True, TurbDADICoupled="transition"` |
| 2    | Decoupled | DADI full 3x3  | `"full"`            | `ANKUseTurbDADI=True, TurbDADICoupled="full"`       |
| 3    | Decoupled | Turb-ANK KSP   | N/A                 | `ANKUseTurbDADI=False`                               |
| 4    | Coupled   | N/A (all in NK)| N/A                 | `ANKCoupledSwitchTol=1.0`                            |
