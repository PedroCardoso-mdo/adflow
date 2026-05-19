# New Options Added for SA-γ-Re̅θt Transition Model

These options do not exist in upstream ADflow. All were added on this branch.

---

## Python Options

### New transition-specific options

| Option | Type | Default | What it does |
|---|---|---|---|
| `"transitionFirstOrderUpwind"` | bool | `True` | First-order upwind for γ and Re̅θt convection. More dissipative but more robust. |
| `"transitionSrcDtRestrict"` | bool | `True` | Enable source-term dt restriction (P&Z Eq. 59). Caps λ_source × dt ≤ 0.9. |
| `"transitionSrcDtLimit"` | float | `0.9` | Threshold for source-term dt restriction (λ_source × dt ≤ this value). |
| `"srcDtDeactivateIters"` | int | `5` | Deactivate source-dt restriction after N consecutive ANK iterations without backtracking (P&Z §IV.B.3). Set to 0 to never deactivate. |
| `"TurbDADICoupled"` | str | `"full"` | DADI coupling mode: `"decoupled"` (3 scalar solves), `"transition"` (SA alone + γ-Re̅θt 2×2), `"full"` (3×3 block). |
| `"turbResScale"` | list/None | `None` (auto) | Residual scaling per equation. Auto-set to `[10000, 10, 10000]` for this model. Override to tune convergence balance. |
| `"transitionDampTheta"` | float | `0.99` | Back-off factor for iterative γ/Re̅θt update damping in DD-ADI (P&Z §3). |
| `"transitionDampMaxIter"` | int | `40` | Max back-off iterations for γ/Re̅θt bounds enforcement in DD-ADI. |

### Turb-ANK KSP physicality options (transition-specific)

| Option | Type | Default | What it does |
|---|---|---|---|
| `"ANKPhysicalLSTolReTheta"` | float | `0.99` | Relative physicality tolerance for Re̅θt in Turb-ANK (replaces `ANKPhysicalLSTolTurb` for Re̅θt). |
| `"omegaMinGamma"` | float | `0.05` | Minimum step factor floor for γ. Prevents collapse in laminar regions where γ→0. |

### Existing ADflow options relevant to turbulent solver path

| Option | Type | Default | What it does |
|---|---|---|---|
| `"ANKUseTurbDADI"` | bool | `True` | `True` = DADI for turbulence. `False` = Turb-ANK KSP (Newton-Krylov). |
| `"ANKNSubiterTurb"` | int | `1` | Inner turbulence iterations per outer ANK step. |
| `"ANKTurbCFLScale"` | float | `1.0` | CFL multiplier for turb equations relative to flow. |
| `"ANKTurbKSPDebug"` | bool | `False` | Print linear residual, KSP iters, step size each Turb-ANK iteration. |
| `"ANKPhysicalLSTolTurb"` | float | `0.99` | Physicality line-search tolerance for ν̃ in Turb-ANK (γ uses absolute bounds instead). |

---

## Examples

### 1. Robust startup (recommended defaults)

```python
solverOptions = {
    # Transition-specific (new)
    "transitionFirstOrderUpwind": True,      # robust convection for γ, Re̅θt
    "transitionSrcDtRestrict": True,         # source limiting ON
    "srcDtDeactivateIters": 5,               # deactivate after 5 clean ANK iters
    "TurbDADICoupled": "full",               # 3×3 coupled DADI
    # turbResScale auto-set to [10000, 10, 10000]

    # Solver path (existing ADflow)
    "ANKUseTurbDADI": True,                  # use DADI for turbulence
}
```

### 2. Accuracy run (restarting from converged solution)

```python
solverOptions = {
    "transitionFirstOrderUpwind": False,  # second-order convection (sharper transition front)
    "transitionSrcDtRestrict": False,     # no source limiting (solution already stable)
    "TurbDADICoupled": "full",
    "ANKUseTurbDADI": True,
}
```

### 3. Debugging convergence — try decoupled DADI

```python
solverOptions = {
    "transitionFirstOrderUpwind": True,
    "transitionSrcDtRestrict": True,
    "TurbDADICoupled": "decoupled",       # simplest: 3 independent scalar solves
    "ANKUseTurbDADI": True,
}
```

### 4. Debugging convergence — try partial coupling

```python
solverOptions = {
    "transitionFirstOrderUpwind": True,
    "transitionSrcDtRestrict": True,
    "TurbDADICoupled": "transition",      # SA alone, γ-Re̅θt coupled as 2×2 block
    "ANKUseTurbDADI": True,
}
```

### 5. Turb-ANK KSP path (Newton-Krylov for turbulence)

```python
solverOptions = {
    "transitionFirstOrderUpwind": True,
    "transitionSrcDtRestrict": True,
    # TurbDADICoupled ignored when ANKUseTurbDADI=False

    "ANKUseTurbDADI": False,              # switch to Turb-ANK KSP
    "ANKNSubiterTurb": 3,                 # more inner Newton iters
    "ANKTurbCFLScale": 0.5,              # lower CFL for turb if unstable
    "ANKTurbKSPDebug": True,             # print convergence info
}
```

### 6. Custom residual scaling

```python
solverOptions = {
    "turbResScale": [5000.0, 1.0, 5000.0],  # [ν̃, γ, Re̅θt] — lower γ scaling
    "TurbDADICoupled": "full",
    "ANKUseTurbDADI": True,
}
```

---

## Solver Path Summary

```
ANKUseTurbDADI = True ──┬── TurbDADICoupled = "decoupled"   → 3 scalar solves
                        ├── TurbDADICoupled = "transition"  → SA scalar + γ-Re̅θt 2×2 block
                        └── TurbDADICoupled = "full"        → 3×3 coupled block (default)

ANKUseTurbDADI = False ──── Turb-ANK KSP (Newton-Krylov, GMRES)
```

---

## Source-Term Eigenvalue Control (P&Z 2020, Eq. 59)

The source-term dt restriction prevents unbounded solution updates by limiting:

```
λ_source × Δt ≤ transitionSrcDtLimit  (default 0.9)
```

where `λ_source` is the **largest positive eigenvalue** of the 3×3 source-term Jacobian:

```
            ⎡ ∂S_ν̃/∂ν̃      ∂S_ν̃/∂γ      0           ⎤
A_source =  ⎢ ∂S_γ/∂ν̃      ∂S_γ/∂γ      ∂S_γ/∂Re̅θt   ⎥
            ⎣ 0            0            ∂S_Re̅θt/∂Re̅θt ⎦
```

### Key points

1. **Block-triangular structure**: A13=A31=A32=0 (P&Z §7.1), so eigenvalues are computed exactly without a cubic solver:
   - λ₃ = A33 (Re̅θt diagonal)
   - λ₁,₂ from 2×2 block [A11,A12; A21,A22] via quadratic formula
   - `λ_source = max(0, λ₁, λ₂, λ₃)`

2. **Independent of `TurbDADICoupled` mode** — coupling mode only affects how DADI solves the system, not eigenvalue computation.

3. **Auto-deactivation**: After `srcDtDeactivateIters` consecutive ANK iterations without backtracking, the restriction turns off. Reactivates on backtracking.

### Examples

#### 7. Conservative eigenvalue control (stiff cases)

```python
solverOptions = {
    "transitionSrcDtRestrict": True,
    "transitionSrcDtLimit": 0.7,             # stricter than default 0.9
    "srcDtDeactivateIters": 10,              # wait longer before deactivating
}
```

#### 8. Debug eigenvalue issues

```python
solverOptions = {
    "transitionSrcDtRestrict": True,
    "srcDtDeactivateIters": 0,               # NEVER deactivate — always restrict
    "ANKTurbKSPDebug": True,                 # print iteration info
}
```

#### 9. Disable source-dt restriction entirely

```python
solverOptions = {
    "transitionSrcDtRestrict": False,        # no eigenvalue computation, no restriction
}
```

---

## γ Physicality Check in Turb-ANK (Redesigned)

In the Turb-ANK KSP path (`ANKUseTurbDADI = False`), γ uses **absolute bound enforcement** instead of a relative tolerance:

- **Full step allowed** if result stays in [gammaLo, gammaHi] (~[1e-10, 2.0])
- **Only reduced** when full step would violate bounds
- **`omegaMinGamma`** (default 0.05) prevents step collapse in laminar regions where γ→0

This differs from ν̃ and Re̅θt which use relative tolerances (`ANKPhysicalLSTolTurb`, `ANKPhysicalLSTolReTheta`).

**Why**: In laminar flow, γ≈0. The old relative check `ratio = γ/update × tol` collapses to near-zero, killing the transition front before it can develop.

---

## Internal State (not user-settable)

These are managed automatically by the solver when `transitionSrcDtRestrict = True`:

- `srcDtRestrictActive`: starts `True`, flips to `False` after `srcDtDeactivateIters` consecutive no-backtrack ANK iterations. Resets to `True` on backtrack.
- `noBacktrackCount`: counter driving the above.
