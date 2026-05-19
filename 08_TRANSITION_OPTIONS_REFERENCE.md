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
| `"transitionSrcDtEigMode"` | str | `"eigenvalue"` | How to compute λ_source: `"gershgorin"` (signed upper bound, cheap) or `"eigenvalue"` (exact 3×3 cubic formula). |
| `"srcDtDeactivateIters"` | int | `5` | Deactivate source-dt restriction after N consecutive ANK iterations without backtracking (P&Z §IV.B.3). Set to 0 to never deactivate. |
| `"TurbDADICoupled"` | str | `"full"` | DADI coupling mode: `"decoupled"` (3 scalar solves), `"transition"` (SA alone + γ-Re̅θt 2×2), `"full"` (3×3 block). |
| `"turbResScale"` | list/None | `None` (auto) | Residual scaling per equation. Auto-set to `[10000, 10, 10000]` for this model. Override to tune convergence balance. |

### Existing ADflow options relevant to turbulent solver path

| Option | Type | Default | What it does |
|---|---|---|---|
| `"ANKUseTurbDADI"` | bool | `True` | `True` = DADI for turbulence. `False` = Turb-ANK KSP (Newton-Krylov). |
| `"ANKNSubiterTurb"` | int | `1` | Inner turbulence iterations per outer ANK step. |
| `"ANKTurbCFLScale"` | float | `1.0` | CFL multiplier for turb equations relative to flow. |
| `"ANKTurbKSPDebug"` | bool | `False` | Print linear residual, KSP iters, step size each Turb-ANK iteration. |
| `"ANKPhysicalLSTolTurb"` | float | `0.99` | Physicality line-search tolerance for Turb-ANK. |

---

## Examples

### 1. Robust startup (recommended defaults)

```python
solverOptions = {
    # Transition-specific (new)
    "transitionFirstOrderUpwind": True,      # robust convection for γ, Re̅θt
    "transitionSrcDtRestrict": True,         # source limiting ON
    "transitionSrcDtEigMode": "eigenvalue",  # exact 3×3 eigenvalue (or "gershgorin")
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
            ⎡ ∂S_ν̃/∂ν̃      ∂S_ν̃/∂γ      ∂S_ν̃/∂Re̅θt  ⎤
A_source =  ⎢ ∂S_γ/∂ν̃      ∂S_γ/∂γ      ∂S_γ/∂Re̅θt   ⎥
            ⎣ ∂S_Re̅θt/∂ν̃  ∂S_Re̅θt/∂γ  ∂S_Re̅θt/∂Re̅θt ⎦
```

### Key points

1. **Always computed on full 3×3 matrix** — independent of `TurbDADICoupled` mode. The coupling mode only affects how DDADI solves the system, not the eigenvalue stability check.

2. **Two computation modes**:
   - `"eigenvalue"` (default): exact eigenvalues via cubic formula (Cardano). More accurate.
   - `"gershgorin"`: signed Gershgorin bound `max(0, max_i[A_ii + Σ|A_ij|])`. Cheaper, always ≥ true eigenvalue.

3. **Auto-deactivation**: After `srcDtDeactivateIters` consecutive ANK iterations without backtracking, the restriction turns off (eigenvalues may stay large but solution is stable). Reactivates on backtracking.

### Examples

#### 7. Conservative eigenvalue control (stiff cases)

```python
solverOptions = {
    "transitionSrcDtRestrict": True,
    "transitionSrcDtEigMode": "eigenvalue",  # exact (default)
    "transitionSrcDtLimit": 0.7,             # stricter than default 0.9
    "srcDtDeactivateIters": 10,              # wait longer before deactivating
}
```

#### 8. Fast eigenvalue (large cases, prioritize speed)

```python
solverOptions = {
    "transitionSrcDtRestrict": True,
    "transitionSrcDtEigMode": "gershgorin",  # cheaper upper bound
    "srcDtDeactivateIters": 3,               # deactivate quickly
}
```

#### 9. Debug eigenvalue issues

```python
solverOptions = {
    "transitionSrcDtRestrict": True,
    "transitionSrcDtEigMode": "eigenvalue",
    "srcDtDeactivateIters": 0,               # NEVER deactivate — always restrict
    "ANKTurbKSPDebug": True,                 # print iteration info
}
```

#### 10. Disable source-dt restriction entirely

```python
solverOptions = {
    "transitionSrcDtRestrict": False,        # no eigenvalue computation, no restriction
}
```

---

## Internal State (not user-settable)

These are managed automatically by the solver when `transitionSrcDtRestrict = True`:

- `srcDtRestrictActive`: starts `True`, flips to `False` after `srcDtDeactivateIters` consecutive no-backtrack ANK iterations. Resets to `True` on backtrack.
- `noBacktrackCount`: counter driving the above.
