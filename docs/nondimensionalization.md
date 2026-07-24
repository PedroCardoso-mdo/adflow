# Non-dimensionalization — ADflow convention (relevant to SA-BCM)

ADflow-wide convention, not SA-BCM-specific, but SA-BCM's transition-trigger terms mix velocity
(vorticity), length (wall distance), and viscosity-ratio quantities, so get this right before
touching any of them.

- **Units are p-ρ non-dimensional, not velocity-based.** Velocity normalizes to `M·√γ` (freestream
  Mach times √(ratio of specific heats)), **not** to 1.
- **Viscosities are ratios to μ∞** (`rlv`, `rev` in the code), not absolute values.
- **`1/Re` is NOT absorbed into viscosity** — it appears explicitly wherever a viscous term is
  formed.

## Where this matters in SA-BCM

- `Rv = ρ d² S / μ` (manuscript Eq. 5, `tterm1`'s numerator): `S` is vorticity magnitude in the
  code's non-dimensional velocity units (`M·√γ` scale), `d` is non-dimensional wall distance,
  `μ` is the non-dimensional (ratio-to-μ∞) laminar viscosity `rlv`. Do not treat `S` as if it
  were dimensional or velocity-normalized-to-1 — a term that looks right in a "textbook"
  non-dimensionalization will be off by a Mach/γ factor here.
- `term2 = νT / (χ2 ν)` (manuscript Eq. 6): both `νT` and `ν` are already viscosity ratios in
  ADflow's convention (`rev`/`rlv`-style quantities), so this ratio is dimensionally consistent
  as written — no additional `1/Re` factor is needed here, but double-check this assumption
  against `sa.F90`'s actual variable definitions (`rev`, `rlv`, not raw `mu`) before editing.

If in doubt, cross-reference the exact variable (`rlv`, `rev`, `ss` for vorticity, etc.) against
its declaration/usage elsewhere in `sa.F90` rather than assuming the manuscript's symbol maps
1:1 to a dimensional quantity.
