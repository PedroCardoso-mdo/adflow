# Adjoint / AD maintenance on this branch

The adjoint of every turbulence model here (SA, SA-BCM, SA-GR, SA-sγ) is
**fully Tapenade-differentiated** — there is no hand-written derivative math.
It is verified to complex-step precision by the ladder in
[`VERIF_00_three_stage_verification.md`](../VERIFICATION/VERIF_00_three_stage_verification.md).
This file is what you need to keep it that way when the primal changes.

## 1. Regenerating (after any change to differentiated code)

```bash
./AD_I.sh        # make -f Makefile_tapenade ad_forward ad_reverse ad_reverse_fast, then make
pip install . --no-deps   # (the /build skill does rm -rf build + install)
```

Output lands in `src/adjoint/output{Forward,Reverse,ReverseFast}/` — never
hand-edit it (CLAUDE.md rule 6). The `autoEdit*.py` post-processors run as part
of each `ad_*` target. Run the `ad_*` steps in the foreground and rebuild the
complex lib too before the CS stage. State **`TAPENADE NEEDED`** in the wrap-up
whenever a differentiated routine or its module variables changed.

## 2. Differentiated heads per model (`src/adjoint/Makefile_tapenade`)

| Model | fullRoutines (forward + reverse) | stateOnly (reverse fast) |
|---|---|---|
| SA / SA-BCM | `sa%saSource`, `sa%saViscous`, `sa%saResScale` | `sa%saSource(w, rlv)`, `saViscous`, `saResScale` |
| SA-GR | `saGammaRetheta%Source(…, uInf, muInf)`, `%Viscous`, `%ResScale` | `Source(w, rlv)`, `Viscous`, `ResScale` |
| SA-sγ | `saGamma%sgSource(…, uInf, muInf)`, `%sgViscous`, `%sgResScale` | `sgSource(w, rlv)`, `sgViscous`, `sgResScale` |
| all | `turbutils%turbAdvection`, `turbutils%saEddyViscosity`, `turbBCRoutines%*` | same |

SA-BCM has no head of its own: its intermittency lives inside `saSource`
(`sa.F90`, gated by `useSABCM`), so the SA heads carry it. The SA-GR/SA-sγ
correlations and `smoothMinMax` (`turbUtils.F90`) are differentiated through
the Source heads.

**Why `uInf`, `muInf` are active inputs:** the P&Z vorticity limiter
`vortLim = U_ref·√(U_ref/(muInf·l))/20` makes the residual depend on the
freestream state. Without them in the head, Tapenade emits `vortlimd = 0` and
`dR/dMach`-type partials miss the limiter in capped cells (dR/dw is
unaffected). Derivation: [`ADFLOW_08`](ADFLOW_08_nondimensionalization.md) §5.

## 3. Master-routine dispatch (`src/adjoint/masterRoutines.F90`)

Each of the five sweeps — `master`, `master_d`, `master_b`,
`master_state_b`, `block_res_state_d` — has one `select case (turbModel)`
calling the model's Source → `turbAdvection` → Viscous → ResScale (reverse
order in `_b`). The advection work array is model-sized: `qq` (SA),
`qqGR(…,3,3)` (SA-GR), `qqSG(…,2,2)` (SA-sγ). A new model needs a case in
**all five** plus its three heads, or its adjoint silently drops terms.

## 4. What is deliberately NOT differentiated

- Everything behind `#ifndef USE_TAPENADE` in `saGammaRetheta.F90` /
  `saGamma.F90`: the DADI LHS / `qq` blocks, `tdia3x3` / `tdia2x2`
  (`autoEdit*` re-inject them as `external`).
- `evalSrcJacBlock` and `computeSrcLambda` / `computeSrcLambdaSaGamma` —
  primal-only (PC and Eq. 59 time-step), not residual.
- `transitionDebug` — output only, no derivative.

Changes confined to these do **not** need a Tapenade regen.

## 5. Tapenade / autoEdit gotchas

- **No in-place `x = smoothMinMax(x, …)`** in differentiated code:
  `autoEditReverseFast.py` strips push/pop pairs, so an overwritten active
  variable loses its saved value in `_fast_b` (Stage 2 of the ladder fails,
  wrong sign). Use distinct targets (`raw` → `clamped` → `local`).
- **No local variable names ending in `_c`** — they collide with Tapenade's
  generated names.
- **Copy module variables into locals** before passing them to a
  differentiated call.
- **A branching routine must `use constants` unrestricted** (no `only:`
  list), or Tapenade's generated code misses symbols.
- **`.pyf` wiring:** an option whose Fortran variable is missing from the
  module block in `src/f2py/adflow.pyf` is a **silent no-op** — grep the
  `.pyf` before trusting a new option.

## 6. Adjoint options specific to this branch

- `frozenTransition` (default False): drops γ and Re̅θt from the
  reverse-mode products (SA-GR only; `master_b`, `master_state_b`, and the
  PC in `adjointUtils.F90` stay consistent). For A/B studies, not production.
- `storePsiHistory` (+ `psiHistoryStep`, `psiHistoryMax`): buffers adjoint
  iterates for convergence diagnostics.
- Transition models (`"gamma"` in `turbulenceModel`) get
  `globalPreconditioner = "field split"` automatically unless the user sets
  it; see [`CORE_01`](../CORE_BASE/CORE_01_architecture.md) for the adjoint
  defaults.

## 7. Verification after a change

Rebuild real + complex, then run the ladder (dot product → `_b` vs `_fast_b`
→ AD vs CS partials → total adjoint vs CS) —
[`VERIF_00`](../VERIFICATION/VERIF_00_three_stage_verification.md). A Stage-2
failure points at §5's in-place rule; a CS failure with Stages 1–2 passing
means the generated code is out of sync with the primal (regen) or an active
input is missing from a head.
