# SA-noft2-Gamma-Retheta (SA-GR) derivative tests

Verifies the SA-GR AD derivatives (Tapenade forward `_d`, reverse `_b`,
reverse-fast `_fast_b`) against complex step, about a converged state on the
standard ADflow tutorial wing (mach 0.15, α 1.8, 8 states). The stage ladder,
invariants and failure table are in
`docs/VERIFICATION/VERIF_00_three_stage_verification.md`; run it with the
`adflow-derivative-verification` user skill. The SA-noft2-Gamma suite
(`reg_sgamma.py`, `*_sgamma.py`, `run_sgamma_tests.sh`) is a copy of this one
with nw = 7.

## File tree

```
tests/reg_tests/
├── reg_sagr.py                      # case config + block/assert helpers (everything SA-GR-specific)
├── run_sagr_tests.sh                # driver: all|real|cs|adjoint|blockette|train|genw
├── test_jacVecProdFWD_sagr.py       # forward mode: AD vs ref / FD / CS
├── test_jacVecProdBWDFast_sagr.py   # _b vs _fast_b + dot products
├── test_adjoint_sagr.py             # adjoint totals vs ref + complex-step totals
├── test_blockette_sagr.py           # block vs blockette residual equivalence
├── dev/generate_sagr_restart.py     # produces the restart CGNS (genw)
├── dev/diag_full_derivatives.py     # full per-DV adjoint-vs-CS table, raw solver output
└── refs/{jacvecfwd,jacvecbwd,adjoint}_sagr_tut_wing.json   # trained references
```

| SA-GR file                       | Mirrors (upstream SA)                                            | Verifies |
|----------------------------------|------------------------------------------------------------------|----------|
| `test_jacVecProdFWD_sagr.py`     | `test_jacVecProdFWD.py`                                          | `computeJacobianVectorProductFwd` (`outputForward/*_d.f90`) |
| `test_jacVecProdBWDFast_sagr.py` | `test_jacVecProdBWDFast.py` + dot products of `test_functionals.py` | `_fast_b` vs `_b`, fwd↔rev transpose consistency |
| `test_adjoint_sagr.py`           | `test_adjoint.py`                                                | adjoint totals (`outputReverse/*_b.f90`) vs complex step |
| `test_blockette_sagr.py`         | —                                                                | block-loop vs blockette SA-GR residual |
| `reg_sagr.py`                    | `reg_default_options.py` + `reg_aeroproblems.py` + `reg_test_utils.py` | fixtures + assert helpers |

## What `reg_sagr.py` holds

- **Case config:** grid/restart/FFD paths, `sagrBaseOptions` (ANK → CANK → NK,
  no wall functions, no ft2, crossflow off), the AeroProblem
  `ap_sagr_tut_wing`, and the aero DVs `alpha, mach, P, T` (mach/P/T drive
  `uInf`/`muInf` and the farfield γ/Re̅θt). The test conditions must match the
  ones the restart was converged at.
- **Block helpers:** `getStateBlocks` / `maskStateVector` split a flattened
  state/residual vector (variable-fastest, `w(i,j,k,1:nw)`) into
  `meanflow` / `nuTilde` / `gamma` / `reThetat` blocks.
- **Assert helpers:**
  - `assert_coupling_blocks_allclose` — norms of all 16 `dR[row]/dw[col]`
    blocks, including the structurally-zero `dR[nuTilde]/dw[reThetat]` as a
    regression guard. Re-seats the real state before each CS column.
  - `assert_transition_xdvdot_allclose` — transition-row norms of
    `dR/d(mach, P, T, alpha)`: the farfield-BC / vorticity-limiter check.
  - `assert_bwdfast_blocks_allclose` — `_b` vs `_fast_b` seeded one
    equation-row block at a time.
  - `assert_coupling_dot_products_allclose` — blockwise transpose tests.
  - `assert_fd_allclose_hsweep` — FD with a step sweep (FD classes are
    `@expectedFailure`; CS is the enforced truth).

## Test classes

- `TestJacVecFwdSAGR` (real, ref-based), `TestJacVecFwdSAGRFD` (real, FD,
  expected-fail), `TestJacVecFwdSAGRCS` (complex, `cmplx_test_*`: records CS
  under the same keys the real class trained, so it compares AD vs CS
  key by key; `cmplx_test_xDvDot_transition_rows` is the decisive check on the
  `uInf`/`muInf` limiter path).
- `TestJacVecBWDFastSAGR` (full and row-block-seeded `_b` vs `_fast_b`, atol
  1e-16, repeated-call determinism), `TestDotProductsSAGR`.
- `TestAdjointSAGR` (real): residuals + `evalFunctionsSens` totals for
  `cl, cd, cmz, drag` vs ref, with the twist/span/shape DVGeo of the SA
  `test_adjoint.py`. Needs `adjointMaxIter = 3000` (set in `reg_sagr.py`):
  the 8-state adjoint takes ~600 iterations to reach 1e-14.
- `TestCmplxStepSAGR` (complex): re-converges the complex solve with a
  1e-40j perturbation per DV and compares against the adjoint refs at
  `rtol = atol = 5e-8`; `mach` and `drag` are reported, not asserted.

## `dev/diag_full_derivatives.py`

Not a testflo test: the full `df/dx` check with the flow/adjoint output
visible, over all aero and geometric DVs including every `shape` component.
`--mode adjoint` (real build) writes a JSON that `--mode cs` (complex build)
reads, printing an abs/rel error table plus the L2 reached by each CS
re-converge, so a loose comparison is attributable to a shallow complex solve.
Narrow the (expensive) sweep with `--shape a,b,c`, `--ncycles`,
`--skip-aero/--skip-geom/--skip-span/--skip-twist`. Run CS at the same `-np`
the refs were trained at (2). For parallel jobs on distinct cores use
mpirun's `--cpu-set`, not a second `--bind-to core`.

## Single test

`testflo test_jacVecProdFWD_sagr.py:TestJacVecFwdSAGR_0_sagr_tut_wing.test_coupling_blocks`
(`testflo --dryrun` lists the exact signatures).
