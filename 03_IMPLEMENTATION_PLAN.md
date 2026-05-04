# Implementation Plan — Per-Task Specs

> Each task below is self-contained. To work a task, jump to its `## T-XX —`
> heading and read only that block plus the files listed in `Context:`.
> Do NOT read the whole document.

---

## How to Read a Task Block

```
## T-XX — title
Goal:        one sentence
AD impact:   yes/no  (yes ⇒ end with STATUS: TAPENADE NEEDED)
Context:     files Claude Code should load before editing
Files:       files Claude Code will edit
Change:      what to do, line-by-line if needed
Compile:     `cd build && make -j` (always)
Done when:   exact bullet criteria
Status row:  which row in 02_IMPLEMENTATION_STATUS.md to flip
End-of-turn: STATUS: READY TO RUN  |  STATUS: TAPENADE NEEDED
```

---

## A1 — Smoke baseline: SA only, transition off

Goal: Confirm the current branch compiles and runs SA without transition,
matching upstream behavior. This is the "we did not break anything" gate
before touching transition code.

AD impact:   no
Context:     `05_TESTING_AND_DEBUG.md` §1
Files:       (none — read-only verification)
Change:
- Run `cd build && make -j`. Fix any compile errors that exist on the branch
  today. Do NOT introduce new code.
- Confirm `python/pyADflow.py` exposes `useTransitionModel` (default False).
  If not, add the option line with default `False` only.
Compile:     standard
Done when:
- `make` exits 0.
- `useTransitionModel` defaults to False; with it False, the binary is
  available for an SA-only run (user runs the smoke script).
Status row:  flip A1 to ✅
End-of-turn: STATUS: READY TO RUN

---

## A2 — Output every non-intermediate transition variable to volume CGNS

Goal: Wire the existing `transitionDebug(:,:,:,:)` block array into ADflow's
volume/surface CGNS output, the same way `eddy` and `cf` are wired. Output
every meaningful state, source, and Jacobian quantity computed in the
`Source` routine — NOT a curated subset. Intermediate scratch variables
(loop-local helpers, sub-expressions used once and discarded) are excluded.

The user opens the resulting CGNS in Tecplot/ParaView and sees every
quantity as a 3D field. This is the standard debug surface for every
subsequent task.

AD impact:   no  (output-side only, not in residual path)

Context:     `04_ARCHITECTURE.md` §3-4

Files:
- `src/modules/paramTurb.F90`            — add `nTransitionDebug` parameter
- `src/modules/inputParam.F90`           — add `storeTransitionDebug` flag (default `.false.`)
- `src/utils/allocMem.F90` (or wherever block arrays allocate; locate via
  `grep -rn 'allocate(blocks' src/`)    — allocate
                                           `transitionDebug(0:ie,0:je,0:ke,nTransitionDebug)`
- `src/turbulence/saGammaRetheta.F90`    — fill block at end of cell loop in
                                           `Source`, gated by `storeTransitionDebug`
- output dispatch — locate via `grep -rn "case ('eddy')" src/` — add one
                    `case` clause per output name. Mirror in surface dispatch.
- `python/pyADflow.py`                    — expose `storeTransitionDebug`,
                                            add output names to the legal
                                            `volumeVariables` /
                                            `surfaceVariables` lists.

Change:

Step 1 — enumerate what to output. In `saGammaRetheta.F90` `Source` routine,
identify every variable that satisfies BOTH:
- it is a named local or block field with physical meaning (state, source,
  flux, switch, correlation value, Jacobian entry); AND
- it is NOT a one-line scratch sub-expression used solely to build the next
  line.

For the current `Source` routine that includes (use this as the working
list — adjust if the routine has additional named quantities at the time
the task is run):

  state-derived         — gammaForSA, ft2, ss (S̃), rTurb (R_T), vortMag,
                          vortMagLim, mu_t (= w·rlv)
  correlation outputs   — reTheta_target, rethetac, flength, fThetat,
                          fOnset, fTurb_val
  BL proxies            — theta_BL, delta_BL, lambda_theta, Re_V_strain
  sources (assembled)   — P_gamma, E_gamma, P_theta, P_SA_raw,
                          P_SA_eff (= γ·P_SA_raw)
  scaling / time        — timeScale
  Jacobian entries      — qq(i,j,k,1,1), (1,2), (1,3),
                          (2,1), (2,2), (2,3),
                          (3,1), (3,2), (3,3)

Set `nTransitionDebug` to the count of unique names from this enumeration
(currently around 26; recount when you make the list — get it right, don't
hardcode 26 if your scan finds 24 or 28).

Step 2 — name them. Use lowercase ASCII names that match the Fortran
identifiers, with no underscores: `gammaforsa`, `ft2`, `stilde`, `rturb`,
`vortmag`, `vortmaglim`, `mut`, `rethetatarget`, `rethetac`, `flength`,
`fthetat`, `fonset`, `fturb`, `thetabl`, `deltabl`, `lambdatheta`,
`revstrain`, `pgamma`, `egamma`, `ptheta`, `psaraw`, `psaeff`,
`timescale`, `qq11`, `qq12`, `qq13`, `qq21`, `qq22`, `qq23`,
`qq31`, `qq32`, `qq33`. Maintain a comment block at the top of
`saGammaRetheta.F90` listing slot → name. Keep slot order stable (do NOT
renumber later — readers' scripts depend on it).

Step 3 — module parameter. In `paramTurb.F90`:

```fortran
! transitionDebug slot map: see comment block at top of saGammaRetheta.F90.
! Do NOT renumber existing slots — append new ones at the end.
integer(kind=intType), parameter :: nTransitionDebug = <count>
```

Step 4 — allocation. In whatever routine sizes the block arrays:

```fortran
allocate(blocks(nn)%transitionDebug(0:ie, 0:je, 0:ke, nTransitionDebug), &
         stat=ierr)
if (ierr /= 0) call terminate('allocMemBlock', &
                              'Allocation failure for transitionDebug')
blocks(nn)%transitionDebug = zero
```

Step 5 — fill block. Inside the existing cell loop in `Source`, after all
named quantities and `qq` are computed:

```fortran
fillDebug: if (storeTransitionDebug) then
    transitionDebug(i,j,k, 1) = gammaForSA
    transitionDebug(i,j,k, 2) = ft2
    transitionDebug(i,j,k, 3) = ss              ! S̃
    transitionDebug(i,j,k, 4) = rTurb
    ! ... one line per slot, in the order of the comment block ...
    transitionDebug(i,j,k, n) = qq(i,j,k,3,3)
end if fillDebug
```

For source terms, if they were never assembled into a single named variable
(e.g. P_γ was applied as additive contributions to `dw`), assemble it into
a local right before the fill block:

```fortran
P_gamma_local = ca1 * flength * fOnset * vortMagLim &
                * (one - rsaGRce1*w(i,j,k,itu2))
! similarly for E_gamma_local, P_theta_local, P_SA_raw_local, P_SA_eff_local
```

Use these local copies in the `transitionDebug(...)` writes. Do NOT change
the existing additive accumulation into `dw`.

Step 6 — output dispatch. For each name in the slot map, add a `case`
clause matching the existing `eddy` pattern in whichever file the dispatch
lives in:

```fortran
case ('fonset')
    do k = 2, kl ; do j = 2, jl ; do i = 2, il
        buffer(i,j,k) = transitionDebug(i,j,k,<slot>)
    end do ; end do ; end do
```

Mirror in the surface dispatch. Also register the standalone state names
`gamma`, `retheta`, `nutilde` if they aren't already exposed through the
generic itu1/itu2/itu3 mechanism.

Step 7 — Python wrapper. In `pyADflow.py`:

```python
'storetransitiondebug': [bool, False],
```

and append every name from the slot map to the legal-output lists.

Step 8 — option-validation guard. If a user puts a transition-debug name in
`volumeVariables` while `storeTransitionDebug=False`, raise:

```
requested '<name>' in volumeVariables but storeTransitionDebug=False —
output would be uninitialized. Set storeTransitionDebug=True or remove
the transition names.
```

Compile:     `cd build && make -j`

Done when:
- `make` exits 0.
- The slot-map comment block exists at the top of `saGammaRetheta.F90` and
  the count matches `nTransitionDebug` in `paramTurb.F90`.
- `grep` for each name finds it in BOTH the volume and surface dispatch.
- `storeTransitionDebug` is a recognized Python option.
- The fill block compiles without warnings.

Status row:  flip A2 to ✅

End-of-turn: STATUS: READY TO RUN

---

## A3 — Smoke run: transition on, γ forced to 1

Goal: Confirm `useTransitionModel=True` runs without NaN when γ is pinned to
1 everywhere. Pinned γ=1 means the model is "neutral" — transition equations
solved but their effect on SA is identity. Result should match SA-only.

AD impact:   yes (touches residual)
Context:     `04_ARCHITECTURE.md` §4
Files:
- `src/modules/inputParam.F90`         — add `forcedGamma` real, default `-1.0`
                                          (sentinel = "not forced")
- `src/turbulence/saGammaRetheta.F90`  — at the point `gammaForSA` is read
                                          (around line 470 in the Source loop),
                                          if `forcedGamma >= zero`, override:
                                          `gammaForSA = forcedGamma`. Same for
                                          the value plugged into source terms
                                          where γ appears. Wrap in a single
                                          `if` to keep the diff small.
- `python/pyADflow.py`                  — expose `forcedGamma`

Change:      Single conditional override. Do NOT modify the transport
             equations themselves — γ still solves its PDE, but the value
             coupled to SA and to the source terms is overridden when
             `forcedGamma >= 0`.

Compile:     standard
Done when:
- `make` exits 0.
- With `forcedGamma=1.0`, the binary is ready for a smoke run.
Status row:  flip A3 to ✅
End-of-turn: STATUS: TAPENADE NEEDED

---

## B1 — Verify timeScale matches ADflow's nondimensionalization convention

Goal: The paper writes `t = 500·μ/(ρ·U²)·(1/Re)` in dimensional form. ADflow
runs fully nondimensional: `rlv = μ/μ_∞`, velocity by U_∞, density by ρ_∞,
length by L_ref. The `1/Re` factor is implicit in many ADflow expressions
because the working variables already absorb it. The current code

```fortran
timeScale = 500.0_realType * rlv(i,j,k) / max(rho * uMag2, xminn)
```

may already be correct under ADflow's convention, or it may need an
explicit `* Reynolds` or `/ Reynolds`. This task DOES NOT apply a blind
fix — it determines the convention from existing code and either documents
the current line as correct or applies the minimal change required to match
the established pattern.

AD impact:   yes  (if a code change is made; no if only a comment is added)

Context:     `01_PAPER_REFERENCE.md` §2.2

Files:
- `src/turbulence/saGammaRetheta.F90` line 430-431 (read first; edit only
  if Step 4 below requires it)

Change:

Step 1 — find the convention. Locate where `rlv` is computed (likely
`computeUtau`, `computeLamViscosity`, or similar in `src/solver/` or
`src/preprocessing/`). Read the comment header / formula and confirm
whether `rlv = μ/μ_∞` (pure ratio, no Re) or `rlv = μ/(μ_∞·Re)` (ratio
with Re absorbed).

Step 2 — find an analogous expression in existing code. The closest
analog is the SA model's own viscous diffusion or destruction terms,
which in dimensional form contain `μ/(ρ·something²)` and in nondim form
should reveal whether ADflow writes `/Reynolds` explicitly or absorbs it.
Look in `saSource` / `saModel.F90` (or wherever SA's viscous terms
assemble) for any line of the form `μ / (ρ·U·L)` and note whether
`Reynolds` (or `RGas`, `Re_inf`, etc.) appears.

Step 3 — compare. The transition `timeScale` assembles a quantity with
units of `length·time/length² = time/length` (or pure time if you
nondimensionalize by L_ref/U_∞). Determine whether the current line
matches the convention of Step 1+2.

Step 4 — three possible outcomes:

  Outcome A — current code is correct under ADflow's convention.
    Action: add a one-line comment above the timeScale assignment:
    ```fortran
    ! Nondim form: rlv already = μ/μ_∞; Re absorbed in nondim NS,
    ! consistent with <reference expression in saSource.F90:NNN>.
    timeScale = 500.0_realType * rlv(i,j,k) / max(rho * uMag2, xminn)
    ```
    No code change. AD impact: no. End STATUS: READY TO RUN.

  Outcome B — current code is missing an explicit factor.
    Action: apply the minimal fix (e.g. `* Reynolds` or `/ Reynolds`)
    matching the pattern from Step 2. Add the same documenting comment.
    AD impact: yes. End STATUS: TAPENADE NEEDED.

  Outcome C — convention is unclear after Steps 1-3.
    Action: do NOT guess. End STATUS: BLOCKED — convention unclear,
    cannot find analogous nondim expression. The user resolves manually.

Compile:     `cd build && make -j`

Done when:
- `make` exits 0.
- A code comment exists above the timeScale line citing the analogous
  expression in another ADflow source file by `file:line`.
- Either no code change (Outcome A) or a minimal documented change
  (Outcome B).

Status row:  flip B1 to ✅ (Outcome A or B) or leave ❌ with note (Outcome C)

End-of-turn: STATUS: READY TO RUN  (Outcome A)
              STATUS: TAPENADE NEEDED (Outcome B)
              STATUS: BLOCKED — <reason>  (Outcome C)

---

## B2 — Verify φ_p overflow safety (Algorithm 1)

Goal: The smooth max/min `phi_p` is used at p=300; naïve `log(exp(p·g1) +
exp(p·g2))/p` overflows. Confirm the implementation in
`saGammaRethetaHelpers.F90` matches the four-branch overflow-safe form in
paper Algorithm 1.

AD impact:   yes (helper is in the residual path)
Context:     `01_PAPER_REFERENCE.md` §6 (Algorithm 1)
Files:       `src/turbulence/saGammaRethetaHelpers.F90` (function `phiP` or `smooth_max_min` — locate)

Change:
- Read the function body. Confirm it has the four branches:
  1. `|p·(g2-g1)| < 1e-15`         → average + correction
  2. `p·(g2-g1) >  20`              → return `g2 - log(2)/p`
  3. `p·(g2-g1) < -20`              → return `g1 - log(2)/p`
  4. otherwise                      → general formula
- If any branch is missing or the threshold differs significantly from 20,
  rewrite to match paper Algorithm 1 exactly (sketch in
  `01_PAPER_REFERENCE.md` §6).
- Do NOT change the function signature.

Compile:     standard
Done when:
- `make` exits 0.
- All four branches present and compile.
Status row:  flip B2 to ✅
End-of-turn: STATUS: TAPENADE NEEDED

---

## B3 — Populate off-diagonal source Jacobian

Goal: The 3×3 source Jacobian `qq(i,j,k,row,col)` currently only has
diagonal entries. Add the six off-diagonal entries listed in paper §7.1 so
DADI's `TurbDADICoupled=True` path is actually coupled.

AD impact:   yes
Context:     `01_PAPER_REFERENCE.md` §7.1, `04_ARCHITECTURE.md` §4 (qq layout)
Files:       `src/turbulence/saGammaRetheta.F90` lines 510-559

Change:
- After the existing diagonal assignments, add:
  - `qq(i,j,k,1,2) = ∂R_ν̃/∂γ`        — derivative of γ·P_SA wrt γ:
    `= P_SA_raw` (i.e. `rsaCb1*(one-ft2)*ss*w(i,j,k,itu1)`).
  - `qq(i,j,k,2,1) = ∂R_γ/∂ν̃`        — comes through `R_T = μ_t/μ` in
    `F_turb` and `E_γ`. Differentiate the destruction:
    `∂E_γ/∂ν̃ = ca2 · Ω_lim · γ · (ce2·γ-1) · ∂F_turb/∂R_T · ∂R_T/∂ν̃`,
    with `∂R_T/∂ν̃ = 1/μ` (i.e. `1.0/rlv(i,j,k)` in nondimensional form).
    For sLM2015 `F_turb = (1-F_onset)·exp(-R_T)`:
    `∂F_turb/∂R_T = -(1-F_onset)·exp(-R_T)`. For LM2015 `F_turb =
    exp(-(R_T/4)^4)`: `∂F_turb/∂R_T = -(R_T^3/64)·exp(-(R_T/4)^4)`.
    Use whichever is active per `useLM2015Fturb`.
  - `qq(i,j,k,2,3) = ∂R_γ/∂Re̅θt`    — through `Re_θc(Re̅θt)` in `F_onset`.
    Compute `∂Re_θc/∂Re̅θt` from the active correlation (sLM2015:
    `0.67 + (24/240)·cos(Re̅θt/240+0.5)`; LM2015: differentiate the
    piecewise polynomial). Then chain through `F_onset` and `P_γ`.
  - `qq(i,j,k,3,2) = ∂R_θt/∂γ`       — through `F_θt(γ)` in `P_θt`.
    Differentiate Eq. 3-4: `∂F_θt/∂γ = -2·(γ-1/ce2)/(1-1/ce2)^2`
    inside the active branch of the outer min/max; if the smooth `phi_p`
    is used, differentiate that. Use a small-difference fallback if the
    analytic derivative gets messy: finite-difference `F_θt` at runtime
    once and store the slope.
  - `qq(i,j,k,1,3)` and `qq(i,j,k,3,1)` — set to zero (paper §7.1: "usually
    ~0"). Add the explicit `= zero` assignments so it is documented.
- Keep all expressions inside the existing cell loop, alongside the
  diagonals. Do not allocate new arrays.

Compile:     standard
Done when:
- `make` exits 0.
- All six off-diagonals are explicitly assigned. (Grep `qq(i,j,k,1,2)`,
  `qq(i,j,k,2,1)`, etc.)
Status row:  flip B3 to ✅
End-of-turn: STATUS: TAPENADE NEEDED

---

## B4 — Initialize row/column scaling factors

Goal: `turbResScale` is declared but never initialized → falls to 1.0 →
scaling is a no-op. Set the three components to (ν̃_max, γ_max, Re̅θt_max) =
(1e3, 10, 1e4) per paper §7.3.

AD impact:   no
Context:     `01_PAPER_REFERENCE.md` §7.3, `04_ARCHITECTURE.md` §3
Files:
- `src/modules/inputParam.F90` lines 293,298
- `src/initFlow/initializeFlow.F90` (or similar early-init point —
  verify by searching for where other scale factors are set)

Change:
- In `inputParam.F90`, set the default values right where `turbResScale` is
  declared:
  ```fortran
  real(kind=realType), dimension(3) :: turbResScale = (/1.0e3_realType, &
                                                        1.0e1_realType, &
                                                        1.0e4_realType/)
  ```
  Match the existing array shape (it may be `dimension(:)` allocatable —
  in that case allocate and fill in `initializeFlow`).
- Expose to Python in `pyADflow.py` as `turbResScale` so user can override.

Compile:     standard
Done when:
- `make` exits 0.
- `turbResScale` has a non-trivial default at first use in
  `saGammaReThetaSolve` (lines 1329-1331).
Status row:  flip B4 to ✅
End-of-turn: STATUS: READY TO RUN

---

## C1 — First-order upwind option for γ, Re̅θt convection

Goal: Paper §IV.A recommends first-order upwind for the transition
transport equations. Add a runtime switch and wire it into the convection
discretization for γ and Re̅θt.

AD impact:   yes (changes residual)
Context:     `01_PAPER_REFERENCE.md` §11, `04_ARCHITECTURE.md` §4
Files:
- `src/modules/inputParam.F90`         — add `transitionFirstOrderUpwind`
                                          logical, default `.true.`
- `src/turbulence/saGammaRetheta.F90`  — locate the convection-flux block
                                          in `Source` (the part that uses
                                          face-value reconstruction for itu2,
                                          itu3). Add a branch:
  ```fortran
  if (transitionFirstOrderUpwind) then
      ! upwind: face value = upstream cell value
      gammaFace  = w(iUp,jUp,kUp, itu2)
      rethetaFace= w(iUp,jUp,kUp, itu3)
  else
      ! existing higher-order reconstruction
      ...
  end if
  ```
  where `(iUp,jUp,kUp)` is determined by the sign of the contravariant
  velocity component on the face.
- `python/pyADflow.py`                  — expose `transitionFirstOrderUpwind`
                                          (default True per paper)

Change:
- Do NOT modify SA convection (itu1) — only γ and Re̅θt.
- Default ON because paper recommends it.

Compile:     standard
Done when:
- `make` exits 0.
- `grep transitionFirstOrderUpwind` finds the option used in the Source
  routine inside an `if`.
Status row:  flip C1 to ✅
End-of-turn: STATUS: TAPENADE NEEDED

---

## C2 — Source-term dt restriction: DADI 3×3 eigenvalue

Goal: For the DADI coupled path, the local time step must be capped by the
largest positive eigenvalue of the 3×3 source Jacobian (paper Eq. 59).
Diagonal-only fallback exists (lines 538, 551). Add the eigenvalue-based
cap that uses the off-diagonals from B3.

AD impact:   no  (it's a time-step cap, outside residual)
Context:     `01_PAPER_REFERENCE.md` §7.1, §7.4
Files:       `src/turbulence/saGammaRetheta.F90` (the DADI solver section,
             ~line 1330 onward)

Change:
- Add a function `max_pos_eigenvalue_3x3(J) result(lambda_max)` near the
  bottom of the file. Use closed-form roots of the 3×3 characteristic
  polynomial OR call LAPACK `dgeev`. Closed-form is preferred to avoid the
  LAPACK call overhead per cell:
  ```fortran
  ! Cubic roots via Cardano. Coefficients of det(J - λI) = 0:
  ! λ^3 + a2·λ^2 + a1·λ + a0 = 0
  ! where a2 = -trace(J)
  !       a1 =  sum of 2x2 principal minors
  !       a0 = -det(J)
  ! Take real parts; lambda_max = max(re(λ_i)).
  ```
  Implementation sketch is enough — Claude Code writes the actual cubic
  solver. Use `realType` throughout. If complex roots arise, take real part
  and largest of the three.
- In the DADI block where `qq(i,j,k,m,m) = max(qq(i,j,k,m,m), -qq(i,j,k,m,m)/0.9)`
  is currently applied, replace with:
  ```fortran
  if (TurbDADICoupled) then
      lambda_max = max_pos_eigenvalue_3x3(qq(i,j,k,1:3,1:3))
      dt_cap = 0.9_realType / max(lambda_max, eps_small)
      ! enforce: shrink dt to dt_cap if smaller than current
      ! ... existing dt update ...
  else
      ! existing diagonal-only path, unchanged
  end if
  ```
- Guard with the `transitionSrcDtRestrict` flag so user can disable if
  desired (default `.true.`).

Compile:     standard
Done when:
- `make` exits 0.
- New function exists and is called only when `TurbDADICoupled=True`.
- Diagonal path is untouched.
Status row:  flip C2 to ✅
End-of-turn: STATUS: READY TO RUN

---

## C3 — Source-term dt restriction: Turb-ANK CFL cap

Goal: The Turb-ANK path needs the same restriction: cap CFL by the source
Jacobian. Paper §7.1.

AD impact:   no
Context:     `04_ARCHITECTURE.md` §1 (Turb-ANK)
Files:       `src/NKSolver/NKSolvers.F90` (or wherever Turb-ANK CFL is
             updated; locate by `grep -rn 'CFL' src/NKSolver/`)

Change:
- Inside the Turb-ANK CFL update, after the standard CFL is computed,
  shrink it by the source-Jacobian cap. Use the diagonal entries of `qq`
  reduced over the local block to a max:
  ```fortran
  qq_diag_max = maxval( abs(reshape([qq(:,:,:,1,1), qq(:,:,:,2,2), qq(:,:,:,3,3)], [1])) )
  CFL_cap = 0.9_realType / max(qq_diag_max, eps_small)
  CFL = min(CFL, CFL_cap)
  ```
  An MPI ALLREDUCE on `qq_diag_max` may be needed for global consistency —
  follow the existing CFL reduction pattern in the file.
- Guard with `transitionSrcDtRestrict`.

Compile:     standard
Done when:
- `make` exits 0.
- The cap is applied only when `useTransitionModel=True`.
Status row:  flip C3 to ✅
End-of-turn: STATUS: READY TO RUN

---

## C4 — 5-iter deactivation switch

Goal: Paper §IV.B: after 5 successive inexact-Newton iterations without
backtracking and with R_d > 1e-5, deactivate the source-term restriction.
Reactivate if backtracking triggers or residual rises.

AD impact:   no
Context:     `01_PAPER_REFERENCE.md` §7.1 (last paragraph), `04_ARCHITECTURE.md` §1
Files:
- `src/modules/inputParam.F90`         — add `srcDtRestrictActive` logical
                                          (module-level state, init `.true.`)
                                          and `noBacktrackCount` integer
                                          (init 0)
- `src/turbulence/saGammaRetheta.F90`  — read `srcDtRestrictActive` to gate
                                          the C2/C3 caps (i.e. wrap them in
                                          `if (srcDtRestrictActive)` blocks)
- `src/NKSolver/NKSolvers.F90` (or wherever the inexact-Newton loop logs
  backtracking) — increment/reset counter and flip
  `srcDtRestrictActive`

Change:
- After each Newton iteration:
  - If backtracking occurred → reset `noBacktrackCount=0`,
    `srcDtRestrictActive=.true.`.
  - Else → increment `noBacktrackCount`.
  - If `noBacktrackCount >= 5` AND last residual drop > 1e-5 →
    `srcDtRestrictActive=.false.`.
  - If residual rises between iterations → reset and reactivate.

Compile:     standard
Done when:
- `make` exits 0.
- Counter logic compiles and the flag actually gates the dt caps in
  C2 and C3.
Status row:  flip C4 to ✅
End-of-turn: STATUS: READY TO RUN

---

## D1-D4 — Solver path smoke tests

Four solver-mode smoke tests. The four paths together cover: simplest
possible path, off-diagonal coupling inside DADI, full coupling at the
ANK level, and Newton-on-turb-but-decoupled-from-flow. If any path NaNs
that the others survive, it isolates the bug to that path's specific
machinery.

**The four paths:**

| Test | Description                                     | What it exercises                                                 |
|------|-------------------------------------------------|-------------------------------------------------------------------|
| D1   | Decoupled flow/turb, DADI **diagonal**          | basic γ-Re̅θt PDE, BCs, halo, γ·P_SA — minimal coupling           |
| D2   | Decoupled flow/turb, DADI **3×3 coupled**        | adds off-diagonal Jacobian (B3) and source-dt cap (C2)            |
| D3   | **Coupled ANK** (flow + turb in one Newton)     | global Jv coupling — γ, ν̃, flow all in one Krylov                |
| D4   | Decoupled flow/turb, **standard ANK** for turb   | turb-block ANK sub-solver with new equations (C3 CFL cap)         |

D3 is the "coupled ANK like solve" — flow and turbulence in the same
Newton-Krylov system, matrix-free. NK as a separate terminal solver is NOT
tested here; it shares D3's coupling structure with stricter tolerances.

**For all four:**

AD impact:   no  (config-only; no source edits)
Context:     `04_ARCHITECTURE.md` §1
Files:       (none modified — config only)
Compile:     `cd build && make -j` should exit 0 immediately (no dirty state).
Done when:   smoke script runs N iterations without NaN in γ, Re̅θt, ν̃.
End-of-turn: STATUS: READY TO RUN
             (or STATUS: BLOCKED — <which variable went NaN at which iter>)

If any path NaNs, do NOT try to fix it inline. The right response is
`STATUS: BLOCKED` with a precise description; the user creates a follow-up
task using `06_TASK_TEMPLATE.md`.

**Per-path option blocks** (exact flag names verified against
`python/pyADflow.py defaultOptions` in the first session of D-block —
update here once confirmed):

### D1 — Decoupled / DADI diagonal

```python
opts = {
    'turbulenceModel':     'SA',
    'useTransitionModel':  True,
    'equationType':        'RANS',
    'useANKSolver':        True,
    'ankCoupledSwitchTol':  0.0,    # never switch to coupled ANK
    'turbDADICoupled':      False,  # diagonal — ignore off-diagonals
    'storeTransitionDebug': True,
    'volumeVariables':      [...standard...] + [...all transition names...],
}
```

### D2 — Decoupled / DADI 3×3 coupled

Same as D1 except:
```python
    'turbDADICoupled':      True,   # use full 3×3 block (off-diags from B3)
```

### D3 — Coupled ANK (flow + turb in one Newton)

```python
opts = {
    'turbulenceModel':     'SA',
    'useTransitionModel':  True,
    'equationType':        'RANS',
    'useANKSolver':        True,
    'ankCoupledSwitchTol':  1.0,    # always coupled — flow+turb together
    'storeTransitionDebug': True,
    'volumeVariables':      [...],
}
```

### D4 — Decoupled / standard ANK for turb

Decoupled at the flow level; turbulence solved by its own ANK Newton-Krylov
sub-solver (not DADI).

```python
opts = {
    'turbulenceModel':       'SA',
    'useTransitionModel':    True,
    'equationType':          'RANS',
    'useANKSolver':          True,
    'ankCoupledSwitchTol':    0.0,    # decoupled at flow level
    'turbSubsolver':         'ANK',   # ← verify exact option name in pyADflow.py
    'storeTransitionDebug':  True,
    'volumeVariables':        [...],
}
```

**First D-task action** (D1 specifically, since options carry over):
`grep -n "ankCoupledSwitchTol\|turbDADICoupled\|turbSubsolver" python/pyADflow.py`
to confirm the exact option names used by this branch. Update the four
blocks above in this file with the verified names before running. If a
flag does not exist (e.g. no `turbSubsolver` option exists yet), output
`STATUS: BLOCKED — option 'turbSubsolver' not exposed in pyADflow.py;
needs prior task to expose it.`

Status rows:  flip D1, D2, D3, D4 to ✅ each as it lands

---

## E1 — Tapenade clean regen — full sweep

Goal: After all source-residual changes (B1, B2, B3, C1, plus A3 for
forcedGamma), the AD-generated files in `src/adjoint/output{Forward,
Reverse,ReverseFast}/` are stale. Confirm the regen commands produce clean
output.

AD impact:   yes (this IS the AD step)
Context:     none — user provides the Tapenade command sequence
Files:       `src/adjoint/output*` (regenerated, not hand-edited)

Change:
- This task is "user runs Tapenade commands, Claude Code stages the
  resulting changes". Do NOT hand-edit AD files.
- After regen, run `cd build && make -j` and resolve any breakage in
  manually-maintained glue code.

Compile:     standard
Done when:
- `make` exits 0 with the regenerated AD files.
- No diff in non-AD source vs HEAD before this task.
Status row:  flip E1 to ✅
End-of-turn: STATUS: READY TO RUN

---

## E2 — Build with AD enabled (adjoint optimization capability gate)

Goal: Confirm ADflow still compiles when built with automatic differentiation
enabled. This is the build configuration that links the Tapenade-generated
forward and reverse derivative routines (regenerated in E1) and produces a
binary capable of computing dCd/dα, dCl/dx_geom, etc. via adjoint — the
entire reason ADflow exists in MDOLab.

The standard build (`make`) tested in A1-D4 produces only the forward RANS
solver. The AD-enabled build is a different target / config flag and has
its own failure modes (Tapenade-generated Fortran that doesn't lint,
missing module paths, type mismatches in hand-written glue between
generated and primal code).

This task is purely a build gate. It does NOT exercise the adjoint at
runtime. The adjoint linearization for γ-Re̅θt is frozen at this stage
(zero seeds), so any optimization run would see through the transition
equations — but the AD routines must still EXIST and COMPILE for the
freeze plumbing to work. E2 confirms that.

If E2 fails, gradient-based optimization is broken on this branch.

AD impact:   no  (build-time only; no source edits beyond glue fixes)

Context:     project README and `Makefile` for AD build flags

Files:       (none modified unless glue code is broken — in which case
              the broken file gets a minimal fix, scoped to the specific
              compile error)

Change:

Step 1 — locate the AD build invocation. Read `Makefile`, any `config/*.mk`
includes, and the project README. Common patterns:
- A separate target: `make -j ad` or `make -j adjoint`
- A flag: `make -j AD=1` or `make -j MODE=adjoint`
- A separate config file: `cp config.LINUX_INTEL_OPENMPI_AD.mk config.mk &&
  make`

Step 2 — invoke the AD build:

```bash
cd build
make clean         # critical: don't reuse standard-build object files
<the AD build command from Step 1>
```

Step 3 — three outcomes:

  Outcome A — AD build exits 0.
    Status row: flip E2 to ✅. End STATUS: READY TO RUN.

  Outcome B — AD build fails with errors in `src/adjoint/output*` files.
    These are Tapenade-generated; do NOT hand-edit them. End STATUS:
    BLOCKED — Tapenade output won't compile, E1 needs to be redone with
    different directives. Quote the first error line in the status.

  Outcome C — AD build fails with errors in hand-written glue code (any
              file outside `src/adjoint/output*/`).
    Apply the minimal fix to the glue file. Common fixes:
    - missing `use` of a new module added by another task
    - new variable in primal code not declared in the corresponding
      hand-edited reverse-mode adjoint file
    - type mismatch where a new variable was added with `realType` but
      the AD-side glue expected `real(kind=8)`
    Do NOT add new functionality, only fix what the compiler complains
    about. End STATUS: READY TO RUN once it compiles.

Compile:     AD-enabled build (per Step 1)

Done when:
- The AD build exits 0.
- No diff in `src/adjoint/output*` files vs after E1 (those are
  regenerated, not hand-edited).
- Any glue-code fix is minimal and scoped to the specific compile error.

Status row:  flip E2 to ✅

End-of-turn: STATUS: READY TO RUN  (Outcomes A, C)
              STATUS: BLOCKED — <reason>  (Outcome B)

---

## After E2

Implementation is feature-complete. User runs validation cases (NLF0416,
S809) on the smoke script, inspects the volume CGNS, compares to paper
figures. Bugs surfaced at that stage become new tasks created via
`06_TASK_TEMPLATE.md`.
