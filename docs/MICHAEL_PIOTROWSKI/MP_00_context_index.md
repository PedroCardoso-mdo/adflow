# Michael Piotrowski — complete publication set — sub-index

Every Piotrowski paper plus his PhD thesis, transcribed and digested. This is
the **author's whole body of work on SA-sLM2015**, of which the 2020 paper
(our physics source of truth) is one piece. Several of these files **correct or
supersede** parts of the 2020 paper.

**Load only what a task needs** (CLAUDE.md rule 9). Full routing lives in the
master index: [`../README.md`](../README.md).

## Structure — two files per source PDF

| Suffix | What it is |
|---|---|
| `_full.md` | **Transcript.** Verbatim `pdftotext -layout` extraction, page-delimited, nothing removed but bare page-number lines. Complete and unedited — but **equations are raw extractions** and may be flattened or split across lines. Confirm any equation against the PDF before implementing. |
| `_digest.md` | **What matters to us.** Hand-written: the findings, equations, algorithms and constants that bear on this branch, with the gaps they close in our own open items. |

Source PDFs: `/home/mdo/Desktop/ARTIGOS SMOTH/`.

## Files

| # | Source | Digest covers | Priority for us |
|---|---|---|---|
| **MP_01** | P&Z 2019, *Investigation of a LCTM in a Newton-Krylov Algorithm* (AIAA SciTech) | The origin paper. Measured failure modes: **non-smooth residual kinks stall the solver at 1e-5 / 3 orders even with source-Δt stepping on**. Exponential-penalty smooth min/max with denormal guard; Gaussian `F_length`; source-Δt rule; matrix-vs-scalar dissipation matters for crossflow. | **High** — our clip sites |
| **MP_02** | P&Z 2020, *Smooth LCTM for the SA Turbulence Model* (AIAA J 58(10)) | The canonical paper. **Use [`../SA_GAMMA_RETHETHA_BASE/SAGR_01_…`](../SA_GAMMA_RETHETHA_BASE/SAGR_01_paper_piotrowski_zingg_2020.md) for equations** (hand-transcribed LaTeX); the digest here is a pointer + a table of what the later papers supersede. | Physics reference |
| **MP_03** | P&Z 2022 Special Session, *Numerical Behaviour … in a Newton-Krylov Flow Solver* | **4 coupling strategies × source-Δt on/off × 3 geometries.** Fully-coupled + source-Δt restriction wins, and the ranking **inverts with case difficulty**. Partial coupling stalls at 4–5 orders on hard 3D cases. | **Highest — our stall** |
| **MP_04** | P&Z 2022, *… in a Discrete-Adjoint ASO Algorithm* | Adjoint verification practice (CS for Jacobian + metric linearization, FD for gradient); residual depends on control points **through off-wall spacing**; streamwise grid requirements tighten with Re; design space is multi-modal. | Medium — adjoint/opt |
| **MP_05** | P&Z 2023, *Compressibility corrections…* (Aeronautical Journal) | `ψ` / `ψ_scf` corrections; **`D_scf` replaced by `F_onset,scf` in the γ equation**; **documented crossflow initialization pathology + the staged-activation fix (crossflow off until rel 1e-5)**; Appendix A = cleanest non-dimensional statement of the whole model. | **High — our crossflow** |
| **MP_06** | Piotrowski, **PhD thesis** (156 pp., U. Toronto) | Superset of all of the above. **§3.1 solver chapter: the explicit `S_r`/`S_c` scaling matrices (Eq. 3.7) and Algorithms 2/3/4 pseudocode.** App. C: transition-length modification. §6.2: his own unsolved problems. | **Highest — solver detail** |
| **MP_07** | Chau, Piotrowski & Duensing 2026 (J. Aircraft) | Cruise-slotted TTBW optimization. **Fully turbulent — no transition model.** Filed for completeness. | Low |

## Quick routing

| Question | Open |
|---|---|
| NK stalls / linear solve saturates / which coupling to use | `MP_03_digest` → `MP_06_digest` §1–3 |
| Row/column scaling of the linear system (`S_r`, `S_c`, Eq. 58) | **`MP_06_digest` §1** (explicit matrix + constants) |
| Physicality check, update damping, line-search design | **`MP_06_digest` §2** (Algorithms 2/3/4) |
| Source-term Δt restriction (Eq. 59) — value, eigenvalue, on/off rule | `MP_06_digest` §3 → `MP_03_digest` |
| Crossflow stalls / wrong transition front | **`MP_05_digest` §2–3** |
| Clipping, kinks, smooth min/max | `MP_01_digest` §1, §3 |
| Transonic / compressibility | `MP_05_digest` §4 |
| Any transition equation or constant | `../SA_GAMMA_RETHETHA_BASE/SAGR_01_…` (source of truth); cross-check non-dimensionalization against `MP_05_full` Appendix A |
| Adjoint verification coverage | `MP_04_digest` §1 |

## Open items in our own docs that these files speak to

| Our open item | Where |
|---|---|
| `S_r` geometric row scaling stalls the NK linear solve (`SAGR_02` §5, §8.3) | **`MP_06_digest` §1** — exponents are `J^{2/3}` (mean flow) vs `J^{−1/3}` (turb/transition), i.e. **opposite signs**; ours are same-signed |
| NK has no physicality/bounds check (`SAGR_02` §7, §8.2) | **`MP_06_digest` §2** — Algorithm 2 pseudocode, incl. reject-and-halve-Δt_ref |
| Deep-NK wall: "no option-level lever remains" (`CORE_02` §Options analysis) | `MP_01_digest` §1 (kinks) + `MP_06_digest` §1 (scaling) — both are code-level, not option-level |
| `transitionCrossflow` plateaus at ~3e-2 (CLAUDE.md rule 3) | **`MP_05_digest` §3** — documented init-dependence + staged activation |
| Hard γ/Re̅θt clips in the residual (`saGammaRetheta.F90:541, 602, 2556, 2592`) | `MP_01_digest` §1 and `MP_06_digest` §2 — the authors damp, never clip, and warn clipping stalls the solver |

## Provenance

Transcribed 2026-08-13 from the seven PDFs in `/home/mdo/Desktop/ARTIGOS SMOTH/`
(17 / 38 / 17 / 20 / 30 / 156 / 35 pages). `S_r` Eq. 3.7 was additionally
verified by reading thesis PDF page 61 directly, not from the text extraction.
