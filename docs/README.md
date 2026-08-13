# Documentation & Knowledge Base — SA-γ-Re̅θt Transition Model

Master index for **all** project documentation. This branch carries two kinds
of docs:

- **Project-specific** (`docs/CORE_BASE/`, `docs/VERIFICATION/`, `docs/ARCHIVE/`,
  `docs/task-log/`) — unique to the Piotrowski & Zingg (2020) γ-Re̅θt
  implementation on this branch.
- **Reference KB** (`docs/ADFLOW_BASE/`, `docs/SA_GAMMA_RETHETHA_BASE/`) — distilled
  extractions of the published papers and the official ADflow docs.

> **Read only the file(s) a task needs — never load the whole KB.** Every file
> is self-contained on purpose (token discipline; see `../CLAUDE.md` rule 9).
> Use the [Task → files routing](#task--files-routing) below to pick the minimum
> set, then open only those.

Upstream ADflow docs (`/README.md`, `doc/`, `tests/README.md`,
`input_files/Readme.md`, `LICENSE.md`) are **not** part of this KB and are left
as-is.

---

## Task → files routing

Read the listed files **in the order shown**; stop as soon as you have the
answer. `→` means "if the first isn't enough, go to the next."

| Your task / question | Read (in order) |
|---|---|
| Implement/extend a transition equation, constant, or algorithm | `SA_GAMMA_RETHETHA_BASE/SAGR_01_paper_piotrowski_zingg_2020.md` → `CORE_BASE/CORE_01_architecture.md` |
| Any equation touching velocity, viscosity, Reynolds number, or time scales | `ADFLOW_BASE/ADFLOW_08_nondimensionalization.md` (**first**) → `SA_GAMMA_RETHETHA_BASE/SAGR_01_paper_piotrowski_zingg_2020.md` |
| Paper symbol ↔ code flag (transition model) | `SA_GAMMA_RETHETHA_BASE/SAGR_01_paper_piotrowski_zingg_2020.md` + `CORE_BASE/CORE_01_architecture.md` (Part 2) |
| Exact **transition** runtime option name / default / enum | `CORE_BASE/CORE_01_architecture.md` (Part 2) |
| Solver architecture, state-vector layout, code/module locations | `CORE_BASE/CORE_01_architecture.md` (Part 1) |
| Flow / ANK / NK theory, equations, defaults | `ADFLOW_BASE/ADFLOW_01_flow_solver_theory.md` |
| What an ANK/NK option does mechanically | `ADFLOW_BASE/ADFLOW_03_concordance.md` → `ADFLOW_BASE/ADFLOW_01_flow_solver_theory.md` |
| Choose R0/R1/R2, coupled vs decoupled, NK switch | `ADFLOW_BASE/ADFLOW_01_flow_solver_theory.md` → `ADFLOW_BASE/ADFLOW_04_debugging_playbook.md` |
| Exact **ADflow** (flow/ANK/NK/adjoint) option name / default / enum | `ADFLOW_BASE/ADFLOW_05_options_and_operations_devguide.md` |
| Paper symbol ↔ code flag (ADflow solver) | `ADFLOW_BASE/ADFLOW_03_concordance.md` |
| Debug a stalling / diverging run | `ADFLOW_BASE/ADFLOW_04_debugging_playbook.md` → `ADFLOW_BASE/ADFLOW_05_options_and_operations_devguide.md` → `ADFLOW_BASE/ADFLOW_01_flow_solver_theory.md` |
| Rotating mesh / rotating frame of reference (rotor, propeller, turbine) | `VERIFICATION/VERIF_02_rotating_frame_audit.md` → `ADFLOW_BASE/ADFLOW_08_nondimensionalization.md` (§5 vortLim, §6 crossflow) |
| Why does the paper converge faster than ADflow / solver-algorithm gaps | `SA_GAMMA_RETHETHA_BASE/SAGR_02_adflow_vs_paper_solver.md` |
| S809/NLF0416 physics validation vs the paper — what was tested, eliminated, still open | `VERIFICATION/VERIF_03_paper_validation_campaign.md` |
| Adjoint / AD theory (general) | `ADFLOW_BASE/ADFLOW_02_adjoint_autodiff_theory.md` |
| Implement or extend the adjoint / AD on **this branch** | `ADFLOW_BASE/ADFLOW_09_adjoint_trace.md` → `ADFLOW_BASE/ADFLOW_02_adjoint_autodiff_theory.md` → `ADFLOW_BASE/ADFLOW_03_concordance.md` |
| Gradients wrong / adjoint won't converge | `ADFLOW_BASE/ADFLOW_04_debugging_playbook.md` → `VERIFICATION/VERIF_01_debugging_derivatives.md` (FD→dot-product→CS checking ladder) → `ADFLOW_BASE/ADFLOW_02_adjoint_autodiff_theory.md` → `ADFLOW_BASE/ADFLOW_09_adjoint_trace.md` |
| Validate partials / plan AD test campaign | `CORE_BASE/CORE_00_current_task.md` → `VERIFICATION/VERIF_01_debugging_derivatives.md` (generic checking ladder) → `VERIFICATION/VERIF_00_three_stage_verification.md` (tests already run, commands, results) → `ARCHIVE/ARCHIVE_03_adjoint_audit_2026-07-07.md` → `ADFLOW_BASE/ADFLOW_07_sst_dev_lessons.md` (verification ladder + watch items) |
| Rerun/reproduce the low-level adjoint verification tests (dot-product, fast_b, 3-way fwd) | `VERIFICATION/VERIF_00_three_stage_verification.md` (§"Canonical way to run" → `tests/reg_tests/run_sagr_tests.sh`) |
| Run / retrain / regen-`w` for the SA-GR partials tests (how-to) | `VERIFICATION/VERIF_00_three_stage_verification.md` §"Canonical way to run" → `tests/reg_tests/dev/README.md` (non-standard `w` generation) |
| How another multi-equation turb model was differentiated (SST precedent) | `ADFLOW_BASE/ADFLOW_07_sst_dev_lessons.md` |
| Why is X implemented this way / was this already discussed | `CORE_BASE/CORE_03_design_decisions.md` (memory of past discussion, not a spec — code/paper win if it disagrees) |
| Converge an SA-GR case (phase ladder, switch tols, LS options, measured limits, falsified levers) | `CORE_BASE/CORE_02_convergence_strategy.md` → `ADFLOW_BASE/ADFLOW_04_debugging_playbook.md` |
| Run a case on the HPC (Deucalion) — envs, job shape, checkpoints, graceful kill | run-tree `HPC_HOWTO.md` (`/home/mdo/Desktop/Run/MDO_PhD/Transition/gama_rethetha/`) + auto-memory HPC rules; **standing rule: all solver runs on HPC, local = analysis only** |
| What's the task in progress right now | `CORE_BASE/CORE_00_current_task.md` |
| What was already finished, and how | `task-log/README.md` (index → per-task case files) |
| Provenance / sources / where a fact came from | this file + `ADFLOW_BASE/ADFLOW_00_context_index.md` |

---

## File index

### Core (`docs/CORE_BASE/`)

| File | What's in it |
|------|--------------|
| [README.md](README.md) | This master index + routing table. |
| [CORE_BASE/CORE_01_architecture.md](CORE_BASE/CORE_01_architecture.md) | Solver architecture, state-vector layout, key code/module locations, user constraints, and the complete reference for every runtime option added for the transition model. |
| [ADFLOW_BASE/ADFLOW_08_nondimensionalization.md](ADFLOW_BASE/ADFLOW_08_nondimensionalization.md) | How ADflow makes the governing equations dimensionless — the **pressure–density (p-ρ) scaling** (velocity normalizes to M·√γ, *not* 1). Read before touching any equation with velocity, viscosity, time scales, or Reynolds number. |
| [CORE_BASE/CORE_04_todo.md](CORE_BASE/CORE_04_todo.md) | Deferred tuning/improvement items (decided, not defects): `turbResScale` calibration, Eq. 59 relaxation options, damping-clip fallback, open dúvida D-A2-3, blockettes, test-infra restructure. Each item links back to its `CORE_BASE/CORE_03_design_decisions.md` analysis. |
| [CORE_BASE/CORE_00_current_task.md](CORE_BASE/CORE_00_current_task.md) | The single task currently in progress (CLAUDE.md: one task per session) — objective, scoped context, working files, checklist. Overwritten each session. |
| [task-log/](task-log/README.md) | Finished-task log, one file per task (a "case"), added over time — index + template in `task-log/README.md`. |
| [CORE_BASE/CORE_02_convergence_strategy.md](CORE_BASE/CORE_02_convergence_strategy.md) | **The validated SA-GR convergence recipe** (ANK→CANK→CSANK→NK ladder, switch tols, non-negotiable options like `turbResScale [1e4, 0.1, 1e-4]` and `ADPC`), measured limits (CSANK floor ~3.5e-8, deep-NK preconditioner wall), the 2026-08-08 corrections (premature-NK-engagement rule; Eisenstat-Walker mitigation falsified), and the index of every acceleration test. |
| [ADFLOW_BASE/ADFLOW_07_sst_dev_lessons.md](ADFLOW_BASE/ADFLOW_07_sst_dev_lessons.md) | Post-mortem of upstream `sst_dev` (SST, PR #331) — how a 2-equation model was differentiated, what broke, and a comparison table vs this branch's SA-GR. Read before any AD-unfreeze work or partials campaign. |
| [CORE_BASE/CORE_03_design_decisions.md](CORE_BASE/CORE_03_design_decisions.md) | **Not a spec — a memory.** Condensed log of resolved code-audit questions (A1-A3): nondim safeguards, convergence-strategy decisions (Eq. 59 forms, damping, `turbResScale`), and SA/SST code-coherence divergences (ft2 default, wall functions, source-Jacobian clips, farfield BC, init values). Read for "why is X implemented this way" — code/paper win if it ever disagrees with this file. |

### Verification (`docs/VERIFICATION/`)

| File | What's in it |
|------|--------------|
| [VERIFICATION/VERIF_00_three_stage_verification.md](VERIFICATION/VERIF_00_three_stage_verification.md) | The 3-stage low-level adjoint verification ladder: (1) reverse↔forward dot-product consistency, (2) reverse vs fast-reverse (`_b` vs `_b_fast`) consistency, (3) 3-way AD/FD/CS forward-mode check. **Start at §"Canonical way to run"**: the ladder is now a registered testflo suite (`tests/reg_tests/test_jacVecProd{FWD,BWDFast}_sagr.py` + `test_adjoint_sagr.py` + `test_blockette_sagr.py` + `reg_sagr.py` + `refs/*.json`) driven by `run_sagr_tests.sh` (`all\|real\|cs\|adjoint\|blockette\|train\|genw`); the suite pins `transitioncrossflow: False` (the runtime default since 2026-07-24). Covers the two invariants not to "fix" (FD `@expectedFailure`; CS state re-seat), the mesh-swap recipe, and results. Old per-script commands are retained as historical (scripts now in `dev/`). |
| [VERIFICATION/VERIF_03_paper_validation_campaign.md](VERIFICATION/VERIF_03_paper_validation_campaign.md) | Consolidated S809/NLF0416 physics-validation record vs Piotrowski & Zingg (2026-08-03/04): exact-paper-grid experiment (+22/+30 counts above bucket ⇒ implementation behavior, not grids/geometry), 10 hypotheses tested & eliminated (incl. λθ clamp, dissipation, C-mesh/α), full clean symbol audit, digitized paper cd/x_tr targets, corner non-uniqueness + late-NK instability findings, defined next step (x_tr slice comparison). |
| [ADFLOW_BASE/ADFLOW_09_adjoint_trace.md](ADFLOW_BASE/ADFLOW_09_adjoint_trace.md) | Paired inventory of adjoint/AD touchpoints (SA vs SA-GR) on this branch: preprocessor guards, Tapenade directives, generated AD files, wiring. Also covers the hand-written adjoint-solver `storePsiHistory` derivative-convergence diagnostic (opt-in, off by default) and its `.pyf` wiring gotcha. |
| [VERIFICATION/VERIF_02_rotating_frame_audit.md](VERIFICATION/VERIF_02_rotating_frame_audit.md) | Rotating-frame consistency fix (2026-07-23): makes the transition model use the **relative** frame — vortLim reference velocity `√(uInf²+\|Ω×r\|²)`, BL velocity scales and helicity on `V_rel`. Verified bit-identical no-op when Ω=0 (residual + full derivative suite). Flags TAPENADE-NEEDED for rotating-case adjoints. |
| [VERIFICATION/VERIF_01_debugging_derivatives.md](VERIFICATION/VERIF_01_debugging_derivatives.md) | Raw source: generic FD → reverse-mode dot-product → complex-step derivative-checking ladder, not covered in `ADFLOW_BASE/02`/`04`. |
| [VERIFICATION/VERIF_04_ank_nk_transition_mods_report.md](VERIFICATION/VERIF_04_ank_nk_transition_mods_report.md) | Historical report: the SA-GR-specific modifications made to the ANK/NK solvers (column scaling, Eq. 59, damping) at the time they were introduced. |
| [VERIFICATION/VERIF_05_ank_nk_complexify_report.md](VERIFICATION/VERIF_05_ank_nk_complexify_report.md) | Historical report: complexification of the ANK/NK transition mods for the CS build. |

### Archive (`docs/ARCHIVE/`)

Historical/point-in-time records, superseded by current docs — trust the
dated snapshot, verify against current code before acting on specifics.
`CORE_BASE/CORE_03_design_decisions.md` and `ADFLOW_BASE/ADFLOW_07_sst_dev_lessons.md` used to live
here too but are actively routed-to, so they moved out instead.

| File | What's in it |
|------|--------------|
| [ARCHIVE/ARCHIVE_03_adjoint_audit_2026-07-07.md](ARCHIVE/ARCHIVE_03_adjoint_audit_2026-07-07.md) | Pre-partials-test visual audit of the SA-GR adjoint wiring: verified-pass table, the `vortlimd = 0` finding (uInf/muInf head activity), and watch items (autoEdit fast_b stripping, limiter kinks). |
| [ARCHIVE/ARCHIVE_00_inventory.md](ARCHIVE/ARCHIVE_00_inventory.md) | Historical (2026-07): implementation-state inventory of the transition model at audit time. |
| [ARCHIVE/ARCHIVE_01_adjoint_wiring.md](ARCHIVE/ARCHIVE_01_adjoint_wiring.md) | Historical (2026-07): deep audit of the SA-GR adjoint wiring that fed the fixes later verified in `VERIFICATION/VERIF_00_three_stage_verification.md`. |
| [ARCHIVE/ARCHIVE_02_test_prep.md](ARCHIVE/ARCHIVE_02_test_prep.md) | Historical (2026-07): preparation notes for the SA-GR test suite (now realized as `tests/reg_tests/README_SAGR.md`). |

### Physics reference (`docs/SA_GAMMA_RETHETHA_BASE/`)

| File | What's in it |
|------|--------------|
| [SA_GAMMA_RETHETHA_BASE/SAGR_01_paper_piotrowski_zingg_2020.md](SA_GAMMA_RETHETHA_BASE/SAGR_01_paper_piotrowski_zingg_2020.md) | Full paper text. **Source of truth for physics** — when code and paper disagree, paper wins. Sole physics reference (distilled summaries were retired after repeated distillation errors; see `CORE_BASE/CORE_03_design_decisions.md` §D3). |
| [SA_GAMMA_RETHETHA_BASE/SAGR_00_context_index.md](SA_GAMMA_RETHETHA_BASE/SAGR_00_context_index.md) | Sub-index for the transition physics KB. |
| [SA_GAMMA_RETHETHA_BASE/SAGR_02_adflow_vs_paper_solver.md](SA_GAMMA_RETHETHA_BASE/SAGR_02_adflow_vs_paper_solver.md) | How ADflow's solver hierarchy (ANK/CANK/CSANK/NK, global-λ step control) differs from the paper's Newton–Krylov–Schur algorithm (Eq. 58 scaling, Alg. 2 per-node damping, Eq. 59 source-Δt + reactivation), why the paper converges in fewer iterations, what's ported, and the open code items. Grounded in the 2026-07-14/15 test campaign. |

### ADflow reference KB (`docs/ADFLOW_BASE/`)

| File | What's in it |
|------|--------------|
| [ADFLOW_00_context_index.md](ADFLOW_BASE/ADFLOW_00_context_index.md) | Sub-index + provenance for the ADflow KB. |
| [ADFLOW_01_flow_solver_theory.md](ADFLOW_BASE/ADFLOW_01_flow_solver_theory.md) | Flow solver, ANK/NK, PTC, R0/R1/R2, defaults — extraction of Yildirim et al. (2019). |
| [ADFLOW_02_adjoint_autodiff_theory.md](ADFLOW_BASE/ADFLOW_02_adjoint_autodiff_theory.md) | Discrete adjoint + automatic differentiation — extraction of Kenway et al. (2019). |
| [ADFLOW_03_concordance.md](ADFLOW_BASE/ADFLOW_03_concordance.md) | Paper-math ↔ code-flag crosswalk + gotchas. |
| [ADFLOW_04_debugging_playbook.md](ADFLOW_BASE/ADFLOW_04_debugging_playbook.md) | Symptom → cause → fix ladder → evidence for stalls, divergence, bad gradients. |
| [ADFLOW_05_options_and_operations_devguide.md](ADFLOW_BASE/ADFLOW_05_options_and_operations_devguide.md) | ADflow options + operations (official docs). |
| [ADFLOW_06_official_solvers_doc.md](ADFLOW_BASE/ADFLOW_06_official_solvers_doc.md) | Raw official-docs source that fed `01`/`04`/`05` — its content is ~fully duplicated in `ADFLOW_05` Part 1; upstream advice, carries an SA-GR override banner. |
| [ADFLOW_07_sst_dev_lessons.md](ADFLOW_BASE/ADFLOW_07_sst_dev_lessons.md) | Branch-specific, not upstream: post-mortem of `sst_dev` (SST, PR #331) — see File index above. |
| [ADFLOW_08_nondimensionalization.md](ADFLOW_BASE/ADFLOW_08_nondimensionalization.md) | Branch-specific, not upstream: ADflow's p-ρ non-dimensional scaling — see File index above. |
| [ADFLOW_09_adjoint_trace.md](ADFLOW_BASE/ADFLOW_09_adjoint_trace.md) | Branch-specific, not upstream: adjoint/AD touchpoint inventory — see Verification section above. |

### Test-suite docs (`tests/reg_tests/`)

| File | What's in it |
|------|--------------|
| [../tests/reg_tests/README_SAGR.md](../tests/reg_tests/README_SAGR.md) | Guide to the SA-GR derivative test suite: file tree, what each test checks, run order, failure ladder. |
| [../tests/reg_tests/dev/README.md](../tests/reg_tests/dev/README.md) | Why `dev/` scripts sit outside testflo; restart regeneration + diagnostic tooling inventory. |

---

## Sources

- Piotrowski & Zingg (2020), *AIAA Journal* 58(10) — γ-Re̅θt transition model
  → `SA_GAMMA_RETHETHA_BASE/`
- Yildirim et al., *J. Comput. Phys.* (2019), doi:10.1016/j.jcp.2019.06.018 → `ADFLOW_01`
- Kenway et al., *Prog. Aerosp. Sci.* (2019), doi:10.1016/j.paerosci.2019.05.002 → `ADFLOW_02`
- ADflow official docs (MDO Lab, retrieved 6 Jul 2026) → `ADFLOW_05`

Project rules and the task roadmap live in [`../CLAUDE.md`](../CLAUDE.md).
