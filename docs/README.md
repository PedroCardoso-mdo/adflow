# Documentation & Knowledge Base — SA-γ-Re̅θt Transition Model

Master index for **all** project documentation. This branch carries two kinds
of docs:

- **Project-specific** (`docs/*.md`) — unique to the Piotrowski & Zingg (2020)
  γ-Re̅θt implementation on this branch.
- **Reference KB** (`docs/ADFLOW_BASE/`, `docs/SA_GAMMA_RETHETHA_BASE/`) — distilled
  extractions of the published papers and the official ADflow docs.

> **Read only the file(s) a task needs — never load the whole KB.** Every file
> is self-contained on purpose (token discipline; see `../CLAUDE.md` rule 8).
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
| Implement/extend a transition equation, constant, or algorithm | `SA_GAMMA_RETHETHA_BASE/Piotrowski_Zingg_2020_SA-sLM2015_clean (1).md` → `architecture.md` |
| Any equation touching velocity, viscosity, Reynolds number, or time scales | `nondimensionalization.md` (**first**) → `SA_GAMMA_RETHETHA_BASE/Piotrowski_Zingg_2020_SA-sLM2015_clean (1).md` |
| Paper symbol ↔ code flag (transition model) | `SA_GAMMA_RETHETHA_BASE/Piotrowski_Zingg_2020_SA-sLM2015_clean (1).md` + `architecture.md` (Part 2) |
| Exact **transition** runtime option name / default / enum | `architecture.md` (Part 2) |
| Solver architecture, state-vector layout, code/module locations | `architecture.md` (Part 1) |
| Flow / ANK / NK theory, equations, defaults | `ADFLOW_BASE/ADFLOW_01_flow_solver_theory.md` |
| What an ANK/NK option does mechanically | `ADFLOW_BASE/ADFLOW_03_concordance.md` → `ADFLOW_BASE/ADFLOW_01_flow_solver_theory.md` |
| Choose R0/R1/R2, coupled vs decoupled, NK switch | `ADFLOW_BASE/ADFLOW_01_flow_solver_theory.md` → `ADFLOW_BASE/ADFLOW_04_debugging_playbook.md` |
| Exact **ADflow** (flow/ANK/NK/adjoint) option name / default / enum | `ADFLOW_BASE/ADFLOW_05_options_and_operations_devguide.md` |
| Paper symbol ↔ code flag (ADflow solver) | `ADFLOW_BASE/ADFLOW_03_concordance.md` |
| Debug a stalling / diverging run | `ADFLOW_BASE/ADFLOW_04_debugging_playbook.md` → `ADFLOW_BASE/ADFLOW_05_options_and_operations_devguide.md` → `ADFLOW_BASE/ADFLOW_01_flow_solver_theory.md` |
| Rotating mesh / rotating frame of reference (rotor, propeller, turbine) | `VERIFICATION/rotating-frame-audit.md` → `nondimensionalization.md` (§5 vortLim, §6 crossflow) |
| Why does the paper converge faster than ADflow / solver-algorithm gaps | `SA_GAMMA_RETHETHA_BASE/adflow-vs-paper-solver.md` |
| Adjoint / AD theory (general) | `ADFLOW_BASE/ADFLOW_02_adjoint_autodiff_theory.md` |
| Implement or extend the adjoint / AD on **this branch** | `VERIFICATION/adjoint-trace.md` → `ADFLOW_BASE/ADFLOW_02_adjoint_autodiff_theory.md` → `ADFLOW_BASE/ADFLOW_03_concordance.md` |
| Gradients wrong / adjoint won't converge | `ADFLOW_BASE/ADFLOW_04_debugging_playbook.md` → `VERIFICATION/debugging_derivatives.md` (FD→dot-product→CS checking ladder) → `ADFLOW_BASE/ADFLOW_02_adjoint_autodiff_theory.md` → `VERIFICATION/adjoint-trace.md` |
| Validate partials / plan AD test campaign | `current-task.md` → `VERIFICATION/debugging_derivatives.md` (generic checking ladder) → `VERIFICATION/three-stage-verification.md` (tests already run, commands, results) → `../audits/adjoint_audit_2026-07-07.md` → `../audits/07_sst_dev_lessons.md` (verification ladder + watch items) |
| Rerun/reproduce the low-level adjoint verification tests (dot-product, fast_b, 3-way fwd) | `VERIFICATION/three-stage-verification.md` (§"Canonical way to run" → `tests/reg_tests/run_sagr_tests.sh`) |
| Run / retrain / regen-`w` for the SA-GR partials tests (how-to) | `VERIFICATION/three-stage-verification.md` §"Canonical way to run" → `tests/reg_tests/dev/README.md` (non-standard `w` generation) |
| How another multi-equation turb model was differentiated (SST precedent) | `../audits/07_sst_dev_lessons.md` |
| Why is X implemented this way / was this already discussed | `../audits/design-decisions.md` (memory of past discussion, not a spec — code/paper win if it disagrees) |
| What's the task in progress right now | `current-task.md` |
| What was already finished, and how | `task-log/README.md` (index → per-task case files) |
| Provenance / sources / where a fact came from | this file + `ADFLOW_BASE/ADFLOW_00_context_index.md` |

---

## File index

### Project-specific (`docs/`)

| File | What's in it |
|------|--------------|
| [README.md](README.md) | This master index + routing table. |
| [architecture.md](architecture.md) | Solver architecture, state-vector layout, key code/module locations, user constraints, and the complete reference for every runtime option added for the transition model. |
| [nondimensionalization.md](nondimensionalization.md) | How ADflow makes the governing equations dimensionless — the **pressure–density (p-ρ) scaling** (velocity normalizes to M·√γ, *not* 1). Read before touching any equation with velocity, viscosity, time scales, or Reynolds number. |
| [TODO.md](TODO.md) | Deferred tuning/improvement items (decided, not defects): `turbResScale` calibration, Eq. 59 relaxation options, damping-clip fallback, open dúvida D-A2-3, blockettes, test-infra restructure. Each item links back to its `audits/design-decisions.md` analysis. |
| [current-task.md](current-task.md) | The single task currently in progress (CLAUDE.md: one task per session) — objective, scoped context, working files, checklist. Overwritten each session. |
| [task-log/](task-log/README.md) | Finished-task log, one file per task (a "case"), added over time — index + template in `task-log/README.md`. |

### Verification (`docs/VERIFICATION/`)

| File | What's in it |
|------|--------------|
| [VERIFICATION/three-stage-verification.md](VERIFICATION/three-stage-verification.md) | The 3-stage low-level adjoint verification ladder: (1) reverse↔forward dot-product consistency, (2) reverse vs fast-reverse (`_b` vs `_b_fast`) consistency, (3) 3-way AD/FD/CS forward-mode check. **Start at §"Canonical way to run"**: the ladder is now a registered testflo suite (`tests/reg_tests/test_jacVecProd{FWD,BWDFast}_sagr.py` + `reg_sagr.py` + `refs/*.json`) driven by `run_sagr_tests.sh` (`all\|real\|cs\|train\|genw`), crossflow always ON. Covers the two invariants not to "fix" (FD `@expectedFailure`; CS state re-seat), the mesh-swap recipe, and results. Old per-script commands are retained as historical (scripts now in `dev/`). |
| [VERIFICATION/adjoint-trace.md](VERIFICATION/adjoint-trace.md) | Paired inventory of adjoint/AD touchpoints (SA vs SA-GR) on this branch: preprocessor guards, Tapenade directives, generated AD files, wiring. |
| [VERIFICATION/rotating-frame-audit.md](VERIFICATION/rotating-frame-audit.md) | Rotating-frame consistency fix (2026-07-23): makes the transition model use the **relative** frame — vortLim reference velocity `√(uInf²+\|Ω×r\|²)`, BL velocity scales and helicity on `V_rel`. Verified bit-identical no-op when Ω=0 (residual + full derivative suite). Flags TAPENADE-NEEDED for rotating-case adjoints. |
| [VERIFICATION/debugging_derivatives.md](VERIFICATION/debugging_derivatives.md) | Raw source: generic FD → reverse-mode dot-product → complex-step derivative-checking ladder, not covered in `ADFLOW_BASE/02`/`04`. |

### Audits (`audits/`)

| File | What's in it |
|------|--------------|
| [../audits/07_sst_dev_lessons.md](../audits/07_sst_dev_lessons.md) | Post-mortem of upstream `sst_dev` (SST, PR #331) — how a 2-equation model was differentiated, what broke, and a comparison table vs this branch's SA-GR. Read before any AD-unfreeze work or partials campaign. |
| [../audits/adjoint_audit_2026-07-07.md](../audits/adjoint_audit_2026-07-07.md) | Pre-partials-test visual audit of the SA-GR adjoint wiring: verified-pass table, the `vortlimd = 0` finding (uInf/muInf head activity), and watch items (autoEdit fast_b stripping, limiter kinks). |
| [../audits/design-decisions.md](../audits/design-decisions.md) | **Not a spec — a memory.** Condensed log of resolved code-audit questions (A1-A3): nondim safeguards, convergence-strategy decisions (Eq. 59 forms, damping, `turbResScale`), and SA/SST code-coherence divergences (ft2 default, wall functions, source-Jacobian clips, farfield BC, init values). Read for "why is X implemented this way" — code/paper win if it ever disagrees with this file. |

### Physics reference (`docs/SA_GAMMA_RETHETHA_BASE/`)

| File | What's in it |
|------|--------------|
| Piotrowski_Zingg_2020_SA-sLM2015_clean (1).md | Full paper text. **Source of truth for physics** — when code and paper disagree, paper wins. Sole physics reference (distilled summaries were retired after repeated distillation errors; see `../audits/design-decisions.md` §D3). |
| [adflow-vs-paper-solver.md](SA_GAMMA_RETHETHA_BASE/adflow-vs-paper-solver.md) | How ADflow's solver hierarchy (ANK/CANK/CSANK/NK, global-λ step control) differs from the paper's Newton–Krylov–Schur algorithm (Eq. 58 scaling, Alg. 2 per-node damping, Eq. 59 source-Δt + reactivation), why the paper converges in fewer iterations, what's ported, and the open code items. Grounded in the 2026-07-14/15 test campaign. |

### ADflow reference KB (`docs/ADFLOW_BASE/`)

| File | What's in it |
|------|--------------|
| [ADFLOW_00_context_index.md](ADFLOW_BASE/ADFLOW_00_context_index.md) | Sub-index + provenance for the ADflow KB. |
| [ADFLOW_01_flow_solver_theory.md](ADFLOW_BASE/ADFLOW_01_flow_solver_theory.md) | Flow solver, ANK/NK, PTC, R0/R1/R2, defaults — extraction of Yildirim et al. (2019). |
| [ADFLOW_02_adjoint_autodiff_theory.md](ADFLOW_BASE/ADFLOW_02_adjoint_autodiff_theory.md) | Discrete adjoint + automatic differentiation — extraction of Kenway et al. (2019). |
| [ADFLOW_03_concordance.md](ADFLOW_BASE/ADFLOW_03_concordance.md) | Paper-math ↔ code-flag crosswalk + gotchas. |
| [ADFLOW_04_debugging_playbook.md](ADFLOW_BASE/ADFLOW_04_debugging_playbook.md) | Symptom → cause → fix ladder → evidence for stalls, divergence, bad gradients. |
| [ADFLOW_05_options_and_operations_devguide.md](ADFLOW_BASE/ADFLOW_05_options_and_operations_devguide.md) | ADflow options + operations (official docs). |
| [adflow_solvers.md](ADFLOW_BASE/adflow_solvers.md) | Raw official-docs source that fed `01`/`04`/`05`; last-resort read for troubleshooting detail not otherwise distilled. |

---

## Sources

- Piotrowski & Zingg (2020), *AIAA Journal* 58(10) — γ-Re̅θt transition model
  → `SA_GAMMA_RETHETHA_BASE/`
- Yildirim et al., *J. Comput. Phys.* (2019), doi:10.1016/j.jcp.2019.06.018 → `ADFLOW_01`
- Kenway et al., *Prog. Aerosp. Sci.* (2019), doi:10.1016/j.paerosci.2019.05.002 → `ADFLOW_02`
- ADflow official docs (MDO Lab, retrieved 6 Jul 2026) → `ADFLOW_05`

Project rules and the task roadmap live in [`../CLAUDE.md`](../CLAUDE.md).
