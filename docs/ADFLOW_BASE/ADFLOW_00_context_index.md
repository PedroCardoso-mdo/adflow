# ADflow KB — sub-index

Reference extraction of the published ADflow papers + official docs. **Load only
the file you need — don't read the whole KB.**

Full task → files routing (including the project-specific docs) lives in the
master index: [`../README.md`](../README.md). Project rules are in `../../CLAUDE.md`.

## Files
| File | What it is |
|---|---|
| `ADFLOW_00_context_index.md` | this sub-index |
| `ADFLOW_01_flow_solver_theory.md` | flow solver, ANK/NK, PTC, R0/R1/R2, paper defaults (Yildirim et al. 2019) — full extraction of Yildirim2019; several defaults differ from current code — code/pyADflow.py wins |
| `ADFLOW_02_adjoint_autodiff_theory.md` | discrete adjoint + AD, variants, appendices — full extraction of Kenway2019 |
| `ADFLOW_03_concordance.md` | paper-math ↔ code-flag crosswalk + gotchas |
| `ADFLOW_04_debugging_playbook.md` | symptom → cause → fix ladder → evidence for stalls, divergence, bad gradients |
| `ADFLOW_05_options_and_operations_devguide.md` | ADflow options + operations (official docs; user-provided) |
| `ADFLOW_06_official_solvers_doc.md` | **Raw source** (official ADflow docs, un-numbered) that fed `01`/`04`/`05`. ~Fully duplicated by `ADFLOW_05` Part 1 (Very Small Step Sizes / High Linear Residuals / EW Stalling are all there too), but this file IS what CLAUDE.md routes to for SA-GR convergence background. *(Corrected 2026-08-12: previously claimed unique troubleshooting detail.)* |
| `../VERIFICATION/VERIF_01_debugging_derivatives.md` | **Raw source**, mostly un-distilled. Holds the generic FD → reverse-mode dot-product → complex-step **derivative-checking ladder** — not covered in `02` or `04`. Referenced by the master index's "Validate partials / plan AD test campaign" and "Gradients wrong / adjoint won't converge" routing rows. |

The full transition-model paper lives in `../SA_GAMMA_RETHETHA_BASE/`.

## Sources
- Yildirim et al., *J. Comput. Phys.* (2019), doi:10.1016/j.jcp.2019.06.018 → `01`
- Kenway et al., *Prog. Aerosp. Sci.* (2019), doi:10.1016/j.paerosci.2019.05.002 → `02`
- ADflow docs (MDO Lab, retrieved 6 Jul 2026) → `05`, `ADFLOW_06_official_solvers_doc.md`, `debugging_derivatives.md`
