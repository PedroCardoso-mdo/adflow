# ADflow KB — sub-index

Reference extraction of the published ADflow papers + official docs. **Load only
the file you need — don't read the whole KB.**

Full task → files routing (including the project-specific docs) lives in the
master index: [`../README.md`](../README.md). Project rules are in `../../CLAUDE.md`.

## Files
| File | What it is |
|---|---|
| `ADFLOW_00_context_index.md` | this sub-index |
| `ADFLOW_01_flow_solver_theory.md` | flow solver, ANK/NK, PTC, R0/R1/R2, defaults — full extraction of Yildirim2019 |
| `ADFLOW_02_adjoint_autodiff_theory.md` | discrete adjoint + AD, variants, appendices — full extraction of Kenway2019 |
| `ADFLOW_03_concordance.md` | paper-math ↔ code-flag crosswalk + gotchas |
| `ADFLOW_04_debugging_playbook.md` | symptom → cause → fix ladder → evidence for stalls, divergence, bad gradients |
| `ADFLOW_05_options_and_operations_devguide.md` | ADflow options + operations (official docs; user-provided) |

The full transition-model paper lives in `../SA_GAMMA_RETHETHA_BASE/`.

## Sources
- Yildirim et al., *J. Comput. Phys.* (2019), doi:10.1016/j.jcp.2019.06.018 → `01`
- Kenway et al., *Prog. Aerosp. Sci.* (2019), doi:10.1016/j.paerosci.2019.05.002 → `02`
- ADflow docs (MDO Lab, retrieved 6 Jul 2026) → `05`
