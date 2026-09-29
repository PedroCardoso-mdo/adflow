# docs/ — index

Load only the file(s) the routing row names, in order.

## Routing

| Task                                                        | Read (in order)                                                                 |
|-------------------------------------------------------------|---------------------------------------------------------------------------------|
| Transition equation / constant / algorithm (SA-GR)          | `SA_GAMMA_RETHETHA_BASE/SAGR_01_paper_piotrowski_zingg_2020.md` → `CORE_BASE/CORE_01_architecture.md` |
| SA-BCM physics                                              | `papers/AIAA20202714_SABCMPartI.md` (Appendix) → `CORE_01` (SA-BCM part)        |
| SA-sγ                                                       | `CORE_01` (SA-sγ part) → external design docs it points to                      |
| Equation with velocity / viscosity / Re / time scale        | `ADFLOW_BASE/ADFLOW_08_nondimensionalization.md` **first** → paper              |
| Runtime option (name / default / meaning)                   | `CORE_BASE/CORE_01_architecture.md`; generic ADflow options: `../doc/options.yaml` |
| Code map, state vector, solver paths, known limitations     | `CORE_BASE/CORE_01_architecture.md`                                             |
| Converge a case / debug a stalling or diverging run         | `CORE_BASE/CORE_02_convergence_strategy.md` → `VERIFICATION/VERIF_06_solver_fidelity_audit.md` |
| A code comment cites an `F<n>` tag                          | `VERIFICATION/VERIF_06_solver_fidelity_audit.md`                                |
| Adjoint / Tapenade wiring, regenerating AD                  | `ADFLOW_BASE/ADFLOW_09_adjoint_trace.md`                                        |
| Verify derivatives / gradients wrong                        | `VERIFICATION/VERIF_00_three_stage_verification.md` (skill `adflow-derivative-verification`) → `ADFLOW_09` |
| Rotating frame (`rotRate ≠ 0`)                              | `VERIFICATION/VERIF_02_rotating_frame_audit.md`                                 |
| Generic ADflow theory (ANK/NK, adjoint)                     | upstream `../doc/*.rst`; Yildirim et al. 2019 (JCP), Kenway et al. 2019 (PAS)   |

## Files

| File                                                   | Holds                                                                  |
|--------------------------------------------------------|------------------------------------------------------------------------|
| `CORE_BASE/CORE_01_architecture.md`                    | Code map, every transition option, design rationale, limitations, SA-sγ and SA-BCM parts |
| `CORE_BASE/CORE_02_convergence_strategy.md`            | Validated solver recipes per model/case, falsified levers, ADflow vs paper solver |
| `ADFLOW_BASE/ADFLOW_08_nondimensionalization.md`       | ADflow's p-ρ scaling                                                   |
| `ADFLOW_BASE/ADFLOW_09_adjoint_trace.md`               | AD heads, master-routine dispatch, Tapenade gotchas                    |
| `VERIFICATION/VERIF_00_three_stage_verification.md`    | Derivative ladder: what each stage proves, how to run, invariants      |
| `VERIFICATION/VERIF_02_rotating_frame_audit.md`        | Relative-frame treatment of the transition terms                       |
| `VERIFICATION/VERIF_06_solver_fidelity_audit.md`       | F-tag → option → default → reason (cited from code)                    |
| `SA_GAMMA_RETHETHA_BASE/SAGR_01_…md`                   | Piotrowski & Zingg 2020, full text — physics source of truth           |
| `papers/`                                              | SA-BCM papers (Mura & Cakmakcioglu 2020; Çakmakçıoğlu et al. 2020)     |
| `MICHAEL_PIOTROWSKI/MP_0{1,3,4,5,6}_*_full.md`         | Transcripts of the later Piotrowski/Zingg papers and thesis (only local copies) |

Test-suite READMEs: `../tests/reg_tests/README_SAGR.md`, `README_BCM.md`, `dev/README.md`.
