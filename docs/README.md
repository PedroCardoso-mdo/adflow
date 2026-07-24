# Documentation & Knowledge Base — SA-BCM Differentiable Transition Model

Master index for all project documentation.

> **Read only the file(s) a task needs — never load the whole KB.**
> (CLAUDE.md rule 8.)

## Task → files routing

| Your task / question | Read (in order) |
|---|---|
| Implement/extend/debug an SA-BCM term (γ, Term1, Term2, KS max, tanh blend) | `papers/AIAA20202714_SABCMPartI.md` (Appendix — authoritative "copy-for-code" formulation) → `architecture.md` |
| Understand the original BCM model / SU2 implementation notation | `papers/AIAA20202706_BCMtransitionmodel.md` (and the "Implementer's reconciliation notes" at the end of the 2020-2714 file — several symbols differ between the two papers) |
| Any equation touching velocity, vorticity, viscosity ratio, or Reynolds number | `nondimensionalization.md` (**first**) → `papers/AIAA20202714_SABCMPartI.md` |
| Exact runtime option name / default / enum | `architecture.md` |
| Solver architecture, code locations, paper-symbol ↔ code-flag mapping | `architecture.md` |
| Adjoint / AD / Tapenade status on this branch | `adjoint-trace.md` |
| Plan or run the derivative-consistency test campaign (FD/CS/adjoint) | `current-task.md` → `tests/reg_tests/README_BCM.md` → `audits/` |
| Why is X implemented this way / was this already discussed | `audits/design-decisions.md` (memory of past discussion, not a spec — code/paper win if it disagrees) |
| What's the task in progress right now | `current-task.md` |
| What was already finished, and how | `task-log/README.md` |
| Deferred tuning/scope items not urgent right now | `TODO.md` |

## File index

| File | Description |
|---|---|
| `papers/AIAA20202714_SABCMPartI.md` | Mura & Cakmakcioglu, AIAA Aviation 2020-2714 — "A Revised One-Equation Transitional Model for External Aerodynamics, Part I". **Primary physics source of truth**; the Appendix is the authoritative, implementation-ready formulation. Includes an "Implementer's reconciliation notes" section reconciling notation/constants against the companion paper below. |
| `papers/AIAA20202706_BCMtransitionmodel.md` | Çakmakçıoğlu, Baş, Mura, Kaynak, AIAA Aviation 2020-2706 — "A Revised One-Equation Transitional Model for External Aerodynamics" (SU2 implementation, companion/shorter paper). Secondary reference; same model, different notation in places (see reconciliation notes above). |
| `architecture.md` | Solver architecture, `saSource` code path, every SA-BCM runtime option (name/default/enum), paper-symbol ↔ code-flag lookup table. |
| `nondimensionalization.md` | ADflow's p-ρ non-dimensionalization convention, relevant wherever SA-BCM mixes velocity/vorticity and viscosity-ratio terms. |
| `adjoint-trace.md` | AD/Tapenade touchpoints for the `use_SABCM` guard: which files are generated, sync status, what triggers a rerun. |
| `current-task.md` | The ONE active task, overwritten each session. |
| `TODO.md` | Deferred, decided-not-urgent items, each linking to its `audits/design-decisions.md` entry. |
| `task-log/README.md` | Index of finished tasks + template. |
| `audits/design-decisions.md` | Consolidated, topic-organized log of resolved audit questions. Explicitly non-normative. |

## Papers were replaced 2026-07-24

The two files in `papers/` above replaced an earlier, unpublished "differentiable reformulation"
manuscript that used to be this KB's stated physics source of truth and carried an explicit
"adjoint-based sensitivity results must be revalidated" caveat (the stated reason the
derivative-test harness in `tests/reg_tests/` exists). That manuscript is no longer in this repo,
so that specific caveat no longer has a source here — **do not assume the underlying bug it
described is either present or fixed; nothing in this repo currently documents what it was.**
What IS established: a line-by-line check of the current `sa.F90`/`blockette.F90` implementation
against the two published AIAA papers now here found no discrepancies (see `adjoint-trace.md`'s
"Paper-vs-code verification" section) — implementation-vs-published-model is faithful. Whether
that also resolves whatever the retired manuscript's caveat was about is unknown; the
derivative-test harness (once actually run, see `current-task.md`) is the real answer to that,
not this note.

## Retired-doc lesson (read once, apply always)

A prior branch's first physics reference was a hand-written "distilled" summary of its paper,
not the paper itself. Over two audit passes it accumulated 5 transcription errors — the code
was correct every time, only the distilled doc was wrong. **Never create a "distilled" or
"summary" physics doc as the working reference.** `CLAUDE.md` rule 9 and every routing-table row
above point directly at the full paper text in `papers/`. `architecture.md`'s paper-symbol ↔
code-flag table is a lookup table, not a re-derivation of physics content — keep it that way.

Project rules and the task roadmap live in [`../CLAUDE.md`](../CLAUDE.md).
