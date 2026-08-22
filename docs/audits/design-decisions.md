# Design Decisions Log — SA-BCM

> **Not a defining/normative file.** This is a *memory* of what was already discussed and
> decided during code audits — not a spec, not a source of truth, and nothing here overrides
> the paper (physics) or the code itself. If this log and the current code/paper ever disagree,
> trust the code/paper and treat this file as stale — update it, don't defer to it. Its only job
> is to save a future reader from re-litigating a question that was already asked and answered,
> and to say *why* a piece of code looks the way it does.

## Repo setup / project scope

### sa.F90 edit-scope rule (kept, with a caveat)

- Question: an earlier project memory forbade any edit to `sa.F90`. SA-BCM hand-edits `sa.F90`
  directly (no separate model file), so the rule as literally stated can't hold on this branch.
- Decision: scope the rule to the `use_SABCM`-gated block only (`sa.F90:294-413`,
  `blockette.F90:1004-1250`) rather than refactoring SA-BCM out into its own file first. The
  extraction option was considered and explicitly deferred — see `../TODO.md`.
- Why: the user confirmed the original "never touch sa.F90" memory was written for a different
  branch (a γ-Re̅θt transition model that *did* live in its own file) and doesn't reflect this
  branch's actual structure.

### No new transport-equation state variables

- Question: does SA-BCM add new `nw` entries like a 2-equation transition model would (γ-Re̅θt
  precedent), requiring coupling-block Jacobian test machinery?
- Decision: no. SA-BCM multiplies the existing SA production term by `tTgamma` and overrides
  `ft2`; `nw` is unchanged. The test harness (`tests/reg_tests/reg_bcm.py`) therefore skips
  state-vector-block splitting entirely and compares SA-BCM's modified Jacobian terms directly
  against plain-SA's (`use_SABCM=False`) unmodified terms instead.
- Why: confirmed by reading `sa.F90:294-413` — no new residual array, no new `itu*` index
  anywhere in the module declarations.

### Physics source of truth: two papers, not one

- Question: which document is the physics reference for `CLAUDE.md` rule 9?
- **Superseded 2026-07-24**: the original decision below named an unpublished internal manuscript
  (`SA-BCM_Differentiable_Reformulation.md`) as primary. That file has been replaced in
  `docs/papers/` by the two actual published AIAA papers — `AIAA20202714_SABCMPartI.md` (Mura &
  Cakmakcioglu; its Appendix is now the authoritative "copy-for-code" formulation) and
  `AIAA20202706_BCMtransitionmodel.md` (Çakmakçıoğlu, Baş, Mura, Kaynak; secondary/companion,
  SU2-notation reference — see that file's "Implementer's reconciliation notes"). The
  `SABCM_Exp=True` ("hard") path is the literal formula from these papers; `SABCM_Exp=False`
  ("smooth", tanh) is a deliberate smoothing **not present in either paper**, kept by explicit
  user choice, confirmed by a line-by-line paper-vs-code check (`docs/adjoint-trace.md`). The
  original decision's premise (below) — an unpublished manuscript's "needs revalidation" caveat —
  no longer has a source in this repo; see `docs/README.md`'s "Papers were replaced" note.

**Original decision (2026-02, retained for history):**
- Decision: `docs/papers/SA-BCM_Differentiable_Reformulation.md` (Cardoso, Marta & Martins) is
  primary — it's what the smoothed/default (`SABCM_Exp=False`) code path should implement.
  `docs/papers/SA-BCM_Transitional_Model.md` (Mura & Cakmakcioglu, AIAA 2020-2714) is secondary
  — it's the original formulation the `SABCM_Exp=True` path implements, and what the
  differentiable reformulation is a reformulation *of*.
- Why: the differentiable-reformulation manuscript is unpublished, author-provided, and carries
  an explicit "adjoint results need revalidation" caveat that is the actual reason this branch's
  work exists — it's not just background reading, it's the task spec.
