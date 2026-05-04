# Task Template

The roadmap in `03_IMPLEMENTATION_PLAN.md` covers all planned tasks. Use
this template only when a NEW task surfaces (e.g. a bug found during
validation after E2). Copy-paste, fill in, and start a fresh Claude Code
session.

---

## Session Prompt

```
@CLAUDE.md
@02_IMPLEMENTATION_STATUS.md

## Task: <T-ID> — <one-line title>

Goal:
<one or two sentences: what changes and why>

AD impact: <yes / no>
   yes ⇒ end with STATUS: TAPENADE NEEDED
   no  ⇒ end with STATUS: READY TO RUN (after make exits 0)

Context (load these and ONLY these):
- <e.g. @01_PAPER_REFERENCE.md §X.Y>
- <e.g. @04_ARCHITECTURE.md §3>

Files to edit:
- <path:line> — <what changes>
- <path>     — <what changes>

Change:
<specific instructions, code sketches if helpful, line numbers if known.
Be tight. Avoid open-ended directives. If the task can't be specified
tightly enough to fit here, it's two tasks, not one.>

Compile: cd build && make -j

Done when:
- make exits 0
- <one or two specific, mechanical checks>

Status row to flip in 02_IMPLEMENTATION_STATUS.md:
- <e.g. add a new row "T-X — <title>" set to ✅>

Do NOT:
- Modify SA model files (except γ·P_SA coupling point if already there)
- Touch files outside the list above
- Verify physics correctness (user does that)
- Create ASCII debug files
```

---

## Filled Example — Hypothetical Post-Validation Bug

```
@CLAUDE.md
@02_IMPLEMENTATION_STATUS.md

## Task: F1 — Cap rTurb in F_onset to prevent overflow at high Re

Goal:
Validation on S809 at Re=2e6 produces R_T values up to 1e8 in the wake,
which makes (R_T/2.5)^3 in the LM2015 F_onset3 overflow to Inf. Cap R_T
at 1e6 inside F_onset only — the destruction term doesn't need it.

AD impact: yes

Context:
- @01_PAPER_REFERENCE.md §4.3 (F_onset definition)

Files to edit:
- src/turbulence/saGammaRetheta.F90 (the LM2015 F_onset branch — locate
  by grep "F_onset3")

Change:
Replace
    rTurbForOnset = rTurb
with
    rTurbForOnset = min(rTurb, 1.0e6_realType)
and use rTurbForOnset only inside F_onset3. Leave R_T untouched
elsewhere (F_turb, debug output, etc.).

Compile: cd build && make -j

Done when:
- make exits 0
- grep finds exactly one usage of `rTurbForOnset` in the F_onset block

Status row:
- Add row "F1 — Cap rTurb in F_onset" set to ✅ in
  02_IMPLEMENTATION_STATUS.md, in a new "Post-validation fixes" section.

Do NOT:
- Cap R_T in F_turb
- Cap R_T in the debug output (slot 3 must show the true value)
- Touch any source term other than F_onset3
```

---

## Session End Checklist

- [ ] `make` exits 0.
- [ ] All "Done when" bullets satisfied.
- [ ] `02_IMPLEMENTATION_STATUS.md` updated.
- [ ] Committed with `T-ID:` prefix.
- [ ] Output ends with one of `STATUS: TAPENADE NEEDED` /
      `STATUS: READY TO RUN` / `STATUS: BLOCKED — <reason>`.
- [ ] User runs `/clear` before next task.

## Triage: Is It One Task or Two?

If you can't write the "Change" section in under ~20 lines without losing
detail, split. Two tight task blocks beat one fuzzy one. Common splits:

- "Add option + use option" → A: add option with default; B: wire it in.
- "Fix bug + regen AD" → A: fix bug (ends TAPENADE NEEDED); user regens
  manually; B: fix any glue breakage post-regen.
- "Implement function + call function" → A: implement and unit-buildable;
  B: call from the residual.
