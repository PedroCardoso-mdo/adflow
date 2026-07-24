# Finished Task Log — SA-BCM

One file per finished task (a "case"), added here over time — not appended into one
ever-growing file. When a task in `../current-task.md` finishes: write
`<YYYY-MM-DD>-<slug>.md` here using the template below, add a row to the index, then reset
`../current-task.md`.

## Index

| Date | Task | File |
|---|---|---|
| 2026-07-24 | Paper verification, Tapenade audit, NK line-search relax option | [2026-07-24-paper-verification-and-nk-relax.md](2026-07-24-paper-verification-and-nk-relax.md) |
| 2026-07-24 | `reg_bcm.py` missing ANK/NK options (useNKSolver silently False); isolated-adjoint drag fix | [2026-07-24-reg-bcm-missing-ank-nk-options.md](2026-07-24-reg-bcm-missing-ank-nk-options.md) |

## Template

```
# <task name> — <YYYY-MM-DD>

**Problem:** what was broken/missing/in question, 1-3 sentences.

**Solution:** what was actually changed (files:lines), 1-5 bullets.

**Files created/touched and why:**
- `path/to/file` — purpose in this task, and why it's structured the way it is (strategy),
  not just what it contains.

**Verification:** what was run to confirm done-ness, per CLAUDE.md's definition of done.

**Follow-ups:** anything punted to `../TODO.md`, with the item name.
```
