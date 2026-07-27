---
name: wrapup
description: End-of-task wrap-up check. Invoke when you (Claude) judge the current task is genuinely finished — propose any needed CLAUDE.md/docs updates, then ask the user whether to commit + push. Not for mid-task pauses.
---

# Task wrap-up

Formerly a Stop hook (`propose_wrapup.sh`) that fired on every stop and
over-triggered on routine pauses. Now a skill: **you** decide when a task is
actually done and invoke this explicitly, instead of it firing automatically
on every turn boundary.

## When to invoke

Only when you believe the task block is fully complete (per CLAUDE.md's
"Definition of done" for the branch) — not after every tool call, not on a
mid-task check-in, not just because you're about to stop talking.

## Steps

1. Check for uncommitted, task-relevant changes:
   ```bash
   git status --porcelain -- src docs tests CLAUDE.md
   ```
   If clean, there's nothing to wrap up — say so briefly and stop.

2. If there are changes, decide whether `CLAUDE.md` or the docs (e.g.
   `docs/README.md`, task logs) need updating to reflect what changed. If so,
   propose the specific edits and apply them.

3. Ask the user whether to commit (and push). Only commit/push after they say
   yes. If they already declined for this same change-set earlier in the
   conversation, don't re-ask — respect that and finish.

## Notes

- This replaces the old always-on Stop hook — don't re-add a Stop hook for
  this purpose without the user asking for it back.
- Keep it to one pass: propose, ask, act on the answer. Don't loop back into
  this skill repeatedly for the same change-set.
