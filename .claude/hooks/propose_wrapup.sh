#!/usr/bin/env bash
#
# Stop hook. When Claude is about to finish and there are uncommitted,
# task-relevant changes (src/, docs/, tests/, CLAUDE.md), block the stop and
# instruct Claude to (1) propose any needed CLAUDE.md / docs updates, then
# (2) ask the user whether to commit + push.
#
# Loop-safe: it records a hash of the current change-set. Once it has prompted
# for a given change-set it stays quiet for that exact set (so "no, don't
# commit" is respected and never re-nags). A new/different change-set re-arms
# it. A clean tree never fires.
set -uo pipefail

REPO="${CLAUDE_PROJECT_DIR:-/home/mdo/MDOLab_3_v2/adflow_sa_gamma_rethetha_paper_solver}"
cd "$REPO" 2>/dev/null || exit 0
git rev-parse --git-dir >/dev/null 2>&1 || exit 0

changes=$(git status --porcelain -- src docs tests CLAUDE.md 2>/dev/null)
[ -z "$changes" ] && exit 0   # clean on task-relevant paths -> allow stop

hash=$(printf '%s' "$changes" | sha1sum | cut -d' ' -f1)
sentinel="$(git rev-parse --git-dir)/claude_wrapup_hash"
if [ -f "$sentinel" ] && [ "$(cat "$sentinel" 2>/dev/null)" = "$hash" ]; then
    exit 0   # already prompted for this exact change-set -> allow stop
fi
printf '%s' "$hash" > "$sentinel"

# Log the wrap-up prompt (change-set snapshot) to the docs activity log.
LOG="$REPO/docs/HOOK_ACTIVITY_LOG.md"
ts=$(date '+%Y-%m-%d %H:%M:%S')
{
    printf -- '- **%s** — wrap-up prompt fired (change-set %s); uncommitted task-relevant paths:\n' "$ts" "${hash:0:12}"
    printf '%s\n' "$changes" | sed 's/^/    /'
} >> "$LOG" 2>/dev/null

cat <<'JSON'
{"decision":"block","reason":"Wrap-up check: there are uncommitted changes under src/, docs/, tests/, or CLAUDE.md. Before finishing this task: (1) decide whether CLAUDE.md and the docs (e.g. docs/current-task.md, docs/task-log/, docs/README.md) need updating to reflect what changed, and if so propose the specific edits and apply them; (2) then ask me whether to commit and push, and only do so after I say yes. If I have already declined committing these same changes, respect that and just finish."}
JSON
exit 0
