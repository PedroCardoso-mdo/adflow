#!/usr/bin/env bash
#
# PreToolUse(Edit|Write|NotebookEdit) guard. Enforces CLAUDE.md hard rule 2
# mechanically, so it no longer needs to live as prose Claude must remember:
#
#   Rule 2  — Do NOT modify the SA model directly. Transition is a modifier;
#             src/turbulence/sa.F90 (and SA-only code) must never change.
#             -> hard DENY.
#
# (Rule 6 — the Tapenade-generated adjoint files — is enforced separately by
#  the "ask" permission globs in .claude/settings.json, so it is deliberately
#  NOT duplicated here; a second gate would just double-prompt.)
#
# Reads the PreToolUse JSON on stdin; emits a permissionDecision. A path that
# matches nothing falls through to exit 0 (no opinion -> normal permissioning).
set -uo pipefail

REPO="${CLAUDE_PROJECT_DIR:-/home/mdo/MDOLab_3_v2/adflow_sa_gamma_rethetha_paper_solver}"

input=$(cat)
fp=$(printf '%s' "$input" | python3 -c 'import sys,json
try: print(json.load(sys.stdin).get("tool_input",{}).get("file_path",""))
except Exception: print("")' 2>/dev/null)

[ -z "$fp" ] && exit 0

# Normalize to a repo-relative-ish tail for matching (works for absolute or
# relative file_path).
emit() {  # $1 = allow|deny|ask, $2 = reason
    python3 - "$1" "$2" <<'PY'
import json, sys
print(json.dumps({"hookSpecificOutput": {
    "hookEventName": "PreToolUse",
    "permissionDecision": sys.argv[1],
    "permissionDecisionReason": sys.argv[2],
}}))
PY
}

log() {  # $1 = short tag
    local LOG="$REPO/docs/HOOK_ACTIVITY_LOG.md"
    local ts; ts=$(date '+%Y-%m-%d %H:%M:%S')
    printf -- '- **%s** — guard: %s (`%s`)\n' "$ts" "$1" "$fp" >> "$LOG" 2>/dev/null
}

# --- Rule 2: the SA model itself is off-limits -------------------------------
case "$fp" in
    */src/turbulence/sa.F90|src/turbulence/sa.F90)
        log "DENIED edit to SA model (rule 2)"
        emit deny "CLAUDE.md rule 2: the SA model (src/turbulence/sa.F90) must never be modified. Transition is a modifier — put the change in src/turbulence/saGammaRetheta.F90 instead."
        exit 0
        ;;
esac

exit 0
