#!/usr/bin/env bash
#
# PreToolUse(Bash) hook. Whenever Claude is about to launch a long-running
# job — a SIMULATION (mpirun python …), a TAPENADE regeneration, or a TEST
# run — announce it and hand the user a ready-to-paste `tail -f` so they can
# follow the run in another terminal.
#
# Identification is by command pattern (see classify() below). The follow
# command needs a log FILE, so the hook extracts a `tee`/`>`/`>>` target from
# the command. Your run scripts print to stdout by default, so if there is no
# redirect the hook tells you how to add one instead of inventing a path.
#
# Fires BEFORE the tool runs (PreToolUse) so you can start tailing as the run
# begins. Never blocks — it only emits an informational systemMessage.
set -uo pipefail

REPO="${CLAUDE_PROJECT_DIR:-/home/mdo/MDOLab_3_v2/adflow_sa_gamma_rethetha_paper_solver}"

input=$(cat)
cmd=$(printf '%s' "$input" | python3 -c 'import sys,json
try: print(json.load(sys.stdin).get("tool_input",{}).get("command",""))
except Exception: print("")' 2>/dev/null)
[ -z "$cmd" ] && exit 0

# A git command can carry any of these words in a commit message / path — never
# a real sim/test/tapenade run. Bail before classifying to avoid false hits
# like: git commit -m "add tapenade note".
printf '%s' "$cmd" | grep -qE '(^|[;&|][[:space:]]*)git([[:space:]]|$)' && exit 0

# --- classify: tapenade -> test -> simulation (first match wins) -------------
# Specific script/makefile names match anywhere (they don't occur in prose);
# bare `tapenade`/`pytest` must be in COMMAND position (start or after ; & | &&)
# so they aren't picked up from inside a quoted argument.
kind=""
if printf '%s' "$cmd" | grep -qE 'Makefile_tapenade|AD_I\.sh|build_tapenade\.sh|(^|[;&|][[:space:]]*)tapenade([[:space:]]|$)'; then
    kind="tapenade regeneration"
elif printf '%s' "$cmd" | grep -qE 'testflo|run_sagr_tests\.sh|(^|[;&|][[:space:]]*)pytest([[:space:]]|$)'; then
    kind="test run"
elif printf '%s' "$cmd" | grep -qE 'mpirun' && printf '%s' "$cmd" | grep -qE '\.py([[:space:]]|$)'; then
    kind="simulation"
fi
[ -z "$kind" ] && exit 0   # not a long run we care about

# --- find the log file the run writes to, if any -----------------------------
# Last `| tee [-a] FILE`, or last `> FILE` / `>> FILE`. Take the last so a
# combined `>out 2>&1 | tee log` yields the tee target.
log=$(printf '%s' "$cmd" | python3 -c '
import sys, re
c = sys.stdin.read()
targets = []
for m in re.finditer(r"tee\s+(?:-a\s+)?([^\s|&;><]+)", c):
    targets.append(m.group(1))
for m in re.finditer(r">>?\s*([^\s|&;><]+)", c):
    targets.append(m.group(1))
print(targets[-1] if targets else "")
' 2>/dev/null)

# --- build the message -------------------------------------------------------
if [ -n "$log" ]; then
    case "$log" in /*) abslog="$log" ;; *) abslog="<run cwd>/$log" ;; esac
    msg="🔭 ${kind} starting — follow it in another terminal with:  tail -f ${abslog}"
else
    msg="🔭 ${kind} starting — this command writes to stdout, not a file, so there is nothing to tail. To follow it, re-run with a log, e.g. append:  2>&1 | tee run.log   then:  tail -f run.log"
fi

python3 - "$msg" <<'PY'
import json, sys
print(json.dumps({"systemMessage": sys.argv[1]}))
PY
exit 0
