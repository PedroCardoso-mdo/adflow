#!/usr/bin/env bash
#
# PostToolUse(Bash) hook. After any successful `make`, reinstall adflow into
# the mach env so the reg tests — which import adflow from site-packages, NOT
# from ./adflow — never run a stale binary. This is the exact failure that
# burned the gammaForSA clamp task: fix compiled into ./adflow but never
# pip-installed, so every test run exercised the old lib.
#
# Reads the PostToolUse JSON on stdin; acts only when the command ran `make`.
set -uo pipefail

REPO="${CLAUDE_PROJECT_DIR:-/home/mdo/MDOLab_3_v2/adflow_sa_gamma_rethetha_paper_solver}"
PIP="${ADFLOW_MACH_PIP:-/home/mdo/packages_v2/mach/bin/pip}"

input=$(cat)
cmd=$(printf '%s' "$input" | python3 -c 'import sys,json
try: print(json.load(sys.stdin).get("tool_input",{}).get("command",""))
except Exception: print("")' 2>/dev/null)

# Only when `make` is in COMMAND position: at the start or right after a shell
# separator (; & | &&). Avoids false positives like: git commit -m "make it".
printf '%s' "$cmd" | grep -qE '(^|[;&|])[[:space:]]*make([[:space:]]|$)' || exit 0

# Which build was it? A complex-step build drives Makefile_CS (produces
# libadflow_cs.so); anything else is the real build. Same `pip install .`
# reinstalls whichever .so was produced — this only labels the log/message so
# HOOK_ACTIVITY_LOG.md tells you which binary just went into mach.
if printf '%s' "$cmd" | grep -qE 'Makefile_CS'; then
    kind='complex-step (CS)'
else
    kind='real'
fi

cd "$REPO" 2>/dev/null || exit 0
LOG="$REPO/docs/HOOK_ACTIVITY_LOG.md"
ts=$(date '+%Y-%m-%d %H:%M:%S')
short=$(printf '%s' "$cmd" | tr '\n' ' ' | cut -c1-100)
if "$PIP" install . --no-deps -q >/tmp/adflow_pip_install.log 2>&1; then
    printf -- '- **%s** — auto-pip: reinstalled adflow (%s build) into mach env after `%s` (site-packages now matches ./adflow)\n' "$ts" "$kind" "$short" >> "$LOG" 2>/dev/null
    printf '{"systemMessage":"✔ auto-installed adflow (%s build) into mach env (post-make) — site-packages now matches ./adflow (logged to docs/HOOK_ACTIVITY_LOG.md)"}\n' "$kind"
else
    printf -- '- **%s** — auto-pip: FAILED after `%s` (see /tmp/adflow_pip_install.log) — site-packages may be STALE\n' "$ts" "$short" >> "$LOG" 2>/dev/null
    echo '{"systemMessage":"⚠ auto pip install after make FAILED (see /tmp/adflow_pip_install.log) — reg tests may run a STALE binary until you reinstall"}'
fi
exit 0
