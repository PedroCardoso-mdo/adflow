#!/usr/bin/env python
"""Compare the two-iteration smoke runs with the cluster runs of the same cases.

    python smoke_compare.py <smoke dir>        -> <smoke dir>/SMOKE_REPORT.md

smoke_reference.json holds, for every smoke case, the first rows of the iteration table of the
cluster log (job, log path, ADflow revision). Iteration 0 is the residual of the free-stream state:
same mesh + same options + same free stream must give the same number to round-off, whatever the
rank count. Iterations 1 and 2 depend on the partition (the preconditioner is per rank and the
linear solves are inexact), so they agree to a few digits only.
"""
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TOL0 = 1e-8          # required relative agreement at iteration 0
COLS = ["totalRes", "Res rho", "C_lift", "C_drag"]


def parseLog(path):
    """Banner, cell count and every row of the ADflow iteration table(s) in a log."""
    info = dict(git=None, cells=None, columns=None, rows=[])
    with open(path, errors="replace") as f:
        for line in f:
            if "| ADflow" in line and "git" in line and info["git"] is None:
                m = re.search(r"git (\w+) \(([^)]+)\)", line)
                info["git"] = m.group(0) + (" DIRTY" if "DIRTY" in line else "")
            elif "Total number of cells" in line:
                info["cells"] = int(line.split(":")[-1])
            elif re.match(r"^#\s+Grid\s+\|\s+Iter", line):
                info["columns"] = [c.strip() for c in line.strip("# \n").split("|") if c.strip()][7:]
            elif re.match(r"^\s+1\s+\d+\s+\d+\s+\S+\s+\S+\s+\S+\s+\S+\s+[-+0-9.E]+", line) and info["columns"]:
                t = line.split()
                try:
                    vals = {c: float(v) for c, v in zip(info["columns"], t[7:])}
                except ValueError:
                    continue
                info["rows"].append(dict(iter=int(t[1]), iterTot=int(t[2]), type=t[3], **vals))
    return info


def rel(a, b):
    return abs(a - b) / abs(b) if b else abs(a - b)


def fmt(x):
    return "—" if x is None else f"{x:.1e}"


def table(header, rows):
    """Markdown table with padded columns (readable as plain text)."""
    w = [max(len(str(r[i])) for r in [header] + rows) for i in range(len(header))]
    line = lambda r: "| " + " | ".join(str(c).ljust(w[i]) for i, c in enumerate(r)) + " |"  # noqa: E731
    return "\n".join([line(header), "| " + " | ".join("-" * n for n in w) + " |"] + [line(r) for r in rows])


def main(smokeDir):
    ref = json.load(open(os.path.join(HERE, "smoke_reference.json")))
    status, agree, detail, problems = [], [], [], []
    localGit = None
    for name, entry in ref.items():
        d = os.path.join(smokeDir, name)
        cl = entry.get("cluster")
        if entry.get("skipped"):
            status.append([name, f"{cl['cells']:,}" if cl else "—", "—", "skipped", "—", entry["skipped"]])
            continue
        if not os.path.exists(os.path.join(d, "run.log")):
            status.append([name, "—", "—", "not run", "—", ""])
            continue
        loc = parseLog(os.path.join(d, "run.log"))
        localGit = localGit or loc["git"]
        rc = open(os.path.join(d, "exit_code")).read().strip() if os.path.exists(os.path.join(d, "exit_code")) else "?"
        res = {}
        if os.path.exists(os.path.join(d, "result.json")):
            res = json.load(open(os.path.join(d, "result.json")))
        nProcs = res.get("provenance", {}).get("nProcs", "?")
        nIter = max([r["iter"] for r in loc["rows"]], default=-1)
        ok = rc == "0" and bool(res) and nIter >= 1
        note = entry.get("note", "")
        status.append([name, f"{loc['cells']:,}" if loc["cells"] else "—", nProcs, "yes" if ok else f"FAILED (rc {rc})",
                       nIter, note])
        if not ok:
            problems.append(f"- `{name}`: did not run cleanly (exit code {rc}, result.json {'present' if res else 'missing'}).")
        if cl is None:
            agree.append([name, "—", "—", "—", "—", "—", "—", "no cluster twin"])
            continue
        if cl["cells"] and loc["cells"] and cl["cells"] != loc["cells"]:
            problems.append(f"- `{name}`: mesh differs (local {loc['cells']} cells, cluster {cl['cells']}).")
        # pick, for each iteration, the local row closest to the cluster row (logs with several
        # tables: two alphas, or the interleaved tables of a multipoint run)
        diffs, nCand = {}, 1
        for cr in cl["rows"][:3]:
            cand = [r for r in loc["rows"] if r["iter"] == cr["iter"] and "totalRes" in r]
            if not cand:
                continue
            nCand = max(nCand, len(cand))
            lr = min(cand, key=lambda r: rel(r["totalRes"], cr["totalRes"]))
            diffs[cr["iter"]] = {c: (lr[c], cr[c], rel(lr[c], cr[c])) for c in COLS if c in lr and c in cr}
        d0 = diffs.get(0, {})
        worst0 = max([v[2] for v in d0.values()], default=None)
        verdict = "—" if worst0 is None else ("OK" if worst0 <= TOL0 else "DIFFERS")
        if verdict == "DIFFERS":
            problems.append(f"- `{name}`: iteration 0 differs from the cluster by {worst0:.1e} (see the detail table).")
        g = lambda i, c: fmt(diffs.get(i, {}).get(c, (None, None, None))[2])  # noqa: E731
        agree.append([name, g(0, "totalRes"), g(0, "Res rho"), g(0, "C_lift") + " / " + g(0, "C_drag"),
                      g(1, "totalRes"), g(2, "totalRes"), verdict,
                      ("job " if cl["job"][:1].isdigit() else "") + f"{cl['job']}, {cl['git'] or 'no banner in the log'}" + (f"; closest of {nCand} local tables" if nCand > 1 else "")])
        if 0 in diffs and "totalRes" in d0:
            detail.append([name, f"{d0['totalRes'][0]:.16E}", f"{d0['totalRes'][1]:.16E}", fmt(d0["totalRes"][2]),
                           "`" + cl["log"].split("TransitionModel/")[-1] + "`"])

    out = ["# Smoke report: two-iteration runs against the cluster runs", "",
           f"Local build: `{localGit}`; every case from free stream, `--nCycles 4` (5 for the SA-BCM wing), "
           "optimisation scripts with `--baselineOnly`.", "",
           "Generated by `examples/transition/smoke_compare.py` from `<case>/run.log` and "
           "`examples/transition/smoke_reference.json`.", "",
           "## 1. Did it run", "",
           table(["case", "cells", "ranks", "ran", "last iteration", "note"], status), "",
           "## 2. Agreement with the cluster (relative differences)", "",
           "Iteration 0 = residual of the free-stream state: must agree to round-off (limit "
           f"{TOL0:g}). Iterations 1-2 depend on the rank count (per-rank preconditioner, inexact linear solves).", "",
           table(["case", "it 0 totalRes", "it 0 Res rho", "it 0 cl / cd", "it 1 totalRes", "it 2 totalRes",
                  "it 0 verdict", "cluster run"], agree), "",
           "## 3. Iteration-0 residual, full precision", "",
           table(["case", "local totalRes", "cluster totalRes", "rel. diff", "cluster log (under TransitionModel/)"], detail), "",
           "## 4. Problems", "", "\n".join(problems) if problems else "None.", ""]
    notes = os.path.join(smokeDir, "SMOKE_NOTES.md")       # hand-written remarks on this particular run
    if os.path.exists(notes):
        out += [open(notes).read()]
    path = os.path.join(smokeDir, "SMOKE_REPORT.md")
    with open(path, "w") as f:
        f.write("\n".join(out))
    print(f"wrote {path}: {sum(1 for r in status if r[3] == 'yes')} ran, "
          f"{sum(1 for r in agree if r[6] == 'OK')} agree at iteration 0, {len(problems)} problem(s)")


if __name__ == "__main__":
    main(sys.argv[1])
