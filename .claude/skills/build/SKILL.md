---
name: build
description: Compile and install ADflow into the mach env, in real and/or complex-step mode. Handles activating the mach venv (needed for the build's Python deps) and selecting the correct PETSc arch (real-debug vs complex-debug). Invoke for "build adflow", "compile and install", "rebuild real/complex/both", or /build.
---

# Build & install ADflow (real / complex)

Compile ADflow and install it into the **mach** env so `site-packages` matches
`./adflow` (the reg tests import from site-packages, not the repo). Two build
modes, runnable individually or together.

## Argument

One of: `real` | `complex` | `both`.
- If the user named a mode, use it.
- If **no mode** was given, ask which (real / complex / both) before building —
  a full build is slow, don't guess.

## Why the setup steps matter

- **Activate the mach venv first.** The Fortran build calls Python (f2py,
  numpy, etc.) from the mach env. Without activation the build uses the wrong
  interpreter and fails or links against the wrong deps. Activation also makes
  `pip`/`python` the mach ones, so `pip install .` lands in the right place.
- **Select the PETSc arch.** Default is `real-debug`. The **complex** build
  must switch to `complex-debug` (real and complex link different PETSc libs).
- The complex build lives in a separate `src_cs/` tree, so **real and complex
  do not share objects** — building `both` needs no clean in between.

## Steps

Run each mode as ONE compound command so it works regardless of prior shell
state (the Bash tool does not persist env between calls). Builds are slow:
give the make a long timeout (up to the 600000 ms max), or run it in the
**background** and monitor if a from-scratch build might exceed 10 min. An
incremental rebuild after editing a few files is usually well under that.

### real
```bash
source /home/mdo/packages_v2/mach/bin/activate && \
export PETSC_ARCH=real-debug && \
make && \
pip install . --no-deps
```
Produces `adflow/libadflow.so`.

### complex
```bash
source /home/mdo/packages_v2/mach/bin/activate && \
export PETSC_ARCH=complex-debug && \
make -f Makefile_CS PETSC_ARCH=complex-debug && \
pip install . --no-deps
```
Produces `adflow/libadflow_cs.so`.

### both
Run **real**, then **complex** (order doesn't matter; they're independent).

## If the build errors: clean and retry once

A compile error is often stale objects (e.g. a changed module interface). On
failure, run the matching clean, then rebuild that mode **once**. If it still
fails after the clean, stop and report the error — don't loop.

- real → `make clean` then the real build again.
- complex → `make -f Makefile_CS clean` then the complex build again.

```bash
# real, clean + rebuild
source /home/mdo/packages_v2/mach/bin/activate && \
export PETSC_ARCH=real-debug && \
make clean && make && pip install . --no-deps
```
```bash
# complex, clean + rebuild
source /home/mdo/packages_v2/mach/bin/activate && \
export PETSC_ARCH=complex-debug && \
make -f Makefile_CS clean && make -f Makefile_CS PETSC_ARCH=complex-debug && \
pip install . --no-deps
```

## Verify

After each install, confirm the mach env can import the freshly installed lib:
```bash
source /home/mdo/packages_v2/mach/bin/activate && \
python -c "import adflow; print('adflow OK:', adflow.__file__)"
```
For complex, also confirm the CS lib is present:
```bash
ls -la adflow/libadflow_cs.so
```

## Notes

- **config.mk must exist** at `config/config.mk` (it does). If it's ever
  missing, the make prints a message to copy one from `config/defaults/` —
  surface that to the user rather than guessing a config.
- The `auto_pip_after_make.sh` PostToolUse hook also reinstalls after any
  `make`, so a redundant install may run — harmless (idempotent). The explicit
  `pip install .` here keeps the skill self-contained if the hook is off.
- **Definition of done** (per CLAUDE.md): the build compiles and installs
  without error. Do not run physics verification here — that's a separate step.
- Report the outcome plainly: which mode(s) built, whether the install and
  import succeeded, and paste any compile error verbatim if it failed.
