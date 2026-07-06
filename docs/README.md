# SA-γ-Re̅θt Transition Model — Documentation

Project-specific documentation for the Piotrowski & Zingg (2020) γ-Re̅θt
transition model implementation in ADflow. Everything here is unique to this
branch; upstream ADflow docs live under `doc/`.

| Document | What's in it |
|----------|--------------|
| [paper-reference.md](paper-reference.md) | All equations, constants, algorithms, and conventions from the paper. **Source of truth for physics** — when code and paper disagree, paper wins. |
| [architecture.md](architecture.md) | Solver architecture, state-vector layout, key code/module locations, user constraints, and the complete reference for every runtime option added for the transition model. |
| [nondimensionalization.md](nondimensionalization.md) | How ADflow makes the governing equations dimensionless — the **pressure–density (p-ρ) scaling** (velocity normalizes to M·√γ, *not* 1), reference scales, and what to watch when porting paper equations. |
| [adjoint-trace.md](adjoint-trace.md) | Paired inventory of adjoint/AD touchpoints (SA vs SA-GR): preprocessor guards, Tapenade directives, generated AD files, and wiring. |

Project rules and the task roadmap live in [`../CLAUDE.md`](../CLAUDE.md).
