# Long-Term Memory

## Conversation History

### 2026-05-29: Initial Setup
- Created documentation structure for SA-BCM timing branch
- Identified >64 rank NaN issue - suspect `rlv` division in sa.F90
- Key files: sa.F90 (lines 296-351), ankProfiling.F90

## Decisions Made

- Documentation stored in `docs/` directory
- Focus on token-efficient .md files for context

## Open Issues

- [ ] NaN with >64 ranks - needs debugging with guards on line 323
- [ ] Verify `rlv` initialization across all ranks

## Learnings

- Timing code (ankProfiling.F90) is safe for >64 ranks
- SA-BCM physics in saSource subroutine
- `use_ANKProfiling` option in pyADflow.py:5862
