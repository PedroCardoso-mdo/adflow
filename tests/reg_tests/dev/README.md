# `dev/`: regenerating the restart states

The transition-model derivative tests linearize about converged states of the tutorial wing that
are not part of the upstream input-file archive (see `../README_transition.md`):

| Script                       | Writes                                    |
|------------------------------|-------------------------------------------|
| `generate_sagr_restart.py`   | `input_files/mdo_tutorial_sagr_dp.cgns`   |
| `generate_sgamma_restart.py` | `input_files/mdo_tutorial_sgamma_dp.cgns` |

Both read the solver options and flow conditions from `../reg_transition.py`, so the state always
matches the tests. `restart_utils.py` holds the common code.

Each script converges the flow from freestream on the real build. This is a full solve: run it on
a cluster, then retrain the reference files.

```bash
mpirun -np 16 python generate_sagr_restart.py
mpirun -np 16 python generate_sgamma_restart.py
```
