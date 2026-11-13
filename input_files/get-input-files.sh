#!/bin/bash
# this file will download the input files for ADflow regression tests
# and extract them to the right place.

DIR=$(dirname $0)
wget -O $DIR/input_files.tar.gz https://websites.umich.edu/~mdolaboratory/repo_files/ADflow/adflow_input_files.tar.gz
tar -xzf $DIR/input_files.tar.gz -C $DIR/../

# ----------------------------------------------------------------------------------------------
# NOTE: transition-model tests
# The restart states used by tests/reg_tests/test_*_sagr.py and test_*_sgamma.py are NOT in the
# archive above:
#   mdo_tutorial_sagr_dp.cgns
#   mdo_tutorial_sgamma_dp.cgns
# They ship with the showcase package (showcase/input_files/): copy or link them into this
# directory. They can also be regenerated with tests/reg_tests/dev/generate_sagr_restart.py and
# tests/reg_tests/dev/generate_sgamma_restart.py (full flow solves).
# See tests/reg_tests/README_transition.md.
# ----------------------------------------------------------------------------------------------
for f in mdo_tutorial_sagr_dp.cgns mdo_tutorial_sgamma_dp.cgns; do
    if [ ! -e "$DIR/$f" ] && [ -e "$DIR/../showcase/input_files/$f" ]; then
        ln -s "../showcase/input_files/$f" "$DIR/$f"
    fi
    [ -e "$DIR/$f" ] || echo "Transition-model restart file missing: input_files/$f (see the note in $0)"
done
