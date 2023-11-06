#!/bin/bash


cd "./test"
LOC_OUT_DIR="./output_0"

EM_COMPS=
RES_PER_COMP=10
FLEX=0
XLS=dss:1:DSS,edc:1:EDC
MAP=
EM_W=1
SHUFFLE=1
N_FRAMES=100

# Copy the modeling script into the job scratch dir.
cp $HOME/mtorc2/scripts/sample/sample_122.py .

# Copy all of the subunit gmms (.mrc).
cp "$HOME"/mtorc2/data/em/comps/122/122_$RES_PER_COMP/* .

python3 sample_122.py "$LOC_OUT_DIR" "$EM_COMPS" "$RES_PER_COMP"  "$FLEX" "$XLS" "$MAP" "$EM_W" "$SHUFFLE" "$N_FRAMES"

cd ..