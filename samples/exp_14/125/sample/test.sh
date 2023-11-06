#!/bin/bash


cd "./test"
LOC_OUT_DIR="./output_0"

EM_COMPS=CRIM,KINASEN,KINASEC
RES_PER_COMP=10
FLEX=1
XLS=dss:1:DSS,edc:1:EDC
MAP=J252/J252_4A/dsfact2_cutoff.05_ng100.txt
EM_W=1
SHUFFLE=0
N_FRAMES=50

# Copy the modeling script into the job scratch dir.
cp $HOME/mtorc2/scripts/sample/sample_122.py .

# Copy all of the subunit gmms (.mrc).
cp "$HOME"/mtorc2/data/em/comps/122/122_$RES_PER_COMP/* .

#python3 sample_122.py "$LOC_OUT_DIR" "$EM_COMPS" "$RES_PER_COMP"  "$FLEX" "$XLS" "$MAP" "$EM_W" "$SHUFFLE" "$N_FRAMES"
mpirun -np 8 python3 sample_122.py "$LOC_OUT_DIR" "$EM_COMPS" "$RES_PER_COMP"  "$FLEX" "$XLS" "$MAP" "$EM_W" "$SHUFFLE" "$N_FRAMES"

cd ..