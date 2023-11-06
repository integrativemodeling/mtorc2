#!/bin/bash


cd ./test
LOC_OUT_DIR="./"

PDB_FILE="$HOME/mtorc2/manuscript/submission_2/models/126_0_2_3/cluster.0/cluster_center_model.pdb"
XLS=dss:1:DSS,edc:1:EDC
PHOS_RES=473
SHUFFLE=1
N_FRAMES=50000

cp $HOME/mtorc2/dev/09_conn_akt_tail_2/scripts/model_akt_tail.py .

python3 model_akt_tail.py "$LOC_OUT_DIR" "$PDB_FILE" "$XLS" "$PHOS_RES" "$SHUFFLE" "$N_FRAMES"

cd ..