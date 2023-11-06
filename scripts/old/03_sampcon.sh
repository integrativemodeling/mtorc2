#!/bin/bash

C=0
ANALYSIS_DIR="$HOME"/mtorc2/exp_8/66/analysis/90004
CLUSTER_DIR="$ANALYSIS_DIR"/cluster_"$C"

A_rmf=$(ls "$CLUSTER_DIR"/A_*.rmf3)
A_txt=$(ls "$CLUSTER_DIR"/A_*.txt)

B_rmf=$(ls "$CLUSTER_DIR"/B_*.rmf3)
B_txt=$(ls "$CLUSTER_DIR"/B_*.txt)

density=$HOME/mtorc2/scripts/density.txt

start=`date +%s`
python3 "$HOME"/imp-sampcon/pyext/src/exhaust.py \
--sysname mtorc2 \
--path "$ANALYSIS_DIR" \
--rmfA "$A_rmf" \
--rmfB "$B_rmf" \
--scoreA "$A_txt" \
--scoreB "$B_txt" \
--density "$density" \
--gnuplot \
--align \
-m cuda -c 4 \
-g 2.0 -gp
end=`date +%s`
runtime=$((end-start))
echo "$runtime"