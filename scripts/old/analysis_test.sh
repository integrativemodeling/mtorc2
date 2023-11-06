#!/bin/bash
module load CBI conda-stage
conda activate imp
module load Sali python3/pyrmsd/4.1.gita558b8a
module load cuda/7

EXP_DIR="exp_9/71"
OUTPUT_NUMS=517999
N_STRUCTURES=100
THRESH=-1
CLUSTER=XL
MODE=gpu
JOB_ID="test"

echo EXP_DIR: "$EXP_DIR"
echo OUTPUT_NUMS: "$OUTPUT_NUMS"
echo N_STRUCTURES: "$N_STRUCTURES"
echo THRESH: "$THRESH"
echo MODE: "$MODE"

python ~/mtorc2/scripts/analysis.py "$EXP_DIR" "$OUTPUT_NUMS" "$N_STRUCTURES" "$MODE" "$JOB_ID" "$THRESH" "$CLUSTER"