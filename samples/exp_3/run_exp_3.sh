#!/bin/bash


conda activate imp
for NUM in {15..24}
do
  JOB_DIR=$HOME/mtorc2/single_traj_experiments/exp_3/$NUM/sample
  cd $JOB_DIR
  nohup python modeling.py & >nohup.out 2>&1
done
