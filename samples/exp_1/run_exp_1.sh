#!/bin/bash


conda activate imp
for NUM in {0..9}
do
  JOB_DIR=$HOME/mtorc2/single_traj_experiments/$NUM/sample
  cd $JOB_DIR
  nohup python modeling.py & >nohup.out 2>&1
done
