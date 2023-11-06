#!/bin/bash


for EXP_DIR in rigid_body_1 exp_8 exp_9 exp_10
do
  REMOTE_PATH="/wynton/group/sali/mhancock/mtorc2/old/$EXP_DIR"
  LOCAL_PATH="/salilab/park1/matthew/mtorc2/old"
  echo "$REMOTE_PATH"

  scp -r mhancock@dt2.wynton.ucsf.edu:"$REMOTE_PATH" "$LOCAL_PATH"
done