#!/bin/bash


for i in {0..7}
do
  cd "$HOME/mtorc2/dev/06_em_arch/J252/$i"
  qsub sample.sh
done