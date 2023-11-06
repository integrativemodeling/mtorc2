#!/bin/bash
for i in {10..33}
do
  cd "$i"
  qsub sample.sh
  cd ..
done