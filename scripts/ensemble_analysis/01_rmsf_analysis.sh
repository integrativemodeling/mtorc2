#!/bin/bash
#$ -S /bin/bash
#$ -o /wynton/home/sali/mhancock/mtorc2/scripts
#$ -N rmsf
#$ -j y
#$ -l h_rt=8:00:00
#$ -l mem_free=5G
#$ -l scratch=50G
#$ -pe smp 25
#$ -l hostname='qb3-id*'

# module load CBI conda-stage
# conda activate imp_217


python /wynton/home/sali/mhancock/mtorc2/scripts/ensemble_analysis/rmsf_analysis.py