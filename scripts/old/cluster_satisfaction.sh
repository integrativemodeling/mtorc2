#$ -S /bin/bash
#$ -N cluster_satisfaction
#$ -o /wynton/home/sali/mhancock/mtorc2/exp_8/66/analysis
#$ -j y
#$ -l h_rt=18:00:00
#$ -l mem_free=5G
#$ -l scratch=10G
#$ -l hostname='qb3-id*'


module load CBI conda-stage
conda activate imp

python ~/mtorc2/scripts/cluster_satisfaction.py

