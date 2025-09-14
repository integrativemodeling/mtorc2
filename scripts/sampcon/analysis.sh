#!/bin/bash
#$ -S /bin/bash
#$ -o /wynton/group/sali/mhancock/mtorc2/analysis/tmp/$JOB_ID.o
#$ -N anlys_136
#$ -j y
#$ -l h_rt=24:00:00
#$ -l mem_free=5G
#$ -l scratch=50G
#$ -pe smp 25
#$ -l hostname='qb3-id*'


# module load Sali
# module load imp
# module load python3/pandas/0.25.3
# module load python3/hdbscan/0.8.33

EXP_NUM=15
TRAJ_NUM=141
SAMPLE_NUM=4410639
ANALYSIS_NUM=0

ANALYS_DIR="/wynton/group/sali/mhancock/mtorc2/analysis/$TRAJ_NUM/$ANALYSIS_NUM"
# rm -r "$ANALYS_DIR/traj"
mkdir -p "$ANALYS_DIR/traj"

SAMPLE_DIR="/wynton/group/sali/mhancock/mtorc2/samples/exp_$EXP_NUM/$TRAJ_NUM/$SAMPLE_NUM/"

## Run analysis.
# CLUSTER_ON=XLs_sum,EM3D_BayesianEM,DR_sum,EV_sum
CLUSTER_ON=DR_sum,EV_sum
EQUIL=.5
MIN_CLUST=8000

python3 ~/mtorc2/scripts/sampcon/analysis.py --analysis_dir "$ANALYS_DIR" --sample_dir "$SAMPLE_DIR" --cluster_on "$CLUSTER_ON" --equil "$EQUIL" --min_clust "$MIN_CLUST"
