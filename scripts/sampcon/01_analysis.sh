#!/bin/bash
#$ -S /bin/bash
#$ -o /wynton/home/sali/mhancock/mtorc2/samples/exp_14/129/analysis
#$ -N anlys_129
#$ -j y
#$ -l h_rt=4:00:00
#$ -l mem_free=5G
#$ -l scratch=50G
#$ -pe smp 25
#$ -l hostname='qb3-id*'


#module load CBI conda-stage
#conda activate imp_217

EXP_NUM=14
TRAJ_NUM=129
SAMPLE_NUM=347466
ANALYSIS_NUM=1

ANALYS_DIR="/wynton/group/sali/mhancock/mtorc2/samples/exp_$EXP_NUM/$TRAJ_NUM/analysis/$ANALYSIS_NUM"
mkdir -p "$ANALYS_DIR/traj"

SAMPLE_DIR="/wynton/group/sali/mhancock/mtorc2/samples/exp_$EXP_NUM/$TRAJ_NUM/sample/output.$SAMPLE_NUM"

## Run analysis.
CLUSTER_ON=XLs_sum,EV_sum
#CLUSTER_ON=XLs_sum,EM3D_BayesianEM
#CLUSTER_ON=XLs_sum,EM3D_BayesianEM,DR_sum,EV_sum
EQUIL=
MIN_CLUST=50000

python ~/mtorc2/scripts/sampcon/analysis.py "$ANALYS_DIR" "$SAMPLE_DIR" "$CLUSTER_ON" "$EQUIL" "$MIN_CLUST"
