#!/bin/bash
#$ -S /bin/bash
#$ -o /wynton/home/sali/mhancock/mtorc2/dev/09_conn_akt_tail_2/analysis
#$ -N tail_analysis
#$ -j y
#$ -l h_rt=2:00:00
#$ -l mem_free=5G
#$ -l scratch=50G
#$ -pe smp 25
#$ -l hostname='qb3-id*'

#
#module load CBI conda-stage
#conda activate imp_217

SAMPLE_NUM=935122
ANALYSIS_NUM=2

ANALYS_DIR="/wynton/group/sali/mhancock/mtorc2/dev/09_conn_akt_tail_2/analysis/$ANALYSIS_NUM"
mkdir -p "$ANALYS_DIR/traj"

SAMPLE_DIR="/wynton/group/sali/mhancock/mtorc2/dev/09_conn_akt_tail_2/sample/output.$SAMPLE_NUM"

## Run analysis.
CLUSTER_ON=DR_sum,EV_sum
#CLUSTER_ON=XLs_sum,EM3D_BayesianEM,DR_sum,EV_sum
EQUIL=
MIN_CLUST=25000

python ~/mtorc2/scripts/sampcon/analysis.py "$ANALYS_DIR" "$SAMPLE_DIR" "$CLUSTER_ON" "$EQUIL" "$MIN_CLUST"
