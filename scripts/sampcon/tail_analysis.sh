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


module load Sali
module load imp
module load python3/pandas/0.25.3
module load python3/hdbscan/0.8.33

SAMPLE_NUM=4291228
ANALYSIS_NUM=1

ANALYS_DIR="/wynton/group/sali/mhancock/mtorc2/dev/09_conn_akt_tail_2/analysis/$ANALYSIS_NUM"
mkdir -p "$ANALYS_DIR/traj"

SAMPLE_DIR="/wynton/group/sali/mhancock/mtorc2/dev/09_conn_akt_tail_2/sample/$SAMPLE_NUM"

## Run analysis.
CLUSTER_ON=DR_sum,EV_sum
EQUIL=.5
MIN_CLUST=16000

python3 ~/mtorc2/scripts/sampcon/analysis.py --analysis_dir "$ANALYS_DIR" --sample_dir "$SAMPLE_DIR" --cluster_on "$CLUSTER_ON" --equil "$EQUIL" --min_clust "$MIN_CLUST"
